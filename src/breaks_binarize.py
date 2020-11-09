import time
import pandas as pd
import pyranges as pr
import argparse

class Timeit():
    def __init__(self, name):
        self.start = None
        self.name = name

    def __enter__(self):
        self.start = time.time()

    def __exit__(self, type, value, traceback):
        print("Finished {} in {:.3f}s".format(self.name, (time.time() - self.start)))


def coverage(intervals, features, feature_name, fun=sum):
    columns_return = ["gene_chrom", "bin_start", "bin_end", "gene_name", "gene_strand", "gene_start", "gene_end", "gene_region_start", "gene_region_end", "attributes", "attributes_details"]
    columns_group = ["gene_chrom", "gene_strand", "gene_name", "gene_start", "gene_end", "gene_region_start", "gene_region_end", "bin_start", "bin_end"]
    columns_preserve = list(set(columns_group + ["attributes", "attributes_details"]) & set(intervals.columns))

    with Timeit("pyranges join"):
        intervals_pr = pr.from_dict({**{'Chromosome': intervals["gene_chrom"], 'Start': intervals["bin_start"]-1, 'End': intervals["bin_end"]+1, 'Strand': intervals["gene_strand"]}, **intervals})
        features_pr = pr.from_dict({**{'Chromosome': features["feature_chrom"], 'Start': features["feature_start"], 'End': features["feature_end"], 'Strand': features["feature_strand"], 'feature_name': features["feature_name"]}, **features})

        overlaps = intervals_pr.join(features_pr, how=False, strandedness="same").as_df()
        overlaps["hit"] = True
        coverage = overlaps.groupby(columns_group, as_index=False).aggregate({'hit': fun})
        coverage_details = overlaps.\
            drop_duplicates(columns_group + ["feature_name" ,"hit"]).\
            groupby(columns_group, as_index=False).\
            agg({'feature_name': ','.join})

        results = intervals[columns_preserve].merge(coverage, how="left", on=columns_group)
        results = results.merge(coverage_details, how="left", on=columns_group)

        results["attributes_details"] = results["attributes"] + "; " + feature_name + "=" + results["feature_name"].fillna("") \
            if "attributes_details" in results \
            else feature_name + "=" + results["feature_name"].fillna("")

        results["attributes"] = results["attributes"] + "; " + feature_name + "=" + results["hit"].fillna(0).astype(str) \
            if "attributes" in results \
            else feature_name + "=" + results["hit"].fillna(0).astype(str)

        return results[columns_return]


def make_windows(annotations, window_size, window_step):
    annotations_bin_keys = {"gene_name", "gene_chrom", "gene_start", "gene_end", "gene_strand", "gene_region_end", "gene_region_start"}
    annotations_bin = {k: [] for k in annotations_bin_keys}
    annotations_bin["bin_start"] = []
    annotations_bin["bin_end"] = []
    for r, row in annotations.iterrows():
        numOfChunks = int((row["gene_region_end"] - row["gene_region_start"] - window_size) / window_step) + 1
        bins = list(range(0, numOfChunks * window_step, window_step))

        # Original strand
        annotations_bin["bin_start"].extend([int(i+row["gene_region_start"]) for i in bins])
        annotations_bin["bin_end"].extend([int(i+window_size+row["gene_region_start"]) for i in bins])
        for k in annotations_bin_keys:
            annotations_bin[k].extend([row[k]]*len(bins))

        # Reverse strand
        annotations_bin["bin_start"].extend([int(i + row["gene_region_start"]) for i in bins])
        annotations_bin["bin_end"].extend([int(i + window_size + row["gene_region_start"]) for i in bins])
        annotations_bin["gene_strand"].extend(["-" if row["gene_strand"] == "+" else "+"]*len(bins))
        for k in (annotations_bin_keys-{"gene_strand"}):
            annotations_bin[k].extend([row[k]]*len(bins))

    return pd.DataFrame.from_dict(annotations_bin)


def main():
    parser = argparse.ArgumentParser(description='Use a sliding window to aggregate breaks in bed file')
    parser.add_argument('input', help='Input .bed file with detected breaks')
    parser.add_argument('annotations', help='Annotation file. If annotation file has gtf or gff extention (possibly .gz) then only transcripts are selected. If .bed file is provided then all annotations from bed file are used')
    parser.add_argument('output', help='Aggregated and annotated file')
    parser.add_argument('-w|--window-size', dest="window_size", default=int(1e5), type=int, help='Window at which to agregate breaks number')
    parser.add_argument('-s|--window-step', dest="window_step", default=int(1e4), type=int, help='Step after each window')
    parser.add_argument('-e|--extend-gene', dest="extend", default=0, type=int, help='Extend each gene both directions by number of base pairs')
    parser.add_argument('-f|--features', dest="features", action="append", default=[], nargs="*", help='Additional features to annotate input file')
    args = parser.parse_args()

    with Timeit("total time"):
        # Read input .bed file
        input = pd.read_csv(args.input, sep="\t", header=0,
                            names=["feature_chrom", "feature_start", "feature_end", "feature_name", "feature_score",
                                   "feature_strand"], dtype={'feature_start': 'int64', 'feature_end': 'int64'})

        # Read annotation file (gff)
        with Timeit('processing GFF file 1'):
            cols = ["gene_chrom", "gene_source", "gene_feature", "gene_start", "gene_end", "gene_score", "gene_strand", "gene_frame", "gene_attribute"]
            dtypes = {'gene_start': 'int64', 'gene_end': 'int64'}
            annotations = pd.read_csv(args.annotations, sep="\t", header=0, names=cols, dtype=dtypes)
            annotations["gene_name"] = annotations["gene_attribute"].str.extract('gene_id "([^"]+)"')
            annotations["gene_length"] = annotations["gene_end"] - annotations["gene_start"]
            annotations = annotations.query("gene_feature=='transcript' & ~gene_chrom.str.contains('_', regex=False)")

            # Select longest transcript to represent gene (TODO: this can be 1s faster if done inline)
            annotations = annotations.sort_values(["gene_strand", "gene_chrom", "gene_name", "gene_length"], ascending=False)
            annotations = annotations.drop_duplicates(["gene_strand", "gene_chrom", "gene_name"], keep="first")

            # Select columns and extend region
            annotations = annotations[["gene_chrom", "gene_name", "gene_start", "gene_end", "gene_strand", "gene_length"]]
            annotations["gene_region_start"] = annotations["gene_start"] - args.extend
            annotations["gene_region_end"] = annotations["gene_end"] + args.extend
            # annotations.to_csv("test/annotations_longest.tsv", index=False, sep="\t")


        #
        # Report genes that were shorter than suggested
        #
        small_genes = annotations.query("gene_region_end - gene_region_start < {window_size}".format(window_size=args.window_size))
        if len(small_genes) > 0:
            small_genes_names = small_genes["gene_name"].values
            small_genes_len = small_genes["gene_region_end"].values - small_genes["gene_region_start"].values
            small_genes_repr = ", ".join("{}|{:0.0f}k".format(n, s / 1000) for n, s in zip(small_genes_names, small_genes_len))
            print("{n} genes were excluded because their regions were smaller than window size ({window_size})".format(window_size=args.window_size, n=len(small_genes_names)))
            # print(small_genes_repr)
        annotations = annotations.query(
            "gene_region_end - gene_region_start >= {window_size}".format(window_size=args.window_size))

        #
        # Split genes into intervals using STEP and WINDOW size
        #
        with Timeit('spliting gene regions to intervals'):
            annotations_bin = make_windows(annotations, args.window_size, args.window_step)

        #
        # Count features
        #
        breaks_bin_agg = coverage(annotations_bin, input, "breaks")
        # Process additional annotations
        for f in args.features:
            print(f)

        # Save results
        breaks_bin_agg.to_csv(args.output, index=False, sep="\t")


if __name__ == "__main__":
    main()
