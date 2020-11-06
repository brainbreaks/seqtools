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
        print("Finished {} in {:.1f}s".format(self.name, (time.time() - self.start)))


def merge_intervals(intervals, features):
    features["feature_chrom_idx"] = pd.Categorical(features["feature_chrom"]).codes
    intervals["gene_chrom_idx"] = pd.Categorical(intervals["gene_chrom"]).codes
    # features = features.sort_values(["feature_chrom_idx", "feature_start", "feature_end"], ignore_index=True)
    # intervals = intervals.sort_values(["gene_chrom_idx", "bin_start", "bin_end"], ignore_index=True)

    # for cat in intervals.catg.unique().tolist():
    #     features[cat] = 0
    left_indexes = []
    right_indexes = []

    # features_it = features.iterrows()
    # intervals_it = intervals.iterrows()
    # features_idx, features_row = next(features_it)
    # intervals_idx, intervals_row = next(intervals_it)

    features_it = iter(range(features.shape[0]))
    intervals_it = iter(range(intervals.shape[0]))
    features_idx = next(features_it)
    intervals_idx = next(intervals_it)
    while True:
        try:
            # r1 = Range(start=features_row.datetime_start, end=features_row.datetime_end)
            # r2 = Range(start=intervals_row.datetime_start, end=intervals_row.datetime_end)
            if intervals.at[intervals_idx, "gene_chrom_idx"] < features.at[features_idx, "feature_chrom_idx"] or intervals.at[intervals_idx, "bin_end"] < features.at[features_idx, "feature_start"]:
                # no overlap. INTERVAL before FEATURE. advance intervals_it
                intervals_idx = next(intervals_it)
            elif features.at[features_idx, "feature_chrom_idx"] < intervals.at[intervals_idx, "gene_chrom_idx"] or features.at[features_idx, "feature_end"] < intervals.at[intervals_idx, "bin_start"]:
                # no overlap. FEATURE before INTERVAL. advance features_it
                features_idx = next(features_it)
            else:
                # overlap. overlap(features_row, intervals_row) must > 0
                # features.loc[features_idx, intervals_row.catg] += overlap(features_row, intervals_row)
                # if overlap(intervals_row, features_row):
                #     left_indexes.append(intervals_idx)
                #     right_indexes.append(features_idx)

                if (features.at[features_idx, "feature_start"] >= intervals.at[intervals_idx, "bin_start"]) & (features.at[features_idx, "feature_start"] <= intervals.at[intervals_idx, "bin_end"]) | \
                    (features.at[features_idx, "feature_end"] >= intervals.at[intervals_idx, "bin_start"]) & (features.at[features_idx, "feature_end"] <= intervals.at[intervals_idx, "bin_end"]) | \
                    (intervals.at[intervals_idx, "bin_start"] >= features.at[features_idx, "feature_start"]) & (intervals.at[intervals_idx, "bin_start"] <= features.at[features_idx, "feature_end"]) | \
                    (intervals.at[intervals_idx, "bin_end"] >= features.at[features_idx, "feature_start"]) & (intervals.at[intervals_idx, "bin_end"] <= features.at[features_idx, "feature_end"]):
                    left_indexes.append(intervals_idx)
                    right_indexes.append(features_idx)

                # determine whether to advance features_it or intervals_it
                if features.at[features_idx, "feature_end"] < intervals.at[intervals_idx, "bin_end"]:
                    # advance features_it
                    features_idx = next(features_it)
                else:
                    # advance intervals_it
                    intervals_idx = next(intervals_it)
        except StopIteration:
            break

    return left_indexes, right_indexes

def map_intervals(intervals, features, feature_name, fun=sum):
    columns_return = ["gene_chrom", "bin_start", "bin_end", "gene_name", "gene_strand", "gene_start", "gene_end", "attributes", "attributes_details"]
    columns_group = ["gene_chrom", "gene_strand", "gene_name", "bin_start", "bin_end"]
    columns_preserve = list({"gene_chrom", "gene_name", "gene_start", "gene_end", "attributes", "attributes_details"} & set(intervals.columns))

    with Timeit("pyranges join"):
        intervals_pr = pr.from_dict({**{'Chromosome': intervals["gene_chrom"], 'Start': intervals["bin_start"], 'End': intervals["bin_end"], 'Strand': intervals["gene_strand"]}, **intervals})
        features_pr = pr.from_dict({**{'Chromosome': features["feature_chrom"], 'Start': features["feature_start"], 'End': features["feature_end"], 'Strand': features["feature_strand"], 'feature_name': features["feature_name"]}, **features})
        results3 = intervals_pr.join(features_pr, how="left", strandedness="same").as_df()
        results3["hit"] = results3["feature_name"].isna() | (results3["feature_name"] == "-1")

        results3_hits = results3.query("hit")[columns_preserve].drop_duplicates(["gene_chrom", "gene_name"])
        results3_agg1 = results3.groupby(columns_group, as_index=False).agg({'hit': fun})
        results3_agg2 = results3.query("hit").groupby(columns_group, as_index=False).agg({'feature_name': ','.join})
        results3_agg = results3_agg1.merge(results3_agg2, how="left", on=columns_group)
        results3_agg = results3_agg.merge(results3_hits, how="inner", on=["gene_chrom", "gene_name"])

        if "attributes" in results3_agg:
            results3_agg["attributes"] = results3_agg["attributes"] + "; " + feature_name + "=" + results3_agg["hit"].astype(str)
        else:
            results3_agg["attributes"] = feature_name + "=" + results3_agg["hit"].astype(str)

        if "attributes_details" in results3_agg:
            results3_agg["attributes_details"] = results3_agg["attributes"] + "; " + feature_name + "=" + results3_agg["feature_name"]
        else:
            results3_agg["attributes_details"] = feature_name + "=" + results3_agg["feature_name"]

        return results3_agg[columns_return]



    # with Timeit("pandas join"):
    #     features = features.sort_values(["feature_chrom", "feature_start", "feature_end"], ignore_index=True)
    #     intervals = intervals.sort_values(["gene_chrom", "bin_start", "bin_end"], ignore_index=True)
    #
    #     results = intervals.merge(features, left_on=["gene_chrom", "gene_strand"], right_on=["feature_chrom", "feature_strand"], how="left")
    #     results["hits"] = (results["feature_start"] >= results["bin_start"]) & (results["feature_start"] <= results["bin_end"]) | \
    #                          (results["feature_end"] >= results["bin_start"]) & (results["feature_end"] <= results["bin_end"]) | \
    #                          (results["bin_start"] >= results["feature_start"]) & (results["bin_start"] <= results["feature_end"]) | \
    #                          (results["bin_end"] >= results["feature_start"]) & (results["bin_end"] <= results["feature_end"])
    #
    #
    #
    #     group_keys = ["gene_chrom", "gene_strand", "bin_start", "bin_end"]
    #     results_keys = ["gene_chrom", "gene_strand", "gene_name", "bin_start", "bin_end", "gene_start", "gene_end"]
    #     if "attributes" in results:
    #         results_keys.append("attributes")
    #     if "attributes_details" in results:
    #         results_keys.append("attributes_details")
    #
    #     time_join = time.time()
    #     results_unique = results[results_keys].drop_duplicates()
    #     results_compute = results.query("hits").groupby(group_keys, as_index=False).agg({'feature_name': ','.join, "hits": fun})
    #     results_agg = results_unique.merge(results_compute, how="left", on=group_keys)
    #     results_agg.fillna({'feature_name': "", 'hits': 0}, inplace=True)
    #     print("Time to map feature {}: {:.1f} seconds".format(feature_name, (time.time() - time_join)))
    #
    #     if "attributes" in results_agg:
    #         results_agg["attributes"] = results_agg["attributes"] + "; " + feature_name + "=" + results_agg["hits"].astype(str)
    #     else:
    #         results_agg["attributes"] = feature_name + "=" + results_agg["hits"].astype(str)
    #
    #     if "attributes_details" in results_agg:
    #         results_agg["attributes_details"] = results_agg["attributes"] + "; " + feature_name + "=" + results_agg["feature_name"]
    #     else:
    #         results_agg["attributes_details"] = feature_name + "=" + results_agg["feature_name"]
    #
    #     return results_agg[["gene_chrom", "bin_start", "bin_end", "gene_name", "gene_strand", "gene_start", "gene_end",
    #                         "attributes", "attributes_details"]]


parser = argparse.ArgumentParser(description='Use a sliding window to aggregate breaks in bed file')
parser.add_argument('input', help='Input .bed file with detected breaks')
parser.add_argument('annotations', help='Annotation file. If annotation file has gtf or gff extention (possibly .gz) then only transcripts are selected. If .bed file is provided then all annotations from bed file are used')
parser.add_argument('output', help='Aggregated and annotated file')
parser.add_argument('-w|--window-size', dest="window_size", default=int(1e5), type=int, help='Window at which to agregate breaks number')
parser.add_argument('-s|--window-step', dest="window_step", default=int(1e4), type=int, help='Step after each window')
parser.add_argument('-e|--extend-gene', dest="extend", default=0, type=int, help='Extend each gene both directions by number of base pairs')
parser.add_argument('-f|--features', dest="features", action="append", default=[], nargs="*", help='Additional features to annotate input file')
args = parser.parse_args()

start_total = time.time()

# Read input .bed file
input = pd.read_csv(args.input, sep="\t", header=0, names=["feature_chrom", "feature_start", "feature_end", "feature_name", "feature_score", "feature_strand"], dtype={'feature_start': 'int64', 'feature_end': 'int64'})

# time_gff = time.time()
# annotations = pr.read_gff(args.annotations)
# annotations = annotations[(annotations.Feature=='transcript') & ~annotations.Chromosome.str.contains('_')]
# annotations_df = annotations.as_df()
# annotations_df["gene_length"] = annotations.lengths()
# annotations_df.sort_values(["gene_name", "gene_length"], ascending=False, inplace=True)
# annotations_df.drop_duplicates(["Strand", "Chromosome", "gene_name"], keep="first", inplace=True)
# print("GFF processing time: {:.1f} seconds".format((time.time() - time_gff)))


# Read annotation file (gff)
with Timeit('processing GFF file'):
    annotations = pd.read_csv(args.annotations, sep="\t", header=0, names=["gene_chrom", "gene_source", "gene_feature", "gene_start", "gene_end", "gene_score", "gene_strand", "gene_frame", "gene_attribute"], dtype={'gene_start': 'int64', 'gene_end': 'int64'})
    annotations = annotations.query("gene_feature=='transcript' & ~gene_chrom.str.contains('_')")
    annotations["gene_name"] = annotations["gene_attribute"].str.extract('gene_id "([^"]+)"')
    annotations["gene_length"] = annotations["gene_end"] - annotations["gene_start"]
    annotations["gene_attributes"] = "gene=" + annotations["gene_name"] + "; length=" + annotations["gene_length"].astype(str)
    annotations.sort_values(["gene_strand", "gene_chrom", "gene_name", "gene_length"], ascending=False, inplace=True)
    annotations.drop_duplicates(["gene_strand", "gene_chrom", "gene_name"], keep="first", inplace=True)
    annotations = annotations[["gene_chrom","gene_name","gene_start","gene_end","gene_strand","gene_length"]]
    annotations["gene_region_start"] = annotations["gene_start"] - args.extend
    annotations["gene_region_end"] = annotations["gene_end"] + args.extend
    #annotations.to_csv("test/annotations_longest.tsv", index=False, sep="\t")

#
# Report genes that were shorter than suggested
#
small_genes = annotations.query("gene_region_end - gene_region_start < {window_size}".format(window_size=args.window_size))
if len(small_genes)>0:
    small_genes_names = small_genes["gene_name"].values
    small_genes_len = small_genes["gene_region_end"].values - small_genes["gene_region_start"].values
    small_genes_repr = ", ".join("{}|{:0.0f}k".format(n, s/1000) for n,s in zip(small_genes_names, small_genes_len))
    print("Some genes were excluded because they were smaller than window size ({window_size})".format(window_size=args.window_size))
    # print(small_genes_repr)
annotations = annotations.query("gene_region_end - gene_region_start >= {window_size}".format(window_size=args.window_size))

#
# Split genes into intervals using STEP and WINDOW size
#
with Timeit('spliting gene regions to intervals'):
    annotations_bin_keys = {"gene_name", "gene_chrom", "gene_start", "gene_end", "gene_strand"}
    annotations_bin = {k: [] for k in annotations_bin_keys}
    annotations_bin["bin_start"] = []
    annotations_bin["bin_end"] = []
    for r, row in annotations.iterrows():
        numOfChunks = int((row["gene_end"] - row["gene_start"] - args.window_size) / args.window_step) + 1
        bins = list(range(0, numOfChunks * args.window_step, args.window_step))

        # Original strand
        annotations_bin["bin_start"].extend([int(i+row["gene_start"]) for i in bins])
        annotations_bin["bin_end"].extend([int(i+args.window_size+row["gene_start"]) for i in bins])
        for k in annotations_bin_keys:
            annotations_bin[k].extend([row[k]]*len(bins))

        # Reverse strand
        annotations_bin["bin_start"].extend([int(i + row["gene_start"]) for i in bins])
        annotations_bin["bin_end"].extend([int(i + args.window_size + row["gene_start"]) for i in bins])
        annotations_bin["gene_strand"].extend(["-" if row["gene_strand"] == "+" else "+"]*len(bins))
        for k in (annotations_bin_keys-{"gene_strand"}):
            annotations_bin[k].extend([row[k]]*len(bins))

    annotations_bin = pd.DataFrame.from_dict(annotations_bin)

#
# Count features
#
breaks_bin_agg = map_intervals(annotations_bin, input, "score")
# Process additional annotations
for f in args.features:
    print(f)

print("Total time: {:.1f} seconds".format((time.time() - start_total)))

# Save results
breaks_bin_agg.to_csv(args.output, index=False, sep="\t")
