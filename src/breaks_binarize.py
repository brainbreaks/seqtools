import time
import os
import re
import pandas as pd
import pyranges as pr
import argparse
import multiprocessing as mp
from glob import glob
from natsort import natsorted


class Timeit():
    def __init__(self, name, display_on_enter=False, pad=0):
        self.start = None
        self.name = name
        self.pad = pad
        self.display_on_enter = display_on_enter

    def __enter__(self):
        if self.display_on_enter:
            print(" "*self.pad + "Starting {}".format(self.name))
        self.start = time.time()

    def __exit__(self, type, value, traceback):
        if self.display_on_enter:
            print(" " * self.pad + "Finished {} in {:.3f}s".format(self.name, (time.time() - self.start)))
        else:
            print(" " * self.pad + "{}: {:.3f}s".format(self.name, (time.time() - self.start)))


def coverage(intervals, features, feature_name, fun=sum, details=True):
    columns_attributes = ["attributes"] + (["attributes_details"] if details else [])
    columns_group = ["bin_start", "bin_end", "bin_strand", "gene_chrom", "gene_name", "gene_strand", "gene_start", "gene_end", "gene_region_start", "gene_region_end"]
    columns_return = columns_group + columns_attributes
    columns_preserve = list(set(columns_return) & set(intervals.columns))

    intervals_pr = pr.from_dict({**{'Chromosome': intervals["gene_chrom"], 'Start': intervals["bin_start"]-1, 'End': intervals["bin_end"]+1, 'Strand': intervals["bin_strand"]}, **intervals})
    features_pr = pr.from_dict({**{'Chromosome': features["feature_chrom"], 'Start': features["feature_start"], 'End': features["feature_end"], 'Strand': features["feature_strand"], 'feature_name': features["feature_name"]}, **features})

    overlaps = intervals_pr.join(features_pr, how=False, strandedness="same").as_df()
    overlaps["hit"] = True
    coverage = overlaps.groupby(columns_group, as_index=False).aggregate({'hit': fun})

    results = intervals[columns_preserve].merge(coverage, how="left", on=columns_group)
    if details:
        coverage_details = overlaps.\
            drop_duplicates(columns_group + ["feature_name", "hit"]).\
            groupby(columns_group, as_index=False).\
            agg({'feature_name': ','.join})
        results = results.merge(coverage_details, how="left", on=columns_group)
        results["attributes_details"] = results["attributes"] + "; " + feature_name + "=" + results["feature_name"].fillna("") \
            if "attributes_details" in results \
            else feature_name + "=" + results["feature_name"].fillna("")

    results["attributes"] = results["attributes"] + "; " + feature_name + "=" + results["hit"].fillna(0).astype(str) \
        if "attributes" in results \
        else feature_name + "=" + results["hit"].fillna(0).astype(str)

    return results[columns_return]


def make_windows(annotations, window_size, step):
    """
    Split gene regions into multiple overlapping windows
    :param annotations: Table with genes and corresponding regions that should be processed
    :param window_size: Window size
    :param step: Step size
    :return:
    """
    annotations_bin_keys = {"gene_name", "gene_chrom", "gene_start", "gene_end", "gene_strand", "gene_region_end", "gene_region_start"}
    annotations_bin = {k: [] for k in annotations_bin_keys}
    annotations_bin["bin_start"] = []
    annotations_bin["bin_end"] = []
    annotations_bin["bin_strand"] = []
    for r, row in annotations.iterrows():
        numOfChunks = int((row["gene_region_end"] - row["gene_region_start"] - window_size) / step) + 1
        bins = list(range(0, numOfChunks * step, step))

        # Original strand
        annotations_bin["bin_start"].extend([int(i+row["gene_region_start"]) for i in bins])
        annotations_bin["bin_end"].extend([int(i+window_size+row["gene_region_start"]) for i in bins])
        annotations_bin["bin_strand"].extend([row["gene_strand"]]*len(bins))
        for k in annotations_bin_keys:
            annotations_bin[k].extend([row[k]]*len(bins))

        # Reverse strand
        annotations_bin["bin_start"].extend([int(i + row["gene_region_start"]) for i in bins])
        annotations_bin["bin_end"].extend([int(i + window_size + row["gene_region_start"]) for i in bins])
        annotations_bin["bin_strand"].extend(["-" if row["gene_strand"] == "+" else "+"]*len(bins))
        for k in annotations_bin_keys:
            annotations_bin[k].extend([row[k]]*len(bins))

    return pd.DataFrame.from_dict(annotations_bin)


def read_breaks_bed(bed_path):
    """
    Read file with breaks positions
    :param bed_path: Path to the breaks file
    :return: Table with breaks
    """
    cols = ["feature_chrom", "feature_start", "feature_end", "feature_name", "feature_score", "feature_strand"]
    dtypes = {'feature_start': 'int64', 'feature_end': 'int64'}

    return pd.read_csv(bed_path, sep="\t", header=0, dtype=dtypes, names=cols)


def read_gff(gff_path):
    """
    Read GFF file and select a longest transcript per each gene
    :param gff_path: Path to GFF file (possibly gziped)
    :return: Gene annotation table
    """
    cols = ["gene_chrom", "gene_source", "gene_feature", "gene_start", "gene_end", "gene_score", "gene_strand", "gene_frame", "gene_attribute"]
    dtypes = {'gene_start': 'int64', 'gene_end': 'int64'}
    annotations = pd.read_csv(gff_path, sep="\t", header=0, names=cols, dtype=dtypes)
    annotations["gene_name"] = annotations["gene_attribute"].str.extract('gene_id "([^"]+)"')
    annotations["gene_length"] = annotations["gene_end"] - annotations["gene_start"]
    annotations = annotations.query("gene_feature=='transcript' & ~gene_chrom.str.contains('_', regex=False)")

    # Select longest transcript to represent gene (TODO: this can be 1s faster if done inline)
    annotations = annotations.sort_values(["gene_strand", "gene_chrom", "gene_name", "gene_length"], ascending=False)
    annotations = annotations.drop_duplicates(["gene_strand", "gene_chrom", "gene_name"], keep="first")

    # Select columns and extend region
    return annotations[["gene_chrom", "gene_name", "gene_start", "gene_end", "gene_strand", "gene_length"]]


def remove_short_genes(annotations, size):
    """
    Remove genes that are shorter than suggested size
    :param annotations: table with gene annotations
    :param size: Minimal required size
    :return: Annotations table with short genes removed
    """
    small_genes = annotations.query("gene_region_end - gene_region_start < {window_size}".format(window_size=size))
    if len(small_genes) > 0:
        small_genes_names = small_genes["gene_name"].values
        small_genes_len = small_genes["gene_region_end"].values - small_genes["gene_region_start"].values
        small_genes_repr = ", ".join(
            "{}|{:0.0f}k".format(n, s / 1000) for n, s in zip(small_genes_names, small_genes_len))
        print("{n} genes were excluded because their regions were smaller than  ({size})bp".format(
            size=size, n=len(small_genes_names)))
        # print(small_genes_repr)
    return annotations.query("gene_region_end - gene_region_start >= {window_size}".format(window_size=size))


def write_aggregated_output(processed, output_path):
    """
    Write aggregated annotations to output path
    :param processed: Processed annotations
    :param output_path: Output path
    :return: Nothing
    """
    processed.to_csv(output_path, index=False, sep="\t")


def process_breaks_table(annotations_bin, breaks, details=False, additional_features_paths=[]):
    """
    Process a single breaks table

    :param annotations_bin: Table with intervals for which to calculate coverage
    :param breaks: Table containing DNA breaks positions
    :param additional_features_paths: Multiple paths to additional annotations
    :return: Processed annotations
    """
    # Calculate breaks coverage over genome annotation (GFF)
    processed = coverage(intervals=annotations_bin, features=breaks, feature_name="breaks", details=details)

    # Process additional annotations
    for f in additional_features_paths:
        print(f)

    return processed


def _call_process_breaks_table(annotations_bin, breaks, details, additional_features_paths, output_path):
    with Timeit(os.path.basename(output_path), pad=2):
        processed = process_breaks_table(annotations_bin=annotations_bin, breaks=breaks, details=details, additional_features_paths=additional_features_paths)
        write_aggregated_output(processed, output_path)


def _handle_error(error):
    pass
    # print(error)
    # raise error


def main(args):
    with Timeit('reading GFF file (only select longest transcript per gene)', display_on_enter=True):
        annotations = read_gff(args.annotations)

        # Select columns and extend region
        annotations = annotations[["gene_chrom", "gene_name", "gene_start", "gene_end", "gene_strand", "gene_length"]]
        annotations["gene_region_start"] = annotations["gene_start"] - args.extend
        annotations["gene_region_end"] = annotations["gene_end"] + args.extend

    with Timeit('Removed short genes'):
        annotations = remove_short_genes(annotations, args.window_size)

    with Timeit('splitting gene regions to intervals', display_on_enter=True):
        annotations_bin = make_windows(annotations, args.window_size, args.window_step)

    multiple_breaks_bed_paths = natsorted([path for input_glob in args.inputs for path in glob(input_glob)])
    threads = min([args.threads, len(multiple_breaks_bed_paths)])
    async_results = []
    with Timeit('calculating coverage for {n} breaks files using {threads} threads'.format(n=len(multiple_breaks_bed_paths), threads=threads), display_on_enter=True):
        pool = mp.Pool(threads)
        for breaks_bed_path in multiple_breaks_bed_paths:
            breaks = read_breaks_bed(breaks_bed_path)

            # Figure out output file name if needed
            output_name = args.output_name
            match_sed = re.fullmatch("/(.+)(?<!\\\)/(.+)/", args.output_name)
            if match_sed:
                pattern, repl = match_sed.groups()
                output_name = re.sub(pattern, repl, os.path.basename(breaks_bed_path)) + ".tsv"
            output_path = os.path.join(args.output_dir, output_name)

            # create output directory
            output_dir = os.path.dirname(output_path)
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
                print("Created {} folder".format(output_dir))

            # Create an asynchronous process for calculating coverage
            r = pool.apply_async(_call_process_breaks_table, kwds={
                'annotations_bin': annotations_bin,
                'breaks': breaks,
                'details': args.details,
                'additional_features_paths': args.additional_features_paths,
                'output_path': output_path}, error_callback=_handle_error)
            async_results.append(r)
        pool.close()
        pool.join()

        for r in async_results:
            r.get()



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Use a sliding window to aggregate breaks in bed file')
    parser.add_argument('inputs', nargs='+', help='Input .bed file with detected breaks. Can also be multiple files or a whildcard expression (e.g.: path/to/*.bed) ')
    parser.add_argument('annotations', help='Annotation file. If annotation file has gtf or gff extention (possibly .gz) then only transcripts are selected. If .bed file is provided then all annotations from bed file are used')
    parser.add_argument('output_dir', help='Directory to which output file will be written')
    parser.add_argument('--output-name', dest="output_name", default="/(.*)\\.[a-z]+/\\1/", help='Name of output file. Can be either a simple name (when applying to single file) or a SED like regular expression (e.g.: /pattern/replace/')
    parser.add_argument('-w|--window-size', dest="window_size", default=int(1e5), type=int, help='Window at which to agregate breaks number')
    parser.add_argument('-s|--window-step', dest="window_step", default=int(1e4), type=int, help='Step after each window')
    parser.add_argument('-e|--extend-gene', dest="extend", default=0, type=int, help='Extend each gene both directions by number of base pairs')
    parser.add_argument('-f|--features', dest="additional_features_paths", action="append", default=[], nargs="*", help='Additional features to annotate input file (need to provide additional -f flag with each feature)')
    parser.add_argument('--details', dest="details", action='store_true', help='Preserve details like break names in the output file (default: no)')
    parser.add_argument('--threads', dest="threads", default=mp.cpu_count(), type=int, help='Number of threads used to calculate coverage (default: {})'.format(mp.cpu_count()))
    args = parser.parse_args()

    main(args)
