# need
# sudo rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/ /usr/bin/ucsc/
#
import pysam
import argparse
import pyranges as pr
import pandas as pd
# from src.functions import *
from glob import glob
import os
import re
import subprocess
import math
from sklearn.neighbors import KernelDensity
from scipy.stats import norm
import numpy as np
from natsort import natsorted
import matplotlib
import tempfile
import multiprocessing as mp


def ranges_kde(breaks_df, chromsizes_df, binwidth=100, bandwidth=500, algorithm="epanechnikov"):
    df = {"Start": [], "End": [], "Strand": [], "Chromosome": [], "Score": []}
    for chr in natsorted(list(set(breaks_df.Chromosome))):
        if not (chromsizes_df.Chromosome == chr).any():
            print("No breaks on chromosome '{}'".format(chr))
            continue

        x_df = breaks_df.loc[breaks_df.Chromosome == chr]

        if x_df.shape[0] < 10:
            print("To few breaks on chromosome '{}' ({} < 10)".format(chr, x_df.shape[0]))
            continue

        csize = chromsizes_df.Size[chromsizes_df.Chromosome == chr].values[0]
        bins = int(math.ceil(csize/binwidth))
        x_df_starts = np.linspace(1, csize, bins).round()
        x_df_model = KernelDensity(kernel=algorithm, bandwidth=bandwidth).fit(x_df["Start"].values[:, np.newaxis])
        x_df_fit = np.exp(x_df_model.score_samples(x_df_starts[:-1, np.newaxis]))
        df["Start"].extend([int(i) for i in x_df_starts[:-1]])
        df["End"].extend([int(i) for i in x_df_starts[1:]])
        df["Strand"].extend(["+"] * (len(x_df_starts) - 1))
        df["Chromosome"].extend([chr] * (len(x_df_starts) - 1))
        df["Score"].extend(x_df_fit * 1e8)

    df_ranges = pr.from_dict(df)

    return df_ranges

def to_bigwig(data, bigwig_path, chromsizes_df, chromsizes_path, name, description="", color="0,0,255", altColor="255,0,0"):
    bigwig_path1 = tempfile.NamedTemporaryFile(delete=False).name
    wig_path = tempfile.NamedTemporaryFile(delete=False).name
    # sizes_path = tempfile.NamedTemporaryFile(delete=False).name
    # chromsizes_df[["Chromosome", "Size"]].to_csv(sizes_path, header=False, index=False, sep="\t")
    data.to_bigwig(path=bigwig_path1, chromosome_sizes=pr.from_dict(chromsizes_df), value_col="Score", rpm=False)

    subprocess.run(["bigWigToWig", bigwig_path1, wig_path], env=os.environ)
    with open(wig_path, 'r') as original: data = original.read()
    with open(wig_path, 'w') as modified:
        modified.write(
            "track type=wiggle_0 name=\"{name}\" visibility=full description=\"This track represents joins to similar strand\" color={color} altColor={altColor}\n".format(
                name=name, url=description, color=color, altColor=altColor) + data)

    subprocess.run(["wigToBigWig", wig_path, chromsizes_path, bigwig_path], env=os.environ)

def main_group(args, group, breaks_df):
    sample = re.sub("\.bed", "", os.path.basename(group))
    print("Processing {}".format(group))
    breaks_group_df = breaks_df.query("group==@group")
    breaks_group_pos_df = breaks_group_df.query("Strand=='+'").copy()
    breaks_group_neg_df = breaks_group_df.query("Strand=='-'").copy()
    chromsizes_df = pd.read_csv(args.chromsizes, sep="\t", header=None, names=["Chromosome", "Size"])
    chromsizes_df["Start"] = 1
    chromsizes_df["End"] = chromsizes_df["Size"]
    # chromsizes_df = breaks_group_df.groupby(["Chromosome"], as_index=False).agg({'Start': min, 'End': max})
    # chromsizes_df["Start"] = 1
    # chromsizes_df["Size"] = chromsizes_df["End"] - chromsizes_df["Start"] + 1

    if args.method == "kde":
        breaks_kde_pos = ranges_kde(breaks_group_pos_df, chromsizes_df=chromsizes_df, binwidth=args.kde_binwidth,
                                    bandwidth=args.kde_bandwidth, algorithm=args.kde_algorithm)
        breaks_kde_neg = ranges_kde(breaks_group_neg_df, chromsizes_df=chromsizes_df, binwidth=args.kde_binwidth,
                                    bandwidth=args.kde_bandwidth, algorithm=args.kde_algorithm)
        kde_max = max([breaks_kde_pos.Score.max(), breaks_kde_neg.Score.max()])
        breaks_kde_pos.Score = breaks_kde_pos.Score / kde_max
        breaks_kde_neg.Score = -breaks_kde_neg.Score / kde_max
        to_bigwig(breaks_kde_pos, os.path.join(args.output_path, sample + "_pos_kde.bw"), chromsizes_df,
                  args.chromsizes, name="{} +".format(os.path.basename(group)))
        to_bigwig(breaks_kde_neg, os.path.join(args.output_path, sample + "_neg_kde.bw"), chromsizes_df,
                  args.chromsizes, name="{} -".format(os.path.basename(group)))
    elif args.method == "window":
        bin_chromosomes = {k: [] for k in ['Chromosome', 'Start', 'End']}
        bin_chromosomes["Strand"] = []
        for r, row in chromsizes_df.iterrows():
            numOfChunks = int((row["End"] - row["Start"] - args.window_size) / args.window_step) + 1
            bins = list(range(0, numOfChunks * args.window_step, args.window_step))

            bin_chromosomes["Start"].extend([int(i + row["Start"]) for i in bins * 2])
            bin_chromosomes["End"].extend([int(i + args.window_size + row["Start"]) for i in bins * 2])
            bin_chromosomes["Strand"].extend(["+"] * len(bins) + ["-"] * len(bins))
            bin_chromosomes["Chromosome"].extend([row["Chromosome"]] * len(bins) * 2)
        bin_chromosomes = pr.from_dict(bin_chromosomes)

        coverage_chromosomes = bin_chromosomes.coverage(pr.PyRanges(breaks_group_df), strandedness="same",
                                                        overlap_col="Score")
        coverage_chromosomes_df = coverage_chromosomes.as_df()
        window_max = coverage_chromosomes_df.Score.max()
        coverage_chromosomes_df.loc[(coverage_chromosomes_df["Strand"] == "-"), "Score"] = -coverage_chromosomes_df.loc[
            (coverage_chromosomes_df["Strand"] == "-"), "Score"]
        coverage_chromosomes_df["End"] = coverage_chromosomes_df["Start"] + math.floor(
            (args.window_size + args.window_step) / 2) - 1
        coverage_chromosomes_df["Start"] = coverage_chromosomes_df["Start"] + math.floor(
            (args.window_size - args.window_step) / 2)
        coverage_chromosomes_df["Score"] = coverage_chromosomes_df["Score"] / window_max

        coverage_chromosomes = pr.PyRanges(coverage_chromosomes_df)
        to_bigwig(coverage_chromosomes[coverage_chromosomes.Strand == "+"],
                  os.path.join(args.output_path, sample + "_pos_window.bw"), chromsizes_df=chromsizes_df,
                  chromsizes_path=args.chromsizes, name=sample + " +")
        to_bigwig(coverage_chromosomes[coverage_chromosomes.Strand == "-"],
                  os.path.join(args.output_path, sample + "_neg_window.bw"), chromsizes_df=chromsizes_df,
                  chromsizes_path=args.chromsizes, name=sample + " -")
    else:
        raise NotImplementedError("Unknown method {}".format(args.method))

def main(args):
    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path, exist_ok=True)

    # Create template with sliding window
    breaks_df = []
    multiple_breaks_bed_paths = [path for input_glob in args.inputs for path in glob(input_glob)]
    for path in multiple_breaks_bed_paths:
        filename = os.path.basename(path)
        bait_chrom = re.sub(".*(chr[^_]+)_.*", r"\1", filename, flags=re.IGNORECASE).lower()
        breaks_i = pd.read_csv(path, sep="\t", header=0, names=["Chromosome", "Start", "End", "Feature", "Score", "Strand"])
        breaks_i["bait_chrom"] = bait_chrom
        breaks_i["filename"] = filename
        breaks_i["group"] = path

        if args.only_bait_chromosome:
            breaks_i = breaks_i.query("Chromosome == bait_chrom").copy()
            breaks_i["group"] = "only_bait_chromosome"

        breaks_df.append(breaks_i)
    breaks_df = pd.concat(breaks_df)


    # Count values in each window (per group)
    groups = sorted(list(set(breaks_df.group)))
    pool = mp.Pool(args.threads)
    for group in groups:
        pool.apply_async(main_group, kwds={'args': args, 'group': group, 'breaks_df': breaks_df})
    pool.close()
    pool.join()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Use a sliding window to aggregate breaks in bed file')
    parser.add_argument('inputs', nargs='+', help='Input .bed files with detected breaks. Can also be multiple files or a whildcard expression (e.g.: path/to/*.bed)')
    parser.add_argument('chromsizes', help='Path to chromsizes file')
    parser.add_argument('output_path', help='Path to folder where all information needed to import to UCSC genome browser will be stored')
    parser.add_argument("method", help="Can be either 'kde' (additional arguments --kde-bandwidth, --kde-bins, --kde-algorithm) or 'window' (additional arguments --sindow-size and --window-step)")
    parser.add_argument('--only-bait-chromosome', dest="only_bait_chromosome", action="store_true", help='Enabling this option will produce a WIG where only bait chromosome is present from each file. The bait chromosome is guessed from file name')
    parser.add_argument('--threads', dest="threads", default=1, type=int, help='Use multiple threads to process different files')
    parser.add_argument('-w', '--window-size', dest="window_size", default=int(1e5), type=int, help='Window at which to agregate breaks number')
    parser.add_argument('-s', '--window-step', dest="window_step", default=int(1e4), type=int, help='Step after each window')
    parser.add_argument('-n', '--kde-bandwidth', dest="kde_bandwidth", default=int(1e3), type=int, help='Number of bins')
    parser.add_argument('-b', '--kde-binwidth', dest="kde_binwidth", default=int(1e5), type=int, help='Number of bins')
    parser.add_argument('-a', '--kde-algorithm', dest="kde_algorithm", default="epanechnikov", help='Algorithm used for kernel density estimation')
    # args = parser.parse_args(["app", "data/breaks/*_DMSO.bed", "data/breaks_window", "window", "--window-size", "100000", "--window-step", "10000"])
    # args = parser.parse_args(["app", "data/breaks/*_DMSO.bed", "data/breaks_kde", "kde", "--kde-bins", "100000", "--kde-bandwidth", "100", "--kde-algorithm", "epanechnikov"])
    args = parser.parse_args()

    main(args)