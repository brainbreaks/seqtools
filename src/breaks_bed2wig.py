import argparse
import pyranges as pr
import pandas as pd
from src.functions import *
from glob import glob
import os
import re


def main(args):
    multiple_breaks_bed_paths = [path for input_glob in args.inputs for path in glob(input_glob)]
    breaks = []
    for path in multiple_breaks_bed_paths:
        chr = re.sub(".*(chr[^_]+)_.*", r"\1", os.path.basename(path), flags=re.IGNORECASE).lower()
        breaks_chr = pd.read_csv(path, sep="\t", header=0, names=["Chromosome", "Start", "End", "Feature", "Score", "Strand"]).\
            query("Chromosome == @chr")
        breaks.append(breaks_chr)
    breaks = pd.concat(breaks)

    window_size = args.window_size
    step = args.window_step
    ann_chromosomes = breaks.groupby(["Chromosome"], as_index=False).agg({'Start': min, 'End': max})
    bin_chromosomes = {k: [] for k in ann_chromosomes.keys()}
    bin_chromosomes["Strand"] = []
    for r, row in ann_chromosomes.iterrows():
        numOfChunks = int((row["End"] - row["Start"] - window_size) / step) + 1
        bins = list(range(0, numOfChunks * step, step))

        bin_chromosomes["Start"].extend([int(i+row["Start"]) for i in bins*2])
        bin_chromosomes["End"].extend([int(i+window_size+row["Start"]) for i in bins*2])
        bin_chromosomes["Strand"].extend(["+"]*len(bins) + ["-"]*len(bins))
        bin_chromosomes["Chromosome"].extend([row["Chromosome"]]*len(bins)*2)

    bin_chromosomes = pr.from_dict(bin_chromosomes)
    coverage_chromosomes = bin_chromosomes.coverage(pr.PyRanges(breaks), strandedness="same", overlap_col="Breaks").\
        as_df().query("Breaks>0").sort_values(["Strand", "Chromosome", "Start", "End"], ignore_index=True)

    for strand, coverage_strand in coverage_chromosomes.groupby(["Strand"]):
        for chr, coverage_chr in coverage_strand.groupby(["Chromosome"]):
            print(strand)
    pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Use a sliding window to aggregate breaks in bed file')
    parser.add_argument('inputs', nargs='+', help='Input .bed files with detected breaks. Can also be multiple files or a whildcard expression (e.g.: path/to/*.bed) ')
    parser.add_argument('output', help='Path to output wig file')
    parser.add_argument('-w|--window-size', dest="window_size", default=int(1e5), type=int, help='Window at which to agregate breaks number')
    parser.add_argument('-s|--window-step', dest="window_step", default=int(1e4), type=int, help='Step after each window')
    args = parser.parse_args()

    main(args)