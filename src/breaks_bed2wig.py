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

def write_coverages_bigwig(coverages, chromsizes_df, bigwig_path, track_name=None):
    if track_name is None:
        track_name = os.path.basename(bigwig_path)

    wig_path = bigwig_path + ".wig"
    chrsizes_path = bigwig_path + ".sizes"
    chromsizes_df[["Chromosome", "Size"]].to_csv(chrsizes_path, header=False, index=False, sep="\t")
    coverages.to_bigwig(path=bigwig_path, chromosome_sizes=pr.PyRanges(chromsizes_df), value_col="Breaks")
    # os.system("bigWigToWig {bigwig} {wig}".format(bigwig=bigwig_path, wig=wig_path))
    print(chromsizes_df)
    subprocess.run(["bigWigToWig", bigwig_path, wig_path], env=os.environ)

    with open(wig_path, 'r') as original: data = original.read()
    with open(wig_path, 'w') as modified:
        modified.write(
            "track type=bigWig name=\"{name}\" description=\"This track represents joins to similar strand\" color=255,0,0\n".format(
                name=track_name, url=os.path.basename(bigwig_path)) + data)

    subprocess.run(["wigToBigWig", wig_path, chrsizes_path, bigwig_path + ".1"], env=os.environ)

def main(args):
    # Create template with sliding window
    breaks_df = []
    multiple_breaks_bed_paths = [path for input_glob in args.inputs for path in glob(input_glob)]
    for path in multiple_breaks_bed_paths:
        filename = os.path.basename(path)
        bait_chrom = re.sub(".*(chr[^_]+)_.*", r"\1", filename, flags=re.IGNORECASE).lower()
        breaks_i = pd.read_csv(path, sep="\t", header=0,
                               names=["Chromosome", "Start", "End", "Feature", "Score", "Strand"])
        breaks_i["bait_chrom"] = bait_chrom
        breaks_i["filename"] = filename
        breaks_i["group"] = filename

        if args.only_bait_chromosome:
            breaks_i = breaks_i.query("Chromosome == bait_chrom")
            breaks_i["group"] = "only_bait_chromosome"

        breaks_df.append(breaks_i)
    breaks_df = pd.concat(breaks_df)


    # Count values in each window (per group)
    for group in sorted(list(set(breaks_df.group))):
        print("Processing {}".format(group))
        breaks_group_df = breaks_df.query("group==@group")
        chromsizes_df = breaks_group_df.groupby(["Chromosome"], as_index=False).agg({'Start': min, 'End': max})

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

        coverage_chromosomes = bin_chromosomes.coverage(pr.PyRanges(breaks_group_df), strandedness="same", overlap_col="Breaks")
        coverage_chromosomes_df = coverage_chromosomes.as_df()
        coverage_chromosomes_df.loc[(coverage_chromosomes_df["Strand"] == "-"), "Breaks"] = -coverage_chromosomes_df.loc[(coverage_chromosomes_df["Strand"] == "-"), "Breaks"]
        coverage_chromosomes_df["Start"] = coverage_chromosomes_df["Start"] + math.floor(args.window_size/2 - args.window_step*2)
        coverage_chromosomes_df["End"] = coverage_chromosomes_df["Start"] + math.floor(args.window_size/2 - args.window_step)
        coverage_chromosomes = pr.PyRanges(coverage_chromosomes_df)

        if not os.path.exists(args.output_path):
            os.makedirs(args.output_path, exist_ok=True)

        # coverage_chromosomes_df = breaks_group_df.copy()
        # coverage_chromosomes_df["Start"] = coverage_chromosomes_df["Start"]-50
        # coverage_chromosomes_df["End"] = coverage_chromosomes_df["End"]+50
        # coverage_chromosomes_df["Breaks"] = 1
        # coverage_chromosomes_df = coverage_chromosomes_df.sort_values(["Chromosome", "Strand", "Start", "End"], ignore_index=True)
        # coverage_chromosomes = pr.PyRanges(coverage_chromosomes_df)

        coverage_chromsizes_df = coverage_chromosomes_df.groupby(["Chromosome"], as_index=False).agg({'Start': min, 'End': max})
        coverage_chromsizes_df["Size"] = coverage_chromsizes_df["End"]
        basename_bigwig_pos, basename_bigwig_neg = "{}_pos.bw".format(group), "{}_neg.bw".format(group)
        path_bigwig_pos, path_bigwig_neg = os.path.join(args.output_path, basename_bigwig_pos), os.path.join(args.output_path, basename_bigwig_neg)
        write_coverages_bigwig(coverage_chromosomes[coverage_chromosomes.Strand == "+"], coverage_chromsizes_df, path_bigwig_pos, track_name="{}:+".format(group))
        write_coverages_bigwig(coverage_chromosomes[coverage_chromosomes.Strand == "-"], coverage_chromsizes_df, path_bigwig_neg, track_name="{}:-".format(group))

        # with open(os.path.join(args.output_path, "{}_custom_tracks.txt".format(group)), 'w') as f:
        #     f.write("#\n# You need to manually replace url to positive and negative strand tracks and \n# add each custom track individually to UCSC genome browser\n#\n")
        #     f.write("track type=bigWig name=\"{name} (pos)\" description=\"This track represents joins to similar strand\" color=255,0,0, bigDataUrl={url}".format(name="{}:{}".format(args.track_name, group), url=basename_pos) + "\n")
        #     f.write("track type=bigWig name=\"{name} (neg)\" description=\"This track represents joins to opposite strand\" color=0,255,0, bigDataUrl={url}".format(name="{}:{}".format(args.track_name, group), url=basename_neg) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Use a sliding window to aggregate breaks in bed file')
    parser.add_argument('inputs', nargs='+', help='Input .bed files with detected breaks. Can also be multiple files or a whildcard expression (e.g.: path/to/*.bed) ')
    parser.add_argument('output_path', help='Path to folder where all information needed to import to UCSC genome browser will be stored')
    parser.add_argument('--track-name', dest="track_name", help='Name of the UCSC track')
    parser.add_argument('--only-bait-chromosome', dest="only_bait_chromosome", action="store_true", help='Enabling this option will produce a WIG where only bait chromosome is present from each file. The bait chromosome is guessed from file name')
    parser.add_argument('-w|--window-size', dest="window_size", default=int(1e5), type=int, help='Window at which to agregate breaks number')
    parser.add_argument('-s|--window-step', dest="window_step", default=int(1e4), type=int, help='Step after each window')
    args = parser.parse_args()
    parser.parse_args(["data/breaks/*_DMSO.bed", "data/mm9/mm9.chrom.sizes", "data/breaks_wig", "-w", "100000", "-s", "10000"])

    if args.track_name is None:
        args.track_name = os.path.basename(args.output_path)

    main(args)