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
import matplotlib.pyplot as plt
import tempfile
import multiprocessing as mp
from scipy.stats.stats import pearsonr

def calculate_crosscorrelation(bed, name, lags=100, bandwidth=100, binwidth=100):
    sizes = bed.groupby(["Chromosome", "Strand"], as_index=False).agg({"Start":min, "End": max}).groupby("Chromosome").agg({"Start":min, "End": max}).to_dict()
    chrom_offsets = bed.groupby("Chromosome", as_index=False).agg({"Start": min, "End": max})
    chrom_offsets["offset_length"] = chrom_offsets.End - chrom_offsets.Start + 1
    chrom_offsets["offset_start"] = chrom_offsets.offset_length.cumsum() - chrom_offsets.offset_length - chrom_offsets.Start
    chrom_offsets["offset_end"] = chrom_offsets.offset_length.cumsum() - chrom_offsets.Start
    chrom_offsets = chrom_offsets[["Chromosome", "offset_start", "offset_end", "offset_length"]]

    bed_offset = bed.merge(chrom_offsets, on="Chromosome")
    bed_offset.loc[:,"Start"] += bed_offset.loc[:,"offset_start"]

    bed_sense = bed_offset.query("Strand=='+'")
    bed_anti = bed_offset.query("Strand=='-'")

    model_sense = KernelDensity(kernel="epanechnikov", bandwidth=bandwidth).fit(bed_sense["Start"].values[:, np.newaxis])
    model_anti = KernelDensity(kernel="epanechnikov", bandwidth=bandwidth).fit(bed_anti["Start"].values[:, np.newaxis])

    prebin = 1e5
    prekde_space = np.arange(0, chrom_offsets["offset_end"].max(), prebin)
    prekde_sense = np.exp(model_sense.score_samples(prekde_space[:, np.newaxis]))
    prekde_anti = np.exp(model_anti.score_samples(prekde_space[:, np.newaxis]))

    q = [0.7, 0.9]
    f1 = (prekde_anti>0) | (prekde_sense>0)
    f = ((np.quantile(prekde_sense[f1], q[0])<prekde_sense) & (prekde_sense<np.quantile(prekde_sense[f1], q[1]))) | ((np.quantile(prekde_anti[f1], q[0])<prekde_anti) & (prekde_anti<np.quantile(prekde_anti[f1], q[1])))

    kde_space = np.array([ss for s in prekde_space[f] for ss in np.arange(s, s+prebin, binwidth)])
    # kde_space = np.arange(binwidth*10, chrom_offsets["offset_end"].max() - binwidth*10, binwidth)
    kde_sense = np.exp(model_sense.score_samples(kde_space[:, np.newaxis]))
    kde_anti = np.exp(model_anti.score_samples(kde_space[:, np.newaxis]))

    q = [0.1, 0.9]
    f1 = (kde_anti>0) | (kde_sense>0)
    f = ((np.quantile(kde_sense[f1], q[0])<kde_sense) & (kde_sense<np.quantile(kde_sense[f1], q[1]))) | ((np.quantile(kde_anti[f1], q[0])<kde_sense) & (kde_sense<np.quantile(kde_anti[f1], q[1])))
    kde_anti = kde_anti[f]
    kde_sense = kde_sense[f]

    lags = np.arange(0, lags)*binwidth
    r = [pearsonr(kde_sense[:len(kde_sense)-i], kde_anti[i:])[1] for i in range(len(lags))]

    fig, ax = plt.subplots()

    plt.title(name)

    ax.plot(lags, r)
    # plt.show()
    print(os.getcwd())
    plt.savefig("{}.png".format(name))
    print("1")


    pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Use a sliding window to aggregate breaks in bed file')
    parser.add_argument('inputs', nargs='+', help='Input .bed files with detected breaks. Can also be multiple files or a whildcard expression (e.g.: path/to/*.bed)')
    parser.add_argument('-n', '--kde-bandwidth', dest="kde_bandwidth", default=int(1e3), type=int, help='Number of bins')
    parser.add_argument('-b', '--kde-binwidth', dest="kde_binwidth", default=int(1e5), type=int, help='Number of bins')
    parser.add_argument('--kde-lags', dest="kde_lags", default=int(1e3), type=int, help='Number of bins')
    parser.add_argument('-a', '--kde-algorithm', dest="kde_algorithm", default="epanechnikov", help='Algorithm used for kernel density estimation')
    # args = parser.parse_args(["app", "data/breaks/*_DMSO.bed", "data/breaks_window", "window", "--window-size", "100000", "--window-step", "10000"])
    args = parser.parse_args()

    # args = parser.parse_args(["app", "data/breaks/*_DMSO.bed", "data/breaks_kde", "kde", "--kde-binwidth", "10000", "--kde-bandwidth", "10000", "--kde-algorithm", "epanechnikov"])

    multiple_breaks_bed_paths = [path for input_glob in args.inputs for path in glob(input_glob)]
    for path in multiple_breaks_bed_paths:
        bed = pd.read_csv(path, sep="\t", header=0, names=["Chromosome", "Start", "End", "Feature", "Score", "Strand"])
        calculate_crosscorrelation(bed, name=os.path.basename(path), lags=args.kde_lags, binwidth=args.kde_binwidth, bandwidth=args.kde_bandwidth)