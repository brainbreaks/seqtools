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
    for chr in set(bed["Chromosome"]):
        kde_space = np.arange(sizes["Start"][chr] + binwidth*10, sizes["End"][chr] - binwidth*10, binwidth)[:, np.newaxis]

        bed_chrom_sense = bed.query("Strand=='+' & Chromosome==@chr")
        model_sense = KernelDensity(kernel="epanechnikov", bandwidth=bandwidth).fit(bed_chrom_sense["Start"].values[:, np.newaxis])
        kde_sense = np.exp(model_sense.score_samples(kde_space))

        bed_chrom_anti = bed.query("Strand=='-' & Chromosome==@chr")
        model_anti = KernelDensity(kernel="epanechnikov", bandwidth=bandwidth).fit(bed_chrom_anti["Start"].values[:, np.newaxis])
        kde_anti = np.exp(model_anti.score_samples(kde_space))

        lags = np.arange(0, lags)*binwidth
        r = [pearsonr(kde_sense[:len(kde_sense)-i], kde_anti[i:])[1] for i in range(len(lags))]

        fig, ax = plt.subplots()
        plt.title(name)
        ax.plot(lags, r)
        plt.show()



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