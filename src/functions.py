import pandas as pd
import time


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