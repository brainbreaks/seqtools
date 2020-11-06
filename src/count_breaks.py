import os
import re
import sys
import time
from pybedtools import BedTool
from pybedtools.featurefuncs import gff2bed, extend_fields
import pandas as pd
import argparse
import tempfile
import warnings
import itertools
warnings.filterwarnings("ignore", message=".*buffering.*")

def filter_strand(feature, strand):
    return feature.strand==strand

def filter_transcript(feature):
    return not re.search("_", feature.fields[0]) and feature.fields[2]=="transcript"

def strand(f, strand="+"):
    f = extend_fields(f, 6)
    f.strand = strand
    return f

def splitDataFrameList(df,target_column,delimiters):
    ''' df = dataframe to split,
    target_column = the column containing the values to split
    separator = the symbol used to perform the split
    returns: a dataframe with each entry for the target column separated, with each element moved into a new row.
    The values in the other columns are duplicated across the newly divided rows.
    '''
    regexPattern = "|".join(map(re.escape,delimiters))
    def splitListToRows(row,row_accumulator,target_column,regexPattern):
        split_row = re.split(regexPattern,row[target_column])
        for s in split_row:
            new_row = row.to_dict()
            new_row[target_column] = s
            row_accumulator.append(new_row)
    new_rows = []
    df.apply(splitListToRows,axis=1,args = (new_rows,target_column,regexPattern))
    new_df = pd.DataFrame(new_rows)
    return new_df


def main():
    parser = argparse.ArgumentParser(description='Use a sliding window to aggregate breaks in bed file')
    parser.add_argument('genome', help='Name of the model used to produce input')
    parser.add_argument('input', help='Input .bed file with detected breaks')
    parser.add_argument('annotations', help='Annotation file. If annotation file has gtf or gff extention (possibly .gz) then only transcripts are selected. If .bed file is provided then all annotations from bed file are used')
    parser.add_argument('output', help='Output .bed file with longest transcripts')
    parser.add_argument('-w|--window-size', dest="window_size", default=int(1e5), type=int, help='Window at which to agregate breaks number')
    parser.add_argument('-s|--window-step', dest="window_step", default=int(1e4), type=int, help='Step after each window')
    parser.add_argument('-f|--features', dest="features", action="append", nargs="*", help='Additional features to annotate input file')

    args = parser.parse_args()
    start = time.time()

    if args.features is None:
        features = []
    else:
        features = list(itertools.chain.from_iterable(args.features))

    output_dir = os.path.dirname(args.output)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print('Processing "{input}" using annotation="{annotation}" window {window}/{step}. Writing output to "{output}"...'.format(
            input=args.input, window=args.window_size, step=args.window_step, output=args.output, annotation=args.annotations))

    # Create temporary files
    tmp = {n: tempfile.NamedTemporaryFile(delete=False).name for n in ["genome_bin_pos", "genome_bin_neg", "genome_bin", "breaks_bin", "results", "all_transcripts", "transcripts"]}

    # Create windows template for sliding window
    genome_bin_pos = BedTool().window_maker(genome=args.genome, w=args.window_size, s=args.window_step).each(strand, "+").saveas(tmp["genome_bin_pos"])
    genome_bin_neg = BedTool().window_maker(genome=args.genome, w=args.window_size, s=args.window_step).each(strand, "-").saveas(tmp["genome_bin_neg"])
    genome_bin = genome_bin_pos.cat(genome_bin_neg, postmerge=False).sort().saveas(tmp["genome_bin"])

    # Read input file
    dna_breaks = BedTool(args.input)

    # Read annotation file
    if re.search(r"\.(gtf|gff)(\.gz)?$", args.annotations):
        annotations = BedTool(args.annotations)
        annotations = annotations.filter(filter_transcript).\
            each(gff2bed, name_field="gene_id").sort().\
            saveas(tmp["all_transcripts"]).\
            groupby(g="1,2,3,6", c="4,5", o="distinct").\
            cut([0,1,2,4,5,3]).\
            saveas(tmp["transcripts"])
    elif re.search(r"\.bed$", args.annotations):
        annotations = BedTool(args.annotations)
    else:
        parser.error("Annotation have to be either in gtf/gff or in bed format")

    bin_breaks = BedTool().intersect(a=genome_bin, b=dna_breaks, wa=True, c=True, s=True). \
        saveas(tmp["breaks_bin"])

    # Map breaks statistics to annotation file
    results = BedTool().map(a=bin_breaks, b=annotations, c="4", o="distinct").cut([0,1,2,7,6,5]).sort().saveas(tmp["results"]) # s=True,
    results_df = splitDataFrameList(results.to_dataframe(), "name", ",")
    results_df = results_df[results_df.name != "."]
    results_df.to_csv(args.output, sep="\t", header=True, index=False)

    # Remove old temporary files
    for f in tmp.values():
        os.remove(f)

    end = time.time()
    print("Total time: {:.1f} minutes".format((end - start)/60))
    #
    # with pd.option_context("display.max_rows", 10, "display.max_columns", 15):
    #     print(results_df)

if __name__ == "__main__":
    main()