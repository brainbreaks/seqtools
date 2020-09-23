import os
import re
import multiprocessing
from pybedtools import BedTool
from pybedtools.featurefuncs import gff2bed, extend_fields
import pandas as pd

def filter_strand(feature, strand):
    return feature.strand==strand

def filter_transcript(feature):
    return not re.search("_", feature.fields[0]) and feature.fields[2]=="transcript"

def filter_covered(feature):
    coverage = float(feature.fields[9]) > 1e-6
    return coverage

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
    if not os.path.exists('../tmp'):
        os.makedirs('../tmp')

    max_size = 1e6
    overlap_size=0

    dna_breaks = BedTool("data/Chr1_APH_no10kb_Merge.bed")
    genome_bin_pos = BedTool().window_maker(genome="mm10", w=max_size, s=max_size-overlap_size).each(strand, "+").saveas("tmp/genome_bin_pos.bed")
    genome_bin_neg = BedTool().window_maker(genome="mm10", w=max_size, s=max_size-overlap_size).each(strand, "-").saveas("tmp/genome_bin_neg.bed")
    genome_bin = genome_bin_pos.cat(genome_bin_neg, postmerge=False).sort().saveas("tmp/genome_bin.bed")
    annotations = BedTool("data/mm10.refGene.gtf.gz")

    transcripts = annotations.filter(filter_transcript).\
        each(gff2bed, name_field="gene_id").sort().\
        saveas("tmp/mm10.refGene.all_transcripts.bed").\
        groupby(g="1,2,3,6", c="4,5", o="distinct").\
        cut([0,1,2,4,5,3]).\
        saveas("tmp/mm10.refGene.transcripts.bed")


    bin_breaks = BedTool().intersect(a=genome_bin, b=dna_breaks, wa=True, c=True, s=True). \
        saveas("tmp/bin_breaks.bed")

    results = BedTool().map(a=bin_breaks, b=transcripts, s=True, c="4", o="distinct").cut([0,1,2,7,6,5]).sort().saveas("tmp/results.bed")
    print(results)

    results_df = splitDataFrameList(results.to_dataframe(), "name", ",")
    results_df.to_csv("results.tsv", sep="\t", header=True, index=False)


    with pd.option_context("display.max_rows", 10, "display.max_columns", 15):
        print(results_df)
    exit()

    # .filter(filter_transcript)

    dna_breaks_pos = dna_breaks.filter(filter_strand, strand="+")
    dna_breaks_neg = dna_breaks.filter(filter_strand, strand="-")

    bin_breaks_pos = BedTool().intersect(a=genome_bin, b=dna_breaks_pos, s=True)
    bin_breaks_neg = BedTool().intersect(a=genome_bin, b=dna_breaks_neg, s=True)


    binann_breaks_pos = BedTool().coverage(a=transcripts, b=bin_breaks_pos).sort().filter(filter_covered).saveas("tmp/results_pos.bed") #
    binann_breaks_neg = BedTool().coverage(a=transcripts, b=bin_breaks_neg).sort().filter(filter_covered).saveas("tmp/results_neg.bed") #
    # print(len(bin_breaks_neg))
    # print(len(bin_breaks_pos))
    print(len(binann_breaks_pos))
    print(len(binann_breaks_neg))
    exit()
    # results = pd.concat(binann_breaks_pos.to_dataframe(), binann_breaks_neg.to_dataframe())
    #results_pos.columns[6:9] = ['unknown', 'binstart', 'binend', 'coverage']

    with pd.option_context("display.max_rows", 10, "display.max_columns", 15):
       # print(binann_breaks_pos.to_dataframe())
        print(binann_breaks_neg.to_dataframe())
        # print(results)
    # print(binann_breaks_pos.to_dataframe())

    exit()

    print("Filter: ")
    print(dna_breaks_pos)
    print("Bin: ")
    print(bin_breaks_pos)
    # .saveas(tx_out_file)

    dna_breaks_pos.window_maker()
    # dna_breaks
    #print(dna_breaks_pos)
    print(dna_breaks_neg)
    print(len(dna_breaks_neg))
    print(len(dna_breaks_pos))

    #pool = multiprocessing.Pool(processes=args.processes)

    # genes = BedTool('hg19.gff')
    #
    # intergenic_snps = snps.subtract(genes)
    # nearby = genes.closest(intergenic_snps, d=True, stream=True)
    #
    # for gene in nearby:
    #     if int(gene[-1]) < 5000:
    #         print(gene.name)

if __name__ == "__main__":
    main()