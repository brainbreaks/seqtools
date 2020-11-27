import glob
import pandas
import re
import os
import subprocess

breaksites_df = pandas.read_csv("data/breaksites_primers.tsv", sep="\t")

macs2_cols = ["junction_chrom", "junction_start", "junction_end", "junction_name", "junction_score", "junction_strand"]
for macs2_input in glob.glob("data/breaks/*.bed"):
    macs2_sample = re.sub(r"\.bed$", "", os.path.basename(macs2_input))
    macs2_outdir = "data/macs2/"
    macs2_input_filtered = re.sub(r"\.bed$", "_filtered.bad", os.path.join("tmp", os.path.basename(macs2_input)))

    print("Processing {}...".format(macs2_sample))

    macs2_input_df = pandas.read_csv(macs2_input, sep="\t", header=None, names=macs2_cols).\
        merge(breaksites_df, how="left", left_on="junction_chrom", right_on="bait_chrom").\
        query("((primer_start - 1500000) > junction_start) | (junction_end > (primer_end + 1500000))")
    macs2_input_df[macs2_cols].to_csv(macs2_input_filtered, sep="\t", index=False, header=False)

    macs2_cmd = "macs2 callpeak -t {input} -f BED -g mm --keep-dup all -n {output} --outdir {outdir} --nomodel --extsize 2000 -q 0.01 --llocal 10000000".format(input=macs2_input_filtered, outdir=macs2_outdir, output=macs2_sample)
    subprocess.check_output(macs2_cmd, shell=True, stderr=subprocess.STDOUT, universal_newlines=True)
