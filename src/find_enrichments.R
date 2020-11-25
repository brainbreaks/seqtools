library(readr)
library(dplyr)
library(stringr)

bait_region_size = 4e6

# breaksites_df = readr::read_tsv("data/breaksites_primers.tsv")
offtargets_df = readr::read_tsv("data/offtargets.tsv") %>% # Doesn't have position for Chr6
  dplyr::filter(offtarget_mismatches<5)
table(offtargets_df$bait_chrom, offtargets_df$offtarget_mismatches)

breaksites_df = offtargets_df %>% # Doesn't have position for Chr6
  dplyr::filter(offtarget_mismatches<5) %>%
  dplyr::select(bait_chrom, bait_position=offtarget_start)

macs2_outdir = "data/macs2"
macs2_cols = c("junction_chrom", "junction_start", "junction_end", "junction_name", "junction_score", "junction_strand")
for(macs2_input in list.files("data/breaks", pattern="*.bed", full.names=T)) {
    bait_chrom = gsub("Chr", "chr", (gsub("_APH_no10kb_Merge\\.bed$", "", basename(macs2_input))))
    macs2_input_filtered = paste0("tmp/", basename(gsub("\\.bed$", "_filtered.bed", macs2_input)))

    writeLines(paste0("Processing ", bait_chrom, "..."))

    macs2_input_df = readr::read_tsv(macs2_input, col_names=macs2_cols) %>%
      dplyr::mutate(bait_chrom=bait_chrom) %>%
      dplyr::inner_join(breaksites_df, by="bait_chrom") %>%
      dplyr::filter(!dplyr::between(junction_start, bait_position - bait_region_size/2, bait_position + bait_region_size/2)) %>%
      dplyr::select(dplyr::starts_with("junction_")) %>%
      dplyr::distinct(junction_chrom, junction_start, .keep_all=T)

    readr::write_tsv(macs2_input_df[macs2_cols], file=macs2_input_filtered, col_names=F)

    macs2_cmd = stringr::str_glue(
      "macs2 callpeak -t {input} -f BED -g mm --keep-dup all -n {output} --outdir {outdir} --nomodel --extsize 2000 -q 0.01 --llocal 10000000",
      input=macs2_input_filtered, outdir=macs2_outdir, output=macs2_sample)
    system(macs2_cmd)

    readr::read_tsv(paste0(macs2_outdir, "/", macs2_sample, "_summits.bed"), col_names=F)
}