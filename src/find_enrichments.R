library(readr)
library(dplyr)
library(stringr)

bait_region_size = 4e6

# breaksites_df = readr::read_tsv("data/breaksites_primers.tsv")
offtargets_df = readr::read_tsv("data/offtargets_predicted.tsv") %>% # Doesn't have position for Chr6
  dplyr::filter(offtarget_mismatches<5)

breaksites_df = offtargets_df %>% # Doesn't have position for Chr6
  dplyr::filter(offtarget_mismatches<5) %>%
  dplyr::select(bait_chrom, bait_position=offtarget_start)

tempdir = "tmp"
enrichment_outdir = "data/enrichment"
junctions_cols = c("junction_chrom", "junction_start", "junction_end", "junction_name", "junction_score", "junction_strand")
peaks_cols = c("control_peak_chrom", "control_peak_start", "control_peak_end", "control_peak_length", "control_peak_summit", "control_peak_pileup", "control_peak_log10pvalue", "control_peak_fc", "control_peak_log10qvalue", "control_peak_name")
peaks_coltypes = cols(
  control_peak_chrom = col_character(),
  control_peak_start = col_double(),
  control_peak_end = col_double(),
  control_peak_length = col_double(),
  control_peak_summit = col_double(),
  control_peak_pileup = col_double(),
  control_peak_log10pvalue = col_double(),
  control_peak_fc = col_double(),
  control_peak_log10qvalue = col_double(),
  control_peak_name = col_character()
)

results = data.frame()
for(breaks_bed in list.files("data/breaks", pattern="*_APH_no10kb_Merge.bed", full.names=T)) {
    bait_chrom = gsub("(chr)([^_]+).*$", "\\L\\1\\2", basename(breaks_bed), perl=T, ignore.case=T)
    breaks_control_bed = gsub("_APH_no10kb_Merge.bed", "_DMSO.bed", breaks_bed)
    breaks_filtered_bed = paste0("tmp/", basename(gsub("\\.bed$", "_filtered.bed", breaks_bed)))
    writeLines(stringr::str_glue("Processing {bait}...", bait=bait_chrom))
    if(!file.exists(breaks_control_bed)) {
        writeLines("Control doesn't exist. Skipping...")
        next()
    }

    output_control_enrichment = paste0(tempdir, "/", bait_chrom, "_control")

    # Run enrichment detection on control
    writeLines(stringr::str_glue("Calculate enrichments for control: '{output}'...", output=output_control_enrichment))

    parser.add_argument('inputs', nargs='+', help='Input .bed files with detected breaks. Can also be multiple files or a whildcard expression (e.g.: path/to/*.bed) ')
    parser.add_argument('chromsizes', help='Chromosome sizes in tab separated format')
    parser.add_argument('output_path', help='Path to folder where all information needed to import to UCSC genome browser will be stored')
    parser.add_argument('--track-name', dest="track_name", help='Name of the UCSC track')
    parser.add_argument('-w|--window-size', dest="window_size", default=int(1e5), type=int, help='Window at which to agregate breaks number')
    parser.add_argument('-s|--window-step', dest="window_step", default=int(1e4), type=int, help='Step after each window')

    cmd_binarize = stringr::str_glue("venv/bin/python3 src/breaks_bed2wig.py ")
    system(cmd_binarize)

    cmd_enrichment_control = stringr::str_glue(
      "macs2 callpeak -t {input} -f BED -g mm --keep-dup all -n {output} --outdir {outdir} --nomodel --extsize 2000 -q 0.01 --llocal 10000000",
      input=breaks_control_bed, outdir=dirname(output_control_enrichment), output=basename(output_control_enrichment))
    system(cmd_enrichment_control)
    control_enrichment_df = readr::read_tsv(paste0(output_control_enrichment, "_peaks.xls"), col_names=names(peaks_coltypes$cols), comment="#", skip=22, col_types=peaks_coltypes) %>%
      dplyr::mutate(bait_chrom=bait_chrom)
    results = dplyr::bind_rows(results, control_enrichment_df)
    #
    #
    # macs2_input_df = readr::read_tsv(macs2_input, col_names=macs2_cols) %>%
    #   dplyr::mutate(bait_chrom=bait_chrom) %>%
    #   dplyr::inner_join(breaksites_df, by="bait_chrom") %>%
    #   dplyr::filter(!dplyr::between(junction_start, bait_position - bait_region_size/2, bait_position + bait_region_size/2)) %>%
    #   dplyr::select(dplyr::starts_with("junction_")) %>%
    #   dplyr::distinct(junction_chrom, junction_start, .keep_all=T)
    #
    # readr::write_tsv(macs2_input_df[macs2_cols], file=macs2_input_filtered, col_names=F)
    #
    # macs2_cmd = stringr::str_glue(
    #   "macs2 callpeak -t {input} -f BED -g mm --keep-dup all -n {output} --outdir {outdir} --nomodel --extsize 2000 -q 0.01 --llocal 10000000",
    #   input=macs2_input_filtered, outdir=macs2_outdir, output=macs2_sample)
    # system(macs2_cmd)
    #
    # readr::read_tsv(paste0(macs2_outdir, "/", macs2_sample, "_summits.bed"), col_names=F)
}


offtargets_ranges = with(offtargets_df, IRanges::IRanges(start=offtarget_start, end=offtarget_end))
control_junctions_ranges = with(results,  IRanges::IRanges(start=control_peak_start, end=control_peak_end))
validated_offtargets_overlapping = validated_offtargets %>%
dplyr::inner_join(as.data.frame(IRanges::mergeByOverlaps(offtargets_ranges, control_junctions_ranges)), by=c("offtarget_start"="control_junctions_ranges.start", "offtarget_end"="validated_offtargets_ranges.end")) %>%
dplyr::mutate(offtarget_validated=T) %>%
dplyr::select(bait_name, bait_chrom, offtarget_is_primary, offtarget_chrom, offtarget_strand, offtarget_start, offtarget_end, offtarget_sequence, offtarget_start, offtarget_end, primer_alignment_start=primers_targets_ranges.start, primer_alignment_end=primers_targets_ranges.end, offtarget_validated)






for(breaks_bed in list.files("data/breaks", pattern="*DMSO.bed", full.names=T)) {
    bait_chrom = gsub("(chr)([^_]+).*$", "\\L\\1\\2", basename(breaks_bed), perl=T, ignore.case=T)


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