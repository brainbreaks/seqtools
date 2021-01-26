library(VennDiagram)
library(dplyr)
library(RColorBrewer)
library(IRanges)
library(GenomicRanges)

macs_call_islands = function(signal_df, background_df, matching_tiles=T, minqvalue, maxgap, minlen) {
  qvalue_path = "data/macs2/bdgcmp_qvalue.bdg"
  islands_path = "data/macs2/bdgcmp_islands.bed"

  # Call islands
  if(matching_tiles) {
    qvalues_df = background_df %>%
      dplyr::select(seqnames, start, end, coverage.control=coverage) %>%
      dplyr::inner_join(signal_df %>% dplyr::select(seqnames, start, end, coverage.sample=coverage), by=c("seqnames", "start", "end")) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(pvalue=pgamma(coverage.sample, shape=coverage.control, rate=1, lower.tail=F)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(qvalue_score=qvalue::qvalue(pvalue)$qvalues) %>%
      dplyr::mutate(qvalue_score=ifelse(qvalue_score==0, 315, -log10(qvalue_score))) %>% # TODO: change 315 to something better
      dplyr::select(qvalue_chrom=seqnames, qvalue_start=start, qvalue_end=end, qvalue_score) %>%
      dplyr::mutate(qvalue_id=1:n())
    readr::write_tsv(qvalues_df, qvalue_path, col_names=F)
  } else {
    sample_bdg_path = "data/macs2/bdgcmp_sample.bdg"
    baseline_bdg_path = "data/macs2/bdgcmp_baseline.bdg"
    readr::write_tsv(signal_df %>% dplyr::select(seqnames, start, end, coverage), sample_bdg_path, col_names=F)
    readr::write_tsv(background_df %>% dplyr::select(seqnames, start, end, coverage), baseline_bdg_path, col_names=F)
    system(stringr::str_glue("macs2 bdgcmp -t {input} -c {noise} -m qpois -o {output}", input=sample_bdg_path, noise=baseline_bdg_path, output=qvalue_path))

    qvalue_cols = cols(qvalue_chrom=col_character(), qvalue_start=col_double(), qvalue_end=col_double(), qvalue_score=col_double())
    qvalues_df = readr::read_tsv(qvalue_path, col_names=names(qvalue_cols$cols), col_types=qvalue_cols) %>%
      dplyr::mutate(qvalue_id=1:n())
  }

  system(stringr::str_glue("macs2 bdgpeakcall -i {qvalue} -c {cutoff} --min-length {format(island_min_length, scientific=F)} --max-gap {format(island_max_gap, scientific=F)} -o {output}",
       qvalue=qvalue_path, output=islands_path, cutoff=minqvalue, island_max_gap=maxgap, island_min_length=minlen))

  # Read results file
  islands_cols = cols(
    island_chrom=col_character(), island_start=col_double(), island_end=col_double(), island_length=col_character(), island_summit_abs=col_double(),
    island_score=col_character(), island_fc=col_double(), island_pvalue_log10=col_double(), island_qvalue_log10=col_double(), island_sammit_offset=col_double()
  )
  islands_df = readr::read_tsv(islands_path, skip=1, col_names=names(islands_cols$cols), col_types=islands_cols) %>%
    dplyr::mutate(island_summit_pos=island_start + island_sammit_offset) %>%
    dplyr::select(island_chrom, island_start, island_end, island_summit_pos)

  if(nrow(islands_df)>0) {
    islands_df = islands_df %>% dplyr::mutate(island_id=1:n())
    qvalues_ranges = GenomicRanges::makeGRangesFromDataFrame(qvalues_df %>% dplyr::mutate(end=qvalue_end, start=qvalue_start, seqnames=qvalue_chrom), keep.extra.columns=T)
    summit_ranges = GenomicRanges::makeGRangesFromDataFrame(islands_df %>% dplyr::mutate(end=island_summit_pos, start=island_summit_pos, seqnames=island_chrom), keep.extra.columns=T)
    sample_ranges = GenomicRanges::makeGRangesFromDataFrame(signal_df %>% dplyr::mutate(junction_chrom=seqnames, junction_start=start, junction_end=end, junction_coverage=coverage), keep.extra.columns=T)
    summit2qvalue.map = as.data.frame(IRanges::mergeByOverlaps(qvalues_ranges, summit_ranges)) %>%
      dplyr::filter(qvalue_chrom==island_chrom) %>%
      dplyr::group_by(island_id) %>%
      dplyr::summarise(island_summit_qvalue=max(qvalue_score), .groups="keep")
    summit2sample.map = as.data.frame(IRanges::mergeByOverlaps(sample_ranges, summit_ranges)) %>%
      dplyr::group_by(island_id) %>%
      dplyr::summarise(island_summit_abs=max(junction_coverage), .groups="keep")

    islands_df = islands_df %>%
      dplyr::select(-dplyr::matches("island_summit_abs|island_summit_qvalue")) %>%
      dplyr::select(-dplyr::matches("island_summit_qvalue")) %>%
      dplyr::inner_join(summit2qvalue.map, by=c("island_id")) %>%
      dplyr::select(-dplyr::matches("island_summit_abs")) %>%
      dplyr::inner_join(summit2sample.map, by=c("island_id"))
  } else {
    islands_df$island_id = c()
    islands_df$island_summit_qvalue = c()
    islands_df$island_summit_abs = c()
    islands_df$island_summit_pos = c()
  }

  list(islands=islands_df, qvalues=qvalues_df)
}