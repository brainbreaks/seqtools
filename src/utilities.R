library(VennDiagram)
library(dplyr)
library(RColorBrewer)
library(IRanges)
library(GenomicRanges)

venn_ranges = function(r1, r2, name1, name2) {
  r1$r1_id = 1:length(r1)
  r1_df = as.data.frame(r1)
  r2$r2_id = 1:length(r2)
  r2_df = as.data.frame(r2)
  r1_r2.map = as.data.frame(IRanges::mergeByOverlaps(r1, r2)) %>% dplyr::select(r1_id, r2_id)

  r1.names = paste("R1", r1_df %>% dplyr::anti_join(r1_r2.map, by="r1_id") %>% .$r1_id)
  r2.names = paste("R2", r2_df %>% dplyr::anti_join(r1_r2.map, by="r2_id") %>% .$r2_id)
  common.names = unique(paste("R1+R2", apply(r1_df %>% dplyr::inner_join(r1_r2.map, by="r1_id") %>% dplyr::select(r1_id, r2_id), 1, paste, collapse="_")))
  r1.names = unique(c(r1.names, common.names))
  r2.names = unique(c(r2.names, common.names))

  # Chart
  p = VennDiagram::venn.diagram(
    x = list(r1.names, r2.names),
    category.names = c(name1 , name2),
    lwd = 2,
    lty = 'blank',
    fill = RColorBrewer::brewer.pal(3, "Pastel2")[1:2],
    cat.cex = 2,
    cat.fontface = "bold",
    filename = NULL
  )

  grid::grid.draw(p)
}

macs_call_peaks = function(signal_df, background_df, matching_tiles=T, minqvalue, maxgap, minlen) {
  qvalue_path = "data/macs2/bdgcmp_qvalue.bdg"
  peaks_path = "data/macs2/bdgcmp_peaks.bed"

  # Call peaks
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

  system(stringr::str_glue("macs2 bdgpeakcall -i {qvalue} -c {cutoff} --min-length {format(peak_min_length, scientific=F)} --max-gap {format(peak_max_gap, scientific=F)} -o {output}",
       qvalue=qvalue_path, output=peaks_path, cutoff=minqvalue, peak_max_gap=maxgap, peak_min_length=minlen))

  # Read results file
  peaks_cols = cols(
    peak_chrom=col_character(), peak_start=col_double(), peak_end=col_double(), peak_length=col_character(), peak_summit_abs=col_double(),
    peak_score=col_character(), peak_fc=col_double(), peak_pvalue_log10=col_double(), peak_qvalue_log10=col_double(), peak_sammit_offset=col_double()
  )
  peaks_df = readr::read_tsv(peaks_path, skip=1, col_names=names(peaks_cols$cols), col_types=peaks_cols) %>%
    dplyr::mutate(peak_summit_pos=peak_start + peak_sammit_offset) %>%
    dplyr::select(peak_chrom, peak_start, peak_end, peak_summit_pos)

  if(nrow(peaks_df)>0) {
    peaks_df = peaks_df %>% dplyr::mutate(peak_id=1:n())
    qvalues_ranges = GenomicRanges::makeGRangesFromDataFrame(qvalues_df %>% dplyr::mutate(end=qvalue_end, start=qvalue_start, seqnames=qvalue_chrom), keep.extra.columns=T)
    summit_ranges = GenomicRanges::makeGRangesFromDataFrame(peaks_df %>% dplyr::mutate(end=peak_summit_pos, start=peak_summit_pos, seqnames=peak_chrom), keep.extra.columns=T)
    sample_ranges = GenomicRanges::makeGRangesFromDataFrame(signal_df %>% dplyr::mutate(break_chrom=seqnames, break_start=start, break_end=end, break_coverage=coverage), keep.extra.columns=T)
    summit2qvalue.map = as.data.frame(IRanges::mergeByOverlaps(qvalues_ranges, summit_ranges)) %>%
      dplyr::filter(qvalue_chrom==peak_chrom) %>%
      dplyr::group_by(peak_id) %>%
      dplyr::summarise(peak_summit_qvalue=max(qvalue_score))
    summit2sample.map = as.data.frame(IRanges::mergeByOverlaps(sample_ranges, summit_ranges)) %>%
      dplyr::filter(break_chrom==break_chrom) %>%
      dplyr::group_by(peak_id) %>%
      dplyr::summarise(peak_summit_abs=max(break_coverage))

    peaks_df = peaks_df %>%
      dplyr::select(-dplyr::matches("peak_summit_abs|peak_summit_qvalue")) %>%
      dplyr::select(-dplyr::matches("peak_summit_qvalue")) %>%
      dplyr::inner_join(summit2qvalue.map, by=c("peak_id")) %>%
      dplyr::select(-dplyr::matches("peak_summit_abs")) %>%
      dplyr::inner_join(summit2sample.map, by=c("peak_id"))
  } else {
    peaks_df$peak_id = c()
    peaks_df$peak_summit_qvalue = c()
    peaks_df$peak_summit_abs = c()
    peaks_df$peak_summit_pos = c()
  }

  list(peaks=peaks_df, qvalues=qvalues_df)
}