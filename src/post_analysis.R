# Title     : TODO
# Objective : TODO
# Created by: s215v
# Created on: 1/25/21



    #pdf("test.pdf", width=10, height=10)
    #random_peaks_df2offtargets = generate_background(sample_peaks_df2offtargets, mm9, n=50)
    #random_peaks_df2offtargets = peaks2offtargets_identity(random_peaks_df2offtargets, targets_df, genome_mm9)
    #x = dplyr::bind_rows(
    #  random_peaks_df2offtargets %>% dplyr::mutate(type="Background", peak_width=paste0(round((peak_end-peak_start)/1e5)/10, "M")) %>% dplyr::select(type, peak_width, bait_alignment_identity),
    #  sample_peaks_df2offtargets %>% dplyr::mutate(type="Signal", peak_width=paste0(round((peak_end-peak_start)/1e5)/10, "M")) %>% dplyr::select(type, peak_width, bait_alignment_identity)
    #) %>% dplyr::filter(bait_alignment_identity<100) %>%
    #  dplyr::group_by(peak_width) %>%
    #  dplyr::mutate(peak_width=paste(peak_width, n())) %>%
    #  dplyr::ungroup()
    #ggplot(x) +
    #  ggridges::geom_density_ridges(aes(fill=type, x=bait_alignment_identity, y=peak_width), alpha=0.8)
    #dev.off()


  #
  # Find breaks close to bait
  #
  #offtargets_all_df %>%
  #  dplyr::filter(bait_chrom==offtarget_chrom) %>%
  #  dplyr::mutate(offtarget_name=stringr::str_glue("{bait}_{mismatches}_{seq}", bait=bait_name, mismatches=offtarget_mismatches, seq=offtarget_sequence), offtarget_score=0) %>%
  #  dplyr::select(offtarget_chrom, offtarget_start, offtarget_end, offtarget_name, offtarget_score, offtarget_strand) %>%
  #  readr::write_tsv("data/macs2/offtargets.bed", col_names=F)
  #junctions_ranges = GenomicRanges::makeGRangesFromDataFrame(junctions_df %>% dplyr::select(-dplyr::matches("is_near_bait|is_very_near_bait|is_very_near_offtarget")) %>% dplyr::mutate(start=junction_start, end=junction_end, seqnames=junction_chrom), keep.extra.columns=T)
  #offtargets_ranges_near = GenomicRanges::makeGRangesFromDataFrame(offtargets_df %>% dplyr::filter(offtarget_mismatches<=3) %>% dplyr::mutate(offtarget_bait_chrom=bait_chrom, start=offtarget_start-2e6, end=offtarget_end+2e6, seqnames=offtarget_chrom), keep.extra.columns=T)
  #offtarget2junctions.map = as.data.frame(IRanges::mergeByOverlaps(offtargets_ranges_near, junctions_ranges)) %>%
  #  dplyr::filter(junction_chrom==offtarget_chrom & bait_chrom==offtarget_bait_chrom) %>%
  #  dplyr::mutate(is_near_bait=ifelse(offtarget_mismatches==0 & offtarget_start-2e6<=junction_start & junction_start<=offtarget_end+2e6, offtarget_chrom, "")) %>%
  #  dplyr::mutate(is_very_near_bait=ifelse(offtarget_mismatches==0 & offtarget_start-1e6<=junction_start & junction_start<=offtarget_end+1e6, offtarget_chrom, "")) %>%
  #  dplyr::mutate(is_very_near_offtarget=ifelse(offtarget_start-1e6<=junction_start & junction_start<=offtarget_end+1e6, offtarget_chrom, "")) %>%
  #  dplyr::group_by(junction_id) %>%
  #  dplyr::summarise(
  #    is_near_bait=paste0(unique(is_near_bait[is_near_bait!=""]), collapse=";"),
  #    is_very_near_bait=paste0(unique(is_very_near_bait[is_very_near_bait!=""]), collapse=";"),
  #    is_very_near_offtarget=paste0(unique(is_very_near_offtarget[is_very_near_offtarget!=""]), collapse=";"))
  #junctions_df = junctions_df %>%
  #  dplyr::select(-dplyr::matches("is_near_bait|is_very_near_bait|is_very_near_offtarget")) %>%
  #  dplyr::left_join(offtarget2junctions.map, by="junction_id") %>%
  #  dplyr::mutate(
  #    is_near_bait=ifelse(is_near_bait!="", is_near_bait, NA_character_),
  #    is_very_near_bait=ifelse(is_very_near_bait!="", is_very_near_bait, NA_character_),
  #    is_very_near_offtarget=ifelse(is_very_near_offtarget!="", is_very_near_offtarget, NA_character_))

calculate_bait_enrichment = function(peaks_df, junctions_df) {
  peaks_ranges = GenomicRanges::makeGRangesFromDataFrame(peaks_df %>% dplyr::mutate(seqnames=peak_chrom, start=peak_start, end=peak_end), keep.extra.columns=T)
  breaks_ranges = GenomicRanges::makeGRangesFromDataFrame(junctions_df %>% dplyr::mutate(seqnames=break_chrom, start=break_start, end=break_end), keep.extra.columns=T)
  peaks2breaks.map = as.data.frame(IRanges::mergeByOverlaps(breaks_ranges, peaks_ranges)) %>%
    dplyr::select(peak_chrom, peak_start, peak_end, break_chrom, break_start, break_end, break_id, peak_id) %>%
    dplyr::filter(break_chrom==peak_chrom) %>%
    dplyr::distinct(break_id, peak_id)

  peak_enrichment = junctions_df %>%
    dplyr::select(-dplyr::matches("^peak_id$")) %>%
    dplyr::left_join(peaks2breaks.map, by="break_id") %>%
    dplyr::left_join(peaks_df, by="peak_id") %>%
    dplyr::filter(!is.na(peak_id)) %>%
    dplyr::mutate(breaks_n=n()) %>%
    dplyr::group_by(break_bait_chrom) %>%
    dplyr::mutate(total_bait_breaks_n=n()) %>%
    dplyr::group_by(peak_id) %>%
    dplyr::mutate(total_peak_breaks_n=n()) %>%
    dplyr::group_by(break_bait_chrom) %>%
    dplyr::mutate(total_bait_breaks_n=n()) %>%
    dplyr::group_by(peak_id, break_bait_chrom) %>%
    dplyr::do((function(y){
      # y = z_clustered_peaks %>% dplyr::filter(break_bait_chrom=="chr1" & peak_id==424)
      yy<<-y
      x11 = nrow(y)
      x12 = y$total_peak_breaks_n[1] - x11
      x21 = y$total_bait_breaks_n[1] - x11
      x22 = y$breaks_n[1] - y$total_bait_breaks_n[1] - x12
      peak_bait_frac = x11 / y$total_peak_breaks_n[1]
      x = matrix(c(x11,x12,x21,x22), ncol=2)
      test = fisher.test(x)
      data.frame(
        peak_id=y$peak_id[1],
        is_near_bait=paste(na.omit(unique(y$is_near_bait)), collapse=""),
        is_very_near_bait=paste(na.omit(unique(y$is_very_near_bait)), collapse=""),
        is_very_near_offtarget=paste(na.omit(unique(y$is_very_near_offtarget)), collapse=""),
        break_bait_chrom=y$break_bait_chrom[1],
        peak_chrom=y$break_chrom[1],
        peak_start=y$peak_start[1],
        peak_end=y$peak_end[1],
        peak_bait_breaks=x11,
        total_peak_breaks=y$total_peak_breaks_n[1],
        peak_bait_frac=peak_bait_frac,
        pvalue=test$p.value,
        odds=test$estimate)
    })(.)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(break_bait_chrom_number = gsub("chr", "", break_bait_chrom)) %>%
    dplyr::mutate(break_bait_chrom_number=dplyr::case_when(
      toupper(break_bait_chrom_number)=="X"~"21",
      toupper(break_bait_chrom_number)=="Y"~"22",
      T~break_bait_chrom_number)) %>%
    dplyr::mutate(break_bait_chrom_number=as.numeric(break_bait_chrom_number))



  peak_enrichment
}


calculate_coverage = function(ranges, mm9, extend) {
  ranges.extended = GenomicRanges::resize(ranges, width=extend, fix="center")
  ranges.extended = GenomicRanges::restrict(ranges.extended, start=0, end=setNames(seqlengths(mm9), seqnames(mm9)))
  coverage_df = as.data.frame(ranges.extended) %>%
    dplyr::group_by(break_bait_chrom) %>%
    dplyr::do((function(z){
      zz<<-z
      z_ranges = GenomicRanges::makeGRangesFromDataFrame(z, keep.extra.columns=T)
      as.data.frame(as(GenomicRanges::coverage(z_ranges), "GRanges")) %>%
        tidyr::crossing(z %>% dplyr::select(break_bait_chrom, scale_factor) %>% dplyr::slice(1)) %>%
        dplyr::select(seqnames, start, end, score, scale_factor)
    })(.)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(score_norm=score*scale_factor) %>%
    dplyr::select(-scale_factor)

  # tempfile()
  total_path  = "data/macs2/total.bdg"
  chr_tmp_path  = "data/macs2/chr.bdg"
  # chr_tmp_path = tempfile()
  # total_path = tempfile()
  i = 0
  for(chr in unique(coverage_df$break_bait_chrom)) {
    # chr = "chr3"
    coverage_df.chr = coverage_df %>% dplyr::filter(break_bait_chrom==chr)

    if(i==0) {
      print("Write total")
      readr::write_tsv(coverage_df.chr %>% dplyr::select(seqnames, start, end, score_norm), total_path, col_names=F)
    } else {
      print(paste0("Write chr: ", chr))
      readr::write_tsv(coverage_df.chr %>% dplyr::mutate(score_norm=-score_norm) %>% dplyr::select(seqnames, start, end, score_norm), chr_tmp_path, col_names=F)
      system(stringr::str_glue("macs2 bdgcmp -m subtract -t {total} -c {chr} -o {total}", total=total_path, chr=chr_tmp_path))
    }

    i = i + 1
  }

  bedgraph_cols = cols(seqnames=col_character(), start=col_double(), end=col_double(), coverage=col_double())
  coverage_norm_df = readr::read_tsv(total_path, col_names=names(bedgraph_cols$cols), col_types=bedgraph_cols) %>%
    dplyr::mutate(coverage_id=1:n())

  coverage_norm_df
}