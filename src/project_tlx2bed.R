library(readr)
library(dplyr)
library(IRanges)
library(tidyr)
library(ggplot2)
library(BSgenome)
library(baseline)
library(foreach)
library(doParallel)
library(qvalue)



ggplot_karyoplot = function(coverage_df=NULL, baseline_df=NULL, peaks_df=NULL, bait_enrichemnt_df=NULL, max_coverage="auto", type="point", qvalue_th=0.01) {
  ggplot.colors = c(enriched="#0000FF", Control="#000000", Sample="#FF0000", "Sample filtered"="#FF0000", Max="#0000FF", Coverage="#000000", RDC="#FFFF00")

  if(max_coverage=="auto") {
    if(!is.null(baseline_df)) { max_coverage = median(baseline_df$coverage_smooth)*12 } else {
      if(!is.null(baseline_df)) { max_coverage = median(coverage_df$coverage)*12 } else {
        max_coverage = 50
      }
    }
  }

  count_annotations = 1
  if(!is.null(peaks_df) & nrow(peaks_df)>0) {
    count_annotations = length(unique(peaks_df$break_exp_condition)) + 3
  }

  p = ggplot() +
    # geom_hex(aes(y=coverage, x=start), data=coverage_df %>% dplyr::filter(dplyr::between(coverage, 2, max_coverage) & seqnames %in% chr), bins=30) +
    #geom_rect(aes(xmin=peak_start, ymin=-5-break_bait_chrom_number, xmax=peak_end, ymax=-4-break_bait_chrom_number, fill="enriched"), data=peak_enrichment.significant %>% dplyr::mutate(seqnames=peak_chrom) %>% dplyr::filter(seqnames %in% chr)) +
    coord_cartesian(ylim=c(ifelse(is.null(bait_enrichemnt_df), -count_annotations*max_coverage*0.1, -25), max_coverage)) +
    scale_color_manual(values=ggplot.colors) +
    scale_fill_manual(values=ggplot.colors) +
    scale_x_continuous(breaks=scale_breaks) +
    facet_wrap(~seqnames, scales="free_x", ncol=1) +
    labs(x="")

  if(!is.null(coverage_df)) {
    coverage_df = coverage_df %>% dplyr::mutate(coverage_limited=ifelse(coverage>max_coverage, max_coverage, coverage))
    if(!("break_exp_condition" %in% colnames(coverage_df))) {
      coverage_df$break_exp_condition = "Coverage"
    }
    coverage_df$peaks_offset = 1
    if(!is.null(peaks_df) & nrow(peaks_df)>0) {
      coverage_df$peaks_offset = length(unique(coverage_df$break_exp_condition)) + 2
    }
    coverage_df$break_exp_condition_number = as.numeric(as.factor(as.character(coverage_df$break_exp_condition)))
    if(type=="line") {
      p = p + geom_line(aes(y=coverage_limited, x=start, color=break_exp_condition), data=coverage_df, alpha=0.5, size=0.5)
    } else {
      p = p + geom_point(aes(y=coverage_limited, x=start, color=break_exp_condition), data=coverage_df, alpha=0.01, size=1)
    }
    p = p + geom_rect(aes(xmin=start, xmax=end, fill=break_exp_condition, alpha=coverage, ymin=-(peaks_offset+1+break_exp_condition_number)*max_coverage*0.1, ymax=-(peaks_offset+break_exp_condition_number)*max_coverage*0.1), data=coverage_df)
  }
  if(!is.null(peaks_df) & nrow(peaks_df)>0) {

    peaks_df$break_exp_condition_number = as.numeric(as.factor(as.character(peaks_df$break_exp_condition)))
    p = p + geom_rect(aes(xmin=peak_start, xmax=peak_end, fill=break_exp_condition, ymin=-(0+break_exp_condition_number)*max_coverage*0.1, ymax=-(1+break_exp_condition_number)*max_coverage*0.1), data=peaks_df)
    if("peak_summit_abs" %in% colnames(peaks_df)) {
      # peaks_df$peak_summit_abs_display = ifelse(peaks_df$peak_summit_abs>=1e3, with(peaks_df, paste0(floor(peak_summit_abs/10^log10(peak_summit_abs)), "e", floor(log10(peak_summit_abs)))), as.character(round(peaks_df$peak_summit_abs)))
      peaks_df$peak_summit_abs_display = round(peaks_df$peak_summit_qvalue)
      p = p + ggrepel::geom_text_repel(aes(x=peak_start, y=-(0.5+break_exp_condition_number)*max_coverage*0.1, label=peak_summit_abs_display), color="#333333", data=peaks_df)
    }
  }
  if(!is.null(baseline_df)) {
    #baseline_df.sign = baseline_df %>% dplyr::filter(qvalue <= qvalue_th)
    #if(nrow(baseline_df.sign)>0) {
    #  p = p + geom_rect(aes(xmin=start, xmax=end, color="mypeak"), ymin=-3, ymax=-2, data=baseline_df.sign)
    #}

    if(!("break_exp_condition" %in% colnames(baseline_df))) { baseline_df$break_exp_condition = "Coverage" }
    baseline_df$break_exp_condition = baseline_df$break_exp_condition
    p = p + geom_line(aes(y=coverage_smooth, x=start, color=break_exp_condition), data=baseline_df, alpha=0.8)
    p = p + geom_line(aes(y=coverage_th01, x=start, color=break_exp_condition), data=baseline_df, linetype="dashed", alpha=0.8)
  }
  if(!is.null(bait_enrichemnt_df)) { p = p + geom_text(aes(x=peak_start, y=-5-break_bait_chrom_number, label=gsub("chr", "", break_bait_chrom)), data=bait_enrichemnt_df, size=3) }

  p
}

macs_call_peaks = function(signal_df, background_df, matching_tiles=T) {
  qvalue_path = "data/macs2/bdgcmp_qvalue.bdg"
  peaks_path = "data/macs2/bdgcmp_peaks.bed"

  # Call peaks
  qvalues_df = NULL
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
       qvalue=qvalue_path, output=peaks_path, cutoff=peaks_minqvalue, peak_max_gap=peaks_maxgap, peak_min_length=peaks_minlen))

  # Read results file
  peaks_cols = cols(
    peak_chrom=col_character(), peak_start=col_double(), peak_end=col_double(), peak_length=col_character(), peak_summit_abs=col_double(),
    peak_score=col_character(), peak_fc=col_double(), peak_pvalue_log10=col_double(), peak_qvalue_log10=col_double(), peak_sammit_offset=col_double()
  )
  peaks_df = readr::read_tsv(peaks_path, skip=1, col_names=names(peaks_cols$cols), col_types=peaks_cols) %>%
    dplyr::mutate(peak_summit_pos=peak_start + peak_sammit_offset) %>% dplyr::select(-peak_summit_abs)

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
  }

  list(peaks=peaks_df, qvalues=qvalues_df)
}

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

fit_baseline = function(coverage_df, binstep, llocal) {
  coverage_df.baseline = coverage_df %>%
    dplyr::mutate(binstep=binstep, llocal=llocal) %>%
    dplyr::mutate(break_exp_condition="Control") %>%
    dplyr::arrange(seqnames, start) %>%
    dplyr::group_by(break_exp_condition, seqnames) %>%
    #dplyr::mutate(coverage_mod=coverage)
    dplyr::mutate(is_q2=coverage>quantile(coverage, 0.2), is_q8=coverage<quantile(coverage, 0.8), is_inner_quantile=is_q2 & is_q8) %>%
    dplyr::mutate(coverage_median=median(coverage[is_inner_quantile], na.rm=T), coverage_sd=sd(coverage[is_inner_quantile], na.rm=T)) %>%
    dplyr::mutate(coverage_min=coverage_median-2*coverage_sd) %>%
    dplyr::mutate(coverage_mod=ifelse(coverage<coverage_min, NA_real_, coverage)) %>%
    dplyr::mutate(is_filled=is.na(coverage_mod), coverage_mod=zoo::na.fill(coverage_mod, "extend"))

  parallelCluster = parallel::makeCluster(19, type="PSOCK", methods=FALSE)
  doParallel::registerDoParallel(parallelCluster)
  coverage_dflist.baseline = coverage_df.baseline %>%
    dplyr::group_split()
  coverage_df.baseline_output = foreach::foreach(z=coverage_dflist.baseline, .combine=rbind) %dopar% {
    x = matrix(as.numeric(z$coverage_mod), nrow=1)
    colnames(x) = z$start
    bc.irls = baseline::baseline(x, method="medianWindow", hws=z$llocal[1]/z$binstep[1], hwm=z$llocal[1]/z$binstep[1]) # !!!!
    z.coverage_smooth = as.numeric(bc.irls@baseline)
    cbind(z[c("seqnames", "start", "end", "coverage", "coverage_mod", "is_filled")], coverage_smooth=z.coverage_smooth)
  }
  stopCluster(parallelCluster)

  class(coverage_df.baseline_output) = c("tbl_df", "tbl", "data.frame")
  coverage_df.baseline_output = coverage_df.baseline_output %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      pvalue=pgamma(coverage, shape=coverage_smooth, rate=1, lower.tail=F),
      coverage_th01=qgamma(0.01, coverage_smooth, 1, lower.tail=F)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(qvalue=qvalue::qvalue(pvalue)$qvalues)


  return(coverage_df.baseline_output)
}


scale_breaks = function(x) {
    breaks = seq(0, max(x), 1e7)
    names(breaks) = paste0(breaks/1e6, "Mb")
    breaks
}


read_breaks_experiments = function(path, chromosomes_map) {
    #
  # Read breaks
  #
  tlx_cols = cols(
    Qname=col_character(), JuncID=col_character(), Rname=col_character(), Junction=col_double(),
    Strand=col_character(), Rstart=col_double(), Rend=col_double(),
    B_Rname=col_character(), B_Rstart=col_double(), B_Rend=col_double(), B_Strand=col_double(),
    B_Qstart=col_double(), B_Qend=col_double(), Qstart=col_double(), Qend=col_double(), Qlen=col_double(),
    B_Cigar=col_character(), Cigar=col_character(), Seq=col_character(), J_Seq=col_character(), Barcode=col_logical(),
    unaligned=col_double(), baitonly=col_double(), uncut=col_double(), misprimed=col_double(), freqcut=col_double(),
    largegap=col_double(), mapqual=col_double(), breaksite=col_double(), sequential=col_double(), repeatseq=col_double(), duplicate=col_double(), Note=col_character()
  )

  # Rename Chr14 to Alt289 (Was Alt289 and only PW255 was Alt287)
  junctions_ann.excluded = "PW121|PW246|PW247|PW248|PW249|JK096|JK097|JK098|JK099|JK100|JK101"
  junctions_ann.excluded = "NO FILTER"
  junctions_ann = data.frame(break_file=list.files(path, full.names=T, recursive=T, pattern="*.tlx")) %>%
     dplyr::filter(!grepl(junctions_ann.excluded, break_file)) %>%
     dplyr::mutate(
       break_bait_chrom1=tolower(basename(dirname(break_file))),
       break_bait_chrom = chromosomes_map[tolower(basename(dirname(break_file)))],
       break_condition = ifelse(grepl("DMSO", break_file), "Control", "Sample"),
       break_tech_sample = gsub("_.*", "", basename(break_file)),
       break_bio_sample=gsub(".*_((Alt|R)[0-9]+).*", "\\1", basename(break_file), ignore.case=T, perl=T)
     ) %>%
    dplyr::group_by(break_bait_chrom, break_bait_chrom1, break_condition) %>%
    dplyr::mutate(break_replicate=1:n()) %>%
    dplyr::ungroup()


  junctions_df = data.frame()
  for(i in 1:nrow(junctions_ann)) {
    tlx_ann = junctions_ann[i,,drop=F]
    tlx = cbind(tlx_ann, readr::read_tsv(tlx_ann$break_file, col_types=tlx_cols))
    bed = tlx %>%
      dplyr::mutate(break_chrom=chromosomes_map[Rname]) %>%
      dplyr::mutate(break_score=".", break_strand=ifelse(Strand=="-1", "-", "+")) %>%
      dplyr::mutate(break_exp_condition=tlx_ann$break_condition, break_file=tlx_ann$break_file, break_bait_chrom=tlx_ann$break_bait_chrom, break_replicate=tlx_ann$break_replicate) %>%
      dplyr::select(break_bait_chrom, break_exp_condition, break_chrom=Rname, break_start=Junction, break_end=Junction, break_name=Qname, break_score, break_strand, break_file)
    junctions_df = rbind(junctions_df, bed)
  }

  junctions_df %>% dplyr::mutate(break_id=1:n())
}

calculate_coverage2 = function(ranges, mm9, extend=5e6, binsize=1000, binstep=100) {
  mm9_original_tiles = unlist(GenomicRanges::tileGenome(mm9, tilewidth=binstep))
  mm9_original_tiles$original_start = GenomicRanges::start(mm9_original_tiles)
  mm9_original_tiles$original_end = GenomicRanges::end(mm9_original_tiles)
  mm9_tiles = GenomicRanges::trim(GenomicRanges::resize(mm9_original_tiles, binsize, fix="center"))

  seqinfo(ranges) = mm9
  ranges.extended = GenomicRanges::trim(GenomicRanges::resize(ranges, width=extend, fix="center"))
  coverage_df = as.data.frame(ranges.extended) %>%
    dplyr::group_by(break_bait_chrom) %>%
    dplyr::do((function(z){
      zz<<-z
      z_ranges = GenomicRanges::makeGRangesFromDataFrame(z, keep.extra.columns=T)
      z_tiles_df = as.data.frame(mm9_tiles) %>%
        dplyr::mutate(coverage=GenomicRanges::countOverlaps(mm9_tiles, z_ranges)) %>%
        tidyr::crossing(z %>% dplyr::select(break_bait_chrom, scale_factor) %>% dplyr::slice(1)) %>%
        dplyr::select(seqnames, start=original_start, end=original_end, coverage, scale_factor)
      z_tiles_df
    })(.)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(coverage_norm=coverage*scale_factor) %>%
    dplyr::select(-scale_factor) %>%
    dplyr::group_by(seqnames, start, end) %>%
    dplyr::summarise(coverage=sum(coverage_norm))

  #write.table(coverage_df %>% dplyr::select(seqnames, start, end, coverage), col.names=F, row.names=F, quote=F, file="data/macs2/coverage.bdg")

  class(coverage_df) = c("tbl_df", "tbl", "data.frame")
  coverage_df
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

venn_ranges = function(r1, r2, name1, name2) {
  r1$r1_id = 1:length(r1)
  r1_df = as.data.frame(r1)
  r2$r2_id = 1:length(r2)
  r2_df = as.data.frame(r2)
  r1_r2.map = as.data.frame(IRanges::mergeByOverlaps(r1, r2)) %>% dplyr::select(r1_id, r2_id)

  r1.names = paste("R1", r1_df %>% dplyr::anti_join(r1_r2.map, by="r1_id") %>% .$r1_id)
  r2.names = paste("R2", r2_df %>% dplyr::anti_join(r1_r2.map, by="r2_id") %>% .$r2_id)
  common.names = paste("R1+R2", r1_df %>% dplyr::inner_join(r1_r2.map, by="r1_id") %>% .$r1)
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

call_peaks = function(sample_df, control_df=NULL, debug=F) {
  #binsize=1e4
  binsize = 1e5
  #binstep=5e3
  binstep = 1e4
  #extend=1e6
  extend = 1e5
  #llocal=5e6
  llocal = 2e6
  peaks_minqvalue=-log10(0.01)
  peaks_maxgap=5e5
  peaks_minlen=200
  debug=F

  do_compare = !is.null(control_df)

  sample_bed_path = "data/breakensembl/Sample.bed"
  control_bed_path = "data/breakensembl/Control.bed"


  sample_df = sample_df %>% dplyr::mutate(break_name=stringr::str_glue("{bait}_{name}", bait=break_bait_chrom, name=break_name))
  if(debug) {
    readr::write_tsv(sample_df %>% dplyr::select(break_chrom, break_start, break_end, break_name, break_score, break_strand), sample_bed_path, col_names=F)
  }

  if(!do_compare) {
    control_df = sample_df
    control_bed_path = sample_bed_path
  } else {
    control_df = control_df %>% dplyr::mutate(break_name=stringr::str_glue("{bait}_{name}", bait=break_bait_chrom, name=break_name))
    if(debug) {
      readr::write_tsv(control_df %>% dplyr::select(break_chrom, break_start, break_end, break_name, break_score, break_strand), control_bed_path, col_names=F)
    }
  }

  single_break_chr = "chr6"

  # slocal background
  options("UCSC.goldenPath.url"="https://hgdownload.cse.ucsc.edu/goldenPath")
  mm9 = GenomeInfoDb::Seqinfo(genome="mm9")[paste0("chr", 1:19)]
  #mm9 = as(tibble::as_tibble(mm9, rownames="s") %>% dplyr::mutate(seqlengths=seqlengths[s==single_break_chr]) %>% tibble::column_to_rownames("s"), "Seqinfo")
  mm9_size = sum(mm9@seqlengths)
  mm9_effective_size=0.74 * mm9_size
  mm9_tiles = unlist(tileGenome(mm9, tilewidth=binsize))
  mm9_tiles$tile_id = 1:length(mm9_tiles)


  #
  # Caclulate control baseline
  #
  control_ranges = GenomicRanges::makeGRangesFromDataFrame(control_df %>% dplyr::filter(break_chrom==break_bait_chrom), end.field="break_end", start.field="break_start", strand.field="break_strand", seqnames.field="break_chrom", keep.extra.columns=T)
  control_tiles_df.coverage = calculate_coverage2(control_ranges, mm9, extend=extend, binsize=binsize, binstep=binstep)
  control_df.baseline = fit_baseline(control_tiles_df.coverage, binstep=binstep, llocal=llocal)


  #
  # Caclulate sample baseline
  #
  sample_ranges = GenomicRanges::makeGRangesFromDataFrame(sample_df %>% dplyr::filter(break_chrom==break_bait_chrom), end.field="break_end", start.field="break_start", strand.field="break_strand", seqnames.field="break_chrom", keep.extra.columns=T)
  #sample_df.coverage = calculate_coverage(sample_ranges, mm9, extend)
  sample_tiles_df.coverage = calculate_coverage2(sample_ranges, mm9, extend=extend, binsize=binsize, binstep=binstep)
  sample_df.baseline = fit_baseline(sample_tiles_df.coverage, binstep=binstep, llocal=llocal)

  #save(sample_df.baseline, control_df.baseline, control_tiles_df.coverage, control_tiles_df.coverage, file="data3.rda")

  #
  # Adjust control baseline to sample library
  #
  control_df.baseline_adjusted = control_df.baseline %>%
    dplyr::select(seqnames, start, end, coverage_smooth.control=coverage_smooth, is_filled.control=is_filled) %>%
    dplyr::inner_join(sample_df.baseline %>% dplyr::select(seqnames, start, end, coverage, coverage_smooth.sample=coverage_smooth, is_filled.sample=is_filled), by=c("seqnames", "start", "end")) %>%
    dplyr::group_by(seqnames) %>%
    dplyr::mutate(
      qc=!is_filled.control & coverage_smooth.control<quantile(coverage_smooth.control[!is_filled.control], 0.5),
      qs=!is_filled.sample & coverage_smooth.sample<quantile(coverage_smooth.sample[!is_filled.sample], 0.5),
      f=mean((coverage_smooth.sample/coverage_smooth.control)[qc&qs]),
      coverage_smooth=f*coverage_smooth.control) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      pvalue=pgamma(coverage, shape=coverage_smooth, rate=1, lower.tail=F),
      coverage_th01=qgamma(0.01, coverage_smooth, 1, lower.tail=F)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(qvalue=qvalue::qvalue(pvalue)$qvalues) %>%
    dplyr::select(seqnames, start, end, coverage, coverage_smooth, pvalue, qvalue, coverage_th01)

  #
  # Max Control/Sample baseline
  #
  mixed_df.baseline = control_df.baseline_adjusted %>%
    dplyr::select(seqnames, start, end, coverage_smooth.control=coverage_smooth) %>%
    dplyr::inner_join(sample_df.baseline %>% dplyr::select(seqnames, start, end, coverage, coverage_smooth.sample=coverage_smooth), by=c("seqnames", "start", "end")) %>%
    dplyr::mutate(coverage_smooth=pmax(coverage_smooth.sample, coverage_smooth.control)) %>%
    dplyr::select(seqnames, start, end, coverage, coverage_smooth) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      pvalue=pgamma(coverage, shape=coverage_smooth, rate=1, lower.tail=F),
      coverage_th01=qgamma(0.01, coverage_smooth, 1, lower.tail=F)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(qvalue=qvalue::qvalue(pvalue)$qvalues) %>%
    dplyr::select(seqnames, start, end, coverage, coverage_smooth, pvalue, qvalue, coverage_th01)


  sample_peaks = macs_call_peaks(
    signal_df=sample_tiles_df.coverage %>% dplyr::select(seqnames, start, end, coverage),
    background_df=mixed_df.baseline %>% dplyr::select(seqnames, start, end, coverage=coverage_smooth))

  control_peaks = macs_call_peaks(
    signal_df=control_tiles_df.coverage %>% dplyr::select(seqnames, start, end, coverage),
    background_df=control_df.baseline %>% dplyr::select(seqnames, start, end, coverage=coverage_smooth))

  sample_peaks_ranges = GenomicRanges::makeGRangesFromDataFrame(sample_peaks$peaks %>% dplyr::select(seqnames=peak_chrom, start=peak_start, end=peak_end, peak_id.sample=peak_id), keep.extra.columns=T)
  control_peaks_ranges = GenomicRanges::makeGRangesFromDataFrame(control_peaks$peaks %>% dplyr::select(seqnames=peak_chrom, start=peak_start, end=peak_end, peak_id.control=peak_id), keep.extra.columns=T)
  control_qvalues_ranges = GenomicRanges::makeGRangesFromDataFrame(control_peaks$qvalues %>% dplyr::select(seqnames=qvalue_chrom, start=qvalue_start, end=qvalue_end, qvalue_score.control=qvalue_score), keep.extra.columns=T)
  sample2control_peaks.map = as.data.frame(IRanges::mergeByOverlaps(sample_peaks_ranges, control_qvalues_ranges)) %>% dplyr::select(peak_id.sample, qvalue_score.control)
  sample_peaks$peaks = sample_peaks$peaks %>%
    dplyr::select(-dplyr::matches("^qvalue_score.control$")) %>%
    dplyr::left_join(sample2control_peaks.map, by=c("peak_id"="peak_id.sample")) %>%
    dplyr::arrange(dplyr::desc(qvalue_score.control)) %>%
    dplyr::distinct(peak_id, .keep_all=T)



  #control_df.coverage = calculate_coverage(sample_ranges, mm9, extend)
  #readr::write_tsv(control_df.coverage %>% dplyr::select(seqnames, start, end, coverage), control_bdg_path, col_names=F)
  #control_ranges.coverage = GenomicRanges::makeGRangesFromDataFrame(control_df.coverage, keep.extra.columns=T)
  #coverage2tile.map = as.data.frame(GenomicRanges::findOverlaps(mm9_tiles, control_ranges.coverage, ignore.strand=TRUE))
  #control_tiles_df.coverage = as.data.frame(mm9_tiles) %>%
  #  dplyr::left_join(coverage2tile.map, by=c("tile_id"="queryHits")) %>%
  #  dplyr::left_join(control_df.coverage %>% dplyr::select(coverage_id, coverage), by=c("subjectHits"="coverage_id")) %>%
  #  dplyr::group_by(seqnames, start, end) %>%
  #  dplyr::summarise(coverage=tidyr::replace_na(max(coverage), 0))
  #readr::write_tsv(control_tiles_df.coverage %>% dplyr::select(seqnames, start, end, coverage), "data/macs2/Control_tile.bdg", col_names=F)

  #peak_enrichment = calculate_bait_enrichment(control_peaks$., control_df)
  #peak_enrichment.significant = peak_enrichment %>%
  #  dplyr::filter(pvalue<0.01 & odds>1) %>%
  #  dplyr::mutate(odds=ifelse(is.finite(odds), odds, max(odds)))

  if(debug) {
    pdf(file="reports/duo_baseline_extend1e5_smooth2e6_bin1e5_step1e4.pdf", width=15, height=8)
    for(chr in paste0("chr", 1:19)) {
      print(chr)
      coverage_df = dplyr::bind_rows(control_tiles_df.coverage %>% dplyr::mutate(break_exp_condition="Control"), sample_tiles_df.coverage %>% dplyr::mutate(break_exp_condition="Sample")) %>%
         dplyr::filter(seqnames %in% chr)
      baseline_df = dplyr::bind_rows(
        control_df.baseline_adjusted %>% dplyr::mutate(break_exp_condition="Control"),
        sample_df.baseline %>% dplyr::mutate(break_exp_condition="Sample"),
        mixed_df.baseline %>% dplyr::mutate(break_exp_condition="Max")) %>%
         dplyr::filter(seqnames %in% chr)
      peaks_df = dplyr::bind_rows(
        control_peaks$peaks %>% dplyr::mutate(break_exp_condition="Control"),
        sample_peaks$peaks %>% dplyr::mutate(break_exp_condition="Sample"),
        sample_peaks$peaks %>% dplyr::filter(qvalue_score.control < peak_summit_qvalue) %>%  dplyr::mutate(break_exp_condition="Sample filtered"),
        rdc %>% dplyr::mutate(break_exp_condition="RDC", peak_start=rdc_start, peak_end=rdc_end, peak_chrom=rdc_chrom)) %>%
         dplyr::mutate(seqnames=peak_chrom) %>%
         dplyr::filter(seqnames %in% chr)
      p = ggplot_karyoplot(
        coverage_df=coverage_df,
        baseline_df=baseline_df,
        peaks_df=peaks_df,
        type="line"
      )
      #p = p + geom_vline(xintercept=c(64e6, 107e6, 117e6))
      #p = p + geom_vline(xintercept=c(6e6, 65.5e6, 97.5e6, 130e6, 135e6, 141e6)) #chr4
      print(p)
    }
    dev.off()
  }

  genome_mm9 = Biostrings::readDNAStringSet("data/mm9/mm9.fa.gz", "fasta")

  targets_df1 = readr::read_tsv("data/breaksites.tsv") %>%
    dplyr::inner_join(targets_df %>% dplyr::select(bait_chrom, offtarget_strand), by=c("breaksite_chrom"="bait_chrom")) %>%
    dplyr::select(bait_chrom=breaksite_chrom, breaksite_pos, offtarget_strand)
  targets_df1$offtarget_primer_sequence = as.character(Biostrings::getSeq(genome_mm9, GenomicRanges::makeGRangesFromDataFrame(
    targets_df1 %>% dplyr::mutate(seqnames=bait_chrom, start=breaksite_pos, end=breaksite_pos+20, strand=offtarget_strand)
  )))

  targets_df1 %>%
    dplyr::inner_join(targets_df %>% dplyr::select(bait_chrom, offtarget_primer_sequence), by="bait_chrom")


  #x = data.frame(
  #  seqnames=c("chr1", "chr11", "chr15", "Chr6", "chr8", "chr9"),
  #  start=c(41964577, 38361129, 61819136, 70900097, 61956495, 25942295),
  #  end=c(41964599, 38361151, 61819158, 70900119, 61956517, 25942317),
  #  strand=c("-","-","-", "-", "-", "-"))
  #x$sequence = as.character(Biostrings::getSeq(genome_mm9, GenomicRanges::makeGRangesFromDataFrame(x)))


  peaks_filtered_df = sample_peaks$peaks
  # %>% dplyr::filter(qvalue_score.control < peak_summit_qvalue)
  peaks_sequences = Biostrings::getSeq(genome_mm9, GenomicRanges::makeGRangesFromDataFrame(peaks_filtered_df %>% dplyr::mutate(seqnames=peak_chrom, start=peak_start, end=peak_end), keep.extra.columns=T))
  peaks_filtered_df$peak_sequence = as.character(peaks_sequences)
  peaks_filtered_df$peak_sequence_rev = as.character(Biostrings::reverseComplement(peaks_sequences))
  peaks_filtered_df = peaks_filtered_df %>%
    dplyr::select(-dplyr::matches("^offtarget_")) %>%
    dplyr::inner_join(targets_df %>% dplyr::select(bait_chrom, offtarget_primer_sequence), by=c("peak_chrom"="bait_chrom"))
  alignments = Biostrings::pairwiseAlignment(peaks_filtered_df$offtarget_primer_sequence, peaks_filtered_df$peak_sequence, type="global-local", gapOpening=2, gapExtension=1)
  alignments_rev = Biostrings::pairwiseAlignment(peaks_filtered_df$offtarget_primer_sequence, peaks_filtered_df$peak_sequence_rev, type="global-local", gapOpening=2, gapExtension=1)
  peaks_filtered_df$bait_sequence=as.character(Biostrings::pattern(alignments))
  peaks_filtered_df$bait_alignment_sequence=as.character(Biostrings::subject(alignments))
  peaks_filtered_df$bait_alignment_score=pmax(Biostrings::score(alignments), Biostrings::score(alignments_rev))
  peaks_filtered_df$bait_alignment_identity=pmax(Biostrings::pid(alignments), Biostrings::pid(alignments_rev))

  peaks_filtered_df %>%
    dplyr::group_by(peak_chrom) %>%
    dplyr::summarise(identity.90=sum(bait_alignment_identity>=90), identity.max=max(bait_alignment_identity))

  #peaks_filtered_df %>%
  #  dplyr::filter(bait_alignment_identity>80) %>%
  #  dplyr::select(offtarget_primer_sequence, bait_sequence, bait_alignment_sequence, bait_alignment_identity, bait_alignment_score) %>%
  #  View()


  pdf("test.pdf", width=8, height=8)
  hist(peaks_filtered_df$bait_alignment_identity)
  plot(sort(peaks_filtered_df$bait_alignment_identity))

  peaks_filtered_df %>%
    dplyr::select(peak_chrom, bait_sequence, bait_alignment_sequence, bait_alignment_identity, bait_alignment_score)  %>%
    ggplot() +
    geom_point(aes(bait_alignment_identity, bait_alignment_score))
    #ggrepel::geom_text_repel(aes(bait_alignment_identity, bait_alignment_score, label=peak_chrom))
  dev.off()


  writeLines(paste0(">", peaks_filtered_df$peak_id, "\n", peaks_filtered_df$peak_sequence), con="peaks.fasta")
  writeLines(paste0(">", offtargets_df$offtarget_id, "\n", offtargets_df$offtarget_sequence), con="targets.fasta")
  system("bwa index -b 100000000 peaks.fasta")
  system("bwa mem -t 30 -k 5 -C -a -T 10 -B 1 -O 100 -E 100 -L 100 -D 0 -c 70000 -y 60000 peaks.fasta targets.fasta > targets.sam")


  targets_df$offtarget_sequence

  offtargets_df %>% dplyr::filter(offtarget_mismatches==0) %>% dplyr::select(bait_chrom, offtarget_start, offtarget_end)
  GenomicRanges::makeGRangesFromDataFrame(as.data.frame(z) %>% dplyr::mutate(seqnames=peak_chrom, start=offtarget_start, end=offtarget_end))

  sample_peaks$peaks %>%
    dplyr::inner_join(offtargets_df %>% dplyr::filter(offtarget_mismatches==0) %>% dplyr::select(bait_chrom, offtarget_start, offtarget_end), by=c("peak_chrom"="bait_chrom")) %>%
    dplyr::rowwise() %>%
    dplyr::do((function(z){
      zz<<-z
      asddsa()
      z.range = GenomicRanges::makeGRangesFromDataFrame(as.data.frame(z) %>% dplyr::mutate(seqnames=peak_chrom, start=peak_start, end=peak_end))
      z.offtarget_range = GenomicRanges::makeGRangesFromDataFrame(as.data.frame(z) %>% dplyr::mutate(seqnames=peak_chrom, start=offtarget_start, end=offtarget_end))
      z.seq = as.character(Biostrings::getSeq(genome_mm9, z.range))
      z.bait_seq = as.character(Biostrings::getSeq(genome_mm9, z.offtarget_range))

      z.seq_sub = substr(z.seq, nchar(z.seq)/2-4e5, nchar(z.seq)/2+4e5)
      writeLines(z.bait_seq, con="seq1.fasta")
      writeLines(z.seq_sub, con="seq2.fasta")

      alignment = pairwiseAlignment(pattern=c(z.seq_sub, z.bait_seq), subject="supersede", type="local")
    })(.))
    dplyr::mutate(sequence=as.character(Biostrings::getSeq(genome_mm9, GRanges(offtarget_chrom, IRanges(start=offtarget_and_pam_start, end=offtarget_and_pam_end), strand=offtarget_strand))))


  jpeg("reports/rdc2peaks_venn.jpg", width=800, height=800)
  peaks_filtered_df = sample_peaks$peaks %>% dplyr::filter(qvalue_score.control < peak_summit_qvalue) %>%  dplyr::mutate(break_exp_condition="Sample filtered")
  rdc_ranges = GenomicRanges::makeGRangesFromDataFrame(rdc %>% dplyr::mutate(seqnames=rdc_chrom, start=rdc_start, end=rdc_end), keep.extra.columns=T)
  peaks_filtered_ranges = GenomicRanges::makeGRangesFromDataFrame(peaks_filtered_df %>% dplyr::mutate(seqnames=peak_chrom, start=peak_start, end=peak_end), keep.extra.columns=T)
  venn_ranges(rdc_ranges, peaks_filtered_ranges, name1="RDC", name2="Peak")
  dev.off()


  if(F) {
    pdf(file="reports/raw_control_vs_sample.pdf", width=20, height=8)
    chr = paste0("chr", 19)
    d = dplyr::bind_rows(control_tiles_df.coverage %>% dplyr::mutate(break_exp_condition="Control"), control_tiles_df.coverage %>% dplyr::mutate(break_exp_condition="Sample")) %>%
      dplyr::mutate(coverage=ifelse(coverage<50, coverage, NA_real_))
    b = dplyr::bind_rows(control_df.baseline %>% dplyr::mutate(break_exp_condition="Control"), sample_df.baseline %>% dplyr::mutate(break_exp_condition="Sample")) %>%
      dplyr::mutate(coverage_smooth=ifelse(coverage_smooth<50, coverage_smooth, NA_real_))
    p = ggplot(d %>% dplyr::filter(seqnames %in% chr)) +
      ggridges::geom_ridgeline(aes(y=break_exp_condition, x=start, height=coverage, fill=break_exp_condition), scale=0.1) +
      ggridges::geom_ridgeline(aes(y=break_exp_condition, x=start, height=coverage_smooth), scale=0.1, data=b %>% dplyr::filter(seqnames %in% chr), color="#FF0000", alpha=0.2, size=2) +
      scale_x_continuous(breaks=scale_breaks) +
      facet_wrap(~seqnames, scales="free_x", ncol=1) +
      theme_bw(base_size=20)
    print(p)
    dev.off()
  }


  #
  # VENN diagram overlap between SAMPLE and CONTROL
  #
  if(debug) {
    control_peaks_ranges = GenomicRanges::makeGRangesFromDataFrame(control_peaks$peaks, end.field="peak_end", start.field="peak_start", seqnames.field="peak_chrom", keep.extra.columns=T)
    sample_peaks_ranges = GenomicRanges::makeGRangesFromDataFrame(sample_peaks$peaks, end.field="peak_end", start.field="peak_start", seqnames.field="peak_chrom", keep.extra.columns=T)
    sample2control.map = as.data.frame(GenomicRanges::findOverlaps(control_peaks_ranges, sample_peaks_ranges))
    x = sample_peaks$peaks %>%
      dplyr::rename(peak_id.sample="peak_id") %>%
      dplyr::full_join(sample2control.map %>% dplyr::rename(peak_id.control="queryHits"), by=c("peak_id.sample"="subjectHits")) %>%
      dplyr::full_join(control_peaks$peaks, by=c("peak_id.control"="peak_id")) %>%
      setNames(gsub("\\.x$", ".sample", colnames(.))) %>%
      setNames(gsub("\\.y$", ".control", colnames(.)))

    names.common_control = unique(x %>% dplyr::filter(!is.na(peak_id.control) & !is.na(peak_id.sample)) %>% .$peak_id.control)
    names.common_sample = unique(x %>% dplyr::filter(!is.na(peak_id.control) & !is.na(peak_id.sample)) %>% .$peak_id.sample)
    if(length(names.common_control) > length(names.common_sample)) { names.common = names.common_control } else { names.common = names.common_sample }
    names.common  = paste("Common", names.common)
    names.sample = c(paste("Sample", unique(x %>% dplyr::filter(is.na(peak_id.control) & !is.na(peak_id.sample)) %>% .$peak_id.sample)), names.common)
    names.control = c(paste("Control", unique(x %>% dplyr::filter(!is.na(peak_id.control) & is.na(peak_id.sample)) %>% .$peak_id.control)), names.common)

    # Chart
    p = VennDiagram::venn.diagram(
      x = list(names.sample, names.control),
      category.names = c("Sample" , "Control"),
      lwd = 2,
      lty = 'blank',
      fill = RColorBrewer::brewer.pal(3, "Pastel2")[1:2],
      cat.cex = 2,
      cat.fontface = "bold",
      filename = NULL
    )
    jpeg("reports/control2sample_venn.jpg", width=400, height=400)
    grid.draw(p)
    dev.off()
  }

}

main = function() {
  chromosomes_map_df = readr::read_tsv("data/mm9_chromosomes_synonyms.tsv")
  chromosomes_map = chromosomes_map_df$chrom
  names(chromosomes_map) = chromosomes_map_df$chrom_synonym

  # Read offtargets
  offtargets_all_df = readr::read_tsv("data/offtargets_predicted.tsv") %>%
    dplyr::mutate(offtarget_id=1:n())
  offtargets_all_df %>%
    dplyr::filter(bait_chrom==offtarget_chrom) %>%
    dplyr::mutate(offtarget_name=stringr::str_glue("{bait}_{mismatches}_{seq}", bait=bait_name, mismatches=offtarget_mismatches, seq=offtarget_sequence), offtarget_score=0) %>%
    dplyr::select(offtarget_chrom, offtarget_start, offtarget_end, offtarget_name, offtarget_score, offtarget_strand) %>%
    readr::write_tsv("data/macs2/offtargets.bed", col_names=F)

  offtargets_df = offtargets_all_df %>%
    dplyr::mutate(bait_chrom=chromosomes_map[tolower(bait_chrom)], offtarget_chrom=chromosomes_map[tolower(offtarget_chrom)]) %>%
    dplyr::filter(offtarget_mismatches<=4) %>%
    dplyr::mutate(
      offtarget_start=ifelse(offtarget_mismatches==0 & bait_chrom=="chr15", 61818880, offtarget_start),
      offtarget_end=ifelse(offtarget_mismatches==0 & bait_chrom=="chr15", 61818881, offtarget_end))
  targets_df = offtargets_df %>%
    dplyr::filter(offtarget_mismatches==0)

  bed_cols = cols(chrom=col_character(), start=col_double(), end=col_double(), name=col_character(), score=col_character(), strand=col_character())

  junctions_df = read_breaks_experiments("data/breaks_tlx", chromosomes_map=chromosomes_map) %>%
    dplyr::filter(grepl("chr[0-9]", break_chrom))


  #
  # Find breaks close to bait
  #
  junctions_ranges = GenomicRanges::makeGRangesFromDataFrame(junctions_df %>% dplyr::select(-dplyr::matches("is_near_bait|is_very_near_bait|is_very_near_offtarget")) %>% dplyr::mutate(start=break_start, end=break_end, seqnames=break_chrom), keep.extra.columns=T)
  offtargets_ranges_near = GenomicRanges::makeGRangesFromDataFrame(offtargets_df %>% dplyr::filter(offtarget_mismatches<=3) %>% dplyr::mutate(offtarget_bait_chrom=bait_chrom, start=offtarget_start-2e6, end=offtarget_end+2e6, seqnames=offtarget_chrom), keep.extra.columns=T)
  offtarget2junctions.map = as.data.frame(IRanges::mergeByOverlaps(offtargets_ranges_near, junctions_ranges)) %>%
    dplyr::filter(break_chrom==offtarget_chrom & break_bait_chrom==offtarget_bait_chrom) %>%
    dplyr::mutate(is_near_bait=ifelse(offtarget_mismatches==0 & offtarget_start-2e6<=break_start & break_start<=offtarget_end+2e6, offtarget_chrom, "")) %>%
    dplyr::mutate(is_very_near_bait=ifelse(offtarget_mismatches==0 & offtarget_start-1e6<=break_start & break_start<=offtarget_end+1e6, offtarget_chrom, "")) %>%
    dplyr::mutate(is_very_near_offtarget=ifelse(offtarget_start-1e6<=break_start & break_start<=offtarget_end+1e6, offtarget_chrom, "")) %>%
    dplyr::group_by(break_id) %>%
    dplyr::summarise(
      is_near_bait=paste0(unique(is_near_bait[is_near_bait!=""]), collapse=";"),
      is_very_near_bait=paste0(unique(is_very_near_bait[is_very_near_bait!=""]), collapse=";"),
      is_very_near_offtarget=paste0(unique(is_very_near_offtarget[is_very_near_offtarget!=""]), collapse=";"))
  junctions_df = junctions_df %>%
    dplyr::select(-dplyr::matches("is_near_bait|is_very_near_bait|is_very_near_offtarget")) %>%
    dplyr::left_join(offtarget2junctions.map, by="break_id") %>%
    dplyr::mutate(
      is_near_bait=ifelse(is_near_bait!="", is_near_bait, NA_character_),
      is_very_near_bait=ifelse(is_very_near_bait!="", is_very_near_bait, NA_character_),
      is_very_near_offtarget=ifelse(is_very_near_offtarget!="", is_very_near_offtarget, NA_character_))


  rdc = readr::read_tsv("data/rdc_pnas.tsv") %>%
    dplyr::mutate(rdc_length=rdc_end-rdc_start) %>%
    tidyr::separate_rows(rdc_gene, sep=" *, *") %>%
    dplyr::filter(rdc_gene != "--") %>%
    dplyr::group_by(rdc_gene) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(rdc_id=1:n())

  #
  # Calculate normalization accross experiment
  #
  scale_factor_df = junctions_df %>%
    dplyr::filter(!grepl("chrX|chrY", break_bait_chrom)) %>%
    dplyr::group_by(break_exp_condition, break_bait_chrom) %>%
    dplyr::summarise(libsize=n()) %>%
    dplyr::mutate(libsize_median=median(libsize)) %>%
    dplyr::mutate(scale_factor=libsize_median/libsize) %>%
    dplyr::select(break_exp_condition, break_bait_chrom, scale_factor)
  junctions_df = junctions_df %>%
    dplyr::select(-dplyr::matches("scale_factor")) %>%
    dplyr::inner_join(scale_factor_df, by=c("break_exp_condition", "break_bait_chrom"))


  # macs2_params = list(extsize=1e5, llocal=1e8, slocal=2e6, qvalue_cutoff=0.05, peak_max_gap=1e5, peak_min_length=20)
  sample_df = junctions_df %>% dplyr::filter(break_exp_condition=="Sample")
  control_df = junctions_df %>% dplyr::filter(break_exp_condition=="Control")
  p = call_peaks(sample_df, control_df)
}