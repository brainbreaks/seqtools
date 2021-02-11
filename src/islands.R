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
library(Gviz)
library(future)
source("src/utilities.R")
source("src/visualization.R")

macs_call_islands = function(signal_df, background_df, params, matching_tiles=T) {
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
      dplyr::mutate(qvalue_score=tidyr::replace_na(qvalue_score, 0)) %>%
      dplyr::select(qvalue_chrom=seqnames, qvalue_start=start, qvalue_end=end, qvalue_score)
    readr::write_tsv(qvalues_df, qvalue_path, col_names=F)
  } else {
    sample_bdg_path = "data/macs2/bdgcmp_sample.bdg"
    baseline_bdg_path = "data/macs2/bdgcmp_baseline.bdg"
    readr::write_tsv(signal_df %>% dplyr::select(seqnames, start, end, coverage), sample_bdg_path, col_names=F)
    readr::write_tsv(background_df %>% dplyr::select(seqnames, start, end, coverage), baseline_bdg_path, col_names=F)
    system(stringr::str_glue("macs2 bdgcmp -t {input} -c {noise} -m qpois -o {output}", input=sample_bdg_path, noise=baseline_bdg_path, output=qvalue_path))

    qvalue_cols = cols(qvalue_chrom=col_character(), qvalue_start=col_double(), qvalue_end=col_double(), qvalue_score=col_double())
    qvalues_df = readr::read_tsv(qvalue_path, col_names=names(qvalue_cols$cols), col_types=qvalue_cols)
  }
  qvalues_df = qvalues_df %>% dplyr::mutate(qvalue_id=1:n())

  system(stringr::str_glue("macs2 bdgpeakcall -i {qvalue} -c {cutoff} --min-length {format(island_min_length, scientific=F)} --max-gap {format(island_max_gap, scientific=F)} -o {output}",
       qvalue=qvalue_path, output=peaks_path, cutoff=-log10(params$minqvalue), island_max_gap=params$maxgap, island_min_length=params$minlen))

  # Read results file
  peaks_cols = cols(
    island_chrom=col_character(), island_start=col_double(), island_end=col_double(), island_length=col_character(), island_summit_abs=col_double(),
    island_score=col_character(), island_fc=col_double(), island_pvalue_log10=col_double(), island_qvalue_log10=col_double(), island_sammit_offset=col_double()
  )
  islands_df = readr::read_tsv(peaks_path, skip=1, col_names=names(peaks_cols$cols), col_types=peaks_cols) %>%
    dplyr::mutate(island_summit_pos=island_start + island_sammit_offset) %>%
    dplyr::select(island_chrom, island_start, island_end, island_summit_pos)

  if(nrow(islands_df)>0) {
    islands_df = islands_df %>% dplyr::mutate(island_id=1:n())
    qvalues_ranges = GenomicRanges::makeGRangesFromDataFrame(qvalues_df %>% dplyr::mutate(end=qvalue_end, start=qvalue_start, seqnames=qvalue_chrom), keep.extra.columns=T)
    summit_ranges = GenomicRanges::makeGRangesFromDataFrame(islands_df %>% dplyr::mutate(end=island_summit_pos, start=island_summit_pos, seqnames=island_chrom), keep.extra.columns=T)
    sample_ranges = GenomicRanges::makeGRangesFromDataFrame(signal_df %>% dplyr::mutate(break_chrom=seqnames, break_start=start, break_end=end, break_coverage=coverage), keep.extra.columns=T)
    summit2qvalue.map = as.data.frame(IRanges::mergeByOverlaps(qvalues_ranges, summit_ranges)) %>%
      dplyr::filter(qvalue_chrom==island_chrom) %>%
      dplyr::group_by(island_id) %>%
      dplyr::summarise(island_summit_qvalue=max(qvalue_score), .groups="keep")
    summit2sample.map = as.data.frame(IRanges::mergeByOverlaps(sample_ranges, summit_ranges)) %>%
      dplyr::filter(break_chrom==break_chrom) %>%
      dplyr::group_by(island_id) %>%
      dplyr::summarise(island_summit_abs=max(break_coverage), .groups="keep")

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

peaks2offtargets_identity = function(islands_df, offtargets_df, genome_mm9) {
  peaks_sequences = Biostrings::getSeq(genome_mm9, GenomicRanges::makeGRangesFromDataFrame(islands_df %>% dplyr::mutate(seqnames=island_chrom, start=island_start, end=island_end), keep.extra.columns=T))
  islands_df$island_sequence = as.character(peaks_sequences)
  islands_df$island_sequence_rev = as.character(Biostrings::reverseComplement(peaks_sequences))
  islands_df = islands_df %>%
    dplyr::select(-dplyr::matches("^offtarget_")) %>%
    dplyr::inner_join(offtargets_df %>% dplyr::select(bait_chrom, offtarget_primer_sequence), by=c("island_chrom"="bait_chrom"))
  alignments = Biostrings::pairwiseAlignment(islands_df$offtarget_primer_sequence, islands_df$island_sequence, type="global-local", gapOpening=10, gapExtension=4)
  alignments_rev = Biostrings::pairwiseAlignment(islands_df$offtarget_primer_sequence, islands_df$island_sequence_rev, type="global-local", gapOpening=10, gapExtension=4)
  alignments = c(alignments, alignments_rev)
  islands_df = rbind(islands_df, islands_df)
  islands_df$bait_alignment_direction = rep(c("sense", "antisense"), each=nrow(islands_df)/2)
  islands_df$bait_alignment_nchar=Biostrings::nchar(alignments)
  islands_df$bait_sequence=as.character(Biostrings::pattern(alignments))
  islands_df$bait_alignment_sequence=as.character(Biostrings::subject(alignments))
  islands_df$bait_alignment_score=Biostrings::score(alignments)
  islands_df$bait_alignment_identity=Biostrings::pid(alignments)
  islands_df = islands_df %>%
    dplyr::arrange(dplyr::desc(bait_alignment_identity)) %>%
    dplyr::filter(bait_alignment_nchar>0) %>%
    dplyr::distinct(island_chrom, island_id, offtarget_primer_sequence, .keep_all=T)

  islands_df
}

fit_baseline = function(coverage_df, params) {
  coverage_df.baseline_1 = coverage_df %>%
    dplyr::mutate(binstep=params$binstep, llocal=params$llocal) %>%
    dplyr::group_by(seqnames) %>%
    dplyr::filter(max(coverage)>0) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(seqnames, start)

  coverage_df.baseline_2 = coverage_df.baseline_1 %>%
    dplyr::group_by(seqnames) %>%
    dplyr::summarize(q1=quantile(coverage, params$baseline_quantiles[1], na.rm=T), q2=quantile(coverage, params$baseline_quantiles[2], na.rm=T), q_str=stringr::str_glue("{chr}: Specified quantiles/coverage => ({p1}/{q1}, {p2}/{q2}) covers {prct}% chromosome", chr=seqnames[1], p1=params$baseline_quantiles[1], p2=params$baseline_quantiles[2], q1=q1, q2=q2, prct=round(mean(dplyr::between(coverage, q1, q2))*100)), .groups="keep")
  writeLines(coverage_df.baseline_2$q_str)

  coverage_df.baseline = coverage_df.baseline_1 %>%
    dplyr::inner_join(coverage_df.baseline_2 %>% dplyr::select(seqnames, q1, q2), by="seqnames") %>%
    dplyr::mutate(is_inner_quantile=dplyr::between(coverage, q1, q2)) %>%
    dplyr::group_by(seqnames) %>%
    dplyr::mutate(coverage_median=mean(coverage[is_inner_quantile], na.rm=T), coverage_sd=sd(coverage[is_inner_quantile], na.rm=T)) %>%
    dplyr::mutate(coverage_min=pmax(coverage_median-2*coverage_sd, 0)) %>%
    dplyr::mutate(coverage_mod=ifelse(coverage<coverage_min, NA_real_, coverage)) %>%
    dplyr::mutate(coverage_mod=zoo::na.fill(coverage_mod, "extend"))

  parallelCluster = parallel::makeCluster(19, type="PSOCK", methods=FALSE)
  doParallel::registerDoParallel(parallelCluster)
  coverage_dflist.baseline = coverage_df.baseline %>% dplyr::group_by(seqnames) %>% dplyr::group_split()
  coverage_df.baseline_output = foreach::foreach(z=coverage_dflist.baseline, .combine=rbind) %dopar% {
    library(dplyr)
    x = matrix(as.numeric(z$coverage_mod), nrow=1)
    colnames(x) = z$start
    bc.irls = baseline::baseline(x, method="medianWindow", hws=z$llocal[1]/z$binstep[1], hwm=z$llocal[1]/z$binstep[1])
    z.coverage_smooth = as.numeric(bc.irls@baseline)
    cbind(z[c("seqnames", "start", "end", "coverage", "coverage_mod")], coverage_smooth=z.coverage_smooth)
      # dplyr::mutate(q1=quantile(coverage_smooth, params$baseline_quantiles[1], na.rm=T), q2=quantile(coverage_smooth, params$baseline_quantiles[2], na.rm=T))
      # dplyr::mutate(is_inner_quantile=dplyr::between(coverage, q1, q2)) %>%
      # dplyr::mutate(coverage_smooth_median=median(coverage[is_inner_quantile], na.rm=T))
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


  # ggplot(coverage_df.baseline_output) +
  #   geom_line(aes(x=start, y=coverage_mod)) +
  #   geom_line(aes(x=start, y=coverage_smooth)) +
  #   geom_line(aes(x=start, y=log10(qvalue), color="red")) +
  #   geom_hline(yintercept=-2) +
  #   facet_wrap(~seqnames)


  return(coverage_df.baseline_output)
}


read_junctions = function(path) {
  #
  # Read breaks
  #
  tlx_cols = cols(
    Qname=readr::col_character(), JuncID=readr::col_character(), Rname=readr::col_character(), Junction=readr::col_double(),
    Strand=readr::col_character(), Rstart=readr::col_double(), Rend=readr::col_double(),
    B_Rname=readr::col_character(), B_Rstart=readr::col_double(), B_Rend=readr::col_double(), B_Strand=readr::col_double(),
    B_Qstart=readr::col_double(), B_Qend=readr::col_double(), Qstart=readr::col_double(), Qend=readr::col_double(), Qlen=readr::col_double(),
    B_Cigar=readr::col_character(), Cigar=readr::col_character(), Seq=readr::col_character(), J_Seq=readr::col_character(), Barcode=readr::col_logical(),
    unaligned=readr::col_double(), baitonly=readr::col_double(), uncut=readr::col_double(), misprimed=readr::col_double(), freqcut=readr::col_double(),
    largegap=readr::col_double(), mapqual=readr::col_double(), breaksite=readr::col_double(), sequential=readr::col_double(), repeatseq=readr::col_double(), duplicate=readr::col_double(), Note=col_character()
  )

  # Rename Chr14 to Alt289 (Was Alt289 and only PW255 was Alt287)
  #junctions_ann.excluded = "PW121|PW246|PW247|PW248|PW249|JK096|JK097|JK098|JK099|JK100|JK101"
  junctions_ann.excluded = "NO FILTER"
  junctions_ann = data.frame(junction_file=Sys.glob(path)) %>%
     dplyr::filter(!grepl(junctions_ann.excluded, junction_file)) %>%
     dplyr::mutate(
       bait_chrom = basename(dirname(junction_file)),
       experimental_condition = ifelse(grepl("DMSO", junction_file), "Control", "Sample"),
       junction_tech_sample = gsub("_.*", "", basename(junction_file)),
       junction_bio_sample=gsub(".*_((Alt|R)[0-9]+).*", "\\1", basename(junction_file), ignore.case=T, perl=T)
     ) %>%
    dplyr::group_by(bait_chrom, experimental_condition) %>%
    dplyr::mutate(junction_replicate=1:n()) %>%
    dplyr::ungroup()

  parallelCluster = parallel::makeCluster(20, type="PSOCK", methods=FALSE)
  doParallel::registerDoParallel(parallelCluster)
  junctions_ann_list = junctions_ann %>%
    dplyr::rowwise() %>%
    dplyr::group_split()
  junctions_df = foreach::foreach(tlx_ann=junctions_ann_list, .combine=rbind) %dopar% {
    tlx = readr::read_tsv(tlx_ann$junction_file, col_types=tlx_cols)
    tlx$junction_chrom = tlx$Rname
    tlx$junction_strand = ifelse(tlx$Strand=="-1", "-", "+")
    tlx$experimental_condition = tlx_ann$experimental_condition
    tlx$junction_file = tlx_ann$junction_file
    tlx$bait_chrom = tlx_ann$bait_chrom
    tlx$junction_replicate = tlx_ann$junction_replicate
    tlx$junction_chrom = tlx$Rname
    tlx$junction_start = tlx$Junction
    tlx$junction_end = tlx$Junction
    tlx$junction_name = tlx$Qname
    cbind(tlx_ann[,c("junction_file", "bait_chrom", "experimental_condition", "junction_replicate")], tlx[,c("junction_chrom", "junction_start", "junction_end", "junction_name", "junction_strand")])
  }
  stopCluster(parallelCluster)

  junctions_df %>% dplyr::mutate(junction_id=1:n())
}


calculate_coverage2 = function(ranges, params) {
  if(is.null(ranges$score)) {
    ranges$score = 1
    writeLines("Data has no score column. Using '1' as a scale factor")
  }

  genome_original_tiles = unlist(GenomicRanges::tileGenome(GenomeInfoDb::seqinfo(ranges), tilewidth=params$binstep))
  genome_original_tiles$original_start = GenomicRanges::start(genome_original_tiles)
  genome_original_tiles$original_end = GenomicRanges::end(genome_original_tiles)
  genome_tiles = suppressWarnings(GenomicRanges::trim(GenomicRanges::resize(genome_original_tiles, params$binsize, fix="center")))

  ranges.extended = suppressWarnings(GenomicRanges::trim(GenomicRanges::resize(ranges, width=params$extend, fix="center")))

  hits = as(findOverlaps(genome_tiles, ranges.extended), "List")
  coverage_df = as.data.frame(genome_tiles) %>%
    dplyr::mutate(coverage = sum(extractList(ranges.extended$score, hits))) %>%
    dplyr::mutate(start=original_start, end=original_end) %>%
    dplyr::group_by(seqnames, start, end) %>%
    dplyr::summarise(coverage=sum(coverage), .groups="keep") %>%
    dplyr::ungroup()

  #write.table(coverage_df %>% dplyr::select(seqnames, start, end, coverage), col.names=F, row.names=F, quote=F, file="data/macs2/coverage.bdg")

  class(coverage_df) = c("tbl_df", "tbl", "data.frame")
  coverage_df
}


call_islands_params = function(binsize=1e5, binstep=1e4, extend=1e5, llocal=2e6, minqvalue=0.01, maxgap=5e5, minlen=200, baseline_quantiles=c(0, 0.8)) {
  list(binsize=binsize, binstep=binstep, extend=extend, llocal=llocal, minqvalue=minqvalue, maxgap=maxgap, minlen=minlen, baseline_quantiles=baseline_quantiles)
}

call_islands = function(sample_ranges, control_ranges, params=call_islands_params(), debug=F) {
  sample_bed_path = "data/breakensembl/Sample.bed"
  control_bed_path = "data/breakensembl/Control.bed"

  if(debug) {
    readr::write_tsv(as.data.frame(sample_ranges), sample_bed_path, col_names=F)
  }

  if(is.null(control_ranges)) {
    control_ranges = sample_ranges
  } else {
    if(debug) {
      readr::write_tsv(as.data.frame(control_ranges), control_bed_path, col_names=F)
    }
  }

  if(is.null(sample_ranges$score)) {
    sample_ranges$score = 1
    writeLines("Sample data has no score column. Using '1' as a scale factor")
  }

  if(is.null(control_ranges$score)) {
    control_ranges$score = 1
    writeLines("Sample data has no score column. Using '1' as a scale factor")
  }

  #
  # Caclulate sample baseline
  #
  sample_tiles_df.coverage = calculate_coverage2(ranges=sample_ranges, params=params)
  sample_df.baseline = fit_baseline(coverage_df=sample_tiles_df.coverage, params=params)

  #
  # Caclulate control baseline
  #
  if(is.null(control_ranges)) {
    control_df.baseline = sample_df.baseline
  } else {
    control_tiles_df.coverage = calculate_coverage2(ranges=control_ranges, params=params)
    control_df.baseline = fit_baseline(coverage_df=control_tiles_df.coverage, params=params)
  }

  #
  # coverage_df.baseline_2 = coverage_df.baseline_1 %>%
  #   dplyr::group_by(seqnames) %>%
  #   dplyr::summarize(q1=quantile(coverage, params$baseline_quantiles[1], na.rm=T), q2=quantile(coverage, params$baseline_quantiles[2], na.rm=T), q_str=stringr::str_glue("{chr}: Specified quantiles/coverage => ({p1}/{q1}, {p2}/{q2}) covers {prct}% chromosome", chr=seqnames[1], p1=mean_quantiles[1], p2=mean_quantiles[2], q1=q1, q2=q2, prct=round(mean(dplyr::between(coverage, q1, q2))*100)), .groups="keep") %>%
  #   dplyr::mutate()
  # writeLines(coverage_df.baseline_2$q_str)
  #
  # coverage_df.baseline = coverage_df.baseline_1 %>%
  #   dplyr::inner_join(coverage_df.baseline_2 %>% dplyr::select(seqnames, q1, q2), by="seqnames") %>%
  #   dplyr::mutate(is_inner_quantile=dplyr::between(coverage, q1, q2)) %>%
  #   dplyr::group_by(seqnames) %>%
  #   dplyr::mutate(coverage_median=median(coverage[is_inner_quantile], na.rm=T), coverage_sd=sd(coverage[is_inner_quantile], na.rm=T)) %>%

  baseline_adjusted_df_1 = control_df.baseline %>%
    dplyr::select(seqnames, start, end, coverage.control=coverage, coverage_smooth.control_vs_control=coverage_smooth) %>%
    dplyr::inner_join(sample_df.baseline %>% dplyr::select(seqnames, start, end, coverage.sample=coverage, coverage_smooth.sample_vs_sample=coverage_smooth), by=c("seqnames", "start", "end"))

  baseline_adjusted_scale_factor = baseline_adjusted_df_1 %>%
    dplyr::group_by(seqnames) %>%
    dplyr::mutate(
      qc_q1=quantile(coverage.control, params$baseline_quantiles[1], na.rm=T), qc_q2=quantile(coverage.control, params$baseline_quantiles[2], na.rm=T),
      qs_q1=quantile(coverage.sample, params$baseline_quantiles[1], na.rm=T), qs_q2=quantile(coverage.sample, params$baseline_quantiles[2], na.rm=T),
      qc=dplyr::between(coverage.control, qc_q1, qc_q2),
      qs=dplyr::between(coverage.sample, qs_q1, qs_q2)) %>%
    dplyr::summarize(
      scale_factor.control=mean(coverage.sample[qs])/mean(coverage.control[qc]),
      q_str=stringr::str_glue("{chr} | control scale factor = {scale_factor} | baseline percentiles = ({p1}, {p2}) | [control] baseline quantiles = ({c_q1}, {c_q2}) covers {c_prct}% chromosome (mean: {c_mean}) | [sample] baseline quantiles => ({s_q1}, {s_q2}) covers {s_prct}% chromosome (mean: {s_mean})",
            chr=seqnames[1], p1=params$baseline_quantiles[1], p2=params$baseline_quantiles[2], scale_factor=round(scale_factor.control, 2),
            c_q1=qc_q1[1], c_q2=qc_q2[2], c_prct=round(mean(qc)*100), c_mean=round(mean(coverage.control[qc]), 2),
            s_q1=qs_q1[1], s_q2=qs_q2[2], s_prct=round(mean(qs)*100), s_mean=round(mean(coverage.sample[qs]), 2)
      ), .groups="keep")
  writeLines(baseline_adjusted_scale_factor$q_str)

  #
  # Adjust control baseline to sample library
  #
  baseline_adjusted_df = baseline_adjusted_df_1 %>%
    dplyr::inner_join(baseline_adjusted_scale_factor %>% dplyr::select(seqnames, scale_factor.control), by="seqnames") %>%
    dplyr::group_by(seqnames) %>%
    dplyr::mutate(
      coverage.sample_scaled=coverage.sample*scale_factor.control,
      coverage_smooth.sample_vs_control=coverage_smooth.control_vs_control*scale_factor.control,
      coverage_smooth.sample_vs_max=pmax(coverage_smooth.sample_vs_control, coverage_smooth.sample_vs_sample),
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      pvalue.control_vs_control=pgamma(coverage.control, shape=coverage_smooth.control_vs_control, rate=1, lower.tail=F),
      coverage_th01.control_vs_control=qgamma(params$minqvalue, coverage_smooth.control_vs_control, 1, lower.tail=F),
      pvalue.sample_vs_sample=pgamma(coverage.sample, shape=coverage_smooth.sample_vs_sample, rate=1, lower.tail=F),
      coverage_th01.sample_vs_sample=qgamma(params$minqvalue, coverage_smooth.sample_vs_sample, 1, lower.tail=F),
      pvalue.sample_vs_control=pgamma(coverage.sample, shape=coverage_smooth.sample_vs_control, rate=1, lower.tail=F),
      coverage_th01.sample_vs_control=qgamma(params$minqvalue, coverage_smooth.sample_vs_sample, 1, lower.tail=F),
      pvalue.sample_vs_max=pgamma(coverage.sample, shape=coverage_smooth.sample_vs_max, rate=1, lower.tail=F),
      coverage_th01.sample_vs_max=qgamma(params$minqvalue, coverage_smooth.sample_vs_sample, 1, lower.tail=F)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      qvalue.control_vs_control=qvalue::qvalue(pvalue.control_vs_control)$qvalues,
      qvalue.sample_vs_sample=qvalue::qvalue(pvalue.sample_vs_sample)$qvalues,
      qvalue.sample_vs_control=qvalue::qvalue(pvalue.sample_vs_control)$qvalues,
      qvalue.sample_vs_max=qvalue::qvalue(pvalue.sample_vs_max)$qvalues) %>%
    dplyr::select(seqnames, start, end, dplyr::matches("\\.(control|sample)"))


  sample_islands = macs_call_islands(
    signal_df=baseline_adjusted_df %>% dplyr::mutate(coverage=coverage.sample),
    background_df=baseline_adjusted_df %>% dplyr::mutate(coverage=coverage_smooth.sample_vs_max),
    params=params)

  control_islands = macs_call_islands(
    signal_df=baseline_adjusted_df %>% dplyr::mutate(coverage=coverage.control),
    background_df=baseline_adjusted_df %>% dplyr::mutate(coverage=coverage_smooth.control_vs_control),
    params=params)

  # ggplot(baseline_adjusted_df) +
  #   geom_line(aes(x=start, y=coverage.sample, color="sample")) +
  #   geom_line(aes(x=start, y=coverage_smooth.sample_vs_max, color="max")) +
  #   coord_cartesian(ylim=c(0, 100))
  #
  # ggplot(sample_islands$qvalues) +
  #   geom_line(aes(x=qvalue_start, y=qvalue_score))


  sample_islands_ranges = GenomicRanges::makeGRangesFromDataFrame(sample_islands$islands %>% dplyr::select(seqnames=island_chrom, start=island_start, end=island_end, island_id.sample=island_id), keep.extra.columns=T)
  control_qvalues_ranges = GenomicRanges::makeGRangesFromDataFrame(control_islands$qvalues %>% dplyr::select(seqnames=qvalue_chrom, start=qvalue_start, end=qvalue_end, island_summit_qvalue.control=qvalue_score), keep.extra.columns=T)
  sample2control_islands.map = as.data.frame(IRanges::mergeByOverlaps(sample_islands_ranges, control_qvalues_ranges)) %>% dplyr::select(island_id.sample, island_summit_qvalue.control)
  sample_islands$islands = sample_islands$islands %>%
    dplyr::select(-dplyr::matches("^island_summit_qvalue.control$")) %>%
    dplyr::left_join(sample2control_islands.map, by=c("island_id"="island_id.sample")) %>%
    dplyr::arrange(dplyr::desc(island_summit_qvalue.control)) %>%
    dplyr::distinct(island_id, .keep_all=T)

  hits = as(findOverlaps(sample_islands_ranges, sample_ranges), "List")
  sample_islands$islands$junctions_count = sum(extractList(sample_ranges$score, hits))
  hits = as(findOverlaps(sample_islands_ranges, control_ranges), "List")
  sample_islands$islands$junctions_count.control = sum(extractList(control_ranges$score, hits))

  list(
    islands=list(sample=sample_islands$islands, control=control_islands$islands),
    coverage=baseline_adjusted_df
  )
}

