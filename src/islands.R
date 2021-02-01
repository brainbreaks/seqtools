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

macs_call_islands = function(signal_df, background_df, matching_tiles=T, minqvalue, maxgap, minlen) {
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
       qvalue=qvalue_path, output=peaks_path, cutoff=minqvalue, island_max_gap=maxgap, island_min_length=minlen))

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

fit_baseline = function(coverage_df, binstep, llocal) {
  coverage_df.baseline = coverage_df %>%
    dplyr::mutate(binstep=binstep, llocal=llocal) %>%
    dplyr::mutate(condition="Control") %>%
    dplyr::arrange(seqnames, start) %>%
    dplyr::group_by(seqnames) %>%
    #dplyr::mutate(coverage_mod=coverage)
    dplyr::mutate(is_q2=coverage>quantile(coverage, 0.2, na.rm=T), is_q8=coverage<quantile(coverage, 0.8, na.rm=T), is_inner_quantile=is_q2 & is_q8) %>%
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
    bc.irls = baseline::baseline(x, method="medianWindow", hws=z$llocal[1]/z$binstep[1], hwm=z$llocal[1]/z$binstep[1])
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


calculate_coverage2 = function(ranges, genome_info, params) {
  if(is.null(ranges$score)) {
    ranges$score = 1
    writeLines("Data has no score column. Using '1' as a scale factor")
  }

  genome_original_tiles = unlist(GenomicRanges::tileGenome(genome_info, tilewidth=params$binstep))
  genome_original_tiles$original_start = GenomicRanges::start(genome_original_tiles)
  genome_original_tiles$original_end = GenomicRanges::end(genome_original_tiles)
  genome_tiles = GenomicRanges::trim(GenomicRanges::resize(genome_original_tiles, params$binsize, fix="center"))

  GenomeInfoDb::seqinfo(ranges) = genome_info
  ranges.extended = GenomicRanges::trim(GenomicRanges::resize(ranges, width=params$extend, fix="center"))


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


call_islands_params = function(binsize=1e5, binstep=1e4, extend=1e5, llocal=2e6, minqvalue=0.01, maxgap=5e5, minlen=200) {
  list(binsize=binsize, binstep=binstep, extend=extend, llocal=llocal, minqvalue=-log10(minqvalue), maxgap=maxgap, minlen=minlen)
}

call_islands = function(sample_ranges, control_ranges, genome_info, params=call_islands_params(), debug=F) {
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
  sample_tiles_df.coverage = calculate_coverage2(sample_ranges, genome_info, params=params)
  sample_df.baseline = fit_baseline(sample_tiles_df.coverage, binstep=params$binstep, llocal=params$llocal)

  #
  # Caclulate control baseline
  #
  if(is.null(control_ranges)) {
    control_df.baseline = sample_df.baseline
  } else {
    control_tiles_df.coverage = calculate_coverage2(control_ranges, genome_info, params=params)
    control_df.baseline = fit_baseline(control_tiles_df.coverage, binstep=params$binstep, llocal=params$llocal)
  }

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


  sample_islands = macs_call_islands(
    signal_df=sample_tiles_df.coverage %>% dplyr::select(seqnames, start, end, coverage),
    background_df=mixed_df.baseline %>% dplyr::select(seqnames, start, end, coverage=coverage_smooth),
    minqvalue=params$minqvalue, maxgap=params$maxgap, minlen=params$minlen)

  control_islands = macs_call_islands(
    signal_df=control_tiles_df.coverage %>% dplyr::select(seqnames, start, end, coverage),
    background_df=control_df.baseline %>% dplyr::select(seqnames, start, end, coverage=coverage_smooth),
    minqvalue=params$minqvalue, maxgap=params$maxgap, minlen=params$minlen)

  sample_islands_ranges = GenomicRanges::makeGRangesFromDataFrame(sample_islands$islands %>% dplyr::select(seqnames=island_chrom, start=island_start, end=island_end, island_id.sample=island_id), keep.extra.columns=T)
  control_islands_ranges = GenomicRanges::makeGRangesFromDataFrame(control_islands$islands %>% dplyr::select(seqnames=island_chrom, start=island_start, end=island_end, island_id.control=island_id), keep.extra.columns=T)
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
    baseline=list(sample=sample_df.baseline, control=control_df.baseline_adjusted, mixed=mixed_df.baseline),
    coverage=list(sample=sample_tiles_df.coverage, control=control_tiles_df.coverage)
  )
}
