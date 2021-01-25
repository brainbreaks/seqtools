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
source("src/utilities.R")
source("src/visualization.R")


peaks2offtargets_identity = function(peaks_df, offtargets_df, genome_mm9) {
  peaks_sequences = Biostrings::getSeq(genome_mm9, GenomicRanges::makeGRangesFromDataFrame(peaks_df %>% dplyr::mutate(seqnames=peak_chrom, start=peak_start, end=peak_end), keep.extra.columns=T))
  peaks_df$peak_sequence = as.character(peaks_sequences)
  peaks_df$peak_sequence_rev = as.character(Biostrings::reverseComplement(peaks_sequences))
  peaks_df = peaks_df %>%
    dplyr::select(-dplyr::matches("^offtarget_")) %>%
    dplyr::inner_join(offtargets_df %>% dplyr::select(bait_chrom, offtarget_primer_sequence), by=c("peak_chrom"="bait_chrom"))
  alignments = Biostrings::pairwiseAlignment(peaks_df$offtarget_primer_sequence, peaks_df$peak_sequence, type="global-local", gapOpening=10, gapExtension=4)
  alignments_rev = Biostrings::pairwiseAlignment(peaks_df$offtarget_primer_sequence, peaks_df$peak_sequence_rev, type="global-local", gapOpening=10, gapExtension=4)
  alignments = c(alignments, alignments_rev)
  peaks_df = rbind(peaks_df, peaks_df)
  peaks_df$bait_alignment_direction = rep(c("sense", "antisense"), each=nrow(peaks_df)/2)
  peaks_df$bait_alignment_nchar=Biostrings::nchar(alignments)
  peaks_df$bait_sequence=as.character(Biostrings::pattern(alignments))
  peaks_df$bait_alignment_sequence=as.character(Biostrings::subject(alignments))
  peaks_df$bait_alignment_score=Biostrings::score(alignments)
  peaks_df$bait_alignment_identity=Biostrings::pid(alignments)
  peaks_df = peaks_df %>%
    dplyr::arrange(dplyr::desc(bait_alignment_identity)) %>%
    dplyr::filter(bait_alignment_nchar>0) %>%
    dplyr::distinct(peak_chrom, peak_id, offtarget_primer_sequence, .keep_all=T)

  peaks_df
}


generate_background = function(peaks_df, mm9, n=5) {
  random_lengths = data.frame(seqnames=peaks_df$peak_chrom, peak_width=round((peaks_df$peak_end - peaks_df$peak_start)/1e4)*1e4)
  random_peaks_df = as.data.frame(mm9) %>%
    tibble::rownames_to_column("seqnames") %>%
    dplyr::inner_join(random_lengths, by="seqnames") %>%
    dplyr::group_by(seqnames) %>%
    dplyr::do((function(z){
      z.start = sample(z$seqlengths[1]-1e6, nrow(z)*n)
      data.frame(seqnames=z$seqnames[1], seqlengths=z$seqlengths[1]-1e6, start=z.start, peak_width=rep(z$peak_width, each=n)) %>%
        dplyr::mutate(start=ifelse(start+peak_width>=seqlengths, start-peak_width, start), end=start+peak_width)
    })(.)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(peak_chrom=seqnames, peak_start=start, peak_end=end) %>%
    dplyr::mutate(peak_id=1:n()) %>%
    dplyr::select(-seqlengths)

  random_peaks_df
}

fit_baseline = function(coverage_df, binstep, llocal) {
  coverage_df.baseline = coverage_df %>%
    dplyr::mutate(binstep=binstep, llocal=llocal) %>%
    dplyr::mutate(condition="Control") %>%
    dplyr::arrange(seqnames, start) %>%
    dplyr::group_by(condition, seqnames) %>%
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
  junctions_ann = data.frame(junction_file=list.files(path, full.names=T, recursive=T, pattern="*.tlx")) %>%
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

call_peaks_params = function(binsize=1e5, binstep=1e4, extend=1e5, llocal=2e6, minqvalue=0.01, maxgap=5e5, minlen=200) {
  list(binsize=binsize, binstep=binstep, extend=extend, llocal=llocal, minqvalue=-log10(minqvalue), maxgap=maxgap, minlen=minlen)
}

call_peaks = function(sample_df, control_df=NULL, params=call_peaks_params(), debug=F, debug_annotations=c("RDC"=rdc)) {
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


  # slocal background
  options("UCSC.goldenPath.url"="https://hgdownload.cse.ucsc.edu/goldenPath")
  mm9 = GenomeInfoDb::Seqinfo(genome="mm9")[paste0("chr", 1:19)]
  #mm9 = as(tibble::as_tibble(mm9, rownames="s") %>% dplyr::mutate(seqlengths=seqlengths[s==single_break_chr]) %>% tibble::column_to_rownames("s"), "Seqinfo")
  mm9_size = sum(mm9@seqlengths)
  mm9_effective_size=0.74 * mm9_size
  mm9_tiles = unlist(tileGenome(mm9, tilewidth=params$binsize))
  mm9_tiles$tile_id = 1:length(mm9_tiles)
  mm9_txdb.gtf = GenomicFeatures::makeTxDbFromGFF('data/mm9/mm9.refGene.gtf.gz', format="gtf")


  #
  # Caclulate control baseline
  #
  control_ranges = GenomicRanges::makeGRangesFromDataFrame(control_df %>% dplyr::filter(break_chrom==break_bait_chrom), end.field="break_end", start.field="break_start", strand.field="break_strand", seqnames.field="break_chrom", keep.extra.columns=T)
  control_tiles_df.coverage = calculate_coverage2(control_ranges, mm9, extend=params$extend, binsize=params$binsize, binstep=params$binstep)
  control_df.baseline = fit_baseline(control_tiles_df.coverage, binstep=params$binstep, llocal=params$llocal)


  #
  # Caclulate sample baseline
  #
  sample_ranges = GenomicRanges::makeGRangesFromDataFrame(sample_df %>% dplyr::filter(break_chrom==break_bait_chrom), end.field="break_end", start.field="break_start", strand.field="break_strand", seqnames.field="break_chrom", keep.extra.columns=T)
  #sample_df.coverage = calculate_coverage(sample_ranges, mm9, params$extend)
  sample_tiles_df.coverage = calculate_coverage2(sample_ranges, mm9, extend=params$extend, binsize=params$binsize, binstep=params$binstep)
  sample_df.baseline = fit_baseline(sample_tiles_df.coverage, binstep=params$binstep, llocal=params$llocal)


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
    background_df=mixed_df.baseline %>% dplyr::select(seqnames, start, end, coverage=coverage_smooth),
    minqvalue=params$minqvalue, maxgap=params$maxgap, minlen=pparams$minlen)

  control_peaks = macs_call_peaks(
    signal_df=control_tiles_df.coverage %>% dplyr::select(seqnames, start, end, coverage),
    background_df=control_df.baseline %>% dplyr::select(seqnames, start, end, coverage=coverage_smooth),
    minqvalue=params$minqvalue, maxgap=params$maxgap, minlen=pparams$minlen)

  sample_peaks_ranges = GenomicRanges::makeGRangesFromDataFrame(sample_peaks$peaks %>% dplyr::select(seqnames=peak_chrom, start=peak_start, end=peak_end, peak_id.sample=peak_id), keep.extra.columns=T)
  control_peaks_ranges = GenomicRanges::makeGRangesFromDataFrame(control_peaks$peaks %>% dplyr::select(seqnames=peak_chrom, start=peak_start, end=peak_end, peak_id.control=peak_id), keep.extra.columns=T)
  control_qvalues_ranges = GenomicRanges::makeGRangesFromDataFrame(control_peaks$qvalues %>% dplyr::select(seqnames=qvalue_chrom, start=qvalue_start, end=qvalue_end, qvalue_score.control=qvalue_score), keep.extra.columns=T)
  sample2control_peaks.map = as.data.frame(IRanges::mergeByOverlaps(sample_peaks_ranges, control_qvalues_ranges)) %>% dplyr::select(peak_id.sample, qvalue_score.control)
  sample_peaks$peaks = sample_peaks$peaks %>%
    dplyr::select(-dplyr::matches("^qvalue_score.control$")) %>%
    dplyr::left_join(sample2control_peaks.map, by=c("peak_id"="peak_id.sample")) %>%
    dplyr::arrange(dplyr::desc(qvalue_score.control)) %>%
    dplyr::distinct(peak_id, .keep_all=T)



  if(debug) {
    genome_mm9 = Biostrings::readDNAStringSet("data/mm9/mm9.fa.gz", "fasta")
    sample_peaks_df2offtargets = peaks2offtargets_identity(sample_peaks$peaks, targets_df, genome_mm9)
    pdf(file="reports/duo_baseline_extend1e5_smooth2e6_bin1e5_step1e4_2.pdf", width=15, height=8)
    for(chr in paste0("chr", 1:19)) {
      print(chr)
      coverage_df = dplyr::bind_rows(control_tiles_df.coverage %>% dplyr::mutate(experimental_condition="Control"), sample_tiles_df.coverage %>% dplyr::mutate(experimental_condition="Sample")) %>%
         dplyr::filter(seqnames %in% chr)
      baseline_df = dplyr::bind_rows(
        control_df.baseline_adjusted %>% dplyr::mutate(experimental_condition="Control"),
        sample_df.baseline %>% dplyr::mutate(experimental_condition="Sample"),
        mixed_df.baseline %>% dplyr::mutate(experimental_condition="Max")) %>%
         dplyr::filter(seqnames %in% chr)
      peaks_df = dplyr::bind_rows(
        control_peaks$peaks %>% dplyr::mutate(experimental_condition="Control"),
        sample_peaks$peaks %>% dplyr::mutate(experimental_condition="Sample"),
        sample_peaks$peaks %>% dplyr::filter(qvalue_score.control < peak_summit_qvalue) %>%  dplyr::mutate(condition="Sample filtered"),
        sample_peaks_df2offtargets %>% dplyr::mutate(seqnames=peak_chrom, peak_summit_qvalue=bait_alignment_identity, score=(bait_alignment_identity-min(bait_alignment_identity))/(100-max(bait_alignment_identity))*max(coverage_df$coverage), condition="Offtarget"),
        rdc %>% dplyr::mutate(experimental_condition="RDC", peak_start=rdc_start, peak_end=rdc_end, peak_chrom=rdc_chrom)) %>%
         dplyr::mutate(seqnames=peak_chrom, score=max(coverage_df$coverage)) %>%
         dplyr::filter(seqnames %in% chr)
      karyoplot.colors=c(enriched="#0000FF", Control="#000000", Sample="#FF0000", Offtarget="#00FF00", "Sample filtered"="#FF0000", Max="#0000FF", Coverage="#000000", RDC="#FFFF00")
      p = ggplot_karyoplot(
        coverage_df=coverage_df,
        baseline_df=baseline_df,
        peaks_df=peaks_df,
        colors=karyoplot.colors
      )
      print(p)
    }
    dev.off()
  }



  jpeg("reports/rdc2peaks_venn.jpg", width=800, height=800)
  peaks_filtered_df = sample_peaks$peaks %>% dplyr::filter(qvalue_score.control < peak_summit_qvalue) %>%  dplyr::mutate(condition="Sample filtered")
  rdc_ranges = GenomicRanges::makeGRangesFromDataFrame(rdc %>% dplyr::mutate(seqnames=rdc_chrom, start=rdc_start, end=rdc_end), keep.extra.columns=T)
  peaks_filtered_ranges = GenomicRanges::makeGRangesFromDataFrame(peaks_filtered_df %>% dplyr::mutate(seqnames=peak_chrom, start=peak_start, end=peak_end), keep.extra.columns=T)
  venn_ranges(rdc_ranges, peaks_filtered_ranges, name1="RDC", name2="Peak")
  dev.off()


  if(F) {
    pdf(file="reports/raw_control_vs_sample.pdf", width=20, height=8)
    chr = paste0("chr", 19)
    d = dplyr::bind_rows(control_tiles_df.coverage %>% dplyr::mutate(experimental_condition="Control"), control_tiles_df.coverage %>% dplyr::mutate(experimental_condition="Sample")) %>%
      dplyr::mutate(coverage=ifelse(coverage<50, coverage, NA_real_))
    b = dplyr::bind_rows(control_df.baseline %>% dplyr::mutate(experimental_condition="Control"), sample_df.baseline %>% dplyr::mutate(experimental_condition="Sample")) %>%
      dplyr::mutate(coverage_smooth=ifelse(coverage_smooth<50, coverage_smooth, NA_real_))
    p = ggplot(d %>% dplyr::filter(seqnames %in% chr)) +
      ggridges::geom_ridgeline(aes(y=experimental_condition, x=start, height=coverage, fill=experimental_condition), scale=0.1) +
      ggridges::geom_ridgeline(aes(y=experimental_condition, x=start, height=coverage_smooth), scale=0.1, data=b %>% dplyr::filter(seqnames %in% chr), color="#FF0000", alpha=0.2, size=2) +
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

  # Read offtargets
  offtargets_df = readr::read_tsv("data/offtargets_predicted.tsv") %>%
    dplyr::inner_join(chromosomes_map_df, by=c("bait_chrom"="chrom_synonym")) %>%
    dplyr::mutate(bait_chrom=unique_chrom) %>%
    dplyr::select(-unique_chrom) %>%
    dplyr::inner_join(chromosomes_map_df, by=c("offtarget_chrom"="chrom_synonym")) %>%
    dplyr::mutate(offtarget_chrom=unique_chrom) %>%
    dplyr::select(-unique_chrom) %>%
    dplyr::mutate(offtarget_id=1:n()) %>%
    dplyr::filter(offtarget_mismatches<=4)
  targets_df = offtargets_df %>%
    dplyr::filter(offtarget_mismatches==0)

  #
  # Calculate normalization accross experiment
  #
  junctions_df = read_junctions("data/breaks_tlx") %>%
    dplyr::inner_join(chromosomes_map_df, by=c("bait_chrom"="chrom_synonym")) %>%
    dplyr::mutate(bait_chrom=unique_chrom) %>%
    dplyr::select(-unique_chrom) %>%
    dplyr::inner_join(chromosomes_map_df, by=c("junction_chrom"="chrom_synonym")) %>%
    dplyr::mutate(junction_chrom=unique_chrom) %>%
    dplyr::select(-unique_chrom)
  scale_factor_df = junctions_df %>%
    dplyr::group_by(condition, bait_chrom) %>%
    dplyr::summarise(libsize=n()) %>%
    dplyr::mutate(libsize_median=median(libsize)) %>%
    dplyr::mutate(scale_factor=libsize_median/libsize) %>%
    dplyr::select(condition, bait_chrom, scale_factor)
  junctions_df = junctions_df %>%
    dplyr::select(-dplyr::matches("scale_factor")) %>%
    dplyr::inner_join(scale_factor_df, by=c("condition", "bait_chrom"))

  rdc_cols = readr::cols(rdc_cluster = col_character(), rdc_chrom = col_character(), rdc_start = col_double(), rdc_end = col_double(), rdc_group = col_double(), rdc_gene = col_character())
  rdc_df = readr::read_tsv("data/rdc_pnas.tsv", col_types=rdc_cols) %>%
    dplyr::mutate(rdc_length=rdc_end-rdc_start) %>%
    tidyr::separate_rows(rdc_gene, sep=" *, *") %>%
    dplyr::filter(rdc_gene != "--") %>%
    dplyr::group_by(rdc_gene) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(rdc_id=1:n())


  sample_df = junctions_df %>% dplyr::filter(experimental_condition=="Sample")
  control_df = junctions_df %>% dplyr::filter(experimental_condition=="Control")
  call_peaks_params = call_peaks_params(binsize=1e5, binstep=1e4, extend=1e5, llocal=2e6, minqvalue=0.01, maxgap=5e5, minlen=200)
  p = call_peaks(sample_df, control_df)
}