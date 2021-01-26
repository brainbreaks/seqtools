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


islands2offtargets_identity = function(islands_df, offtargets_df, genome_mm9) {
  islands_sequences = Biostrings::getSeq(genome_mm9, GenomicRanges::makeGRangesFromDataFrame(islands_df %>% dplyr::mutate(seqnames=island_chrom, start=island_start, end=island_end), keep.extra.columns=T))
  islands_df$island_sequence = as.character(islands_sequences)
  islands_df$island_sequence_rev = as.character(Biostrings::reverseComplement(islands_sequences))
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
  coverage_df.baseline_output = foreach::foreach(z=coverage_dflist.baseline, .combine=rbind) %do% {
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

calculate_coverage2 = function(ranges, genome_info, params) {
  genome_original_tiles = unlist(GenomicRanges::tileGenome(genome_info, tilewidth=params$binstep))
  genome_original_tiles$original_start = GenomicRanges::start(genome_original_tiles)
  genome_original_tiles$original_end = GenomicRanges::end(genome_original_tiles)
  genome_tiles = GenomicRanges::trim(GenomicRanges::resize(genome_original_tiles, params$binsize, fix="center"))

  GenomeInfoDb::seqinfo(ranges) = genome_info
  ranges.extended = GenomicRanges::trim(GenomicRanges::resize(ranges, width=params$extend, fix="center"))


  hits = as(findOverlaps(genome_tiles, ranges.extended), "List")
  coverage_df = as.data.frame(genome_tiles) %>%
    dplyr::mutate(coverage = sum(extractList(ranges.extended$scale_factor, hits))) %>%
    dplyr::mutate(start=original_start, end=original_end) %>%
    dplyr::group_by(seqnames, start, end) %>%
    dplyr::summarise(coverage=sum(coverage))

  #write.table(coverage_df %>% dplyr::select(seqnames, start, end, coverage), col.names=F, row.names=F, quote=F, file="data/macs2/coverage.bdg")

  class(coverage_df) = c("tbl_df", "tbl", "data.frame")
  coverage_df
}

call_islands_params = function(binsize=1e5, binstep=1e4, extend=1e5, llocal=2e6, minqvalue=0.01, maxgap=5e5, minlen=200) {
  list(binsize=binsize, binstep=binstep, extend=extend, llocal=llocal, minqvalue=-log10(minqvalue), maxgap=maxgap, minlen=minlen)
}

call_islands = function(sample_ranges, control_ranges, genome_info, params=call_islands_params(), debug=F) {
  # Prepare tiles
  genome_tiles = unlist(tileGenome(genome_info, tilewidth=params$binsize))
  genome_tiles$tile_id = 1:length(genome_tiles)

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

  #
  # Caclulate sample baseline
  #
  sample_tiles_df.coverage = calculate_coverage2(sample_ranges, genome_info, params=params)
  sample_df.baseline = fit_baseline(sample_tiles_df.coverage, binstep=params$binstep, llocal=params$llocal)

  #
  # Caclulate control baseline
  #
  if(is.null(control_df)) {
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
  control_qvalues_ranges = GenomicRanges::makeGRangesFromDataFrame(control_islands$qvalues %>% dplyr::select(seqnames=qvalue_chrom, start=qvalue_start, end=qvalue_end, qvalue_score.control=qvalue_score), keep.extra.columns=T)
  sample2control_islands.map = as.data.frame(IRanges::mergeByOverlaps(sample_islands_ranges, control_qvalues_ranges)) %>% dplyr::select(island_id.sample, qvalue_score.control)
  sample_islands$islands = sample_islands$islands %>%
    dplyr::select(-dplyr::matches("^qvalue_score.control$")) %>%
    dplyr::left_join(sample2control_islands.map, by=c("island_id"="island_id.sample")) %>%
    dplyr::arrange(dplyr::desc(qvalue_score.control)) %>%
    dplyr::distinct(island_id, .keep_all=T)

  list(
    islands=list(sample=sample_islands$islands, control=control_islands$islands),
    baseline=list(sample=sample_df.baseline, control=control_df.baseline_adjusted, mixed=mixed_df.baseline),
    coverage=list(sample=sample_tiles_df.coverage, control=control_tiles_df.coverage)
  )
}

main = function() {
  read_mm9_job = parallel::mcparallel(Biostrings::readDNAStringSet("data/mm9/mm9.fa.gz", "fasta"))

  chromsizes_cols = readr::cols(seqnames=col_character(), seqlengths=col_double())
  genome_info = with(readr::read_tsv("data/mm9/mm9.chrom.sizes", col_names=names(chromsizes_cols$cols), col_types=chromsizes_cols),
             GenomeInfoDb::Seqinfo(seqnames, seqlengths, isCircular=rep(F, length(seqnames)), genome=rep("mm9", length(seqnames))))
  genome_info = genome_info[paste0("chr", c(1:19, "X", "Y"))]



  chromosomes_map_df = readr::read_tsv("data/mm9_chromosomes_synonyms.tsv")

  # Read offtargets
  offtargets_cols = readr::cols(bait_name=col_character(), bait_chrom=col_character(),
    offtarget_chrom=col_character(), offtarget_strand=col_character(), offtarget_start=col_double(), offtarget_end=col_double(),
    offtarget_mismatches=col_double(), offtarget_primer_sequence=col_character(), offtarget_sequence=col_character())
  offtargets_df = readr::read_tsv("data/offtargets_predicted.tsv", col_types=offtargets_cols) %>%
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
    dplyr::select(-unique_chrom) %>%
    dplyr::mutate(junction_name=stringr::str_glue("{bait}_{name}", bait=bait_chrom, name=junction_name))
  scale_factor_df = junctions_df %>%
    dplyr::group_by(experimental_condition, bait_chrom) %>%
    dplyr::summarise(libsize=n()) %>%
    dplyr::mutate(libsize_median=median(libsize)) %>%
    dplyr::mutate(scale_factor=libsize_median/libsize) %>%
    dplyr::select(experimental_condition, bait_chrom, scale_factor)
  junctions_df = junctions_df %>%
    dplyr::select(-dplyr::matches("scale_factor")) %>%
    dplyr::inner_join(scale_factor_df, by=c("experimental_condition", "bait_chrom"))

  rdc_cols = readr::cols(rdc_cluster=col_character(), rdc_chrom=col_character(), rdc_start=col_double(), rdc_end=col_double(), rdc_group=col_double(), rdc_gene=col_character())
  rdc_df = readr::read_tsv("data/rdc_pnas.tsv", col_types=rdc_cols) %>%
    dplyr::inner_join(chromosomes_map_df, by=c("rdc_chrom"="chrom_synonym")) %>%
    dplyr::mutate(rdc_chrom=unique_chrom) %>%
    dplyr::select(-unique_chrom) %>%
    dplyr::mutate(rdc_length=rdc_end-rdc_start) %>%
    dplyr::mutate(rdc_id=1:n())

  sample_df = junctions_df %>% dplyr::filter(experimental_condition=="Sample" & junction_chrom==bait_chrom)
  sample_ranges = GenomicRanges::makeGRangesFromDataFrame(sample_df %>% dplyr::mutate(end=junction_end, start=junction_start, seqnames=junction_chrom), keep.extra.columns=T)
  control_df = junctions_df %>% dplyr::filter(experimental_condition=="Control" & junction_chrom==bait_chrom)
  control_ranges = GenomicRanges::makeGRangesFromDataFrame(control_df %>% dplyr::mutate(end=junction_end, start=junction_start, seqnames=junction_chrom), keep.extra.columns=T)
  params = call_islands_params(binsize=1e5, binstep=1e4, extend=1e5, llocal=2e6, minqvalue=0.01, maxgap=5e5, minlen=200)
  results = call_islands(sample_ranges, control_ranges, genome_info, params)


  # mm9_txdb.gtf = GenomicFeatures::makeTxDbFromGFF('data/mm9/mm9.refGene.gtf.gz', format="gtf")



  if(debug) {
    genome_mm9 = mccollect(list(read_mm9_job))[[1]]
    sample_islands_df2offtargets = islands2offtargets_identity(results$islands$sample, targets_df, genome_mm9)
    pdf(file="reports/duo_baseline_extend1e5_smooth2e6_bin1e5_step1e4_2.pdf", width=15, height=8)
    pdf(file="test.pdf", width=15, height=8)
    for(chr in paste0("chr", 1:19)) {
      print(chr)
      coverage_df = dplyr::bind_rows(
        results$coverage$control %>% dplyr::mutate(experimental_condition="Control"),
        results$coverage$sample %>% dplyr::mutate(experimental_condition="Sample")) %>%
        dplyr::filter(seqnames %in% chr)
      baseline_df = dplyr::bind_rows(
        results$baseline$control %>% dplyr::mutate(experimental_condition="Control"),
        results$baseline$sample %>% dplyr::mutate(experimental_condition="Sample"),
        results$baseline$mixed %>% dplyr::mutate(experimental_condition="Max")) %>%
         dplyr::filter(seqnames %in% chr)
      islands_df = dplyr::bind_rows(
        results$islands$control %>% dplyr::mutate(experimental_condition="Control"),
        results$islands$sample %>% dplyr::mutate(experimental_condition="Sample"),
        sample_islands_df2offtargets %>% dplyr::mutate(seqnames=island_chrom, island_summit_qvalue=bait_alignment_identity, score=(bait_alignment_identity-min(bait_alignment_identity))/(100-max(bait_alignment_identity))*max(coverage_df$coverage), condition="Offtarget"),
        rdc_df %>% dplyr::mutate(experimental_condition="RDC", island_start=rdc_start, island_end=rdc_end, island_chrom=rdc_chrom)) %>%
         dplyr::mutate(seqnames=island_chrom, score=max(results$coverage$sample$coverage)) %>%
         dplyr::filter(seqnames %in% chr)
      karyoplot.colors=c(enriched="#0000FF", Control="#000000", Sample="#FF0000", Offtarget="#00FF00", "Sample filtered"="#FF0000", Max="#0000FF", Coverage="#000000", RDC="#FFFF00")
      p = ggplot_karyoplot(
        coverage_df=coverage_df,
        baseline_df=baseline_df,
        islands_df=islands_df,
        colors=karyoplot.colors
      )
      print(p)
    }
    dev.off()
  }



  jpeg("reports/rdc2islands_venn.jpg", width=800, height=800)
  islands_filtered_df = sample_islands$islands %>% dplyr::filter(qvalue_score.control < island_summit_qvalue) %>%  dplyr::mutate(condition="Sample filtered")
  rdc_ranges = GenomicRanges::makeGRangesFromDataFrame(rdc %>% dplyr::mutate(seqnames=rdc_chrom, start=rdc_start, end=rdc_end), keep.extra.columns=T)
  islands_filtered_ranges = GenomicRanges::makeGRangesFromDataFrame(islands_filtered_df %>% dplyr::mutate(seqnames=island_chrom, start=island_start, end=island_end), keep.extra.columns=T)
  venn_ranges(rdc_ranges, islands_filtered_ranges, name1="RDC", name2="island")
  dev.off()


  #
  # VENN diagram overlap between SAMPLE and CONTROL
  #
  if(debug) {
    control_islands_ranges = GenomicRanges::makeGRangesFromDataFrame(control_islands$islands, end.field="island_end", start.field="island_start", seqnames.field="island_chrom", keep.extra.columns=T)
    sample_islands_ranges = GenomicRanges::makeGRangesFromDataFrame(sample_islands$islands, end.field="island_end", start.field="island_start", seqnames.field="island_chrom", keep.extra.columns=T)
    sample2control.map = as.data.frame(GenomicRanges::findOverlaps(control_islands_ranges, sample_islands_ranges))
    x = sample_islands$islands %>%
      dplyr::rename(island_id.sample="island_id") %>%
      dplyr::full_join(sample2control.map %>% dplyr::rename(island_id.control="queryHits"), by=c("island_id.sample"="subjectHits")) %>%
      dplyr::full_join(control_islands$islands, by=c("island_id.control"="island_id")) %>%
      setNames(gsub("\\.x$", ".sample", colnames(.))) %>%
      setNames(gsub("\\.y$", ".control", colnames(.)))

    names.common_control = unique(x %>% dplyr::filter(!is.na(island_id.control) & !is.na(island_id.sample)) %>% .$island_id.control)
    names.common_sample = unique(x %>% dplyr::filter(!is.na(island_id.control) & !is.na(island_id.sample)) %>% .$island_id.sample)
    if(length(names.common_control) > length(names.common_sample)) { names.common = names.common_control } else { names.common = names.common_sample }
    names.common  = paste("Common", names.common)
    names.sample = c(paste("Sample", unique(x %>% dplyr::filter(is.na(island_id.control) & !is.na(island_id.sample)) %>% .$island_id.sample)), names.common)
    names.control = c(paste("Control", unique(x %>% dplyr::filter(!is.na(island_id.control) & is.na(island_id.sample)) %>% .$island_id.control)), names.common)

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