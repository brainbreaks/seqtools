library(dplyr)
library(reshape2)
library(readr)
library(GenomicRanges)
library(GenomicFeatures)
library(ggplot2)
library(tidyr)
library(baseline)
library(smoother)
library(dbscan)
library(ggpmisc)
library(rtracklayer)
library(ggbeeswarm)
library(plotly)
library(Rtsne)
library(randomForest)
source("src/visualization.R")
source("src/islands.R")


liftOverRanges = function(ranges, chain) {
  ranges$ranges_id = 1:length(ranges)
  as.data.frame(unlist(rtracklayer::liftOver(ranges, chain))) %>%
    dplyr::group_by(ranges_id) %>%
    dplyr::mutate(start=min(start), end=max(end)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(ranges_id, .keep_all=T) %>%
    dplyr::select(-ranges_id)
}

repliseq_read = function(path) {
  repliseq_df = as.data.frame(t(readr::read_tsv(path, col_names=F))) %>%
    dplyr::rename(repliseq_chrom="V1", repliseq_start="V2", repliseq_end="V3") %>%
    dplyr::mutate(repliseq_start=as.numeric(repliseq_start), repliseq_end=as.numeric(repliseq_end)) %>%
    reshape2::melt(id.vars=c("repliseq_chrom", "repliseq_start", "repliseq_end"), variable.name="repliseq_fraction", value.name="repliseq_value") %>%
    dplyr::mutate(repliseq_fraction=as.numeric(gsub("V", "", repliseq_fraction))-3, repliseq_value=as.numeric(repliseq_value))
  repliseq_df.keep = repliseq_df %>%
    dplyr::arrange(repliseq_start) %>%
    dplyr::group_by(repliseq_chrom, repliseq_start, repliseq_end) %>%
    dplyr::summarize(repliseq_isna=sum(is.na(repliseq_value))>10) %>%
    dplyr::group_by(repliseq_chrom) %>%
    dplyr::mutate(first_value_position=min(which(!repliseq_isna)), last_value_position=max(which(!repliseq_isna)), keep_value=dplyr::between(1:n(), first_value_position, last_value_position)) %>%
    dplyr::filter(keep_value) %>%
    dplyr::ungroup() %>%
    dplyr::select(repliseq_chrom, repliseq_start, repliseq_end)
  repliseq_df.f = repliseq_df %>%
    dplyr::inner_join(repliseq_df.keep, by=c("repliseq_chrom", "repliseq_start", "repliseq_end"))

  repliseq_df.f
}

repliseq_summarize = function(repliseq_df, window=5) {
  th.repliseq_value_norm = 0.1
  repliseq_df = repliseq_df %>%
    dplyr::group_by(repliseq_chrom, repliseq_start, repliseq_end) %>%
    dplyr::mutate(repliseq_value_norm=((repliseq_value)/max(repliseq_value))) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(repliseq_chrom, repliseq_start, repliseq_end) %>%
    do((function(z) {
      zz<<-z

      i.max = z$repliseq_fraction[which.max(z$repliseq_value_norm)]
      lb = which(z$repliseq_value_norm < th.repliseq_value_norm & z$repliseq_fraction < i.max)
      lb = z$repliseq_fraction[ifelse(length(lb), lb[length(lb)]+1, 1)]
      ub = which(z$repliseq_value_norm < th.repliseq_value_norm & z$repliseq_fraction > i.max)
      ub = z$repliseq_fraction[ifelse(length(ub), ub[1]-1, nrow(z))]

      z$repliseq_value_in_scope = dplyr::between(z$repliseq_fraction, lb, ub)
      z
    })(.)) %>%
    dplyr::ungroup()

  repliseq_time_df = repliseq_df %>%
    dplyr::mutate(repliseqTime_chrom=repliseq_chrom, repliseqTime_start=repliseq_start, repliseqTime_end=repliseq_end) %>%
    dplyr::group_by(repliseqTime_chrom, repliseqTime_start, repliseqTime_end) %>%
    dplyr::mutate(repliseqTime_avg=weighted.mean(repliseq_fraction[repliseq_value_in_scope], repliseq_value_norm[repliseq_value_in_scope], na.rm=T)) %>%
    dplyr::mutate(lb=dplyr::between(repliseq_fraction, floor(repliseqTime_avg[1]-4), ceiling(repliseqTime_avg[1]-2)), repliseqTime_min=ifelse(any(lb), weighted.mean(repliseq_fraction[lb], repliseq_value_norm[lb]), NA_real_)) %>%
    dplyr::mutate(ub=dplyr::between(repliseq_fraction, floor(repliseqTime_avg[1]+2), ceiling(repliseqTime_avg[1]+4)), repliseqTime_max=ifelse(any(ub), weighted.mean(repliseq_fraction[ub], repliseq_value_norm[ub]), NA_real_)) %>%
    dplyr::summarize(repliseqTime_avg=repliseqTime_avg[1], repliseqTime_min=repliseqTime_min[1], repliseqTime_max=repliseqTime_max[1], repliseqTime_lb=min(which(repliseq_value_in_scope)), repliseqTime_ub=max(which(repliseq_value_in_scope))) %>%
    dplyr::group_by(repliseqTime_chrom) %>%
    dplyr::mutate(repliseqTime_avg=smoother::smth.gaussian(repliseqTime_avg, window=window), repliseqTime_min=smoother::smth.gaussian(repliseqTime_min, window=window), repliseqTime_max=smoother::smth.gaussian(repliseqTime_max, window=window)) %>%
    dplyr::mutate(repliseqTime_avg=zoo::na.fill(repliseqTime_avg, "extend"), repliseqTime_min=zoo::na.fill(repliseqTime_min, "extend"), repliseqTime_max=zoo::na.fill(repliseqTime_max, "extend")) %>%
    dplyr::ungroup()

  repliseq_time_df
}

repliseq_preprocess = function() {
  # Load repliseq data
  repliseq_ESC_df = repliseq_read("data/zhao_bmc_repliseq_2020/GSE137764_mESC_Gaussiansmooth_scaled_autosome.mat") %>% dplyr::mutate(repliseq_celltype="esc")
  repliseq_NPC_df = repliseq_read("data/zhao_bmc_repliseq_2020/GSE137764_mNPC_Gaussiansmooth_scaled_autosome.mat") %>% dplyr::mutate(repliseq_celltype="npc")
  repliseq_ESC2NPC_df = dplyr::bind_rows(repliseq_ESC_df, repliseq_NPC_df)
  readr::write_tsv(repliseq_ESC2NPC_df,"data/zhao_bmc_repliseq_2020/preprocessed/repliseq.tsv")

  repliseqTime_ESC_df = repliseq_summarize(repliseq_ESC_df) %>% dplyr::mutate(repliseqTime_celltype="esc")
  repliseqTime_NPC_df = repliseq_summarize(repliseq_NPC_df) %>% dplyr::mutate(repliseqTime_celltype="npc")
  repliseqTime_ESC2NPC_df = dplyr::bind_rows(repliseqTime_ESC_df, repliseqTime_NPC_df)
  readr::write_tsv(repliseqTime_ESC2NPC_df,"data/zhao_bmc_repliseq_2020/preprocessed/repliseqTime.tsv")


  repliseq_merged_df = repliseq_ESC2NPC_df %>%
    dplyr::group_by(repliseq_chrom, repliseq_start, repliseq_end, repliseq_fraction) %>%
    dplyr::summarize(repliseq_value=mean(repliseq_value, na.rm=T)) %>%
    dplyr::arrange(repliseq_chrom, repliseq_start, repliseq_end, repliseq_fraction)
  repliseqTime_merged_df = repliseq_summarize(repliseq_merged_df) %>% dplyr::mutate(repliseqTime_celltype="esc/npc")
  readr::write_tsv(repliseqTime_merged_df,"data/zhao_bmc_repliseq_2020/preprocessed/repliseqTime_merged.tsv")
}

test = function() {
    dplyr::select(repliseqTime_chrom, repliseqTime_start, repliseqTime_end, repliseqTime_min.esc=repliseqTime_min, repliseqTime_max.esc=repliseqTime_max, repliseqTime_avg.esc=repliseqTime_avg) %>%
    dplyr::inner_join(repliseqTime_NPC_df %>% dplyr::select(repliseqTime_chrom, repliseqTime_start, repliseqTime_end, repliseqTime_min.npc=repliseqTime_min, repliseqTime_max.npc=repliseqTime_max, repliseqTime_avg.npc=repliseqTime_avg), by=c("repliseqTime_chrom", "repliseqTime_start", "repliseqTime_end")) %>%
    reshape2::melt(id.vars=c("repliseqTime_chrom", "repliseqTime_start", "repliseqTime_end"), value.name="repliseqTime_time", variable.name="repliseqTime_source") %>%
    dplyr::rename(repliseq_value.esc="repliseq_value") %>%
    dplyr::inner_join(repliseq_NPC_df %>% dplyr::rename(repliseq_value.npc="repliseq_value"), by=c("repliseq_chrom", "repliseq_start", "repliseq_end", "repliseq_fraction"))
  repliseq_cor = repliseqTime_NPC2ESC2RDC_df %>%
    dplyr::group_by(rdc_cluster, rdc_gene, rdc_cluster_details) %>%
    dplyr::summarize(R=cor(repliseq_time_avg.npc, repliseq_time_avg.esc)) %>%
    dplyr::arrange(R)

  ggplot(repliseq_cor) +
    geom_bar(aes(y=R, x=reorder(rdc_cluster_details, R)), stat="identity") +
    coord_flip()

  ggplot(repliseqTime_NPC2ESC2RDC_df %>% dplyr::filter(repliseq_chrom=="chr3")) +
    geom_point(aes(repliseq_time_avg.esc, repliseq_time_avg.npc, color=rdc_cluster_details)) +
    facet_wrap(~repliseq_chrom)
}

main = function() {
  m9_m10 = rtracklayer::import.chain("data/mm9/mm9ToMm10.over.chain")
  m10_m9 = rtracklayer::import.chain("data/mm10/mm10ToMm9.over.chain")

  chromosomes_map_df = readr::read_tsv("data/mm9_chromosomes_synonyms.tsv")

  # Load genome info
  chromsizes_cols = readr::cols(seqnames=col_character(), seqlengths=col_double())
  genome_info = with(readr::read_tsv("data/mm10/mm10.chrom.sizes", col_names=names(chromsizes_cols$cols), col_types=chromsizes_cols),
           GenomeInfoDb::Seqinfo(seqnames, seqlengths, isCircular=rep(F, length(seqnames)), genome=rep("mm10", length(seqnames))))
  genome_info = genome_info[paste0("chr", c(1:19, "X", "Y"))]

  genome_txdb = GenomicFeatures::makeTxDbFromGFF('data/mm10/mm10.refGene.gtf.gz', format="gtf")
  genes_df = as.data.frame(GenomicFeatures::genes(genome_txdb)) %>%
    dplyr::mutate(gene_chrom=as.character(seqnames), gene_start=start, gene_end=end, gene_strand=strand, gene_length=gene_end-gene_start)
  genes_ranges = GenomicRanges::makeGRangesFromDataFrame(genes_df, keep.extra.columns=T)

  genes_reduced_ranges = genes_ranges
  strand(genes_reduced_ranges) = "*"
  genes_reduced_ranges = GenomicRanges::reduce(genes_reduced_ranges)
  genes_reduced_ranges$gene_cluster = 1:length(genes_reduced_ranges)
  genes_df = as.data.frame(IRanges::mergeByOverlaps(genes_ranges, genes_reduced_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::group_by(gene_chrom, gene_cluster) %>%
    dplyr::summarize(gene_cluster_size=n(), gene_dominated_length=max(gene_end-gene_start), gene_id=paste0(gene_id, collapse=","), gene_start=min(gene_start), gene_end=max(gene_end), gene_strand="*", gene_length=gene_end-gene_start) %>%
    dplyr::select(-gene_dominated_length)

  chromsizes_cols = readr::cols(seqnames=col_character(), seqlengths=col_double())
  genome_info_mm9 = with(readr::read_tsv("data/mm9/mm9.chrom.sizes", col_names=names(chromsizes_cols$cols), col_types=chromsizes_cols),
             GenomeInfoDb::Seqinfo(seqnames, seqlengths, isCircular=rep(F, length(seqnames)), genome=rep("mm9", length(seqnames))))
  genome_info_mm9 = genome_info_mm9[paste0("chr", c(1:19, "X", "Y"))]

  #
  # Calculate normalization accross experiment
  #
  junctions_df = read_junctions("data/breaks_tlx/*/*.tlx") %>%
    dplyr::inner_join(chromosomes_map_df, by=c("bait_chrom"="chrom_synonym")) %>%
    dplyr::mutate(bait_chrom=unique_chrom) %>%
    dplyr::select(-unique_chrom) %>%
    dplyr::inner_join(chromosomes_map_df, by=c("junction_chrom"="chrom_synonym")) %>%
    dplyr::mutate(junction_chrom=unique_chrom) %>%
    dplyr::select(-unique_chrom) %>%
    dplyr::mutate(junction_name=stringr::str_glue("{bait}_{name}", bait=bait_chrom, name=junction_name)) %>%
    dplyr::filter(junction_chrom %in% seqlevels(genome_info_mm9)) %>%
    dplyr::filter(junction_chrom==bait_chrom)

  sample_df = junctions_df %>% dplyr::filter(experimental_condition=="Sample")
  sample_ranges = GenomicRanges::makeGRangesFromDataFrame(sample_df %>% dplyr::mutate(end=junction_end, start=junction_start, seqnames=junction_chrom), seqinfo=genome_info_mm9, keep.extra.columns=T)
  control_df = junctions_df %>% dplyr::filter(experimental_condition=="Control")
  control_ranges = GenomicRanges::makeGRangesFromDataFrame(control_df %>% dplyr::mutate(end=junction_end, start=junction_start, seqnames=junction_chrom), seqinfo=genome_info_mm9,  keep.extra.columns=T)
  params = call_islands_params(binsize=1e5, binstep=1e4, extend=1e5, llocal=2e6, minqvalue=0.01, maxgap=5e5, minlen=200)
  breaks_df = call_islands(sample_ranges, control_ranges, params, debug=F)$coverage
  breaks_ranges = GenomicRanges::makeGRangesFromDataFrame(breaks_df %>% dplyr::mutate(junction_chrom=seqnames, junction_start=start, junction_end=end, qvalue=qvalue.sample_vs_max), keep.extra.columns=T)
  breaks_df = liftOverRanges(breaks_ranges, m9_m10) %>%
    dplyr::mutate(junction_chrom=seqnames, junction_start=start, junction_end=end) %>%
    dplyr::select(-(seqnames:strand))

  # Load TAD
  tad_df = readr::read_tsv("data/HiC_bonev2017/topological_associated_domains.tsv") %>%
    dplyr::mutate(tad_celltype_num=as.numeric(as.factor(tad_celltype)))

  # Load CTCF
  ctcf_df = readr::read_tsv("data/CTCFBSDB.tsv") %>%
    dplyr::filter(ctcfdb_species=="Mouse")

  # Load RPKM (data fraom Peggy)
  rpkm_df = readr::read_tsv("data/rpkm/rpkm_data.tsv") %>%
    dplyr::mutate(rpkm_group=as.numeric(factor(rpkm_celltype))*10 + as.numeric(factor(rpkm_replicate))) %>%
    dplyr::group_by(rpkm_gene, rpkm_celltype, rpkm_chrom, rpkm_start, rpkm_end) %>%
    dplyr::summarize(rpkm_value=mean(rpkm_value)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(rpkm_celltype_num=as.numeric(as.factor(rpkm_celltype)))
  rpkm_ranges = GenomicRanges::makeGRangesFromDataFrame(rpkm_df %>% dplyr::mutate(seqnames=rpkm_chrom, start=rpkm_start, end=rpkm_end), keep.extra.columns=T)
  rpkm_df = liftOverRanges(rpkm_ranges, m9_m10) %>%
    dplyr::mutate(rpkm_chrom=seqnames, rpkm_start=start, rpkm_end=end) %>%
    dplyr::select(dplyr::matches("^rpkm_"))



  # Load transcription
  # Induction of recurrent break cluster genes in neuralprogenitor cells differentiated from embryonic stemcells in culture
  transcription_df = readr::read_tsv("data/tena2020/tena2020_mm10_transcription.tsv") %>%
    dplyr::group_by(tena2020transcription_gene, tena2020transcription_celltype, tena2020transcription_chrom, tena2020transcription_start, tena2020transcription_end) %>%
    dplyr::summarize(tena2020transcription_rpkm=mean(tena2020transcription_rpkm)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(tena2020transcription_celltype_num=as.numeric(as.factor(tena2020transcription_celltype)))
  # transcription_ranges = GenomicRanges::makeGRangesFromDataFrame(transcription_df %>% dplyr::mutate(seqnames=tena2020transcription_chrom, start=tena2020transcription_start, end=tena2020transcription_end), keep.extra.columns=T)
  # transcription_df = liftOverRanges(transcription_ranges, m9_m10) %>%
  #   dplyr::mutate(tena2020transcription_chrom=seqnames, tena2020transcription_start=start, tena2020transcription_end=end) %>%
  #   dplyr::select(dplyr::matches("^tena2020transcription_"))

  # Load RDC (Convert to mm10)
  rdc_cols = readr::cols(rdc_cluster=col_character(), rdc_chrom=col_character(), rdc_start=col_double(), rdc_end=col_double(), rdc_group=col_double(), rdc_gene=col_character())
  rdc_df = readr::read_tsv("data/rdc_pnas.tsv", col_types=rdc_cols) %>%
    dplyr::inner_join(chromosomes_map_df, by=c("rdc_chrom"="chrom_synonym")) %>%
    dplyr::mutate(rdc_chrom=unique_chrom) %>%
    dplyr::select(-unique_chrom) %>%
    dplyr::mutate(rdc_length=rdc_end-rdc_start) %>%
    dplyr::filter(rdc_length>2e5 & grepl("[a-zA-Z]", rdc_gene))
  rdc_ranges = GenomicRanges::makeGRangesFromDataFrame(rdc_df %>% dplyr::mutate(seqnames=rdc_chrom, start=rdc_start, end=rdc_end), keep.extra.columns=T)
  rdc_df = liftOverRanges(rdc_ranges, m9_m10) %>%
    dplyr::mutate(rdc_start=start, rdc_end=end, rdc_chrom=seqnames) %>%
    dplyr::select(dplyr::matches("^rdc_")) %>%
    dplyr::mutate(signal="RDC")

  # RDC genes
  re_genes = paste(rdc_df %>%
    tidyr::separate_rows(rdc_gene, sep=", *") %>%
    dplyr::distinct(rdc_gene) %>%
    dplyr::filter(!grepl("-", rdc_gene) & rdc_gene!="") %>%
    .$rdc_gene, collapse="|")

  # Load GC
  gc_df = readr::read_tsv("data/mm10/mm10_tile1000.freq", col_types=cols(gc_chrom=col_character(), gc_start=col_double(), gc_end=col_double(), gc_freq=col_double()))


  # Load repliseq
  # High-resolution Repli-Seq defines the temporal choreography of initiation, elongation and termination of replication in mammalian cells
  repliseq_df = readr::read_tsv("data/zhao_bmc_repliseq_2020/preprocessed/repliseq.tsv")
  repliseqTime_df = readr::read_tsv("data/zhao_bmc_repliseq_2020/preprocessed/repliseqTime.tsv")

  # Calculate IZ
  starts = unique(repliseqTime_df$repliseqTime_start)
  repliseqIZ_df = repliseqTime_df %>%
    reshape2::melt(id.vars=c("repliseqTime_celltype", "repliseqTime_chrom", "repliseqTime_start", "repliseqTime_end"), value.name="repliseqTime_time", variable.name="repliseqTime_source") %>%
    tidyr::extract(repliseqTime_source, "repliseqTime_source", "repliseqTime_([^.]+)") %>%
    dplyr::filter(repliseqTime_source %in% c("avg")) %>%
    dplyr::rename(iz_chrom="repliseqTime_chrom", iz_source="repliseqTime_source", iz_celltype="repliseqTime_celltype") %>%
    dplyr::group_by(iz_chrom, iz_celltype, iz_source) %>%
    dplyr::do((function(z){
      zz<<-z
      data.frame(iz_start=which(ggpmisc:::find_peaks(17-z$repliseqTime_time, ignore_threshold = 1/16, span=9, strict=F))) %>%
      dplyr::mutate(iz_start=z$repliseqTime_start[iz_start])
    })(.)) %>%
    dplyr::arrange(iz_start) %>%
    dplyr::group_by(iz_chrom) %>%
    dplyr::mutate(iz_cluster=dbscan::dbscan(as.matrix(iz_start), eps=50000*5, minPts=1)$cluster) %>%
    dplyr::group_by(iz_chrom, iz_cluster, iz_source) %>%
    dplyr::mutate(
      iz_start=starts[which.min(abs(starts-mean(iz_start)))],
      iz_source=paste(unique(iz_source), collapse=",")) %>%
    dplyr::inner_join(repliseqTime_df %>% dplyr::select(repliseqTime_chrom, repliseqTime_start, repliseqTime_celltype, repliseqTime_avg), by=c("iz_chrom"="repliseqTime_chrom", "iz_start"="repliseqTime_start", "iz_celltype"="repliseqTime_celltype")) %>%
    dplyr::distinct(iz_chrom, iz_start, iz_source, iz_celltype, .keep_all=T)
  repliseqIZmerged_df = repliseqIZ_df %>%
    reshape2::dcast(iz_chrom+iz_start+iz_source~iz_celltype, value.var="repliseqTime_avg")

  if(F) {
    pdf(file="reports/repliseq_termination.pdf", width=80, height=40)
    celltype_pal = c("npc"="#E41A1C", "esc"="#E97F02")
    repliseq_gg = repliseq_df %>% dplyr::filter(repliseq_celltype=="npc", grepl("chr6$", repliseq_chrom) & dplyr::between(repliseq_start, 6e7, 9e7)) %>% dplyr::mutate(chrom=repliseq_chrom)
    repliseqTime_gg = repliseqTime_df %>% dplyr::filter(grepl("chr6$", repliseqTime_chrom) & dplyr::between(repliseqTime_start, 6e7, 9e7)) %>% dplyr::mutate(chrom=repliseqTime_chrom)
    repliseqIZ_gg = repliseqIZ_df %>% dplyr::filter(iz_source=="avg" & grepl("chr6$", iz_chrom) & dplyr::between(iz_start, 6e7, 9e7)) %>% dplyr::mutate(chrom=iz_chrom)
    repliseqIZmerged_gg = repliseqIZmerged_df %>% dplyr::filter(iz_source=="avg" & grepl("chr6$", iz_chrom) & dplyr::between(iz_start, 6e7, 9e7)) %>% dplyr::mutate(chrom=iz_chrom)
    transcription_gg = transcription_df %>% dplyr::filter(grepl("chr6$", tena2020transcription_chrom) & dplyr::between(tena2020transcription_start, 6e7, 9e7)) %>% dplyr::mutate(chrom=tena2020transcription_chrom)

    ggplot() +
      geom_tile(aes(y=repliseq_fraction, x=repliseq_start, fill=repliseq_value), data=repliseq_gg) +
      geom_line(aes(x=repliseqTime_start, y=repliseqTime_avg, color=repliseqTime_celltype), data=repliseqTime_gg, size=1) +
      geom_segment(aes(x=iz_start, xend=iz_start, y=pmin(npc, esc), yend=pmax(npc, esc), color="npc"), size=0.5, linetype="dashed", data=repliseqIZmerged_gg) +
      geom_point(aes(x=iz_start, y=repliseqTime_avg, color=iz_celltype), size=3,  data=repliseqIZ_gg) +
      geom_rect(aes(xmin=tena2020transcription_start, xmax=tena2020transcription_end, ymin=-2*(tena2020transcription_celltype_num-1)-2, ymax=-2*(tena2020transcription_celltype_num-1), alpha=tena2020transcription_rpkm), fill="#FF8C00", data=transcription_gg) +
      scale_fill_gradientn(colours=c("#666666", "#CCCCCC", "#FFFFFF", "#00FFFF", "#000066"), values=c(0, 0.1, 0.3, 0.5, 1)) +
      scale_color_manual(values=celltype_pal)+
      facet_wrap(~chrom, ncol=1, scale="free_x")
    dev.off()
  }

  # Create random regions
  random_lengths = rdc_df %>% dplyr::select(rdc_chrom, rdc_length)
  random_df = random_lengths %>%
    dplyr::inner_join(as.data.frame(genome_info) %>% tibble::rownames_to_column("rdc_chrom"), by=c("rdc_chrom")) %>%
    dplyr::group_by(rdc_chrom) %>%
    dplyr::do((function(z){
      zz<<-z
      n = 5
      z.start = sample(z$seqlengths[1], nrow(z)*n)
      data.frame(rdc_chrom=z$rdc_chrom[1], seqlengths=z$seqlengths[1], start=z.start, rdc_length=rep(z$rdc_length, each=n)) %>%
        dplyr::mutate(rdc_start=ifelse(start+rdc_length>=seqlengths, start-rdc_length, start), rdc_end=start+z$rdc_length, signal="Random positions", rdc_gene="") %>%
        dplyr::mutate(rdc_cluster=paste0(rdc_chrom, "_", 1:n()))
    })(.)) %>%
    dplyr::ungroup()

  random2_df = genes_df %>%
    dplyr::filter(!grepl(re_genes, gene_id) & gene_length>2e5 & gene_chrom %in% paste0("chr", 1:19)) %>%
    dplyr::mutate(rdc_chrom=gene_chrom, rdc_cluster=paste0(rdc_chrom, "_", 1:n()), rdc_start=gene_start, rdc_end=gene_end, rdc_length=gene_length, signal="Non-RDC genes", rdc_gene=gene_id)
  rdc_df = dplyr::bind_rows(rdc_df, random2_df, random_df) %>% dplyr::select(rdc_chrom, rdc_start, rdc_end, rdc_length, rdc_cluster, rdc_gene, signal)
  rdcMargins_df = data.frame()
  for(m in c(0, 5e5, 1e6, 5e6)) {
    rdcMargins_df = dplyr::bind_rows(rdcMargins_df, rdc_df %>% dplyr::mutate(margin=m))
  }

  #
  # Merge different datasets
  #
  repliseqTime_ranges = GenomicRanges::makeGRangesFromDataFrame(repliseqTime_df %>% dplyr::mutate(seqnames=repliseqTime_chrom, start=repliseqTime_start, end=repliseqTime_start), keep.extra.columns=T)
  repliseqIZ_ranges = GenomicRanges::makeGRangesFromDataFrame(repliseqIZ_df %>% dplyr::mutate(seqnames=iz_chrom, start=iz_start, end=iz_start), keep.extra.columns=T)
  repliseqIZmerged_ranges = GenomicRanges::makeGRangesFromDataFrame(repliseqIZmerged_df %>% dplyr::mutate(seqnames=iz_chrom, start=iz_start, end=iz_start), keep.extra.columns=T)
  repliseq_ranges = GenomicRanges::makeGRangesFromDataFrame(repliseq_df %>% dplyr::mutate(seqnames=repliseq_chrom, start=repliseq_start, end=repliseq_start), keep.extra.columns=T)
  rdc_ranges = GenomicRanges::makeGRangesFromDataFrame(rdcMargins_df %>% dplyr::mutate(seqnames=rdc_chrom, start=rdc_start-margin, end=rdc_end+margin), keep.extra.columns=T)
  genes_ranges = GenomicRanges::makeGRangesFromDataFrame(genes_df %>% dplyr::mutate(seqnames=gene_chrom, start=gene_start, end=gene_end), keep.extra.columns=T)
  transcription_ranges = GenomicRanges::makeGRangesFromDataFrame(transcription_df %>% dplyr::mutate(seqnames=tena2020transcription_chrom, start=tena2020transcription_start, end=tena2020transcription_end), keep.extra.columns=T)
  breaks_ranges = GenomicRanges::makeGRangesFromDataFrame(breaks_df %>% dplyr::mutate(seqnames=junction_chrom, start=junction_start, end=junction_end, qvalue=qvalue.sample_vs_max), keep.extra.columns=T)
  ctcf_ranges = GenomicRanges::makeGRangesFromDataFrame(ctcf_df %>% dplyr::mutate(seqnames=ctcfdb_chrom, start=ctcfdb_start, end=ctcfdb_end), keep.extra.columns=T)
  tad_ranges = GenomicRanges::makeGRangesFromDataFrame(tad_df %>% dplyr::mutate(seqnames=tad_chrom, start=tad_start, end=tad_end), keep.extra.columns=T)
  gc_ranges = GenomicRanges::makeGRangesFromDataFrame(gc_df %>% dplyr::mutate(seqnames=gc_chrom, start=gc_start, end=gc_end), keep.extra.columns=T)
  # rpkm_ranges = GenomicRanges::makeGRangesFromDataFrame(rpkm_df %>% dplyr::mutate(seqnames=rpkm_chrom, start=rpkm_start, end=rpkm_end), keep.extra.columns=T)

  #
  # Merge different data frames
  #
  repliseq2rdc_df = as.data.frame(IRanges::mergeByOverlaps(repliseq_ranges, rdc_ranges)) %>% dplyr::select(-dplyr::matches("_ranges\\."))
  repliseqTime2rdc_df = as.data.frame(IRanges::mergeByOverlaps(repliseqTime_ranges, rdc_ranges)) %>% dplyr::select(-dplyr::matches("_ranges\\."))
  genes2rdc_df = as.data.frame(IRanges::mergeByOverlaps(genes_ranges, rdc_ranges)) %>% dplyr::select(-dplyr::matches("_ranges\\."))
  repliseqIZ2rdc_df = as.data.frame(IRanges::mergeByOverlaps(repliseqIZ_ranges, rdc_ranges)) %>% dplyr::select(-dplyr::matches("_ranges\\."))
  repliseqIZmerged2rdc_df = as.data.frame(IRanges::mergeByOverlaps(repliseqIZmerged_ranges, rdc_ranges)) %>% dplyr::select(-dplyr::matches("_ranges\\."))
  transcription2rdc_df = as.data.frame(IRanges::mergeByOverlaps(transcription_ranges, rdc_ranges)) %>% dplyr::select(-dplyr::matches("_ranges\\."))
  # rpkm2rdc_df = as.data.frame(IRanges::mergeByOverlaps(rpkm_ranges, rdc_ranges)) %>% dplyr::select(-dplyr::matches("_ranges\\."))
  breaks2rdc_df = as.data.frame(IRanges::mergeByOverlaps(breaks_ranges, rdc_ranges)) %>% dplyr::select(-dplyr::matches("_ranges\\."))
  ctcf2rdc_df = as.data.frame(IRanges::mergeByOverlaps(ctcf_ranges, rdc_ranges)) %>% dplyr::select(-dplyr::matches("_ranges\\."))
  tad2rdc_df = as.data.frame(IRanges::mergeByOverlaps(tad_ranges, rdc_ranges)) %>% dplyr::select(-dplyr::matches("_ranges\\."))

  repliseqIZ2rdc_ranges = GenomicRanges::makeGRangesFromDataFrame(repliseqIZ2rdc_df %>% dplyr::mutate(seqnames=rdc_chrom, start=rdc_start, end=rdc_end), keep.extra.columns=T)
  repliseqIZ2rdc2transcription_df = as.data.frame(IRanges::mergeByOverlaps(repliseqIZ2rdc_ranges, transcription_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::filter(iz_celltype==tena2020transcription_celltype) %>%
    dplyr::arrange(dplyr::desc(tena2020transcription_rpkm)) %>%
    dplyr::distinct(rdc_cluster, iz_celltype, iz_start, margin, signal, .keep_all=T)

  breaks2transcription_df = as.data.frame(IRanges::mergeByOverlaps(breaks_ranges, transcription_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::group_by(tena2020transcription_celltype, tena2020transcription_gene) %>%
    dplyr::summarize(tena2020transcription_rpkm=max(tena2020transcription_rpkm), qvalue=max(qvalue.sample_vs_max))
  breaks2transcription_df %>%
    dplyr::group_by(tena2020transcription_celltype) %>%
    dplyr::summarize(R=cor(qvalue, tena2020transcription_rpkm, method="spearman", use="pairwise.complete.obs"))


  #
  # T-SNE
  #
  long_genes_ranges = GenomicRanges::makeGRangesFromDataFrame(genes_df %>% dplyr::mutate(seqnames=gene_chrom, start=gene_start-2e5, end=gene_end+2e5) %>% dplyr::filter(gene_length>2e5), keep.extra.columns=T)
  long_genes_ranges$gene_rdc = countOverlaps(long_genes_ranges, rdc_ranges[rdc_ranges$signal=="RDC"]) > 1
  breaks2genes_df = as.data.frame(IRanges::mergeByOverlaps(long_genes_ranges, breaks_ranges)) %>%
    dplyr::group_by(gene_id, gene_chrom, gene_start, gene_end, gene_strand, gene_length, gene_rdc) %>%
    dplyr::summarize(hit=any(qvalue<1e-5))
  breaks2genes_ranges = GenomicRanges::makeGRangesFromDataFrame(breaks2genes_df %>% dplyr::mutate(seqnames=gene_chrom, start=gene_start-2e5, end=gene_end+2e5) %>% dplyr::filter(gene_length>2e5), keep.extra.columns=T)
  repliseq2genes_cluster_df = as.data.frame(IRanges::mergeByOverlaps(breaks2genes_ranges, repliseqTime_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::filter(repliseqTime_celltype=="npc") %>%
    dplyr::group_by(repliseqTime_celltype, gene_id, hit, gene_rdc) %>%
    dplyr::do((function(z){
      zz <<- z

      n = 50
      as.data.frame(matrix(approx(z$repliseqTime_start, z$repliseqTime_avg, n=n)$y, nrow=1, dimnames=list(row.names=1, col.names=paste0("part_", 1:n)) ))
    })(.)) %>%
    dplyr::ungroup()

  repliseq2genes_cluster_matrix = jitter(as.matrix(repliseq2genes_cluster_df %>% dplyr::select(dplyr::matches("^part_"))))
  repliseq2genes_cluster_tsne = Rtsne(repliseq2genes_cluster_matrix, perplexity=8)
  x = repliseq2genes_cluster_df %>% dplyr::mutate(x=repliseq2genes_cluster_tsne$Y[,1], y=repliseq2genes_cluster_tsne$Y[,2])
  # repliseq2genes_cluster_tsne = prcomp(repliseq2genes_cluster_matrix, scale=F)
  # x = repliseq2genes_cluster_df %>% dplyr::mutate(x=repliseq2genes_cluster_tsne$x[,1], y=repliseq2genes_cluster_tsne$x[,2])
  g = ggplot(x) +
    geom_point(aes(x=x, y=y, gene_id=gene_id, color=gene_rdc))
    # ggrepel::geom_text_repel(aes(x=x, y=y, label=gene_id))
  g
  # p = plotly::ggplotly(g)
  # htmlwidgets::saveWidget(plotly::as_widget(p), "index3.html")

  # x =repliseq2genes_cluster_df %>% dplyr::filter(grepl("Rnf220|Pitpnc1|Ctif|Atg7|Msi2|Pknox2|Srgap3", gene_id))
  # x =repliseq2genes_cluster_df %>% dplyr::filter(grepl("Sik3|Shroom3|Frmd5|Nav2|Megf11|Slc39a11|Trerf1", gene_id))
  # x =repliseq2genes_cluster_df %>% dplyr::filter(grepl("Cphx3|Tmem254c|Duxbl1|Cphx2", gene_id))
  # x =repliseq2genes_cluster_df %>% dplyr::filter(grepl("Ntm|Bmpr1b|Arhgap10|Stk39|Crim1A930003O13Rik|Slc22a22", gene_id))
  # x =repliseq2genes_cluster_df %>% dplyr::filter(grepl("Fam205a1|Mir6419|Myo16|Csmd2|Nckap5|Sh3rf3|Phactr1", gene_id))
  # x =repliseq2genes_cluster_df %>% dplyr::filter(grepl("Syne1Rbfox1|Adamts19|Cdh2Zfp385b|Rnf148|Cdkn2b|Pvt1|Klf12|Col22a1|Pik3c2g", gene_id))
  # x =repliseq2genes_cluster_df %>% dplyr::filter(grepl("Evi2a|Ppm1e|Pitpnc1", gene_id))
  # x = repliseq2genes_cluster_df %>% dplyr::filter(grepl("Dtna|Adamts3|Ntm|Mctp1|Ptprr|Rnf150|Dlgap1", gene_id))
  # g = ggplot(x) +
  #   geom_line(aes(x=repliseqTime_start, y=repliseqTime_avg)) +
  #   geom_vline(aes(xintercept=gene_start), linetype="dashed", data=x %>% dplyr::distinct(gene_id, .keep_all=T)) +
  #   geom_vline(aes(xintercept=gene_end), linetype="dashed", data=x %>% dplyr::distinct(gene_id, .keep_all=T)) +
  #   facet_wrap(~gene_id, scales="free_x")
  # g
  # p = plotly::ggplotly(g)
  # htmlwidgets::saveWidget(plotly::as_widget(p), "index4.html")




  #
  # RANDOM FOREST DATA
  # ==================================
  long_genes_ranges = GenomicRanges::makeGRangesFromDataFrame(genes_df %>% dplyr::mutate(seqnames=gene_chrom, start=gene_start-2e5, end=gene_end+2e5) %>% dplyr::filter(gene_length>2e5), keep.extra.columns=T)
  long_genes_ranges$gene_rdc = countOverlaps(long_genes_ranges, rdc_ranges[rdc_ranges$signal=="RDC"]) > 1

  repliseq2genes_df = as.data.frame(IRanges::mergeByOverlaps(long_genes_ranges, repliseqTime_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::group_by(repliseqTime_celltype, gene_id, gene_chrom, gene_start, gene_end, gene_length, gene_rdc) %>%
    dplyr::summarize(repliseqTime_slope=max(diff(abs(repliseqTime_avg)), na.rm=T), repliseqTime_max=max(repliseqTime_avg), repliseqTime_min=min(repliseqTime_avg), repliseqTime_diff=repliseqTime_max-repliseqTime_min)
  repliseq2genes_ranges = GenomicRanges::makeGRangesFromDataFrame(repliseq2genes_df %>% dplyr::mutate(seqnames=gene_chrom, start=gene_start, end=gene_end), keep.extra.columns=T)

  transcription2repliseq2genes_df = as.data.frame(IRanges::mergeByOverlaps(repliseq2genes_ranges, transcription_ranges)) %>%
    dplyr::filter(tena2020transcription_celltype==repliseqTime_celltype) %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::group_by(repliseqTime_celltype, gene_id, gene_chrom, gene_start, gene_end, gene_length, gene_rdc, repliseqTime_max, repliseqTime_min, repliseqTime_diff, repliseqTime_slope) %>%
    dplyr::summarize(tena2020transcription_rpkm=max(tena2020transcription_rpkm))
  transcription2repliseq2genes_ranges = GenomicRanges::makeGRangesFromDataFrame(transcription2repliseq2genes_df %>% dplyr::mutate(seqnames=gene_chrom, start=gene_start, end=gene_end), keep.extra.columns=T)

  breaks2transcription2repliseq2genes_df = as.data.frame(IRanges::mergeByOverlaps(breaks_ranges, transcription2repliseq2genes_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::filter(repliseqTime_celltype=="npc") %>%
    dplyr::mutate(qvalue=qvalue.sample_vs_max, qvalue_norm=-log10(qvalue), qvalue_norm=ifelse(qvalue_norm>20, 20, qvalue_norm)) %>%
    dplyr::group_by(repliseqTime_celltype, gene_id, gene_chrom, gene_start, gene_end, gene_length, gene_rdc, tena2020transcription_rpkm, repliseqTime_max, repliseqTime_min, repliseqTime_diff, repliseqTime_slope) %>%
    dplyr::summarize(qvalue_norm=max(qvalue_norm, na.rm=T))
  breaks2transcription2repliseq2genes_ranges = GenomicRanges::makeGRangesFromDataFrame(breaks2transcription2repliseq2genes_df %>% dplyr::mutate(seqnames=gene_chrom, start=gene_start, end=gene_end), keep.extra.columns=T)

  breaks2transcription2repliseq2genes2tad_df = as.data.frame(IRanges::mergeByOverlaps(breaks2transcription2repliseq2genes_ranges, tad_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::mutate(tad_celltype==repliseqTime_celltype) %>%
    dplyr::group_by(repliseqTime_celltype, gene_id, gene_chrom, gene_start, gene_end, gene_length, gene_rdc, tena2020transcription_rpkm, repliseqTime_max, repliseqTime_min, repliseqTime_diff, repliseqTime_slope, qvalue_norm) %>%
    dplyr::summarize(A_TAD_all=all(tad_compartment=="A"), B_TAD_all=all(tad_compartment=="B"), A_TAD_any=any(tad_compartment=="A"), B_TAD_any=any(tad_compartment=="B")) %>%
    dplyr::ungroup()
  breaks2transcription2repliseq2genes2tad_ranges = GenomicRanges::makeGRangesFromDataFrame(breaks2transcription2repliseq2genes2tad_df %>% dplyr::mutate(seqnames=gene_chrom, start=gene_start-5e5, end=gene_end+5e5), keep.extra.columns=T)

  breaks2transcription2repliseq2genes2tad2IZ_df = as.data.frame(IRanges::mergeByOverlaps(breaks2transcription2repliseq2genes2tad_ranges, repliseqIZmerged_ranges)) %>%
    dplyr::filter(!is.na(npc)) %>%
    # dplyr::mutate(iz_celltype==repliseqTime_celltype) %>%
    dplyr::group_by(gene_id, gene_chrom, gene_start, gene_end, gene_length, gene_rdc, tena2020transcription_rpkm, repliseqTime_max, repliseqTime_min, repliseqTime_diff, repliseqTime_slope, qvalue_norm, A_TAD_all, B_TAD_all, A_TAD_any, B_TAD_any, .drop=F) %>%
    dplyr::summarize(IZ_count=n(), IZ_dist=min(abs(gene_start+gene_length/2 - iz_start)/gene_length), IZ_outside=sum(iz_start>gene_end | iz_start<gene_start), IZ_inside=sum(iz_start>gene_end | iz_start<gene_start)) %>%
    dplyr::ungroup()

  breaks2transcription2repliseq2genes2tad2IZ_ranges = GenomicRanges::makeGRangesFromDataFrame(breaks2transcription2repliseq2genes2tad2IZ_df %>% dplyr::mutate(seqnames=gene_chrom, start=gene_start-5e5, end=gene_end+5e5), keep.extra.columns=T)
  breaks2transcription2repliseq2genes2tad2IZ2gc_df = as.data.frame(IRanges::mergeByOverlaps(breaks2transcription2repliseq2genes2tad2IZ_ranges, gc_ranges)) %>%
    dplyr::group_by(gene_id, gene_chrom, gene_start, gene_end, gene_length, gene_rdc, tena2020transcription_rpkm, repliseqTime_max, repliseqTime_min, repliseqTime_diff, repliseqTime_slope, qvalue_norm, A_TAD_all, B_TAD_all, A_TAD_any, B_TAD_any, IZ_count, IZ_dist, IZ_outside, IZ_inside, .drop=F) %>%
    dplyr::summarize(gc_freq=max(gc_freq, na.rm=T)) %>%
    dplyr::ungroup()


  load("data/methylation/data_all_npc.rda")
  breaks2transcription2repliseq2genes2tad2IZ2gc_ranges = GenomicRanges::makeGRangesFromDataFrame(breaks2transcription2repliseq2genes2tad2IZ2gc_df %>% dplyr::mutate(seqnames=gene_chrom, start=gene_start-5e5, end=gene_end+5e5), keep.extra.columns=T)
  methylation_ranges = GenomicRanges::makeGRangesFromDataFrame(methylation_npc_df %>% dplyr::mutate(seqnames=methylation_chrom, start=methylation_start, end=methylation_end), keep.extra.columns=T)
  breaks2transcription2repliseq2genes2tad2IZ2gc2meth_df = as.data.frame(IRanges::mergeByOverlaps(breaks2transcription2repliseq2genes2tad2IZ2gc_ranges, methylation_ranges)) %>%
    dplyr::group_by(gene_id, gene_chrom, gene_start, gene_end, gene_length, gene_rdc, tena2020transcription_rpkm, repliseqTime_max, repliseqTime_min, repliseqTime_diff, repliseqTime_slope, qvalue_norm, A_TAD_all, B_TAD_all, A_TAD_any, B_TAD_any, IZ_count, IZ_dist, IZ_outside, IZ_inside, gc_freq, methylation_celltype, methylation_type, .drop=F) %>%
    dplyr::summarize(methylation_maxscore=max(methylation_score, na.rm=T)) %>%
    dplyr::ungroup() %>%
    reshape2::dcast(...~methylation_type, value.var="methylation_maxscore")


  #
  # Random forest
  #
  # prediction_cols = c("repliseqTime_avg", "A_TAD_all", "B_TAD_all", "A_TAD_any", "B_TAD_any", "tena2020transcription_rpkm")
  # prediction_cols = c("repliseqTime_avg", "tena2020transcription_rpkm", "gene_length")
  # prediction_cols = c("repliseqTime_max", "repliseqTime_min", "repliseqTime_diff", "repliseqTime_slope", "tena2020transcription_rpkm", "gene_length", "A_TAD_all", "B_TAD_all", "A_TAD_any", "B_TAD_any", "IZ_count", "IZ_dist", "IZ_outside", "IZ_inside", "gc_freq", "H3K27me3", "H3K4me3", "H3K9me3")
  prediction_cols = c("repliseqTime_diff", "repliseqTime_slope", "tena2020transcription_rpkm", "gene_length", "IZ_dist", "gc_freq", "H3K27me3", "H3K4me3", "H3K9me3")

  all_df = breaks2transcription2repliseq2genes2tad2IZ2gc2meth_df
  all_df = all_df[rowSums(is.na(all_df[,setdiff(prediction_cols, col)]))==0,]

  perf_results = data.frame()
  varimp = data.frame()
  for(col in c("all", prediction_cols)) {
    # all_df$y = factor(ifelse(all_df$qvalue_norm>=3, "breaks", "no breaks"))
    all_df$y = factor(ifelse(all_df$gene_rdc, "breaks", "no breaks"))
    for(i in 1:50) {
      training_size = round(sum(all_df$y=="breaks")*0.8)
      training_hit_i = sample(which(all_df$y=="breaks"), training_size)
      training_non_i = sample(setdiff(1:nrow(all_df), training_hit_i), training_size)
      data_training = all_df[c(training_hit_i, training_non_i),]

      testing_hit_i = setdiff(which(all_df$y=="breaks"), training_hit_i)
      testing_non_i = sample(setdiff(which(all_df$y!="breaks"), training_non_i), length(testing_hit_i))
      data_testing = all_df[c(testing_hit_i, testing_non_i),]
      model = randomForest(y=data_training$y, x=data_training[,setdiff(prediction_cols, col)], importance=T)

      # library(party)
      # model = party::ctree(y ~ repliseqTime_slope, data=data_training, controls=cforest_control(mtry=2, mincriterion=0))
      # plot(model)

      y_ = predict(model, data_testing[,setdiff(prediction_cols, col)], type = "prob")
      pred <- ROCR::prediction(as.matrix(y_[,2]), as.matrix(as.numeric(data_testing$y)))
      perf <- ROCR::performance(pred, "tpr", "fpr")
      auc <- ROCR::performance(pred, "auc")@y.values[[1]][1]
      perf_results = dplyr::bind_rows(perf_results, data.frame(TPR=perf@y.values[[1]], FPR=perf@x.values[[1]], auc=auc, missing_column=col, i=i))

      if(col=="all") {
        varimp_i = as.data.frame(randomForest::importance(model)) %>%
          tibble::rownames_to_column("variable") %>%
          dplyr::mutate(i=i, MeanDecreaseAccuracy=MeanDecreaseAccuracy/max(MeanDecreaseAccuracy), MeanDecreaseGini=MeanDecreaseGini/max(MeanDecreaseGini)) %>%
          dplyr::select(variable, i, MeanDecreaseAccuracy, MeanDecreaseGini)

        varimp_i$R = NA_real_
        for(col2 in prediction_cols) {
          cor_i = cor(data_testing[[col2]], data_testing$qvalue_norm, method="spearman", use="pairwise.complete.obs")
          varimp_i$R[varimp_i$variable==col2] = cor_i
        }

        varimp = rbind(varimp, varimp_i %>% reshape2::melt(id.vars=c("i", "variable"), variable.name="metric"))

      }
    }
  }

  varimp.sorted = varimp %>%
    dplyr::filter(metric=="MeanDecreaseGini") %>%
    dplyr::group_by(variable, metric) %>%
    dplyr::summarize(value.sd=sd(value), value.mean=mean(value)) %>%
    dplyr::arrange(value.mean) %>%
    dplyr::distinct(dplyr::desc(variable)) %>%
    .$variable
  varimp = varimp %>%
    dplyr::mutate(variable=factor(variable, varimp.sorted)) %>%
    dplyr::arrange(variable)

  perf_results.sum = perf_results %>%
    dplyr::group_by(missing_column, FPR) %>%
    dplyr::summarize(TPR=mean(TPR, na.rm=T)) %>%
    dplyr::group_by(missing_column) %>%
    dplyr::mutate(TPR=smoother::smth.gaussian(TPR, window=10)) %>%
    dplyr::ungroup()

  pdf(file="reports/random_forest_new3.pdf", width=8, height=6)
  ggplot(varimp) +
    geom_hline(yintercept=0, color="#000000") +
    geom_boxplot(aes(x=variable, y=value)) +
    coord_flip() +
    facet_wrap(~metric)

  # randomForest::varImpPlot(models[["all"]])
  perf_results.ggplot = perf_results.sum
  ggplot(perf_results.ggplot) +
    geom_line(aes(x=FPR, y=TPR, color=missing_column), size=ifelse(perf_results.ggplot$missing_column=="all", 2, 0.5)) +
    geom_abline(slope=1) +
    coord_cartesian(xlim=c(0, 1), ylim=c(0,1))
  dev.off()



  model = randomForest(y=all_df$y, x=all_df[,prediction_cols], importance=T)
  y_ = predict(model, all_df[,prediction_cols], type = "prob")
  table(all_df$y, y_[,1]>y_[,2])

  transcription2repliseqTime_df = as.data.frame(IRanges::mergeByOverlaps(transcription_ranges, repliseqTime_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::mutate(tena2020transcription_length=tena2020transcription_end-tena2020transcription_start) %>%
    dplyr::left_join(transcription2rdc_df %>% dplyr::filter(margin==0) %>% dplyr::distinct(tena2020transcription_celltype, rdc_cluster, tena2020transcription_gene, signal), by=c("tena2020transcription_celltype", "tena2020transcription_gene")) %>%
    dplyr::inner_join(tad2rdc_df %>% dplyr::filter(margin==0) %>% dplyr::distinct(tad_celltype, tad_compartment, rdc_cluster, signal), by=c("tena2020transcription_celltype"="tad_celltype", "rdc_cluster", "signal")) %>%
    dplyr::filter(repliseqTime_celltype==tena2020transcription_celltype & !is.na(signal) & tena2020transcription_length>1e5) %>%
    dplyr::group_by(tena2020transcription_celltype, signal, tena2020transcription_gene, tena2020transcription_rpkm, tena2020transcription_length) %>%
    dplyr::summarize(repliseqTime_avg=max(repliseqTime_avg), tad_A=ifelse(any(tad_compartment=="A"), "A-TAD", "No A-TAD"), tad_B=ifelse(any(tad_compartment=="B"), "B-TAD", "No B-TAD")) %>%
    dplyr::ungroup()

  transcription2repliseqTime_df %>%
    dplyr::group_by(tena2020transcription_celltype, signal) %>%
    dplyr::summarize(cor=cor(tena2020transcription_rpkm, repliseqTime_avg, method="spearman", use="pairwise.complete.obs"))

  pdf(file="reports/transcription_vs_rpkm.pdf", width=8, height=6)
  transcription2repliseqTime_df %>%
    dplyr::filter(tena2020transcription_celltype=="npc") %>%
    # dplyr::filter(rpkm_value<5) %>%
    dplyr::mutate(tena2020transcription_rpkm=cut(tena2020transcription_rpkm, 10), tena2020transcription_length=cut(tena2020transcription_length, 5)) %>%
    ggplot() +
    geom_boxplot(aes(x=tena2020transcription_rpkm, y=repliseqTime_avg, fill=signal))
    # geom_beeswarm(aes(x=tena2020transcription_rpkm, y=repliseqTime_avg, group=signal), dodge.width=0.8, size=0.5) +
    # facet_wrap(~tad_A)
  dev.off()


  if(F)
  {
    pdf(file="reports/ctcf_counts.pdf", width=10, height=8)
    x = ctcf2rdc_df %>%
      dplyr::group_by(rdc_cluster, margin, signal) %>%
      dplyr::summarize(CTCF=1e6*n()/(rdc_length+margin*2)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(margin) %>%
      dplyr::mutate(margin_text=paste0(round(margin/1e6, 1), "M")) %>%
      dplyr::mutate(margin_text=factor(margin_text, unique(margin_text)))
    ggplot(x) +
      geom_boxplot(aes(x=margin_text, y=CTCF, fill=signal)) +
      labs(y="Experimenta CTCF count per million", x="") +
      theme_bw(base_size = 30) +
      theme(legend.position="bottom", axis.title.x=element_blank(), legend.text=element_text(size=20),legend.title=element_text(size=20))
    dev.off()

    pdf(file="reports/tad_counts.pdf", width=10, height=8)
    x = tad2rdc_df %>%
      dplyr::group_by(tad_celltype, tad_compartment, rdc_cluster, margin, signal) %>%
      dplyr::summarize(TAD=1e6*n()/(rdc_length+margin*2)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(margin) %>%
      dplyr::mutate(margin_text=paste0(round(margin/1e6, 1), "M")) %>%
      dplyr::mutate(margin_text=factor(margin_text, unique(margin_text)))
    ggplot(x) +
      geom_boxplot(aes(x=margin_text, y=TAD, fill=signal)) +
      labs(y="Experimenta TAD count per million", x="") +
      theme_bw(base_size = 30) +
      theme(legend.position="bottom", axis.title.x=element_blank(), legend.text=element_text(size=20),legend.title=element_text(size=20)) +
      facet_wrap(~tad_celltype+tad_compartment, scales="free")
    dev.off()

    pdf(file="reports/iz_counts.pdf", width=10, height=8)
    x = repliseqIZ2rdc2transcription_df %>%
      dplyr::filter(signal=="RDC") %>%
      dplyr::group_by(rdc_cluster, signal, margin) %>%
      dplyr::mutate(transcribed=paste0(unique(iz_celltype[tena2020transcription_rpkm>=0.1]), collapse="/"))  %>%
      dplyr::mutate(transcribed=ifelse(transcribed=="", "none", transcribed)) %>%
      dplyr::group_by(rdc_cluster, iz_celltype, signal, margin, transcribed) %>%
      dplyr::summarize(IZ=1e6*n()/rdc_length[1]) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(margin) %>%
      dplyr::mutate(margin_text=paste0(round(margin/1e6, 1), "M")) %>%
      dplyr::mutate(margin_text=factor(margin_text, unique(margin_text)))
    ggplot(x) +
      geom_boxplot(aes(x=margin_text, y=IZ, fill=transcribed), outlier.alpha=0) +
      geom_beeswarm(aes(x=margin_text, y=IZ, group=transcribed), dodge.width=0.8, size=0.5) +
      labs(y="Repliseq peaks (per Mb)", x="") +
      theme_bw(base_size = 30) +
      facet_wrap(~iz_celltype, scales="free_y") +
      theme(legend.position="bottom", axis.title.x=element_blank(), legend.text=element_text(size=20),legend.title=element_text(size=20))


    y = x %>%
      reshape2::dcast(rdc_cluster+signal+margin_text+transcribed~iz_celltype, value.var="IZ")
    ggplot(y) +
      geom_point(aes(x=esc, y=npc, color=transcribed)) +
      geom_abline(slope=1) +
      theme_bw(base_size = 30) +
      labs(x="Repliseq peaks in ESC (per Mb)", y="Repliseq peaks in NPC (per Mb)") +
      facet_wrap(~margin_text)


    x = repliseqIZ2rdc_df %>%
      dplyr::group_by(rdc_cluster, iz_celltype, signal, margin) %>%
      dplyr::summarize(IZ=1e6*n()/rdc_length) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(margin) %>%
      dplyr::mutate(margin_text=paste0(round(margin/1e6, 1), "M")) %>%
      dplyr::mutate(margin_text=factor(margin_text, unique(margin_text)))
    ggplot(x) +
      geom_boxplot(aes(x=margin_text, y=IZ, fill=signal)) +
      labs(y="Repliseq peaks (per Mb)", x="") +
      theme_bw(base_size = 30) +
      facet_wrap(~iz_celltype) +
      theme(legend.position="bottom", axis.title.x=element_blank(), legend.text=element_text(size=20),legend.title=element_text(size=20))
    dev.off()
  }

  pdf(file="reports/repliseq_to_rdc_new4.pdf", width=25, height=60)
  re_genes = "Ctnna"
  # re_genes = "Lsamp|Ctnna|Csmd1"
  # re_genes = "1700048|O20Rik|Ackr2|Astn2|Auts2|Bai3|Ccser1|Cdh13|Celf4|Csmd1|Ctnna2|Ctnnd2|Dcc|Dock1|Fhit|Gpc6|Grid2|Grik2|Hdac9|Lrp1b|Magi2|Naaladl2|Nbea|Nkain2|Nrg3|Nrxn1|Pard3b|Pcdh9|Prkg1|Rbfox1|Rbms3|Sdk1|Tenm4|Tpgs2|Wwox"
  # re_genes = "[A-Za-z]*"

  repliseqIZ2rdc_gg = repliseqIZ2rdc_df %>% dplyr::filter(margin==1e6 & signal=="RDC" & grepl(re_genes, rdc_gene))
  repliseqIZmerged2rdc_gg = repliseqIZmerged2rdc_df %>% dplyr::filter(margin==1e6 & signal=="RDC" & grepl(re_genes, rdc_gene))
  repliseq2rdc_gg = repliseq2rdc_df %>% dplyr::filter(margin==1e6 & signal=="RDC" & grepl(re_genes, rdc_gene)) %>% dplyr::mutate(repliseq_celltype=="npc")
  repliseqTime2rdc_gg = repliseqTime2rdc_df %>% dplyr::filter(margin==1e6 & signal=="RDC" & grepl(re_genes, rdc_gene))
  rdc_gg = rdcMargins_df %>% dplyr::filter(margin==1e6 & signal=="RDC" & grepl(re_genes, rdc_gene))
  genes2rdc_gg = genes2rdc_df %>% dplyr::filter(margin==1e6 & signal=="RDC" & grepl(re_genes, gene_id))
  transcription2rdc_gg = transcription2rdc_df %>% dplyr::filter(margin==1e6 & signal=="RDC" & grepl(re_genes, rdc_gene)) %>%
    dplyr::mutate(tena2020transcription_rpkm=ifelse(tena2020transcription_rpkm>=20, 20, tena2020transcription_rpkm))
  # rpkm2rdc_gg = rpkm2rdc_df %>% dplyr::filter(margin==1e6 & signal=="RDC" & grepl(re_genes, rdc_gene))
  ctcf2rdc_gg = ctcf2rdc_df %>% dplyr::filter(margin==1e6 & signal=="RDC" & grepl(re_genes, rdc_gene))
  tad2rdc_gg = tad2rdc_df %>% dplyr::filter(margin==1e6 & signal=="RDC" & grepl(re_genes, rdc_gene))
  breaks2rdc_gg = breaks2rdc_df %>% dplyr::filter(margin==1e6 & signal=="RDC" & grepl(re_genes, rdc_gene)) %>%
    dplyr::mutate(qvalue_norm = -log10(qvalue), ifelse(!is.finite(qvalue_norm), max(qvalue_norm[is.finite(qvalue_norm)]), qvalue_norm)) %>%
    dplyr::group_by(rdc_cluster) %>%
    dplyr::mutate(qvalue_max=max(qvalue_norm), qvalue_norm=qvalue_norm/qvalue_max) %>%
    dplyr::filter(junction_end-junction_start < 15e3) %>%
    dplyr::ungroup()

  celltype_pal = c("npc"="#E41A1C", "esc"="#E97F02")
  # , alpha=breaks2rdc_gg$qvalue_max/max(breaks2rdc_gg$qvalue_max, na.rm=T)
  ggplot() +
    geom_line(aes(x=junction_start, y=17+qvalue_norm*5), data=breaks2rdc_gg, size=1, alpha=breaks2rdc_gg$qvalue_max)  +
    geom_tile(aes(y=repliseq_fraction, x=repliseq_start, fill=repliseq_value), data=repliseq2rdc_gg) +
    geom_line(aes(x=repliseqTime_start, y=repliseqTime_avg, color=repliseqTime_celltype), data=repliseqTime2rdc_gg, size=1) +
    geom_ribbon(aes(x=repliseqTime_start, ymin=repliseqTime_min, ymax=repliseqTime_max, color=repliseqTime_celltype), data=repliseqTime2rdc_gg, alpha=0.1) +
    geom_segment(aes(x=iz_start, xend=iz_start, y=pmin(npc, esc), yend=pmax(npc, esc), color="npc"), size=0.2, linetype="twodash", data=repliseqIZmerged2rdc_gg) +
    geom_point(aes(x=iz_start, y=repliseqTime_avg, color=iz_celltype), size=3,  data=repliseqIZ2rdc_gg) +
    geom_rect(aes(xmin=rdc_start, xmax=rdc_end, ymin=-2, ymax=0), data=rdc_gg) +
    geom_text(aes(x=(rdc_end+rdc_start)/2, y=-1), label="rdc", data=rdc_gg) +
    # geom_rect(aes(xmin=ctcfdb_start, xmax=ctcfdb_end, ymin=-4, ymax=-2), data=ctcf2rdc_gg, color="#E7298A") +
    geom_rect(aes(xmin=tad_start, xmax=tad_end, ymin=-5+tad_celltype_num, ymax=-4+tad_celltype_num), fill=c("A"="#FF0000","B"="#0000FF")[tad2rdc_gg$tad_compartment], color="#FFFFFF00", data=tad2rdc_gg, color="#E7298A") +
    geom_rect(aes(xmin=gene_start, xmax=gene_end, ymin=-6, ymax=-4), fill="#FF8C00", color="#000000", data=genes2rdc_gg) +
    geom_text(aes(x=(gene_start+gene_end)/2, y=-5, label=gene_id), data=genes2rdc_gg, size=5) +

    geom_rect(aes(xmin=tena2020transcription_start, xmax=tena2020transcription_end, ymin=-9+tena2020transcription_celltype_num, ymax=-8+tena2020transcription_celltype_num, alpha=tena2020transcription_rpkm), fill="#377EB8", data=transcription2rdc_gg) +
    geom_text(aes(x=rdc_start, y=-8.5+tena2020transcription_celltype_num, label=paste(tena2020transcription_celltype)), data=transcription2rdc_gg %>% dplyr::distinct(rdc_cluster, tena2020transcription_celltype, .keep_all=T)) +

    # geom_rect(aes(xmin=rpkm_start, xmax=rpkm_end, ymin=-11+rpkm_celltype_num, ymax=-10+rpkm_celltype_num, alpha=rpkm_value), fill="#377EB8", data=rpkm2rdc_gg) +
    # geom_text(aes(x=rdc_start, y=-10.5+rpkm_celltype_num, label=paste(rpkm_celltype, "(RPKM)")), data=rpkm2rdc_gg %>% dplyr::distinct(rdc_cluster, rpkm_celltype, .keep_all=T)) +

    scale_fill_gradientn(colours=c("#666666", "#CCCCCC", "#FFFFFF", "#00FFFF", "#000066"), values=c(0, 0.1, 0.3, 0.5, 1)) +
    scale_color_manual(values=celltype_pal) +
    scale_x_continuous(breaks=scale_1mb) +
    coord_cartesian(ylim=c(-9, 22)) +
    facet_wrap(rdc_cluster~., scales="free_x", ncol=5, strip.position="left")


  # +
  #     facet_wrap(celltype~., scales="free_x", ncol=5, strip.position="left")
  # p = plotly::ggplotly(g)
  # htmlwidgets::saveWidget(plotly::as_widget(p), "index.html")
  dev.off()
}