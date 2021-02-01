library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
source("src/visualization.R")
source("src/islands.R")

main = function() {
  chromsizes_cols = readr::cols(seqnames=col_character(), seqlengths=col_double())
  genome_info = with(readr::read_tsv("data/mm10/mm10.chrom.sizes", col_names=names(chromsizes_cols$cols), col_types=chromsizes_cols),
             GenomeInfoDb::Seqinfo(seqnames, seqlengths, isCircular=rep(F, length(seqnames)), genome=rep("mm9", length(seqnames))))
  genome_info = genome_info[paste0("chr", c(1:19, "X", "Y"))]

  sample_annotations = readr::read_tsv("data/vivien/sample_annotations.tsv")
  junctions_df = read_junctions("data/vivien/chr*/*.tlx") %>%
    dplyr::mutate(library=gsub("_.*", "", basename(junction_file))) %>%
    dplyr::inner_join(sample_annotations, by=c("library", "bait_chrom")) %>%
    dplyr::mutate(experimental_condition=ifelse(treatment=="DMSO", "Control", "Sample")) %>%
    dplyr::select(-treatment)
  junctions_bait_df = junctions_df %>%
    dplyr::filter(bait_chrom==junction_chrom)


  scale_factor_df = junctions_bait_df %>%
    dplyr::group_by(junction_file) %>%
    dplyr::summarise(libsize=n(), .groups="keep") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(libsize_median=median(libsize)) %>%
    dplyr::mutate(score=libsize_median/libsize) %>%
    dplyr::select(junction_file, score)
  junctions_bait_df = junctions_bait_df %>%
    dplyr::select(-dplyr::matches("score")) %>%
    dplyr::inner_join(scale_factor_df, by="junction_file")

  params = call_islands_params(binsize=1e5, binstep=1e4, extend=1e5, llocal=2e6, minqvalue=0.01, maxgap=5e5, minlen=200)
  results = list()
  for(gene in c("none", "csmd1", "ctnna2")) {
    sample_df = junctions_bait_df %>% dplyr::filter(experimental_condition=="Sample" & (is.na(gene_silenced) & gene=="none" | gene_silenced==gene))
    sample_ranges = GenomicRanges::makeGRangesFromDataFrame(sample_df %>% dplyr::mutate(end=junction_end, start=junction_start, seqnames=junction_chrom), seqinfo=genome_info, keep.extra.columns=T)
    control_df = junctions_bait_df %>% dplyr::filter(experimental_condition=="Control" & (is.na(gene_silenced) & gene=="none" | gene_silenced==gene))
    control_ranges = GenomicRanges::makeGRangesFromDataFrame(control_df %>% dplyr::mutate(end=junction_end, start=junction_start, seqnames=junction_chrom), seqinfo=genome_info,keep.extra.columns=T)
    results[[gene]] = call_islands(
      sample_ranges=unreated_sample_ranges,
      control_ranges=unreated_control_ranges,
      genome_info=genome_info,
      params=params,
      debug=F)
  }

    # sample_islands_df2offtargets = islands2offtargets_identity(results$islands$sample, targets_df, genome_mm9)
    pdf(file="reports/duo_baseline_extend1e5_smooth2e6_bin1e5_step1e4_3.pdf", width=15, height=8)
    for(chr in paste0("chr", 1:19)) {
      coverage_df = dplyr::bind_rows(
        results$coverage$control %>% dplyr::mutate(experimental_condition="Control"),
        results$coverage$sample %>% dplyr::mutate(experimental_condition="Sample")) %>%
        dplyr::filter(seqnames %in% chr)
      baseline_df = dplyr::bind_rows(
        # results$baseline$control %>% dplyr::mutate(experimental_condition="Control"),
        # results$baseline$sample %>% dplyr::mutate(experimental_condition="Sample"),
        results$baseline$mixed %>% dplyr::mutate(experimental_condition="Max")) %>%
         dplyr::filter(seqnames %in% chr)
      islands_df = dplyr::bind_rows(
        results$islands$control %>% dplyr::mutate(experimental_condition="Control"),
        results$islands$sample %>% dplyr::mutate(experimental_condition="Sample")) %>%
        dplyr::filter(island_chrom %in% chr)
      karyoplot.colors=c(enriched="#0000FF", Control="#000000", Sample="#FF0000", Offtarget="#00FF00", "Sample filtered"="#FF0000", Max="#0000FF", Coverage="#000000", RDC="#FFFF00")
      p = ggplot_karyoplot(
        coverage_df=coverage_df,
        baseline_df=baseline_df,
        islands_df=islands_df,
        colors=karyoplot.colors
      )
      print(p)

      results$islands$sample %>% dplyr::filter(dplyr::between(island_start, 25e6, 27e6))
    }
    dev.off()

  head(p$sample)
}