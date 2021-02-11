library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(GenomicFeatures)
library(GenomicRanges)
library(gridExtra)
source("src/visualization.R")
source("src/islands.R")

main = function() {
  genome_txdb = GenomicFeatures::makeTxDbFromGFF('data/mm9/mm9.refGene.gtf.gz', format="gtf")
  genes_df = as.data.frame(genes(genome_txdb)) %>%
    dplyr::mutate(seqnames=as.character(seqnames))

  chromsizes_cols = readr::cols(seqnames=col_character(), seqlengths=col_double())
  genome_info = with(readr::read_tsv("data/mm10/mm10.chrom.sizes", col_names=names(chromsizes_cols$cols), col_types=chromsizes_cols),
             GenomeInfoDb::Seqinfo(seqnames, seqlengths, isCircular=rep(F, length(seqnames)), genome=rep("mm9", length(seqnames))))
  genome_info = genome_info[paste0("chr", c(1:19, "X", "Y"))]

  sample_annotations = readr::read_tsv("data/vivien/sample_annotations.tsv")
  junctions_df = read_junctions("data/vivien/chr*/*.tlx") %>%
    dplyr::mutate(library=gsub("_.*", "", basename(junction_file))) %>%
    dplyr::inner_join(sample_annotations, by=c("library", "bait_chrom"))
  junctions_bait_df = junctions_df %>%
    dplyr::filter(bait_chrom==junction_chrom)

  # scale_factor_df = junctions_bait_df %>%
  #   dplyr::group_by(bait_chrom, gene_silenced, experimental_condition) %>%
  #   dplyr::mutate(replicates=length(unique(junction_file)))  %>%
  #   dplyr::group_by(bait_chrom, gene_silenced) %>%
  #   dplyr::mutate(replicates_factor=replicates/min(replicates)) %>%
  #   dplyr::group_by(junction_file) %>%
  #   dplyr::summarise(libsize=n(), replicates_factor=replicates_factor[1], .groups="keep") %>%
  #   dplyr::ungroup() %>%
  #   dplyr::mutate(libsize_median=median(libsize)) %>%
  #   dplyr::mutate(score=replicates_factor*(libsize_median/libsize)) %>%
  #   dplyr::select(junction_file, score)
  # junctions_bait_df = junctions_bait_df %>%
  #   dplyr::select(-dplyr::matches("score")) %>%
  #   dplyr::inner_join(scale_factor_df, by="junction_file")

      # sample_annotations %>%
      #   dplyr::filter(gene_effected==gene & treatment=="APH" & description=="parental cells")
  #

  params = call_islands_params(binsize=1e5, binstep=1e4, extend=1e5, llocal=2e6, minqvalue=1e-5, maxgap=5e5, minlen=200)
  results = list()
  conditions = junctions_df %>%
    dplyr::group_by(gene_effected, description) %>%
    dplyr::summarize(n.sample=sum(treatment=="APH"), n_files.sample=length(unique(junction_file[treatment=="APH"])), n.control=sum(treatment=="DMSO"), n_files.control=length(unique(junction_file[treatment=="DMSO"])), proportion=ifelse(n.sample>n.control, n.sample/n.control, -n.control/n.sample)) %>%
    dplyr::mutate(gene_id=gene_effected) %>%
    dplyr::inner_join(genes_df, by="gene_id")
  conditions.f = conditions %>% dplyr::filter(n.sample>0 & n.control>0 & abs(proportion) < 4)
  for(d in 1:nrow(conditions.f)) {
    # control_df = junctions_bait_df %>% dplyr::filter(junction_chrom==gene_df$seqnames & gene_effected==gene & treatment=="APH" & description=="parental cells")
    # sample_df = junctions_bait_df %>% dplyr::filter(junction_chrom==gene_df$seqnames & gene_effected==gene & treatment=="APH" & description==d)
    control_df = junctions_bait_df %>% dplyr::inner_join(conditions.f[d,,drop=F], by=c("junction_chrom"="seqnames", "description", "gene_effected")) %>% dplyr::filter(treatment=="DMSO")
    sample_df = junctions_bait_df %>% dplyr::inner_join(conditions.f[d,,drop=F], by=c("junction_chrom"="seqnames", "description", "gene_effected")) %>% dplyr::filter(treatment=="APH")
    if(!nrow(control_df) || !nrow(sample_df)) {
      next
    }

    genome_info.f = genome_info[conditions.f$seqnames[d]]
    control_ranges = GenomicRanges::makeGRangesFromDataFrame(control_df %>% dplyr::mutate(end=junction_end, start=junction_start, seqnames=junction_chrom), seqinfo=genome_info.f, keep.extra.columns=T)
    sample_ranges = GenomicRanges::makeGRangesFromDataFrame(sample_df %>% dplyr::mutate(end=junction_end, start=junction_start, seqnames=junction_chrom), seqinfo=genome_info.f, keep.extra.columns=T)

    results[[conditions.f$gene_id[d]]][[conditions.f$description[d]]] = call_islands(
      sample_ranges=sample_ranges,
      control_ranges=control_ranges,
      params=params,
      debug=F)
  }

    # sample_islands_df2offtargets = islands2offtargets_identity(results$islands$sample, targets_df, genome_mm9)
    pdf(file="reports/viven_experiment_diff.pdf", width=15, height=8)
    for(gene in c("Ctnna2", "Csmd1")) {
      gene_df = genes_df %>%
        dplyr::filter(gene_id==gene) %>%
        dplyr::mutate(island_chrom=seqnames, island_start=start, island_end=end, experimental_condition="Gene")

      coverage_df = do.call(rbind, lapply(names(results[[gene]]), function(d) {
        results[[gene]][[d]]$coverage %>%
          dplyr::mutate(experimental_condition=d) %>%
          dplyr::group_by(experimental_condition, seqnames) %>%
          dplyr::mutate(coverage.control=ifelse(coverage.control==0, NA_real_, coverage.control)) %>%
          dplyr::mutate(coverage.control=zoo::na.fill(coverage.control, "extend")) %>%
          dplyr::mutate(coverage=coverage.sample_scaled-coverage.control)
      }))
      coverage_df = coverage_df%>%
        # dplyr::mutate(coverage=ifelse(is.nan(coverage) | coverage < 1, 0, log10(coverage))) %>%
        dplyr::group_by(experimental_condition) %>%
        dplyr::mutate(coverage.median=mean(coverage[dplyr::between(coverage, quantile(coverage,0.3), quantile(coverage,0.7))])) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(coverage_adj=max(coverage.median)/coverage.median) %>%
        dplyr::mutate(experimental_condition=as.character(experimental_condition), experimental_condition=factor(experimental_condition, c("allelic deletion", "parental cells", "enhancer deleted", "promoter deleted"))) %>%
        dplyr::filter(seqnames == gene_df$seqnames & dplyr::between(start, gene_df$start-3e7, gene_df$end+3e7))
      coverage_df.f = coverage_df %>% dplyr::filter(seqnames == gene_df$seqnames & dplyr::between(start, gene_df$start-1e7, gene_df$end+1e7))

    manual_colors=c("allelic deletion"="#386CB0", "parental cells"="#F0027F", "promoter deleted"="#BBBBBB", "enhancer deleted"="#888888")
    p = ggplot() +
      scale_fill_manual(values=manual_colors) +
      scale_color_manual(values=manual_colors) +
      theme_classic(base_size=16) +
      scale_x_continuous(breaks=scale_breaks) +
      guides(fill=guide_legend(title="Treatment"), color=guide_legend(title="")) +
      labs(x="", y="Coverage (adjusted)")
    p1 = p +
      geom_area(aes(x=start, y=coverage*coverage_adj, fill=experimental_condition), data=coverage_df %>% dplyr::filter(!grepl("parental|allelic", experimental_condition)), position=ggplot2::position_identity()) +
      geom_line(aes(x=start, y=coverage*coverage_adj, color=experimental_condition), data=coverage_df %>% dplyr::filter(grepl("parental|allelic", experimental_condition)), size=0.3) +
      theme(legend.position="none")
    p2 = p +
      geom_area(aes(x=start, y=coverage*coverage_adj, fill=experimental_condition), data=coverage_df.f %>% dplyr::filter(!grepl("parental|allelic", experimental_condition)), position=ggplot2::position_identity()) +
      geom_line(aes(x=start, y=coverage*coverage_adj, color=experimental_condition), data=coverage_df.f %>% dplyr::filter(grepl("parental|allelic", experimental_condition)), size=0.3) +
      coord_cartesian(ylim=c(-10, ifelse(gene=="Ctnnna2", 100, 300))) +
      geom_rect(aes(xmin=start, xmax=end), ymin=-10, ymax=-0, fill="#FF7F00", data=gene_df) +
      geom_text(aes(x=start, label=gene_id), y=-2.5, hjust=0, data=gene_df) +
      theme(legend.direction="vertical", legend.box="horizontal", legend.position="bottom")

    gridExtra::grid.arrange(p1, p2, heights=c(3, 7))


      # baseline_df = dplyr::bind_rows(
      #   results[[gene]][["promoter deleted"]]$coverage %>% dplyr::mutate(coverage_smooth=coverage_smooth.sample_vs_max, experimental_condition="Promoter deleted", coverage_th01=coverage_th01.sample_vs_max),
      #   # results[[gene]]$coverage %>% dplyr::mutate(coverage_smooth=coverage_smooth.control_vs_control, experimental_condition="Control", coverage_th01=coverage_th01.control_vs_control)
      # ) %>% dplyr::filter(seqnames %in% chr)
      # islands_df = dplyr::bind_rows(
      #   results[[gene]][["promoter deleted"]]$islands$sample %>% dplyr::mutate(experimental_condition="Promoter deleted"),
      #   results[[gene]][["enhancer deleted"]]$islands$sample %>% dplyr::mutate(experimental_condition="Enhancer deleted"),
      #   results[[gene]][["allelic deletion"]]$islands$sample %>% dplyr::mutate(experimental_condition="Allelic deletion"),
      #   results[[gene]][["allelic deletion"]]$islands$control %>% dplyr::mutate(experimental_condition="Founder cells"),
      #   gene_df %>% dplyr::mutate(island_start=start, island_end=island_end, island_chrom=seqnames)) %>%
      #   dplyr::mutate(seqnames=island_chrom) %>%
      #   dplyr::filter(seqnames == gene_df$seqnames & dplyr::between(island_start, gene_df$start-1e7, gene_df$end+1e7))
      # karyoplot.colors=c("Allelic deletion"="#000000", "Promoter deleted"="#0000FF", "Enhancer deleted"="#FF0000", Gene="#00FF00", "Founder cells"="#AA3344")
      # p = ggplot_karyoplot(
      #   coverage_df=coverage_df,
      #   # baseline_df=baseline_df,
      #   islands_df=islands_df,
      #   colors=karyoplot.colors,
      #   max_coverage = 1000
      # )
      # print(p)

    }
    dev.off()

  head(p$sample)
}