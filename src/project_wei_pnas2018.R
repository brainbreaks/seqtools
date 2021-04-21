library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(doParallel)
source("src/visualization.R")
source("src/islands.R")

main = function() {
  read_mm9_job = parallel::mcparallel(Biostrings::readDNAStringSet("data/mm9/mm9.fa.gz", "fasta"))

  chromsizes_cols = readr::cols(seqnames=col_character(), seqlengths=col_double())
  genome_info = with(readr::read_tsv("data/mm9/mm9.chrom.sizes", col_names=names(chromsizes_cols$cols), col_types=chromsizes_cols),
             GenomeInfoDb::Seqinfo(seqnames, seqlengths, isCircular=rep(F, length(seqnames)), genome=rep("mm9", length(seqnames))))
  genome_info = genome_info[paste0("chr", c(1:19, "X", "Y"))]

  chromosomes_map_df = readr::read_tsv("data/mm9_chromosomes_synonyms.tsv")

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
    dplyr::mutate(junction_name=stringr::str_glue("{bait}_{name}", bait=bait_chrom, name=junction_name))
  scale_factor_df = junctions_df %>%
    dplyr::group_by(experimental_condition, bait_chrom) %>%
    dplyr::summarise(libsize=n()) %>%
    dplyr::mutate(libsize_median=median(libsize)) %>%
    dplyr::mutate(scale_factor=libsize_median/libsize) %>%
    dplyr::select(experimental_condition, bait_chrom, scale_factor)
  junctions_df = junctions_df %>%
    dplyr::select(-dplyr::matches("scale_factor")) %>%
    dplyr::inner_join(scale_factor_df, by=c("experimental_condition", "bait_chrom")) %>%
    dplyr::filter(junction_chrom %in% seqlevels(genome_info)) %>%
    dplyr::filter(junction_chrom==bait_chrom)

  sample_df = junctions_df %>% dplyr::filter(experimental_condition=="Sample")
  sample_ranges = GenomicRanges::makeGRangesFromDataFrame(sample_df %>% dplyr::mutate(end=junction_end, start=junction_start, seqnames=junction_chrom), seqinfo=genome_info, keep.extra.columns=T)
  control_df = junctions_df %>% dplyr::filter(experimental_condition=="Control")
  control_ranges = GenomicRanges::makeGRangesFromDataFrame(control_df %>% dplyr::mutate(end=junction_end, start=junction_start, seqnames=junction_chrom), seqinfo=genome_info,  keep.extra.columns=T)
  params = call_islands_params(binsize=1e5, binstep=1e4, extend=1e5, llocal=2e6, minqvalue=0.01, maxgap=5e5, minlen=200)
  results = call_islands(sample_ranges, control_ranges, params, debug=F)



  if(debug) {
    rdc_cols = readr::cols(rdc_cluster=col_character(), rdc_chrom=col_character(), rdc_start=col_double(), rdc_end=col_double(), rdc_group=col_double(), rdc_gene=col_character())
    rdc_df = readr::read_tsv("data/rdc_pnas.tsv", col_types=rdc_cols) %>%
      dplyr::inner_join(chromosomes_map_df, by=c("rdc_chrom"="chrom_synonym")) %>%
      dplyr::mutate(rdc_chrom=unique_chrom) %>%
      dplyr::select(-unique_chrom) %>%
      dplyr::mutate(rdc_length=rdc_end-rdc_start) %>%
      dplyr::mutate(rdc_id=1:n())

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


    genome_mm9 = mccollect(list(read_mm9_job))[[1]]
    sample_islands_df2offtargets = islands2offtargets_identity(results$islands$sample, targets_df, genome_mm9)
    pdf(file="reports/duo_baseline_extend1e5_smooth2e6_bin1e5_step1e4_3.pdf", width=15, height=8)
    for(chr in paste0("chr", 1:19)) {
      print(chr)


        geom_line(aes(x=start, y=coverage*coverage_adj, color=experimental_condition), data=coverage_df %>% dplyr::filter(grepl("parental|allelic", experimental_condition)), size=0.3) +
        scale_fill_manual(values=manual_colors) +
        scale_color_manual(values=manual_colors) +
        theme_classic(base_size=16) +
        scale_x_continuous(breaks=scale_breaks) +
        guides(fill=guide_legend(title="Treatment"), color=guide_legend(title="")) +
        labs(x="", y="Coverage (adjusted)")

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
        # sample_islands_df2offtargets %>% dplyr::mutate(seqnames=island_chrom, island_summit_qvalue=bait_alignment_identity, score=(bait_alignment_identity-min(bait_alignment_identity))/(100-max(bait_alignment_identity))*max(coverage_df$coverage), condition="Offtarget"),
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


  #
  # VENN diagram overlap RDC and sample islands
  #
  jpeg("reports/rdc2islands_venn_2.jpg", width=800, height=800)
  islands_filtered_df = results$islands$sample %>% dplyr::filter(island_summit_qvalue.control < island_summit_qvalue) %>%  dplyr::mutate(condition="Sample filtered")
  rdc_ranges = GenomicRanges::makeGRangesFromDataFrame(rdc_df %>% dplyr::mutate(seqnames=rdc_chrom, start=rdc_start, end=rdc_end), keep.extra.columns=T)
  islands_filtered_ranges = GenomicRanges::makeGRangesFromDataFrame(islands_filtered_df %>% dplyr::mutate(seqnames=island_chrom, start=island_start, end=island_end), keep.extra.columns=T)
  venn_ranges(rdc_ranges, islands_filtered_ranges, name1="RDC", name2="island")
  dev.off()


  #
  # VENN diagram overlap between SAMPLE and CONTROL
  #
  jpeg("reports/control_vs_sample_venn.jpg", width=800, height=800)
  control_islands_ranges = GenomicRanges::makeGRangesFromDataFrame(results$islands$control, end.field="island_end", start.field="island_start", seqnames.field="island_chrom", keep.extra.columns=T)
  sample_islands_ranges = GenomicRanges::makeGRangesFromDataFrame(results$islands$sample, end.field="island_end", start.field="island_start", seqnames.field="island_chrom", keep.extra.columns=T)
  venn_ranges(control_islands_ranges, sample_islands_ranges, name1="Control", name2="Sample")
  dev.off()
}


test_robustness = function() {
  chromsizes_cols = readr::cols(seqnames=col_character(), seqlengths=col_double())
  genome_info = with(readr::read_tsv("data/mm9/mm9.chrom.sizes", col_names=names(chromsizes_cols$cols), col_types=chromsizes_cols),
             GenomeInfoDb::Seqinfo(seqnames, seqlengths, isCircular=rep(F, length(seqnames)), genome=rep("mm9", length(seqnames))))[paste0("chr", c(1:19))]
  chromosomes_map_df = readr::read_tsv("data/mm9_chromosomes_synonyms.tsv")

  junctions_df = read_junctions("data/breaks_tlx/*/*.tlx") %>%
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
  junctions_df.f = junctions_df %>%
    dplyr::select(-dplyr::matches("scale_factor")) %>%
    dplyr::inner_join(scale_factor_df, by=c("experimental_condition", "bait_chrom")) %>%
    dplyr::filter(junction_chrom %in% seqlevels(genome_info))
    #dplyr::filter(junction_chrom=="chr2" & junction_chrom==bait_chrom)
    #dplyr::mutate(bait_chrom="chr2")
  genome_info.f = genome_info[unique(junctions_df.f$junction_chrom)]

  i = 0
  peaks_robustness = data.frame()
  for(prop in rev(seq(0.3, 1, 0.05))) {
    if(prop==1) {
      rep_all = 1
    } else {
      rep_all = 1:10
    }
    for(rep in rep_all) {
      i = i + 1
      print(paste(i, ":", prop, " / ", rep))
      sample_df = junctions_df.f %>% dplyr::filter(experimental_condition=="Sample" & junction_chrom==bait_chrom) %>% dplyr::sample_frac(prop)
      control_df = junctions_df.f %>% dplyr::filter(experimental_condition=="Control" & junction_chrom==bait_chrom) %>% dplyr::sample_frac(prop)

      sample_ranges = GenomicRanges::makeGRangesFromDataFrame(sample_df %>% dplyr::mutate(end=junction_end, start=junction_start, seqnames=junction_chrom), keep.extra.columns=T, seqinfo=genome_info.f)
      control_ranges = GenomicRanges::makeGRangesFromDataFrame(control_df %>% dplyr::mutate(end=junction_end, start=junction_start, seqnames=junction_chrom), keep.extra.columns=T, seqinfo=genome_info.f)
      params = call_islands_params(binsize=1e5, binstep=1e4, extend=1e5, llocal=2e6, minqvalue=0.01, maxgap=5e5, minlen=200)
      results = call_islands(sample_ranges, control_ranges, genome_info.f, params)

      counts = data.frame(robustness_rep=rep, robustness_prop=prop, robustness_count.sample=nrow(sample_df), robustness_count.control=nrow(control_df))
      peaks_robustness = rbind(peaks_robustness, results$islands$sample %>% tidyr::crossing(counts))
    }
  }

  peaks_robustness.sum = peaks_robustness %>%
    tidyr::crossing(data.frame(robustness_qvalue=2:20)) %>%
    dplyr::group_by(robustness_prop, robustness_rep, robustness_qvalue) %>%
    dplyr::summarise(n=sum(island_summit_qvalue>robustness_qvalue)-1)

 pdf("reports/robustness3.pdf", width=16, height=16)
 ggplot(peaks_robustness.sum) +
   geom_point(aes(x=robustness_prop, y=n), alpha=0.2) +
   geom_smooth(aes(x=robustness_prop, y=n)) +
   facet_wrap(~robustness_qvalue, scales="free_y")
  dev.off()
