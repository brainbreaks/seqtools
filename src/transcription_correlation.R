library(ggplot2)
library(readr)
library(GenomicRanges)

chromosomes_map_df = readr::read_tsv("data/mm9_chromosomes_synonyms.tsv")


rdc_cols = readr::cols(rdc_cluster=col_character(), rdc_chrom=col_character(), rdc_start=col_double(), rdc_end=col_double(), rdc_group=col_double(), rdc_gene=col_character())
rdc_df = readr::read_tsv("data/rdc_pnas.tsv", col_types=rdc_cols) %>%
  dplyr::inner_join(chromosomes_map_df, by=c("rdc_chrom"="chrom_synonym")) %>%
  dplyr::mutate(rdc_chrom=unique_chrom) %>%
  dplyr::select(-unique_chrom) %>%
  dplyr::mutate(rdc_length=rdc_end-rdc_start) %>%
  dplyr::filter(rdc_length>2e5 & grepl("[a-zA-Z]", rdc_gene))

transcription_df = readr::read_tsv("data/tena2020/tena2020_mm9_transcription.tsv")
transcription_avg_df = transcription_df %>%
  dplyr::group_by(tena2020transcription_gene, tena2020transcription_celltype, tena2020transcription_chrom, tena2020transcription_start, tena2020transcription_end) %>%
  dplyr::summarize(tena2020transcription_rpkm=mean(tena2020transcription_rpkm))

transcription_df %>%
  dplyr::mutate(tena2020transcription_replicate=paste0("Replicate", tena2020transcription_replicate)) %>%
  reshape2::dcast(tena2020transcription_gene+tena2020transcription_celltype ~ tena2020transcription_replicate, value.var="tena2020transcription_rpkm") %>%
  ggplot() +
    geom_point(aes(x=Replicate1, y=Replicate2)) +
    facet_wrap(~tena2020transcription_celltype, scales="free")

rdc_ranges = GenomicRanges::makeGRangesFromDataFrame(rdc_df %>% dplyr::mutate(seqnames=rdc_chrom, start=rdc_start, end=rdc_end), keep.extra.columns=T)
transcription_avg_ranges = GenomicRanges::makeGRangesFromDataFrame(transcription_avg_df %>% dplyr::mutate(seqnames=tena2020transcription_chrom, start=tena2020transcription_start, end=tena2020transcription_end), keep.extra.columns=T)
transcription2rdc_df = as.data.frame(IRanges::mergeByOverlaps(transcription_avg_ranges, rdc_ranges))

transcription2rdc_df.f = transcription2rdc_df %>%
  dplyr::group_by(tena2020transcription_gene) %>%
  dplyr::filter(all(tena2020transcription_rpkm[tena2020transcription_celltype=="NPC"]>0.25 & tena2020transcription_rpkm[tena2020transcription_celltype=="ESC"]<=0.25))

ggplot(transcription2rdc_df.f) +
  geom_boxplot(aes(x=factor(rdc_group), y=tena2020transcription_rpkm, fill=tena2020transcription_celltype)) +
  coord_cartesian(ylim=c(0,3))

transcription2rdc_df.f %>% dplyr::distinct(tena2020transcription_gene, rdc_ranges.rdc_gene)
