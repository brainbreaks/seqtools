# Title     : TODO
# Objective : TODO
# Created by: sandr
# Created on: 4/14/2021

library(dplyr)
library(readr)

repliseq_df = readr::read_tsv("data/zhao_bmc_repliseq_2020/preprocessed/repliseq.tsv") %>%
  dplyr::filter(repliseq_celltype=="npc")

repliseqTime_df = data.frame()
for(i in c(4:16)) {
  if(i == 16) {
    repliseqApprox_df = repliseq_df
  } else {
    i.seq = seq(0, 17, length.out=i)
    i.seq = i.seq[c(-1, -length(i.seq))]

    repliseqApprox_df = repliseq_df %>%
      # dplyr::filter(repliseq_chrom=="chr1" & dplyr::between(repliseq_start, 3.6e7, 3.6e7)) %>%
      dplyr::arrange(repliseq_celltype, repliseq_chrom, repliseq_start, repliseq_end) %>%
      dplyr::group_by(repliseq_celltype, repliseq_chrom, repliseq_start, repliseq_end) %>%
      dplyr::do((function(z){
        zz<<-z
        if(any(is.na(z$repliseq_value))) {
          z.approx = approx(rep(1, 16), xout=i.seq)
          z.result = dplyr::bind_cols(z %>% dplyr::select(-repliseq_fraction, -repliseq_value) %>% dplyr::slice(1), data.frame(repliseq_fraction=z.approx$x, repliseq_value=NA_real_))
        } else {
          z.approx = approx(z$repliseq_value, xout=i.seq)
          z.result = dplyr::bind_cols(z %>% dplyr::select(-repliseq_fraction, -repliseq_value) %>% dplyr::slice(1), data.frame(repliseq_fraction=z.approx$x, repliseq_value=z.approx$y))
        }
        #
        # plot(y=z$repliseq_value, x=z$repliseq_fraction, type="l")
        # lines(y=z.result$repliseq_value, x=z.result$repliseq_fraction, col="#ff0000")

        z.result
      })(.)) %>%
      dplyr::ungroup()
  }

  repliseqTime_df = rbind(repliseqTime_df, repliseq_summarize(repliseqApprox_df, window=3) %>% dplyr::mutate(repliseqTime_nfractions=as.character(i)))
}
repliseqTime_df = repliseqTime_df %>%
  dplyr::mutate(repliseqTime_nfractions_int=dplyr::case_when(
    dplyr::between(as.numeric(repliseqTime_nfractions), 4, 5) ~ "4-5",
    dplyr::between(as.numeric(repliseqTime_nfractions), 6, 7) ~ "6-7",
    dplyr::between(as.numeric(repliseqTime_nfractions), 8, 9) ~ "8-9",
    dplyr::between(as.numeric(repliseqTime_nfractions), 10, 11) ~ "10-11",
    dplyr::between(as.numeric(repliseqTime_nfractions), 12, 13) ~ "12-13",
    dplyr::between(as.numeric(repliseqTime_nfractions), 14, 15) ~ "14-15",
    dplyr::between(as.numeric(repliseqTime_nfractions), 16, 16) ~ "16"
  ), repliseqTime_nfractions_offset=as.numeric(repliseqTime_nfractions)-as.numeric(gsub("-.*", "", repliseqTime_nfractions_int)))

chromsizes_cols = readr::cols(seqnames=col_character(), seqlengths=col_double())
genome_info = with(readr::read_tsv("data/mm9/mm9.chrom.sizes", col_names=names(chromsizes_cols$cols), col_types=chromsizes_cols),
   GenomeInfoDb::Seqinfo(seqnames, seqlengths, isCircular=rep(F, length(seqnames)), genome=rep("mm10", length(seqnames))))
genome_info = genome_info[paste0("chr", c(1:19, "X", "Y"))]
repliseqTime_ranges = GenomicRanges::makeGRangesFromDataFrame(repliseqTime_df %>% dplyr::mutate(seqnames=repliseqTime_chrom, start=repliseqTime_start, end=repliseqTime_start), keep.extra.columns=T, seqinfo=genome_info)
tiles_df = as.data.frame(unlist(GenomicRanges::tileGenome(GenomeInfoDb::seqinfo(repliseqTime_ranges), tilewidth=1e6))) %>%
  dplyr::mutate(tile_start=start, tile_end=end, tile_chrom=seqnames)
tiles_ranges = GenomicRanges::makeGRangesFromDataFrame(tiles_df %>% dplyr::mutate(seqnames=tile_chrom, start=tile_start, tile_end=tile_end), keep.extra.columns=T, seqinfo=genome_info)

repliseqTime2tiles_df = as.data.frame(IRanges::mergeByOverlaps(repliseqTime_ranges, tiles_ranges)) %>% dplyr::select(-dplyr::matches("_ranges\\."))
repliseqTime2slopes_df = repliseqTime2tiles_df %>%
  dplyr::group_by(repliseqTime_nfractions, repliseqTime_nfractions_int, repliseqTime_nfractions_offset, tile_chrom, tile_start, tile_end) %>%
  dplyr::summarize(tile_slope=max(diff(abs(repliseqTime_avg)))) %>%
  dplyr::filter(is.finite(tile_slope)) %>%
  dplyr::mutate(repliseqTime_nfractions=as.numeric(repliseqTime_nfractions))

repliseqTime2slopes_R = repliseqTime2slopes_df %>%
  dplyr::inner_join(repliseqTime2slopes_df, by=c("tile_chrom", "tile_start", "tile_end")) %>%
  dplyr::group_by(repliseqTime_nfractions.x, repliseqTime_nfractions.y) %>%
  dplyr::summarize(R=cor(tile_slope.x, tile_slope.y))
repliseqTime2slopes_R.short = repliseqTime2slopes_R %>%
  dplyr::filter(repliseqTime_nfractions.x==16)

pdf(file="reports/evaluate_repliseq_fraction_number.pdf", width=11, height=6)
ggplot(repliseqTime2slopes_R.short) +
  geom_line(aes(x=repliseqTime_nfractions.y, y=R)) +
  geom_point(aes(x=repliseqTime_nfractions.y, y=R)) +
  scale_x_continuous(breaks=seq(0, 16, by=1))

repliseqTime_ggplot = repliseqTime_df %>%
  dplyr::filter(repliseqTime_chrom=="chr1" & dplyr::between(repliseqTime_start, 3.6e7, 4.605e7)) %>%
  dplyr::mutate(fractions=ifelse(repliseqTime_nfractions == "16", "16", "less"))
repliseqTime_ggplot = dplyr::bind_rows(repliseqTime_ggplot, repliseqTime_ggplot %>% dplyr::filter(repliseqTime_nfractions==16) %>% dplyr::mutate(repliseqTime_nfractions_int="4-5"))
repliseqTime_ggplot = dplyr::bind_rows(repliseqTime_ggplot, repliseqTime_ggplot %>% dplyr::filter(repliseqTime_nfractions==16) %>% dplyr::mutate(repliseqTime_nfractions_int="6-7"))
repliseqTime_ggplot = dplyr::bind_rows(repliseqTime_ggplot, repliseqTime_ggplot %>% dplyr::filter(repliseqTime_nfractions==16) %>% dplyr::mutate(repliseqTime_nfractions_int="8-9"))
repliseqTime_ggplot = dplyr::bind_rows(repliseqTime_ggplot, repliseqTime_ggplot %>% dplyr::filter(repliseqTime_nfractions==16) %>% dplyr::mutate(repliseqTime_nfractions_int="10-11"))
repliseqTime_ggplot = repliseqTime_ggplot %>%
  dplyr::filter(repliseqTime_nfractions_int %in% c("4-5", "6-7", "8-9", "10-11")) %>%
  dplyr::mutate(repliseqTime_nfractions_offset=ifelse(repliseqTime_nfractions=="16", "16", paste0("+", repliseqTime_nfractions_offset)))
ggplot(repliseqTime_ggplot) +
  geom_line(aes(x=repliseqTime_start, y=repliseqTime_avg, group=repliseqTime_nfractions, color=repliseqTime_nfractions_offset, size=fractions, alpha=fractions)) +
  facet_wrap(~repliseqTime_chrom) +
  labs(y="Fraction", x="chr1: 3.6M-46M", color="Fractions (check plot name)") +
  scale_size_manual(values=c("16"=3, "less"=0.5)) +
  scale_alpha_manual(values=c("16"=0.5, "less"=1)) +
  scale_color_manual(values=c("16"="#666666", "+0"="#FF0000", "+1"="#0000FF")) +
  facet_wrap(~repliseqTime_nfractions_int) +
  guides(alpha=F, size=F) +
  theme_bw()
dev.off()
