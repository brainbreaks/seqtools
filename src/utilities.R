library(VennDiagram)
library(dplyr)
library(RColorBrewer)
library(IRanges)
library(GenomicRanges)

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