# Title     : TODO
# Objective : TODO
# Created by: sandr
# Created on: 3/24/2021
library(Biostrings)
library(GenomicRanges)

genome = Biostrings::readDNAStringSet("data/mm10/mm10.fa.gz", "fasta")[paste0("chr", c(1:19, "X", "Y"))]

tile_width = 1e3
tiles <- unlist(tileGenome(seqinfo(genome), tilewidth=1e3))
dna = Biostrings::getSeq(genome, tiles)
freq = alphabetFrequency(dna)
gc = rowSums(freq[,c("G","C")])/ rowSums(freq[,c("G","C", "A", "T")])
result = as.data.frame(tiles) %>%
  dplyr::mutate(gc_freq=gc) %>%
  dplyr::rename(gc_chrom="seqnames", gc_start="start",  gc_end="end") %>%
  dplyr::select(dplyr::matches("gc_"))

readr::write_tsv(result, file="data/mm10/mm10_tile1000.freq")
