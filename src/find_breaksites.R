library(Rsamtools)
library(dplyr)
library(readr)

#
# Find primers position in genome
#
# TODO: Find multiple alignments
breaksites_primers = readr::read_tsv("data/breaksites_primers.tsv")  %>% 
  reshape2::melt(var="Bait", measure=c("sgRNA_and_PAM", "Bio_primer", "Red_Primer"), variable.name="Primer", value.name="seq") %>%
  dplyr::mutate(
    seq=gsub("(5('|-)Bio-| )", "", seq), 
    seq_name=paste0(Bait, "_plus_", Primer))
breaksites_fasta = with(breaksites_primers, paste0(">", seq_name, "\n", seq))
writeLines(breaksites_fasta, "data/breaksites_primers.fa")
system("bowtie2 --local -x data/mm9 -f data/breaksites_primers.fa -S data/breaksites_primers.sam --no-unal")
system("samtools view -S -b data/breaksites_primers.sam > data/breaksites_primers.bam")
file.remove(c("data/breaksites_primers.fa", "data/breaksites_primers.sam"))



#
# Analyze genome position and give final output
#
primers_alignment_bam = Rsamtools::scanBam("data/breaksites_primers.bam")
primers_alignment_table = lapply(primers_alignment_bam, as.data.frame)[[1]] %>%
  tidyr::extract(qname, into=c("Bait", "Primer"), regex="(.*)_plus_(.*)") %>%
  dplyr::select(Bait, Primer, alignment_chrom=rname, alignment_strand=strand, alignment_pos=pos)

primers_alignment_table %>%
  dplyr::group_by(Bait, Primer) %>%
  dplyr::filter(n()>1)

breaksites_compare = breaksites_primers %>%
  dplyr::select(Bait, Chromosome, Primer, start, end, strand, Comment) %>%
  dplyr::left_join(primers_alignment_table, by=c("Bait", "Primer"))

#
# Not all breaks should be integrated
# 
breaksites_compare %>%
  dplyr::mutate(aligned=ifelse(is.na(alignment_pos), "N", ifelse(strand != alignment_strand, "(!)", "Y")), Comment=tidyr::replace_na(Comment, "")) %>%
  # dplyr::group_by(Bait) %>%
  # dplyr::filter(any(is.na(alignment_pos))) %>%
  reshape2::dcast(Bait+start+end+Comment ~ Primer, value.var="aligned") %>%
  dplyr::select(-Comment, Comment)

breaksites_compare %>% dplyr::filter(!is.na(alignment_pos) & strand != alignment_strand)

ggplot(breaksites_compare) +
  geom_point(aes(x=alignment_pos, y=start))



primers_alignment_table %>%