library(Rsamtools)
library(GenomicRanges)
library(IRanges)
library(dplyr)
library(readr)
library(ggpmisc)
library(SamSeq)


analyze.primers_in_mm9 = function() {
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

  #system("bwa index mm9.fa")

  system("bwa aln -n 6 -o 6 -k 2 data/mm9.fa data/breaksites_primers.fa > data/breaksites_primers.sai")
  system("bwa samse data/mm9.fa data/breaksites_primers.sai data/breaksites_primers.fa > data/breaksites_primers.sam")
  system("samtools view -S -b data/breaksites_primers.sam > data/breaksites_primers.bam")

  #system("bwa mem -a -T 19 -k 2 data/mm9.fa data/breaksites_primers.fa > data/breaksites_primers.sam")
  #system("samtools view -S -b data/breaksites_primers.sam > data/breaksites_primers.bam")

  #system("bowtie2 --local -x data/mm9 -f data/breaksites_primers.fa -S data/breaksites_primers.sam --no-unal")
  #system("samtools view -S -b data/breaksites_primers.sam > data/breaksites_primers.bam")
  #file.remove(c("data/breaksites_primers.fa", "data/breaksites_primers.sam"))


  #
  # Analyze genome position and give final output
  #
  primers_alignment_bam = Rsamtools::scanBam("data/breaksites_primers.bam")
  primers_alignment_table = lapply(primers_alignment_bam, as.data.frame)[[1]] %>%
    dplyr::rowwise() %>%
    dplyr::do((function(z) cbind(z, as.data.frame(t(SamSeq::samFlags(z$flag)))))(.)) %>%
    dplyr::filter(!READ_UNMAPPED) %>%
    tidyr::extract(qname, into=c("Bait", "Primer"), regex="(.*)_plus_(.*)") %>%
    dplyr::select(Bait, Primer, alignment_chrom=rname, alignment_strand=strand, alignment_pos=pos, alignment_seq=seq, alignment_flag=flag, alignment_cigar=cigar)

  breaksites_compare = breaksites_primers %>%
    dplyr::select(Bait, Chromosome, Primer, start, end, strand, Comment) %>%
    dplyr::left_join(primers_alignment_table, by=c("Bait", "Primer"))
  breaksites_compare.sum = breaksites_compare %>%
    dplyr::mutate(
      Bait=paste0(Bait, strand),
      Comment=tidyr::replace_na(Comment, ""),
      aligned=dplyr::case_when(
        is.na(alignment_pos) ~ "N/A",
        strand != alignment_strand ~ paste0("(", alignment_strand, ")"),
        alignment_pos != start ~ paste0(round((alignment_pos-start)/1e6,1), "Mbp"),
        T ~ "")) %>%
    # dplyr::group_by(Bait) %>%
    # dplyr::filter(any(is.na(alignment_pos))) %>%
    reshape2::dcast(Bait+Comment ~ Primer, value.var="aligned") %>%
    dplyr::select(-Comment, Comment)


  pdf("reports/primers_in_mm9.pdf", width=10, height=6)
  ggplot(breaksites_compare) +
    geom_bar(aes(x=Bait, y=(alignment_pos-start)/1e6, fill=Primer), stat="identity", position="dodge") +
    coord_flip() +
    labs(y="Distance to old position (Mbp)") +
    ggpmisc::geom_table(aes(x=x, y=y, label=tb), data=tibble(y=7.5, x=2, tb = list(breaksites_compare.sum)), size=2)

    ggplot(breaksites_compare) +
      geom_point(aes(x=alignment_pos, y=start, color=Primer, Bait=Bait))
  dev.off()

  p = ggplot(breaksites_compare) +
    geom_point(aes(x=alignment_pos, y=start, color=Primer, Bait=Bait))

  plotly::ggplotly(p)
}