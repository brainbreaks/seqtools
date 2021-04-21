library(Rsamtools)
library(GenomicRanges)
library(IRanges)
library(dplyr)
library(readr)
library(ggpmisc)
library(SamSeq)
library(biomaRt)
library(stringi)
library(stringr)
library(ggplot2)
library(Biostrings)
library(BSgenome)

analyze.primers_in_mm10 = function() {
  genome_txdb = GenomicFeatures::makeTxDbFromGFF('data/mm10/mm10.refGene.gtf.gz', format="gtf")
  genes_df = as.data.frame(GenomicFeatures::genes(genome_txdb)) %>% dplyr::mutate(seqnames=as.character(seqnames))
  genes_ranges = GenomicRanges::makeGRangesFromDataFrame(genes_df %>% dplyr::mutate(gene_chrom=seqnames, gene_start=start, gene_end=end), keep.extra.columns=T)


  chromosomes_map_df = readr::read_tsv("data/mm9_chromosomes_synonyms.tsv")
  chromosomes_map = chromosomes_map_df$chrom
  names(chromosomes_map) = chromosomes_map_df$chrom_synonym

  #
  # Find primers position in genome
  #
  primers = readr::read_tsv("data/Lorenzo/annotations.tsv")  %>%
    reshape2::melt(var=c("library", "treatment_gene"), measure=c("treatment_gene_sgRNA", "treatment_bait_sgRNA"), variable.name="primer_name", value.name="primer_sequence") %>%
    dplyr::mutate(primer_sequence_name=paste0(library, "_plus_", primer_name)) %>%
    dplyr::mutate(primer_len=nchar(primer_sequence)) %>%
    dplyr::filter(primer_len>5) %>%
    dplyr::filter(primer_name=="treatment_gene_sgRNA")
  writeLines(with(primers, paste0(">", primer_sequence_name, "\n", primer_sequence)), "tmp/primers.fa")


  # system("tail -q -n +24 data/Lorenzo/callpeak/*_peaks.xls >> data/Lorenzo/macs.bed")
  # system("bedtools slop -i data/Lorenzo/macs.bed -g data/mm10/mm10.chrom.sizes -b 5000 > data/Lorenzo/macs_extended.bed; mv data/Lorenzo/macs_extended.bed data/Lorenzo/macs.bed")
  # islands_cols = cols(
  #   island_chrom=col_character(), island_start=col_double(), island_end=col_double(), island_length=col_double(), island_summit=col_double(), island_pileup=col_double(),
  #   island_pvalue_log10=col_character(), island_fc=col_double(), island_qvalue_log10=col_double(), island_name=col_character()
  # )
  # islands_df = readr::read_tsv("data/Lorenzo/macs.bed", col_names=names(islands_cols$cols), col_types=islands_cols)
  islands_df = do.call(rbind, lapply(Sys.glob("data/Lorenzo/callpeak/*_peaks.xls"), readr::read_tsv, col_names=names(islands_cols$cols), col_types=islands_cols, skip=23)) %>%
    dplyr::mutate(island_start=island_start-5000, island_end=island_end+5000) %>%
    dplyr::mutate(library=gsub("(^L|_.*)", "", island_name))

  islands_df %>%
    dplyr::mutate(island_strand="+") %>%
    dplyr::select(island_chrom, island_start, island_end, island_name, island_qvalue_log10, island_strand) %>%
    readr::write_tsv(file="data/Lorenzo/macs.bed", col_names=F)

  system("bedtools getfasta -nameOnly -fi data/mm10/mm10.fa -bed data/Lorenzo/macs.bed -fo data/Lorenzo/macs.fa")
  genome = Biostrings::readDNAStringSet("data/Lorenzo/macs.fa", "fasta")


  system("bwa index -b 100000000 data/Lorenzo/macs.fa")
  system("bwa mem -t 30 -k 5 -C -a -T 10 -B 1 -O 100 -E 100 -L 100 -D 0 -c 70000 -y 60000 data/Lorenzo/macs.fa tmp/primers.fa > tmp/primers_bwa.sam")
  system("samtools sort -@ 30 tmp/primers_bwa.sam > tmp/primers_bwa.bam")

  system("bowtie-build --threads 30 data/Lorenzo/macs.fa data/Lorenzo/macs")
  system("bowtie --threads 30 --tryhard --all -v 3 --seedlen 5 -f --sam --no-unal data/Lorenzo/macs tmp/primers.fa > tmp/primers_bowtie3.sam")
  system("samtools sort tmp/primers_bowtie3.sam > tmp/primers_bowtie3.bam")
  system("samtools view -S -b tmp/primers_bowtie3.sam > tmp/primers_bowtie3.bam")




  # Analyze genome position and give final output
  primers_alignments_param = ScanBamParam(tag=c("AS", "NM"), what=c("qname", "rname", "strand", "flag", "pos", "qwidth",  "cigar", "mapq", "qual"))
  primers_alignments = data.frame()
  for(bam_file in c("tmp/primers_bwa.bam", "tmp/primers_bowtie3.bam")) {
    primers_alignments.d = lapply(Rsamtools::scanBam("tmp/primers_bwa.bam", param=primers_alignments_param), function(z) {
      zz<<-z
      d = data.frame(z[names(z)!="tag"])
      d$primer_alignment_mismatches=z$tag$NM
      cbind(d, do.call(rbind, lapply(d$flag, SamSeq::samFlags))) })[[1]]
    if(!nrow(primers_alignments.d)) next

    primers_alignments.d = primers_alignments.d %>%
        dplyr::mutate(primer_sequence_name=qname, primer_alignment_end=pos+qwidth-1) %>%
        dplyr::mutate(primer_alignment_primary=!NOT_PRIMARY_ALIGNMENT) %>%
        dplyr::filter(!READ_UNMAPPED & rname!="chrM") %>%
        dplyr::select(primer_sequence_name, primer_alignment_seqname=rname, primer_alignment_strand=strand, primer_alignment_start=pos, primer_alignment_end, primer_alignment_flag=flag, primer_alignment_cigar=cigar, primer_alignment_len=qwidth, primer_alignment_mismatches, primer_alignment_primary)
    primers_alignments = rbind(primers_alignments, primers_alignments.d)
  }

  primers_alignments = primers_alignments %>%
    dplyr::distinct(primer_sequence_name, primer_alignment_seqname, primer_alignment_strand, primer_alignment_start, primer_alignment_end, .keep_all=T) %>%
    dplyr::mutate(
      primer_alignment_and_pam_start=ifelse(primer_alignment_strand=="+", primer_alignment_start, primer_alignment_start-3),
      primer_alignment_and_pam_end=ifelse(primer_alignment_strand=="+", primer_alignment_end+3, primer_alignment_end),
      primer_alignment_and_pam_sequence=as.character(Biostrings::getSeq(genome, GRanges(primer_alignment_seqname, IRanges(start=primer_alignment_and_pam_start, end=primer_alignment_and_pam_end), strand=primer_alignment_strand))),
      primer_alignment_sequence=substr(primer_alignment_and_pam_sequence, 1, nchar(primer_alignment_and_pam_sequence)-3),
      primer_alignment_pam_sequence=substr(primer_alignment_and_pam_sequence, nchar(primer_alignment_and_pam_sequence)-3+1, nchar(primer_alignment_and_pam_sequence)),
      primer_alignment_N_total=nchar(gsub("[^N]", "", primer_alignment_sequence)),
      primer_alignment_N_max=max(nchar(stringr::str_split(primer_alignment_sequence, "[^N]")[[1]])),
      primer_alignment_has_pam=grepl("^[A-Z]GG$", toupper(primer_alignment_pam_sequence))) %>%
      dplyr::inner_join(primers, by="primer_sequence_name") %>%
      dplyr::inner_join(islands_df, by=c("primer_alignment_seqname"="island_name", "library")) %>%
      dplyr::inner_join(genes_df %>% dplyr::select(treatment_gene_chrom=seqnames, gene_id), by=c("treatment_gene"="gene_id"))


  primers_alignments %>%
    dplyr::select(gene=treatment_gene, gene_chrom=treatment_gene_chrom, library=library, alignment_length=primer_alignment_len, alignment_mismatches=primer_alignment_mismatches, alignment_cigar=primer_alignment_cigar, alignment_strand=primer_alignment_strand, island_name=primer_alignment_seqname, island_chrom, island_start, island_end) %>%
    dplyr::arrange(gene, alignment_mismatches, library) %>%
    readr::write_tsv("data/Lorenzo/offtargets.tsv")

  table(gene=paste(primers_alignments$treatment_gene, primers_alignments$treatment_gene_chrom), island=primers_alignments$island_chrom, primers_alignments$primer_alignment_mismatches)

  primers_alignments.f = primers_alignments %>% dplyr::filter(island_chrom)

  primers_alignments_ranges = GenomicRanges::makeGRangesFromDataFrame(primers_alignments %>% dplyr::mutate(seqnames=island_chrom, start=primer_alignment_start, end=primer_alignment_end), keep.extra.columns=T)
  primers_alignments2genes = as.data.frame(IRanges::mergeByOverlaps(primers_alignments_ranges, genes_ranges))

  primers_alignments
  primers_alignments2genes = primers_alignments
  genes_ranges

  x = primers_alignments %>% dplyr::select(library, treatment_bait_chrom, primer_alignment_sequence, primer_alignment_pam_sequence, primer_alignment_mismatches)
  table(x$primer_alignment_mismatches, x$library)
}