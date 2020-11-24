library(Rsamtools)
library(GenomicRanges)
library(IRanges)
library(dplyr)
library(readr)
library(ggpmisc)
library(SamSeq)
library(biomaRt)
library(stringi)
library(ggplot2)


analyze.primers_in_mm9 = function() {
  validated_offtargets = readr::read_tsv("data/offtargets.tsv")

  #
  # Find primers position in genome
  #
  # TODO: Find multiple alignments
  primers = readr::read_tsv("data/breaksites_primers.tsv")  %>%
    dplyr::mutate(sgRNA_and_PAM=gsub(" ", "", sgRNA_and_PAM)) %>%
    dplyr::mutate(sgRNA=substr(sgRNA_and_PAM, 1, nchar(sgRNA_and_PAM)-3)) %>%
    dplyr::mutate(Bio=gsub("(5('|-)Bio-| )", "", Bio)) %>%
    reshape2::melt(var="bait_name", measure=c("sgRNA_and_PAM", "sgRNA", "Bio", "Red"), variable.name="primer_name", value.name="primer_sequence") %>%
    dplyr::mutate(primer_sequence_name=paste0(bait_name, "_plus_", primer_name)) %>%
    dplyr::mutate(primer_len=nchar(primer_sequence)) %>%
    dplyr::filter(primer_name=="sgRNA")
  writeLines(with(primers, paste0(">", primer_sequence_name, "\n", primer_sequence)), "tmp/primers.fa")


  #system("bwa mem -t 30 -k 4 -C -a -T 1 -B 1 -O 100 -E 100 -L 100 -D 0.3 -c 10000 -y 8000 data/mm9/mm9.fa.gz tmp/primers.fa > tmp/primers.sam")
  #system("bwa mem -t 30 -k 4 -C -a -T 1 -B 1 -O 100 -E 100 -L 100 -D 0.1 -c 12000 -y 9000 data/mm9/mm9.fa.gz tmp/primers.fa > tmp/primers.sam")
  #system("bwa mem -t 30 -k 4 -C -a -T 1 -B 1 -O 100 -E 100 -L 100 -D 0.05 -c 17000 -y 14000 data/mm9/mm9.fa.gz tmp/primers.fa > tmp/primers.sam")
  if(!file.exists("data/mm9/mm9.fa.gz.bwt")) {
    system("bwa index -b 100000000 data/mm9/mm9.fa.gz")
  }
  system("bwa mem -t 30 -k 5 -C -a -T 10 -B 1 -O 100 -E 100 -L 100 -D 0 -c 70000 -y 60000 data/mm9/mm9.fa.gz tmp/primers.fa > tmp/primers.sam")
  system("samtools sort -@ 30 tmp/primers.sam > tmp/primers_bwa.bam")

  if(!file.exists("data/mm9/bowtie1/mm9.1.ebwt")) {
    dir.create("data/mm9/bowtie1")
    system("bowtie-build --threads 30 data/mm9/mm9.fa.gz data/mm9/bowtie1/mm9")
  }
  system("bowtie --threads 30 --tryhard --all -v 3 --seedlen 5 -f --sam --no-unal data/mm9/bowtie1/mm9 tmp/primers.fa > tmp/primers.sam")
  system("samtools sort -@ 30 tmp/primers.sam > tmp/primers_bowtie.bam")

  # Analyze genome position and give final output
  primers_alignments_param = ScanBamParam(tag=c("AS", "NM"), what=c("qname", "rname", "strand", "flag", "pos", "qwidth",  "cigar", "mapq", "qual"))
  primers_alignments = data.frame()
  for(bam in c("tmp/primers_bwa.bam", "tmp/primers_bowtie.bam")) {
    primers_alignments.d = lapply(Rsamtools::scanBam("tmp/primers_bwa.bam", param=primers_alignments_param), function(z) {
      d = data.frame(z[names(z)!="tag"])
      d$primer_alignment_mismatches=z$tag$NM
      cbind(d, do.call(rbind, lapply(d$flag, SamSeq::samFlags))) })[[1]] %>%
      dplyr::mutate(primer_sequence_name=qname, primer_alignment_end=pos+qwidth) %>%
      dplyr::mutate(primer_alignment_primary=!NOT_PRIMARY_ALIGNMENT) %>%
      dplyr::filter(!READ_UNMAPPED) %>%
      dplyr::select(primer_sequence_name, primer_alignment_chrom=rname, primer_alignment_strand=strand, primer_alignment_start=pos, primer_alignment_end, primer_alignment_flag=flag, primer_alignment_cigar=cigar, primer_alignment_len=qwidth, primer_alignment_mismatches, primer_alignment_primary)
    primers_alignments = rbind(primers_alignments, primers_alignments.d)
  }
  primers_alignments = primers_alignments %>%
    dplyr::distinct(primer_alignment_strand, primer_alignment_chrom, primer_alignment_start, primer_alignment_end, .keep_all=T) %>%
    dplyr::mutate(primer_alignment_id=paste0(primer_sequence_name, "_", 1:n()))

  # Extract sequences that were aligned to primers
  primers_alignments %>%
    dplyr::mutate(
      score=0,
      primer_alignment_start=ifelse(primer_alignment_strand=="+", primer_alignment_start-1, primer_alignment_start-3-1),
      primer_alignment_end=ifelse(primer_alignment_strand=="+", primer_alignment_end+3-1, primer_alignment_end-1)) %>%
    dplyr::select(primer_alignment_chrom, primer_alignment_start, primer_alignment_end, primer_alignment_id, score, primer_alignment_strand) %>%
    readr::write_tsv(file=paste0("tmp/primers_alignments_pos.bed"), col_names=F)
  system("bedtools getfasta -fi data/mm9/mm9.fa -bed tmp/primers_alignments_pos.bed -bedOut -s > tmp/primers_alignments_seq.bed")

  primers_alignments_sequences = readr::read_tsv("tmp/primers_alignments_seq.bed", col_names=c("primer_alignment_chrom", "primer_alignment_start", "primer_alignment_end", "primer_alignment_id", "primer_alignment_score", "primer_alignment_strand", "primer_alignment_and_pam_sequence")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      primer_alignment_sequence=substr(primer_alignment_and_pam_sequence, 1, nchar(primer_alignment_and_pam_sequence)-3),
      primer_alignment_pam_sequence=substr(primer_alignment_and_pam_sequence, nchar(primer_alignment_and_pam_sequence)-3+1, nchar(primer_alignment_and_pam_sequence)),
      primer_alignment_N_total=nchar(gsub("[^N]", "", primer_alignment_sequence)),
      primer_alignment_N_max=max(nchar(stringr::str_split(primer_alignment_sequence, "[^N]")[[1]])),
      primer_alignment_has_pam=grepl("^[A-Z]GG$", toupper(primer_alignment_pam_sequence))
    ) %>%
    dplyr::filter(primer_alignment_N_max<3) %>%
    dplyr::select(primer_alignment_id, primer_alignment_and_pam_sequence, primer_alignment_sequence, primer_alignment_pam_sequence, primer_alignment_N_total, primer_alignment_N_max, primer_alignment_has_pam)

  validated_offtargets_ranges = with(validated_offtargets, IRanges::IRanges(start=offtarget_start, end=offtarget_end))
  primers_targets_ranges = with(primers_alignments,  IRanges::IRanges(start=primer_alignment_start, end=primer_alignment_end))
  validated_offtargets_overlapping = validated_offtargets %>%
    dplyr::inner_join(as.data.frame(IRanges::mergeByOverlaps(primers_targets_ranges, validated_offtargets_ranges)), by=c("offtarget_start"="validated_offtargets_ranges.start", "offtarget_end"="validated_offtargets_ranges.end")) %>%
    dplyr::mutate(offtarget_validated=T) %>%
    dplyr::select(bait_name, bait_chrom, offtarget_is_primary, offtarget_chrom, offtarget_strand, offtarget_start, offtarget_end, primer_alignment_mismatches, offtarget_sequence, offtarget_start, offtarget_end, primer_alignment_start=primers_targets_ranges.start, primer_alignment_end=primers_targets_ranges.end, offtarget_validated)

  primers_targets = primers %>%
    dplyr::inner_join(primers_alignments, by="primer_sequence_name") %>%
    dplyr::inner_join(primers_alignments_sequences, by="primer_alignment_id") %>%
    dplyr::filter(primer_alignment_len==primer_len) %>%
    dplyr::left_join(validated_offtargets_overlapping, by=c("bait_name", "bait_chrom", "primer_alignment_strand"="offtarget_strand", "primer_alignment_chrom"="offtarget_chrom", "primer_alignment_start", "primer_alignment_end")) %>%
    dplyr::mutate(offtarget_validated=tidyr::replace_na(offtarget_validated, F)) %>%
    dplyr::select(-primer_sequence_name, -primer_alignment_id) %>%
    dplyr::mutate(primer_alignment_mismatches_rate=primer_alignment_mismatches/primer_len)

  x = as.matrix(strsplit(primers_targets$primer_sequence,""))
  y = as.matrix(strsplit(primers_targets$primer_alignment_sequence,""))
  primers_targets$primer_alignment_mismatches_check = sapply(1:nrow(primers_targets), function(i) sum(x[[i]]!=y[[i]]))
  primers_targets$primer_alignment_mismatches = primers_targets$primer_alignment_mismatches_check
  #primers_targets = primers_targets %>% dplyr::filter(primer_alignment_mismatches_check==primer_alignment_mismatches)
  primers_targets = primers_targets %>%
    dplyr::filter(primer_alignment_mismatches_check<=10 & primer_alignment_has_pam & !grepl("random|chrY", primer_alignment_chrom)) %>%
    dplyr::filter(primer_name=="sgRNA" & primer_alignment_has_pam & !grepl("random|chrY", primer_alignment_chrom))
  primers_targets %>%
    dplyr::mutate(offtarget_validated=ifelse(offtarget_validated, 1, 0)) %>%
    dplyr::select(bait_name, bait_chrom, offtarget_chrom=primer_alignment_chrom, offtarget_strand=primer_alignment_strand, offtarget_start=primer_alignment_start, offtarget_end=primer_alignment_end, offtarget_mismatches=primer_alignment_mismatches, offtarget_primer_sequence=primer_sequence, offtarget_sequence=primer_alignment_sequence) %>%
    readr::write_tsv(file="data/offtargets_predicted.tsv")

  pdf("reports/primers_offtargets.pdf", width=10, height=6)
  ggplot(primers_targets) +
    geom_histogram(aes(primer_alignment_mismatches, color=primer_name), bins=10) +
    labs(x="Error rate") +
    theme_bw()

  primers_targets.ggplot = primers_targets %>%
    dplyr::mutate(primer_alignment_chrom=factor(primer_alignment_chrom)) %>%
    dplyr::filter(primer_alignment_mismatches<=5)

  primers_targets.ggplot_sum = primers_targets.ggplot %>%
    dplyr::group_by(bait_name, primer_alignment_mismatches) %>%
    dplyr::summarise(mismatches=n()) %>%
    reshape2::dcast(bait_name ~ primer_alignment_mismatches, value.var="mismatches") %>%
    replace(is.na(.), 0)

  ggplot(primers_targets.ggplot) +
    geom_segment(aes(x=primer_alignment_start/1e6, xend=primer_alignment_start/1e6, y=as.numeric(primer_alignment_chrom)-0.2*(1-primer_alignment_mismatches_rate), yend=as.numeric(primer_alignment_chrom)+0.2*(1-primer_alignment_mismatches_rate), color=bait_chrom, size=(1-primer_alignment_mismatches_rate)*0.5)) +
    scale_size_identity() +
    labs(color="Bait chromosome", y="Chromosome", x="Chromosome position") +
    scale_y_continuous(labels=levels(primers_targets.ggplot$primer_alignment_chrom), breaks=1:nlevels(primers_targets.ggplot$primer_alignment_chrom)) +
    theme_bw() +
    ggpmisc::geom_table(aes(x=x, y=y, label=tb), data=tibble(y=2, x=200, tb=list(primers_targets.ggplot_sum)), size=2)
  dev.off()

}