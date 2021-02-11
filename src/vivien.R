library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(plotly)
library(IRanges)


breaksites = readr::read_tsv("data/offtargets_predicted.tsv") %>%
   dplyr::filter(offtarget_mismatches==0) %>%
   dplyr::mutate(offtarget_start=offtarget_start-5e6, offtarget_end=offtarget_end+5e6) %>%
   dplyr::mutate(offtarget_id=1:n())

# breaksites = readr::read_tsv("breaksites.tsv") %>%
#    dplyr::mutate(breaksite_start=breaksite_pos-1e6, breaksite_end=breaksite_pos+1e6) %>%
#    dplyr::mutate(breaksite_id=1:n())

groseq_cols=cols(groseq_chrom=col_character(), groseq_start=col_double(), groseq_end=col_double(), groseq_strand=col_character(), groseq_transcript=col_character(), groseq_gene=col_character(), groseq_reads=col_double(), groseq_length=col_double(), groseq_rpkm=col_double())
groseq = readr::read_tsv("data/Xrcc4_p53_double.transcript.longest.cnt", col_names=names(groseq_cols$cols), col_types=groseq_cols) %>%
   dplyr::mutate(groseq_id=1:n()) %>%
   dplyr::mutate(groseq_gene=dplyr::case_when(
        groseq_gene=="Large" ~ "Large1",
        groseq_gene=="Naalad2" ~ "Naaladl2",
        T~groseq_gene))



breaks = data.frame()
breaks_cols = cols(break_chrom=col_character(), break_start=col_double(), break_end=col_double(), break_name=col_character(), break_score=col_character(), break_strand=col_character())
for(breaks_path in list.files("data/breaks", pattern="*.bed", full.names=T)) {
    chr = gsub("(chr)([^_]+).*$", "\\L\\1\\2", basename(breaks_path), perl=T, ignore.case=T)
    condition = ifelse(grepl("DMSO", breaks_path), "Control", "Sample")
    sample = gsub("\\.bed", "", basename(breaks_path))


    x = readr::read_tsv(breaks_path, col_names=names(breaks_cols$cols), col_types=breaks_cols) %>%
      dplyr::mutate(break_bait_chrom=chr, break_exp_condition=condition, break_sample=sample, breaks_file=breaks_path)
    breaks = rbind(breaks, x)
}
breaks = breaks %>% dplyr::mutate(break_id=1:n())


#
# Load RDC
#
rdc_all= readr::read_tsv("data/rdc_pnas.tsv")
rdc = rdc_all %>%
  dplyr::mutate(rdc_length=rdc_end-rdc_start) %>%
  tidyr::separate_rows(rdc_gene, sep=" *, *") %>%
  dplyr::filter(rdc_gene != "--") %>%
  dplyr::group_by(rdc_gene) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(rdc_id=1:n())


#
# Remove breaks close to bait
#
breaksites_ranges = with(breaksites, IRanges::IRanges(start=offtarget_start, end=offtarget_end, offtarget_id=offtarget_id, offtarget_chrom=offtarget_chrom))
breaks_ranges = with(breaks,  IRanges::IRanges(start=break_start, end=break_end, break_id=break_id, break_chrom=break_chrom))
breaksites2break = as.data.frame(IRanges::mergeByOverlaps(breaksites_ranges, breaks_ranges)) %>%
  dplyr::filter(break_chrom==offtarget_chrom) %>%
  dplyr::distinct(break_id)
breaks.f = breaks %>% dplyr::filter(!(break_id %in% breaksites2break$break_id))
libsizes.f = breaks.f %>%
  dplyr::group_by(break_sample, break_chrom) %>%
  dplyr::summarise(library_size=n())


########################################################################################################################
########################################################################################################################
scale_factor_df = breaks %>%
    dplyr::group_by(breaks_file, break_sample, break_chrom) %>%
    dplyr::summarise(libsize=n(), .groups="keep") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(libsize_median=median(libsize)) %>%
    dplyr::mutate(score=libsize_median/libsize) %>%
    dplyr::select(breaks_file, score)

mm9_txdb.gtf = GenomicFeatures::makeTxDbFromGFF('data/mm9/mm9.refGene.gtf.gz', format="gtf")
mm9_txdb.tlen = as.data.frame(GenomicFeatures::transcripts(mm9_txdb.gtf, columns=c("tx_id", "tx_name", "gene_id"))) %>%
  dplyr::arrange(dplyr::desc(width)) %>%
  dplyr::distinct(gene_id, .keep_all=T)
mm9_txdb.tlen$gene_id = unlist(mm9_txdb.tlen$gene_id)
mm9_genes_df = as.data.frame(GenomicFeatures::transcripts(mm9_txdb.gtf)) %>%
  dplyr::inner_join(mm9_txdb.tlen %>% dplyr::select(tx_id, gene_id), by="tx_id") %>%
  dplyr::filter(seqnames %in% paste0("chr", 1:19)) %>%
  dplyr::inner_join(groseq, by=c("gene_id"="groseq_gene")) %>%
  dplyr::left_join(rdc, by=c("gene_id"="rdc_gene")) %>%
  dplyr::mutate(is_rdc=ifelse(is.na(rdc_cluster), "Not RDC", "RDC"))
# mm9_genes_df = as.data.frame(GenomicFeatures::genes(mm9_txdb.gtf)) %>%
#   dplyr::inner_join(groseq, by=c("gene_id"="groseq_gene"))

mm9_genes_ranges = GenomicRanges::makeGRangesFromDataFrame(mm9_genes_df)
breaks_bait = breaks.f
#%>% dplyr::filter(break_chrom==break_bait_chrom)
breaks_bait = breaks_bait %>% dplyr::select(-dplyr::matches("score")) %>% dplyr::inner_join(scale_factor_df, by="breaks_file")
breaks_bait_ranges = GenomicRanges::makeGRangesFromDataFrame(breaks_bait %>% dplyr::mutate(seqnames=break_chrom, start=break_start, end=break_end), keep.extra.columns=T)

hits = as(GenomicRanges::findOverlaps(mm9_genes_ranges, breaks_bait_ranges, ignore.strand=T), "List")
mm9_genes_df$junctions_count = sum(extractList(breaks_bait_ranges$score, hits))
mm9_genes_df$junctions_density = 1e6*(mm9_genes_df$junctions_count/mm9_genes_df$width)
# mm9_genes_df$is_rdc = ifelse(GenomicRanges::countOverlaps(mm9_genes_ranges, GenomicRanges::makeGRangesFromDataFrame(rdc_all %>% dplyr::mutate(seqnames=rdc_chrom, start=rdc_start, end=rdc_end), keep.extra.columns=T), ignore.strand=T), "Not RDC", "RDC")
mm9_genes_df = mm9_genes_df  %>%
  tidyr::crossing(data.frame(width_th=1:6 * 1e5)) %>%
  dplyr::filter(width >= width_th)

library(BSgenome)
genome = Biostrings::readDNAStringSet("data/mm9/mm9.fa.gz", "fasta")


data_df = list()
data_df[["rdc"]] = mm9_genes_df %>% dplyr::filter(groseq_chrom=="chr6" & gene_id %in% c("Sox6", "Npas6", "Nlgn1", "Auts2", "Exoc4", "Sdk1", "Ctnna2", "Grid2", "Anks1b", "Ptprt", "Pard3d", "Cdh13") & groseq_rpkm > 0.05 & is_rdc=="RDC" & width_th==2e5)
data_df[["none"]] = mm9_genes_df %>%
  dplyr::filter(groseq_rpkm > 0.05 & width_th==1e5 & is_rdc!="RDC") %>%
  tidyr::crossing(data_df[["rdc"]]  %>% dplyr::select(gene_id.rdc=gene_id, width.rdc=width)) %>%
  dplyr::arrange(abs(width-width.rdc)) %>%
  dplyr::distinct(gene_id, .keep_all=T) %>%
  dplyr::distinct(gene_id.rdc, .keep_all=T)
data_df[["none"]] %>% dplyr::select(gene_id, gene_id.rdc)
ranges = lapply(data_df, function(d) {
    r = GenomicRanges::makeGRangesFromDataFrame(d)
    names(r) = d$gene_id
    r
})
dna = lapply(ranges, function(r) Biostrings::getSeq(genome, r))


combinations = rbind(
    data.frame(t(combn(1:length(dna[["rdc"]]), 2))) %>% dplyr::mutate(type="rdc"),
    data.frame(t(combn(1:length(dna[["none"]]), 2))) %>% dplyr::mutate(type="none")
)

writeLines(c(data_df[["none"]]$gene_id, data_df[["rdc"]]$gene_id))

alignments = data.frame()
for(i in 1:nrow(combinations)) {
    writeLines(stringr::str_glue("{a}/{n}", a=which(i==1:nrow(combinations)), n=nrow(combinations)))
    i1 = combinations[i, "X1"]
    i2 = combinations[i, "X2"]
    type = combinations[i, "type"]
    writeXStringSet(dna[[type]][c(i1,i2)], 'data/rdc.fa')
    writeXStringSet(dna[[type]][i1], 'data/rdc1.fa')
    writeXStringSet(dna[[type]][i2], 'data/rdc2.fa')
    system("yass data/rdc1.fa data/rdc2.fa -d 3 -r 2 -o data/rdc.aln -E 1e-5 -O 100")
    a = readr::read_tsv("data/rdc.aln")
    a$type = type
    a$gene_id1 = data_df[[type]]$gene_id[i1]
    a$gene_id2 = data_df[[type]]$gene_id[i2]
    alignments = rbind(alignments, a)
}

alignments %>% dplyr::arrange(dplyr::desc(`bit score`)) %>% View()

alignments.sum = alignments %>%
  dplyr::group_by(type, gene_id1, gene_id2) %>%
  dplyr::summarize(bitscore=quantile(`bit score`, 1))

ggplot(alignments.sum) +
  ggridges::geom_density_ridges(aes(x=bitscore, y=type))
# plotly::ggplotly()

p = ggplot(mm9_genes_df %>% dplyr::filter(groseq_rpkm > 0.05), aes(x=groseq_rpkm, y=junctions_density)) +
  geom_point(aes(color=is_rdc, gene=gene_id, size=width, chrom=rdc_chrom), alpha=0.4) +
  geom_smooth(aes(color=is_rdc), method="lm") +
  facet_wrap(~is_rdc) +
  scale_size_continuous(breaks=seq(0, 1e6, 1e4)) +
  labs(y="Breaks (1/M)", x="RPKM") +
  # coord_cartesian(ylim=c(0,2e5)) +
  facet_wrap(~width_th, scales="free")
p
pp = plotly::ggplotly(p)
htmlwidgets::saveWidget(plotly::as_widget(pp), "rdc_vs_none_density2rpkm_notarget_allchrom.html")

p = ggplot(mm9_genes_df) +
  geom_density(aes(x=groseq_rpkm, color=is_rdc)) +
  labs(x="RPKM") +
  facet_wrap(~width_th, scales="free")
pp = plotly::ggplotly(p)
htmlwidgets::saveWidget(plotly::as_widget(pp), "rdc_vs_none_density_notarget_allchrom.html")


rdc2groseq = rdc %>% dplyr::inner_join(groseq, by=c("rdc_gene"="groseq_gene"))
rdc2groseq_ranges = GenomicRanges::makeGRangesFromDataFrame(rdc2groseq %>% dplyr::mutate(seqnames=rdc_chrom, start=rdc_start, end=rdc_end))

hits = as(findOverlaps(rdc2groseq_ranges, breaks_bait_ranges, ignore.strand=T), "List")
rdc2groseq$junctions_count = sum(extractList(breaks_bait_ranges$score, hits))
rdc2groseq_grouped = rdc2groseq %>%
  dplyr::mutate(junctions_density=(junctions_count/groseq_length)*1e6) %>%
  dplyr::arrange(dplyr::desc(groseq_length)) %>%
  dplyr::group_by(rdc_cluster) %>%
  dplyr::mutate(groseq_rpkm=min(groseq_rpkm)) %>%
  dplyr::slice(1)

ggplot(rdc2groseq_grouped, aes(x=groseq_rpkm, y=junctions_density)) +
  geom_point(aes(size=rdc_length, gene=rdc_gene))

########################################################################################################################
########################################################################################################################



#
# Remove genes shorter than 800kb
#
groseq.f = groseq %>%
  dplyr::filter(groseq_rpkm>0.05) # %>%
  #dplyr::filter(dplyr::between(groseq_length, 10e3, 50e3)) # %>%
  #dplyr::filter(groseq_length>=0e3)
groseq_ranges.f = with(groseq.f, IRanges::IRanges(start=groseq_start, end=groseq_end, groseq_id=groseq_id, groseq_chrom=groseq_chrom))
breaks_ranges.f = with(breaks.f, IRanges::IRanges(start=break_start, end=break_end, break_id=break_id, break_chrom=break_chrom))
groseq2break.f.map = as.data.frame(IRanges::mergeByOverlaps(groseq_ranges.f, breaks_ranges.f)) %>%
  dplyr::filter(break_chrom==groseq_chrom) %>%
  dplyr::select(groseq_id=groseq_id, break_id=break_id)

groseq2break.f = groseq.f %>%
  dplyr::inner_join(groseq2break.f.map, by="groseq_id") %>%
  dplyr::inner_join(breaks.f, by="break_id") %>%
  dplyr::group_by(break_sample, break_bait_chrom, break_exp_condition, groseq_chrom, groseq_gene, groseq_strand, groseq_length, groseq_rpkm) %>%
  dplyr::summarise(breaks_count=n()) %>%
  dplyr::inner_join(libsizes.f, by=c("break_sample", "groseq_chrom"="break_chrom")) %>%
  dplyr::left_join(rdc, by=c("groseq_gene"="rdc_gene")) %>%
  dplyr::mutate(y=(breaks_count/groseq_length) / library_size, x=groseq_length, rdc_group=tidyr::replace_na(rdc_group, "N/A"), rdc_cluster=tidyr::replace_na(rdc_cluster, "N/A"), rdc_length=tidyr::replace_na(rdc_length, 0)) %>%
  dplyr::filter(break_bait_chrom==groseq_chrom)
groseq2break_fc.f = groseq2break.f %>%
  dplyr::group_by(groseq_gene) %>%
  dplyr::summarise(
    s=which(break_exp_condition=="Sample"),
    c=which(break_exp_condition=="Control"),
    breaks_density_fc=log2((breaks_count[s]/breaks_count[c]) * (library_size[c]/library_size[s])))
groseq2break.f = groseq2break.f %>%
  dplyr::filter(break_exp_condition=="Sample") %>%
  dplyr::inner_join(groseq2break_fc.f, by="groseq_gene")


ggplot(groseq2break.f %>% dplyr::filter(groseq_length<80e4), aes(y=breaks_density_fc, x=groseq_length/1000)) +
  geom_point(aes(size=groseq_rpkm, color=rdc_group, name=groseq_gene, rdc_cluster=rdc_cluster, rdc_group=rdc_group), alpha=0.5) +
  geom_smooth() +
  geom_hline(yintercept=0)

ggplot(groseq2break.f, aes(y=(breaks_count/groseq_length) / library_size, x=groseq_length/1000)) +
  geom_point(aes(size=groseq_rpkm, color=rdc_group, name=groseq_gene, rdc_cluster=rdc_cluster, rdc_group=rdc_group)) +
  geom_vline(xintercept=800) +
  geom_vline(xintercept=80) +
  #ggrepel::geom_text_repel(aes(label=groseq_gene)) +
  labs(x="gene length", y="junction density") +
  coord_cartesian(ylim=c(0, 2e-7))


#p = plotly::ggplotly()
#htmlwidgets::saveWidget(plotly::as_widget(p), "viven-0.10.html")

plotly::ggplotly()


validated_offtargets_overlapping = validated_offtargets %>%
dplyr::inner_join(as.data.frame(IRanges::mergeByOverlaps(offtargets_ranges, control_junctions_ranges)), by=c("offtarget_start"="control_junctions_ranges.start", "offtarget_end"="validated_offtargets_ranges.end")) %>%
dplyr::mutate(offtarget_validated=T) %>%
dplyr::select(bait_name, bait_chrom, offtarget_is_primary, offtarget_chrom, offtarget_strand, offtarget_start, offtarget_end, offtarget_sequence, offtarget_start, offtarget_end, primer_alignment_start=primers_targets_ranges.start, primer_alignment_end=primers_targets_ranges.end, offtarget_validated)



validated_offtargets_overlapping = validated_offtargets %>%
dplyr::inner_join(as.data.frame(IRanges::mergeByOverlaps(offtargets_ranges, control_junctions_ranges)), by=c("offtarget_start"="control_junctions_ranges.start", "offtarget_end"="validated_offtargets_ranges.end")) %>%
dplyr::mutate(offtarget_validated=T) %>%
dplyr::select(bait_name, bait_chrom, offtarget_is_primary, offtarget_chrom, offtarget_strand, offtarget_start, offtarget_end, offtarget_sequence, offtarget_start, offtarget_end, primer_alignment_start=primers_targets_ranges.start, primer_alignment_end=primers_targets_ranges.end, offtarget_validated)




