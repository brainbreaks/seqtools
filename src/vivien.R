library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(plotly)


breaksites = readr::read_tsv("data/offtargets_predicted.tsv") %>%
   dplyr::filter(offtarget_mismatches==0) %>%
   dplyr::mutate(offtarget_id=1:n())

breaksites = readr::read_tsv("breaksites.tsv") %>%
   dplyr::mutate(breaksite_start=breaksite_pos-1e6, breaksite_end=breaksite_pos+1e6) %>%
   dplyr::mutate(breaksite_id=1:n())

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
      dplyr::mutate(break_bait_chrom=chr, break_exp_condition=condition, break_sample=sample)
    breaks = rbind(breaks, x)
}
breaks = breaks %>% dplyr::mutate(break_id=1:n())


#
# Load RDC
#
rdc = readr::read_tsv("data/rdc_pnas.tsv") %>%
  dplyr::mutate(rdc_length=rdc_end-rdc_start) %>%
  tidyr::separate_rows(rdc_gene, sep=" *, *") %>%
  dplyr::filter(rdc_gene != "--") %>%
  dplyr::group_by(rdc_gene) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(rdc_id=1:n())


breaks_ranges = with(breaks,  IRanges::IRanges(start=break_start, end=break_end, break_id=break_id, break_chrom=break_chrom))
breaksites2break = as.data.frame(IRanges::mergeByOverlaps(breaksites_ranges, breaks_ranges)) %>%
  dplyr::filter(break_chrom==breaksite_chrom) %>%
  dplyr::distinct(break_id)
breaks.f = breaks %>% dplyr::filter(!(break_id %in% breaksites2break$break_id))
libsizes.f = breaks.f %>%
  dplyr::group_by(break_sample, break_chrom) %>%
  dplyr::summarise(library_size=n())

#
# Remove breaks close to bait
#
breaksites_ranges = with(breaksites, IRanges::IRanges(start=breaksite_start, end=breaksite_end, breaksite_id=breaksite_id, breaksite_chrom=breaksite_chrom))
breaks_ranges = with(breaks,  IRanges::IRanges(start=break_start, end=break_end, break_id=break_id, break_chrom=break_chrom))
breaksites2break = as.data.frame(IRanges::mergeByOverlaps(breaksites_ranges, breaks_ranges)) %>%
  dplyr::filter(break_chrom==breaksite_chrom) %>%
  dplyr::distinct(break_id)
breaks.f = breaks %>% dplyr::filter(!(break_id %in% breaksites2break$break_id))
libsizes.f = breaks.f %>%
  dplyr::group_by(break_sample, break_chrom) %>%
  dplyr::summarise(library_size=n())


#
# Remove genes shorter than 800kb
#
groseq.f = groseq %>%
  dplyr::filter(groseq_rpkm>0.1) %>%
  #dplyr::filter(dplyr::between(groseq_length, 10e3, 50e3)) # %>%
  dplyr::filter(groseq_length>=800e3)
groseq_ranges.f = with(groseq.f, IRanges::IRanges(start=groseq_start, end=groseq_end, groseq_id=groseq_id, groseq_chrom=groseq_chrom))
breaks_ranges.f = with(breaks.f, IRanges::IRanges(start=break_start, end=break_end, break_id=break_id, break_chrom=break_chrom))
groseq2break.f.map = as.data.frame(IRanges::mergeByOverlaps(groseq_ranges.f, breaks_ranges.f)) %>%
  dplyr::filter(break_chrom==groseq_chrom) %>%
  dplyr::select(groseq_id=groseq_id, break_id=break_id)

groseq2break.f = groseq.f %>%
  dplyr::inner_join(groseq2break.f.map, by=c("groseq_id")) %>%
  dplyr::inner_join(breaks.f, by=c("break_id")) %>%
  dplyr::group_by(break_sample, break_bait_chrom, break_exp_condition, groseq_chrom, groseq_gene, groseq_strand, groseq_length, groseq_rpkm) %>%
  dplyr::summarise(breaks_count=n()) %>%
  dplyr::inner_join(libsizes.f, by=c("break_sample", "groseq_chrom"="break_chrom")) %>%
  dplyr::left_join(rdc, by=c("groseq_gene"="rdc_gene")) %>%
  dplyr::mutate(y=(breaks_count/groseq_length) / library_size, x=groseq_length, rdc_group=tidyr::replace_na(rdc_group, "N/A"), rdc_cluster=tidyr::replace_na(rdc_cluster, "N/A"), rdc_length=tidyr::replace_na(rdc_length, 0)) %>%
  dplyr::filter(break_exp_condition=="Sample" & break_bait_chrom==groseq_chrom)


ggplot(groseq2break.f, aes(x=x, y=y)) +
  geom_point(aes(size=rdc_length, color=rdc_group, name=groseq_gene, rdc_cluster=rdc_cluster, rdc_group=rdc_group)) +
  ggrepel::geom_text_repel(aes(label=groseq_gene)) +
  labs(x="gene length", y="junction density") +
  coord_cartesian(ylim=c(0, 0.0000001))


p = plotly::ggplotly()
htmlwidgets::saveWidget(plotly::as_widget(p), "viven.html")
plotly::ggplotly()



validated_offtargets_overlapping = validated_offtargets %>%
dplyr::inner_join(as.data.frame(IRanges::mergeByOverlaps(offtargets_ranges, control_junctions_ranges)), by=c("offtarget_start"="control_junctions_ranges.start", "offtarget_end"="validated_offtargets_ranges.end")) %>%
dplyr::mutate(offtarget_validated=T) %>%
dplyr::select(bait_name, bait_chrom, offtarget_is_primary, offtarget_chrom, offtarget_strand, offtarget_start, offtarget_end, offtarget_sequence, offtarget_start, offtarget_end, primer_alignment_start=primers_targets_ranges.start, primer_alignment_end=primers_targets_ranges.end, offtarget_validated)



validated_offtargets_overlapping = validated_offtargets %>%
dplyr::inner_join(as.data.frame(IRanges::mergeByOverlaps(offtargets_ranges, control_junctions_ranges)), by=c("offtarget_start"="control_junctions_ranges.start", "offtarget_end"="validated_offtargets_ranges.end")) %>%
dplyr::mutate(offtarget_validated=T) %>%
dplyr::select(bait_name, bait_chrom, offtarget_is_primary, offtarget_chrom, offtarget_strand, offtarget_start, offtarget_end, offtarget_sequence, offtarget_start, offtarget_end, primer_alignment_start=primers_targets_ranges.start, primer_alignment_end=primers_targets_ranges.end, offtarget_validated)




