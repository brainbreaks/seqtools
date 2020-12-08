library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(plotly)
library(GenomeInfoDb)
library(BSgenome)
library(Gviz)
# BiocManager::install("BSgenome.Mmusculus.UCSC.mm9")
# BiocManager::install("TxDb.Mmusculus.UCSC.mm9.knownGene")


plot_breaks = function(breaks2groseq, gene_name, mm9, mm9_txdb.gtf) {
  breaks_display = breaks2groseq %>%
    dplyr::filter(groseq_gene==gene_name) %>%
    dplyr::mutate(group=break_strand)

  breaks_ranges_display.control_plus = GenomicRanges::makeGRangesFromDataFrame(breaks_display %>% dplyr::filter(break_exp_condition=="Control" & break_strand=="+"), start.field="break_start", end.field="break_end", seqnames.field="break_chrom", seqinfo=mm9, strand.field="break_strand", keep.extra.columns=T)
  breaks_ranges_display.control_minus = GenomicRanges::makeGRangesFromDataFrame(breaks_display %>% dplyr::filter(break_exp_condition=="Control" & break_strand=="-"), start.field="break_start", end.field="break_end", seqnames.field="break_chrom", seqinfo=mm9, strand.field="break_strand", keep.extra.columns=T)
  breaks_ranges_display.sample_plus = GenomicRanges::makeGRangesFromDataFrame(breaks_display %>% dplyr::filter(break_exp_condition=="Sample" & break_strand=="+"), start.field="break_start", end.field="break_end", seqnames.field="break_chrom", seqinfo=mm9, strand.field="break_strand", keep.extra.columns=T)
  breaks_ranges_display.sample_minus = GenomicRanges::makeGRangesFromDataFrame(breaks_display %>% dplyr::filter(break_exp_condition=="Sample" & break_strand=="-"), start.field="break_start", end.field="break_end", seqnames.field="break_chrom", seqinfo=mm9, strand.field="break_strand", keep.extra.columns=T)

  breaks_count = breaks_display %>%
    dplyr::group_by(break_exp_condition) %>%
    dplyr::summarize(n=n()) %>%
    tibble::column_to_rownames("break_exp_condition")

  breaks_display.density = breaks_display %>%
    dplyr::group_by(break_exp_condition, groseq_gene, groseq_start, groseq_end, break_chrom, break_strand, display_start, display_end) %>%
    dplyr::do((function(z) {
      zz<<-z
        from = min(z$display_start)
        to = max(z$display_end)
        n = round((to - from)/10000)
        d = density(z$break_start, from=from, to=to, bw=2e4, n=n)
        d.bin = d$x[2L] - d$x[1L]
        d.count = zoo::rollsum(d$y, 2)*d.bin*nrow(z)/2

      data.frame(density_start=d$x[-length(d$x)], density_end=d$x[-1], density_value=d.count*ifelse(z$break_strand[1]=="+",1,-1))

        #x = seq(d$x[1], to=d$x[length(d$x)-1], by=10000)
        #y = stats::approx(d$x[-length(d$x)], d.count, xout=x)$y
        #
        #l = length(x)
        #data.frame(density_start=x[-l], density_end=x[-1], density_value=nrow(z)*ifelse(z$break_strand[1]=="+",1,-1) * y[-1])
    })(.)) %>%
    reshape2::dcast(break_chrom+break_exp_condition+density_start+density_end ~ break_strand, value.var="density_value")


  breaks_ranges_display.control_density = GenomicRanges::makeGRangesFromDataFrame(breaks_display.density %>% dplyr::filter(break_exp_condition=="Control"), start.field="density_start", end.field="density_end", seqnames.field="break_chrom", seqinfo=mm9, strand.field="break_strand", keep.extra.columns=T)
  breaks_ranges_display.sample_density = GenomicRanges::makeGRangesFromDataFrame(breaks_display.density %>% dplyr::filter(break_exp_condition=="Sample"), start.field="density_start", end.field="density_end", seqnames.field="break_chrom", seqinfo=mm9, strand.field="break_strand", keep.extra.columns=T)

  max_density = max(abs(c(breaks_display.density$`+`, breaks_display.density$`-`)))


  axisTrack = GenomeAxisTrack()
  breaksControlDensityTrack = Gviz::DataTrack(range=breaks_ranges_display.control_density, data=t(cbind(breaks_ranges_display.control_density$`-`, breaks_ranges_display.control_density$`+`)), type="histogram", name="Control\njunctions\ndensity", groups=c("+", "-"), col=c("red","blue"), ylim=c(-max_density, max_density), legend=F)
  breaksSampleDensityTrack = Gviz::DataTrack(range=breaks_ranges_display.sample_density, data=t(cbind(breaks_ranges_display.sample_density$`-`, breaks_ranges_display.sample_density$`+`)), type="histogram", name="Sample\njunctions\ndensity", groups=c("+", "-"), col=c("red","blue"), ylim=c(-max_density, max_density))
  genesTrack = GeneRegionTrack(mm9_txdb.gtf, genome="mm9", chromosome=breaks_display$groseq_chrom[1], name="refGene", transcriptAnnotation='gene', showAxis=T, col="black")
  breaksSamplePlusTrack = Gviz::AnnotationTrack(range=breaks_ranges_display.sample_plus, name="A + (Sample)", col="red", stacking = "dense", alpha=0.5, alpha.title=1, stackHeight=1, shape="box", collapse=F, shape="box", col.line="transparent")
  breaksSampleMinusTrack = Gviz::AnnotationTrack(range=breaks_ranges_display.sample_minus, name="A - (Sample)", col="blue", stacking = "dense", alpha=0.5, alpha.title=1, stackHeight=1, shape="box", collapse=F, shape="box", col.line="transparent")
  breaksControlPlusTrack = Gviz::AnnotationTrack(range=breaks_ranges_display.control_plus, name="A + (Control)", col="red", stacking = "dense", alpha=0.5, alpha.title=1, stackHeight=1, shape="box", collapse=F, shape="box", col.line="transparent")
  breaksControlMinusTrack = Gviz::AnnotationTrack(range=breaks_ranges_display.control_minus, name="A - (Control)", col="blue", stacking = "dense", alpha=0.5, alpha.title=1, stackHeight=1, shape="box", collapse=F, shape="box", col.line="transparent")
  plotTracks(
    main=stringr::str_glue("{gene} (junctions: {ctrl} in control, {smpl} in sample)", gene=gene_name, ctrl=breaks_count["Control","n"], smpl=breaks_count["Sample","n"]),
    c(axisTrack, genesTrack, breaksControlDensityTrack, breaksSampleDensityTrack, breaksSamplePlusTrack, breaksSampleMinusTrack),
    from=min(breaks_display$display_start), to=max(breaks_display$display_end),
    col.axis=col.axis, fontcolor.title="black", background.title="transparent", fontsize=16, title.width = 1.6, cex.group=1.2)
  #
  #availableDisplayPars(genesTrack)
  #displayPars(breaksSampleDensityTrack)
  #displayPars(breaksSampleDensityTrack)[grepl("font", names(displayPars(breaksSampleDensityTrack)))]
}



breaksites = readr::read_tsv("data/offtargets_predicted.tsv") %>%
   dplyr::filter(offtarget_mismatches==0) %>%
   dplyr::mutate(offtarget_id=1:n())

breaksites = readr::read_tsv("data/breaksites.tsv") %>%
   dplyr::mutate(breaksite_start=breaksite_pos-1e6, breaksite_end=breaksite_pos+1e6) %>%
   dplyr::mutate(breaksite_id=1:n())

groseq_cols=cols(groseq_chrom=col_character(), groseq_start=col_double(), groseq_end=col_double(), groseq_strand=col_character(), groseq_transcript=col_character(), groseq_gene=col_character(), groseq_reads=col_double(), groseq_length=col_double(), groseq_rpkm=col_double())
groseq = readr::read_tsv("data/Xrcc4_p53_double.transcript.longest.cnt", col_names=names(groseq_cols$cols), col_types=groseq_cols) %>%
   dplyr::mutate(groseq_id=1:n()) %>%
   dplyr::mutate(groseq_gene=dplyr::case_when(
        groseq_gene=="Large" ~ "Large1",
        groseq_gene=="Naalad2" ~ "Naaladl2",
        T~groseq_gene)) %>%
  dplyr::mutate(display_start=groseq_start-2e5, display_end=groseq_end+2e5)

breaks = data.frame()
breaks_cols = cols(break_chrom=col_character(), break_start=col_double(), break_end=col_double(), break_name=col_character(), break_score=col_character(), break_strand=col_character())
for(breaks_path in list.files("data/breaks", pattern="*.bed", full.names=T)) {
    chr = gsub("(chr)([^_]+).*$", "\\L\\1\\2", basename(breaks_path), perl=T, ignore.case=T)
    condition = ifelse(grepl("DMSO", breaks_path), "Control", "Sample")
    sample = gsub("\\.bed", "", basename(breaks_path))


    x = readr::read_tsv(breaks_path, col_names=names(breaks_cols$cols), col_types=breaks_cols) %>%
      dplyr::mutate(break_bait_chrom=chr, break_exp_condition=condition, break_sample=sample, break_end=break_start+1)
    breaks = rbind(breaks, x)
}
breaks = breaks %>% dplyr::mutate(break_id=1:n())


groseq_ranges = with(groseq, IRanges::IRanges(start=display_start, end=display_end, groseq_id=groseq_id, groseq_chrom=groseq_chrom))
breaks_ranges = with(breaks, IRanges::IRanges(start=break_start, end=break_end, break_id=break_id, break_chrom=break_chrom))
groseq2break.map = as.data.frame(IRanges::mergeByOverlaps(groseq_ranges, breaks_ranges)) %>%
  dplyr::filter(break_chrom==groseq_chrom) %>%
  dplyr::select(groseq_id=groseq_id, break_id=break_id)

rdc = readr::read_tsv("data/rdc_pnas.tsv") %>%
  dplyr::mutate(rdc_length=rdc_end-rdc_start) %>%
  tidyr::separate_rows(rdc_gene, sep=" *, *") %>%
  dplyr::filter(rdc_gene != "--") %>%
  dplyr::group_by(rdc_gene) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(rdc_id=1:n())

libsizes = breaks %>%
  dplyr::group_by(break_sample, break_chrom) %>%
  dplyr::summarise(library_size=n())

breaks2groseq = breaks %>%
  dplyr::filter(break_bait_chrom==break_chrom) %>%
  dplyr::inner_join(libsizes, by=c("break_sample", "break_chrom")) %>%
  dplyr::inner_join(groseq2break.map, by="break_id") %>%
  dplyr::inner_join(groseq, by="groseq_id") %>%
  dplyr::left_join(rdc, by=c("groseq_gene"="rdc_gene")) %>%
  dplyr::group_by(break_sample, groseq_gene) %>%
  dplyr::mutate(gene_breaks_count=sum(dplyr::between(break_start, groseq_start, groseq_end))) %>%
  dplyr::ungroup()

genes_of_interest = breaks2groseq %>%
  dplyr::filter(!is.na(rdc_group) & rdc_group==1 & gene_breaks_count>50) %>%
  dplyr::distinct(groseq_gene) %>%
  .$groseq_gene

#
# Visualize with Gviz
#
mm9 = GenomeInfoDb::Seqinfo(genome="mm9")
mm9_txdb.gtf = GenomicFeatures::makeTxDbFromGFF('data/mm9/mm9.refGene.gtf.gz', format="gtf")
pdf("reports/breaks_shift_examples.pdf", width=10, height=6)
for(gene_name in genes_of_interest) {
  plot_breaks(breaks2groseq=breaks2groseq, gene_name=gene_name, mm9=mm9, mm9_txdb.gtf=mm9_txdb.gtf)
}
dev.off()