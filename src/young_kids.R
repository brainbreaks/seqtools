library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(plotly)
library(GenomeInfoDb)
library(BSgenome)
library(Gviz)
library(karyoploteR)
library(regioneR)
library(zoo)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
# BiocManager::install("BSgenome.Mmusculus.UCSC.mm9")
# BiocManager::install("TxDb.Mmusculus.UCSC.mm9.knownGene")


plot_breaks = function(breaks2groseq, gene_name, mm9, mm9_txdb.gtf) {
  breaks_display = breaks2groseq %>%
    dplyr::filter(groseq_gene==gene_name & break_bait_chrom==break_chrom) %>%
    dplyr::mutate(group=break_strand)

  breaks_ranges_display.control_plus = GenomicRanges::makeGRangesFromDataFrame(breaks_display %>% dplyr::filter(break_exp_condition=="Control" & break_strand=="+"), start.field="break_start", end.field="break_end", seqnames.field="break_chrom", seqinfo=mm9, strand.field="break_strand", keep.extra.columns=T)
  breaks_ranges_display.control_minus = GenomicRanges::makeGRangesFromDataFrame(breaks_display %>% dplyr::filter(break_exp_condition=="Control" & break_strand=="-"), start.field="break_start", end.field="break_end", seqnames.field="break_chrom", seqinfo=mm9, strand.field="break_strand", keep.extra.columns=T)
  breaks_ranges_display.sample_plus = GenomicRanges::makeGRangesFromDataFrame(breaks_display %>% dplyr::filter(break_exp_condition=="Sample" & break_strand=="+"), start.field="break_start", end.field="break_end", seqnames.field="break_chrom", seqinfo=mm9, strand.field="break_strand", keep.extra.columns=T)
  breaks_ranges_display.sample_minus = GenomicRanges::makeGRangesFromDataFrame(breaks_display %>% dplyr::filter(break_exp_condition=="Sample" & break_strand=="-"), start.field="break_start", end.field="break_end", seqnames.field="break_chrom", seqinfo=mm9, strand.field="break_strand", keep.extra.columns=T)

  breaks_display.density = breaks_display %>%
    dplyr::group_by(break_exp_condition, groseq_gene, groseq_start, groseq_end, break_chrom, break_strand, display_start, display_end) %>%
    dplyr::do((function(z) {
      zz<<-z
        d = density(z$break_start, from=min(z$display_start), to=max(z$display_end), bw=2e4)
        x = seq(min(z$display_start), to=max(z$display_end), by=10000)
        y = stats::approx(d$x, d$y, xout=x)$y
        l = length(x)
        data.frame(density_start=x[-l], density_end=x[-1], density_value=nrow(z)*ifelse(z$break_strand[1]=="+",1,-1) * y[-1])
    })(.)) %>%
    reshape2::dcast(break_chrom+break_exp_condition+density_start+density_end ~ break_strand, value.var="density_value")
  breaks_ranges_display.control_density = GenomicRanges::makeGRangesFromDataFrame(breaks_display.density %>% dplyr::filter(break_exp_condition=="Control"), start.field="density_start", end.field="density_end", seqnames.field="break_chrom", seqinfo=mm9, strand.field="break_strand", keep.extra.columns=T)
  breaks_ranges_display.sample_density = GenomicRanges::makeGRangesFromDataFrame(breaks_display.density %>% dplyr::filter(break_exp_condition=="Sample"), start.field="density_start", end.field="density_end", seqnames.field="break_chrom", seqinfo=mm9, strand.field="break_strand", keep.extra.columns=T)

  max_density = max(abs(c(breaks_display.density$`+`, breaks_display.density$`-`)))

  axisTrack = GenomeAxisTrack()
  breaksControlDensityTrack = Gviz::DataTrack(range=breaks_ranges_display.control_density, data=t(cbind(breaks_ranges_display.control_density$`-`, breaks_ranges_display.control_density$`+`)), type="histogram", name="Breaks density (Control)", groups=c("+", "-"), col=c("red","blue"), ylim=c(-max_density, max_density))
  breaksSampleDensityTrack = Gviz::DataTrack(range=breaks_ranges_display.sample_density, data=t(cbind(breaks_ranges_display.sample_density$`-`, breaks_ranges_display.sample_density$`+`)), type="histogram", name="Breaks density (Control)", groups=c("+", "-"), col=c("red","blue"), ylim=c(-max_density, max_density))
  genesTrack = GeneRegionTrack(mm9_txdb.gtf, genome="mm9", chromosome=breaks_display$groseq_chrom[1], name="refGene", transcriptAnnotation='gene') # collapseTranscripts = "meta
  breaksSamplePlusTrack = Gviz::AnnotationTrack(range=breaks_ranges_display.sample_plus, name="+ (Sample)", col="red", stacking = "dense", alpha=0.5, alpha.title=1, stackHeight=1, shape="box", collapse=F, shape="box", col.line="transparent")
  breaksSampleMinusTrack = Gviz::AnnotationTrack(range=breaks_ranges_display.sample_minus, name="- (Sample)", col="blue", stacking = "dense", alpha=0.5, alpha.title=1, stackHeight=1, shape="box", collapse=F, shape="box", col.line="transparent")
  breaksControlPlusTrack = Gviz::AnnotationTrack(range=breaks_ranges_display.control_plus, name="+ (Control)", col="red", stacking = "dense", alpha=0.5, alpha.title=1, stackHeight=1, shape="box", collapse=F, shape="box", col.line="transparent")
  breaksControlMinusTrack = Gviz::AnnotationTrack(range=breaks_ranges_display.control_minus, name="- (Control)", col="blue", stacking = "dense", alpha=0.5, alpha.title=1, stackHeight=1, shape="box", collapse=F, shape="box", col.line="transparent")
  plotTracks(c(axisTrack, genesTrack, breaksControlDensityTrack, breaksSampleDensityTrack, breaksSamplePlusTrack, breaksSampleMinusTrack), from=min(breaks_display$display_start), to=max(breaks_display$display_end)) #, from=77.55e6, to=77.63e
}



mm9 = GenomeInfoDb::Seqinfo(genome="mm9")

# mm9_txdb.ucsc = makeTxDbFromUCSC(genome="mm9", tablename="refGene",)
mm9_txdb.gtf = GenomicFeatures::makeTxDbFromGFF('data/mm9/mm9.refGene.gtf.gz', format="gtf")

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


breaks2groseq = breaks %>%
  dplyr::inner_join(groseq2break.map, by="break_id") %>%
  dplyr::inner_join(groseq, by="groseq_id")




#
# Visualize with Gviz
#

mm9 = GenomeInfoDb::Seqinfo(genome="mm9")
mm9_txdb.gtf = GenomicFeatures::makeTxDbFromGFF('data/mm9/mm9.refGene.gtf.gz', format="gtf")
plot_breaks(breaks2groseq=breaks2groseq, gene_name="Ctnna2", mm9=mm9, mm9_txdb.gtf=mm9_txdb.gtf)
plot_breaks(breaks2groseq=breaks2groseq, gene_name="Csmd3", mm9=mm9, mm9_txdb.gtf=mm9_txdb.gtf)




# x = breaks_display %>%
#   dplyr::mutate(fixed_strand="+", lo=ifelse(break_strand=="+", 0, 2), hi=ifelse(break_strand=="+", 1, 1)) %>%
#   reshape2::melt(measure.vars=c("lo", "hi")) %>%
#   reshape2::dcast(break_chrom+break_strand+break_start+break_end ~ break_strand+variable, value.var="value", fun.aggregate=function(z) ifelse(length(z)==0, NA_real_, z))
# xx = GenomicRanges::makeGRangesFromDataFrame(x, start.field="break_start", end.field="break_end", seqnames.field="break_chrom", seqinfo=mm9, strand.field="fixed_strand", keep.extra.columns=T)
# breaksTrack = Gviz::DataTrack(range=xx, type="S", data=t(x[, c("-_lo", "-_hi", "+_lo", "+_hi")]), name="+/-", groups=c("-","-", "+","+"), ifelse(x$break_strand=="+", "red", "blue"), alpha=0.5, alpha.title=1)
# plotTracks(c(axisTrack, genesTrack, breaksDensityTrack,  breaksXTrack), from=min(breaks_display$display_start), to=max(breaks_display$display_end)) #, from=77.55e6, to=77.63e6



gcTrack = Gviz::UcscTrack(genome="mm9", chromosome=breaks_display$groseq_chrom[1],
                       track = "GC Percent", table = "gc5Base",
                       from=min(breaks_display$groseq_start), to=max(breaks_display$groseq_end),
                       trackType = "DataTrack",
                       start = "start", end = "end", data = "score",
                       type = "hist", window = -1, windowSize = 1500,
                       fill.histogram = "black", col.histogram = "black",
                       ylim = c(30, 70), name = "GC Percent")

gcContent <- UcscTrack(genome = "mm9", chromosome = breaks_display$groseq_chrom[1],
                       track = "GC Percent", table = "gc5Base",
                       from = from, to = to,
                       trackType = "DataTrack",
                       start = "start", end = "end", data = "score",
                       type = "hist", window = -1, windowSize = 1500,
                       fill.histogram = "black", col.histogram = "black",
                       ylim = c(30, 70), name = "GC Percent")

plotTracks(c(axisTrack, genesTrack, breaksDensityTrack, breaksPlusTrack, breaksMinusTrack), from=min(breaks_display$groseq_start), to=max(breaks_display$groseq_end)) #, from=77.55e6, to=77.63e6

ranges(genesTrack)$gene

# mcols(genes(mm9_txdb))$gene_id
# mcols(genes(mm9_txdb))