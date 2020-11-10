library(dplyr)
library(ggplot2)
library(scales)
library(reshape2)
library(ggrepel)
library(readr)

Sys.setenv("plotly_username" = "sandrejev")
Sys.setenv("plotly_api_key" = "AxuixdXZeEivct7XGi2g")


#
# Thresholds
#
th.ccf_pvalue = 0.05
th.ccf_minlag = 0.1
th.ccf_maxlag = 0.1
th.ccf_correlation = 0.5
th.min_bin_breaks = 5


# TODO: why there is no Y chromosome. 
# TODO: Is it true that in results chrXM => chrY and chrXF => chrX
breaksites = readr::read_tsv("data/breaksites.tsv") %>%
  dplyr::mutate(breaksite_full=paste0(breaksite_bait, breaksite_priming_direction))

#
# Read result files
#
results_chr_map = c("chrxf"="chrX", "chrxm"="chrY")
results =  data.frame()
for(results_path in list.files("results_parallel3/", full.names=T)) {
  writeLines(paste0("Reading file ", basename(results_path), "..."))
  chr = tolower(gsub(".*(chr[^_]+)_.*", "\\1", basename(results_path), ignore.case=T))
  chr = ifelse(chr %in% names(results_chr_map), results_chr_map[chr], chr)
  results.chr = readr::read_tsv(results_path, col_types=cols(
    bin_start = col_double(),
    bin_end = col_double(),
    bin_strand = col_character(),
    gene_chrom = col_character(),
    gene_name = col_character(),
    gene_strand = col_character(),
    gene_start = col_double(),
    gene_end = col_double(),
    gene_region_start = col_double(),
    gene_region_end = col_double(),
    attributes = col_character()
  )) %>% dplyr::mutate(breaksite_chrom=chr)
  results = dplyr::bind_rows(results, results.chr)
}

results = results  %>%
  dplyr::filter(gene_name!=".") %>%
  tidyr::extract(attributes, into=c("breaks"), regex="breaks=([^;]+)") %>%
  dplyr::mutate(gene_full=paste0(gene_chrom, ":", gene_name)) %>%
  dplyr::mutate(gene_strand_name=ifelse(gene_strand=="+", "positive", "negative")) %>%
  dplyr::mutate(breaks=as.numeric(breaks), breaks_sign=ifelse(gene_strand=="+", breaks, -breaks)) %>%
  dplyr::group_by(breaksite_chrom, gene_full, gene_start, gene_end) %>% 
  dplyr::mutate(gene_length=max(gene_end)-min(gene_start)) %>%
  dplyr::mutate(breaks_max=max(breaks)) %>%
  dplyr::filter(breaks_max>=th.min_bin_breaks)

#
# Display relation between distance to the breaksite and amount of breaks
#
results.dist_sum = results %>%
  dplyr::filter(breaksite_chrom==gene_chrom) %>%
  dplyr::inner_join(breaksites, by=c("breaksite_chrom")) %>%
  dplyr::group_by(breaksite_full, breaksite_bait, breaksite_priming_direction, gene_chrom, gene_name, gene_full, gene_strand, gene_start, gene_end, gene_length, bin_strand) %>%
  dplyr::summarise(
    breaks_max=max(breaks), breaks_max_norm=breaks_max/gene_length[1], 
    breaks_sum=sum(breaks), breaks_sum_norm=breaks_sum/gene_length[1], 
    breaksite_dist=pmin(abs(gene_start[1]-breaksite_pos[1]), abs(gene_end[1]-breaksite_pos[1])),
    breaks_few=breaks_max<10) %>%
  dplyr::group_by(breaksite_full, breaks_few, bin_strand, gene_strand) %>%
  dplyr::do((function(z) {
    zz<<-z
    if(z$breaks_few[1]) {
      z = z[sample(1:nrow(z), pmin(100, nrow(z))),]
    }
    z
  })(.)) %>%
  dplyr::mutate(strand_full=paste0("g", gene_strand, "/b", bin_strand)) %>%
  dplyr::ungroup() %>%
  dplyr::select(breaksite_full, gene_name, gene_length, breaks_max, breaksite_dist, strand_full)

g = ggplot(results.dist_sum, aes(x=log10(breaksite_dist), y=log10(breaks_max), color=strand_full)) +
  geom_point(aes(size=gene_length, gene_name=gene_name), alpha=0.3) +
  # geom_hline(yintercept=th.breaks_many) + 
  # ggrepel::geom_text_repel(aes(label=gene_name), data=. %>% dplyr::filter(breaks>th.breaks_many)) +
  geom_smooth() +
  facet_wrap(~breaksite_full, scale="free_y")
ggplotly(g)
p = plotly::ggplotly(g)
htmlwidgets::saveWidget(plotly::as_widget(p), "index.html")
# plotly::plotly_POST(g, "breaks2breaksite_dist", sharing="public")

#
# Calculate cross-correlation
#
results_ccf = results %>%
  dplyr::filter(breaks_max>=th.min_bin_breaks) %>% 
  dplyr::arrange(breaksite_chrom, gene_full) %>%
  dplyr::group_by(breaksite_chrom, gene_full, gene_name, gene_chrom, gene_start, gene_end, gene_region_start, gene_region_end) %>%
  dplyr::do((function(z){
    zz<<-z
    z.pos = z[z$gene_strand=="+",]
    z.neg = z[z$gene_strand=="-",]
    
    ts.positive = ts(z.pos$breaks)
    ts.negative = ts(z.neg$breaks)
    ts.negative_rev = ts(rev(z.neg$breaks))
    
    ccf.sense = ccf(ts.positive, ts.negative, plot=F)
    if(all(is.na(ccf.sense$acf))) {
      acf.sense = NA
      lag.sense = NA
      pvalue.sense = NA
    } else {
      f = which.max(abs(ccf.sense$acf))
      acf.sense = ccf.sense$acf[f]
      lag.sense = ccf.sense$lag[f]/nrow(z.pos)
      pvalue.sense = cor.test(z.pos$breaks, z.neg$breaks)$p.value
    }
    
    
    ccf.antisense = ccf(ts.positive, ts.negative_rev, plot=F)
    if(all(is.na(ccf.antisense$acf))) {
      acf.antisense = NA
      lag.antisense = NA
      pvalue.antisense = NA
    } else {
      f = which.max(abs(ccf.antisense$acf))
      acf.antisense = ccf.antisense$acf[f]
      lag.antisense = ccf.antisense$lag[f]/nrow(z.pos)
      pvalue.antisense = cor.test(z.pos$breaks, rev(z.neg$breaks))$p.value
    }
    
    data.frame(sense=c("sense", "antisense"), lag=c(lag.sense, lag.antisense), cor=c(acf.sense, acf.antisense), pvalue=c(pvalue.sense, pvalue.antisense))
  })(.)) %>%
  dplyr::group_by(gene_full) %>%
  dplyr::mutate(is_hit=pvalue <= th.ccf_pvalue  & abs(lag) < th.ccf_maxlag & cor>=th.ccf_correlation) %>%
  dplyr::mutate(is_hit=is_hit & !all(is_hit)) %>%
  dplyr::ungroup() 

# results_hits = results_ccf %>%
#   dplyr::inner_join(results %>% dplyr::distinct(gene_full, breaksite_chrom, gene_length), by=c("gene_full", "breaksite_chrom"))  %>% 
#   dplyr::filter(is_hit & sense=="antisense" & gene_length>600000) %>%
#   dplyr::select(-gene_length)
# 
# 
# results_hits = results_ccf %>% 
#   dplyr::filter(sense=="antisense" & grepl(":(Dgkb|Csmd3|Plcb1|Inpp4b|Lrrtm4|Grik2|Ctnna3|Ctnnd2|Csmd2|Pard3b|Creb5|Ccser1|Dgki|Exoc6b|Magi1|Erc1|Grin2b|Prickle2|Ctnna2|Grid2|Cacna1c)$", gene_full)) %>%
#   dplyr::mutate(is_hit=grepl(":(Ccser1|Dgkb|Ctnna3|Inpp4b)$", gene_full))


re_genes = ":(Npas3|Csmd3| Lsamp|Prkg1||Ctnna2|Ptn|Pard3b|Maml2)$"
re_genes = ":(Npas3|Lsamp|Nrxn1|Ptn|Nfia|Ctnna2|Sdk1|Grid2|Csmd1|Pard3b|Prkg1|Maml2|Csmd3|Tcf4|Lrp1b|Cdh13|Grik2|Nrxn3|Gpc6|Ctnnd2|Rbfox1|Cadm2|Nlgn1|Auts2)$"
results_hits = results_ccf %>% 
  dplyr::filter(sense=="antisense" & grepl(re_genes, gene_full))


pdf("results.pdf", paper="a4r")
for(chr in unique(results_hits$breaksite_chrom)) {
  results.chr = results %>%
    dplyr::inner_join(results_hits %>% dplyr::select(breaksite_chrom, gene_full, gene_start, gene_end, gene_region_start, gene_region_end, lag, cor, is_hit), by=c("gene_full", "breaksite_chrom", "gene_start", "gene_end", "gene_region_start", "gene_region_end")) %>%
    dplyr::filter(breaksite_chrom == gene_chrom | breaks_max>=th.min_bin_breaks) %>%
    dplyr::mutate(facet=paste0(gene_full, "(", round(gene_start/1e6), "Mb + ", round(gene_length/1e6,2), "Mb)")) %>%
    dplyr::mutate(facet_cor=paste0(" l:", round(lag*100), "%, R:.", round(cor*100,0)))
  # results.chr = results.chr %>% dplyr::filter(breaksite_chrom==chr)
  
  ggplot(results.chr) +
    geom_hline(aes(yintercept=0), color="#999999") +
    # geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf), alpha=0.05, fill="#FF0000", data=. %>% dplyr::filter(is_hit) %>% dplyr::distinct(facet)) +
    geom_vline(aes(xintercept=gene_start), data=. %>% dplyr::distinct(facet, .keep_all=T)) +
    geom_vline(aes(xintercept=gene_end), data=. %>% dplyr::distinct(facet, .keep_all=T)) +
    geom_line(aes(x=bin_start+(bin_end-bin_start)/2, y=breaks_sign, group=paste(gene_strand, breaksite_chrom), color=breaksite_chrom, size=ifelse(breaksite_chrom==gene_chrom, "yes", "no"), alpha=0.3)) +
    geom_line(aes(x=bin_start+(bin_end-bin_start)/2, y=breaks_sign, group=paste(gene_strand, breaksite_chrom), color=breaksite_chrom, size=ifelse(breaksite_chrom==gene_chrom, "yes", "no"), alpha=0.3)) +
    geom_blank(aes(y=breaks)) + 
    geom_blank(aes(y=-breaks)) + 
    facet_wrap(~facet, scales="free_x", ncol=7) + 
    scale_size_manual(values=c("yes"=2, "no"=0.5)) + 
    scale_x_continuous(labels=scales::label_number_si(accuracy=0.2)) + 
    labs(size="In induced chromosome")
}
dev.off()
plotly::ggplotly(p)




# for(strand in unique(results$gene_strand_name)) {
#   r.con = file("results.wig", "w")
#   wig = results %>%
#     dplyr::filter(breaksite_chrom==gene_chrom & gene_strand_name==strand) %>%
#     dplyr::arrange(gene_chrom, gene_name, bin_start, bin_end)
#   wig.regions = wig %>% dplyr::distinct(gene_chrom, gene_name, gene_strand_name)
#   for(r in 1:nrow(wig.regions)) {
#     r.chrom = wig.regions$gene_chrom[r]
#     r.strand = wig.regions$gene_strand_name[r]
#     r.gene = wig.regions$gene_name[r]
#     r.wig = wig %>% dplyr::filter(gene_chrom==r.chrom & gene_name==r.gene & gene_strand_name==r.strand)
#     if(nrow(r.wig)<5) next
#     r.step = unique(diff(r.wig$bin_start)); stopifnot(length(r.step) == 1)
#     
#     r.wig = r.wig %>% 
#       dplyr::mutate(which_first=ifelse(any(breaks_sign!=0), min(which(breaks_sign!=0)), Inf)) %>%
#       dplyr::mutate(keep=1:n() >= which_first) %>%
#       dplyr::filter(keep)
#     if(nrow(r.wig)==0) next
#     
#     writeLines(paste0("fixedStep chrom=", r.chrom, " start=", min(r.wig$bin_start)," step=", r.step, " span=", r.step), con=r.con)
#     writeLines(apply(r.wig %>% dplyr::select(bin_start, breaks), 1, paste, collapse="\t"), con=r.con)
#   }   
#   close(r.con)
# }

#fixedStep chrom=chr3 start=400601 step=100 span=5  