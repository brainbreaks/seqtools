
# Grid2 - has only few hypomethylation sites
# Erc1 - has only few hypomethylation sites

library(dplyr)
library(ggplot2)
library(scales)
library(reshape2)
library(ggrepel)


th.ccf_pvalue = 0.05
th.ccf_minlag = 0.1
th.ccf_maxlag = 0.1
th.ccf_correlation = 0.5
th.min_bin_breaks = 5

results =  data.frame()
for(results_path in list.files("results", full.names=T)) {
  chrom_induction = gsub(".(tsv|bed)$", "", basename(results_path))
  results.chr = readr::read_tsv(results_path)  %>% 
    dplyr::mutate(chrom_induction=chrom_induction, gene_full=paste0(gene_chrom, ":", gene_name), gene_pos=(gene_end+gene_start)/2) %>%
    tidyr::extract(attributes, into=c("score"), regex="score=([^;]+)") %>%
    dplyr::mutate(score=as.numeric(score), score_sign=ifelse(gene_strand=="+", score, -score)) %>%
    dplyr::mutate(gene_strand_name=ifelse(gene_strand=="+", "positive", "negative")) %>%
    dplyr::filter(gene_name!=".") %>%
    dplyr::group_by(gene_full) %>% 
    dplyr::mutate(score_max=max(score), gene_length=max(gene_end)-min(gene_start)) %>%
    dplyr::filter(any(score_max>th.min_bin_breaks) & n()>1) %>% 
    dplyr::group_by(gene_full, gene_strand_name) %>%
    dplyr::mutate(peak_relative_loc=which.max(score)/n(), score_var=var(score)) %>%
    dplyr::ungroup()
  results = dplyr::bind_rows(results, results.chr)
}

results %>%
  reshape2::dcast(gene_full + gene_strand + bin_start + bin_end ~ chrom_induction, value.var="score") %>%
  ggplot() +
  geom_point(aes(x=chr8, y=chr7)) + 
  coord_cartesian(xlim=c(1,20), ylim=c(1,20))


results_ccf = results %>%
  dplyr::arrange(chrom_induction, gene_full, gene_pos) %>%
  dplyr::group_by(chrom_induction, gene_full, gene_name, gene_length) %>%
  dplyr::do((function(z){
    zz<<-z
    # asda()
    # if(z$gene_name=="Setd5") {
    #   zz<<-z
    #   asdas()
    # }
    z.pos = z[z$gene_strand=="+",]
    z.neg = z[z$gene_strand=="-",]
    
    ts.positive = ts(z.pos$score)
    ts.negative = ts(z.neg$score)
    ts.negative_rev = ts(rev(z.neg$score))
    
    ccf.sense = ccf(ts.positive, ts.negative, plot=F)
    if(all(is.na(ccf.sense$acf))) {
      acf.sense = NA
      lag.sense = NA
      pvalue.sense = NA
    } else {
      f = which.max(abs(ccf.sense$acf))
      acf.sense = ccf.sense$acf[f]
      lag.sense = ccf.sense$lag[f]/nrow(z.pos)
      pvalue.sense = cor.test(z.pos$score, z.neg$score)$p.value
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
      pvalue.antisense = cor.test(z.pos$score, rev(z.neg$score))$p.value
    }
    
    data.frame(sense=c("sense", "antisense"), lag=c(lag.sense, lag.antisense), cor=c(acf.sense, acf.antisense), pvalue=c(pvalue.sense, pvalue.antisense))
  })(.)) %>%
  dplyr::group_by(gene_full) %>%
  dplyr::mutate(is_hit=pvalue <= th.ccf_pvalue  & abs(lag) < th.ccf_maxlag & cor>=th.ccf_correlation) %>%
  dplyr::mutate(is_hit=is_hit & !all(is_hit)) %>%
  dplyr::ungroup() 
results_hits = results_ccf %>% 
  dplyr::filter(is_hit & sense=="antisense")

results_hits = results_ccf %>% 
  dplyr::filter(sense=="antisense" & gene_name %in% c("Creb5", "Ccser1", "Dgki", "Exoc6b", "Magi1", "Erc1", "Grin2b", "Prickle2", "Ctnna2", "Grid2", "Cacna1c"))%>%
  dplyr::mutate(is_hit=gene_name %in% c("Ccser1", "Ctnna2", "Erc1", "Grid2"))

p = results %>% 
  dplyr::inner_join(results_hits %>% dplyr::select(-gene_length), by=c("gene_full", "gene_name")) %>%
  dplyr::mutate(facet=paste0(gene_full, "(", round(gene_start/1e6), "Mb + ", round(gene_length/1e6,2), "Mb) ", " l:", round(lag*100), "%, R:.", round(cor*100,0))) %>%
  ggplot() +
    geom_hline(aes(yintercept=0), color="#999999") +
    geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf), alpha=0.05, fill="#FF0000", data=. %>% dplyr::filter(is_hit) %>% dplyr::distinct(facet)) +
    geom_line(aes(x=gene_pos, y=score_sign, color=gene_strand)) +
    geom_blank(aes(y=score)) + 
    geom_blank(aes(y=-score)) + 
    facet_wrap(~facet, scales="free") + 
    scale_x_continuous(labels=scales::label_number_si(accuracy=0.2))
p
plotly::ggplotly(p)
  