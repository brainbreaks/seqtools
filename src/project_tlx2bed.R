library(readr)
library(dplyr)
library(IRanges)
library(tidyr)
library(ggplot2)
library(BSgenome)
library(baseline)

seqlengths <- c(chr1=60, chr2=20, chr3=25)



scale_breaks = function(x) {
    breaks = seq(0, max(x), 1e7)
    names(breaks) = paste0(breaks/1e6, "Mb")
    breaks
}

calculate_coverage2 = function(ranges, mm9, extend=5e6, binsize=1000, binstep=100) {
  mm9_tiles = GenomicRanges::trim(GenomicRanges::resize(unlist(GenomicRanges::tileGenome(mm9, tilewidth=binstep)), binsize, fix="center"))

  seqinfo(ranges) = mm9
  ranges.extended = GenomicRanges::trim(GenomicRanges::resize(ranges, width=extend, fix="center"))
  coverage_df = as.data.frame(ranges.extended) %>%
    dplyr::group_by(break_bait_chrom) %>%
    dplyr::do((function(z){
      zz<<-z
      z_ranges = GenomicRanges::makeGRangesFromDataFrame(z, keep.extra.columns=T)
      z_tiles_df = as.data.frame(mm9_tiles) %>%
        dplyr::mutate(coverage=GenomicRanges::countOverlaps(mm9_tiles, z_ranges)) %>%
        tidyr::crossing(z %>% dplyr::select(break_bait_chrom, scale_factor) %>% dplyr::slice(1)) %>%
        dplyr::select(seqnames, start, end, coverage, scale_factor)
      z_tiles_df
    })(.)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(coverage_norm=coverage*scale_factor) %>%
    dplyr::select(-scale_factor) %>%
    dplyr::group_by(seqnames, start, end) %>%
    dplyr::summarise(coverage=sum(coverage_norm))
  coverage_df
}

calculate_coverage = function(ranges, mm9, extend) {
  ranges.extended = GenomicRanges::resize(ranges, width=extend, fix="center")
  ranges.extended = GenomicRanges::restrict(ranges.extended, start=0, end=setNames(seqlengths(mm9), seqnames(mm9)))
  coverage_df = as.data.frame(ranges.extended) %>%
    dplyr::group_by(break_bait_chrom) %>%
    dplyr::do((function(z){
      zz<<-z
      z_ranges = GenomicRanges::makeGRangesFromDataFrame(z, keep.extra.columns=T)
      as.data.frame(as(GenomicRanges::coverage(z_ranges), "GRanges")) %>%
        tidyr::crossing(z %>% dplyr::select(break_bait_chrom, scale_factor) %>% dplyr::slice(1)) %>%
        dplyr::select(seqnames, start, end, score, scale_factor)
    })(.)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(score_norm=score*scale_factor) %>%
    dplyr::select(-scale_factor)

  # tempfile()
  total_path  = "data/macs2/total.bdg"
  chr_tmp_path  = "data/macs2/chr.bdg"
  # chr_tmp_path = tempfile()
  # total_path = tempfile()
  i = 0
  for(chr in unique(coverage_df$break_bait_chrom)) {
    # chr = "chr3"
    coverage_df.chr = coverage_df %>% dplyr::filter(break_bait_chrom==chr)

    if(i==0) {
      print("Write total")
      readr::write_tsv(coverage_df.chr %>% dplyr::select(seqnames, start, end, score_norm), total_path, col_names=F)
    } else {
      print(paste0("Write chr: ", chr))
      readr::write_tsv(coverage_df.chr %>% dplyr::mutate(score_norm=-score_norm) %>% dplyr::select(seqnames, start, end, score_norm), chr_tmp_path, col_names=F)
      system(stringr::str_glue("macs2 bdgcmp -m subtract -t {total} -c {chr} -o {total}", total=total_path, chr=chr_tmp_path))
    }

    i = i + 1
  }

  bedgraph_cols = cols(seqnames=col_character(), start=col_double(), end=col_double(), coverage=col_double())
  coverage_norm_df = readr::read_tsv(total_path, col_names=names(bedgraph_cols$cols), col_types=bedgraph_cols) %>%
    dplyr::mutate(coverage_id=1:n())

  coverage_norm_df
}

call_peaks = function(sample_df, control_df=NULL, debug=F) {
  binsize = 1e3
  binstep = 5e2
  extend = 1e5
  llocal = 6e6
  peaks_minqvalue=-log10(0.01)
  peaks_maxgap=1e5
  peaks_minlen=200
  debug=F

  do_compare = !is.null(control_df)

  sample_bed_path = "data/breakensembl/Sample.bed"
  control_bed_path = "data/breakensembl/Control.bed"
  baseline_bdg_path = "data/macs2/baseline.bdg"
  sample_bdg_path = "data/macs2/Sample.bdg"
  control_bdg_path = "data/macs2/Control.bdg"
  qvalue_path = "data/macs2/qvalue.bdg"
  peaks_path = "data/macs2/peaks.bed"


  sample_df = sample_df %>% dplyr::mutate(break_name=stringr::str_glue("{bait}_{name}", bait=break_bait_chrom, name=break_name))
  if(debug) {
    readr::write_tsv(sample_df %>% dplyr::select(break_chrom, break_start, break_end, break_name, break_score, break_strand), sample_bed_path, col_names=F)
  }

  if(!do_compare) {
    control_df = sample_df
    control_bed_path = sample_bed_path
  } else {
    control_df = control_df %>% dplyr::mutate(break_name=stringr::str_glue("{bait}_{name}", bait=break_bait_chrom, name=break_name))
    if(debug) {
      readr::write_tsv(control_df %>% dplyr::select(break_chrom, break_start, break_end, break_name, break_score, break_strand), control_bed_path, col_names=F)
    }
  }

  # slocal background
  options("UCSC.goldenPath.url"="https://hgdownload.cse.ucsc.edu/goldenPath")
  mm9 = GenomeInfoDb::Seqinfo(genome="mm9")[paste0("chr", 1:19)]
  mm9_size = sum(mm9@seqlengths)
  mm9_effective_size=0.74 * mm9_size
  mm9_tiles = unlist(tileGenome(mm9, tilewidth=binsize))
  mm9_tiles$tile_id = 1:length(mm9_tiles)

  control_ranges = GenomicRanges::makeGRangesFromDataFrame(control_df, end.field="break_end", start.field="break_start", strand.field="break_strand", seqnames.field="break_chrom", keep.extra.columns=T)
  control_tiles_df.coverage = calculate_coverage2(control_ranges, mm9, extend=extend, binsize=binsize, binstep=binstep)
  readr::write_tsv(control_tiles_df.coverage %>% dplyr::select(seqnames, start, end, coverage), control_bdg_path, col_names=F)

  #control_df.coverage = calculate_coverage(sample_ranges, mm9, extend)
  #readr::write_tsv(control_df.coverage %>% dplyr::select(seqnames, start, end, coverage), control_bdg_path, col_names=F)
  #control_ranges.coverage = GenomicRanges::makeGRangesFromDataFrame(control_df.coverage, keep.extra.columns=T)
  #coverage2tile.map = as.data.frame(GenomicRanges::findOverlaps(mm9_tiles, control_ranges.coverage, ignore.strand=TRUE))
  #control_tiles_df.coverage = as.data.frame(mm9_tiles) %>%
  #  dplyr::left_join(coverage2tile.map, by=c("tile_id"="queryHits")) %>%
  #  dplyr::left_join(control_df.coverage %>% dplyr::select(coverage_id, coverage), by=c("subjectHits"="coverage_id")) %>%
  #  dplyr::group_by(seqnames, start, end) %>%
  #  dplyr::summarise(coverage=tidyr::replace_na(max(coverage), 0))
  #readr::write_tsv(control_tiles_df.coverage %>% dplyr::select(seqnames, start, end, coverage), "data/macs2/Control_tile.bdg", col_names=F)

  sample_ranges = GenomicRanges::makeGRangesFromDataFrame(sample_df, end.field="break_end", start.field="break_start", strand.field="break_strand", seqnames.field="break_chrom", keep.extra.columns=T)
  sample_df.coverage = calculate_coverage(sample_ranges, mm9, extend)
  readr::write_tsv(sample_df.coverage %>% dplyr::select(seqnames, start, end, coverage), sample_bdg_path, col_names=F)

  fit_baseline = function(z) {
    zz<<-z
    x = matrix(as.numeric(z$coverage_mod), nrow=1)
    colnames(x) = z$start
    bc.irls = baseline::baseline(x, method="medianWindow", hws=5e6/binstep, hwm=5e6/binstep) # !!!!
    z.coverage_smooth = as.numeric(bc.irls@baseline)
    # bc.irls = baseline::baseline(x, method="medianWindow", hws=5e6/binstep, hwm=4e6/binstep) # !!!!
    # z.coverage_smooth2 = as.numeric(bc.irls@baseline)

    cbind(z[c("seqnames", "start", "end", "coverage", "coverage_mod")], coverage_smooth=z.coverage_smooth)
  }

  # cluster = multidplyr::new_cluster(20)
  # multidplyr::cluster_assign(cluster, fit_baseline=fit_baseline)
  control_df.baseline = control_tiles_df.coverage %>%
    dplyr::mutate(break_exp_condition="Control") %>%
    dplyr::arrange(seqnames, start) %>%
    dplyr::group_by(break_exp_condition, seqnames) %>%
    dplyr::mutate(is_q2=coverage>quantile(coverage, 0.2), is_q8=coverage<quantile(coverage, 0.8), is_inner_quantile=is_q2 & is_q8) %>%
    dplyr::mutate(coverage_median=median(coverage[is_inner_quantile], na.rm=T), coverage_sd=sd(coverage[is_inner_quantile], na.rm=T)) %>%
    dplyr::mutate(coverage_min=coverage_median-2*coverage_sd) %>%
    dplyr::mutate(coverage_mod=ifelse(coverage<coverage_min, NA_real_, coverage)) %>%
    dplyr::mutate(is_filled=is.na(coverage_mod), coverage_mod=zoo::na.fill(coverage_mod, "extend")) %>%
    dplyr::do(fit_baseline(.)) %>%
    dplyr::ungroup()
  class(control_tiles_df.coverage) = c("tbl_df", "tbl", "data.frame")
  readr::write_tsv(control_df.baseline %>% dplyr::select(seqnames, start, end, coverage_smooth), baseline_bdg_path, col_names=F)

  control_df.significant = control_df.baseline %>%
    dplyr::rowwise() %>%
    dplyr::mutate(pvalue=pgamma(coverage, coverage_smooth, 1, lower.tail=F))

  # Call peaks
  readr::write_tsv(control_df.baseline %>% dplyr::select(seqnames, start, end, coverage), sample_bdg_path, col_names=F)
  readr::write_tsv(control_df.baseline %>% dplyr::select(seqnames, start, end, coverage_smooth), baseline_bdg_path, col_names=F)
  system(stringr::str_glue("macs2 bdgcmp -t {input} -c {noise} -m qpois -o {output}", input=sample_bdg_path, noise=baseline_bdg_path, output=qvalue_path))
  system(stringr::str_glue("macs2 bdgpeakcall -i {qvalue} -c {cutoff} --min-length {format(peak_min_length, scientific=F)} --max-gap {format(peak_max_gap, scientific=F)} -o {output}",
       qvalue=qvalue_path, output=peaks_path, cutoff=peaks_minqvalue, peak_max_gap=peaks_maxgap, peak_min_length=peaks_minlen))

  # Read results file
  peaks_cols = cols(
    peak_chrom=col_character(), peak_start=col_double(), peak_end=col_double(), peak_length=col_character(), peak_abs_summit=col_double(),
    peak_score=col_character(), peak_fc=col_double(), peak_pvalue_log10=col_double(), peak_qvalue_log10=col_double(), peak_sammit_offset=col_double()
  )
  qvalue_cols = cols(qvalue_chrom=col_character(), qvalue_start=col_double(), qvalue_end=col_double(), qvalue_score=col_double())
  qvalues_df = readr::read_tsv(qvalue_path, col_names=names(qvalue_cols$cols), col_types=qvalue_cols) %>%
    dplyr::mutate(qvalue_id=1:n())
  peaks_df = readr::read_tsv(peaks_path, skip=1, col_names=names(peaks_cols$cols), col_types=peaks_cols) %>%
    dplyr::mutate(peak_summit_pos=peak_start + peak_sammit_offset) %>%
    dplyr::mutate(peak_id=1:n())


  if(debug) {
    pdf(file="reports/baseline_5e6_6.pdf", width=15, height=4)
      control_tiles_df.coverage_sample = control_tiles_df.coverage %>% dplyr::sample_frac(0.2)
      for(chr in paste0("chr", 1:19)) {
        print(chr)
        ggplot.colors = c("pileup"="#333333", "smooth"="#FF0000", "smooth2"="#00AA00", mypeak="#CCCCCC", macspeaks="#FF0000")
        p = ggplot() +
          # geom_hex(aes(y=coverage, x=start), data=control_tiles_df.coverage %>% dplyr::filter(dplyr::between(coverage, 2, 50) & seqnames %in% chr), bins=30) +
          geom_point(aes(y=coverage, x=start), data=control_tiles_df.coverage_sample %>% dplyr::filter(dplyr::between(coverage, 0, 50) & seqnames %in% chr), alpha=0.05, size=0.5) +
          geom_rect(aes(xmin=peak_start, ymin=-2, xmax=peak_end, ymax=-1, color="macspeaks"), data=peaks_df %>% dplyr::mutate(seqnames=peak_chrom) %>% dplyr::filter(seqnames %in% chr)) +
          geom_rect(aes(xmin=start, ymin=-3, xmax=end, ymax=-2, color="mypeak"), data=control_df.significant %>% dplyr::filter(seqnames %in% chr & pvalue < 0.01)) +
          geom_line(aes(y=coverage_smooth, x=start, color="smooth"), data=control_df.baseline %>% dplyr::filter(seqnames %in% chr), alpha=0.8) +
          # geom_line(aes(y=coverage_smooth2, x=start, color="smooth2"), data=control_df.baseline %>% dplyr::filter(seqnames %in% chr), alpha=0.8) +
          coord_cartesian(ylim=c(-5, 50)) +
          scale_color_manual(values=ggplot.colors) +
          scale_fill_manual(values=ggplot.colors) +
          scale_x_continuous(breaks=scale_breaks) +
          facet_wrap(~seqnames, scales="free_x", ncol=1)
        print(p)
      }
    dev.off()
  }






  cx.raw = control_df.baseline %>%
    dplyr::filter(!is_filled & coverage_smooth<1000) %>%
    dplyr::mutate(coverage_lo=floor(coverage_smooth), coverage_hi=ceiling(coverage_smooth)) %>%
    dplyr::group_by(coverage_lo, coverage_hi) %>%
    dplyr::filter(n()>1000) %>%
    dplyr::filter(coverage<=quantile(coverage, 0.9)) %>%
    dplyr::ungroup()

  x = x.raw %>%
    dplyr::group_by(coverage_start, coverage_end) %>%
    dplyr::do((function(z){
      zz<<-z
      my.mle = MASS::fitdistr(z$coverage, densfun="gamma")
      dplyr::bind_cols(z[1, c("coverage_start", "coverage_end")], as.data.frame(t(my.mle$estimate)), n=nrow(z))
    })(.))

  pdf(file="reports/coverage_distribution_across_baseline.pdf", width=10, height=20)
  ggplot(x.raw) +
    ggridges::geom_density_ridges(aes(x=coverage, y=factor(coverage_lo), fill="raw"))

  ggplot(x.raw) +
    geom_density(aes(x=coverage, y=..count..)) +
    facet_wrap(~factor(coverage_lo), ncol=1)
  dev.off()

  ggplot(x %>% dplyr::mutate(value=shape/rate)) +
    geom_line(aes(x=coverage_start, y=value))


  ggplot(x %>% reshape2::melt(measure.vars=c("shape", "rate"))) +
    geom_line(aes(x=coverage_start, y=value)) +
    facet_wrap(~variable, scale="free", ncol=1)

  x = control_df.baseline %>%
    dplyr::filter(!is_filled & coverage_smooth<100) %>%
    dplyr::mutate(coverage_smooth_cut=cut(coverage_smooth, 100)) %>%
    dplyr::group_by(coverage_smooth_cut) %>%
    dplyr::summarise(q=quantile(coverage_smooth - coverage, 0.9), coverage_smooth_min=min(coverage_smooth)) %>%
    dplyr::filter(q>0)
  ggplot(x) +
    geom_smooth(aes(coverage_smooth_min, q)) +
    coord_cartesian(ylim=c(0, 10))

  ggplot(control_df.baseline %>% dplyr::filter(!is_filled & coverage_smooth<100 & abs(coverage_smooth - coverage)<100)) +
    geom_density_2d_filled(aes(coverage_smooth, coverage_smooth - coverage), bins=100)

  plot(control_df.baseline$coverage, control_df.baseline$coverage_smooth - control_df.baseline$coverage)
  control_df.baseline$coverage_smooth - control_df.baseline$coverage
  control_df.baseline$score
  #
  # Y1 <- 100
  # set.seed(101)
  # x.pois<-rpois(Y1, 20)
  # hist(x.pois, breaks=100,freq=FALSE)
  # lines(density(x.pois, bw=0.8), col="red")
  # library(MASS)
  # (my.mle<-fitdistr(x.pois, densfun="poisson"))
  # ##      lambda
  # ##   20.6700000
  # ##  ( 0.4546427)
  # BIC(my.mle)


  pdf(file="reports/control_estimated_noise.pdf", width=20, height=20)
  # my.mle = MASS::fitdistr(x$score, densfun="gamma")

  # x.poi<-rpois(n=1e6,lambda=7.433379128)
  # x.poi<-rgamma(n=1e5, shape=3.258554685, rate=0.438367973)
  # plot(density(x.poi, bw=1), xlim=c(0,70))

  x = control_df.baseline %>%
    dplyr::mutate(residual=coverage-coverage_smooth) %>%
    dplyr::filter(coverage>0) %>%
    dplyr:::filter(dplyr::between(residual, quantile(residual, 0.01), quantile(residual, 0.99)))

  ggplot(x) +
    geom_density(aes(x=coverage))

  MASS::fitdistr(x$coverage, densfun="gamma")
  pgamma(25, 3.430902785, 0.433324396, lower.tail=F)

  shape = 3.258554685
  rate = 0.438367973
  mean = mean(x$coverage)
  sd = sd(x$coverage)
  shape1 = mean^2/sd^2
  rate1 = mean/sd^2
  shape1 = 3.0
  rate1 = 0.441
  ggplot(x) +
    geom_density(aes(x=coverage)) +
    geom_line(aes(x=coverage, y=density, color="gamma"), data=with(density(rgamma(n=1e5, shape=shape, rate=rate), from=min(x$coverage, na.rm=T), to=max(x$coverage, na.rm=T)), data.frame(coverage=x, density=y))) +
    geom_line(aes(x=coverage, y=density, color="gamma1"), data=with(density(rgamma(n=1e5, shape=shape1, rate=rate1), from=min(x$coverage, na.rm=T), to=max(x$coverage, na.rm=T)), data.frame(coverage=x, density=y)))

  ggplot(x) +
    geom_density(aes(x=residual)) +
    coord_cartesian(xlim=c(-23, 23))

  dev.off()

  if(F) {
    pdf(file="reports/raw_control_vs_sample.pdf", width=10, height=80)
      chr = paste0("chr", 1:19)
      d = dplyr::bind_rows(control_df.coverage %>% dplyr::mutate(break_exp_condition="Control"), sample_df.coverage %>% dplyr::mutate(break_exp_condition="Sample"))
      p = ggplot(d %>% dplyr::filter(seqnames %in% chr)) +
        ggridges::geom_ridgeline(aes(y=break_exp_condition, x=start, height=coverage, fill=break_exp_condition), scale=0.1) +
        ggridges::geom_ridgeline(aes(y=break_exp_condition, x=start, height=coverage), scale=0.1, data=control_df.baseline %>% dplyr::filter(seqnames %in% chr), color="#FF0000") +
        coord_cartesian(ylim=c(0, 50)) +
        scale_x_continuous(breaks=scale_breaks) +
        facet_wrap(~seqnames, scales="free_x", ncol=1)
      print(p)
    dev.off()
  }


  # Call peaks
  system(stringr::str_glue("macs2 bdgcmp -t {input} -c {noise} -m qpois -o {output}", input=sample_bdg_path, noise=baseline_bdg_path, output=qvalue_path))
  system(stringr::str_glue("macs2 bdgpeakcall -i {qvalue} -c {cutoff} --min-length {format(peak_min_length, scientific=F)} --max-gap {format(peak_max_gap, scientific=F)} -o {output}",
       qvalue=qvalue_path, output=peaks_path, cutoff=peaks_minqvalue, peak_max_gap=peaks_maxgap, peak_min_length=peaks_minlen))

  # Read results file
  peaks_cols = cols(
    peak_chrom=col_character(), peak_start=col_double(), peak_end=col_double(), peak_length=col_character(), peak_abs_summit=col_double(),
    peak_score=col_character(), peak_fc=col_double(), peak_pvalue_log10=col_double(), peak_qvalue_log10=col_double(), peak_sammit_offset=col_double()
  )
  qvalue_cols = cols(qvalue_chrom=col_character(), qvalue_start=col_double(), qvalue_end=col_double(), qvalue_score=col_double())
  qvalues_df = readr::read_tsv(qvalue_path, col_names=names(qvalue_cols$cols), col_types=qvalue_cols) %>%
    dplyr::mutate(qvalue_id=1:n())

  peaks_df = readr::read_tsv(peaks_path, skip=1, col_names=names(peaks_cols$cols), col_types=peaks_cols) %>%
    dplyr::mutate(peak_summit_pos=peak_start + peak_sammit_offset) %>%
    dplyr::mutate(peak_id=1:n())

  qvalues_ranges = GenomicRanges::makeGRangesFromDataFrame(qvalues_df, end.field="qvalue_end", start.field="qvalue_start", seqnames.field="qvalue_chrom", keep.extra.columns=T)
  summit_ranges = GenomicRanges::makeGRangesFromDataFrame(peaks_df, end.field="peak_summit_pos", start.field="peak_summit_pos", seqnames.field="peak_chrom", keep.extra.columns=T)
  summit2qvalue.map = as.data.frame(IRanges::mergeByOverlaps(qvalues_ranges, summit_ranges)) %>%
    dplyr::filter(qvalues_ranges.seqnames==summit_ranges.seqnames) %>%
    dplyr::select(peak_id, peak_summit_qvalue=qvalues_ranges.qvalue_score)
  peaks_df = peaks_df %>%
    dplyr::select(-dplyr::matches("peak_summit_qvalue")) %>%
    dplyr::inner_join(summit2qvalue.map, by=c("peak_id"))

  qvalues %>%
    dplyr::inner_join()
  peaks$peak_start + peaks$peak_sammit_offset

  # Call peaks
  macs2_peaks_path = "data/macs2/Sample_peaks.bed"
  system(stringr::str_glue("macs2 bdgpeakcall -i {qvalue} -c {cutoff} --min-length {format(peak_min_length, scientific=F)} --max-gap {format(peak_max_gap, scientific=F)} -o {output}", qvalue=macs2_qvalue_path, output=macs2_peaks_path, cutoff=-log10(macs2_params$qvalue_cutoff), peak_max_gap=macs2_params$peak_max_gap, peak_min_length=macs2_params$peak_min_length))

  export(bins %>% dplyr::select(seqnames, start, end, coverage_smooth), format="bedGraph")


  macs2_params$extsize / macs2_params$slocal
  readr::write_tsv(bins, file="data/macs2/coverage.bdg", col_names=F)


  # llocal background
  macs2_llocal_path = "data/macs2/Control_llocal.bdg"
  system(stringr::str_glue("macs2 pileup -i {input} -f BED -B --extsize {format(extsize, scientific=F)} -o {output}", input=control_bed_path, output=macs2_llocal_path, extsize=extend/2))
  macs2_llocal_norm_path = "data/macs2/Control_llocal_norm.bdg"
  system(stringr::str_glue("macs2 bdgopt -i {input} -m multiply -p {norm} -o {output}", input=macs2_llocal_path, output=macs2_llocal_norm_path, norm=macs2_params$extsize / macs2_params$llocal))



  # llocal background
  macs2_llocal_path = "data/macs2/Control_llocal.bdg"
  system(stringr::str_glue("macs2 pileup -i {input} -f BED -B --extsize {format(extsize, scientific=F)} -o {output}", input=control_bed_path, output=macs2_llocal_path, extsize=macs2_params$llocal/2))
  macs2_llocal_norm_path = "data/macs2/Control_llocal_norm.bdg"
  system(stringr::str_glue("macs2 bdgopt -i {input} -m multiply -p {norm} -o {output}", input=macs2_llocal_path, output=macs2_llocal_norm_path, norm=macs2_params$extsize / macs2_params$llocal))

  #T

  # glocal background
  #BiocManager::install("BSgenome.Mmusculus.UCSC.mm9")
  library(BSgenome)
  options("UCSC.goldenPath.url"="https://hgdownload.cse.ucsc.edu/goldenPath")
  mm9 = GenomeInfoDb::Seqinfo(genome="mm9")
  mm9_size = sum(mm9@seqlengths)
  mm9_effective_size=0.74 * mm9_size
  macs2_glocal = nrow(z.control)*macs2_params$extsize/mm9_effective_size

  # Combined noise
  macs2_noise_path = "data/macs2/Control_noise.bdg"
  system(stringr::str_glue("macs2 bdgcmp -m max -t {slocal} -c {llocal} -o {noise}", slocal=macs2_slocal_norm_path, llocal=macs2_llocal_norm_path, noise=macs2_noise_path))
  if(do_compare) {
    system(stringr::str_glue("macs2 bdgcmp -m max -t {noise} -c {dlocal} -o {noise}", dlocal=macs2_dlocal_path, noise=macs2_noise_path))
  }
  system(stringr::str_glue("macs2 bdgopt -i {noise} -m max -p {glocal} -o {noise}", glocal=macs2_glocal, noise=macs2_noise_path))
  system(stringr::str_glue("macs2 bdgopt -i {noise} -m multiply -p {factor} -o {noise}", noise=macs2_noise_path, factor=macs2_control_scale))

  # Pileup
  macs2_pileup_path = "data/macs2/Sample_pileup.bdg"
  system(stringr::str_glue("macs2 pileup -i {input} -f BED -B --extsize {format(extsize, scientific=F)} -o {output}", input=sample_bed_path, output=macs2_pileup_path, extsize=macs2_params$extsize/2))

  # Compare
  macs2_qvalue_path = "data/macs2/Sample_qvalue.bdg"
  system(stringr::str_glue("macs2 bdgcmp -t {input} -c {noise} -m qpois -o {output}", input=macs2_pileup_path, noise=macs2_noise_path, output=macs2_qvalue_path))

  # Call peaks
  macs2_peaks_path = "data/macs2/Sample_peaks.bed"
  system(stringr::str_glue("macs2 bdgpeakcall -i {qvalue} -c {cutoff} --min-length {format(peak_min_length, scientific=F)} --max-gap {format(peak_max_gap, scientific=F)} -o {output}", qvalue=macs2_qvalue_path, output=macs2_peaks_path, cutoff=-log10(macs2_params$qvalue_cutoff), peak_max_gap=macs2_params$peak_max_gap, peak_min_length=macs2_params$peak_min_length))

  # Read results file
  peaks_cols = cols(
    peak_chrom=col_character(), peak_start=col_double(), peak_end=col_double(), peak_length=col_character(), peak_abs_summit=col_double(),
    peak_score=col_character(), peak_fc=col_double(), peak_pvalue_log10=col_double(), peak_qvalue_log10=col_double(), peak_sammit_offset=col_character()
  )
  peaks = readr::read_tsv(macs2_peaks_path, skip=1, col_names=names(peaks_cols$cols), col_types=peaks_cols) %>%
    dplyr::mutate(peak_id=1:n())



  z_ranges = with(z.sample, IRanges::IRanges(start=break_start, end=break_end, break_id=break_id, break_chrom=break_chrom, break_bait_chrom=break_bait_chrom, break_exp_condition=break_exp_condition))
  peaks_ranges = with(peaks, IRanges::IRanges(start=peak_start, end=peak_end, peak_id=peak_id, peak_chrom=peak_chrom))
  peaks2z.map = as.data.frame(IRanges::mergeByOverlaps(z_ranges, peaks_ranges)) %>%
    dplyr::filter(break_chrom==peak_chrom) %>%
    dplyr::distinct(break_id, peak_id)
  z_clustered = z.sample %>%
    dplyr::left_join(peaks2z.map, by="break_id")

  z_clustered_peaks = z_clustered %>%
    dplyr::filter(!is.na(peak_id)) %>%
    dplyr::inner_join(peaks, by="peak_id") %>%
    dplyr::mutate(breaks_n=n()) %>%
    dplyr::group_by(break_bait_chrom) %>%
    dplyr::mutate(total_bait_breaks_n=n()) %>%
    dplyr::group_by(peak_id) %>%
    dplyr::mutate(total_peak_breaks_n=n()) %>%
    dplyr::group_by(break_bait_chrom) %>%
    dplyr::mutate(total_bait_breaks_n=n()) %>%
    dplyr::group_by(peak_id, break_bait_chrom) %>%
    dplyr::do((function(y){
      # y = z_clustered_peaks %>% dplyr::filter(break_bait_chrom=="chr1" & peak_id==424)
      yy<<-y
      x11 = nrow(y)
      x12 = y$total_peak_breaks_n[1] - x11
      x21 = y$total_bait_breaks_n[1] - x11
      x22 = y$breaks_n[1] - y$total_bait_breaks_n[1] - x12
      peak_bait_frac = x11 / y$total_peak_breaks_n[1]
      x = matrix(c(x11,x12,x21,x22), ncol=2)
      test = fisher.test(x)
      data.frame(
        peak_id=y$peak_id[1],
        break_bait_chrom=y$break_bait_chrom[1],
        peak_chrom=y$peak_chrom[1],
        peak_start=y$peak_start[1],
        peak_end=y$peak_end[1],
        peak_bait_breaks=x11,
        total_peak_breaks=y$total_peak_breaks_n[1],
        peak_bait_frac=peak_bait_frac,
        pvalue=test$p.value,
        odds=test$estimate)
    })(.))


  # z_clustered_peaks_ranges = with(z_clustered_peaks, IRanges::IRanges(start=peak_start-1e6, end=peak_end+1e6, peak_id=peak_id, peak_chrom=peak_chrom, peak_bait_chrom=break_bait_chrom, peak_bait_frac=peak_bait_frac))
  # offtargets_ranges_near = with(offtargets_df, IRanges::IRanges(start=offtarget_start, end=offtarget_end, offtarget_id=offtarget_id, offtarget_chrom=offtarget_chrom, offtarget_bait_chrom=bait_chrom))
  # z_clustered_peaks_ranges_map = as.data.frame(IRanges::mergeByOverlaps(offtargets_ranges_near, z_clustered_peaks_ranges)) %>%
  #   dplyr::filter(peak_chrom==offtarget_chrom & peak_bait_frac>0.1) %>%
  #   dplyr::mutate(peak_bait_frac=round(peak_bait_frac, 2)) %>%
  #   dplyr::select(peak_id, peak_bait_chrom, offtarget_bait_chrom, offtarget_chrom, peak_chrom, peak_bait_frac)
  #   dplyr::distinct(peak_id, offtarget_bait_chrom) %>%
  #   dplyr::mutate(is_offtarget=T)
  # junctions_df = junctions_df %>%
  #   dplyr::select(-dplyr::matches("is_near_bait")) %>%
  #   dplyr::left_join(offtarget2junctions.map, by="break_id") %>%
  #   dplyr::mutate(is_near_bait=tidyr::replace_na(is_near_bait, F))


  z_clustered_peaks %>%
    dplyr::filter(peak_bait_breaks>20 & peak_bait_frac>0.2 & odds>1 & pvalue<0.05) %>%
    dplyr::group_by(peak_chrom, peak_start, peak_end, peak_id) %>%
    dplyr::arrange(dplyr::desc(peak_bait_frac)) %>%
    dplyr::summarise(peak_bait_content=paste(paste0(break_bait_chrom, "=", round(peak_bait_frac, 2)), collapse=", "), score="0", peak_strand="+") %>%
    dplyr::select(peak_chrom, peak_start, peak_end, peak_bait_content, score, peak_strand) %>%
    readr::write_tsv("data/macs2/enriched_peaks.bed", col_names=F)

  list(peaks=peaks, peaks_bait_stats=z_clustered_peaks)
}

read_breaks_experiments = function(path, chromosomes_map) {
    #
  # Read breaks
  #
  tlx_cols = cols(
    Qname=col_character(), JuncID=col_character(), Rname=col_character(), Junction=col_double(),
    Strand=col_character(), Rstart=col_double(), Rend=col_double(),
    B_Rname=col_character(), B_Rstart=col_double(), B_Rend=col_double(), B_Strand=col_double(),
    B_Qstart=col_double(), B_Qend=col_double(), Qstart=col_double(), Qend=col_double(), Qlen=col_double(),
    B_Cigar=col_character(), Cigar=col_character(), Seq=col_character(), J_Seq=col_character(), Barcode=col_logical(),
    unaligned=col_double(), baitonly=col_double(), uncut=col_double(), misprimed=col_double(), freqcut=col_double(),
    largegap=col_double(), mapqual=col_double(), breaksite=col_double(), sequential=col_double(), repeatseq=col_double(), duplicate=col_double(), Note=col_character()
  )

  # Rename Chr14 to Alt289 (Was Alt289 and only PW255 was Alt287)
  junctions_ann.excluded = "PW121|PW246|PW247|PW248|PW249|JK096|JK097|JK098|JK099|JK100|JK101"
  junctions_ann.excluded = "NO FILTER"
  junctions_ann = data.frame(break_file=list.files(path, full.names=T, recursive=T, pattern="*.tlx")) %>%
     dplyr::filter(!grepl(junctions_ann.excluded, break_file)) %>%
     dplyr::mutate(
       break_bait_chrom1=tolower(basename(dirname(break_file))),
       break_bait_chrom = chromosomes_map[tolower(basename(dirname(break_file)))],
       break_condition = ifelse(grepl("DMSO", break_file), "Control", "Sample"),
       break_tech_sample = gsub("_.*", "", basename(break_file)),
       break_bio_sample=gsub(".*_((Alt|R)[0-9]+).*", "\\1", basename(break_file), ignore.case=T, perl=T)
     ) %>%
    dplyr::group_by(break_bait_chrom, break_bait_chrom1, break_condition) %>%
    dplyr::mutate(break_replicate=1:n()) %>%
    dplyr::ungroup()


  junctions_df = data.frame()
  for(i in 1:nrow(junctions_ann)) {
    tlx_ann = junctions_ann[i,,drop=F]
    tlx = cbind(tlx_ann, readr::read_tsv(tlx_ann$break_file, col_types=tlx_cols))
    bed = tlx %>%
      dplyr::mutate(break_chrom=chromosomes_map[Rname]) %>%
      dplyr::mutate(break_score=".", break_strand=ifelse(Strand=="-1", "-", "+")) %>%
      dplyr::mutate(break_exp_condition=tlx_ann$break_condition, break_file=tlx_ann$break_file, break_bait_chrom=tlx_ann$break_bait_chrom, break_replicate=tlx_ann$break_replicate) %>%
      dplyr::select(break_bait_chrom, break_exp_condition, break_chrom=Rname, break_start=Junction, break_end=Junction, break_name=Qname, break_score, break_strand, break_file)
    junctions_df = rbind(junctions_df, bed)
  }

  junctions_df %>% dplyr::mutate(break_id=1:n())
}

main = function() {
  chromosomes_map_df = readr::read_tsv("data/mm9_chromosomes_synonyms.tsv")
  chromosomes_map = chromosomes_map_df$chrom
  names(chromosomes_map) = chromosomes_map_df$chrom_synonym

  # Read offtargets
  offtargets_all_df = readr::read_tsv("data/offtargets_predicted.tsv") %>%
    dplyr::mutate(offtarget_id=1:n())
  offtargets_all_df %>%
    dplyr::filter(bait_chrom==offtarget_chrom) %>%
    dplyr::mutate(offtarget_name=stringr::str_glue("{bait}_{mismatches}_{seq}", bait=bait_name, mismatches=offtarget_mismatches, seq=offtarget_sequence), offtarget_score=0) %>%
    dplyr::select(offtarget_chrom, offtarget_start, offtarget_end, offtarget_name, offtarget_score, offtarget_strand) %>%
    readr::write_tsv("data/macs2/offtargets.bed", col_names=F)

  offtargets_df = offtargets_all_df %>%
    dplyr::mutate(bait_chrom=chromosomes_map[tolower(bait_chrom)], offtarget_chrom=chromosomes_map[tolower(offtarget_chrom)]) %>%
    dplyr::filter(offtarget_mismatches==0) %>%
    dplyr::mutate(
      offtarget_start=ifelse(bait_chrom=="chr15", 61818880, offtarget_start),
      offtarget_end=ifelse(bait_chrom=="chr15", 61818881, offtarget_end))

  bed_cols = cols(chrom=col_character(), start=col_double(), end=col_double(), name=col_character(), score=col_character(), strand=col_character())

  junctions_df = read_breaks_experiments("data/breaks_tlx", chromosomes_map=chromosomes_map) %>%
    dplyr::filter(grepl("chr[0-9]", break_chrom))


  #
  # Find breaks close to bait
  #
  junctions_ranges = with(junctions_df, IRanges::IRanges(start=break_start, end=break_end, break_id=break_id, break_chrom=break_chrom, break_bait_chrom=break_bait_chrom, break_exp_condition=break_exp_condition))
  offtargets_ranges_near = with(offtargets_df %>% dplyr::filter(offtarget_mismatches==0), IRanges::IRanges(start=offtarget_start-2e6, end=offtarget_end+2e6, offtarget_id=offtarget_id, offtarget_chrom=offtarget_chrom, offtarget_bait_chrom=bait_chrom))
  offtarget2junctions.map = as.data.frame(IRanges::mergeByOverlaps(offtargets_ranges_near, junctions_ranges)) %>%
    dplyr::filter(break_chrom==offtarget_chrom & break_bait_chrom==offtarget_bait_chrom) %>%
    dplyr::mutate(is_near_bait=T) %>%
    dplyr::distinct(break_id, is_near_bait)
  junctions_df = junctions_df %>%
    dplyr::select(-dplyr::matches("is_near_bait")) %>%
    dplyr::left_join(offtarget2junctions.map, by="break_id") %>%
    dplyr::mutate(is_near_bait=tidyr::replace_na(is_near_bait, F))

  #
  # Calculate normalization accross experiment
  #
  scale_factor_df = junctions_df %>%
    dplyr::filter(!grepl("chrX|chrY", break_bait_chrom)) %>%
    dplyr::group_by(break_exp_condition, break_bait_chrom) %>%
    dplyr::summarise(libsize=n()) %>%
    dplyr::mutate(libsize_median=median(libsize)) %>%
    dplyr::mutate(scale_factor=libsize_median/libsize) %>%
    dplyr::select(break_exp_condition, break_bait_chrom, scale_factor)
  junctions_df = junctions_df %>%
    dplyr::select(-dplyr::matches("scale_factor")) %>%
    dplyr::inner_join(scale_factor_df, by=c("break_exp_condition", "break_bait_chrom"))


  # macs2_params = list(extsize=1e5, llocal=1e8, slocal=2e6, qvalue_cutoff=0.05, peak_max_gap=1e5, peak_min_length=20)
  sample_df = junctions_df %>% dplyr::filter(break_exp_condition=="Sample")
  control_df = junctions_df %>% dplyr::filter(break_exp_condition=="Control")
  p = call_peaks(sample_df, control_df)

  #
  # Run MACS2
  #
  z = junctions_df %>%
    dplyr::mutate(break_name=stringr::str_glue("{bait}_{name}", bait=break_bait_chrom, name=break_name)) %>%
    dplyr::filter(!grepl("chrX|chrY", break_bait_chrom)) %>%
    dplyr::filter(!is_near_bait & break_bait_chrom != break_chrom)
  z.sample = z %>% dplyr::filter(break_exp_condition=="Sample")
  z.control = z %>% dplyr::filter(break_exp_condition=="Control")


  #macs2_params = list(extsize=1e5, qvalue=0.01, llocal=1e7, slocal=2e6, shift=5e4, qvalue_cutoff=0.05, peak_max_gap=1e5, peak_min_length=1)
  macs2_params = list(extsize=1e5, llocal=1e8, slocal=2e6, qvalue_cutoff=0.05, peak_max_gap=1e5, peak_min_length=1)
  call_peaks(z.sample, macs2_params, z.control)





  pdf("reports/breaks_distance_to_bait.pdf", width=16, height=10)
  junctions_df %>%
  dplyr::group_by(break_file, break_bait_chrom, break_exp_condition) %>%
  dplyr::summarise(near_bait_count=sum(is_near_bait), far_bait_count=sum(!is_near_bait)) %>%
  reshape2::melt(measure.vars=c("near_bait_count", "far_bait_count")) %>%
  ggplot() +
    geom_bar(aes(x=break_file, fill=variable, y=value), stat="identity", position="stack") +
    coord_flip() +
    facet_wrap(~break_bait_chrom, scale="free_y", ncol=2)
  dev.off()

  #
  # Create a bed file where all breaks from outside bait region are put together
  #

      bed = list(
        "all" = z,
        "plus" = z %>% dplyr::filter(break_strand=="+"),
        "minus" = z %>% dplyr::filter(break_strand=="-")
      )
      bed_path = sapply(names(bed), function(s) stringr::str_glue("data/breakensembl/{exp}_{strand}.bed", exp=bed[[s]]$break_exp_condition[1], strand=s))
      o = sapply(names(bed), function(s) readr::write_tsv(bed[[s]] %>% dplyr::select(break_chrom, break_start, break_end, break_name, break_score, break_strand), bed_path[[s]], col_names=F))

# chr7:240 (R=0.5)        chr6:320 (R=0.1)        chr10:790 (R=0.8)       chr16: 300 (R=0.1)      chr9: 300  (R=0.3)      chr8: 450  (R=0.7)     chr3: 520 (R=0.9)


    #
    # Select right data
    #
    z = junctions_df %>%
      dplyr::mutate(break_name=stringr::str_glue("{bait}_{name}", bait=break_bait_chrom, name=break_name)) %>%
      dplyr::filter(!grepl("chrX|chrY", break_bait_chrom)) %>%
      dplyr::filter(!is_near_bait & break_bait_chrom != break_chrom)
    z.sample = z %>% dplyr::filter(break_exp_condition=="Sample")
    z.control = z %>% dplyr::filter(break_exp_condition=="Control")

    macs2_params = list(extsize=1e4, qvalue=0.01, llocal=1e7, slocal=2e6, shift=5e4, qvalue_cutoff=0.05, peak_max_gap=1e5, peak_min_length=1)
    call_peaks(z.sample, z.control, macs2_params)

    macs2_control_scale = nrow(junctions_df %>% dplyr::filter(break_exp_condition=="Sample"))/nrow(junctions_df %>% dplyr::filter(break_exp_condition=="Control"))

    #
    # Save bed files
    #
    sample_bed_path = "data/breakensembl/Sample.bed"
    readr::write_tsv(z.sample %>% dplyr::select(break_chrom, break_start, break_end, break_name, break_score, break_strand), sample_bed_path, col_names=F)
    control_bed_path = "data/breakensembl/Control.bed"
    readr::write_tsv(z.control %>% dplyr::select(break_chrom, break_start, break_end, break_name, break_score, break_strand), control_bed_path, col_names=F)

    macs2_params = list(extsize=1e4, qvalue=0.01, llocal=1e7, slocal=2e6, shift=5e4, qvalue_cutoff=0.05, peak_max_gap=1e5, peak_min_length=1)



    macs2_output_path = sapply(names(bed), function(s) stringr::str_glue("data/macs2/{exp}_{strand}_peaks.xls", exp=bed[[s]]$break_exp_condition[1], strand=s))
    macs2_cmd = stringr::str_glue("macs2 callpeak -t {input} -f BED -g mm --keep-dup all -n {output} -B --outdir {outdir} --nomodel --extsize {format(extsize, scientific=F)} -q {format(qvalue, scientific=F)} --llocal {format(llocal, scientific=F)} {pipe}",
      input="data/breakensembl/Sample.bed",
      outdir="data/macs2",
      output="test",
      extsize=macs2_params$extsize,
      qvalue=macs2_params$qvalue,
      llocal=macs2_params$llocal,
      pipe="", pipe2="2>/dev/null"
    )
    system(macs2_cmd)

    macs2_output = readr::read_tsv(macs2_output_path[["all"]], col_names=names(macs2_cols$cols), comment="#", skip=23, col_types=macs2_cols) %>%
      dplyr::mutate(peak_method="macs2") %>%
      dplyr::select(peak_method, peak_chrom, peak_start, peak_end, peak_name, peak_score)

    macs2_output %>%
        dplyr::mutate(peak_strand="+") %>%
        dplyr::select(peak_chrom, peak_start, peak_end, peak_name, peak_score, peak_strand) %>%
        readr::write_tsv(gsub("_peaks.xls", "_peaks.bed", macs2_output_path[["all"]]), col_names=F)



    islands = dplyr::bind_rows(islands, macs2_output)
}







old = function()
{

    ##################################
    # Create temporary files
    ##################################
    path_bed = tempfile()
    readr::write_tsv(z %>% dplyr::select(break_chrom, break_start, break_end, break_name, break_score, break_strand), path_bed, col_names=F)

    path_bed_far = tempfile()
    readr::write_tsv(z %>% dplyr::filter(!is_near_bait) %>% dplyr::select(break_chrom, break_start, break_end, break_name, break_score, break_strand), path_bed_far, col_names=F)

junctions_df.f = junctions_df
junctions_g = junctions_df %>%
  dplyr::arrange(break_start, break_end) %>%
  dplyr::group_by(break_bait_chrom, break_exp_condition, break_chrom)
all_islands = junctions_g %>% dplyr::do((function(z, g){
    gg <<- g
    g = c("break_bait_chrom", "break_exp_condition", "break_chrom")
    z = junctions_df %>% dplyr::filter(break_bait_chrom=="chr15" & break_chrom=="chr15" & break_exp_condition=="Sample")
    zz<<-z

    writeLines(paste0(g, "=", z[1,g], collapse=" "))

    ##################################
    # Create temporary files
    ##################################
    path_bed = tempfile()
    readr::write_tsv(z %>% dplyr::select(break_chrom, break_start, break_end, break_name, break_score, break_strand), path_bed, col_names=F)

    path_bed_far = tempfile()
    readr::write_tsv(z %>% dplyr::filter(!is_near_bait) %>% dplyr::select(break_chrom, break_start, break_end, break_name, break_score, break_strand), path_bed_far, col_names=F)

    islands = data.frame()

    ##################################
    # DBSCAN
    ##################################
    if(F) {
      writeLines("DBSCAN")
      dbscan_params = list(minPts=20, eps_cl=1e5)
      res.optics = dbscan::optics(matrix(z$break_start), minPts=dbscan_params$minPts, eps=2e7)
      res.dbscan = dbscan::extractDBSCAN(res.optics, eps_cl=dbscan_params$eps_cl)
      #plot(res.dbscan)
      res.dbscan = dbscan::dbscan(matrix(z$break_start), minPts=dbscan_params$minPts, eps=dbscan_params$eps_cl)
      dbscan_output = z %>%
        dplyr::mutate(cluster=res.dbscan$cluster) %>%
        dplyr::filter(cluster>0) %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(peak_method="dbscan", peak_start=min(break_start), peak_end=max(break_end), peak_score=0) %>%
        dplyr::select(-cluster)
      islands = rbind(islands, dbscan_output)
    }

    if(F) {
      ##################################
      # SICER
      ##################################
      writeLines("SICER")
      sicer_params = list(window=10000, gap=30000, e_value=0.1, fdr=0.01)
      sicer_control_results.cols = cols("peak_chrom"=col_character(), peak_start=col_double(), peak_end=col_double(), peak_score = col_double())
      sicer_cmd = stringr::str_glue("sh SICER1.1/SICER/SICER-rb.sh {input_dir} {input} {output_dir} {genome} {format(r_threshold, scientific=F)} {format(window_size, scientific=F)} {format(fragment_size, scientific=F)} {format(fraction, scientific=F)} {format(gap_size, scientific=F)} {format(e_value, scientific=F)}",
                       input_dir=dirname(path_bed_far), input=basename(path_bed_far), output_dir="data/sicer1",
                       genome="mm9", fraction=0.74, r_threshold=5, window_size=sicer_params["window"], fragment_size=20, gap_size=sicer_params["gap"], e_value=sicer_params["e_value"]
      )
      system(sicer_cmd)
      sicer_output_path = stringr::str_glue("data/sicer1/{sample}-W{format(window, scientific=F)}-G{format(gap, scientific=F)}-E{format(e_value, scientific=F)}.scoreisland", window=sicer_params["window"], gap=sicer_params["gap"], e_value=sicer_params["e_value"], sample=gsub("\\.bed$", "", basename(path_bed_far)))
      if(file.exists(sicer_output_path) && file.size(sicer_output_path)>0) {
        sicer_output = readr::read_tsv(sicer_output_path, col_names=names(sicer_control_results.cols$cols), col_types=sicer_control_results.cols) %>%
          dplyr::mutate(peak_method="sicer") %>%
          dplyr::select(peak_method, peak_start, peak_end, peak_score)
        islands = dplyr::bind_rows(islands, sicer_output)
      }
    }

    ##################################
    # MACS2
    ##################################
    writeLines("MACS2")
    macs2_output_path = tempfile(tmpdir="data/macs2")
    macs2_params = list(extsize=20, qvalue=0.01, llocal=3e4, shift=5e4)
    #macs2_params$llocal = macs2_params$extsize*5
    #macs2_params$shift = macs2_params$extsize/2
    macs2_cmd = stringr::str_glue("macs2 callpeak -t {input} -f BED -g mm --keep-dup all -n {output} --outdir {outdir} --nomodel --extsize {extsize} -q {qvalue} --shift {format(shift, scientific=F)} --llocal {format(llocal, scientific=F)}  2>/dev/null",
      input=path_bed_far, outdir=dirname(macs2_output_path), output=basename(macs2_output_path), extsize=macs2_params$extsize, qvalue=macs2_params$qvalue, llocal=macs2_params$llocal, shift=macs2_params$shift)
    system(macs2_cmd)

    macs2_pileup_cmd = stringr::str_glue("macs2 pileup -i {input} -f BED --extsize {format(extsize, scientific=F)} --outdir {outdir} -o {output}",
      input=path_bed, outdir=dirname(macs2_output_path), output=paste0(basename(macs2_output_path), ".bdg"), extsize=1e4)
    system(macs2_pileup_cmd)

    readr::write_tsv(junctions_df %>% dplyr::filter(break_bait_chrom=="chr15" & break_exp_condition=="Sample") %>% dplyr::select(break_chrom, break_start, break_end, break_name, break_score, break_strand), path_bed, col_names=F)
    macs2_predictd_cmd = stringr::str_glue("macs2 predictd -i {input} -g mm -m 1 50000 --bw 50",
      input=path_bed_far)
    system(macs2_predictd_cmd)



    macs2_output = readr::read_tsv(paste0(macs2_output_path, "_peaks.xls"), col_names=names(macs2_cols$cols), comment="#", skip=23, col_types=macs2_cols) %>%
      dplyr::mutate(peak_method="macs2") %>%
      dplyr::select(peak_method, peak_chrom, peak_start, peak_end, peak_score)
    islands = dplyr::bind_rows(islands, macs2_output)

    print(macs2_output_path)
    breaks = readr::read_tsv(path_bed_far, col_names=names(bed_cols$cols), col_types=bed_cols)
    breaks %>% dplyr::filter(dplyr::between(start, 57585942, 57670940))
    readr::read_tsv(paste0(macs2_output_path, "_peaks.xls"), col_names=names(macs2_cols$cols), comment="#", skip=23, col_types=macs2_cols) %>%
        dplyr::mutate(peak_strand="+") %>%
        dplyr::select(peak_chrom, peak_start, peak_end,  peak_name, peak_score, peak_strand) %>%
        readr::write_tsv(paste0(macs2_output_path, "_peaks.bed"), col_names=F)


    #z_ranges.cluster = IRanges::IRanges(start=z.cluster$cluster_start, end=z.cluster$cluster_end)
    #plotRanges(z_ranges.cluster, xlim=c(0, max(z.cluster$cluster_end)))
    islands %>% dplyr::bind_cols(z[1,g, drop=F])
  })(., group_vars(junctions_g)))

plot(density(all_islands.f$peak_end-all_islands.f$peak_start))

all_islands.f = all_islands %>% dplyr::filter(break_bait_chrom!=break_chrom) %>% dplyr::mutate(i=1:n())
all_islands_cross = all_islands.f %>%
  dplyr::inner_join(all_islands.f, by=setdiff(c("peak_method", group_vars(junctions_g)), "break_bait_chrom")) %>%
  #dplyr::filter(i.y>i.x) %>%
  dplyr::filter(break_bait_chrom.x!=break_bait_chrom.y) %>%
  dplyr::filter(dplyr::between(peak_start.x, peak_start.y, peak_end.y) | dplyr::between(peak_start.y, peak_start.x, peak_end.x))

  # +

ggplot() +
  ggridges::geom_density_ridges(aes(y=break_chrom, x=peak_start.x, color=break_bait_chrom.x), bandwidth=1e5, data=all_islands_cross, alpha=0.1) +
  geom_point(aes(y=bait_chrom, x=offtarget_start, color=bait_chrom), data=offtargets_df, size=5, alpha=0.5) +
  facet_wrap(~break_exp_condition)
}




  #
  # sample_sizes = sample_df.coverage %>%
  #   dplyr::group_by(seqnames) %>%
  #   dplyr::summarise(sample_start=min(start[score>=quantile(score, 0.05)], na.rm=T), sample_end=max(end[score>quantile(score, 0.05)], na.rm=T)) %>%
  #   dplyr::select(seqnames, sample_start, sample_end)
  #
  # #
  # # Calculate baseline
  # #
  # control_ranges.bins = unlist(tileGenome(mm9, tilewidth=binsize))
  # control_ranges.bins$coverage = GenomicRanges::countOverlaps(control_ranges.bins, control_ranges.extended)
  # baseline_df.extended = data.table::as.data.table(control_ranges.bins) %>%
  #   # dplyr::inner_join(sample_sizes, by="seqnames") %>%
  #   # dplyr::filter(start>=sample_start & end<=sample_end) %>%
  #   dplyr::arrange(seqnames, start) %>%
  #   dplyr::group_by(seqnames) %>%
  #   dplyr::mutate(is_q2=coverage>quantile(coverage, 0.2), is_q8=coverage<quantile(coverage, 0.8), is_inner_quantile=is_q2 & is_q8) %>%
  #   dplyr::mutate(coverage_median=median(coverage[is_inner_quantile], na.rm=T), coverage_sd=sd(coverage[is_inner_quantile], na.rm=T)) %>%
  #   dplyr::mutate(coverage_min=coverage_median-2*coverage_sd) %>%
  #   dplyr::mutate(coverage_mod=ifelse(coverage<coverage_min, NA_real_, coverage)) %>%
  #   dplyr::mutate(coverage_mod=zoo::na.fill(coverage_mod, "extend")) %>%
  #   dplyr::do((function(z){
  #     zz<<-z
  #     # asdsad()
  #     x = matrix(as.numeric(z$coverage_mod), nrow=1)
  #     colnames(x) = z$start
  #     bc.irls = baseline::baseline(x, method="medianWindow", hws=1e6/binsize, hwm=1e6/binsize) # !!!!
  #     # bc.irls = baseline::baseline(x, method="medianWindow", hws=2e6/binsize, hwm=5e6/binsize) # !!!!
  #     # bc.irls = baseline(x, method="rollingBall", wm=1e6/binsize, ws=2e6/binsize)
  #     # bc.irls = baseline::baseline(x, method="irls", lambda1=7, lambda2=13) # 1e2
  #     # bc.irls = baseline::baseline(x, method="irls" lambda1=5, lambda2=10)
  #     # bc.irls = baseline::baseline(x, method="irls", lambda1=4, lambda2=9) # For binsize = 1e4 / extend = 1e5 !!!!!
  #     # bc.irls = baseline::baseline(x, method="irls", lambda1=3, lambda2=7) # For binsize = 1e4 / extend = 1e5 !!!!!
  #     z$coverage_smooth = as.numeric(bc.irls@baseline)
  #     z
  #   })(.)) %>%
  #   dplyr::ungroup()
