library(ggplot2)

scale_breaks = function(x) {
    breaks = seq(0, max(x), 1e7)
    names(breaks) = paste0(breaks/1e6, "Mb")
    breaks
}

ggplot_karyoplot = function(coverage_df=NULL, baseline_df=NULL, peaks_df=NULL, max_coverage="auto", colors=NULL) {
  if(max_coverage=="auto") {
    if(!is.null(baseline_df)) { max_coverage = median(baseline_df$coverage_smooth)*12 } else {
      if(!is.null(baseline_df)) { max_coverage = median(coverage_df$coverage)*12 } else {
        max_coverage = 50
      }
    }
  }

  count_annotations = 0
  if(!is.null(peaks_df) & nrow(peaks_df)>0) {
    count_annotations = count_annotations + length(unique(peaks_df$condition))
  }
  if(!is.null(coverage_df) & nrow(coverage_df)>0) {
    count_annotations = count_annotations + length(unique(coverage_df$condition))
  }

  p = ggplot() +
    coord_cartesian(ylim=c(-(count_annotations+2)*0.05*max_coverage ,max_coverage)) +
    scale_color_manual(values=colors) +
    scale_fill_manual(values=colors) +
    scale_x_continuous(breaks=scale_breaks) +
    facet_wrap(~seqnames, scales="free_x", ncol=1) +
    labs(x="")

  if(!is.null(coverage_df)) {
    coverage_df = coverage_df %>% dplyr::mutate(coverage_limited=ifelse(coverage>max_coverage, max_coverage, coverage))
    if(!("condition" %in% colnames(coverage_df))) {
      coverage_df$condition = "Coverage"
    }
    coverage_df$peaks_offset = 1
    if(!is.null(peaks_df) & nrow(peaks_df)>0) {
      coverage_df$peaks_offset = length(unique(coverage_df$condition)) + 3
    }
    coverage_df$condition_number = as.numeric(as.factor(as.character(coverage_df$condition)))
    p = p + geom_line(aes(y=coverage_limited, x=start, color=condition), data=coverage_df, alpha=0.5, size=0.5)
    p = p + geom_rect(aes(xmin=start, xmax=end, fill=condition, alpha=coverage, ymin=-(peaks_offset+1+condition_number)*max_coverage*0.05, ymax=-(peaks_offset+condition_number)*max_coverage*0.05), data=coverage_df)
  }
  if(!is.null(peaks_df) & nrow(peaks_df)>0) {
    if("score" %in% colnames(peaks_df)) {
      peaks_df$score = max(coverage_df$coverage)
    }
    peaks_df$score = ifelse(is.na(peaks_df$score), max(coverage_df$coverage), peaks_df$score)

    peaks_df$condition_number = as.numeric(as.factor(as.character(peaks_df$condition)))
    p = p + geom_rect(aes(xmin=peak_start, xmax=peak_end, fill=condition, ymin=-(0+condition_number)*max_coverage*0.05, ymax=-(1+condition_number)*max_coverage*0.05, alpha=score), data=peaks_df)
    if("peak_summit_abs" %in% colnames(peaks_df)) {
      # peaks_df$peak_summit_abs_display = ifelse(peaks_df$peak_summit_abs>=1e3, with(peaks_df, paste0(floor(peak_summit_abs/10^log10(peak_summit_abs)), "e", floor(log10(peak_summit_abs)))), as.character(round(peaks_df$peak_summit_abs)))
      peaks_df$peak_summit_abs_display = round(peaks_df$peak_summit_qvalue)
      p = p + ggrepel::geom_text_repel(aes(x=peak_start, y=-(0.5+condition_number)*max_coverage*0.05, label=peak_summit_abs_display), color="#333333", data=peaks_df)
    }
  }
  if(!is.null(baseline_df)) {
    if(!("condition" %in% colnames(baseline_df))) { baseline_df$condition = "Coverage" }
    baseline_df$condition = baseline_df$condition
    p = p + geom_line(aes(y=coverage_smooth, x=start, color=condition), data=baseline_df, alpha=0.8)
    p = p + geom_line(aes(y=coverage_th01, x=start, color=condition), data=baseline_df, linetype="dashed", alpha=0.8)
  }

  p
}