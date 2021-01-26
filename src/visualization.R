library(ggplot2)

scale_breaks = function(x) {
    breaks = seq(0, max(x), 1e7)
    names(breaks) = paste0(breaks/1e6, "Mb")
    breaks
}

ggplot_karyoplot = function(coverage_df=NULL, baseline_df=NULL, islands_df=NULL, max_coverage="auto", colors=NULL) {
  if(max_coverage=="auto") {
    if(!is.null(baseline_df)) { max_coverage = median(baseline_df$coverage_smooth)*12 } else {
      if(!is.null(baseline_df)) { max_coverage = median(coverage_df$coverage)*12 } else {
        max_coverage = 50
      }
    }
  }

  count_annotations = 0
  if(!is.null(islands_df) & nrow(islands_df)>0) {
    count_annotations = count_annotations + length(unique(islands_df$experimental_condition))
  }
  if(!is.null(coverage_df) & nrow(coverage_df)>0) {
    count_annotations = count_annotations + length(unique(coverage_df$experimental_condition))
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
    if(!("experimental_condition" %in% colnames(coverage_df))) {
      coverage_df$experimental_condition = "Coverage"
    }
    coverage_df$islands_offset = 1
    if(!is.null(islands_df) & nrow(islands_df)>0) {
      coverage_df$islands_offset = length(unique(coverage_df$experimental_condition)) + 3
    }
    coverage_df$experimental_condition_number = as.numeric(as.factor(as.character(coverage_df$experimental_condition)))
    p = p + geom_line(aes(y=coverage_limited, x=start, color=experimental_condition), data=coverage_df, alpha=0.5, size=0.5)
    p = p + geom_rect(aes(xmin=start, xmax=end, fill=experimental_condition, alpha=coverage, ymin=-(islands_offset+1+experimental_condition_number)*max_coverage*0.05, ymax=-(islands_offset+experimental_condition_number)*max_coverage*0.05), data=coverage_df)
  }
  if(!is.null(islands_df) & nrow(islands_df)>0) {
    if("score" %in% colnames(islands_df)) {
      islands_df$score = max(coverage_df$coverage)
    }
    islands_df$score = ifelse(is.na(islands_df$score), max(coverage_df$coverage), islands_df$score)

    islands_df$experimental_condition_number = as.numeric(as.factor(as.character(islands_df$experimental_condition)))
    p = p + geom_rect(aes(xmin=island_start, xmax=island_end, fill=experimental_condition, ymin=-(0+experimental_condition_number)*max_coverage*0.05, ymax=-(1+experimental_condition_number)*max_coverage*0.05, alpha=score), data=islands_df)
    if("island_summit_abs" %in% colnames(islands_df)) {
      # islands_df$island_summit_abs_display = ifelse(islands_df$island_summit_abs>=1e3, with(islands_df, paste0(floor(island_summit_abs/10^log10(island_summit_abs)), "e", floor(log10(island_summit_abs)))), as.character(round(islands_df$island_summit_abs)))
      islands_df$island_summit_abs_display = round(islands_df$island_summit_qvalue)
      p = p + ggrepel::geom_text_repel(aes(x=island_start, y=-(0.5+experimental_condition_number)*max_coverage*0.05, label=island_summit_abs_display), color="#333333", data=islands_df)
    }
  }
  if(!is.null(baseline_df)) {
    if(!("experimental_condition" %in% colnames(baseline_df))) { baseline_df$experimental_condition = "Coverage" }
    baseline_df$experimental_condition = baseline_df$experimental_condition
    p = p + geom_line(aes(y=coverage_smooth, x=start, color=experimental_condition), data=baseline_df, alpha=0.8)
    p = p + geom_line(aes(y=coverage_th01, x=start, color=experimental_condition), data=baseline_df, linetype="dashed", alpha=0.8)
  }

  p
}