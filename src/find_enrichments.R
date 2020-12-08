library(readr)
library(dplyr)
library(dbscan)
library(stringr)
library(ggplot2)
library(ROCR)

plotRanges <- function(x, xlim=x, main=deparse(substitute(x)), col="black", sep=0.5, ...)
{
    height <- 1
    if (is(xlim, "IntegerRanges"))
    xlim <- c(min(start(xlim)), max(end(xlim)))
    bins <- disjointBins(IRanges(start(x), end(x) + 1))
    plot.new()
    plot.window(xlim, c(0, max(bins)*(height + sep)))
    ybottom <- bins * (sep + height) - height
    rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col=col, ...)
    title(main)
    axis(1)
}

bait_region_size = 4e6

# breaksites_df = readr::read_tsv("data/breaksites_primers.tsv")
offtargets_df = readr::read_tsv("data/offtargets_predicted.tsv") %>% # Doesn't have position for Chr6
  dplyr::filter(offtarget_mismatches<5)

breaksites_df = offtargets_df %>% # Doesn't have position for Chr6
  dplyr::filter(offtarget_mismatches<5) %>%
  dplyr::select(bait_chrom, bait_position=offtarget_start)

tempdir = "tmp"
enrichment_outdir = "data/enrichment"
junctions_cols = c("junction_chrom", "junction_start", "junction_end", "junction_name", "junction_score", "junction_strand")


dir.create("data/epic2", recursive=T, showWarnings=F)
dir.create("data/sicer1", recursive=T, showWarnings=F)
dir.create("data/sicer2", recursive=T, showWarnings=F)

sicer_params = list(window=30000, gap=90000, e_value=0.1, fdr=0.01)
epic2_control_results.cols = cols("peak_chrom"=col_character(), peak_start=col_double(), peak_end=col_double(), peak_breaks_count=col_double(), peak_score = col_double(), peak_strand=col_character())
sicer_control_results.cols = cols("peak_chrom"=col_character(), peak_start=col_double(), peak_end=col_double(), peak_score = col_double())
epic2_control_results = data.frame()
sicer_control_results = data.frame()
sicer_diff_results = data.frame()
for(breaks_bed in list.files("data/breaks", pattern="*_APH_no10kb_Merge.bed", full.names=T)) {
    x.bait_chrom = gsub("(chr)([^_]+).*$", "\\L\\1\\2", basename(breaks_bed), perl=T, ignore.case=T)
    x.sample = gsub("\\.bed$", "", basename(breaks_control_bed))
    breaks_control_bed = gsub("_APH_no10kb_Merge.bed", "_DMSO.bed", breaks_bed)
    breaks_filtered_bed = paste0("tmp/", basename(gsub("\\.bed$", "_filtered.bed", breaks_bed)))
    breaks_control_filtered_bed = paste0("tmp/", basename(gsub("\\.bed$", "_filtered.bed", breaks_control_bed)))
    epic2_control_bed = stringr::str_glue("data/epic2/{chrom}_both", chrom=x.bait_chrom)
    sicer_control_bed = stringr::str_glue("data/sicer1/{sample}-W{format(window, scientific=F)}-G{format(gap, scientific=F)}-E{format(e_value, scientific=F)}.scoreisland", window=sicer_params["window"], gap=sicer_params["gap"], e_value=sicer_params["e_value"], sample=gsub("\\.bed$", "", basename(breaks_control_filtered_bed)))
    sicer_diff_bed = stringr::str_glue("data/sicer1/{sample}-W{format(window, scientific=F)}-G{format(gap, scientific=F)}-islands-summary", window=sicer_params["window"], gap=sicer_params["gap"], fdr=sicer_params["fdr"], sample=gsub("\\.bed$", "", basename(breaks_filtered_bed)))

    writeLines(stringr::str_glue("Processing {bait}...", bait=x.bait_chrom))
    if(!file.exists(breaks_control_bed)) {
        writeLines("Control doesn't exist. Skipping...")
        next()
    }

    bed_cols = cols(chr=col_character(), start=col_double(), end=col_double(), name=col_character(), score=col_character(), strand=col_character())

    x.offtargets = offtargets_df %>% dplyr::filter(offtarget_mismatches==0 & bait_chrom==x.bait_chrom)
    readr::read_tsv(breaks_control_bed, col_names=names(bed_cols$cols), col_types=bed_cols) %>%
      dplyr::filter(!(start >= x.offtargets$offtarget_start-2e6 & end<=x.offtargets$offtarget_start+2e6 & x.offtargets$offtarget_chrom==chr)) %>%
      dplyr::mutate(start=ifelse(strand=="-", start-20, start), end=ifelse(strand=="+", end+20, end)) %>%
      dplyr::mutate(score=1) %>%
      readr::write_tsv(file=breaks_control_filtered_bed, col_names=F)

    readr::read_tsv(breaks_bed, col_names=names(bed_cols$cols), col_types=bed_cols) %>%
      dplyr::filter(!(start >= x.offtargets$offtarget_start-2e6 & end<=x.offtargets$offtarget_start+2e6 & x.offtargets$offtarget_chrom==chr)) %>%
      dplyr::mutate(start=ifelse(strand=="-", start-20, start), end=ifelse(strand=="+", end+20, end)) %>%
      dplyr::mutate(score=1) %>%
      readr::write_tsv(file=breaks_filtered_bed, col_names=F)

    breaks = dplyr::bind_rows(
      readr::read_tsv(breaks_bed, col_names=names(bed_cols$cols), col_types=bed_cols) %>% dplyr::mutate(break_exp_condition="Sample"),
      readr::read_tsv(breaks_control_bed, col_names=names(bed_cols$cols), col_types=bed_cols) %>% dplyr::mutate(break_exp_condition="Control")) %>%
      dplyr::mutate(end=start+1)
    breaks_clusters = breaks %>%
        dplyr::group_by(break_exp_condition, chr) %>%
        dplyr::do((function(z){
            #z = breaks_control %>% dplyr::filter(chr=="chr6")
            zz<<-z

            res.optics = dbscan::optics(matrix(z$start), minPts=50, eps=1e4)
            res.xi = dbscan::extractXi(res.optics, xi=0.1)
            #plot(res.xi)
          plot(res.xi$coredist, type="l", ylim=c(0, 100))
            if(is.null(res.xi$clusters_xi)) return(data.frame())

            z.clusters = data.frame(res.xi$clusters_xi) %>%
              dplyr::mutate(start=pmin(start, end), end=pmax(start, end)) %>%
              dplyr::mutate(optics_start=start, optics_end=end, start=z$start[optics_start], end=z$end[optics_end])

            z.clusters_ranges = with(z.clusters, IRanges::IRanges(start=start, end=end, cluster_id=cluster_id))
            z.clusters_reduced = as.data.frame(IRanges::reduce(z.clusters_ranges, min.gapwidth=1, drop.empty.ranges=F, with.revmap=T)) %>% dplyr::mutate(reduced_cluster_id=1:n())
            z.clusters_reduced_map = data.frame(reduced_cluster_id=c(rep(z.clusters_reduced$reduced_cluster_id, sapply(z.clusters_reduced$revmap, length)), 0), cluster_id=c(unlist(z.clusters_reduced$revmap), 0))
            z.clusters_reduced_sizes = z.clusters_reduced_map %>%
              dplyr::inner_join(z.clusters, by="cluster_id")  %>%
              dplyr::group_by(reduced_cluster_id) %>%
              dplyr::summarise(optics_start=min(optics_start), optics_end=max(optics_end))

            z.clusters_reduced = z.clusters_reduced %>%
              dplyr::inner_join(z.clusters_reduced_sizes, by="reduced_cluster_id") %>%
              dplyr::mutate(chr=z$chr[1], cluster_size=optics_end-optics_start+1) %>%
              dplyr::select(-revmap)

            z.clusters_reduced

            ## Debug merged
            #res.xi$cluster = data.frame(cluster_id=as.numeric(res.xi$cluster)) %>% dplyr::inner_join(z.clusters_reduced_map, by="cluster_id") %>% .$reduced_cluster_id
            #res.xi$clusters_xi = z.clusters_reduced %>%
            #  dplyr::select(start=optics_start, end=optics_end, cluster_id=reduced_cluster_id) %>%
            #  data.frame()
            #plot(res.xi)
            #
            #dbscan::hullplot(matrix(z$start), res.xi)
            #
            #plotRanges(z.clusters_ranges)
            #plotRanges(IRanges::reduce(z.clusters_ranges, min.gapwidth=1000, drop.empty.ranges=F, with.revmap=T))
            #
            #res.xi = dbscan::extractXi(res.optics, xi = 0.05)
            #plot(res.xi)

        })(.))

    optics_control_bed = stringr::str_glue("data/optics/{sample}_optics.bed", sample=x.sample)
    breaks_clusters_bed = breaks_clusters  %>%
      dplyr::mutate(cluster_id=stringr::str_glue("{sample}_{chr}_{cluster}", sample=x.sample, chr=chr, cluster=1:n()), score=0) %>%
      dplyr::select(chr, start, end, cluster_id, score)
    readr::write_tsv(breaks_clusters_bed, file=optics_control_bed, col_names=F)

    optics = IRanges(breaks_clusters$start, breaks_clusters$end, optics_chr=breaks_clusters$chr)
    sicer = IRanges(sicer_control_output$peak_start, sicer_control_output$peak_end, sicer_chr=sicer_control_output$peak_chrom)

    table(breaks_clusters$chr)
    table(sicer_control_output$peak_chrom)
    as.data.frame(IRanges::mergeByOverlaps(optics, sicer)) %>%
      dplyr::filter(optics_chr==sicer_chr)

    next
    output_control_enrichment = paste0(tempdir, "/", x.bait_chrom, "_control")

    # Run enrichment detection on control
    writeLines(stringr::str_glue("Calculate enrichments for control: '{output}'...", output=output_control_enrichment))

    # #Species=mm9, redundancy threshold=5, window size=30000, fragment size=1, effective genome fraction=0.74, gap size=90000, E-value=0.1
    # system(stringr::str_glue("venv/bin/epic2 --genome mm9 --bin-size 30000 --fragment-size 1 --gaps-allowed 3 -t {input} -c {control} -o {output}", input=breaks_filtered_bed, control=breaks_control_filtered_bed, output=output_control_bed))
    # system(stringr::str_glue("venv/bin/epic2 --genome mm9 --bin-size 30000 --fragment-size 1 --gaps-allowed 3 -t {input} -o {output}", input=breaks_filtered_bed, output=output_control_bed))
    # system(stringr::str_glue("venv/bin/sicer -s mm9 --redundancy_threshold 5 --window_size 30000 --fragment_size 1 -egf 0.74 --gap_size 90000 --e_value 0.1 --cpu 24  -t {input} -c {control} -o {output_dir}", input=breaks_filtered_bed, control=breaks_control_filtered_bed, output_dir="data/sicer2"))
    # system(stringr::str_glue("venv/bin/sicer -s mm9 --redundancy_threshold 5 --window_size 30000 --fragment_size 1 -egf 0.74 --gap_size 90000 --e_value 0.1 --cpu 24  -t {input} -o {output_dir}", input=breaks_filtered_bed, output_dir="data/sicer2"))
    #
    # system(stringr::str_glue("venv/bin/recognicer -s mm9 --redundancy_threshold 5 --window_size 100 --fragment_size 1 -egf 0.74 --gap_size 100 --e_value 0.1 --cpu 24  -t {input} -o {output_dir}", input=breaks_filtered_bed, output_dir="data/sicer2"))
    #
    #
    # # Species=mm9, redundancy threshold=5, window size=30000, fragment size=1, effective genome fraction=0.74, gap size=90000, FDR=0.01
    # system(stringr::str_glue("sh SICER1.1/SICER/SICER.sh {input_dir} {input} {control} {output_dir} {genome} {r_threshold} {window_size} {fragment_size} {fraction} {gap_size} {fdr}",
    #                          input_dir=dirname(breaks_filtered_bed),
    #                          input=basename(breaks_filtered_bed),
    #                          control=basename(breaks_control_filtered_bed),
    #                          output_dir="data/sicer1",
    #                          genome="mm9",
    #                          fraction=0.74,
    #                          r_threshold=5,
    #                          window_size=30000,
    #                          fragment_size=1,
    #                          gap_size=90000,
    #                          fdr=0.01
    # ))
    #
    # system("rm data/sicer1/*")
    #
    # # ORIGINAL FROM WEI

    # curl https://bootstrap.pypa.io/get-pip.py --output get-pip.py
    # python -m pip install scipy

    #system(stringr::str_glue("sh SICER1.1/SICER/SICER.sh {input_dir} {input} {control} {output_dir} {genome} {format(r_threshold, scientific=F)} {format(window_size, scientific=F)} {format(fragment_size, scientific=F)} {format(fraction, scientific=F)} {format(gap_size, scientific=F)} {format(fdr, scientific=F)}",
    #                 input_dir=dirname(breaks_control_filtered_bed),
    #                 input=basename(breaks_filtered_bed),
    #                 control=basename(breaks_control_filtered_bed),
    #                 output_dir="data/sicer1",
    #                 genome="mm9",
    #                 fraction=0.74,
    #                 r_threshold=5,
    #                 window_size=sicer_params["window"],
    #                 fragment_size=1,
    #                 gap_size=sicer_params["gap"],
    #                 fdr=sicer_params["fdr"]
    #))
    #
    #sicer_diff_bed.cols = cols(peak_chrom=col_character(), peak_start=col_double(), peak_end=col_double(), peak_sample_breaks_count=col_double(), peak_control_breaks_count=col_double(), peak_pvalue=col_double(), peak_fc=col_double(), peak_fdr=col_double())
    #sicer_diff_output = readr::read_tsv(sicer_diff_bed, col_names=names(sicer_diff_bed.cols$cols), col_types=sicer_diff_bed.cols) %>%
    #  dplyr::mutate(bait_chrom=x.bait_chrom)
    #sicer_diff_results = rbind(sicer_diff_results, sicer_diff_output)



    system(stringr::str_glue("sh SICER1.1/SICER/SICER-rb.sh {input_dir} {input} {output_dir} {genome} {format(r_threshold, scientific=F)} {format(window_size, scientific=F)} {format(fragment_size, scientific=F)} {format(fraction, scientific=F)} {format(gap_size, scientific=F)} {format(e_value, scientific=F)}",
                     input_dir=dirname(breaks_control_filtered_bed),
                     input=basename(breaks_control_filtered_bed),
                     output_dir="data/sicer1",
                     genome="mm9",
                     fraction=0.74,
                     r_threshold=5,
                     window_size=sicer_params["window"],
                     fragment_size=1,
                     gap_size=sicer_params["gap"],
                     e_value=sicer_params["e_value"]
    ))
    sicer_control_output = readr::read_tsv(sicer_control_bed, col_names=names(sicer_control_results.cols$cols), col_types=sicer_control_results.cols) %>%
      dplyr::mutate(bait_chrom=x.bait_chrom)
    readr::write_tsv(sicer_control_output %>% dplyr::mutate(peak_id=paste0("peak", 1:n())) %>% dplyr::select(peak_chrom, peak_start, peak_end, peak_id, peak_score), col_names=F, file=paste0(sicer_control_bed, ".bed"))
    sicer_control_results = rbind(sicer_control_results, sicer_control_output)
    #
    # # Working for Chrom3 and Chrom1 (but needs tuning everytime)
    # system(stringr::str_glue("sh SICER1.1/SICER/SICER-rb.sh {input_dir} {input} {output_dir} {genome} {format(r_threshold, scientific=F)} {format(window_size, scientific=F)} {fragment_size} {fraction} {format(gap_size, scientific=F)} {format(e_value, scientific=F)}",
    #                  input_dir=dirname(breaks_filtered_bed),
    #                  input=basename(breaks_control_filtered_bed),
    #                  output_dir="data/sicer1",
    #                  genome="mm9",
    #                  fraction=0.74,
    #                  r_threshold=100000,
    #                  window_size=1000,
    #                  fragment_size=20,
    #                  gap_size=3000,
    #                  e_value=1e-15
    # ))
    #
    #
    # # Equivalent with SICER2, but e-value in sicer 2 has to be integer (BUG!)
    # system(stringr::str_glue("venv/bin/sicer -s mm9 --redundancy_threshold 100000 --window_size 1000 --fragment_size 20 -egf 0.74 --gap_size 3000 --e_value 1 --cpu 24  -t {input} -o {output_dir}", input=breaks_control_filtered_bed, output_dir="data/sicer2"))
    # system(stringr::str_glue("sh SICER1.1/SICER/SICER-rb.sh {input_dir} {input} {output_dir} {genome} {format(r_threshold, scientific=F)} {format(window_size, scientific=F)} {fragment_size} {fraction} {format(gap_size, scientific=F)} {format(e_value, scientific=F)}",
    #                  input_dir=dirname(breaks_filtered_bed),
    #                  input=basename(breaks_control_filtered_bed),
    #                  output_dir="data/sicer1",
    #                  genome="mm9",
    #                  fraction=0.74,
    #                  r_threshold=100000,
    #                  window_size=1000,
    #                  fragment_size=20,
    #                  gap_size=3000,
    #                  e_value=1
    # ))

    # # Species=mm9, redundancy threshold=5, window size=30000, fragment size=1, effective genome fraction=0.74, gap size=90000, FDR=0.01
    # system(stringr::str_glue("epic2 --genome mm9 --bin-size 10000 -kd --fragment-size 20 --gaps-allowed 3 -e 1 -t {input} -o {output}", input=breaks_control_filtered_bed, output=epic2_control_bed))
    system(stringr::str_glue("epic2 --genome mm9 --bin-size 10000 -kd --fragment-size 20 --gaps-allowed 3 -e 1 -t {input} -o {output}", input=breaks_control_filtered_bed, output=epic2_control_bed))
    epic2_control_output = readr::read_tsv(epic2_control_bed, col_names=names(epic2_control_results.cols$cols), col_types=epic2_control_results.cols, skip=1) %>%
      dplyr::mutate(bait_chrom=x.bait_chrom)
    epic2_control_results = rbind(epic2_control_results, epic2_control_output)


    #cmd_enrichment_control = stringr::str_glue(
    #  "macs2 callpeak -t {input} -f BED -g mm --keep-dup all -n {output} --outdir {outdir} --nomodel --extsize 2000 -q 0.01 --llocal 10000000",
    #  input=breaks_control_bed, outdir=dirname(output_control_enrichment), output=basename(output_control_enrichment))
    #system(cmd_enrichment_control)
    #control_enrichment_df = readr::read_tsv(paste0(output_control_enrichment, "_peaks.xls"), col_names=names(peaks_coltypes$cols), comment="#", skip=22, col_types=peaks_coltypes) %>%
    #  dplyr::mutate(bait_chrom=bait_chrom)
    #results = dplyr::bind_rows(results, control_enrichment_df)
    #
    #
    # macs2_input_df = readr::read_tsv(macs2_input, col_names=macs2_cols) %>%
    #   dplyr::mutate(bait_chrom=bait_chrom) %>%
    #   dplyr::inner_join(breaksites_df, by="bait_chrom") %>%
    #   dplyr::filter(!dplyr::between(junction_start, bait_position - bait_region_size/2, bait_position + bait_region_size/2)) %>%
    #   dplyr::select(dplyr::starts_with("junction_")) %>%
    #   dplyr::distinct(junction_chrom, junction_start, .keep_all=T)
    #
    # readr::write_tsv(macs2_input_df[macs2_cols], file=macs2_input_filtered, col_names=F)
    #
    # macs2_cmd = stringr::str_glue(
    #   "macs2 callpeak -t {input} -f BED -g mm --keep-dup all -n {output} --outdir {outdir} --nomodel --extsize 2000 -q 0.01 --llocal 10000000",
    #   input=macs2_input_filtered, outdir=macs2_outdir, output=macs2_sample)
    # system(macs2_cmd)
    #
    # readr::read_tsv(paste0(macs2_outdir, "/", macs2_sample, "_summits.bed"), col_names=F)
}

x = sicer_control_results %>%
    dplyr::filter(peak_score>100) %>%
    dplyr::left_join(offtargets_df %>% dplyr::filter(offtarget_mismatches>0), by=c("bait_chrom", "peak_chrom"="offtarget_chrom")) %>%
    dplyr::mutate(peak_is_predicted=!is.na(offtarget_start) & peak_start>=offtarget_start & offtarget_start<=peak_end) %>%
    dplyr::arrange(dplyr::desc(peak_is_predicted)) %>%
    dplyr::distinct(bait_chrom, peak_chrom, peak_start, peak_end, .keep_all=T)
table(predicted=x$peak_is_predicted)
table(predicted=x$peak_is_predicted, mismatches=x$offtarget_mismatches, x$bait_chrom==x$peak_chrom)

x = epic2_control_results %>%
    # dplyr::filter(peak_score>200) %>%
    dplyr::left_join(offtargets_df %>% dplyr::filter(offtarget_mismatches>0), by=c("bait_chrom", "peak_chrom"="offtarget_chrom")) %>%
    dplyr::mutate(peak_is_predicted=!is.na(offtarget_start) & peak_start<=offtarget_start & offtarget_start<=peak_end) %>%
    dplyr::arrange(dplyr::desc(peak_is_predicted)) %>%
    dplyr::distinct(bait_chrom, peak_chrom, peak_start, peak_end, .keep_all=T)
table(predicted=x$peak_is_predicted)

y = offtargets_df %>%
  dplyr::filter(offtarget_mismatches>0) %>%
  dplyr::left_join(sicer_control_results %>% dplyr::filter(peak_score>100), by=c("bait_chrom", "offtarget_chrom"="peak_chrom")) %>%
  dplyr::mutate(offtarget_is_confirmed=!is.na(peak_score) & peak_start<=offtarget_start & offtarget_start<=peak_end) %>%
  dplyr::arrange(dplyr::desc(offtarget_is_confirmed)) %>%
  dplyr::distinct(bait_chrom, offtarget_chrom, offtarget_start, offtarget_start, .keep_all=T)
table(confirmed=y$offtarget_is_confirmed)


pred <- ROCR::prediction(x$peak_score, x$peak_is_predicted)
perf <- ROCR::performance(pred, measure = "prec", x.measure="rec")
perf <- ROCR::performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf, colorize = TRUE); abline(a=0, b=1)

x = sicer_diff_results %>%
  dplyr::filter(peak_fc > 1)
sicer_diff_results %>%
    dplyr::filter(peak_fc > 1) %>%
    ggplot() +
      geom_point(aes(x=peak_fc, y=-log10(peak_pvalue), color=peak_fdr<0.001))

ggplot(x) +
  geom_boxplot(aes(x=peak_is_predicted, y=log10(peak_score)))
  # geom_violin(aes(x=peak_is_predicted, y=log10(peak_score)))
ggplot(y) +
  geom_violin(aes(x=offtarget_is_confirmed, y=offtarget_mismatches))




offtargets_ranges = with(offtargets_df, IRanges::IRanges(start=offtarget_start, end=offtarget_end))
control_junctions_ranges = with(results,  IRanges::IRanges(start=control_peak_start, end=control_peak_end))
validated_offtargets_overlapping = validated_offtargets %>%
dplyr::inner_join(as.data.frame(IRanges::mergeByOverlaps(offtargets_ranges, control_junctions_ranges)), by=c("offtarget_start"="control_junctions_ranges.start", "offtarget_end"="validated_offtargets_ranges.end")) %>%
dplyr::mutate(offtarget_validated=T) %>%
dplyr::select(bait_name, bait_chrom, offtarget_is_primary, offtarget_chrom, offtarget_strand, offtarget_start, offtarget_end, offtarget_sequence, offtarget_start, offtarget_end, primer_alignment_start=primers_targets_ranges.start, primer_alignment_end=primers_targets_ranges.end, offtarget_validated)






for(breaks_bed in list.files("data/breaks", pattern="*DMSO.bed", full.names=T)) {
    bait_chrom = gsub("(chr)([^_]+).*$", "\\L\\1\\2", basename(breaks_bed), perl=T, ignore.case=T)


    macs2_input_filtered = paste0("tmp/", basename(gsub("\\.bed$", "_filtered.bed", macs2_input)))

    writeLines(paste0("Processing ", bait_chrom, "..."))

    macs2_input_df = readr::read_tsv(macs2_input, col_names=macs2_cols) %>%
      dplyr::mutate(bait_chrom=bait_chrom) %>%
      dplyr::inner_join(breaksites_df, by="bait_chrom") %>%
      dplyr::filter(!dplyr::between(junction_start, bait_position - bait_region_size/2, bait_position + bait_region_size/2)) %>%
      dplyr::select(dplyr::starts_with("junction_")) %>%
      dplyr::distinct(junction_chrom, junction_start, .keep_all=T)

    readr::write_tsv(macs2_input_df[macs2_cols], file=macs2_input_filtered, col_names=F)

    macs2_cmd = stringr::str_glue(
      "macs2 callpeak -t {input} -f BED -g mm --keep-dup all -n {output} --outdir {outdir} --nomodel --extsize 2000 -q 0.01 --llocal 10000000",
      input=macs2_input_filtered, outdir=macs2_outdir, output=macs2_sample)
    system(macs2_cmd)

    readr::read_tsv(paste0(macs2_outdir, "/", macs2_sample, "_summits.bed"), col_names=F)
}