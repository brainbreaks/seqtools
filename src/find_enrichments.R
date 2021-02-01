library(readr)
library(dplyr)
library(dbscan)
library(stringr)
library(ggplot2)
library(ROCR)
library(IRanges)
library(Gviz)

breaksDensity = function(z, from=1, to=NULL, width=10000, bw=1e3) {
    if(is.null(to)) to = max(z$end)
    z.n = round((to - from)/width)
    z.d = density(z$start, from=from, to=to, bw=bw, n=z.n)
    d.bin = z.d$x[2L] - z.d$x[1L]
    d.count = zoo::rollsum(z.d$y, 2)*d.bin
    d.count = d.count*nrow(z)/2
    data.frame(density_chrom=z$chr[1], density_start=z.d$x[-length(z.d$x)], density_end=z.d$x[-1], density_value=d.count)
}

read.breaks = function(dir) {
  breaks = data.frame()
  breaks_cols = cols(break_chrom=col_character(), break_start=col_double(), break_end=col_double(), break_name=col_character(), break_score=col_character(), break_strand=col_character())
  for(breaks_path in list.files(dir, pattern="*.bed", full.names=T)) {
      chr = gsub("(chr)([^_]+).*$", "\\L\\1\\2", basename(breaks_path), perl=T, ignore.case=T)
      condition = ifelse(grepl("DMSO", breaks_path), "Control", "Sample")
      sample = gsub("\\.bed", "", basename(breaks_path))


      x = readr::read_tsv(breaks_path, col_names=names(breaks_cols$cols), col_types=breaks_cols) %>%
        dplyr::mutate(break_name_copy=break_name) %>%
        dplyr::mutate(bait_chrom=chr, experimental_condition=condition, break_sample=sample, break_end=break_start+1) %>%
        tidyr::extract(break_name_copy, "break_replicate", "M01407:320:000000000-ARNUY:1:([^:]).*")
      breaks = rbind(breaks, x)
  }
  breaks = breaks %>% dplyr::mutate(break_id=1:n())

  breaks
}

control2sample_correlation = function() {
  breaks = data.frame()
  breaks_cols = cols(break_chrom=col_character(), break_start=col_double(), break_end=col_double(), break_name=col_character(), break_score=col_character(), break_strand=col_character())
  for(breaks_path in list.files("data/breaks", pattern="*.bed", full.names=T)) {
      chr = gsub("(chr)([^_]+).*$", "\\L\\1\\2", basename(breaks_path), perl=T, ignore.case=T)
      experimental_condition = ifelse(grepl("DMSO", breaks_path), "Control", "Sample")
      sample = gsub("\\.bed", "", basename(breaks_path))


      x = readr::read_tsv(breaks_path, col_names=names(breaks_cols$cols), col_types=breaks_cols) %>%
        dplyr::mutate(bait_chrom=chr, experimental_condition=condition, break_sample=sample, break_end=break_start+1)
      breaks = rbind(breaks, x)
  }
  breaks = breaks %>% dplyr::mutate(break_id=1:n())
  breaks = breaks %>% dplyr::filter(break_chrom==bait_chrom)


  offtargets_df = readr::read_tsv("data/offtargets_predicted.tsv") %>%
    dplyr::mutate(offtarget_id=1:n()) %>%
    dplyr::filter(offtarget_mismatches==0)

  breaks_ranges = with(breaks, IRanges::IRanges(start=break_start, end=break_end, break_id=break_id, break_chrom=break_chrom, break_bait_chrom=bait_chrom))
  offtargets_ranges = with(offtargets_df, IRanges::IRanges(start=offtarget_start-3e6, end=offtarget_end+3e6, offtarget_id=offtarget_id, offtarget_chrom=offtarget_chrom, offtarget_bait_chrom=bait_chrom))
  offtarget2break.map = as.data.frame(IRanges::mergeByOverlaps(offtargets_ranges, breaks_ranges)) %>%
    dplyr::filter(break_chrom==offtarget_chrom) %>%
    dplyr::select(offtarget_id, break_id)


  breaks.densities %>%
    dplyr::filter(dplyr::between(break_start, 57209692 - 3e6, 57209692+3e6) & bait_chrom=="chr13" & break_chrom=="chr13" )
  breaks.densities %>% dplyr::filter(bait_chrom=="chr13" & density_chrom=="chr13") %>% dplyr::arrange(dplyr::desc(Sample)) %>% dplyr::slice(1)

  breaks.densities = breaks %>%
    dplyr::filter(!(break_chrom %in% c("chrM", "chrX", "chrY"))) %>%
    dplyr::anti_join(offtarget2break.map, by="break_id") %>%
    dplyr::group_by(bait_chrom, break_chrom) %>%
    dplyr::mutate(break_pos_max=max(c(break_start, break_end))) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(bait_chrom, experimental_condition, break_chrom) %>%
    dplyr::do((function(z){
      #z = breaks.densities %>% dplyr::filter(experimental_condition=="Sample" & bait_chrom=="chr13" & density_chrom=="chr13")
      zz<<-z
      z = z %>% dplyr::mutate(chr=break_chrom, start=break_start, end=break_end)
      d = breaksDensity(z, to=z$break_pos_max[1], width=2e5) %>%
        dplyr::mutate(bait_chrom=z$bait_chrom[1], experimental_condition=z$experimental_condition[1])
      d
    })(.)) %>%
    reshape2::dcast(bait_chrom + density_chrom + density_start ~ experimental_condition, value.var="density_value") %>%
    dplyr::filter(Control>0 | Sample>0)
  breaks.densities.R = breaks.densities %>%
    dplyr::group_by(density_chrom) %>%
    dplyr::summarise(R=cor(Control, Sample, use="pairwise.complete.obs"), x=max(Control, na.rm=T)*0.8, y=max(Sample, na.rm=T)*0.8)

  png("reports/correlation_between_control_and_sample_all.png", width=2048, height=2048)
  ggplot(breaks.densities) +
    geom_point(aes(x=Control, y=Sample, color=bait_chrom), alpha=0.1) +
    geom_text(aes(x=x, y=y, label=round(R, 2)), data=breaks.densities.R) +
    facet_wrap(~density_chrom, scales="free")
  dev.off()

}

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


find_clusters = function() {
  bait_region_size = 2e6

  offtargets_df = readr::read_tsv("data/offtargets_predicted.tsv") %>% # Doesn't have position for Chr6
    dplyr::mutate(offtarget_id=1:n()) %>%
    dplyr::filter(offtarget_mismatches==0)

  breaks = read.breaks("data/breaks")
      # %>% dplyr::filter(break_chrom==bait_chrom)

  breaks_ranges = with(breaks,  IRanges::IRanges(start=break_start, end=break_end, break_id=break_id, break_bait_chrom=bait_chrom, break_chrom=break_chrom))
  offtargets_ranges = with(offtargets_df, IRanges::IRanges(start=offtarget_start-bait_region_size, end=offtarget_end+bait_region_size, offtarget_id=offtarget_id, offtarget_bait_chrom=bait_chrom, offtarget_chrom=offtarget_chrom))
  breaks2offtargets = as.data.frame(IRanges::mergeByOverlaps(offtargets_ranges, breaks_ranges)) %>%
    dplyr::filter(offtarget_bait_chrom==break_bait_chrom & break_chrom==offtarget_chrom) %>%
    dplyr::select(break_id, offtarget_id)
  breaks = breaks %>%
    dplyr::mutate(bait_overlap=break_id %in% breaks2offtargets$break_id)

  breaks %>%
    dplyr::group_by(bait_chrom, break_chrom) %>%
    dplyr::do((function(z){
      zz<<-z
      asdasd()
      # TODO: estimate minpoints
      res.optics = dbscan::optics(matrix(z$break_start), minPts=10, eps=1e8)
      plot(res.optics)
      res.dbscan = dbscan::extractDBSCAN(res.optics, eps_cl=1e5)
      plot(res.dbscan)
      res.dbscan = dbscan::extractXi(res.optics, xi=0.2)
      plot(res.dbscan)

      z$cluster = res.dbscan$cluster
      z.cluster = z %>%
        dplyr::filter(cluster>0) %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(cluster_start=min(break_start), cluster_end=max(break_end), cluster_width=cluster_end-cluster_start+1) %>%
        dplyr::arrange(dplyr::desc(cluster_width))
      z_ranges.cluster = IRanges::IRanges(start=z.cluster$cluster_start, end=z.cluster$cluster_end)

      plotRanges(z_ranges.cluster, xlim=c(0, max(z.cluster$cluster_end)))
    })(.))

  breaks1 %>%
    dplyr::group_by(bait_chrom, experimental_condition) %>%
    dplyr::summarise(n=n(), n_bait=sum(bait_overlap), n_nonbait=sum(!bait_overlap)) %>%
    data.frame()

  breaks1 %>%
    dplyr::group_by(break_chrom, bait_chrom, experimental_condition) %>%
    dplyr::summarise(n_bait=sum(bait_overlap), n_nonbait=sum(!bait_overlap))  %>%
    reshape2::melt(measure.vars=c("n_bait", "n_nonbait")) %>%
    ggplot() +
      geom_bar(aes(y=value, x=))





  breaks %>%
    dplyr::group_by(bait_chrom, break_chrom) %>%
    dplyr::do((function(z){
      #zz<<-z
      #asdasd()
      # TODO: estimate minpoints
      res.optics = dbscan::optics(matrix(z$break_start), minPts=10, eps=1e8)
      plot(res.optics)
      res.dbscan = dbscan::extractDBSCAN(res.optics, eps_cl=1e5)
    })(.))


              z$cluster = res.dbscan$cluster[res.dbscan$order]
              zf = z %>% dplyr::filter(cluster>0) %>%
                dplyr::group_by(cluster) %>%
                dplyr::summarise(start=min(start), end=max(end), chr="chr6")
              zf.clusters.control = GenomicRanges::makeGRangesFromDataFrame(zf, keep.extra.columns=T)
              rbind(zf.clusters.control, zf.clusters.sample)

  table(breaks$experimental_condition , breaks$break_replicate)
}









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

    bait_neighbour_distance = 3e6
    x.offtargets = offtargets_df %>% dplyr::filter(offtarget_mismatches==0 & bait_chrom==x.bait_chrom)
    breaks_control_filtered = readr::read_tsv(breaks_control_bed, col_names=names(bed_cols$cols), col_types=bed_cols) %>%
      dplyr::filter(!(start >= x.offtargets$offtarget_start-bait_neighbour_distance & end<=x.offtargets$offtarget_start+bait_neighbour_distance & x.offtargets$offtarget_chrom==chr)) %>%
      dplyr::mutate(start=ifelse(strand=="-", start-20, start), end=ifelse(strand=="+", end+20, end)) %>%
      dplyr::mutate(score=1)
    readr::write_tsv(breaks_control_filtered, file=breaks_control_filtered_bed, col_names=F)

    breaks_filtered = readr::read_tsv(breaks_bed, col_names=names(bed_cols$cols), col_types=bed_cols) %>%
      dplyr::filter(!(start >= x.offtargets$offtarget_start-bait_neighbour_distance & end<=x.offtargets$offtarget_start+bait_neighbour_distance & x.offtargets$offtarget_chrom==chr)) %>%
      dplyr::mutate(start=ifelse(strand=="-", start-20, start), end=ifelse(strand=="+", end+20, end)) %>%
      dplyr::mutate(score=1)
    readr::write_tsv(breaks_filtered, file=breaks_filtered_bed, col_names=F)

    dplyr::bind_rows(
      breaks_control_filtered %>% dplyr::mutate(experimental_condition="Control"),
      breaks_filtered %>% dplyr::mutate(experimental_condition="Sample")) %>%
        dplyr::group_by(chr) %>%
        dplyr::summarise(
          control=sum(experimental_condition=="Control"),
          control_filtered=sum(experimental_condition=="Control" & !(start >= x.offtargets$offtarget_start-3e6 & end<=x.offtargets$offtarget_start+3e6 & x.offtargets$offtarget_chrom==chr)),
          sample=sum(experimental_condition=="Sample"),
          sample_filtered=sum(experimental_condition=="Sample" & !(start >= x.offtargets$offtarget_start-3e6 & end<=x.offtargets$offtarget_start+3e6 & x.offtargets$offtarget_chrom==chr))) %>%
      View()

  x = readr::read_tsv(breaks_bed, col_names=names(bed_cols$cols), col_types=bed_cols) %>%
    dplyr::mutate(experimental_condition="Sample") %>%
    tidyr::extract(name, "replicate", "M01407:320:000000000-ARNUY:1:([^:]).*")

  table(x$replicate)



table(stringr::str_extract_all("M01407:320:000000000-ARNUY:1:", x$name))
    breaks = dplyr::bind_rows(
      readr::read_tsv(breaks_bed, col_names=names(bed_cols$cols), col_types=bed_cols) %>% dplyr::mutate(experimental_condition="Sample"),
      readr::read_tsv(breaks_control_bed, col_names=names(bed_cols$cols), col_types=bed_cols) %>% dplyr::mutate(experimental_condition="Control")) %>%
      dplyr::mutate(end=start+1)

    breaks$name


    breaks_clusters = breaks %>%
        dplyr::group_by(experimental_condition, chr) %>%
        dplyr::do((function(z){
            if(nrow(z) < 20) return(data.frame())

            #z.extend = 100
            #z.span = GenomicRanges::makeGRangesFromDataFrame(z %>% dplyr::mutate(start=start - ifelse(strand=="-", z.extend, 0), end=end + ifelse(strand=="+", z.extend, 0), id=1:nrow(z)), seqnames.field="chr", keep.extra.columns=T)
            #l = length(z.span)
            #z.dist_pairwise = data.frame(
            #  id1=rep(1:l, each=l),
            #  id2=rep(1:l, times=l),
            #  distance=GenomicRanges::distance(rep(z.span, each=l), rep(z.span, times=l), ignore.strand=T))
            #
            #z.dist_wide = z.dist_pairwise %>%
            #  reshape2::dcast(id1 ~ id2, value.var="distance") %>%
            #  tibble::column_to_rownames("id1")
            #z.dist = as.dist(z.dist_wide)
            #res.optics = dbscan::optics(z.dist, minPts=20, eps=1e8)

            z.max = max(z = breaks %>% dplyr::filter(chr=="chr6") %>% .$start)
            z = breaks %>% dplyr::filter(chr=="chr6" & experimental_condition=="Control")
            res.optics = dbscan::optics(matrix(z$start), minPts=10, eps=1e8)
            res.dbscan = dbscan::extractDBSCAN(res.optics, eps_cl=1e5)

            z$cluster = res.dbscan$cluster[res.dbscan$order]
            zf = z %>% dplyr::filter(cluster>0) %>%
              dplyr::group_by(cluster) %>%
              dplyr::summarise(start=min(start), end=max(end), chr="chr6")
            zf.clusters.control = GenomicRanges::makeGRangesFromDataFrame(zf, keep.extra.columns=T)
            rbind(zf.clusters.control, zf.clusters.sample)

            plotRanges(zf.clusters, xlim=c(0, z.max))


            z.clusters_ranges_reduced = IRanges::reduce(z.clusters_ranges, min.gapwidth=1, drop.empty.ranges=F, with.revmap=T)
            plotRanges(z.clusters_ranges)


          x = z %>%
            dplyr::arrange(start, cluster)

          zoo::rollapply(x$cluster, 2, function(zz) zz[1]==0 | zz[2]==0 | zz[2]>zz[1])
            plot(res.dbscan$cluster[res.dbscan$order])

            plot(res.xi)
            plot(res.optics)

            res.xi = dbscan::extractDBSCAN(res.optics,  eps_cl=1e5)
            res.hdbscan = dbscan::hdbscan(matrix(z$start), minPts=50, gen_hdbscan_tree=T, gen_simplified_tree=F)
            res.jpclust = dbscan::jpclust(matrix(z$start), k=5, kt=2)
            res.ssnclust = dbscan::sNNclust(matrix(z$start), k=5, eps=1e6, minPts=5)
          table(res.ssnclust$cluster)

            plot(res.hdbscan, show_flat=T)
            res.optics = dbscan::optics(matrix(z$start), minPts=50, eps=1e8)
            res.optics$cluster = res.hdbscan$cluster
            res.optics$cluster = extractFOSC(res.hdbscan$hc, minPts=50)$cluster
            res.optics$cluster = res.jpclust$cluster
            res.optics$cluster = res.ssnclust$cluster
            plot(res.optics)

            plot(res.optics)
            plot(res.xi)
            #plot(res.xi$coredist, type="l", ylim=c(0, 100))
            if(is.null(res.xi$clusters_xi)) return(data.frame())

            z.clusters = data.frame(res.xi$clusters_xi) %>%
              dplyr::mutate(start=pmin(start, end), end=pmax(start, end)) %>%
              dplyr::mutate(optics_start=start, optics_end=end, start=z$start[optics_start], end=z$end[optics_end])

            z.clusters_ranges = with(z.clusters, IRanges::IRanges(start=start, end=end, cluster_id=cluster_id))
            z.clusters_ranges_reduced = IRanges::reduce(z.clusters_ranges, min.gapwidth=1, drop.empty.ranges=F, with.revmap=T)
            z.clusters_reduced = as.data.frame(z.clusters_ranges_reduced) %>% dplyr::mutate(reduced_cluster_id=1:n())
            z.clusters_reduced_map = data.frame(reduced_cluster_id=c(rep(z.clusters_reduced$reduced_cluster_id, sapply(z.clusters_reduced$revmap, length)), 0), cluster_id=c(unlist(z.clusters_reduced$revmap), 0))
            z.clusters_reduced_sizes = z.clusters_reduced_map %>%
              dplyr::inner_join(z.clusters, by="cluster_id")  %>%
              dplyr::group_by(reduced_cluster_id) %>%
              dplyr::summarise(optics_start=min(optics_start), optics_end=max(optics_end))

            z.clusters_reduced = z.clusters_reduced %>%
              dplyr::inner_join(z.clusters_reduced_sizes, by="reduced_cluster_id") %>%
              dplyr::mutate(experimental_condition=z$experimental_condition[1], chr=z$chr[1], cluster_size=optics_end-optics_start+1) %>%
              dplyr::select(-revmap)

          axisTrack = GenomeAxisTrack()
          density_range = GenomicRanges::makeGRangesFromDataFrame(breaksDensity(z), start.field="density_start", end.field="density_end", seqnames.field="density_chrom", keep.extra.columns=T)
          breaksDensityTrack = Gviz::DataTrack(range=density_range, data=breaksDensity(z)$density_value, type="l", name="Sample\njunctions\ndensity")
          breaksClusterTrack = Gviz::AnnotationTrack(range=GenomicRanges::makeGRangesFromDataFrame(z.clusters %>% dplyr::mutate(chromosome=z$chr[1]), seqnames.field="chromosome", start.field="start", end.field="end"), name="Cluster", col="blue", stacking = "dense", alpha=0.5, alpha.title=1, stackHeight=1, shape="box", collapse=F, shape="box", col.line="transparent")
          plotTracks(c(axisTrack, breaksDensityTrack, breaksClusterTrack), from=1, to=max(z$end))


            z.clusters_reduced

            # Debug merged
            res.xi$cluster = data.frame(cluster_id=as.numeric(res.xi$cluster)) %>% dplyr::inner_join(z.clusters_reduced_map, by="cluster_id") %>% .$reduced_cluster_id
            res.xi$clusters_xi = z.clusters_reduced %>%
              dplyr::select(start=optics_start, end=optics_end, cluster_id=reduced_cluster_id) %>%
              data.frame()
            plot(res.xi)

            dbscan::hullplot(matrix(z$start), res.xi)

            plotRanges(z.clusters_ranges)
            plotRanges(IRanges::reduce(z.clusters_ranges, min.gapwidth=1000, drop.empty.ranges=F, with.revmap=T))

            res.xi = dbscan::extractXi(res.optics, xi = 0.05)
            plot(res.xi)

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