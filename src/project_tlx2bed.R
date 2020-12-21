library(readr)
library(dplyr)
library(IRanges)
library(tidyr)
library(ggplot2)

chromosomes_map_df = readr::read_tsv("data/mm9_chromosomes_synonyms.tsv")
chromosomes_map = chromosomes_map_df$chrom
names(chromosomes_map) = chromosomes_map_df$chrom_synonym


# Read offtargets
offtargets_df = readr::read_tsv("data/offtargets_predicted.tsv") %>%
  dplyr::mutate(offtarget_id=1:n()) %>%
  dplyr::mutate(bait_chrom=chromosomes_map[tolower(bait_chrom)], offtarget_chrom=chromosomes_map[tolower(offtarget_chrom)]) %>%
  dplyr::filter(offtarget_mismatches==0) %>%
  dplyr::mutate(
    offtarget_start=ifelse(bait_chrom=="chr15", 61818880, offtarget_start),
    offtarget_end=ifelse(bait_chrom=="chr15", 61818881, offtarget_end))

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

macs2_cols = cols(
  peak_chrom=col_character(), peak_start=col_double(), peak_end=col_double(), peak_length=col_double(), peak_abs_summit=col_double(),
  peak_pileup=col_double(), peak_pvalue=col_double(), peak_fc=col_double(), peak_score=col_double(), peak_name=col_character()
)

bed_cols = cols(chrom=col_character(), start=col_double(), end=col_double(), name=col_character(), score=col_character(), strand=col_character())

# Rename Chr14 to Alt289 (Was Alt289 and only PW255 was Alt287)
junctions_ann.excluded = "PW121|PW246|PW247|PW248|PW249|JK096|JK097|JK098|JK099|JK100|JK101"
junctions_ann.excluded = "NO FILTER"
junctions_ann = data.frame(break_file=list.files("data/breaks_tlx", full.names=T, recursive=T, pattern="*.tlx")) %>%
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
fragment_size = 20
junctions_df = junctions_df %>%
  dplyr::mutate(break_id=1:n()) %>%
  dplyr::mutate(
      break_start=ifelse(break_strand=="+", break_start, break_start-fragment_size),
      break_end=ifelse(break_strand=="+", break_end+fragment_size, break_end))


#
# Find breaks close to bait
#
junctions_ranges = with(junctions_df, IRanges::IRanges(start=break_start, end=break_end, break_id=break_id, break_chrom=break_chrom, break_bait_chrom=break_bait_chrom, break_exp_condition=break_exp_condition))
offtargets_ranges = with(offtargets_df %>% dplyr::filter(offtarget_mismatches==0), IRanges::IRanges(start=offtarget_start-3e6, end=offtarget_end+3e6, offtarget_id=offtarget_id, offtarget_chrom=offtarget_chrom, offtarget_bait_chrom=bait_chrom))
offtarget2junctions.map = as.data.frame(IRanges::mergeByOverlaps(offtargets_ranges, junctions_ranges)) %>%
  dplyr::filter(break_chrom==offtarget_chrom & break_bait_chrom==offtarget_bait_chrom) %>%
  dplyr::mutate(is_near_bait=T) %>%
  dplyr::distinct(break_id, is_near_bait)
junctions_df = junctions_df %>%
  dplyr::select(-dplyr::matches("is_near_bait")) %>%
  dplyr::left_join(offtarget2junctions.map, by="break_id") %>%
  dplyr::mutate(is_near_bait=tidyr::replace_na(is_near_bait, F))

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
junctions_df %>%
  dplyr::filter(!is_near_bait & break_bait_chrom != break_chrom) %>%
  dplyr::group_by(break_exp_condition) %>%
  dplyr::do((function(z){
    zz<<-z

    bed = list(
      "all" = z,
      "plus" = z %>% dplyr::filter(break_strand=="+"),
      "minus" = z %>% dplyr::filter(break_strand=="-")
    )
    bed_path = sapply(names(bed), function(s) stringr::str_glue("data/breakensembl/{exp}_{strand}.bed", exp=bed[[s]]$break_exp_condition[1], strand=s))
    o = sapply(names(bed), function(s) readr::write_tsv(bed[[s]] %>% dplyr::select(break_chrom, break_start, break_end, break_name, break_score, break_strand), bed_path[[s]], col_names=F))

    macs2_pileup_cmd = sapply(names(bed), function(s) stringr::str_glue("macs2 pileup -i {input} -f BED --extsize {format(extsize, scientific=F)} -o {output}",
      input=ensemble_bed, output=paste0(bed_path[[s]], ".bdg"), extsize=1e4))
    o = lapply(macs2_pileup_cmd, system)

    macs2_output_path = tempfile(tmpdir="data/macs2")
    macs2_params = list(extsize=20, qvalue=0.01, llocal=3e4, shift=5e4)
    macs2_cmd = stringr::str_glue("macs2 callpeak -t {input} -f BED -g mm --keep-dup all -n {output} -B --outdir {outdir} --nomodel --extsize {extsize} -q {qvalue} --llocal {format(llocal, scientific=F)} {pipe}",
      input=bed_path[["all"]],
      outdir=dirname(macs2_output_path),
      output=basename(macs2_output_path),
      extsize=macs2_params$extsize,
      qvalue=macs2_params$qvalue,
      llocal=macs2_params$llocal,
      pipe="", pipe2="2>/dev/null"
    )
    system(macs2_cmd)

    data.frame()
  })(.))

x = junctions_df %>%
  dplyr::filter(!is_near_bait & break_bait_chrom != break_chrom)
table(x$break_chrom)


    ##################################
    # Create temporary files
    ##################################
    path_bed = tempfile()
    readr::write_tsv(z %>% dplyr::select(break_chrom, break_start, break_end, break_name, break_score, break_strand), path_bed, col_names=F)

    path_bed_far = tempfile()
    readr::write_tsv(z %>% dplyr::filter(!is_near_bait) %>% dplyr::select(break_chrom, break_start, break_end, break_name, break_score, break_strand), path_bed_far, col_names=F)

junctions_df.f = junctions_df %>% dplyr::filter(!grepl("chrX|chrY", break_bait_chrom))
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
