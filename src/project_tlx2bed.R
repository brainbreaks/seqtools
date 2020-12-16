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

# Rename Chr14 to Alt289 (Was Alt289 and only PW255 was Alt287)
junctions_ann.excluded = "PW121|PW246|PW247|PW248|PW249|JK096|JK097|JK098|JK099|JK100|JK101"
junctions_ann.excluded = ""
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
junctions_df = junctions_df %>% dplyr::mutate(break_id=1:n())


#
# Find breaks close to bait
#
junctions_ranges = with(junctions_df, IRanges::IRanges(start=break_start, end=break_end, break_id=break_id, break_chrom=break_chrom, break_bait_chrom=break_bait_chrom))
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




islands = data.frame()
junctions_df %>%
  dplyr::group_by(break_bait_chrom, break_chrom) %>%
  dplyr::do((function(z){
    zz<<-z
    asdasd()

    # DBSCAN
    dbscan_params = list(minPts=20, eps_cl=1e5)
    res.optics = dbscan::optics(matrix(z$break_start), minPts=dbscan_params$minPts, eps=1e6)
    res.dbscan = dbscan::extractDBSCAN(res.optics, eps_cl=dbscan_params$eps_cl)
    plot(res.dbscan, ylim=c(0,1e5))

    dbscan_output = z %>%
      dplyr::mutate(cluster=res.dbscan$cluster) %>%
      dplyr::filter(cluster>0) %>%
      dplyr::group_by(cluster) %>%
      dplyr::summarise(peak_method="dbscan", bait_chrom=z$break_bait_chrom[1], peak_chrom=z$break_chrom[1], peak_start=min(break_start), peak_end=max(break_end), peak_score=0) %>%
      dplyr::select(-cluster)
    islands = rbind(islands, dbscan_output)

    path_bed = tempfile()
    readr::write_tsv(z %>% dplyr::select(break_chrom, break_start, break_end, break_name, break_score, break_strand), path_bed, col_names=F)

    path_bed_far = tempfile()
    readr::write_tsv(z %>% dplyr::filter(!is_near_bait) %>% dplyr::select(break_chrom, break_start, break_end, break_name, break_score, break_strand), path_bed_far, col_names=F)


    # SICER
    sicer_params = list(window=30000, gap=90000, e_value=0.1, fdr=0.01)
    sicer_control_results.cols = cols("peak_chrom"=col_character(), peak_start=col_double(), peak_end=col_double(), peak_score = col_double())
    system(stringr::str_glue("sh SICER1.1/SICER/SICER-rb.sh {input_dir} {input} {output_dir} {genome} {format(r_threshold, scientific=F)} {format(window_size, scientific=F)} {format(fragment_size, scientific=F)} {format(fraction, scientific=F)} {format(gap_size, scientific=F)} {format(e_value, scientific=F)}",
                     input_dir=dirname(path_bed_far), input=basename(path_bed_far), output_dir="data/sicer1",
                     genome="mm9", fraction=0.74, r_threshold=5, window_size=sicer_params["window"], fragment_size=1, gap_size=sicer_params["gap"], e_value=sicer_params["e_value"]
    ))
    sicer_output_path = stringr::str_glue("data/sicer1/{sample}-W{format(window, scientific=F)}-G{format(gap, scientific=F)}-E{format(e_value, scientific=F)}.scoreisland", window=sicer_params["window"], gap=sicer_params["gap"], e_value=sicer_params["e_value"], sample=gsub("\\.bed$", "", basename(path_bed_far)))
    sicer_output = readr::read_tsv(sicer_output_path, col_names=names(sicer_control_results.cols$cols), col_types=sicer_control_results.cols) %>%
      dplyr::mutate(bait_chrom=z$break_bait_chrom[1], peak_method="sicer") %>%
      dplyr::select(peak_method, bait_chrom, peak_chrom, peak_start, peak_end, peak_score)

    MACS2

    #z_ranges.cluster = IRanges::IRanges(start=z.cluster$cluster_start, end=z.cluster$cluster_end)
    #plotRanges(z_ranges.cluster, xlim=c(0, max(z.cluster$cluster_end)))
  })(.))

sicer =

dim(junctions)

junctions = junctions %>%
  dplyr::group_by()