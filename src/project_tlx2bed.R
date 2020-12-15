
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
junctions_ann = data.frame(break_file=list.files("data/breaks_tlx", full.names=T, recursive=T, pattern="*.tlx")) %>%
  dplyr::filter(!grepl(junctions_ann.excluded, break_file)) %>%
   dplyr::mutate(
     break_bait_chrom = tolower(basename(dirname(break_file))),
     break_condition = ifelse(grepl("DMSO", break_file), "Control", "Sample"),
     break_tech_sample = gsub("_.*", "", basename(break_file)),
     break_bio_sample=gsub(".*_((Alt|R)[0-9]+).*", "\\1", basename(break_file), ignore.case=T, perl=T)
   ) %>%
  dplyr::group_by(break_bait_chrom, break_condition) %>%
  dplyr::mutate(break_replicate=1:n()) %>%
  dplyr::ungroup()

junctions = data.frame()
for(i in 1:nrow(junctions_ann)) {
  tlx_ann = junctions_ann[i,,drop=F]
  tlx = cbind(tlx_ann, readr::read_tsv(tlx_ann$break_file, col_types=tlx_cols))
  bed = tlx %>%
    dplyr::mutate(break_score=".", break_strand=ifelse(Strand=="-1", "-", "+")) %>%
    dplyr::mutate(break_exp_condition=tlx.condition, break_file=tlx_path, break_bait_chrom=tlx.bait_chrom) %>%
    dplyr::select(break_bait_chrom, break_exp_condition, break_chrom=Rname, break_start=Junction, break_end=Junction, break_name=Qname, break_score, break_strand, break_file)
  junctions = rbind(junctions, bed)
}



  junctions %>%
    dplyr::group_by(break_exp_condition, bait_chrom, break_chrom) %>%
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

dim(junctions)

junctions = junctions %>%
  dplyr::group_by()