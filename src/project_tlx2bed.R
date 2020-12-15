
tlx_cols = cols(
  Qname=col_character(), JuncID=col_character(), Rname=col_character(), Junction=col_double(),
  Strand=col_character(), Rstart=col_double(), Rend=col_double(),
  B_Rname=col_character(), B_Rstart=col_double(), B_Rend=col_double(), B_Strand=col_double(),
  B_Qstart=col_double(), B_Qend=col_double(), Qstart=col_double(), Qend=col_double(), Qlen=col_double(),
  B_Cigar=col_character(), Cigar=col_character(), Seq=col_character(), J_Seq=col_character(), Barcode=col_logical(),
  unaligned=col_double(), baitonly=col_double(), uncut=col_double(), misprimed=col_double(), freqcut=col_double(),
  largegap=col_double(), mapqual=col_double(), breaksite=col_double(), sequential=col_double(), repeatseq=col_double(), duplicate=col_double(), Note=col_character()
)

junctions = data.frame()
junctions_ann = data.frame(breaks_file=list.files("data/breaks_tlx", full.names=T, recursive=T, pattern="*.tlx")) %>%
   dplyr::mutate(
     break_bait_chrom = tolower(basename(dirname(breaks_file))),
     break_condition = ifelse(grepl("DMSO", breaks_file), "Control", "Sample"),
     break_tech_sample = gsub("_.*", "", basename(breaks_file)),
     break_bio_sample=gsub(".*(DMSO|APH)_([^_]+).*", "\\2", basename(breaks_file), ignore.case=T, perl=T) #,
     #break_smth=gsub(".*intersect_(.*)\\.tlx", "\\1", breaks_file, ignore.case=T)
   ) %>%
  dplyr::group_by(break_bait_chrom, break_condition, break_bio_sample) %>%
  dplyr::mutate(break_tech_replicate=match(break_tech_sample, unique(break_tech_sample))) %>%
  dplyr::group_by(break_bait_chrom, break_condition) %>%
  dplyr::mutate(break_bio_replicate=match(break_bio_sample, unique(break_bio_sample))) %>%
  dplyr::ungroup()
table(junctions_ann$break_bio_sample, junctions_ann$break_tech_replicate)
table(bio=junctions_ann$break_bio_replicate, tech=junctions_ann$break_tech_replicate)

junctions_ann%>%
  dplyr::group_by(break_bait_chrom, break_condition) %>%
  dplyr::filter(n() > 4) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(break_bait_chrom, break_condition, break_bio_sample, break_tech_sample) %>%
  View()
for(tlx_path in list.files("data/breaks_tlx", full.names=T, recursive=T, pattern="*.tlx")) {
  print(tlx_path)
  tlx.bait_chrom = tolower(basename(dirname(tlx_path)))
  tlx.condition = ifelse(grepl("DMSO", tlx_path), "Control", "Sample")

  tlx = readr::read_tsv(tlx_path, col_types=tlx_cols)
  bed = tlx %>%
    dplyr::mutate(Score=".", break_strand=ifelse(Strand=="-1", "-", "+")) %>%
    dplyr::mutate(break_exp_condition=tlx.condition, break_file=tlx_path, break_bait_chrom=tlx.bait_chrom) %>%
    dplyr::select(break_bait_chrom, break_exp_condition, break_chrom=Rname, break_start=Junction, break_end=Junction, break_name=Qname, break_score=Score, break_strand, break_file)
  junctions = rbind(junctions, bed)
}

junctions = junctions %>%
  dplyr::group_by()