# ============================================================
# Multi-organ causal network: hub-and-spoke layout
#
# Design:
#   Each mortality clock is shown as a separate hub.
#   Disease endpoints are arranged around each clock as satellites.
#
# Layers:
#   1. LDSC genetic correlation:
#        analysis_group == "Disease_endpoint" & P < 0.05 / 527
#        dotted line, no arrow
#
#   2. MR Clock2DE:
#        solid arrow Clock -> Disease
#
#   3. MR DE2Clock:
#        solid arrow Disease -> Clock
#
# Edge color:
#   Red  = positive association, rg > 0 or OR > 1
#   Blue = negative association, rg < 0 or OR < 1
#
# Important:
#   MR edge endpoints are shortened so arrowheads are not hidden under nodes.
#   Bidirectional MR pairs are plotted as two curved arrows in opposite directions.
#   Disease labels are abbreviated but not truncated with "...".
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(scales)
  library(grid)
})

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------
result_dir <- "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/Result"
out_dir <- result_dir

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ldsc_candidates <- c(
  file.path(result_dir, "LDSC_gc_mortality_clocks_all_targets.tsv"),
  file.path(result_dir, "LDSC_gc_mortality_clocks_all_targets(3).tsv"),
  "/mnt/data/LDSC_gc_mortality_clocks_all_targets(3).tsv"
)

clock2de_candidates <- c(
  file.path(result_dir, "2SampleMR_MortalityClock2DE_ivw_sig_by_DE.tsv"),
  file.path(result_dir, "2SampleMR_MortalityClock2DE_ivw_sig_by_DE(1).tsv"),
  "/mnt/data/2SampleMR_MortalityClock2DE_ivw_sig_by_DE(1).tsv"
)

de2clock_candidates <- c(
  file.path(result_dir, "2SampleMR_DE2clock_ivw_sig_by_DE.tsv"),
  file.path(result_dir, "2SampleMR_DE2clock_ivw_sig_by_DE(1).tsv"),
  "/mnt/data/2SampleMR_DE2clock_ivw_sig_by_DE(1).tsv"
)

pick_existing_file <- function(candidates, pattern = NULL) {
  hits <- candidates[file.exists(candidates)]
  
  if (length(hits) > 0) {
    return(hits[1])
  }
  
  if (!is.null(pattern)) {
    pattern_hits <- Sys.glob(pattern)
    if (length(pattern_hits) > 0) {
      return(sort(pattern_hits, decreasing = TRUE)[1])
    }
  }
  
  stop(
    "Could not find any candidate file:\n",
    paste(candidates, collapse = "\n")
  )
}

ldsc_file <- pick_existing_file(
  ldsc_candidates,
  file.path(result_dir, "LDSC_gc_mortality_clocks_all_targets*.tsv")
)

clock2de_file <- pick_existing_file(
  clock2de_candidates,
  file.path(result_dir, "2SampleMR_MortalityClock2DE_ivw_sig_by_DE*.tsv")
)

de2clock_file <- pick_existing_file(
  de2clock_candidates,
  file.path(result_dir, "2SampleMR_DE2clock_ivw_sig_by_DE*.tsv")
)

message("Using files:")
message("  LDSC      : ", ldsc_file)
message("  Clock2DE  : ", clock2de_file)
message("  DE2Clock  : ", de2clock_file)

# ------------------------------------------------------------
# 2. Plot controls
# ------------------------------------------------------------
ldsc_p_threshold <- 0.05 / 527

# Set these to finite values such as 5 or 10 to reduce density.
top_n_ldsc_per_clock <- Inf
top_n_clock2de_per_clock <- Inf
top_n_de2clock_per_clock <- Inf

# Association values are now encoded by color only.
show_edge_labels <- FALSE

# Disease labels are abbreviated and wrapped, but not truncated with "...".
show_disease_labels <- TRUE
disease_label_wrap_width <- 18

# Layout controls
n_cluster_cols <- 4
cluster_margin_x <- c(0.11, 0.89)
cluster_margin_y <- c(0.92, 0.08)

# Satellite rings
max_single_ring_nodes <- 10
max_two_ring_nodes <- 24

# Edge endpoint offsets prevent arrowheads from being hidden by nodes.
hub_edge_offset <- 0.015
disease_edge_offset <- 0.010
ldsc_edge_offset_hub <- 0.012
ldsc_edge_offset_disease <- 0.008

# Figure size
fig_width <- 18
fig_height <- 22
fig_dpi <- 500

# y-radius adjustment makes radial clusters look closer to circles
aspect_y_adjust <- fig_width / fig_height

out_prefix <- file.path(out_dir, "multi_organ_causal_network_radial_clusters_arrow_fixed")

# ------------------------------------------------------------
# 3. General helper functions
# ------------------------------------------------------------
clean_code <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_replace_all("\\s+", "_") %>%
    str_to_upper()
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

get_first_col <- function(tbl, candidates, default = NA_character_) {
  hits <- intersect(candidates, names(tbl))
  
  if (length(hits) == 0) {
    return(rep(default, nrow(tbl)))
  }
  
  tbl[[hits[1]]]
}

rescale_safe <- function(x, to = c(0, 1)) {
  x <- as.numeric(x)
  
  if (all(is.na(x))) {
    return(rep(mean(to), length(x)))
  }
  
  rng <- range(x, na.rm = TRUE)
  
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) {
    return(rep(mean(to), length(x)))
  }
  
  scales::rescale(x, to = to, from = rng)
}

wrap_label <- function(x, width = 18) {
  map_chr(as.character(x), ~ stringr::str_wrap(.x, width = width))
}

# Shorten disease names using abbreviations, without using "..."
abbreviate_disease_label <- function(x) {
  x <- str_squish(as.character(x))
  
  x <- str_replace_all(x, regex("Alzheimer's disease", ignore_case = TRUE), "AD")
  x <- str_replace_all(x, regex("Autism spectrum disorder", ignore_case = TRUE), "ASD")
  x <- str_replace_all(x, regex("Alcohol use disorder", ignore_case = TRUE), "AUD")
  x <- str_replace_all(x, regex("Bipolar disorder", ignore_case = TRUE), "BIP")
  x <- str_replace_all(x, regex("Schizophrenia", ignore_case = TRUE), "SCZ")
  
  x <- str_replace_all(x, regex("Cardiovascular disease", ignore_case = TRUE), "CVD")
  x <- str_replace_all(x, regex("Coronary heart disease", ignore_case = TRUE), "CHD")
  x <- str_replace_all(x, regex("Ischemic heart disease", ignore_case = TRUE), "IHD")
  x <- str_replace_all(x, regex("Coronary atherosclerosis", ignore_case = TRUE), "Coronary athero.")
  x <- str_replace_all(x, regex("Atrial fibrillation", ignore_case = TRUE), "AF")
  x <- str_replace_all(x, regex("Myocardial infarction", ignore_case = TRUE), "MI")
  x <- str_replace_all(x, regex("Non-ischemic cardiomyopathy", ignore_case = TRUE), "Non-ischemic CMP")
  x <- str_replace_all(x, regex("Antihypertensive use", ignore_case = TRUE), "Anti-HTN use")
  x <- str_replace_all(x, regex("Hypertension", ignore_case = TRUE), "HTN")
  x <- str_replace_all(x, regex("Essential HTN", ignore_case = TRUE), "Essential HTN")
  x <- str_replace_all(x, regex("Coronary revascularization", ignore_case = TRUE), "Coronary revasc.")
  
  x <- str_replace_all(x, regex("Hypercholesterolemia", ignore_case = TRUE), "Hyperchol.")
  x <- str_replace_all(x, regex("Hyperlipidemia", ignore_case = TRUE), "Hyperlipid.")
  x <- str_replace_all(x, regex("Lipoprotein disorder", ignore_case = TRUE), "Lipoprotein dis.")
  x <- str_replace_all(x, regex("Metabolic disorder", ignore_case = TRUE), "Metabolic dis.")
  x <- str_replace_all(x, regex("Type 1 diabetes", ignore_case = TRUE), "T1D")
  x <- str_replace_all(x, regex("Type 2 diabetes", ignore_case = TRUE), "T2D")
  x <- str_replace_all(x, regex("Diabetes with complications", ignore_case = TRUE), "DM w/ complications")
  x <- str_replace_all(x, regex("Diabetic retinopathy", ignore_case = TRUE), "DM retinopathy")
  x <- str_replace_all(x, regex("Diabetic ketoacidosis", ignore_case = TRUE), "DKA")
  x <- str_replace_all(x, regex("Diabetic hypoglycemia", ignore_case = TRUE), "DM hypoglycemia")
  x <- str_replace_all(x, regex("Gestational diabetes", ignore_case = TRUE), "Gestational DM")
  x <- str_replace_all(x, regex("Insulin-treated diabetes", ignore_case = TRUE), "Insulin-treated DM")
  x <- str_replace_all(x, regex("Non-toxic thyroid disorder", ignore_case = TRUE), "Non-toxic thyroid dis.")
  
  x <- str_replace_all(x, regex("Autoimmune disease", ignore_case = TRUE), "Autoimmune dis.")
  x <- str_replace_all(x, regex("Non-thyroid autoimmune disease", ignore_case = TRUE), "Non-thyroid autoimmune dis.")
  x <- str_replace_all(x, regex("Strict non-thyroid autoimmune disease", ignore_case = TRUE), "Strict non-thyroid autoimmune dis.")
  x <- str_replace_all(x, regex("Rheumatoid arthritis", ignore_case = TRUE), "RA")
  x <- str_replace_all(x, regex("Seropositive RA", ignore_case = TRUE), "Seropositive RA")
  x <- str_replace_all(x, regex("Systemic connective tissue disease", ignore_case = TRUE), "Systemic connective tissue dis.")
  x <- str_replace_all(x, regex("Other systemic connective disease", ignore_case = TRUE), "Other systemic connective dis.")
  x <- str_replace_all(x, regex("Papulosquamous disorder", ignore_case = TRUE), "Papulosquamous dis.")
  x <- str_replace_all(x, regex("Oral lichen planus", ignore_case = TRUE), "Oral lichen planus")
  
  x <- str_replace_all(x, regex("Allergic asthma", ignore_case = TRUE), "Allergic asthma")
  x <- str_replace_all(x, regex("Asthma with infections", ignore_case = TRUE), "Asthma w/ infections")
  
  x <- str_replace_all(x, regex("Renal failure", ignore_case = TRUE), "Renal failure")
  x <- str_replace_all(x, regex("Abdominal hernia", ignore_case = TRUE), "Abdominal hernia")
  x <- str_replace_all(x, regex("Intestinal infections", ignore_case = TRUE), "Intestinal infections")
  x <- str_replace_all(x, regex("Crohn's second-line therapy", ignore_case = TRUE), "Crohn's 2nd-line tx")
  x <- str_replace_all(x, regex("Oral disease / reimbursement", ignore_case = TRUE), "Oral disease")
  x <- str_replace_all(x, regex("Any mental disorder", ignore_case = TRUE), "Any mental dis.")
  x <- str_replace_all(x, regex("Endometriosis", ignore_case = TRUE), "Endometriosis")
  
  x <- str_replace_all(x, regex("disorder", ignore_case = TRUE), "dis.")
  x <- str_replace_all(x, regex("disease", ignore_case = TRUE), "dis.")
  x <- str_replace_all(x, regex("complications", ignore_case = TRUE), "comp.")
  x <- str_replace_all(x, regex("treatment", ignore_case = TRUE), "tx")
  x <- str_replace_all(x, regex("therapy", ignore_case = TRUE), "tx")
  x <- str_replace_all(x, regex("reimbursement", ignore_case = TRUE), "reimb.")
  x <- str_replace_all(x, regex("cardiomyopathy", ignore_case = TRUE), "CMP")
  
  str_squish(x)
}

shorten_and_wrap_disease_label <- function(x, width = 18) {
  x %>%
    abbreviate_disease_label() %>%
    wrap_label(width = width)
}

# ------------------------------------------------------------
# 4. Clock helper functions
# ------------------------------------------------------------
parse_clock_organ <- function(clock_id) {
  x <- str_to_lower(as.character(clock_id))
  
  case_when(
    str_detect(x, "^brain") ~ "Brain",
    str_detect(x, "^eye") ~ "Eye",
    str_detect(x, "^pulmonary|^lung") ~ "Pulmonary",
    str_detect(x, "^heart") ~ "Heart",
    str_detect(x, "^hepatic|^liver") ~ "Hepatic",
    str_detect(x, "^renal|^kidney") ~ "Renal",
    str_detect(x, "^pancreas") ~ "Pancreas",
    str_detect(x, "^spleen") ~ "Spleen",
    str_detect(x, "^immune") ~ "Immune",
    str_detect(x, "^endocrine") ~ "Endocrine",
    str_detect(x, "^digestive") ~ "Digestive",
    str_detect(x, "^metabolic") ~ "Metabolic",
    str_detect(x, "^adipose") ~ "Adipose",
    str_detect(x, "^skin") ~ "Skin",
    str_detect(x, "^reproductive") ~ "Reproductive",
    TRUE ~ "Other"
  )
}

parse_clock_display_organ <- function(clock_id) {
  x <- str_to_lower(as.character(clock_id))
  
  case_when(
    str_detect(x, "^reproductive_male") ~ "Reproductive male",
    str_detect(x, "^reproductive_female") ~ "Reproductive female",
    TRUE ~ parse_clock_organ(clock_id)
  )
}

parse_clock_modality <- function(clock_id) {
  x <- str_to_lower(as.character(clock_id))
  
  case_when(
    str_detect(x, "_mri$|_mri_") ~ "MRI",
    str_detect(x, "proteomics") ~ "Proteomics",
    str_detect(x, "metabolomics") ~ "Metabolomics",
    TRUE ~ "Other"
  )
}

pretty_clock_label <- function(clock_id) {
  organ <- parse_clock_display_organ(clock_id)
  modality <- parse_clock_modality(clock_id)
  
  modality_short <- case_when(
    modality == "MRI" ~ "MRI",
    modality == "Proteomics" ~ "Prot",
    modality == "Metabolomics" ~ "Met",
    TRUE ~ modality
  )
  
  paste0(organ, "\n", modality_short)
}

# ------------------------------------------------------------
# 5. Disease manual map
# ------------------------------------------------------------
disease_map <- tribble(
  ~disease_code, ~disease_name, ~disease_class,
  
  # Cardiometabolic / endocrine
  "E4_OBESITYNAS", "Obesity", "Cardiometabolic",
  "E4_HYPERCHOL", "Hypercholesterolemia", "Cardiometabolic",
  "E4_HYPERLIPNAS", "Hyperlipidemia", "Cardiometabolic",
  "E4_LIPOPROT", "Lipoprotein disorder", "Cardiometabolic",
  "E4_METABOLIA", "Metabolic disorder", "Cardiometabolic",
  "E4_DM1NASCOMP", "T1D complications", "Endocrine / Diabetes",
  "E4_DM1OPTH", "T1D eye complications", "Endocrine / Diabetes",
  "E4_DM2NASCOMP", "T2D complications", "Endocrine / Diabetes",
  "E4_NONTOXIC_THYROID", "Non-toxic thyroid disorder", "Endocrine / Diabetes",
  "DM_RETINOPATHY_EXMORE", "Diabetic retinopathy", "Endocrine / Diabetes",
  "DM_HYPOGLYC", "Diabetic hypoglycemia", "Endocrine / Diabetes",
  "DM_KETOACIDOSIS", "Diabetic ketoacidosis", "Endocrine / Diabetes",
  "DM_SEVERAL_COMPLICATIONS", "Diabetes with complications", "Endocrine / Diabetes",
  "GEST_DIABETES", "Gestational diabetes", "Endocrine / Diabetes",
  "KELA_DIAB_INSUL_EXMORE", "Insulin-treated diabetes", "Endocrine / Diabetes",
  "T1D_WIDE", "Type 1 diabetes", "Endocrine / Diabetes",
  "T2D", "Type 2 diabetes", "Endocrine / Diabetes",
  "THYROTOXICOSIS", "Thyrotoxicosis", "Endocrine / Diabetes",
  
  # Cardiovascular
  "I9_CHD", "Coronary heart disease", "Cardiovascular",
  "I9_CORATHER", "Coronary atherosclerosis", "Cardiovascular",
  "I9_IHD", "Ischemic heart disease", "Cardiovascular",
  "I9_HYPTENS", "Hypertension", "Cardiovascular",
  "I9_HYPTENSESS", "Essential hypertension", "Cardiovascular",
  "I9_AF_REIMB", "Atrial fibrillation", "Cardiovascular",
  "I9_NONISCHCARDMYOP", "Non-ischemic cardiomyopathy", "Cardiovascular",
  "I9_MI_STRICT", "Myocardial infarction", "Cardiovascular",
  "I9_ANGIO", "Angina", "Cardiovascular",
  "I9_REVASC", "Coronary revascularization", "Cardiovascular",
  "FG_CVD", "Cardiovascular disease", "Cardiovascular",
  "FG_DOAAC", "DOAC use", "Cardiovascular",
  "RX_ANTIHYP", "Antihypertensive use", "Cardiovascular",
  "RX_STATIN", "Statin use", "Cardiovascular",
  
  # Renal / eye
  "N14_RENFAIL", "Renal failure", "Renal",
  "H7_RETINOPATHYDIAB", "Diabetic retinopathy", "Eye",
  
  # Immune / inflammatory
  "AUTOIMMUNE", "Autoimmune disease", "Immune / Inflammatory",
  "AUTOIMMUNE_NONTHYROID", "Non-thyroid autoimmune disease", "Immune / Inflammatory",
  "AUTOIMMUNE_NONTHYROID_STRICT", "Strict non-thyroid autoimmune disease", "Immune / Inflammatory",
  "M13_GOUT", "Gout", "Immune / Inflammatory",
  "M13_RHEUMA", "Rheumatoid arthritis", "Immune / Inflammatory",
  "M13_POLYARTHROPATHIES", "Polyarthropathies", "Immune / Inflammatory",
  "RHEUMA_SEROPOS_WIDE", "Seropositive RA", "Immune / Inflammatory",
  "RHEUMA_SEROPOS_OTH", "Other seropositive RA", "Immune / Inflammatory",
  "M13_SYSTCONNECT", "Systemic connective tissue disease", "Immune / Inflammatory",
  "OTHER_SYSTCON_FG", "Other systemic connective disease", "Immune / Inflammatory",
  "L12_DERMATITISNAS", "Dermatitis", "Immune / Inflammatory",
  "L12_PAPULOSQUAMOUS", "Papulosquamous disorder", "Immune / Inflammatory",
  "K11_ORAL_LICHEN_PLANUS_WIDE", "Oral lichen planus", "Immune / Inflammatory",
  
  # Respiratory / allergy
  "ALLERG_ASTHMA", "Allergic asthma", "Respiratory / Allergy",
  "ASTHMA_INFECTIONS", "Asthma with infections", "Respiratory / Allergy",
  
  # Neuropsychiatric
  "AD", "Alzheimer's disease", "Neuropsychiatric",
  "ADHD", "ADHD", "Neuropsychiatric",
  "ASD", "Autism spectrum disorder", "Neuropsychiatric",
  "AUD", "Alcohol use disorder", "Neuropsychiatric",
  "BIP", "Bipolar disorder", "Neuropsychiatric",
  "SCZ", "Schizophrenia", "Neuropsychiatric",
  "G6_AD_WIDE", "Alzheimer's disease", "Neuropsychiatric",
  "KRA_PSY_ANYMENTAL", "Any mental disorder", "Neuropsychiatric",
  
  # GI / infection / medication
  "ABDOM_HERNIA", "Abdominal hernia", "Digestive / GI",
  "AB1_INTESTINAL_INFECTIONS", "Intestinal infections", "Digestive / GI",
  "RX_CROHN_2NDLINE", "Crohn's second-line therapy", "Digestive / GI",
  "K11_REIMB_202", "Oral disease / reimbursement", "Digestive / GI",
  
  # Reproductive
  "N14_ENDOMETRIOSIS", "Endometriosis", "Reproductive"
) %>%
  mutate(disease_code = clean_code(disease_code)) %>%
  distinct(disease_code, .keep_all = TRUE)

infer_disease_class_one <- function(code, source = NA_character_) {
  code <- clean_code(code)
  
  hit <- disease_map %>%
    filter(disease_code == code)
  
  if (nrow(hit) > 0) {
    return(hit$disease_class[1])
  }
  
  case_when(
    str_detect(code, "^RX_") ~ "Medication",
    str_detect(code, "^I9_|^FG_CVD|^FG_DOAAC") ~ "Cardiovascular",
    str_detect(code, "^E4_|^DM_|^T1D|^T2D|GEST_DIABETES|THYROID|THYRO") ~ "Endocrine / Diabetes",
    str_detect(code, "AUTOIMMUNE|RHEUMA|GOUT|DERMAT|LICHEN|PAPULOS|SYSTCON") ~ "Immune / Inflammatory",
    str_detect(code, "ASTHMA|ALLERG") ~ "Respiratory / Allergy",
    str_detect(code, "RENAL|RENFAIL") ~ "Renal",
    str_detect(code, "RETINOPATHY|^H7_") ~ "Eye",
    code %in% c("AD", "ADHD", "ASD", "AUD", "BIP", "SCZ") |
      str_detect(code, "PSY|MENTAL|G6_AD") ~ "Neuropsychiatric",
    str_detect(code, "HERNIA|CROHN|INTESTINAL|^K11_") ~ "Digestive / GI",
    str_detect(code, "ENDOMETRIOSIS") ~ "Reproductive",
    !is.na(source) & source %in% c("PGC", "FinnGen") ~ as.character(source),
    TRUE ~ "Other"
  )
}

infer_disease_name_one <- function(code, fallback = NA_character_) {
  code <- clean_code(code)
  
  hit <- disease_map %>%
    filter(disease_code == code)
  
  if (nrow(hit) > 0) {
    return(hit$disease_name[1])
  }
  
  fallback <- as.character(fallback)
  fallback <- str_squish(fallback)
  
  if (!is.na(fallback) &&
      fallback != "" &&
      clean_code(fallback) != code &&
      !str_detect(str_to_lower(fallback), "^na$|^unknown$|unavailable")) {
    return(str_squish(str_replace_all(fallback, "_", " ")))
  }
  
  code %>%
    str_replace_all("_", " ") %>%
    str_to_lower() %>%
    str_to_title()
}

# ------------------------------------------------------------
# 6. Read files
# ------------------------------------------------------------
ldsc_raw <- readr::read_tsv(ldsc_file, show_col_types = FALSE)
clock2de_raw <- readr::read_tsv(clock2de_file, show_col_types = FALSE)
de2clock_raw <- readr::read_tsv(de2clock_file, show_col_types = FALSE)

# ------------------------------------------------------------
# 7. Standardize LDSC edges
# ------------------------------------------------------------
ldsc_p_vec <- safe_num(get_first_col(ldsc_raw, c("P", "p", "pval", "p_value")))
ldsc_rg_vec <- safe_num(get_first_col(ldsc_raw, c("gc_mean", "rg", "rg_mean", "rg_estimate", "rg_value")))

ldsc_clock_vec <- get_first_col(
  ldsc_raw,
  c("mortality_clock", "mortality_clock_from_file", "clock_id", "clock_folder")
)

ldsc_target_vec <- get_first_col(
  ldsc_raw,
  c("target_id", "target", "outcome", "disease_code")
)

ldsc_target_name_vec <- get_first_col(
  ldsc_raw,
  c("target_display", "target_name", "trait_name", "phenotype", "description")
)

ldsc_source_vec <- get_first_col(
  ldsc_raw,
  c("target_source", "source")
)

ldsc_analysis_group_vec <- get_first_col(
  ldsc_raw,
  c("analysis_group")
)

ldsc_edges <- ldsc_raw %>%
  mutate(
    P_num = ldsc_p_vec,
    rg = ldsc_rg_vec,
    clock_id = clean_code(ldsc_clock_vec),
    disease_code = clean_code(ldsc_target_vec),
    disease_name_raw = as.character(ldsc_target_name_vec),
    target_source = as.character(ldsc_source_vec),
    analysis_group = as.character(ldsc_analysis_group_vec)
  ) %>%
  filter(
    analysis_group == "Disease_endpoint",
    !is.na(P_num),
    !is.na(rg),
    P_num < ldsc_p_threshold
  ) %>%
  mutate(
    disease_name = map2_chr(disease_code, disease_name_raw, infer_disease_name_one),
    disease_class = map2_chr(disease_code, target_source, infer_disease_class_one),
    relation_type = "LDSC",
    from_id = clock_id,
    to_id = disease_code,
    effect_value = rg,
    effect_label = sprintf("rg=%+.2f", rg),
    p_value = P_num,
    assoc_sign = if_else(effect_value >= 0, "Positive", "Negative"),
    edge_strength = pmax(-log10(pmax(p_value, 1e-300)), 1)
  ) %>%
  select(
    relation_type,
    from_id,
    to_id,
    clock_id,
    disease_code,
    disease_name,
    disease_class,
    effect_value,
    effect_label,
    p_value,
    assoc_sign,
    edge_strength,
    target_source
  )

if (is.finite(top_n_ldsc_per_clock)) {
  ldsc_edges <- ldsc_edges %>%
    group_by(clock_id) %>%
    arrange(desc(abs(effect_value)), p_value, .by_group = TRUE) %>%
    slice_head(n = top_n_ldsc_per_clock) %>%
    ungroup()
}

# ------------------------------------------------------------
# 8. Standardize MR Clock -> Disease edges
# ------------------------------------------------------------
clock2de_or_vec <- safe_num(get_first_col(clock2de_raw, c("or", "OR", "odds_ratio")))
clock2de_b_vec <- safe_num(get_first_col(clock2de_raw, c("b", "beta")))
clock2de_p_vec <- safe_num(get_first_col(clock2de_raw, c("pval", "P", "p", "p_value")))

clock2de_clock_vec <- get_first_col(
  clock2de_raw,
  c("clock_id", "clock_folder", "exposure")
)

clock2de_disease_vec <- get_first_col(
  clock2de_raw,
  c("outcome_code", "outcome", "disease_code")
)

clock2de_source_vec <- get_first_col(
  clock2de_raw,
  c("target_source", "source")
)

clock2de_edges <- clock2de_raw %>%
  mutate(
    p_value = clock2de_p_vec,
    beta = clock2de_b_vec,
    or_raw = clock2de_or_vec,
    or_value = if_else(!is.na(or_raw), or_raw, exp(beta)),
    clock_id = clean_code(clock2de_clock_vec),
    disease_code = clean_code(clock2de_disease_vec),
    target_source = as.character(clock2de_source_vec)
  ) %>%
  filter(!is.na(p_value), !is.na(or_value), !is.na(clock_id), !is.na(disease_code)) %>%
  mutate(
    disease_name = map2_chr(disease_code, disease_code, infer_disease_name_one),
    disease_class = map2_chr(disease_code, target_source, infer_disease_class_one),
    relation_type = "Clock_to_Disease",
    from_id = clock_id,
    to_id = disease_code,
    effect_value = or_value,
    effect_label = sprintf("OR=%.2f", effect_value),
    p_value = p_value,
    assoc_sign = if_else(effect_value >= 1, "Positive", "Negative"),
    edge_strength = pmax(-log10(pmax(p_value, 1e-300)), 1)
  ) %>%
  select(
    relation_type,
    from_id,
    to_id,
    clock_id,
    disease_code,
    disease_name,
    disease_class,
    effect_value,
    effect_label,
    p_value,
    assoc_sign,
    edge_strength,
    target_source
  )

if (is.finite(top_n_clock2de_per_clock)) {
  clock2de_edges <- clock2de_edges %>%
    group_by(clock_id) %>%
    arrange(desc(abs(log(effect_value))), p_value, .by_group = TRUE) %>%
    slice_head(n = top_n_clock2de_per_clock) %>%
    ungroup()
}

# ------------------------------------------------------------
# 9. Standardize MR Disease -> Clock edges
# ------------------------------------------------------------
de2clock_or_vec <- safe_num(get_first_col(de2clock_raw, c("or", "OR", "odds_ratio")))
de2clock_b_vec <- safe_num(get_first_col(de2clock_raw, c("b", "beta")))
de2clock_p_vec <- safe_num(get_first_col(de2clock_raw, c("pval", "P", "p", "p_value")))

de2clock_clock_vec <- get_first_col(
  de2clock_raw,
  c("clock_id", "clock_folder", "outcome")
)

de2clock_disease_vec <- get_first_col(
  de2clock_raw,
  c("disease_code", "exposure", "exposure_original")
)

de2clock_source_vec <- get_first_col(
  de2clock_raw,
  c("target_source", "source")
)

de2clock_edges <- de2clock_raw %>%
  mutate(
    p_value = de2clock_p_vec,
    beta = de2clock_b_vec,
    or_raw = de2clock_or_vec,
    or_value = if_else(!is.na(or_raw), or_raw, exp(beta)),
    clock_id = clean_code(de2clock_clock_vec),
    disease_code = clean_code(de2clock_disease_vec),
    target_source = as.character(de2clock_source_vec)
  ) %>%
  filter(!is.na(p_value), !is.na(or_value), !is.na(clock_id), !is.na(disease_code)) %>%
  mutate(
    disease_name = map2_chr(disease_code, disease_code, infer_disease_name_one),
    disease_class = map2_chr(disease_code, target_source, infer_disease_class_one),
    relation_type = "Disease_to_Clock",
    from_id = disease_code,
    to_id = clock_id,
    effect_value = or_value,
    effect_label = sprintf("OR=%.2f", effect_value),
    p_value = p_value,
    assoc_sign = if_else(effect_value >= 1, "Positive", "Negative"),
    edge_strength = pmax(-log10(pmax(p_value, 1e-300)), 1)
  ) %>%
  select(
    relation_type,
    from_id,
    to_id,
    clock_id,
    disease_code,
    disease_name,
    disease_class,
    effect_value,
    effect_label,
    p_value,
    assoc_sign,
    edge_strength,
    target_source
  )

if (is.finite(top_n_de2clock_per_clock)) {
  de2clock_edges <- de2clock_edges %>%
    group_by(clock_id) %>%
    arrange(desc(abs(log(effect_value))), p_value, .by_group = TRUE) %>%
    slice_head(n = top_n_de2clock_per_clock) %>%
    ungroup()
}

# ------------------------------------------------------------
# 10. Combine all edge layers
# ------------------------------------------------------------
all_edges <- bind_rows(
  ldsc_edges,
  clock2de_edges,
  de2clock_edges
) %>%
  mutate(
    clock_organ = parse_clock_organ(clock_id),
    clock_display_organ = parse_clock_display_organ(clock_id),
    clock_modality = parse_clock_modality(clock_id),
    clock_label = pretty_clock_label(clock_id),
    pair_id = paste(clock_id, disease_code, sep = "__")
  )

message("Edge counts:")
print(all_edges %>% count(relation_type, name = "n_edges"))
message("LDSC threshold: ", signif(ldsc_p_threshold, 3))
message("Total edges: ", nrow(all_edges))

if (nrow(all_edges) == 0) {
  stop("No edges were available after filtering.")
}

# ------------------------------------------------------------
# 11. Create radial cluster layout
# ------------------------------------------------------------
organ_order <- c(
  "Brain", "Eye", "Pulmonary", "Heart", "Hepatic", "Renal", "Pancreas",
  "Spleen", "Immune", "Endocrine", "Digestive", "Metabolic",
  "Adipose", "Skin", "Reproductive", "Other"
)

modality_order <- c("MRI", "Proteomics", "Metabolomics", "Other")

clock_clusters <- all_edges %>%
  distinct(clock_id, clock_label, clock_organ, clock_display_organ, clock_modality) %>%
  mutate(
    clock_organ_factor = factor(clock_organ, levels = organ_order),
    clock_modality_factor = factor(clock_modality, levels = modality_order),
    n_edges_clock = map_int(clock_id, ~ sum(all_edges$clock_id == .x)),
    n_diseases_clock = map_int(clock_id, ~ n_distinct(all_edges$disease_code[all_edges$clock_id == .x]))
  ) %>%
  arrange(clock_organ_factor, clock_modality_factor, clock_label) %>%
  mutate(
    cluster_index = row_number(),
    cluster_col = ((cluster_index - 1) %% n_cluster_cols) + 1,
    cluster_row = floor((cluster_index - 1) / n_cluster_cols) + 1
  )

n_cluster_rows <- max(clock_clusters$cluster_row)

x_centers <- seq(cluster_margin_x[1], cluster_margin_x[2], length.out = n_cluster_cols)
y_centers <- seq(cluster_margin_y[1], cluster_margin_y[2], length.out = n_cluster_rows)

clock_clusters <- clock_clusters %>%
  mutate(
    cx = x_centers[cluster_col],
    cy = y_centers[cluster_row],
    hub_size = 5.1 + 0.35 * sqrt(n_edges_clock)
  )

pair_nodes <- all_edges %>%
  group_by(clock_id, disease_code) %>%
  summarise(
    disease_name = first(disease_name),
    disease_class = first(disease_class),
    n_edges_pair = n(),
    max_edge_strength = max(edge_strength, na.rm = TRUE),
    has_ldsc = any(relation_type == "LDSC"),
    has_clock_to_disease = any(relation_type == "Clock_to_Disease"),
    has_disease_to_clock = any(relation_type == "Disease_to_Clock"),
    has_bidirectional_mr = has_clock_to_disease & has_disease_to_clock,
    .groups = "drop"
  ) %>%
  left_join(
    clock_clusters %>%
      select(clock_id, cx, cy, clock_organ, clock_display_organ, clock_modality, clock_label),
    by = "clock_id"
  ) %>%
  group_by(clock_id) %>%
  arrange(desc(has_bidirectional_mr), disease_class, desc(n_edges_pair), disease_name, .by_group = TRUE) %>%
  mutate(
    node_index = row_number(),
    n_nodes_clock = n()
  ) %>%
  ungroup()

assign_radial_positions <- function(tbl) {
  n <- nrow(tbl)
  
  if (n == 0) {
    return(tbl)
  }
  
  if (n <= max_single_ring_nodes) {
    ring_capacity <- c(n)
    ring_radius <- c(0.060)
  } else if (n <= max_two_ring_nodes) {
    inner_n <- min(10, n)
    outer_n <- n - inner_n
    ring_capacity <- c(inner_n, outer_n)
    ring_radius <- c(0.050, 0.083)
  } else {
    inner_n <- min(10, n)
    middle_n <- min(16, n - inner_n)
    outer_n <- n - inner_n - middle_n
    ring_capacity <- c(inner_n, middle_n, outer_n)
    ring_radius <- c(0.046, 0.075, 0.105)
  }
  
  ring_id <- rep(seq_along(ring_capacity), ring_capacity)
  pos_in_ring <- unlist(map(ring_capacity, seq_len))
  
  out <- tbl %>%
    mutate(
      ring_id = ring_id[seq_len(n)],
      pos_in_ring = pos_in_ring[seq_len(n)]
    )
  
  out <- out %>%
    group_by(ring_id) %>%
    mutate(
      n_in_ring = n(),
      angle = pi / 2 - 2 * pi * (row_number() - 1) / n_in_ring,
      radius = ring_radius[ring_id[1]]
    ) %>%
    ungroup() %>%
    mutate(
      disease_x = cx + radius * cos(angle),
      disease_y = cy + radius * aspect_y_adjust * sin(angle),
      disease_node_size = 2.05 + 0.22 * sqrt(n_edges_pair),
      disease_label = shorten_and_wrap_disease_label(disease_name, width = disease_label_wrap_width),
      label_x = disease_x + if_else(cos(angle) >= 0, 0.011, -0.011),
      label_y = disease_y + 0.004 * sign(sin(angle)),
      label_hjust = if_else(cos(angle) >= 0, 0, 1)
    )
  
  out
}

pair_nodes_layout <- pair_nodes %>%
  group_by(clock_id) %>%
  group_modify(~ assign_radial_positions(.x)) %>%
  ungroup()

# Hub nodes
clock_nodes <- clock_clusters %>%
  transmute(
    node_id = clock_id,
    node_type = "Mortality clock",
    clock_id,
    label = clock_label,
    x = cx,
    y = cy,
    clock_organ,
    clock_display_organ,
    clock_modality,
    node_size = hub_size,
    n_edges_clock,
    n_diseases_clock
  )

# Disease satellite nodes are intentionally duplicated per clock
disease_nodes <- pair_nodes_layout %>%
  transmute(
    node_id = paste(clock_id, disease_code, sep = "__"),
    node_type = "Disease endpoint",
    clock_id,
    disease_code,
    label = disease_label,
    disease_name,
    disease_class,
    x = disease_x,
    y = disease_y,
    label_x,
    label_y,
    label_hjust,
    node_size = disease_node_size,
    n_edges_pair,
    has_ldsc,
    has_clock_to_disease,
    has_disease_to_clock,
    has_bidirectional_mr
  )

# ------------------------------------------------------------
# 12. Attach coordinates to edge table
# ------------------------------------------------------------
edge_xy <- all_edges %>%
  left_join(
    clock_nodes %>%
      select(clock_id, clock_x = x, clock_y = y),
    by = "clock_id"
  ) %>%
  left_join(
    disease_nodes %>%
      select(clock_id, disease_code, disease_x = x, disease_y = y),
    by = c("clock_id", "disease_code")
  ) %>%
  filter(
    !is.na(clock_x),
    !is.na(clock_y),
    !is.na(disease_x),
    !is.na(disease_y)
  ) %>%
  mutate(
    dx_clock_to_disease = disease_x - clock_x,
    dy_clock_to_disease = disease_y - clock_y,
    edge_len = sqrt(dx_clock_to_disease^2 + dy_clock_to_disease^2),
    edge_len = if_else(edge_len == 0 | is.na(edge_len), 1e-6, edge_len),
    ux = dx_clock_to_disease / edge_len,
    uy = dy_clock_to_disease / edge_len,
    
    edge_width = rescale_safe(edge_strength, to = c(0.38, 1.38)),
    edge_color = case_when(
      assoc_sign == "Positive" ~ "#C53A3A",
      assoc_sign == "Negative" ~ "#2C6FB7",
      TRUE ~ "#777777"
    )
  )

# LDSC: dotted, no arrow, clock-disease shortened segment
edges_ldsc <- edge_xy %>%
  filter(relation_type == "LDSC") %>%
  mutate(
    x = clock_x + ux * ldsc_edge_offset_hub,
    y = clock_y + uy * ldsc_edge_offset_hub,
    xend = disease_x - ux * ldsc_edge_offset_disease,
    yend = disease_y - uy * ldsc_edge_offset_disease,
    label_x = x + 0.50 * (xend - x),
    label_y = y + 0.50 * (yend - y),
    curve = 0.00
  )

# Clock -> Disease MR: solid arrow with visible arrowhead near disease node
edges_clock2de <- edge_xy %>%
  filter(relation_type == "Clock_to_Disease") %>%
  mutate(
    x = clock_x + ux * hub_edge_offset,
    y = clock_y + uy * hub_edge_offset,
    xend = disease_x - ux * disease_edge_offset,
    yend = disease_y - uy * disease_edge_offset,
    label_x = x + 0.55 * (xend - x),
    label_y = y + 0.55 * (yend - y),
    curve = 0.18
  )

# Disease -> Clock MR: solid arrow with visible arrowhead near clock node
edges_de2clock <- edge_xy %>%
  filter(relation_type == "Disease_to_Clock") %>%
  mutate(
    x = disease_x - ux * disease_edge_offset,
    y = disease_y - uy * disease_edge_offset,
    xend = clock_x + ux * hub_edge_offset,
    yend = clock_y + uy * hub_edge_offset,
    label_x = x + 0.45 * (xend - x),
    label_y = y + 0.45 * (yend - y),
    curve = -0.18
  )

edge_labels <- bind_rows(edges_ldsc, edges_clock2de, edges_de2clock) %>%
  group_by(relation_type) %>%
  arrange(desc(edge_strength), .by_group = TRUE) %>%
  ungroup()

# ------------------------------------------------------------
# 13. Colors
# ------------------------------------------------------------
organ_palette <- c(
  "Brain" = "#4169A1",
  "Eye" = "#74A9CF",
  "Pulmonary" = "#5AA05A",
  "Heart" = "#2E6B45",
  "Hepatic" = "#B08D00",
  "Renal" = "#3B5BA9",
  "Pancreas" = "#8E5BAE",
  "Spleen" = "#6A3D9A",
  "Immune" = "#1B9E77",
  "Endocrine" = "#9BAA4F",
  "Digestive" = "#A65E2E",
  "Metabolic" = "#C17D11",
  "Adipose" = "#E6AB02",
  "Skin" = "#C77CFF",
  "Reproductive" = "#E7298A",
  "Other" = "#999999"
)

edge_sign_palette <- c(
  "Positive" = "#C53A3A",
  "Negative" = "#2C6FB7"
)

# ------------------------------------------------------------
# 14. Subtitle / caption
# ------------------------------------------------------------
n_bidirectional_pairs <- disease_nodes %>%
  summarise(n = sum(has_bidirectional_mr, na.rm = TRUE)) %>%
  pull(n)

subtitle_text <- paste0(
  "Each mortality clock is shown as a separate hub. ",
  "LDSC dotted no-arrow edges: P < 0.05/527 = ",
  signif(ldsc_p_threshold, 3),
  "; MR solid arrows show causal direction. ",
  "LDSC n=", nrow(edges_ldsc),
  ", Clock-to-Disease MR n=", nrow(edges_clock2de),
  ", Disease-to-Clock MR n=", nrow(edges_de2clock),
  ", bidirectional MR pairs n=", n_bidirectional_pairs,
  "."
)

caption_text <- paste0(
  "Red edges indicate positive associations: rg > 0 for LDSC or OR > 1 for MR. ",
  "Blue edges indicate negative associations: rg < 0 for LDSC or OR < 1 for MR. ",
  "Dotted edges are non-directional LDSC genetic correlations. ",
  "Solid arrows are MR directions; bidirectional causal relationships are shown as two opposite curved arrows."
)

# ------------------------------------------------------------
# 15. Plot
# ------------------------------------------------------------
p_net <- ggplot() +
  
  # Light cluster background circles
  geom_point(
    data = clock_clusters,
    aes(x = cx, y = cy),
    shape = 21,
    size = 44,
    fill = "#F8F8F8",
    color = "#E5E5E5",
    stroke = 0.35,
    alpha = 0.85
  ) +
  
  # LDSC dotted, no arrow
  geom_curve(
    data = edges_ldsc,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      color = assoc_sign,
      linewidth = edge_width
    ),
    curvature = 0.00,
    linetype = "22",
    alpha = 0.70,
    lineend = "round"
  ) +
  
  # Clock -> Disease MR arrows
  geom_curve(
    data = edges_clock2de,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      color = assoc_sign,
      linewidth = edge_width
    ),
    curvature = 0.18,
    alpha = 0.82,
    arrow = arrow(length = unit(0.13, "inches"), type = "closed"),
    lineend = "round"
  ) +
  
  # Disease -> Clock MR arrows
  geom_curve(
    data = edges_de2clock,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      color = assoc_sign,
      linewidth = edge_width
    ),
    curvature = -0.18,
    alpha = 0.82,
    arrow = arrow(length = unit(0.13, "inches"), type = "closed"),
    lineend = "round"
  ) +
  
  # Edge labels are off by default
  {
    if (show_edge_labels && nrow(edge_labels) > 0) {
      geom_text(
        data = edge_labels,
        aes(
          x = label_x,
          y = label_y,
          label = effect_label,
          color = assoc_sign
        ),
        size = 2.0,
        alpha = 0.95,
        check_overlap = TRUE
      )
    }
  } +
  
  # Disease satellite nodes
  geom_point(
    data = disease_nodes,
    aes(
      x = x,
      y = y,
      size = node_size
    ),
    shape = 22,
    fill = "white",
    color = "grey20",
    stroke = 0.38,
    alpha = 0.98
  ) +
  
  # Mortality clock hubs
  geom_point(
    data = clock_nodes,
    aes(
      x = x,
      y = y,
      fill = clock_organ,
      size = node_size
    ),
    shape = 21,
    color = "grey10",
    stroke = 0.55,
    alpha = 0.98
  ) +
  
  # Clock labels
  geom_text(
    data = clock_nodes,
    aes(
      x = x,
      y = y + 0.044,
      label = label
    ),
    size = 3.0,
    lineheight = 0.88,
    fontface = "bold",
    color = "grey15"
  ) +
  
  # Disease endpoint labels
  {
    if (show_disease_labels) {
      geom_text(
        data = disease_nodes,
        aes(
          x = label_x,
          y = label_y,
          label = label,
          hjust = label_hjust
        ),
        size = 1.85,
        lineheight = 0.88,
        color = "grey18",
        check_overlap = TRUE
      )
    }
  } +
  
  # Manual legend: LDSC dotted
  annotate(
    "segment",
    x = 0.035,
    xend = 0.085,
    y = 0.035,
    yend = 0.035,
    color = "grey35",
    linewidth = 0.75,
    linetype = "22"
  ) +
  annotate(
    "text",
    x = 0.095,
    y = 0.035,
    label = "LDSC rg",
    hjust = 0,
    vjust = 0.5,
    size = 3.0,
    color = "grey25"
  ) +
  
  # Manual legend: MR arrow
  geom_segment(
    data = tibble(
      x = 0.035,
      xend = 0.085,
      y = 0.017,
      yend = 0.017
    ),
    aes(x = x, xend = xend, y = y, yend = yend),
    color = "grey35",
    linewidth = 0.75,
    arrow = arrow(length = unit(0.10, "inches"), type = "closed"),
    inherit.aes = FALSE
  ) +
  annotate(
    "text",
    x = 0.095,
    y = 0.017,
    label = "MR OR",
    hjust = 0,
    vjust = 0.5,
    size = 3.0,
    color = "grey25"
  ) +
  
  # Manual legend: color sign
  annotate(
    "segment",
    x = 0.035,
    xend = 0.085,
    y = 0.065,
    yend = 0.065,
    color = "#C53A3A",
    linewidth = 0.85
  ) +
  annotate(
    "text",
    x = 0.095,
    y = 0.065,
    label = "Positive rg / OR>1",
    hjust = 0,
    vjust = 0.5,
    size = 3.0,
    color = "grey25"
  ) +
  annotate(
    "segment",
    x = 0.035,
    xend = 0.085,
    y = 0.052,
    yend = 0.052,
    color = "#2C6FB7",
    linewidth = 0.85
  ) +
  annotate(
    "text",
    x = 0.095,
    y = 0.052,
    label = "Negative rg / OR<1",
    hjust = 0,
    vjust = 0.5,
    size = 3.0,
    color = "grey25"
  ) +
  
  scale_color_manual(
    values = edge_sign_palette,
    drop = FALSE,
    name = "Association sign"
  ) +
  scale_fill_manual(
    values = organ_palette,
    drop = FALSE,
    name = "Mortality-clock organ"
  ) +
  scale_size_identity() +
  scale_linewidth_identity() +
  
  coord_cartesian(
    xlim = c(-0.02, 1.02),
    ylim = c(0.00, 1.00),
    clip = "off"
  ) +
  
  labs(
    title = "Multi-organ genetic and causal network linking mortality clocks and disease endpoints",
    subtitle = subtitle_text,
    caption = caption_text
  ) +
  
  theme_void(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(
      face = "bold",
      size = 16.5,
      color = "grey10",
      hjust = 0
    ),
    plot.subtitle = element_text(
      size = 10.2,
      color = "grey25",
      hjust = 0,
      lineheight = 1.05
    ),
    plot.caption = element_text(
      size = 9.0,
      color = "grey35",
      hjust = 0,
      lineheight = 1.05
    ),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 9.5),
    legend.text = element_text(size = 8.2),
    plot.margin = margin(14, 16, 12, 16)
  )

print(p_net)

# ------------------------------------------------------------
# 16. Save figure
# ------------------------------------------------------------
ggsave(
  filename = paste0(out_prefix, ".pdf"),
  plot = p_net,
  width = fig_width,
  height = fig_height,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  filename = paste0(out_prefix, ".png"),
  plot = p_net,
  width = fig_width,
  height = fig_height,
  units = "in",
  dpi = fig_dpi,
  bg = "white"
)

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    filename = paste0(out_prefix, ".svg"),
    plot = p_net,
    width = fig_width,
    height = fig_height,
    units = "in",
    device = svglite::svglite,
    bg = "white"
  )
}

# ------------------------------------------------------------
# 17. Save edge/node tables
# ------------------------------------------------------------
readr::write_tsv(
  edge_xy,
  paste0(out_prefix, "_edges.tsv")
)

readr::write_tsv(
  clock_nodes,
  paste0(out_prefix, "_clock_nodes.tsv")
)

readr::write_tsv(
  disease_nodes,
  paste0(out_prefix, "_disease_satellite_nodes.tsv")
)

edge_summary <- edge_xy %>%
  count(
    relation_type,
    assoc_sign,
    clock_organ,
    clock_modality,
    disease_class,
    name = "n_edges",
    sort = TRUE
  )

readr::write_tsv(
  edge_summary,
  paste0(out_prefix, "_edge_summary.tsv")
)

readr::write_tsv(
  edges_ldsc,
  paste0(out_prefix, "_LDSC_edges.tsv")
)

readr::write_tsv(
  edges_clock2de,
  paste0(out_prefix, "_MR_Clock_to_Disease_edges.tsv")
)

readr::write_tsv(
  edges_de2clock,
  paste0(out_prefix, "_MR_Disease_to_Clock_edges.tsv")
)

bidirectional_pairs <- disease_nodes %>%
  filter(has_bidirectional_mr) %>%
  arrange(clock_id, disease_code)

readr::write_tsv(
  bidirectional_pairs,
  paste0(out_prefix, "_bidirectional_MR_pairs.tsv")
)

message("Done.")
message("Saved:")
message("  ", paste0(out_prefix, ".pdf"))
message("  ", paste0(out_prefix, ".png"))
if (requireNamespace("svglite", quietly = TRUE)) {
  message("  ", paste0(out_prefix, ".svg"))
}
message("  ", paste0(out_prefix, "_edges.tsv"))
message("  ", paste0(out_prefix, "_clock_nodes.tsv"))
message("  ", paste0(out_prefix, "_disease_satellite_nodes.tsv"))
message("  ", paste0(out_prefix, "_edge_summary.tsv"))
message("  ", paste0(out_prefix, "_LDSC_edges.tsv"))
message("  ", paste0(out_prefix, "_MR_Clock_to_Disease_edges.tsv"))
message("  ", paste0(out_prefix, "_MR_Disease_to_Clock_edges.tsv"))
message("  ", paste0(out_prefix, "_bidirectional_MR_pairs.tsv"))