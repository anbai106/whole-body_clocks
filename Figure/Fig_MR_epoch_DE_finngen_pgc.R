# ============================================================
# Fig. 2E-style bidirectional MR forest plot
# Mortality EPOCH clocks vs FinnGen / PGC disease endpoints
#
# Directions:
#   1) EPOCH -> disease
#   2) Disease -> EPOCH
#
# Input files are already significant after correction by
# the number of disease endpoints.
#
# One row = one significant MR association
# x-axis = MR OR
# horizontal line = 95% CI
# disease endpoint name is annotated after the CI bar
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(grid)
})

# ----------------------------
# 1. Input / output
# ----------------------------
result_dir <- "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/Result"

clock2de_file <- file.path(
  result_dir,
  "2SampleMR_MortalityClock2DE_ivw_sig_by_DE.tsv"
)

de2clock_file <- file.path(
  result_dir,
  "2SampleMR_DE2clock_ivw_sig_by_DE.tsv"
)

if (!file.exists(clock2de_file)) {
  hits <- Sys.glob(file.path(result_dir, "2SampleMR_MortalityClock2DE_ivw_sig_by_DE*.tsv"))
  if (length(hits) == 0) stop("Cannot find Clock2DE MR file.")
  clock2de_file <- sort(hits, decreasing = TRUE)[1]
}

if (!file.exists(de2clock_file)) {
  hits <- Sys.glob(file.path(result_dir, "2SampleMR_DE2clock_ivw_sig_by_DE*.tsv"))
  if (length(hits) == 0) stop("Cannot find DE2Clock MR file.")
  de2clock_file <- sort(hits, decreasing = TRUE)[1]
}

out_prefix <- file.path(
  result_dir,
  "Fig2E_bidirectional_MR_EPOCH_DE_by_direction_forest"
)

message("Using files:")
message("  EPOCH -> disease: ", clock2de_file)
message("  Disease -> EPOCH: ", de2clock_file)

clock2de_raw <- readr::read_tsv(clock2de_file, show_col_types = FALSE)
de2clock_raw <- readr::read_tsv(de2clock_file, show_col_types = FALSE)

# ----------------------------
# 2. Helper functions
# ----------------------------
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

parse_bool <- function(x) {
  if (is.logical(x)) return(x)
  x_chr <- str_to_lower(as.character(x))
  case_when(
    x_chr %in% c("true", "t", "1", "yes", "y") ~ TRUE,
    x_chr %in% c("false", "f", "0", "no", "n") ~ FALSE,
    is.na(x_chr) ~ NA,
    TRUE ~ NA
  )
}

format_pval <- function(p) {
  case_when(
    is.na(p) ~ NA_character_,
    p < 1e-300 ~ "<1e-300",
    p < 1e-3 ~ formatC(p, format = "e", digits = 2),
    TRUE ~ sprintf("%.3f", p)
  )
}

clean_text <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("_", " ") %>%
    str_replace_all("EXALLC|INCLAVO|EXMORE", "") %>%
    str_squish()
}

clean_clock_id <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_replace("^FinnGen::", "") %>%
    str_replace("^PGC::", "") %>%
    str_replace("_clock_acceleration_z$", "") %>%
    str_replace("_mortality_clock$", "") %>%
    str_replace("_mortality$", "") %>%
    str_replace("_clock$", "") %>%
    str_squish()
}

clean_disease_code <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_replace("^FinnGen::", "") %>%
    str_replace("^PGC::", "") %>%
    str_squish()
}

shorten_endpoint <- function(x) {
  x %>%
    str_replace_all("Alzheimer's disease", "AD") %>%
    str_replace_all("Autism spectrum disorder", "ASD") %>%
    str_replace_all("Alcohol use disorder", "AUD") %>%
    str_replace_all("Bipolar disorder", "BIP") %>%
    str_replace_all("Schizophrenia", "SCZ") %>%
    str_replace_all("Attention deficit hyperactivity disorder", "ADHD") %>%
    str_replace_all("Cardiovascular disease", "CVD") %>%
    str_replace_all("Coronary heart disease", "CHD") %>%
    str_replace_all("Ischemic heart disease", "IHD") %>%
    str_replace_all("Coronary atherosclerosis", "Coronary athero.") %>%
    str_replace_all("Atrial fibrillation", "AF") %>%
    str_replace_all("Myocardial infarction", "MI") %>%
    str_replace_all("Hypertension", "HTN") %>%
    str_replace_all("Essential hypertension", "Essential HTN") %>%
    str_replace_all("Antihypertensive use", "Anti-HTN use") %>%
    str_replace_all("Hypercholesterolemia", "Hyperchol.") %>%
    str_replace_all("Hyperlipidemia", "Hyperlipid.") %>%
    str_replace_all("Lipoprotein disorder", "Lipoprotein dis.") %>%
    str_replace_all("Metabolic disorder", "Metabolic dis.") %>%
    str_replace_all("Type 1 diabetes", "T1D") %>%
    str_replace_all("Type 2 diabetes", "T2D") %>%
    str_replace_all("Rheumatoid arthritis", "RA") %>%
    str_replace_all("Autoimmune disease", "Autoimmune dis.") %>%
    str_replace_all("Non-thyroid autoimmune disease", "Non-thyroid autoimmune dis.") %>%
    str_replace_all("Strict non-thyroid autoimmune disease", "Strict non-thyroid autoimmune dis.") %>%
    str_replace_all("Systemic connective tissue disease", "Systemic connective tissue dis.") %>%
    str_replace_all("Other systemic connective disease", "Other systemic connective dis.") %>%
    str_replace_all("Diabetic retinopathy", "DM retinopathy") %>%
    str_replace_all("Diabetes with complications", "DM w/ comp.") %>%
    str_replace_all("Diabetic hypoglycemia", "DM hypoglycemia") %>%
    str_replace_all("Diabetic ketoacidosis", "DKA") %>%
    str_replace_all("Non-toxic thyroid disorder", "Non-toxic thyroid dis.") %>%
    str_replace_all("Papulosquamous disorder", "Papulosquamous dis.") %>%
    str_replace_all("Oral lichen planus", "Oral lichen planus") %>%
    str_replace_all("Abdominal hernia", "Abdominal hernia") %>%
    str_replace_all("Renal failure", "Renal failure") %>%
    str_replace_all("Endometriosis", "Endometriosis") %>%
    str_replace_all("disease", "dis.") %>%
    str_replace_all("disorder", "dis.") %>%
    str_replace_all("complications", "comp.") %>%
    str_replace_all("treatment", "tx") %>%
    str_replace_all("therapy", "tx") %>%
    str_replace_all("reimbursement", "reimb.") %>%
    str_squish()
}

# ----------------------------
# 3. Mortality EPOCH clock labels
# ----------------------------
organ_levels <- c(
  "Brain",
  "Eye",
  "Pulmonary",
  "Heart",
  "Liver / hepatic",
  "Kidney / renal",
  "Pancreas",
  "Spleen",
  "Immune",
  "Endocrine",
  "Digestive",
  "Metabolic",
  "Adipose",
  "Skin",
  "Reproductive"
)

clock_key <- tribble(
  ~clock_id,                          ~clock_block,        ~clock_label,                         ~clock_suborder,
  "brain_mri",                       "Brain",             "Brain MRI",                          1,
  "Brain_proteomics",                "Brain",             "Brain protein",                      2,
  
  "Eye_proteomics",                  "Eye",               "Eye protein",                        1,
  
  "Pulmonary_proteomics",            "Pulmonary",         "Pulmonary protein",                  1,
  
  "heart_mri",                       "Heart",             "Heart MRI",                          1,
  "Heart_proteomics",                "Heart",             "Heart protein",                      2,
  
  "liver_mri",                       "Liver / hepatic",   "Liver MRI",                          1,
  "Hepatic_proteomics",              "Liver / hepatic",   "Hepatic protein",                    2,
  "Hepatic_metabolomics",            "Liver / hepatic",   "Hepatic metabolite",                 3,
  
  "kidney_mri",                      "Kidney / renal",    "Kidney MRI",                         1,
  "Renal_proteomics",                "Kidney / renal",    "Renal protein",                      2,
  
  "pancreas_mri",                    "Pancreas",          "Pancreas MRI",                       1,
  
  "spleen_mri",                      "Spleen",            "Spleen MRI",                         1,
  
  "Immune_proteomics",               "Immune",            "Immune protein",                     1,
  "Immune_metabolomics",             "Immune",            "Immune metabolite",                  2,
  
  "Endocrine_proteomics",            "Endocrine",         "Endocrine protein",                  1,
  "Endocrine_metabolomics",          "Endocrine",         "Endocrine metabolite",               2,
  
  "Digestive_metabolomics",          "Digestive",         "Digestive metabolite",               1,
  
  "Metabolic_metabolomics",          "Metabolic",         "Metabolic metabolite",               1,
  
  "adipose_mri",                     "Adipose",           "Adipose MRI",                        1,
  
  "Skin_proteomics",                 "Skin",              "Skin protein",                       1,
  
  "Reproductive_female_proteomics",  "Reproductive",      "Reproductive female protein",        1,
  "Reproductive_male_proteomics",    "Reproductive",      "Reproductive male protein",          2
) %>%
  mutate(
    clock_block = factor(clock_block, levels = organ_levels),
    clock_modality = case_when(
      str_detect(str_to_lower(clock_id), "mri") ~ "MRI",
      str_detect(str_to_lower(clock_id), "proteomics") ~ "Proteomics",
      str_detect(str_to_lower(clock_id), "metabolomics") ~ "Metabolomics",
      TRUE ~ "Other"
    ),
    clock_id_join = str_to_lower(clock_id)
  ) %>%
  arrange(clock_block, clock_suborder) %>%
  mutate(clock_rank = row_number())

# ----------------------------
# 4. Disease endpoint chapter assignment
# ----------------------------
assign_disease_chapter <- function(target_source, disease_code) {
  
  disease_code <- str_to_upper(as.character(disease_code))
  
  case_when(
    target_source == "PGC" ~ "PGC brain/psychiatric",
    
    str_detect(disease_code, "^(AD_LO|G6|NEURODEG|MIGRAINE|TRAUMBRAIN)") ~ "Neurologic",
    str_detect(disease_code, "^(F5|KRA|AUD|ALCOHOL|ANTIDEPRESSANTS|SLEEP|ADHD|ASD|BIP|SCZ)") ~ "Psychiatric / sleep",
    
    str_detect(disease_code, "^H7") ~ "Eye",
    str_detect(disease_code, "^(H8|FE)") ~ "Ear / sensory",
    
    str_detect(disease_code, "^(I9|FG|CARDIAC|OTHER_SYSTCON|RX_ANTIHYP|RX_STATIN)") ~ "Cardiovascular",
    
    str_detect(disease_code, "^(J10|ASTHMA|COPD|BRONCHITIS|CPAP|NIV|PULM|VOCALCORD|INFLUENZA)") ~ "Respiratory",
    
    str_detect(disease_code, "^(K11|ABDOM|APPEND|DENTAL|TEMPOROMANDIB|RX_CROHN)") ~ "Digestive / oral",
    
    str_detect(disease_code, "^(E4|DM|T1D|T2D|KELA_DIAB|HYPOTHY|THYRO|GOUT|AUTOIMMUNE|M13_GOUT)") ~ "Endocrine / immune",
    
    str_detect(disease_code, "^(D3|BLEEDING)") ~ "Blood",
    
    str_detect(disease_code, "^(C3|C$|CD2)") ~ "Cancer / benign neoplasm",
    
    str_detect(disease_code, "^(L12|ALLERG|NONALLERG|POLLEN|DRY)") ~ "Skin / allergy",
    
    str_detect(disease_code, "^(M13|RHEU|RHEUMA|SPONDYLO|JOINTPAIN|PAIN|PRIM|FALLS|ST19|VWXY20|RX_CODEINE|RX_PARACETAMOL|RX_GLUCO)") ~ "MSK / injury",
    
    str_detect(disease_code, "^(N14|R18|RX_INFERTILITY|Z21)") ~ "GU / reproductive",
    
    str_detect(disease_code, "^(O15|GEST)") ~ "Pregnancy",
    
    str_detect(disease_code, "^(AB1)") ~ "Infectious",
    
    TRUE ~ "Other"
  )
}

chapter_levels_all <- c(
  "PGC brain/psychiatric",
  "Neurologic",
  "Psychiatric / sleep",
  "Eye",
  "Ear / sensory",
  "Cardiovascular",
  "Respiratory",
  "Digestive / oral",
  "Endocrine / immune",
  "Blood",
  "Cancer / benign neoplasm",
  "Skin / allergy",
  "MSK / injury",
  "GU / reproductive",
  "Pregnancy",
  "Infectious",
  "Other"
)

# ----------------------------
# 5. Standardization function
# ----------------------------
standardize_mr <- function(tbl, direction = c("clock2de", "de2clock")) {
  
  direction <- match.arg(direction)
  
  if (direction == "clock2de") {
    direction_label <- "EPOCH -> disease"
    direction_rank <- 1
    
    clock_candidates <- c(
      "clock_id",
      "clock_folder",
      "exposure",
      "exposure_id",
      "id.exposure"
    )
    
    disease_candidates <- c(
      "outcome_code",
      "disease_code",
      "outcome",
      "outcome_id",
      "id.outcome"
    )
    
    disease_name_candidates <- c(
      "outcome_name",
      "outcome_display",
      "disease_name",
      "target_display",
      "phenotype",
      "trait",
      "outcome"
    )
    
  } else {
    direction_label <- "Disease -> EPOCH"
    direction_rank <- 2
    
    clock_candidates <- c(
      "clock_id",
      "clock_folder",
      "outcome",
      "outcome_id",
      "id.outcome"
    )
    
    disease_candidates <- c(
      "disease_code",
      "exposure_original",
      "exposure",
      "exposure_id",
      "id.exposure"
    )
    
    disease_name_candidates <- c(
      "exposure_name",
      "exposure_display",
      "disease_name",
      "target_display",
      "phenotype",
      "trait",
      "exposure_original",
      "exposure"
    )
  }
  
  clock_raw <- get_first_col(tbl, clock_candidates)
  disease_raw <- get_first_col(tbl, disease_candidates)
  disease_name_raw <- get_first_col(tbl, disease_name_candidates, default = disease_raw)
  
  source_raw <- get_first_col(tbl, c("target_source", "source"), default = "FinnGen")
  
  p_raw <- safe_num(get_first_col(tbl, c("pval", "P", "p", "p_value")))
  beta_raw <- safe_num(get_first_col(tbl, c("b", "beta")))
  se_raw <- safe_num(get_first_col(tbl, c("se", "SE", "beta_se")))
  
  or_raw <- safe_num(get_first_col(tbl, c("or", "OR", "odds_ratio")))
  
  lci_raw <- safe_num(get_first_col(
    tbl,
    c("or_lci95", "OR_lci95", "lo_ci", "lower_ci", "lci", "ci_low", "or_lower")
  ))
  
  uci_raw <- safe_num(get_first_col(
    tbl,
    c("or_uci95", "OR_uci95", "up_ci", "upper_ci", "uci", "ci_high", "or_upper")
  ))
  
  p_bon_raw <- safe_num(get_first_col(tbl, c("P_bon_n_disease", "p_bon_n_disease", "p_threshold")))
  
  sig_raw <- get_first_col(tbl, c("sig_by_disease", "significant_by_disease"), default = TRUE)
  sig_by_disease <- parse_bool(sig_raw)
  sig_by_disease[is.na(sig_by_disease)] <- TRUE
  
  out <- tbl %>%
    mutate(
      MR_direction_label = direction_label,
      direction_rank = direction_rank,
      
      clock_id = clean_clock_id(clock_raw),
      clock_id_join = str_to_lower(clock_id),
      
      disease_code = clean_disease_code(disease_raw),
      target_source = as.character(source_raw),
      
      disease_name_raw = disease_name_raw,
      endpoint_label = clean_text(disease_name_raw),
      endpoint_label = if_else(
        is.na(endpoint_label) | endpoint_label == "" | endpoint_label == "NA",
        clean_text(disease_code),
        endpoint_label
      ),
      
      p_value = p_raw,
      beta = beta_raw,
      se = se_raw,
      
      or_value = or_raw,
      lci_raw = lci_raw,
      uci_raw = uci_raw,
      
      # If OR is missing, compute from beta.
      or_value = if_else(is.na(or_value) & !is.na(beta), exp(beta), or_value),
      
      # CI handling:
      # If CI columns are already OR-scale, use them.
      # If they are beta-scale or missing, compute from beta +/- 1.96*SE.
      or_lci95 = case_when(
        !is.na(lci_raw) & !is.na(uci_raw) & lci_raw > 0 & uci_raw > 0 ~ lci_raw,
        !is.na(lci_raw) & !is.na(uci_raw) & (lci_raw <= 0 | uci_raw <= 0) ~ exp(lci_raw),
        !is.na(beta) & !is.na(se) ~ exp(beta - 1.96 * se),
        TRUE ~ NA_real_
      ),
      
      or_uci95 = case_when(
        !is.na(lci_raw) & !is.na(uci_raw) & lci_raw > 0 & uci_raw > 0 ~ uci_raw,
        !is.na(lci_raw) & !is.na(uci_raw) & (lci_raw <= 0 | uci_raw <= 0) ~ exp(uci_raw),
        !is.na(beta) & !is.na(se) ~ exp(beta + 1.96 * se),
        TRUE ~ NA_real_
      ),
      
      P_bon_n_disease = p_bon_raw,
      sig_by_disease = sig_by_disease,
      
      disease_chapter = assign_disease_chapter(target_source, disease_code)
    )
  
  out
}

# ----------------------------
# 6. Standardize and combine MR directions
# ----------------------------
clock2de <- standardize_mr(clock2de_raw, "clock2de")
de2clock <- standardize_mr(de2clock_raw, "de2clock")

mr_sig <- bind_rows(clock2de, de2clock) %>%
  filter(
    !is.na(p_value),
    !is.na(or_value),
    !is.na(or_lci95),
    !is.na(or_uci95)
  ) %>%
  filter(
    sig_by_disease |
      (!is.na(P_bon_n_disease) & p_value < P_bon_n_disease)
  ) %>%
  left_join(
    clock_key %>%
      select(
        clock_id_join,
        clock_block,
        clock_label,
        clock_suborder,
        clock_modality,
        clock_rank
      ),
    by = "clock_id_join"
  ) %>%
  mutate(
    clock_label = if_else(is.na(clock_label), clean_text(clock_id), clock_label),
    clock_modality = if_else(
      is.na(clock_modality),
      case_when(
        str_detect(str_to_lower(clock_id), "mri") ~ "MRI",
        str_detect(str_to_lower(clock_id), "proteomics") ~ "Proteomics",
        str_detect(str_to_lower(clock_id), "metabolomics") ~ "Metabolomics",
        TRUE ~ "Other"
      ),
      as.character(clock_modality)
    ),
    clock_modality = factor(clock_modality, levels = c("MRI", "Proteomics", "Metabolomics", "Other")),
    disease_chapter = factor(disease_chapter, levels = chapter_levels_all),
    MR_direction_label = factor(
      MR_direction_label,
      levels = c("EPOCH -> disease", "Disease -> EPOCH")
    ),
    neg_log10_P = -log10(p_value),
    endpoint_short = shorten_endpoint(endpoint_label),
    p_label = format_pval(p_value),
    or_label = sprintf("OR=%.2f", or_value),
    assoc_sign = if_else(or_value >= 1, "Positive", "Negative")
  ) %>%
  arrange(
    MR_direction_label,
    disease_chapter,
    endpoint_short,
    clock_label,
    p_value
  )

if (nrow(mr_sig) == 0) {
  stop("No significant-by-disease MR results found after filtering.")
}

message("Significant MR rows plotted: ", nrow(mr_sig))
message("  EPOCH -> disease: ", sum(mr_sig$MR_direction_label == "EPOCH -> disease"))
message("  Disease -> EPOCH: ", sum(mr_sig$MR_direction_label == "Disease -> EPOCH"))
message("Unique disease endpoints: ", n_distinct(mr_sig$disease_code))
message("Unique mortality EPOCH clocks: ", n_distinct(mr_sig$clock_id))

readr::write_tsv(
  mr_sig,
  paste0(out_prefix, "_significant_by_DE_MR_results.tsv")
)

summary_tbl <- mr_sig %>%
  count(MR_direction_label, disease_chapter, clock_modality, name = "n_signals") %>%
  arrange(MR_direction_label, disease_chapter, desc(n_signals))

readr::write_tsv(
  summary_tbl,
  paste0(out_prefix, "_summary_by_direction_disease_system_modality.tsv")
)

# ----------------------------
# 7. Plotting table
# ----------------------------
plot_df <- mr_sig %>%
  mutate(
    # y-axis shows the EPOCH clock; disease is annotated after the CI bar.
    row_clock_label = clock_label,
    
    # unique row ID prevents duplicate clock labels from collapsing.
    row_id = paste0(
      sprintf("%03d", row_number()),
      "__",
      row_clock_label
    ),
    
    disease_annotation = endpoint_short,
    
    disease_chapter = fct_drop(disease_chapter)
  ) %>%
  arrange(MR_direction_label, disease_chapter, endpoint_short, clock_label, p_value) %>%
  mutate(
    row_id = factor(row_id, levels = rev(unique(row_id))),
    
    # place disease labels just after the CI bar
    label_x = pmax(or_uci95, or_value, na.rm = TRUE) * 1.06
  )

# ----------------------------
# 8. Plot limits
# ----------------------------
x_min <- min(plot_df$or_lci95, plot_df$or_value, na.rm = TRUE)
x_max <- max(plot_df$label_x, plot_df$or_uci95, plot_df$or_value, na.rm = TRUE)

x_lim <- c(
  max(0.10, x_min * 0.80),
  x_max * 1.08
)

# ----------------------------
# 9. Colors and shapes
# ----------------------------
chapter_palette <- c(
  "PGC brain/psychiatric" = "#3B6EA5",
  "Neurologic" = "#4C78A8",
  "Psychiatric / sleep" = "#2F5D8C",
  "Eye" = "#74A9CF",
  "Ear / sensory" = "#9ECAE1",
  "Cardiovascular" = "#2E6B45",
  "Respiratory" = "#5AA05A",
  "Digestive / oral" = "#8E63A9",
  "Endocrine / immune" = "#9BAA4F",
  "Blood" = "#B2182B",
  "Cancer / benign neoplasm" = "#8B1A1A",
  "Skin / allergy" = "#C77CFF",
  "MSK / injury" = "#E39CB1",
  "GU / reproductive" = "#00008B",
  "Pregnancy" = "#E7298A",
  "Infectious" = "#A65E2E",
  "Other" = "#999999"
)

shape_values <- c(
  "MRI" = 21,
  "Proteomics" = 24,
  "Metabolomics" = 22,
  "Other" = 23
)

# ----------------------------
# 10. Build bidirectional MR forest figure
# ----------------------------
p <- ggplot(
  plot_df,
  aes(
    x = or_value,
    y = row_id,
    color = disease_chapter,
    fill = disease_chapter,
    shape = clock_modality
  )
) +
  
  # Null line, OR = 1.
  geom_vline(
    xintercept = 1,
    color = "#B2182B",
    linewidth = 0.50,
    linetype = "22"
  ) +
  
  # Error bars from MR output.
  geom_segment(
    aes(
      x = or_lci95,
      xend = or_uci95,
      yend = row_id
    ),
    linewidth = 0.70,
    alpha = 0.80,
    lineend = "round"
  ) +
  
  # Point estimate.
  geom_point(
    size = 3.4,
    stroke = 0.38,
    alpha = 0.98
  ) +
  
  # Disease endpoint name after the dot + error bar.
  geom_text(
    aes(
      x = label_x,
      label = disease_annotation
    ),
    hjust = 0,
    size = 2.55,
    lineheight = 0.88,
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = chapter_palette, drop = TRUE) +
  scale_fill_manual(values = chapter_palette, drop = TRUE) +
  scale_shape_manual(values = shape_values, drop = TRUE) +
  
  scale_x_log10(
    limits = x_lim,
    breaks = c(0.25, 0.5, 0.67, 0.8, 1, 1.25, 1.5, 2, 3, 4),
    labels = function(x) sprintf("%.2g", x)
  ) +
  
  scale_y_discrete(
    drop = TRUE,
    labels = function(x) str_replace(x, "^\\d+__", "")
  ) +
  
  facet_grid(
    rows = vars(MR_direction_label, disease_chapter),
    scales = "free_y",
    space = "free_y",
    switch = "y",
    drop = TRUE
  ) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    tag = "E",
    title = "Bidirectional Mendelian randomization links mortality EPOCH clocks and disease endpoints",
    subtitle = paste0(
      "Significant IVW MR associations corrected by the number of disease endpoints. ",
      "Panels distinguish EPOCH-to-disease and disease-to-EPOCH directions. ",
      "Points show MR OR; horizontal bars show 95% CI; dashed line marks OR = 1."
    ),
    x = "MR effect estimate, OR scale",
    y = "Mortality EPOCH clock",
    color = "Disease system",
    fill = "Disease system",
    shape = "Mortality EPOCH modality",
    caption = paste0(
      "Disease endpoint names are annotated to the right of each point and confidence interval. ",
      "EPOCH -> disease tests genetically predicted mortality EPOCH effects on disease endpoints; ",
      "Disease -> EPOCH tests genetically predicted disease liability effects on mortality EPOCH clocks. ",
      "Rows plotted: ", nrow(plot_df), "."
    )
  ) +
  
  guides(
    color = guide_legend(
      override.aes = list(size = 4),
      title.position = "top",
      title.hjust = 0
    ),
    fill = "none",
    shape = guide_legend(
      title.position = "top",
      title.hjust = 0
    )
  ) +
  
  theme_minimal(base_size = 9) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    
    panel.grid.major.y = element_line(color = "grey91", linewidth = 0.22),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey88", linewidth = 0.30),
    
    strip.placement = "outside",
    strip.background = element_rect(
      fill = "grey95",
      color = "grey78",
      linewidth = 0.25
    ),
    strip.text.y.left = element_text(
      angle = 0,
      face = "bold",
      size = 7.4,
      color = "grey20",
      margin = margin(l = 2, r = 4)
    ),
    
    axis.text.y = element_text(
      size = 5.9,
      color = "grey20",
      lineheight = 0.90
    ),
    axis.text.x = element_text(
      size = 8,
      color = "grey20"
    ),
    axis.title.x = element_text(
      size = 9,
      face = "bold",
      margin = margin(t = 6)
    ),
    axis.title.y = element_text(
      size = 9,
      face = "bold",
      margin = margin(r = 6)
    ),
    
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(face = "bold", size = 8.5),
    legend.text = element_text(size = 7.4),
    
    plot.tag = element_text(face = "bold", size = 22),
    plot.tag.position = c(0.006, 0.995),
    
    plot.title = element_text(face = "bold", size = 14.5),
    plot.subtitle = element_text(size = 8.7, color = "grey25", lineheight = 1.05),
    plot.caption = element_text(size = 7.2, color = "grey35", hjust = 0, lineheight = 1.05),
    
    # extra right margin for disease-name annotations
    plot.margin = margin(8, 185, 8, 8)
  )

print(p)

# ----------------------------
# 11. Save outputs
# ----------------------------
fig_height <- max(12.5, 0.18 * nrow(plot_df) + 3.2)
fig_width <- 17.2

ggsave(
  filename = paste0(out_prefix, ".pdf"),
  plot = p,
  width = fig_width,
  height = fig_height,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  filename = paste0(out_prefix, ".png"),
  plot = p,
  width = fig_width,
  height = fig_height,
  units = "in",
  dpi = 500,
  bg = "white"
)

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    filename = paste0(out_prefix, ".svg"),
    plot = p,
    width = fig_width,
    height = fig_height,
    units = "in",
    device = svglite::svglite,
    bg = "white"
  )
}

# ----------------------------
# 12. Also save direction-specific tables
# ----------------------------
readr::write_tsv(
  plot_df %>% filter(MR_direction_label == "EPOCH -> disease"),
  paste0(out_prefix, "_EPOCH_to_disease_rows.tsv")
)

readr::write_tsv(
  plot_df %>% filter(MR_direction_label == "Disease -> EPOCH"),
  paste0(out_prefix, "_disease_to_EPOCH_rows.tsv")
)

message("Done.")
message("Saved:")
message("  ", paste0(out_prefix, ".pdf"))
message("  ", paste0(out_prefix, ".png"))
if (requireNamespace("svglite", quietly = TRUE)) {
  message("  ", paste0(out_prefix, ".svg"))
}
message("  ", paste0(out_prefix, "_significant_by_DE_MR_results.tsv"))
message("  ", paste0(out_prefix, "_summary_by_direction_disease_system_modality.tsv"))
message("  ", paste0(out_prefix, "_EPOCH_to_disease_rows.tsv"))
message("  ", paste0(out_prefix, "_disease_to_EPOCH_rows.tsv"))