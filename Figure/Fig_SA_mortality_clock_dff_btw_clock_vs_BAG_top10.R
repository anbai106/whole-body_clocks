# ============================================================
# Panel D: Clock vs BAG effect-size difference for disease onset
#
# Fair selection:
#   Option 1:
#     clock_p < 0.05 / diseases / clocks OR
#     bag_p   < 0.05 / diseases / clocks
#
#   Option 2:
#     all valid clock-BAG disease pairs
#
# Effect size:
#   joint_log_HR_diff_clock_minus_BAG =
#     beta_clock_joint - beta_BAG_joint
#
# Main plotted effect:
#   joint_HR_ratio_clock_vs_BAG =
#     exp(beta_clock_joint - beta_BAG_joint)
#
# Matching statistical test:
#   joint_p_diff
#
# Plot:
#   Top 10 disease endpoints per mortality clock ranked by:
#     abs(joint_log_HR_diff_clock_minus_BAG)
#
# Labeling:
#   y-axis: ICD10 code + mortality clock only
#   inside plot: shortened disease name placed to the right of HR-ratio point
#   far right: number of cases / non-cases for each line
#
# Updates in this version:
#   1. Added manual ICD10 fallback mapping for codes not found by icd.data.
#   2. Removed white rectangle background under disease-name annotations by
#      replacing geom_label() with geom_text().
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(scales)
  library(grid)
})

# ----------------------------
# 1. Input / output
# ----------------------------
infile <- "/Users/hao/Downloads/clock_vs_BAG_all_rows_all_status.xlsx"

out_dir <- "/Users/hao/Downloads"

out_prefix <- file.path(
  out_dir,
  "panelD_clock_vs_BAG_FAIR_clock_or_BAG_top10_per_clock_short_disease_labels_cases"
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 2. Plot and analysis controls
# ----------------------------
diff_p_cutoff <- 0.05
top_n_per_clock <- 10

# Fair comparison selection mode:
#   "clock_or_bag" = selected if clock_p OR bag_p passes Bonferroni threshold
#   "all_valid"    = all valid clock-BAG disease pairs, regardless of marginal P
fair_selection_mode <- "clock_or_bag"

# Maximum length of disease-name labels inside the figure.
max_disease_label_chars <- 30

# ICD display:
#   FALSE: J440, E119, I10
#   TRUE:  J44.0, E11.9, I10
format_icd_with_dot <- FALSE

# Optional external ICD10 mapping file.
# Expected columns:
#   icd_code, disease_name
# Codes may be with or without dots, e.g. I10, E11.9, E119.
icd10_name_file <- NA_character_

# Optional external short-name file.
# Expected columns:
#   icd_code, disease_short_name
# This is optional and only used if you want manual display-name overrides.
icd10_short_name_file <- NA_character_

# Use icd.data if available, or install if allowed.
use_icd_data_package <- TRUE
auto_install_icd_data <- TRUE

# ----------------------------
# 3. Read data
# ----------------------------
df_raw <- readxl::read_excel(infile, sheet = 1)

required_cols <- c(
  "disease_id",
  "organ",
  "modality",
  "mortality_clock",
  "bag",
  "status",
  "N",
  "N_case",
  "N_noncase",
  "clock_p",
  "bag_p",
  "clock_hr",
  "bag_hr",
  "joint_p_diff"
)

missing_cols <- setdiff(required_cols, names(df_raw))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns from input file: ",
    paste(missing_cols, collapse = ", ")
  )
}

optional_cols <- c(
  "joint_beta_diff_clock_minus_bag",
  "joint_se_diff",
  "joint_z_diff",
  "clock_joint_beta",
  "bag_joint_beta",
  "clock_joint_hr",
  "bag_joint_hr",
  "clock_cindex",
  "bag_cindex",
  "base_cindex",
  "delta_cindex_clock_minus_bag"
)

for (cc in optional_cols) {
  if (!cc %in% names(df_raw)) {
    df_raw[[cc]] <- NA_real_
  }
}

candidate_name_cols <- c(
  "disease_full_name",
  "disease_name",
  "endpoint_full_name",
  "endpoint_name",
  "trait_name",
  "trait",
  "phenotype",
  "description",
  "outcome_name"
)

name_col <- intersect(candidate_name_cols, names(df_raw))

if (length(name_col) > 0) {
  name_col <- name_col[1]
  message("Using disease full-name column from file: ", name_col)
  df_raw$disease_full_name_from_file <- as.character(df_raw[[name_col]])
} else {
  message("No explicit disease-name column found. ICD10 mapping will be used.")
  df_raw$disease_full_name_from_file <- NA_character_
}

num_cols <- intersect(
  c(
    "N", "N_case", "N_noncase",
    "clock_beta", "clock_se", "clock_hr", "clock_ci_lo", "clock_ci_hi", "clock_p",
    "bag_beta", "bag_se", "bag_hr", "bag_ci_lo", "bag_ci_hi", "bag_p",
    "clock_joint_beta", "clock_joint_hr", "clock_joint_p",
    "bag_joint_beta", "bag_joint_hr", "bag_joint_p",
    "joint_beta_diff_clock_minus_bag",
    "joint_se_diff",
    "joint_z_diff",
    "joint_p_diff",
    "base_cindex", "clock_cindex", "bag_cindex", "both_cindex",
    "delta_cindex_clock_minus_bag",
    "delta_cindex_clock_minus_base",
    "delta_cindex_bag_minus_base"
  ),
  names(df_raw)
)

df <- df_raw %>%
  mutate(across(all_of(num_cols), ~ suppressWarnings(as.numeric(.x)))) %>%
  filter(status == "ok") %>%
  filter(!is.na(clock_p), !is.na(bag_p), !is.na(joint_p_diff))

# ----------------------------
# 4. Organ ordering and clock labels
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

standardize_clock_block <- function(organ) {
  case_when(
    str_detect(organ, regex("brain", ignore_case = TRUE)) ~ "Brain",
    str_detect(organ, regex("eye", ignore_case = TRUE)) ~ "Eye",
    str_detect(organ, regex("pulmonary|lung", ignore_case = TRUE)) ~ "Pulmonary",
    str_detect(organ, regex("heart", ignore_case = TRUE)) ~ "Heart",
    str_detect(organ, regex("liver|hepatic", ignore_case = TRUE)) ~ "Liver / hepatic",
    str_detect(organ, regex("kidney|renal", ignore_case = TRUE)) ~ "Kidney / renal",
    str_detect(organ, regex("pancreas", ignore_case = TRUE)) ~ "Pancreas",
    str_detect(organ, regex("spleen", ignore_case = TRUE)) ~ "Spleen",
    str_detect(organ, regex("immune", ignore_case = TRUE)) ~ "Immune",
    str_detect(organ, regex("endocrine", ignore_case = TRUE)) ~ "Endocrine",
    str_detect(organ, regex("digestive", ignore_case = TRUE)) ~ "Digestive",
    str_detect(organ, regex("metabolic", ignore_case = TRUE)) ~ "Metabolic",
    str_detect(organ, regex("adipose", ignore_case = TRUE)) ~ "Adipose",
    str_detect(organ, regex("skin", ignore_case = TRUE)) ~ "Skin",
    str_detect(organ, regex("reproductive", ignore_case = TRUE)) ~ "Reproductive",
    TRUE ~ organ
  )
}

pretty_modality <- function(modality) {
  case_when(
    str_detect(modality, regex("^MRI$|mri", ignore_case = TRUE)) ~ "MRI",
    str_detect(modality, regex("proteomics|protein", ignore_case = TRUE)) ~ "protein",
    str_detect(modality, regex("metabolomics|metabolite", ignore_case = TRUE)) ~ "metabolite",
    TRUE ~ as.character(modality)
  )
}

modality_group_fun <- function(modality) {
  case_when(
    str_detect(modality, regex("^MRI$|mri", ignore_case = TRUE)) ~ "MRI",
    str_detect(modality, regex("proteomics|protein", ignore_case = TRUE)) ~ "Proteomics",
    str_detect(modality, regex("metabolomics|metabolite", ignore_case = TRUE)) ~ "Metabolomics",
    TRUE ~ as.character(modality)
  )
}

modality_rank_fun <- function(modality) {
  case_when(
    str_detect(modality, regex("^MRI$|mri", ignore_case = TRUE)) ~ 1,
    str_detect(modality, regex("proteomics|protein", ignore_case = TRUE)) ~ 2,
    str_detect(modality, regex("metabolomics|metabolite", ignore_case = TRUE)) ~ 3,
    TRUE ~ 9
  )
}

df <- df %>%
  mutate(
    clock_block = standardize_clock_block(organ),
    clock_block = factor(clock_block, levels = organ_levels),
    modality_pretty = pretty_modality(modality),
    modality_group = modality_group_fun(modality),
    modality_group = factor(modality_group, levels = c("MRI", "Proteomics", "Metabolomics")),
    modality_rank = modality_rank_fun(modality),
    clock_label = paste(organ, modality_pretty),
    clock_label = str_replace_all(clock_label, "_", " "),
    clock_label = str_squish(clock_label)
  )

clock_block_palette <- c(
  "Brain" = "#4169A1",
  "Eye" = "#74A9CF",
  "Pulmonary" = "#5AA05A",
  "Heart" = "#2E6B45",
  "Liver / hepatic" = "#B08D00",
  "Kidney / renal" = "#0000A8",
  "Pancreas" = "#8E5BAE",
  "Spleen" = "#6A3D9A",
  "Immune" = "#1B9E77",
  "Endocrine" = "#9BAA4F",
  "Digestive" = "#8E5BAE",
  "Metabolic" = "#A6761D",
  "Adipose" = "#E6AB02",
  "Skin" = "#C77CFF",
  "Reproductive" = "#E7298A"
)

# ----------------------------
# 5. ICD10 helper functions
# ----------------------------
clean_icd_code <- function(x) {
  x %>%
    as.character() %>%
    str_to_upper() %>%
    str_replace_all("\\.", "") %>%
    str_replace_all("\\s+", "") %>%
    str_squish()
}

format_icd_code <- function(x) {
  x <- clean_icd_code(x)

  if (!format_icd_with_dot) {
    return(x)
  }

  ifelse(
    str_detect(x, "^[A-Z][0-9]{2}[A-Z0-9]+$"),
    paste0(str_sub(x, 1, 3), ".", str_sub(x, 4)),
    x
  )
}

extract_icd_code <- function(x) {
  x_clean <- x %>%
    as.character() %>%
    str_to_upper() %>%
    str_replace_all("\\.", "")

  str_extract(x_clean, "^[A-Z][0-9]{2}[A-Z0-9]*")
}

clean_name_candidate <- function(x) {
  x <- str_squish(as.character(x))

  x <- if_else(
    is.na(x) |
      x == "" |
      str_detect(str_to_lower(x), "disease name unavailable") |
      str_detect(str_to_lower(x), "^unavailable$") |
      str_detect(str_to_lower(x), "^unknown$") |
      str_detect(str_to_lower(x), "^na$"),
    NA_character_,
    x
  )

  x
}

clean_disease_name_from_id <- function(disease_id, icd_clean) {
  out <- disease_id %>%
    as.character() %>%
    str_remove("_clock_disease_free$") %>%
    str_remove("_disease_free$") %>%
    str_remove("_diagnosis$") %>%
    str_remove("_incident$") %>%
    str_remove("_prevalent$") %>%
    str_replace_all("_", " ") %>%
    str_squish()

  icd_pretty <- format_icd_code(icd_clean)

  out <- ifelse(
    !is.na(icd_clean),
    str_remove(out, paste0("^", icd_clean, "\\s*")),
    out
  )

  out <- ifelse(
    !is.na(icd_pretty),
    str_remove(out, paste0("^", fixed(icd_pretty), "\\s*")),
    out
  )

  out <- str_squish(out)
  out <- ifelse(out == "" | out == icd_clean | out == icd_pretty, NA_character_, out)
  out <- clean_name_candidate(out)

  out
}

# ----------------------------
# 6. Manual ICD10 fallback mapping
# ----------------------------
# This is used only as a fallback if the package / file mapping does not
# provide a useful name, and as a source of short figure labels.

manual_icd10_fallback_map <- tribble(
  ~icd_clean, ~disease_full_manual, ~disease_short_manual,

  "A099", "Gastroenteritis and colitis of infectious origin, unspecified", "Infectious gastroenteritis",
  "A419", "Sepsis, unspecified organism", "Sepsis",
  "B182", "Chronic viral hepatitis C", "Chronic hepatitis C",
  "B348", "Other viral infections of unspecified site", "Viral infection",
  "B968", "Other specified bacterial agents as the cause of diseases classified elsewhere", "Bacterial infection",

  "C220", "Liver cell carcinoma", "Liver cancer",
  "C64",  "Malignant neoplasm of kidney, except renal pelvis", "Kidney cancer",
  "C819", "Hodgkin lymphoma, unspecified", "Hodgkin lymphoma",
  "C880", "Waldenstrom macroglobulinemia", "Waldenstrom macroglobulinemia",

  "D45",  "Polycythemia vera", "Polycythemia vera",
  "D471", "Chronic myeloproliferative disease", "Myeloproliferative disease",
  "D474", "Osteomyelofibrosis", "Osteomyelofibrosis",
  "D509", "Iron deficiency anemia, unspecified", "Iron deficiency anemia",
  "D539", "Nutritional anemia, unspecified", "Nutritional anemia",
  "D638", "Anemia in other chronic diseases classified elsewhere", "Anemia of chronic disease",
  "D649", "Anemia, unspecified", "Anemia",
  "D696", "Thrombocytopenia, unspecified", "Thrombocytopenia",
  "D70",  "Agranulocytosis", "Agranulocytosis",

  "E102", "Type 1 diabetes mellitus with kidney complications", "T1D with kidney complications",
  "E103", "Type 1 diabetes mellitus with ophthalmic complications", "T1D with eye complications",
  "E104", "Type 1 diabetes mellitus with neurological complications", "T1D with neurologic complications",
  "E109", "Type 1 diabetes mellitus without complications", "T1D",
  "E111", "Type 2 diabetes mellitus with ketoacidosis", "T2D with ketoacidosis",
  "E116", "Type 2 diabetes mellitus with other specified complications", "T2D with complications",
  "E119", "Type 2 diabetes mellitus without complications", "T2D",
  "E143", "Unspecified diabetes mellitus with ophthalmic complications", "Diabetes with eye complications",
  "E211", "Secondary hyperparathyroidism", "Secondary hyperparathyroidism",
  "E271", "Primary adrenocortical insufficiency", "Adrenal insufficiency",
  "E834", "Disorders of magnesium metabolism", "Magnesium metabolism disorder",
  "E871", "Hypo-osmolality and hyponatremia", "Hyponatremia",
  "E892", "Postprocedural hypoparathyroidism", "Postprocedural endocrine disorder",

  "F101", "Mental and behavioral disorders due to use of alcohol, harmful use", "Alcohol harmful use",

  "G312", "Degeneration of nervous system due to alcohol", "Alcohol-related neurodegeneration",
  "G35",  "Multiple sclerosis", "Multiple sclerosis",
  "G590", "Diabetic mononeuropathy", "Diabetic mononeuropathy",
  "G990", "Autonomic neuropathy in diseases classified elsewhere", "Autonomic neuropathy",

  "I052", "Rheumatic mitral stenosis with insufficiency", "Rheumatic mitral stenosis",
  "I10",  "Essential hypertension", "Hypertension",
  "I20",  "Angina pectoris", "Angina",
  "I21",  "Acute myocardial infarction", "Acute MI",
  "I25",  "Chronic ischemic heart disease", "Ischemic heart disease",
  "I251", "Atherosclerotic heart disease of native coronary artery", "Coronary atherosclerosis",
  "I279", "Pulmonary heart disease, unspecified", "Pulmonary heart disease",
  "I421", "Obstructive hypertrophic cardiomyopathy", "Hypertrophic cardiomyopathy",
  "I429", "Cardiomyopathy, unspecified", "Cardiomyopathy",
  "I442", "Atrioventricular block, complete", "Complete AV block",
  "I447", "Left bundle-branch block, unspecified", "Left bundle branch block",
  "I460", "Cardiac arrest with successful resuscitation", "Cardiac arrest",
  "I48",  "Atrial fibrillation and flutter", "Atrial fibrillation/flutter",
  "I500", "Congestive heart failure", "Congestive heart failure",
  "I501", "Left ventricular failure", "Left ventricular failure",
  "I509", "Heart failure, unspecified", "Heart failure",
  "I518", "Other ill-defined heart diseases", "Ill-defined heart disease",
  "I7021", "Atherosclerosis of native arteries of extremities", "Peripheral atherosclerosis",
  "I7080", "Atherosclerosis of other arteries", "Atherosclerosis",
  "I745", "Embolism and thrombosis of iliac artery", "Iliac artery thrombosis",
  "I780", "Hereditary hemorrhagic telangiectasia", "Hereditary telangiectasia",
  "I850", "Esophageal varices with bleeding", "Bleeding esophageal varices",
  "I859", "Esophageal varices without bleeding", "Esophageal varices",
  "I864", "Gastric varices", "Gastric varices",
  "I982", "Esophageal varices in diseases classified elsewhere", "Secondary esophageal varices",

  "J150", "Pneumonia due to Klebsiella pneumoniae", "Klebsiella pneumonia",
  "J181", "Lobar pneumonia, unspecified organism", "Lobar pneumonia",
  "J189", "Pneumonia, unspecified organism", "Pneumonia",
  "J432", "Centrilobular emphysema", "Centrilobular emphysema",
  "J440", "Chronic obstructive pulmonary disease with acute lower respiratory infection", "COPD with acute LRI",
  "J449", "Chronic obstructive pulmonary disease, unspecified", "COPD",
  "J848", "Other specified interstitial pulmonary diseases", "Interstitial lung disease",

  "K210", "Gastro-esophageal reflux disease with esophagitis", "GERD with esophagitis",
  "K500", "Crohn disease of small intestine", "Crohn disease",
  "K573", "Diverticular disease of large intestine without perforation or abscess", "Diverticular disease",
  "K658", "Other peritonitis", "Peritonitis",
  "K701", "Alcoholic hepatitis", "Alcoholic hepatitis",
  "K704", "Alcoholic hepatic failure", "Alcoholic liver failure",
  "K709", "Alcoholic liver disease, unspecified", "Alcoholic liver disease",
  "K743", "Primary biliary cirrhosis", "Primary biliary cirrhosis",
  "K745", "Biliary cirrhosis, unspecified", "Biliary cirrhosis",
  "K754", "Autoimmune hepatitis", "Autoimmune hepatitis",
  "K758", "Other specified inflammatory liver diseases", "Inflammatory liver disease",
  "K767", "Hepatorenal syndrome", "Hepatorenal syndrome",
  "K802", "Calculus of gallbladder without cholecystitis", "Gallstones",
  "K912", "Postsurgical malabsorption", "Postsurgical malabsorption",

  "M0590", "Felty syndrome, unspecified site", "Seropositive RA",
  "M0599", "Seropositive rheumatoid arthritis, unspecified site", "Seropositive RA",
  "M0691", "Rheumatoid arthritis, unspecified, shoulder", "Rheumatoid arthritis",
  "M1309", "Polyarthritis, unspecified", "Polyarthritis",
  "M313", "Wegener granulomatosis", "Granulomatosis with polyangiitis",
  "M341", "CR(E)ST syndrome", "CREST syndrome",
  "M349", "Systemic sclerosis, unspecified", "Systemic sclerosis",
  "M351", "Sicca syndrome", "Sicca syndrome",
  "M7986", "Other specified soft tissue disorders", "Soft tissue disorder",

  "N028", "Recurrent and persistent hematuria", "Recurrent hematuria",
  "N049", "Nephrotic syndrome, unspecified", "Nephrotic syndrome",
  "N083", "Glomerular disorders in diabetes mellitus", "Diabetic kidney disease",
  "N180", "End-stage renal disease", "End-stage renal disease",
  "N184", "Chronic kidney disease, stage 4", "CKD stage 4",
  "N185", "Chronic kidney disease, stage 5", "CKD stage 5",
  "N189", "Chronic kidney disease, unspecified", "CKD",
  "N258", "Other disorders resulting from impaired renal tubular function", "Renal tubular dysfunction",
  "N811", "Cystocele", "Cystocele",

  "Q612", "Polycystic kidney, adult type", "Adult polycystic kidney",

  "R000", "Tachycardia, unspecified", "Tachycardia",
  "R268", "Other abnormalities of gait and mobility", "Gait abnormality",
  "R296", "Tendency to fall, not elsewhere classified", "Falls tendency",
  "R35",  "Polyuria", "Polyuria",
  "R410", "Disorientation, unspecified", "Disorientation",
  "R509", "Fever, unspecified", "Fever",
  "R601", "Generalized edema", "Generalized edema",
  "R69",  "Unknown and unspecified causes of morbidity", "Unspecified morbidity",

  "T179", "Foreign body in respiratory tract, part unspecified", "Airway foreign body",
  "T824", "Mechanical complication of vascular dialysis catheter", "Dialysis catheter complication",
  "T861", "Kidney transplant failure and rejection", "Kidney transplant failure",

  "W010", "Fall on same level from slipping, tripping and stumbling", "Slip/trip fall",
  "W190", "Unspecified fall", "Unspecified fall",

  "Y495", "Adverse effects of psychotropic drugs", "Psychotropic drug adverse effect",
  "Y545", "Adverse effects of mineralocorticoids and their antagonists", "Mineralocorticoid adverse effect",
  "Y575", "Adverse effects of contrast media", "Contrast media adverse effect",
  "Y830", "Surgical operation with transplant of whole organ", "Transplant surgery complication",
  "Y841", "Kidney dialysis as the cause of abnormal reaction", "Dialysis complication",
  "Y95",  "Nosocomial condition", "Nosocomial condition",

  "Z490", "Preparatory care for dialysis", "Dialysis preparatory care",
  "Z491", "Extracorporeal dialysis", "Extracorporeal dialysis",
  "Z515", "Palliative care", "Palliative care",
  "Z852", "Personal history of malignant neoplasm of other respiratory and intrathoracic organs", "History of thoracic cancer",
  "Z903", "Acquired absence of part of stomach", "Partial gastrectomy status",
  "Z922", "Personal history of long-term drug therapy", "Long-term medication history",
  "Z940", "Kidney transplant status", "Kidney transplant status",
  "Z944", "Liver transplant status", "Liver transplant status",
  "Z950", "Presence of cardiac pacemaker", "Pacemaker status",
  "Z992", "Dependence on renal dialysis", "Dialysis dependence"
) %>%
  mutate(
    icd_clean = clean_icd_code(icd_clean),
    disease_full_manual = str_squish(disease_full_manual),
    disease_short_manual = str_squish(disease_short_manual)
  ) %>%
  distinct(icd_clean, .keep_all = TRUE)

# ----------------------------
# 7. Load ICD10 mapping
# ----------------------------
guess_code_col <- function(nms) {
  nms_low <- str_to_lower(nms)
  hits <- nms[str_detect(nms_low, "^(code|icd|icd_code|icd10|icd10_code)$")]
  if (length(hits) == 0) hits <- nms[str_detect(nms_low, "code|icd")]
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

guess_desc_col <- function(nms) {
  nms_low <- str_to_lower(nms)
  hits <- nms[str_detect(nms_low, "long.*desc|long.*description|description|desc|name|label|title")]
  if (length(hits) == 0) hits <- nms[str_detect(nms_low, "short.*desc")]
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

standardize_icd_map <- function(tbl, source_name = "unknown") {
  tbl <- as_tibble(tbl)

  code_col <- guess_code_col(names(tbl))
  desc_col <- guess_desc_col(names(tbl))

  if (is.na(code_col) || is.na(desc_col)) {
    warning("Could not identify code/name columns in ICD mapping source: ", source_name)
    return(tibble(icd_clean = character(), disease_name_icd = character(), source = character()))
  }

  tbl %>%
    transmute(
      icd_clean = clean_icd_code(.data[[code_col]]),
      disease_name_icd = clean_name_candidate(.data[[desc_col]])
    ) %>%
    filter(!is.na(icd_clean), icd_clean != "", !is.na(disease_name_icd), disease_name_icd != "") %>%
    mutate(
      disease_name_icd = str_squish(disease_name_icd),
      source = source_name
    ) %>%
    distinct(icd_clean, .keep_all = TRUE)
}

read_external_icd10_map <- function(path) {
  if (is.na(path) || !file.exists(path)) {
    return(tibble(icd_clean = character(), disease_name_icd = character(), source = character()))
  }

  ext <- readr::read_tsv(path, show_col_types = FALSE)

  if (!all(c("icd_code", "disease_name") %in% names(ext))) {
    stop("External ICD10 mapping file must contain columns: icd_code and disease_name")
  }

  ext %>%
    transmute(
      icd_clean = clean_icd_code(icd_code),
      disease_name_icd = clean_name_candidate(disease_name)
    ) %>%
    filter(!is.na(icd_clean), icd_clean != "", !is.na(disease_name_icd), disease_name_icd != "") %>%
    mutate(source = "external_file") %>%
    distinct(icd_clean, .keep_all = TRUE)
}

read_external_short_name_map <- function(path) {
  if (is.na(path) || !file.exists(path)) {
    return(tibble(icd_clean = character(), disease_short_external = character()))
  }

  ext <- readr::read_tsv(path, show_col_types = FALSE)

  if (!all(c("icd_code", "disease_short_name") %in% names(ext))) {
    stop("External short-name file must contain columns: icd_code and disease_short_name")
  }

  ext %>%
    transmute(
      icd_clean = clean_icd_code(icd_code),
      disease_short_external = clean_name_candidate(disease_short_name)
    ) %>%
    filter(!is.na(icd_clean), icd_clean != "", !is.na(disease_short_external), disease_short_external != "") %>%
    distinct(icd_clean, .keep_all = TRUE)
}

load_icd_data_map <- function(auto_install = FALSE) {

  if (!requireNamespace("icd.data", quietly = TRUE)) {
    if (auto_install) {
      try(
        install.packages("icd.data", repos = "https://cloud.r-project.org"),
        silent = TRUE
      )
    }
  }

  if (!requireNamespace("icd.data", quietly = TRUE)) {
    warning("Package icd.data is not installed. ICD10 names will rely on external file, manual map, or disease_id parsing.")
    return(tibble(icd_clean = character(), disease_name_icd = character(), source = character()))
  }

  pkg_data <- utils::data(package = "icd.data")$results

  if (is.null(pkg_data) || nrow(pkg_data) == 0) {
    warning("No datasets found in icd.data.")
    return(tibble(icd_clean = character(), disease_name_icd = character(), source = character()))
  }

  data_names <- pkg_data[, "Item"] %>% as.character()

  icd10_names <- data_names[str_detect(data_names, "^icd10cm[0-9]{4}$")]

  if (length(icd10_names) == 0) {
    icd10_names <- data_names[str_detect(str_to_lower(data_names), "icd10")]
  }

  if (length(icd10_names) == 0) {
    warning("No ICD10 dataset found in icd.data.")
    return(tibble(icd_clean = character(), disease_name_icd = character(), source = character()))
  }

  chosen <- sort(icd10_names, decreasing = TRUE)[1]
  message("Using ICD10 mapping from icd.data object: ", chosen)

  e <- new.env()
  utils::data(list = chosen, package = "icd.data", envir = e)
  obj <- get(chosen, envir = e)

  if (is.vector(obj) && !is.null(names(obj)) && !is.data.frame(obj)) {
    out <- tibble(
      icd_clean = clean_icd_code(names(obj)),
      disease_name_icd = clean_name_candidate(obj),
      source = paste0("icd.data::", chosen)
    ) %>%
      filter(!is.na(icd_clean), icd_clean != "", !is.na(disease_name_icd), disease_name_icd != "") %>%
      distinct(icd_clean, .keep_all = TRUE)

    return(out)
  }

  standardize_icd_map(obj, source_name = paste0("icd.data::", chosen))
}

manual_icd10_full_map <- manual_icd10_fallback_map %>%
  transmute(
    icd_clean,
    disease_name_icd = disease_full_manual,
    source = "manual_fallback"
  )

manual_icd10_short_map <- manual_icd10_fallback_map %>%
  transmute(
    icd_clean,
    disease_short_manual
  )

icd10_map <- bind_rows(
  read_external_icd10_map(icd10_name_file),
  if (use_icd_data_package) {
    load_icd_data_map(auto_install = auto_install_icd_data)
  } else {
    tibble(icd_clean = character(), disease_name_icd = character(), source = character())
  },
  manual_icd10_full_map
) %>%
  filter(!is.na(icd_clean), icd_clean != "") %>%
  distinct(icd_clean, .keep_all = TRUE)

icd10_map_icd3 <- icd10_map %>%
  mutate(icd3 = str_sub(icd_clean, 1, 3)) %>%
  group_by(icd3) %>%
  summarise(
    disease_name_icd3 = disease_name_icd[which.min(nchar(icd_clean))][1],
    .groups = "drop"
  )

icd10_short_map <- bind_rows(
  read_external_short_name_map(icd10_short_name_file),
  manual_icd10_short_map %>%
    transmute(
      icd_clean,
      disease_short_external = disease_short_manual
    )
) %>%
  filter(!is.na(icd_clean), icd_clean != "", !is.na(disease_short_external), disease_short_external != "") %>%
  distinct(icd_clean, .keep_all = TRUE)

readr::write_tsv(
  icd10_map,
  paste0(out_prefix, "_ICD10_mapping_used.tsv")
)

readr::write_tsv(
  manual_icd10_fallback_map,
  paste0(out_prefix, "_manual_ICD10_fallback_map.tsv")
)

# ----------------------------
# 8. Disease-name shortening function
# ----------------------------
shorten_disease_name <- function(x) {
  x <- clean_name_candidate(x)

  x <- str_replace_all(x, regex("chronic obstructive pulmonary disease", ignore_case = TRUE), "COPD")
  x <- str_replace_all(x, regex("acute lower respiratory infection", ignore_case = TRUE), "acute LRI")
  x <- str_replace_all(x, regex("lower respiratory infection", ignore_case = TRUE), "LRI")
  x <- str_replace_all(x, regex("gastro-esophageal reflux disease", ignore_case = TRUE), "GERD")
  x <- str_replace_all(x, regex("gastroesophageal reflux disease", ignore_case = TRUE), "GERD")
  x <- str_replace_all(x, regex("type 1 diabetes mellitus", ignore_case = TRUE), "T1D")
  x <- str_replace_all(x, regex("type 2 diabetes mellitus", ignore_case = TRUE), "T2D")
  x <- str_replace_all(x, regex("diabetes mellitus", ignore_case = TRUE), "diabetes")
  x <- str_replace_all(x, regex("chronic kidney disease", ignore_case = TRUE), "CKD")
  x <- str_replace_all(x, regex("end-stage renal disease", ignore_case = TRUE), "ESRD")
  x <- str_replace_all(x, regex("rheumatoid arthritis", ignore_case = TRUE), "RA")
  x <- str_replace_all(x, regex("malignant neoplasm of", ignore_case = TRUE), "cancer of")
  x <- str_replace_all(x, regex("atherosclerotic heart disease of native coronary artery", ignore_case = TRUE), "coronary atherosclerosis")
  x <- str_replace_all(x, regex("atrial fibrillation and flutter", ignore_case = TRUE), "atrial fibrillation/flutter")
  x <- str_replace_all(x, regex("pneumonia, unspecified organism", ignore_case = TRUE), "pneumonia")
  x <- str_replace_all(x, regex("heart failure, unspecified", ignore_case = TRUE), "heart failure")
  x <- str_replace_all(x, regex("essential \\(primary\\) hypertension", ignore_case = TRUE), "hypertension")
  x <- str_replace_all(x, regex("essential hypertension", ignore_case = TRUE), "hypertension")
  x <- str_replace_all(x, regex("acute myocardial infarction", ignore_case = TRUE), "acute MI")
  x <- str_replace_all(x, regex("myocardial infarction", ignore_case = TRUE), "MI")
  x <- str_replace_all(x, regex("atherosclerosis of native arteries of extremities", ignore_case = TRUE), "peripheral atherosclerosis")
  x <- str_replace_all(x, regex("conductive and sensorineural hearing loss", ignore_case = TRUE), "hearing loss")
  x <- str_replace_all(x, regex("personal history of malignant neoplasm", ignore_case = TRUE), "history of cancer")
  x <- str_replace_all(x, regex("adverse effects? of", ignore_case = TRUE), "adverse effect:")

  x <- str_replace_all(x, regex(", unspecified", ignore_case = TRUE), "")
  x <- str_replace_all(x, regex("unspecified ", ignore_case = TRUE), "")
  x <- str_replace_all(x, regex("not elsewhere classified", ignore_case = TRUE), "NEC")
  x <- str_replace_all(x, regex("without mention of", ignore_case = TRUE), "without")
  x <- str_replace_all(x, regex("with other specified complications", ignore_case = TRUE), "with complications")
  x <- str_replace_all(x, regex("without complications", ignore_case = TRUE), "")
  x <- str_replace_all(x, regex("other specified", ignore_case = TRUE), "other")
  x <- str_replace_all(x, regex("other ", ignore_case = TRUE), "")

  x <- str_squish(x)

  x <- if_else(
    is.na(x) | x == "",
    NA_character_,
    paste0(str_to_upper(str_sub(x, 1, 1)), str_sub(x, 2))
  )

  x
}

# ----------------------------
# 9. Disease labels
# ----------------------------
df <- df %>%
  mutate(
    icd_clean = extract_icd_code(disease_id),
    icd_clean = if_else(is.na(icd_clean), clean_icd_code(disease_id), icd_clean),
    icd_pretty = format_icd_code(icd_clean),
    icd3 = str_sub(icd_clean, 1, 3),
    disease_name_from_id = clean_disease_name_from_id(disease_id, icd_clean)
  ) %>%
  left_join(
    icd10_map %>%
      select(icd_clean, disease_name_icd, source),
    by = "icd_clean"
  ) %>%
  left_join(
    icd10_map_icd3,
    by = "icd3"
  ) %>%
  left_join(
    icd10_short_map,
    by = "icd_clean"
  ) %>%
  mutate(
    disease_full_name = coalesce(
      clean_name_candidate(disease_full_name_from_file),
      clean_name_candidate(disease_name_from_id),
      clean_name_candidate(disease_name_icd),
      clean_name_candidate(disease_name_icd3)
    ),

    disease_short_name = coalesce(
      clean_name_candidate(disease_short_external),
      shorten_disease_name(disease_full_name),
      icd_pretty
    ),

    disease_short_name = if_else(
      is.na(disease_short_name) | disease_short_name == "",
      icd_pretty,
      disease_short_name
    ),

    disease_short_name = str_trunc(
      disease_short_name,
      width = max_disease_label_chars,
      side = "right"
    ),

    disease_label_icd_only = icd_pretty
  )

# ----------------------------
# 10. Construct standardized Clock - BAG effect size
# ----------------------------
df <- df %>%
  mutate(
    joint_log_HR_diff_clock_minus_BAG = case_when(
      !is.na(joint_beta_diff_clock_minus_bag) ~ joint_beta_diff_clock_minus_bag,

      is.na(joint_beta_diff_clock_minus_bag) &
        !is.na(clock_joint_beta) &
        !is.na(bag_joint_beta) ~ clock_joint_beta - bag_joint_beta,

      is.na(joint_beta_diff_clock_minus_bag) &
        !is.na(clock_joint_hr) &
        !is.na(bag_joint_hr) &
        clock_joint_hr > 0 &
        bag_joint_hr > 0 ~ log(clock_joint_hr) - log(bag_joint_hr),

      TRUE ~ NA_real_
    ),

    abs_joint_log_HR_diff_clock_minus_BAG = abs(joint_log_HR_diff_clock_minus_BAG),
    joint_HR_ratio_clock_vs_BAG = exp(joint_log_HR_diff_clock_minus_BAG),

    diff_direction = case_when(
      joint_log_HR_diff_clock_minus_BAG > 0 ~ "Mortality clock stronger",
      joint_log_HR_diff_clock_minus_BAG < 0 ~ "BAG stronger",
      TRUE ~ "No difference"
    )
  )

if (all(is.na(df$joint_log_HR_diff_clock_minus_BAG))) {
  stop(
    "Could not construct joint_log_HR_diff_clock_minus_BAG. ",
    "Need one of: joint_beta_diff_clock_minus_bag, ",
    "clock_joint_beta + bag_joint_beta, or clock_joint_hr + bag_joint_hr."
  )
}

# ----------------------------
# 11. Fair clock-BAG comparison threshold
# ----------------------------
n_unique_diseases <- n_distinct(df$disease_id)
n_mortality_clocks <- n_distinct(df$mortality_clock)

p_pair_bonf <- 0.05 / n_unique_diseases / n_mortality_clocks

threshold_tbl <- tibble(
  selection_rule = case_when(
    fair_selection_mode == "clock_or_bag" ~
      "clock_p < 0.05 / diseases / clocks OR bag_p < 0.05 / diseases / clocks",
    fair_selection_mode == "all_valid" ~
      "all valid clock-BAG disease pairs",
    TRUE ~ fair_selection_mode
  ),
  n_unique_diseases = n_unique_diseases,
  n_mortality_clocks = n_mortality_clocks,
  n_tests = n_unique_diseases * n_mortality_clocks,
  pairwise_bonferroni_threshold = p_pair_bonf,
  effect_size = "joint_log_HR_diff_clock_minus_BAG",
  effect_size_display = "joint_HR_ratio_clock_vs_BAG = exp(joint_log_HR_diff_clock_minus_BAG)",
  top_n_per_mortality_clock_for_plot = top_n_per_clock,
  ranking_for_plot = "absolute value of joint_log_HR_diff_clock_minus_BAG",
  difference_test_p = "joint_p_diff",
  difference_p_cutoff = diff_p_cutoff,
  fair_selection_mode = fair_selection_mode
)

print(threshold_tbl)

readr::write_tsv(
  threshold_tbl,
  paste0(out_prefix, "_fair_selection_thresholds.tsv")
)

# ----------------------------
# 12. Select fair comparison disease results
# ----------------------------
df_fair_base <- df %>%
  mutate(
    clock_sig_bonf = clock_p < p_pair_bonf,
    bag_sig_bonf = bag_p < p_pair_bonf,

    selected_by = case_when(
      clock_sig_bonf & bag_sig_bonf ~ "Both clock and BAG significant",
      clock_sig_bonf & !bag_sig_bonf ~ "Clock only significant",
      !clock_sig_bonf & bag_sig_bonf ~ "BAG only significant",
      TRUE ~ "Neither significant"
    ),

    selected_by = factor(
      selected_by,
      levels = c(
        "Both clock and BAG significant",
        "Clock only significant",
        "BAG only significant",
        "Neither significant"
      )
    )
  )

if (fair_selection_mode == "clock_or_bag") {

  sig_all <- df_fair_base %>%
    filter(clock_sig_bonf | bag_sig_bonf)

} else if (fair_selection_mode == "all_valid") {

  sig_all <- df_fair_base %>%
    filter(
      !is.na(clock_p),
      !is.na(bag_p),
      !is.na(joint_p_diff),
      !is.na(joint_log_HR_diff_clock_minus_BAG)
    )

} else {

  stop("Unknown fair_selection_mode. Use 'clock_or_bag' or 'all_valid'.")
}

sig_all <- sig_all %>%
  filter(!is.na(joint_log_HR_diff_clock_minus_BAG), !is.na(joint_p_diff)) %>%
  mutate(
    diff_sig_nominal = joint_p_diff < diff_p_cutoff,
    diff_p_fdr = p.adjust(joint_p_diff, method = "BH"),
    diff_p_bonf = p.adjust(joint_p_diff, method = "bonferroni"),
    diff_sig_fdr = diff_p_fdr < 0.05,
    diff_sig_bonf = diff_p_bonf < 0.05,

    clock_p_label = case_when(
      clock_p < 1e-300 ~ "clock p<1e-300",
      clock_p < 1e-3 ~ paste0("clock p=", formatC(clock_p, format = "e", digits = 1)),
      TRUE ~ paste0("clock p=", formatC(clock_p, format = "f", digits = 3))
    ),

    bag_p_label = case_when(
      bag_p < 1e-300 ~ "BAG p<1e-300",
      bag_p < 1e-3 ~ paste0("BAG p=", formatC(bag_p, format = "e", digits = 1)),
      TRUE ~ paste0("BAG p=", formatC(bag_p, format = "f", digits = 3))
    ),

    diff_p_label = case_when(
      joint_p_diff < 1e-300 ~ "joint p<1e-300",
      joint_p_diff < 1e-3 ~ paste0("joint p=", formatC(joint_p_diff, format = "e", digits = 1)),
      TRUE ~ paste0("joint p=", formatC(joint_p_diff, format = "f", digits = 3))
    ),

    case_noncase_label = paste0(
      "Cases/non-cases: ",
      comma(N_case),
      "/",
      comma(N_noncase)
    ),

    case_noncase_label_short = paste0(
      comma(N_case),
      "/",
      comma(N_noncase)
    ),

    overpower_class = case_when(
      joint_p_diff >= diff_p_cutoff ~ "No significant Clock-BAG difference",

      joint_p_diff < diff_p_cutoff &
        joint_log_HR_diff_clock_minus_BAG > 0 ~ "Mortality clock overpowered",

      joint_p_diff < diff_p_cutoff &
        joint_log_HR_diff_clock_minus_BAG < 0 ~ "BAG overpowered",

      TRUE ~ "No significant Clock-BAG difference"
    ),

    overpower_class = factor(
      overpower_class,
      levels = c(
        "No significant Clock-BAG difference",
        "Mortality clock overpowered",
        "BAG overpowered"
      )
    )
  ) %>%
  arrange(
    desc(abs_joint_log_HR_diff_clock_minus_BAG),
    joint_p_diff,
    pmin(clock_p, bag_p, na.rm = TRUE)
  ) %>%
  mutate(global_abs_diff_rank = row_number())

message("Fair selection mode: ", fair_selection_mode)
message("Number of selected fair-comparison disease associations: ", nrow(sig_all))
message("Number selected by both clock and BAG: ", sum(sig_all$selected_by == "Both clock and BAG significant", na.rm = TRUE))
message("Number selected by clock only: ", sum(sig_all$selected_by == "Clock only significant", na.rm = TRUE))
message("Number selected by BAG only: ", sum(sig_all$selected_by == "BAG only significant", na.rm = TRUE))
message("Number with nominally significant Clock-vs-BAG difference: ", sum(sig_all$diff_sig_nominal, na.rm = TRUE))
message("Number with FDR-significant Clock-vs-BAG difference: ", sum(sig_all$diff_sig_fdr, na.rm = TRUE))
message("Number with Bonferroni-significant Clock-vs-BAG difference: ", sum(sig_all$diff_sig_bonf, na.rm = TRUE))

if (nrow(sig_all) == 0) {
  stop("No disease associations were selected under the current fair_selection_mode.")
}

# Save all selected rows, not just plotted rows
sig_all_out <- sig_all %>%
  select(
    global_abs_diff_rank,
    selected_by,
    clock_sig_bonf,
    bag_sig_bonf,
    overpower_class,
    disease_id,
    icd_clean,
    icd_pretty,
    disease_full_name,
    disease_short_name,
    source,
    mortality_clock,
    bag,
    organ,
    modality,
    modality_group,
    clock_block,
    clock_label,
    N,
    N_case,
    N_noncase,
    case_noncase_label,
    clock_hr,
    bag_hr,
    clock_p,
    bag_p,
    clock_p_label,
    bag_p_label,
    joint_log_HR_diff_clock_minus_BAG,
    abs_joint_log_HR_diff_clock_minus_BAG,
    joint_HR_ratio_clock_vs_BAG,
    joint_beta_diff_clock_minus_bag,
    joint_se_diff,
    joint_z_diff,
    joint_p_diff,
    diff_p_label,
    diff_p_fdr,
    diff_p_bonf,
    diff_sig_nominal,
    diff_sig_fdr,
    diff_sig_bonf,
    diff_direction,
    clock_cindex,
    bag_cindex,
    delta_cindex_clock_minus_bag,
    everything()
  )

readr::write_tsv(
  sig_all_out,
  paste0(out_prefix, "_fair_selected_with_short_disease_labels.tsv")
)

unmapped_icd <- sig_all_out %>%
  filter(is.na(disease_full_name) | disease_full_name == "" | disease_short_name == icd_pretty) %>%
  distinct(icd_clean, icd_pretty, disease_id) %>%
  arrange(icd_clean)

readr::write_tsv(
  unmapped_icd,
  paste0(out_prefix, "_unmapped_ICD10_codes.tsv")
)

selection_summary <- sig_all %>%
  count(selected_by, overpower_class, name = "n_results") %>%
  group_by(selected_by) %>%
  mutate(
    total_selected_by_group = sum(n_results),
    percent_within_selected_by_group = n_results / total_selected_by_group
  ) %>%
  ungroup()

readr::write_tsv(
  selection_summary,
  paste0(out_prefix, "_fair_selection_summary.tsv")
)

modality_summary <- sig_all %>%
  count(modality_group, selected_by, overpower_class, name = "n_results") %>%
  group_by(modality_group) %>%
  mutate(
    total_modality_results = sum(n_results),
    percent_modality_results = n_results / total_modality_results
  ) %>%
  ungroup()

readr::write_tsv(
  modality_summary,
  paste0(out_prefix, "_fair_modality_summary.tsv")
)

# ----------------------------
# 13. Select top 10 absolute differences per mortality clock
# ----------------------------
clock_order_tbl <- sig_all %>%
  distinct(mortality_clock, clock_block, clock_label, modality_rank) %>%
  arrange(clock_block, modality_rank, clock_label) %>%
  mutate(clock_plot_rank = row_number())

sig_plot <- sig_all %>%
  group_by(mortality_clock) %>%
  arrange(
    desc(abs_joint_log_HR_diff_clock_minus_BAG),
    joint_p_diff,
    pmin(clock_p, bag_p, na.rm = TRUE),
    .by_group = TRUE
  ) %>%
  slice_head(n = top_n_per_clock) %>%
  mutate(rank_within_clock_abs_diff = row_number()) %>%
  ungroup() %>%
  left_join(
    clock_order_tbl %>%
      select(mortality_clock, clock_plot_rank),
    by = "mortality_clock"
  ) %>%
  arrange(
    clock_plot_rank,
    rank_within_clock_abs_diff,
    desc(abs_joint_log_HR_diff_clock_minus_BAG),
    joint_p_diff,
    pmin(clock_p, bag_p, na.rm = TRUE)
  ) %>%
  mutate(
    plot_row_number = row_number(),
    y = rev(plot_row_number),

    row_label_plot = paste0(
      icd_pretty,
      " | ",
      as.character(clock_label)
    ),

    row_label_plot = str_squish(row_label_plot)
  )

message("Number of rows shown in plot: ", nrow(sig_plot))
message("Maximum plotted rows per mortality clock: ", top_n_per_clock)

readr::write_tsv(
  sig_plot,
  paste0(out_prefix, "_fair_top10_per_clock_plot_rows_short_disease_labels.tsv")
)

# ----------------------------
# 14. Main plot: adjusted HR ratio, Clock / BAG
# ----------------------------
sig_plot_ratio <- sig_plot %>%
  filter(!is.na(joint_HR_ratio_clock_vs_BAG), joint_HR_ratio_clock_vs_BAG > 0)

if (nrow(sig_plot_ratio) == 0) {
  stop("No valid HR-ratio rows to plot.")
}

ratio_vals <- sig_plot_ratio$joint_HR_ratio_clock_vs_BAG
ratio_min <- min(ratio_vals, na.rm = TRUE)
ratio_max <- max(ratio_vals, na.rm = TRUE)

ratio_lim_data <- c(
  min(ratio_min * 0.80, 1 / ratio_max * 0.80),
  max(ratio_max * 1.20, 1 / ratio_min * 1.20)
)

sig_plot_ratio <- sig_plot_ratio %>%
  mutate(
    disease_label_x = case_when(
      joint_HR_ratio_clock_vs_BAG >= 1 ~ joint_HR_ratio_clock_vs_BAG * 1.08,
      joint_HR_ratio_clock_vs_BAG < 1 ~ 1.08,
      TRUE ~ 1.08
    ),
    disease_label_hjust = 0
  )

ratio_lim <- ratio_lim_data
ratio_lim[2] <- max(
  ratio_lim_data[2],
  max(sig_plot_ratio$disease_label_x, na.rm = TRUE) * 1.20
)

ratio_breaks <- c(0.1, 0.2, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 3, 4, 5, 8, 10, 15)
ratio_breaks <- ratio_breaks[ratio_breaks >= ratio_lim[1] & ratio_breaks <= ratio_lim[2]]

case_label_x_ratio <- ratio_lim[2] * 1.06
case_header_y <- max(sig_plot_ratio$y, na.rm = TRUE) + 1.1

plot_height <- max(8.0, 2.8 + 0.135 * nrow(sig_plot_ratio))

subtitle_txt <- paste0(
  "Disease associations were selected using a fair comparison rule: ",
  if_else(
    fair_selection_mode == "clock_or_bag",
    paste0(
      "clock P or BAG P < 0.05 / ",
      n_unique_diseases,
      " diseases / ",
      n_mortality_clocks,
      " clocks = ",
      signif(p_pair_bonf, 3)
    ),
    "all valid clock-BAG disease pairs"
  ),
  ". The figure shows the top ",
  top_n_per_clock,
  " absolute Clock - BAG differences per mortality clock."
)

caption_txt <- paste0(
  "Adjusted HR ratio is exp(beta_clock - beta_BAG) from the joint Cox model. ",
  "Values >1 indicate stronger conditional mortality-clock effects; values <1 indicate stronger BAG effects. ",
  "Triangles indicate joint_p_diff < ",
  diff_p_cutoff,
  ". Y-axis shows ICD10 code + mortality clock; shortened disease names are placed to the right of each HR-ratio point. ",
  "The far-right column shows cases/non-cases."
)

p_ratio <- ggplot(sig_plot_ratio, aes(y = y)) +

  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "#D62728",
    linewidth = 0.45
  ) +

  geom_segment(
    aes(
      x = 1,
      xend = joint_HR_ratio_clock_vs_BAG,
      yend = y,
      color = clock_block
    ),
    linewidth = 0.65,
    alpha = 0.70
  ) +

  geom_point(
    aes(
      x = joint_HR_ratio_clock_vs_BAG,
      fill = clock_block,
      shape = diff_sig_nominal
    ),
    color = "grey10",
    stroke = 0.30,
    size = 2.35,
    alpha = 0.98
  ) +

  # Disease-name annotation without rectangle background.
  geom_text(
    aes(
      x = disease_label_x,
      y = y,
      label = disease_short_name,
      color = clock_block
    ),
    hjust = 0,
    size = 1.75,
    lineheight = 0.90,
    show.legend = FALSE
  ) +

  geom_text(
    data = sig_plot_ratio,
    aes(
      x = case_label_x_ratio,
      y = y,
      label = case_noncase_label_short
    ),
    hjust = 0,
    size = 1.75,
    color = "grey15",
    inherit.aes = FALSE
  ) +

  annotate(
    "text",
    x = case_label_x_ratio,
    y = case_header_y,
    label = "Cases/non-cases",
    hjust = 0,
    vjust = 0,
    size = 2.25,
    fontface = "bold",
    color = "grey10"
  ) +

  scale_y_continuous(
    breaks = sig_plot_ratio$y,
    labels = sig_plot_ratio$row_label_plot,
    expand = expansion(add = c(0.6, 1.5))
  ) +

  scale_x_log10(
    breaks = ratio_breaks,
    labels = ratio_breaks,
    expand = expansion(mult = c(0.03, 0.06))
  ) +

  coord_cartesian(
    xlim = ratio_lim,
    clip = "off"
  ) +

  scale_color_manual(
    values = clock_block_palette,
    drop = FALSE,
    name = "Mortality-clock organ"
  ) +

  scale_fill_manual(
    values = clock_block_palette,
    drop = FALSE,
    name = "Mortality-clock organ"
  ) +

  scale_shape_manual(
    values = c(
      "FALSE" = 21,
      "TRUE" = 24
    ),
    labels = c(
      "FALSE" = paste0("joint_p_diff >= ", diff_p_cutoff),
      "TRUE" = paste0("joint_p_diff < ", diff_p_cutoff)
    ),
    name = "Clock vs BAG difference"
  ) +

  labs(
    tag = "D",
    title = "Mortality clocks versus matched BAGs across fair-selected disease endpoints",
    subtitle = subtitle_txt,
    x = "Adjusted HR ratio, mortality clock / BAG",
    y = NULL,
    caption = caption_txt
  ) +

  guides(
    color = guide_legend(
      override.aes = list(size = 3),
      ncol = 1
    ),
    fill = "none",
    shape = guide_legend(
      override.aes = list(size = 3)
    )
  ) +

  theme_minimal(base_size = 9) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),

    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.22),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey88", linewidth = 0.22),

    axis.text.y = element_text(
      size = 5.4,
      color = "grey15",
      lineheight = 0.90
    ),
    axis.text.x = element_text(
      size = 8.0,
      color = "grey15"
    ),
    axis.title.x = element_text(
      size = 9.6,
      face = "bold"
    ),

    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(face = "bold", size = 8.5),
    legend.text = element_text(size = 7.2),

    plot.tag = element_text(face = "bold", size = 24),
    plot.tag.position = c(0.005, 0.995),

    plot.title = element_text(face = "bold", size = 13.5),
    plot.subtitle = element_text(
      size = 8.6,
      color = "grey25",
      lineheight = 1.05
    ),
    plot.caption = element_text(
      size = 7.2,
      color = "grey35",
      hjust = 0,
      lineheight = 1.05
    ),

    plot.margin = margin(8, 120, 8, 8)
  )

print(p_ratio)

ggsave(
  filename = paste0(out_prefix, "_HR_ratio_top10_per_clock_short_disease_labels_cases.pdf"),
  plot = p_ratio,
  width = 16.8,
  height = plot_height,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = paste0(out_prefix, "_HR_ratio_top10_per_clock_short_disease_labels_cases.png"),
  plot = p_ratio,
  width = 16.8,
  height = plot_height,
  units = "in",
  dpi = 500,
  bg = "white"
)

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    filename = paste0(out_prefix, "_HR_ratio_top10_per_clock_short_disease_labels_cases.svg"),
    plot = p_ratio,
    width = 16.8,
    height = plot_height,
    units = "in",
    device = svglite::svglite,
    bg = "white"
  )
}

# ----------------------------
# 15. Secondary plot: joint log-HR difference
# ----------------------------
sig_plot_loghr <- sig_plot %>%
  mutate(
    abs_x = abs(joint_log_HR_diff_clock_minus_BAG)
  )

x_vals <- sig_plot_loghr$joint_log_HR_diff_clock_minus_BAG
x_abs <- max(abs(x_vals), na.rm = TRUE)
x_abs <- max(x_abs, 0.05)

x_lim_data <- c(-1.15 * x_abs, 1.15 * x_abs)
x_range_data <- diff(x_lim_data)

sig_plot_loghr <- sig_plot_loghr %>%
  mutate(
    disease_label_x = case_when(
      joint_log_HR_diff_clock_minus_BAG >= 0 ~ joint_log_HR_diff_clock_minus_BAG + 0.035 * x_range_data,
      joint_log_HR_diff_clock_minus_BAG < 0 ~ 0 + 0.035 * x_range_data,
      TRUE ~ 0 + 0.035 * x_range_data
    )
  )

x_lim <- x_lim_data
x_lim[2] <- max(
  x_lim_data[2],
  max(sig_plot_loghr$disease_label_x, na.rm = TRUE) + 0.15 * x_range_data
)

x_range <- diff(x_lim)

case_label_x_loghr <- x_lim[2] + 0.06 * x_range
case_header_y_loghr <- max(sig_plot_loghr$y, na.rm = TRUE) + 1.1

p_loghr <- ggplot(sig_plot_loghr, aes(y = y)) +

  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "#D62728",
    linewidth = 0.45
  ) +

  geom_segment(
    aes(
      x = 0,
      xend = joint_log_HR_diff_clock_minus_BAG,
      yend = y,
      color = clock_block
    ),
    linewidth = 0.65,
    alpha = 0.70
  ) +

  geom_point(
    aes(
      x = joint_log_HR_diff_clock_minus_BAG,
      fill = clock_block,
      shape = diff_sig_nominal
    ),
    color = "grey10",
    stroke = 0.30,
    size = 2.35,
    alpha = 0.98
  ) +

  # Disease-name annotation without rectangle background.
  geom_text(
    aes(
      x = disease_label_x,
      y = y,
      label = disease_short_name,
      color = clock_block
    ),
    hjust = 0,
    size = 1.75,
    lineheight = 0.90,
    show.legend = FALSE
  ) +

  geom_text(
    data = sig_plot_loghr,
    aes(
      x = case_label_x_loghr,
      y = y,
      label = case_noncase_label_short
    ),
    hjust = 0,
    size = 1.75,
    color = "grey15",
    inherit.aes = FALSE
  ) +

  annotate(
    "text",
    x = case_label_x_loghr,
    y = case_header_y_loghr,
    label = "Cases/non-cases",
    hjust = 0,
    vjust = 0,
    size = 2.25,
    fontface = "bold",
    color = "grey10"
  ) +

  scale_y_continuous(
    breaks = sig_plot_loghr$y,
    labels = sig_plot_loghr$row_label_plot,
    expand = expansion(add = c(0.6, 1.5))
  ) +

  scale_x_continuous(
    labels = label_number(accuracy = 0.01),
    expand = expansion(mult = c(0.03, 0.06))
  ) +

  coord_cartesian(
    xlim = x_lim,
    clip = "off"
  ) +

  scale_color_manual(
    values = clock_block_palette,
    drop = FALSE,
    name = "Mortality-clock organ"
  ) +

  scale_fill_manual(
    values = clock_block_palette,
    drop = FALSE,
    name = "Mortality-clock organ"
  ) +

  scale_shape_manual(
    values = c(
      "FALSE" = 21,
      "TRUE" = 24
    ),
    labels = c(
      "FALSE" = paste0("joint_p_diff >= ", diff_p_cutoff),
      "TRUE" = paste0("joint_p_diff < ", diff_p_cutoff)
    ),
    name = "Clock vs BAG difference"
  ) +

  labs(
    tag = "D",
    title = "Mortality clocks versus matched BAGs across fair-selected disease endpoints",
    subtitle = subtitle_txt,
    x = "Joint log-HR difference, mortality clock - BAG",
    y = NULL,
    caption = paste0(
      "Values >0 indicate stronger mortality-clock effects; values <0 indicate stronger BAG effects. ",
      "Y-axis shows ICD10 code + mortality clock; shortened disease names are placed to the right of each effect-size point. ",
      "The far-right column shows cases/non-cases."
    )
  ) +

  guides(
    color = guide_legend(
      override.aes = list(size = 3),
      ncol = 1
    ),
    fill = "none",
    shape = guide_legend(
      override.aes = list(size = 3)
    )
  ) +

  theme_minimal(base_size = 9) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.22),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey88", linewidth = 0.22),

    axis.text.y = element_text(size = 5.4, color = "grey15", lineheight = 0.90),
    axis.text.x = element_text(size = 8.0, color = "grey15"),
    axis.title.x = element_text(size = 9.6, face = "bold"),

    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(face = "bold", size = 8.5),
    legend.text = element_text(size = 7.2),

    plot.tag = element_text(face = "bold", size = 24),
    plot.tag.position = c(0.005, 0.995),

    plot.title = element_text(face = "bold", size = 13.5),
    plot.subtitle = element_text(size = 8.6, color = "grey25", lineheight = 1.05),
    plot.caption = element_text(size = 7.2, color = "grey35", hjust = 0, lineheight = 1.05),

    plot.margin = margin(8, 120, 8, 8)
  )

print(p_loghr)

ggsave(
  filename = paste0(out_prefix, "_joint_logHR_diff_top10_per_clock_short_disease_labels_cases.pdf"),
  plot = p_loghr,
  width = 16.8,
  height = plot_height,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = paste0(out_prefix, "_joint_logHR_diff_top10_per_clock_short_disease_labels_cases.png"),
  plot = p_loghr,
  width = 16.8,
  height = plot_height,
  units = "in",
  dpi = 500,
  bg = "white"
)

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    filename = paste0(out_prefix, "_joint_logHR_diff_top10_per_clock_short_disease_labels_cases.svg"),
    plot = p_loghr,
    width = 16.8,
    height = plot_height,
    units = "in",
    device = svglite::svglite,
    bg = "white"
  )
}

# ----------------------------
# 16. Fair comparison summary tables
# ----------------------------
summary_by_selection <- sig_all %>%
  group_by(selected_by, overpower_class) %>%
  summarise(
    n_pairs = n(),
    n_unique_diseases = n_distinct(disease_id),
    n_unique_clocks = n_distinct(mortality_clock),
    median_joint_log_HR_diff = median(joint_log_HR_diff_clock_minus_BAG, na.rm = TRUE),
    median_abs_joint_log_HR_diff = median(abs_joint_log_HR_diff_clock_minus_BAG, na.rm = TRUE),
    median_HR_ratio = median(joint_HR_ratio_clock_vs_BAG, na.rm = TRUE),
    .groups = "drop"
  )

summary_by_modality <- sig_all %>%
  group_by(modality_group, overpower_class) %>%
  summarise(
    n_pairs = n(),
    n_unique_diseases = n_distinct(disease_id),
    n_unique_clocks = n_distinct(mortality_clock),
    median_joint_log_HR_diff = median(joint_log_HR_diff_clock_minus_BAG, na.rm = TRUE),
    median_abs_joint_log_HR_diff = median(abs_joint_log_HR_diff_clock_minus_BAG, na.rm = TRUE),
    median_HR_ratio = median(joint_HR_ratio_clock_vs_BAG, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(modality_group) %>%
  mutate(
    total_modality_pairs = sum(n_pairs),
    percent_modality_pairs = n_pairs / total_modality_pairs
  ) %>%
  ungroup()

summary_by_clock_all <- sig_all %>%
  group_by(clock_block, clock_label, mortality_clock, modality_group) %>%
  summarise(
    n_selected_pairs = n(),
    n_clock_sig_only = sum(selected_by == "Clock only significant", na.rm = TRUE),
    n_bag_sig_only = sum(selected_by == "BAG only significant", na.rm = TRUE),
    n_both_sig = sum(selected_by == "Both clock and BAG significant", na.rm = TRUE),
    n_clock_overpowered = sum(overpower_class == "Mortality clock overpowered", na.rm = TRUE),
    n_bag_overpowered = sum(overpower_class == "BAG overpowered", na.rm = TRUE),
    n_no_diff = sum(overpower_class == "No significant Clock-BAG difference", na.rm = TRUE),
    median_joint_log_HR_diff = median(joint_log_HR_diff_clock_minus_BAG, na.rm = TRUE),
    median_abs_joint_log_HR_diff = median(abs_joint_log_HR_diff_clock_minus_BAG, na.rm = TRUE),
    median_HR_ratio = median(joint_HR_ratio_clock_vs_BAG, na.rm = TRUE),
    min_joint_p_diff = min(joint_p_diff, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(clock_block, desc(median_abs_joint_log_HR_diff))

summary_by_clock_top10 <- sig_plot %>%
  group_by(clock_block, clock_label, mortality_clock, modality_group) %>%
  summarise(
    n_plotted_pairs = n(),
    total_cases_plotted = sum(N_case, na.rm = TRUE),
    total_noncases_plotted = sum(N_noncase, na.rm = TRUE),
    n_clock_overpowered = sum(overpower_class == "Mortality clock overpowered", na.rm = TRUE),
    n_bag_overpowered = sum(overpower_class == "BAG overpowered", na.rm = TRUE),
    n_no_diff = sum(overpower_class == "No significant Clock-BAG difference", na.rm = TRUE),
    max_abs_joint_log_HR_diff = max(abs_joint_log_HR_diff_clock_minus_BAG, na.rm = TRUE),
    median_abs_joint_log_HR_diff = median(abs_joint_log_HR_diff_clock_minus_BAG, na.rm = TRUE),
    min_joint_p_diff = min(joint_p_diff, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(clock_block, desc(max_abs_joint_log_HR_diff))

readr::write_tsv(
  summary_by_selection,
  paste0(out_prefix, "_summary_by_fair_selection.tsv")
)

readr::write_tsv(
  summary_by_modality,
  paste0(out_prefix, "_summary_by_modality.tsv")
)

readr::write_tsv(
  summary_by_clock_all,
  paste0(out_prefix, "_summary_by_clock_all_selected.tsv")
)

readr::write_tsv(
  summary_by_clock_top10,
  paste0(out_prefix, "_summary_by_clock_top10_plotted.tsv")
)

message("Saved files:")
message("  ", paste0(out_prefix, "_fair_selection_thresholds.tsv"))
message("  ", paste0(out_prefix, "_ICD10_mapping_used.tsv"))
message("  ", paste0(out_prefix, "_manual_ICD10_fallback_map.tsv"))
message("  ", paste0(out_prefix, "_unmapped_ICD10_codes.tsv"))
message("  ", paste0(out_prefix, "_fair_selected_with_short_disease_labels.tsv"))
message("  ", paste0(out_prefix, "_fair_selection_summary.tsv"))
message("  ", paste0(out_prefix, "_fair_modality_summary.tsv"))
message("  ", paste0(out_prefix, "_fair_top10_per_clock_plot_rows_short_disease_labels.tsv"))
message("  ", paste0(out_prefix, "_HR_ratio_top10_per_clock_short_disease_labels_cases.pdf/png/svg"))
message("  ", paste0(out_prefix, "_joint_logHR_diff_top10_per_clock_short_disease_labels_cases.pdf/png/svg"))
message("  ", paste0(out_prefix, "_summary_by_fair_selection.tsv"))
message("  ", paste0(out_prefix, "_summary_by_modality.tsv"))
message("  ", paste0(out_prefix, "_summary_by_clock_all_selected.tsv"))
message("  ", paste0(out_prefix, "_summary_by_clock_top10_plotted.tsv"))