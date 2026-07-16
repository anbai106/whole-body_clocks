#!/usr/bin/env Rscript

# ============================================================
# Mortality EPOCH-only disease-endpoint architecture
#
# Revised version:
#   1. Plot only 22 mortality EPOCH clocks.
#   2. Use same Bonferroni threshold as Clock-vs-BAG:
#        clock_p < 0.05 / disease endpoints / mortality EPOCH clocks
#   3. Reassign Z-code and organ-informative external/complication
#      endpoints to dominant organ domains when possible.
#   4. Replace Panel B heatmap with endpoint-organ-domain word-cloud boxes.
#   5. Show shortened disease names, not ICD10 codes, inside each box.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(scales)
  library(grid)
  library(patchwork)
})

# ============================================================
# 1. Input / output
# ============================================================

infile_candidates <- c(
  "/Users/hao/Downloads/EPOCH_result/clock_vs_BAG_all_rows_all_status.xlsx",
  "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/Result/clock_vs_BAG_survival_summary/clock_vs_BAG_all_rows_all_status.xlsx")

infile_candidates <- infile_candidates[file.exists(infile_candidates)]

if (length(infile_candidates) == 0) {
  stop("Could not find input file. Please edit infile_candidates or set infile manually.")
}

infile <- infile_candidates[1]

out_dir <- file.path(
  dirname(infile),
  "mortality_EPOCH_only_whole_body_risk_architecture_wordcloud"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_prefix <- file.path(
  out_dir,
  "mortality_EPOCH_only_whole_body_risk_architecture_wordcloud"
)

message("Input file:  ", infile)
message("Output dir:  ", out_dir)

# ============================================================
# 2. Analysis controls
# ============================================================

alpha_bonf <- 0.05
discard_metabolic_metabolomics <- FALSE

format_icd_with_dot <- FALSE
max_disease_label_chars <- 34
neglog10_cap <- 50

MAX_WORDS_PER_ENDPOINT_ORGAN_DOMAIN <- 24
MIN_ENDPOINT_PER_EDGE <- 1

# If TRUE, tries to install icd.data if missing.
# For CUBIC, I recommend FALSE.
use_icd_data_package <- TRUE
auto_install_icd_data <- FALSE

# Optional user-provided ICD mapping.
# Expected columns: icd_code, disease_name
external_icd10_name_file <- NA_character_

# ============================================================
# 3. Read and validate data
# ============================================================

df_raw <- readxl::read_excel(infile, sheet = 1)

required_cols <- c(
  "disease_id",
  "organ",
  "modality",
  "mortality_clock",
  "status",
  "N",
  "N_case",
  "N_noncase",
  "clock_p"
)

missing_cols <- setdiff(required_cols, names(df_raw))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns from input file: ",
    paste(missing_cols, collapse = ", ")
  )
}

if (!"clock_hr" %in% names(df_raw)) {
  if ("clock_beta" %in% names(df_raw)) {
    df_raw$clock_hr <- exp(as.numeric(df_raw$clock_beta))
  } else {
    df_raw$clock_hr <- NA_real_
  }
}

optional_numeric_cols <- c(
  "N",
  "N_case",
  "N_noncase",
  "clock_beta",
  "clock_se",
  "clock_hr",
  "clock_ci_lo",
  "clock_ci_hi",
  "clock_p",
  "clock_cindex",
  "base_cindex",
  "delta_cindex_clock_minus_base"
)

for (cc in optional_numeric_cols) {
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
  message("Using disease-name column from input: ", name_col)
  df_raw$disease_full_name_from_file <- as.character(df_raw[[name_col]])
} else {
  message("No explicit disease-name column found; ICD10/name parsing will be used.")
  df_raw$disease_full_name_from_file <- NA_character_
}

num_cols <- intersect(optional_numeric_cols, names(df_raw))

df <- df_raw %>%
  mutate(across(all_of(num_cols), ~ suppressWarnings(as.numeric(.x)))) %>%
  filter(status == "ok") %>%
  filter(!is.na(clock_p), is.finite(clock_p), clock_p > 0)

if (nrow(df) == 0) {
  stop("No status == 'ok' mortality EPOCH rows with finite clock_p.")
}

# ============================================================
# 4. Organ, modality, and clock labels
# ============================================================

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
  "Reproductive",
  "Other"
)

endpoint_organ_levels <- c(
  "Brain / neurologic",
  "Eye / ENT",
  "Pulmonary",
  "Heart / cardiovascular",
  "Liver / hepatic",
  "Kidney / renal",
  "Digestive",
  "Endocrine / metabolic",
  "Immune / blood",
  "Musculoskeletal",
  "Skin",
  "Reproductive",
  "Cancer / neoplasm",
  "Infectious",
  "Injury / external",
  "Clinical symptoms",
  "Other"
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
    TRUE ~ "Other"
  )
}

standardize_organ_display <- function(organ) {
  x <- as.character(organ)
  x <- str_replace_all(x, "_", " ")
  x <- str_squish(x)
  
  case_when(
    str_detect(x, regex("^brain$", ignore_case = TRUE)) ~ "Brain",
    str_detect(x, regex("^eye$", ignore_case = TRUE)) ~ "Eye",
    str_detect(x, regex("^pulmonary|^lung", ignore_case = TRUE)) ~ "Pulmonary",
    str_detect(x, regex("^heart$", ignore_case = TRUE)) ~ "Heart",
    str_detect(x, regex("^liver$|^hepatic$", ignore_case = TRUE)) ~ "Hepatic",
    str_detect(x, regex("^kidney$|^renal$", ignore_case = TRUE)) ~ "Renal",
    str_detect(x, regex("^pancreas$", ignore_case = TRUE)) ~ "Pancreas",
    str_detect(x, regex("^spleen$", ignore_case = TRUE)) ~ "Spleen",
    str_detect(x, regex("^immune$", ignore_case = TRUE)) ~ "Immune",
    str_detect(x, regex("^endocrine$", ignore_case = TRUE)) ~ "Endocrine",
    str_detect(x, regex("^digestive$", ignore_case = TRUE)) ~ "Digestive",
    str_detect(x, regex("^metabolic$", ignore_case = TRUE)) ~ "Metabolic",
    str_detect(x, regex("^adipose$", ignore_case = TRUE)) ~ "Adipose",
    str_detect(x, regex("^skin$", ignore_case = TRUE)) ~ "Skin",
    str_detect(x, regex("reproductive female", ignore_case = TRUE)) ~ "Reproductive female",
    str_detect(x, regex("reproductive male", ignore_case = TRUE)) ~ "Reproductive male",
    str_detect(x, regex("reproductive", ignore_case = TRUE)) ~ "Reproductive",
    TRUE ~ str_to_sentence(str_to_lower(x))
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

modality_shape_values <- c(
  "MRI" = 21,
  "Proteomics" = 24,
  "Metabolomics" = 22
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
  "Reproductive" = "#E7298A",
  "Other" = "#999999"
)

endpoint_domain_palette <- c(
  "Brain / neurologic" = "#4169A1",
  "Eye / ENT" = "#74A9CF",
  "Pulmonary" = "#5AA05A",
  "Heart / cardiovascular" = "#2E6B45",
  "Liver / hepatic" = "#B08D00",
  "Kidney / renal" = "#0000A8",
  "Digestive" = "#8E5BAE",
  "Endocrine / metabolic" = "#A6761D",
  "Immune / blood" = "#1B9E77",
  "Musculoskeletal" = "#6A3D9A",
  "Skin" = "#C77CFF",
  "Reproductive" = "#E7298A",
  "Cancer / neoplasm" = "#B2182B",
  "Infectious" = "#1B7837",
  "Injury / external" = "#666666",
  "Clinical symptoms" = "#8C8C8C",
  "Other" = "#999999"
)

df <- df %>%
  mutate(
    clock_block = standardize_clock_block(organ),
    clock_block = factor(clock_block, levels = organ_levels),
    organ_display = standardize_organ_display(organ),
    modality_pretty = pretty_modality(modality),
    modality_group = modality_group_fun(modality),
    modality_group = factor(modality_group, levels = c("MRI", "Proteomics", "Metabolomics")),
    modality_rank = modality_rank_fun(modality),
    clock_label = paste(organ_display, modality_pretty),
    clock_label = str_replace_all(clock_label, "_", " "),
    clock_label = str_squish(clock_label)
  )

if (discard_metabolic_metabolomics) {
  df <- df %>%
    filter(!(as.character(clock_block) == "Metabolic" & as.character(modality_group) == "Metabolomics"))
}

# ============================================================
# 5. ICD10 helpers
# ============================================================

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
      str_detect(str_to_lower(x), "^na$") |
      str_detect(str_to_lower(x), "^nan$"),
    NA_character_,
    x
  )
  
  x
}

is_icd_like <- function(x) {
  str_detect(str_to_upper(str_squish(as.character(x))), "^[A-Z][0-9]{2}[A-Z0-9]*$")
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

# ============================================================
# 6. ICD10 mapping
# ============================================================

manual_icd10_fallback_map <- tribble(
  ~icd_clean, ~disease_full_manual, ~disease_short_manual,
  
  "A099", "Gastroenteritis and colitis of infectious origin, unspecified", "Infectious gastroenteritis",
  "A419", "Sepsis, unspecified organism", "Sepsis",
  "B182", "Chronic viral hepatitis C", "Chronic hepatitis C",
  "B348", "Other viral infections of unspecified site", "Viral infection",
  "B968", "Other specified bacterial agents as the cause of diseases classified elsewhere", "Bacterial infection",
  
  "C220", "Liver cell carcinoma", "Liver cancer",
  "C64", "Malignant neoplasm of kidney, except renal pelvis", "Kidney cancer",
  "C819", "Hodgkin lymphoma, unspecified", "Hodgkin lymphoma",
  "C880", "Waldenstrom macroglobulinemia", "Waldenstrom macroglobulinemia",
  
  "D45", "Polycythemia vera", "Polycythemia vera",
  "D471", "Chronic myeloproliferative disease", "Myeloproliferative disease",
  "D474", "Osteomyelofibrosis", "Osteomyelofibrosis",
  "D509", "Iron deficiency anemia, unspecified", "Iron deficiency anemia",
  "D539", "Nutritional anemia, unspecified", "Nutritional anemia",
  "D638", "Anemia in other chronic diseases classified elsewhere", "Anemia of chronic disease",
  "D649", "Anemia, unspecified", "Anemia",
  "D696", "Thrombocytopenia, unspecified", "Thrombocytopenia",
  "D70", "Agranulocytosis", "Agranulocytosis",
  
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
  "G35", "Multiple sclerosis", "Multiple sclerosis",
  "G590", "Diabetic mononeuropathy", "Diabetic mononeuropathy",
  "G990", "Autonomic neuropathy in diseases classified elsewhere", "Autonomic neuropathy",
  
  "I052", "Rheumatic mitral stenosis with insufficiency", "Rheumatic mitral stenosis",
  "I10", "Essential hypertension", "Hypertension",
  "I20", "Angina pectoris", "Angina",
  "I21", "Acute myocardial infarction", "Acute MI",
  "I25", "Chronic ischemic heart disease", "Ischemic heart disease",
  "I251", "Atherosclerotic heart disease of native coronary artery", "Coronary atherosclerosis",
  "I279", "Pulmonary heart disease, unspecified", "Pulmonary heart disease",
  "I421", "Obstructive hypertrophic cardiomyopathy", "Hypertrophic cardiomyopathy",
  "I429", "Cardiomyopathy, unspecified", "Cardiomyopathy",
  "I442", "Atrioventricular block, complete", "Complete AV block",
  "I447", "Left bundle-branch block, unspecified", "Left bundle branch block",
  "I460", "Cardiac arrest with successful resuscitation", "Cardiac arrest",
  "I48", "Atrial fibrillation and flutter", "Atrial fibrillation/flutter",
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
  "R35", "Polyuria", "Polyuria",
  "R410", "Disorientation, unspecified", "Disorientation",
  "R509", "Fever, unspecified", "Fever",
  "R601", "Generalized edema", "Generalized edema",
  "R69", "Unknown and unspecified causes of morbidity", "Unspecified morbidity",
  
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
  "Y95", "Nosocomial condition", "Nosocomial condition",
  
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
    mutate(source = source_name) %>%
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
    warning("Package icd.data is not installed. ICD10 names will rely on disease_id parsing and manual fallback map.")
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
      disease_name_icd = clean_name_candidate(disease_name),
      source = "external_file"
    ) %>%
    filter(!is.na(icd_clean), icd_clean != "", !is.na(disease_name_icd), disease_name_icd != "") %>%
    distinct(icd_clean, .keep_all = TRUE)
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
  read_external_icd10_map(external_icd10_name_file),
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

readr::write_tsv(
  icd10_map,
  paste0(out_prefix, "_ICD10_mapping_used.tsv")
)

readr::write_tsv(
  manual_icd10_fallback_map,
  paste0(out_prefix, "_manual_ICD10_fallback_map.tsv")
)

# ============================================================
# 7. Disease-name shortening and domain reassignment
# ============================================================

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
  x <- str_replace_all(x, regex("acute myocardial infarction", ignore_case = TRUE), "acute MI")
  x <- str_replace_all(x, regex("myocardial infarction", ignore_case = TRUE), "MI")
  x <- str_replace_all(x, regex("heart failure, unspecified", ignore_case = TRUE), "heart failure")
  x <- str_replace_all(x, regex("essential \\(primary\\) hypertension", ignore_case = TRUE), "hypertension")
  x <- str_replace_all(x, regex("essential hypertension", ignore_case = TRUE), "hypertension")
  x <- str_replace_all(x, regex("malignant neoplasm of", ignore_case = TRUE), "cancer of")
  x <- str_replace_all(x, regex("personal history of malignant neoplasm", ignore_case = TRUE), "history of cancer")
  x <- str_replace_all(x, regex("dependence on renal dialysis", ignore_case = TRUE), "dialysis dependence")
  x <- str_replace_all(x, regex("extracorporeal dialysis", ignore_case = TRUE), "dialysis")
  x <- str_replace_all(x, regex("preparatory care for dialysis", ignore_case = TRUE), "dialysis preparation")
  x <- str_replace_all(x, regex("kidney transplant status", ignore_case = TRUE), "kidney transplant")
  x <- str_replace_all(x, regex("liver transplant status", ignore_case = TRUE), "liver transplant")
  x <- str_replace_all(x, regex("presence of cardiac pacemaker", ignore_case = TRUE), "pacemaker")
  x <- str_replace_all(x, regex("acquired absence of part of stomach", ignore_case = TRUE), "partial gastrectomy")
  x <- str_replace_all(x, regex("palliative care", ignore_case = TRUE), "palliative care")
  x <- str_replace_all(x, regex("personal history of long-term drug therapy", ignore_case = TRUE), "long-term medication history")
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

classify_icd_organ_domain_initial <- function(icd_clean) {
  x <- clean_icd_code(icd_clean)
  letter <- str_sub(x, 1, 1)
  num <- suppressWarnings(as.integer(str_extract(x, "(?<=^[A-Z])[0-9]{2}")))
  
  case_when(
    letter %in% c("A", "B") ~ "Infectious",
    letter == "C" ~ "Cancer / neoplasm",
    letter == "D" & !is.na(num) & num <= 49 ~ "Cancer / neoplasm",
    letter == "D" & !is.na(num) & num >= 50 ~ "Immune / blood",
    letter == "E" ~ "Endocrine / metabolic",
    letter == "F" ~ "Brain / neurologic",
    letter == "G" ~ "Brain / neurologic",
    letter == "H" ~ "Eye / ENT",
    letter == "I" ~ "Heart / cardiovascular",
    letter == "J" ~ "Pulmonary",
    letter == "K" ~ "Digestive",
    letter == "L" ~ "Skin",
    letter == "M" ~ "Musculoskeletal",
    letter == "N" ~ "Kidney / renal",
    letter == "O" ~ "Reproductive",
    letter == "P" ~ "Other",
    letter == "Q" ~ "Other",
    letter == "R" ~ "Clinical symptoms",
    letter %in% c("S", "T", "V", "W", "X", "Y") ~ "Injury / external",
    letter == "Z" ~ "Other",
    TRUE ~ "Other"
  )
}

assign_special_code_to_organ_domain <- function(icd_clean, disease_text) {
  x <- clean_icd_code(icd_clean)
  txt <- str_to_lower(str_squish(as.character(disease_text)))
  
  case_when(
    # Z-code renal care/status
    str_detect(x, "^Z49") |
      str_detect(x, "^Z940") |
      str_detect(x, "^Z992") |
      str_detect(txt, "dialysis|kidney|renal|nephro") ~ "Kidney / renal",
    
    # Z-code liver status
    str_detect(x, "^Z944") |
      str_detect(txt, "liver|hepatic") ~ "Liver / hepatic",
    
    # Z-code cardiac device/status
    str_detect(x, "^Z950") |
      str_detect(txt, "pacemaker|cardiac|heart|coronary") ~ "Heart / cardiovascular",
    
    # Z-code digestive surgical status
    str_detect(x, "^Z903") |
      str_detect(txt, "stomach|gastrectomy|gastric|intestinal|bowel|digestive") ~ "Digestive",
    
    # Z-code pulmonary/thoracic cancer history
    str_detect(x, "^Z852") |
      str_detect(txt, "thoracic|intrathoracic|lung|respiratory") ~ "Pulmonary",
    
    # Z-code oncology history without single organ
    str_detect(txt, "history of cancer|malignant neoplasm|neoplasm|cancer") ~ "Cancer / neoplasm",
    
    # Z-code reproductive status
    str_detect(txt, "pregnancy|obstetric|maternal|prostate|ovary|uterus|breast|reproductive") ~ "Reproductive",
    
    # External/complication organ reassignment
    str_detect(x, "^W") |
      str_detect(txt, "fall|slip|trip|fracture|gait") ~ "Musculoskeletal",
    
    str_detect(x, "^Y841") |
      str_detect(x, "^T824") |
      str_detect(x, "^T861") ~ "Kidney / renal",
    
    str_detect(x, "^Y545") |
      str_detect(txt, "mineralocorticoid|adrenal|endocrine") ~ "Endocrine / metabolic",
    
    str_detect(x, "^Y495") |
      str_detect(txt, "psychotropic|nervous system|brain") ~ "Brain / neurologic",
    
    str_detect(x, "^Y95") |
      str_detect(txt, "nosocomial|infection") ~ "Infectious",
    
    str_detect(x, "^T179") |
      str_detect(txt, "respiratory tract|airway") ~ "Pulmonary",
    
    TRUE ~ NA_character_
  )
}

# ============================================================
# 8. Disease labels and organ-domain assignment
# ============================================================

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
    manual_icd10_short_map,
    by = "icd_clean"
  ) %>%
  mutate(
    disease_full_name = coalesce(
      clean_name_candidate(disease_full_name_from_file),
      clean_name_candidate(disease_name_from_id),
      clean_name_candidate(disease_name_icd),
      clean_name_candidate(disease_name_icd3),
      clean_name_candidate(icd_pretty)
    ),
    
    disease_short_name_raw = coalesce(
      clean_name_candidate(disease_short_manual),
      shorten_disease_name(disease_full_name),
      clean_name_candidate(disease_full_name),
      icd_pretty
    ),
    
    disease_short_name = if_else(
      is_icd_like(disease_short_name_raw) & !is.na(disease_full_name) & !is_icd_like(disease_full_name),
      shorten_disease_name(disease_full_name),
      disease_short_name_raw
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
    
    endpoint_organ_domain_initial = classify_icd_organ_domain_initial(icd_clean),
    
    endpoint_organ_domain_special = assign_special_code_to_organ_domain(
      icd_clean,
      paste(disease_full_name, disease_short_name, disease_id)
    ),
    
    endpoint_organ_domain = coalesce(
      endpoint_organ_domain_special,
      endpoint_organ_domain_initial,
      "Other"
    ),
    
    endpoint_organ_domain = factor(endpoint_organ_domain, levels = endpoint_organ_levels),
    
    neg_log10_p = -log10(clock_p),
    neg_log10_p_capped = pmin(neg_log10_p, neglog10_cap),
    
    clock_log_hr = case_when(
      !is.na(clock_hr) & clock_hr > 0 ~ log(clock_hr),
      !is.na(clock_beta) ~ as.numeric(clock_beta),
      TRUE ~ NA_real_
    ),
    
    clock_hr_for_plot = case_when(
      !is.na(clock_hr) & clock_hr > 0 ~ clock_hr,
      !is.na(clock_log_hr) ~ exp(clock_log_hr),
      TRUE ~ NA_real_
    )
  )

unmapped_name_tbl <- df %>%
  filter(is_icd_like(disease_short_name)) %>%
  distinct(
    disease_id,
    icd_clean,
    icd_pretty,
    disease_full_name,
    disease_short_name
  ) %>%
  arrange(icd_clean)

readr::write_tsv(
  unmapped_name_tbl,
  paste0(out_prefix, "_unmapped_ICD10_codes_still_shown_as_codes.tsv")
)

# ============================================================
# 9. Same Bonferroni threshold as fair Clock-vs-BAG analysis
# ============================================================

n_unique_diseases <- n_distinct(df$disease_id)
n_mortality_clocks <- n_distinct(df$mortality_clock)

p_clock_bonf <- alpha_bonf / n_unique_diseases / n_mortality_clocks

threshold_tbl <- tibble(
  selection_rule = "clock_p < 0.05 / disease endpoints / mortality EPOCH clocks",
  n_unique_diseases = n_unique_diseases,
  n_mortality_clocks = n_mortality_clocks,
  n_tests = n_unique_diseases * n_mortality_clocks,
  bonferroni_threshold = p_clock_bonf,
  alpha = alpha_bonf,
  analysis = "Mortality EPOCH only",
  endpoint_label_rule = "Word-cloud boxes use shortened disease names; ICD codes are used only if no mapping is available.",
  special_code_rule = "Z-code and organ-informative external/complication endpoints are reassigned to organ domains when possible."
)

print(threshold_tbl)

readr::write_tsv(
  threshold_tbl,
  paste0(out_prefix, "_Bonferroni_threshold.tsv")
)

df_all <- df %>%
  mutate(
    clock_sig_bonf = clock_p < p_clock_bonf,
    case_noncase_label_short = paste0(
      comma(N_case),
      "/",
      comma(N_noncase)
    )
  )

sig_all <- df_all %>%
  filter(clock_sig_bonf) %>%
  arrange(clock_p)

message("Bonferroni threshold: ", signif(p_clock_bonf, 4))
message("Selected mortality EPOCH associations: ", nrow(sig_all))
message("Unique disease endpoints selected: ", n_distinct(sig_all$disease_id))
message("Mortality EPOCH clocks with signal: ", n_distinct(sig_all$mortality_clock), " / ", n_mortality_clocks)

if (nrow(sig_all) == 0) {
  stop("No mortality EPOCH disease associations passed Bonferroni correction.")
}

readr::write_tsv(
  sig_all %>%
    select(
      mortality_clock,
      clock_label,
      clock_block,
      organ,
      organ_display,
      modality,
      modality_group,
      disease_id,
      icd_clean,
      icd_pretty,
      disease_full_name,
      disease_short_name,
      endpoint_organ_domain_initial,
      endpoint_organ_domain_special,
      endpoint_organ_domain,
      source,
      N,
      N_case,
      N_noncase,
      clock_hr_for_plot,
      clock_hr,
      clock_log_hr,
      clock_p,
      neg_log10_p,
      everything()
    ),
  paste0(out_prefix, "_Bonferroni_significant_mortality_EPOCH_only_rows.tsv")
)

special_code_assignment_tbl <- df_all %>%
  filter(
    str_starts(clean_icd_code(icd_clean), "Z") |
      str_detect(clean_icd_code(icd_clean), "^[V-Y]")
  ) %>%
  distinct(
    disease_id,
    icd_clean,
    icd_pretty,
    disease_full_name,
    disease_short_name,
    endpoint_organ_domain_initial,
    endpoint_organ_domain_special,
    endpoint_organ_domain
  ) %>%
  arrange(endpoint_organ_domain, icd_clean)

readr::write_tsv(
  special_code_assignment_tbl,
  paste0(out_prefix, "_Z_and_external_code_organ_domain_assignment.tsv")
)

# ============================================================
# 10. Clock-level and endpoint-domain summaries
# ============================================================

clock_order_tbl <- df_all %>%
  distinct(mortality_clock, clock_block, clock_label, modality_group, modality_rank) %>%
  mutate(
    clock_block_chr = as.character(clock_block),
    clock_block_order = match(clock_block_chr, organ_levels),
    clock_block_order = if_else(is.na(clock_block_order), 999L, clock_block_order)
  ) %>%
  arrange(clock_block_order, modality_rank, clock_label, mortality_clock) %>%
  mutate(clock_plot_order = row_number())

clock_summary <- df_all %>%
  group_by(mortality_clock, clock_block, clock_label, modality_group) %>%
  summarise(
    n_tested_endpoints = n_distinct(disease_id),
    n_sig_endpoints = sum(clock_sig_bonf, na.rm = TRUE),
    n_sig_endpoint_organ_domains = n_distinct(endpoint_organ_domain[clock_sig_bonf], na.rm = TRUE),
    min_clock_p = min(clock_p, na.rm = TRUE),
    max_neg_log10_p = max(neg_log10_p, na.rm = TRUE),
    median_hr_sig = median(clock_hr_for_plot[clock_sig_bonf], na.rm = TRUE),
    max_hr_sig = max(clock_hr_for_plot[clock_sig_bonf], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    clock_order_tbl %>%
      select(mortality_clock, clock_plot_order),
    by = "mortality_clock"
  ) %>%
  mutate(
    clock_block_chr = as.character(clock_block),
    clock_label_plot = paste0(clock_label, " EPOCH"),
    clock_label_short = str_replace_all(clock_label_plot, " protein", "\nprotein"),
    clock_label_short = str_replace_all(clock_label_short, " metabolite", "\nmetabolite"),
    clock_label_short = str_replace_all(clock_label_short, " MRI", "\nMRI"),
    n_sig_label = paste0(n_sig_endpoints, " endpoints")
  ) %>%
  arrange(clock_plot_order)

summary_by_clock_out <- paste0(out_prefix, "_summary_by_mortality_EPOCH_clock.tsv")

readr::write_tsv(
  clock_summary,
  summary_by_clock_out
)

endpoint_domain_summary <- sig_all %>%
  count(endpoint_organ_domain, name = "n_significant_clock_endpoint_pairs") %>%
  mutate(endpoint_organ_domain = as.character(endpoint_organ_domain)) %>%
  arrange(desc(n_significant_clock_endpoint_pairs))

readr::write_tsv(
  endpoint_domain_summary,
  paste0(out_prefix, "_summary_by_endpoint_organ_domain.tsv")
)

# ============================================================
# 11. Human anatomy-style coordinates
# ============================================================

body_coords <- tribble(
  ~clock_block_chr,       ~x0,   ~y0,
  "Brain",                0.00,  5.95,
  "Eye",                  0.32,  5.70,
  "Pulmonary",           -0.48,  4.65,
  "Heart",                0.42,  4.45,
  "Liver / hepatic",     -0.45,  3.75,
  "Kidney / renal",       0.48,  3.35,
  "Pancreas",             0.05,  3.70,
  "Spleen",              -0.82,  3.72,
  "Immune",              -0.92,  2.82,
  "Endocrine",            0.02,  5.10,
  "Digestive",            0.00,  2.82,
  "Metabolic",            0.00,  2.15,
  "Adipose",              0.00,  1.65,
  "Skin",                 0.95,  3.20,
  "Reproductive",         0.00,  1.15,
  "Other",               -1.20,  0.75
)

clock_summary <- clock_summary %>%
  left_join(body_coords, by = "clock_block_chr") %>%
  group_by(clock_block_chr) %>%
  arrange(clock_plot_order, .by_group = TRUE) %>%
  mutate(
    n_in_block = n(),
    index_in_block = row_number(),
    x = x0 + (index_in_block - (n_in_block + 1) / 2) * 0.19,
    y = y0 + if_else(index_in_block %% 2 == 0, 0.08, -0.08)
  ) %>%
  ungroup()

domain_nodes <- sig_all %>%
  count(endpoint_organ_domain, name = "n_pairs") %>%
  mutate(endpoint_organ_domain = as.character(endpoint_organ_domain)) %>%
  filter(!is.na(endpoint_organ_domain), endpoint_organ_domain != "") %>%
  mutate(
    endpoint_organ_domain = factor(endpoint_organ_domain, levels = endpoint_organ_levels)
  ) %>%
  arrange(endpoint_organ_domain) %>%
  mutate(
    x_domain = 3.45,
    y_domain = seq(5.85, 0.85, length.out = n())
  )

edge_df <- sig_all %>%
  group_by(mortality_clock, clock_block, endpoint_organ_domain) %>%
  summarise(
    n_endpoint = n_distinct(disease_id),
    n_pairs = n(),
    min_p = min(clock_p, na.rm = TRUE),
    max_neg_log10_p = max(neg_log10_p, na.rm = TRUE),
    median_hr = median(clock_hr_for_plot, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_endpoint >= MIN_ENDPOINT_PER_EDGE) %>%
  left_join(
    clock_summary %>%
      select(mortality_clock, clock_label_plot, clock_block_chr, x, y),
    by = "mortality_clock"
  ) %>%
  mutate(endpoint_organ_domain = as.character(endpoint_organ_domain)) %>%
  left_join(
    domain_nodes %>%
      mutate(endpoint_organ_domain = as.character(endpoint_organ_domain)),
    by = "endpoint_organ_domain"
  ) %>%
  filter(!is.na(x), !is.na(y), !is.na(x_domain), !is.na(y_domain)) %>%
  mutate(
    edge_alpha = pmin(0.75, 0.18 + 0.06 * log1p(max_neg_log10_p)),
    edge_width = pmin(2.2, 0.18 + 0.22 * sqrt(n_endpoint))
  )

readr::write_tsv(
  edge_df,
  paste0(out_prefix, "_whole_body_network_edges.tsv")
)

# ============================================================
# 12. Panel A: whole-body risk architecture
# ============================================================

body_col <- "#F3F3F3"
body_outline <- "#D5D5D5"

p_architecture <- ggplot() +
  
  annotate(
    "point",
    x = 0,
    y = 5.75,
    size = 36,
    shape = 21,
    fill = body_col,
    color = body_outline,
    stroke = 0.6
  ) +
  annotate(
    "segment",
    x = 0,
    xend = 0,
    y = 4.95,
    yend = 2.05,
    linewidth = 34,
    lineend = "round",
    color = body_col
  ) +
  annotate(
    "segment",
    x = -0.20,
    xend = -1.22,
    y = 4.55,
    yend = 3.15,
    linewidth = 14,
    lineend = "round",
    color = body_col
  ) +
  annotate(
    "segment",
    x = 0.20,
    xend = 1.22,
    y = 4.55,
    yend = 3.15,
    linewidth = 14,
    lineend = "round",
    color = body_col
  ) +
  annotate(
    "segment",
    x = -0.22,
    xend = -0.62,
    y = 2.05,
    yend = 0.55,
    linewidth = 16,
    lineend = "round",
    color = body_col
  ) +
  annotate(
    "segment",
    x = 0.22,
    xend = 0.62,
    y = 2.05,
    yend = 0.55,
    linewidth = 16,
    lineend = "round",
    color = body_col
  ) +
  
  geom_curve(
    data = edge_df,
    aes(
      x = x,
      y = y,
      xend = x_domain,
      yend = y_domain,
      color = clock_block,
      linewidth = n_endpoint,
      alpha = max_neg_log10_p
    ),
    curvature = 0.15,
    lineend = "round"
  ) +
  
  geom_point(
    data = domain_nodes,
    aes(
      x = x_domain,
      y = y_domain,
      size = n_pairs
    ),
    shape = 21,
    fill = "white",
    color = "grey20",
    stroke = 0.45
  ) +
  
  geom_label(
    data = domain_nodes,
    aes(
      x = x_domain + 0.16,
      y = y_domain,
      label = paste0(endpoint_organ_domain, " (", n_pairs, ")")
    ),
    hjust = 0,
    size = 2.55,
    label.size = 0.15,
    label.padding = unit(0.10, "lines"),
    fill = alpha("white", 0.90),
    color = "grey10"
  ) +
  
  geom_point(
    data = clock_summary,
    aes(
      x = x,
      y = y,
      fill = clock_block,
      size = n_sig_endpoints,
      shape = modality_group
    ),
    color = "grey10",
    stroke = 0.35,
    alpha = 0.96
  ) +
  
  geom_label(
    data = clock_summary,
    aes(
      x = x,
      y = y + 0.31,
      label = paste0(clock_label_short, "\n", n_sig_endpoints, " sig.")
    ),
    size = 2.15,
    lineheight = 0.86,
    label.size = 0.12,
    label.padding = unit(0.08, "lines"),
    fill = alpha("white", 0.82),
    color = "grey10"
  ) +
  
  annotate(
    "text",
    x = 0,
    y = 6.48,
    label = "Mortality EPOCH clocks",
    size = 4.2,
    fontface = "bold",
    color = "grey10"
  ) +
  
  annotate(
    "text",
    x = 3.45,
    y = 6.48,
    label = "Disease endpoint organ domains",
    size = 4.2,
    fontface = "bold",
    color = "grey10"
  ) +
  
  scale_color_manual(
    values = clock_block_palette,
    drop = FALSE,
    name = "EPOCH organ/system"
  ) +
  scale_fill_manual(
    values = clock_block_palette,
    drop = FALSE,
    name = "EPOCH organ/system"
  ) +
  scale_shape_manual(
    values = modality_shape_values,
    drop = FALSE,
    name = "Modality"
  ) +
  scale_size_area(
    max_size = 10,
    name = "No. significant\nendpoints"
  ) +
  scale_linewidth_continuous(
    range = c(0.15, 1.35),
    name = "Endpoints per\norgan-domain edge"
  ) +
  scale_alpha_continuous(
    range = c(0.10, 0.55),
    guide = "none"
  ) +
  
  coord_cartesian(
    xlim = c(-1.75, 5.35),
    ylim = c(0.15, 6.65),
    clip = "off"
  ) +
  
  labs(
    tag = "A",
    title = "Brain-body mortality EPOCH risk architecture",
    subtitle = paste0(
      "Bonferroni-significant disease associations: clock P < ",
      signif(p_clock_bonf, 3),
      " = 0.05 / ",
      n_unique_diseases,
      " disease endpoints / ",
      n_mortality_clocks,
      " mortality EPOCH clocks"
    ),
    x = NULL,
    y = NULL
  ) +
  
  guides(
    color = guide_legend(override.aes = list(linewidth = 1.5, alpha = 1)),
    fill = "none",
    linewidth = guide_legend(order = 3),
    size = guide_legend(order = 2),
    shape = guide_legend(order = 1)
  ) +
  
  theme_void(base_size = 9) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.tag = element_text(face = "bold", size = 22),
    plot.tag.position = c(0.01, 0.99),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 8.5, color = "grey30", lineheight = 1.05),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(face = "bold", size = 8),
    legend.text = element_text(size = 7),
    plot.margin = margin(8, 16, 8, 8)
  )

# ============================================================
# 13. Panel B: word-cloud boxes with shortened disease names
# ============================================================

dominant_clock_for_endpoint <- sig_all %>%
  group_by(endpoint_organ_domain, disease_id, disease_short_name) %>%
  arrange(clock_p, desc(abs(clock_log_hr)), .by_group = TRUE) %>%
  summarise(
    dominant_clock_block = as.character(clock_block[1]),
    dominant_clock_label = clock_label[1],
    dominant_mortality_clock = mortality_clock[1],
    min_p = min(clock_p, na.rm = TRUE),
    max_neg_log10_p = max(neg_log10_p, na.rm = TRUE),
    n_clock_endpoint_pairs = n(),
    n_unique_clocks = n_distinct(mortality_clock),
    median_hr = median(clock_hr_for_plot, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    endpoint_organ_domain = as.character(endpoint_organ_domain),
    word = disease_short_name,
    word = if_else(is.na(word) | word == "", "Unnamed endpoint", word),
    word_weight_raw = log1p(n_clock_endpoint_pairs) + 0.25 * pmin(max_neg_log10_p, neglog10_cap)
  ) %>%
  group_by(endpoint_organ_domain) %>%
  arrange(desc(word_weight_raw), min_p, .by_group = TRUE) %>%
  slice_head(n = MAX_WORDS_PER_ENDPOINT_ORGAN_DOMAIN) %>%
  mutate(
    word_rank = row_number(),
    n_words_domain = n()
  ) %>%
  ungroup()

safe_rescale <- function(x, to = c(2.6, 7.6)) {
  if (length(x) == 0 || all(is.na(x))) {
    return(rep(mean(to), length(x)))
  }
  
  rng <- range(x, na.rm = TRUE)
  
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) {
    return(rep(mean(to), length(x)))
  }
  
  scales::rescale(x, to = to, from = rng)
}

dominant_clock_for_endpoint <- dominant_clock_for_endpoint %>%
  mutate(
    word_size = safe_rescale(word_weight_raw, to = c(2.6, 7.6)),
    word_alpha = safe_rescale(pmin(max_neg_log10_p, neglog10_cap), to = c(0.65, 1.00))
  )

if (nrow(dominant_clock_for_endpoint) == 0) {
  stop("No endpoint words available for word cloud.")
}

make_domain_word_layout <- function(n) {
  base_pos <- tibble(
    x = c(
      0.00, -0.42, 0.42, -0.05, 0.05,
      -0.68, 0.68, -0.28, 0.28, -0.55, 0.55,
      -0.82, 0.82, -0.12, 0.12, -0.40, 0.40,
      -0.72, 0.72, -0.22, 0.22, -0.58, 0.58, 0.00
    ),
    y = c(
      0.00, 0.22, -0.22, -0.34, 0.34,
      0.00, 0.00, 0.48, -0.48, -0.42, 0.42,
      0.28, -0.28, 0.66, -0.66, -0.70, 0.70,
      -0.58, 0.58, 0.82, -0.82, 0.76, -0.76, 0.92
    )
  )
  
  if (n <= nrow(base_pos)) {
    return(base_pos[seq_len(n), , drop = FALSE])
  }
  
  extra_n <- n - nrow(base_pos)
  theta <- seq(0, 4.2 * pi, length.out = extra_n)
  r <- seq(0.25, 0.95, length.out = extra_n)
  
  extra <- tibble(
    x = r * cos(theta),
    y = 0.90 * r * sin(theta)
  )
  
  bind_rows(base_pos, extra)
}

word_df <- dominant_clock_for_endpoint %>%
  group_by(endpoint_organ_domain) %>%
  group_modify(~ {
    pos <- make_domain_word_layout(nrow(.x))
    bind_cols(.x, pos)
  }) %>%
  ungroup() %>%
  mutate(
    endpoint_organ_domain = factor(endpoint_organ_domain, levels = endpoint_organ_levels)
  )

domain_box_df <- word_df %>%
  distinct(endpoint_organ_domain) %>%
  mutate(
    xmin = -1.05,
    xmax = 1.05,
    ymin = -1.05,
    ymax = 1.05
  )

readr::write_tsv(
  word_df %>%
    select(
      endpoint_organ_domain,
      word,
      disease_id,
      disease_short_name,
      dominant_clock_block,
      dominant_clock_label,
      dominant_mortality_clock,
      n_clock_endpoint_pairs,
      n_unique_clocks,
      min_p,
      max_neg_log10_p,
      median_hr,
      x,
      y,
      word_size,
      word_alpha
    ),
  paste0(out_prefix, "_endpoint_organ_domain_wordcloud_terms_short_names.tsv")
)

p_wordcloud <- ggplot() +
  
  geom_rect(
    data = domain_box_df,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax
    ),
    fill = "white",
    color = "grey78",
    linewidth = 0.35
  ) +
  
  geom_text(
    data = word_df,
    aes(
      x = x,
      y = y,
      label = word,
      size = word_size,
      color = dominant_clock_block,
      alpha = word_alpha
    ),
    fontface = "bold",
    lineheight = 0.86,
    check_overlap = TRUE
  ) +
  
  facet_wrap(
    ~ endpoint_organ_domain,
    ncol = 3,
    drop = TRUE
  ) +
  
  scale_color_manual(
    values = clock_block_palette,
    drop = FALSE,
    name = "Dominant mortality\nEPOCH organ/system"
  ) +
  
  scale_size_identity() +
  scale_alpha_identity() +
  
  coord_cartesian(
    xlim = c(-1.08, 1.08),
    ylim = c(-1.08, 1.08),
    clip = "off"
  ) +
  
  labs(
    tag = "B",
    title = "Disease endpoint word clouds by organ domain",
    subtitle = "Words show shortened disease names; color denotes the dominant mortality EPOCH organ/system."
  ) +
  
  theme_void(base_size = 9) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(face = "bold", size = 8.2, color = "grey10"),
    strip.background = element_rect(fill = "grey96", color = "grey80", linewidth = 0.25),
    panel.spacing = unit(0.55, "lines"),
    plot.tag = element_text(face = "bold", size = 22),
    plot.tag.position = c(0.01, 0.99),
    plot.title = element_text(face = "bold", size = 13.5),
    plot.subtitle = element_text(size = 8.2, color = "grey30", lineheight = 1.05),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 8),
    legend.text = element_text(size = 7),
    plot.margin = margin(8, 8, 8, 8)
  )

# ============================================================
# 14. Panel C: significant endpoints per mortality EPOCH clock
# ============================================================

clock_levels <- clock_summary %>%
  arrange(clock_plot_order) %>%
  pull(clock_label_plot)

p_clock_bar <- clock_summary %>%
  mutate(
    clock_label_plot = factor(clock_label_plot, levels = rev(clock_levels))
  ) %>%
  ggplot(
    aes(
      x = n_sig_endpoints,
      y = clock_label_plot,
      fill = clock_block
    )
  ) +
  geom_col(width = 0.72, color = "white", linewidth = 0.20) +
  geom_text(
    aes(label = n_sig_endpoints),
    hjust = -0.12,
    size = 2.4,
    color = "grey10",
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = clock_block_palette,
    drop = FALSE,
    guide = "none"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.01, 0.15)),
    labels = comma
  ) +
  labs(
    tag = "C",
    title = "Significant disease endpoints per mortality EPOCH clock",
    x = "No. Bonferroni-significant disease endpoints",
    y = NULL
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 7.2, color = "grey15"),
    axis.text.x = element_text(size = 7.5, color = "grey15"),
    axis.title.x = element_text(face = "bold", size = 9),
    plot.tag = element_text(face = "bold", size = 20),
    plot.title = element_text(face = "bold", size = 12.5),
    plot.margin = margin(8, 8, 8, 8)
  )

# ============================================================
# 15. Combined figure
# ============================================================

p_combine <- (
  p_architecture | p_wordcloud
) / p_clock_bar +
  patchwork::plot_layout(
    widths = c(1.08, 1.12),
    heights = c(1.00, 0.42),
    guides = "collect"
  ) +
  patchwork::plot_annotation(
    title = "Mortality EPOCH reveals a brain-body, whole-body risk architecture for future disease onset",
    subtitle = paste0(
      "Only mortality EPOCH associations are shown. Significance is defined using the same Bonferroni correction: ",
      "P < 0.05 / disease endpoints / mortality EPOCH clocks = ",
      signif(p_clock_bonf, 3),
      "."
    ),
    caption = paste0(
      "Nodes represent organ- or system-specific mortality EPOCH clocks positioned on a simplified human anatomy layout. ",
      "Curves connect each clock to organ-assigned disease endpoint domains when at least one Bonferroni-significant association is present. ",
      "Panel B shows shortened disease endpoint names grouped by endpoint organ domain; word color denotes the dominant mortality EPOCH organ/system. ",
      "Z-code and organ-informative external/complication endpoints were reassigned to organ domains when possible; otherwise they were assigned to Other."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 17, hjust = 0),
      plot.subtitle = element_text(size = 9.5, color = "grey25", hjust = 0),
      plot.caption = element_text(size = 7.5, color = "grey35", hjust = 0, lineheight = 1.05),
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

print(p_combine)

# ============================================================
# 16. Save outputs
# ============================================================

combined_pdf <- paste0(out_prefix, "_combined_whole_body_architecture_wordcloud_short_names.pdf")
combined_png <- paste0(out_prefix, "_combined_whole_body_architecture_wordcloud_short_names.png")
combined_rds <- paste0(out_prefix, "_combined_whole_body_architecture_wordcloud_short_names_p_combine.rds")

ggsave(
  filename = combined_pdf,
  plot = p_combine,
  width = 19.0,
  height = 15.8,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  filename = combined_png,
  plot = p_combine,
  width = 19.0,
  height = 15.8,
  units = "in",
  dpi = 500,
  bg = "white"
)

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    filename = paste0(out_prefix, "_combined_whole_body_architecture_wordcloud_short_names.svg"),
    plot = p_combine,
    width = 19.0,
    height = 15.8,
    units = "in",
    device = svglite::svglite,
    bg = "white"
  )
}

saveRDS(
  p_combine,
  file = combined_rds
)

ggsave(
  filename = paste0(out_prefix, "_panelA_anatomy_network.pdf"),
  plot = p_architecture,
  width = 10.5,
  height = 9.5,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  filename = paste0(out_prefix, "_panelB_endpoint_organ_domain_wordcloud_short_names.pdf"),
  plot = p_wordcloud,
  width = 10.7,
  height = 9.5,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  filename = paste0(out_prefix, "_panelC_clock_count_barplot.pdf"),
  plot = p_clock_bar,
  width = 9.5,
  height = 6.4,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

readr::write_tsv(
  clock_summary,
  paste0(out_prefix, "_clock_summary.tsv")
)

message("============================================================")
message("Finished mortality EPOCH-only whole-body architecture plot with shortened disease-name word clouds.")
message("Input file:")
message("  ", infile)
message("Output directory:")
message("  ", out_dir)
message("Bonferroni threshold:")
message("  ", signif(p_clock_bonf, 5))
message("Main figure:")
message("  ", combined_pdf)
message("  ", combined_png)
message("R object:")
message("  ", combined_rds)
message("Still-unmapped ICD10 codes:")
message("  ", paste0(out_prefix, "_unmapped_ICD10_codes_still_shown_as_codes.tsv"))
message("Special code assignments:")
message("  ", paste0(out_prefix, "_Z_and_external_code_organ_domain_assignment.tsv"))
message("Word-cloud terms:")
message("  ", paste0(out_prefix, "_endpoint_organ_domain_wordcloud_terms_short_names.tsv"))
message("Selected rows:")
message("  ", paste0(out_prefix, "_Bonferroni_significant_mortality_EPOCH_only_rows.tsv"))
message("Clock summary:")
message("  ", summary_by_clock_out)
message("============================================================")