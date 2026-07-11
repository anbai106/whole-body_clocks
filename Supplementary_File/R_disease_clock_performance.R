#!/usr/bin/env Rscript

# ============================================================
# Collect model performance summary for 5 disease EPOCH clocks
# Diseases: dementia, COPD, asthma, MI, stroke
#
# Hard-coded 23 organ/modality clock folders per disease.
# Does NOT scan the full base directory.
#
# Input:
#   <base_dir>/<hard-coded disease clock folder>/<prefix>_performance.json
#
# Output:
#   /Users/hao/Dropbox/2026_EPOCH/Supplementary_File/
#     Supplementary_Table_Disease_EPOCH_Model_Performance_5_diseases.xlsx
#
# Excel sheets:
#   1) Summary_5_disease_clocks
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(jsonlite)
  library(openxlsx)
  library(glue)
})

# ============================================================
# 1. Paths
# ============================================================

base_dir_candidates <- c(
  "/gpfs/fs001/cbica/home/wenju/Reproducibile_paper/WholeBodyClock",
  "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock",
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock",
  getwd()
)

base_dir <- Sys.getenv("WHOLEBODYCLOCK_BASE_DIR", unset = NA_character_)

if (is.na(base_dir) || base_dir == "") {
  base_dir <- base_dir_candidates[dir.exists(base_dir_candidates)][1]
}

if (is.na(base_dir) || !dir.exists(base_dir)) {
  stop(
    "Could not find base_dir. Please set it manually, for example:\n",
    "base_dir <- '/gpfs/fs001/cbica/home/wenju/Reproducibile_paper/WholeBodyClock'"
  )
}

base_dir <- normalizePath(base_dir, mustWork = TRUE)

out_dir <- "/Users/hao/Dropbox/2026_EPOCH/Supplementary_File"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_xlsx <- file.path(
  out_dir,
  "Supplementary_Table_Disease_EPOCH_Model_Performance_5_diseases.xlsx"
)

message("Base directory:")
message("  ", base_dir)
message("Output Excel:")
message("  ", out_xlsx)

# ============================================================
# 2. Hard-coded disease list
# ============================================================

disease_tbl <- tribble(
  ~disease_key, ~disease_label,
  "dementia",   "Dementia",
  "copd",       "COPD",
  "asthma",     "Asthma",
  "mi",         "MI",
  "stroke",     "Stroke"
)

# ============================================================
# 3. Hard-coded 23 organ/modality clock templates
# ============================================================

clock_template <- tribble(
  ~modality,       ~modality_key,    ~organ_key,              ~organ_label,              ~folder_organ_name,        ~prefix_organ_name,
  "MRI",           "mri",            "adipose",               "Adipose",                 "adipose",                "adipose",
  "MRI",           "mri",            "brain",                 "Brain",                   "brain",                  "brain",
  "MRI",           "mri",            "heart",                 "Heart",                   "heart",                  "heart",
  "MRI",           "mri",            "kidney",                "Kidney",                  "kidney",                 "kidney",
  "MRI",           "mri",            "liver",                 "Liver",                   "liver",                  "liver",
  "MRI",           "mri",            "pancreas",              "Pancreas",                "pancreas",               "pancreas",
  "MRI",           "mri",            "spleen",                "Spleen",                  "spleen",                 "spleen",

  "Proteomics",    "proteomics",     "brain",                 "Brain",                   "Brain",                  "brain",
  "Proteomics",    "proteomics",     "endocrine",             "Endocrine",               "Endocrine",              "endocrine",
  "Proteomics",    "proteomics",     "eye",                   "Eye",                     "Eye",                    "eye",
  "Proteomics",    "proteomics",     "heart",                 "Heart",                   "Heart",                  "heart",
  "Proteomics",    "proteomics",     "hepatic",               "Hepatic",                 "Hepatic",                "hepatic",
  "Proteomics",    "proteomics",     "immune",                "Immune",                  "Immune",                 "immune",
  "Proteomics",    "proteomics",     "pulmonary",             "Pulmonary",               "Pulmonary",              "pulmonary",
  "Proteomics",    "proteomics",     "renal",                 "Renal",                   "Renal",                  "renal",
  "Proteomics",    "proteomics",     "reproductive_female",   "Reproductive female",     "Reproductive_female",    "reproductive_female",
  "Proteomics",    "proteomics",     "reproductive_male",     "Reproductive male",       "Reproductive_male",      "reproductive_male",
  "Proteomics",    "proteomics",     "skin",                  "Skin",                    "Skin",                   "skin",

  "Metabolomics",  "metabolomics",   "digestive",             "Digestive",               "Digestive",              "digestive",
  "Metabolomics",  "metabolomics",   "endocrine",             "Endocrine",               "Endocrine",              "endocrine",
  "Metabolomics",  "metabolomics",   "hepatic",               "Hepatic",                 "Hepatic",                "hepatic",
  "Metabolomics",  "metabolomics",   "immune",                "Immune",                  "Immune",                 "immune",
  "Metabolomics",  "metabolomics",   "metabolic",             "Metabolic",               "Metabolic",              "metabolic"
)

# ============================================================
# 4. Build hard-coded 5 x 23 manifest
# ============================================================

clock_manifest <- tidyr::crossing(
  disease_tbl,
  clock_template
) %>%
  mutate(
    folder = paste0(
      folder_organ_name,
      "_",
      modality_key,
      "_",
      disease_key,
      "_clock"
    ),
    prefix = paste0(
      prefix_organ_name,
      "_",
      modality_key,
      "_",
      disease_key,
      "_clock"
    ),
    clock_id = paste(disease_key, organ_key, modality_key, sep = "__"),
    disease_clock_label = paste(disease_label, "EPOCH"),
    clock_label = paste(disease_label, organ_label, modality),
    organ_modality_label = paste(organ_label, modality),
    json_file = file.path(base_dir, folder, paste0(prefix, "_performance.json")),
    json_file_relative = file.path(folder, paste0(prefix, "_performance.json")),
    json_exists = file.exists(json_file)
  ) %>%
  arrange(
    factor(disease_key, levels = disease_tbl$disease_key),
    factor(modality, levels = c("MRI", "Proteomics", "Metabolomics")),
    organ_label
  )

message("Hard-coded expected disease-clock JSON files: ", nrow(clock_manifest))
message("Existing JSON files: ", sum(clock_manifest$json_exists), " / ", nrow(clock_manifest))

missing_files <- clock_manifest %>% filter(!json_exists)

if (nrow(missing_files) > 0) {
  warning(
    "Some hard-coded JSON files are missing:\n",
    paste0("  ", missing_files$json_file_relative, collapse = "\n")
  )
}

# ============================================================
# 5. Helper functions
# ============================================================

safe_scalar <- function(x) {
  if (is.null(x)) return(NA)
  if (length(x) == 0) return(NA)
  if (is.atomic(x) && length(x) == 1) return(x)
  return(NA)
}

safe_numeric_scalar <- function(x) {
  suppressWarnings(as.numeric(safe_scalar(x)))
}

safe_character_scalar <- function(x) {
  out <- safe_scalar(x)
  if (all(is.na(out))) return(NA_character_)
  as.character(out)
}

safe_logical_scalar <- function(x) {
  out <- safe_scalar(x)
  if (all(is.na(out))) return(NA)
  as.logical(out)
}

format_p_display <- function(p) {
  case_when(
    is.na(p) ~ NA_character_,
    p == 0 ~ "<1e-300",
    p < 0.001 ~ formatC(p, format = "e", digits = 2),
    TRUE ~ as.character(signif(p, 3))
  )
}

read_json_safe <- function(path) {
  tryCatch(
    jsonlite::fromJSON(path, simplifyDataFrame = FALSE),
    error = function(e) {
      warning("Failed to read JSON: ", path, "\n", e$message)
      NULL
    }
  )
}

list_rows_to_tibble <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(tibble())
  }

  purrr::map_dfr(x, function(row_i) {
    if (is.null(row_i)) return(tibble())
    as_tibble(as.list(row_i))
  })
}

extract_test_model_cindex <- function(model_tbl, model_prefix) {
  if (nrow(model_tbl) == 0) return(NA_real_)
  if (!all(c("split", "model", "cindex") %in% colnames(model_tbl))) return(NA_real_)

  out <- model_tbl %>%
    filter(
      split == "test",
      stringr::str_starts(model, model_prefix)
    ) %>%
    slice(1)

  if (nrow(out) == 0) return(NA_real_)
  suppressWarnings(as.numeric(out$cindex[[1]]))
}

extract_test_model_features <- function(model_tbl, model_prefix) {
  if (nrow(model_tbl) == 0) return(NA_real_)
  if (!all(c("split", "model", "n_features") %in% colnames(model_tbl))) return(NA_real_)

  out <- model_tbl %>%
    filter(
      split == "test",
      stringr::str_starts(model, model_prefix)
    ) %>%
    slice(1)

  if (nrow(out) == 0) return(NA_real_)
  suppressWarnings(as.numeric(out$n_features[[1]]))
}

extract_test_model_label <- function(model_tbl, model_prefix) {
  if (nrow(model_tbl) == 0) return(NA_character_)
  if (!all(c("split", "model", "model_label") %in% colnames(model_tbl))) return(NA_character_)

  out <- model_tbl %>%
    filter(
      split == "test",
      stringr::str_starts(model, model_prefix)
    ) %>%
    slice(1)

  if (nrow(out) == 0) return(NA_character_)
  as.character(out$model_label[[1]])
}

extract_delta_summary <- function(delta_tbl) {
  if (nrow(delta_tbl) == 0) {
    return(tibble(
      delta_comparison = NA_character_,
      delta_cindex = NA_real_,
      delta_cindex_ci_lower = NA_real_,
      delta_cindex_ci_upper = NA_real_,
      delta_p_two_sided = NA_real_,
      delta_p_one_sided_le_0 = NA_real_,
      delta_n_bootstrap_requested = NA_real_,
      delta_n_bootstrap_successful = NA_real_,
      delta_interpretation = NA_character_
    ))
  }

  d <- delta_tbl %>% slice(1)

  tibble(
    delta_comparison = if ("comparison" %in% colnames(d)) as.character(d$comparison[[1]]) else NA_character_,
    delta_cindex = if ("delta_cindex" %in% colnames(d)) as.numeric(d$delta_cindex[[1]]) else NA_real_,
    delta_cindex_ci_lower = if ("delta_cindex_ci_lower" %in% colnames(d)) as.numeric(d$delta_cindex_ci_lower[[1]]) else NA_real_,
    delta_cindex_ci_upper = if ("delta_cindex_ci_upper" %in% colnames(d)) as.numeric(d$delta_cindex_ci_upper[[1]]) else NA_real_,
    delta_p_two_sided = if ("empirical_p_two_sided_delta_not_equal_0" %in% colnames(d)) {
      as.numeric(d$empirical_p_two_sided_delta_not_equal_0[[1]])
    } else {
      NA_real_
    },
    delta_p_one_sided_le_0 = if ("empirical_p_one_sided_delta_le_0" %in% colnames(d)) {
      as.numeric(d$empirical_p_one_sided_delta_le_0[[1]])
    } else {
      NA_real_
    },
    delta_n_bootstrap_requested = if ("n_bootstrap_requested" %in% colnames(d)) {
      as.numeric(d$n_bootstrap_requested[[1]])
    } else {
      NA_real_
    },
    delta_n_bootstrap_successful = if ("n_bootstrap_successful" %in% colnames(d)) {
      as.numeric(d$n_bootstrap_successful[[1]])
    } else {
      NA_real_
    },
    delta_interpretation = if ("interpretation" %in% colnames(d)) as.character(d$interpretation[[1]]) else NA_character_
  )
}

coalesce_numeric <- function(...) {
  vals <- list(...)
  for (v in vals) {
    v_num <- suppressWarnings(as.numeric(v))
    if (length(v_num) > 0 && is.finite(v_num[[1]])) {
      return(v_num[[1]])
    }
  }
  NA_real_
}

# ============================================================
# 6. Parse one JSON into summary-row only
# ============================================================

parse_one_summary <- function(meta_row) {
  meta <- meta_row

  if (!isTRUE(meta$json_exists)) {
    return(meta %>%
      transmute(
        disease_key,
        disease_label,
        disease_clock_label,
        clock_id,
        modality,
        organ_label,
        organ_modality_label,
        clock_label,
        folder,
        prefix,
        json_file_relative,
        json_exists,
        parse_status = "missing_json",

        n_total = NA_real_,
        n_events_total = NA_real_,
        n_censored_total = NA_real_,
        median_followup_years = NA_real_,

        n_train = NA_real_,
        n_events_train = NA_real_,
        n_validation = NA_real_,
        n_events_validation = NA_real_,
        n_test = NA_real_,
        n_events_test = NA_real_,

        cindex_train = NA_real_,
        cindex_validation = NA_real_,
        cindex_trainval = NA_real_,
        cindex_test = NA_real_,

        cindex_test_M0_age_sex = NA_real_,
        cindex_test_M1_covariate_baseline = NA_real_,
        cindex_test_M2_features_only = NA_real_,
        cindex_test_M3_full_model = NA_real_,

        delta_cindex_test_M3_vs_M1 = NA_real_,
        delta_cindex_test_M3_vs_M1_ci_lower = NA_real_,
        delta_cindex_test_M3_vs_M1_ci_upper = NA_real_,
        delta_cindex_test_M3_vs_M1_p_two_sided = NA_real_,
        delta_cindex_test_M3_vs_M1_p_two_sided_display = NA_character_,
        delta_cindex_test_M3_vs_M1_p_one_sided_le_0 = NA_real_,
        delta_n_bootstrap_requested = NA_real_,
        delta_n_bootstrap_successful = NA_real_,
        delta_significant = NA,

        train_minus_test_cindex = NA_real_,
        validation_minus_test_cindex = NA_real_,

        m0_n_features = NA_real_,
        m1_n_features = NA_real_,
        m2_n_features = NA_real_,
        m3_n_features = NA_real_,
        m2_model_label = NA_character_,
        m3_model_label = NA_character_,

        best_l1_ratio = NA_real_,
        best_alpha = NA_real_,
        best_validation_cindex_during_tuning = NA_real_,
        used_penalty_factor = NA,

        n_original_features = NA_real_,
        n_numeric_cols_kept = NA_real_,
        n_categorical_cols_kept = NA_real_,
        n_nonzero_coefficients = NA_real_,
        n_residualization_covariates = NA_real_,

        admin_censor_date = NA_character_,
        time_zero = NA_character_,
        event_date = NA_character_,
        note = NA_character_
      ))
  }

  message("Reading: ", meta$json_file_relative)

  js <- read_json_safe(meta$json_file)

  if (is.null(js)) {
    return(meta %>%
      transmute(
        disease_key,
        disease_label,
        disease_clock_label,
        clock_id,
        modality,
        organ_label,
        organ_modality_label,
        clock_label,
        folder,
        prefix,
        json_file_relative,
        json_exists,
        parse_status = "json_read_failed"
      ))
  }

  model_tbl <- list_rows_to_tibble(js$incremental_value_model_comparison)
  delta_tbl <- list_rows_to_tibble(js$incremental_value_delta_cindex)
  delta_summary <- extract_delta_summary(delta_tbl)

  cindex_test_M0 <- extract_test_model_cindex(model_tbl, "M0")
  cindex_test_M1 <- extract_test_model_cindex(model_tbl, "M1")
  cindex_test_M2 <- extract_test_model_cindex(model_tbl, "M2")
  cindex_test_M3 <- extract_test_model_cindex(model_tbl, "M3")

  n_features_test_M0 <- extract_test_model_features(model_tbl, "M0")
  n_features_test_M1 <- extract_test_model_features(model_tbl, "M1")
  n_features_test_M2 <- extract_test_model_features(model_tbl, "M2")
  n_features_test_M3 <- extract_test_model_features(model_tbl, "M3")

  label_M2 <- extract_test_model_label(model_tbl, "M2")
  label_M3 <- extract_test_model_label(model_tbl, "M3")

  n_original_features <- coalesce_numeric(
    safe_numeric_scalar(js$n_original_organ_features),
    safe_numeric_scalar(js$n_original_brain_features),
    safe_numeric_scalar(js$n_original_protein_features),
    safe_numeric_scalar(js$n_original_metabolite_features),
    safe_numeric_scalar(js$n_original_adipose_features),
    safe_numeric_scalar(js$n_original_heart_features),
    safe_numeric_scalar(js$n_original_kidney_features),
    safe_numeric_scalar(js$n_original_liver_features),
    safe_numeric_scalar(js$n_original_pancreas_features),
    safe_numeric_scalar(js$n_original_spleen_features),
    safe_numeric_scalar(js$n_original_features)
  )

  meta %>%
    transmute(
      disease_key,
      disease_label,
      disease_clock_label,
      clock_id,
      modality,
      organ_label,
      organ_modality_label,
      clock_label,
      folder,
      prefix,
      json_file_relative,
      json_exists,
      parse_status = "ok",

      n_total = safe_numeric_scalar(js$n_total),
      n_events_total = safe_numeric_scalar(js$n_events_total),
      n_censored_total = safe_numeric_scalar(js$n_censored_total),
      median_followup_years = safe_numeric_scalar(js$median_followup_years),

      n_train = safe_numeric_scalar(js$n_train),
      n_events_train = safe_numeric_scalar(js$n_events_train),
      n_validation = safe_numeric_scalar(js$n_validation),
      n_events_validation = safe_numeric_scalar(js$n_events_validation),
      n_test = safe_numeric_scalar(js$n_test),
      n_events_test = safe_numeric_scalar(js$n_events_test),

      cindex_train = safe_numeric_scalar(js$cindex_train),
      cindex_validation = safe_numeric_scalar(js$cindex_validation),
      cindex_trainval = safe_numeric_scalar(js$cindex_trainval),
      cindex_test = safe_numeric_scalar(js$cindex_test),

      cindex_test_M0_age_sex = cindex_test_M0,
      cindex_test_M1_covariate_baseline = cindex_test_M1,
      cindex_test_M2_features_only = cindex_test_M2,
      cindex_test_M3_full_model = cindex_test_M3,

      delta_cindex_test_M3_vs_M1 = delta_summary$delta_cindex,
      delta_cindex_test_M3_vs_M1_ci_lower = delta_summary$delta_cindex_ci_lower,
      delta_cindex_test_M3_vs_M1_ci_upper = delta_summary$delta_cindex_ci_upper,
      delta_cindex_test_M3_vs_M1_p_two_sided = delta_summary$delta_p_two_sided,
      delta_cindex_test_M3_vs_M1_p_two_sided_display = format_p_display(delta_summary$delta_p_two_sided),
      delta_cindex_test_M3_vs_M1_p_one_sided_le_0 = delta_summary$delta_p_one_sided_le_0,
      delta_n_bootstrap_requested = delta_summary$delta_n_bootstrap_requested,
      delta_n_bootstrap_successful = delta_summary$delta_n_bootstrap_successful,
      delta_significant = case_when(
        is.finite(delta_summary$delta_cindex_ci_lower) &
          delta_summary$delta_cindex_ci_lower > 0 ~ TRUE,
        is.finite(delta_summary$delta_p_two_sided) &
          delta_summary$delta_p_two_sided < 0.05 &
          delta_summary$delta_cindex > 0 ~ TRUE,
        TRUE ~ FALSE
      ),

      train_minus_test_cindex = cindex_train - cindex_test,
      validation_minus_test_cindex = cindex_validation - cindex_test,

      m0_n_features = n_features_test_M0,
      m1_n_features = n_features_test_M1,
      m2_n_features = n_features_test_M2,
      m3_n_features = n_features_test_M3,
      m2_model_label = label_M2,
      m3_model_label = label_M3,

      best_l1_ratio = safe_numeric_scalar(js$best_l1_ratio),
      best_alpha = safe_numeric_scalar(js$best_alpha),
      best_validation_cindex_during_tuning = safe_numeric_scalar(js$best_validation_cindex_during_tuning),
      used_penalty_factor = safe_logical_scalar(js$used_penalty_factor),

      n_original_features = n_original_features,
      n_numeric_cols_kept = safe_numeric_scalar(js$n_numeric_cols_kept),
      n_categorical_cols_kept = safe_numeric_scalar(js$n_categorical_cols_kept),
      n_nonzero_coefficients = safe_numeric_scalar(js$n_nonzero_coefficients),
      n_residualization_covariates = safe_numeric_scalar(js$n_residualization_covariates),

      admin_censor_date = safe_character_scalar(js$admin_censor_date),
      time_zero = safe_character_scalar(js$time_zero),
      event_date = safe_character_scalar(js$event_date),
      note = safe_character_scalar(js$note)
    )
}

# ============================================================
# 7. Parse all hard-coded disease-clock JSONs
# ============================================================

summary_tbl <- purrr::map_dfr(
  seq_len(nrow(clock_manifest)),
  ~ parse_one_summary(clock_manifest[.x, ])
) %>%
  arrange(
    factor(disease_key, levels = disease_tbl$disease_key),
    factor(modality, levels = c("MRI", "Proteomics", "Metabolomics")),
    organ_label
  )

message("Rows in summary table: ", nrow(summary_tbl))
message("Rows with parse_status == ok: ", sum(summary_tbl$parse_status == "ok", na.rm = TRUE))

# ============================================================
# 8. Save Excel workbook only
# ============================================================

wb <- openxlsx::createWorkbook()

header_style <- openxlsx::createStyle(
  fontColour = "white",
  fgFill = "#1F4E78",
  halign = "center",
  valign = "center",
  textDecoration = "bold",
  border = "Bottom"
)

body_style <- openxlsx::createStyle(
  valign = "top",
  wrapText = TRUE
)

num_style <- openxlsx::createStyle(
  numFmt = "0.000"
)

p_style <- openxlsx::createStyle(
  numFmt = "0.00E+00"
)

sheet_name <- "Summary_5_disease_clocks"

openxlsx::addWorksheet(wb, sheet_name)
openxlsx::writeData(wb, sheet_name, summary_tbl, withFilter = TRUE)

openxlsx::addStyle(
  wb,
  sheet = sheet_name,
  style = header_style,
  rows = 1,
  cols = 1:ncol(summary_tbl),
  gridExpand = TRUE
)

openxlsx::addStyle(
  wb,
  sheet = sheet_name,
  style = body_style,
  rows = 2:(nrow(summary_tbl) + 1),
  cols = 1:ncol(summary_tbl),
  gridExpand = TRUE,
  stack = TRUE
)

openxlsx::freezePane(
  wb,
  sheet = sheet_name,
  firstActiveRow = 2,
  firstActiveCol = 6
)

openxlsx::setColWidths(
  wb,
  sheet = sheet_name,
  cols = 1:ncol(summary_tbl),
  widths = "auto"
)

summary_cols <- colnames(summary_tbl)

numeric_cols <- which(
  stringr::str_detect(
    summary_cols,
    "cindex|followup|alpha|ratio|delta|n_|events|features|coefficients"
  )
)

if (length(numeric_cols) > 0) {
  openxlsx::addStyle(
    wb,
    sheet = sheet_name,
    style = num_style,
    rows = 2:(nrow(summary_tbl) + 1),
    cols = numeric_cols,
    gridExpand = TRUE,
    stack = TRUE
  )
}

p_cols <- which(
  stringr::str_detect(summary_cols, "p_|p_two_sided|p_one_sided")
)

if (length(p_cols) > 0) {
  openxlsx::addStyle(
    wb,
    sheet = sheet_name,
    style = p_style,
    rows = 2:(nrow(summary_tbl) + 1),
    cols = p_cols,
    gridExpand = TRUE,
    stack = TRUE
  )
}

openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)

message("Saved Excel workbook:")
message("  ", out_xlsx)

# ============================================================
# 9. Console summary
# ============================================================

message("\n===== Disease-clock summary =====")

print(
  summary_tbl %>%
    group_by(disease_label, parse_status) %>%
    summarise(n = n(), .groups = "drop")
)

message("\n===== Test-set performance overview =====")

print(
  summary_tbl %>%
    filter(parse_status == "ok") %>%
    select(
      disease_label,
      modality,
      organ_label,
      n_total,
      n_events_total,
      n_test,
      n_events_test,
      cindex_test,
      cindex_test_M1_covariate_baseline,
      cindex_test_M3_full_model,
      delta_cindex_test_M3_vs_M1,
      delta_cindex_test_M3_vs_M1_ci_lower,
      delta_cindex_test_M3_vs_M1_ci_upper,
      delta_cindex_test_M3_vs_M1_p_two_sided_display,
      delta_significant
    ) %>%
    arrange(
      factor(disease_label, levels = disease_tbl$disease_label),
      factor(modality, levels = c("MRI", "Proteomics", "Metabolomics")),
      desc(cindex_test)
    )
)

message("\nDone.")