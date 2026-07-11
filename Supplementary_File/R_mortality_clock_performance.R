#!/usr/bin/env Rscript

# ============================================================
# Collect model performance statistics for 23 mortality EPOCH clocks
# Hard-coded clock folders to avoid scanning the full base directory
#
# Input:
#   <base_dir>/<hard-coded mortality clock folder>/<prefix>_performance.json
#
# Output:
#   /Users/hao/Dropbox/2026_EPOCH/Supplementary_File/
#     Supplementary_Table_Mortality_EPOCH_Model_Performance_23_clocks.xlsx
#
# Excel sheets:
#   1) README
#   2) Summary_23_clocks
#   3) Model_comparison
#   4) Delta_Cindex
#   5) Top_level_scalars
#   6) Residualization_covariates
#   7) JSON_manifest
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
  "Supplementary_Table_Mortality_EPOCH_Model_Performance_23_clocks.xlsx"
)

out_summary_tsv <- file.path(
  out_dir,
  "Supplementary_Table_Mortality_EPOCH_Model_Performance_23_clocks_summary.tsv"
)

out_model_tsv <- file.path(
  out_dir,
  "Supplementary_Table_Mortality_EPOCH_Model_Performance_23_clocks_model_comparison.tsv"
)

out_delta_tsv <- file.path(
  out_dir,
  "Supplementary_Table_Mortality_EPOCH_Model_Performance_23_clocks_delta_cindex.tsv"
)

message("Base directory:")
message("  ", base_dir)
message("Output Excel:")
message("  ", out_xlsx)

# ============================================================
# 2. Hard-coded 23 mortality clock manifest
# ============================================================

clock_manifest <- tribble(
  ~modality,       ~modality_key,    ~organ_key,              ~organ_label,              ~folder,                                                ~prefix,
  "MRI",           "mri",            "adipose",               "Adipose",                 "adipose_mri_mortality_clock",                         "adipose_mri_mortality_clock",
  "MRI",           "mri",            "brain",                 "Brain",                   "brain_mri_mortality_clock",                           "brain_mri_mortality_clock",
  "MRI",           "mri",            "heart",                 "Heart",                   "heart_mri_mortality_clock",                           "heart_mri_mortality_clock",
  "MRI",           "mri",            "kidney",                "Kidney",                  "kidney_mri_mortality_clock",                          "kidney_mri_mortality_clock",
  "MRI",           "mri",            "liver",                 "Liver",                   "liver_mri_mortality_clock",                           "liver_mri_mortality_clock",
  "MRI",           "mri",            "pancreas",              "Pancreas",                "pancreas_mri_mortality_clock",                        "pancreas_mri_mortality_clock",
  "MRI",           "mri",            "spleen",                "Spleen",                  "spleen_mri_mortality_clock",                          "spleen_mri_mortality_clock",
  
  "Proteomics",    "proteomics",     "brain",                 "Brain",                   "Brain_proteomics_mortality_clock",                    "brain_proteomics_mortality_clock",
  "Proteomics",    "proteomics",     "endocrine",             "Endocrine",               "Endocrine_proteomics_mortality_clock",                "endocrine_proteomics_mortality_clock",
  "Proteomics",    "proteomics",     "eye",                   "Eye",                     "Eye_proteomics_mortality_clock",                      "eye_proteomics_mortality_clock",
  "Proteomics",    "proteomics",     "heart",                 "Heart",                   "Heart_proteomics_mortality_clock",                    "heart_proteomics_mortality_clock",
  "Proteomics",    "proteomics",     "hepatic",               "Hepatic",                 "Hepatic_proteomics_mortality_clock",                  "hepatic_proteomics_mortality_clock",
  "Proteomics",    "proteomics",     "immune",                "Immune",                  "Immune_proteomics_mortality_clock",                   "immune_proteomics_mortality_clock",
  "Proteomics",    "proteomics",     "pulmonary",             "Pulmonary",               "Pulmonary_proteomics_mortality_clock",                "pulmonary_proteomics_mortality_clock",
  "Proteomics",    "proteomics",     "renal",                 "Renal",                   "Renal_proteomics_mortality_clock",                    "renal_proteomics_mortality_clock",
  "Proteomics",    "proteomics",     "reproductive_female",   "Reproductive female",     "Reproductive_female_proteomics_mortality_clock",      "reproductive_female_proteomics_mortality_clock",
  "Proteomics",    "proteomics",     "reproductive_male",     "Reproductive male",       "Reproductive_male_proteomics_mortality_clock",        "reproductive_male_proteomics_mortality_clock",
  "Proteomics",    "proteomics",     "skin",                  "Skin",                    "Skin_proteomics_mortality_clock",                     "skin_proteomics_mortality_clock",
  
  "Metabolomics",  "metabolomics",   "digestive",             "Digestive",               "Digestive_metabolomics_mortality_clock",              "digestive_metabolomics_mortality_clock",
  "Metabolomics",  "metabolomics",   "endocrine",             "Endocrine",               "Endocrine_metabolomics_mortality_clock",              "endocrine_metabolomics_mortality_clock",
  "Metabolomics",  "metabolomics",   "hepatic",               "Hepatic",                 "Hepatic_metabolomics_mortality_clock",                "hepatic_metabolomics_mortality_clock",
  "Metabolomics",  "metabolomics",   "immune",                "Immune",                  "Immune_metabolomics_mortality_clock",                 "immune_metabolomics_mortality_clock",
  "Metabolomics",  "metabolomics",   "metabolic",             "Metabolic",               "Metabolic_metabolomics_mortality_clock",              "metabolic_metabolomics_mortality_clock"
) %>%
  mutate(
    clock_id = paste(organ_key, modality_key, sep = "__"),
    clock_label = paste(organ_label, modality),
    json_file = file.path(base_dir, folder, paste0(prefix, "_performance.json")),
    json_file_relative = file.path(folder, paste0(prefix, "_performance.json")),
    exists = file.exists(json_file)
  ) %>%
  arrange(
    factor(modality, levels = c("MRI", "Proteomics", "Metabolomics")),
    organ_label
  )

message("Hard-coded clocks: ", nrow(clock_manifest))
message("Detected existing JSON files: ", sum(clock_manifest$exists), " / ", nrow(clock_manifest))

missing_files <- clock_manifest %>% filter(!exists)

if (nrow(missing_files) > 0) {
  warning(
    "Some hard-coded JSON files are missing:\n",
    paste0("  ", missing_files$json_file_relative, collapse = "\n")
  )
}

# ============================================================
# 3. Helper functions
# ============================================================

is_scalar_value <- function(x) {
  is.atomic(x) && length(x) == 1
}

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

get_scalar_tbl <- function(js) {
  scalar_names <- names(js)[purrr::map_lgl(js, is_scalar_value)]
  
  if (length(scalar_names) == 0) {
    return(tibble())
  }
  
  as_tibble(as.list(js[scalar_names]))
}

extract_top_level_numeric <- function(scalar_tbl, pattern) {
  hit <- colnames(scalar_tbl)[stringr::str_detect(colnames(scalar_tbl), pattern)]
  if (length(hit) == 0) return(NA_real_)
  suppressWarnings(as.numeric(scalar_tbl[[hit[[1]]]][[1]]))
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
# 4. Parse each of the 23 hard-coded JSON files
# ============================================================

parsed_list <- purrr::map(seq_len(nrow(clock_manifest)), function(i) {
  meta <- clock_manifest[i, ]
  
  if (!isTRUE(meta$exists)) {
    warning("Skipping missing JSON: ", meta$json_file_relative)
    return(NULL)
  }
  
  message("Reading: ", meta$json_file_relative)
  
  js <- read_json_safe(meta$json_file)
  
  if (is.null(js)) {
    return(NULL)
  }
  
  scalar_tbl <- get_scalar_tbl(js) %>%
    mutate(clock_id = meta$clock_id, .before = 1)
  
  model_tbl <- list_rows_to_tibble(js$incremental_value_model_comparison) %>%
    mutate(
      clock_id = meta$clock_id,
      modality = meta$modality,
      organ_label = meta$organ_label,
      clock_label = meta$clock_label,
      folder = meta$folder,
      prefix = meta$prefix,
      json_file = meta$json_file,
      .before = 1
    )
  
  delta_tbl <- list_rows_to_tibble(js$incremental_value_delta_cindex) %>%
    mutate(
      clock_id = meta$clock_id,
      modality = meta$modality,
      organ_label = meta$organ_label,
      clock_label = meta$clock_label,
      folder = meta$folder,
      prefix = meta$prefix,
      json_file = meta$json_file,
      .before = 1
    )
  
  residual_cov_tbl <- tibble()
  
  if (!is.null(js$residualization_covariates)) {
    residual_cov_tbl <- tibble(
      clock_id = meta$clock_id,
      modality = meta$modality,
      organ_label = meta$organ_label,
      clock_label = meta$clock_label,
      folder = meta$folder,
      residualization_covariate = as.character(js$residualization_covariates)
    )
  }
  
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
  
  # Fallback: top-level direct fields.
  if (!is.finite(cindex_test_M1)) {
    cindex_test_M1 <- extract_top_level_numeric(
      scalar_tbl,
      "^cindex_test_M1_covariate_baseline$"
    )
  }
  
  if (!is.finite(cindex_test_M3)) {
    cindex_test_M3 <- extract_top_level_numeric(
      scalar_tbl,
      "^cindex_test_M3_full_covariates_plus_"
    )
  }
  
  n_original_features <- coalesce_numeric(
    safe_numeric_scalar(js$n_original_brain_features),
    safe_numeric_scalar(js$n_original_protein_features),
    safe_numeric_scalar(js$n_original_metabolite_features),
    safe_numeric_scalar(js$n_original_features)
  )
  
  summary_tbl_i <- meta %>%
    transmute(
      clock_id,
      modality,
      organ_label,
      clock_label,
      folder,
      prefix,
      json_file,
      json_file_relative,
      
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
      delta_interpretation = delta_summary$delta_interpretation,
      
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
  
  list(
    summary = summary_tbl_i,
    model = model_tbl,
    delta = delta_tbl,
    scalar = scalar_tbl,
    residual_covariates = residual_cov_tbl
  )
})

parsed_list <- purrr::compact(parsed_list)

if (length(parsed_list) == 0) {
  stop("No hard-coded JSON files could be parsed.")
}

summary_tbl <- purrr::map_dfr(parsed_list, "summary") %>%
  arrange(
    factor(modality, levels = c("MRI", "Proteomics", "Metabolomics")),
    organ_label
  )

model_comparison_tbl <- purrr::map_dfr(parsed_list, "model") %>%
  mutate(
    split = factor(
      split,
      levels = c("train", "validation", "test", "trainval")
    )
  ) %>%
  arrange(
    factor(modality, levels = c("MRI", "Proteomics", "Metabolomics")),
    organ_label,
    model,
    split
  )

delta_tbl <- purrr::map_dfr(parsed_list, "delta") %>%
  arrange(
    factor(modality, levels = c("MRI", "Proteomics", "Metabolomics")),
    organ_label
  ) %>%
  mutate(
    p_two_sided_display = if ("empirical_p_two_sided_delta_not_equal_0" %in% colnames(.)) {
      format_p_display(as.numeric(empirical_p_two_sided_delta_not_equal_0))
    } else {
      NA_character_
    }
  )

scalar_tbl <- purrr::map_dfr(parsed_list, "scalar") %>%
  left_join(
    clock_manifest %>%
      select(clock_id, modality, organ_label, clock_label, folder, prefix, json_file, json_file_relative),
    by = "clock_id"
  ) %>%
  relocate(modality, organ_label, clock_label, folder, prefix, json_file, json_file_relative, .after = clock_id)

residual_cov_tbl <- purrr::map_dfr(parsed_list, "residual_covariates") %>%
  arrange(
    factor(modality, levels = c("MRI", "Proteomics", "Metabolomics")),
    organ_label,
    residualization_covariate
  )

manifest_tbl <- clock_manifest %>%
  select(
    modality,
    organ_label,
    clock_label,
    clock_id,
    folder,
    prefix,
    json_file_relative,
    json_file,
    exists
  )

# ============================================================
# 5. Write TSV backups
# ============================================================

readr::write_tsv(summary_tbl, out_summary_tsv)
readr::write_tsv(model_comparison_tbl, out_model_tsv)
readr::write_tsv(delta_tbl, out_delta_tsv)

message("Saved TSV backups:")
message("  ", out_summary_tsv)
message("  ", out_model_tsv)
message("  ", out_delta_tsv)

# ============================================================
# 6. Write Excel workbook
# ============================================================

readme_tbl <- tibble(
  field = c(
    "Purpose",
    "Input files",
    "Hard-coded clocks",
    "Detected JSON files",
    "Summary_23_clocks",
    "Model_comparison",
    "Delta_Cindex",
    "Top_level_scalars",
    "Residualization_covariates",
    "Primary discrimination metric",
    "Primary incremental-value metric",
    "Output generated on"
  ),
  description = c(
    "Collect model performance statistics for all 23 mortality EPOCH clocks.",
    "Only the hard-coded 23 performance JSON files are read; the script does not scan the full base directory.",
    as.character(nrow(clock_manifest)),
    as.character(sum(clock_manifest$exists)),
    "One row per mortality clock with sample sizes, event counts, C-index values, incremental ΔC-index, and tuning/feature summary fields.",
    "Long-format model-comparison table extracted from incremental_value_model_comparison in each JSON.",
    "Long-format incremental ΔC-index table extracted from incremental_value_delta_cindex in each JSON.",
    "All scalar top-level JSON fields retained for reproducibility.",
    "One row per residualization covariate per clock.",
    "Held-out test C-index for the mortality EPOCH model.",
    "Test-set ΔC-index comparing M3 full model versus M1 covariate baseline.",
    as.character(Sys.time())
  )
)

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

write_sheet <- function(wb, sheet_name, df, freeze_col = 1) {
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(wb, sheet_name, df, withFilter = TRUE)
  
  if (nrow(df) >= 1 && ncol(df) >= 1) {
    openxlsx::addStyle(
      wb,
      sheet = sheet_name,
      style = header_style,
      rows = 1,
      cols = 1:ncol(df),
      gridExpand = TRUE
    )
    
    if (nrow(df) >= 1) {
      openxlsx::addStyle(
        wb,
        sheet = sheet_name,
        style = body_style,
        rows = 2:(nrow(df) + 1),
        cols = 1:ncol(df),
        gridExpand = TRUE,
        stack = TRUE
      )
    }
    
    openxlsx::freezePane(
      wb,
      sheet = sheet_name,
      firstActiveRow = 2,
      firstActiveCol = freeze_col + 1
    )
    
    openxlsx::setColWidths(
      wb,
      sheet = sheet_name,
      cols = 1:ncol(df),
      widths = "auto"
    )
  }
}

write_sheet(wb, "README", readme_tbl, freeze_col = 1)
write_sheet(wb, "Summary_23_clocks", summary_tbl, freeze_col = 4)
write_sheet(wb, "Model_comparison", model_comparison_tbl, freeze_col = 5)
write_sheet(wb, "Delta_Cindex", delta_tbl, freeze_col = 5)
write_sheet(wb, "Top_level_scalars", scalar_tbl, freeze_col = 5)
write_sheet(wb, "Residualization_covariates", residual_cov_tbl, freeze_col = 4)
write_sheet(wb, "JSON_manifest", manifest_tbl, freeze_col = 4)

# Apply numeric formatting to selected summary columns.
summary_cols <- colnames(summary_tbl)

numeric_cols_summary <- which(
  stringr::str_detect(
    summary_cols,
    "cindex|followup|alpha|ratio|delta|n_|events|features|coefficients"
  )
)

if (length(numeric_cols_summary) > 0) {
  openxlsx::addStyle(
    wb,
    sheet = "Summary_23_clocks",
    style = num_style,
    rows = 2:(nrow(summary_tbl) + 1),
    cols = numeric_cols_summary,
    gridExpand = TRUE,
    stack = TRUE
  )
}

p_cols_summary <- which(
  stringr::str_detect(summary_cols, "p_|p_two_sided|p_one_sided")
)

if (length(p_cols_summary) > 0) {
  openxlsx::addStyle(
    wb,
    sheet = "Summary_23_clocks",
    style = p_style,
    rows = 2:(nrow(summary_tbl) + 1),
    cols = p_cols_summary,
    gridExpand = TRUE,
    stack = TRUE
  )
}

openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)

message("Saved Excel workbook:")
message("  ", out_xlsx)

# ============================================================
# 7. Console summary
# ============================================================

message("\n===== Summary of collected clocks =====")

print(
  summary_tbl %>%
    select(
      modality,
      organ_label,
      n_total,
      n_events_total,
      median_followup_years,
      cindex_train,
      cindex_validation,
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
      factor(modality, levels = c("MRI", "Proteomics", "Metabolomics")),
      desc(cindex_test)
    )
)

message("\nHard-coded clocks: ", nrow(clock_manifest))
message("Existing JSON files: ", sum(clock_manifest$exists), " / ", nrow(clock_manifest))
message("Rows in Summary_23_clocks: ", nrow(summary_tbl))
message("Rows in Model_comparison: ", nrow(model_comparison_tbl))
message("Rows in Delta_Cindex: ", nrow(delta_tbl))
message("\nDone.")