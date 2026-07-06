# ============================================================
# Disease L'EPOCH incremental value + year-scale QC visualization
#
# Revised version:
#   - Does NOT read *_predictions.tsv
#   - Gets N / cases from:
#       all_<disease>_lepoch_summary.tsv
#       or small *_performance.json files
#   - Gets M3-M1 delta and CI from:
#       all_<disease>_lepoch_summary.tsv
#       all_<disease>_lepoch_incremental_value_summary.tsv
#       or per-clock *_incremental_value_delta_cindex.tsv
#   - Gets scale issues from:
#       <disease>_lepoch_year_scale_qc/<disease>_year_scale_qc_summary.tsv
#
# Forest plot:
#   - filled circle = significant + stable scale
#   - open triangle = significant but scale issue
#   - open circle = non-significant + stable scale
#   - cross = non-significant + scale issue
#   - x-axis = -abs(max CI) to +abs(max CI)
#   - labels to right show N, cases, and CI
#
# Diseases:
#   asthma, dementia, copd, mi, stroke
#
# Outputs:
#   <base_dir>/all_disease_lepoch_incremental_value_scale_qc/
# ============================================================

.libPaths('/gpfs/fs001/cbica/home/wenju/R/x86_64-pc-linux-gnu-library/4.3')

suppressPackageStartupMessages({
  library(tidyverse)
  library(jsonlite)
  library(glue)
  library(scales)
  library(patchwork)
})

# ============================================================
# 1. Settings
# ============================================================

possible_base_dirs <- c(
  "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock",
  "/gpfs/fs001/cbica/home/wenju/Reproducibile_paper/WholeBodyClock",
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock",
  getwd()
)

base_dir <- possible_base_dirs[file.exists(possible_base_dirs)][1]
if (is.na(base_dir) || is.null(base_dir)) {
  stop("Could not detect base_dir. Please set base_dir manually.")
}

diseases <- c("asthma", "dementia", "copd", "mi", "stroke")

disease_labels <- c(
  asthma = "Asthma",
  dementia = "Dementia",
  copd = "COPD",
  mi = "MI",
  stroke = "Stroke"
)

alpha_sig <- 0.05

outdir <- file.path(
  base_dir,
  "all_disease_lepoch_incremental_value_scale_qc"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

message("Base directory: ", base_dir)
message("Output directory: ", outdir)

# ============================================================
# 2. Theme and palettes
# ============================================================

theme_clock <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(color = "#17202A"),
      plot.title = element_text(
        face = "bold",
        size = base_size + 4,
        margin = margin(b = 6)
      ),
      plot.subtitle = element_text(
        size = base_size,
        color = "#566573",
        margin = margin(b = 8)
      ),
      plot.caption = element_text(
        size = base_size - 3,
        color = "#7B7D7D"
      ),
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(color = "#2C3E50"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "#E5E8E8", linewidth = 0.35),
      panel.grid.major.y = element_line(color = "#EAECEE", linewidth = 0.35),
      strip.text = element_text(face = "bold", size = base_size + 1),
      strip.background = element_rect(fill = "#F4F6F7", color = "#D5D8DC"),
      legend.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.box = "vertical",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

# Van-Gogh-inspired, high-contrast, colorblind-conscious palette.
vangogh23 <- c(
  "#1F4E79", "#2E86AB", "#3A86FF", "#00A5CF", "#148F77",
  "#2A9D8F", "#4DAA57", "#84A59D", "#6A4C93", "#8E44AD",
  "#7B2CBF", "#4361EE", "#F2CC8F", "#E9C46A", "#F4A261",
  "#DDA15E", "#B9770E", "#FB8500", "#E76F51", "#C1121F",
  "#7D6608", "#566573", "#17202A"
)

# Shape 2 = open / non-filled triangle.
shape_values <- c(
  "Good: significant + stable scale" = 16,
  "Significant but scale issue" = 2,
  "Non-significant + stable scale" = 1,
  "Non-significant + scale issue" = 4,
  "Missing M3-M1 or QC" = 3
)

status_cols <- c(
  "Good: significant + stable scale" = "#1B9E77",
  "Significant but scale issue" = "#D95F02",
  "Non-significant + stable scale" = "#7570B3",
  "Non-significant + scale issue" = "#E7298A",
  "Missing M3-M1 or QC" = "#7F8C8D"
)

# ============================================================
# 3. Helper functions
# ============================================================

safe_read_tsv <- function(path) {
  if (is.na(path) || is.null(path) || !file.exists(path)) {
    return(tibble())
  }
  
  readr::read_tsv(
    path,
    show_col_types = FALSE,
    progress = FALSE
  )
}

safe_read_json <- function(path) {
  if (is.na(path) || is.null(path) || !file.exists(path)) {
    return(list())
  }
  
  tryCatch(
    jsonlite::fromJSON(path),
    error = function(e) list()
  )
}

json_get_num <- function(x, field) {
  if (is.null(x[[field]])) return(NA_real_)
  out <- suppressWarnings(as.numeric(x[[field]]))
  if (length(out) == 0) return(NA_real_)
  out[[1]]
}

first_existing_col <- function(df, candidates) {
  hits <- candidates[candidates %in% colnames(df)]
  if (length(hits) == 0) return(NA_character_)
  hits[[1]]
}

pick_num_col <- function(df, candidates) {
  col <- first_existing_col(df, candidates)
  if (is.na(col)) return(rep(NA_real_, nrow(df)))
  suppressWarnings(as.numeric(df[[col]]))
}

pick_chr_col <- function(df, candidates) {
  col <- first_existing_col(df, candidates)
  if (is.na(col)) return(rep(NA_character_, nrow(df)))
  as.character(df[[col]])
}

safe_first_num <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  x[[1]]
}

infer_sig_positive <- function(delta, lo, p, existing_sig, alpha = 0.05) {
  existing_sig <- as.character(existing_sig)
  
  case_when(
    is.finite(delta) & is.finite(lo) ~ delta > 0 & lo > 0,
    is.finite(delta) & is.finite(p) ~ delta > 0 & p < alpha,
    is.finite(delta) &
      !is.na(existing_sig) &
      stringr::str_to_lower(existing_sig) %in% c("significant", "true", "yes") ~ delta > 0,
    TRUE ~ FALSE
  )
}

format_ci_label <- function(lo, hi) {
  case_when(
    is.finite(lo) & is.finite(hi) ~ sprintf("CI %.3f, %.3f", lo, hi),
    TRUE ~ "CI NA"
  )
}

format_n_cases_label <- function(n_total, n_events) {
  case_when(
    is.finite(n_total) & is.finite(n_events) ~ paste0(
      "N=", scales::comma(round(n_total)),
      " (cases=", scales::comma(round(n_events)), ")"
    ),
    is.finite(n_events) ~ paste0("cases=", scales::comma(round(n_events))),
    TRUE ~ "N/cases=NA"
  )
}

# ============================================================
# 4. Manifest for 23 disease clocks
# ============================================================

make_manifest_one_disease <- function(disease) {
  mri_organs <- c(
    "brain",
    "heart",
    "adipose",
    "kidney",
    "liver",
    "pancreas",
    "spleen"
  )
  
  proteomics_organs <- c(
    "Pulmonary",
    "Heart",
    "Brain",
    "Eye",
    "Hepatic",
    "Renal",
    "Endocrine",
    "Immune",
    "Skin",
    "Reproductive_female",
    "Reproductive_male"
  )
  
  metabolomics_organs <- c(
    "Endocrine",
    "Digestive",
    "Hepatic",
    "Immune",
    "Metabolic"
  )
  
  bind_rows(
    tibble(
      disease = disease,
      modality = "MRI",
      modality_key = "mri",
      organ_folder_name = mri_organs,
      organ_key = stringr::str_to_lower(mri_organs),
      folder = paste0(stringr::str_to_lower(mri_organs), "_mri_", disease, "_clock"),
      prefix = paste0(stringr::str_to_lower(mri_organs), "_mri_", disease, "_clock")
    ),
    tibble(
      disease = disease,
      modality = "Proteomics",
      modality_key = "proteomics",
      organ_folder_name = proteomics_organs,
      organ_key = stringr::str_to_lower(proteomics_organs),
      folder = paste0(proteomics_organs, "_proteomics_", disease, "_clock"),
      prefix = paste0(stringr::str_to_lower(proteomics_organs), "_proteomics_", disease, "_clock")
    ),
    tibble(
      disease = disease,
      modality = "Metabolomics",
      modality_key = "metabolomics",
      organ_folder_name = metabolomics_organs,
      organ_key = stringr::str_to_lower(metabolomics_organs),
      folder = paste0(metabolomics_organs, "_metabolomics_", disease, "_clock"),
      prefix = paste0(stringr::str_to_lower(metabolomics_organs), "_metabolomics_", disease, "_clock")
    )
  ) %>%
    mutate(
      organ_label = organ_folder_name %>%
        stringr::str_replace_all("_", " ") %>%
        stringr::str_to_sentence(),
      clock_label = paste(organ_label, modality),
      clock_id = paste(organ_key, modality_key, sep = "__"),
      disease_label = unname(disease_labels[disease]),
      clock_dir = file.path(base_dir, folder),
      
      performance_file = file.path(clock_dir, paste0(prefix, "_performance.json")),
      model_comparison_file = file.path(clock_dir, paste0(prefix, "_model_comparison.tsv")),
      delta_file = file.path(clock_dir, paste0(prefix, "_incremental_value_delta_cindex.tsv"))
    )
}

clock_manifest <- purrr::map_dfr(diseases, make_manifest_one_disease)

clock_order_ref <- make_manifest_one_disease(diseases[[1]]) %>%
  distinct(clock_id, clock_label, modality) %>%
  mutate(
    modality = factor(modality, levels = c("MRI", "Proteomics", "Metabolomics"))
  ) %>%
  arrange(modality, clock_label)

clock_levels <- clock_order_ref$clock_label

clock_cols <- setNames(
  vangogh23[seq_along(clock_levels)],
  clock_levels
)

manifest_check_out <- file.path(outdir, "all_disease_lepoch_manifest_file_check.tsv")

manifest_check <- clock_manifest %>%
  mutate(
    clock_dir_exists = file.exists(clock_dir),
    performance_file_exists = file.exists(performance_file),
    model_comparison_file_exists = file.exists(model_comparison_file),
    delta_file_exists = file.exists(delta_file)
  )

readr::write_tsv(manifest_check, manifest_check_out)

message("Manifest file check written to:")
message("  ", manifest_check_out)

# ============================================================
# 5. Read year-scale QC summaries
# ============================================================

read_qc_for_disease <- function(disease) {
  qc_dir <- file.path(base_dir, paste0(disease, "_lepoch_year_scale_qc"))
  
  qc_file <- file.path(
    qc_dir,
    paste0(disease, "_year_scale_qc_summary.tsv")
  )
  
  qc_problematic_file <- file.path(
    qc_dir,
    paste0(disease, "_problematic_year_scale_clocks.tsv")
  )
  
  qc <- safe_read_tsv(qc_file)
  
  if (nrow(qc) == 0) {
    warning("Missing or empty QC summary for ", disease, ": ", qc_file)
    
    return(
      make_manifest_one_disease(disease) %>%
        transmute(
          disease,
          folder,
          qc_file = qc_file,
          qc_problematic_file = qc_problematic_file,
          scale_qc_status = NA_character_,
          scale_qc_reason = NA_character_,
          age_beta = NA_real_,
          scale_issue = NA
        )
    )
  }
  
  qc %>%
    transmute(
      disease = disease,
      folder = as.character(clock_folder),
      qc_file = qc_file,
      qc_problematic_file = qc_problematic_file,
      scale_qc_status = as.character(final_year_scale_qc_status),
      scale_qc_reason = as.character(final_year_scale_qc_reason),
      age_beta = suppressWarnings(as.numeric(age_beta)),
      scale_issue = stringr::str_detect(scale_qc_status, "^(WARN|FAIL)")
    )
}

qc_all <- purrr::map_dfr(diseases, read_qc_for_disease)

qc_out <- file.path(outdir, "all_disease_lepoch_year_scale_qc_merged.tsv")
readr::write_tsv(qc_all, qc_out)

message("Merged QC table:")
message("  ", qc_out)

# ============================================================
# 6. Read sample sizes from small performance JSON files
# ============================================================

read_performance_counts <- function(performance_file) {
  perf <- safe_read_json(performance_file)
  
  tibble(
    n_total_json = json_get_num(perf, "n_total"),
    n_events_json = json_get_num(perf, "n_events_total")
  )
}

performance_counts <- clock_manifest %>%
  mutate(row_id = row_number()) %>%
  split(.$row_id) %>%
  purrr::map_dfr(function(meta_row) {
    meta_row <- meta_row[1, ]
    
    cnt <- read_performance_counts(meta_row$performance_file)
    
    bind_cols(
      meta_row %>%
        select(disease, clock_id, folder, performance_file),
      cnt
    )
  })

performance_counts_out <- file.path(
  outdir,
  "all_disease_lepoch_performance_json_counts.tsv"
)

readr::write_tsv(performance_counts, performance_counts_out)

message("Performance JSON count table:")
message("  ", performance_counts_out)

# ============================================================
# 7. Read incremental-value summaries
# ============================================================

standardize_combined_summary <- function(tbl, disease) {
  if (nrow(tbl) == 0) return(tibble())
  
  tibble(
    disease = disease,
    clock_id = pick_chr_col(tbl, c("clock_id")),
    
    n_total_summary = pick_num_col(
      tbl,
      c("n_total", "N", "n")
    ),
    
    n_events_summary = pick_num_col(
      tbl,
      c("n_events_total", "events", "n_events", "n_cases", "cases")
    ),
    
    cindex_train = pick_num_col(tbl, c("cindex_train")),
    cindex_validation = pick_num_col(tbl, c("cindex_validation")),
    cindex_test = pick_num_col(tbl, c("cindex_test")),
    
    cindex_test_M1_covariate_baseline = pick_num_col(
      tbl,
      c(
        "cindex_test_M1_covariate_baseline",
        "cindex_test_M1",
        "m1_cindex_test"
      )
    ),
    
    cindex_test_M3_full_model = pick_num_col(
      tbl,
      c(
        "cindex_test_M3_full_model",
        "cindex_test_M3",
        "m3_cindex_test"
      )
    ),
    
    delta_cindex = pick_num_col(
      tbl,
      c(
        "delta_cindex_test_M3_vs_M1",
        "delta_cindex",
        "delta_cindex_test_M3_vs_M1_value"
      )
    ),
    
    delta_cindex_ci_lower = pick_num_col(
      tbl,
      c(
        "delta_cindex_test_M3_vs_M1_ci_lower",
        "delta_cindex_ci_lower",
        "ci_lower"
      )
    ),
    
    delta_cindex_ci_upper = pick_num_col(
      tbl,
      c(
        "delta_cindex_test_M3_vs_M1_ci_upper",
        "delta_cindex_ci_upper",
        "ci_upper"
      )
    ),
    
    delta_cindex_p = pick_num_col(
      tbl,
      c(
        "delta_cindex_test_M3_vs_M1_p_two_sided",
        "empirical_p_two_sided_delta_not_equal_0",
        "p",
        "p_value"
      )
    ),
    
    delta_significant_existing = pick_chr_col(
      tbl,
      c(
        "delta_significant",
        "significance",
        "delta_cindex_significant_positive"
      )
    )
  )
}

read_incremental_from_delta_file <- function(meta_row) {
  meta <- as.list(meta_row)
  
  delta_tbl <- safe_read_tsv(meta$delta_file)
  perf <- safe_read_json(meta$performance_file)
  
  delta <- NA_real_
  lo <- NA_real_
  hi <- NA_real_
  p <- NA_real_
  
  if (nrow(delta_tbl) > 0) {
    delta <- safe_first_num(
      pick_num_col(
        delta_tbl,
        c(
          "delta_cindex",
          "delta_cindex_test_M3_vs_M1"
        )
      )
    )
    
    lo <- safe_first_num(
      pick_num_col(
        delta_tbl,
        c(
          "delta_cindex_ci_lower",
          "delta_cindex_test_M3_vs_M1_ci_lower"
        )
      )
    )
    
    hi <- safe_first_num(
      pick_num_col(
        delta_tbl,
        c(
          "delta_cindex_ci_upper",
          "delta_cindex_test_M3_vs_M1_ci_upper"
        )
      )
    )
    
    p <- safe_first_num(
      pick_num_col(
        delta_tbl,
        c(
          "empirical_p_two_sided_delta_not_equal_0",
          "delta_cindex_test_M3_vs_M1_p_two_sided",
          "p",
          "p_value"
        )
      )
    )
  }
  
  cindex_train <- json_get_num(perf, "cindex_train")
  cindex_validation <- json_get_num(perf, "cindex_validation")
  cindex_test <- json_get_num(perf, "cindex_test")
  m1 <- json_get_num(perf, "cindex_test_M1_covariate_baseline")
  
  n_total_json <- json_get_num(perf, "n_total")
  n_events_json <- json_get_num(perf, "n_events_total")
  
  m3_fields <- c(
    paste0(
      "cindex_test_M3_full_covariates_plus_",
      meta$organ_key,
      "_",
      meta$modality_key
    ),
    paste0(
      "cindex_test_M3_full_model_plus_",
      meta$organ_key,
      "_",
      meta$modality_key
    ),
    "cindex_test_M3_full_model",
    "cindex_test_M3"
  )
  
  m3 <- NA_real_
  
  for (field in m3_fields) {
    val <- json_get_num(perf, field)
    if (is.finite(val)) {
      m3 <- val
      break
    }
  }
  
  if ((!is.finite(m1) || !is.finite(m3)) && file.exists(meta$model_comparison_file)) {
    mc <- safe_read_tsv(meta$model_comparison_file)
    
    if (nrow(mc) > 0 && all(c("split", "model", "cindex") %in% colnames(mc))) {
      m1_row <- mc %>%
        filter(split == "test", stringr::str_detect(model, "^M1")) %>%
        slice(1)
      
      m3_row <- mc %>%
        filter(split == "test", stringr::str_detect(model, "^M3")) %>%
        slice(1)
      
      if (nrow(m1_row) > 0) {
        m1 <- suppressWarnings(as.numeric(m1_row$cindex[[1]]))
      }
      
      if (nrow(m3_row) > 0) {
        m3 <- suppressWarnings(as.numeric(m3_row$cindex[[1]]))
      }
    }
  }
  
  if (!is.finite(delta) && is.finite(m1) && is.finite(m3)) {
    delta <- m3 - m1
  }
  
  tibble(
    disease = meta$disease,
    clock_id = meta$clock_id,
    n_total_summary = n_total_json,
    n_events_summary = n_events_json,
    cindex_train = cindex_train,
    cindex_validation = cindex_validation,
    cindex_test = cindex_test,
    cindex_test_M1_covariate_baseline = m1,
    cindex_test_M3_full_model = m3,
    delta_cindex = delta,
    delta_cindex_ci_lower = lo,
    delta_cindex_ci_upper = hi,
    delta_cindex_p = p,
    delta_significant_existing = NA_character_
  )
}

read_incremental_for_disease <- function(disease) {
  disease_outdir <- file.path(
    base_dir,
    paste0("all_", disease, "_lepoch_model_performance")
  )
  
  summary_file <- file.path(
    disease_outdir,
    paste0("all_", disease, "_lepoch_summary.tsv")
  )
  
  incremental_file <- file.path(
    disease_outdir,
    paste0("all_", disease, "_lepoch_incremental_value_summary.tsv")
  )
  
  if (file.exists(summary_file)) {
    message("Reading combined model-performance summary for ", disease, ":")
    message("  ", summary_file)
    
    return(
      standardize_combined_summary(
        safe_read_tsv(summary_file),
        disease
      )
    )
  }
  
  if (file.exists(incremental_file)) {
    message("Reading combined incremental-value summary for ", disease, ":")
    message("  ", incremental_file)
    
    return(
      standardize_combined_summary(
        safe_read_tsv(incremental_file),
        disease
      )
    )
  }
  
  warning(
    "No combined summary found for ",
    disease,
    ". Falling back to per-clock *_incremental_value_delta_cindex.tsv."
  )
  
  make_manifest_one_disease(disease) %>%
    split(seq_len(nrow(.))) %>%
    purrr::map_dfr(~ read_incremental_from_delta_file(.x))
}

incremental_all <- purrr::map_dfr(
  diseases,
  read_incremental_for_disease
)

incremental_out <- file.path(
  outdir,
  "all_disease_lepoch_incremental_value_merged.tsv"
)

readr::write_tsv(incremental_all, incremental_out)

message("Merged incremental-value table:")
message("  ", incremental_out)

# ============================================================
# 8. Merge, classify, and save analysis table
# ============================================================

plot_tbl <- clock_manifest %>%
  left_join(
    incremental_all,
    by = c("disease", "clock_id")
  ) %>%
  left_join(
    qc_all,
    by = c("disease", "folder")
  ) %>%
  left_join(
    performance_counts %>%
      select(
        disease,
        clock_id,
        folder,
        n_total_json,
        n_events_json
      ),
    by = c("disease", "clock_id", "folder")
  ) %>%
  mutate(
    n_total_final = case_when(
      is.finite(n_total_summary) ~ n_total_summary,
      is.finite(n_total_json) ~ n_total_json,
      TRUE ~ NA_real_
    ),
    n_events_final = case_when(
      is.finite(n_events_summary) ~ n_events_summary,
      is.finite(n_events_json) ~ n_events_json,
      TRUE ~ NA_real_
    ),
    
    disease_label = factor(
      disease_labels[disease],
      levels = disease_labels[diseases]
    ),
    
    clock_label = factor(
      clock_label,
      levels = rev(clock_levels)
    ),
    
    modality = factor(
      modality,
      levels = c("MRI", "Proteomics", "Metabolomics")
    ),
    
    scale_issue = case_when(
      is.na(scale_issue) ~ NA,
      TRUE ~ scale_issue
    ),
    
    delta_significant_positive = infer_sig_positive(
      delta = delta_cindex,
      lo = delta_cindex_ci_lower,
      p = delta_cindex_p,
      existing_sig = delta_significant_existing,
      alpha = alpha_sig
    ),
    
    plot_status = case_when(
      !is.finite(delta_cindex) | is.na(scale_issue) ~ "Missing M3-M1 or QC",
      delta_significant_positive & !scale_issue ~ "Good: significant + stable scale",
      delta_significant_positive & scale_issue ~ "Significant but scale issue",
      !delta_significant_positive & !scale_issue ~ "Non-significant + stable scale",
      !delta_significant_positive & scale_issue ~ "Non-significant + scale issue",
      TRUE ~ "Missing M3-M1 or QC"
    ),
    
    plot_status = factor(
      plot_status,
      levels = names(shape_values)
    ),
    
    main_text_keep = plot_status == "Good: significant + stable scale",
    
    ci_text = format_ci_label(
      delta_cindex_ci_lower,
      delta_cindex_ci_upper
    ),
    
    n_cases_text = format_n_cases_label(
      n_total_final,
      n_events_final
    ),
    
    forest_label = paste0(
      n_cases_text,
      "; ",
      ci_text
    ),
    
    delta_label = ifelse(
      is.finite(delta_cindex),
      sprintf("%+.3f", delta_cindex),
      "NA"
    )
  )

plot_table_out <- file.path(
  outdir,
  "all_disease_lepoch_incremental_value_scale_qc_plot_table.tsv"
)

readr::write_tsv(plot_tbl, plot_table_out)

status_summary <- plot_tbl %>%
  count(disease, disease_label, plot_status, name = "n_clocks") %>%
  arrange(disease_label, plot_status)

status_summary_out <- file.path(
  outdir,
  "all_disease_lepoch_incremental_value_scale_qc_status_summary.tsv"
)

readr::write_tsv(status_summary, status_summary_out)

main_text_tbl <- plot_tbl %>%
  filter(main_text_keep) %>%
  arrange(disease_label, modality, desc(delta_cindex)) %>%
  select(
    disease,
    disease_label,
    modality,
    organ_label,
    clock_label,
    clock_id,
    folder,
    n_total_final,
    n_events_final,
    cindex_train,
    cindex_validation,
    cindex_test,
    cindex_test_M1_covariate_baseline,
    cindex_test_M3_full_model,
    delta_cindex,
    delta_cindex_ci_lower,
    delta_cindex_ci_upper,
    delta_cindex_p,
    scale_qc_status,
    scale_qc_reason,
    performance_file,
    delta_file
  )

main_text_out <- file.path(
  outdir,
  "all_disease_lepoch_main_text_good_clocks.tsv"
)

readr::write_tsv(main_text_tbl, main_text_out)

message("\nStatus summary:")
print(status_summary)

message("\nMain-text keep clocks:")
print(main_text_tbl, n = Inf)

# ============================================================
# 9. Figure 1: faceted forest plot with annotations
# ============================================================

# Axis scale is determined only by CI values:
#   xlim = -abs(max CI) to +abs(max CI)
x_vals_ci <- c(
  plot_tbl$delta_cindex_ci_lower,
  plot_tbl$delta_cindex_ci_upper
)

x_vals_ci <- x_vals_ci[is.finite(x_vals_ci)]

if (length(x_vals_ci) == 0) {
  x_vals_ci <- plot_tbl$delta_cindex
  x_vals_ci <- x_vals_ci[is.finite(x_vals_ci)]
}

if (length(x_vals_ci) == 0) {
  stop("No finite ΔC-index or CI values found.")
}

x_abs_ci <- max(abs(x_vals_ci), na.rm = TRUE)

# Guard against degenerate scales.
x_abs_ci <- max(x_abs_ci, 0.005)

# Labels are drawn outside plotting region with clip = "off".
# They do not determine the x-axis scale.
label_pad <- max(0.0015, x_abs_ci * 0.035)

plot_tbl_forest <- plot_tbl %>%
  mutate(
    ci_or_point_right = case_when(
      is.finite(delta_cindex_ci_upper) ~ delta_cindex_ci_upper,
      is.finite(delta_cindex) ~ delta_cindex,
      TRUE ~ NA_real_
    ),
    label_x = ci_or_point_right + label_pad
  )

p_forest <- ggplot(plot_tbl_forest, aes(y = clock_label)) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "#7B7D7D",
    linewidth = 0.45
  ) +
  geom_segment(
    data = plot_tbl_forest %>%
      filter(
        is.finite(delta_cindex_ci_lower),
        is.finite(delta_cindex_ci_upper)
      ),
    aes(
      x = delta_cindex_ci_lower,
      xend = delta_cindex_ci_upper,
      yend = clock_label,
      color = clock_label
    ),
    linewidth = 0.72,
    alpha = 0.78
  ) +
  geom_point(
    aes(
      x = delta_cindex,
      color = clock_label,
      shape = plot_status
    ),
    size = 3.35,
    stroke = 1.2,
    alpha = 0.98,
    na.rm = TRUE
  ) +
  geom_text(
    data = plot_tbl_forest %>%
      filter(is.finite(label_x)),
    aes(
      x = label_x,
      label = forest_label
    ),
    hjust = 0,
    size = 2.05,
    color = "#17202A",
    na.rm = TRUE
  ) +
  facet_wrap(
    ~ disease_label,
    nrow = 1
  ) +
  scale_color_manual(
    values = clock_cols,
    guide = "none"
  ) +
  scale_shape_manual(
    values = shape_values,
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 5),
    labels = function(x) sprintf("%+.3f", x),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  coord_cartesian(
    xlim = c(-x_abs_ci, x_abs_ci),
    clip = "off"
  ) +
  labs(
    title = "Incremental predictive value and year-scale QC across disease L’EPOCH clocks",
    subtitle = paste0(
      "Effect size is test-set ΔC-index = M3 full model − M1 covariate baseline; ",
      "horizontal bars show 95% CI. ",
      "The x-axis is fixed to −abs(max CI) to +abs(max CI)."
    ),
    x = "Test-set ΔC-index beyond covariates, M3 − M1",
    y = NULL,
    shape = "Clock class",
    caption = paste0(
      "Good clocks are significant for M3−M1 and pass year-scale QC. ",
      "Open triangles indicate significant clocks with scale issues. ",
      "Right-side labels show N, cases, and ΔC-index 95% CI. ",
      "Colors identify the 23 clocks and are kept consistent across diseases."
    )
  ) +
  theme_clock(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 8.2),
    axis.text.x = element_text(size = 8.4),
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    panel.spacing.x = unit(2.2, "lines"),
    plot.margin = margin(t = 8, r = 160, b = 8, l = 8)
  ) +
  guides(
    shape = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(
        color = "#17202A",
        size = 3.8,
        stroke = 1.2
      )
    )
  )

out_forest_pdf <- file.path(
  outdir,
  "all_disease_lepoch_incremental_value_scale_qc_forest.pdf"
)

out_forest_png <- file.path(
  outdir,
  "all_disease_lepoch_incremental_value_scale_qc_forest.png"
)

if (capabilities("cairo")) {
  ggsave(
    filename = out_forest_pdf,
    plot = p_forest,
    width = 27,
    height = 12.5,
    device = cairo_pdf
  )
} else {
  ggsave(
    filename = out_forest_pdf,
    plot = p_forest,
    width = 27,
    height = 12.5
  )
}

ggsave(
  filename = out_forest_png,
  plot = p_forest,
  width = 27,
  height = 12.5,
  dpi = 320
)

# ============================================================
# 10. Figure 2: disease-by-clock matrix
# ============================================================

p_matrix <- plot_tbl %>%
  mutate(
    delta_abs = ifelse(is.finite(delta_cindex), abs(delta_cindex), NA_real_)
  ) %>%
  ggplot(aes(x = disease_label, y = clock_label)) +
  geom_tile(
    fill = "#F8F9F9",
    color = "white",
    linewidth = 0.35
  ) +
  geom_point(
    aes(
      color = clock_label,
      shape = plot_status,
      size = delta_abs
    ),
    stroke = 1.1,
    alpha = 0.96,
    na.rm = TRUE
  ) +
  geom_text(
    aes(
      label = ifelse(is.finite(delta_cindex), sprintf("%+.3f", delta_cindex), "")
    ),
    nudge_x = 0.29,
    size = 2.25,
    color = "#17202A",
    na.rm = TRUE
  ) +
  scale_color_manual(
    values = clock_cols,
    guide = "none"
  ) +
  scale_shape_manual(
    values = shape_values,
    drop = FALSE
  ) +
  scale_size_continuous(
    range = c(1.8, 5.2),
    name = "|ΔC-index|",
    labels = number_format(accuracy = 0.001)
  ) +
  labs(
    title = "Disease-by-clock map of incremental value and year-scale QC",
    subtitle = paste0(
      "Rows are the 23 disease clocks and columns are the five diseases. ",
      "Numbers show test-set ΔC-index, M3 − M1. ",
      "Shape indicates whether the clock is significant and scale-stable."
    ),
    x = NULL,
    y = NULL,
    shape = "Clock class"
  ) +
  theme_clock(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 8.2),
    axis.text.x = element_text(face = "bold"),
    legend.position = "bottom",
    legend.box = "vertical"
  ) +
  guides(
    shape = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(color = "#17202A", size = 3.6)
    ),
    size = guide_legend(nrow = 1)
  )

out_matrix_pdf <- file.path(
  outdir,
  "all_disease_lepoch_incremental_value_scale_qc_matrix.pdf"
)

out_matrix_png <- file.path(
  outdir,
  "all_disease_lepoch_incremental_value_scale_qc_matrix.png"
)

if (capabilities("cairo")) {
  ggsave(
    filename = out_matrix_pdf,
    plot = p_matrix,
    width = 14.5,
    height = 12.5,
    device = cairo_pdf
  )
} else {
  ggsave(
    filename = out_matrix_pdf,
    plot = p_matrix,
    width = 14.5,
    height = 12.5
  )
}

ggsave(
  filename = out_matrix_png,
  plot = p_matrix,
  width = 14.5,
  height = 12.5,
  dpi = 320
)

# ============================================================
# 11. Figure 3: combined figure
# ============================================================

combined_fig <- p_forest / p_matrix +
  plot_layout(heights = c(1.1, 1.0)) +
  plot_annotation(
    title = "Cross-disease L’EPOCH clock prioritization",
    subtitle = "Combining incremental predictive value beyond covariates with year-scale QC across 115 disease clocks.",
    caption = "Main-text candidates are clocks with significant positive M3−M1 ΔC-index and stable year-scale QC.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 20, color = "#17202A"),
      plot.subtitle = element_text(size = 11, color = "#566573"),
      plot.caption = element_text(size = 9, color = "#7B7D7D")
    )
  )

out_combined_pdf <- file.path(
  outdir,
  "all_disease_lepoch_incremental_value_scale_qc_combined.pdf"
)

out_combined_png <- file.path(
  outdir,
  "all_disease_lepoch_incremental_value_scale_qc_combined.png"
)

if (capabilities("cairo")) {
  ggsave(
    filename = out_combined_pdf,
    plot = combined_fig,
    width = 27,
    height = 24,
    device = cairo_pdf
  )
} else {
  ggsave(
    filename = out_combined_pdf,
    plot = combined_fig,
    width = 27,
    height = 24
  )
}

ggsave(
  filename = out_combined_png,
  plot = combined_fig,
  width = 27,
  height = 24,
  dpi = 320
)

# ============================================================
# 12. Final messages
# ============================================================

message("\n============================================================")
message("Finished disease-clock incremental value + scale-QC plotting.")
message("============================================================")
message("Output directory:")
message("  ", outdir)

message("\nKey tables:")
message("  Manifest file check:")
message("    ", manifest_check_out)
message("  Merged QC table:")
message("    ", qc_out)
message("  Performance JSON counts:")
message("    ", performance_counts_out)
message("  Merged incremental-value table:")
message("    ", incremental_out)
message("  Plot table:")
message("    ", plot_table_out)
message("  Status summary:")
message("    ", status_summary_out)
message("  Main-text good clocks:")
message("    ", main_text_out)

message("\nFigures:")
message("  Forest with symmetric CI-based x-axis and annotation:")
message("    ", out_forest_pdf)
message("    ", out_forest_png)
message("  Matrix:")
message("    ", out_matrix_pdf)
message("    ", out_matrix_png)
message("  Combined:")
message("    ", out_combined_pdf)
message("    ", out_combined_png)

message("\nInterpretation guide:")
message("  - Good: significant + stable scale = recommended main-text clocks.")
message("  - Open triangle = significant M3-M1 but scale issue.")
message("  - Open circle = non-significant but stable scale.")
message("  - Cross = non-significant and scale issue.")
message("  - Forest x-axis is fixed to -abs(max CI) to +abs(max CI).")
message("  - Right-side forest labels show N, cases, and ΔC-index 95% CI.")
message("============================================================")