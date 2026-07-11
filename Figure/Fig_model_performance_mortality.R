# ============================================================
# Improved combined overview of all mortality clocks
# MRI + Proteomics + Metabolomics
# Independent RStudio-ready script
#
# Main goals:
#   A. Split-wise discrimination using bar plots
#   B. Train-test optimism only
#   C. Incremental value beyond covariates (M3 - M1)
#   D. Practical model-quality map
#   E. Test-set mortality separation by training-defined risk quartile
#
# This script reads the original per-clock outputs directly
# (predictions TSV, performance JSON, model-comparison TSV,
#  incremental-value TSV), so it can be run independently.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(jsonlite)
  library(patchwork)
  library(scales)
  library(glue)
})

# ============================================================
# 1. Settings
# ============================================================

possible_base_dirs <- c(
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock",
  "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock",
  getwd()
)
base_dir <- possible_base_dirs[file.exists(possible_base_dirs)][1]
if (is.na(base_dir) || is.null(base_dir)) {
  stop("Please set base_dir manually.")
}

skip_missing <- TRUE
combined_outdir <- file.path(base_dir, "all_mortality_clock_model_performance")
dir.create(combined_outdir, recursive = TRUE, showWarnings = FALSE)

# If TRUE, panel E is limited to 15 years to match the individual figures better.
km_time_limit <- 15

message("Base directory: ", base_dir)
message("Combined output directory: ", combined_outdir)

# ============================================================
# 2. Clock manifest
# ============================================================

make_manifest <- function(base_dir) {
  mri_organs <- c("brain", "heart", "adipose", "kidney", "liver", "pancreas", "spleen")
  proteomics_organs <- c(
    "Reproductive_female", "Pulmonary", "Heart", "Brain", "Eye", "Hepatic",
    "Renal", "Reproductive_male", "Endocrine", "Immune", "Skin"
  )
  metabolomics_organs <- c("Endocrine", "Digestive", "Hepatic", "Immune", "Metabolic")
  
  bind_rows(
    tibble(
      modality = "MRI",
      modality_key = "mri",
      organ_folder_name = mri_organs,
      organ_key = stringr::str_to_lower(mri_organs),
      folder = paste0(stringr::str_to_lower(mri_organs), "_mri_mortality_clock"),
      prefix = paste0(stringr::str_to_lower(mri_organs), "_mri_mortality_clock")
    ),
    tibble(
      modality = "Proteomics",
      modality_key = "proteomics",
      organ_folder_name = proteomics_organs,
      organ_key = stringr::str_to_lower(proteomics_organs),
      folder = paste0(proteomics_organs, "_proteomics_mortality_clock"),
      prefix = paste0(stringr::str_to_lower(proteomics_organs), "_proteomics_mortality_clock")
    ),
    tibble(
      modality = "Metabolomics",
      modality_key = "metabolomics",
      organ_folder_name = metabolomics_organs,
      organ_key = stringr::str_to_lower(metabolomics_organs),
      folder = paste0(metabolomics_organs, "_metabolomics_mortality_clock"),
      prefix = paste0(stringr::str_to_lower(metabolomics_organs), "_metabolomics_mortality_clock")
    )
  ) %>%
    mutate(
      organ_label = organ_folder_name %>%
        stringr::str_replace_all("_", " ") %>%
        stringr::str_to_sentence(),
      clock_label = paste(organ_label, modality),
      clock_id = paste(organ_key, modality_key, sep = "__"),
      clock_dir = file.path(base_dir, folder),
      prediction_file = file.path(clock_dir, paste0(prefix, "_predictions.tsv")),
      performance_file = file.path(clock_dir, paste0(prefix, "_performance.json")),
      model_comparison_file = file.path(clock_dir, paste0(prefix, "_model_comparison.tsv")),
      delta_file = file.path(clock_dir, paste0(prefix, "_incremental_value_delta_cindex.tsv"))
    )
}

clock_manifest <- make_manifest(base_dir)
message("Expected clocks: ", nrow(clock_manifest))

# ============================================================
# 3. General helpers
# ============================================================

theme_clock <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(color = "#17202A"),
      plot.title = element_text(face = "bold", size = base_size + 3, margin = margin(b = 6)),
      plot.subtitle = element_text(size = base_size, color = "#566573", margin = margin(b = 8)),
      plot.caption = element_text(size = base_size - 3, color = "#7B7D7D"),
      axis.title = element_text(face = "bold", size = base_size - 0.2),
      axis.text = element_text(color = "#2C3E50"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "#EFF3F6", linewidth = 0.35),
      panel.grid.major.y = element_line(color = "#EFF3F6", linewidth = 0.35),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.position = "bottom",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

split_cols <- c("train" = "#2E86AB", "validation" = "#F18F01", "test" = "#6A994E")
modality_cols <- c("MRI" = "#2E86AB", "Proteomics" = "#8E44AD", "Metabolomics" = "#D35400")
significance_cols <- c("Significant" = "#2C7BB6", "Not significant" = "#9AA0A6", "Missing" = "#D0D3D4")
quartile_cols <- c(
  "Q1 lowest risk" = "#1B9E77",
  "Q2" = "#7570B3",
  "Q3" = "#E6AB02",
  "Q4 highest risk" = "#D95F02"
)

json_get_num <- function(x, field) {
  if (is.null(x[[field]])) return(NA_real_)
  out <- suppressWarnings(as.numeric(x[[field]]))
  if (length(out) == 0) return(NA_real_)
  out[[1]]
}

normalize_event_column <- function(x) {
  x_chr <- as.character(x)
  dplyr::case_when(
    x_chr %in% c("TRUE", "True", "true", "1", "1.0") ~ TRUE,
    x_chr %in% c("FALSE", "False", "false", "0", "0.0") ~ FALSE,
    TRUE ~ as.logical(x)
  )
}

calc_cindex <- function(df) {
  if (nrow(df) < 5 || sum(df$event, na.rm = TRUE) < 2) return(NA_real_)
  fit <- survival::concordance(
    survival::Surv(time_years, event) ~ I(-risk_score),
    data = df
  )
  as.numeric(fit$concordance)
}

first_existing_col <- function(df, candidates) {
  hits <- candidates[candidates %in% colnames(df)]
  if (length(hits) == 0) return(NA_character_)
  hits[[1]]
}

read_delta_from_files <- function(meta, perf) {
  delta <- json_get_num(perf, "delta_cindex_test_M3_vs_M1")
  lo <- json_get_num(perf, "delta_cindex_test_M3_vs_M1_ci_lower")
  hi <- json_get_num(perf, "delta_cindex_test_M3_vs_M1_ci_upper")
  p <- json_get_num(perf, "delta_cindex_test_M3_vs_M1_p_two_sided")
  
  if ((!is.finite(delta) || !is.finite(lo) || !is.finite(hi)) && file.exists(meta$delta_file)) {
    delta_tbl <- readr::read_tsv(meta$delta_file, show_col_types = FALSE, progress = FALSE)
    if (nrow(delta_tbl) > 0) {
      if ("delta_cindex" %in% colnames(delta_tbl)) delta <- as.numeric(delta_tbl$delta_cindex[[1]])
      if ("delta_cindex_ci_lower" %in% colnames(delta_tbl)) lo <- as.numeric(delta_tbl$delta_cindex_ci_lower[[1]])
      if ("delta_cindex_ci_upper" %in% colnames(delta_tbl)) hi <- as.numeric(delta_tbl$delta_cindex_ci_upper[[1]])
      if ("empirical_p_two_sided_delta_not_equal_0" %in% colnames(delta_tbl)) {
        p <- as.numeric(delta_tbl$empirical_p_two_sided_delta_not_equal_0[[1]])
      }
    }
  }
  
  tibble(
    delta_cindex_test_M3_vs_M1 = delta,
    delta_cindex_test_M3_vs_M1_ci_lower = lo,
    delta_cindex_test_M3_vs_M1_ci_upper = hi,
    delta_cindex_test_M3_vs_M1_p_two_sided = p
  )
}

tidy_survfit_mortality <- function(fit) {
  ss <- summary(fit)
  tibble(
    time = ss$time,
    surv = ss$surv,
    surv_lower = ss$lower,
    surv_upper = ss$upper,
    strata = if (is.null(ss$strata)) "All" else as.character(ss$strata)
  ) %>%
    mutate(
      strata = stringr::str_replace(strata, "^risk_quartile=", ""),
      cum_mortality = 1 - surv,
      cum_lower = 1 - surv_upper,
      cum_upper = 1 - surv_lower
    )
}

# ============================================================
# 4. Per-clock reader
# ============================================================

read_one_clock <- function(meta_row, skip_missing = TRUE) {
  meta <- as.list(meta_row)
  
  if (!file.exists(meta$clock_dir) || !file.exists(meta$prediction_file) || !file.exists(meta$performance_file)) {
    msg <- paste0(
      "Missing input for ", meta$clock_label, ":\n",
      "  clock_dir: ", meta$clock_dir, " exists=", file.exists(meta$clock_dir), "\n",
      "  predictions: ", meta$prediction_file, " exists=", file.exists(meta$prediction_file), "\n",
      "  performance: ", meta$performance_file, " exists=", file.exists(meta$performance_file)
    )
    if (skip_missing) {
      warning(msg)
      return(NULL)
    } else {
      stop(msg)
    }
  }
  
  pred <- readr::read_tsv(meta$prediction_file, show_col_types = FALSE, progress = FALSE)
  perf <- jsonlite::fromJSON(meta$performance_file)
  
  feature_token <- paste0(meta$organ_key, "_", meta$modality_key)
  risk_score_col <- paste0(feature_token, "_mortality_risk_score")
  
  if (!risk_score_col %in% colnames(pred)) {
    fallback_cols <- grep("_mortality_risk_score$", colnames(pred), value = TRUE)
    if (length(fallback_cols) == 1) {
      risk_score_col <- fallback_cols[[1]]
    } else {
      stop("Could not identify mortality risk-score column for ", meta$clock_label)
    }
  }
  
  required_cols <- c("participant_id", "time_years", "event", "split", risk_score_col)
  missing_cols <- setdiff(required_cols, colnames(pred))
  if (length(missing_cols) > 0) {
    stop("Missing required columns for ", meta$clock_label, ": ", paste(missing_cols, collapse = ", "))
  }
  
  pred <- pred %>%
    mutate(
      event = normalize_event_column(event),
      split = factor(split, levels = c("train", "validation", "test")),
      risk_score = .data[[risk_score_col]]
    ) %>%
    filter(!is.na(time_years), !is.na(event), !is.na(split), !is.na(risk_score))
  
  cindex_tbl <- pred %>%
    group_by(split) %>%
    summarise(cindex_calc = calc_cindex(cur_data()), .groups = "drop") %>%
    mutate(
      cindex_json = c(
        json_get_num(perf, "cindex_train"),
        json_get_num(perf, "cindex_validation"),
        json_get_num(perf, "cindex_test")
      ),
      cindex = dplyr::coalesce(cindex_json, cindex_calc),
      clock_id = meta$clock_id,
      modality = meta$modality,
      modality_key = meta$modality_key,
      organ_key = meta$organ_key,
      organ_label = meta$organ_label,
      clock_label = meta$clock_label,
      folder = meta$folder,
      prefix = meta$prefix
    )
  
  cindex_train <- cindex_tbl %>% filter(split == "train") %>% pull(cindex)
  cindex_validation <- cindex_tbl %>% filter(split == "validation") %>% pull(cindex)
  cindex_test <- cindex_tbl %>% filter(split == "test") %>% pull(cindex)
  
  delta_tbl <- read_delta_from_files(meta, perf)
  
  m3_json_field <- paste0("cindex_test_M3_full_covariates_plus_", meta$organ_key, "_", meta$modality_key)
  m3_cindex <- json_get_num(perf, m3_json_field)
  if (!is.finite(m3_cindex) && file.exists(meta$model_comparison_file)) {
    mc <- readr::read_tsv(meta$model_comparison_file, show_col_types = FALSE, progress = FALSE)
    if (all(c("split", "model", "cindex") %in% colnames(mc))) {
      m3_row <- mc %>% filter(split == "test", stringr::str_detect(model, "^M3_"))
      if (nrow(m3_row) > 0) m3_cindex <- as.numeric(m3_row$cindex[[1]])
    }
  }
  
  summary_row <- tibble(
    clock_id = meta$clock_id,
    modality = meta$modality,
    modality_key = meta$modality_key,
    organ_key = meta$organ_key,
    organ_label = meta$organ_label,
    clock_label = meta$clock_label,
    folder = meta$folder,
    prefix = meta$prefix,
    n_total = nrow(pred),
    n_events_total = sum(pred$event, na.rm = TRUE),
    median_followup_years = median(pred$time_years, na.rm = TRUE),
    cindex_train = cindex_train,
    cindex_validation = cindex_validation,
    cindex_test = cindex_test,
    train_minus_test = cindex_train - cindex_test,
    validation_minus_test = cindex_validation - cindex_test,
    cindex_test_M1_covariate_baseline = json_get_num(perf, "cindex_test_M1_covariate_baseline"),
    cindex_test_M3_full_model = m3_cindex
  ) %>%
    bind_cols(delta_tbl) %>%
    mutate(
      delta_significant = case_when(
        is.na(delta_cindex_test_M3_vs_M1) ~ "Missing",
        is.finite(delta_cindex_test_M3_vs_M1_ci_lower) & delta_cindex_test_M3_vs_M1_ci_lower > 0 ~ "Significant",
        is.finite(delta_cindex_test_M3_vs_M1_p_two_sided) & delta_cindex_test_M3_vs_M1_p_two_sided < 0.05 & delta_cindex_test_M3_vs_M1 > 0 ~ "Significant",
        TRUE ~ "Not significant"
      )
    )
  
  # Risk quartiles are defined from the training split and then applied to the test set.
  train_risk <- pred %>% filter(split == "train") %>% pull(risk_score)
  risk_breaks <- quantile(train_risk, probs = seq(0, 1, by = 0.25), na.rm = TRUE, type = 8)
  risk_breaks[1] <- -Inf
  risk_breaks[length(risk_breaks)] <- Inf
  
  if (length(unique(risk_breaks)) < length(risk_breaks)) {
    pred <- pred %>%
      group_by(split) %>%
      mutate(
        risk_quartile = ntile(risk_score, 4),
        risk_quartile = factor(risk_quartile, levels = 1:4, labels = c("Q1 lowest risk", "Q2", "Q3", "Q4 highest risk"))
      ) %>%
      ungroup()
  } else {
    pred <- pred %>%
      mutate(
        risk_quartile = cut(
          risk_score,
          breaks = risk_breaks,
          include.lowest = TRUE,
          labels = c("Q1 lowest risk", "Q2", "Q3", "Q4 highest risk")
        )
      )
  }
  
  test_for_km <- pred %>%
    filter(split == "test", !is.na(risk_quartile))
  
  if (nrow(test_for_km) > 0 && sum(test_for_km$event, na.rm = TRUE) > 0 && dplyr::n_distinct(test_for_km$risk_quartile) >= 2) {
    km_fit <- survival::survfit(survival::Surv(time_years, event) ~ risk_quartile, data = test_for_km)
    km_tbl <- tidy_survfit_mortality(km_fit) %>%
      transmute(
        clock_id = meta$clock_id,
        modality = meta$modality,
        organ_key = meta$organ_key,
        organ_label = meta$organ_label,
        clock_label = meta$clock_label,
        time = time,
        risk_quartile = factor(strata, levels = c("Q1 lowest risk", "Q2", "Q3", "Q4 highest risk")),
        cum_mortality = cum_mortality,
        cum_lower = cum_lower,
        cum_upper = cum_upper
      )
  } else {
    km_tbl <- tibble(
      clock_id = character(), modality = character(), organ_key = character(), organ_label = character(),
      clock_label = character(), time = numeric(), risk_quartile = factor(levels = c("Q1 lowest risk", "Q2", "Q3", "Q4 highest risk")),
      cum_mortality = numeric(), cum_lower = numeric(), cum_upper = numeric()
    )
  }
  
  list(
    summary = summary_row,
    cindex = cindex_tbl,
    km = km_tbl
  )
}

# ============================================================
# 5. Read all clocks
# ============================================================

results <- purrr::map(
  seq_len(nrow(clock_manifest)),
  ~ read_one_clock(clock_manifest[.x, ], skip_missing = skip_missing)
)
results <- purrr::compact(results)

if (length(results) == 0) {
  stop("No mortality clocks were successfully processed.")
}

combined_summary <- purrr::map_dfr(results, "summary")
combined_cindex <- purrr::map_dfr(results, "cindex")
combined_km <- purrr::map_dfr(results, "km")

# Save combined tables.
out_combined_summary <- file.path(combined_outdir, "all_mortality_clock_summary.tsv")
out_combined_cindex <- file.path(combined_outdir, "all_mortality_clock_cindex_summary.tsv")
out_combined_incremental <- file.path(combined_outdir, "all_mortality_clock_incremental_value_summary.tsv")
out_combined_test_km <- file.path(combined_outdir, "all_mortality_clock_test_km_curves.tsv")

readr::write_tsv(combined_summary, out_combined_summary)
readr::write_tsv(combined_cindex, out_combined_cindex)
readr::write_tsv(
  combined_summary %>%
    select(
      clock_id, modality, organ_label, clock_label,
      cindex_test_M1_covariate_baseline,
      cindex_test_M3_full_model,
      delta_cindex_test_M3_vs_M1,
      delta_cindex_test_M3_vs_M1_ci_lower,
      delta_cindex_test_M3_vs_M1_ci_upper,
      delta_cindex_test_M3_vs_M1_p_two_sided,
      delta_significant
    ),
  out_combined_incremental
)
if (nrow(combined_km) > 0) {
  readr::write_tsv(combined_km, out_combined_test_km)
}

# ============================================================
# 6. Ordering and plotting tables
# ============================================================

clock_order_tbl <- combined_summary %>%
  mutate(modality = factor(modality, levels = c("MRI", "Proteomics", "Metabolomics"))) %>%
  arrange(modality, desc(cindex_test)) %>%
  mutate(clock_plot = paste0(clock_label, "  "))

clock_levels <- rev(clock_order_tbl$clock_plot)

combined_cindex_plot_tbl <- combined_cindex %>%
  left_join(clock_order_tbl %>% select(clock_id, clock_plot), by = "clock_id") %>%
  mutate(
    clock_plot = factor(clock_plot, levels = clock_levels),
    split = factor(split, levels = c("train", "validation", "test")),
    modality = factor(modality, levels = c("MRI", "Proteomics", "Metabolomics"))
  )

combined_summary_plot_tbl <- combined_summary %>%
  left_join(clock_order_tbl %>% select(clock_id, clock_plot), by = "clock_id") %>%
  mutate(
    clock_plot = factor(clock_plot, levels = clock_levels),
    modality = factor(modality, levels = c("MRI", "Proteomics", "Metabolomics")),
    delta_significant = factor(delta_significant, levels = c("Significant", "Not significant", "Missing")),
    train_optimism_flag = case_when(
      train_minus_test > 0.05 ~ "High",
      train_minus_test > 0.02 ~ "Mild",
      TRUE ~ "Low/none"
    ),
    train_optimism_flag = factor(train_optimism_flag, levels = c("Low/none", "Mild", "High"))
  )

combined_km_plot_tbl <- combined_km %>%
  left_join(clock_order_tbl %>% select(clock_id, clock_plot), by = "clock_id") %>%
  mutate(
    clock_plot = factor(clock_plot, levels = clock_levels),
    modality = factor(modality, levels = c("MRI", "Proteomics", "Metabolomics")),
    risk_quartile = factor(risk_quartile, levels = c("Q1 lowest risk", "Q2", "Q3", "Q4 highest risk"))
  )

# ============================================================
# 7. Panel A: split-wise discrimination (bar plot)
# ============================================================

p_cindex_bar <- combined_cindex_plot_tbl %>%
  ggplot(aes(x = clock_plot, y = cindex, fill = split)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "#BFC9CA", linewidth = 0.45) +
  geom_col(position = position_dodge(width = 0.76), width = 0.70) +
  coord_flip() +
  scale_fill_manual(values = split_cols, drop = FALSE) +
  labs(
    title = "A. Discrimination across train, validation, and test",
    subtitle = "Bar plot of C-index from the original clock outputs",
    x = NULL,
    y = "C-index",
    fill = "Split"
  ) +
  theme_clock(base_size = 11) +
  theme(axis.text.y = element_text(size = 8.4))

# ============================================================
# 8. Panel B: train-test optimism only
# ============================================================

p_optimism <- combined_summary_plot_tbl %>%
  ggplot(aes(x = train_minus_test, y = clock_plot, color = modality)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#7B7D7D", linewidth = 0.45) +
  geom_vline(xintercept = 0.02, linetype = "dotted", color = "#B9770E", linewidth = 0.45) +
  geom_vline(xintercept = 0.05, linetype = "dotted", color = "#922B21", linewidth = 0.45) +
  geom_segment(aes(x = 0, xend = train_minus_test, y = clock_plot, yend = clock_plot), alpha = 0.65, linewidth = 0.55) +
  geom_point(size = 2.2, alpha = 0.95) +
  scale_color_manual(values = modality_cols) +
  labs(
    title = "B. Train-test optimism",
    subtitle = "Dotted lines mark +0.02 and +0.05 C-index gaps",
    x = "Train - test C-index gap",
    y = NULL,
    color = "Modality"
  ) +
  theme_clock(base_size = 11) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# ============================================================
# 9. Panel C: incremental value beyond covariates
# ============================================================

p_delta <- combined_summary_plot_tbl %>%
  ggplot(aes(x = delta_cindex_test_M3_vs_M1, y = clock_plot, color = modality)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#7B7D7D", linewidth = 0.45) +
  geom_errorbarh(
    aes(xmin = delta_cindex_test_M3_vs_M1_ci_lower, xmax = delta_cindex_test_M3_vs_M1_ci_upper),
    height = 0.18,
    linewidth = 0.55,
    na.rm = TRUE
  ) +
  geom_point(aes(size = cindex_test), alpha = 0.95, na.rm = TRUE) +
  scale_color_manual(values = modality_cols, drop = FALSE) +
  scale_size_continuous(range = c(1.8, 4.2), limits = range(combined_summary_plot_tbl$cindex_test, na.rm = TRUE)) +
  labs(
    title = "C. Incremental value beyond covariates",
    subtitle = "Test-set ΔC-index = M3 full model - M1 covariate baseline; horizontal bars show bootstrap 95% CI",
    x = "ΔC-index on test set",
    y = NULL,
    color = "Modality",
    size = "Test C-index"
  ) +
  theme_clock(base_size = 11) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# ============================================================
# 10. Panel D: practical quality map
# ============================================================

base_quality <- combined_summary_plot_tbl %>%
  ggplot(aes(x = train_minus_test, y = delta_cindex_test_M3_vs_M1, color = modality)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#BDC3C7", linewidth = 0.45) +
  geom_vline(xintercept = 0.05, linetype = "dotted", color = "#922B21", linewidth = 0.45) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#BDC3C7", linewidth = 0.45) +
  geom_point(aes(size = cindex_test, shape = delta_significant), alpha = 0.92) +
  scale_color_manual(values = modality_cols) +
  scale_shape_manual(values = c("Significant" = 16, "Not significant" = 1, "Missing" = 4), drop = FALSE) +
  scale_size_continuous(range = c(2.3, 5.2)) +
  labs(
    title = "D. Practical model-quality map",
    subtitle = "Ideal clocks appear near low optimism and positive ΔC-index",
    x = "Train - test C-index gap",
    y = "Test-set ΔC-index, M3 - M1",
    color = "Modality",
    shape = "Incremental value",
    size = "Test C-index"
  ) +
  theme_clock(base_size = 11)

if (requireNamespace("ggrepel", quietly = TRUE)) {
  p_quality <- base_quality +
    ggrepel::geom_text_repel(
      aes(label = organ_label),
      size = 2.8,
      max.overlaps = 40,
      min.segment.length = 0,
      box.padding = 0.25,
      show.legend = FALSE
    )
} else {
  warning("Package ggrepel is not installed. Panel D will use geom_text with check_overlap instead.")
  p_quality <- base_quality +
    geom_text(aes(label = organ_label), size = 2.5, check_overlap = TRUE, vjust = -0.7, show.legend = FALSE)
}

# ============================================================
# 11. Panel E: test-set mortality separation for all clocks
# ============================================================

if (nrow(combined_km_plot_tbl) > 0) {
  km_tbl_plot <- combined_km_plot_tbl %>%
    filter(is.finite(time), is.finite(cum_mortality)) %>%
    mutate(
      time = pmin(time, km_time_limit),
      cum_mortality = pmax(0, cum_mortality),
      cum_lower = pmax(0, cum_lower),
      cum_upper = pmax(cum_lower, cum_upper)
    )
  
  km_ymax <- max(km_tbl_plot$cum_upper, na.rm = TRUE)
  if (!is.finite(km_ymax) || km_ymax <= 0) km_ymax <- 0.30
  km_ymax <- min(max(0.10, km_ymax * 1.04), 0.40)
  
  p_km <- km_tbl_plot %>%
    ggplot(aes(x = time, y = cum_mortality, color = risk_quartile, fill = risk_quartile)) +
    geom_ribbon(aes(ymin = cum_lower, ymax = cum_upper), alpha = 0.12, color = NA) +
    geom_step(linewidth = 0.82) +
    facet_wrap(~ clock_plot, ncol = 4) +
    scale_color_manual(values = quartile_cols, drop = FALSE) +
    scale_fill_manual(values = quartile_cols, drop = FALSE) +
    scale_x_continuous(limits = c(0, km_time_limit), expand = expansion(mult = c(0, 0.02))) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, km_ymax), expand = expansion(mult = c(0, 0.02))) +
    labs(
      title = "E. Test-set mortality separation",
      subtitle = "Kaplan-Meier cumulative mortality by training-defined risk quartile",
      x = "Years after clock assessment",
      y = "Cumulative mortality",
      color = NULL,
      fill = NULL
    ) +
    theme_clock(base_size = 11) +
    theme(
      strip.text = element_text(size = 9, face = "bold"),
      legend.position = "bottom"
    )
} else {
  p_km <- ggplot() +
    annotate("text", x = 0, y = 0, label = "No test-set KM data available.", size = 5, fontface = "bold") +
    xlim(-1, 1) + ylim(-1, 1) +
    theme_void() +
    labs(title = "E. Test-set mortality separation")
}

# ============================================================
# 12. Assemble and save the combined overview figure
# ============================================================

combined_fig <-
  (p_cindex_bar | p_optimism | p_delta) /
  p_quality /
  p_km +
  plot_layout(heights = c(2.15, 1.35, 3.15)) +
  plot_annotation(
    title = "Performance summary of MRI, proteomics, and metabolomics mortality clocks",
    subtitle = glue(
      "Processed {nrow(combined_summary)} of {nrow(clock_manifest)} expected clocks. ",
      "Panels summarize discrimination, overfitting risk, incremental value beyond covariates, ",
      "and test-set mortality separation."
    ),
    caption = paste0(
      "M1 = covariate baseline. M3 = covariates plus imaging/protein/metabolite features. ",
      "The script reads the original per-clock outputs directly. Missing or still-running clocks are skipped automatically."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 20, color = "#17202A"),
      plot.subtitle = element_text(size = 11, color = "#566573"),
      plot.caption = element_text(size = 9, color = "#7B7D7D")
    )
  )

out_combined_fig_pdf <- file.path(combined_outdir, "all_mortality_clock_model_performance_overview_improved.pdf")
out_combined_fig_png <- file.path(combined_outdir, "all_mortality_clock_model_performance_overview_improved.png")

if (capabilities("cairo")) {
  ggsave(filename = out_combined_fig_pdf, plot = combined_fig, width = 20, height = 24, device = cairo_pdf)
} else {
  ggsave(filename = out_combined_fig_pdf, plot = combined_fig, width = 20, height = 24)
}

ggsave(filename = out_combined_fig_png, plot = combined_fig, width = 20, height = 24, dpi = 320)

message("\nSaved combined overview figure:")
message("  ", out_combined_fig_pdf)
message("  ", out_combined_fig_png)
message("\nSaved combined summary tables:")
message("  ", out_combined_summary)
message("  ", out_combined_cindex)
message("  ", out_combined_incremental)
if (file.exists(out_combined_test_km)) message("  ", out_combined_test_km)

message("\n===== Combined mortality clock summary =====")
print(
  combined_summary %>%
    arrange(factor(modality, levels = c("MRI", "Proteomics", "Metabolomics")), desc(cindex_test)) %>%
    select(
      modality,
      organ_label,
      n_total,
      n_events_total,
      median_followup_years,
      cindex_train,
      cindex_validation,
      cindex_test,
      train_minus_test,
      validation_minus_test,
      delta_cindex_test_M3_vs_M1,
      delta_cindex_test_M3_vs_M1_ci_lower,
      delta_cindex_test_M3_vs_M1_ci_upper,
      delta_significant
    )
)

message("\nInterpretation guide:")
message("  - Similar train, validation, and test C-index values suggest limited overfitting.")
message("  - Large positive train-test gaps, especially >0.05, suggest overfitting or split instability.")
message("  - Positive M3-M1 delta indicates that imaging/protein/metabolite features add value beyond covariates.")
message("  - A bootstrap CI for M3-M1 entirely above zero provides stronger evidence of incremental value.")
message("  - Clear Q1-to-Q4 separation in panel E supports useful risk stratification on the held-out test set.")
