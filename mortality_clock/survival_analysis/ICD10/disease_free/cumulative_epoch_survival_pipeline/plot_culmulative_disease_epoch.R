#!/usr/bin/env Rscript

# =============================================================================
# Plot cumulative C-index for one disease endpoint
# Proteomics + Metabolomics mortality EPOCH cumulative survival analysis
#
# This script plots how Harrell's C-index changes as the 11 proteomics and
# 4 metabolomics mortality EPOCH clocks are added cumulatively to the
# clinical baseline model.
#
# C-index options:
#   C_INDEX_TYPE=cv        : 5-fold cross-validated C-index, default
#   C_INDEX_TYPE=apparent  : apparent/in-sample Harrell C-index
#
# Example:
#   DISEASE_CODE=N183 DISEASE_NAME="CKD" C_INDEX_TYPE=cv Rscript plot_culmulative_disease_epoch.R
#   DISEASE_CODE=N183 DISEASE_NAME="CKD" C_INDEX_TYPE=apparent Rscript plot_culmulative_disease_epoch.R
# =============================================================================

options(stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# User-editable settings
# -----------------------------------------------------------------------------

DISEASE_CODE <- Sys.getenv("DISEASE_CODE", unset = "I500")

DISEASE_NAME <- Sys.getenv(
  "DISEASE_NAME",
  unset = "CHF"
)

INPUT_DIR <- Sys.getenv(
  "INPUT_DIR",
  unset = "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/mortality_clock/SA/output_cumulative_EPOCH_PM/disease_free_cv"
)

INPUT_FILE <- Sys.getenv(
  "INPUT_FILE",
  unset = file.path(
    INPUT_DIR,
    paste0("cox_cumulative_EPOCH_PM_", DISEASE_CODE, ".tsv")
  )
)

OUTPUT_DIR <- Sys.getenv(
  "OUTPUT_DIR",
  unset = file.path(INPUT_DIR, "figures_cumulative_cindex")
)

# Options:
#   "cv"       = plot 5-fold cross-validated Harrell C-index
#   "apparent" = plot apparent/in-sample Harrell C-index
C_INDEX_TYPE <- tolower(Sys.getenv("C_INDEX_TYPE", unset = "cv"))

# Sequential likelihood-ratio P-value thresholds for point annotations.
# These P values are from the apparent Cox models, even when C_INDEX_TYPE="cv".
P_CUT_1 <- 0.05
P_CUT_2 <- 0.01
P_CUT_3 <- 0.001

SHOW_CINDEX_LABELS <- TRUE
SHOW_INCREMENT_LABELS <- TRUE
SHOW_CLOCK_LABELS <- TRUE

PDF_WIDTH <- 14
PDF_HEIGHT <- 7.5
PNG_DPI <- 320

# -----------------------------------------------------------------------------
# Packages
# -----------------------------------------------------------------------------

required_packages <- c(
  "readr",
  "dplyr",
  "stringr",
  "ggplot2",
  "scales",
  "tibble"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Please install these R packages first: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(tibble)
})

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

normalize_icd <- function(x) {
  x <- toupper(trimws(as.character(x)))
  gsub("\\.", "", x)
}

to_numeric_safe <- function(x) {
  suppressWarnings(as.numeric(x))
}

to_logical_safe <- function(x) {
  x_chr <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    x_chr %in% c("true", "t", "1", "yes", "y") ~ TRUE,
    x_chr %in% c("false", "f", "0", "no", "n") ~ FALSE,
    TRUE ~ NA
  )
}

format_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 1e-300 ~ "P < 1e-300",
    p < 0.001 ~ paste0("P = ", format(p, scientific = TRUE, digits = 2)),
    TRUE ~ paste0("P = ", format(round(p, 3), nsmall = 3))
  )
}

significance_symbol <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < P_CUT_3 ~ "***",
    p < P_CUT_2 ~ "**",
    p < P_CUT_1 ~ "*",
    TRUE ~ "ns"
  )
}

infer_modality <- function(added_clock, added_modality) {
  out <- as.character(added_modality)

  missing_modality <- is.na(out) | out == "" | out == "NA" | out == "BASE"

  out[missing_modality & grepl("proteomics", added_clock, ignore.case = TRUE)] <- "Proteomics"
  out[missing_modality & grepl("metabolomics", added_clock, ignore.case = TRUE)] <- "Metabolomics"
  out[missing_modality & added_clock == "BASE"] <- "Baseline"
  out[is.na(out) | out == "" | out == "NA"] <- "Unknown"

  out
}

pretty_clock_label <- function(x) {
  x <- as.character(x)
  x <- gsub("_proteomics$", "\nProt", x, ignore.case = TRUE)
  x <- gsub("_metabolomics$", "\nMet", x, ignore.case = TRUE)
  x <- gsub("_", " ", x)
  x
}

safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9_-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

# -----------------------------------------------------------------------------
# Validate C-index type
# -----------------------------------------------------------------------------

if (!C_INDEX_TYPE %in% c("cv", "apparent")) {
  stop(
    "C_INDEX_TYPE must be one of: cv, apparent. Current value: ",
    C_INDEX_TYPE,
    call. = FALSE
  )
}

metric_label <- dplyr::case_when(
  C_INDEX_TYPE == "cv" ~ "5-fold cross-validated Harrell C-index",
  C_INDEX_TYPE == "apparent" ~ "Apparent Harrell C-index",
  TRUE ~ "Harrell C-index"
)

metric_tag <- safe_filename(C_INDEX_TYPE)

# -----------------------------------------------------------------------------
# Read and validate input
# -----------------------------------------------------------------------------

if (!file.exists(INPUT_FILE)) {
  stop("Input file does not exist: ", INPUT_FILE, call. = FALSE)
}

dat <- readr::read_tsv(
  INPUT_FILE,
  show_col_types = FALSE,
  progress = FALSE,
  na = c("", "NA", "NaN", "nan")
)

required_columns <- c(
  "disease_id",
  "N",
  "N_case",
  "N_noncase",
  "event_rate",
  "median_followup_years",
  "cumulative_step",
  "added_pair_id",
  "added_clock",
  "added_organ",
  "added_modality",
  "clocks_in_model",
  "n_clocks",
  "n_model_parameters",
  "events_per_parameter",
  "low_EPV_flag",
  "rank_significant",
  "individual_clock_hr",
  "individual_clock_p",
  "added_clock_hr_in_cumulative_model",
  "added_clock_ci_lo",
  "added_clock_ci_hi",
  "added_clock_p",
  "c_index",
  "base_c_index",
  "delta_c_index_vs_base",
  "delta_c_index_vs_previous",
  "lr_chi2_vs_base",
  "lr_p_vs_base",
  "sequential_lr_chi2_vs_previous",
  "sequential_lr_p_vs_previous",
  "status"
)

if (C_INDEX_TYPE == "cv") {
  required_columns <- c(
    required_columns,
    "cv_folds",
    "cv_seed",
    "cv_c_index",
    "cv_base_c_index",
    "delta_cv_c_index_vs_base",
    "delta_cv_c_index_vs_previous",
    "cv_status"
  )
}

missing_columns <- setdiff(required_columns, names(dat))
if (length(missing_columns) > 0L) {
  stop(
    "Input file is missing required columns for C_INDEX_TYPE='",
    C_INDEX_TYPE,
    "': ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

disease_code_clean <- normalize_icd(DISEASE_CODE)

disease_dat_all <- dat %>%
  mutate(
    disease_id_clean = normalize_icd(disease_id),
    cumulative_step = as.integer(cumulative_step)
  ) %>%
  filter(disease_id_clean == disease_code_clean) %>%
  arrange(cumulative_step)

if (nrow(disease_dat_all) == 0L) {
  stop(
    "Disease code ",
    disease_code_clean,
    " was not found in the cumulative-results file.",
    call. = FALSE
  )
}

if (anyDuplicated(disease_dat_all$cumulative_step)) {
  stop(
    "More than one row exists for at least one cumulative step for ",
    disease_code_clean,
    ".",
    call. = FALSE
  )
}

# -----------------------------------------------------------------------------
# Clean and prepare plotting data
# -----------------------------------------------------------------------------

if (C_INDEX_TYPE == "apparent") {
  failed_models <- disease_dat_all %>%
    filter(
      tolower(as.character(status)) != "ok" |
        is.na(c_index) |
        !is.finite(to_numeric_safe(c_index))
    )

  plot_dat <- disease_dat_all %>%
    filter(
      tolower(as.character(status)) == "ok",
      is.finite(to_numeric_safe(c_index))
    ) %>%
    mutate(
      c_index_plot = to_numeric_safe(c_index),
      base_c_index_plot = to_numeric_safe(base_c_index),
      delta_vs_base_plot = to_numeric_safe(delta_c_index_vs_base),
      delta_vs_previous_plot = to_numeric_safe(delta_c_index_vs_previous),
      plot_status = as.character(status)
    )
} else {
  failed_models <- disease_dat_all %>%
    mutate(
      cv_status_tmp = as.character(cv_status),
      cv_c_index_tmp = to_numeric_safe(cv_c_index)
    ) %>%
    filter(
      tolower(cv_status_tmp) != "ok" |
        is.na(cv_c_index_tmp) |
        !is.finite(cv_c_index_tmp)
    )

  plot_dat <- disease_dat_all %>%
    mutate(
      cv_status_tmp = as.character(cv_status),
      cv_c_index_tmp = to_numeric_safe(cv_c_index)
    ) %>%
    filter(
      tolower(cv_status_tmp) == "ok",
      is.finite(cv_c_index_tmp)
    ) %>%
    mutate(
      c_index_plot = to_numeric_safe(cv_c_index),
      base_c_index_plot = to_numeric_safe(cv_base_c_index),
      delta_vs_base_plot = to_numeric_safe(delta_cv_c_index_vs_base),
      delta_vs_previous_plot = to_numeric_safe(delta_cv_c_index_vs_previous),
      plot_status = as.character(cv_status)
    )
}

if (nrow(failed_models) > 0L) {
  warning(
    "Omitting ",
    nrow(failed_models),
    " failed or unavailable cumulative model(s) for C_INDEX_TYPE='",
    C_INDEX_TYPE,
    "': ",
    paste(
      paste0(
        "step ",
        failed_models$cumulative_step,
        " (",
        failed_models$added_clock,
        ")"
      ),
      collapse = ", "
    ),
    call. = FALSE
  )
}

plot_dat <- plot_dat %>%
  mutate(
    is_baseline = cumulative_step == 0,
    added_clock = as.character(added_clock),
    added_pair_id = as.character(added_pair_id),
    added_modality = infer_modality(added_clock, added_modality),
    added_modality = if_else(is_baseline, "Baseline", added_modality),
    rank_significant_bool = to_logical_safe(rank_significant),
    low_EPV_bool = to_logical_safe(low_EPV_flag),

    apparent_c_index = to_numeric_safe(c_index),
    apparent_base_c_index = to_numeric_safe(base_c_index),
    apparent_delta_vs_base = to_numeric_safe(delta_c_index_vs_base),
    apparent_delta_vs_previous = to_numeric_safe(delta_c_index_vs_previous),

    cv_c_index_numeric = if ("cv_c_index" %in% names(.)) to_numeric_safe(cv_c_index) else NA_real_,
    cv_base_c_index_numeric = if ("cv_base_c_index" %in% names(.)) to_numeric_safe(cv_base_c_index) else NA_real_,
    cv_delta_vs_base_numeric = if ("delta_cv_c_index_vs_base" %in% names(.)) to_numeric_safe(delta_cv_c_index_vs_base) else NA_real_,
    cv_delta_vs_previous_numeric = if ("delta_cv_c_index_vs_previous" %in% names(.)) to_numeric_safe(delta_cv_c_index_vs_previous) else NA_real_,

    sequential_p = to_numeric_safe(sequential_lr_p_vs_previous),
    sequential_chi2 = to_numeric_safe(sequential_lr_chi2_vs_previous),

    model_label = if_else(
      is_baseline,
      "Clinical\nbaseline",
      paste0("+", pretty_clock_label(added_clock))
    ),
    model_label = factor(model_label, levels = model_label),

    cindex_label = sprintf("%.3f", c_index_plot),
    increment_label = if_else(
      is_baseline | is.na(delta_vs_previous_plot),
      "",
      sprintf("%+.4f", delta_vs_previous_plot)
    ),
    delta_base_label = sprintf("%+.4f", delta_vs_base_plot),
    sequential_sig = significance_symbol(sequential_p),
    sequential_p_label = format_p(sequential_p),

    point_type = case_when(
      is_baseline ~ "Baseline",
      low_EPV_bool %in% TRUE ~ "Low EPV",
      TRUE ~ "Adequate EPV"
    )
  )

if (nrow(plot_dat) < 2L) {
  stop(
    "Fewer than two successful cumulative models are available for ",
    disease_code_clean,
    " using C_INDEX_TYPE='",
    C_INDEX_TYPE,
    "'.",
    call. = FALSE
  )
}

baseline_c <- plot_dat %>%
  filter(is_baseline) %>%
  pull(c_index_plot)

if (length(baseline_c) != 1L) {
  stop("Exactly one successful baseline row is required.", call. = FALSE)
}

best_row <- plot_dat %>%
  filter(c_index_plot == max(c_index_plot, na.rm = TRUE)) %>%
  slice(1)

final_row <- plot_dat %>%
  arrange(desc(cumulative_step)) %>%
  slice(1)

# -----------------------------------------------------------------------------
# Stable color palette by modality
# -----------------------------------------------------------------------------

modality_palette <- c(
  "Baseline" = "#444444",
  "Proteomics" = "#2F5D8C",
  "Metabolomics" = "#C98245",
  "Unknown" = "#777777"
)

available_modalities <- unique(plot_dat$added_modality)
missing_palette <- setdiff(available_modalities, names(modality_palette))
if (length(missing_palette) > 0L) {
  extra_cols <- rep("#777777", length(missing_palette))
  names(extra_cols) <- missing_palette
  modality_palette <- c(modality_palette, extra_cols)
}

# -----------------------------------------------------------------------------
# Useful numbers for subtitle and caption
# -----------------------------------------------------------------------------

N_value <- plot_dat$N[1]
N_case_value <- plot_dat$N_case[1]
N_noncase_value <- plot_dat$N_noncase[1]
event_rate_value <- plot_dat$event_rate[1]
median_followup_value <- plot_dat$median_followup_years[1]

rank_p_threshold_value <- if ("rank_p_threshold" %in% names(plot_dat)) {
  plot_dat$rank_p_threshold[1]
} else {
  NA
}

analysis_clock_set_value <- if ("analysis_clock_set" %in% names(plot_dat)) {
  plot_dat$analysis_clock_set[1]
} else {
  NA
}

complete_case_mode_value <- if ("complete_case_mode" %in% names(plot_dat)) {
  plot_dat$complete_case_mode[1]
} else {
  NA
}

cv_folds_value <- if ("cv_folds" %in% names(plot_dat)) {
  suppressWarnings(unique(plot_dat$cv_folds[!is.na(plot_dat$cv_folds)]))
} else {
  NA
}
cv_folds_value <- if (length(cv_folds_value) > 0L) cv_folds_value[1] else NA

low_epv_steps <- plot_dat %>%
  filter(!is_baseline, low_EPV_bool %in% TRUE) %>%
  pull(cumulative_step)

low_epv_caption <- if (length(low_epv_steps) > 0L) {
  paste0(
    " Low EPV warning at steps: ",
    paste(low_epv_steps, collapse = ", "),
    "."
  )
} else {
  ""
}

metric_caption <- if (C_INDEX_TYPE == "cv") {
  paste0(
    "Y-axis shows ",
    metric_label,
    ifelse(!is.na(cv_folds_value), paste0(" (", cv_folds_value, " folds)"), ""),
    ". Step-wise significance symbols are based on sequential likelihood-ratio tests from the corresponding apparent Cox models."
  )
} else {
  paste0(
    "Y-axis shows ",
    metric_label,
    ". Sequential LRT significance: * P<0.05, ** P<0.01, *** P<0.001."
  )
}

# Dynamic y-axis limits with room for labels.
c_min <- min(plot_dat$c_index_plot, na.rm = TRUE)
c_max <- max(plot_dat$c_index_plot, na.rm = TRUE)
c_range <- max(c_max - c_min, 0.01)

y_lower <- max(0, c_min - 0.30 * c_range)
y_upper <- min(1, c_max + 0.42 * c_range)

# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------

p <- ggplot(
  plot_dat,
  aes(
    x = model_label,
    y = c_index_plot,
    group = 1
  )
) +
  geom_ribbon(
    aes(
      ymin = baseline_c,
      ymax = c_index_plot,
      group = 1
    ),
    fill = "#DCE7EF",
    alpha = 0.45
  ) +
  geom_hline(
    yintercept = baseline_c,
    linetype = "dashed",
    linewidth = 0.65,
    color = "#555555"
  ) +
  geom_line(
    linewidth = 0.9,
    color = "#555555"
  ) +
  geom_point(
    aes(
      color = added_modality,
      shape = point_type
    ),
    size = 4.0,
    stroke = 0.8
  ) +
  scale_color_manual(
    values = modality_palette,
    breaks = intersect(
      c("Baseline", "Proteomics", "Metabolomics", "Unknown"),
      available_modalities
    ),
    name = "Added clock modality"
  ) +
  scale_shape_manual(
    values = c(
      "Baseline" = 15,
      "Adequate EPV" = 16,
      "Low EPV" = 17
    ),
    name = "Model status"
  ) +
  scale_y_continuous(
    limits = c(y_lower, y_upper),
    breaks = scales::pretty_breaks(n = 6),
    labels = scales::label_number(accuracy = 0.001),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Sequential cumulative model",
    y = metric_label,
    title = paste0(
      DISEASE_NAME,
      " (",
      disease_code_clean,
      "): cumulative prediction by omics mortality EPOCH clocks"
    ),
    subtitle = paste0(
      "Clinical baseline = age, sex, smoking, BMI, diastolic BP, systolic BP; ",
      "clocks are added by disease-specific HR among significant signals"
    ),
    caption = paste0(
      "N = ",
      scales::comma(N_value),
      "; cases/noncases = ",
      scales::comma(N_case_value),
      "/",
      scales::comma(N_noncase_value),
      "; event rate = ",
      scales::percent(event_rate_value, accuracy = 0.1),
      "; median follow-up = ",
      sprintf("%.2f", median_followup_value),
      " years. Baseline C = ",
      sprintf("%.3f", baseline_c),
      "; maximum C = ",
      sprintf("%.3f", best_row$c_index_plot),
      " at step ",
      best_row$cumulative_step,
      " (ΔC = ",
      sprintf("%+.4f", best_row$c_index_plot - baseline_c),
      "). ",
      metric_caption,
      low_epv_caption
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 9.5,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 10,
      color = "black"
    ),
    axis.title = element_text(
      size = 12,
      face = "bold",
      color = "black"
    ),
    plot.title = element_text(
      size = 15,
      face = "bold",
      color = "black"
    ),
    plot.subtitle = element_text(
      size = 10.5,
      color = "black"
    ),
    plot.caption = element_text(
      size = 8.2,
      hjust = 0,
      color = "black",
      margin = margin(t = 10)
    ),
    panel.grid.major.y = element_line(
      color = "grey88",
      linewidth = 0.35
    ),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", color = "black"),
    legend.text = element_text(color = "black"),
    plot.margin = margin(12, 20, 12, 12)
  )

if (SHOW_CINDEX_LABELS) {
  p <- p +
    geom_text(
      aes(label = cindex_label),
      vjust = -1.15,
      size = 3.2,
      color = "black"
    )
}

if (SHOW_INCREMENT_LABELS) {
  p <- p +
    geom_text(
      aes(label = increment_label),
      vjust = 1.55,
      size = 2.8,
      color = "grey30"
    )
}

if (SHOW_CLOCK_LABELS) {
  p <- p +
    geom_text(
      data = plot_dat %>% filter(!is_baseline),
      aes(
        label = sequential_sig,
        y = c_index_plot + 0.16 * c_range
      ),
      size = 3.3,
      fontface = "bold",
      color = "black",
      inherit.aes = TRUE
    )
}

p <- p +
  annotate(
    "text",
    x = 1,
    y = baseline_c,
    label = paste0("  Baseline C = ", sprintf("%.3f", baseline_c)),
    hjust = 0,
    vjust = -0.75,
    size = 3.2,
    color = "#444444"
  ) +
  annotate(
    "text",
    x = best_row$cumulative_step + 1,
    y = best_row$c_index_plot,
    label = paste0(
      "  Best C = ",
      sprintf("%.3f", best_row$c_index_plot),
      "\n  ΔC = ",
      sprintf("%+.4f", best_row$c_index_plot - baseline_c)
    ),
    hjust = 0,
    vjust = -0.15,
    size = 3.0,
    color = "#222222"
  )

# -----------------------------------------------------------------------------
# Export
# -----------------------------------------------------------------------------

safe_disease <- safe_filename(paste0(disease_code_clean, "_", DISEASE_NAME))

pdf_path <- file.path(
  OUTPUT_DIR,
  paste0("Cumulative_Cindex_EPOCH_PM_", safe_disease, "_", metric_tag, ".pdf")
)

png_path <- file.path(
  OUTPUT_DIR,
  paste0("Cumulative_Cindex_EPOCH_PM_", safe_disease, "_", metric_tag, ".png")
)

plot_data_path <- file.path(
  OUTPUT_DIR,
  paste0("Cumulative_Cindex_EPOCH_PM_", safe_disease, "_", metric_tag, "_plot_data.tsv")
)

summary_path <- file.path(
  OUTPUT_DIR,
  paste0("Cumulative_Cindex_EPOCH_PM_", safe_disease, "_", metric_tag, "_summary.txt")
)

ggsave(
  filename = pdf_path,
  plot = p,
  width = PDF_WIDTH,
  height = PDF_HEIGHT,
  units = "in",
  device = cairo_pdf,
  limitsize = FALSE
)

ggsave(
  filename = png_path,
  plot = p,
  width = PDF_WIDTH,
  height = PDF_HEIGHT,
  units = "in",
  dpi = PNG_DPI,
  bg = "white",
  limitsize = FALSE
)

plot_export <- plot_dat %>%
  transmute(
    disease_id = disease_code_clean,
    disease_name = DISEASE_NAME,
    c_index_type = C_INDEX_TYPE,
    cumulative_step,
    added_clock,
    added_organ,
    added_modality,
    n_clocks,
    clocks_in_model,
    N,
    N_case,
    N_noncase,
    event_rate,
    median_followup_years,

    plotted_c_index = c_index_plot,
    plotted_base_c_index = base_c_index_plot,
    plotted_delta_c_index_vs_base = delta_vs_base_plot,
    plotted_delta_c_index_vs_previous = delta_vs_previous_plot,

    apparent_c_index,
    apparent_base_c_index,
    apparent_delta_c_index_vs_base = apparent_delta_vs_base,
    apparent_delta_c_index_vs_previous = apparent_delta_vs_previous,

    cv_c_index = cv_c_index_numeric,
    cv_base_c_index = cv_base_c_index_numeric,
    delta_cv_c_index_vs_base = cv_delta_vs_base_numeric,
    delta_cv_c_index_vs_previous = cv_delta_vs_previous_numeric,

    sequential_lr_chi2_vs_previous = sequential_chi2,
    sequential_lr_p_vs_previous = sequential_p,
    sequential_significance = sequential_sig,
    individual_clock_hr,
    individual_clock_p,
    added_clock_hr_in_cumulative_model,
    added_clock_ci_lo,
    added_clock_ci_hi,
    added_clock_p,
    events_per_parameter,
    low_EPV_flag = low_EPV_bool,
    rank_significant = rank_significant_bool,
    apparent_status = status,
    plot_status
  )

readr::write_tsv(plot_export, plot_data_path)

summary_lines <- c(
  paste0("Input file: ", INPUT_FILE),
  paste0("Disease: ", DISEASE_NAME, " (", disease_code_clean, ")"),
  paste0("C-index type plotted: ", C_INDEX_TYPE),
  paste0("Metric label: ", metric_label),
  paste0("Analysis clock set: ", analysis_clock_set_value),
  paste0("Complete-case mode: ", complete_case_mode_value),
  paste0("Rank P threshold: ", rank_p_threshold_value),
  paste0("N: ", N_value),
  paste0("Cases: ", N_case_value),
  paste0("Noncases: ", N_noncase_value),
  paste0("Event rate: ", sprintf("%.6f", event_rate_value)),
  paste0("Median follow-up years: ", sprintf("%.6f", median_followup_value)),
  paste0("Successful models plotted: ", nrow(plot_dat)),
  paste0("Failed/unavailable models omitted: ", nrow(failed_models)),
  paste0("Baseline plotted C-index: ", sprintf("%.6f", baseline_c)),
  paste0(
    "Final plotted C-index: ",
    sprintf("%.6f", final_row$c_index_plot),
    " at step ",
    final_row$cumulative_step
  ),
  paste0(
    "Final plotted improvement over baseline: ",
    sprintf("%.6f", final_row$c_index_plot - baseline_c)
  ),
  paste0(
    "Maximum plotted C-index: ",
    sprintf("%.6f", best_row$c_index_plot),
    " at step ",
    best_row$cumulative_step,
    " after adding ",
    best_row$added_clock
  ),
  paste0(
    "Maximum plotted improvement over baseline: ",
    sprintf("%.6f", best_row$c_index_plot - baseline_c)
  ),
  paste0(
    "Low EPV steps: ",
    ifelse(length(low_epv_steps) == 0L, "None", paste(low_epv_steps, collapse = ", "))
  ),
  paste0("PDF: ", pdf_path),
  paste0("PNG: ", png_path),
  paste0("Plot data: ", plot_data_path)
)

writeLines(summary_lines, summary_path)

message("Disease: ", DISEASE_NAME, " (", disease_code_clean, ")")
message("Input: ", INPUT_FILE)
message("C-index type: ", C_INDEX_TYPE)
message("Baseline plotted C-index: ", sprintf("%.4f", baseline_c))
message(
  "Final plotted C-index: ",
  sprintf("%.4f", final_row$c_index_plot),
  " | final plotted Delta C: ",
  sprintf("%+.4f", final_row$c_index_plot - baseline_c)
)
message(
  "Best plotted C-index: ",
  sprintf("%.4f", best_row$c_index_plot),
  " at step ",
  best_row$cumulative_step,
  " | best plotted Delta C: ",
  sprintf("%+.4f", best_row$c_index_plot - baseline_c)
)
message("Wrote PDF: ", pdf_path)
message("Wrote PNG: ", png_path)
message("Wrote plotting data: ", plot_data_path)
message("Wrote summary: ", summary_path)