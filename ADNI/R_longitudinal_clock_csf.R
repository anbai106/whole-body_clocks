# ============================================================
# Plot baseline CSF biomarkers versus longitudinal AD L'EPOCH slope_per_year
#
# Corrected version:
#   The analysis dataset does NOT contain a column named change_metric.
#   Instead, slope_per_year is a direct column in the analysis dataset.
#
# Main visualization:
#   x-axis: baseline CSF biomarker
#   y-axis: longitudinal AD L'EPOCH slope_per_year
#
# Primary clock:
#   adni_brain_mri_ad_lepoch_acceleration_z
#
# Statistical annotation from combined association table:
#   slope_per_year ~ baseline_CSF
#                    + baseline_clock_value
#                    + conversion_group
#                    + Age + Sex + ICV + APOE4
#                    + followup_span_years
#
# Update:
#   For p_raw only, exclude y-axis outlier points:
#      abs(slope_per_year) < 1
#
# Outputs:
#   1. Raw scatter:
#      baseline CSF biomarker vs slope_per_year
#      with abs(slope_per_year) < 1
#
#   2. Pathology-direction scatter:
#      worse CSF pathology vs slope_per_year
#
#   3. Adjusted residual scatter:
#      residualized baseline CSF biomarker vs residualized slope_per_year
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(ggplot2)
  library(scales)
  library(grid)
})

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

base_dir <- "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/adni_lepoch/longitudinal_change_vs_baseline_csf"

analysis_file <- file.path(
  base_dir,
  "adni_brain_mri_ad_lepoch_longitudinal_change_vs_baseline_csf_analysis_dataset.tsv"
)

combined_assoc_file <- file.path(
  base_dir,
  "adni_brain_mri_ad_lepoch_longitudinal_change_vs_baseline_csf_combined_associations.tsv"
)

out_dir <- file.path(base_dir, "plots_slope_per_year_vs_baseline_csf")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 2. Main options
# ------------------------------------------------------------

primary_clock_col <- "adni_brain_mri_ad_lepoch_acceleration_z"
primary_change_metric <- "slope_per_year"

# Raw scatter y-axis outlier filter.
# Applied only to p_raw.
raw_y_abs_limit <- 1

csf_vars <- c(
  "Abeta_CSF",
  "Tau_CSF",
  "PTau_CSF"
)

csf_labels <- c(
  "Abeta_CSF" = "CSF Aβ42",
  "Tau_CSF" = "CSF total tau",
  "PTau_CSF" = "CSF p-tau"
)

csf_pathology_labels <- c(
  "Abeta_CSF" = "Worse CSF Aβ pathology (-Aβ42)",
  "Tau_CSF" = "Worse CSF tau pathology",
  "PTau_CSF" = "Worse CSF p-tau pathology"
)

group_levels <- c(
  "Non-event & censored",
  "CN-MCI",
  "CN-AD"
)

group_labels <- c(
  "Non-event & censored" = "Censored / non-event",
  "CN-MCI" = "CN to MCI",
  "CN-AD" = "CN to AD"
)

group_palette <- c(
  "Censored / non-event" = "#8C6D31",
  "CN to MCI" = "#2A6F9E",
  "CN to AD" = "#A23E48"
)

# ------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------

format_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "P = NA",
    p < 2.2e-16 ~ "P < 2.2e-16",
    p < 0.001 ~ paste0("P = ", formatC(p, format = "e", digits = 2)),
    TRUE ~ paste0("P = ", signif(p, 3))
  )
}

format_q <- function(q) {
  dplyr::case_when(
    is.na(q) ~ "FDR = NA",
    q < 2.2e-16 ~ "FDR < 2.2e-16",
    q < 0.001 ~ paste0("FDR = ", formatC(q, format = "e", digits = 2)),
    TRUE ~ paste0("FDR = ", signif(q, 3))
  )
}

format_num <- function(x, digits = 3) {
  dplyr::case_when(
    is.na(x) ~ "NA",
    abs(x) < 0.001 ~ formatC(x, format = "e", digits = 2),
    TRUE ~ sprintf(paste0("%.", digits, "f"), x)
  )
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

make_annotation_positions <- function(dat, x_col, y_col) {
  dat %>%
    group_by(csf_var, biomarker_display) %>%
    summarise(
      x_min = min(.data[[x_col]], na.rm = TRUE),
      x_max = max(.data[[x_col]], na.rm = TRUE),
      y_min = min(.data[[y_col]], na.rm = TRUE),
      y_max = max(.data[[y_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      x_range = if_else(
        is.finite(x_max - x_min) & (x_max - x_min) > 0,
        x_max - x_min,
        1
      ),
      y_range = if_else(
        is.finite(y_max - y_min) & (y_max - y_min) > 0,
        y_max - y_min,
        1
      ),
      x_label = x_min + 0.04 * x_range,
      y_label = y_max - 0.04 * y_range
    )
}

# ------------------------------------------------------------
# 4. Read data
# ------------------------------------------------------------

if (!file.exists(analysis_file)) {
  stop("Cannot find analysis file: ", analysis_file)
}

if (!file.exists(combined_assoc_file)) {
  stop("Cannot find combined association file: ", combined_assoc_file)
}

analysis_df <- readr::read_tsv(analysis_file, show_col_types = FALSE)
assoc_df <- readr::read_tsv(combined_assoc_file, show_col_types = FALSE)

message("Loaded analysis file: ", analysis_file)
message("Loaded association file: ", combined_assoc_file)
message("Analysis columns:")
print(colnames(analysis_df))

# ------------------------------------------------------------
# 5. Check required columns
# ------------------------------------------------------------

required_analysis_cols <- c(
  "PTID",
  "clock_col",
  "metric_status",
  "conversion_group_3level",
  "baseline_clock_value",
  "followup_span_years",
  "_age_model",
  "_sex_male_model",
  "_icv_model",
  "_apoe4_model",
  primary_change_metric,
  "csf_var",
  "csf_value"
)

missing_analysis_cols <- setdiff(required_analysis_cols, colnames(analysis_df))

if (length(missing_analysis_cols) > 0) {
  stop(
    "Missing required columns in analysis file: ",
    paste(missing_analysis_cols, collapse = ", "),
    "\n\nImportant note: this corrected script does NOT require change_metric in the analysis file."
  )
}

required_assoc_cols <- c(
  "analysis_type",
  "clock_col",
  "change_metric",
  "csf_var",
  "n",
  "status",
  "std_beta_predictor",
  "partial_r_predictor",
  "p_raw_predictor"
)

missing_assoc_cols <- setdiff(required_assoc_cols, colnames(assoc_df))

if (length(missing_assoc_cols) > 0) {
  stop(
    "Missing required columns in combined association table: ",
    paste(missing_assoc_cols, collapse = ", ")
  )
}

if (!"p_bh_all_tests" %in% colnames(assoc_df)) {
  assoc_df$p_bh_all_tests <- NA_real_
}

if (!"std_beta_pathology_direction" %in% colnames(assoc_df)) {
  assoc_df <- assoc_df %>%
    mutate(
      std_beta_pathology_direction = case_when(
        csf_var == "Abeta_CSF" ~ -1 * std_beta_predictor,
        TRUE ~ std_beta_predictor
      )
    )
}

if (!"partial_r_pathology_direction" %in% colnames(assoc_df)) {
  assoc_df <- assoc_df %>%
    mutate(
      partial_r_pathology_direction = case_when(
        csf_var == "Abeta_CSF" ~ -1 * partial_r_predictor,
        TRUE ~ partial_r_predictor
      )
    )
}

# ------------------------------------------------------------
# 6. Clean analysis dataset
# ------------------------------------------------------------

plot_df <- analysis_df %>%
  filter(
    clock_col == primary_clock_col,
    metric_status == "ok",
    csf_var %in% csf_vars
  ) %>%
  mutate(
    PTID = as.character(PTID),
    conversion_group_3level = factor(
      as.character(conversion_group_3level),
      levels = group_levels
    ),
    group_display = factor(
      group_labels[as.character(conversion_group_3level)],
      levels = unname(group_labels)
    ),
    biomarker_display = factor(
      csf_labels[csf_var],
      levels = unname(csf_labels)
    ),
    pathology_display = factor(
      csf_pathology_labels[csf_var],
      levels = unname(csf_pathology_labels)
    ),
    csf_value = safe_numeric(csf_value),
    slope_per_year = safe_numeric(.data[[primary_change_metric]]),
    baseline_clock_value = safe_numeric(baseline_clock_value),
    followup_span_years = safe_numeric(followup_span_years),
    age_model = safe_numeric(`_age_model`),
    sex_male_model = safe_numeric(`_sex_male_model`),
    icv_model = safe_numeric(`_icv_model`),
    apoe4_model = safe_numeric(`_apoe4_model`)
  )

if ("csf_pathology_value" %in% colnames(plot_df)) {
  plot_df <- plot_df %>%
    mutate(
      csf_pathology_value = safe_numeric(csf_pathology_value)
    )
} else {
  plot_df <- plot_df %>%
    mutate(
      csf_pathology_value = case_when(
        csf_var == "Abeta_CSF" ~ -1 * csf_value,
        TRUE ~ csf_value
      )
    )
}

plot_df <- plot_df %>%
  filter(
    is.finite(csf_value),
    is.finite(csf_pathology_value),
    is.finite(slope_per_year),
    is.finite(baseline_clock_value),
    is.finite(followup_span_years),
    is.finite(age_model),
    is.finite(sex_male_model),
    is.finite(icv_model),
    is.finite(apoe4_model),
    !is.na(conversion_group_3level)
  ) %>%
  droplevels()

if (nrow(plot_df) == 0) {
  stop(
    "No rows available for plotting after filtering.\n",
    "Check that clock_col == ", primary_clock_col,
    ", metric_status == ok, and CSF variables are present."
  )
}

# Raw scatter only: remove extreme y-axis outliers.
plot_df_raw <- plot_df %>%
  filter(abs(slope_per_year) < raw_y_abs_limit) %>%
  droplevels()

message("Rows available for plotting, full plot_df: ", nrow(plot_df))
message("Subjects available for plotting, full plot_df: ", dplyr::n_distinct(plot_df$PTID))
message("Rows available for p_raw after abs(slope_per_year) < ", raw_y_abs_limit, ": ", nrow(plot_df_raw))
message("Subjects available for p_raw after y-axis outlier filtering: ", dplyr::n_distinct(plot_df_raw$PTID))
message("Rows excluded from p_raw only: ", nrow(plot_df) - nrow(plot_df_raw))

# ------------------------------------------------------------
# 7. Clean combined association annotations
# ------------------------------------------------------------

assoc_plot <- assoc_df %>%
  filter(
    analysis_type == "combined_group_adjusted",
    clock_col == primary_clock_col,
    change_metric == primary_change_metric,
    csf_var %in% csf_vars
  ) %>%
  mutate(
    biomarker_display = factor(
      csf_labels[csf_var],
      levels = unname(csf_labels)
    ),
    pathology_display = factor(
      csf_pathology_labels[csf_var],
      levels = unname(csf_pathology_labels)
    ),
    label_raw = paste0(
      "N = ", n,
      "\nstd β = ", format_num(std_beta_predictor),
      "\npartial r = ", format_num(partial_r_predictor),
      "\n", format_p(p_raw_predictor),
      "\n", format_q(p_bh_all_tests)
    ),
    label_pathology = paste0(
      "N = ", n,
      "\nstd βpath = ", format_num(std_beta_pathology_direction),
      "\npartial rpath = ", format_num(partial_r_pathology_direction),
      "\n", format_p(p_raw_predictor),
      "\n", format_q(p_bh_all_tests)
    )
  )

if (nrow(assoc_plot) == 0) {
  stop(
    "No matching rows in association table.\n",
    "Expected rows with analysis_type == combined_group_adjusted, clock_col == ",
    primary_clock_col,
    ", and change_metric == ",
    primary_change_metric
  )
}

readr::write_tsv(
  plot_df,
  file.path(out_dir, "adni_lepoch_slope_per_year_vs_baseline_csf_plot_data.tsv")
)

readr::write_tsv(
  plot_df_raw,
  file.path(out_dir, "adni_lepoch_slope_per_year_vs_baseline_csf_raw_scatter_y_filtered_plot_data.tsv")
)

readr::write_tsv(
  assoc_plot,
  file.path(out_dir, "adni_lepoch_slope_per_year_vs_baseline_csf_annotation_table.tsv")
)

# ------------------------------------------------------------
# 8. Plot 1: Raw scatter
#    Uses plot_df_raw only:
#       abs(slope_per_year) < raw_y_abs_limit
# ------------------------------------------------------------

raw_annot_tbl <- make_annotation_positions(
  dat = plot_df_raw,
  x_col = "csf_value",
  y_col = "slope_per_year"
) %>%
  left_join(
    assoc_plot %>%
      select(csf_var, biomarker_display, label_raw),
    by = c("csf_var", "biomarker_display")
  )

p_raw <- ggplot(
  plot_df_raw,
  aes(
    x = csf_value,
    y = slope_per_year
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.35,
    color = "grey70"
  ) +
  geom_point(
    aes(color = group_display),
    alpha = 0.62,
    size = 1.9,
    stroke = 0
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 1.0,
    color = "black",
    alpha = 0.16
  ) +
  geom_label(
    data = raw_annot_tbl,
    aes(
      x = x_label,
      y = y_label,
      label = label_raw
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.0,
    label.size = 0.25,
    label.padding = unit(0.17, "lines"),
    fill = "white",
    alpha = 0.90,
    color = "black"
  ) +
  facet_wrap(
    ~ biomarker_display,
    scales = "free_x",
    nrow = 1
  ) +
  scale_color_manual(values = group_palette, name = NULL) +
  scale_x_continuous(
    name = "Baseline CSF biomarker",
    labels = number_format(accuracy = 0.1)
  ) +
  scale_y_continuous(
    name = "AD L’EPOCH acceleration_z slope per year",
    labels = number_format(accuracy = 0.001)
  ) +
  labs(
    title = "Baseline CSF biomarkers versus longitudinal AD L’EPOCH slope",
    subtitle = paste0(
      "Black line shows the unadjusted visual trend after excluding y-axis outliers with |slope per year| ≥ ",
      raw_y_abs_limit,
      "; labels report the full-sample adjusted model."
    ),
    caption = "Adjusted model: slope_per_year ~ baseline CSF + baseline clock + conversion group + age + sex + ICV + APOE4 + follow-up span."
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10.5),
    plot.caption = element_text(size = 8.5, hjust = 0),
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 12, 8, 8)
  )

ggsave(
  file.path(out_dir, "adni_lepoch_slope_per_year_vs_baseline_csf_raw_scatter.pdf"),
  p_raw,
  width = 12.5,
  height = 5.6
)

ggsave(
  file.path(out_dir, "adni_lepoch_slope_per_year_vs_baseline_csf_raw_scatter.png"),
  p_raw,
  width = 12.5,
  height = 5.6,
  dpi = 300
)

# ------------------------------------------------------------
# 9. Plot 2: Pathology-direction scatter
#    Uses full plot_df, unchanged.
# ------------------------------------------------------------

path_annot_tbl <- make_annotation_positions(
  dat = plot_df,
  x_col = "csf_pathology_value",
  y_col = "slope_per_year"
) %>%
  left_join(
    assoc_plot %>%
      select(csf_var, biomarker_display, label_pathology),
    by = c("csf_var", "biomarker_display")
  )

p_pathology <- ggplot(
  plot_df,
  aes(
    x = csf_pathology_value,
    y = slope_per_year
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.35,
    color = "grey70"
  ) +
  geom_point(
    aes(color = group_display),
    alpha = 0.62,
    size = 1.9,
    stroke = 0
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 1.0,
    color = "black",
    alpha = 0.16
  ) +
  geom_label(
    data = path_annot_tbl,
    aes(
      x = x_label,
      y = y_label,
      label = label_pathology
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.0,
    label.size = 0.25,
    label.padding = unit(0.17, "lines"),
    fill = "white",
    alpha = 0.90,
    color = "black"
  ) +
  facet_wrap(
    ~ pathology_display,
    scales = "free_x",
    nrow = 1
  ) +
  scale_color_manual(values = group_palette, name = NULL) +
  scale_x_continuous(
    name = "Baseline CSF pathology direction",
    labels = number_format(accuracy = 0.1)
  ) +
  scale_y_continuous(
    name = "AD L’EPOCH acceleration_z slope per year",
    labels = number_format(accuracy = 0.001)
  ) +
  labs(
    title = "Baseline CSF pathology versus longitudinal AD L’EPOCH slope",
    subtitle = "Higher x-axis values indicate worse CSF pathology for all three biomarkers.",
    caption = "For Aβ42, the sign is reversed so that higher values indicate worse amyloid pathology. Black line shows the unadjusted trend; labels report the adjusted model."
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10.5),
    plot.caption = element_text(size = 8.5, hjust = 0),
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 12, 8, 8)
  )

ggsave(
  file.path(out_dir, "adni_lepoch_slope_per_year_vs_baseline_csf_pathology_direction_scatter.pdf"),
  p_pathology,
  width = 12.5,
  height = 5.6
)

ggsave(
  file.path(out_dir, "adni_lepoch_slope_per_year_vs_baseline_csf_pathology_direction_scatter.png"),
  p_pathology,
  width = 12.5,
  height = 5.6,
  dpi = 300
)

# ------------------------------------------------------------
# 10. Plot 3: Adjusted residual scatter
#     Uses full plot_df, unchanged.
# ------------------------------------------------------------

make_adjusted_residuals_one <- function(dat, csf_var_i) {
  d <- dat %>%
    filter(csf_var == csf_var_i) %>%
    mutate(
      conversion_group_3level = droplevels(conversion_group_3level)
    )
  
  if (nrow(d) < 8) {
    return(tibble())
  }
  
  candidate_covariates <- c(
    "baseline_clock_value",
    "conversion_group_3level",
    "age_model",
    "sex_male_model",
    "icv_model",
    "apoe4_model",
    "followup_span_years"
  )
  
  keep_covariates <- candidate_covariates[
    map_lgl(
      d[candidate_covariates],
      ~ dplyr::n_distinct(.x, na.rm = TRUE) > 1
    )
  ]
  
  rhs <- paste(keep_covariates, collapse = " + ")
  
  y_formula <- as.formula(paste0("slope_per_year ~ ", rhs))
  x_formula <- as.formula(paste0("csf_value ~ ", rhs))
  
  fit_y <- lm(y_formula, data = d)
  fit_x <- lm(x_formula, data = d)
  
  d %>%
    mutate(
      slope_per_year_adjusted_residual = residuals(fit_y),
      csf_value_adjusted_residual = residuals(fit_x),
      residualization_covariates = rhs
    )
}

residual_df <- map_dfr(
  csf_vars,
  ~ make_adjusted_residuals_one(plot_df, .x)
)

readr::write_tsv(
  residual_df,
  file.path(out_dir, "adni_lepoch_slope_per_year_vs_baseline_csf_adjusted_residual_plot_data.tsv")
)

resid_annot_tbl <- make_annotation_positions(
  dat = residual_df,
  x_col = "csf_value_adjusted_residual",
  y_col = "slope_per_year_adjusted_residual"
) %>%
  left_join(
    assoc_plot %>%
      select(csf_var, biomarker_display, label_raw),
    by = c("csf_var", "biomarker_display")
  )

p_residual <- ggplot(
  residual_df,
  aes(
    x = csf_value_adjusted_residual,
    y = slope_per_year_adjusted_residual
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.35,
    color = "grey70"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.35,
    color = "grey70"
  ) +
  geom_point(
    aes(color = group_display),
    alpha = 0.62,
    size = 1.9,
    stroke = 0
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 1.0,
    color = "black",
    alpha = 0.16
  ) +
  geom_label(
    data = resid_annot_tbl,
    aes(
      x = x_label,
      y = y_label,
      label = label_raw
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.0,
    label.size = 0.25,
    label.padding = unit(0.17, "lines"),
    fill = "white",
    alpha = 0.90,
    color = "black"
  ) +
  facet_wrap(
    ~ biomarker_display,
    scales = "free",
    nrow = 1
  ) +
  scale_color_manual(values = group_palette, name = NULL) +
  scale_x_continuous(
    name = "Adjusted residual of baseline CSF biomarker",
    labels = number_format(accuracy = 0.1)
  ) +
  scale_y_continuous(
    name = "Adjusted residual of AD L’EPOCH slope per year",
    labels = number_format(accuracy = 0.001)
  ) +
  labs(
    title = "Adjusted partial association between baseline CSF biomarkers and AD L’EPOCH slope",
    subtitle = "Both axes are residualized for baseline clock, conversion group, age, sex, ICV, APOE4, and follow-up span.",
    caption = "Black line shows the full-sample partial association. Points are colored by future conversion group."
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10.5),
    plot.caption = element_text(size = 8.5, hjust = 0),
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 12, 8, 8)
  )

ggsave(
  file.path(out_dir, "adni_lepoch_slope_per_year_vs_baseline_csf_adjusted_residual_scatter.pdf"),
  p_residual,
  width = 12.5,
  height = 5.6
)

ggsave(
  file.path(out_dir, "adni_lepoch_slope_per_year_vs_baseline_csf_adjusted_residual_scatter.png"),
  p_residual,
  width = 12.5,
  height = 5.6,
  dpi = 300
)

# ------------------------------------------------------------
# 11. Print summary
# ------------------------------------------------------------

message("============================================================")
message("Done.")
message("Plots saved to: ", out_dir)
message("Primary plotted clock column: ", primary_clock_col)
message("Primary plotted longitudinal metric: ", primary_change_metric)
message("Raw scatter y-axis filter: abs(slope_per_year) < ", raw_y_abs_limit)
message("============================================================")

message("Combined model annotations:")
print(
  assoc_plot %>%
    select(
      csf_var,
      n,
      std_beta_predictor,
      partial_r_predictor,
      std_beta_pathology_direction,
      partial_r_pathology_direction,
      p_raw_predictor,
      p_bh_all_tests
    )
)

message("Plot data sample-size summary, full plot_df:")
print(
  plot_df %>%
    group_by(csf_var, biomarker_display, group_display) %>%
    summarise(
      n = n(),
      n_subjects = n_distinct(PTID),
      mean_slope_per_year = mean(slope_per_year, na.rm = TRUE),
      sd_slope_per_year = sd(slope_per_year, na.rm = TRUE),
      .groups = "drop"
    )
)

message("Plot data sample-size summary, p_raw only after y-axis filtering:")
print(
  plot_df_raw %>%
    group_by(csf_var, biomarker_display, group_display) %>%
    summarise(
      n = n(),
      n_subjects = n_distinct(PTID),
      mean_slope_per_year = mean(slope_per_year, na.rm = TRUE),
      sd_slope_per_year = sd(slope_per_year, na.rm = TRUE),
      .groups = "drop"
    )
)