#!/usr/bin/env Rscript

# ==============================================================================
# Publication-ready 4-panel figure (E-F-G-H)
# Color-harmonized to the survival-analysis A-B-C-D figure
# Brain MRI mortality EPOCH vs 9 AI-derived disease subtype scores
# Fixed 3-year logistic-regression sensitivity analysis
#
# This figure is designed to echo the survival-analysis A-B-C-D panels:
#
#   E = adjusted odds-ratio forest plot
#       (echoes survival Panel A: adjusted HR forest plot)
#
#   F = incremental repeated out-of-fold ROC-AUC vs covariates only
#       (echoes survival Panel B: incremental C-index)
#
#   G = paired-bootstrap ROC-AUC difference, EPOCH minus each subtype
#       (echoes survival Panel C: paired-bootstrap C-index difference)
#
#   H = covariates vs EPOCH vs all 9 subtypes vs all 9 + EPOCH
#       (echoes survival Panel D: combined-panel comparison)
#
# Fixed-horizon endpoint:
#   Case    = incident inpatient ICD-10 F* or G* diagnosis within 3 years after MRI
#   Control = F/G-disease-free and observed through the full 3-year horizon
#
# Inputs expected in fixed_3y/:
#   single_marker_adjusted.tsv
#   cv_model_metrics.tsv
#   bootstrap_discrimination_comparisons.tsv
#   combined_adjusted_models.tsv
#
# Outputs:
#   Individual panel PDF/PNG files
#   Combined E-F-G-H PDF/PNG
#   Compact plotting summary TSV
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. User settings
# ------------------------------------------------------------------------------

input_dir_candidates <- c(
  paste0(
    "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/",
    "brain_mri_mortality_clock/",
    "comparative_with_9_disease_subtypes_any_FG/",
    "logistic_regression/fixed_3y"
  ),
  paste0(
    "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/",
    "brain_mri_mortality_clock/",
    "comparative_with_9_disease_subtypes_any_FG/",
    "logistic_regression/fixed_3y"
  )
)

existing_input_dirs <- input_dir_candidates[
  dir.exists(input_dir_candidates)
]

if (length(existing_input_dirs) == 0) {
  stop(
    "Could not find the fixed_3y result directory. Checked:\n",
    paste(input_dir_candidates, collapse = "\n")
  )
}

input_dir <- existing_input_dirs[[1]]

output_dir <- file.path(
  input_dir,
  "plots"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

prefix <- "brain_mri_mortality_EPOCH_vs_9_AI_subtypes_fixed3y_logistic"

single_width <- 8.3
single_height <- 5.8

combined_width <- 13.4
combined_height <- 10.5


# ------------------------------------------------------------------------------
# 2. Packages
# ------------------------------------------------------------------------------

required_packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "ggplot2",
  "forcats",
  "scales",
  "stringr",
  "patchwork"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall with:\ninstall.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))"
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(scales)
  library(stringr)
  library(patchwork)
})


# ------------------------------------------------------------------------------
# 3. Read inputs
# ------------------------------------------------------------------------------

single_marker_file <- file.path(
  input_dir,
  "single_marker_adjusted.tsv"
)

cv_file <- file.path(
  input_dir,
  "cv_model_metrics.tsv"
)

bootstrap_file <- file.path(
  input_dir,
  "bootstrap_discrimination_comparisons.tsv"
)

combined_file <- file.path(
  input_dir,
  "combined_adjusted_models.tsv"
)

required_files <- c(
  single_marker_file,
  cv_file,
  bootstrap_file,
  combined_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    "Missing required input files:\n",
    paste(missing_files, collapse = "\n")
  )
}

single_tbl <- readr::read_tsv(
  single_marker_file,
  show_col_types = FALSE,
  progress = FALSE
)

cv_tbl <- readr::read_tsv(
  cv_file,
  show_col_types = FALSE,
  progress = FALSE
)

bootstrap_tbl <- readr::read_tsv(
  bootstrap_file,
  show_col_types = FALSE,
  progress = FALSE
)

combined_tbl <- readr::read_tsv(
  combined_file,
  show_col_types = FALSE,
  progress = FALSE
)


# ------------------------------------------------------------------------------
# 4. Validate inputs
# ------------------------------------------------------------------------------

check_columns <- function(df, required_cols, table_name) {
  missing_cols <- setdiff(
    required_cols,
    names(df)
  )
  
  if (length(missing_cols) > 0) {
    stop(
      table_name,
      " is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
}

check_columns(
  single_tbl,
  c(
    "marker",
    "N",
    "N_case",
    "N_control",
    "OR_per_1SD",
    "OR_CI_low",
    "OR_CI_high",
    "p_marker",
    "lrt_p_vs_base"
  ),
  "single_marker_adjusted.tsv"
)

check_columns(
  cv_tbl,
  c(
    "model",
    "N",
    "N_case",
    "N_control",
    "ROC_AUC_repeated_OOF",
    "AUPRC_repeated_OOF",
    "Brier_repeated_OOF",
    "log_loss_repeated_OOF",
    "cv_folds",
    "cv_repeats"
  ),
  "cv_model_metrics.tsv"
)

check_columns(
  bootstrap_tbl,
  c(
    "comparison",
    "context",
    "model_A",
    "model_B",
    "delta_ROC_AUC",
    "ROC_AUC_CI_low",
    "ROC_AUC_CI_high",
    "ROC_AUC_bootstrap_p",
    "delta_AUPRC",
    "AUPRC_CI_low",
    "AUPRC_CI_high",
    "AUPRC_bootstrap_p"
  ),
  "bootstrap_discrimination_comparisons.tsv"
)

check_columns(
  combined_tbl,
  c(
    "N",
    "N_case",
    "N_control",
    "lrt_p_all9_vs_covariates",
    "lrt_p_EPOCH_beyond_all9",
    "lrt_p_all9_beyond_EPOCH",
    "EPOCH_conditional_OR_given_all9",
    "EPOCH_conditional_CI_low",
    "EPOCH_conditional_CI_high",
    "EPOCH_conditional_p_given_all9"
  ),
  "combined_adjusted_models.tsv"
)

if (nrow(combined_tbl) != 1) {
  stop(
    "combined_adjusted_models.tsv should contain exactly one summary row."
  )
}


# ------------------------------------------------------------------------------
# 5. Labels, ordering, colors, helpers
# ------------------------------------------------------------------------------

marker_order <- c(
  "Brain_MRI_mortality_EPOCH",
  "AD1",
  "AD2",
  "ASD1",
  "ASD2",
  "ASD3",
  "LLD1",
  "LLD2",
  "SCZ1",
  "SCZ2"
)

pretty_marker <- c(
  "Brain_MRI_mortality_EPOCH" = "Mortality EPOCH",
  "AD1" = "AD1",
  "AD2" = "AD2",
  "ASD1" = "ASD1",
  "ASD2" = "ASD2",
  "ASD3" = "ASD3",
  "LLD1" = "LLD1",
  "LLD2" = "LLD2",
  "SCZ1" = "SCZ1",
  "SCZ2" = "SCZ2"
)

adjusted_model_order <- c(
  "EPOCH_adjusted",
  "AD1_adjusted",
  "AD2_adjusted",
  "ASD1_adjusted",
  "ASD2_adjusted",
  "ASD3_adjusted",
  "LLD1_adjusted",
  "LLD2_adjusted",
  "SCZ1_adjusted",
  "SCZ2_adjusted"
)

adjusted_model_to_marker <- c(
  "EPOCH_adjusted" = "Brain_MRI_mortality_EPOCH",
  "AD1_adjusted" = "AD1",
  "AD2_adjusted" = "AD2",
  "ASD1_adjusted" = "ASD1",
  "ASD2_adjusted" = "ASD2",
  "ASD3_adjusted" = "ASD3",
  "LLD1_adjusted" = "LLD1",
  "LLD2_adjusted" = "LLD2",
  "SCZ1_adjusted" = "SCZ1",
  "SCZ2_adjusted" = "SCZ2"
)

# ------------------------------------------------------------------------------
# Color palette copied from the survival A-B-C-D figure so that E-F-G-H uses
# the same visual language throughout the manuscript.
# ------------------------------------------------------------------------------

col_epoch     <- "#B55239"  # mortality EPOCH: warm rust/red
col_subtype   <- "#5E6D7A"  # AI disease subtypes: blue-grey
col_reference <- "#8A8A8A"  # covariate-only reference: neutral grey
col_combined  <- "#6A5A89"  # all 9 subtypes + EPOCH: muted purple

method_palette <- c(
  "Mortality EPOCH" = col_epoch,
  "AI disease subtype" = col_subtype
)

combined_palette <- c(
  "Reference" = col_reference,
  "Mortality EPOCH" = col_epoch,
  "AI disease subtype" = col_subtype,
  "Combined" = col_combined
)

# Panel F mirrors survival Panel B:
# EPOCH = filled circle; subtype markers = open circles.
method_shape <- c(
  "Mortality EPOCH" = 16,
  "AI disease subtype" = 1
)

significance_shape <- c(
  "P < 0.05" = 16,
  "P >= 0.05" = 1
)

format_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "P = NA",
    p < 0.001 ~ paste0(
      "P = ",
      formatC(
        p,
        format = "e",
        digits = 2
      )
    ),
    TRUE ~ paste0(
      "P = ",
      formatC(
        p,
        format = "f",
        digits = 3
      )
    )
  )
}

format_delta <- function(x, digits = 3) {
  sprintf(
    paste0("%+.", digits, "f"),
    x
  )
}

theme_epoch <- function(base_size = 11.5) {
  theme_classic(
    base_size = base_size
  ) +
    theme(
      axis.text = element_text(
        color = "black"
      ),
      axis.title = element_text(
        face = "bold",
        color = "black"
      ),
      plot.title = element_text(
        face = "bold",
        size = base_size + 1.5
      ),
      plot.subtitle = element_text(
        size = base_size - 0.6
      ),
      plot.caption = element_text(
        size = base_size - 2.0,
        hjust = 0,
        lineheight = 1.0
      ),
      panel.grid.major.y = element_line(
        color = "grey92",
        linewidth = 0.35
      ),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.margin = margin(
        8,
        12,
        10,
        8
      )
    )
}


# ------------------------------------------------------------------------------
# 6. Common sample information
# ------------------------------------------------------------------------------

common_n <- single_tbl$N[[1]]
common_cases <- single_tbl$N_case[[1]]
common_controls <- single_tbl$N_control[[1]]

base_auc_row <- cv_tbl |>
  filter(
    model == "Covariates_only"
  )

if (nrow(base_auc_row) != 1) {
  stop(
    "Expected exactly one Covariates_only row in cv_model_metrics.tsv."
  )
}

base_auc <- base_auc_row$ROC_AUC_repeated_OOF[[1]]

cv_folds <- base_auc_row$cv_folds[[1]]
cv_repeats <- base_auc_row$cv_repeats[[1]]


# ------------------------------------------------------------------------------
# 7. Panel E: adjusted OR forest plot
#    Echoes survival Panel A
# ------------------------------------------------------------------------------

panel_e_df <- single_tbl |>
  filter(
    marker %in% marker_order
  ) |>
  mutate(
    marker = factor(
      marker,
      levels = marker_order
    ),
    label = unname(
      pretty_marker[
        as.character(marker)
      ]
    ),
    method_group = if_else(
      as.character(marker) == "Brain_MRI_mortality_EPOCH",
      "Mortality EPOCH",
      "AI disease subtype"
    ),
    significance = if_else(
      p_marker < 0.05,
      "P < 0.05",
      "P >= 0.05"
    )
  ) |>
  arrange(
    marker
  ) |>
  mutate(
    label = factor(
      label,
      levels = rev(
        unname(
          pretty_marker[
            marker_order
          ]
        )
      )
    )
  )

if (nrow(panel_e_df) != 10) {
  warning(
    "Expected 10 biomarker rows for Panel E; found ",
    nrow(panel_e_df),
    "."
  )
}

p_or <- ggplot(
  panel_e_df,
  aes(
    x = OR_per_1SD,
    y = label
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey55"
  ) +
  geom_errorbarh(
    aes(
      xmin = OR_CI_low,
      xmax = OR_CI_high,
      color = method_group
    ),
    height = 0.16,
    linewidth = 0.65
  ) +
  geom_point(
    aes(
      color = method_group,
      shape = significance
    ),
    size = 2.9,
    stroke = 0.9
  ) +
  scale_color_manual(
    values = method_palette
  ) +
  scale_shape_manual(
    values = significance_shape
  ) +
  scale_x_continuous(
    name = "Odds ratio per 1-SD higher score",
    breaks = pretty_breaks(
      n = 5
    )
  ) +
  scale_y_discrete(
    name = NULL
  ) +
  labs(
    title = "E  Association with 3-year incident F/G disease",
    subtitle = "Logistic models adjusted for age at imaging, sex, smoking, and BMI",
    caption = "Filled points indicate nominal P < 0.05."
  ) +
  theme_epoch()


# ------------------------------------------------------------------------------
# 8. Panel F: incremental repeated OOF ROC-AUC vs covariates
#    Echoes survival Panel B
# ------------------------------------------------------------------------------

panel_f_df <- cv_tbl |>
  filter(
    model %in% adjusted_model_order
  ) |>
  mutate(
    marker = unname(
      adjusted_model_to_marker[
        model
      ]
    ),
    label = unname(
      pretty_marker[
        marker
      ]
    ),
    method_group = if_else(
      marker == "Brain_MRI_mortality_EPOCH",
      "Mortality EPOCH",
      "AI disease subtype"
    ),
    delta_auc_vs_covariates = ROC_AUC_repeated_OOF - base_auc,
    marker = factor(
      marker,
      levels = marker_order
    ),
    label = factor(
      label,
      levels = rev(
        unname(
          pretty_marker[
            marker_order
          ]
        )
      )
    )
  ) |>
  arrange(
    marker
  )

if (nrow(panel_f_df) != 10) {
  warning(
    "Expected 10 adjusted biomarker models for Panel F; found ",
    nrow(panel_f_df),
    "."
  )
}

delta_auc_range <- range(
  panel_f_df$delta_auc_vs_covariates,
  na.rm = TRUE
)

delta_auc_pad <- max(
  0.002,
  diff(
    delta_auc_range
  ) * 0.15
)

p_delta_auc <- ggplot(
  panel_f_df,
  aes(
    x = delta_auc_vs_covariates,
    y = label
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey55"
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = delta_auc_vs_covariates,
      yend = label,
      color = method_group
    ),
    linewidth = 1.0,
    alpha = 0.80
  ) +
  geom_point(
    aes(
      color = method_group,
      shape = method_group
    ),
    size = 3.0,
    stroke = 0.9
  ) +
  scale_color_manual(
    values = method_palette
  ) +
  scale_shape_manual(
    values = method_shape
  ) +
  scale_x_continuous(
    name = expression(
      Delta*" ROC-AUC vs covariate-only model"
    ),
    labels = label_number(
      accuracy = 0.001
    ),
    limits = c(
      min(
        -0.005,
        delta_auc_range[1] - delta_auc_pad
      ),
      delta_auc_range[2] + delta_auc_pad
    )
  ) +
  scale_y_discrete(
    name = NULL
  ) +
  labs(
    title = "F  Incremental discrimination",
    subtitle = paste0(
      cv_folds,
      "-fold CV x ",
      cv_repeats,
      " repeats; covariate-only ROC-AUC = ",
      sprintf(
        "%.3f",
        base_auc
      )
    ),
    caption = "Positive values indicate higher repeated out-of-fold ROC-AUC than covariates alone."
  ) +
  theme_epoch()


# ------------------------------------------------------------------------------
# 9. Panel G: direct paired-bootstrap EPOCH-vs-subtype ROC-AUC comparison
#    Echoes survival Panel C
# ------------------------------------------------------------------------------

panel_g_df <- bootstrap_tbl |>
  filter(
    context == "adjusted",
    comparison %in% paste0(
      "EPOCH_vs_",
      c(
        "AD1",
        "AD2",
        "ASD1",
        "ASD2",
        "ASD3",
        "LLD1",
        "LLD2",
        "SCZ1",
        "SCZ2"
      )
    )
  ) |>
  mutate(
    subtype = str_remove(
      comparison,
      "^EPOCH_vs_"
    ),
    subtype = factor(
      subtype,
      levels = c(
        "AD1",
        "AD2",
        "ASD1",
        "ASD2",
        "ASD3",
        "LLD1",
        "LLD2",
        "SCZ1",
        "SCZ2"
      )
    ),
    subtype = fct_rev(
      subtype
    )
  )

if (nrow(panel_g_df) != 9) {
  warning(
    "Expected 9 adjusted EPOCH-vs-subtype bootstrap rows for Panel G; found ",
    nrow(panel_g_df),
    "."
  )
}

g_ci_range <- range(
  c(
    panel_g_df$ROC_AUC_CI_low,
    panel_g_df$ROC_AUC_CI_high
  ),
  na.rm = TRUE
)

g_ci_pad <- max(
  0.002,
  diff(
    g_ci_range
  ) * 0.10
)

n_boot <- unique(
  panel_g_df$successful_bootstrap_replicates
)

n_boot_label <- if (
  length(
    n_boot
  ) == 1
) {
  as.character(
    n_boot
  )
} else {
  "paired"
}

p_boot_auc <- ggplot(
  panel_g_df,
  aes(
    x = delta_ROC_AUC,
    y = subtype
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey55"
  ) +
  geom_errorbarh(
    aes(
      xmin = ROC_AUC_CI_low,
      xmax = ROC_AUC_CI_high
    ),
    height = 0.16,
    linewidth = 0.7,
    color = col_epoch
  ) +
  geom_point(
    size = 3.1,
    color = col_epoch
  ) +
  scale_x_continuous(
    name = expression(
      Delta*" ROC-AUC (EPOCH - disease subtype)"
    ),
    labels = label_number(
      accuracy = 0.001
    ),
    limits = c(
      min(
        g_ci_range[1] - g_ci_pad,
        -0.005
      ),
      g_ci_range[2] + g_ci_pad
    )
  ) +
  scale_y_discrete(
    name = NULL
  ) +
  labs(
    title = "G  Direct discrimination comparison",
    subtitle = paste0(
      "Paired bootstrap 95% confidence intervals (",
      n_boot_label,
      " replicates)"
    ),
    caption = paste0(
      "Positive values favor EPOCH. ",
      "A 95% CI crossing zero does not establish significant ROC-AUC superiority."
    )
  ) +
  theme_epoch()


# ------------------------------------------------------------------------------
# 10. Panel H: EPOCH versus combined 9-subtype panel
#     Echoes survival Panel D
#
# NOTE:
# This uses repeated out-of-fold ROC-AUC, not apparent/in-sample AUC.
# A dot plot is used rather than bars because the AUCs are close together and
# a truncated bar baseline would be visually misleading.
# ------------------------------------------------------------------------------

get_cv_value <- function(model_name, column_name) {
  tmp <- cv_tbl |>
    filter(
      model == model_name
    )
  
  if (nrow(tmp) != 1) {
    stop(
      "Expected exactly one row for model: ",
      model_name
    )
  }
  
  tmp[[column_name]][[1]]
}

auc_cov <- get_cv_value(
  "Covariates_only",
  "ROC_AUC_repeated_OOF"
)

auc_epoch <- get_cv_value(
  "EPOCH_adjusted",
  "ROC_AUC_repeated_OOF"
)

auc_all9 <- get_cv_value(
  "All9_adjusted",
  "ROC_AUC_repeated_OOF"
)

auc_all10 <- get_cv_value(
  "All9_plus_EPOCH_adjusted",
  "ROC_AUC_repeated_OOF"
)

boot_epoch_vs_all9 <- bootstrap_tbl |>
  filter(
    comparison == "EPOCH_vs_all9",
    context == "adjusted"
  )

if (nrow(boot_epoch_vs_all9) != 1) {
  stop(
    "Expected exactly one adjusted EPOCH_vs_all9 bootstrap row."
  )
}

boot_increment <- bootstrap_tbl |>
  filter(
    comparison == "EPOCH_increment_beyond_all9",
    context == "adjusted"
  )

if (nrow(boot_increment) != 1) {
  stop(
    "Expected exactly one adjusted EPOCH_increment_beyond_all9 bootstrap row."
  )
}

epoch_vs_all9_delta <- boot_epoch_vs_all9$delta_ROC_AUC[[1]]
epoch_vs_all9_lo <- boot_epoch_vs_all9$ROC_AUC_CI_low[[1]]
epoch_vs_all9_hi <- boot_epoch_vs_all9$ROC_AUC_CI_high[[1]]
epoch_vs_all9_p <- boot_epoch_vs_all9$ROC_AUC_bootstrap_p[[1]]

increment_delta <- boot_increment$delta_ROC_AUC[[1]]
increment_lo <- boot_increment$ROC_AUC_CI_low[[1]]
increment_hi <- boot_increment$ROC_AUC_CI_high[[1]]
increment_p <- boot_increment$ROC_AUC_bootstrap_p[[1]]

lrt_epoch_beyond_all9 <- combined_tbl$lrt_p_EPOCH_beyond_all9[[1]]
lrt_all9_beyond_epoch <- combined_tbl$lrt_p_all9_beyond_EPOCH[[1]]

panel_h_df <- tibble(
  model = factor(
    c(
      "Covariates only",
      "Covariates + EPOCH",
      "Covariates + 9 subtypes",
      "Covariates + 9 subtypes + EPOCH"
    ),
    levels = c(
      "Covariates only",
      "Covariates + EPOCH",
      "Covariates + 9 subtypes",
      "Covariates + 9 subtypes + EPOCH"
    )
  ),
  auc = c(
    auc_cov,
    auc_epoch,
    auc_all9,
    auc_all10
  ),
  model_group = c(
    "Reference",
    "Mortality EPOCH",
    "AI disease subtype",
    "Combined"
  )
)

h_min <- min(
  panel_h_df$auc,
  na.rm = TRUE
)

h_max <- max(
  panel_h_df$auc,
  na.rm = TRUE
)

h_pad <- max(
  0.006,
  diff(
    range(
      panel_h_df$auc,
      na.rm = TRUE
    )
  ) * 0.20
)

p_combined_auc <- ggplot(
  panel_h_df,
  aes(
    x = model,
    y = auc,
    color = model_group
  )
) +
  geom_hline(
    yintercept = auc_cov,
    linetype = "dashed",
    linewidth = 0.5,
    color = col_reference
  ) +
  geom_segment(
    aes(
      x = model,
      xend = model,
      y = auc_cov,
      yend = auc,
      color = model_group
    ),
    linewidth = 1.15,
    alpha = 0.78
  ) +
  geom_point(
    size = 4.0
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.3f",
        auc
      )
    ),
    color = "black",
    vjust = -1.0,
    size = 3.6,
    fontface = "bold"
  ) +
  scale_color_manual(
    values = combined_palette
  ) +
  scale_y_continuous(
    name = "Repeated out-of-fold ROC-AUC",
    labels = label_number(
      accuracy = 0.001
    ),
    limits = c(
      h_min - h_pad,
      h_max + h_pad
    ),
    breaks = pretty_breaks(
      n = 5
    )
  ) +
  scale_x_discrete(
    name = NULL,
    labels = c(
      "Covariates only" = "Covariates\nonly",
      "Covariates + EPOCH" = "Covariates\n+ EPOCH",
      "Covariates + 9 subtypes" = "Covariates\n+ 9 subtypes",
      "Covariates + 9 subtypes + EPOCH" = "Covariates\n+ 9 subtypes\n+ EPOCH"
    )
  ) +
  labs(
    title = "H  EPOCH versus the combined subtype panel",
    subtitle = paste0(
      "EPOCH vs all 9: ",
      "ΔAUC = ",
      format_delta(
        epoch_vs_all9_delta,
        3
      ),
      " [",
      sprintf(
        "%.3f",
        epoch_vs_all9_lo
      ),
      ", ",
      sprintf(
        "%.3f",
        epoch_vs_all9_hi
      ),
      "], ",
      format_p(
        epoch_vs_all9_p
      )
    ),
    caption = paste0(
      "Adding EPOCH to all 9 subtypes: ",
      "ΔAUC = ",
      format_delta(
        increment_delta,
        3
      ),
      " [",
      sprintf(
        "%.3f",
        increment_lo
      ),
      ", ",
      sprintf(
        "%.3f",
        increment_hi
      ),
      "], ",
      format_p(
        increment_p
      ),
      ". Logistic LRT: EPOCH beyond all 9, ",
      format_p(
        lrt_epoch_beyond_all9
      ),
      "; all 9 beyond EPOCH, ",
      format_p(
        lrt_all9_beyond_epoch
      ),
      "."
    )
  ) +
  theme_epoch() +
  theme(
    axis.text.x = element_text(
      size = 9.2,
      lineheight = 0.95
    ),
    plot.subtitle = element_text(
      size = 10.3
    ),
    plot.caption = element_text(
      size = 8.8,
      hjust = 0
    )
  )


# ------------------------------------------------------------------------------
# 11. Save individual panels
# ------------------------------------------------------------------------------

save_panel <- function(
    plot_object,
    stem,
    width = single_width,
    height = single_height
) {
  ggsave(
    filename = file.path(
      output_dir,
      paste0(
        prefix,
        "_",
        stem,
        ".pdf"
      )
    ),
    plot = plot_object,
    width = width,
    height = height,
    device = cairo_pdf
  )
  
  ggsave(
    filename = file.path(
      output_dir,
      paste0(
        prefix,
        "_",
        stem,
        ".png"
      )
    ),
    plot = plot_object,
    width = width,
    height = height,
    dpi = 400,
    bg = "white"
  )
}

save_panel(
  p_or,
  "E_OR_forest"
)

save_panel(
  p_delta_auc,
  "F_delta_ROC_AUC_vs_covariates"
)

save_panel(
  p_boot_auc,
  "G_bootstrap_EPOCH_minus_subtypes_ROC_AUC"
)

save_panel(
  p_combined_auc,
  "H_combined_9subtypes_vs_EPOCH_ROC_AUC"
)


# ------------------------------------------------------------------------------
# 12. Combined E-F-G-H figure
# ------------------------------------------------------------------------------

combined_figure <- (
  p_or |
    p_delta_auc
) /
  (
    p_boot_auc |
      p_combined_auc
  ) +
  plot_annotation(
    title = paste0(
      "Brain MRI mortality EPOCH versus AI-derived disease subtypes ",
      "for 3-year incident mental or nervous-system disorders"
    ),
    subtitle = paste0(
      "Fixed-horizon logistic sensitivity analysis; held-out EPOCH test split; ",
      "N = ",
      format(
        common_n,
        big.mark = ","
      ),
      " (",
      common_cases,
      " cases, ",
      common_controls,
      " controls)"
    ),
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 15
      ),
      plot.subtitle = element_text(
        size = 11
      )
    )
  )

ggsave(
  filename = file.path(
    output_dir,
    paste0(
      prefix,
      "_combined_EFGH.pdf"
    )
  ),
  plot = combined_figure,
  width = combined_width,
  height = combined_height,
  device = cairo_pdf
)

ggsave(
  filename = file.path(
    output_dir,
    paste0(
      prefix,
      "_combined_EFGH.png"
    )
  ),
  plot = combined_figure,
  width = combined_width,
  height = combined_height,
  dpi = 400,
  bg = "white"
)


# ------------------------------------------------------------------------------
# 13. Save compact plotting summary
# ------------------------------------------------------------------------------

single_summary <- single_tbl |>
  filter(
    marker %in% marker_order
  ) |>
  select(
    marker,
    N,
    N_case,
    N_control,
    OR_per_1SD,
    OR_CI_low,
    OR_CI_high,
    p_marker,
    lrt_p_vs_base
  )

cv_summary <- panel_f_df |>
  transmute(
    marker = as.character(
      marker
    ),
    adjusted_ROC_AUC = ROC_AUC_repeated_OOF,
    adjusted_AUPRC = AUPRC_repeated_OOF,
    delta_ROC_AUC_vs_covariates = delta_auc_vs_covariates
  )

bootstrap_summary <- panel_g_df |>
  transmute(
    marker = as.character(
      subtype
    ),
    EPOCH_minus_marker_delta_ROC_AUC = delta_ROC_AUC,
    EPOCH_minus_marker_ROC_AUC_CI_low = ROC_AUC_CI_low,
    EPOCH_minus_marker_ROC_AUC_CI_high = ROC_AUC_CI_high,
    EPOCH_minus_marker_ROC_AUC_bootstrap_p = ROC_AUC_bootstrap_p
  )

plot_summary <- single_summary |>
  left_join(
    cv_summary,
    by = "marker"
  ) |>
  left_join(
    bootstrap_summary,
    by = "marker"
  )

write_tsv(
  plot_summary,
  file.path(
    output_dir,
    paste0(
      prefix,
      "_plotting_summary.tsv"
    )
  )
)


# ------------------------------------------------------------------------------
# 14. Print key results
# ------------------------------------------------------------------------------

message("============================================================")
message("Fixed 3-year logistic E-F-G-H plotting complete.")
message("Input directory: ", input_dir)
message("Output directory: ", output_dir)
message("")
message(
  "Common sample: N = ",
  common_n,
  "; cases = ",
  common_cases,
  "; controls = ",
  common_controls
)
message(
  "Covariate-only repeated OOF ROC-AUC = ",
  sprintf(
    "%.4f",
    auc_cov
  )
)
message(
  "Covariates + EPOCH ROC-AUC = ",
  sprintf(
    "%.4f",
    auc_epoch
  )
)
message(
  "Covariates + all 9 subtypes ROC-AUC = ",
  sprintf(
    "%.4f",
    auc_all9
  )
)
message(
  "Covariates + all 9 + EPOCH ROC-AUC = ",
  sprintf(
    "%.4f",
    auc_all10
  )
)
message(
  "EPOCH vs all 9 subtypes: delta AUC = ",
  sprintf(
    "%+.4f",
    epoch_vs_all9_delta
  ),
  "; 95% CI ",
  sprintf(
    "%.4f",
    epoch_vs_all9_lo
  ),
  " to ",
  sprintf(
    "%.4f",
    epoch_vs_all9_hi
  ),
  "; ",
  format_p(
    epoch_vs_all9_p
  )
)
message(
  "Logistic LRT EPOCH beyond all 9 subtypes: ",
  format_p(
    lrt_epoch_beyond_all9
  )
)
message(
  "Logistic LRT all 9 subtypes beyond EPOCH: ",
  format_p(
    lrt_all9_beyond_epoch
  )
)
message("============================================================")