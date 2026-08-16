#!/usr/bin/env Rscript

# =============================================================================
# Multi-panel figure: influence of mortality time-to-event training horizon
# on the Endocrine metabolomics EPOCH mortality clock
#
# Main figure panels
#   A. Pairwise Pearson correlations of baseline EPOCH acceleration scores
#      (held-out mortality-clock test set)
#   B. Mortality discrimination (Uno C-index) at common 5- and 10-year
#      evaluation horizons for 5y-, 10y-, and full-follow-up-trained clocks
#   C. Incident disease associations (HR per 1 SD higher EPOCH) for five major
#      diseases, comparing the three mortality-training horizons
#   D. Paired bootstrap differences in disease-prediction C-index between the
#      three horizon-specific clocks (covariates + EPOCH models)
#
# Usage
#   Rscript plot_endocrine_epoch_time_horizon_main_figure.R
#
# or
#   Rscript plot_endocrine_epoch_time_horizon_main_figure.R \
#     /path/to/Endocrine_metabolomics_mortality_horizon_clocks \
#     /path/to/output_figures
#
# Required R packages:
#   ggplot2, dplyr, patchwork
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

root_dir <- if (length(args) >= 1) {
  normalizePath(args[1], mustWork = FALSE)
} else {
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/Endocrine_metabolomics_mortality_horizon_clocks"
}

out_dir <- if (length(args) >= 2) {
  normalizePath(args[2], mustWork = FALSE)
} else {
  file.path(root_dir, "figures_time_horizon")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("ggplot2", "dplyr", "patchwork")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script."
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

read_tsv_checked <- function(path) {
  if (!file.exists(path)) {
    stop("Required input file does not exist:\n", path)
  }
  read.delim(
    path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )
}

assert_columns <- function(df, cols, file_label) {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    stop(
      file_label, " is missing required columns: ",
      paste(missing, collapse = ", ")
    )
  }
}

horizon_order <- c("5y", "10y", "full")
horizon_labels <- c(
  "5y" = "5-year",
  "10y" = "10-year",
  "full" = "Full follow-up"
)

# Okabe-Ito-inspired, color-blind-friendly palette.
horizon_colors <- c(
  "5y" = "#0072B2",
  "10y" = "#E69F00",
  "full" = "#009E73"
)

pair_order <- c("5y_minus_10y", "5y_minus_full", "10y_minus_full")
pair_labels <- c(
  "5y_minus_10y" = "5y - 10y",
  "5y_minus_full" = "5y - Full",
  "10y_minus_full" = "10y - Full"
)
pair_colors <- c(
  "5y_minus_10y" = "#0072B2",
  "5y_minus_full" = "#CC79A7",
  "10y_minus_full" = "#009E73"
)

disease_order <- c(
  "all_cause_dementia",
  "asthma",
  "myocardial_infarction",
  "copd",
  "stroke"
)
disease_labels <- c(
  "all_cause_dementia" = "All-cause dementia",
  "asthma" = "Asthma",
  "myocardial_infarction" = "Myocardial infarction",
  "copd" = "COPD",
  "stroke" = "Stroke"
)

base_theme <- theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 11.5, hjust = 0),
    plot.subtitle = element_text(size = 9.5, margin = margin(b = 5)),
    axis.title = element_text(size = 10.5),
    axis.text = element_text(size = 9.5),
    legend.title = element_text(size = 9.5),
    legend.text = element_text(size = 9),
    legend.position = "bottom",
    plot.margin = margin(6, 8, 6, 6)
  )

# -----------------------------------------------------------------------------
# Input files
# -----------------------------------------------------------------------------

score_corr_file <- file.path(
  root_dir,
  "endocrine_metabolomics_mortality_horizon_clocks_cross_horizon_score_correlations.tsv"
)

mortality_eval_file <- file.path(
  root_dir,
  "endocrine_metabolomics_mortality_horizon_clocks_common_test_mortality_evaluation.tsv"
)

disease_dir <- file.path(root_dir, "disease_onset_test")
if (!dir.exists(disease_dir)) {
  stop("Disease-onset results directory does not exist:\n", disease_dir)
}

# -----------------------------------------------------------------------------
# Panel A: EPOCH acceleration correlations across training horizons
# -----------------------------------------------------------------------------

score_corr <- read_tsv_checked(score_corr_file)
assert_columns(
  score_corr,
  c("score_type", "subset", "horizon_a", "horizon_b", "pearson_r"),
  basename(score_corr_file)
)

corr_pairs <- score_corr %>%
  filter(score_type == "acceleration_z", subset == "test") %>%
  select(horizon_a, horizon_b, pearson_r)

if (nrow(corr_pairs) != 3) {
  warning(
    "Expected 3 pairwise acceleration_z/test correlations, found ",
    nrow(corr_pairs), "."
  )
}

corr_diag <- data.frame(
  horizon_a = horizon_order,
  horizon_b = horizon_order,
  pearson_r = 1,
  stringsAsFactors = FALSE
)

corr_mirror <- corr_pairs %>%
  transmute(
    horizon_a = horizon_b,
    horizon_b = horizon_a,
    pearson_r = pearson_r
  )

corr_plot_df <- bind_rows(corr_pairs, corr_mirror, corr_diag) %>%
  mutate(
    horizon_a = factor(horizon_a, levels = horizon_order, labels = horizon_labels[horizon_order]),
    horizon_b = factor(horizon_b, levels = horizon_order, labels = horizon_labels[horizon_order])
  )

corr_min <- suppressWarnings(min(corr_plot_df$pearson_r, na.rm = TRUE))
fill_min <- min(0.95, floor(corr_min * 100) / 100)
fill_min <- max(0, fill_min)

pA <- ggplot(corr_plot_df, aes(x = horizon_a, y = horizon_b, fill = pearson_r)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.3f", pearson_r)), size = 3.7, fontface = "bold") +
  scale_fill_gradient(
    low = "#F7FBFF",
    high = "#08519C",
    limits = c(fill_min, 1),
    breaks = unique(round(c(fill_min, (fill_min + 1) / 2, 1), 2)),
    name = "Pearson r"
  ) +
  coord_equal() +
  labs(
    title = "EPOCH acceleration stability",
    subtitle = "Held-out test set",
    x = "Mortality-training horizon",
    y = "Mortality-training horizon"
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "right"
  )

# -----------------------------------------------------------------------------
# Panel B: mortality discrimination at COMMON evaluation horizons
# -----------------------------------------------------------------------------

mortality_eval <- read_tsv_checked(mortality_eval_file)
assert_columns(
  mortality_eval,
  c("model_horizon", "evaluation_horizon_years", "uno_c"),
  basename(mortality_eval_file)
)

mortality_plot_df <- mortality_eval %>%
  filter(evaluation_horizon_years %in% c(5, 10)) %>%
  mutate(
    model_horizon = factor(model_horizon, levels = horizon_order),
    training_label = factor(
      as.character(model_horizon),
      levels = horizon_order,
      labels = horizon_labels[horizon_order]
    ),
    evaluation_label = factor(
      evaluation_horizon_years,
      levels = c(5, 10),
      labels = c("5-year mortality", "10-year mortality")
    )
  )

if (nrow(mortality_plot_df) != 6) {
  warning("Expected 6 mortality-evaluation rows (3 clocks x 2 horizons); found ", nrow(mortality_plot_df), ".")
}

mort_range <- range(mortality_plot_df$uno_c, na.rm = TRUE)
mort_pad <- max(0.004, diff(mort_range) * 0.30)
mort_ylim <- c(mort_range[1] - mort_pad, mort_range[2] + mort_pad)

pB <- ggplot(
  mortality_plot_df,
  aes(x = training_label, y = uno_c, group = 1)
) +
  geom_line(color = "grey55", linewidth = 0.65) +
  geom_point(
    aes(color = model_horizon),
    size = 3.2
  ) +
  geom_text(
    aes(label = sprintf("%.3f", uno_c)),
    vjust = -0.9,
    size = 3.0,
    color = "black"
  ) +
  facet_wrap(~evaluation_label, nrow = 1) +
  scale_color_manual(
    values = horizon_colors,
    breaks = horizon_order,
    labels = horizon_labels[horizon_order],
    name = "Training horizon"
  ) +
  coord_cartesian(ylim = mort_ylim, clip = "off") +
  labs(
    title = "Mortality discrimination",
    subtitle = "All clocks evaluated on the same held-out participants and outcome horizon",
    x = "Mortality-training horizon",
    y = "Uno C-index"
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 9.5),
    legend.position = "none"
  )

# -----------------------------------------------------------------------------
# Read and combine five disease-association outputs
# -----------------------------------------------------------------------------

disease_result_list <- list()
disease_delta_list <- list()

for (endpoint in disease_order) {
  main_file <- file.path(
    disease_dir,
    paste0("cox_endocrine_horizon_clocks_", endpoint, "_test.tsv")
  )
  delta_file <- file.path(
    disease_dir,
    paste0("cox_endocrine_horizon_clocks_", endpoint, "_test_paired_delta_cindex.tsv")
  )

  main_df <- read_tsv_checked(main_file)
  delta_df <- read_tsv_checked(delta_file)

  assert_columns(
    main_df,
    c(
      "endpoint", "endpoint_label", "clock_training_horizon", "status",
      "clock_hr_per_1sd", "clock_ci_lo", "clock_ci_hi", "clock_p",
      "full_model_cindex", "covariate_only_cindex",
      "delta_cindex_full_minus_covariates"
    ),
    basename(main_file)
  )

  assert_columns(
    delta_df,
    c(
      "endpoint", "endpoint_label", "prediction_type", "comparison",
      "delta_cindex_a_minus_b", "delta_cindex_ci_lower",
      "delta_cindex_ci_upper"
    ),
    basename(delta_file)
  )

  disease_result_list[[endpoint]] <- main_df
  disease_delta_list[[endpoint]] <- delta_df
}

disease_results <- bind_rows(disease_result_list) %>%
  filter(status == "ok") %>%
  mutate(
    endpoint = factor(endpoint, levels = disease_order),
    disease_label = factor(
      as.character(endpoint),
      levels = disease_order,
      labels = disease_labels[disease_order]
    ),
    clock_training_horizon = factor(clock_training_horizon, levels = horizon_order)
  )

disease_deltas <- bind_rows(disease_delta_list) %>%
  filter(prediction_type == "covariates_plus_epoch") %>%
  mutate(
    endpoint = factor(endpoint, levels = disease_order),
    disease_label = factor(
      as.character(endpoint),
      levels = disease_order,
      labels = disease_labels[disease_order]
    ),
    comparison = factor(comparison, levels = pair_order)
  )

# Save compact combined tables used by the figure.
write.table(
  disease_results,
  file = file.path(out_dir, "combined_endocrine_horizon_disease_results.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)
write.table(
  disease_deltas,
  file = file.path(out_dir, "combined_endocrine_horizon_disease_delta_cindex.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)

# -----------------------------------------------------------------------------
# Panel C: incident disease HRs across training horizons
# -----------------------------------------------------------------------------

pd_hr <- position_dodge(width = 0.68)

pC <- ggplot(
  disease_results,
  aes(
    x = disease_label,
    y = clock_hr_per_1sd,
    color = clock_training_horizon,
    group = clock_training_horizon
  )
) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.55, color = "grey45") +
  geom_errorbar(
    aes(ymin = clock_ci_lo, ymax = clock_ci_hi),
    width = 0.18,
    linewidth = 0.75,
    position = pd_hr
  ) +
  geom_point(size = 2.7, position = pd_hr) +
  coord_flip() +
  scale_y_log10() +
  scale_x_discrete(limits = rev(disease_labels[disease_order])) +
  scale_color_manual(
    values = horizon_colors,
    breaks = horizon_order,
    labels = horizon_labels[horizon_order],
    name = "Training horizon"
  ) +
  labs(
    title = "Incident disease associations",
    subtitle = "Hazard ratio per 1 SD higher baseline Endocrine EPOCH",
    x = NULL,
    y = "Hazard ratio (95% CI; log scale)"
  ) +
  base_theme +
  theme(legend.position = "bottom")

# -----------------------------------------------------------------------------
# Panel D: paired disease C-index differences across horizon-specific clocks
# -----------------------------------------------------------------------------

pd_delta <- position_dodge(width = 0.68)

pD <- ggplot(
  disease_deltas,
  aes(
    x = disease_label,
    y = delta_cindex_a_minus_b,
    color = comparison,
    group = comparison
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.55, color = "grey40") +
  geom_errorbar(
    aes(
      ymin = delta_cindex_ci_lower,
      ymax = delta_cindex_ci_upper
    ),
    width = 0.18,
    linewidth = 0.75,
    position = pd_delta
  ) +
  geom_point(size = 2.6, position = pd_delta) +
  coord_flip() +
  scale_x_discrete(limits = rev(disease_labels[disease_order])) +
  scale_color_manual(
    values = pair_colors,
    breaks = pair_order,
    labels = pair_labels[pair_order],
    name = "Clock comparison"
  ) +
  labs(
    title = "Differences in disease discrimination",
    subtitle = "Paired bootstrap: covariates + EPOCH models",
    x = NULL,
    y = expression(Delta * " C-index (model A - model B), 95% CI")
  ) +
  base_theme +
  theme(legend.position = "bottom")

# -----------------------------------------------------------------------------
# Assemble 2 x 2 main figure
# -----------------------------------------------------------------------------

main_figure <- (pA | pB) / (pC | pD) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 15),
      plot.tag.position = c(0.01, 0.99)
    )
  ) &
  theme(plot.tag = element_text(face = "bold"))

pdf_file <- file.path(out_dir, "Endocrine_EPOCH_time_horizon_main_figure.pdf")
png_file <- file.path(out_dir, "Endocrine_EPOCH_time_horizon_main_figure.png")
svg_file <- file.path(out_dir, "Endocrine_EPOCH_time_horizon_main_figure.svg")

ggsave(
  filename = pdf_file,
  plot = main_figure,
  width = 13.5,
  height = 10.2,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = png_file,
  plot = main_figure,
  width = 13.5,
  height = 10.2,
  units = "in",
  dpi = 400,
  bg = "white"
)

# SVG is optional because some clusters do not have svglite installed.
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    filename = svg_file,
    plot = main_figure,
    width = 13.5,
    height = 10.2,
    units = "in",
    device = svglite::svglite
  )
} else {
  message("Package 'svglite' not installed; skipping SVG export.")
}

# Also save individual panels for easy manuscript assembly/revision.
individual_plots <- list(
  A_epoch_acceleration_correlations = pA,
  B_common_horizon_mortality_cindex = pB,
  C_incident_disease_HR = pC,
  D_disease_delta_cindex = pD
)

for (nm in names(individual_plots)) {
  ggsave(
    filename = file.path(out_dir, paste0(nm, ".pdf")),
    plot = individual_plots[[nm]],
    width = 7.2,
    height = 5.4,
    units = "in",
    device = cairo_pdf
  )
  ggsave(
    filename = file.path(out_dir, paste0(nm, ".png")),
    plot = individual_plots[[nm]],
    width = 7.2,
    height = 5.4,
    units = "in",
    dpi = 400,
    bg = "white"
  )
}

# -----------------------------------------------------------------------------
# Optional secondary figure: disease C-index and incremental C-index
# This is useful if reviewers ask for the absolute predictive-performance view.
# -----------------------------------------------------------------------------

pS1 <- ggplot(
  disease_results,
  aes(
    x = clock_training_horizon,
    y = full_model_cindex,
    color = clock_training_horizon,
    group = 1
  )
) +
  geom_line(color = "grey55", linewidth = 0.65) +
  geom_point(size = 2.8) +
  facet_wrap(~disease_label, scales = "free_y", ncol = 3) +
  scale_color_manual(values = horizon_colors, guide = "none") +
  scale_x_discrete(labels = horizon_labels[horizon_order]) +
  labs(
    title = "Disease discrimination across EPOCH training horizons",
    x = "Mortality-training horizon",
    y = "Full-model Harrell C-index"
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 9)
  )

pS2 <- ggplot(
  disease_results,
  aes(
    x = clock_training_horizon,
    y = delta_cindex_full_minus_covariates,
    color = clock_training_horizon,
    group = 1
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_line(color = "grey55", linewidth = 0.65) +
  geom_point(size = 2.8) +
  facet_wrap(~disease_label, scales = "free_y", ncol = 3) +
  scale_color_manual(values = horizon_colors, guide = "none") +
  scale_x_discrete(labels = horizon_labels[horizon_order]) +
  labs(
    title = "Incremental disease discrimination beyond common covariates",
    x = "Mortality-training horizon",
    y = expression(Delta * " C-index vs covariate-only model")
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 9)
  )

secondary_figure <- pS1 / pS2 +
  plot_annotation(tag_levels = "A")

ggsave(
  filename = file.path(out_dir, "Endocrine_EPOCH_time_horizon_disease_cindex_supplement.pdf"),
  plot = secondary_figure,
  width = 11.5,
  height = 9.5,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(out_dir, "Endocrine_EPOCH_time_horizon_disease_cindex_supplement.png"),
  plot = secondary_figure,
  width = 11.5,
  height = 9.5,
  units = "in",
  dpi = 400,
  bg = "white"
)

message("============================================================")
message("Finished plotting Endocrine EPOCH time-horizon results")
message("Input root: ", root_dir)
message("Output directory: ", out_dir)
message("Main figure PDF: ", pdf_file)
message("Main figure PNG: ", png_file)
if (file.exists(svg_file)) message("Main figure SVG: ", svg_file)
message("============================================================")