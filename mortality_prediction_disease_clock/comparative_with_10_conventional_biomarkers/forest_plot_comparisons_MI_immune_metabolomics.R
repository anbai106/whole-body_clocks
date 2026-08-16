#!/usr/bin/env Rscript

# ==============================================================================
# Combined two-panel comparison:
#
#   Panel A: Hazard ratios (95% CI) for all-cause mortality
#   Panel B: Incremental discrimination (Delta C-index)
#
# MI immune-metabolomics EPOCH vs 10 conventional mortality biomarkers
#
# All predictors come from the same strict apple-to-apple analysis sample.
#
# NOTE:
# The summary TSV contains confidence intervals for HRs but does not contain
# confidence intervals for Delta C-index. Therefore Panel B shows point estimates
# of Delta C-index only. If bootstrap CIs are generated later, they can be added
# to Panel B as horizontal error bars.
#
# Outputs:
#   mi_immune_metabolomics_EPOCH_vs_10_biomarkers_combined_AB.pdf
#   mi_immune_metabolomics_EPOCH_vs_10_biomarkers_combined_AB.png
#   mi_immune_metabolomics_EPOCH_vs_10_biomarkers_combined_plot_data.tsv
# ==============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------------------------

input_file <- paste0(
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/",
  "comparative_with_10_conventional_biomarkers/",
  "2_mortality_comparison_mi_immune_metabolomics/",
  "mi_immune_metabolomics_EPOCH_vs_10_biomarkers_mortality_summary.tsv"
)

output_dir <- dirname(input_file)

pdf_file <- file.path(
  output_dir,
  "mi_immune_metabolomics_EPOCH_vs_10_biomarkers_combined_AB.pdf"
)

png_file <- file.path(
  output_dir,
  "mi_immune_metabolomics_EPOCH_vs_10_biomarkers_combined_AB.png"
)

plot_data_file <- file.path(
  output_dir,
  "mi_immune_metabolomics_EPOCH_vs_10_biomarkers_combined_plot_data.tsv"
)

# ------------------------------------------------------------------------------
# 2. Read data
# ------------------------------------------------------------------------------

if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file)
}

dat <- read_tsv(
  input_file,
  show_col_types = FALSE,
  progress = FALSE
)

required_cols <- c(
  "predictor_order",
  "predictor_type",
  "predictor",
  "n_analysis_rows",
  "n_deaths",
  "hr_per_1sd",
  "hr_ci_lower",
  "hr_ci_upper",
  "p_value",
  "cindex_baseline_covariates",
  "cindex_baseline_plus_predictor",
  "delta_cindex_vs_baseline"
)

missing_cols <- setdiff(required_cols, colnames(dat))

if (length(missing_cols) > 0) {
  stop(
    "Missing required column(s): ",
    paste(missing_cols, collapse = ", ")
  )
}

# ------------------------------------------------------------------------------
# 3. Prepare plotting data
# ------------------------------------------------------------------------------

format_p <- function(p) {
  if (is.na(p)) {
    return("NA")
  }
  
  if (p < 0.001) {
    return(format(p, scientific = TRUE, digits = 2))
  }
  
  sprintf("%.3f", p)
}

plot_dat <- dat %>%
  mutate(
    predictor_order = as.numeric(predictor_order),
    hr_per_1sd = as.numeric(hr_per_1sd),
    hr_ci_lower = as.numeric(hr_ci_lower),
    hr_ci_upper = as.numeric(hr_ci_upper),
    p_value = as.numeric(p_value),
    cindex_baseline_covariates = as.numeric(cindex_baseline_covariates),
    cindex_baseline_plus_predictor = as.numeric(cindex_baseline_plus_predictor),
    delta_cindex_vs_baseline = as.numeric(delta_cindex_vs_baseline),
    
    # Standardize the EPOCH display label regardless of the exact predictor
    # string written by the Python survival script.
    predictor_plot = case_when(
      predictor_type == "EPOCH" ~
        "MI immune-metabolomics EPOCH",
      predictor == "C-reactive protein (log1p)" ~
        "C-reactive protein",
      TRUE ~ predictor
    ),
    
    group = if_else(
      predictor_type == "EPOCH",
      "EPOCH",
      "Conventional biomarker"
    ),
    
    hr_ci_text = sprintf(
      "%.2f (%.2f\u2013%.2f)",
      hr_per_1sd,
      hr_ci_lower,
      hr_ci_upper
    ),
    
    delta_cindex_text = sprintf(
      "%+.4f",
      delta_cindex_vs_baseline
    ),
    
    p_text = vapply(
      p_value,
      format_p,
      FUN.VALUE = character(1)
    )
  ) %>%
  arrange(predictor_order)

# Preserve identical y-axis order across both panels, with EPOCH at the top.
predictor_levels <- rev(plot_dat$predictor_plot)

plot_dat$predictor_plot <- factor(
  plot_dat$predictor_plot,
  levels = predictor_levels
)

write_tsv(plot_dat, plot_data_file)

# ------------------------------------------------------------------------------
# 4. Study information
# ------------------------------------------------------------------------------

n_values <- unique(plot_dat$n_analysis_rows)
death_values <- unique(plot_dat$n_deaths)
baseline_c_values <- unique(plot_dat$cindex_baseline_covariates)

if (length(n_values) != 1) {
  warning("More than one analysis N is present in the input.")
}

if (length(death_values) != 1) {
  warning("More than one death count is present in the input.")
}

if (length(baseline_c_values) != 1) {
  warning("More than one baseline C-index is present in the input.")
}

analysis_n <- n_values[1]
analysis_deaths <- death_values[1]
baseline_cindex <- baseline_c_values[1]

subtitle_text <- paste0(
  "All-cause mortality; N = ",
  format(analysis_n, big.mark = ","),
  "; deaths = ",
  format(analysis_deaths, big.mark = ",")
)

# ------------------------------------------------------------------------------
# 5. Van Gogh-inspired palette and shapes
# ------------------------------------------------------------------------------

epoch_color <- "#E0A526"
biomarker_color <- "#345995"

group_colors <- c(
  "EPOCH" = epoch_color,
  "Conventional biomarker" = biomarker_color
)

group_shapes <- c(
  "EPOCH" = 18,
  "Conventional biomarker" = 16
)

# ------------------------------------------------------------------------------
# 6. Shared theme
# ------------------------------------------------------------------------------

common_theme <- theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 12.5,
      hjust = 0
    ),
    plot.subtitle = element_text(
      size = 10.2,
      hjust = 0,
      margin = margin(b = 8)
    ),
    axis.text.y = element_text(
      size = 10.2,
      color = "black"
    ),
    axis.text.x = element_text(
      size = 9.5,
      color = "black"
    ),
    axis.title.x = element_text(
      size = 10.5,
      margin = margin(t = 8)
    ),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 9.5),
    legend.key.width = unit(1.1, "lines"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 10, 8, 8)
  )

# ------------------------------------------------------------------------------
# 7. Panel A: HR forest plot
# ------------------------------------------------------------------------------

hr_lower <- min(plot_dat$hr_ci_lower, na.rm = TRUE)
hr_upper <- max(plot_dat$hr_ci_upper, na.rm = TRUE)

hr_xlim_lower <- exp(log(hr_lower) - 0.08)
hr_xlim_upper <- exp(log(hr_upper) + 0.08)

candidate_hr_breaks <- c(
  0.50, 0.60, 0.70, 0.75, 0.80, 0.90,
  1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.75, 2.00
)

hr_breaks <- candidate_hr_breaks[
  candidate_hr_breaks >= hr_xlim_lower &
    candidate_hr_breaks <= hr_xlim_upper
]

if (!any(abs(hr_breaks - 1) < 1e-12)) {
  hr_breaks <- sort(unique(c(hr_breaks, 1)))
}

panel_a <- ggplot(
  plot_dat,
  aes(
    x = hr_per_1sd,
    y = predictor_plot,
    color = group,
    shape = group
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.55,
    color = "grey45"
  ) +
  geom_errorbarh(
    aes(
      xmin = hr_ci_lower,
      xmax = hr_ci_upper
    ),
    height = 0.16,
    linewidth = 0.75,
    show.legend = FALSE
  ) +
  geom_point(
    size = 3.5,
    stroke = 0.7
  ) +
  scale_color_manual(
    values = group_colors,
    breaks = c("EPOCH", "Conventional biomarker")
  ) +
  scale_shape_manual(
    values = group_shapes,
    breaks = c("EPOCH", "Conventional biomarker")
  ) +
  scale_x_log10(
    limits = c(hr_xlim_lower, hr_xlim_upper),
    breaks = hr_breaks,
    labels = label_number(accuracy = 0.01)
  ) +
  labs(
    title = "Hazard ratios",
    subtitle = subtitle_text,
    x = "Hazard ratio for all-cause mortality (per 1 SD)",
    y = NULL,
    color = NULL,
    shape = NULL
  ) +
  common_theme

# ------------------------------------------------------------------------------
# 8. Panel B: Delta C-index forest-style plot
# ------------------------------------------------------------------------------

max_delta <- max(plot_dat$delta_cindex_vs_baseline, na.rm = TRUE)
min_delta <- min(plot_dat$delta_cindex_vs_baseline, na.rm = TRUE)

delta_range <- max_delta - min_delta

if (!is.finite(delta_range) || delta_range == 0) {
  delta_range <- max(abs(max_delta), 0.001)
}

delta_pad <- max(0.0015, 0.08 * delta_range)
delta_lower <- min(0, min_delta) - delta_pad
delta_upper <- max(0, max_delta) + delta_pad

panel_b <- ggplot(
  plot_dat,
  aes(
    x = delta_cindex_vs_baseline,
    y = predictor_plot,
    color = group,
    shape = group
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.55,
    color = "grey45"
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = delta_cindex_vs_baseline,
      y = predictor_plot,
      yend = predictor_plot
    ),
    linewidth = 0.65,
    alpha = 0.65,
    show.legend = FALSE
  ) +
  geom_point(
    size = 3.5,
    stroke = 0.7
  ) +
  scale_color_manual(
    values = group_colors,
    breaks = c("EPOCH", "Conventional biomarker")
  ) +
  scale_shape_manual(
    values = group_shapes,
    breaks = c("EPOCH", "Conventional biomarker")
  ) +
  scale_x_continuous(
    limits = c(delta_lower, delta_upper),
    labels = label_number(
      accuracy = 0.005,
      trim = TRUE
    ),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  labs(
    title = "Incremental discrimination",
    subtitle = paste0(
      "Baseline C-index = ",
      sprintf("%.3f", baseline_cindex)
    ),
    x = expression(Delta * " C-index vs baseline covariates"),
    y = NULL,
    color = NULL,
    shape = NULL
  ) +
  common_theme +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()
  )

# ------------------------------------------------------------------------------
# 9. Combine Panels A-B
# ------------------------------------------------------------------------------

combined_plot <- (
  panel_a + panel_b +
    plot_layout(
      widths = c(1.38, 1.0),
      guides = "collect"
    )
) +
  plot_annotation(
    title = paste0(
      "MI immune-metabolomics EPOCH versus conventional biomarkers ",
      "for all-cause mortality"
    ),
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 14,
        hjust = 0
      ),
      plot.tag = element_text(
        face = "bold",
        size = 14
      )
    )
  ) &
  theme(
    legend.position = "bottom"
  )

# ------------------------------------------------------------------------------
# 10. Save combined figure
# ------------------------------------------------------------------------------

ggsave(
  filename = pdf_file,
  plot = combined_plot,
  width = 12.4,
  height = 6.3,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = png_file,
  plot = combined_plot,
  width = 12.4,
  height = 6.3,
  units = "in",
  dpi = 600,
  bg = "white"
)

# ------------------------------------------------------------------------------
# 11. Console summary
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("MI immune-metabolomics combined Panels A-B completed\n")
cat("============================================================\n")
cat("Input:\n  ", input_file, "\n", sep = "")
cat("Analysis N: ", format(analysis_n, big.mark = ","), "\n", sep = "")
cat("Deaths: ", format(analysis_deaths, big.mark = ","), "\n", sep = "")
cat(
  "Baseline C-index: ",
  sprintf("%.4f", baseline_cindex),
  "\n",
  sep = ""
)

cat("\nVan Gogh-inspired coding:\n")
cat("  EPOCH: golden-yellow diamond\n")
cat("  Conventional biomarkers: deep-blue circles\n")

cat("\nPDF:\n  ", pdf_file, "\n", sep = "")
cat("\nPNG:\n  ", png_file, "\n", sep = "")
cat("\nPlot data:\n  ", plot_data_file, "\n", sep = "")

cat("\nPredictor results:\n")

print(
  plot_dat %>%
    select(
      predictor,
      hr_per_1sd,
      hr_ci_lower,
      hr_ci_upper,
      p_value,
      cindex_baseline_plus_predictor,
      delta_cindex_vs_baseline
    ),
  n = Inf
)

cat("\nNOTE: Panel B displays Delta C-index point estimates only.\n")
cat("Bootstrap confidence intervals can be added when available.\n")
cat("============================================================\n")