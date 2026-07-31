#!/usr/bin/env Rscript

# ============================================================
# Plot MCI-to-AD cumulative CSF + AD EPOCH results
#
# Panel A: Hazard ratios from the final cumulative Cox model
#          for CSF amyloid, CSF total tau, CSF p-tau,
#          and AD EPOCH only. Age and sex are omitted.
#
# Panel B: Five-fold cross-validated cumulative C-index across
#          nested models, with the incremental p-value for adding
#          AD EPOCH beyond demographics + CSF amyloid + CSF tau.
#
# Color palette:
#   Okabe-Ito colorblind-friendly palette.
#
# Required packages:
#   readr, dplyr, stringr, ggplot2, patchwork, scales
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# ----------------------------
# User-editable paths
# ----------------------------

input_dir <- "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/adni_lepoch/results_brain_mri_ad_lepoch_mci2ad_cumulative_csf_epoch"

prefix <- "adni_mci2ad_cumulative_csf_epoch"

coefficients_file <- file.path(input_dir, paste0(prefix, "_coefficients.tsv"))
performance_file  <- file.path(input_dir, paste0(prefix, "_performance.tsv"))
comparisons_file  <- file.path(input_dir, paste0(prefix, "_comparisons.tsv"))

output_dir <- file.path(input_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_stem <- file.path(
  output_dir,
  paste0(prefix, "_HR_and_CV_Cindex_colorblind")
)

final_model_name <- "M3_plus_AD_EPOCH"
previous_model_name <- "M2_plus_tau"
analysis_type_keep <- "common_complete_case"

# ----------------------------
# Colorblind-friendly palette
# ----------------------------

# Okabe-Ito palette. These colors remain distinguishable for common forms
# of color-vision deficiency and also reproduce well in print.
cb_palette <- c(
  "CSF amyloid-beta" = "#0072B2",          # blue
  "CSF total tau" = "#E69F00",             # orange
  "CSF phosphorylated tau" = "#CC79A7",    # reddish purple
  "AD EPOCH" = "#009E73"                   # bluish green
)

model_palette <- c(
  "Age + sex" = "#999999",                 # neutral gray
  "+ CSF amyloid-beta" = "#0072B2",        # blue
  "+ CSF total tau + p-tau" = "#E69F00",   # orange
  "+ AD EPOCH" = "#009E73"                 # bluish green
)

line_color <- "#4D4D4D"
reference_color <- "#666666"
annotation_fill <- "#F7F7F7"

# ----------------------------
# Helpers
# ----------------------------

stop_if_missing <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    stop("Missing or empty required file: ", path, call. = FALSE)
  }
}

format_p <- function(p) {
  vapply(p, function(x) {
    if (length(x) == 0 || is.na(x)) return("P = NA")
    if (x < 0.001) return("P < 0.001")
    paste0("P = ", formatC(x, format = "f", digits = 3))
  }, character(1))
}

term_to_label <- function(term) {
  case_when(
    term == "Abeta_CSF" ~ "CSF amyloid-beta",
    term == "Tau_CSF" ~ "CSF total tau",
    term == "PTau_CSF" ~ "CSF phosphorylated tau",
    term == "adni_brain_mri_ad_lepoch_risk_score" ~ "AD EPOCH",
    TRUE ~ term
  )
}

model_to_label <- function(model) {
  case_when(
    model == "M0_age_sex" ~ "Age + sex",
    model == "M1_plus_amyloid" ~ "+ CSF amyloid-beta",
    model == "M2_plus_tau" ~ "+ CSF total tau + p-tau",
    model == "M3_plus_AD_EPOCH" ~ "+ AD EPOCH",
    TRUE ~ model
  )
}

# ----------------------------
# Read inputs
# ----------------------------

invisible(lapply(
  c(coefficients_file, performance_file, comparisons_file),
  stop_if_missing
))

coef_df <- read_tsv(coefficients_file, show_col_types = FALSE)
perf_df <- read_tsv(performance_file, show_col_types = FALSE)
comp_df <- read_tsv(comparisons_file, show_col_types = FALSE)

# ----------------------------
# Panel A: biomarker HR forest plot
# ----------------------------

required_coef_cols <- c(
  "model", "term", "hazard_ratio",
  "hazard_ratio_ci_lower", "hazard_ratio_ci_upper", "p_value"
)
missing_coef_cols <- setdiff(required_coef_cols, names(coef_df))
if (length(missing_coef_cols) > 0) {
  stop(
    "Coefficient file is missing required columns: ",
    paste(missing_coef_cols, collapse = ", "),
    call. = FALSE
  )
}

if ("analysis_type" %in% names(coef_df)) {
  coef_df <- coef_df %>% filter(analysis_type == analysis_type_keep)
}

# Age and sex are intentionally excluded from Panel A.
terms_keep <- c(
  "Abeta_CSF",
  "Tau_CSF",
  "PTau_CSF",
  "adni_brain_mri_ad_lepoch_risk_score"
)

panel_a_df <- coef_df %>%
  filter(model == final_model_name, term %in% terms_keep) %>%
  mutate(
    label = term_to_label(term),
    p_label = format_p(p_value),
    hr_label = paste0(
      formatC(hazard_ratio, format = "f", digits = 2),
      " (",
      formatC(hazard_ratio_ci_lower, format = "f", digits = 2),
      "-",
      formatC(hazard_ratio_ci_upper, format = "f", digits = 2),
      ")"
    )
  )

if (nrow(panel_a_df) == 0) {
  stop(
    "No biomarker coefficients found for final model '", final_model_name,
    "'. Available models: ", paste(unique(coef_df$model), collapse = ", "),
    call. = FALSE
  )
}

label_order <- c(
  "CSF amyloid-beta",
  "CSF total tau",
  "CSF phosphorylated tau",
  "AD EPOCH"
)

panel_a_df <- panel_a_df %>%
  mutate(label = factor(label, levels = rev(label_order))) %>%
  arrange(label)

x_min <- min(panel_a_df$hazard_ratio_ci_lower, na.rm = TRUE)
x_max <- max(panel_a_df$hazard_ratio_ci_upper, na.rm = TRUE)
annotation_x <- x_max * 1.10
plot_x_max <- x_max * 1.48

panel_a <- ggplot(panel_a_df, aes(x = hazard_ratio, y = label, color = label)) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.65,
    color = reference_color
  ) +
  geom_segment(
    aes(
      x = hazard_ratio_ci_lower,
      xend = hazard_ratio_ci_upper,
      y = label,
      yend = label
    ),
    linewidth = 1.05,
    lineend = "round"
  ) +
  geom_point(size = 3.5) +
  geom_text(
    aes(
      x = annotation_x,
      label = paste0("HR ", hr_label, "; ", p_label)
    ),
    hjust = 0,
    size = 3.4,
    color = "#222222",
    show.legend = FALSE
  ) +
  scale_color_manual(values = cb_palette, drop = FALSE) +
  scale_x_log10(
    limits = c(max(x_min * 0.82, 0.05), plot_x_max),
    breaks = c(0.5, 0.75, 1, 1.5, 2, 3),
    labels = label_number(accuracy = 0.01)
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "A",
    subtitle = "Biomarker hazard ratios in the full demographic + CSF + AD EPOCH model",
    x = "Hazard ratio per 1-SD increase",
    y = NULL
  ) +
  guides(color = "none") +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0),
    plot.subtitle = element_text(size = 10.5),
    axis.text.y = element_text(size = 10.5, color = "#222222"),
    axis.text.x = element_text(color = "#222222"),
    axis.title.x = element_text(color = "#222222"),
    plot.margin = margin(7, 100, 7, 7)
  )

# ----------------------------
# Panel B: cumulative CV C-index
# ----------------------------

required_perf_cols <- c("model_step", "model", "cv_cindex")
missing_perf_cols <- setdiff(required_perf_cols, names(perf_df))
if (length(missing_perf_cols) > 0) {
  stop(
    "Performance file is missing required columns: ",
    paste(missing_perf_cols, collapse = ", "),
    call. = FALSE
  )
}

if ("analysis_type" %in% names(perf_df)) {
  perf_df <- perf_df %>% filter(analysis_type == analysis_type_keep)
}
if ("status" %in% names(perf_df)) {
  perf_df <- perf_df %>% filter(status == "success")
}

panel_b_df <- perf_df %>%
  filter(model %in% c(
    "M0_age_sex",
    "M1_plus_amyloid",
    "M2_plus_tau",
    "M3_plus_AD_EPOCH"
  )) %>%
  arrange(model_step) %>%
  mutate(
    model_label = model_to_label(model),
    model_label = factor(model_label, levels = model_to_label(model))
  )

if (nrow(panel_b_df) != 4) {
  warning(
    "Expected 4 cumulative models but found ", nrow(panel_b_df),
    ". Available models: ", paste(unique(perf_df$model), collapse = ", ")
  )
}

if ("analysis_type" %in% names(comp_df)) {
  comp_df <- comp_df %>% filter(analysis_type == analysis_type_keep)
}

final_comparison <- comp_df %>%
  filter(
    new_model == final_model_name,
    previous_model == previous_model_name
  ) %>%
  slice(1)

if (nrow(final_comparison) == 0) {
  warning(
    "Could not find the M3 versus M2 comparison. The panel will be plotted without a p-value annotation."
  )
  comparison_label <- NULL
} else {
  delta_value <- final_comparison$delta_cindex[[1]]
  bootstrap_p <- final_comparison$bootstrap_p_two_sided[[1]]
  lr_p <- final_comparison$likelihood_ratio_p[[1]]
  comparison_label <- paste0(
    "AD EPOCH vs demographic + CSF model\n",
    "Delta C = ", formatC(delta_value, format = "f", digits = 3),
    "; bootstrap ", format_p(bootstrap_p),
    "; LR ", format_p(lr_p)
  )
}

cindex_min <- min(panel_b_df$cv_cindex, na.rm = TRUE)
cindex_max <- max(panel_b_df$cv_cindex, na.rm = TRUE)
y_lower <- max(0.45, floor((cindex_min - 0.025) * 100) / 100)
y_upper <- min(1.00, ceiling((cindex_max + 0.040) * 100) / 100)

panel_b <- ggplot(
  panel_b_df,
  aes(x = model_label, y = cv_cindex, group = 1)
) +
  geom_line(
    linewidth = 0.9,
    color = line_color
  ) +
  geom_point(
    aes(color = model_label),
    size = 3.8
  ) +
  geom_text(
    aes(label = formatC(cv_cindex, format = "f", digits = 3)),
    vjust = -1.0,
    size = 3.5,
    color = "#222222"
  ) +
  scale_color_manual(values = model_palette, drop = FALSE) +
  scale_y_continuous(
    limits = c(y_lower, y_upper),
    breaks = pretty_breaks(n = 5),
    labels = label_number(accuracy = 0.01)
  ) +
  labs(
    title = "B",
    subtitle = "Cumulative five-fold cross-validated discrimination",
    x = NULL,
    y = "Cross-validated Harrell C-index"
  ) +
  guides(color = "none") +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0),
    plot.subtitle = element_text(size = 10.5),
    axis.text.x = element_text(
      angle = 24,
      hjust = 1,
      vjust = 1,
      color = "#222222"
    ),
    axis.text.y = element_text(color = "#222222"),
    axis.title.y = element_text(color = "#222222"),
    plot.margin = margin(7, 7, 7, 7)
  )

if (!is.null(comparison_label)) {
  final_x <- which(levels(panel_b_df$model_label) == "+ AD EPOCH")
  if (length(final_x) == 1) {
    panel_b <- panel_b +
      annotate(
        "label",
        x = final_x,
        y = y_upper - 0.004,
        label = comparison_label,
        hjust = 1,
        vjust = 1,
        size = 3.05,
        label.size = 0.25,
        color = "#222222",
        fill = annotation_fill
      )
  }
}

# ----------------------------
# Combine and export
# ----------------------------

combined <- panel_a + panel_b +
  plot_layout(widths = c(1.08, 1.0))

ggsave(
  paste0(output_stem, ".pdf"),
  combined,
  width = 13.2,
  height = 5.4,
  units = "in",
  device = cairo_pdf
)

ggsave(
  paste0(output_stem, ".svg"),
  combined,
  width = 13.2,
  height = 5.4,
  units = "in"
)

ggsave(
  paste0(output_stem, ".png"),
  combined,
  width = 13.2,
  height = 5.4,
  units = "in",
  dpi = 400
)

write_tsv(panel_a_df, paste0(output_stem, "_panelA_data.tsv"))
write_tsv(panel_b_df, paste0(output_stem, "_panelB_data.tsv"))

message("Finished.")
message("PDF: ", paste0(output_stem, ".pdf"))
message("SVG: ", paste0(output_stem, ".svg"))
message("PNG: ", paste0(output_stem, ".png"))