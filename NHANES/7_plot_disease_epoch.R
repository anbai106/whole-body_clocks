#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# ============================================================
# STEP III forest plot:
# NHANES mortality EPOCH acceleration vs disease status
#
# Input:
#   disease_association_summary_all.tsv
# or:
#   logistic_results_healthy_control.tsv
#   logistic_results_disease_negative.tsv
#
# Main output:
#   NHANES_EPOCH_disease_forest_combined.pdf/png ## p_healthy main result
#   NHANES_EPOCH_disease_forest_primary_healthy_control.pdf/png
#   NHANES_EPOCH_disease_forest_sensitivity_disease_negative.pdf/png
#   forest_plot_results_with_bonferroni.tsv
# ============================================================

# -----------------------------
# Parse command-line arguments
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx == length(args)) return(default)
  args[idx + 1]
}

INDIR <- get_arg(
  "--indir",
  "/Users/hao/Dropbox/NHANES/output_dir/model2_nondisease_mortality_epoch/step3_disease_associations"
)

OUTDIR <- get_arg(
  "--outdir",
  file.path(INDIR, "figures")
)

BASE_FAMILY <- get_arg("--base_family", "Times New Roman")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

message("Input directory:  ", INDIR)
message("Output directory: ", OUTDIR)

# -----------------------------
# Colors and labels
# -----------------------------
category_cols <- c(
  cardiovascular = "#D55E00",
  cardiometabolic = "#E69F00",
  metabolic = "#F39C12",
  pulmonary = "#3E95B5",
  brain_vascular = "#0072B2",
  musculoskeletal = "#009E73",
  cancer = "#CC79A7",
  renal = "#6A3D9A",
  hepatic = "#8C564B",
  endocrine = "#7F7F7F",
  global = "#4D4D4D",
  other = "#999999"
)

category_order <- c(
  "cardiovascular",
  "cardiometabolic",
  "metabolic",
  "brain_vascular",
  "pulmonary",
  "renal",
  "hepatic",
  "endocrine",
  "musculoskeletal",
  "cancer",
  "global",
  "other"
)

disease_label_map <- c(
  diabetes = "Diabetes",
  hypertension = "Hypertension",
  high_cholesterol = "High cholesterol",
  asthma_ever = "Asthma ever",
  asthma_current = "Current asthma",
  arthritis = "Arthritis",
  congestive_heart_failure = "Congestive heart failure",
  coronary_heart_disease = "Coronary heart disease",
  angina = "Angina",
  myocardial_infarction = "Myocardial infarction",
  stroke = "Stroke",
  emphysema = "Emphysema",
  chronic_bronchitis = "Chronic bronchitis",
  copd_composite = "COPD composite",
  liver_condition = "Liver condition",
  thyroid_condition = "Thyroid condition",
  cancer_any = "Cancer",
  kidney_disease = "Kidney disease",
  cardiovascular_disease_composite = "CVD composite",
  heart_disease_composite = "Heart disease composite",
  major_disease_composite = "Major disease composite"
)

mode_label_map <- c(
  healthy_control = "Primary: disease cases vs strict healthy controls",
  disease_negative = "Sensitivity: disease cases vs disease-negative controls"
)

theme_epoch_forest <- function(base_size = 11) {
  theme_bw(base_size = base_size, base_family = BASE_FAMILY) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.border = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 0.35),
      axis.line.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "grey88", linewidth = 0.35),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(face = "bold", size = base_size + 2),
      plot.subtitle = element_text(size = base_size - 1, color = "#2C3E50"),
      legend.title = element_text(face = "bold"),
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey95", color = "grey80"),
      strip.text = element_text(face = "bold"),
      plot.margin = margin(8, 24, 8, 8)
    )
}

save_plot_pdf_png <- function(plot, basename, width, height) {
  pdf_file <- file.path(OUTDIR, paste0(basename, ".pdf"))
  png_file <- file.path(OUTDIR, paste0(basename, ".png"))

  pdf_device <- if (capabilities("cairo")) cairo_pdf else pdf

  ggsave(
    filename = pdf_file,
    plot = plot,
    width = width,
    height = height,
    device = pdf_device
  )

  ggsave(
    filename = png_file,
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )

  message("Saved: ", pdf_file)
  message("Saved: ", png_file)
}

format_p <- function(p) {
  out <- rep(NA_character_, length(p))
  out[is.na(p)] <- "NA"
  out[!is.na(p) & p == 0] <- "<1e-300"
  out[!is.na(p) & p > 0 & p < 1e-3] <- sprintf("%.1e", p[!is.na(p) & p > 0 & p < 1e-3])
  out[!is.na(p) & p >= 1e-3] <- sprintf("%.3f", p[!is.na(p) & p >= 1e-3])
  out
}

clean_category <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "other"
  x[!x %in% names(category_cols)] <- "other"
  x
}

# -----------------------------
# Read input files
# -----------------------------
summary_file <- file.path(INDIR, "disease_association_summary_all.tsv")
healthy_file <- file.path(INDIR, "logistic_results_healthy_control.tsv")
negative_file <- file.path(INDIR, "logistic_results_disease_negative.tsv")

if (file.exists(summary_file)) {
  dt <- fread(summary_file)
} else {
  if (!file.exists(healthy_file)) stop("Missing file: ", healthy_file)
  if (!file.exists(negative_file)) stop("Missing file: ", negative_file)

  healthy <- fread(healthy_file)
  negative <- fread(negative_file)
  dt <- rbindlist(list(healthy, negative), fill = TRUE)
}

required_cols <- c(
  "disease",
  "control_mode",
  "n",
  "n_case",
  "n_control",
  "or_per_1sd",
  "or_lower_95",
  "or_upper_95",
  "p"
)

missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# -----------------------------
# Clean and recompute Bonferroni correction
# -----------------------------
num_cols <- c(
  "n", "n_case", "n_control",
  "or_per_1sd", "or_lower_95", "or_upper_95",
  "p"
)

for (v in intersect(num_cols, names(dt))) {
  dt[, (v) := suppressWarnings(as.numeric(get(v)))]
}

if (!"category" %in% names(dt)) {
  dt[, category := "other"]
}

if (!"description" %in% names(dt)) {
  dt[, description := disease]
}

dt[, category := clean_category(category)]

dt <- dt[
  control_mode %in% c("healthy_control", "disease_negative") &
    is.finite(or_per_1sd) &
    is.finite(or_lower_95) &
    is.finite(or_upper_95) &
    or_per_1sd > 0 &
    or_lower_95 > 0 &
    or_upper_95 > 0 &
    !is.na(p)
]

# Exclude global composite from the main forest by default because it is not
# disease-specific. It remains saved if present in the raw result table.
dt[, is_global_composite := disease %in% c("major_disease_composite")]
plot_dt <- dt[is_global_composite == FALSE]

# Bonferroni correction based on number of disease endpoints tested within each comparison mode.
plot_dt[, n_diseases_tested := uniqueN(disease), by = control_mode]
plot_dt[, bonferroni_threshold := 0.05 / n_diseases_tested]
plot_dt[, p_bonferroni := pmin(p * n_diseases_tested, 1)]
plot_dt[, bonferroni_significant := p < bonferroni_threshold]

plot_dt[, p_label := format_p(p)]
plot_dt[, p_bonferroni_label := format_p(p_bonferroni)]
plot_dt[, or_label := sprintf("%.2f [%.2f, %.2f]", or_per_1sd, or_lower_95, or_upper_95)]
plot_dt[, n_label := paste0("N=", scales::comma(n), "; cases=", scales::comma(n_case))]
plot_dt[, sig_label := ifelse(bonferroni_significant, "Bonferroni significant", "Nominal / not significant")]

plot_dt[, disease_label := ifelse(
  disease %in% names(disease_label_map),
  disease_label_map[disease],
  gsub("_", " ", disease)
)]

plot_dt[, category := factor(category, levels = category_order)]
plot_dt[, control_mode_label := factor(
  mode_label_map[control_mode],
  levels = mode_label_map[c("healthy_control", "disease_negative")]
)]

# Disease order is defined from the primary healthy-control analysis.
primary_order <- plot_dt[
  control_mode == "healthy_control"
][
  order(category, -or_per_1sd, disease_label),
  disease
]

remaining_order <- setdiff(unique(plot_dt$disease), primary_order)
disease_order <- c(primary_order, remaining_order)

label_order <- unique(plot_dt[
  match(disease, disease_order)
][order(match(disease, disease_order)), disease_label])

plot_dt[, disease_label_factor := factor(disease_label, levels = rev(label_order))]

# Save table with Bonferroni correction.
fwrite(
  plot_dt,
  file.path(OUTDIR, "forest_plot_results_with_bonferroni.tsv"),
  sep = "\t"
)

# -----------------------------
# Forest plot function
# -----------------------------
plot_forest_one_mode <- function(d, title, subtitle, show_right_text = TRUE) {
  d <- copy(d)

  if (nrow(d) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, label = "No valid rows", size = 5) +
        theme_void() +
        labs(title = title)
    )
  }

  x_min <- min(d$or_lower_95, na.rm = TRUE)
  x_max <- max(d$or_upper_95, na.rm = TRUE)

  x_min <- max(0.5, floor(x_min * 10) / 10)
  x_max <- min(12, x_max * 1.25)
  if (x_max < 5) x_max <- 5

  text_x <- x_max * 1.05
  d[, text_x := text_x]

  d_sig <- d[bonferroni_significant == TRUE]
  d_nonsig <- d[bonferroni_significant == FALSE]

  p <- ggplot(d, aes(y = disease_label_factor)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey45", linewidth = 0.45) +
    geom_errorbarh(
      aes(xmin = or_lower_95, xmax = or_upper_95, color = category),
      height = 0.18,
      linewidth = 0.65,
      alpha = 0.95
    ) +
    geom_point(
      data = d_nonsig,
      aes(x = or_per_1sd, color = category),
      shape = 1,
      size = 2.8,
      stroke = 0.9
    ) +
    geom_point(
      data = d_sig,
      aes(x = or_per_1sd, fill = category),
      shape = 21,
      size = 3.2,
      color = "black",
      stroke = 0.25
    ) +
    scale_color_manual(values = category_cols, drop = FALSE) +
    scale_fill_manual(values = category_cols, drop = FALSE) +
    scale_x_log10(
      breaks = c(0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 7.5, 10),
      labels = c("0.5", "0.75", "1", "1.5", "2", "3", "4", "5", "7.5", "10")
    ) +
    coord_cartesian(xlim = c(x_min, x_max), clip = "off") +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Odds ratio per 1-SD mortality EPOCH acceleration, log scale",
      y = NULL,
      color = "Disease domain",
      fill = "Disease domain"
    ) +
    theme_epoch_forest(base_size = 11)

  if (show_right_text) {
    p <- p +
      geom_text(
        aes(x = text_x, label = or_label),
        hjust = 0,
        size = 3.0,
        family = BASE_FAMILY
      ) +
      annotate(
        "text",
        x = text_x,
        y = length(levels(d$disease_label_factor)) + 0.8,
        label = "OR [95% CI]",
        hjust = 0,
        fontface = "bold",
        size = 3.1,
        family = BASE_FAMILY
      ) +
      theme(plot.margin = margin(8, 90, 8, 8))
  }

  p
}

# -----------------------------
# Individual plots
# -----------------------------
healthy_dt <- plot_dt[control_mode == "healthy_control"]
negative_dt <- plot_dt[control_mode == "disease_negative"]

n_healthy_tests <- unique(healthy_dt$n_diseases_tested)
n_negative_tests <- unique(negative_dt$n_diseases_tested)

healthy_subtitle <- paste0(
  "Cases vs strict healthy controls; Bonferroni threshold = 0.05 / ",
  n_healthy_tests,
  " = ",
  signif(0.05 / n_healthy_tests, 3)
)

negative_subtitle <- paste0(
  "Cases vs disease-negative controls; Bonferroni threshold = 0.05 / ",
  n_negative_tests,
  " = ",
  signif(0.05 / n_negative_tests, 3)
)

p_healthy <- plot_forest_one_mode(
  healthy_dt,
  title = "Mortality EPOCH acceleration is associated with baseline disease burden",
  subtitle = healthy_subtitle,
  show_right_text = TRUE
)

p_negative <- plot_forest_one_mode(
  negative_dt,
  title = "Sensitivity analysis using disease-negative controls",
  subtitle = negative_subtitle,
  show_right_text = TRUE
)

# -----------------------------
# Combined two-panel figure
# -----------------------------
combined_fig <- p_healthy / p_negative +
  plot_annotation(
    title = "NHANES mortality EPOCH acceleration and disease status",
    subtitle = paste0(
      "Odds ratios are per 1-SD increase in mortality_epoch_acceleration_z. ",
      "Filled points indicate Bonferroni-significant associations."
    ),
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16, family = BASE_FAMILY),
      plot.subtitle = element_text(size = 11, family = BASE_FAMILY)
    )
  )

# -----------------------------
# Save figures
# -----------------------------
save_plot_pdf_png(
  p_healthy,
  "NHANES_EPOCH_disease_forest_primary_healthy_control",
  width = 11,
  height = 8.5
)

save_plot_pdf_png(
  p_negative,
  "NHANES_EPOCH_disease_forest_sensitivity_disease_negative",
  width = 11,
  height = 8.5
)

save_plot_pdf_png(
  combined_fig,
  "NHANES_EPOCH_disease_forest_combined",
  width = 11,
  height = 16
)

# -----------------------------
# Compact significant table
# -----------------------------
sig_table <- plot_dt[
  bonferroni_significant == TRUE
][
  order(control_mode, p, disease)
][
  ,
  .(
    control_mode,
    disease,
    disease_label,
    category,
    n,
    n_case,
    n_control,
    or_per_1sd,
    or_lower_95,
    or_upper_95,
    p,
    p_bonferroni,
    bonferroni_threshold
  )
]

fwrite(
  sig_table,
  file.path(OUTDIR, "bonferroni_significant_disease_associations.tsv"),
  sep = "\t"
)

# -----------------------------
# Console summary
# -----------------------------
message("\nDone.")
message("Bonferroni correction was computed within each comparison mode:")
message("  healthy_control:   0.05 / ", n_healthy_tests, " = ", signif(0.05 / n_healthy_tests, 4))
message("  disease_negative:  0.05 / ", n_negative_tests, " = ", signif(0.05 / n_negative_tests, 4))

message("\nMain outputs:")
message("  ", file.path(OUTDIR, "NHANES_EPOCH_disease_forest_combined.pdf"))
message("  ", file.path(OUTDIR, "NHANES_EPOCH_disease_forest_combined.png"))
message("  ", file.path(OUTDIR, "NHANES_EPOCH_disease_forest_primary_healthy_control.pdf"))
message("  ", file.path(OUTDIR, "NHANES_EPOCH_disease_forest_sensitivity_disease_negative.pdf"))
message("  ", file.path(OUTDIR, "forest_plot_results_with_bonferroni.tsv"))
message("  ", file.path(OUTDIR, "bonferroni_significant_disease_associations.tsv"))