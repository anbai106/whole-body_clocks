#!/usr/bin/env Rscript

# ==============================================================================
# Forest plot: brain-proteomics EPOCH vs 53 underlying brain-enriched proteins
# for incident Alzheimer's disease onset (G309)
#
# Primary display:
#   Hazard ratio per 1-SD higher predictor, with 95% CI
#
# The survival models were already fitted in Python on a common sample.
# This script does NOT refit the Cox models.
#
# It also creates a companion plot showing the change in C-index relative to
# the covariate-only model, because HR magnitude alone is not a direct measure
# of predictive discrimination.
#
# Usage:
#   Rscript plot_G309_EPOCH_vs_brain_proteins.R
#
# Optional:
#   Rscript plot_G309_EPOCH_vs_brain_proteins.R input.tsv output_directory
# ==============================================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# Input / output
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

default_input <- paste0(
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/",
  "mortality_clock/SA/output_EPOCH_vs_underlying_proteins/G309/",
  "brain_proteomics_EPOCH_vs_underlying_proteins_G309.tsv"
)

input_tsv <- if (length(args) >= 1) args[1] else default_input

default_out_dir <- dirname(input_tsv)
out_dir <- if (length(args) >= 2) args[2] else default_out_dir

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

forest_pdf <- file.path(
  out_dir,
  "brain_proteomics_EPOCH_vs_underlying_proteins_G309_forest.pdf"
)

forest_png <- file.path(
  out_dir,
  "brain_proteomics_EPOCH_vs_underlying_proteins_G309_forest.png"
)

cindex_pdf <- file.path(
  out_dir,
  "brain_proteomics_EPOCH_vs_underlying_proteins_G309_delta_cindex.pdf"
)

cindex_png <- file.path(
  out_dir,
  "brain_proteomics_EPOCH_vs_underlying_proteins_G309_delta_cindex.png"
)

plot_data_tsv <- file.path(
  out_dir,
  "brain_proteomics_EPOCH_vs_underlying_proteins_G309_plot_data.tsv"
)

sig_protein_tsv <- file.path(
  out_dir,
  "brain_proteomics_EPOCH_vs_underlying_proteins_G309_FDR_significant_proteins.tsv"
)

if (!file.exists(input_tsv)) {
  stop("Input TSV not found: ", input_tsv)
}

# ------------------------------------------------------------------------------
# Read results
# ------------------------------------------------------------------------------

dat <- read.delim(
  input_tsv,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  na.strings = c("NA", "NaN", ""),
  check.names = FALSE
)

required_cols <- c(
  "predictor_type",
  "predictor",
  "status",
  "N",
  "N_case",
  "N_noncase",
  "beta_per_1SD",
  "hr_per_1SD",
  "ci95_lower",
  "ci95_upper",
  "p_value",
  "p_fdr_bh_all_predictors",
  "p_fdr_bh_proteins_only",
  "base_cindex",
  "predictor_cindex",
  "delta_cindex_vs_base"
)

missing_cols <- setdiff(required_cols, names(dat))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns in input TSV: ",
    paste(missing_cols, collapse = ", ")
  )
}

dat <- dat[dat$status == "ok", , drop = FALSE]

if (nrow(dat) == 0) {
  stop("No rows with status == 'ok'.")
}

numeric_cols <- c(
  "N",
  "N_case",
  "N_noncase",
  "beta_per_1SD",
  "hr_per_1SD",
  "ci95_lower",
  "ci95_upper",
  "p_value",
  "p_fdr_bh_all_predictors",
  "p_fdr_bh_proteins_only",
  "base_cindex",
  "predictor_cindex",
  "delta_cindex_vs_base"
)

for (x in numeric_cols) {
  dat[[x]] <- suppressWarnings(as.numeric(dat[[x]]))
}

# ------------------------------------------------------------------------------
# Labels and significance
# ------------------------------------------------------------------------------

epoch_name <- "Brain proteomic EPOCH"

dat$display_name <- dat$predictor

dat$display_name[dat$predictor_type == "EPOCH"] <- epoch_name

# For the 53 individual proteins, use the protein-only BH correction.
# EPOCH is a prespecified benchmark and is therefore shown separately.
dat$protein_FDR_significant <- (
  dat$predictor_type == "protein" &
  !is.na(dat$p_fdr_bh_proteins_only) &
  dat$p_fdr_bh_proteins_only < 0.05
)

dat$plot_group <- "Other protein"
dat$plot_group[dat$protein_FDR_significant] <- "FDR-significant protein"
dat$plot_group[dat$predictor_type == "EPOCH"] <- "EPOCH"

dat$plot_group <- factor(
  dat$plot_group,
  levels = c("EPOCH", "FDR-significant protein", "Other protein")
)

# Make a compact significance label for optional inspection.
dat$FDR_label <- ""
dat$FDR_label[dat$protein_FDR_significant] <- "BH-FDR < 0.05"
dat$FDR_label[dat$predictor_type == "EPOCH"] <- "Prespecified EPOCH"

# ------------------------------------------------------------------------------
# Ordering for forest plot
# ------------------------------------------------------------------------------

# EPOCH first, followed by proteins ordered from largest to smallest Cox beta.
epoch_rows <- dat[dat$predictor_type == "EPOCH", , drop = FALSE]
protein_rows <- dat[dat$predictor_type == "protein", , drop = FALSE]

protein_rows <- protein_rows[
  order(protein_rows$beta_per_1SD, decreasing = TRUE, na.last = TRUE),
  ,
  drop = FALSE
]

plot_dat <- rbind(epoch_rows, protein_rows)

top_to_bottom <- plot_dat$display_name

# ggplot places the first factor level at the bottom, hence rev().
plot_dat$display_name <- factor(
  plot_dat$display_name,
  levels = rev(top_to_bottom)
)

# ------------------------------------------------------------------------------
# Summary table for significant proteins
# ------------------------------------------------------------------------------

sig_proteins <- plot_dat[
  plot_dat$predictor_type == "protein" &
  plot_dat$protein_FDR_significant,
  ,
  drop = FALSE
]

sig_proteins <- sig_proteins[
  order(sig_proteins$p_fdr_bh_proteins_only),
  ,
  drop = FALSE
]

write.table(
  sig_proteins,
  sig_protein_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

write.table(
  plot_dat,
  plot_data_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

# ------------------------------------------------------------------------------
# Forest plot
# ------------------------------------------------------------------------------

group_colors <- c(
  "EPOCH" = "#B2182B",
  "FDR-significant protein" = "#2166AC",
  "Other protein" = "#BDBDBD"
)

group_shapes <- c(
  "EPOCH" = 18,
  "FDR-significant protein" = 16,
  "Other protein" = 16
)

n_total <- unique(plot_dat$N[!is.na(plot_dat$N)])
n_case <- unique(plot_dat$N_case[!is.na(plot_dat$N_case)])

n_total_text <- if (length(n_total) == 1) {
  format(n_total, big.mark = ",", scientific = FALSE)
} else {
  "common sample"
}

n_case_text <- if (length(n_case) == 1) {
  format(n_case, big.mark = ",", scientific = FALSE)
} else {
  "incident cases"
}

forest_subtitle <- paste0(
  "Incident G309; N = ", n_total_text,
  " (", n_case_text, " cases). ",
  "Adjusted for age, sex, BMI and smoking."
)

p_forest <- ggplot(
  plot_dat,
  aes(
    y = display_name,
    color = plot_group
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey45"
  ) +
  geom_segment(
    aes(
      x = ci95_lower,
      xend = ci95_upper,
      yend = display_name
    ),
    linewidth = 0.6
  ) +
  geom_point(
    aes(
      x = hr_per_1SD,
      shape = plot_group
    ),
    size = 2.7
  ) +
  scale_x_log10(
    breaks = c(0.7, 0.8, 1.0, 1.25, 1.5, 2.0),
    labels = c("0.70", "0.80", "1.00", "1.25", "1.50", "2.00"),
    expand = expansion(mult = c(0.04, 0.08))
  ) +
  scale_color_manual(
    values = group_colors,
    drop = FALSE
  ) +
  scale_shape_manual(
    values = group_shapes,
    drop = FALSE
  ) +
  labs(
    title = "Brain proteomic EPOCH versus underlying brain-enriched proteins",
    subtitle = forest_subtitle,
    x = "Hazard ratio per 1-SD higher predictor",
    y = NULL,
    color = NULL,
    shape = NULL,
    caption = paste0(
      "Protein significance: Benjamini-Hochberg FDR < 0.05 across 53 proteins. ",
      "Hazard ratios are from separate covariate-adjusted Cox models."
    )
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9),
    axis.text.y = element_text(
      size = 7.3,
      color = "black"
    ),
    axis.text.x = element_text(color = "black"),
    axis.title.x = element_text(size = 10),
    legend.position = "top",
    legend.justification = "left",
    legend.box = "horizontal",
    legend.text = element_text(size = 8.5),
    plot.caption = element_text(
      size = 7.5,
      hjust = 0
    ),
    plot.margin = margin(
      t = 8,
      r = 12,
      b = 8,
      l = 8
    )
  )

ggsave(
  forest_pdf,
  p_forest,
  width = 8.0,
  height = 12.5,
  units = "in"
)

ggsave(
  forest_png,
  p_forest,
  width = 8.0,
  height = 12.5,
  units = "in",
  dpi = 400
)

# ------------------------------------------------------------------------------
# Companion plot: incremental C-index
# ------------------------------------------------------------------------------

# This is useful because HR magnitude is not equivalent to predictive power.
# Keep the same row order as the forest plot.

p_cindex <- ggplot(
  plot_dat,
  aes(
    y = display_name,
    x = delta_cindex_vs_base,
    color = plot_group
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey45"
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = delta_cindex_vs_base,
      yend = display_name
    ),
    linewidth = 0.6
  ) +
  geom_point(
    aes(shape = plot_group),
    size = 2.7
  ) +
  scale_color_manual(
    values = group_colors,
    drop = FALSE
  ) +
  scale_shape_manual(
    values = group_shapes,
    drop = FALSE
  ) +
  labs(
    title = "Incremental discrimination beyond conventional covariates",
    subtitle = forest_subtitle,
    x = expression(Delta * " C-index versus covariate-only model"),
    y = NULL,
    color = NULL,
    shape = NULL,
    caption = paste0(
      "Covariate-only model includes age, sex, BMI and smoking. ",
      "Positive values indicate improved discrimination."
    )
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9),
    axis.text.y = element_text(
      size = 7.3,
      color = "black"
    ),
    axis.text.x = element_text(color = "black"),
    axis.title.x = element_text(size = 10),
    legend.position = "top",
    legend.justification = "left",
    legend.box = "horizontal",
    legend.text = element_text(size = 8.5),
    plot.caption = element_text(
      size = 7.5,
      hjust = 0
    ),
    plot.margin = margin(
      t = 8,
      r = 12,
      b = 8,
      l = 8
    )
  )

ggsave(
  cindex_pdf,
  p_cindex,
  width = 8.0,
  height = 12.5,
  units = "in"
)

ggsave(
  cindex_png,
  p_cindex,
  width = 8.0,
  height = 12.5,
  units = "in",
  dpi = 400
)

# ------------------------------------------------------------------------------
# Console summary
# ------------------------------------------------------------------------------

epoch <- plot_dat[plot_dat$predictor_type == "EPOCH", , drop = FALSE]

cat("\n")
cat("====================================================================\n")
cat("G309 EPOCH vs underlying brain-enriched proteins\n")
cat("====================================================================\n")

if (nrow(epoch) == 1) {
  cat(
    sprintf(
      "EPOCH: HR = %.3f (95%% CI %.3f-%.3f), P = %.3g\n",
      epoch$hr_per_1SD,
      epoch$ci95_lower,
      epoch$ci95_upper,
      epoch$p_value
    )
  )
  cat(
    sprintf(
      "EPOCH: C-index = %.4f; delta C-index = %+.4f\n",
      epoch$predictor_cindex,
      epoch$delta_cindex_vs_base
    )
  )
}

cat(
  sprintf(
    "FDR-significant proteins: %d / %d\n",
    nrow(sig_proteins),
    nrow(protein_rows)
  )
)

if (nrow(protein_rows) > 0) {
  best_cindex_idx <- which.max(protein_rows$predictor_cindex)
  best <- protein_rows[best_cindex_idx, , drop = FALSE]

  cat(
    sprintf(
      "Best individual protein by C-index: %s, C-index = %.4f; delta C-index = %+.4f\n",
      best$predictor,
      best$predictor_cindex,
      best$delta_cindex_vs_base
    )
  )

  if (nrow(epoch) == 1) {
    n_below_epoch <- sum(
      protein_rows$predictor_cindex < epoch$predictor_cindex,
      na.rm = TRUE
    )

    cat(
      sprintf(
        "Proteins with lower C-index than EPOCH: %d / %d\n",
        n_below_epoch,
        nrow(protein_rows)
      )
    )
  }
}

cat("\nOutputs:\n")
cat("  Forest PDF:       ", forest_pdf, "\n", sep = "")
cat("  Forest PNG:       ", forest_png, "\n", sep = "")
cat("  Delta-C PDF:      ", cindex_pdf, "\n", sep = "")
cat("  Delta-C PNG:      ", cindex_png, "\n", sep = "")
cat("  Plot data TSV:    ", plot_data_tsv, "\n", sep = "")
cat("  Significant TSV:  ", sig_protein_tsv, "\n", sep = "")
cat("====================================================================\n")