#!/usr/bin/env Rscript

# =============================================================================
# Brain proteomics BAG-EPOCH residual resilience figure
#
# Panel A:
#   Relationship between standardized Brain_ProtBAG and brain proteomics
#   mortality EPOCH acceleration. The solid line is the fitted OLS expectation:
#
#       EPOCH_z = alpha + beta * BAG_z
#
#   Vertical dashed lines mark low/high BAG thresholds, and diagonal dashed
#   lines mark low/high EPOCH|BAG residual thresholds. The four p20 phenotypes
#   are highlighted.
#
# Panel B:
#   Residualized coordinate system: BAG_z versus standardized EPOCH|BAG
#   discordance residual. This directly visualizes the four phenotypes:
#
#       CFA: low BAG  + low residual
#       CRA: high BAG + low residual
#       LVA: low BAG  + high residual
#       CUA: high BAG + high residual
#
# Default threshold is p20. Set THRESHOLD_TAG=p10 or p25 to reproduce the
# alternative definitions already present in the participant-level TSV.
#
# Expected input columns:
#   participant_id
#   Brain_ProtBAG
#   Brain_ProtBAG_z
#   brain_proteomics_mortality_clock_acceleration_z
#   brain_proteomics_EPOCH_z_predicted_from_BAG
#   EPOCH_BAG_discordance_residual
#   EPOCH_BAG_discordance_residual_z
#   CFA_p10 / CRA_p10 / LVA_p10 / CUA_p10 / aging_phenotype_p10
#   CFA_p20 / CRA_p20 / LVA_p20 / CUA_p20 / aging_phenotype_p20
#   CFA_p25 / CRA_p25 / LVA_p25 / CUA_p25 / aging_phenotype_p25
#
# Outputs:
#   Panel_A_BAG_vs_EPOCH_<threshold>.pdf/png/svg
#   Panel_B_BAG_vs_residual_<threshold>.pdf/png/svg
#   Figure_BAG_EPOCH_resilience_AB_<threshold>.pdf/png/svg
#
# =============================================================================

options(stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# User settings
# -----------------------------------------------------------------------------
ROOT <- Sys.getenv(
  "ROOT",
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/Brain_proteomics_mortality_clock"
)

INPUT_TSV <- Sys.getenv(
  "INPUT_TSV",
  file.path(ROOT, "brain_proteomics_BAG_EPOCH_discordance_resilience.tsv")
)

SUMMARY_TSV <- Sys.getenv(
  "SUMMARY_TSV",
  file.path(ROOT, "brain_proteomics_BAG_EPOCH_discordance_summary.tsv")
)

OUT_DIR <- Sys.getenv(
  "OUT_DIR",
  file.path(ROOT, "figures_BAG_EPOCH_resilience")
)

# Primary categorical definition. Alternatives: p10, p25.
THRESHOLD_TAG <- Sys.getenv("THRESHOLD_TAG", "p20")

# Plot dimensions.
PANEL_WIDTH_IN <- 5.2
PANEL_HEIGHT_IN <- 4.5
COMBINED_WIDTH_IN <- 10.4
COMBINED_HEIGHT_IN <- 4.6
PNG_DPI <- 400

# Point settings. To keep PDF/SVG files lightweight and easy to edit in
# Inkscape, only a deterministic random subsample of non-extreme (grey)
# participants is plotted. All CFA/CRA/LVA/CUA participants are retained.
BACKGROUND_N <- as.integer(Sys.getenv("BACKGROUND_N", "1500"))
BACKGROUND_SEED <- as.integer(Sys.getenv("BACKGROUND_SEED", "20260814"))
BACKGROUND_ALPHA <- 0.18
BACKGROUND_SIZE <- 0.60
PHENOTYPE_ALPHA <- 0.72
PHENOTYPE_SIZE <- 1.00

if (!is.finite(BACKGROUND_N) || BACKGROUND_N < 0) {
  stop("BACKGROUND_N must be a non-negative integer.")
}
if (!is.finite(BACKGROUND_SEED)) {
  stop("BACKGROUND_SEED must be a finite integer.")
}

# -----------------------------------------------------------------------------
# Packages
# -----------------------------------------------------------------------------
required_packages <- c("ggplot2", "patchwork")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nInstall with: install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))"
  )
}

library(ggplot2)
library(patchwork)

# -----------------------------------------------------------------------------
# Constants / expected column names
# -----------------------------------------------------------------------------
ID_COL <- "participant_id"
BAG_RAW_COL <- "Brain_ProtBAG"
BAG_Z_COL <- "Brain_ProtBAG_z"
EPOCH_Z_COL <- "brain_proteomics_mortality_clock_acceleration_z"
EPOCH_PRED_COL <- "brain_proteomics_EPOCH_z_predicted_from_BAG"
RESID_COL <- "EPOCH_BAG_discordance_residual"
RESID_Z_COL <- "EPOCH_BAG_discordance_residual_z"

PHENO_SHORT <- c(
  "Concordant_favorable_ager"   = "CFA",
  "Candidate_resilient_ager"    = "CRA",
  "Latent_vulnerability_ager"   = "LVA",
  "Concordant_unfavorable_ager" = "CUA"
)

PHENO_LONG <- c(
  "CFA" = "Concordant favorable agers",
  "CRA" = "Candidate resilient agers",
  "LVA" = "Latent vulnerability agers",
  "CUA" = "Concordant unfavorable agers"
)

# Fixed, intuitive phenotype colors.
PHENO_COLORS <- c(
  "CFA" = "#2A9D8F",
  "CRA" = "#457B9D",
  "LVA" = "#E9A23B",
  "CUA" = "#D1495B"
)

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------
read_tsv_base <- function(path) {
  if (!file.exists(path)) {
    stop("Input file not found: ", path)
  }
  read.delim(
    path,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    na.strings = c("NA", "", "NaN", "nan"),
    quote = "",
    comment.char = ""
  )
}

as_logical_robust <- function(x) {
  if (is.logical(x)) return(x)
  z <- trimws(tolower(as.character(x)))
  z %in% c("true", "t", "1", "yes", "y")
}

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  cor(x[ok], y[ok], method = "pearson")
}

format_n <- function(x) {
  format(as.integer(x), big.mark = ",", scientific = FALSE)
}

theme_epoch <- function() {
  theme_classic(base_size = 10) +
    theme(
      text = element_text(family = "Arial", colour = "black"),
      axis.text = element_text(colour = "black", size = 9),
      axis.title = element_text(colour = "black", size = 10),
      plot.title = element_text(face = "bold", colour = "black", size = 11),
      plot.subtitle = element_text(colour = "black", size = 9),
      legend.title = element_blank(),
      legend.text = element_text(colour = "black", size = 8.5),
      legend.key.height = grid::unit(0.42, "cm"),
      plot.margin = margin(7, 8, 7, 7)
    )
}

save_plot_all <- function(plot, stem, width, height) {
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  
  ggsave(
    filename = file.path(OUT_DIR, paste0(stem, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf,
    bg = "white"
  )
  
  ggsave(
    filename = file.path(OUT_DIR, paste0(stem, ".png")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = PNG_DPI,
    bg = "white"
  )
  
  # SVG requires the svglite package. Skip gracefully if unavailable.
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(
      filename = file.path(OUT_DIR, paste0(stem, ".svg")),
      plot = plot,
      width = width,
      height = height,
      units = "in",
      device = svglite::svglite,
      bg = "white"
    )
  } else {
    message("Package 'svglite' not installed; SVG output skipped.")
  }
}

# -----------------------------------------------------------------------------
# Read participant-level data
# -----------------------------------------------------------------------------
message("Reading participant-level TSV: ", INPUT_TSV)
dat <- read_tsv_base(INPUT_TSV)

group_col <- paste0("aging_phenotype_", THRESHOLD_TAG)
cfa_col <- paste0("CFA_", THRESHOLD_TAG)
cra_col <- paste0("CRA_", THRESHOLD_TAG)
lva_col <- paste0("LVA_", THRESHOLD_TAG)
cua_col <- paste0("CUA_", THRESHOLD_TAG)

required_cols <- c(
  ID_COL,
  BAG_Z_COL,
  EPOCH_Z_COL,
  RESID_COL,
  RESID_Z_COL,
  cfa_col,
  cra_col,
  lva_col,
  cua_col
)

missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop(
    "Participant TSV is missing required column(s): ",
    paste(missing_cols, collapse = ", ")
  )
}

# Numeric cleanup.
for (cc in c(BAG_Z_COL, EPOCH_Z_COL, RESID_COL, RESID_Z_COL)) {
  dat[[cc]] <- suppressWarnings(as.numeric(dat[[cc]]))
}

# Complete plotting population.
keep <- is.finite(dat[[BAG_Z_COL]]) &
  is.finite(dat[[EPOCH_Z_COL]]) &
  is.finite(dat[[RESID_COL]]) &
  is.finite(dat[[RESID_Z_COL]])
dat <- dat[keep, , drop = FALSE]

if (nrow(dat) < 20) {
  stop("Too few complete participants for plotting.")
}

# -----------------------------------------------------------------------------
# Derive the four-phenotype label directly from Boolean columns.
# This remains robust even if the group column name/labels change.
# -----------------------------------------------------------------------------
dat[[cfa_col]] <- as_logical_robust(dat[[cfa_col]])
dat[[cra_col]] <- as_logical_robust(dat[[cra_col]])
dat[[lva_col]] <- as_logical_robust(dat[[lva_col]])
dat[[cua_col]] <- as_logical_robust(dat[[cua_col]])

dat$phenotype4 <- "Other"
dat$phenotype4[dat[[cfa_col]]] <- "CFA"
dat$phenotype4[dat[[cra_col]]] <- "CRA"
dat$phenotype4[dat[[lva_col]]] <- "LVA"
dat$phenotype4[dat[[cua_col]]] <- "CUA"

dat$phenotype4 <- factor(
  dat$phenotype4,
  levels = c("Other", "CFA", "CRA", "LVA", "CUA")
)

extreme <- dat[dat$phenotype4 != "Other", , drop = FALSE]

# Display-only background subsample. Statistical quantities, thresholds, model
# estimates, and phenotype counts continue to use the full analysis population.
background_pool <- dat[dat$phenotype4 == "Other", , drop = FALSE]
set.seed(BACKGROUND_SEED)
if (BACKGROUND_N == 0L) {
  background <- background_pool[0, , drop = FALSE]
} else if (nrow(background_pool) > BACKGROUND_N) {
  background <- background_pool[
    sample.int(nrow(background_pool), size = BACKGROUND_N, replace = FALSE),
    ,
    drop = FALSE
  ]
} else {
  background <- background_pool
}

message(
  "Plotting ", format_n(nrow(background)), " grey background participants out of ",
  format_n(nrow(background_pool)), " non-extreme participants; all ",
  format_n(nrow(extreme)), " four-phenotype participants are retained."
)

# -----------------------------------------------------------------------------
# Recover threshold fraction from tag: p10 -> 0.10, p20 -> 0.20, p2p5 -> 0.025.
# -----------------------------------------------------------------------------
tag_to_fraction <- function(tag) {
  if (!grepl("^p", tag)) stop("THRESHOLD_TAG must look like p10, p20, p25, etc.")
  txt <- sub("^p", "", tag)
  txt <- sub("p", ".", txt, fixed = TRUE)
  as.numeric(txt) / 100
}

tail_fraction <- tag_to_fraction(THRESHOLD_TAG)
if (!is.finite(tail_fraction) || tail_fraction <= 0 || tail_fraction >= 0.5) {
  stop("Invalid THRESHOLD_TAG: ", THRESHOLD_TAG)
}

# Calculate the actual cutoffs from the participant-level data.
bag_low <- as.numeric(quantile(dat[[BAG_Z_COL]], probs = tail_fraction, na.rm = TRUE))
bag_high <- as.numeric(quantile(dat[[BAG_Z_COL]], probs = 1 - tail_fraction, na.rm = TRUE))
resid_z_low <- as.numeric(quantile(dat[[RESID_Z_COL]], probs = tail_fraction, na.rm = TRUE))
resid_z_high <- as.numeric(quantile(dat[[RESID_Z_COL]], probs = 1 - tail_fraction, na.rm = TRUE))
resid_raw_low <- as.numeric(quantile(dat[[RESID_COL]], probs = tail_fraction, na.rm = TRUE))
resid_raw_high <- as.numeric(quantile(dat[[RESID_COL]], probs = 1 - tail_fraction, na.rm = TRUE))

# -----------------------------------------------------------------------------
# Obtain OLS model parameters.
# Prefer summary TSV when available; otherwise reconstruct from participant data.
# -----------------------------------------------------------------------------
intercept <- NA_real_
beta_bag <- NA_real_
r_squared <- NA_real_
pearson_r <- safe_cor(dat[[BAG_Z_COL]], dat[[EPOCH_Z_COL]])
resid_bag_r <- safe_cor(dat[[BAG_Z_COL]], dat[[RESID_Z_COL]])

if (file.exists(SUMMARY_TSV)) {
  message("Reading summary TSV: ", SUMMARY_TSV)
  sm <- read_tsv_base(SUMMARY_TSV)
  if (nrow(sm) > 0) {
    if ("OLS_intercept" %in% names(sm)) {
      intercept <- suppressWarnings(as.numeric(sm$OLS_intercept[1]))
    }
    if ("OLS_beta_BAG_z" %in% names(sm)) {
      beta_bag <- suppressWarnings(as.numeric(sm$OLS_beta_BAG_z[1]))
    }
    if ("OLS_R_squared" %in% names(sm)) {
      r_squared <- suppressWarnings(as.numeric(sm$OLS_R_squared[1]))
    }
    if ("Pearson_r_BAGz_vs_EPOCHz" %in% names(sm)) {
      pearson_r <- suppressWarnings(as.numeric(sm$Pearson_r_BAGz_vs_EPOCHz[1]))
    }
    if ("Pearson_r_BAGz_vs_residual" %in% names(sm)) {
      resid_bag_r <- suppressWarnings(as.numeric(sm$Pearson_r_BAGz_vs_residual[1]))
    }
  }
}

# Fallback: refit OLS from the plotting population.
if (!is.finite(intercept) || !is.finite(beta_bag) || !is.finite(r_squared)) {
  fit <- lm(
    dat[[EPOCH_Z_COL]] ~ dat[[BAG_Z_COL]]
  )
  cf <- coef(fit)
  intercept <- unname(cf[1])
  beta_bag <- unname(cf[2])
  r_squared <- summary(fit)$r.squared
}

# Counts for legend.
counts <- table(factor(extreme$phenotype4, levels = c("CFA", "CRA", "LVA", "CUA")))
legend_labels <- c(
  "CFA" = paste0("CFA: Concordant favorable (n=", format_n(counts["CFA"]), ")"),
  "CRA" = paste0("CRA: Candidate resilient (n=", format_n(counts["CRA"]), ")"),
  "LVA" = paste0("LVA: Latent vulnerability (n=", format_n(counts["LVA"]), ")"),
  "CUA" = paste0("CUA: Concordant unfavorable (n=", format_n(counts["CUA"]), ")")
)

# -----------------------------------------------------------------------------
# PANEL A
# BAG versus EPOCH with OLS line and residual boundaries.
#
# A residual cutoff r corresponds to a line parallel to the fitted regression:
#     EPOCH = intercept + beta * BAG + r
# -----------------------------------------------------------------------------
panel_a <- ggplot(dat, aes(x = .data[[BAG_Z_COL]], y = .data[[EPOCH_Z_COL]])) +
  geom_point(
    data = background,
    colour = "grey55",
    alpha = BACKGROUND_ALPHA,
    size = BACKGROUND_SIZE,
    stroke = 0
  ) +
  geom_abline(
    intercept = intercept + resid_raw_low,
    slope = beta_bag,
    linetype = "dashed",
    linewidth = 0.55,
    colour = "grey35"
  ) +
  geom_abline(
    intercept = intercept + resid_raw_high,
    slope = beta_bag,
    linetype = "dashed",
    linewidth = 0.55,
    colour = "grey35"
  ) +
  geom_vline(
    xintercept = c(bag_low, bag_high),
    linetype = "dotted",
    linewidth = 0.55,
    colour = "grey35"
  ) +
  geom_point(
    data = extreme,
    aes(colour = phenotype4),
    alpha = PHENOTYPE_ALPHA,
    size = PHENOTYPE_SIZE,
    stroke = 0
  ) +
  geom_abline(
    intercept = intercept,
    slope = beta_bag,
    linewidth = 0.85,
    colour = "black"
  ) +
  scale_colour_manual(
    values = PHENO_COLORS,
    breaks = c("CFA", "CRA", "LVA", "CUA"),
    labels = legend_labels,
    drop = FALSE
  ) +
  labs(
    title = "BAG and mortality EPOCH capture overlapping but distinct variation",
    subtitle = paste0(
      "Solid line: EPOCH expected from BAG; dashed diagonal lines: ",
      THRESHOLD_TAG, " residual boundaries; grey points are a display subsample"
    ),
    x = "Brain proteomic BAG (z-score)",
    y = "Brain proteomic mortality EPOCH acceleration (z-score)"
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    hjust = -0.08,
    vjust = 1.25,
    size = 3.1,
    family = "Arial",
    colour = "black",
    label = sprintf(
      "N = %s\nr = %.3f\nR\u00B2 = %.3f",
      format_n(nrow(dat)),
      pearson_r,
      r_squared
    )
  ) +
  theme_epoch() +
  theme(legend.position = "bottom")

# -----------------------------------------------------------------------------
# PANEL B
# Residualized BAG-EPOCH plane with four interpretable corner phenotypes.
# -----------------------------------------------------------------------------
x_range <- range(dat[[BAG_Z_COL]], finite = TRUE)
y_range <- range(dat[[RESID_Z_COL]], finite = TRUE)
x_span <- diff(x_range)
y_span <- diff(y_range)

# Label coordinates placed inside the threshold-defined corners.
x_left <- bag_low - 0.06 * x_span
x_right <- bag_high + 0.06 * x_span
y_bottom <- resid_z_low - 0.05 * y_span
y_top <- resid_z_high + 0.05 * y_span

panel_b <- ggplot(dat, aes(x = .data[[BAG_Z_COL]], y = .data[[RESID_Z_COL]])) +
  geom_point(
    data = background,
    colour = "grey55",
    alpha = BACKGROUND_ALPHA,
    size = BACKGROUND_SIZE,
    stroke = 0
  ) +
  geom_vline(
    xintercept = c(bag_low, bag_high),
    linetype = "dashed",
    linewidth = 0.6,
    colour = "grey30"
  ) +
  geom_hline(
    yintercept = c(resid_z_low, resid_z_high),
    linetype = "dashed",
    linewidth = 0.6,
    colour = "grey30"
  ) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.35,
    colour = "grey65"
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.35,
    colour = "grey65"
  ) +
  geom_point(
    data = extreme,
    aes(colour = phenotype4),
    alpha = PHENOTYPE_ALPHA,
    size = PHENOTYPE_SIZE,
    stroke = 0
  ) +
  scale_colour_manual(
    values = PHENO_COLORS,
    breaks = c("CFA", "CRA", "LVA", "CUA"),
    labels = legend_labels,
    drop = FALSE
  ) +
  annotate(
    "label",
    x = x_left,
    y = y_bottom,
    label = paste0("CFA\nn=", format_n(counts["CFA"])),
    hjust = 1,
    vjust = 1,
    size = 2.75,
    family = "Arial",
    label.size = 0.18,
    fill = "white",
    colour = PHENO_COLORS["CFA"]
  ) +
  annotate(
    "label",
    x = x_right,
    y = y_bottom,
    label = paste0("CRA\nn=", format_n(counts["CRA"])),
    hjust = 0,
    vjust = 1,
    size = 2.75,
    family = "Arial",
    label.size = 0.18,
    fill = "white",
    colour = PHENO_COLORS["CRA"]
  ) +
  annotate(
    "label",
    x = x_left,
    y = y_top,
    label = paste0("LVA\nn=", format_n(counts["LVA"])),
    hjust = 1,
    vjust = 0,
    size = 2.75,
    family = "Arial",
    label.size = 0.18,
    fill = "white",
    colour = PHENO_COLORS["LVA"]
  ) +
  annotate(
    "label",
    x = x_right,
    y = y_top,
    label = paste0("CUA\nn=", format_n(counts["CUA"])),
    hjust = 0,
    vjust = 0,
    size = 2.75,
    family = "Arial",
    label.size = 0.18,
    fill = "white",
    colour = PHENO_COLORS["CUA"]
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    hjust = -0.08,
    vjust = 1.25,
    size = 3.1,
    family = "Arial",
    colour = "black",
    label = sprintf(
      "r(BAG, residual) = %.3g",
      resid_bag_r
    )
  ) +
  labs(
    title = "Residual discordance identifies four aging-vulnerability states",
    subtitle = paste0(
      "Primary definition: ", THRESHOLD_TAG,
      " tails; grey points are a display subsample of non-extreme participants"
    ),
    x = "Brain proteomic BAG (z-score)",
    y = "EPOCH\u2013BAG discordance residual (z-score)"
  ) +
  theme_epoch() +
  theme(legend.position = "bottom")

# -----------------------------------------------------------------------------
# Combined A+B figure
# -----------------------------------------------------------------------------
combined <- panel_a + panel_b +
  plot_layout(guides = "collect", widths = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "bottom",
    plot.tag = element_text(
      family = "Arial",
      face = "bold",
      size = 13,
      colour = "black"
    )
  )

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------
stem_a <- paste0("Panel_A_BAG_vs_EPOCH_", THRESHOLD_TAG)
stem_b <- paste0("Panel_B_BAG_vs_residual_", THRESHOLD_TAG)
stem_ab <- paste0("Figure_BAG_EPOCH_resilience_AB_", THRESHOLD_TAG)

save_plot_all(panel_a, stem_a, PANEL_WIDTH_IN, PANEL_HEIGHT_IN)
save_plot_all(panel_b, stem_b, PANEL_WIDTH_IN, PANEL_HEIGHT_IN)
save_plot_all(combined, stem_ab, COMBINED_WIDTH_IN, COMBINED_HEIGHT_IN)

# -----------------------------------------------------------------------------
# Console report
# -----------------------------------------------------------------------------
cat("\nResidual-based resilience figure complete\n")
cat("  Input N:", format_n(nrow(dat)), "\n")
cat("  Threshold:", THRESHOLD_TAG, "(", tail_fraction * 100, "% tails )\n")
cat(sprintf("  OLS: EPOCH_z = %.6f + %.6f * BAG_z + residual\n", intercept, beta_bag))
cat(sprintf("  Pearson r(BAG, EPOCH) = %.6f\n", pearson_r))
cat(sprintf("  R^2 = %.6f\n", r_squared))
cat(sprintf("  r(BAG, residual) = %.6g\n", resid_bag_r))
cat("  CFA:", format_n(counts["CFA"]), "\n")
cat("  CRA:", format_n(counts["CRA"]), "\n")
cat("  LVA:", format_n(counts["LVA"]), "\n")
cat("  CUA:", format_n(counts["CUA"]), "\n")
cat("  Grey background points plotted:", format_n(nrow(background)),
    "of", format_n(nrow(background_pool)), "non-extreme participants\n")
cat("  Background sampling seed:", BACKGROUND_SEED, "\n")
cat("  Output directory:", OUT_DIR, "\n")