#!/usr/bin/env Rscript

# ==============================================================================
# Publication-style plots for:
# Brain proteomics mortality EPOCH vs 10 conventional biomarkers
# predicting incident Alzheimer's disease (G309) and heart failure (I500)
#
# Default input root:
# /Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/
# Brain_proteomics_mortality_clock/
# comparative_with_10_conventional_biomarkers_disease_onset
#
# Expected files under all/G309 and all/I500:
#   *_individual_predictor_summary.tsv
#   *_model_panel_summary.tsv
#   *_paired_delta_cindex.tsv                    [optional but recommended]
#   *_epoch_beyond_10_biomarker_panel.tsv       [optional]
#
# Main 2x2 figure:
#   A. HR forest plot: Alzheimer's disease
#   B. HR forest plot: Heart failure
#   C. Model-level C-index comparison
#   D. Incremental C-index of EPOCH beyond the complete 10-biomarker panel
#
# Also saves a supplementary individual-predictor delta-C figure.
# ============================================================================== 

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)

ROOT <- if (length(args) >= 1) args[[1]] else
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/Brain_proteomics_mortality_clock/comparative_with_10_conventional_biomarkers_disease_onset"

ANALYSIS_SPLIT <- if (length(args) >= 2) args[[2]] else "all"

OUTDIR <- if (length(args) >= 3) args[[3]] else
  file.path(ROOT, "figures_brain_EPOCH_vs_10_biomarkers")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

message("Input root: ", ROOT)
message("Analysis split: ", ANALYSIS_SPLIT)
message("Output directory: ", OUTDIR)

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

find_one <- function(dir, pattern, required = TRUE) {
  if (!dir.exists(dir)) {
    if (required) stop("Directory not found: ", dir)
    return(NA_character_)
  }
  hits <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(hits) == 0) {
    if (required) stop("No file matching /", pattern, "/ in: ", dir)
    return(NA_character_)
  }
  if (length(hits) > 1) {
    warning("Multiple files matched /", pattern, "/ in ", dir,
            "; using: ", basename(hits[[1]]))
  }
  hits[[1]]
}

read_tsv_safe <- function(path) {
  readr::read_tsv(path, show_col_types = FALSE, na = c("NA", "NaN", ""))
}

fmt_p <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 0.001, formatC(p, format = "e", digits = 2),
                sprintf("%.3f", p)))
}

# Color-blind-friendly palette.
COL_EPOCH <- "#0072B2"
COL_BIO   <- "#666666"
COL_BASE  <- "#7F7F7F"
COL_PANEL <- "#D55E00"
COL_BOTH  <- "#009E73"
COL_GRID  <- "#D9D9D9"

base_theme <- theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 11.5),
    plot.subtitle = element_text(size = 9.5, colour = "grey30"),
    axis.title = element_text(size = 10.5),
    axis.text = element_text(size = 9.2, colour = "black"),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(7, 10, 7, 7)
  )

# ------------------------------------------------------------------------------
# Locate and read files
# ------------------------------------------------------------------------------

codes <- c("G309", "I500")
disease_labels <- c(G309 = "Alzheimer's disease", I500 = "Heart failure")

ind_list <- list()
model_list <- list()
pair_list <- list()
inc_list <- list()

for (code in codes) {
  ddir <- file.path(ROOT, ANALYSIS_SPLIT, code)

  ind_file <- find_one(ddir, "individual_predictor_summary\\.tsv$")
  model_file <- find_one(ddir, "model_panel_summary\\.tsv$")
  pair_file <- find_one(ddir, "paired_delta_cindex\\.tsv$", required = FALSE)
  inc_file <- find_one(ddir, "epoch_beyond_10_biomarker_panel\\.tsv$", required = FALSE)

  message(code, " individual predictor file: ", basename(ind_file))
  message(code, " model panel file: ", basename(model_file))

  ind_list[[code]] <- read_tsv_safe(ind_file)
  model_list[[code]] <- read_tsv_safe(model_file)

  if (!is.na(pair_file)) {
    message(code, " paired C-index file: ", basename(pair_file))
    pair_list[[code]] <- read_tsv_safe(pair_file)
  }
  if (!is.na(inc_file)) {
    inc_list[[code]] <- read_tsv_safe(inc_file)
  }
}

ind <- bind_rows(ind_list, .id = "source_code") %>%
  mutate(
    disease_code = ifelse(is.na(disease_code), source_code, disease_code),
    disease_label = ifelse(is.na(disease_label), disease_labels[disease_code], disease_label),
    predictor_type = ifelse(predictor_type == "EPOCH", "Brain EPOCH", "Conventional biomarker"),
    predictor_order = as.numeric(predictor_order),
    hr_per_1sd = as.numeric(hr_per_1sd),
    hr_ci_lower = as.numeric(hr_ci_lower),
    hr_ci_upper = as.numeric(hr_ci_upper),
    p_value = as.numeric(p_value),
    cindex_baseline_covariates = as.numeric(cindex_baseline_covariates),
    cindex_baseline_plus_predictor = as.numeric(cindex_baseline_plus_predictor),
    delta_cindex_vs_baseline = as.numeric(delta_cindex_vs_baseline),
    n_analysis_rows = as.numeric(n_analysis_rows),
    n_cases = as.numeric(n_cases)
  )

models <- bind_rows(model_list, .id = "source_code") %>%
  mutate(
    disease_code = ifelse(is.na(disease_code), source_code, disease_code),
    disease_label = ifelse(is.na(disease_label), disease_labels[disease_code], disease_label),
    cindex = as.numeric(cindex),
    delta_cindex_vs_baseline = as.numeric(delta_cindex_vs_baseline),
    penalizer = as.numeric(penalizer)
  )

pairs <- if (length(pair_list) > 0) bind_rows(pair_list, .id = "source_code") else tibble()
incremental <- if (length(inc_list) > 0) bind_rows(inc_list, .id = "source_code") else tibble()

# Save combined tables for convenience.
write_tsv(ind, file.path(OUTDIR, "combined_individual_predictor_results.tsv"), na = "NA")
write_tsv(models, file.path(OUTDIR, "combined_model_panel_results.tsv"), na = "NA")
if (nrow(pairs) > 0) write_tsv(pairs, file.path(OUTDIR, "combined_paired_delta_cindex.tsv"), na = "NA")
if (nrow(incremental) > 0) write_tsv(incremental, file.path(OUTDIR, "combined_epoch_beyond_panel.tsv"), na = "NA")

# ------------------------------------------------------------------------------
# Forest plot function
# ------------------------------------------------------------------------------

make_forest <- function(df, target_code, panel_title) {
  d <- df %>%
    filter(.data$disease_code == .env$target_code) %>%
    arrange(predictor_order)

  # Defensive QC: each disease should contribute one row per predictor.
  # The previous version used `filter(.data$disease_code == disease_code)`,
  # where dplyr data masking caused the RHS to resolve to the column itself.
  # That retained BOTH diseases and produced duplicate factor levels.
  if (nrow(d) == 0) {
    stop("No individual-predictor rows found for disease code: ", target_code)
  }
  if (anyDuplicated(d$predictor)) {
    dup <- unique(d$predictor[duplicated(d$predictor)])
    stop(
      "Duplicate predictor rows remain within disease ", target_code, ": ",
      paste(dup, collapse = "; "),
      ". Check the input summary files."
    )
  }

  predictor_levels <- rev(unique(as.character(d$predictor)))

  d <- d %>%
    mutate(
      predictor_plot = factor(as.character(predictor), levels = predictor_levels),
      type_plot = factor(predictor_type, levels = c("Brain EPOCH", "Conventional biomarker"))
    )

  n_txt <- unique(d$n_analysis_rows)
  case_txt <- unique(d$n_cases)
  subtitle_txt <- paste0("N = ", comma(n_txt[[1]]), "; incident cases = ", comma(case_txt[[1]]),
                         "; HR per 1-SD higher predictor")

  ggplot(d, aes(x = predictor_plot, y = hr_per_1sd, colour = type_plot)) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.45, colour = "grey45") +
    geom_errorbar(aes(ymin = hr_ci_lower, ymax = hr_ci_upper), width = 0.18, linewidth = 0.65) +
    geom_point(size = 2.6) +
    coord_flip() +
    scale_y_log10(
      breaks = c(0.5, 0.75, 1, 1.5, 2, 2.5),
      limits = c(0.45, 2.65),
      labels = label_number(accuracy = 0.01)
    ) +
    scale_colour_manual(values = c("Brain EPOCH" = COL_EPOCH,
                                   "Conventional biomarker" = COL_BIO)) +
    labs(
      title = panel_title,
      subtitle = subtitle_txt,
      x = NULL,
      y = "Hazard ratio (95% CI)"
    ) +
    base_theme +
    theme(legend.position = "none")
}

pA <- make_forest(ind, "G309", "A  Alzheimer's disease")
pB <- make_forest(ind, "I500", "B  Heart failure")

# ------------------------------------------------------------------------------
# Panel C: model-level C-index comparison
# ------------------------------------------------------------------------------

model_label_map <- c(
  "Baseline covariates" = "Clinical covariates",
  "Baseline covariates + Brain EPOCH" = "+ Brain EPOCH",
  "Baseline covariates + 10-biomarker panel" = "+ 10-biomarker panel",
  "Baseline covariates + 10-biomarker panel + Brain EPOCH" = "+ 10-biomarker panel + EPOCH"
)

model_order <- unname(model_label_map)

models_plot <- models %>%
  mutate(
    model_short = recode(model, !!!model_label_map),
    model_short = factor(model_short, levels = rev(model_order)),
    disease_label = factor(disease_label, levels = c("Alzheimer's disease", "Heart failure")),
    c_label = sprintf("%.3f", cindex)
  )

# Data-driven x range, padded but not anchored at zero because C-index is concentrated near 1.
xmin_c <- max(0.5, min(models_plot$cindex, na.rm = TRUE) - 0.025)
xmax_c <- min(1.0, max(models_plot$cindex, na.rm = TRUE) + 0.035)

pC <- ggplot(models_plot, aes(x = cindex, y = model_short)) +
  geom_point(aes(colour = model_short), size = 3) +
  geom_text(aes(label = c_label), hjust = -0.45, size = 3.0, colour = "black") +
  facet_wrap(~ disease_label, nrow = 1) +
  scale_colour_manual(values = c(
    "Clinical covariates" = COL_BASE,
    "+ Brain EPOCH" = COL_EPOCH,
    "+ 10-biomarker panel" = COL_PANEL,
    "+ 10-biomarker panel + EPOCH" = COL_BOTH
  )) +
  scale_x_continuous(limits = c(xmin_c, xmax_c), labels = label_number(accuracy = 0.01)) +
  labs(
    title = "C  Disease discrimination across model specifications",
    subtitle = "Harrell C-index in the full analysis sample",
    x = "C-index",
    y = NULL
  ) +
  base_theme +
  theme(legend.position = "none")

# ------------------------------------------------------------------------------
# Panel D: EPOCH added beyond the complete 10-biomarker panel
# Prefer paired-bootstrap results if available.
# ------------------------------------------------------------------------------

if (nrow(pairs) > 0 && "comparison_type" %in% names(pairs) &&
    any(pairs$comparison_type == "incremental_EPOCH_beyond_10_biomarker_panel", na.rm = TRUE)) {
  dD <- pairs %>%
    filter(comparison_type == "incremental_EPOCH_beyond_10_biomarker_panel") %>%
    mutate(
      disease_label = ifelse(is.na(disease_label), disease_labels[disease_code], disease_label),
      delta = as.numeric(delta_cindex_a_minus_b),
      lo = as.numeric(delta_cindex_ci_lower),
      hi = as.numeric(delta_cindex_ci_upper),
      p = as.numeric(empirical_p_two_sided),
      label = paste0("ΔC = ", sprintf("%.3f", delta), "; P = ", fmt_p(p))
    )
  d_subtitle <- "EPOCH added to clinical covariates + all 10 biomarkers; paired bootstrap 95% CI"
} else {
  warning("Paired delta-C files not found. Panel D will show point estimates without bootstrap CIs.")
  dD <- models %>%
    select(disease_code, disease_label, model, cindex) %>%
    pivot_wider(names_from = model, values_from = cindex) %>%
    transmute(
      disease_code,
      disease_label,
      delta = `Baseline covariates + 10-biomarker panel + Brain EPOCH` -
              `Baseline covariates + 10-biomarker panel`,
      lo = NA_real_, hi = NA_real_, p = NA_real_,
      label = paste0("ΔC = ", sprintf("%.3f", delta))
    )
  d_subtitle <- "EPOCH added to clinical covariates + all 10 biomarkers; point estimates only"
}

dD <- dD %>%
  mutate(disease_label = factor(disease_label, levels = rev(c("Alzheimer's disease", "Heart failure"))))

# Make enough room for labels to the right.
d_lim <- max(abs(c(dD$lo, dD$hi, dD$delta)), na.rm = TRUE)
if (!is.finite(d_lim)) d_lim <- max(abs(dD$delta), na.rm = TRUE)
d_lim <- max(d_lim * 1.55, 0.01)

pD <- ggplot(dD, aes(y = disease_label, x = delta)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, colour = "grey45") +
  geom_segment(
    data = dD %>% filter(is.finite(lo), is.finite(hi)),
    aes(x = lo, xend = hi, y = disease_label, yend = disease_label),
    linewidth = 0.75, colour = COL_EPOCH
  ) +
  geom_point(size = 3.2, colour = COL_EPOCH) +
  geom_text(aes(label = label), hjust = -0.08, vjust = -0.8, size = 3.0) +
  scale_x_continuous(limits = c(min(-0.005, -0.08 * d_lim), d_lim),
                     labels = label_number(accuracy = 0.001)) +
  labs(
    title = "D  Incremental discrimination of Brain EPOCH beyond 10 biomarkers",
    subtitle = d_subtitle,
    x = expression(Delta*"C-index (panel + EPOCH minus panel)"),
    y = NULL
  ) +
  base_theme +
  theme(legend.position = "none")

# ------------------------------------------------------------------------------
# Assemble and save main figure
# ------------------------------------------------------------------------------

main_fig <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title = "Brain proteomics mortality EPOCH provides disease-predictive information beyond conventional biomarkers",
    subtitle = paste0(
      "Incident Alzheimer's disease and heart failure; ", ANALYSIS_SPLIT,
      " mortality-EPOCH sample. All predictor comparisons within each disease use the same complete-case population."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10.5, colour = "grey30"),
      plot.margin = margin(8, 8, 8, 8)
    )
  )

main_pdf <- file.path(OUTDIR, "Brain_proteomics_EPOCH_vs_10_biomarkers_main_figure.pdf")
main_png <- file.path(OUTDIR, "Brain_proteomics_EPOCH_vs_10_biomarkers_main_figure.png")

ggsave(main_pdf, main_fig, width = 13.2, height = 10.2, units = "in", device = "pdf")
ggsave(main_png, main_fig, width = 13.2, height = 10.2, units = "in", dpi = 450)

if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite(file.path(OUTDIR, "Brain_proteomics_EPOCH_vs_10_biomarkers_main_figure.svg"),
                   width = 13.2, height = 10.2)
  print(main_fig)
  dev.off()
}

# Save individual panels.
for (nm in names(list(A = pA, B = pB, C = pC, D = pD))) {
  pp <- list(A = pA, B = pB, C = pC, D = pD)[[nm]]
  ggsave(file.path(OUTDIR, paste0("panel_", nm, ".pdf")), pp,
         width = 6.5, height = 5.0, units = "in", device = "pdf")
  ggsave(file.path(OUTDIR, paste0("panel_", nm, ".png")), pp,
         width = 6.5, height = 5.0, units = "in", dpi = 450)
}

# ------------------------------------------------------------------------------
# Supplementary: individual predictor delta C vs clinical covariates
# ------------------------------------------------------------------------------

supp <- ind %>%
  mutate(
    predictor_plot = factor(predictor, levels = rev(unique(ind %>% arrange(predictor_order) %>% pull(predictor)))),
    type_plot = predictor_type
  )

p_supp <- ggplot(supp, aes(x = delta_cindex_vs_baseline, y = predictor_plot, colour = type_plot)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.45, colour = "grey45") +
  geom_point(size = 2.8) +
  facet_wrap(~ disease_label, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = c("Brain EPOCH" = COL_EPOCH,
                                 "Conventional biomarker" = COL_BIO)) +
  labs(
    title = "Incremental C-index of Brain EPOCH and conventional biomarkers",
    subtitle = "Each point is the C-index change after adding one standardized predictor to the same clinical-covariate model",
    x = expression(Delta*"C-index vs clinical covariates"),
    y = NULL
  ) +
  base_theme +
  theme(legend.position = "bottom")

ggsave(file.path(OUTDIR, "Supplementary_individual_predictor_delta_cindex.pdf"),
       p_supp, width = 11.5, height = 6.0, units = "in", device = "pdf")
ggsave(file.path(OUTDIR, "Supplementary_individual_predictor_delta_cindex.png"),
       p_supp, width = 11.5, height = 6.0, units = "in", dpi = 450)

message("Finished. Main figure:")
message("  ", main_pdf)
message("  ", main_png)