#!/usr/bin/env Rscript

# =============================================================================
# STEP 6: MHAS-only panels I, J, and K for mortality EPOCH
# FINAL REVISED VERSION
# =============================================================================
#
# What this script fixes:
#   1) Panel I is explicitly MHAS-only.
#      It uses Train / Validation / Test labels, not NHANES calendar-year labels.
#
#   2) Panel I bars are drawn correctly even when the y-axis is zoomed.
#      Instead of geom_col() from 0, the bars are drawn from the visible y-axis
#      baseline to the C-index value using geom_rect().
#
#   3) Panel J keeps the cumulative mortality/KM plot by EPOCH quartile.
#
#   4) Panel K uses only updated Step 4 baseline disease-burden associations:
#        analysis_type == "baseline_prevalent"
#      It includes:
#        - disease-domain legend
#        - colored points and CI bars by disease domain
#        - OR (95% CI) text at the far right side of the plot
#
# Required inputs:
#   /Users/hao/Dropbox/MHAS/step2_mortality_epoch_model/
#     mhas_mortality_epoch_predictions.tsv
#
#   /Users/hao/Dropbox/MHAS/step3_mortality_epoch_validation/
#     mhas_step3_cindex_bootstrap_ci.tsv
#
#   /Users/hao/Dropbox/MHAS/step4_epoch_disease_associations/
#     mhas_step4_disease_association_summary.tsv
#
# Outputs:
#   /Users/hao/Dropbox/MHAS/step6_panels_ijk/
#     Step6_panel_i_model_performance.pdf/png
#     Step6_panel_j_km_mortality_quartile.pdf/png
#     Step6_panel_k_baseline_burden_forest.pdf/png
#     Step6_panels_IJK_combined.pdf/png, if patchwork is installed
#     Step6_panel_i_plot_data.tsv
#     Step6_panel_j_km_plot_data.tsv
#     Step6_panel_k_plot_data.tsv
#     Step6_panels_ijk_audit.txt
#
# Run:
#   cd /Users/hao/Project/whole-body_clocks/MHAS
#   Rscript 6_step6_plot_mhas_panels_IJK_MHAS_final.R
#
# Optional:
#   MHAS_BASE_DIR=/Users/hao/Dropbox/MHAS \
#   PANEL_I_SCORE=lp_total \
#   KM_SPLIT=test \
#   DISEASE_SPLIT=all \
#   TOP_DISEASE_N=12 \
#   Rscript 6_step6_plot_mhas_panels_IJK_MHAS_final.R
#
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(survival)
})

HAS_PATCHWORK <- requireNamespace("patchwork", quietly = TRUE)

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------

BASE_DIR <- Sys.getenv("MHAS_BASE_DIR", unset = "/Users/hao/Dropbox/MHAS")

STEP2_DIR <- file.path(BASE_DIR, "step2_mortality_epoch_model")
STEP3_DIR <- file.path(BASE_DIR, "step3_mortality_epoch_validation")
STEP4_DIR <- file.path(BASE_DIR, "step4_epoch_disease_associations")

OUT_DIR <- Sys.getenv("STEP6_OUT_DIR", unset = file.path(BASE_DIR, "step6_panels_ijk"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Panel I score:
#   lp_total = full mortality EPOCH Cox linear predictor; recommended for discrimination.
#   mortality_epoch_acceleration_z = age/sex-residualized EPOCH acceleration.
#   clinical_baseline_lp = clinical comparator, if available.
PANEL_I_SCORE <- Sys.getenv("PANEL_I_SCORE", unset = "lp_total")

# Panel J:
#   test = internal hold-out split for mortality stratification.
#   all = descriptive all-sample curve.
KM_SPLIT <- Sys.getenv("KM_SPLIT", unset = "test")

# Panel K:
#   all is recommended because baseline disease-burden analysis is cross-sectional.
DISEASE_SPLIT <- Sys.getenv("DISEASE_SPLIT", unset = "all")

TOP_DISEASE_N <- as.integer(Sys.getenv("TOP_DISEASE_N", unset = "12"))

# Optional KM x-axis truncation. Empty/NA uses all available follow-up.
KM_MAX_YEARS_TXT <- Sys.getenv("KM_MAX_YEARS", unset = "")
KM_MAX_YEARS <- suppressWarnings(as.numeric(KM_MAX_YEARS_TXT))
if (length(KM_MAX_YEARS) == 0 || is.na(KM_MAX_YEARS)) KM_MAX_YEARS <- NA_real_

# Whether to show N and deaths under Train/Validation/Test in panel I.
SHOW_PANEL_I_COUNTS <- tolower(Sys.getenv("SHOW_PANEL_I_COUNTS", unset = "TRUE")) %in% c("true", "1", "yes", "y")

# -----------------------------------------------------------------------------
# Colors matched to existing panels F, G, and H
# -----------------------------------------------------------------------------

# Panel F-like Train / Validation / Test colors.
SPLIT_COLORS <- c(
  "Train" = "#2E91A8",       # teal-blue
  "Validation" = "#E69F00",  # orange
  "Test" = "#76A95C"         # green
)

# Panel G-like EPOCH quartile colors.
QUARTILE_COLORS <- c(
  "Q1 lowest" = "#4DBD73",   # green
  "Q2" = "#6B9FD8",          # blue
  "Q3" = "#D6A84A",          # gold
  "Q4 highest" = "#F06A4A"   # orange-red
)

# Panel H-like disease-domain colors.
DOMAIN_COLORS <- c(
  "Cardiovascular" = "#D55E00",
  "Metabolic" = "#0072B2",
  "Respiratory" = "#009E73",
  "Neurological" = "#9467BD",
  "Cancer" = "#CC79A7",
  "Musculoskeletal" = "#8C564B",
  "Other" = "#666666"
)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

message2 <- function(...) {
  message(sprintf(...))
}

read_tsv_required <- function(path) {
  if (!file.exists(path)) {
    stop("Required file does not exist: ", path)
  }
  data.table::fread(path)
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

theme_epoch <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = element_text(size = base_size, color = "grey30", hjust = 0),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      strip.background = element_rect(fill = "grey95", color = "grey70"),
      strip.text = element_text(face = "bold"),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.25),
      panel.grid.major.x = element_blank()
    )
}

save_plot <- function(p, basename, width = 6, height = 4.5) {
  pdf_path <- file.path(OUT_DIR, paste0(basename, ".pdf"))
  png_path <- file.path(OUT_DIR, paste0(basename, ".png"))

  ggsave(pdf_path, p, width = width, height = height, device = cairo_pdf)
  ggsave(png_path, p, width = width, height = height, dpi = 300)

  message2("Saved: %s", pdf_path)
  message2("Saved: %s", png_path)
}

pretty_endpoint <- function(x) {
  x %>%
    str_replace_all("_incl_", "_including_") %>%
    str_replace_all("_", " ") %>%
    str_replace_all("respiratory disease including asthma", "respiratory disease/asthma") %>%
    str_to_sentence()
}

format_p <- function(p) {
  p <- safe_num(p)
  dplyr::case_when(
    is.na(p) ~ "NA",
    p < 1e-300 ~ "<1e-300",
    p < 0.001 ~ sprintf("%.2e", p),
    TRUE ~ sprintf("%.3f", p)
  )
}

format_or_ci <- function(or, lo, hi) {
  sprintf("%.2f (%.2f, %.2f)", safe_num(or), safe_num(lo), safe_num(hi))
}

prefer_adjusted <- function(dt) {
  if (!"adjusted" %in% names(dt)) return(dt)
  adjusted_bool <- as.logical(dt$adjusted)
  adjusted_bool[is.na(adjusted_bool)] <- FALSE
  if (any(adjusted_bool)) {
    return(dt[adjusted_bool, ])
  }
  dt[!adjusted_bool, ]
}

require_baseline_prevalent <- function(dt, table_name) {
  if (!"analysis_type" %in% names(dt)) {
    stop(
      table_name, " does not contain analysis_type. ",
      "This looks like old Step 4 output. Re-run the updated baseline-prevalent Step 4 script."
    )
  }

  out <- dt[as.character(dt$analysis_type) == "baseline_prevalent", ]

  if (nrow(out) == 0) {
    stop(
      table_name, " contains no analysis_type == baseline_prevalent rows. ",
      "This looks like old incident/future-disease output."
    )
  }

  out
}

assign_disease_domain <- function(endpoint) {
  e <- tolower(as.character(endpoint))

  dplyr::case_when(
    str_detect(e, "heart|hypertension|stroke|mi|attack") ~ "Cardiovascular",
    str_detect(e, "diabetes") ~ "Metabolic",
    str_detect(e, "respiratory|asthma|copd") ~ "Respiratory",
    str_detect(e, "dementia|brain|neuro") ~ "Neurological",
    str_detect(e, "cancer") ~ "Cancer",
    str_detect(e, "arthritis|bone|musculo") ~ "Musculoskeletal",
    TRUE ~ "Other"
  )
}

tidy_survfit_cuminc <- function(fit) {
  s <- summary(fit)

  if (is.null(s$strata)) {
    strata <- rep("All", length(s$time))
  } else {
    strata <- as.character(s$strata)
    strata <- str_replace(strata, "^mortality_epoch_quartile=", "")
  }

  data.frame(
    time = s$time,
    n_risk = s$n.risk,
    n_event = s$n.event,
    survival = s$surv,
    lower = s$lower,
    upper = s$upper,
    quartile = strata,
    cumulative_mortality = 1 - s$surv,
    cumulative_mortality_lower = 1 - s$upper,
    cumulative_mortality_upper = 1 - s$lower
  )
}

split_simple_label <- function(split) {
  dplyr::recode(
    as.character(split),
    "train" = "Train",
    "validation" = "Validation",
    "test" = "Test",
    .default = as.character(split)
  )
}

split_display_label <- function(split, n, deaths) {
  simple <- split_simple_label(split)

  if (!SHOW_PANEL_I_COUNTS || is.na(n) || is.na(deaths)) {
    return(simple)
  }

  paste0(
    simple,
    "\nN=", format(as.integer(n), big.mark = ","),
    "\nDeaths=", format(as.integer(deaths), big.mark = ",")
  )
}

# -----------------------------------------------------------------------------
# Input files
# -----------------------------------------------------------------------------

pred_file <- file.path(STEP2_DIR, "mhas_mortality_epoch_predictions.tsv")
cindex_file <- file.path(STEP3_DIR, "mhas_step3_cindex_bootstrap_ci.tsv")
disease_file <- file.path(STEP4_DIR, "mhas_step4_disease_association_summary.tsv")

pred_dt <- read_tsv_required(pred_file)
cindex_dt <- read_tsv_required(cindex_file)
disease_dt <- read_tsv_required(disease_file)
disease_dt <- require_baseline_prevalent(disease_dt, "mhas_step4_disease_association_summary.tsv")

# -----------------------------------------------------------------------------
# Panel I: MHAS model performance
# -----------------------------------------------------------------------------

make_panel_i <- function(cindex_dt, score_use = PANEL_I_SCORE) {
  required <- c("split", "score", "cindex", "ci_lower", "ci_upper")
  missing <- setdiff(required, names(cindex_dt))

  if (length(missing) > 0) {
    stop("C-index table missing required columns: ", paste(missing, collapse = ", "))
  }

  score_label <- dplyr::case_when(
    score_use == "lp_total" ~ "Mortality EPOCH LP",
    score_use == "mortality_epoch_acceleration_z" ~ "Mortality EPOCH acceleration",
    score_use == "clinical_baseline_lp" ~ "Clinical baseline",
    TRUE ~ score_use
  )

  dt <- as.data.frame(cindex_dt) %>%
    filter(
      split %in% c("train", "validation", "test"),
      score == score_use
    ) %>%
    mutate(
      split_order = as.integer(factor(split, levels = c("train", "validation", "test"))),
      split_simple = split_simple_label(split),
      cindex = safe_num(cindex),
      ci_lower = safe_num(ci_lower),
      ci_upper = safe_num(ci_upper),
      n = if ("n" %in% names(.)) safe_num(n) else NA_real_,
      deaths = if ("deaths" %in% names(.)) safe_num(deaths) else NA_real_
    ) %>%
    filter(!is.na(cindex)) %>%
    arrange(split_order)

  if (nrow(dt) == 0) {
    stop("No C-index rows found for PANEL_I_SCORE=", score_use)
  }

  # Make explicit x positions and visible bar baseline.
  y_min <- max(0.50, floor((min(dt$ci_lower, na.rm = TRUE) - 0.03) * 20) / 20)
  y_max <- min(1.00, ceiling((max(dt$ci_upper, na.rm = TRUE) + 0.03) * 20) / 20)

  # Ensure there is visible height even if cindex is close to y_min.
  if (any(dt$cindex <= y_min)) {
    y_min <- max(0.40, min(dt$cindex, na.rm = TRUE) - 0.05)
  }

  dt <- dt %>%
    mutate(
      x = split_order,
      xmin = x - 0.34,
      xmax = x + 0.34,
      y_base = y_min,
      split_label = mapply(split_display_label, split, n, deaths),
      split_simple = factor(split_simple, levels = c("Train", "Validation", "Test"))
    )

  data.table::fwrite(dt, file.path(OUT_DIR, "Step6_panel_i_plot_data.tsv"), sep = "\t")

  p <- ggplot(dt) +
    geom_rect(
      aes(xmin = xmin, xmax = xmax, ymin = y_base, ymax = cindex, fill = split_simple),
      color = "white",
      linewidth = 0.25
    ) +
    geom_errorbar(
      aes(x = x, ymin = ci_lower, ymax = ci_upper),
      width = 0.16,
      linewidth = 0.55,
      color = "black"
    ) +
    geom_text(
      aes(x = x, y = cindex, label = sprintf("%.3f", cindex)),
      vjust = -0.55,
      size = 3.4,
      color = "black"
    ) +
    scale_fill_manual(values = SPLIT_COLORS, guide = "none") +
    scale_x_continuous(
      breaks = dt$x,
      labels = dt$split_label,
      expand = expansion(mult = c(0.08, 0.08))
    ) +
    scale_y_continuous(
      breaks = pretty_breaks(n = 4),
      expand = expansion(mult = c(0, 0.05))
    ) +
    coord_cartesian(ylim = c(y_min, y_max), clip = "off") +
    labs(
      title = "I",
      subtitle = paste0("Model performance: ", score_label),
      x = NULL,
      y = "Harrell C-index"
    ) +
    theme_epoch(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 18),
      axis.text.x = element_text(size = 8.4, lineheight = 0.95),
      plot.margin = margin(t = 8, r = 10, b = 8, l = 8)
    )

  p
}

panel_i <- make_panel_i(cindex_dt, PANEL_I_SCORE)
save_plot(panel_i, "Step6_panel_i_model_performance", width = 4.8, height = 4.6)

# -----------------------------------------------------------------------------
# Panel J: KM/cumulative mortality by EPOCH quartile
# -----------------------------------------------------------------------------

make_panel_j <- function(pred_dt, split_use = KM_SPLIT, max_years = KM_MAX_YEARS) {
  required <- c("split", "followup_years", "event_death", "mortality_epoch_quartile")
  missing <- setdiff(required, names(pred_dt))

  if (length(missing) > 0) {
    stop("Prediction table missing required columns: ", paste(missing, collapse = ", "))
  }

  dt <- as.data.frame(pred_dt) %>%
    mutate(
      followup_years = safe_num(followup_years),
      event_death = as.integer(safe_num(event_death)),
      mortality_epoch_quartile = as.character(mortality_epoch_quartile)
    ) %>%
    filter(!is.na(followup_years), !is.na(event_death), !is.na(mortality_epoch_quartile))

  if (split_use != "all") {
    dt <- dt %>% filter(split == split_use)
  }

  dt <- dt %>%
    filter(mortality_epoch_quartile %in% c("Q1_lowest", "Q2", "Q3", "Q4_highest")) %>%
    mutate(
      mortality_epoch_quartile = factor(
        mortality_epoch_quartile,
        levels = c("Q1_lowest", "Q2", "Q3", "Q4_highest"),
        labels = c("Q1 lowest", "Q2", "Q3", "Q4 highest")
      )
    )

  if (!is.na(max_years)) {
    dt <- dt %>%
      mutate(
        event_death = ifelse(followup_years > max_years, 0L, event_death),
        followup_years = pmin(followup_years, max_years)
      )
  }

  if (nrow(dt) == 0 || length(unique(dt$event_death)) < 2) {
    stop("KM input has no valid rows or no mortality events for split=", split_use)
  }

  fit <- survival::survfit(
    survival::Surv(followup_years, event_death) ~ mortality_epoch_quartile,
    data = dt
  )

  km <- tidy_survfit_cuminc(fit)
  km$quartile <- factor(km$quartile, levels = c("Q1 lowest", "Q2", "Q3", "Q4 highest"))

  data.table::fwrite(km, file.path(OUT_DIR, "Step6_panel_j_km_plot_data.tsv"), sep = "\t")

  split_label <- ifelse(split_use == "all", "all participants", paste0(split_use, " split"))

  p <- ggplot(km, aes(x = time, y = cumulative_mortality, color = quartile)) +
    geom_step(linewidth = 1.0) +
    scale_color_manual(values = QUARTILE_COLORS, drop = FALSE) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      title = "J",
      subtitle = paste0("Cumulative mortality by EPOCH quartile: ", split_label),
      x = "Years after baseline exam",
      y = "Cumulative mortality",
      color = "EPOCH quartile"
    ) +
    theme_epoch(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 18),
      legend.position = "top",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 9)
    )

  p
}

panel_j <- make_panel_j(pred_dt, KM_SPLIT, KM_MAX_YEARS)
save_plot(panel_j, "Step6_panel_j_km_mortality_quartile", width = 6.9, height = 4.2)

# -----------------------------------------------------------------------------
# Panel K: Baseline disease-burden forest plot
# -----------------------------------------------------------------------------

make_panel_k <- function(disease_dt, split_use = DISEASE_SPLIT, top_n = TOP_DISEASE_N) {
  required <- c("analysis_type", "endpoint", "split", "score", "or_per_sd",
                "ci_lower", "ci_upper", "p", "status")
  missing <- setdiff(required, names(disease_dt))

  if (length(missing) > 0) {
    stop("Disease table missing required columns: ", paste(missing, collapse = ", "))
  }

  dt <- as.data.frame(disease_dt) %>%
    filter(
      analysis_type == "baseline_prevalent",
      split == split_use,
      score == "mortality_epoch_acceleration_z",
      str_detect(as.character(status), "^ok")
    )

  if ("adjusted" %in% names(dt)) {
    dt <- prefer_adjusted(dt)
  }

  dt <- dt %>%
    mutate(
      or_per_sd = safe_num(or_per_sd),
      ci_lower = safe_num(ci_lower),
      ci_upper = safe_num(ci_upper),
      p = safe_num(p),
      fdr_bh = if ("fdr_bh" %in% names(.)) safe_num(fdr_bh) else NA_real_,
      endpoint_label = pretty_endpoint(endpoint),
      disease_domain = assign_disease_domain(endpoint),
      or_ci_label = format_or_ci(or_per_sd, ci_lower, ci_upper)
    ) %>%
    filter(!is.na(or_per_sd), !is.na(ci_lower), !is.na(ci_upper), !is.na(p)) %>%
    arrange(p, desc(or_per_sd)) %>%
    slice_head(n = top_n) %>%
    arrange(or_per_sd) %>%
    mutate(
      endpoint_label = factor(endpoint_label, levels = endpoint_label),
      disease_domain = factor(
        disease_domain,
        levels = c("Cardiovascular", "Metabolic", "Respiratory", "Neurological",
                   "Cancer", "Musculoskeletal", "Other")
      )
    )

  if (nrow(dt) == 0) {
    stop("No valid baseline-prevalent disease rows found for split=", split_use)
  }

  data.table::fwrite(dt, file.path(OUT_DIR, "Step6_panel_k_plot_data.tsv"), sep = "\t")

  x_data_min <- min(dt$ci_lower, na.rm = TRUE)
  x_data_max <- max(dt$ci_upper, na.rm = TRUE)

  x_min <- max(0.05, x_data_min * 0.90)
  x_range <- max(0.1, x_data_max - x_min)
  x_annot <- x_data_max + x_range * 0.18
  x_max <- x_annot + x_range * 0.70

  p <- ggplot(dt, aes(x = or_per_sd, y = endpoint_label, color = disease_domain)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.45, color = "grey45") +
    geom_errorbarh(
      aes(xmin = ci_lower, xmax = ci_upper),
      height = 0.18,
      linewidth = 0.65
    ) +
    geom_point(size = 2.9) +
    geom_text(
      aes(x = x_annot, label = or_ci_label),
      hjust = 0,
      size = 2.9,
      color = "black"
    ) +
    annotate(
      "text",
      x = x_annot,
      y = nrow(dt) + 0.65,
      label = "OR (95% CI)",
      hjust = 0,
      fontface = "bold",
      size = 3.0,
      color = "black"
    ) +
    scale_color_manual(values = DOMAIN_COLORS, drop = TRUE, name = "Disease domain") +
    scale_x_continuous(
      limits = c(x_min, x_max),
      breaks = pretty_breaks(n = 4)
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = "K",
      subtitle = paste0("Baseline disease-burden associations: ", split_use, " split"),
      x = "Odds ratio per 1-SD EPOCH acceleration",
      y = NULL
    ) +
    theme_epoch(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 18),
      legend.position = "right",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8.5),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.25),
      plot.margin = margin(t = 8, r = 80, b = 8, l = 8)
    )

  p
}

panel_k <- make_panel_k(disease_dt, DISEASE_SPLIT, TOP_DISEASE_N)
save_plot(
  panel_k,
  "Step6_panel_k_baseline_burden_forest",
  width = 8.8,
  height = max(4.4, TOP_DISEASE_N * 0.35 + 1.5)
)

# -----------------------------------------------------------------------------
# Combined I-J-K figure
# -----------------------------------------------------------------------------

if (HAS_PATCHWORK) {
  combined <- (panel_i | panel_j) / panel_k +
    patchwork::plot_layout(heights = c(1.0, 1.18)) +
    patchwork::plot_annotation(
      title = "MHAS mortality EPOCH validation and baseline disease burden",
      subtitle = "Panels I-K: MHAS discrimination, mortality stratification, and baseline disease-burden associations"
    )

  save_plot(combined, "Step6_panels_IJK_combined", width = 12.8, height = 9.2)
} else {
  message("Package 'patchwork' not installed; skipping combined I-J-K figure.")
  message("Install with: install.packages('patchwork')")
}

# -----------------------------------------------------------------------------
# Audit
# -----------------------------------------------------------------------------

audit_lines <- c(
  "MHAS STEP 6 panels I-J-K audit",
  "================================",
  paste0("BASE_DIR: ", BASE_DIR),
  paste0("OUT_DIR: ", OUT_DIR),
  paste0("PANEL_I_SCORE: ", PANEL_I_SCORE),
  paste0("KM_SPLIT: ", KM_SPLIT),
  paste0("DISEASE_SPLIT: ", DISEASE_SPLIT),
  paste0("TOP_DISEASE_N: ", TOP_DISEASE_N),
  paste0("KM_MAX_YEARS: ", ifelse(is.na(KM_MAX_YEARS), "all available follow-up", KM_MAX_YEARS)),
  paste0("SHOW_PANEL_I_COUNTS: ", SHOW_PANEL_I_COUNTS),
  "",
  "Critical MHAS correction:",
  "Panel I uses MHAS train/validation/test labels only.",
  "No NHANES temporal labels such as 1999-2010/2011-2014/2015-2018 are used.",
  "Panel I bars are drawn with geom_rect from the visible y-axis baseline to the C-index value.",
  "",
  "Panel design:",
  "Panel I: bar plot with Train/Validation/Test colors matching panel F.",
  "Panel J: KM/cumulative mortality by EPOCH quartile matching panel G.",
  "Panel K: baseline disease-burden forest plot with disease-domain legend and far-right OR (95% CI) annotations.",
  "",
  "Input files:",
  paste0("Predictions: ", pred_file),
  paste0("C-index: ", cindex_file),
  paste0("Disease baseline burden: ", disease_file),
  "",
  "Output figures:",
  "Step6_panel_i_model_performance.pdf/png",
  "Step6_panel_j_km_mortality_quartile.pdf/png",
  "Step6_panel_k_baseline_burden_forest.pdf/png",
  "Step6_panels_IJK_combined.pdf/png, if patchwork is installed",
  "",
  "Analysis included:",
  "Baseline prevalent disease-burden associations only.",
  "",
  "Analysis excluded:",
  "Old incident/future disease analyses, old downstream-disease filenames, and NHANES temporal split labels."
)

writeLines(audit_lines, file.path(OUT_DIR, "Step6_panels_ijk_audit.txt"))
message2("Saved audit: %s", file.path(OUT_DIR, "Step6_panels_ijk_audit.txt"))

message("STEP 6 panels I-J-K finished successfully.")
