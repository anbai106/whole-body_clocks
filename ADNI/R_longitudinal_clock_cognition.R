#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(scales)
})

# ============================================================
# Paths
# ============================================================

BASE_DIR <- "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/adni_lepoch"

COG_DIR <- file.path(
  BASE_DIR,
  "baseline_associations",
  "cognition_biomarker_comparison"
)

SLOPE_DATA_FILE <- file.path(COG_DIR, "analysis_dataset_ad_epoch_slope.tsv")

BASELINE_EPOCH_FILE <- file.path(
  BASE_DIR,
  "results_brain_mri_ad_lepoch",
  "adni_brain_mri_ad_lepoch_predictions.tsv"
)

OUTDIR <- file.path(
  BASE_DIR,
  "baseline_associations",
  "ad_epoch_slope_cognition_figure_slope_filtered"
)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# User settings
# ============================================================

cognitive_cols <- c(
  "Animal_Fluency",
  "RAVLT",
  "TMT_A",
  "TMT_B"
)

cognitive_labels <- c(
  "Animal_Fluency" = "Animal Fluency",
  "RAVLT" = "RAVLT",
  "TMT_A" = "TMT-A",
  "TMT_B" = "TMT-B"
)

covariates <- c(
  "Age",
  "Sex",
  "Education_Years",
  "APOE4_Alleles",
  "DLICV",
  "SITE"
)

slope_col <- "AD_EPOCH_slope_per_year"

# Main requested filter:
# Keep participants with -1 < AD EPOCH slope per year < 1.
slope_filter_min <- -1
slope_filter_max <- 1

window_days <- 365

# ============================================================
# Helper functions
# ============================================================

p_format <- function(p) {
  ifelse(
    is.na(p), "NA",
    ifelse(
      p < 1e-4,
      formatC(p, format = "e", digits = 2),
      sprintf("%.4f", p)
    )
  )
}

num_format <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", sprintf(paste0("%.", digits, "f"), x))
}

safe_z <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - m) / s
}

partial_r_from_lm <- function(fit, term) {
  sm <- summary(fit)
  
  if (!term %in% rownames(sm$coefficients)) {
    return(NA_real_)
  }
  
  tval <- sm$coefficients[term, "t value"]
  df <- fit$df.residual
  
  if (!is.finite(tval) || !is.finite(df) || df <= 0) {
    return(NA_real_)
  }
  
  as.numeric(tval / sqrt(tval^2 + df))
}

get_model_stats <- function(df, cog_col, slope_col, covariates) {
  needed <- c(cog_col, slope_col, covariates)
  needed <- needed[needed %in% names(df)]
  
  dat <- df %>%
    dplyr::select(dplyr::all_of(needed)) %>%
    mutate(
      cognitive_z = safe_z(.data[[cog_col]]),
      slope_z = safe_z(.data[[slope_col]])
    )
  
  model_vars <- c("slope_z", "cognitive_z", covariates[covariates %in% names(dat)])
  
  dat <- dat %>%
    dplyr::select(dplyr::all_of(model_vars)) %>%
    tidyr::drop_na()
  
  if (nrow(dat) < 30) {
    return(data.frame(
      cognitive_score = cog_col,
      n = nrow(dat),
      std_beta = NA_real_,
      se = NA_real_,
      t = NA_real_,
      p = NA_real_,
      partial_r = NA_real_,
      r2 = NA_real_,
      status = "insufficient_n"
    ))
  }
  
  form <- as.formula(
    paste(
      "slope_z ~ cognitive_z +",
      paste(covariates[covariates %in% names(dat)], collapse = " + ")
    )
  )
  
  fit <- lm(form, data = dat)
  sm <- summary(fit)
  
  if (!"cognitive_z" %in% rownames(sm$coefficients)) {
    return(data.frame(
      cognitive_score = cog_col,
      n = nrow(dat),
      std_beta = NA_real_,
      se = NA_real_,
      t = NA_real_,
      p = NA_real_,
      partial_r = NA_real_,
      r2 = NA_real_,
      status = "missing_term"
    ))
  }
  
  data.frame(
    cognitive_score = cog_col,
    n = nrow(dat),
    std_beta = sm$coefficients["cognitive_z", "Estimate"],
    se = sm$coefficients["cognitive_z", "Std. Error"],
    t = sm$coefficients["cognitive_z", "t value"],
    p = sm$coefficients["cognitive_z", "Pr(>|t|)"],
    partial_r = partial_r_from_lm(fit, "cognitive_z"),
    r2 = sm$r.squared,
    status = "ok"
  )
}

# ============================================================
# Read data
# ============================================================

if (!file.exists(SLOPE_DATA_FILE)) {
  stop("Cannot find slope analysis dataset: ", SLOPE_DATA_FILE)
}

df_raw <- fread(SLOPE_DATA_FILE) %>% as.data.frame()

missing_cols <- setdiff(c("PTID", slope_col, cognitive_cols, covariates), names(df_raw))
if (length(missing_cols) > 0) {
  stop("Missing required columns in slope dataset: ", paste(missing_cols, collapse = ", "))
}

# ============================================================
# Apply AD EPOCH slope-per-year filter
# ============================================================

df <- df_raw %>%
  mutate(
    AD_EPOCH_slope_per_year_raw = as.numeric(.data[[slope_col]])
  ) %>%
  filter(
    is.finite(AD_EPOCH_slope_per_year_raw),
    AD_EPOCH_slope_per_year_raw > slope_filter_min,
    AD_EPOCH_slope_per_year_raw < slope_filter_max
  ) %>%
  mutate(
    !!slope_col := AD_EPOCH_slope_per_year_raw
  )

filter_qc <- data.frame(
  n_before_slope_filter = nrow(df_raw),
  n_after_slope_filter = nrow(df),
  n_removed_by_slope_filter = nrow(df_raw) - nrow(df),
  slope_filter_min_exclusive = slope_filter_min,
  slope_filter_max_exclusive = slope_filter_max
)

fwrite(
  filter_qc,
  file.path(OUTDIR, "ad_epoch_slope_filter_qc.tsv"),
  sep = "\t"
)

message("Slope filter QC:")
print(filter_qc)

# ============================================================
# Add CN-to-MCI / CN-to-AD / censored group
# ============================================================

if (file.exists(BASELINE_EPOCH_FILE)) {
  base_pred <- fread(BASELINE_EPOCH_FILE) %>% as.data.frame()
  
  if (all(c("PTID", "event", "event_or_censor_dx") %in% names(base_pred))) {
    group_df <- base_pred %>%
      dplyr::select(PTID, event, event_or_censor_dx) %>%
      distinct() %>%
      mutate(
        event = as.character(event),
        event_bool = event %in% c("True", "TRUE", "true", "1", "T"),
        outcome_group = case_when(
          event_bool & event_or_censor_dx == "MCI" ~ "CN to MCI",
          event_bool & event_or_censor_dx == "AD"  ~ "CN to AD",
          TRUE ~ "Censored / non-event"
        )
      ) %>%
      dplyr::select(PTID, outcome_group)
    
    df <- df %>%
      left_join(group_df, by = "PTID")
  }
}

if (!"outcome_group" %in% names(df)) {
  df$outcome_group <- "Censored / non-event"
}

df$outcome_group <- factor(
  df$outcome_group,
  levels = c("Censored / non-event", "CN to MCI", "CN to AD")
)

# Save filtered analysis dataset
fwrite(
  df,
  file.path(OUTDIR, "analysis_dataset_ad_epoch_slope_filtered_minus1_to_1.tsv"),
  sep = "\t"
)

# ============================================================
# Re-run adjusted statistics after slope filtering
# ============================================================

stats_tbl <- bind_rows(
  lapply(cognitive_cols, function(cog) {
    get_model_stats(df, cog, slope_col, covariates)
  })
) %>%
  mutate(
    fdr = p.adjust(p, method = "BH"),
    cognitive_label = cognitive_labels[cognitive_score],
    annotation = paste0(
      "N = ", n,
      "\nstd \u03b2 = ", num_format(std_beta, 3),
      "\npartial r = ", num_format(partial_r, 3),
      "\nP = ", p_format(p),
      "\nFDR = ", p_format(fdr)
    ),
    slope_filter = paste0(
      slope_filter_min,
      " < ",
      slope_col,
      " < ",
      slope_filter_max
    )
  )

fwrite(
  stats_tbl,
  file.path(OUTDIR, "ad_epoch_slope_vs_cognition_adjusted_stats_filtered_minus1_to_1.tsv"),
  sep = "\t"
)

# ============================================================
# Long plotting table
# ============================================================

plot_df <- df %>%
  dplyr::select(
    PTID,
    outcome_group,
    dplyr::all_of(slope_col),
    dplyr::all_of(cognitive_cols)
  ) %>%
  pivot_longer(
    cols = dplyr::all_of(cognitive_cols),
    names_to = "cognitive_score",
    values_to = "cognitive_value"
  ) %>%
  mutate(
    cognitive_label = cognitive_labels[cognitive_score],
    cognitive_label = factor(
      cognitive_label,
      levels = cognitive_labels[cognitive_cols]
    ),
    AD_EPOCH_slope_per_year = .data[[slope_col]]
  ) %>%
  filter(
    is.finite(cognitive_value),
    is.finite(AD_EPOCH_slope_per_year),
    AD_EPOCH_slope_per_year > slope_filter_min,
    AD_EPOCH_slope_per_year < slope_filter_max
  )

# Annotation positions per facet
annot_df <- plot_df %>%
  group_by(cognitive_score, cognitive_label) %>%
  summarise(
    x_pos = quantile(cognitive_value, 0.03, na.rm = TRUE),
    y_pos = min(
      quantile(AD_EPOCH_slope_per_year, 0.97, na.rm = TRUE),
      slope_filter_max - 0.05
    ),
    .groups = "drop"
  ) %>%
  left_join(
    stats_tbl %>% dplyr::select(cognitive_score, annotation),
    by = "cognitive_score"
  )

# ============================================================
# Plot
# ============================================================

p <- ggplot(
  plot_df,
  aes(
    x = cognitive_value,
    y = AD_EPOCH_slope_per_year
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.35,
    color = "grey65"
  ) +
  geom_point(
    aes(color = outcome_group),
    alpha = 0.75,
    size = 1.8
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "black",
    linewidth = 0.8
  ) +
  geom_label(
    data = annot_df,
    aes(
      x = x_pos,
      y = y_pos,
      label = annotation
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.1,
    label.size = 0.25,
    fill = "white",
    alpha = 0.92
  ) +
  facet_wrap(
    ~ cognitive_label,
    scales = "free_x",
    nrow = 1,
    strip.position = "top"
  ) +
  coord_cartesian(
    ylim = c(slope_filter_min, slope_filter_max)
  ) +
  scale_color_manual(
    values = c(
      "Censored / non-event" = "#B08A44",
      "CN to MCI" = "#4C9AC7",
      "CN to AD" = "#C76E7C"
    ),
    drop = FALSE
  ) +
  labs(
    x = "Cognitive score",
    y = "AD EPOCH slope per year",
    color = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, color = "grey88"),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    strip.text = element_text(face = "bold", size = 10),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    legend.position = c(0.54, 0.88),
    legend.background = element_rect(fill = "white", color = "grey60", linewidth = 0.25),
    legend.text = element_text(size = 9),
    plot.margin = margin(6, 8, 6, 8)
  )

# ============================================================
# Save
# ============================================================

pdf_file <- file.path(
  OUTDIR,
  "ad_epoch_slope_per_year_vs_cognition_filtered_minus1_to_1.pdf"
)

png_file <- file.path(
  OUTDIR,
  "ad_epoch_slope_per_year_vs_cognition_filtered_minus1_to_1.png"
)

ggsave(pdf_file, p, width = 12.5, height = 4.2, device = cairo_pdf)
ggsave(png_file, p, width = 12.5, height = 4.2, dpi = 300)

message("Saved filtered figure:")
message(pdf_file)
message(png_file)

message("Saved filtered adjusted statistics:")
message(file.path(
  OUTDIR,
  "ad_epoch_slope_vs_cognition_adjusted_stats_filtered_minus1_to_1.tsv"
))

message("Saved slope-filter QC:")
message(file.path(OUTDIR, "ad_epoch_slope_filter_qc.tsv"))