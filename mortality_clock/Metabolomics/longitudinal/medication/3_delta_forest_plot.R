# ============================================================
# Plot medication-cluster association with metabolomics
# delta mortality clocks.
#
# Reads outputs from:
#   Project/whole-body_clocks/mortality_clock/Metabolomics/longitudinal/medication/2_run_lm_across_clusters_and_delta_clocks.py
#
# Main plots:
#   1) Forest plot of cluster beta vs No/minimal medication
#   2) Adjusted mean delta clock age by medication cluster
#   3) Heatmap of beta effects and Bonferroni/nominal significance
#
# Significance encoding in forest plot:
#   Solid circle        = Bonferroni-significant, p_bonferroni < 0.05
#   Light filled circle = nominal only, raw P < 0.05 but Bonferroni P >= 0.05
#   Empty circle        = non-significant, raw P >= 0.05
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(stringr)
  library(forcats)
})

# -----------------------------
# 1. Paths
# -----------------------------

results_dir <- "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/mortality_clock/SA/delta_metabolomics_algorithmic_disease_onset/medication_cluster_delta_clock_lm_results"

out_dir <- file.path(results_dir, "plots_bonferroni")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

effects_tsv <- file.path(
  results_dir,
  "medication_cluster_delta_clock_lm_cluster_effects.tsv"
)

means_tsv <- file.path(
  results_dir,
  "medication_cluster_delta_clock_lm_adjusted_means.tsv"
)

counts_tsv <- file.path(
  results_dir,
  "medication_cluster_delta_clock_lm_cluster_counts.tsv"
)

summaries_tsv <- file.path(
  results_dir,
  "medication_cluster_delta_clock_lm_model_summaries.tsv"
)

# -----------------------------
# 2. Settings
# -----------------------------

organ_order <- c("Endocrine", "Digestive", "Hepatic", "Immune")

cluster_order <- c(
  "No/minimal medication",
  "Cardiometabolic medication cluster",
  "Respiratory medication cluster",
  "Psychiatric/pain medication cluster",
  "High polypharmacy cluster"
)

effect_cluster_order <- c(
  "Cardiometabolic medication cluster",
  "Respiratory medication cluster",
  "Psychiatric/pain medication cluster",
  "High polypharmacy cluster"
)

cluster_short <- c(
  "No/minimal medication" = "No/minimal",
  "Cardiometabolic medication cluster" = "Cardiometabolic",
  "Respiratory medication cluster" = "Respiratory",
  "Psychiatric/pain medication cluster" = "Psych/pain",
  "High polypharmacy cluster" = "High polypharmacy"
)

# Van Gogh-inspired, colorblind-friendly palette
cluster_palette <- c(
  "No/minimal medication" = "#6B6B6B",
  "Cardiometabolic medication cluster" = "#C9A227",
  "Respiratory medication cluster" = "#4E6FAE",
  "Psychiatric/pain medication cluster" = "#3F7F5F",
  "High polypharmacy cluster" = "#B65E16"
)

# Lighter versions for nominal-only signals
cluster_palette_light <- c(
  "No/minimal medication" = "#BDBDBD",
  "Cardiometabolic medication cluster" = "#E7D48A",
  "Respiratory medication cluster" = "#A9BCE3",
  "Psychiatric/pain medication cluster" = "#A7C7AE",
  "High polypharmacy cluster" = "#E2AA72"
)

plot_bg <- "#FFF9E8"
strip_bg <- "#F3E7C8"

sig_levels <- c(
  "Bonferroni-significant",
  "Nominal only",
  "Not significant"
)

# -----------------------------
# 3. Helper functions
# -----------------------------

fmt_p <- function(p) {
  ifelse(
    is.na(p), "NA",
    ifelse(
      p < 1e-300, "<1e-300",
      ifelse(
        p < 0.001,
        formatC(p, format = "e", digits = 2),
        formatC(p, format = "f", digits = 3)
      )
    )
  )
}

fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

# -----------------------------
# 4. Read inputs
# -----------------------------

if (!file.exists(effects_tsv)) {
  stop("Missing effects TSV: ", effects_tsv)
}
if (!file.exists(means_tsv)) {
  stop("Missing adjusted means TSV: ", means_tsv)
}

effects <- fread(effects_tsv)
means <- fread(means_tsv)

if (file.exists(counts_tsv)) {
  counts <- fread(counts_tsv)
} else {
  counts <- NULL
}

if (file.exists(summaries_tsv)) {
  summaries <- fread(summaries_tsv)
} else {
  summaries <- NULL
}

# -----------------------------
# 5. Clean fields and define Bonferroni/nominal significance
# -----------------------------

if (!"p_bonferroni" %in% names(effects)) {
  effects$p_bonferroni <- NA_real_
}

if (!"p_fdr_bh" %in% names(effects)) {
  effects$p_fdr_bh <- NA_real_
}

effect_num_cols <- c(
  "N", "N_reference", "N_exposure", "beta", "se", "ci_lo", "ci_hi",
  "p", "p_fdr_bh", "p_bonferroni"
)

for (cc in intersect(effect_num_cols, names(effects))) {
  effects[[cc]] <- safe_num(effects[[cc]])
}

# If p_bonferroni is missing, compute it across all status == ok tests.
valid_p_idx <- which(effects$status == "ok" & !is.na(effects$p))
if (length(valid_p_idx) > 0 && all(is.na(effects$p_bonferroni[valid_p_idx]))) {
  effects$p_bonferroni[valid_p_idx] <- p.adjust(effects$p[valid_p_idx], method = "bonferroni")
}

mean_num_cols <- c(
  "N_cluster", "adjusted_mean", "adjusted_mean_se",
  "adjusted_ci_lo", "adjusted_ci_hi"
)

for (cc in intersect(mean_num_cols, names(means))) {
  means[[cc]] <- safe_num(means[[cc]])
}

effects_plot <- effects %>%
  filter(status == "ok") %>%
  mutate(
    organ_label = factor(organ_label, levels = organ_order),
    exposure_cluster = factor(exposure_cluster, levels = effect_cluster_order),
    exposure_cluster_short = recode(as.character(exposure_cluster), !!!cluster_short),
    exposure_cluster_short = factor(
      exposure_cluster_short,
      levels = cluster_short[effect_cluster_order]
    ),
    
    bonf_sig = !is.na(p_bonferroni) & p_bonferroni < 0.05,
    nominal_only = !bonf_sig & !is.na(p) & p < 0.05,
    nonsig = is.na(p) | p >= 0.05,
    
    sig_group = case_when(
      bonf_sig ~ "Bonferroni-significant",
      nominal_only ~ "Nominal only",
      TRUE ~ "Not significant"
    ),
    sig_group = factor(sig_group, levels = sig_levels),
    
    annotation = paste0(
      "N=", comma(N_exposure),
      "; beta=", fmt_num(beta, 2),
      "; P=", fmt_p(p),
      "; Bonf P=", fmt_p(p_bonferroni)
    )
  )

means_plot <- means %>%
  mutate(
    organ_label = factor(organ_label, levels = organ_order),
    medication_cluster = factor(medication_cluster, levels = cluster_order),
    medication_cluster_short = recode(as.character(medication_cluster), !!!cluster_short),
    medication_cluster_short = factor(
      medication_cluster_short,
      levels = cluster_short[cluster_order]
    )
  )

# -----------------------------
# 6. Save cleaned tables
# -----------------------------

fwrite(
  effects_plot,
  file.path(out_dir, "medication_cluster_delta_clock_effects_for_plotting_bonferroni.tsv"),
  sep = "\t"
)

fwrite(
  means_plot,
  file.path(out_dir, "medication_cluster_delta_clock_adjusted_means_for_plotting.tsv"),
  sep = "\t"
)

fwrite(
  effects_plot %>% filter(bonf_sig),
  file.path(out_dir, "bonferroni_significant_medication_cluster_delta_clock_effects.tsv"),
  sep = "\t"
)

fwrite(
  effects_plot %>% filter(nominal_only),
  file.path(out_dir, "nominal_only_medication_cluster_delta_clock_effects.tsv"),
  sep = "\t"
)

# -----------------------------
# 7. Plot 1: Forest plot of beta effects
# -----------------------------

if (nrow(effects_plot) > 0) {
  
  x_min <- min(effects_plot$ci_lo, na.rm = TRUE)
  x_max <- max(effects_plot$ci_hi, na.rm = TRUE)
  x_pad <- 0.18 * (x_max - x_min)
  
  effects_plot <- effects_plot %>%
    mutate(
      x_text = x_max + 0.05 * (x_max - x_min)
    )
  
  p_forest <- ggplot(
    effects_plot,
    aes(y = exposure_cluster_short)
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "#4D4D4D",
      linewidth = 0.7
    ) +
    geom_segment(
      aes(
        x = ci_lo,
        xend = ci_hi,
        yend = exposure_cluster_short,
        color = exposure_cluster
      ),
      linewidth = 1.0,
      alpha = 0.95
    ) +
    
    # Bonferroni-significant: solid circle
    geom_point(
      data = effects_plot %>% filter(sig_group == "Bonferroni-significant"),
      aes(
        x = beta,
        color = exposure_cluster
      ),
      shape = 16,
      size = 3.9,
      stroke = 1.0,
      show.legend = FALSE
    ) +
    
    # Nominal-only: lighter filled circle
    geom_point(
      data = effects_plot %>% filter(sig_group == "Nominal only"),
      aes(
        x = beta,
        color = exposure_cluster,
        fill = exposure_cluster
      ),
      shape = 21,
      size = 3.7,
      stroke = 1.1,
      alpha = 0.55,
      show.legend = FALSE
    ) +
    
    # Non-significant: empty circle
    geom_point(
      data = effects_plot %>% filter(sig_group == "Not significant"),
      aes(
        x = beta,
        color = exposure_cluster
      ),
      shape = 1,
      size = 3.7,
      stroke = 1.1,
      show.legend = FALSE
    ) +
    
    # Clean significance legend using shape only.
    # drop = FALSE keeps all 3 categories even if one is absent in the data.
    geom_point(
      aes(
        x = beta,
        shape = sig_group
      ),
      color = "#333333",
      fill = "#BDBDBD",
      size = 3.6,
      stroke = 1.0,
      alpha = 0,
      show.legend = TRUE
    ) +
    
    geom_text(
      aes(
        x = x_text,
        label = annotation,
        color = exposure_cluster
      ),
      hjust = 0,
      size = 3.0,
      show.legend = FALSE
    ) +
    facet_wrap(
      ~ organ_label,
      ncol = 2,
      scales = "free_y"
    ) +
    scale_color_manual(
      values = cluster_palette,
      name = "Medication cluster"
    ) +
    scale_fill_manual(
      values = cluster_palette_light,
      name = "Medication cluster"
    ) +
    scale_shape_manual(
      values = c(
        "Bonferroni-significant" = 16,
        "Nominal only" = 21,
        "Not significant" = 1
      ),
      breaks = sig_levels,
      drop = FALSE,
      name = "Significance"
    ) +
    guides(
      shape = guide_legend(
        override.aes = list(
          alpha = 1,
          color = "#333333",
          fill = "#BDBDBD",
          size = 3.8
        )
      )
    ) +
    coord_cartesian(
      xlim = c(x_min - x_pad, x_max + 2.4 * x_pad),
      clip = "off"
    ) +
    labs(
      x = paste0(
        "Adjusted difference in delta clock age years\n",
        "versus no/minimal medication"
      ),
      y = NULL,
      title = "Medication clusters and longitudinal metabolomics mortality-clock change",
      subtitle = paste0(
        "Solid = Bonferroni-significant; light-filled = nominal only; empty = not significant. ",
        "Models adjusted for baseline clock acceleration, chronological age, sex, smoking, BMI, BP, and age change."
      )
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.background = element_rect(fill = plot_bg, color = NA),
      panel.background = element_rect(fill = plot_bg, color = NA),
      legend.background = element_rect(fill = plot_bg, color = NA),
      legend.box.background = element_rect(fill = plot_bg, color = NA),
      strip.background = element_rect(fill = strip_bg, color = NA),
      strip.text = element_text(face = "bold", size = 12, hjust = 0),
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10.5),
      axis.text.y = element_text(face = "bold", size = 10),
      axis.title.x = element_text(size = 12),
      legend.position = "bottom",
      legend.box = "vertical",
      plot.margin = margin(10, 205, 10, 10)
    )
  
  ggsave(
    file.path(out_dir, "medication_cluster_delta_clock_beta_forest_plot_bonferroni.pdf"),
    p_forest,
    width = 14.8,
    height = 7.8
  )
  
  ggsave(
    file.path(out_dir, "medication_cluster_delta_clock_beta_forest_plot_bonferroni.png"),
    p_forest,
    width = 14.8,
    height = 7.8,
    dpi = 300
  )
  
  ggsave(
    file.path(out_dir, "medication_cluster_delta_clock_beta_forest_plot_bonferroni.svg"),
    p_forest,
    width = 14.8,
    height = 7.8
  )
}

# -----------------------------
# 8. Plot 2: Adjusted mean delta clock age
# -----------------------------

if (nrow(means_plot) > 0) {
  
  p_means <- ggplot(
    means_plot,
    aes(
      x = medication_cluster_short,
      y = adjusted_mean,
      color = medication_cluster,
      fill = medication_cluster
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "#4D4D4D",
      linewidth = 0.6
    ) +
    geom_errorbar(
      aes(
        ymin = adjusted_ci_lo,
        ymax = adjusted_ci_hi
      ),
      width = 0.14,
      linewidth = 0.85,
      alpha = 0.95
    ) +
    geom_point(
      shape = 21,
      size = 3.3,
      stroke = 0.9
    ) +
    geom_text(
      aes(
        label = paste0("N=", comma(N_cluster))
      ),
      vjust = -1.25,
      size = 2.8,
      show.legend = FALSE
    ) +
    facet_wrap(
      ~ organ_label,
      ncol = 2,
      scales = "free_y"
    ) +
    scale_color_manual(values = cluster_palette, guide = "none") +
    scale_fill_manual(values = cluster_palette_light, guide = "none") +
    labs(
      x = NULL,
      y = "Adjusted mean delta clock age years",
      title = "Adjusted delta clock-age profiles by medication cluster",
      subtitle = paste0(
        "Predicted at mean covariate values from the organ-specific linear models"
      )
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.background = element_rect(fill = plot_bg, color = NA),
      panel.background = element_rect(fill = plot_bg, color = NA),
      strip.background = element_rect(fill = strip_bg, color = NA),
      strip.text = element_text(face = "bold", size = 12, hjust = 0),
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 11),
      axis.text.x = element_text(angle = 35, hjust = 1, face = "bold", size = 9.5),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 12),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  ggsave(
    file.path(out_dir, "medication_cluster_delta_clock_adjusted_means.pdf"),
    p_means,
    width = 11.5,
    height = 7.8
  )
  
  ggsave(
    file.path(out_dir, "medication_cluster_delta_clock_adjusted_means.png"),
    p_means,
    width = 11.5,
    height = 7.8,
    dpi = 300
  )
  
  ggsave(
    file.path(out_dir, "medication_cluster_delta_clock_adjusted_means.svg"),
    p_means,
    width = 11.5,
    height = 7.8
  )
}

# -----------------------------
# 9. Plot 3: Heatmap of beta effects
# -----------------------------

if (nrow(effects_plot) > 0) {
  
  heat_tbl <- effects_plot %>%
    mutate(
      beta_label = paste0(
        fmt_num(beta, 2),
        "\nBonf=", fmt_p(p_bonferroni)
      ),
      sig_symbol = case_when(
        bonf_sig ~ "\u25CF",        # solid circle
        nominal_only ~ "\u25D0",    # half-filled circle approximation
        TRUE ~ "\u25CB"             # empty circle
      ),
      heat_label = paste0(beta_label, "\n", sig_symbol)
    )
  
  p_heat <- ggplot(
    heat_tbl,
    aes(
      x = organ_label,
      y = exposure_cluster_short,
      fill = beta
    )
  ) +
    geom_tile(color = plot_bg, linewidth = 0.9) +
    geom_text(
      aes(label = heat_label),
      size = 3.0,
      lineheight = 0.88,
      color = "#1F1F1F"
    ) +
    scale_fill_gradient2(
      low = "#4E6FAE",
      mid = "#FFF9E8",
      high = "#B65E16",
      midpoint = 0,
      name = "Beta"
    ) +
    labs(
      x = NULL,
      y = NULL,
      title = "Medication-cluster effects on delta metabolomics mortality-clock age",
      subtitle = "\u25CF Bonferroni-significant; \u25D0 nominal only; \u25CB not significant"
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.background = element_rect(fill = plot_bg, color = NA),
      panel.background = element_rect(fill = plot_bg, color = NA),
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 11),
      axis.text.x = element_text(face = "bold", size = 11),
      axis.text.y = element_text(face = "bold", size = 10),
      legend.position = "right"
    )
  
  ggsave(
    file.path(out_dir, "medication_cluster_delta_clock_beta_heatmap_bonferroni.pdf"),
    p_heat,
    width = 9.2,
    height = 5.2
  )
  
  ggsave(
    file.path(out_dir, "medication_cluster_delta_clock_beta_heatmap_bonferroni.png"),
    p_heat,
    width = 9.2,
    height = 5.2,
    dpi = 300
  )
  
  ggsave(
    file.path(out_dir, "medication_cluster_delta_clock_beta_heatmap_bonferroni.svg"),
    p_heat,
    width = 9.2,
    height = 5.2
  )
}

# -----------------------------
# 10. Console summary
# -----------------------------

cat("\n============================================================\n")
cat("Finished Bonferroni-coded medication-cluster delta-clock plots.\n")
cat("Input results directory:\n", results_dir, "\n\n")
cat("Output plot directory:\n", out_dir, "\n\n")

cat("Bonferroni-significant signals:\n")
print(
  effects_plot %>%
    filter(bonf_sig) %>%
    arrange(p_bonferroni) %>%
    select(
      organ_label,
      exposure_cluster,
      N,
      N_exposure,
      beta,
      ci_lo,
      ci_hi,
      p,
      p_bonferroni
    )
)

cat("\nNominal-only signals:\n")
print(
  effects_plot %>%
    filter(nominal_only) %>%
    arrange(p) %>%
    select(
      organ_label,
      exposure_cluster,
      N,
      N_exposure,
      beta,
      ci_lo,
      ci_hi,
      p,
      p_bonferroni
    )
)

cat("============================================================\n")