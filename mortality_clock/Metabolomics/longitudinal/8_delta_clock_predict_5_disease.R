# ============================================================
# Forest plot: longitudinal metabolomics delta mortality clocks
# versus algorithmically defined disease endpoints
#
# Van Gogh-inspired, colorblind-friendly palette
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

results_dir <- "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/mortality_clock/SA/delta_metabolomics_algorithmic_disease_onset"

out_dir <- file.path(results_dir, "forest_plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 2. Settings
# -----------------------------

organ_order <- c("Endocrine", "Digestive", "Hepatic", "Immune")

endpoint_order <- c(
  "All-cause dementia",
  "Asthma",
  "COPD",
  "Myocardial infarction",
  "Stroke"
)

# Van Gogh-inspired but not too bright:
# Sunflower ochre, iris blue, burnt sienna, cypress green
organ_palette <- c(
  "Endocrine" = "#C9A227",
  "Digestive" = "#4E6FAE",
  "Hepatic" = "#B65E16",
  "Immune" = "#3F7F5F"
)

organ_palette_light <- c(
  "Endocrine" = "#E7D48A",
  "Digestive" = "#A9BCE3",
  "Hepatic" = "#E2AA72",
  "Immune" = "#A7C7AE"
)

plot_bg <- "#FFF9E8"
strip_bg <- "#F3E7C8"
grid_col <- "#D8C99B"

bonf_p <- 0.05 / 4 / 5
bonf_label <- paste0("Bonferroni P < ", signif(bonf_p, 2))

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
# 4. Read and combine results
# -----------------------------

result_files <- list.files(
  results_dir,
  pattern = "^cox_delta_metabolomics_.*\\.tsv$",
  full.names = TRUE
)

result_files <- result_files[
  !grepl("analysis_dataset|complete_cases|run_summary", result_files)
]

if (length(result_files) == 0) {
  stop("No cox_delta_metabolomics_*.tsv files found in: ", results_dir)
}

res <- rbindlist(
  lapply(result_files, function(f) {
    x <- fread(f, sep = "\t", fill = TRUE)
    x[, source_file := basename(f)]
    x
  }),
  fill = TRUE
)

# -----------------------------
# 5. Clean numeric fields
# -----------------------------

num_cols <- c(
  "N",
  "N_case",
  "N_noncase",
  "N_prevalent_excluded",
  "event_rate",
  "delta_hr",
  "delta_ci_lo",
  "delta_ci_hi",
  "delta_p",
  "baseline_hr",
  "baseline_ci_lo",
  "baseline_ci_hi",
  "baseline_p",
  "reduced_cindex",
  "full_cindex",
  "delta_cindex_full_minus_reduced",
  "lrt_chisq_delta_vs_reduced",
  "lrt_p_delta_vs_reduced"
)

for (cc in intersect(num_cols, names(res))) {
  res[[cc]] <- safe_num(res[[cc]])
}

plot_tbl <- res %>%
  filter(status == "ok") %>%
  filter(
    !is.na(delta_hr),
    !is.na(delta_ci_lo),
    !is.na(delta_ci_hi),
    !is.na(delta_p)
  ) %>%
  mutate(
    endpoint_label = factor(endpoint_label, levels = endpoint_order),
    organ_label = factor(organ_label, levels = rev(organ_order)),
    
    bonferroni_significant = delta_p < bonf_p,
    nominal_significant = delta_p < 0.05 & delta_p >= bonf_p,
    
    significance_group = case_when(
      bonferroni_significant ~ bonf_label,
      nominal_significant ~ "Nominal P < 0.05",
      TRUE ~ "Not significant"
    ),
    
    significance_group = factor(
      significance_group,
      levels = c(
        bonf_label,
        "Nominal P < 0.05",
        "Not significant"
      )
    ),
    
    sig_size_key = ifelse(bonferroni_significant, "Bonferroni", "Other"),
    sig_line_key = ifelse(bonferroni_significant, "Bonferroni", "Other"),
    
    hr_label = paste0(
      "HR ", fmt_num(delta_hr, 2),
      " (", fmt_num(delta_ci_lo, 2),
      "-", fmt_num(delta_ci_hi, 2), ")"
    ),
    
    p_label = paste0(
      "P = ", fmt_p(delta_p),
      "; LRT P = ", fmt_p(lrt_p_delta_vs_reduced)
    ),
    
    n_label = paste0(
      "N = ", comma(N),
      "; cases = ", comma(N_case)
    ),
    
    cindex_label = paste0(
      "\u0394C = ",
      fmt_num(delta_cindex_full_minus_reduced, 4)
    ),
    
    annotation_label = paste0(
      n_label,
      "\n",
      hr_label,
      "\n",
      p_label
    )
  )

if (nrow(plot_tbl) == 0) {
  stop("No valid status == 'ok' rows with HR/CI/P values were found.")
}

# Save cleaned combined table
fwrite(
  plot_tbl,
  file.path(out_dir, "combined_delta_metabolomics_disease_onset_forest_plot_table.tsv"),
  sep = "\t"
)

# Save Bonferroni-significant rows
fwrite(
  plot_tbl %>% filter(bonferroni_significant),
  file.path(out_dir, "bonferroni_significant_delta_metabolomics_disease_onset.tsv"),
  sep = "\t"
)

# -----------------------------
# 6. Annotation positioning
# -----------------------------

x_min_data <- min(plot_tbl$delta_ci_lo, na.rm = TRUE)
x_max_data <- max(plot_tbl$delta_ci_hi, na.rm = TRUE)

x_min <- min(0.75, x_min_data * 0.85)
x_max <- max(2.10, x_max_data * 1.60)

plot_tbl <- plot_tbl %>%
  mutate(
    x_text = x_max_data * 1.08
  )

# -----------------------------
# 7. Main forest plot: all endpoints
# -----------------------------

shape_values <- setNames(
  c(21, 21, 1),
  c(bonf_label, "Nominal P < 0.05", "Not significant")
)

size_values <- c(
  "Bonferroni" = 3.9,
  "Other" = 3.1
)

linewidth_values <- c(
  "Bonferroni" = 1.35,
  "Other" = 0.90
)

p_all <- ggplot(
  plot_tbl,
  aes(y = organ_label)
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "#5B5B5B",
    linewidth = 0.75
  ) +
  geom_segment(
    aes(
      x = delta_ci_lo,
      xend = delta_ci_hi,
      yend = organ_label,
      color = organ_label,
      linewidth = sig_line_key
    ),
    alpha = 0.98
  ) +
  geom_point(
    aes(
      x = delta_hr,
      color = organ_label,
      fill = organ_label,
      shape = significance_group,
      size = sig_size_key
    ),
    stroke = 1.05
  ) +
  geom_text(
    aes(
      x = x_text,
      label = annotation_label,
      color = organ_label
    ),
    hjust = 0,
    size = 3.05,
    lineheight = 0.90,
    show.legend = FALSE
  ) +
  facet_wrap(
    ~ endpoint_label,
    ncol = 1,
    scales = "free_y"
  ) +
  scale_color_manual(values = organ_palette, name = "Delta clock") +
  scale_fill_manual(values = organ_palette_light, name = "Delta clock") +
  scale_shape_manual(
    values = shape_values,
    name = "Significance"
  ) +
  scale_size_manual(
    values = size_values,
    guide = "none"
  ) +
  scale_linewidth_manual(
    values = linewidth_values,
    guide = "none"
  ) +
  scale_x_log10(
    limits = c(x_min, x_max),
    breaks = c(0.75, 0.90, 1.00, 1.10, 1.25, 1.50, 2.00),
    labels = number_format(accuracy = 0.01)
  ) +
  coord_cartesian(clip = "off") +
  labs(
    x = paste0(
      "Hazard ratio per 1 SD higher \u0394 metabolomics mortality-clock biomarker\n",
      "Bonferroni threshold: P < 0.05 / 4 / 5 = ",
      signif(bonf_p, 3)
    ),
    y = NULL,
    title = "Longitudinal metabolomics mortality-clock change and incident disease onset",
    subtitle = paste0(
      "Cox models adjusted for age at instance 1, sex, smoking, BMI, blood pressure, ",
      "and baseline mortality-clock acceleration"
    )
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.background = element_rect(fill = plot_bg, color = NA),
    panel.background = element_rect(fill = plot_bg, color = NA),
    legend.background = element_rect(fill = plot_bg, color = NA),
    legend.box.background = element_rect(fill = plot_bg, color = NA),
    
    plot.title = element_text(face = "bold", size = 15, color = "#2B2B2B"),
    plot.subtitle = element_text(size = 11, color = "#3A3A3A"),
    
    strip.text = element_text(face = "bold", size = 12, hjust = 0, color = "#2B2B2B"),
    strip.background = element_rect(fill = strip_bg, color = NA),
    
    axis.text.y = element_text(face = "bold", size = 11, color = "#2B2B2B"),
    axis.text.x = element_text(size = 10, color = "#2B2B2B"),
    axis.title.x = element_text(size = 12, color = "#2B2B2B"),
    
    axis.line = element_line(color = "#2B2B2B", linewidth = 0.6),
    axis.ticks = element_line(color = "#2B2B2B"),
    
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    
    plot.margin = margin(10, 165, 10, 10)
  )

ggsave(
  file.path(out_dir, "delta_metabolomics_4clocks_5endpoints_forest_plot_vangogh.pdf"),
  p_all,
  width = 13.5,
  height = 11.0
)

ggsave(
  file.path(out_dir, "delta_metabolomics_4clocks_5endpoints_forest_plot_vangogh.png"),
  p_all,
  width = 13.5,
  height = 11.0,
  dpi = 300
)

ggsave(
  file.path(out_dir, "delta_metabolomics_4clocks_5endpoints_forest_plot_vangogh.svg"),
  p_all,
  width = 13.5,
  height = 11.0
)

# -----------------------------
# 8. Compact heatmap-style summary
# -----------------------------

heat_tbl <- plot_tbl %>%
  mutate(
    endpoint_label = factor(endpoint_label, levels = endpoint_order),
    organ_label = factor(organ_label, levels = organ_order),
    neg_log10_p = -log10(delta_p),
    label = paste0(
      "HR=", fmt_num(delta_hr, 2),
      "\nP=", fmt_p(delta_p)
    )
  )

p_heat <- ggplot(
  heat_tbl,
  aes(x = endpoint_label, y = organ_label, fill = neg_log10_p)
) +
  geom_tile(color = plot_bg, linewidth = 0.9) +
  geom_text(aes(label = label), size = 3.0, lineheight = 0.85, color = "#1F1F1F") +
  geom_point(
    data = heat_tbl %>% filter(bonferroni_significant),
    aes(x = endpoint_label, y = organ_label),
    inherit.aes = FALSE,
    shape = 8,
    size = 3.2,
    color = "#1F1F1F"
  ) +
  scale_fill_gradientn(
    colors = c("#F7EBC2", "#E7C85C", "#C98A00", "#4E6FAE", "#273A68"),
    name = expression(-log[10](P))
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Delta metabolomics mortality-clock associations across disease endpoints",
    subtitle = paste0(
      "Star marks Bonferroni-significant associations, P < ",
      signif(bonf_p, 3)
    )
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.background = element_rect(fill = plot_bg, color = NA),
    panel.background = element_rect(fill = plot_bg, color = NA),
    legend.background = element_rect(fill = plot_bg, color = NA),
    
    plot.title = element_text(face = "bold", color = "#2B2B2B"),
    plot.subtitle = element_text(size = 11, color = "#3A3A3A"),
    
    axis.text.x = element_text(angle = 35, hjust = 1, face = "bold", color = "#2B2B2B"),
    axis.text.y = element_text(face = "bold", color = "#2B2B2B"),
    
    axis.line = element_line(color = "#2B2B2B", linewidth = 0.6),
    axis.ticks = element_line(color = "#2B2B2B"),
    
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

ggsave(
  file.path(out_dir, "delta_metabolomics_4clocks_5endpoints_heatmap_summary_vangogh.pdf"),
  p_heat,
  width = 9.5,
  height = 4.8
)

ggsave(
  file.path(out_dir, "delta_metabolomics_4clocks_5endpoints_heatmap_summary_vangogh.png"),
  p_heat,
  width = 9.5,
  height = 4.8,
  dpi = 300
)

ggsave(
  file.path(out_dir, "delta_metabolomics_4clocks_5endpoints_heatmap_summary_vangogh.svg"),
  p_heat,
  width = 9.5,
  height = 4.8
)

# -----------------------------
# 9. Console summary
# -----------------------------

cat("\n============================================================\n")
cat("Finished plotting Van Gogh-style delta metabolomics disease-onset forest plots.\n")
cat("Input directory:\n", results_dir, "\n\n")
cat("Output directory:\n", out_dir, "\n\n")
cat("Bonferroni threshold:", bonf_p, "\n\n")
cat("Bonferroni-significant signals:\n")

print(
  plot_tbl %>%
    filter(bonferroni_significant) %>%
    arrange(delta_p) %>%
    select(
      endpoint_label,
      organ_label,
      N,
      N_case,
      delta_hr,
      delta_ci_lo,
      delta_ci_hi,
      delta_p,
      lrt_p_delta_vs_reduced,
      delta_cindex_full_minus_reduced
    )
)

cat("============================================================\n")