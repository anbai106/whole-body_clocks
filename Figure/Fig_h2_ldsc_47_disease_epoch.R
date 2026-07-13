#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(forcats)
  library(scales)
})

# ============================================================
# Paths
# ============================================================

INPUT_FILE <- "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/Result/LDSC_h2_intercept_47_disease_clocks.tsv"

OUTDIR <- "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/Figure/ldsc_h2_47_disease_epoch_clocks_organ_colors"

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Settings
# ============================================================

disease_order <- c("Asthma", "COPD", "Dementia", "MI", "Stroke")
modality_order <- c("MRI", "Proteomics", "Metabolomics")

# Better colorblind-friendly, Van Gogh / Okabe-Ito-inspired organ palette.
# The same organ/system uses the same color across all disease panels.
organ_colors <- c(
  "Brain"                 = "#0072B2",  # deep blue
  "Heart"                 = "#D55E00",  # vermilion
  "Hepatic"               = "#009E73",  # green
  "Immune"                = "#CC79A7",  # reddish purple
  "Endocrine"             = "#E69F00",  # sunflower orange
  "Digestive"             = "#56B4E9",  # sky blue
  "Metabolic"             = "#F0E442",  # yellow
  "Pulmonary"             = "#7B61A8",  # iris purple
  "Spleen"                = "#8B5A2B",  # earthy brown
  "Reproductive female"   = "#E78AC3",  # soft magenta
  "Reproductive male"     = "#4D4D4D",  # dark grey
  "Other"                 = "#999999"
)

# Filled plotting symbols.
# These allow fill = organ and shape = modality at the same time.
modality_shapes <- c(
  "MRI"           = 21,  # filled circle
  "Proteomics"   = 24,  # filled triangle
  "Metabolomics" = 22   # filled square
)

# ============================================================
# Helper functions
# ============================================================

clean_modality <- function(x) {
  x <- as.character(x)
  
  case_when(
    str_detect(x, regex("^mri$", ignore_case = TRUE)) ~ "MRI",
    str_detect(x, regex("proteomics", ignore_case = TRUE)) ~ "Proteomics",
    str_detect(x, regex("metabolomics", ignore_case = TRUE)) ~ "Metabolomics",
    TRUE ~ x
  )
}

clean_disease <- function(x) {
  x <- as.character(x)
  
  case_when(
    str_detect(x, regex("^asthma$", ignore_case = TRUE)) ~ "Asthma",
    str_detect(x, regex("^copd$", ignore_case = TRUE)) ~ "COPD",
    str_detect(x, regex("^dementia$", ignore_case = TRUE)) ~ "Dementia",
    str_detect(x, regex("^mi$", ignore_case = TRUE)) ~ "MI",
    str_detect(x, regex("^stroke$", ignore_case = TRUE)) ~ "Stroke",
    TRUE ~ x
  )
}

clean_organ <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "_", " ")
  x <- str_squish(x)
  
  case_when(
    str_detect(x, regex("brain", ignore_case = TRUE)) ~ "Brain",
    str_detect(x, regex("heart", ignore_case = TRUE)) ~ "Heart",
    str_detect(x, regex("hepatic|liver", ignore_case = TRUE)) ~ "Hepatic",
    str_detect(x, regex("immune", ignore_case = TRUE)) ~ "Immune",
    str_detect(x, regex("endocrine", ignore_case = TRUE)) ~ "Endocrine",
    str_detect(x, regex("digestive", ignore_case = TRUE)) ~ "Digestive",
    str_detect(x, regex("metabolic", ignore_case = TRUE)) ~ "Metabolic",
    str_detect(x, regex("pulmonary|lung", ignore_case = TRUE)) ~ "Pulmonary",
    str_detect(x, regex("spleen", ignore_case = TRUE)) ~ "Spleen",
    str_detect(x, regex("reproductive female", ignore_case = TRUE)) ~ "Reproductive female",
    str_detect(x, regex("reproductive male", ignore_case = TRUE)) ~ "Reproductive male",
    TRUE ~ "Other"
  )
}

short_modality <- function(x) {
  case_when(
    x == "MRI" ~ "MRI",
    x == "Proteomics" ~ "Prot.",
    x == "Metabolomics" ~ "Metab.",
    TRUE ~ as.character(x)
  )
}

make_clock_label <- function(organ, modality) {
  paste0(organ, "\n", short_modality(modality))
}

# ============================================================
# Read and QC data
# ============================================================

if (!file.exists(INPUT_FILE)) {
  stop("Cannot find input file: ", INPUT_FILE)
}

df <- fread(INPUT_FILE) %>% as.data.frame()

required_cols <- c(
  "clock_folder",
  "disease_label",
  "modality",
  "organ_label",
  "h2_mean",
  "h2_std",
  "h2_p_calc",
  "intercept",
  "lambda_gc",
  "mean_chi2",
  "h2_significant_bonferroni_47"
)

missing_cols <- setdiff(required_cols, names(df))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

plot_df <- df %>%
  mutate(
    h2_mean = as.numeric(h2_mean),
    h2_std = as.numeric(h2_std),
    h2_p_calc = as.numeric(h2_p_calc),
    intercept = as.numeric(intercept),
    lambda_gc = as.numeric(lambda_gc),
    mean_chi2 = as.numeric(mean_chi2),
    
    disease_label = clean_disease(disease_label),
    modality = clean_modality(modality),
    organ_color_group = clean_organ(organ_label),
    
    disease_label = factor(disease_label, levels = disease_order),
    modality = factor(modality, levels = modality_order),
    
    h2_lower = pmax(0, h2_mean - h2_std),
    h2_upper = h2_mean + h2_std,
    
    h2_sig_label = case_when(
      h2_p_calc < 0.05 / 47 ~ "Bonferroni",
      h2_p_calc < 0.05 ~ "Nominal",
      TRUE ~ "NS"
    ),
    
    clock_plot = make_clock_label(organ_color_group, modality)
  ) %>%
  filter(
    !is.na(disease_label),
    !is.na(modality),
    is.finite(h2_mean),
    is.finite(h2_std)
  )

# Add colors for any unexpected organ labels.
missing_color_organs <- setdiff(unique(plot_df$organ_color_group), names(organ_colors))

if (length(missing_color_organs) > 0) {
  extra_colors <- rep("#999999", length(missing_color_organs))
  names(extra_colors) <- missing_color_organs
  organ_colors <- c(organ_colors, extra_colors)
}

organ_order <- names(organ_colors)[names(organ_colors) %in% unique(plot_df$organ_color_group)]

plot_df <- plot_df %>%
  mutate(
    organ_color_group = factor(organ_color_group, levels = organ_order)
  )

# Order clocks within each disease panel by modality and descending h2.
plot_df <- plot_df %>%
  group_by(disease_label) %>%
  arrange(modality, desc(h2_mean), .by_group = TRUE) %>%
  mutate(
    clock_plot_ordered = factor(clock_plot, levels = unique(clock_plot))
  ) %>%
  ungroup()

# ============================================================
# Save plotting data and summaries
# ============================================================

fwrite(
  plot_df,
  file.path(OUTDIR, "LDSC_h2_47_disease_epoch_clocks_plot_data_organ_colors.tsv"),
  sep = "\t"
)

disease_summary <- plot_df %>%
  group_by(disease_label) %>%
  summarise(
    n_clocks = n(),
    mean_h2 = mean(h2_mean, na.rm = TRUE),
    median_h2 = median(h2_mean, na.rm = TRUE),
    min_h2 = min(h2_mean, na.rm = TRUE),
    max_h2 = max(h2_mean, na.rm = TRUE),
    mean_intercept = mean(intercept, na.rm = TRUE),
    mean_lambda_gc = mean(lambda_gc, na.rm = TRUE),
    mean_chi2 = mean(mean_chi2, na.rm = TRUE),
    n_bonferroni = sum(h2_p_calc < 0.05 / 47, na.rm = TRUE),
    .groups = "drop"
  )

organ_summary <- plot_df %>%
  group_by(organ_color_group) %>%
  summarise(
    n_clocks = n(),
    mean_h2 = mean(h2_mean, na.rm = TRUE),
    median_h2 = median(h2_mean, na.rm = TRUE),
    min_h2 = min(h2_mean, na.rm = TRUE),
    max_h2 = max(h2_mean, na.rm = TRUE),
    n_bonferroni = sum(h2_p_calc < 0.05 / 47, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_h2))

modality_summary <- plot_df %>%
  group_by(modality) %>%
  summarise(
    n_clocks = n(),
    mean_h2 = mean(h2_mean, na.rm = TRUE),
    median_h2 = median(h2_mean, na.rm = TRUE),
    min_h2 = min(h2_mean, na.rm = TRUE),
    max_h2 = max(h2_mean, na.rm = TRUE),
    mean_intercept = mean(intercept, na.rm = TRUE),
    mean_lambda_gc = mean(lambda_gc, na.rm = TRUE),
    mean_chi2 = mean(mean_chi2, na.rm = TRUE),
    n_bonferroni = sum(h2_p_calc < 0.05 / 47, na.rm = TRUE),
    .groups = "drop"
  )

fwrite(
  disease_summary,
  file.path(OUTDIR, "LDSC_h2_47_disease_epoch_clocks_disease_summary.tsv"),
  sep = "\t"
)

fwrite(
  organ_summary,
  file.path(OUTDIR, "LDSC_h2_47_disease_epoch_clocks_organ_summary.tsv"),
  sep = "\t"
)

fwrite(
  modality_summary,
  file.path(OUTDIR, "LDSC_h2_47_disease_epoch_clocks_modality_summary.tsv"),
  sep = "\t"
)

# ============================================================
# Main figure: disease panels, organ colors
# ============================================================

p_main <- ggplot(
  plot_df,
  aes(
    x = clock_plot_ordered,
    y = h2_mean,
    fill = organ_color_group,
    color = organ_color_group,
    shape = modality
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.35,
    color = "grey45"
  ) +
  geom_errorbar(
    aes(ymin = h2_lower, ymax = h2_upper),
    width = 0.18,
    linewidth = 0.50,
    alpha = 0.85
  ) +
  geom_point(
    size = 3.1,
    stroke = 0.70,
    alpha = 0.96
  ) +
  facet_wrap(
    ~ disease_label,
    scales = "free_x",
    nrow = 1
  ) +
  scale_fill_manual(
    values = organ_colors,
    drop = FALSE,
    name = "Organ/system"
  ) +
  scale_color_manual(
    values = organ_colors,
    drop = FALSE,
    name = "Organ/system"
  ) +
  scale_shape_manual(
    values = modality_shapes,
    drop = FALSE,
    name = "Modality"
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0.02, 0.12)),
    labels = number_format(accuracy = 0.01)
  ) +
  labs(
    x = NULL,
    y = expression("LDSC SNP heritability (" * h^2 * ")")
  ) +
  guides(
    color = guide_legend(
      title = "Organ/system",
      override.aes = list(shape = 21, size = 3.6, alpha = 1),
      nrow = 2,
      byrow = TRUE
    ),
    fill = "none",
    shape = guide_legend(
      title = "Modality",
      override.aes = list(fill = "grey80", color = "black", size = 3.6),
      nrow = 1
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.25, color = "grey88"),
    
    strip.background = element_rect(fill = "#F7F3E8", color = "black", linewidth = 0.8),
    strip.text = element_text(face = "bold", size = 11, color = "black"),
    
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8.2),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(face = "bold", size = 12),
    
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 9.5),
    legend.key.size = unit(0.45, "cm"),
    
    plot.margin = margin(6, 8, 6, 8)
  )

pdf_main <- file.path(
  OUTDIR,
  "LDSC_h2_47_disease_epoch_clocks_disease_panels_same_organ_colors.pdf"
)

png_main <- file.path(
  OUTDIR,
  "LDSC_h2_47_disease_epoch_clocks_disease_panels_same_organ_colors.png"
)

ggsave(pdf_main, p_main, width = 16.2, height = 6.0, device = cairo_pdf)
ggsave(png_main, p_main, width = 16.2, height = 6.0, dpi = 300)

# ============================================================
# Alternative ranked figure: all 47 clocks
# ============================================================

rank_df <- plot_df %>%
  arrange(h2_mean) %>%
  mutate(
    clock_rank_label = paste0(
      disease_label,
      " | ",
      organ_color_group,
      " ",
      modality
    ),
    clock_rank_label = factor(clock_rank_label, levels = clock_rank_label)
  )

p_ranked <- ggplot(
  rank_df,
  aes(
    x = h2_mean,
    y = clock_rank_label,
    fill = organ_color_group,
    color = organ_color_group,
    shape = modality
  )
) +
  geom_errorbarh(
    aes(xmin = h2_lower, xmax = h2_upper),
    height = 0.18,
    linewidth = 0.45,
    alpha = 0.80
  ) +
  geom_point(
    size = 2.9,
    stroke = 0.70,
    alpha = 0.96
  ) +
  scale_fill_manual(
    values = organ_colors,
    drop = FALSE,
    name = "Organ/system"
  ) +
  scale_color_manual(
    values = organ_colors,
    drop = FALSE,
    name = "Organ/system"
  ) +
  scale_shape_manual(
    values = modality_shapes,
    drop = FALSE,
    name = "Modality"
  ) +
  scale_x_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0.02, 0.08)),
    labels = number_format(accuracy = 0.01)
  ) +
  labs(
    x = expression("LDSC SNP heritability (" * h^2 * ")"),
    y = NULL
  ) +
  guides(
    color = guide_legend(
      title = "Organ/system",
      override.aes = list(shape = 21, size = 3.5, alpha = 1),
      nrow = 3,
      byrow = TRUE
    ),
    fill = "none",
    shape = guide_legend(
      title = "Modality",
      override.aes = list(fill = "grey80", color = "black", size = 3.5),
      nrow = 1
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.25, color = "grey88"),
    
    axis.text.y = element_text(size = 7.4),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(face = "bold", size = 12),
    
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 9.2),
    legend.key.size = unit(0.42, "cm"),
    
    plot.margin = margin(6, 8, 6, 8)
  )

pdf_ranked <- file.path(
  OUTDIR,
  "LDSC_h2_47_disease_epoch_clocks_ranked_same_organ_colors.pdf"
)

png_ranked <- file.path(
  OUTDIR,
  "LDSC_h2_47_disease_epoch_clocks_ranked_same_organ_colors.png"
)

ggsave(pdf_ranked, p_ranked, width = 9.0, height = 11.8, device = cairo_pdf)
ggsave(png_ranked, p_ranked, width = 9.0, height = 11.8, dpi = 300)

# ============================================================
# Additional compact disease-by-organ summary plot
# ============================================================

summary_plot_df <- plot_df %>%
  group_by(disease_label, organ_color_group) %>%
  summarise(
    mean_h2 = mean(h2_mean, na.rm = TRUE),
    se_h2 = sd(h2_mean, na.rm = TRUE) / sqrt(n()),
    n_clocks = n(),
    .groups = "drop"
  ) %>%
  mutate(
    h2_lower = pmax(0, mean_h2 - se_h2),
    h2_upper = mean_h2 + se_h2
  )

p_summary <- ggplot(
  summary_plot_df,
  aes(
    x = organ_color_group,
    y = mean_h2,
    fill = organ_color_group,
    color = organ_color_group
  )
) +
  geom_col(
    width = 0.72,
    alpha = 0.88
  ) +
  geom_errorbar(
    aes(ymin = h2_lower, ymax = h2_upper),
    width = 0.18,
    linewidth = 0.45,
    color = "grey25"
  ) +
  facet_wrap(
    ~ disease_label,
    scales = "free_x",
    nrow = 1
  ) +
  scale_fill_manual(values = organ_colors, drop = FALSE) +
  scale_color_manual(values = organ_colors, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0.02, 0.12)),
    labels = number_format(accuracy = 0.01)
  ) +
  labs(
    x = NULL,
    y = expression("Mean LDSC SNP heritability (" * h^2 * ")")
  ) +
  guides(
    fill = "none",
    color = "none"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.25, color = "grey88"),
    
    strip.background = element_rect(fill = "#F7F3E8", color = "black", linewidth = 0.8),
    strip.text = element_text(face = "bold", size = 11),
    
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8.5),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(face = "bold", size = 12),
    
    plot.margin = margin(6, 8, 6, 8)
  )

pdf_summary <- file.path(
  OUTDIR,
  "LDSC_h2_47_disease_epoch_clocks_mean_by_organ_disease_panels.pdf"
)

png_summary <- file.path(
  OUTDIR,
  "LDSC_h2_47_disease_epoch_clocks_mean_by_organ_disease_panels.png"
)

ggsave(pdf_summary, p_summary, width = 15.5, height = 5.0, device = cairo_pdf)
ggsave(png_summary, p_summary, width = 15.5, height = 5.0, dpi = 300)

# ============================================================
# Console output
# ============================================================

message("Finished plotting LDSC h2 for 47 disease EPOCH clocks using same organ colors.")
message("Input file:")
message(INPUT_FILE)
message("")
message("Main disease-panel figure:")
message(pdf_main)
message(png_main)
message("")
message("Ranked all-clock figure:")
message(pdf_ranked)
message(png_ranked)
message("")
message("Mean organ-by-disease summary figure:")
message(pdf_summary)
message(png_summary)
message("")
message("Output directory:")
message(OUTDIR)