# ============================================================
# APOE e4/e4 vs e2/e2: dementia L'EPOCH vs matched BAG plots
#
# Revised changes:
#   1) p_beta includes horizontal grid lines.
#   2) Significance threshold for both L'EPOCH and BAG:
#        P < 0.05 / 11 / 2
#      Circle = dementia L'EPOCH
#      Triangle = matched BAG
#      Filled symbol = significant
#      Empty symbol = non-significant
#   3) Van Gogh-inspired, colorblind-friendly palette for the
#      11 organ/modality clock-BAG pairs.
#   4) Fold labels are directly added to p_beta for pairs where
#      at least one of L'EPOCH or BAG is significant.
#      Fold label also denotes beta direction.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(forcats)
  library(scales)
  library(ggrepel)
})

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

result_dir <- "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/apoe_status_ukbb/apoe_e4e4_vs_e2e2_dementia_clock_results"

out_dir <- file.path(result_dir, "plots_bonferroni_half_alpha_filled_symbols_fold_labels")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Two-sided Bonferroni threshold requested:
# P < 0.05 / 11 / 2
sig_alpha <- 0.05 / 11 / 2

# ------------------------------------------------------------
# 2. Read all per-clock result files
# ------------------------------------------------------------

result_files <- list.files(
  result_dir,
  pattern = "dementia_.*_clock_acceleration_z_association_and_fold\\.tsv$",
  full.names = TRUE
)

stopifnot(length(result_files) == 11)

res <- purrr::map_dfr(result_files, readr::read_tsv, show_col_types = FALSE)

# ------------------------------------------------------------
# 3. Clean labels
# ------------------------------------------------------------

clean_lepoch_label <- function(x) {
  x %>%
    str_remove("^dementia_") %>%
    str_remove("_clock_acceleration_z$") %>%
    str_replace("_proteomics$", "_Prot") %>%
    str_replace("_metabolomics$", "_Met") %>%
    str_replace_all("_", " ") %>%
    str_to_title()
}

clean_bag_label <- function(x) {
  x %>%
    str_replace("_ProtBAG$", " ProtBAG") %>%
    str_replace("_MetBAG$", " MetBAG") %>%
    str_replace_all("_", " ")
}

beta_sign_label <- function(x) {
  case_when(
    is.na(x) ~ "NA",
    x > 0 ~ "+",
    x < 0 ~ "-",
    TRUE ~ "0"
  )
}

res2 <- res %>%
  mutate(
    clock_label = clean_lepoch_label(lepoch_clock),
    bag_label = clean_bag_label(matched_bag),
    
    # Significance for both dementia L'EPOCH and matched BAG.
    lepoch_sig = lepoch_std_p < sig_alpha,
    bag_sig = bag_std_p < sig_alpha,
    any_pair_sig = lepoch_sig | bag_sig,
    
    lepoch_direction = case_when(
      lepoch_std_beta > 0 ~ "Higher in APOE e4/e4",
      lepoch_std_beta < 0 ~ "Lower in APOE e4/e4",
      TRUE ~ "No difference"
    ),
    
    bag_direction = case_when(
      bag_std_beta > 0 ~ "Higher in APOE e4/e4",
      bag_std_beta < 0 ~ "Lower in APOE e4/e4",
      TRUE ~ "No difference"
    ),
    
    same_direction_recomputed = case_when(
      is.na(lepoch_std_beta) | is.na(bag_std_beta) ~ NA,
      lepoch_std_beta == 0 | bag_std_beta == 0 ~ FALSE,
      sign(lepoch_std_beta) == sign(bag_std_beta) ~ TRUE,
      TRUE ~ FALSE
    ),
    
    abs_lepoch_beta_std = abs(lepoch_std_beta),
    abs_bag_beta_std = abs(bag_std_beta),
    
    # L'EPOCH / BAG fold, as originally computed.
    lepoch_vs_bag_abs_effect_fold = abs_lepoch_beta_std / abs_bag_beta_std,
    
    # For direct annotation, report whichever marker has the larger absolute effect.
    stronger_marker_recomputed = case_when(
      is.na(abs_lepoch_beta_std) | is.na(abs_bag_beta_std) ~ "NA",
      abs_lepoch_beta_std > abs_bag_beta_std ~ "Clock",
      abs_bag_beta_std > abs_lepoch_beta_std ~ "BAG",
      TRUE ~ "Equal"
    ),
    
    stronger_fold_value = case_when(
      stronger_marker_recomputed == "Clock" ~ abs_lepoch_beta_std / abs_bag_beta_std,
      stronger_marker_recomputed == "BAG" ~ abs_bag_beta_std / abs_lepoch_beta_std,
      stronger_marker_recomputed == "Equal" ~ 1,
      TRUE ~ NA_real_
    ),
    
    lepoch_beta_sign = beta_sign_label(lepoch_std_beta),
    bag_beta_sign = beta_sign_label(bag_std_beta),
    
    direction_label = paste0("β ", lepoch_beta_sign, "/", bag_beta_sign),
    
    # Only show fold when at least one member of the pair is significant.
    # Label format examples:
    #   Clock 24.0x (β +/+)
    #   BAG 2.1x (β -/-)
    #   Clock 31.6x (β +/-)
    fold_label_for_plot = case_when(
      !any_pair_sig ~ "",
      is.na(stronger_fold_value) ~ "",
      stronger_marker_recomputed == "Equal" ~ paste0("Equal 1.0x (", direction_label, ")"),
      TRUE ~ paste0(
        stronger_marker_recomputed,
        " ",
        sprintf("%.1fx", stronger_fold_value),
        " (",
        direction_label,
        ")"
      )
    ),
    
    direction_status = case_when(
      same_direction_recomputed ~ "Same direction",
      !same_direction_recomputed ~ "Opposite direction",
      TRUE ~ "Unknown"
    )
  )

# ------------------------------------------------------------
# 4. Order rows by absolute L'EPOCH effect size
# ------------------------------------------------------------

clock_order <- res2 %>%
  arrange(abs(lepoch_std_beta)) %>%
  pull(clock_label)

res2 <- res2 %>%
  mutate(clock_label = factor(clock_label, levels = clock_order))

# ------------------------------------------------------------
# 5. Van Gogh-inspired colorblind-friendly palette
#    Same color is used for L'EPOCH and BAG within each pair.
# ------------------------------------------------------------

clock_levels <- levels(res2$clock_label)

vg_palette <- c(
  "#1F3B73",  # starry night deep blue
  "#2A6F9E",  # ultramarine blue
  "#0F766E",  # cypress teal
  "#4F772D",  # olive green
  "#DDAA33",  # sunflower yellow
  "#C65D35",  # burnt orange
  "#A23E48",  # red ochre
  "#6D597A",  # iris purple
  "#8C6D31",  # wheat brown
  "#3D5A80",  # blue grey
  "#B08968"   # warm beige
)

names(vg_palette) <- clock_levels

# ------------------------------------------------------------
# 6. Export clean summary table
# ------------------------------------------------------------

summary_tbl <- res2 %>%
  transmute(
    clock_label,
    lepoch_clock,
    matched_bag,
    
    lepoch_n = lepoch_std_n,
    lepoch_n_e2e2 = lepoch_std_n_e2e2,
    lepoch_n_e4e4 = lepoch_std_n_e4e4,
    
    lepoch_beta_std = lepoch_std_beta,
    lepoch_se_std = lepoch_std_se,
    lepoch_p = lepoch_std_p,
    lepoch_sig_threshold = sig_alpha,
    lepoch_sig,
    
    bag_beta_std = bag_std_beta,
    bag_se_std = bag_std_se,
    bag_p = bag_std_p,
    bag_sig_threshold = sig_alpha,
    bag_sig,
    
    any_pair_sig,
    same_direction = same_direction_recomputed,
    direction_status,
    
    lepoch_vs_bag_abs_effect_fold,
    stronger_marker_recomputed,
    stronger_fold_value,
    fold_label_for_plot
  ) %>%
  arrange(desc(abs(lepoch_beta_std)))

readr::write_tsv(
  summary_tbl,
  file.path(out_dir, "clean_apoe_e4e4_vs_e2e2_lepoch_vs_bag_summary_revised.tsv")
)

# ------------------------------------------------------------
# 7. Long-format table for beta comparison plot
# ------------------------------------------------------------

plot_df <- res2 %>%
  select(
    clock_label,
    lepoch_std_beta,
    lepoch_std_se,
    lepoch_std_p,
    lepoch_sig,
    bag_std_beta,
    bag_std_se,
    bag_std_p,
    bag_sig
  ) %>%
  pivot_longer(
    cols = c(lepoch_std_beta, bag_std_beta),
    names_to = "marker",
    values_to = "beta_std"
  ) %>%
  mutate(
    marker = recode(
      marker,
      "lepoch_std_beta" = "Dementia L'EPOCH",
      "bag_std_beta" = "Matched BAG"
    ),
    
    se_std = if_else(marker == "Dementia L'EPOCH", lepoch_std_se, bag_std_se),
    p_value = if_else(marker == "Dementia L'EPOCH", lepoch_std_p, bag_std_p),
    significant = if_else(marker == "Dementia L'EPOCH", lepoch_sig, bag_sig),
    
    ci_low = beta_std - 1.96 * se_std,
    ci_high = beta_std + 1.96 * se_std,
    
    marker = factor(marker, levels = c("Dementia L'EPOCH", "Matched BAG"))
  )

# Annotation table for fold labels on p_beta
fold_annot_df <- res2 %>%
  filter(any_pair_sig) %>%
  mutate(
    clock_label = factor(clock_label, levels = levels(res2$clock_label))
  )

# Dynamic x-axis range with room for fold labels on the right
x_min <- min(plot_df$ci_low, na.rm = TRUE)
x_max <- max(plot_df$ci_high, na.rm = TRUE)
x_range <- x_max - x_min
fold_x <- x_max + 0.08 * x_range

# ------------------------------------------------------------
# 8. Plot 1: standardized beta comparison
# ------------------------------------------------------------

p_beta <- ggplot(
  plot_df,
  aes(x = beta_std, y = clock_label)
) +
  # Horizontal grid lines, similar to previous figures
  geom_hline(
    yintercept = seq_along(levels(plot_df$clock_label)),
    color = "grey90",
    linewidth = 0.35
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.55,
    linetype = "dashed",
    color = "grey45"
  ) +
  
  # 95% CI lines
  geom_errorbarh(
    aes(
      xmin = ci_low,
      xmax = ci_high,
      color = clock_label
    ),
    height = 0.20,
    position = position_dodge(width = 0.65),
    linewidth = 0.50,
    alpha = 0.90
  ) +
  
  # Non-significant points: empty circle/triangle
  geom_point(
    data = plot_df %>% filter(!significant),
    aes(
      color = clock_label,
      shape = marker
    ),
    fill = "white",
    position = position_dodge(width = 0.65),
    size = 3.1,
    stroke = 0.95
  ) +
  
  # Significant points: filled circle/triangle
  geom_point(
    data = plot_df %>% filter(significant),
    aes(
      color = clock_label,
      fill = clock_label,
      shape = marker
    ),
    position = position_dodge(width = 0.65),
    size = 3.1,
    stroke = 0.95
  ) +
  
  # Fold labels for pairs with at least one significant marker
  geom_text(
    data = fold_annot_df,
    aes(
      x = fold_x,
      y = clock_label,
      label = fold_label_for_plot,
      color = clock_label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.2,
    fontface = "bold"
  ) +
  
  scale_shape_manual(
    values = c(
      "Dementia L'EPOCH" = 21,
      "Matched BAG" = 24
    )
  ) +
  scale_color_manual(values = vg_palette) +
  scale_fill_manual(values = vg_palette, guide = "none") +
  scale_x_continuous(
    name = "Adjusted standardized beta: APOE ε4/ε4 vs ε2/ε2",
    limits = c(x_min - 0.06 * x_range, fold_x + 0.38 * x_range),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_discrete(name = NULL) +
  coord_cartesian(clip = "off") +
  labs(
    title = "APOE ε4/ε4 vs ε2/ε2 effects on dementia L'EPOCH clocks and matched BAGs",
    subtitle = paste0(
      "Circle = dementia L'EPOCH; triangle = matched BAG. ",
      "Filled symbols indicate P < 0.05/11/2 = ",
      signif(sig_alpha, 3),
      "."
    ),
    color = "Organ/modality pair",
    shape = "Marker",
    caption = paste0(
      "Fold labels are shown only when at least one marker in the pair is significant. ",
      "Example: Clock 24.0x (β +/+) means the clock has a 24-fold larger absolute standardized beta than the matched BAG, with both betas positive."
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10.5),
    plot.caption = element_text(size = 8.8, hjust = 0),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 8.5),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 11, face = "bold"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 120, 8, 8)
  )

ggsave(
  filename = file.path(out_dir, "apoe_e4e4_vs_e2e2_standardized_beta_lepoch_vs_bag_revised.pdf"),
  plot = p_beta,
  width = 12.5,
  height = 6.8
)

ggsave(
  filename = file.path(out_dir, "apoe_e4e4_vs_e2e2_standardized_beta_lepoch_vs_bag_revised.png"),
  plot = p_beta,
  width = 12.5,
  height = 6.8,
  dpi = 300
)

# ------------------------------------------------------------
# 9. Plot 2: fold comparison for significant pairs
#    This plot uses the same fold labels shown in p_beta.
# ------------------------------------------------------------

fold_df_sig <- res2 %>%
  filter(any_pair_sig) %>%
  mutate(
    fold_plot_value = stronger_fold_value,
    clock_label = fct_reorder(clock_label, fold_plot_value),
    direction_status = factor(
      direction_status,
      levels = c("Same direction", "Opposite direction", "Unknown")
    )
  )

p_fold_sig <- ggplot(
  fold_df_sig,
  aes(
    x = fold_plot_value,
    y = clock_label,
    fill = direction_status
  )
) +
  geom_vline(
    xintercept = 1,
    linewidth = 0.55,
    linetype = "dashed",
    color = "grey45"
  ) +
  geom_col(width = 0.70, color = "grey25", linewidth = 0.25) +
  geom_text(
    aes(label = fold_label_for_plot),
    hjust = -0.05,
    size = 3.3,
    fontface = "bold"
  ) +
  scale_x_continuous(
    name = "Larger absolute standardized beta fold",
    trans = "log10",
    breaks = c(1, 2, 5, 10, 25, 50),
    labels = c("1", "2", "5", "10", "25", "50"),
    expand = expansion(mult = c(0.04, 0.34))
  ) +
  scale_y_discrete(name = NULL) +
  scale_fill_manual(
    values = c(
      "Same direction" = "#2A6F9E",
      "Opposite direction" = "#C65D35",
      "Unknown" = "grey65"
    )
  ) +
  labs(
    title = "Relative APOE effect size for significant L'EPOCH–BAG pairs",
    subtitle = paste0(
      "Included pairs have at least one marker with P < 0.05/11/2 = ",
      signif(sig_alpha, 3),
      ". Labels identify whether the clock or BAG has the larger absolute standardized beta."
    ),
    fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10.5),
    legend.position = "top",
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 11, face = "bold"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 80, 8, 8)
  ) +
  coord_cartesian(clip = "off")

ggsave(
  filename = file.path(out_dir, "apoe_e4e4_vs_e2e2_significant_pairs_effect_size_fold_revised.pdf"),
  plot = p_fold_sig,
  width = 10.5,
  height = 6.2
)

ggsave(
  filename = file.path(out_dir, "apoe_e4e4_vs_e2e2_significant_pairs_effect_size_fold_revised.png"),
  plot = p_fold_sig,
  width = 10.5,
  height = 6.2,
  dpi = 300
)

# ------------------------------------------------------------
# 10. Plot 3: signed L'EPOCH vs BAG effect-size scatter
# ------------------------------------------------------------

scatter_df <- res2 %>%
  mutate(
    clock_label = fct_reorder(clock_label, abs(lepoch_std_beta)),
    pair_sig_status = case_when(
      lepoch_sig & bag_sig ~ "Both significant",
      lepoch_sig & !bag_sig ~ "Clock significant only",
      !lepoch_sig & bag_sig ~ "BAG significant only",
      TRUE ~ "Neither significant"
    ),
    pair_sig_status = factor(
      pair_sig_status,
      levels = c(
        "Both significant",
        "Clock significant only",
        "BAG significant only",
        "Neither significant"
      )
    )
  )

p_scatter <- ggplot(
  scatter_df,
  aes(
    x = bag_std_beta,
    y = lepoch_std_beta,
    label = clock_label
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.4,
    linetype = "dashed",
    color = "grey55"
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.4,
    linetype = "dashed",
    color = "grey55"
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linewidth = 0.5,
    linetype = "dotted",
    color = "grey40"
  ) +
  geom_point(
    aes(
      color = clock_label,
      shape = pair_sig_status
    ),
    size = 3.2,
    stroke = 0.9
  ) +
  ggrepel::geom_text_repel(
    aes(color = clock_label),
    size = 3.2,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.25,
    segment.size = 0.25,
    show.legend = FALSE
  ) +
  scale_color_manual(values = vg_palette, guide = "none") +
  scale_shape_manual(
    values = c(
      "Both significant" = 16,
      "Clock significant only" = 21,
      "BAG significant only" = 24,
      "Neither significant" = 1
    )
  ) +
  coord_equal() +
  labs(
    title = "Signed APOE effect sizes for matched BAGs and dementia L'EPOCH clocks",
    subtitle = paste0(
      "Significance threshold for both markers: P < 0.05/11/2 = ",
      signif(sig_alpha, 3),
      "."
    ),
    x = "Matched BAG standardized beta",
    y = "Dementia L'EPOCH standardized beta",
    shape = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10.5),
    legend.position = "top",
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.35),
    panel.grid.major.x = element_line(color = "grey94", linewidth = 0.25),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(out_dir, "apoe_e4e4_vs_e2e2_signed_effect_scatter_revised.pdf"),
  plot = p_scatter,
  width = 7.8,
  height = 6.8
)

ggsave(
  filename = file.path(out_dir, "apoe_e4e4_vs_e2e2_signed_effect_scatter_revised.png"),
  plot = p_scatter,
  width = 7.8,
  height = 6.8,
  dpi = 300
)

# ------------------------------------------------------------
# 11. Print concise summaries
# ------------------------------------------------------------

message("Significance threshold for both L'EPOCH and BAG: P < ", signif(sig_alpha, 5))

message("Number of significant L'EPOCH clocks: ",
        sum(res2$lepoch_sig, na.rm = TRUE), " / 11")

message("Number of significant matched BAGs: ",
        sum(res2$bag_sig, na.rm = TRUE), " / 11")

message("Number of pairs with at least one significant marker: ",
        sum(res2$any_pair_sig, na.rm = TRUE), " / 11")

message("Outputs saved to: ", out_dir)