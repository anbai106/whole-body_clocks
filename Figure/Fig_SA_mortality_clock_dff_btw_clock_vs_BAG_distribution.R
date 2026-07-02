# ============================================================
# Three-panel modality summary:
# MRI, Proteomics, Metabolomics
#
# Fair comparison question:
#   Among disease endpoints where either the mortality clock OR the matched BAG
#   is significant, how often does the mortality clock outperform its matched BAG?
#
# Fair selection:
#   clock_p < 0.05 / number of unique diseases / number of mortality clocks
#      OR
#   bag_p   < 0.05 / number of unique diseases / number of mortality clocks
#
# Alternative:
#   Set fair_selection_mode <- "all_valid" to count all valid clock-BAG pairs.
#
# Classification:
#   1. No significant Clock-BAG difference:
#        joint_p_diff >= 0.05
#   2. Mortality clock overpowered:
#        joint_p_diff < 0.05 and beta_clock - beta_BAG > 0
#   3. BAG overpowered:
#        joint_p_diff < 0.05 and beta_clock - beta_BAG < 0
#
# Effect size:
#   joint_log_HR_diff_clock_minus_BAG =
#     beta_clock_joint - beta_BAG_joint
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(scales)
})

# ----------------------------
# 1. Input / output
# ----------------------------
infile <- "/Users/hao/Downloads/clock_vs_BAG_all_rows_all_status.xlsx"

out_dir <- "/Users/hao/Downloads"

out_prefix <- file.path(
  out_dir,
  "panel_modality_overpower_summary_FAIR_clock_or_BAG"
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 2. Analysis controls
# ----------------------------
diff_p_cutoff <- 0.05

# Fair comparison selection mode:
#   "clock_or_bag" = selected if clock_p OR bag_p passes Bonferroni threshold
#   "all_valid"    = all valid clock-BAG disease pairs, regardless of marginal P
fair_selection_mode <- "clock_or_bag"

# ----------------------------
# 3. Read data
# ----------------------------
df_raw <- readxl::read_excel(infile, sheet = 1)

required_cols <- c(
  "disease_id",
  "organ",
  "modality",
  "mortality_clock",
  "bag",
  "status",
  "clock_p",
  "bag_p",
  "joint_p_diff"
)

missing_cols <- setdiff(required_cols, names(df_raw))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns from input file: ",
    paste(missing_cols, collapse = ", ")
  )
}

optional_cols <- c(
  "joint_beta_diff_clock_minus_bag",
  "joint_se_diff",
  "joint_z_diff",
  "clock_joint_beta",
  "bag_joint_beta",
  "clock_joint_hr",
  "bag_joint_hr",
  "N",
  "N_case",
  "N_noncase",
  "clock_hr",
  "bag_hr",
  "clock_cindex",
  "bag_cindex",
  "delta_cindex_clock_minus_bag"
)

for (cc in optional_cols) {
  if (!cc %in% names(df_raw)) {
    df_raw[[cc]] <- NA_real_
  }
}

num_cols <- intersect(
  c(
    "N", "N_case", "N_noncase",
    "clock_beta", "clock_se", "clock_hr", "clock_ci_lo", "clock_ci_hi", "clock_p",
    "bag_beta", "bag_se", "bag_hr", "bag_ci_lo", "bag_ci_hi", "bag_p",
    "clock_joint_beta", "clock_joint_hr", "clock_joint_p",
    "bag_joint_beta", "bag_joint_hr", "bag_joint_p",
    "joint_beta_diff_clock_minus_bag",
    "joint_se_diff",
    "joint_z_diff",
    "joint_p_diff",
    "base_cindex", "clock_cindex", "bag_cindex", "both_cindex",
    "delta_cindex_clock_minus_bag",
    "delta_cindex_clock_minus_base",
    "delta_cindex_bag_minus_base"
  ),
  names(df_raw)
)

df <- df_raw %>%
  mutate(across(all_of(num_cols), ~ suppressWarnings(as.numeric(.x)))) %>%
  filter(status == "ok") %>%
  filter(!is.na(clock_p), !is.na(bag_p), !is.na(joint_p_diff))

# ----------------------------
# 4. Standardize modality and organ labels
# ----------------------------
modality_levels <- c("MRI", "Proteomics", "Metabolomics")

standardize_modality <- function(modality) {
  case_when(
    str_detect(modality, regex("^MRI$|mri", ignore_case = TRUE)) ~ "MRI",
    str_detect(modality, regex("proteomics|protein", ignore_case = TRUE)) ~ "Proteomics",
    str_detect(modality, regex("metabolomics|metabolite", ignore_case = TRUE)) ~ "Metabolomics",
    TRUE ~ as.character(modality)
  )
}

organ_levels <- c(
  "Brain",
  "Eye",
  "Pulmonary",
  "Heart",
  "Liver / hepatic",
  "Kidney / renal",
  "Pancreas",
  "Spleen",
  "Immune",
  "Endocrine",
  "Digestive",
  "Metabolic",
  "Adipose",
  "Skin",
  "Reproductive"
)

standardize_clock_block <- function(organ) {
  case_when(
    str_detect(organ, regex("brain", ignore_case = TRUE)) ~ "Brain",
    str_detect(organ, regex("eye", ignore_case = TRUE)) ~ "Eye",
    str_detect(organ, regex("pulmonary|lung", ignore_case = TRUE)) ~ "Pulmonary",
    str_detect(organ, regex("heart", ignore_case = TRUE)) ~ "Heart",
    str_detect(organ, regex("liver|hepatic", ignore_case = TRUE)) ~ "Liver / hepatic",
    str_detect(organ, regex("kidney|renal", ignore_case = TRUE)) ~ "Kidney / renal",
    str_detect(organ, regex("pancreas", ignore_case = TRUE)) ~ "Pancreas",
    str_detect(organ, regex("spleen", ignore_case = TRUE)) ~ "Spleen",
    str_detect(organ, regex("immune", ignore_case = TRUE)) ~ "Immune",
    str_detect(organ, regex("endocrine", ignore_case = TRUE)) ~ "Endocrine",
    str_detect(organ, regex("digestive", ignore_case = TRUE)) ~ "Digestive",
    str_detect(organ, regex("metabolic", ignore_case = TRUE)) ~ "Metabolic",
    str_detect(organ, regex("adipose", ignore_case = TRUE)) ~ "Adipose",
    str_detect(organ, regex("skin", ignore_case = TRUE)) ~ "Skin",
    str_detect(organ, regex("reproductive", ignore_case = TRUE)) ~ "Reproductive",
    TRUE ~ organ
  )
}

df <- df %>%
  mutate(
    modality_group = standardize_modality(modality),
    modality_group = factor(modality_group, levels = modality_levels),
    clock_block = standardize_clock_block(organ),
    clock_block = factor(clock_block, levels = organ_levels),
    clock_label = paste(organ, modality_group),
    clock_label = str_replace_all(clock_label, "_", " "),
    clock_label = str_squish(clock_label)
  ) %>%
  filter(!is.na(modality_group))

# ----------------------------
# 5. Construct standardized Clock - BAG effect size
# ----------------------------
df <- df %>%
  mutate(
    joint_log_HR_diff_clock_minus_BAG = case_when(
      !is.na(joint_beta_diff_clock_minus_bag) ~ joint_beta_diff_clock_minus_bag,

      is.na(joint_beta_diff_clock_minus_bag) &
        !is.na(clock_joint_beta) &
        !is.na(bag_joint_beta) ~ clock_joint_beta - bag_joint_beta,

      is.na(joint_beta_diff_clock_minus_bag) &
        !is.na(clock_joint_hr) &
        !is.na(bag_joint_hr) &
        clock_joint_hr > 0 &
        bag_joint_hr > 0 ~ log(clock_joint_hr) - log(bag_joint_hr),

      TRUE ~ NA_real_
    ),

    abs_joint_log_HR_diff_clock_minus_BAG = abs(joint_log_HR_diff_clock_minus_BAG),
    joint_HR_ratio_clock_vs_BAG = exp(joint_log_HR_diff_clock_minus_BAG)
  )

if (all(is.na(df$joint_log_HR_diff_clock_minus_BAG))) {
  stop(
    "Could not construct joint_log_HR_diff_clock_minus_BAG. ",
    "Need one of: joint_beta_diff_clock_minus_bag, ",
    "clock_joint_beta + bag_joint_beta, or clock_joint_hr + bag_joint_hr."
  )
}

# ----------------------------
# 6. Fair comparison threshold
# ----------------------------
n_unique_diseases <- n_distinct(df$disease_id)
n_mortality_clocks <- n_distinct(df$mortality_clock)

p_pair_bonf <- 0.05 / n_unique_diseases / n_mortality_clocks

threshold_tbl <- tibble(
  selection_rule = case_when(
    fair_selection_mode == "clock_or_bag" ~
      "clock_p < 0.05 / diseases / clocks OR bag_p < 0.05 / diseases / clocks",
    fair_selection_mode == "all_valid" ~
      "all valid clock-BAG disease pairs",
    TRUE ~ fair_selection_mode
  ),
  n_unique_diseases = n_unique_diseases,
  n_mortality_clocks = n_mortality_clocks,
  n_tests = n_unique_diseases * n_mortality_clocks,
  pairwise_bonferroni_threshold = p_pair_bonf,
  difference_test_p = "joint_p_diff",
  difference_p_cutoff = diff_p_cutoff,
  counted_population = fair_selection_mode
)

print(threshold_tbl)

readr::write_tsv(
  threshold_tbl,
  paste0(out_prefix, "_thresholds.tsv")
)

# ----------------------------
# 7. Fair selection and Clock vs BAG classification
# ----------------------------
df_fair_base <- df %>%
  mutate(
    clock_sig_bonf = clock_p < p_pair_bonf,
    bag_sig_bonf = bag_p < p_pair_bonf,

    selected_by = case_when(
      clock_sig_bonf & bag_sig_bonf ~ "Both clock and BAG significant",
      clock_sig_bonf & !bag_sig_bonf ~ "Clock only significant",
      !clock_sig_bonf & bag_sig_bonf ~ "BAG only significant",
      TRUE ~ "Neither significant"
    ),

    selected_by = factor(
      selected_by,
      levels = c(
        "Both clock and BAG significant",
        "Clock only significant",
        "BAG only significant",
        "Neither significant"
      )
    )
  )

if (fair_selection_mode == "clock_or_bag") {

  analysis_df <- df_fair_base %>%
    filter(clock_sig_bonf | bag_sig_bonf)

} else if (fair_selection_mode == "all_valid") {

  analysis_df <- df_fair_base %>%
    filter(
      !is.na(clock_p),
      !is.na(bag_p),
      !is.na(joint_p_diff),
      !is.na(joint_log_HR_diff_clock_minus_BAG)
    )

} else {

  stop("Unknown fair_selection_mode. Use 'clock_or_bag' or 'all_valid'.")
}

class_levels <- c(
  "No significant\nClock-BAG difference",
  "Mortality clock\noverpowered",
  "BAG\noverpowered"
)

analysis_df <- analysis_df %>%
  mutate(
    overpower_class = case_when(
      is.na(joint_p_diff) | is.na(joint_log_HR_diff_clock_minus_BAG) ~
        "No significant\nClock-BAG difference",

      joint_p_diff >= diff_p_cutoff ~
        "No significant\nClock-BAG difference",

      joint_p_diff < diff_p_cutoff &
        joint_log_HR_diff_clock_minus_BAG > 0 ~
        "Mortality clock\noverpowered",

      joint_p_diff < diff_p_cutoff &
        joint_log_HR_diff_clock_minus_BAG < 0 ~
        "BAG\noverpowered",

      TRUE ~ "No significant\nClock-BAG difference"
    ),

    overpower_class = factor(overpower_class, levels = class_levels),

    direction = case_when(
      overpower_class == "Mortality clock\noverpowered" ~ "Mortality clock stronger",
      overpower_class == "BAG\noverpowered" ~ "BAG stronger",
      TRUE ~ "No significant difference"
    )
  )

message("Fair selection mode: ", fair_selection_mode)
message("Number of selected fair-comparison disease associations: ", nrow(analysis_df))
message("Number selected by both clock and BAG: ", sum(analysis_df$selected_by == "Both clock and BAG significant", na.rm = TRUE))
message("Number selected by clock only: ", sum(analysis_df$selected_by == "Clock only significant", na.rm = TRUE))
message("Number selected by BAG only: ", sum(analysis_df$selected_by == "BAG only significant", na.rm = TRUE))
message("Number with no significant Clock-BAG difference: ", sum(analysis_df$overpower_class == "No significant\nClock-BAG difference", na.rm = TRUE))
message("Number with mortality-clock overpower: ", sum(analysis_df$overpower_class == "Mortality clock\noverpowered", na.rm = TRUE))
message("Number with BAG overpower: ", sum(analysis_df$overpower_class == "BAG\noverpowered", na.rm = TRUE))

# Save classified row-level data
readr::write_tsv(
  analysis_df,
  paste0(out_prefix, "_classified_clock_BAG_results_FAIR.tsv")
)

# ----------------------------
# 8. Modality-level counts
# ----------------------------
modality_counts <- analysis_df %>%
  count(modality_group, overpower_class, name = "n_results") %>%
  complete(
    modality_group = factor(modality_levels, levels = modality_levels),
    overpower_class = factor(class_levels, levels = class_levels),
    fill = list(n_results = 0)
  ) %>%
  group_by(modality_group) %>%
  mutate(
    total_results = sum(n_results),
    percent_results = if_else(total_results > 0, n_results / total_results, 0)
  ) %>%
  ungroup() %>%
  mutate(
    bar_label = paste0(
      comma(n_results),
      "\n(",
      percent(percent_results, accuracy = 0.1),
      ")"
    )
  )

readr::write_tsv(
  modality_counts,
  paste0(out_prefix, "_modality_counts_FAIR.tsv")
)

# ----------------------------
# 9. Clock-level counts
# ----------------------------
clock_counts <- analysis_df %>%
  count(
    modality_group,
    clock_block,
    clock_label,
    mortality_clock,
    selected_by,
    overpower_class,
    name = "n_results"
  ) %>%
  complete(
    nesting(modality_group, clock_block, clock_label, mortality_clock, selected_by),
    overpower_class = factor(class_levels, levels = class_levels),
    fill = list(n_results = 0)
  ) %>%
  group_by(modality_group, clock_label, mortality_clock) %>%
  mutate(
    total_results = sum(n_results),
    percent_results = if_else(total_results > 0, n_results / total_results, 0)
  ) %>%
  ungroup() %>%
  arrange(modality_group, clock_block, clock_label, selected_by, overpower_class)

readr::write_tsv(
  clock_counts,
  paste0(out_prefix, "_clock_level_counts_FAIR.tsv")
)

# A collapsed version for the stacked clock-level bar plot
clock_counts_collapsed <- analysis_df %>%
  count(
    modality_group,
    clock_block,
    clock_label,
    mortality_clock,
    overpower_class,
    name = "n_results"
  ) %>%
  complete(
    nesting(modality_group, clock_block, clock_label, mortality_clock),
    overpower_class = factor(class_levels, levels = class_levels),
    fill = list(n_results = 0)
  ) %>%
  group_by(modality_group, clock_label, mortality_clock) %>%
  mutate(
    total_results = sum(n_results),
    percent_results = if_else(total_results > 0, n_results / total_results, 0)
  ) %>%
  ungroup() %>%
  arrange(modality_group, clock_block, clock_label, overpower_class)

readr::write_tsv(
  clock_counts_collapsed,
  paste0(out_prefix, "_clock_level_counts_collapsed_FAIR.tsv")
)

# ----------------------------
# 10. Van Gogh-inspired colors
# ----------------------------
# Starry-night blue-gray, sunflower gold, cypress green.
class_palette <- c(
  "No significant\nClock-BAG difference" = "#5B6F95",
  "Mortality clock\noverpowered" = "#E3A21A",
  "BAG\noverpowered" = "#2F6B4F"
)

# Light parchment panel background inspired by Van Gogh paper/canvas tones.
plot_bg <- "#FBF4DF"
panel_bg <- "#FFF9E8"
grid_col <- "#D8CFAF"
ink_col <- "#1F2A44"

# ----------------------------
# 11. Three-panel modality bar plot
# ----------------------------
subtitle_txt <- paste0(
  "Counts use a fair selection rule: ",
  if_else(
    fair_selection_mode == "clock_or_bag",
    paste0(
      "clock P or BAG P < 0.05 / ",
      n_unique_diseases,
      " diseases / ",
      n_mortality_clocks,
      " clocks = ",
      signif(p_pair_bonf, 3)
    ),
    "all valid clock-BAG disease pairs"
  ),
  ". Clock-BAG difference threshold: joint_p_diff < ",
  diff_p_cutoff,
  "."
)

p_modality <- ggplot(
  modality_counts,
  aes(
    x = overpower_class,
    y = n_results,
    fill = overpower_class
  )
) +

  geom_col(
    width = 0.68,
    color = ink_col,
    linewidth = 0.28
  ) +

  geom_text(
    aes(label = bar_label),
    vjust = -0.22,
    size = 3.0,
    lineheight = 0.90,
    fontface = "bold",
    color = ink_col
  ) +

  facet_wrap(
    ~ modality_group,
    nrow = 1,
    scales = "free_y"
  ) +

  scale_fill_manual(
    values = class_palette,
    drop = FALSE
  ) +

  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.18))
  ) +

  labs(
    tag = "E",
    title = "Fair comparison of mortality clocks and matched BAGs by modality",
    subtitle = subtitle_txt,
    x = NULL,
    y = "Number of clock-disease results",
    caption = paste0(
      "Mortality clock overpowered: joint_p_diff < ",
      diff_p_cutoff,
      " and beta_clock - beta_BAG > 0. ",
      "BAG overpowered: joint_p_diff < ",
      diff_p_cutoff,
      " and beta_clock - beta_BAG < 0."
    )
  ) +

  guides(fill = "none") +

  theme_minimal(base_size = 10) +
  theme(
    plot.background = element_rect(fill = plot_bg, color = NA),
    panel.background = element_rect(fill = panel_bg, color = NA),

    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = grid_col, linewidth = 0.28),

    strip.text = element_text(
      face = "bold",
      size = 12,
      color = ink_col
    ),
    strip.background = element_rect(
      fill = "#E8D99B",
      color = "#B69E4A",
      linewidth = 0.40
    ),

    axis.text.x = element_text(
      size = 8.6,
      color = ink_col,
      face = "bold",
      lineheight = 0.92
    ),
    axis.text.y = element_text(
      size = 8.5,
      color = ink_col
    ),
    axis.title.y = element_text(
      size = 10.2,
      face = "bold",
      color = ink_col
    ),

    plot.tag = element_text(face = "bold", size = 24, color = ink_col),
    plot.tag.position = c(0.005, 0.995),

    plot.title = element_text(
      face = "bold",
      size = 14,
      color = ink_col
    ),
    plot.subtitle = element_text(
      size = 8.5,
      color = "#3A3A35",
      lineheight = 1.05
    ),
    plot.caption = element_text(
      size = 7.5,
      color = "#3A3A35",
      hjust = 0,
      lineheight = 1.05
    ),

    plot.margin = margin(8, 8, 8, 8)
  )

print(p_modality)

ggsave(
  filename = paste0(out_prefix, "_3panel_modality_counts_FAIR_vangogh.pdf"),
  plot = p_modality,
  width = 11.69,
  height = 4.8,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = paste0(out_prefix, "_3panel_modality_counts_FAIR_vangogh.png"),
  plot = p_modality,
  width = 11.69,
  height = 4.8,
  units = "in",
  dpi = 600,
  bg = plot_bg
)

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    filename = paste0(out_prefix, "_3panel_modality_counts_FAIR_vangogh.svg"),
    plot = p_modality,
    width = 11.69,
    height = 4.8,
    units = "in",
    device = svglite::svglite,
    bg = plot_bg
  )
}

# ----------------------------
# 12. Optional: clock-level stacked bar plot by modality
# ----------------------------
clock_order_tbl <- analysis_df %>%
  distinct(modality_group, clock_block, clock_label, mortality_clock) %>%
  arrange(modality_group, clock_block, clock_label) %>%
  mutate(clock_order = row_number())

clock_counts_plot <- clock_counts_collapsed %>%
  left_join(
    clock_order_tbl %>% select(mortality_clock, clock_order),
    by = "mortality_clock"
  ) %>%
  mutate(
    clock_label_plot = fct_reorder(clock_label, -clock_order)
  )

# Choose label color dynamically for readability:
# sunflower bars get dark labels; blue/green bars get cream labels.
clock_counts_plot <- clock_counts_plot %>%
  mutate(
    label_color = case_when(
      overpower_class == "Mortality clock\noverpowered" ~ ink_col,
      TRUE ~ "#FFF9E8"
    )
  )

p_clock <- ggplot(
  clock_counts_plot,
  aes(
    x = n_results,
    y = clock_label_plot,
    fill = overpower_class
  )
) +

  geom_col(
    width = 0.72,
    color = panel_bg,
    linewidth = 0.15
  ) +

  geom_text(
    data = clock_counts_plot %>% filter(n_results > 0),
    aes(
      label = n_results,
      color = label_color
    ),
    position = position_stack(vjust = 0.5),
    size = 2.45,
    fontface = "bold",
    show.legend = FALSE
  ) +

  facet_grid(
    modality_group ~ .,
    scales = "free_y",
    space = "free_y"
  ) +

  scale_fill_manual(
    values = class_palette,
    drop = FALSE,
    name = NULL
  ) +

  scale_color_identity() +

  scale_x_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.05))
  ) +

  labs(
    tag = "F",
    title = "Clock-level distribution of Clock-BAG differences by modality",
    subtitle = subtitle_txt,
    x = "Number of clock-disease results",
    y = NULL,
    caption = "Stacked bars show the same fair-selection classification, stratified by each mortality clock."
  ) +

  theme_minimal(base_size = 9) +
  theme(
    plot.background = element_rect(fill = plot_bg, color = NA),
    panel.background = element_rect(fill = panel_bg, color = NA),

    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = grid_col, linewidth = 0.28),

    strip.text.y = element_text(
      face = "bold",
      size = 10,
      angle = 0,
      color = ink_col
    ),
    strip.background = element_rect(
      fill = "#E8D99B",
      color = "#B69E4A",
      linewidth = 0.40
    ),

    axis.text.y = element_text(
      size = 7.2,
      color = ink_col
    ),
    axis.text.x = element_text(
      size = 8.0,
      color = ink_col
    ),
    axis.title.x = element_text(
      size = 9.5,
      face = "bold",
      color = ink_col
    ),

    legend.position = "top",
    legend.text = element_text(size = 8.0, color = ink_col),

    plot.tag = element_text(face = "bold", size = 22, color = ink_col),
    plot.tag.position = c(0.005, 0.995),

    plot.title = element_text(
      face = "bold",
      size = 13.2,
      color = ink_col
    ),
    plot.subtitle = element_text(
      size = 8.2,
      color = "#3A3A35",
      lineheight = 1.05
    ),
    plot.caption = element_text(
      size = 7.3,
      color = "#3A3A35",
      hjust = 0,
      lineheight = 1.05
    ),

    plot.margin = margin(8, 8, 8, 8)
  )

print(p_clock)

ggsave(
  filename = paste0(out_prefix, "_clock_level_stacked_counts_FAIR_vangogh.pdf"),
  plot = p_clock,
  width = 8.27,
  height = 11.69,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = paste0(out_prefix, "_clock_level_stacked_counts_FAIR_vangogh.png"),
  plot = p_clock,
  width = 8.27,
  height = 11.69,
  units = "in",
  dpi = 600,
  bg = plot_bg
)

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    filename = paste0(out_prefix, "_clock_level_stacked_counts_FAIR_vangogh.svg"),
    plot = p_clock,
    width = 8.27,
    height = 11.69,
    units = "in",
    device = svglite::svglite,
    bg = plot_bg
  )
}

# ----------------------------
# 13. Selection-basis summary plot
# ----------------------------
# This extra plot helps diagnose whether the comparison is driven by
# clock-only, BAG-only, or both-significant endpoints.

selection_counts <- analysis_df %>%
  count(modality_group, selected_by, name = "n_results") %>%
  complete(
    modality_group = factor(modality_levels, levels = modality_levels),
    selected_by = factor(
      c(
        "Both clock and BAG significant",
        "Clock only significant",
        "BAG only significant",
        "Neither significant"
      ),
      levels = c(
        "Both clock and BAG significant",
        "Clock only significant",
        "BAG only significant",
        "Neither significant"
      )
    ),
    fill = list(n_results = 0)
  ) %>%
  group_by(modality_group) %>%
  mutate(
    total_results = sum(n_results),
    percent_results = if_else(total_results > 0, n_results / total_results, 0),
    bar_label = paste0(comma(n_results), "\n", percent(percent_results, accuracy = 0.1))
  ) %>%
  ungroup()

selection_palette <- c(
  "Both clock and BAG significant" = "#2F6B4F",
  "Clock only significant" = "#E3A21A",
  "BAG only significant" = "#5B6F95",
  "Neither significant" = "#B9A77E"
)

p_selection <- ggplot(
  selection_counts,
  aes(
    x = selected_by,
    y = n_results,
    fill = selected_by
  )
) +
  geom_col(
    width = 0.68,
    color = ink_col,
    linewidth = 0.28
  ) +
  geom_text(
    aes(label = bar_label),
    vjust = -0.22,
    size = 2.8,
    lineheight = 0.90,
    fontface = "bold",
    color = ink_col
  ) +
  facet_wrap(
    ~ modality_group,
    nrow = 1,
    scales = "free_y"
  ) +
  scale_fill_manual(values = selection_palette, drop = FALSE) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.18))
  ) +
  labs(
    tag = "G",
    title = "Selection basis under the fair clock-or-BAG rule",
    subtitle = paste0(
      "Bonferroni threshold: P < ",
      signif(p_pair_bonf, 3),
      " for either mortality clock or matched BAG."
    ),
    x = NULL,
    y = "Number of selected results",
    caption = "This diagnostic panel shows whether selected endpoints were driven by clock significance, BAG significance, or both."
  ) +
  guides(fill = "none") +
  theme_minimal(base_size = 10) +
  theme(
    plot.background = element_rect(fill = plot_bg, color = NA),
    panel.background = element_rect(fill = panel_bg, color = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = grid_col, linewidth = 0.28),
    strip.text = element_text(face = "bold", size = 12, color = ink_col),
    strip.background = element_rect(fill = "#E8D99B", color = "#B69E4A", linewidth = 0.40),
    axis.text.x = element_text(size = 7.4, color = ink_col, face = "bold", lineheight = 0.90),
    axis.text.y = element_text(size = 8.5, color = ink_col),
    axis.title.y = element_text(size = 10.2, face = "bold", color = ink_col),
    plot.tag = element_text(face = "bold", size = 24, color = ink_col),
    plot.tag.position = c(0.005, 0.995),
    plot.title = element_text(face = "bold", size = 14, color = ink_col),
    plot.subtitle = element_text(size = 8.5, color = "#3A3A35", lineheight = 1.05),
    plot.caption = element_text(size = 7.5, color = "#3A3A35", hjust = 0, lineheight = 1.05),
    plot.margin = margin(8, 8, 8, 8)
  )

print(p_selection)

ggsave(
  filename = paste0(out_prefix, "_selection_basis_counts_FAIR_vangogh.pdf"),
  plot = p_selection,
  width = 11.69,
  height = 4.8,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = paste0(out_prefix, "_selection_basis_counts_FAIR_vangogh.png"),
  plot = p_selection,
  width = 11.69,
  height = 4.8,
  units = "in",
  dpi = 600,
  bg = plot_bg
)

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    filename = paste0(out_prefix, "_selection_basis_counts_FAIR_vangogh.svg"),
    plot = p_selection,
    width = 11.69,
    height = 4.8,
    units = "in",
    device = svglite::svglite,
    bg = plot_bg
  )
}

message("Saved files:")
message("  ", paste0(out_prefix, "_thresholds.tsv"))
message("  ", paste0(out_prefix, "_classified_clock_BAG_results_FAIR.tsv"))
message("  ", paste0(out_prefix, "_modality_counts_FAIR.tsv"))
message("  ", paste0(out_prefix, "_clock_level_counts_FAIR.tsv"))
message("  ", paste0(out_prefix, "_clock_level_counts_collapsed_FAIR.tsv"))
message("  ", paste0(out_prefix, "_3panel_modality_counts_FAIR_vangogh.pdf/png/svg"))
message("  ", paste0(out_prefix, "_clock_level_stacked_counts_FAIR_vangogh.pdf/png/svg"))
message("  ", paste0(out_prefix, "_selection_basis_counts_FAIR_vangogh.pdf/png/svg"))