#!/usr/bin/env Rscript

# ============================================================
# Fig: Disease L'EPOCH clocks predict future all-cause mortality
#
# Robustly collects:
#   <base_dir>/<clock_folder>/survival_analysis_mortality/
#      *_mortality_survival_summary.tsv
#
# Main effect size:
#   HR per 1-SD higher disease-clock acceleration-z
#
# Revised figure rules:
#   i)  Do NOT annotate P-values.
#   ii) Use solid filled circles only for Bonferroni-significant
#       signals: P < 0.05 / 47 / 5.
#   iii) Annotate number of mortality cases/deaths instead.
#
# Outputs:
#   all_disease_clock_mortality_survival_summary_merged.tsv
#   all_disease_clock_mortality_hr_forest.pdf/png
#   all_disease_clock_mortality_hr_heatmap.pdf/png
#   all_disease_clock_mortality_hr_ranked_lollipop.pdf/png
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(forcats)
  library(grid)
})

# ============================================================
# 1. User settings
# ============================================================

base_dir_candidates <- c(
  "/gpfs/fs001/cbica/home/wenju/Reproducibile_paper/WholeBodyClock",
  "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock",
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock"
)

base_dir <- Sys.getenv("WHOLEBODYCLOCK_BASE_DIR", unset = NA_character_)

if (is.na(base_dir) || base_dir == "") {
  base_dir <- base_dir_candidates[dir.exists(base_dir_candidates)][1]
}

if (is.na(base_dir) || !dir.exists(base_dir)) {
  stop(
    "Could not find base_dir. Please set it manually, e.g.\n",
    "base_dir <- '/gpfs/fs001/cbica/home/wenju/Reproducibile_paper/WholeBodyClock'"
  )
}

base_dir <- normalizePath(base_dir, mustWork = TRUE)

survival_dir_name <- "survival_analysis_mortality"

out_dir <- Sys.getenv(
  "MORTALITY_FIGURE_OUT_DIR",
  unset = file.path(base_dir, "all_disease_clock_mortality_survival_plots")
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

disease_order <- c("asthma", "dementia", "copd", "mi", "stroke")

disease_labels <- c(
  asthma   = "Asthma",
  dementia = "Dementia",
  copd     = "COPD",
  mi       = "MI",
  stroke   = "Stroke"
)

modality_order <- c("MRI", "Proteomics", "Metabolomics")

modality_cols <- c(
  MRI          = "#3B6EA8",
  Proteomics   = "#C65A1E",
  Metabolomics = "#2A9D8F",
  Other        = "#666666"
)

disease_cols <- c(
  Asthma   = "#4C78A8",
  Dementia = "#B279A2",
  COPD     = "#59A14F",
  MI       = "#E15759",
  Stroke   = "#F28E2B"
)

# Bonferroni threshold requested by user.
n_clocks_for_bonferroni <- 47
n_endpoint_families_for_bonferroni <- 5
bonferroni_p_threshold <- 0.05 / n_clocks_for_bonferroni / n_endpoint_families_for_bonferroni

message("Base directory:")
message("  ", base_dir)
message("Output directory:")
message("  ", out_dir)
message("Bonferroni threshold:")
message("  0.05 / 47 / 5 = ", signif(bonferroni_p_threshold, 4))

# ============================================================
# 2. Helper functions
# ============================================================

format_hr_ci <- function(hr, lo, hi) {
  case_when(
    is.na(hr) ~ "NA",
    is.na(lo) | is.na(hi) ~ sprintf("%.2f", hr),
    TRUE ~ sprintf("%.2f (%.2f-%.2f)", hr, lo, hi)
  )
}

safe_read_summary <- function(path) {
  tryCatch(
    {
      readr::read_tsv(path, show_col_types = FALSE) %>%
        mutate(source_file = path)
    },
    error = function(e) {
      message("Failed to read: ", path)
      message("  ", e$message)
      tibble()
    }
  )
}

theme_mortality <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 3, hjust = 0),
      plot.subtitle = element_text(size = base_size, color = "grey25", hjust = 0),
      plot.caption = element_text(size = base_size - 2, color = "grey35", hjust = 0),
      axis.title = element_text(face = "bold", color = "grey15"),
      axis.text = element_text(color = "grey15"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "grey85", linewidth = 0.25),
      strip.text = element_text(face = "bold", size = base_size, color = "grey10"),
      strip.background = element_rect(fill = "grey92", color = NA),
      legend.title = element_text(face = "bold"),
      legend.position = "top",
      plot.margin = margin(12, 55, 12, 12)
    )
}

# ============================================================
# 3. Fast find summary files
# ============================================================

clock_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)

candidate_survival_dirs <- file.path(clock_dirs, survival_dir_name)
candidate_survival_dirs <- candidate_survival_dirs[dir.exists(candidate_survival_dirs)]

message("Detected survival_analysis_mortality directories: ", length(candidate_survival_dirs))

if (length(candidate_survival_dirs) == 0) {
  stop(
    "No survival_analysis_mortality directories found under:\n",
    "  ", base_dir, "\n\n",
    "Expected structure:\n",
    "  <base_dir>/<clock_folder>/survival_analysis_mortality/"
  )
}

summary_files <- unlist(
  lapply(
    candidate_survival_dirs,
    function(d) {
      list.files(
        d,
        pattern = "_mortality_survival_summary\\.tsv$",
        full.names = TRUE,
        recursive = FALSE
      )
    }
  ),
  use.names = FALSE
)

summary_files <- unique(summary_files)
summary_files <- summary_files[file.exists(summary_files)]

summary_manifest <- tibble(
  summary_file = summary_files,
  clock_folder = basename(dirname(dirname(summary_files))),
  survival_dir = basename(dirname(summary_files))
) %>%
  arrange(clock_folder, summary_file)

manifest_out <- file.path(out_dir, "mortality_summary_file_manifest.tsv")
readr::write_tsv(summary_manifest, manifest_out)

message("Summary-file manifest saved:")
message("  ", manifest_out)

if (length(summary_files) == 0) {
  stop(
    "No mortality survival summary files found under detected survival directories.\n\n",
    "Expected pattern:\n",
    "  <base_dir>/<clock_folder>/survival_analysis_mortality/*_mortality_survival_summary.tsv\n\n",
    "Manifest written to:\n",
    "  ", manifest_out
  )
}

message("Found ", length(summary_files), " mortality survival summary files.")

if (length(summary_files) != 47) {
  warning(
    "Expected 47 stable/significant disease-clock mortality summary files, but found ",
    length(summary_files), ". The script will continue with all detected files."
  )
}

message("First detected summary files:")
print(head(summary_files, 10))

# ============================================================
# 4. Read and validate summary files
# ============================================================

res_raw <- purrr::map_dfr(summary_files, safe_read_summary)

required_cols <- c(
  "disease",
  "clock_label",
  "folder",
  "modality",
  "score_col",
  "n_analysis_rows",
  "n_deaths",
  "clock_hr_per_1sd",
  "clock_hr_ci_lower",
  "clock_hr_ci_upper",
  "clock_p",
  "cindex_covariates",
  "cindex_covariates_plus_clock",
  "delta_cindex_clock_vs_covariates"
)

missing_cols <- setdiff(required_cols, colnames(res_raw))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# ============================================================
# 5. Clean and annotate table
# ============================================================

res <- res_raw %>%
  mutate(
    disease = tolower(as.character(disease)),
    disease_label = dplyr::recode(
      disease,
      !!!disease_labels,
      .default = str_to_title(disease)
    ),
    disease = factor(disease, levels = disease_order),
    disease_label = factor(
      disease_label,
      levels = unname(disease_labels[disease_order])
    ),
    
    modality = as.character(modality),
    modality = if_else(modality %in% modality_order, modality, "Other"),
    modality = factor(modality, levels = c(modality_order, "Other")),
    
    clock_label = as.character(clock_label),
    folder = as.character(folder),
    score_col = as.character(score_col),
    
    n_analysis_rows = as.numeric(n_analysis_rows),
    n_deaths = as.numeric(n_deaths),
    event_rate = as.numeric(event_rate),
    median_followup_years = as.numeric(median_followup_years),
    
    clock_hr_per_1sd = as.numeric(clock_hr_per_1sd),
    clock_hr_ci_lower = as.numeric(clock_hr_ci_lower),
    clock_hr_ci_upper = as.numeric(clock_hr_ci_upper),
    clock_coef = as.numeric(clock_coef),
    clock_coef_se = as.numeric(clock_coef_se),
    clock_p = as.numeric(clock_p),
    
    cindex_covariates = as.numeric(cindex_covariates),
    cindex_covariates_plus_clock = as.numeric(cindex_covariates_plus_clock),
    delta_cindex_clock_vs_covariates = as.numeric(delta_cindex_clock_vs_covariates),
    
    clock_q_bh = p.adjust(clock_p, method = "BH"),
    
    bonferroni_significant = !is.na(clock_p) & clock_p < bonferroni_p_threshold,
    bonferroni_status = if_else(
      bonferroni_significant,
      "Bonferroni significant",
      "Not Bonferroni significant"
    ),
    bonferroni_status = factor(
      bonferroni_status,
      levels = c("Bonferroni significant", "Not Bonferroni significant")
    ),
    
    hr_label = format_hr_ci(clock_hr_per_1sd, clock_hr_ci_lower, clock_hr_ci_upper),
    cases_label = paste0("Cases=", comma(n_deaths, accuracy = 1)),
    clock_display = str_squish(clock_label),
    
    # Requested: annotate cases, not P-values.
    annotation_label = cases_label,
    delta_cindex_label = sprintf("%+.3f", delta_cindex_clock_vs_covariates)
  ) %>%
  filter(
    !is.na(clock_hr_per_1sd),
    !is.na(clock_hr_ci_lower),
    !is.na(clock_hr_ci_upper),
    clock_hr_per_1sd > 0,
    clock_hr_ci_lower > 0,
    clock_hr_ci_upper > 0
  ) %>%
  arrange(disease, desc(clock_hr_per_1sd))

merged_out <- file.path(out_dir, "all_disease_clock_mortality_survival_summary_merged.tsv")
readr::write_tsv(res, merged_out)

message("Merged summary saved:")
message("  ", merged_out)

# ============================================================
# 6. Main Figure 1: Disease-faceted HR forest plot
# ============================================================

forest_tbl <- res %>%
  group_by(disease_label) %>%
  arrange(desc(clock_hr_per_1sd), .by_group = TRUE) %>%
  mutate(
    row_id = row_number(),
    plot_id = paste0(
      as.character(disease_label),
      "___",
      sprintf("%02d", row_id),
      "___",
      clock_display
    )
  ) %>%
  ungroup() %>%
  arrange(disease_label, row_id) %>%
  mutate(
    plot_id = factor(plot_id, levels = rev(unique(plot_id))),
    label_x = clock_hr_ci_upper * 1.035
  )

x_min <- min(c(1, forest_tbl$clock_hr_ci_lower), na.rm = TRUE)
x_max <- max(c(1, forest_tbl$clock_hr_ci_upper), na.rm = TRUE)

x_min <- max(0.50, x_min * 0.90)
x_max_label <- x_max * 1.45

hr_breaks <- c(
  0.5, 0.6, 0.7, 0.8, 0.9,
  1.0, 1.1, 1.2, 1.3, 1.5,
  1.75, 2.0, 2.5, 3.0, 4.0, 5.0
)

hr_breaks <- hr_breaks[hr_breaks >= x_min & hr_breaks <= x_max_label]

forest_height <- max(8.5, 2.4 + 0.24 * nrow(forest_tbl))

p_forest <- ggplot(forest_tbl, aes(y = plot_id)) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey35"
  ) +
  geom_segment(
    aes(
      x = clock_hr_ci_lower,
      xend = clock_hr_ci_upper,
      yend = plot_id,
      color = modality
    ),
    linewidth = 1.05,
    alpha = 0.88,
    lineend = "round"
  ) +
  
  # Open circles: not Bonferroni significant.
  geom_point(
    data = forest_tbl %>% filter(!bonferroni_significant),
    aes(
      x = clock_hr_per_1sd,
      y = plot_id,
      color = modality
    ),
    shape = 21,
    fill = "white",
    size = 3.4,
    stroke = 0.85,
    inherit.aes = FALSE
  ) +
  
  # Solid filled circles: Bonferroni significant.
  geom_point(
    data = forest_tbl %>% filter(bonferroni_significant),
    aes(
      x = clock_hr_per_1sd,
      y = plot_id,
      fill = modality
    ),
    shape = 21,
    color = "grey10",
    size = 3.7,
    stroke = 0.35,
    inherit.aes = FALSE
  ) +
  
  # Annotate number of mortality cases/deaths, not P-values.
  geom_text(
    aes(
      x = label_x,
      label = annotation_label
    ),
    hjust = 0,
    size = 2.55,
    color = "grey10"
  ) +
  facet_grid(
    disease_label ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_y_discrete(
    labels = function(x) str_replace(x, "^.*___", "")
  ) +
  scale_x_log10(
    breaks = hr_breaks,
    labels = function(x) sprintf("%.2f", x),
    limits = c(x_min, x_max_label),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_color_manual(values = modality_cols, drop = FALSE) +
  scale_fill_manual(values = modality_cols, drop = FALSE) +
  labs(
    title = "Disease L'EPOCH clocks predict future all-cause mortality",
    subtitle = paste0(
      "Effect sizes are hazard ratios per 1-SD higher disease-clock acceleration score. ",
      "Models were adjusted for age, sex, ethnicity, assessment center, smoking, BMI, and blood pressure."
    ),
    x = "Hazard ratio for mortality per 1-SD clock acceleration",
    y = NULL,
    color = "Modality",
    fill = "Modality",
    caption = paste0(
      "Dashed line indicates HR = 1. Error bars show 95% CI. ",
      "Filled circles indicate Bonferroni-significant associations, P < 0.05/47/5 = ",
      signif(bonferroni_p_threshold, 3),
      ". Open circles are not Bonferroni significant. Right-side labels show mortality cases/deaths."
    )
  ) +
  coord_cartesian(clip = "off") +
  theme_mortality(base_size = 11) +
  theme(
    legend.box = "vertical",
    legend.margin = margin(0, 0, 0, 0),
    strip.text.y = element_text(angle = 0),
    axis.text.y = element_text(size = 8.6),
    panel.spacing.y = unit(0.55, "lines")
  )

forest_pdf <- file.path(out_dir, "all_disease_clock_mortality_hr_forest.pdf")
forest_png <- file.path(out_dir, "all_disease_clock_mortality_hr_forest.png")

ggsave(forest_pdf, p_forest, width = 13.2, height = forest_height, units = "in")
ggsave(forest_png, p_forest, width = 13.2, height = forest_height, units = "in", dpi = 320)

message("Forest plot saved:")
message("  ", forest_pdf)
message("  ", forest_png)

# ============================================================
# 7. Main Figure 2: HR heatmap
# ============================================================

heat_tbl <- res %>%
  mutate(
    disease_label = factor(disease_label, levels = unname(disease_labels[disease_order])),
    clock_row = clock_display
  ) %>%
  group_by(clock_row) %>%
  mutate(mean_log_hr = mean(log(clock_hr_per_1sd), na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(mean_log_hr), clock_row) %>%
  mutate(
    clock_row = factor(clock_row, levels = rev(unique(clock_row))),
    heat_label = paste0(
      sprintf("%.2f", clock_hr_per_1sd),
      "\n",
      "Cases=", comma(n_deaths, accuracy = 1)
    )
  )

heat_height <- max(8, 2.2 + 0.30 * n_distinct(heat_tbl$clock_row))

p_heat <- ggplot(
  heat_tbl,
  aes(
    x = disease_label,
    y = clock_row,
    fill = clock_hr_per_1sd
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.55,
    width = 0.95,
    height = 0.95
  ) +
  
  # Small solid dot inside heatmap cell for Bonferroni-significant signals.
  geom_point(
    data = heat_tbl %>% filter(bonferroni_significant),
    aes(x = disease_label, y = clock_row),
    inherit.aes = FALSE,
    shape = 21,
    fill = "grey10",
    color = "grey10",
    size = 1.7,
    position = position_nudge(x = 0.34, y = 0.28)
  ) +
  
  geom_text(
    aes(label = heat_label),
    size = 2.45,
    lineheight = 0.85,
    color = "grey10"
  ) +
  scale_fill_gradient2(
    low = "#3B6EA8",
    mid = "#F7F7F7",
    high = "#B33A3A",
    midpoint = 1,
    na.value = "grey92",
    name = "HR"
  ) +
  labs(
    title = "Mortality hazard ratios across disease L'EPOCH clocks",
    subtitle = "Each cell shows HR and number of mortality cases. Solid dots indicate Bonferroni-significant signals.",
    x = "Disease clock target",
    y = NULL,
    caption = paste0(
      "Solid dots indicate P < 0.05/47/5 = ",
      signif(bonferroni_p_threshold, 3),
      "."
    )
  ) +
  theme_mortality(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(size = 8.5),
    legend.position = "right"
  )

heat_pdf <- file.path(out_dir, "all_disease_clock_mortality_hr_heatmap.pdf")
heat_png <- file.path(out_dir, "all_disease_clock_mortality_hr_heatmap.png")

ggsave(heat_pdf, p_heat, width = 10.2, height = heat_height, units = "in")
ggsave(heat_png, p_heat, width = 10.2, height = heat_height, units = "in", dpi = 320)

message("Heatmap saved:")
message("  ", heat_pdf)
message("  ", heat_png)

# ============================================================
# 8. Main Figure 3: Ranked lollipop plot
# ============================================================

rank_tbl <- res %>%
  arrange(clock_hr_per_1sd) %>%
  mutate(
    rank_label = paste0(disease_label, ": ", clock_display),
    rank_label = factor(rank_label, levels = rank_label),
    label_text = paste0("Cases=", comma(n_deaths, accuracy = 1)),
    label_x = clock_hr_ci_upper * 1.025
  )

x_min_rank <- min(c(1, rank_tbl$clock_hr_ci_lower), na.rm = TRUE)
x_max_rank <- max(c(rank_tbl$clock_hr_ci_upper), na.rm = TRUE)

x_min_rank <- max(0.50, x_min_rank * 0.90)
x_max_rank_label <- x_max_rank * 1.35

rank_breaks <- c(
  0.5, 0.6, 0.7, 0.8, 0.9,
  1.0, 1.1, 1.2, 1.3, 1.5,
  1.75, 2.0, 2.5, 3.0, 4.0, 5.0
)

rank_breaks <- rank_breaks[rank_breaks >= x_min_rank & rank_breaks <= x_max_rank_label]

rank_height <- max(8.5, 2.2 + 0.23 * nrow(rank_tbl))

p_rank <- ggplot(rank_tbl, aes(y = rank_label)) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey35"
  ) +
  geom_segment(
    aes(
      x = 1,
      xend = clock_hr_per_1sd,
      yend = rank_label,
      color = disease_label
    ),
    linewidth = 1.05,
    alpha = 0.70,
    lineend = "round"
  ) +
  geom_segment(
    aes(
      x = clock_hr_ci_lower,
      xend = clock_hr_ci_upper,
      yend = rank_label
    ),
    linewidth = 0.55,
    color = "grey25",
    alpha = 0.70
  ) +
  
  # Open circles: not Bonferroni significant.
  geom_point(
    data = rank_tbl %>% filter(!bonferroni_significant),
    aes(
      x = clock_hr_per_1sd,
      y = rank_label,
      color = disease_label
    ),
    shape = 21,
    fill = "white",
    size = 3.1,
    stroke = 0.85,
    inherit.aes = FALSE
  ) +
  
  # Solid filled circles: Bonferroni significant.
  geom_point(
    data = rank_tbl %>% filter(bonferroni_significant),
    aes(
      x = clock_hr_per_1sd,
      y = rank_label,
      fill = disease_label
    ),
    shape = 21,
    color = "grey10",
    size = 3.4,
    stroke = 0.35,
    inherit.aes = FALSE
  ) +
  
  geom_text(
    aes(
      x = label_x,
      label = label_text
    ),
    hjust = 0,
    size = 2.4,
    color = "grey10"
  ) +
  scale_x_log10(
    breaks = rank_breaks,
    labels = function(x) sprintf("%.2f", x),
    limits = c(x_min_rank, x_max_rank_label),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_color_manual(values = disease_cols, drop = FALSE) +
  scale_fill_manual(values = disease_cols, drop = FALSE) +
  labs(
    title = "Ranked mortality associations of disease L'EPOCH clocks",
    subtitle = "Clocks are ranked by HR per 1-SD higher disease-clock acceleration score.",
    x = "Hazard ratio for mortality per 1-SD clock acceleration",
    y = NULL,
    color = "Disease target",
    fill = "Disease target",
    caption = paste0(
      "Filled circles indicate Bonferroni-significant associations, P < 0.05/47/5 = ",
      signif(bonferroni_p_threshold, 3),
      ". Open circles are not Bonferroni significant. Right-side labels show mortality cases/deaths."
    )
  ) +
  coord_cartesian(clip = "off") +
  theme_mortality(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 7.3),
    legend.box = "vertical"
  )

rank_pdf <- file.path(out_dir, "all_disease_clock_mortality_hr_ranked_lollipop.pdf")
rank_png <- file.path(out_dir, "all_disease_clock_mortality_hr_ranked_lollipop.png")

ggsave(rank_pdf, p_rank, width = 13.2, height = rank_height, units = "in")
ggsave(rank_png, p_rank, width = 13.2, height = rank_height, units = "in", dpi = 320)

message("Ranked lollipop plot saved:")
message("  ", rank_pdf)
message("  ", rank_png)

# ============================================================
# 9. Manuscript-ready result table
# ============================================================

manuscript_tbl <- res %>%
  transmute(
    disease = as.character(disease_label),
    clock = clock_display,
    modality = as.character(modality),
    n = n_analysis_rows,
    deaths = n_deaths,
    event_rate = event_rate,
    median_followup_years = median_followup_years,
    HR_per_1SD = clock_hr_per_1sd,
    HR_95CI_lower = clock_hr_ci_lower,
    HR_95CI_upper = clock_hr_ci_upper,
    P = clock_p,
    BH_FDR = clock_q_bh,
    Bonferroni_threshold = bonferroni_p_threshold,
    Bonferroni_significant = bonferroni_significant,
    C_index_covariates = cindex_covariates,
    C_index_covariates_plus_clock = cindex_covariates_plus_clock,
    Delta_C_index = delta_cindex_clock_vs_covariates,
    folder = folder,
    score_col = score_col,
    source_file = source_file
  ) %>%
  arrange(disease, desc(HR_per_1SD))

manuscript_out <- file.path(out_dir, "all_disease_clock_mortality_survival_manuscript_table.tsv")
readr::write_tsv(manuscript_tbl, manuscript_out)

message("Manuscript-ready table saved:")
message("  ", manuscript_out)

# ============================================================
# 10. Save session info
# ============================================================

sink(file.path(out_dir, "plotting_session_info.txt"))
sessionInfo()
sink()

message("Done.")