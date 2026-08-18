#!/usr/bin/env Rscript

# ==============================================================================
# Publication-ready comparison plots:
# Brain MRI mortality EPOCH vs 9 AI-derived disease subtype scores
#
# Endpoint:
#   First incident inpatient ICD-10 F* or G* diagnosis after brain MRI
#
# Primary analysis:
#   Held-out EPOCH test split, strict common complete-case sample
#
# Required inputs:
#   marker_models_common_sample.tsv
#   epoch_vs_subtype_pairwise.tsv
#   combined_subtypes_vs_epoch.tsv
#   bootstrap_cindex_common_sample.tsv
#
# Outputs:
#   A. HR forest plot
#   B. Incremental C-index vs covariate-only model
#   C. Paired bootstrap C-index difference: EPOCH - each subtype
#   D. Base vs EPOCH vs all-9-subtypes vs all-9+EPOCH
#   Combined 4-panel PDF/PNG
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. User settings
# ------------------------------------------------------------------------------

input_dir <- paste0(
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/",
  "brain_mri_mortality_clock/",
  "comparative_with_9_disease_subtypes_any_FG"
)

output_dir <- file.path(input_dir, "plots")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

prefix <- "brain_mri_mortality_EPOCH_vs_9_AI_subtypes_any_FG"

# ------------------------------------------------------------------------------
# 2. Packages
# ------------------------------------------------------------------------------

required_packages <- c(
  "readr", "dplyr", "tidyr", "ggplot2",
  "forcats", "scales", "stringr", "patchwork"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall with:\ninstall.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))"
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(scales)
  library(stringr)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# 3. Read inputs
# ------------------------------------------------------------------------------

marker_file <- file.path(input_dir, "marker_models_common_sample.tsv")
pairwise_file <- file.path(input_dir, "epoch_vs_subtype_pairwise.tsv")
combined_file <- file.path(input_dir, "combined_subtypes_vs_epoch.tsv")
bootstrap_file <- file.path(input_dir, "bootstrap_cindex_common_sample.tsv")

required_files <- c(marker_file, pairwise_file, combined_file, bootstrap_file)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required input files:\n",
    paste(missing_files, collapse = "\n")
  )
}

marker_tbl <- read_tsv(marker_file, show_col_types = FALSE, progress = FALSE)
pairwise_tbl <- read_tsv(pairwise_file, show_col_types = FALSE, progress = FALSE)
combined_tbl <- read_tsv(combined_file, show_col_types = FALSE, progress = FALSE)
bootstrap_tbl <- read_tsv(bootstrap_file, show_col_types = FALSE, progress = FALSE)

# ------------------------------------------------------------------------------
# 4. Validation
# ------------------------------------------------------------------------------

check_columns <- function(df, required_cols, table_name) {
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      table_name,
      " is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
}

check_columns(
  marker_tbl,
  c(
    "marker", "N", "N_case", "HR_per_1SD", "HR_CI_low", "HR_CI_high",
    "p_marker", "base_cindex", "marker_cindex", "delta_cindex_vs_base"
  ),
  "marker_models_common_sample.tsv"
)

check_columns(
  pairwise_tbl,
  c(
    "subtype", "EPOCH_standalone_cindex", "subtype_standalone_cindex",
    "delta_cindex_EPOCH_minus_subtype", "lrt_p_EPOCH_beyond_subtype",
    "lrt_p_subtype_beyond_EPOCH"
  ),
  "epoch_vs_subtype_pairwise.tsv"
)

check_columns(
  bootstrap_tbl,
  c(
    "comparison", "observed_difference", "bootstrap_CI_low",
    "bootstrap_CI_high", "bootstrap_p_two_sided"
  ),
  "bootstrap_cindex_common_sample.tsv"
)

check_columns(
  combined_tbl,
  c(
    "base_cindex", "EPOCH_only_cindex", "all_9_subtypes_cindex",
    "all_9_subtypes_plus_EPOCH_cindex", "delta_cindex_EPOCH_beyond_all9",
    "lrt_p_EPOCH_beyond_all9"
  ),
  "combined_subtypes_vs_epoch.tsv"
)

if (nrow(combined_tbl) != 1) {
  stop("combined_subtypes_vs_epoch.tsv should contain exactly one row.")
}

# ------------------------------------------------------------------------------
# 5. Labels/helpers
# ------------------------------------------------------------------------------

marker_order <- c(
  "Brain_MRI_mortality_EPOCH",
  "AD1", "AD2",
  "ASD1", "ASD2", "ASD3",
  "LLD1", "LLD2",
  "SCZ1", "SCZ2"
)

pretty_marker <- c(
  "Brain_MRI_mortality_EPOCH" = "Mortality EPOCH",
  "AD1" = "AD1",
  "AD2" = "AD2",
  "ASD1" = "ASD1",
  "ASD2" = "ASD2",
  "ASD3" = "ASD3",
  "LLD1" = "LLD1",
  "LLD2" = "LLD2",
  "SCZ1" = "SCZ1",
  "SCZ2" = "SCZ2"
)

format_p <- function(p) {
  case_when(
    is.na(p) ~ "P = NA",
    p < 0.001 ~ paste0("P = ", formatC(p, format = "e", digits = 2)),
    TRUE ~ paste0("P = ", formatC(p, format = "f", digits = 3))
  )
}

theme_epoch <- function(base_size = 11.5) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold", color = "black"),
      plot.title = element_text(face = "bold", size = base_size + 1.5),
      plot.subtitle = element_text(size = base_size - 0.6),
      plot.caption = element_text(size = base_size - 2.0, hjust = 0),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.margin = margin(8, 12, 8, 8)
    )
}

# EPOCH is emphasized; subtype scores share a common visual identity.
method_palette <- c(
  "Mortality EPOCH" = "#B55239",
  "AI disease subtype" = "#5E6D7A"
)

significance_shape <- c(
  "P < 0.05" = 16,
  "P >= 0.05" = 1
)

# ------------------------------------------------------------------------------
# 6. Panel A: hazard-ratio forest plot
# ------------------------------------------------------------------------------

panel_a_df <- marker_tbl |>
  mutate(
    marker = factor(marker, levels = marker_order),
    label = unname(pretty_marker[as.character(marker)]),
    method_group = if_else(
      as.character(marker) == "Brain_MRI_mortality_EPOCH",
      "Mortality EPOCH",
      "AI disease subtype"
    ),
    significance = if_else(p_marker < 0.05, "P < 0.05", "P >= 0.05")
  ) |>
  arrange(marker) |>
  mutate(
    label = factor(
      label,
      levels = rev(unname(pretty_marker[marker_order]))
    )
  )

p_hr <- ggplot(panel_a_df, aes(x = HR_per_1SD, y = label)) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey55"
  ) +
  geom_errorbarh(
    aes(
      xmin = HR_CI_low,
      xmax = HR_CI_high,
      color = method_group
    ),
    height = 0.16,
    linewidth = 0.65
  ) +
  geom_point(
    aes(
      color = method_group,
      shape = significance
    ),
    size = 2.9,
    stroke = 0.9
  ) +
  scale_color_manual(values = method_palette) +
  scale_shape_manual(values = significance_shape) +
  scale_x_continuous(
    name = "Hazard ratio per 1-SD higher score",
    breaks = pretty_breaks(n = 5)
  ) +
  scale_y_discrete(name = NULL) +
  labs(
    title = "A  Association with incident F/G disease",
    subtitle = "Cox models adjusted for age at imaging, sex, smoking, and BMI",
    caption = "Filled points indicate nominal P < 0.05."
  ) +
  theme_epoch()

# ------------------------------------------------------------------------------
# 7. Panel B: incremental C-index versus covariate-only model
# ------------------------------------------------------------------------------

panel_b_df <- marker_tbl |>
  mutate(
    marker = factor(marker, levels = marker_order),
    label = unname(pretty_marker[as.character(marker)]),
    method_group = if_else(
      as.character(marker) == "Brain_MRI_mortality_EPOCH",
      "Mortality EPOCH",
      "AI disease subtype"
    ),
    significance = if_else(p_marker < 0.05, "P < 0.05", "P >= 0.05")
  ) |>
  arrange(marker) |>
  mutate(
    label = factor(
      label,
      levels = rev(unname(pretty_marker[marker_order]))
    )
  )

delta_range <- range(panel_b_df$delta_cindex_vs_base, na.rm = TRUE)
delta_pad <- max(0.001, diff(delta_range) * 0.18)

p_delta <- ggplot(panel_b_df, aes(x = delta_cindex_vs_base, y = label)) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey55"
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = delta_cindex_vs_base,
      yend = label,
      color = method_group
    ),
    linewidth = 1.0,
    alpha = 0.80
  ) +
  geom_point(
    aes(
      color = method_group,
      shape = significance
    ),
    size = 3.0,
    stroke = 0.9
  ) +
  scale_color_manual(values = method_palette) +
  scale_shape_manual(values = significance_shape) +
  scale_x_continuous(
    name = expression(Delta*" C-index vs covariate-only model"),
    labels = label_number(accuracy = 0.001),
    limits = c(
      min(-0.0015, delta_range[1] - delta_pad),
      delta_range[2] + delta_pad
    )
  ) +
  scale_y_discrete(name = NULL) +
  labs(
    title = "B  Incremental discrimination",
    subtitle = paste0(
      "Covariate-only C-index = ",
      sprintf("%.3f", marker_tbl$base_cindex[[1]])
    ),
    caption = "Positive values indicate higher C-index than the covariate-only model."
  ) +
  theme_epoch()

# ------------------------------------------------------------------------------
# 8. Panel C: direct paired-bootstrap EPOCH-vs-subtype C-index comparison
# ------------------------------------------------------------------------------

panel_c_df <- bootstrap_tbl |>
  filter(str_detect(comparison, "^EPOCH - ")) |>
  mutate(
    subtype = str_remove(comparison, "^EPOCH - "),
    subtype = factor(
      subtype,
      levels = c(
        "AD1", "AD2",
        "ASD1", "ASD2", "ASD3",
        "LLD1", "LLD2",
        "SCZ1", "SCZ2"
      )
    ),
    subtype = fct_rev(subtype)
  )

if (nrow(panel_c_df) != 9) {
  warning(
    "Expected 9 EPOCH-minus-subtype bootstrap rows, found ",
    nrow(panel_c_df),
    "."
  )
}

ci_range <- range(
  c(panel_c_df$bootstrap_CI_low, panel_c_df$bootstrap_CI_high),
  na.rm = TRUE
)
ci_pad <- max(0.001, diff(ci_range) * 0.10)

p_boot <- ggplot(
  panel_c_df,
  aes(x = observed_difference, y = subtype)
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey55"
  ) +
  geom_errorbarh(
    aes(
      xmin = bootstrap_CI_low,
      xmax = bootstrap_CI_high
    ),
    height = 0.16,
    linewidth = 0.7,
    color = method_palette[["Mortality EPOCH"]]
  ) +
  geom_point(
    size = 3.1,
    color = method_palette[["Mortality EPOCH"]]
  ) +
  scale_x_continuous(
    name = expression(Delta*" C-index (EPOCH - disease subtype)"),
    labels = label_number(accuracy = 0.001),
    limits = c(
      min(ci_range[1] - ci_pad, -0.001),
      ci_range[2] + ci_pad
    )
  ) +
  scale_y_discrete(name = NULL) +
  labs(
    title = "C  Direct discrimination comparison",
    subtitle = "Paired bootstrap 95% confidence intervals (200 replicates)",
    caption = paste0(
      "Positive values favor EPOCH. A 95% CI crossing zero does not establish ",
      "significant C-index superiority."
    )
  ) +
  theme_epoch()

# ------------------------------------------------------------------------------
# 9. Panel D: EPOCH versus combined nine-subtype panel
# ------------------------------------------------------------------------------

c0 <- combined_tbl$base_cindex[[1]]
ce <- combined_tbl$EPOCH_only_cindex[[1]]
c9 <- combined_tbl$all_9_subtypes_cindex[[1]]
c10 <- combined_tbl$all_9_subtypes_plus_EPOCH_cindex[[1]]

p_epoch_beyond_all9 <- combined_tbl$lrt_p_EPOCH_beyond_all9[[1]]
delta_epoch_beyond_all9 <- combined_tbl$delta_cindex_EPOCH_beyond_all9[[1]]

panel_d_df <- tibble(
  model = factor(
    c(
      "Covariates only",
      "Covariates + EPOCH",
      "Covariates + 9 subtypes",
      "Covariates + 9 subtypes + EPOCH"
    ),
    levels = c(
      "Covariates only",
      "Covariates + EPOCH",
      "Covariates + 9 subtypes",
      "Covariates + 9 subtypes + EPOCH"
    )
  ),
  cindex = c(c0, ce, c9, c10),
  model_group = c(
    "Reference",
    "Mortality EPOCH",
    "AI disease subtype",
    "Combined"
  )
)

combined_palette <- c(
  "Reference" = "#8A8A8A",
  "Mortality EPOCH" = "#B55239",
  "AI disease subtype" = "#5E6D7A",
  "Combined" = "#6A5A89"
)

d_min <- min(panel_d_df$cindex, na.rm = TRUE)
d_max <- max(panel_d_df$cindex, na.rm = TRUE)

p_combined <- ggplot(
  panel_d_df,
  aes(x = model, y = cindex, fill = model_group)
) +
  geom_col(
    width = 0.68,
    color = "grey25",
    linewidth = 0.25
  ) +
  geom_text(
    aes(label = sprintf("%.3f", cindex)),
    vjust = -0.45,
    size = 3.5,
    fontface = "bold"
  ) +
  scale_fill_manual(values = combined_palette) +
  scale_y_continuous(
    name = "C-index",
    labels = label_number(accuracy = 0.001),
    limits = c(
      max(0, d_min - 0.015),
      d_max + 0.012
    ),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_x_discrete(
    name = NULL,
    labels = c(
      "Covariates only" = "Covariates\nonly",
      "Covariates + EPOCH" = "Covariates\n+ EPOCH",
      "Covariates + 9 subtypes" = "Covariates\n+ 9 subtypes",
      "Covariates + 9 subtypes + EPOCH" = "Covariates\n+ 9 subtypes\n+ EPOCH"
    )
  ) +
  labs(
    title = "D  EPOCH versus the combined subtype panel",
    subtitle = paste0(
      "Adding EPOCH beyond all 9 subtypes: ",
      "Delta C = ", sprintf("%.3f", delta_epoch_beyond_all9),
      "; ", format_p(p_epoch_beyond_all9),
      " (likelihood-ratio test)"
    ),
    caption = paste0(
      "The likelihood-ratio test evaluates incremental Cox-model fit, ",
      "not a paired test of C-index superiority."
    )
  ) +
  theme_epoch() +
  theme(
    axis.text.x = element_text(size = 9.2, lineheight = 0.95)
  )

# ------------------------------------------------------------------------------
# 10. Save individual panels
# ------------------------------------------------------------------------------

save_panel <- function(plot_object, stem, width = 8.2, height = 5.7) {
  ggsave(
    filename = file.path(output_dir, paste0(prefix, "_", stem, ".pdf")),
    plot = plot_object,
    width = width,
    height = height,
    device = cairo_pdf
  )

  ggsave(
    filename = file.path(output_dir, paste0(prefix, "_", stem, ".png")),
    plot = plot_object,
    width = width,
    height = height,
    dpi = 400,
    bg = "white"
  )
}

save_panel(p_hr, "A_HR_forest")
save_panel(p_delta, "B_delta_Cindex_vs_base")
save_panel(p_boot, "C_bootstrap_EPOCH_minus_subtypes")
save_panel(p_combined, "D_combined_9subtypes_vs_EPOCH")

# ------------------------------------------------------------------------------
# 11. Combined 4-panel figure
# ------------------------------------------------------------------------------

combined_figure <- (
  p_hr | p_delta
) / (
  p_boot | p_combined
) +
  plot_annotation(
    title = paste0(
      "Brain MRI mortality EPOCH versus AI-derived disease subtypes ",
      "for incident mental or nervous-system disorders"
    ),
    subtitle = paste0(
      "Held-out EPOCH test split; all predictors evaluated in one strict ",
      "common complete-case sample"
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 11)
    )
  )

ggsave(
  filename = file.path(output_dir, paste0(prefix, "_combined_4panel.pdf")),
  plot = combined_figure,
  width = 13.2,
  height = 10.2,
  device = cairo_pdf
)

ggsave(
  filename = file.path(output_dir, paste0(prefix, "_combined_4panel.png")),
  plot = combined_figure,
  width = 13.2,
  height = 10.2,
  dpi = 400,
  bg = "white"
)

# ------------------------------------------------------------------------------
# 12. Save compact plotting summary
# ------------------------------------------------------------------------------

plot_summary <- marker_tbl |>
  select(
    marker,
    N,
    N_case,
    HR_per_1SD,
    HR_CI_low,
    HR_CI_high,
    p_marker,
    base_cindex,
    marker_cindex,
    delta_cindex_vs_base
  ) |>
  left_join(
    bootstrap_tbl |>
      filter(str_detect(comparison, "^EPOCH - ")) |>
      transmute(
        marker = str_remove(comparison, "^EPOCH - "),
        EPOCH_minus_subtype_delta_cindex = observed_difference,
        EPOCH_minus_subtype_bootstrap_CI_low = bootstrap_CI_low,
        EPOCH_minus_subtype_bootstrap_CI_high = bootstrap_CI_high,
        EPOCH_minus_subtype_bootstrap_p = bootstrap_p_two_sided
      ),
    by = "marker"
  )

write_tsv(
  plot_summary,
  file.path(output_dir, paste0(prefix, "_plotting_summary.tsv"))
)

# ------------------------------------------------------------------------------
# 13. Print key values
# ------------------------------------------------------------------------------

message("============================================================")
message("Plotting complete.")
message("Input directory: ", input_dir)
message("Output directory: ", output_dir)
message("")
message("Common-sample N = ", marker_tbl$N[[1]])
message("Incident F/G events = ", marker_tbl$N_case[[1]])
message("Covariate-only C-index = ", sprintf("%.4f", marker_tbl$base_cindex[[1]]))
message(
  "EPOCH C-index = ",
  sprintf(
    "%.4f",
    marker_tbl |>
      filter(marker == "Brain_MRI_mortality_EPOCH") |>
      pull(marker_cindex)
  )
)
message("All 9 subtypes C-index = ", sprintf("%.4f", c9))
message("All 9 subtypes + EPOCH C-index = ", sprintf("%.4f", c10))
message("LRT for EPOCH beyond all 9 subtypes: ", format_p(p_epoch_beyond_all9))
message("============================================================")