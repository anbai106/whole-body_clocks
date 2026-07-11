#!/usr/bin/env Rscript

# ============================================================
# Combined mortality-clock acceleration-year distribution plot
# for 23 mortality clocks
# MRI + proteomics + metabolomics
#
# Input:
#   <base_dir>/<clock_folder>/<prefix>_predictions.tsv
#
# Plotted value:
#   <organ>_<modality>_mortality_clock_acceleration_years
#
# Output:
#   all_mortality_clock_acceleration_years_distribution_data.tsv
#   all_mortality_clock_acceleration_years_distribution_summary.tsv
#   all_mortality_clock_acceleration_years_distribution_combined.pdf/png
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
  library(glue)
})

# ============================================================
# 1. Settings
# ============================================================

possible_base_dirs <- c(
  "/gpfs/fs001/cbica/home/wenju/Reproducibile_paper/WholeBodyClock",
  "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock",
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock",
  getwd()
)

base_dir <- Sys.getenv("WHOLEBODYCLOCK_BASE_DIR", unset = NA_character_)

if (is.na(base_dir) || base_dir == "") {
  base_dir <- possible_base_dirs[file.exists(possible_base_dirs)][1]
}

if (is.na(base_dir) || is.null(base_dir) || !dir.exists(base_dir)) {
  stop(
    "Could not find base_dir. Please set it manually, e.g.\n",
    "base_dir <- '/gpfs/fs001/cbica/home/wenju/Reproducibile_paper/WholeBodyClock'"
  )
}

base_dir <- normalizePath(base_dir, mustWork = TRUE)

combined_outdir <- file.path(
  base_dir,
  "all_mortality_clock_acceleration_year_distributions"
)

dir.create(combined_outdir, recursive = TRUE, showWarnings = FALSE)

skip_missing <- TRUE

# Figure options
ncol_panels <- 4
use_free_x <- TRUE
use_free_y <- TRUE

split_cols <- c(
  "train" = "#2E86AB",
  "validation" = "#F18F01",
  "test" = "#6A994E"
)

modality_order <- c("MRI", "Proteomics", "Metabolomics")

message("Base directory: ", base_dir)
message("Output directory: ", combined_outdir)

# ============================================================
# 2. Clock manifest: 23 mortality clocks
# ============================================================

make_manifest <- function(base_dir) {
  mri_organs <- c(
    "brain", "heart", "adipose", "kidney", "liver", "pancreas", "spleen"
  )

  proteomics_organs <- c(
    "Reproductive_female", "Pulmonary", "Heart", "Brain", "Eye", "Hepatic",
    "Renal", "Reproductive_male", "Endocrine", "Immune", "Skin"
  )

  metabolomics_organs <- c(
    "Endocrine", "Digestive", "Hepatic", "Immune", "Metabolic"
  )

  bind_rows(
    tibble(
      modality = "MRI",
      modality_key = "mri",
      organ_folder_name = mri_organs,
      organ_key = stringr::str_to_lower(mri_organs),
      folder = paste0(stringr::str_to_lower(mri_organs), "_mri_mortality_clock"),
      prefix = paste0(stringr::str_to_lower(mri_organs), "_mri_mortality_clock")
    ),
    tibble(
      modality = "Proteomics",
      modality_key = "proteomics",
      organ_folder_name = proteomics_organs,
      organ_key = stringr::str_to_lower(proteomics_organs),
      folder = paste0(proteomics_organs, "_proteomics_mortality_clock"),
      prefix = paste0(stringr::str_to_lower(proteomics_organs), "_proteomics_mortality_clock")
    ),
    tibble(
      modality = "Metabolomics",
      modality_key = "metabolomics",
      organ_folder_name = metabolomics_organs,
      organ_key = stringr::str_to_lower(metabolomics_organs),
      folder = paste0(metabolomics_organs, "_metabolomics_mortality_clock"),
      prefix = paste0(stringr::str_to_lower(metabolomics_organs), "_metabolomics_mortality_clock")
    )
  ) %>%
    mutate(
      organ_label = organ_folder_name %>%
        stringr::str_replace_all("_", " ") %>%
        stringr::str_to_sentence(),
      clock_label = paste(organ_label, modality),
      clock_id = paste(organ_key, modality_key, sep = "__"),
      clock_dir = file.path(base_dir, folder),
      prediction_file = file.path(clock_dir, paste0(prefix, "_predictions.tsv"))
    )
}

clock_manifest <- make_manifest(base_dir)

message("Expected clocks: ", nrow(clock_manifest))

# ============================================================
# 3. Helpers
# ============================================================

theme_combined_density <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(color = "#17202A"),
      plot.title = element_text(face = "bold", size = base_size + 5),
      plot.subtitle = element_text(size = base_size + 1, color = "#566573"),
      plot.caption = element_text(size = base_size - 2, color = "#7B7D7D", hjust = 0),

      axis.title.x = element_text(face = "bold", size = 9),
      axis.title.y = element_text(face = "bold", size = 9),

      axis.text.x = element_text(size = 13, color = "#2C3E50"),
      axis.text.y = element_text(size = 13, color = "#2C3E50"),

      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "#EAECEE", linewidth = 0.30),

      strip.text = element_text(face = "bold", size = 8.5),
      strip.background = element_rect(fill = "grey95", color = NA),

      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 8),

      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.spacing = unit(0.85, "lines")
    )
}

normalize_split <- function(x) {
  x_chr <- as.character(x)

  case_when(
    x_chr %in% c("train", "training", "Train", "TRAIN") ~ "train",
    x_chr %in% c("validation", "valid", "val", "Validation", "VALIDATION") ~ "validation",
    x_chr %in% c("test", "Test", "TEST") ~ "test",
    TRUE ~ x_chr
  )
}

read_prediction_one_clock <- function(meta_row) {
  meta <- as.list(meta_row)

  if (!file.exists(meta$prediction_file)) {
    msg <- paste0(
      "Missing prediction file for ", meta$clock_label, ":\n  ",
      meta$prediction_file
    )

    if (skip_missing) {
      warning(msg)
      return(NULL)
    } else {
      stop(msg)
    }
  }

  pred <- readr::read_tsv(
    meta$prediction_file,
    show_col_types = FALSE,
    progress = FALSE
  )

  feature_token <- paste0(meta$organ_key, "_", meta$modality_key)

  expected_accel_year_col <- paste0(
    feature_token,
    "_mortality_clock_acceleration_years"
  )

  accel_year_col <- expected_accel_year_col

  # Robust fallback if exact column is not present.
  if (!accel_year_col %in% colnames(pred)) {
    fallback_cols <- grep(
      "_mortality_clock_acceleration_years$",
      colnames(pred),
      value = TRUE
    )

    if (length(fallback_cols) == 1) {
      accel_year_col <- fallback_cols[[1]]
    } else {
      warning(
        "Could not identify mortality-clock acceleration-years column for ",
        meta$clock_label,
        ". Expected: ",
        expected_accel_year_col,
        ". Available acceleration-years columns: ",
        paste(fallback_cols, collapse = ", ")
      )
      return(NULL)
    }
  }

  required_cols <- c("split", accel_year_col)
  missing_cols <- setdiff(required_cols, colnames(pred))

  if (length(missing_cols) > 0) {
    warning(
      "Skipping ", meta$clock_label,
      " because missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
    return(NULL)
  }

  participant_col <- dplyr::case_when(
    "participant_id" %in% colnames(pred) ~ "participant_id",
    "eid" %in% colnames(pred) ~ "eid",
    "IID" %in% colnames(pred) ~ "IID",
    "id" %in% colnames(pred) ~ "id",
    TRUE ~ NA_character_
  )

  out <- pred %>%
    mutate(
      participant_id = if (!is.na(participant_col)) {
        as.character(.data[[participant_col]])
      } else {
        as.character(row_number())
      },
      split = normalize_split(.data[["split"]]),
      split = factor(split, levels = c("train", "validation", "test")),
      clock_acceleration_years = as.numeric(.data[[accel_year_col]])
    ) %>%
    filter(
      !is.na(split),
      is.finite(clock_acceleration_years)
    ) %>%
    transmute(
      participant_id,
      split,
      clock_acceleration_years,
      clock_id = meta$clock_id,
      modality = meta$modality,
      modality_key = meta$modality_key,
      organ_key = meta$organ_key,
      organ_label = meta$organ_label,
      clock_label = meta$clock_label,
      folder = meta$folder,
      prefix = meta$prefix,
      acceleration_years_col = accel_year_col,
      prediction_file = meta$prediction_file
    )

  if (nrow(out) == 0) {
    warning(
      "No usable acceleration-years rows after filtering for ",
      meta$clock_label
    )
    return(NULL)
  }

  message(
    "Loaded ", meta$clock_label,
    ": N = ", nrow(out),
    "; acceleration-years column = ", accel_year_col
  )

  out
}

# ============================================================
# 4. Read all prediction files
# ============================================================

accel_df <- purrr::map_dfr(
  seq_len(nrow(clock_manifest)),
  ~ read_prediction_one_clock(clock_manifest[.x, ])
)

if (nrow(accel_df) == 0) {
  stop("No acceleration-years data were loaded. Please check prediction-file paths.")
}

accel_df <- accel_df %>%
  mutate(
    modality = factor(modality, levels = modality_order),
    split = factor(split, levels = c("train", "validation", "test"))
  )

# ============================================================
# 5. Summary table
# ============================================================

summary_tbl <- accel_df %>%
  group_by(modality, organ_label, clock_label, clock_id, split) %>%
  summarise(
    n = n(),
    mean = mean(clock_acceleration_years, na.rm = TRUE),
    sd = sd(clock_acceleration_years, na.rm = TRUE),
    median = median(clock_acceleration_years, na.rm = TRUE),
    q1 = quantile(clock_acceleration_years, 0.25, na.rm = TRUE),
    q3 = quantile(clock_acceleration_years, 0.75, na.rm = TRUE),
    min = min(clock_acceleration_years, na.rm = TRUE),
    max = max(clock_acceleration_years, na.rm = TRUE),
    .groups = "drop"
  )

clock_order_tbl <- accel_df %>%
  distinct(modality, organ_label, clock_label, clock_id) %>%
  mutate(
    modality = factor(modality, levels = modality_order)
  ) %>%
  arrange(modality, organ_label) %>%
  mutate(clock_label_plot = clock_label)

accel_df <- accel_df %>%
  left_join(
    clock_order_tbl %>%
      select(clock_id, clock_label_plot),
    by = "clock_id"
  ) %>%
  mutate(
    clock_label_plot = factor(
      clock_label_plot,
      levels = clock_order_tbl$clock_label_plot
    )
  )

summary_tbl <- summary_tbl %>%
  left_join(
    clock_order_tbl %>%
      select(clock_id, clock_label_plot),
    by = "clock_id"
  )

# Save tables.
accel_data_out <- file.path(
  combined_outdir,
  "all_mortality_clock_acceleration_years_distribution_data.tsv"
)

summary_out <- file.path(
  combined_outdir,
  "all_mortality_clock_acceleration_years_distribution_summary.tsv"
)

manifest_out <- file.path(
  combined_outdir,
  "all_mortality_clock_acceleration_years_distribution_manifest.tsv"
)

readr::write_tsv(accel_df, accel_data_out)
readr::write_tsv(summary_tbl, summary_out)
readr::write_tsv(clock_manifest, manifest_out)

message("Saved combined acceleration-years data:")
message("  ", accel_data_out)
message("  ", summary_out)
message("  ", manifest_out)

# ============================================================
# 6. Combined density-panel figure
# ============================================================

facet_scales <- case_when(
  use_free_x & use_free_y ~ "free",
  use_free_x & !use_free_y ~ "free_x",
  !use_free_x & use_free_y ~ "free_y",
  TRUE ~ "fixed"
)

p_combined_density <- ggplot(
  accel_df,
  aes(
    x = clock_acceleration_years,
    color = split,
    fill = split
  )
) +
  geom_density(
    alpha = 0.18,
    linewidth = 0.70,
    adjust = 1.05,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~ clock_label_plot,
    ncol = ncol_panels,
    scales = facet_scales
  ) +
  scale_color_manual(values = split_cols, drop = FALSE) +
  scale_fill_manual(values = split_cols, drop = FALSE) +
  labs(
    title = "Mortality EPOCH acceleration-year distributions across train, validation, and test sets",
    subtitle = glue(
      "Combined density plots for {n_distinct(accel_df$clock_id)} mortality clocks. ",
      "Broadly overlapping distributions indicate comparable acceleration-year ranges across data splits."
    ),
    x = "Mortality EPOCH acceleration years",
    y = "Density",
    color = NULL,
    fill = NULL,
    caption = "Each panel shows one organ/system-specific mortality EPOCH clock. Acceleration years are read directly from each clock's prediction TSV file."
  ) +
  theme_combined_density(base_size = 10)

combined_height <- max(
  12,
  ceiling(n_distinct(accel_df$clock_id) / ncol_panels) * 2.25
)

combined_pdf <- file.path(
  combined_outdir,
  "all_mortality_clock_acceleration_years_distribution_combined.pdf"
)

combined_png <- file.path(
  combined_outdir,
  "all_mortality_clock_acceleration_years_distribution_combined.png"
)

if (capabilities("cairo")) {
  ggsave(
    filename = combined_pdf,
    plot = p_combined_density,
    width = 15,
    height = combined_height,
    units = "in",
    device = cairo_pdf
  )
} else {
  ggsave(
    filename = combined_pdf,
    plot = p_combined_density,
    width = 15,
    height = combined_height,
    units = "in"
  )
}

ggsave(
  filename = combined_png,
  plot = p_combined_density,
  width = 15,
  height = combined_height,
  units = "in",
  dpi = 350
)

message("Saved combined acceleration-year density figure:")
message("  ", combined_pdf)
message("  ", combined_png)

# ============================================================
# 7. Optional modality-separated figures
# ============================================================

for (mod_i in modality_order) {
  df_i <- accel_df %>% filter(modality == mod_i)

  if (nrow(df_i) == 0) next

  p_i <- ggplot(
    df_i,
    aes(
      x = clock_acceleration_years,
      color = split,
      fill = split
    )
  ) +
    geom_density(
      alpha = 0.18,
      linewidth = 0.75,
      adjust = 1.05,
      na.rm = TRUE
    ) +
    facet_wrap(
      ~ clock_label_plot,
      ncol = 3,
      scales = facet_scales
    ) +
    scale_color_manual(values = split_cols, drop = FALSE) +
    scale_fill_manual(values = split_cols, drop = FALSE) +
    labs(
      title = paste0(mod_i, " mortality EPOCH acceleration-year distributions"),
      subtitle = "Train, validation, and test split density curves",
      x = "Mortality EPOCH acceleration years",
      y = "Density",
      color = NULL,
      fill = NULL
    ) +
    theme_combined_density(base_size = 10)

  mod_file <- stringr::str_to_lower(mod_i)

  pdf_i <- file.path(
    combined_outdir,
    paste0("mortality_clock_acceleration_years_distribution_", mod_file, ".pdf")
  )

  png_i <- file.path(
    combined_outdir,
    paste0("mortality_clock_acceleration_years_distribution_", mod_file, ".png")
  )

  height_i <- max(4.5, ceiling(n_distinct(df_i$clock_id) / 3) * 2.4)

  if (capabilities("cairo")) {
    ggsave(
      filename = pdf_i,
      plot = p_i,
      width = 11,
      height = height_i,
      units = "in",
      device = cairo_pdf
    )
  } else {
    ggsave(
      filename = pdf_i,
      plot = p_i,
      width = 11,
      height = height_i,
      units = "in"
    )
  }

  ggsave(
    filename = png_i,
    plot = p_i,
    width = 11,
    height = height_i,
    units = "in",
    dpi = 350
  )

  message("Saved ", mod_i, " acceleration-year density figure:")
  message("  ", pdf_i)
  message("  ", png_i)
}

# ============================================================
# 8. Print summary
# ============================================================

message("\n===== Loaded clocks =====")
print(
  accel_df %>%
    distinct(
      modality,
      organ_label,
      clock_label,
      clock_id,
      acceleration_years_col
    ) %>%
    arrange(modality, organ_label)
)

message("\n===== Split sample-size summary =====")
print(
  summary_tbl %>%
    select(
      modality,
      organ_label,
      clock_label,
      split,
      n,
      mean,
      sd,
      median
    ) %>%
    arrange(modality, organ_label, split)
)

message("\nDone.")