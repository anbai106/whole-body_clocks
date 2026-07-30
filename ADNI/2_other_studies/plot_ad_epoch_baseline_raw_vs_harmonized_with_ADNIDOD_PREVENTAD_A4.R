#!/usr/bin/env Rscript

# ==============================================================================
# Baseline-only AD EPOCH acceleration-years distributions:
#   - ADNI training baseline as the top reference row
#   - External studies shown by STUDY and SITE
#   - Raw features shown as unfilled violins
#   - Harmonized features shown as filled violins
#   - Van Gogh-inspired colors assigned by STUDY
#   - Aligned baseline-age range panel
#   - Quantification of between-SITE distribution shift before/after harmonization
#
# RStudio use:
#   1. Open this script.
#   2. Review the paths in "User settings".
#   3. Click Source.
# ==============================================================================

# ------------------------------------------------------------------------------
# User settings
# ------------------------------------------------------------------------------

raw_file <- paste0(
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/adni_lepoch/",
  "results_external_longitudinal_ad_epoch_raw/",
  "external_5_studies_adni_brain_mri_ad_epoch_raw_scan_level_predictions.tsv"
)

harmonized_file <- paste0(
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/adni_lepoch/",
  "results_external_longitudinal_ad_epoch_harmonized/",
  "external_5_studies_adni_brain_mri_ad_epoch_harmonized_scan_level_predictions.tsv"
)

adni_reference_file <- paste0(
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/adni_lepoch/",
  "results_brain_mri_ad_lepoch/",
  "adni_brain_mri_ad_lepoch_predictions.tsv"
)

outdir <- paste0(
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/adni_lepoch/",
  "results_external_longitudinal_ad_epoch_comparison"
)

prefix <- "external_5_studies_baseline_ad_epoch_with_adni_reference"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

message("Raw input: ", raw_file)
message("Harmonized input: ", harmonized_file)
message("ADNI reference input: ", adni_reference_file)
message("Output directory: ", outdir)

# ------------------------------------------------------------------------------
# Packages
# ------------------------------------------------------------------------------

required_packages <- c(
  "readr", "dplyr", "tidyr", "ggplot2",
  "forcats", "patchwork", "scales", "stringr"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
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
  library(patchwork)
  library(scales)
  library(stringr)
})

# ------------------------------------------------------------------------------
# Visual settings
# ------------------------------------------------------------------------------

# Van Gogh-inspired colors assigned by STUDY.
van_gogh_study_colors <- c(
  "ADNI"      = "#D9A400",  # sunflower ochre
  "ADNI_DOD"  = "#355C9A",  # cobalt blue
  "AIBL"      = "#2A9D8F",  # blue-green
  "BLSA"      = "#7A5195",  # violet
  "OASIS"     = "#D97706",  # burnt orange
  "PreventAD" = "#6B8E23"   # olive green
)

reference_fill <- "#F2C14E"
reference_outline <- "#8A6500"

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

detect_column <- function(
  df,
  preferred,
  regex,
  label,
  required = TRUE
) {
  direct <- preferred[preferred %in% names(df)]

  if (length(direct) >= 1) {
    return(direct[[1]])
  }

  candidates <- grep(
    regex,
    names(df),
    value = TRUE,
    ignore.case = TRUE
  )

  if (length(candidates) == 1) {
    return(candidates[[1]])
  }

  if (!required) {
    return(NA_character_)
  }

  if (length(candidates) > 1) {
    stop(
      "Multiple candidate columns found for ", label, ": ",
      paste(candidates, collapse = ", ")
    )
  }

  stop(
    "Could not identify ", label, ". Available columns include:\n",
    paste(head(names(df), 100), collapse = ", ")
  )
}

clean_character <- function(x, missing_label = "Unknown") {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "nan", "None", "null")] <- NA_character_
  x <- trimws(x)
  tidyr::replace_na(x, missing_label)
}

visit_to_month <- function(x) {
  x <- tolower(trimws(as.character(x)))
  out <- rep(NA_real_, length(x))

  baseline_mask <- x %in% c(
    "bl", "base", "baseline", "m00", "m0",
    "screen", "screening", "sc"
  )
  out[baseline_mask] <- 0

  remaining <- which(is.na(out) & !is.na(x))
  if (length(remaining) > 0) {
    matches <- stringr::str_match(
      x[remaining],
      "m(?:onth)?\\s*0*([0-9]+)"
    )
    valid <- !is.na(matches[, 2])
    out[remaining[valid]] <- as.numeric(matches[valid, 2])
  }

  remaining <- which(is.na(out) & !is.na(x))
  if (length(remaining) > 0) {
    matches <- stringr::str_match(
      x[remaining],
      "([0-9]+)"
    )
    valid <- !is.na(matches[, 2])
    out[remaining[valid]] <- as.numeric(matches[valid, 2])
  }

  out
}

weighted_sd <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  x <- x[keep]
  w <- w[keep]

  if (length(x) < 2 || sum(w) <= 0) {
    return(NA_real_)
  }

  mu <- weighted.mean(x, w)
  sqrt(sum(w * (x - mu)^2) / sum(w))
}

read_external_file <- function(path, feature_type) {
  if (!file.exists(path)) {
    stop("Input file does not exist: ", path)
  }

  df <- readr::read_tsv(
    path,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "unique"
  )

  accel_col <- detect_column(
    df,
    preferred = c(
      "adni_brain_mri_ad_epoch_acceleration_years",
      "adni_brain_mri_ad_lepoch_acceleration_years"
    ),
    regex = "acceleration[_\\.]*years$",
    label = "acceleration-years column"
  )

  study_col <- detect_column(
    df,
    preferred = c("external_Study", "Study", "STUDY"),
    regex = "(^|_)study$",
    label = "STUDY column"
  )

  site_col <- detect_column(
    df,
    preferred = c("external_SITE", "SITE", "Site"),
    regex = "(^|_)site$",
    label = "SITE column"
  )

  age_col <- detect_column(
    df,
    preferred = c("Age", "age_at_scan_used_for_model", "AGE"),
    regex = "(^|_)age($|_at_scan)",
    label = "age column"
  )

  id_col <- detect_column(
    df,
    preferred = c("PTID", "participant_id", "IID", "eid"),
    regex = "(^ptid$|participant.*id|^iid$|^eid$)",
    label = "participant ID column"
  )

  date_col <- detect_column(
    df,
    preferred = c("Date", "scan_date", "MRI_Date"),
    regex = "(^|_)date$",
    label = "scan date column",
    required = FALSE
  )

  visit_col <- detect_column(
    df,
    preferred = c("Visit_Code", "VISCODE", "visit_code", "Visit"),
    regex = "(visit|viscode)",
    label = "visit column",
    required = FALSE
  )

  baseline_flag_col <- detect_column(
    df,
    preferred = c(
      "is_external_baseline_scan",
      "is_baseline_scan"
    ),
    regex = "baseline.*scan",
    label = "baseline scan indicator",
    required = FALSE
  )

  scan_number_col <- detect_column(
    df,
    preferred = c(
      "longitudinal_scan_number",
      "scan_number"
    ),
    regex = "scan.*number",
    label = "longitudinal scan number",
    required = FALSE
  )

  years_col <- detect_column(
    df,
    preferred = c(
      "years_since_external_baseline",
      "Delta_Baseline",
      "years_since_baseline"
    ),
    regex = "(years.*baseline|delta_baseline)",
    label = "years since baseline",
    required = FALSE
  )

  df |>
    transmute(
      participant_id = as.character(.data[[id_col]]),
      study = clean_character(.data[[study_col]]),
      site = clean_character(.data[[site_col]]),
      age = suppressWarnings(as.numeric(.data[[age_col]])),
      acceleration_years = suppressWarnings(
        as.numeric(.data[[accel_col]])
      ),
      scan_date = if (!is.na(date_col)) {
        suppressWarnings(as.Date(.data[[date_col]]))
      } else {
        as.Date(NA)
      },
      visit_code = if (!is.na(visit_col)) {
        as.character(.data[[visit_col]])
      } else {
        NA_character_
      },
      baseline_flag = if (!is.na(baseline_flag_col)) {
        as.logical(.data[[baseline_flag_col]])
      } else {
        NA
      },
      scan_number = if (!is.na(scan_number_col)) {
        suppressWarnings(as.numeric(.data[[scan_number_col]]))
      } else {
        NA_real_
      },
      years_since_baseline = if (!is.na(years_col)) {
        suppressWarnings(as.numeric(.data[[years_col]]))
      } else {
        NA_real_
      },
      feature_type = feature_type,
      dataset_role = "External"
    ) |>
    filter(
      !is.na(participant_id),
      !is.na(acceleration_years),
      !is.na(age)
    )
}

select_one_external_baseline <- function(df) {
  df |>
    mutate(
      visit_month = visit_to_month(visit_code),
      explicit_baseline_visit = visit_month == 0,
      baseline_priority = case_when(
        baseline_flag %in% TRUE ~ 1,
        scan_number == 1 ~ 2,
        !is.na(years_since_baseline) ~ 3,
        explicit_baseline_visit %in% TRUE ~ 4,
        !is.na(scan_date) ~ 5,
        TRUE ~ 6
      ),
      baseline_distance = case_when(
        baseline_priority == 1 ~ 0,
        baseline_priority == 2 ~ 0,
        baseline_priority == 3 ~ abs(years_since_baseline),
        baseline_priority == 4 ~ abs(visit_month),
        baseline_priority == 5 ~ as.numeric(scan_date),
        TRUE ~ row_number()
      )
    ) |>
    arrange(
      participant_id,
      baseline_priority,
      baseline_distance,
      scan_date,
      visit_month
    ) |>
    group_by(participant_id) |>
    slice_head(n = 1) |>
    ungroup() |>
    select(
      -visit_month,
      -explicit_baseline_visit,
      -baseline_priority,
      -baseline_distance
    )
}

read_adni_reference <- function(path) {
  if (!file.exists(path)) {
    stop("ADNI reference file does not exist: ", path)
  }

  df <- readr::read_tsv(
    path,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "unique"
  )

  accel_col <- detect_column(
    df,
    preferred = c(
      "adni_brain_mri_ad_lepoch_acceleration_years",
      "adni_brain_mri_ad_epoch_acceleration_years"
    ),
    regex = "acceleration[_\\.]*years$",
    label = "ADNI acceleration-years column"
  )

  id_col <- detect_column(
    df,
    preferred = c("PTID", "participant_id"),
    regex = "(^ptid$|participant.*id)",
    label = "ADNI participant ID"
  )

  age_col <- detect_column(
    df,
    preferred = c("Age", "AGE"),
    regex = "^age$",
    label = "ADNI age"
  )

  df |>
    transmute(
      participant_id = as.character(.data[[id_col]]),
      study = "ADNI",
      site = "Training reference",
      age = suppressWarnings(as.numeric(.data[[age_col]])),
      acceleration_years = suppressWarnings(
        as.numeric(.data[[accel_col]])
      ),
      scan_date = as.Date(NA),
      visit_code = NA_character_,
      baseline_flag = TRUE,
      scan_number = 1,
      years_since_baseline = 0,
      feature_type = "ADNI reference",
      dataset_role = "Reference"
    ) |>
    filter(
      !is.na(participant_id),
      !is.na(acceleration_years),
      !is.na(age)
    ) |>
    distinct(participant_id, .keep_all = TRUE)
}

# ------------------------------------------------------------------------------
# Read and prepare data
# ------------------------------------------------------------------------------

raw_baseline <- read_external_file(
  raw_file,
  "Raw features"
) |>
  select_one_external_baseline()

harmonized_baseline <- read_external_file(
  harmonized_file,
  "Harmonized features"
) |>
  select_one_external_baseline()

adni_reference <- read_adni_reference(adni_reference_file)

# Match raw and harmonized external participants.
matched_ids <- inner_join(
  raw_baseline |> distinct(participant_id),
  harmonized_baseline |> distinct(participant_id),
  by = "participant_id"
)

raw_matched <- raw_baseline |>
  semi_join(matched_ids, by = "participant_id")

harmonized_matched <- harmonized_baseline |>
  semi_join(matched_ids, by = "participant_id")

if (nrow(raw_matched) == 0 || nrow(harmonized_matched) == 0) {
  stop("No matched external baseline participants remained.")
}

# Check STUDY/SITE consistency between raw and harmonized files.
metadata_check <- full_join(
  raw_matched |>
    select(
      participant_id,
      raw_study = study,
      raw_site = site
    ),
  harmonized_matched |>
    select(
      participant_id,
      harmonized_study = study,
      harmonized_site = site
    ),
  by = "participant_id"
)

metadata_mismatches <- metadata_check |>
  filter(
    raw_study != harmonized_study |
      raw_site != harmonized_site
  )

if (nrow(metadata_mismatches) > 0) {
  warning(
    nrow(metadata_mismatches),
    " matched participants have differing STUDY/SITE labels across files. ",
    "A mismatch table has been written for review."
  )

  readr::write_tsv(
    metadata_mismatches,
    file.path(
      outdir,
      paste0(prefix, "_study_site_metadata_mismatches.tsv")
    )
  )
}

external_combined <- bind_rows(
  raw_matched,
  harmonized_matched
)

plot_data <- bind_rows(
  adni_reference,
  external_combined
)

# ------------------------------------------------------------------------------
# Summaries
# ------------------------------------------------------------------------------

external_site_summary <- external_combined |>
  group_by(feature_type, study, site) |>
  summarise(
    n_participants = n_distinct(participant_id),
    median_age = median(age, na.rm = TRUE),
    age_p05 = quantile(age, 0.05, na.rm = TRUE),
    age_q25 = quantile(age, 0.25, na.rm = TRUE),
    age_q75 = quantile(age, 0.75, na.rm = TRUE),
    age_p95 = quantile(age, 0.95, na.rm = TRUE),
    mean_acceleration_years = mean(
      acceleration_years,
      na.rm = TRUE
    ),
    sd_acceleration_years = sd(
      acceleration_years,
      na.rm = TRUE
    ),
    median_acceleration_years = median(
      acceleration_years,
      na.rm = TRUE
    ),
    acceleration_q25 = quantile(
      acceleration_years,
      0.25,
      na.rm = TRUE
    ),
    acceleration_q75 = quantile(
      acceleration_years,
      0.75,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

adni_summary <- adni_reference |>
  summarise(
    feature_type = "ADNI reference",
    study = "ADNI",
    site = "Training reference",
    n_participants = n_distinct(participant_id),
    median_age = median(age, na.rm = TRUE),
    age_p05 = quantile(age, 0.05, na.rm = TRUE),
    age_q25 = quantile(age, 0.25, na.rm = TRUE),
    age_q75 = quantile(age, 0.75, na.rm = TRUE),
    age_p95 = quantile(age, 0.95, na.rm = TRUE),
    mean_acceleration_years = mean(
      acceleration_years,
      na.rm = TRUE
    ),
    sd_acceleration_years = sd(
      acceleration_years,
      na.rm = TRUE
    ),
    median_acceleration_years = median(
      acceleration_years,
      na.rm = TRUE
    ),
    acceleration_q25 = quantile(
      acceleration_years,
      0.25,
      na.rm = TRUE
    ),
    acceleration_q75 = quantile(
      acceleration_years,
      0.75,
      na.rm = TRUE
    )
  )

site_summary <- bind_rows(
  adni_summary,
  external_site_summary
)

# Age summaries use one row per participant/site. For external participants,
# raw metadata are used because the matched sample is identical.
age_summary <- bind_rows(
  adni_reference |>
    select(
      participant_id,
      study,
      site,
      age
    ),
  raw_matched |>
    select(
      participant_id,
      study,
      site,
      age
    )
) |>
  group_by(study, site) |>
  summarise(
    n_participants = n_distinct(participant_id),
    median_age = median(age, na.rm = TRUE),
    age_p05 = quantile(age, 0.05, na.rm = TRUE),
    age_q25 = quantile(age, 0.25, na.rm = TRUE),
    age_q75 = quantile(age, 0.75, na.rm = TRUE),
    age_p95 = quantile(age, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

# Domain-shift metrics compare only external sites.
site_mean_metrics <- external_site_summary |>
  group_by(feature_type) |>
  summarise(
    n_sites = n(),
    unweighted_sd_of_site_means = sd(
      mean_acceleration_years,
      na.rm = TRUE
    ),
    weighted_sd_of_site_means = weighted_sd(
      mean_acceleration_years,
      n_participants
    ),
    range_of_site_means = diff(
      range(mean_acceleration_years, na.rm = TRUE)
    ),
    mean_absolute_site_mean = mean(
      abs(mean_acceleration_years),
      na.rm = TRUE
    ),
    .groups = "drop"
  )

raw_metrics <- site_mean_metrics |>
  filter(feature_type == "Raw features")

harmonized_metrics <- site_mean_metrics |>
  filter(feature_type == "Harmonized features")

domain_shift_metrics <- tibble(
  metric = c(
    "Unweighted SD of external SITE means",
    "Participant-weighted SD of external SITE means",
    "Range of external SITE means",
    "Mean absolute external SITE mean"
  ),
  raw = c(
    raw_metrics$unweighted_sd_of_site_means,
    raw_metrics$weighted_sd_of_site_means,
    raw_metrics$range_of_site_means,
    raw_metrics$mean_absolute_site_mean
  ),
  harmonized = c(
    harmonized_metrics$unweighted_sd_of_site_means,
    harmonized_metrics$weighted_sd_of_site_means,
    harmonized_metrics$range_of_site_means,
    harmonized_metrics$mean_absolute_site_mean
  )
) |>
  mutate(
    absolute_reduction = raw - harmonized,
    percent_reduction = 100 * (raw - harmonized) / raw
  )

readr::write_tsv(
  site_summary,
  file.path(
    outdir,
    paste0(prefix, "_baseline_site_summary.tsv")
  )
)

readr::write_tsv(
  domain_shift_metrics,
  file.path(
    outdir,
    paste0(prefix, "_domain_shift_metrics.tsv")
  )
)

readr::write_tsv(
  tibble(
    n_adni_reference_participants = nrow(adni_reference),
    n_raw_external_baseline = nrow(raw_baseline),
    n_harmonized_external_baseline = nrow(harmonized_baseline),
    n_matched_external_baseline_participants = nrow(matched_ids),
    n_study_site_metadata_mismatches = nrow(metadata_mismatches)
  ),
  file.path(
    outdir,
    paste0(prefix, "_sample_summary.tsv")
  )
)

# ------------------------------------------------------------------------------
# Shared row order: ADNI reference on top, followed by STUDY and SITE
# ------------------------------------------------------------------------------

external_rows <- age_summary |>
  filter(study != "ADNI") |>
  arrange(study, site) |>
  transmute(row_label = paste0(study, "  |  ", site)) |>
  pull(row_label)

desired_top_to_bottom <- c(
  "ADNI  |  Training reference",
  external_rows
)

factor_levels_bottom_to_top <- rev(desired_top_to_bottom)

plot_data <- plot_data |>
  mutate(
    study_site = factor(
      paste0(study, "  |  ", site),
      levels = factor_levels_bottom_to_top
    ),
    feature_type = factor(
      feature_type,
      levels = c(
        "Raw features",
        "Harmonized features",
        "ADNI reference"
      )
    )
  )

site_summary <- site_summary |>
  mutate(
    study_site = factor(
      paste0(study, "  |  ", site),
      levels = factor_levels_bottom_to_top
    ),
    feature_type = factor(
      feature_type,
      levels = c(
        "Raw features",
        "Harmonized features",
        "ADNI reference"
      )
    )
  )

age_summary <- age_summary |>
  mutate(
    study_site = factor(
      paste0(study, "  |  ", site),
      levels = factor_levels_bottom_to_top
    )
  )

# Add fallback colors if an unexpected study appears.
missing_studies <- setdiff(
  unique(plot_data$study),
  names(van_gogh_study_colors)
)

if (length(missing_studies) > 0) {
  fallback <- scales::hue_pal()(length(missing_studies))
  names(fallback) <- missing_studies
  van_gogh_study_colors <- c(
    van_gogh_study_colors,
    fallback
  )
}

# ------------------------------------------------------------------------------
# Distribution panel
# ------------------------------------------------------------------------------

raw_data <- plot_data |>
  filter(feature_type == "Raw features")

harmonized_data <- plot_data |>
  filter(feature_type == "Harmonized features")

adni_data <- plot_data |>
  filter(feature_type == "ADNI reference")

raw_means <- site_summary |>
  filter(feature_type == "Raw features")

harmonized_means <- site_summary |>
  filter(feature_type == "Harmonized features")

adni_mean <- site_summary |>
  filter(feature_type == "ADNI reference")

distribution_plot <- ggplot() +
  # ADNI training-reference distribution at the top.
  geom_violin(
    data = adni_data,
    aes(
      x = acceleration_years,
      y = study_site,
      group = study_site
    ),
    fill = reference_fill,
    color = reference_outline,
    width = 0.72,
    scale = "width",
    trim = TRUE,
    alpha = 0.78,
    linewidth = 0.6
  ) +
  geom_boxplot(
    data = adni_data,
    aes(
      x = acceleration_years,
      y = study_site,
      group = study_site
    ),
    width = 0.14,
    fill = "white",
    color = reference_outline,
    outlier.shape = NA,
    linewidth = 0.45
  ) +
  # Raw features: unfilled/outline-only violin.
  geom_violin(
    data = raw_data,
    aes(
      x = acceleration_years,
      y = study_site,
      color = study,
      group = interaction(study_site, feature_type)
    ),
    fill = NA,
    position = position_nudge(y = 0.17),
    width = 0.58,
    scale = "width",
    trim = TRUE,
    linewidth = 0.75
  ) +
  geom_boxplot(
    data = raw_data,
    aes(
      x = acceleration_years,
      y = study_site,
      color = study,
      group = interaction(study_site, feature_type)
    ),
    fill = "white",
    position = position_nudge(y = 0.17),
    width = 0.10,
    outlier.shape = NA,
    linewidth = 0.42
  ) +
  # Harmonized features: filled violin using STUDY color.
  geom_violin(
    data = harmonized_data,
    aes(
      x = acceleration_years,
      y = study_site,
      fill = study,
      color = study,
      group = interaction(study_site, feature_type)
    ),
    position = position_nudge(y = -0.17),
    width = 0.58,
    scale = "width",
    trim = TRUE,
    alpha = 0.62,
    linewidth = 0.55
  ) +
  geom_boxplot(
    data = harmonized_data,
    aes(
      x = acceleration_years,
      y = study_site,
      color = study,
      group = interaction(study_site, feature_type)
    ),
    fill = "white",
    position = position_nudge(y = -0.17),
    width = 0.10,
    outlier.shape = NA,
    linewidth = 0.42
  ) +
  # Mean markers provide the feature-type legend.
  geom_point(
    data = raw_means,
    aes(
      x = mean_acceleration_years,
      y = study_site,
      shape = feature_type
    ),
    position = position_nudge(y = 0.17),
    size = 2.6,
    stroke = 0.85,
    fill = "white",
    color = "grey15"
  ) +
  geom_point(
    data = harmonized_means,
    aes(
      x = mean_acceleration_years,
      y = study_site,
      shape = feature_type
    ),
    position = position_nudge(y = -0.17),
    size = 2.8,
    stroke = 0.65,
    color = "grey15",
    fill = reference_fill
  ) +
  geom_point(
    data = adni_mean,
    aes(
      x = mean_acceleration_years,
      y = study_site,
      shape = feature_type
    ),
    size = 3.0,
    stroke = 0.75,
    color = reference_outline,
    fill = reference_fill
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey35"
  ) +
  scale_y_discrete(
    limits = factor_levels_bottom_to_top,
    drop = FALSE
  ) +
  scale_color_manual(
    values = van_gogh_study_colors,
    name = "Study"
  ) +
  scale_fill_manual(
    values = van_gogh_study_colors,
    name = "Study"
  ) +
  scale_shape_manual(
    values = c(
      "Raw features" = 21,
      "Harmonized features" = 23,
      "ADNI reference" = 24
    ),
    name = "Feature representation"
  ) +
  guides(
    fill = guide_legend(
      order = 1,
      override.aes = list(alpha = 0.70)
    ),
    color = guide_legend(
      order = 1,
      override.aes = list(fill = NA, linewidth = 1.0)
    ),
    shape = guide_legend(
      order = 2,
      override.aes = list(
        size = 3.4,
        color = "grey15"
      )
    )
  ) +
  labs(
    x = "AD EPOCH acceleration (years)",
    y = NULL,
    title = "Baseline acceleration distributions",
    subtitle = paste0(
      "ADNI training data are shown as the top reference; ",
      "raw features are outline-only and harmonized features are filled"
    )
  ) +
  theme_minimal(base_size = 11.5) +
  theme(
    panel.grid.major.y = element_line(
      color = "grey90",
      linewidth = 0.35
    ),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(
      size = 8.4,
      color = "grey15"
    ),
    axis.text.x = element_text(color = "grey20"),
    plot.title = element_text(
      face = "bold",
      size = 13
    ),
    plot.subtitle = element_text(
      size = 9.4,
      color = "grey35"
    ),
    legend.position = "top",
    legend.box = "vertical",
    legend.justification = "left",
    plot.margin = margin(8, 8, 8, 8)
  )

# ------------------------------------------------------------------------------
# Age-range panel
# ------------------------------------------------------------------------------

age_plot <- ggplot(
  age_summary,
  aes(
    y = study_site,
    color = study
  )
) +
  geom_segment(
    aes(
      x = age_p05,
      xend = age_p95,
      yend = study_site
    ),
    linewidth = 2.7,
    alpha = 0.38,
    lineend = "round"
  ) +
  geom_segment(
    aes(
      x = age_q25,
      xend = age_q75,
      yend = study_site
    ),
    linewidth = 5.4,
    alpha = 0.94,
    lineend = "round"
  ) +
  geom_point(
    aes(x = median_age),
    shape = 21,
    fill = "#F7E8A4",
    size = 3.2,
    stroke = 0.9
  ) +
  geom_text(
    aes(
      x = age_p95,
      y = study_site,
      label = paste0(
        "n=",
        scales::comma(n_participants)
      )
    ),
    hjust = -0.15,
    size = 3.0,
    color = "grey25",
    inherit.aes = FALSE
  ) +
  scale_y_discrete(
    limits = factor_levels_bottom_to_top,
    drop = FALSE
  ) +
  scale_color_manual(
    values = van_gogh_study_colors,
    name = "Study"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.03, 0.20))
  ) +
  labs(
    x = "Baseline age (years)",
    y = NULL,
    title = "Age range",
    subtitle = "Rows are locked to the acceleration panel; thin: 5th–95th percentile; thick: IQR; point: median"
  ) +
  theme_minimal(base_size = 11.5) +
  theme(
    panel.grid.major.y = element_line(
      color = "grey90",
      linewidth = 0.35
    ),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(color = "grey20"),
    plot.title = element_text(
      face = "bold",
      size = 13
    ),
    plot.subtitle = element_text(
      size = 9.4,
      color = "grey35"
    ),
    legend.position = "none",
    plot.margin = margin(8, 25, 8, 2)
  )

# ------------------------------------------------------------------------------
# Domain-shift annotation
# ------------------------------------------------------------------------------

weighted_row <- domain_shift_metrics |>
  filter(
    metric ==
      "Participant-weighted SD of external SITE means"
  )

range_row <- domain_shift_metrics |>
  filter(
    metric ==
      "Range of external SITE means"
  )

annotation_text <- paste0(
  "Between-SITE dispersion of baseline means\n",
  "Weighted SD: ",
  sprintf("%.2f", weighted_row$raw),
  " → ",
  sprintf("%.2f", weighted_row$harmonized),
  " years (",
  sprintf("%.1f", weighted_row$percent_reduction),
  "% reduction)\n",
  "Range: ",
  sprintf("%.2f", range_row$raw),
  " → ",
  sprintf("%.2f", range_row$harmonized),
  " years (",
  sprintf("%.1f", range_row$percent_reduction),
  "% reduction)"
)

annotation_plot <- ggplot() +
  annotate(
    "rect",
    xmin = 0,
    xmax = 1,
    ymin = 0,
    ymax = 1,
    fill = "#F3E6B3",
    color = "#355C9A",
    linewidth = 0.8
  ) +
  annotate(
    "text",
    x = 0.04,
    y = 0.90,
    label = "Does harmonization reduce external domain shift?",
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = 4.2,
    color = "#1B4965"
  ) +
  annotate(
    "text",
    x = 0.04,
    y = 0.70,
    label = annotation_text,
    hjust = 0,
    vjust = 1,
    size = 3.65,
    lineheight = 1.15,
    color = "grey15"
  ) +
  annotate(
    "text",
    x = 0.04,
    y = 0.14,
    label = paste0(
      "ADNI provides the original training-reference distribution. ",
      "A reduction in between-SITE dispersion after harmonization is ",
      "consistent with reduced scanner/cohort shift, while within-SITE ",
      "biological heterogeneity remains visible."
    ),
    hjust = 0,
    vjust = 0,
    size = 3.25,
    lineheight = 1.10,
    color = "grey30"
  ) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1),
    clip = "off"
  ) +
  theme_void() +
  theme(
    plot.margin = margin(4, 8, 8, 8)
  )

# ------------------------------------------------------------------------------
# Assemble and save figure
# ------------------------------------------------------------------------------

top_panels <- distribution_plot + age_plot +
  plot_layout(
    widths = c(2.55, 1.05),
    guides = "collect"
  )

final_plot <- top_panels / annotation_plot +
  plot_layout(
    heights = c(4.9, 1.05)
  ) +
  plot_annotation(
    title = "Harmonization and cross-study transport of the AD EPOCH clock",
    subtitle = paste0(
      "Baseline-only matched comparison across ",
      length(unique(external_combined$study)),
      " external studies and ",
      length(unique(external_combined$site)),
      " external scanner/site groups, with ADNI as the training reference"
    ),
    caption = paste0(
      "Acceleration years are calculated using the ADNI-trained model. ",
      "The external raw and harmonized analyses use the same matched baseline participants."
    ),
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 16,
        color = "#1B4965"
      ),
      plot.subtitle = element_text(
        size = 11,
        color = "grey30"
      ),
      plot.caption = element_text(
        size = 8.5,
        color = "grey40"
      )
    )
  ) &
  theme(
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    panel.background = element_rect(
      fill = "white",
      color = NA
    )
  )

figure_height <- max(
  9.2,
  0.50 * length(desired_top_to_bottom) + 4.5
)

pdf_file <- file.path(
  outdir,
  paste0(
    prefix,
    "_distribution_age_range.pdf"
  )
)

png_file <- file.path(
  outdir,
  paste0(
    prefix,
    "_distribution_age_range.png"
  )
)

ggsave(
  filename = pdf_file,
  plot = final_plot,
  width = 15.8,
  height = figure_height,
  units = "in",
  device = cairo_pdf,
  limitsize = FALSE
)

ggsave(
  filename = png_file,
  plot = final_plot,
  width = 15.8,
  height = figure_height,
  units = "in",
  dpi = 350,
  bg = "white",
  limitsize = FALSE
)

message("============================================================")
message("Figure complete.")
message("ADNI reference participants: ", nrow(adni_reference))
message("Matched external baseline participants: ", nrow(matched_ids))
message("PDF: ", pdf_file)
message("PNG: ", png_file)
message("============================================================")