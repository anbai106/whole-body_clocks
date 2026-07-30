#!/usr/bin/env Rscript

# ==============================================================================
# ADNI-DOD demographic summary for the brain MRI AD EPOCH analytic sample
# ==============================================================================

prediction_file <- paste0(
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/adni_lepoch/",
  "results_external_longitudinal_ad_epoch_raw/",
  "external_5_studies_adni_brain_mri_ad_epoch_raw_scan_level_predictions.tsv"
)

istaging_file <- paste0(
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/adni_lepoch/",
  "external_5_studies_istaging.tsv"
)

out_dir <- paste0(
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/adni_lepoch/",
  "results_external_longitudinal_ad_epoch_comparison/ADNI_DOD_demographics"
)

prefix <- "ADNI_DOD_brain_MRI_AD_EPOCH"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("readr", "dplyr", "stringr", "tidyr", "tibble")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(tibble)
})

detect_column <- function(df, preferred, regex, label, required = TRUE) {
  direct <- preferred[preferred %in% names(df)]
  if (length(direct) >= 1) return(direct[[1]])

  candidates <- grep(regex, names(df), value = TRUE, ignore.case = TRUE)
  if (length(candidates) == 1) return(candidates[[1]])
  if (!required) return(NA_character_)
  if (length(candidates) > 1) {
    stop("Multiple candidate columns found for ", label, ": ", paste(candidates, collapse = ", "))
  }
  stop("Could not identify ", label, ". Available columns include:\n",
       paste(head(names(df), 200), collapse = ", "))
}

clean_character <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "nan", "None", "null", "<NA>")] <- NA_character_
  x
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(gsub(",", "", clean_character(x), fixed = TRUE)))
}

normalize_study <- function(x) {
  x <- toupper(clean_character(x))
  x <- gsub("[^A-Z0-9]", "", x)
  ifelse(x == "ADNIDOD", "ADNI-DOD", x)
}

normalize_sex <- function(x) {
  x <- toupper(clean_character(x))
  case_when(
    x %in% c("F", "FEMALE", "WOMAN", "WOMEN", "2") ~ "Female",
    x %in% c("M", "MALE", "MAN", "MEN", "1") ~ "Male",
    TRUE ~ NA_character_
  )
}

parse_date_flexibly <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, c("POSIXct", "POSIXlt"))) return(as.Date(x))
  if (is.numeric(x) || is.integer(x)) return(rep(as.Date(NA_character_), length(x)))

  x <- clean_character(x)
  out <- rep(as.Date(NA_character_), length(x))
  formats <- c("%Y-%m-%d", "%Y-%m-%d %H:%M:%S", "%m/%d/%Y", "%m/%d/%Y %H:%M:%S", "%d-%b-%Y", "%d/%m/%Y")
  for (fmt in formats) {
    idx <- which(is.na(out) & !is.na(x))
    if (length(idx) == 0) break
    parsed <- suppressWarnings(as.Date(x[idx], format = fmt))
    out[idx[!is.na(parsed)]] <- parsed[!is.na(parsed)]
  }
  out
}

format_mean_sd <- function(mean_value, sd_value, digits = 2) {
  paste0(sprintf(paste0("%.", digits, "f"), mean_value), " ± ",
         sprintf(paste0("%.", digits, "f"), sd_value))
}

format_n_pct <- function(n_value, pct_value, digits = 1) {
  paste0(n_value, " (", sprintf(paste0("%.", digits, "f"), pct_value), "%)")
}

if (!file.exists(prediction_file)) stop("Prediction file does not exist: ", prediction_file)

prediction_df <- read_tsv(prediction_file, show_col_types = FALSE, progress = FALSE, name_repair = "unique")

study_col <- detect_column(prediction_df, c("external_Study", "Study", "STUDY"), "(^|_)study$", "study")
id_col <- detect_column(prediction_df, c("PTID", "participant_id", "IID", "RID"), "(^ptid$|participant.*id|^iid$|^rid$)", "participant ID")
age_col <- detect_column(prediction_df, c("Age", "age_at_scan_used_for_model", "AGE"), "(^age$|age.*scan)", "age")
sex_col <- detect_column(prediction_df, c("Sex", "SEX", "Gender", "PTGENDER"), "(^sex$|gender)", "sex", FALSE)
date_col <- detect_column(prediction_df, c("Date", "scan_date", "MRI_Date"), "(^date$|scan.*date|mri.*date)", "scan date", FALSE)
years_col <- detect_column(prediction_df, c("years_since_external_baseline", "years_since_baseline", "Delta_Baseline"), "(years.*baseline|delta_baseline)", "years since baseline", FALSE)
scan_number_col <- detect_column(prediction_df, c("longitudinal_scan_number", "scan_number"), "scan.*number", "scan number", FALSE)
baseline_flag_col <- detect_column(prediction_df, c("is_external_baseline_scan", "is_baseline_scan"), "baseline.*scan", "baseline scan flag", FALSE)

adni_dod_scans <- prediction_df |>
  transmute(
    source_row = row_number(),
    participant_id = as.character(.data[[id_col]]),
    study = normalize_study(.data[[study_col]]),
    age = safe_numeric(.data[[age_col]]),
    sex = if (!is.na(sex_col)) normalize_sex(.data[[sex_col]]) else NA_character_,
    scan_date = if (!is.na(date_col)) parse_date_flexibly(.data[[date_col]]) else as.Date(NA_character_),
    years_since_baseline = if (!is.na(years_col)) safe_numeric(.data[[years_col]]) else NA_real_,
    scan_number = if (!is.na(scan_number_col)) safe_numeric(.data[[scan_number_col]]) else NA_real_,
    baseline_flag = if (!is.na(baseline_flag_col)) as.logical(.data[[baseline_flag_col]]) else NA
  ) |>
  filter(study == "ADNI-DOD", !is.na(participant_id), participant_id != "", is.finite(age))

if (nrow(adni_dod_scans) == 0) stop("No ADNI-DOD scans were found after study filtering.")

if (all(is.na(adni_dod_scans$sex)) && file.exists(istaging_file)) {
  istaging_df <- read_tsv(istaging_file, show_col_types = FALSE, progress = FALSE, name_repair = "unique")
  istudy_col <- detect_column(istaging_df, c("Study", "STUDY"), "(^|_)study$", "iSTAGING study")
  iid_col <- detect_column(istaging_df, c("PTID", "participant_id", "IID", "RID"), "(^ptid$|participant.*id|^iid$|^rid$)", "iSTAGING participant ID")
  isex_col <- detect_column(istaging_df, c("Sex", "SEX", "Gender", "PTGENDER"), "(^sex$|gender)", "iSTAGING sex")

  sex_lookup <- istaging_df |>
    transmute(
      participant_id = as.character(.data[[iid_col]]),
      study = normalize_study(.data[[istudy_col]]),
      sex_fallback = normalize_sex(.data[[isex_col]])
    ) |>
    filter(study == "ADNI-DOD", !is.na(participant_id)) |>
    group_by(participant_id) |>
    summarise(
      sex_fallback = first(sex_fallback[!is.na(sex_fallback)], default = NA_character_),
      .groups = "drop"
    )

  adni_dod_scans <- adni_dod_scans |>
    left_join(sex_lookup, by = "participant_id") |>
    mutate(sex = coalesce(sex, sex_fallback)) |>
    select(-sex_fallback)
}

baseline_demographics <- adni_dod_scans |>
  mutate(
    baseline_priority = case_when(
      baseline_flag %in% TRUE ~ 1,
      scan_number == 1 ~ 2,
      is.finite(years_since_baseline) ~ 3,
      !is.na(scan_date) ~ 4,
      TRUE ~ 5
    ),
    baseline_distance = case_when(
      baseline_priority %in% c(1, 2) ~ 0,
      baseline_priority == 3 ~ abs(years_since_baseline),
      baseline_priority == 4 ~ as.numeric(scan_date),
      TRUE ~ source_row
    )
  ) |>
  arrange(participant_id, baseline_priority, baseline_distance, scan_date, age, source_row) |>
  group_by(participant_id) |>
  slice_head(n = 1) |>
  ungroup() |>
  transmute(
    participant_id,
    baseline_age = age,
    sex,
    baseline_scan_date = scan_date,
    baseline_years_since_baseline = years_since_baseline,
    baseline_selection_priority = baseline_priority
  )

followup_by_participant <- adni_dod_scans |>
  group_by(participant_id) |>
  summarise(
    n_scans = n(),
    min_scan_age = min(age, na.rm = TRUE),
    max_scan_age = max(age, na.rm = TRUE),
    age_based_followup_years = max_scan_age - min_scan_age,
    min_years_since_baseline = if (any(is.finite(years_since_baseline))) min(years_since_baseline, na.rm = TRUE) else NA_real_,
    max_years_since_baseline = if (any(is.finite(years_since_baseline))) max(years_since_baseline, na.rm = TRUE) else NA_real_,
    explicit_followup_years = if (is.finite(min_years_since_baseline) && is.finite(max_years_since_baseline)) max_years_since_baseline - min_years_since_baseline else NA_real_,
    first_scan_date = if (any(!is.na(scan_date))) min(scan_date, na.rm = TRUE) else as.Date(NA_character_),
    last_scan_date = if (any(!is.na(scan_date))) max(scan_date, na.rm = TRUE) else as.Date(NA_character_),
    date_based_followup_years = if (!is.na(first_scan_date) && !is.na(last_scan_date)) as.numeric(last_scan_date - first_scan_date) / 365.25 else NA_real_,
    .groups = "drop"
  ) |>
  mutate(
    followup_years = coalesce(explicit_followup_years, date_based_followup_years, age_based_followup_years),
    followup_source = case_when(
      is.finite(explicit_followup_years) ~ "years_since_baseline",
      is.finite(date_based_followup_years) ~ "scan dates",
      is.finite(age_based_followup_years) ~ "scan age difference",
      TRUE ~ "unavailable"
    )
  )

participant_level <- baseline_demographics |>
  left_join(followup_by_participant, by = "participant_id", relationship = "one-to-one")

n_total <- n_distinct(participant_level$participant_id)
n_with_sex <- sum(!is.na(participant_level$sex))
female_n <- sum(participant_level$sex == "Female", na.rm = TRUE)
female_pct <- if (n_with_sex > 0) 100 * female_n / n_with_sex else NA_real_

followup_values_longitudinal <- participant_level$followup_years[
  is.finite(participant_level$followup_years) & participant_level$n_scans >= 2
]

summary_numeric <- tibble(
  study = "ADNI-DOD",
  country = "USA",
  n_participants = n_total,
  mean_age = mean(participant_level$baseline_age, na.rm = TRUE),
  sd_age = sd(participant_level$baseline_age, na.rm = TRUE),
  age_min = min(participant_level$baseline_age, na.rm = TRUE),
  age_max = max(participant_level$baseline_age, na.rm = TRUE),
  n_with_sex = n_with_sex,
  female_n = female_n,
  female_pct = female_pct,
  male_n = sum(participant_level$sex == "Male", na.rm = TRUE),
  n_missing_sex = sum(is.na(participant_level$sex)),
  n_with_at_least_2_scans = sum(participant_level$n_scans >= 2, na.rm = TRUE),
  median_n_scans = median(participant_level$n_scans, na.rm = TRUE),
  maximum_n_scans = max(participant_level$n_scans, na.rm = TRUE),
  followup_min_longitudinal = if (length(followup_values_longitudinal) > 0) min(followup_values_longitudinal) else NA_real_,
  followup_max_longitudinal = if (length(followup_values_longitudinal) > 0) max(followup_values_longitudinal) else NA_real_,
  median_followup_longitudinal = if (length(followup_values_longitudinal) > 0) median(followup_values_longitudinal) else NA_real_
)

summary_manuscript <- summary_numeric |>
  transmute(
    `Data type` = "Individual",
    `BAG/Omics` = "Brain MRI",
    Study = study,
    Country = country,
    N = n_participants,
    Age = format_mean_sd(mean_age, sd_age, 2),
    `Sex (Female)` = format_n_pct(female_n, female_pct, 1),
    `Event follow-up years` = if (
      is.finite(followup_min_longitudinal) && is.finite(followup_max_longitudinal)
    ) {
      paste0(sprintf("%.1f", followup_min_longitudinal), "–", sprintf("%.1f", followup_max_longitudinal))
    } else {
      NA_character_
    }
  )

write_tsv(participant_level, file.path(out_dir, paste0(prefix, "_participant_level_demographics.tsv")))
write_tsv(summary_numeric, file.path(out_dir, paste0(prefix, "_demographic_summary_numeric.tsv")))
write_tsv(summary_manuscript, file.path(out_dir, paste0(prefix, "_demographic_summary_manuscript.tsv")))

message("============================================================")
message("ADNI-DOD demographic analysis complete.")
message("")
message("Manuscript-ready row:")
print(summary_manuscript)
message("")
message("Detailed numeric summary:")
print(summary_numeric)
message("")
message("Sex counts:")
print(participant_level |> count(sex, name = "n_participants", .drop = FALSE))
message("")
message("Follow-up source counts:")
print(participant_level |> count(followup_source, name = "n_participants", .drop = FALSE))
message("")
message("Outputs saved to: ", out_dir)
message("============================================================")