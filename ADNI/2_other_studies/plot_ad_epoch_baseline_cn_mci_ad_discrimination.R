#!/usr/bin/env Rscript

# ==============================================================================
# Maximized baseline CN/MCI/AD discrimination using harmonized AD EPOCH
#
# Studies:
#   AIBL, BLSA, OASIS
#
# DIAGNOSIS BASELINE
# ------------------
# Derived from the full iSTAGING file:
#   1. Keep rows with an available DX_Binary.
#   2. Within Study + PTID, choose the earliest valid Date.
#   3. If no Date is available, choose the youngest Age.
#
# EPOCH MATCHING
# --------------
# To maximize the usable baseline sample, do NOT require the prediction file's
# internally selected baseline scan. Instead, for each diagnosed participant:
#
#   1. consider every scored harmonized EPOCH scan;
#   2. select the scan closest to the diagnosis-baseline Date;
#   3. if dates cannot be compared, select the scan closest in Age;
#   4. if neither Date nor Age can be compared, use the earliest scored scan.
#
# This retains participants whose original MRI baseline failed ROI coverage but
# who have another scored MRI close to the diagnosis baseline.
#
# IMPORTANT:
# The final sample still cannot exceed participants with at least one scored
# harmonized EPOCH scan.
#
# FIGURE
# ------
#   AIBL:  CN, MCI, AD
#   BLSA:  CN, MCI, AD
#   OASIS: CN, AD only
#
# Cohen's d is shown only when a two-sided Welch t-test has nominal p < 0.05.
# ==============================================================================

# ------------------------------------------------------------------------------
# User settings
# ------------------------------------------------------------------------------

harmonized_file <- paste0(
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/adni_lepoch/",
  "results_external_longitudinal_ad_epoch_harmonized/",
  "external_5_studies_adni_brain_mri_ad_epoch_harmonized_scan_level_predictions.tsv"
)

istaging_file <- paste0(
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/adni_lepoch/",
  "external_5_studies_istaging.tsv"
)

outdir <- paste0(
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/adni_lepoch/",
  "results_external_longitudinal_ad_epoch_comparison"
)

prefix <- paste0(
  "external_AIBL_BLSA_OASIS_",
  "maximized_diagnosis_baseline_harmonized_epoch"
)

study_order <- c("AIBL", "BLSA", "OASIS")
diagnosis_order <- c("CN", "MCI", "AD")

min_group_n_for_test <- 2
annotation_p_threshold <- 0.05

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

message("Harmonized prediction input: ", harmonized_file)
message("Full iSTAGING input: ", istaging_file)
message("Output directory: ", outdir)

# ------------------------------------------------------------------------------
# Packages
# ------------------------------------------------------------------------------

required_packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "ggplot2",
  "stringr",
  "scales",
  "patchwork"
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
  library(stringr)
  library(scales)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# Colors
# ------------------------------------------------------------------------------

diagnosis_colors <- c(
  "CN"  = "#355C9A",
  "MCI" = "#E3A018",
  "AD"  = "#B55239"
)

# ------------------------------------------------------------------------------
# General helpers
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
      "Multiple candidate columns found for ",
      label,
      ": ",
      paste(candidates, collapse = ", ")
    )
  }
  
  stop(
    "Could not identify ",
    label,
    ". Available columns include:\n",
    paste(head(names(df), 160), collapse = ", ")
  )
}

clean_optional_character <- function(x) {
  x <- trimws(as.character(x))
  
  x[
    x %in% c(
      "",
      "NA",
      "NaN",
      "nan",
      "None",
      "null",
      "<NA>"
    )
  ] <- NA_character_
  
  x
}

clean_character <- function(
    x,
    missing_label = "Unknown"
) {
  tidyr::replace_na(
    clean_optional_character(x),
    missing_label
  )
}

normalize_study <- function(x) {
  x <- clean_character(x)
  upper <- toupper(x)
  
  case_when(
    str_detect(upper, "^AIBL") ~ "AIBL",
    str_detect(upper, "^BLSA") ~ "BLSA",
    str_detect(upper, "^OASIS") ~ "OASIS",
    TRUE ~ x
  )
}

normalize_diagnosis <- function(x) {
  upper <- toupper(
    clean_optional_character(x)
  )
  
  case_when(
    upper %in% c(
      "CN",
      "NC",
      "NORMAL",
      "COGNITIVELY NORMAL",
      "COGNITIVE NORMAL",
      "CONTROL",
      "HEALTHY CONTROL",
      "HC",
      "0"
    ) ~ "CN",
    
    upper %in% c(
      "MCI",
      "LMCI",
      "EMCI",
      "EARLY MCI",
      "MILD COGNITIVE IMPAIRMENT",
      "1"
    ) ~ "MCI",
    
    upper %in% c(
      "AD",
      "DEMENTIA",
      "ALZHEIMER",
      "ALZHEIMER'S DISEASE",
      "ALZHEIMERS DISEASE",
      "ALZHEIMER DISEASE",
      "2"
    ) ~ "AD",
    
    str_detect(
      upper,
      "(^|[^A-Z])MCI([^A-Z]|$)"
    ) ~ "MCI",
    
    str_detect(
      upper,
      "ALZHEIMER|DEMENTIA"
    ) ~ "AD",
    
    str_detect(
      upper,
      "COGNITIVELY NORMAL|COGNITIVE NORMAL"
    ) ~ "CN",
    
    TRUE ~ NA_character_
  )
}

parse_date_flexibly <- function(x) {
  x <- clean_optional_character(x)
  
  output <- rep(
    as.Date(NA),
    length(x)
  )
  
  formats <- c(
    "%Y-%m-%d",
    "%Y-%m-%d %H:%M:%S",
    "%m/%d/%Y",
    "%m/%d/%Y %H:%M:%S",
    "%d-%b-%Y",
    "%d/%m/%Y"
  )
  
  for (format_string in formats) {
    unresolved <- is.na(output) & !is.na(x)
    
    if (!any(unresolved)) {
      break
    }
    
    indices <- which(unresolved)
    
    parsed <- suppressWarnings(
      as.Date(
        x[indices],
        format = format_string
      )
    )
    
    valid <- !is.na(parsed)
    
    if (any(valid)) {
      output[
        indices[valid]
      ] <- parsed[valid]
    }
  }
  
  output
}

# ------------------------------------------------------------------------------
# Baseline diagnosis from full iSTAGING
# ------------------------------------------------------------------------------

select_earliest_diagnosed_record <- function(df) {
  df |>
    mutate(
      diagnosis_date_missing = is.na(diagnosis_date),
      
      diagnosis_baseline_source = case_when(
        !diagnosis_date_missing ~
          "Earliest diagnosed Date",
        
        diagnosis_date_missing &
          !is.na(diagnosis_age) ~
          "Youngest diagnosed Age",
        
        TRUE ~
          "Diagnosed row without valid Date or Age"
      )
    ) |>
    arrange(
      study,
      participant_id,
      diagnosis_date_missing,
      diagnosis_date,
      diagnosis_age
    ) |>
    group_by(
      study,
      participant_id
    ) |>
    slice_head(
      n = 1
    ) |>
    ungroup()
}

# ------------------------------------------------------------------------------
# Match every diagnosis baseline to the closest scored EPOCH scan
# ------------------------------------------------------------------------------

select_closest_prediction <- function(
    diagnosis_baseline,
    prediction_scans
) {
  candidate_pairs <- prediction_scans |>
    inner_join(
      diagnosis_baseline,
      by = c(
        "study",
        "participant_id"
      ),
      relationship = "many-to-many"
    ) |>
    mutate(
      date_difference_days = ifelse(
        !is.na(prediction_date) &
          !is.na(diagnosis_date),
        abs(
          as.numeric(
            prediction_date -
              diagnosis_date
          )
        ),
        NA_real_
      ),
      
      age_difference_years = ifelse(
        !is.na(prediction_age) &
          !is.na(diagnosis_age),
        abs(
          prediction_age -
            diagnosis_age
        ),
        NA_real_
      ),
      
      match_priority = case_when(
        is.finite(date_difference_days) ~ 1,
        is.finite(age_difference_years) ~ 2,
        !is.na(prediction_date) ~ 3,
        !is.na(prediction_age) ~ 4,
        TRUE ~ 5
      ),
      
      match_distance = case_when(
        match_priority == 1 ~
          date_difference_days,
        
        match_priority == 2 ~
          age_difference_years,
        
        match_priority == 3 ~
          as.numeric(prediction_date),
        
        match_priority == 4 ~
          prediction_age,
        
        TRUE ~
          prediction_source_row
      ),
      
      epoch_match_source = case_when(
        match_priority == 1 ~
          "Closest scored scan by Date",
        
        match_priority == 2 ~
          "Closest scored scan by Age",
        
        match_priority == 3 ~
          "Earliest scored scan by Date",
        
        match_priority == 4 ~
          "Youngest scored scan by Age",
        
        TRUE ~
          "First available scored scan"
      )
    ) |>
    arrange(
      study,
      participant_id,
      match_priority,
      match_distance,
      prediction_date,
      prediction_age,
      prediction_source_row
    ) |>
    group_by(
      study,
      participant_id
    ) |>
    slice_head(
      n = 1
    ) |>
    ungroup()
  
  candidate_pairs
}

# ------------------------------------------------------------------------------
# Statistical helpers
# ------------------------------------------------------------------------------

calculate_pairwise_effect <- function(
    data,
    group_1,
    group_2
) {
  x <- data |>
    filter(
      diagnosis == group_1
    ) |>
    pull(
      acceleration_years
    )
  
  y <- data |>
    filter(
      diagnosis == group_2
    ) |>
    pull(
      acceleration_years
    )
  
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  
  n_x <- length(x)
  n_y <- length(y)
  
  result <- tibble(
    comparison = paste(
      group_1,
      "vs",
      group_2
    ),
    group_1 = group_1,
    group_2 = group_2,
    n_group_1 = n_x,
    n_group_2 = n_y,
    mean_group_1 = ifelse(
      n_x > 0,
      mean(x),
      NA_real_
    ),
    mean_group_2 = ifelse(
      n_y > 0,
      mean(y),
      NA_real_
    ),
    mean_difference = ifelse(
      n_x > 0 && n_y > 0,
      mean(y) - mean(x),
      NA_real_
    ),
    cohens_d = NA_real_,
    t_statistic = NA_real_,
    degrees_freedom = NA_real_,
    p_value = NA_real_,
    significant = FALSE
  )
  
  if (
    n_x < min_group_n_for_test ||
    n_y < min_group_n_for_test
  ) {
    return(result)
  }
  
  pooled_variance <- (
    (n_x - 1) * stats::var(x) +
      (n_y - 1) * stats::var(y)
  ) / (
    n_x + n_y - 2
  )
  
  if (
    !is.finite(pooled_variance) ||
    pooled_variance <= 0
  ) {
    return(result)
  }
  
  d_value <- (
    mean(y) - mean(x)
  ) / sqrt(
    pooled_variance
  )
  
  test <- tryCatch(
    stats::t.test(
      x,
      y,
      alternative = "two.sided",
      var.equal = FALSE
    ),
    error = function(e) NULL
  )
  
  result$cohens_d <- d_value
  
  if (!is.null(test)) {
    result$t_statistic <- unname(
      test$statistic
    )
    
    result$degrees_freedom <- unname(
      test$parameter
    )
    
    result$p_value <- unname(
      test$p.value
    )
    
    result$significant <- (
      is.finite(result$p_value) &&
        result$p_value <
        annotation_p_threshold
    )
  }
  
  result
}

format_p_value <- function(p) {
  if (is.na(p)) {
    return("p=NA")
  }
  
  if (p < 0.001) {
    return("p<0.001")
  }
  
  paste0(
    "p=",
    sprintf("%.3f", p)
  )
}

# ------------------------------------------------------------------------------
# Validate input files
# ------------------------------------------------------------------------------

if (!file.exists(istaging_file)) {
  stop(
    "iSTAGING input file does not exist: ",
    istaging_file
  )
}

if (!file.exists(harmonized_file)) {
  stop(
    "Harmonized prediction file does not exist: ",
    harmonized_file
  )
}

# ------------------------------------------------------------------------------
# Read and derive baseline diagnosis from full iSTAGING
# ------------------------------------------------------------------------------

istaging_df <- readr::read_tsv(
  istaging_file,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "unique"
)

istaging_study_col <- detect_column(
  istaging_df,
  preferred = c(
    "Study",
    "STUDY"
  ),
  regex = "(^|_)study$",
  label = "iSTAGING Study"
)

istaging_id_col <- detect_column(
  istaging_df,
  preferred = c(
    "PTID",
    "participant_id",
    "IID"
  ),
  regex = "(^ptid$|participant.*id|^iid$)",
  label = "iSTAGING participant ID"
)

istaging_dx_col <- detect_column(
  istaging_df,
  preferred = c(
    "DX_Binary",
    "Dx_binary",
    "dx_binary"
  ),
  regex = "(^|_)dx[_\\.]*binary$",
  label = "iSTAGING DX_Binary"
)

istaging_date_col <- detect_column(
  istaging_df,
  preferred = c(
    "Date",
    "scan_date",
    "MRI_Date"
  ),
  regex = "(^|_)date$",
  label = "iSTAGING Date",
  required = FALSE
)

istaging_age_col <- detect_column(
  istaging_df,
  preferred = c(
    "Age",
    "AGE"
  ),
  regex = "^age$",
  label = "iSTAGING Age"
)

diagnosis_rows <- istaging_df |>
  transmute(
    study = normalize_study(
      .data[[istaging_study_col]]
    ),
    
    participant_id = as.character(
      .data[[istaging_id_col]]
    ),
    
    diagnosis_original = clean_optional_character(
      .data[[istaging_dx_col]]
    ),
    
    diagnosis = normalize_diagnosis(
      .data[[istaging_dx_col]]
    ),
    
    diagnosis_date = if (!is.na(
      istaging_date_col
    )) {
      parse_date_flexibly(
        .data[[istaging_date_col]]
      )
    } else {
      as.Date(NA)
    },
    
    diagnosis_age = suppressWarnings(
      as.numeric(
        .data[[istaging_age_col]]
      )
    )
  ) |>
  filter(
    study %in% study_order,
    !is.na(participant_id),
    participant_id != "",
    !is.na(diagnosis_original)
  )

baseline_diagnosis_all_labels <- diagnosis_rows |>
  select_earliest_diagnosed_record()

baseline_diagnosis <- baseline_diagnosis_all_labels |>
  filter(
    !is.na(diagnosis)
  )

diagnosis_counts_full_istaging <- baseline_diagnosis |>
  count(
    study,
    diagnosis,
    name = "n_diagnosed_participants"
  ) |>
  tidyr::complete(
    study = study_order,
    diagnosis = diagnosis_order,
    fill = list(
      n_diagnosed_participants = 0
    )
  ) |>
  arrange(
    study,
    diagnosis
  )

# ------------------------------------------------------------------------------
# Read all scored harmonized EPOCH scans
# ------------------------------------------------------------------------------

prediction_df <- readr::read_tsv(
  harmonized_file,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "unique"
)

prediction_acceleration_col <- detect_column(
  prediction_df,
  preferred = c(
    "adni_brain_mri_ad_epoch_acceleration_years",
    "adni_brain_mri_ad_lepoch_acceleration_years"
  ),
  regex = "acceleration[_\\.]*years$",
  label = "prediction acceleration-years"
)

prediction_study_col <- detect_column(
  prediction_df,
  preferred = c(
    "external_Study",
    "Study",
    "STUDY"
  ),
  regex = "(^|_)study$",
  label = "prediction Study"
)

prediction_id_col <- detect_column(
  prediction_df,
  preferred = c(
    "PTID",
    "participant_id",
    "IID"
  ),
  regex = "(^ptid$|participant.*id|^iid$)",
  label = "prediction participant ID"
)

prediction_site_col <- detect_column(
  prediction_df,
  preferred = c(
    "external_SITE",
    "SITE",
    "Site"
  ),
  regex = "(^|_)site$",
  label = "prediction SITE",
  required = FALSE
)

prediction_date_col <- detect_column(
  prediction_df,
  preferred = c(
    "Date",
    "scan_date",
    "MRI_Date"
  ),
  regex = "(^|_)date$",
  label = "prediction Date",
  required = FALSE
)

prediction_age_col <- detect_column(
  prediction_df,
  preferred = c(
    "Age",
    "age_at_scan_used_for_model",
    "AGE"
  ),
  regex = "(^|_)age($|_at_scan)",
  label = "prediction Age"
)

prediction_rows <- prediction_df |>
  transmute(
    prediction_source_row = row_number(),
    
    study = normalize_study(
      .data[[prediction_study_col]]
    ),
    
    participant_id = as.character(
      .data[[prediction_id_col]]
    ),
    
    site = if (!is.na(
      prediction_site_col
    )) {
      clean_character(
        .data[[prediction_site_col]]
      )
    } else {
      "Unknown"
    },
    
    prediction_date = if (!is.na(
      prediction_date_col
    )) {
      parse_date_flexibly(
        .data[[prediction_date_col]]
      )
    } else {
      as.Date(NA)
    },
    
    prediction_age = suppressWarnings(
      as.numeric(
        .data[[prediction_age_col]]
      )
    ),
    
    acceleration_years = suppressWarnings(
      as.numeric(
        .data[[prediction_acceleration_col]]
      )
    )
  ) |>
  filter(
    study %in% study_order,
    !is.na(participant_id),
    participant_id != "",
    !is.na(acceleration_years)
  )

prediction_subject_counts <- prediction_rows |>
  group_by(
    study
  ) |>
  summarise(
    n_scored_scans = n(),
    n_scored_participants = n_distinct(
      participant_id
    ),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# Maximize diagnosed participants with a harmonized EPOCH scan
# ------------------------------------------------------------------------------

matched_baseline <- select_closest_prediction(
  diagnosis_baseline = baseline_diagnosis,
  prediction_scans = prediction_rows
) |>
  mutate(
    study = factor(
      study,
      levels = study_order
    ),
    
    diagnosis = factor(
      diagnosis,
      levels = diagnosis_order
    )
  )

if (nrow(matched_baseline) == 0) {
  stop(
    "No diagnosed iSTAGING participants matched to harmonized predictions."
  )
}

# OASIS: display CN and AD only, as requested previously.
oasis_mci_rows <- matched_baseline |>
  filter(
    as.character(study) == "OASIS",
    diagnosis == "MCI"
  )

if (nrow(oasis_mci_rows) > 0) {
  readr::write_tsv(
    oasis_mci_rows,
    file.path(
      outdir,
      paste0(
        prefix,
        "_excluded_OASIS_MCI_rows.tsv"
      )
    )
  )
}

plot_data <- matched_baseline |>
  filter(
    !(
      as.character(study) == "OASIS" &
        diagnosis == "MCI"
    )
  )

# ------------------------------------------------------------------------------
# Matching audit
# ------------------------------------------------------------------------------

diagnosis_subjects <- baseline_diagnosis |>
  distinct(
    study,
    participant_id,
    diagnosis
  )

prediction_subjects <- prediction_rows |>
  distinct(
    study,
    participant_id
  )

matching_summary <- diagnosis_subjects |>
  mutate(
    study = factor(
      study,
      levels = study_order
    )
  ) |>
  group_by(
    study
  ) |>
  summarise(
    n_diagnosed_istaging_participants = n_distinct(
      participant_id
    ),
    
    .groups = "drop"
  ) |>
  left_join(
    prediction_subject_counts,
    by = "study"
  ) |>
  left_join(
    matched_baseline |>
      group_by(
        study
      ) |>
      summarise(
        n_diagnosed_with_scored_epoch = n_distinct(
          participant_id
        ),
        
        .groups = "drop"
      ),
    by = "study"
  ) |>
  mutate(
    n_diagnosed_without_scored_epoch = (
      n_diagnosed_istaging_participants -
        n_diagnosed_with_scored_epoch
    ),
    
    percent_diagnosed_retained = (
      100 *
        n_diagnosed_with_scored_epoch /
        n_diagnosed_istaging_participants
    )
  )

match_source_summary <- matched_baseline |>
  count(
    study,
    epoch_match_source,
    name = "n_participants"
  ) |>
  arrange(
    study,
    epoch_match_source
  )

match_distance_summary <- matched_baseline |>
  group_by(
    study,
    epoch_match_source
  ) |>
  summarise(
    n = n(),
    
    median_date_difference_days = median(
      date_difference_days,
      na.rm = TRUE
    ),
    
    q95_date_difference_days = quantile(
      date_difference_days,
      0.95,
      na.rm = TRUE
    ),
    
    median_age_difference_years = median(
      age_difference_years,
      na.rm = TRUE
    ),
    
    q95_age_difference_years = quantile(
      age_difference_years,
      0.95,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# Final counts and descriptive summaries
# ------------------------------------------------------------------------------

diagnosis_counts <- plot_data |>
  count(
    study,
    diagnosis,
    name = "n"
  ) |>
  tidyr::complete(
    study = factor(
      study_order,
      levels = study_order
    ),
    
    diagnosis = factor(
      diagnosis_order,
      levels = diagnosis_order
    ),
    
    fill = list(
      n = 0
    )
  ) |>
  arrange(
    study,
    diagnosis
  )

study_summary <- diagnosis_counts |>
  tidyr::pivot_wider(
    names_from = diagnosis,
    values_from = n,
    values_fill = 0
  ) |>
  mutate(
    total_plotted = CN + MCI + AD
  ) |>
  left_join(
    matching_summary,
    by = "study"
  )

site_summary <- plot_data |>
  group_by(
    study,
    site,
    diagnosis
  ) |>
  summarise(
    n = n_distinct(
      participant_id
    ),
    .groups = "drop"
  ) |>
  arrange(
    study,
    site,
    diagnosis
  )

descriptive_summary <- plot_data |>
  group_by(
    study,
    diagnosis
  ) |>
  summarise(
    n = n_distinct(
      participant_id
    ),
    
    mean_prediction_age = mean(
      prediction_age,
      na.rm = TRUE
    ),
    
    sd_prediction_age = sd(
      prediction_age,
      na.rm = TRUE
    ),
    
    median_prediction_age = median(
      prediction_age,
      na.rm = TRUE
    ),
    
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
    
    q25_acceleration_years = quantile(
      acceleration_years,
      0.25,
      na.rm = TRUE
    ),
    
    q75_acceleration_years = quantile(
      acceleration_years,
      0.75,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# Pairwise Cohen's d and Welch tests
# ------------------------------------------------------------------------------

pairwise_results <- bind_rows(
  calculate_pairwise_effect(
    plot_data |>
      filter(
        as.character(study) == "AIBL"
      ),
    "CN",
    "MCI"
  ) |>
    mutate(
      study = "AIBL",
      .before = 1
    ),
  
  calculate_pairwise_effect(
    plot_data |>
      filter(
        as.character(study) == "AIBL"
      ),
    "MCI",
    "AD"
  ) |>
    mutate(
      study = "AIBL",
      .before = 1
    ),
  
  calculate_pairwise_effect(
    plot_data |>
      filter(
        as.character(study) == "BLSA"
      ),
    "CN",
    "MCI"
  ) |>
    mutate(
      study = "BLSA",
      .before = 1
    ),
  
  calculate_pairwise_effect(
    plot_data |>
      filter(
        as.character(study) == "BLSA"
      ),
    "MCI",
    "AD"
  ) |>
    mutate(
      study = "BLSA",
      .before = 1
    ),
  
  calculate_pairwise_effect(
    plot_data |>
      filter(
        as.character(study) == "OASIS"
      ),
    "CN",
    "AD"
  ) |>
    mutate(
      study = "OASIS",
      .before = 1
    )
) |>
  mutate(
    study = factor(
      study,
      levels = study_order
    )
  )

# ------------------------------------------------------------------------------
# Write all QC and result tables
# ------------------------------------------------------------------------------

readr::write_tsv(
  diagnosis_counts_full_istaging,
  file.path(
    outdir,
    paste0(
      prefix,
      "_full_istaging_diagnosis_counts.tsv"
    )
  )
)

readr::write_tsv(
  prediction_subject_counts,
  file.path(
    outdir,
    paste0(
      prefix,
      "_prediction_subject_counts.tsv"
    )
  )
)

readr::write_tsv(
  matching_summary,
  file.path(
    outdir,
    paste0(
      prefix,
      "_matching_summary.tsv"
    )
  )
)

readr::write_tsv(
  match_source_summary,
  file.path(
    outdir,
    paste0(
      prefix,
      "_match_source_summary.tsv"
    )
  )
)

readr::write_tsv(
  match_distance_summary,
  file.path(
    outdir,
    paste0(
      prefix,
      "_match_distance_summary.tsv"
    )
  )
)

readr::write_tsv(
  matched_baseline,
  file.path(
    outdir,
    paste0(
      prefix,
      "_selected_maximized_baseline_rows.tsv"
    )
  )
)

readr::write_tsv(
  diagnosis_counts,
  file.path(
    outdir,
    paste0(
      prefix,
      "_diagnosis_counts.tsv"
    )
  )
)

readr::write_tsv(
  study_summary,
  file.path(
    outdir,
    paste0(
      prefix,
      "_study_summary.tsv"
    )
  )
)

readr::write_tsv(
  site_summary,
  file.path(
    outdir,
    paste0(
      prefix,
      "_site_by_diagnosis_counts.tsv"
    )
  )
)

readr::write_tsv(
  descriptive_summary,
  file.path(
    outdir,
    paste0(
      prefix,
      "_descriptive_summary.tsv"
    )
  )
)

readr::write_tsv(
  pairwise_results,
  file.path(
    outdir,
    paste0(
      prefix,
      "_pairwise_cohens_d_welch_t_test.tsv"
    )
  )
)

# ------------------------------------------------------------------------------
# Study panels
# ------------------------------------------------------------------------------

make_study_panel <- function(
    study_name,
    displayed_diagnoses,
    show_y_title = FALSE
) {
  dat <- plot_data |>
    filter(
      as.character(study) == study_name,
      diagnosis %in% displayed_diagnoses
    ) |>
    mutate(
      diagnosis_plot = factor(
        as.character(diagnosis),
        levels = displayed_diagnoses
      )
    )
  
  if (nrow(dat) == 0) {
    stop(
      "No data available for study: ",
      study_name
    )
  }
  
  counts <- dat |>
    count(
      diagnosis_plot,
      name = "n"
    ) |>
    tidyr::complete(
      diagnosis_plot = factor(
        displayed_diagnoses,
        levels = displayed_diagnoses
      ),
      fill = list(
        n = 0
      )
    )
  
  count_text <- paste0(
    as.character(
      counts$diagnosis_plot
    ),
    " n=",
    scales::comma(
      counts$n
    ),
    collapse = " | "
  )
  
  panel_title <- paste0(
    study_name,
    "\n",
    count_text,
    " | total n=",
    scales::comma(
      nrow(dat)
    )
  )
  
  annotation <- pairwise_results |>
    filter(
      as.character(study) == study_name,
      significant,
      is.finite(cohens_d),
      is.finite(p_value)
    ) |>
    mutate(
      label = paste0(
        comparison,
        ": d=",
        sprintf(
          "%.2f",
          cohens_d
        ),
        ", ",
        vapply(
          p_value,
          format_p_value,
          FUN.VALUE = character(1)
        )
      )
    )
  
  y_min <- min(
    dat$acceleration_years,
    na.rm = TRUE
  )
  
  y_max <- max(
    dat$acceleration_years,
    na.rm = TRUE
  )
  
  y_range <- max(
    y_max - y_min,
    1
  )
  
  if (nrow(annotation) > 0) {
    annotation_df <- tibble(
      x = mean(
        seq_along(
          displayed_diagnoses
        )
      ),
      
      y = y_max +
        0.08 * y_range,
      
      label = paste(
        annotation$label,
        collapse = "\n"
      )
    )
  } else {
    annotation_df <- tibble(
      x = numeric(0),
      y = numeric(0),
      label = character(0)
    )
  }
  
  ggplot(
    dat,
    aes(
      x = diagnosis_plot,
      y = acceleration_years,
      fill = diagnosis_plot,
      color = diagnosis_plot
    )
  ) +
    geom_violin(
      trim = TRUE,
      scale = "width",
      width = 0.84,
      alpha = 0.66,
      linewidth = 0.65
    ) +
    geom_boxplot(
      width = 0.16,
      outlier.shape = NA,
      fill = "white",
      linewidth = 0.48
    ) +
    geom_jitter(
      width = 0.10,
      height = 0,
      size = 0.60,
      alpha = 0.17,
      show.legend = FALSE
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 23,
      size = 3.0,
      stroke = 0.70,
      fill = "white",
      color = "grey15"
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.45,
      color = "grey35"
    ) +
    geom_text(
      data = annotation_df,
      aes(
        x = x,
        y = y,
        label = label
      ),
      inherit.aes = FALSE,
      hjust = 0.5,
      vjust = 0,
      size = 3.0,
      lineheight = 1.08,
      color = "#1B4965",
      fontface = "bold"
    ) +
    scale_x_discrete(
      limits = displayed_diagnoses,
      drop = FALSE
    ) +
    scale_fill_manual(
      values = diagnosis_colors,
      limits = diagnosis_order,
      drop = FALSE
    ) +
    scale_color_manual(
      values = diagnosis_colors,
      limits = diagnosis_order,
      drop = FALSE
    ) +
    scale_y_continuous(
      expand = expansion(
        mult = c(
          0.06,
          0.22
        )
      )
    ) +
    labs(
      x = NULL,
      
      y = if (show_y_title) {
        "Harmonized AD EPOCH acceleration (years)"
      } else {
        NULL
      },
      
      title = panel_title
    ) +
    theme_minimal(
      base_size = 11.5
    ) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      
      plot.title = element_text(
        face = "bold",
        size = 10.0,
        color = "#1B4965",
        hjust = 0.5,
        margin = margin(
          b = 8
        )
      ),
      
      axis.text.x = element_text(
        face = "bold",
        size = 10.2
      ),
      
      axis.text.y = element_text(
        color = "grey20"
      ),
      
      axis.title.y = element_text(
        face = "bold"
      ),
      
      legend.position = "none",
      
      plot.margin = margin(
        8,
        8,
        8,
        8
      )
    )
}

p_aibl <- make_study_panel(
  "AIBL",
  c(
    "CN",
    "MCI",
    "AD"
  ),
  show_y_title = TRUE
)

p_blsa <- make_study_panel(
  "BLSA",
  c(
    "CN",
    "MCI",
    "AD"
  ),
  show_y_title = FALSE
)

p_oasis <- make_study_panel(
  "OASIS",
  c(
    "CN",
    "AD"
  ),
  show_y_title = FALSE
)

final_plot <- (
  p_aibl |
    p_blsa |
    p_oasis
) +
  patchwork::plot_layout(
    nrow = 1,
    widths = c(
      1,
      1,
      0.82
    )
  ) +
  patchwork::plot_annotation(
    title = paste0(
      "CN, MCI, and AD distributions using maximized ",
      "diagnosis-baseline matching"
    ),
    
    subtitle = paste0(
      "Each earliest diagnosed iSTAGING participant is matched to the ",
      "closest available harmonized EPOCH scan by Date or Age"
    ),
    
    caption = paste0(
      "OASIS is displayed as CN versus AD. Cohen's d is shown only ",
      "when the corresponding two-sided Welch t-test has nominal p<0.05. ",
      "The selected EPOCH scan may be later than the original MRI baseline ",
      "when the original scan did not pass harmonized ROI coverage QC."
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

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------

pdf_file <- file.path(
  outdir,
  paste0(
    prefix,
    "_one_row.pdf"
  )
)

png_file <- file.path(
  outdir,
  paste0(
    prefix,
    "_one_row.png"
  )
)

ggsave(
  filename = pdf_file,
  plot = final_plot,
  width = 16.0,
  height = 6.6,
  units = "in",
  device = cairo_pdf,
  limitsize = FALSE
)

ggsave(
  filename = png_file,
  plot = final_plot,
  width = 16.0,
  height = 6.6,
  units = "in",
  dpi = 350,
  bg = "white",
  limitsize = FALSE
)

message("============================================================")
message("Analysis complete.")
message("Included studies: ", paste(study_order, collapse = ", "))
message("")
message("Full iSTAGING diagnosed participants:")
print(diagnosis_counts_full_istaging)
message("")
message("Scored harmonized EPOCH availability:")
print(prediction_subject_counts)
message("")
message("Maximum diagnosis-to-EPOCH matching:")
print(matching_summary)
message("")
message("Final plotted diagnosis counts:")
print(study_summary)
message("")
message("EPOCH scan matching sources:")
print(match_source_summary)
message("")
message("Significant pairwise comparisons annotated:")
print(
  pairwise_results |>
    filter(
      significant
    )
)
message("")
message("PDF: ", pdf_file)
message("PNG: ", png_file)
message("============================================================")