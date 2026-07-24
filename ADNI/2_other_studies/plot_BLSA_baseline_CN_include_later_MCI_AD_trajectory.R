#!/usr/bin/env Rscript

# ==============================================================================
# Longitudinal harmonized AD EPOCH trajectories in BLSA participants
# who are CN at baseline, including later MCI/AD sessions
#
# Cohort definition
# -----------------
# 1. Restrict the full iSTAGING dataset to BLSA.
# 2. Keep participants whose earliest mapped diagnosis is CN.
# 3. Retain all later scored harmonized EPOCH scans, including scans occurring
#    after later MCI or AD diagnoses.
#
# Baseline definition
# -------------------
# The baseline diagnosis is the earliest mapped diagnosis record:
#   - earliest valid Date;
#   - youngest Age when Date is unavailable.
#
# Scan-level diagnosis assignment
# -------------------------------
# Each scored EPOCH scan is assigned the closest available longitudinal
# diagnosis record for that participant:
#   - closest diagnosis Date when both dates are available;
#   - otherwise closest diagnosis Age;
#   - otherwise the earliest available diagnosis record.
#
# The scan dots are colored by assigned diagnosis:
#   CN, MCI, or AD.
#
# Main outputs
# ------------
# 1. Individual trajectories with scan points colored by diagnosis.
# 2. Population-level trend across all baseline-CN participants.
# 3. Diagnosis-stratified population trajectories.
# 4. Within-subject change from the first qualified EPOCH scan.
# 5. Subject-level slope distribution by eventual diagnosis group.
# 6. Detailed QC and summary tables.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. User settings
# ------------------------------------------------------------------------------

istaging_file <- paste0(
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/adni_lepoch/",
  "external_5_studies_istaging.tsv"
)

harmonized_prediction_file <- paste0(
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/adni_lepoch/",
  "results_external_longitudinal_ad_epoch_harmonized/",
  "external_5_studies_adni_brain_mri_ad_epoch_harmonized_scan_level_predictions.tsv"
)

out_dir <- paste0(
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/adni_lepoch/",
  "results_external_longitudinal_ad_epoch_comparison/",
  "BLSA_baseline_CN_longitudinal_trajectory"
)

prefix <- "BLSA_baseline_CN_harmonized_AD_EPOCH"

bin_width <- 0.5
minimum_scans_per_subject <- 2

# Small tolerance for scans occurring very close to the selected CN baseline.
baseline_time_tolerance_years <- 0.05

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

message("iSTAGING input: ", istaging_file)
message("Harmonized predictions: ", harmonized_prediction_file)
message("Output directory: ", out_dir)

# ------------------------------------------------------------------------------
# 2. Packages
# ------------------------------------------------------------------------------

required_packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "ggplot2",
  "stringr",
  "scales",
  "purrr",
  "tibble"
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
  library(purrr)
  library(tibble)
})

# ------------------------------------------------------------------------------
# 3. Visual settings
# ------------------------------------------------------------------------------

diagnosis_palette <- c(
  "CN"  = "#355C9A",
  "MCI" = "#E3A018",
  "AD"  = "#B55239"
)

eventual_group_palette <- c(
  "Remained CN" = "#355C9A",
  "Later MCI"   = "#E3A018",
  "Later AD"    = "#B55239"
)

baseline_fill <- "#F7E8A4"

# ------------------------------------------------------------------------------
# 4. Helper functions
# ------------------------------------------------------------------------------

detect_column <- function(
    df,
    preferred,
    regex,
    label,
    required = TRUE
) {
  direct <- preferred[
    preferred %in% names(df)
  ]
  
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

normalize_study <- function(x) {
  x <- toupper(clean_optional_character(x))
  
  case_when(
    str_detect(x, "^BLSA") ~ "BLSA",
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

format_p <- function(p) {
  case_when(
    is.na(p) ~ "P = NA",
    p < 2.2e-16 ~ "P < 2.2e-16",
    p < 0.001 ~ paste0(
      "P = ",
      formatC(
        p,
        format = "e",
        digits = 2
      )
    ),
    TRUE ~ paste0(
      "P = ",
      signif(p, 3)
    )
  )
}

format_beta <- function(x) {
  case_when(
    is.na(x) ~ "NA",
    abs(x) < 0.001 ~ formatC(
      x,
      format = "e",
      digits = 2
    ),
    TRUE ~ sprintf(
      "%.3f",
      x
    )
  )
}

select_earliest_mapped_diagnosis <- function(df) {
  df |>
    mutate(
      date_missing = is.na(diagnosis_date),
      
      baseline_selection_source = case_when(
        !date_missing ~
          "Earliest mapped diagnosis Date",
        
        date_missing &
          !is.na(diagnosis_age) ~
          "Youngest mapped diagnosis Age",
        
        TRUE ~
          "Mapped diagnosis without Date or Age"
      )
    ) |>
    arrange(
      participant_id,
      date_missing,
      diagnosis_date,
      diagnosis_age
    ) |>
    group_by(
      participant_id
    ) |>
    slice_head(
      n = 1
    ) |>
    ungroup()
}

fit_population_trend <- function(d) {
  d <- d |>
    filter(
      is.finite(
        years_since_cn_baseline
      ),
      is.finite(
        acceleration_years
      )
    )
  
  n_subjects <- n_distinct(
    d$participant_id
  )
  
  n_scans <- nrow(d)
  
  if (
    n_subjects < 2 ||
    n_scans < 4 ||
    n_distinct(
      d$years_since_cn_baseline
    ) < 2
  ) {
    return(
      tibble(
        n_subjects = n_subjects,
        n_scans = n_scans,
        model = "insufficient data",
        beta_slope_per_year = NA_real_,
        se = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        p_value = NA_real_,
        std_beta = NA_real_
      )
    )
  }
  
  sd_x <- sd(
    d$years_since_cn_baseline,
    na.rm = TRUE
  )
  
  sd_y <- sd(
    d$acceleration_years,
    na.rm = TRUE
  )
  
  if (
    requireNamespace(
      "nlme",
      quietly = TRUE
    )
  ) {
    fit_lme <- tryCatch(
      nlme::lme(
        acceleration_years ~
          years_since_cn_baseline,
        random = ~ 1 |
          participant_id,
        data = d,
        method = "REML",
        na.action = na.omit,
        control = nlme::lmeControl(
          opt = "optim",
          msMaxIter = 200,
          returnObject = TRUE
        )
      ),
      error = function(e) NULL
    )
    
    if (!is.null(fit_lme)) {
      tt <- summary(
        fit_lme
      )$tTable
      
      row_name <- "years_since_cn_baseline"
      
      if (row_name %in% rownames(tt)) {
        beta <- unname(
          tt[
            row_name,
            "Value"
          ]
        )
        
        se <- unname(
          tt[
            row_name,
            "Std.Error"
          ]
        )
        
        p_value <- unname(
          tt[
            row_name,
            "p-value"
          ]
        )
        
        model_df <- unname(
          tt[
            row_name,
            "DF"
          ]
        )
        
        ci_low <- beta -
          qt(
            0.975,
            df = model_df
          ) * se
        
        ci_high <- beta +
          qt(
            0.975,
            df = model_df
          ) * se
        
        std_beta <- ifelse(
          is.finite(sd_x) &&
            is.finite(sd_y) &&
            sd_y > 0,
          beta * sd_x / sd_y,
          NA_real_
        )
        
        return(
          tibble(
            n_subjects = n_subjects,
            n_scans = n_scans,
            model = paste0(
              "linear mixed model: ",
              "random intercept for participant"
            ),
            beta_slope_per_year = beta,
            se = se,
            ci_low = ci_low,
            ci_high = ci_high,
            p_value = p_value,
            std_beta = std_beta
          )
        )
      }
    }
  }
  
  fit_lm <- lm(
    acceleration_years ~
      years_since_cn_baseline,
    data = d
  )
  
  tt <- summary(
    fit_lm
  )$coefficients
  
  row_name <- "years_since_cn_baseline"
  
  beta <- unname(
    tt[
      row_name,
      "Estimate"
    ]
  )
  
  se <- unname(
    tt[
      row_name,
      "Std. Error"
    ]
  )
  
  p_value <- unname(
    tt[
      row_name,
      "Pr(>|t|)"
    ]
  )
  
  model_df <- fit_lm$df.residual
  
  ci_low <- beta -
    qt(
      0.975,
      df = model_df
    ) * se
  
  ci_high <- beta +
    qt(
      0.975,
      df = model_df
    ) * se
  
  std_beta <- ifelse(
    is.finite(sd_x) &&
      is.finite(sd_y) &&
      sd_y > 0,
    beta * sd_x / sd_y,
    NA_real_
  )
  
  tibble(
    n_subjects = n_subjects,
    n_scans = n_scans,
    model = "ordinary linear model fallback",
    beta_slope_per_year = beta,
    se = se,
    ci_low = ci_low,
    ci_high = ci_high,
    p_value = p_value,
    std_beta = std_beta
  )
}

# ------------------------------------------------------------------------------
# 5. Validate inputs
# ------------------------------------------------------------------------------

if (!file.exists(istaging_file)) {
  stop(
    "iSTAGING file does not exist: ",
    istaging_file
  )
}

if (!file.exists(
  harmonized_prediction_file
)) {
  stop(
    "Harmonized prediction file does not exist: ",
    harmonized_prediction_file
  )
}

# ------------------------------------------------------------------------------
# 6. Read full BLSA longitudinal diagnosis records
# ------------------------------------------------------------------------------

istaging_df <- readr::read_tsv(
  istaging_file,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "unique"
)

study_col <- detect_column(
  istaging_df,
  preferred = c(
    "Study",
    "STUDY"
  ),
  regex = "(^|_)study$",
  label = "iSTAGING Study"
)

id_col <- detect_column(
  istaging_df,
  preferred = c(
    "PTID",
    "participant_id",
    "IID"
  ),
  regex = "(^ptid$|participant.*id|^iid$)",
  label = "iSTAGING participant ID"
)

dx_col <- detect_column(
  istaging_df,
  preferred = c(
    "DX_Binary",
    "Dx_binary",
    "dx_binary"
  ),
  regex = "(^|_)dx[_\\.]*binary$",
  label = "iSTAGING DX_Binary"
)

date_col <- detect_column(
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

age_col <- detect_column(
  istaging_df,
  preferred = c(
    "Age",
    "AGE"
  ),
  regex = "^age$",
  label = "iSTAGING Age"
)

blsa_diagnosis_long <- istaging_df |>
  transmute(
    participant_id = as.character(
      .data[[id_col]]
    ),
    
    study = normalize_study(
      .data[[study_col]]
    ),
    
    diagnosis_original = clean_optional_character(
      .data[[dx_col]]
    ),
    
    diagnosis = normalize_diagnosis(
      .data[[dx_col]]
    ),
    
    diagnosis_date = if (!is.na(
      date_col
    )) {
      parse_date_flexibly(
        .data[[date_col]]
      )
    } else {
      as.Date(NA)
    },
    
    diagnosis_age = suppressWarnings(
      as.numeric(
        .data[[age_col]]
      )
    )
  ) |>
  filter(
    study == "BLSA",
    !is.na(participant_id),
    participant_id != "",
    !is.na(diagnosis)
  ) |>
  distinct(
    participant_id,
    diagnosis,
    diagnosis_date,
    diagnosis_age,
    .keep_all = TRUE
  )

# ------------------------------------------------------------------------------
# 7. Define baseline-CN participants
# ------------------------------------------------------------------------------

baseline_diagnosis <- blsa_diagnosis_long |>
  select_earliest_mapped_diagnosis()

baseline_cn_ids <- baseline_diagnosis |>
  filter(
    diagnosis == "CN"
  ) |>
  select(
    participant_id
  )

cn_baseline <- baseline_diagnosis |>
  filter(
    diagnosis == "CN"
  ) |>
  transmute(
    participant_id,
    baseline_diagnosis = diagnosis,
    diagnosis_date,
    diagnosis_age,
    baseline_selection_source
  )

# Eventual diagnosis grouping for subject-level summaries.
eventual_group <- blsa_diagnosis_long |>
  semi_join(
    baseline_cn_ids,
    by = "participant_id"
  ) |>
  group_by(
    participant_id
  ) |>
  summarise(
    ever_mci = any(
      diagnosis == "MCI",
      na.rm = TRUE
    ),
    
    ever_ad = any(
      diagnosis == "AD",
      na.rm = TRUE
    ),
    
    eventual_diagnosis_group = case_when(
      ever_ad ~ "Later AD",
      ever_mci ~ "Later MCI",
      TRUE ~ "Remained CN"
    ),
    
    .groups = "drop"
  ) |>
  mutate(
    eventual_diagnosis_group = factor(
      eventual_diagnosis_group,
      levels = c(
        "Remained CN",
        "Later MCI",
        "Later AD"
      )
    )
  )

message(
  "BLSA participants with mapped diagnosis: ",
  n_distinct(
    blsa_diagnosis_long$participant_id
  )
)

message(
  "Participants whose earliest mapped diagnosis is CN: ",
  nrow(
    baseline_cn_ids
  )
)

# ------------------------------------------------------------------------------
# 8. Read all scored harmonized BLSA EPOCH scans
# ------------------------------------------------------------------------------

prediction_df <- readr::read_tsv(
  harmonized_prediction_file,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "unique"
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

acceleration_years_col <- detect_column(
  prediction_df,
  preferred = c(
    "adni_brain_mri_ad_epoch_acceleration_years",
    "adni_brain_mri_ad_lepoch_acceleration_years"
  ),
  regex = "acceleration[_\\.]*years$",
  label = "AD EPOCH acceleration-years"
)

acceleration_z_col <- detect_column(
  prediction_df,
  preferred = c(
    "adni_brain_mri_ad_epoch_acceleration_z",
    "adni_brain_mri_ad_lepoch_acceleration_z"
  ),
  regex = "acceleration[_\\.]*z$",
  label = "AD EPOCH acceleration-z",
  required = FALSE
)

risk_score_col <- detect_column(
  prediction_df,
  preferred = c(
    "adni_brain_mri_ad_epoch_risk_score",
    "adni_brain_mri_ad_lepoch_risk_score"
  ),
  regex = "risk[_\\.]*score$",
  label = "AD EPOCH risk score",
  required = FALSE
)

blsa_epoch_scans <- prediction_df |>
  transmute(
    prediction_source_row = row_number(),
    
    participant_id = as.character(
      .data[[prediction_id_col]]
    ),
    
    study = normalize_study(
      .data[[prediction_study_col]]
    ),
    
    site = if (!is.na(
      prediction_site_col
    )) {
      clean_optional_character(
        .data[[prediction_site_col]]
      )
    } else {
      NA_character_
    },
    
    scan_date = if (!is.na(
      prediction_date_col
    )) {
      parse_date_flexibly(
        .data[[prediction_date_col]]
      )
    } else {
      as.Date(NA)
    },
    
    scan_age = suppressWarnings(
      as.numeric(
        .data[[prediction_age_col]]
      )
    ),
    
    acceleration_years = suppressWarnings(
      as.numeric(
        .data[[acceleration_years_col]]
      )
    ),
    
    acceleration_z = if (!is.na(
      acceleration_z_col
    )) {
      suppressWarnings(
        as.numeric(
          .data[[acceleration_z_col]]
        )
      )
    } else {
      NA_real_
    },
    
    risk_score = if (!is.na(
      risk_score_col
    )) {
      suppressWarnings(
        as.numeric(
          .data[[risk_score_col]]
        )
      )
    } else {
      NA_real_
    }
  ) |>
  filter(
    study == "BLSA",
    !is.na(participant_id),
    participant_id != "",
    !is.na(acceleration_years)
  )

# ------------------------------------------------------------------------------
# 9. Align scored scans to baseline and assign nearest diagnosis
# ------------------------------------------------------------------------------

baseline_aligned_scans <- blsa_epoch_scans |>
  semi_join(
    baseline_cn_ids,
    by = "participant_id"
  ) |>
  inner_join(
    cn_baseline,
    by = "participant_id",
    relationship = "many-to-one"
  ) |>
  mutate(
    years_since_cn_baseline_by_date = ifelse(
      !is.na(scan_date) &
        !is.na(diagnosis_date),
      as.numeric(
        scan_date -
          diagnosis_date
      ) / 365.25,
      NA_real_
    ),
    
    years_since_cn_baseline_by_age = ifelse(
      !is.na(scan_age) &
        !is.na(diagnosis_age),
      scan_age -
        diagnosis_age,
      NA_real_
    ),
    
    time_source = case_when(
      is.finite(
        years_since_cn_baseline_by_date
      ) ~ "Date",
      
      is.finite(
        years_since_cn_baseline_by_age
      ) ~ "Age",
      
      TRUE ~ "Unavailable"
    ),
    
    years_since_cn_baseline = coalesce(
      years_since_cn_baseline_by_date,
      years_since_cn_baseline_by_age
    )
  ) |>
  filter(
    is.finite(
      years_since_cn_baseline
    ),
    years_since_cn_baseline >=
      -baseline_time_tolerance_years
  ) |>
  mutate(
    years_since_cn_baseline = ifelse(
      abs(
        years_since_cn_baseline
      ) <=
        baseline_time_tolerance_years,
      0,
      years_since_cn_baseline
    )
  )

# Candidate scan-diagnosis pairs.
scan_diagnosis_candidates <- baseline_aligned_scans |>
  select(
    prediction_source_row,
    participant_id,
    scan_date,
    scan_age
  ) |>
  inner_join(
    blsa_diagnosis_long |>
      semi_join(
        baseline_cn_ids,
        by = "participant_id"
      ),
    by = "participant_id",
    relationship = "many-to-many"
  ) |>
  mutate(
    diagnosis_date_difference_days = ifelse(
      !is.na(scan_date) &
        !is.na(diagnosis_date),
      abs(
        as.numeric(
          scan_date -
            diagnosis_date
        )
      ),
      NA_real_
    ),
    
    diagnosis_age_difference_years = ifelse(
      !is.na(scan_age) &
        !is.na(diagnosis_age),
      abs(
        scan_age -
          diagnosis_age
      ),
      NA_real_
    ),
    
    diagnosis_match_priority = case_when(
      is.finite(
        diagnosis_date_difference_days
      ) ~ 1,
      
      is.finite(
        diagnosis_age_difference_years
      ) ~ 2,
      
      !is.na(
        diagnosis_date
      ) ~ 3,
      
      !is.na(
        diagnosis_age
      ) ~ 4,
      
      TRUE ~ 5
    ),
    
    diagnosis_match_distance = case_when(
      diagnosis_match_priority == 1 ~
        diagnosis_date_difference_days,
      
      diagnosis_match_priority == 2 ~
        diagnosis_age_difference_years,
      
      diagnosis_match_priority == 3 ~
        as.numeric(
          diagnosis_date
        ),
      
      diagnosis_match_priority == 4 ~
        diagnosis_age,
      
      TRUE ~ Inf
    ),
    
    diagnosis_match_source = case_when(
      diagnosis_match_priority == 1 ~
        "Closest diagnosis by Date",
      
      diagnosis_match_priority == 2 ~
        "Closest diagnosis by Age",
      
      diagnosis_match_priority == 3 ~
        "Earliest diagnosis Date",
      
      diagnosis_match_priority == 4 ~
        "Youngest diagnosis Age",
      
      TRUE ~
        "Unavailable"
    )
  ) |>
  arrange(
    prediction_source_row,
    diagnosis_match_priority,
    diagnosis_match_distance,
    diagnosis_date,
    diagnosis_age
  ) |>
  group_by(
    prediction_source_row
  ) |>
  slice_head(
    n = 1
  ) |>
  ungroup() |>
  transmute(
    prediction_source_row,
    scan_diagnosis = diagnosis,
    matched_diagnosis_date = diagnosis_date,
    matched_diagnosis_age = diagnosis_age,
    diagnosis_match_source,
    diagnosis_date_difference_days,
    diagnosis_age_difference_years
  )

aligned_scans <- baseline_aligned_scans |>
  left_join(
    scan_diagnosis_candidates,
    by = "prediction_source_row",
    relationship = "one-to-one"
  ) |>
  left_join(
    eventual_group,
    by = "participant_id",
    relationship = "many-to-one"
  ) |>
  mutate(
    scan_diagnosis = factor(
      scan_diagnosis,
      levels = c(
        "CN",
        "MCI",
        "AD"
      )
    )
  ) |>
  arrange(
    participant_id,
    years_since_cn_baseline,
    scan_date,
    scan_age
  )

# Remove duplicated participant-time records.
aligned_scans <- aligned_scans |>
  group_by(
    participant_id,
    years_since_cn_baseline
  ) |>
  slice_head(
    n = 1
  ) |>
  ungroup()

# Keep participants with at least two qualified scans and nonzero follow-up.
longitudinal_subjects <- aligned_scans |>
  group_by(
    participant_id
  ) |>
  summarise(
    eventual_diagnosis_group = first(
      eventual_diagnosis_group
    ),
    
    n_scans = n(),
    
    followup_span_years = max(
      years_since_cn_baseline,
      na.rm = TRUE
    ) -
      min(
        years_since_cn_baseline,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) |>
  filter(
    n_scans >=
      minimum_scans_per_subject,
    followup_span_years > 0
  )

traj_df <- aligned_scans |>
  semi_join(
    longitudinal_subjects,
    by = "participant_id"
  ) |>
  group_by(
    participant_id
  ) |>
  arrange(
    years_since_cn_baseline,
    .by_group = TRUE
  ) |>
  mutate(
    scan_relation = ifelse(
      row_number() == 1,
      "selected_EPOCH_baseline",
      "followup"
    )
  ) |>
  ungroup() |>
  mutate(
    scan_relation = factor(
      scan_relation,
      levels = c(
        "selected_EPOCH_baseline",
        "followup"
      )
    )
  )

if (nrow(traj_df) == 0) {
  stop(
    "No baseline-CN BLSA participants had at least ",
    minimum_scans_per_subject,
    " qualified longitudinal EPOCH scans."
  )
}

# ------------------------------------------------------------------------------
# 10. QC and summaries
# ------------------------------------------------------------------------------

cohort_flow <- tibble(
  stage = c(
    "BLSA participants with mapped diagnosis",
    "Participants whose earliest mapped diagnosis is CN",
    "Baseline-CN participants with >=1 scored EPOCH scan",
    paste0(
      "Baseline-CN participants with >=",
      minimum_scans_per_subject,
      " longitudinal EPOCH scans"
    ),
    "Participants remaining CN",
    "Participants later diagnosed with MCI",
    "Participants later diagnosed with AD"
  ),
  
  n_participants = c(
    n_distinct(
      blsa_diagnosis_long$participant_id
    ),
    
    nrow(
      baseline_cn_ids
    ),
    
    n_distinct(
      aligned_scans$participant_id
    ),
    
    n_distinct(
      traj_df$participant_id
    ),
    
    sum(
      longitudinal_subjects$eventual_diagnosis_group ==
        "Remained CN"
    ),
    
    sum(
      longitudinal_subjects$eventual_diagnosis_group ==
        "Later MCI"
    ),
    
    sum(
      longitudinal_subjects$eventual_diagnosis_group ==
        "Later AD"
    )
  )
)

scan_summary <- traj_df |>
  group_by(
    scan_diagnosis
  ) |>
  summarise(
    n_scans = n(),
    
    n_subjects = n_distinct(
      participant_id
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
    
    .groups = "drop"
  )

subject_summary <- traj_df |>
  group_by(
    participant_id
  ) |>
  arrange(
    years_since_cn_baseline,
    .by_group = TRUE
  ) |>
  summarise(
    eventual_diagnosis_group = first(
      eventual_diagnosis_group
    ),
    
    n_scans = n(),
    
    n_cn_scans = sum(
      scan_diagnosis == "CN",
      na.rm = TRUE
    ),
    
    n_mci_scans = sum(
      scan_diagnosis == "MCI",
      na.rm = TRUE
    ),
    
    n_ad_scans = sum(
      scan_diagnosis == "AD",
      na.rm = TRUE
    ),
    
    followup_span_years = max(
      years_since_cn_baseline
    ) -
      min(
        years_since_cn_baseline
      ),
    
    baseline_acceleration_years = first(
      acceleration_years
    ),
    
    last_acceleration_years = last(
      acceleration_years
    ),
    
    change_last_minus_baseline = (
      last_acceleration_years -
        baseline_acceleration_years
    ),
    
    .groups = "drop"
  )

diagnosis_match_summary <- traj_df |>
  count(
    diagnosis_match_source,
    scan_diagnosis,
    name = "n_scans"
  )

# ------------------------------------------------------------------------------
# 11. Separate longitudinal trend models by eventual diagnosis group
# ------------------------------------------------------------------------------

# Fit one longitudinal model separately for:
#   1. participants who remained CN;
#   2. participants who later developed MCI;
#   3. participants who later developed AD.
#
# The group is participant-level and is based on all available mapped
# longitudinal diagnoses after the initial CN diagnosis.

group_trend_tbl <- traj_df |>
  filter(
    !is.na(eventual_diagnosis_group)
  ) |>
  group_by(
    eventual_diagnosis_group
  ) |>
  group_modify(
    ~ fit_population_trend(.x)
  ) |>
  ungroup() |>
  mutate(
    eventual_diagnosis_group = factor(
      eventual_diagnosis_group,
      levels = c(
        "Remained CN",
        "Later MCI",
        "Later AD"
      )
    ),
    
    trend_label = paste0(
      as.character(
        eventual_diagnosis_group
      ),
      "\nN = ",
      n_subjects,
      "; scans = ",
      n_scans,
      "\n\u03b2 = ",
      format_beta(
        beta_slope_per_year
      ),
      " y/y",
      "\n95% CI: ",
      format_beta(
        ci_low
      ),
      " to ",
      format_beta(
        ci_high
      ),
      "\nstd \u03b2 = ",
      format_beta(
        std_beta
      ),
      "; ",
      format_p(
        p_value
      )
    )
  )

# Also retain the pooled model as a secondary output for completeness.
overall_trend_tbl <- fit_population_trend(
  traj_df
) |>
  mutate(
    analysis_group = "All baseline-CN participants",
    .before = 1
  )

# Create nonoverlapping annotation positions for the three group-specific fits.
plot_ranges <- traj_df |>
  summarise(
    x_min = min(
      years_since_cn_baseline,
      na.rm = TRUE
    ),
    
    x_max = max(
      years_since_cn_baseline,
      na.rm = TRUE
    ),
    
    y_min = min(
      acceleration_years,
      na.rm = TRUE
    ),
    
    y_max = max(
      acceleration_years,
      na.rm = TRUE
    )
  ) |>
  mutate(
    x_range = pmax(
      x_max - x_min,
      1
    ),
    
    y_range = pmax(
      y_max - y_min,
      1
    )
  )

group_annotation_positions <- group_trend_tbl |>
  arrange(
    eventual_diagnosis_group
  ) |>
  mutate(
    annotation_index = row_number()
  ) |>
  mutate(
    x = plot_ranges$x_min +
      0.03 *
      plot_ranges$x_range,
    
    y = plot_ranges$y_max -
      (
        0.03 +
          0.20 *
          (
            annotation_index -
              1
          )
      ) *
      plot_ranges$y_range
  )

# ------------------------------------------------------------------------------
# 12. Plot 1: individual trajectories with separate fitted lines by outcome
# ------------------------------------------------------------------------------

relation_shapes <- c(
  "selected_EPOCH_baseline" = 21,
  "followup" = 24
)

p_spaghetti <- ggplot(
  traj_df,
  aes(
    x = years_since_cn_baseline,
    y = acceleration_years,
    group = participant_id
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey55"
  ) +
  
  # Individual longitudinal participant trajectories.
  geom_line(
    color = "grey60",
    alpha = 0.24,
    linewidth = 0.48
  ) +
  
  # Session-level diagnosis is shown by point fill.
  geom_point(
    aes(
      fill = scan_diagnosis,
      shape = scan_relation
    ),
    color = "grey20",
    alpha = 0.82,
    size = 2.05,
    stroke = 0.38
  ) +
  
  # Separate linear fitting lines for the three participant-level outcomes.
  geom_smooth(
    aes(
      color = eventual_diagnosis_group,
      group = eventual_diagnosis_group
    ),
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    linewidth = 1.35,
    alpha = 1
  ) +
  
  # Group-specific slope/model annotations.
  geom_label(
    data = group_annotation_positions,
    aes(
      x = x,
      y = y,
      label = trend_label,
      color = eventual_diagnosis_group
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 2.85,
    lineheight = 0.98,
    label.size = 0.25,
    label.padding = grid::unit(
      0.16,
      "lines"
    ),
    fill = "white",
    alpha = 0.92,
    show.legend = FALSE
  ) +
  
  scale_fill_manual(
    values = diagnosis_palette,
    limits = c(
      "CN",
      "MCI",
      "AD"
    ),
    drop = FALSE,
    name = "Diagnosis at MRI session"
  ) +
  
  scale_color_manual(
    values = eventual_group_palette,
    limits = c(
      "Remained CN",
      "Later MCI",
      "Later AD"
    ),
    drop = FALSE,
    name = "Participant outcome"
  ) +
  
  scale_shape_manual(
    values = relation_shapes,
    drop = FALSE,
    name = "Scan relation",
    labels = c(
      "Selected EPOCH baseline",
      "Follow-up"
    )
  ) +
  
  scale_x_continuous(
    name = "Years since selected CN diagnosis baseline",
    breaks = pretty_breaks(
      n = 6
    )
  ) +
  
  scale_y_continuous(
    name = "Harmonized AD EPOCH acceleration (years)",
    labels = number_format(
      accuracy = 0.1
    ),
    expand = expansion(
      mult = c(
        0.08,
        0.25
      )
    )
  ) +
  
  labs(
    title = paste0(
      "Longitudinal AD EPOCH trajectories in ",
      "baseline-CN BLSA participants"
    ),
    
    subtitle = paste0(
      "Separate linear fits are shown for participants who remained CN, ",
      "later developed MCI, or later developed AD; ",
      "participants are required to have at least ",
      minimum_scans_per_subject,
      " qualified scans"
    ),
    
    caption = paste0(
      "Grey lines represent individual participants. Point fill indicates ",
      "the diagnosis assigned to each MRI session. Colored fitted lines and ",
      "annotations correspond to participant-level eventual diagnosis. ",
      "\u03b2 is the estimated change in AD EPOCH acceleration years per year."
    )
  ) +
  
  coord_cartesian(
    clip = "off"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 13
    ),
    
    plot.subtitle = element_text(
      size = 10.3
    ),
    
    plot.caption = element_text(
      size = 8.5,
      hjust = 0
    ),
    
    legend.position = "top",
    legend.box = "vertical",
    legend.justification = "left",
    
    axis.text = element_text(
      color = "black"
    ),
    
    axis.title = element_text(
      face = "bold"
    ),
    
    panel.grid.major.y = element_line(
      color = "grey90",
      linewidth = 0.35
    ),
    
    panel.grid.minor = element_blank(),
    
    plot.margin = margin(
      8,
      12,
      8,
      8
    )
  )

# ------------------------------------------------------------------------------
# 13. Plot 2: diagnosis-stratified population trajectories
# ------------------------------------------------------------------------------

population_by_dx <- traj_df |>
  mutate(
    time_bin = round(
      years_since_cn_baseline /
        bin_width
    ) *
      bin_width
  ) |>
  group_by(
    scan_diagnosis,
    time_bin
  ) |>
  summarise(
    n_scans = n(),
    
    n_subjects = n_distinct(
      participant_id
    ),
    
    mean_acceleration_years = mean(
      acceleration_years,
      na.rm = TRUE
    ),
    
    sd_acceleration_years = sd(
      acceleration_years,
      na.rm = TRUE
    ),
    
    se_acceleration_years = (
      sd_acceleration_years /
        sqrt(
          n_scans
        )
    ),
    
    ci_low = (
      mean_acceleration_years -
        1.96 *
        se_acceleration_years
    ),
    
    ci_high = (
      mean_acceleration_years +
        1.96 *
        se_acceleration_years
    ),
    
    .groups = "drop"
  ) |>
  filter(
    n_scans >= 2
  )

p_population_dx <- ggplot(
  population_by_dx,
  aes(
    x = time_bin,
    y = mean_acceleration_years,
    color = scan_diagnosis,
    fill = scan_diagnosis
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey55"
  ) +
  geom_ribbon(
    aes(
      ymin = ci_low,
      ymax = ci_high
    ),
    alpha = 0.15,
    color = NA
  ) +
  geom_line(
    linewidth = 1.05
  ) +
  geom_point(
    aes(
      size = n_subjects
    ),
    alpha = 0.90
  ) +
  scale_color_manual(
    values = diagnosis_palette,
    drop = FALSE
  ) +
  scale_fill_manual(
    values = diagnosis_palette,
    drop = FALSE
  ) +
  scale_size_continuous(
    name = "N subjects",
    range = c(
      1.8,
      4.4
    )
  ) +
  scale_x_continuous(
    name = "Years since selected CN diagnosis baseline",
    breaks = pretty_breaks(
      n = 6
    )
  ) +
  scale_y_continuous(
    name = "Mean harmonized AD EPOCH acceleration (years)",
    labels = number_format(
      accuracy = 0.1
    )
  ) +
  labs(
    title = "Population trajectories by scan diagnosis",
    
    subtitle = paste0(
      "Mean AD EPOCH acceleration in ",
      bin_width,
      "-year bins"
    ),
    
    color = "Scan diagnosis",
    fill = "Scan diagnosis",
    
    caption = paste0(
      "Shaded regions are approximate 95% confidence intervals. ",
      "Bins with fewer than two scans are omitted."
    )
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 13
    ),
    
    plot.subtitle = element_text(
      size = 10.5
    ),
    
    plot.caption = element_text(
      size = 8.5,
      hjust = 0
    ),
    
    legend.position = "top",
    
    axis.text = element_text(
      color = "black"
    ),
    
    axis.title = element_text(
      face = "bold"
    ),
    
    panel.grid.major.y = element_line(
      color = "grey90",
      linewidth = 0.35
    ),
    
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------------------------------------
# 14. Plot 3: within-subject delta from first qualified EPOCH scan
# ------------------------------------------------------------------------------

baseline_epoch_tbl <- traj_df |>
  group_by(
    participant_id
  ) |>
  arrange(
    years_since_cn_baseline,
    .by_group = TRUE
  ) |>
  summarise(
    selected_epoch_baseline_year = first(
      years_since_cn_baseline
    ),
    
    baseline_acceleration_years = first(
      acceleration_years
    ),
    
    .groups = "drop"
  )

delta_df <- traj_df |>
  left_join(
    baseline_epoch_tbl,
    by = "participant_id"
  ) |>
  mutate(
    years_since_selected_epoch_baseline = (
      years_since_cn_baseline -
        selected_epoch_baseline_year
    ),
    
    delta_acceleration_years = (
      acceleration_years -
        baseline_acceleration_years
    )
  )

delta_population_df <- delta_df |>
  mutate(
    time_bin = round(
      years_since_selected_epoch_baseline /
        bin_width
    ) *
      bin_width
  ) |>
  group_by(
    scan_diagnosis,
    time_bin
  ) |>
  summarise(
    n_scans = n(),
    
    n_subjects = n_distinct(
      participant_id
    ),
    
    mean_delta_acceleration_years = mean(
      delta_acceleration_years,
      na.rm = TRUE
    ),
    
    sd_delta_acceleration_years = sd(
      delta_acceleration_years,
      na.rm = TRUE
    ),
    
    se_delta_acceleration_years = (
      sd_delta_acceleration_years /
        sqrt(
          n_scans
        )
    ),
    
    ci_low = (
      mean_delta_acceleration_years -
        1.96 *
        se_delta_acceleration_years
    ),
    
    ci_high = (
      mean_delta_acceleration_years +
        1.96 *
        se_delta_acceleration_years
    ),
    
    .groups = "drop"
  ) |>
  filter(
    n_scans >= 2
  )

p_delta <- ggplot(
  delta_population_df,
  aes(
    x = time_bin,
    y = mean_delta_acceleration_years,
    color = scan_diagnosis,
    fill = scan_diagnosis
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey55"
  ) +
  geom_ribbon(
    aes(
      ymin = ci_low,
      ymax = ci_high
    ),
    alpha = 0.15,
    color = NA
  ) +
  geom_line(
    linewidth = 1.05
  ) +
  geom_point(
    aes(
      size = n_subjects
    ),
    alpha = 0.90
  ) +
  scale_color_manual(
    values = diagnosis_palette,
    drop = FALSE
  ) +
  scale_fill_manual(
    values = diagnosis_palette,
    drop = FALSE
  ) +
  scale_size_continuous(
    name = "N subjects",
    range = c(
      1.8,
      4.4
    )
  ) +
  scale_x_continuous(
    name = "Years since first qualified EPOCH scan",
    breaks = pretty_breaks(
      n = 6
    )
  ) +
  scale_y_continuous(
    name = "Mean change in AD EPOCH acceleration (years)",
    labels = number_format(
      accuracy = 0.1
    )
  ) +
  labs(
    title = "Within-subject change by scan diagnosis",
    
    subtitle = paste0(
      "Each participant is centered at the first qualified ",
      "harmonized EPOCH scan"
    ),
    
    color = "Scan diagnosis",
    fill = "Scan diagnosis",
    
    caption = paste0(
      "Positive values indicate increasing AD EPOCH acceleration ",
      "relative to the participant's selected EPOCH baseline."
    )
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 13
    ),
    
    plot.subtitle = element_text(
      size = 10.5
    ),
    
    plot.caption = element_text(
      size = 8.5,
      hjust = 0
    ),
    
    legend.position = "top",
    
    axis.text = element_text(
      color = "black"
    ),
    
    axis.title = element_text(
      face = "bold"
    ),
    
    panel.grid.major.y = element_line(
      color = "grey90",
      linewidth = 0.35
    ),
    
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------------------------------------
# 15. Plot 4: subject-level slopes by eventual diagnosis group
# ------------------------------------------------------------------------------

slope_tbl <- traj_df |>
  group_by(
    participant_id
  ) |>
  filter(
    n() >=
      minimum_scans_per_subject
  ) |>
  summarise(
    eventual_diagnosis_group = first(
      eventual_diagnosis_group
    ),
    
    n_scans = n(),
    
    followup_span_years = max(
      years_since_cn_baseline
    ) -
      min(
        years_since_cn_baseline
      ),
    
    slope_acceleration_years_per_year = {
      fit <- lm(
        acceleration_years ~
          years_since_cn_baseline
      )
      
      unname(
        coef(fit)[["years_since_cn_baseline"]]
      )
    },
    
    .groups = "drop"
  ) |>
  filter(
    is.finite(
      slope_acceleration_years_per_year
    )
  )

slope_summary <- slope_tbl |>
  group_by(
    eventual_diagnosis_group
  ) |>
  summarise(
    n_subjects = n(),
    
    mean_slope = mean(
      slope_acceleration_years_per_year,
      na.rm = TRUE
    ),
    
    sd_slope = sd(
      slope_acceleration_years_per_year,
      na.rm = TRUE
    ),
    
    median_slope = median(
      slope_acceleration_years_per_year,
      na.rm = TRUE
    ),
    
    q1 = quantile(
      slope_acceleration_years_per_year,
      0.25,
      na.rm = TRUE
    ),
    
    q3 = quantile(
      slope_acceleration_years_per_year,
      0.75,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

slope_labels <- slope_summary |>
  mutate(
    label = paste0(
      as.character(
        eventual_diagnosis_group
      ),
      "\n(n = ",
      n_subjects,
      ")"
    )
  ) |>
  select(
    eventual_diagnosis_group,
    label
  ) |>
  tibble::deframe()

p_slope <- ggplot(
  slope_tbl,
  aes(
    x = eventual_diagnosis_group,
    y = slope_acceleration_years_per_year,
    fill = eventual_diagnosis_group,
    color = eventual_diagnosis_group
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey55"
  ) +
  geom_violin(
    width = 0.80,
    alpha = 0.28,
    trim = FALSE,
    linewidth = 0.40
  ) +
  geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.90,
    linewidth = 0.48
  ) +
  geom_jitter(
    width = 0.10,
    alpha = 0.58,
    size = 1.7
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 23,
    size = 3.2,
    fill = baseline_fill,
    color = "black"
  ) +
  scale_fill_manual(
    values = eventual_group_palette,
    guide = "none"
  ) +
  scale_color_manual(
    values = eventual_group_palette,
    guide = "none"
  ) +
  scale_x_discrete(
    name = NULL,
    labels = slope_labels
  ) +
  scale_y_continuous(
    name = paste0(
      "Subject-level slope\n",
      "(acceleration years per year)"
    ),
    labels = number_format(
      accuracy = 0.1
    )
  ) +
  labs(
    title = "Subject-level AD EPOCH slopes by eventual diagnosis",
    
    subtitle = paste0(
      "All participants were CN at their earliest mapped diagnosis"
    ),
    
    caption = paste0(
      "Positive slopes indicate increasing AD EPOCH acceleration ",
      "during longitudinal follow-up."
    )
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 13
    ),
    
    plot.subtitle = element_text(
      size = 10.5
    ),
    
    plot.caption = element_text(
      size = 8.5,
      hjust = 0
    ),
    
    axis.text = element_text(
      color = "black"
    ),
    
    axis.title = element_text(
      face = "bold"
    ),
    
    panel.grid.major.y = element_line(
      color = "grey90",
      linewidth = 0.35
    ),
    
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------------------------------------
# 16. Save tables
# ------------------------------------------------------------------------------

readr::write_tsv(
  baseline_diagnosis,
  file.path(
    out_dir,
    paste0(
      prefix,
      "_earliest_mapped_diagnosis.tsv"
    )
  )
)

readr::write_tsv(
  cn_baseline,
  file.path(
    out_dir,
    paste0(
      prefix,
      "_selected_CN_baselines.tsv"
    )
  )
)

readr::write_tsv(
  eventual_group,
  file.path(
    out_dir,
    paste0(
      prefix,
      "_eventual_diagnosis_groups.tsv"
    )
  )
)

readr::write_tsv(
  aligned_scans,
  file.path(
    out_dir,
    paste0(
      prefix,
      "_all_aligned_scored_scans_with_diagnosis.tsv"
    )
  )
)

readr::write_tsv(
  traj_df,
  file.path(
    out_dir,
    paste0(
      prefix,
      "_longitudinal_scans_min2.tsv"
    )
  )
)

readr::write_tsv(
  cohort_flow,
  file.path(
    out_dir,
    paste0(
      prefix,
      "_cohort_flow.tsv"
    )
  )
)

readr::write_tsv(
  scan_summary,
  file.path(
    out_dir,
    paste0(
      prefix,
      "_scan_summary_by_diagnosis.tsv"
    )
  )
)

readr::write_tsv(
  subject_summary,
  file.path(
    out_dir,
    paste0(
      prefix,
      "_subject_summary.tsv"
    )
  )
)

readr::write_tsv(
  diagnosis_match_summary,
  file.path(
    out_dir,
    paste0(
      prefix,
      "_diagnosis_match_summary.tsv"
    )
  )
)

readr::write_tsv(
  group_trend_tbl,
  file.path(
    out_dir,
    paste0(
      prefix,
      "_trend_models_by_eventual_diagnosis.tsv"
    )
  )
)

readr::write_tsv(
  overall_trend_tbl,
  file.path(
    out_dir,
    paste0(
      prefix,
      "_overall_population_trend_model.tsv"
    )
  )
)

readr::write_tsv(
  population_by_dx,
  file.path(
    out_dir,
    paste0(
      prefix,
      "_population_timebin_by_diagnosis.tsv"
    )
  )
)

readr::write_tsv(
  delta_df,
  file.path(
    out_dir,
    paste0(
      prefix,
      "_delta_from_selected_epoch_baseline.tsv"
    )
  )
)

readr::write_tsv(
  delta_population_df,
  file.path(
    out_dir,
    paste0(
      prefix,
      "_population_delta_timebin_by_diagnosis.tsv"
    )
  )
)

readr::write_tsv(
  slope_tbl,
  file.path(
    out_dir,
    paste0(
      prefix,
      "_subject_level_slopes.tsv"
    )
  )
)

readr::write_tsv(
  slope_summary,
  file.path(
    out_dir,
    paste0(
      prefix,
      "_subject_level_slope_summary.tsv"
    )
  )
)

# ------------------------------------------------------------------------------
# 17. Save figures
# ------------------------------------------------------------------------------

ggsave(
  file.path(
    out_dir,
    paste0(
      prefix,
      "_individual_trajectories_diagnosis_colored.pdf"
    )
  ),
  p_spaghetti,
  width = 9.2,
  height = 6.0
)

ggsave(
  file.path(
    out_dir,
    paste0(
      prefix,
      "_individual_trajectories_diagnosis_colored.png"
    )
  ),
  p_spaghetti,
  width = 9.2,
  height = 6.0,
  dpi = 350,
  bg = "white"
)

ggsave(
  file.path(
    out_dir,
    paste0(
      prefix,
      "_population_trajectory_by_scan_diagnosis.pdf"
    )
  ),
  p_population_dx,
  width = 8.7,
  height = 5.7
)

ggsave(
  file.path(
    out_dir,
    paste0(
      prefix,
      "_population_trajectory_by_scan_diagnosis.png"
    )
  ),
  p_population_dx,
  width = 8.7,
  height = 5.7,
  dpi = 350,
  bg = "white"
)

ggsave(
  file.path(
    out_dir,
    paste0(
      prefix,
      "_population_delta_by_scan_diagnosis.pdf"
    )
  ),
  p_delta,
  width = 8.7,
  height = 5.7
)

ggsave(
  file.path(
    out_dir,
    paste0(
      prefix,
      "_population_delta_by_scan_diagnosis.png"
    )
  ),
  p_delta,
  width = 8.7,
  height = 5.7,
  dpi = 350,
  bg = "white"
)

ggsave(
  file.path(
    out_dir,
    paste0(
      prefix,
      "_subject_level_slope_by_eventual_diagnosis.pdf"
    )
  ),
  p_slope,
  width = 8.0,
  height = 5.7
)

ggsave(
  file.path(
    out_dir,
    paste0(
      prefix,
      "_subject_level_slope_by_eventual_diagnosis.png"
    )
  ),
  p_slope,
  width = 8.0,
  height = 5.7,
  dpi = 350,
  bg = "white"
)

# ------------------------------------------------------------------------------
# 18. Print summaries
# ------------------------------------------------------------------------------

message("============================================================")
message("Analysis complete.")
message("")
message("Cohort flow:")
print(cohort_flow)
message("")
message("Scan summary by assigned diagnosis:")
print(scan_summary)
message("")
message("Separate trend models by eventual diagnosis:")
print(group_trend_tbl)
message("")
message("Overall pooled trend model:")
print(overall_trend_tbl)
message("")
message("Subject-level slope summary:")
print(slope_summary)
message("")
message("Outputs saved to: ", out_dir)
message("============================================================")