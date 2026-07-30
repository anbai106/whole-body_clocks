#!/usr/bin/env Rscript

# ==============================================================================
# Longitudinal raw-MUSE AD EPOCH trajectories in A4
#
# Primary analyses
# ----------------
# 1. Compare longitudinal AD EPOCH trajectories between randomized treatment
#    arms (Placebo versus active drug).
# 2. Among participants who are amyloid-negative at their earliest available
#    amyloid assessment, compare participants who later convert to amyloid
#    positive with those who remain amyloid negative.
# 3. Produce an optional joint treatment-by-amyloid-trajectory figure when
#    sufficient participants are present in the corresponding cells.
#
# Data sources
# ------------
# - A4 AD EPOCH scan-level predictions generated from raw MUSE features.
# - Derived_Data/SUBJINFO.csv for treatment assignment and participant metadata.
# - External_Data/imaging_SUVR_amyloid.csv for longitudinal amyloid measurements.
#
# Important scientific note
# -------------------------
# A4 randomized participants were selected on the basis of elevated amyloid,
# whereas amyloid-negative participants may primarily arise from the LEARN
# cohort. Therefore, the treatment-arm and amyloid-conversion analyses are
# implemented as complementary analyses. The script reports the observed
# treatment-by-amyloid group counts and only fits a joint interaction model
# when the necessary combinations are represented.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. User settings
# ------------------------------------------------------------------------------

a4_prediction_file <- paste0(
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/adni_lepoch/",
  "results_a4_longitudinal_ad_epoch/",
  "a4_adni_brain_mri_ad_epoch_scan_level_predictions.tsv"
)

subjinfo_file <- paste0(
  "/Users/hao/cubic-projects/MULTI/download/A4/Clinical/Derived_Data/SUBJINFO.csv"
)

amyloid_file <- paste0(
  "/Users/hao/cubic-projects/MULTI/download/A4/Clinical/External_Data/imaging_SUVR_amyloid.csv"
)

out_dir <- paste0(
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/adni_lepoch/",
  "results_a4_longitudinal_ad_epoch/",
  "A4_treatment_amyloid_longitudinal_trajectory"
)

prefix <- "A4_raw_MUSE_AD_EPOCH"

minimum_scans_per_subject <- 2
minimum_followup_years <- 0
bin_width_years <- 0.5

# A4 amyloid-elevation rule used at screening:
#   - SUVR >= 1.15: elevated/positive
#   - SUVR 1.10 to <1.15: elevated only when the central visual read was positive
#   - SUVR < 1.10: not elevated
#
# Using the published florbetapir conversion:
#   Centiloid = 183 * SUVR - 177
# these correspond approximately to:
#   - automatic-positive threshold: 33.45 Centiloids
#   - lower borderline threshold:   24.30 Centiloids
#
# A4 therefore did not use a single Centiloid-only cutoff for all scans.
amyloid_positive_suvr_threshold <- 1.15
amyloid_borderline_suvr_lower <- 1.10
amyloid_positive_centiloid_threshold <- 33.45
amyloid_borderline_centiloid_lower <- 24.30

# Maximum temporal distance for assigning an amyloid assessment to an EPOCH MRI.
# This assignment is used for scan-level visualization, not for defining eventual
# participant-level conversion.
maximum_amyloid_match_distance_years <- 1.5

# If TRUE, the conversion-specific plots include only participants whose earliest
# amyloid assessment was negative.
restrict_conversion_analysis_to_baseline_negative <- TRUE

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("A4 EPOCH predictions: ", a4_prediction_file)
message("A4 participant information: ", subjinfo_file)
message("A4 amyloid imaging: ", amyloid_file)
message("Output directory: ", out_dir)

# ------------------------------------------------------------------------------
# 2. Packages
# ------------------------------------------------------------------------------

required_packages <- c(
  "readr", "dplyr", "tidyr", "ggplot2", "stringr",
  "scales", "purrr", "tibble", "forcats"
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
  library(stringr)
  library(scales)
  library(purrr)
  library(tibble)
  library(forcats)
})

# ------------------------------------------------------------------------------
# 3. Visual settings
# ------------------------------------------------------------------------------

treatment_palette <- c(
  "Placebo" = "#355C9A",
  "Drug" = "#B55239",
  "Not randomized/Unknown" = "#777777"
)

amyloid_palette <- c(
  "Remained amyloid negative" = "#355C9A",
  "Converted to amyloid positive" = "#D97706",
  "Amyloid positive at baseline" = "#B55239",
  "Amyloid status unavailable" = "#777777"
)

scan_amyloid_palette <- c(
  "Amyloid negative" = "#355C9A",
  "Amyloid positive" = "#B55239",
  "Amyloid status unavailable" = "#777777"
)

# ------------------------------------------------------------------------------
# 4. Helper functions
# ------------------------------------------------------------------------------

detect_column <- function(df, preferred, regex, label, required = TRUE) {
  direct <- preferred[preferred %in% names(df)]
  if (length(direct) >= 1) return(direct[[1]])
  
  candidates <- grep(regex, names(df), value = TRUE, ignore.case = TRUE)
  
  if (length(candidates) == 1) return(candidates[[1]])
  if (!required) return(NA_character_)
  
  if (length(candidates) > 1) {
    stop(
      "Multiple candidate columns found for ", label, ": ",
      paste(candidates, collapse = ", ")
    )
  }
  
  stop(
    "Could not identify ", label, ". Available columns include:\n",
    paste(head(names(df), 200), collapse = ", ")
  )
}

clean_character <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "nan", "None", "null", "<NA>")] <- NA_character_
  x
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(gsub(",", "", clean_character(x), fixed = TRUE)))
}

normalize_visit <- function(x) {
  x <- clean_character(x)
  numeric_visit <- suppressWarnings(as.integer(x))
  ifelse(!is.na(numeric_visit), as.character(numeric_visit), toupper(x))
}

normalize_treatment <- function(x) {
  x_upper <- toupper(clean_character(x))
  
  case_when(
    str_detect(x_upper, "PLACEBO") ~ "Placebo",
    str_detect(x_upper, "SOLANEZUMAB|ACTIVE|DRUG|TREAT") ~ "Drug",
    TRUE ~ "Not randomized/Unknown"
  )
}

normalize_binary_amyloid <- function(x) {
  x_upper <- toupper(clean_character(x))
  
  case_when(
    x_upper %in% c(
      "1", "POS", "POSITIVE", "AMYLOID POSITIVE", "A+",
      "ABNORMAL", "ELEVATED", "YES", "Y", "TRUE"
    ) ~ TRUE,
    x_upper %in% c(
      "0", "NEG", "NEGATIVE", "AMYLOID NEGATIVE", "A-",
      "NORMAL", "NOT ELEVATED", "NO", "N", "FALSE"
    ) ~ FALSE,
    TRUE ~ NA
  )
}


detect_first_amyloid_column <- function(
    df,
    preferred,
    regex,
    label,
    required = FALSE
) {
  direct <- preferred[preferred %in% names(df)]
  
  if (length(direct) >= 1) {
    message(label, " column selected by preferred name: ", direct[[1]])
    return(direct[[1]])
  }
  
  candidates <- grep(
    regex,
    names(df),
    value = TRUE,
    ignore.case = TRUE
  )
  
  if (length(candidates) >= 1) {
    # Prefer global/composite/summary measures and whole-cerebellum references.
    score <- rep(0, length(candidates))
    score <- score + 4 * str_detect(
      toupper(candidates),
      "GLOBAL|COMPOSITE|SUMMARY|CORTICAL"
    )
    score <- score + 3 * str_detect(
      toupper(candidates),
      "WHOLE.*CEREB|WCEREB|CEREBELLUM"
    )
    score <- score + 2 * str_detect(
      toupper(candidates),
      "CENTILOID|SUVR"
    )
    score <- score - 4 * str_detect(
      toupper(candidates),
      "REGION|ROI|LEFT|RIGHT|FRONTAL|PARIETAL|TEMPORAL"
    )
    
    selected <- candidates[order(score, decreasing = TRUE)][[1]]
    
    if (length(candidates) > 1) {
      warning(
        "Multiple candidate columns found for ", label, ": ",
        paste(candidates, collapse = ", "),
        ". Selected: ", selected,
        ". Review the column audit output."
      )
    } else {
      message(label, " column selected by regex: ", selected)
    }
    
    return(selected)
  }
  
  if (required) {
    stop(
      "Could not identify ", label, ". Available columns include:\n",
      paste(names(df), collapse = ", ")
    )
  }
  
  NA_character_
}

normalize_visual_read <- function(x) {
  x_upper <- toupper(clean_character(x))
  
  case_when(
    x_upper %in% c(
      "1", "POS", "POSITIVE", "AMYLOID POSITIVE",
      "ABNORMAL", "ELEVATED", "YES", "Y", "TRUE"
    ) ~ TRUE,
    x_upper %in% c(
      "0", "NEG", "NEGATIVE", "AMYLOID NEGATIVE",
      "NORMAL", "NOT ELEVATED", "NO", "N", "FALSE"
    ) ~ FALSE,
    str_detect(x_upper, "POS|ELEVATED|ABNORMAL") ~ TRUE,
    str_detect(x_upper, "NEG|NOT ELEVATED|NORMAL") ~ FALSE,
    TRUE ~ NA
  )
}

format_p <- function(p) {
  case_when(
    is.na(p) ~ "P = NA",
    p < 2.2e-16 ~ "P < 2.2e-16",
    p < 0.001 ~ paste0("P = ", formatC(p, format = "e", digits = 2)),
    TRUE ~ paste0("P = ", signif(p, 3))
  )
}

format_beta <- function(x) {
  case_when(
    is.na(x) ~ "NA",
    abs(x) < 0.001 ~ formatC(x, format = "e", digits = 2),
    TRUE ~ sprintf("%.3f", x)
  )
}

fit_group_slope <- function(d, time_col = "years_since_epoch_baseline") {
  d <- d |>
    filter(
      is.finite(.data[[time_col]]),
      is.finite(acceleration_years)
    )
  
  n_subjects <- n_distinct(d$participant_id)
  n_scans <- nrow(d)
  
  if (
    n_subjects < 2 ||
    n_scans < 4 ||
    n_distinct(d[[time_col]]) < 2
  ) {
    return(tibble(
      n_subjects = n_subjects,
      n_scans = n_scans,
      model = "insufficient data",
      beta_slope_per_year = NA_real_,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_
    ))
  }
  
  formula_obj <- as.formula(
    paste0("acceleration_years ~ ", time_col)
  )
  
  if (requireNamespace("nlme", quietly = TRUE)) {
    fit <- tryCatch(
      nlme::lme(
        fixed = formula_obj,
        random = ~ 1 | participant_id,
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
    
    if (!is.null(fit)) {
      tt <- summary(fit)$tTable
      if (time_col %in% rownames(tt)) {
        beta <- unname(tt[time_col, "Value"])
        se <- unname(tt[time_col, "Std.Error"])
        p_value <- unname(tt[time_col, "p-value"])
        model_df <- unname(tt[time_col, "DF"])
        
        return(tibble(
          n_subjects = n_subjects,
          n_scans = n_scans,
          model = "linear mixed model with participant random intercept",
          beta_slope_per_year = beta,
          se = se,
          ci_low = beta - qt(0.975, df = model_df) * se,
          ci_high = beta + qt(0.975, df = model_df) * se,
          p_value = p_value
        ))
      }
    }
  }
  
  fit <- lm(formula_obj, data = d)
  tt <- summary(fit)$coefficients
  beta <- unname(tt[time_col, "Estimate"])
  se <- unname(tt[time_col, "Std. Error"])
  p_value <- unname(tt[time_col, "Pr(>|t|)"])
  
  tibble(
    n_subjects = n_subjects,
    n_scans = n_scans,
    model = "ordinary linear model fallback",
    beta_slope_per_year = beta,
    se = se,
    ci_low = beta - qt(0.975, df = fit$df.residual) * se,
    ci_high = beta + qt(0.975, df = fit$df.residual) * se,
    p_value = p_value
  )
}

fit_interaction_model <- function(d, group_col, output_label) {
  d <- d |>
    filter(
      is.finite(years_since_epoch_baseline),
      is.finite(acceleration_years),
      !is.na(.data[[group_col]])
    )
  
  group_counts <- d |>
    distinct(participant_id, .data[[group_col]]) |>
    count(.data[[group_col]], name = "n_subjects")
  
  if (
    nrow(group_counts) < 2 ||
    any(group_counts$n_subjects < 2) ||
    n_distinct(d$years_since_epoch_baseline) < 2
  ) {
    return(tibble(
      analysis = output_label,
      model = "insufficient data",
      term = NA_character_,
      estimate = NA_real_,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_
    ))
  }
  
  formula_obj <- as.formula(
    paste0(
      "acceleration_years ~ years_since_epoch_baseline * ",
      group_col
    )
  )
  
  if (requireNamespace("nlme", quietly = TRUE)) {
    fit <- tryCatch(
      nlme::lme(
        fixed = formula_obj,
        random = ~ 1 | participant_id,
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
    
    if (!is.null(fit)) {
      tt <- summary(fit)$tTable
      interaction_rows <- grep(
        "years_since_epoch_baseline.*:",
        rownames(tt),
        value = TRUE
      )
      
      if (length(interaction_rows) > 0) {
        return(map_dfr(interaction_rows, function(term_name) {
          beta <- unname(tt[term_name, "Value"])
          se <- unname(tt[term_name, "Std.Error"])
          model_df <- unname(tt[term_name, "DF"])
          
          tibble(
            analysis = output_label,
            model = "linear mixed model with participant random intercept",
            term = term_name,
            estimate = beta,
            se = se,
            ci_low = beta - qt(0.975, df = model_df) * se,
            ci_high = beta + qt(0.975, df = model_df) * se,
            p_value = unname(tt[term_name, "p-value"])
          )
        }))
      }
    }
  }
  
  fit <- lm(formula_obj, data = d)
  tt <- summary(fit)$coefficients
  interaction_rows <- grep(
    "years_since_epoch_baseline.*:",
    rownames(tt),
    value = TRUE
  )
  
  if (length(interaction_rows) == 0) {
    interaction_rows <- grep(
      ":years_since_epoch_baseline",
      rownames(tt),
      value = TRUE
    )
  }
  
  map_dfr(interaction_rows, function(term_name) {
    beta <- unname(tt[term_name, "Estimate"])
    se <- unname(tt[term_name, "Std. Error"])
    
    tibble(
      analysis = output_label,
      model = "ordinary linear model fallback",
      term = term_name,
      estimate = beta,
      se = se,
      ci_low = beta - qt(0.975, df = fit$df.residual) * se,
      ci_high = beta + qt(0.975, df = fit$df.residual) * se,
      p_value = unname(tt[term_name, "Pr(>|t|)"])
    )
  })
}

# ------------------------------------------------------------------------------
# 5. Validate and read inputs
# ------------------------------------------------------------------------------

for (path in c(a4_prediction_file, subjinfo_file, amyloid_file)) {
  if (!file.exists(path)) {
    stop("Required input file does not exist: ", path)
  }
}

prediction_df <- readr::read_tsv(
  a4_prediction_file,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "unique"
)

subjinfo_df <- readr::read_csv(
  subjinfo_file,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "unique"
)

amyloid_df <- readr::read_csv(
  amyloid_file,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "unique"
)

# ------------------------------------------------------------------------------
# 6. Prepare participant-level treatment metadata
# ------------------------------------------------------------------------------

subj_id_col <- detect_column(
  subjinfo_df,
  preferred = c("BID", "participant_id", "PTID"),
  regex = "(^bid$|participant.*id|^ptid$)",
  label = "SUBJINFO participant ID"
)

tx_col <- detect_column(
  subjinfo_df,
  preferred = c("TX", "TRT", "TREATMENT", "ARM"),
  regex = "(^tx$|treat|trt|arm)",
  label = "SUBJINFO treatment assignment"
)

subj_age_col <- detect_column(
  subjinfo_df,
  preferred = c("AGEYR", "Age", "AGE"),
  regex = "(^ageyr$|^age$)",
  label = "SUBJINFO baseline age",
  required = FALSE
)

subj_sex_col <- detect_column(
  subjinfo_df,
  preferred = c("SEX", "Sex"),
  regex = "^sex$",
  label = "SUBJINFO sex",
  required = FALSE
)

participant_metadata <- subjinfo_df |>
  transmute(
    participant_id = as.character(.data[[subj_id_col]]),
    treatment_original = clean_character(.data[[tx_col]]),
    treatment_group = normalize_treatment(.data[[tx_col]]),
    baseline_age_subjinfo = if (!is.na(subj_age_col)) {
      safe_numeric(.data[[subj_age_col]])
    } else {
      NA_real_
    },
    sex_subjinfo = if (!is.na(subj_sex_col)) {
      clean_character(.data[[subj_sex_col]])
    } else {
      NA_character_
    }
  ) |>
  filter(!is.na(participant_id)) |>
  distinct(participant_id, .keep_all = TRUE)

# ------------------------------------------------------------------------------
# 7. Prepare A4 longitudinal EPOCH scans
# ------------------------------------------------------------------------------

pred_id_col <- detect_column(
  prediction_df,
  preferred = c("BID", "participant_id", "PTID"),
  regex = "(^bid$|participant.*id|^ptid$)",
  label = "prediction participant ID"
)

pred_visit_col <- detect_column(
  prediction_df,
  preferred = c(
    "Visit_Code", "MUSE_VISIT", "VISCODE",
    "visit_code", "Visit"
  ),
  regex = "(visit|viscode)",
  label = "prediction visit",
  required = FALSE
)

pred_age_col <- detect_column(
  prediction_df,
  preferred = c("Age", "age_at_scan_used_for_model", "AGE"),
  regex = "(^|_)age($|_at_scan)",
  label = "prediction age"
)

pred_years_col <- detect_column(
  prediction_df,
  preferred = c(
    "years_since_external_baseline",
    "years_since_baseline"
  ),
  regex = "years.*baseline",
  label = "prediction years since baseline",
  required = FALSE
)

pred_days_col <- detect_column(
  prediction_df,
  preferred = c(
    "Date_DAYS_CONSENT",
    "scan_days_from_consent",
    "days_from_consent"
  ),
  regex = "(date.*days.*consent|scan.*days.*consent|days.*consent)",
  label = "prediction days from consent",
  required = FALSE
)

acceleration_years_col <- detect_column(
  prediction_df,
  preferred = c(
    "adni_brain_mri_ad_epoch_acceleration_years",
    "adni_brain_mri_ad_lepoch_acceleration_years"
  ),
  regex = "acceleration[_\\.]*years$",
  label = "AD EPOCH acceleration years"
)

acceleration_z_col <- detect_column(
  prediction_df,
  preferred = c(
    "adni_brain_mri_ad_epoch_acceleration_z",
    "adni_brain_mri_ad_lepoch_acceleration_z"
  ),
  regex = "acceleration[_\\.]*z$",
  label = "AD EPOCH acceleration z",
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

epoch_scans <- prediction_df |>
  transmute(
    prediction_source_row = row_number(),
    participant_id = as.character(.data[[pred_id_col]]),
    visit_code = if (!is.na(pred_visit_col)) {
      normalize_visit(.data[[pred_visit_col]])
    } else {
      NA_character_
    },
    scan_age = safe_numeric(.data[[pred_age_col]]),
    years_from_consent = case_when(
      !is.na(pred_days_col) ~ safe_numeric(.data[[pred_days_col]]) / 365.25,
      !is.na(pred_years_col) ~ safe_numeric(.data[[pred_years_col]]),
      TRUE ~ NA_real_
    ),
    acceleration_years = safe_numeric(.data[[acceleration_years_col]]),
    acceleration_z = if (!is.na(acceleration_z_col)) {
      safe_numeric(.data[[acceleration_z_col]])
    } else {
      NA_real_
    },
    risk_score = if (!is.na(risk_score_col)) {
      safe_numeric(.data[[risk_score_col]])
    } else {
      NA_real_
    }
  ) |>
  filter(
    !is.na(participant_id),
    is.finite(acceleration_years)
  ) |>
  left_join(
    participant_metadata,
    by = "participant_id",
    relationship = "many-to-one"
  ) |>
  group_by(participant_id) |>
  arrange(
    years_from_consent,
    scan_age,
    visit_code,
    .by_group = TRUE
  ) |>
  mutate(
    epoch_scan_number = row_number(),
    baseline_epoch_age = first(scan_age),
    baseline_epoch_year_from_consent = first(years_from_consent),
    years_since_epoch_baseline = case_when(
      is.finite(years_from_consent) &
        is.finite(baseline_epoch_year_from_consent) ~
        years_from_consent - baseline_epoch_year_from_consent,
      is.finite(scan_age) &
        is.finite(baseline_epoch_age) ~
        scan_age - baseline_epoch_age,
      TRUE ~ NA_real_
    )
  ) |>
  ungroup()

# Remove duplicated participant-time records.
epoch_scans <- epoch_scans |>
  arrange(
    participant_id,
    years_since_epoch_baseline,
    prediction_source_row
  ) |>
  group_by(participant_id, years_since_epoch_baseline) |>
  slice_head(n = 1) |>
  ungroup()

# ------------------------------------------------------------------------------
# 8. Prepare longitudinal amyloid measurements
# ------------------------------------------------------------------------------

amy_id_col <- detect_column(
  amyloid_df,
  preferred = c("BID", "participant_id", "PTID"),
  regex = "(^bid$|participant.*id|^ptid$)",
  label = "amyloid participant ID"
)

amy_visit_col <- detect_column(
  amyloid_df,
  preferred = c("VISCODE", "VISITCD", "Visit_Code", "visit_code"),
  regex = "(viscode|visit)",
  label = "amyloid visit code",
  required = FALSE
)

amy_days_col <- detect_column(
  amyloid_df,
  preferred = c(
    "Date_DAYS_CONSENT",
    "EXAMDATE_DAYS_CONSENT",
    "SCANDATE_DAYS_CONSENT"
  ),
  regex = "(date.*days.*consent|exam.*days.*consent|scan.*days.*consent)",
  label = "amyloid assessment days from consent",
  required = FALSE
)

amy_centiloid_col <- detect_first_amyloid_column(
  amyloid_df,
  preferred = c(
    "AMYLCENT", "CENTILOID", "Centiloid", "CENTILOIDS",
    "CL", "Composite_Centiloid", "GLOBAL_CENTILOID",
    "SUMMARY_CENTILOID", "CENTILOID_WHOLECEREBELLUM"
  ),
  regex = "(amy.*cent|centiloid|^cl$)",
  label = "amyloid Centiloid",
  required = FALSE
)

amy_suvr_col <- detect_first_amyloid_column(
  amyloid_df,
  preferred = c(
    "SUVRCER", "SUVR", "COMPOSITE_SUVR",
    "Global_SUVR", "GLOBAL_SUVR", "SUMMARY_SUVR",
    "CORTICAL_SUMMARY_SUVR", "CORTICAL_COMPOSITE_SUVR",
    "WHOLECEREBELLUM_SUVR", "WHOLE_CEREBELLUM_SUVR",
    "FLORBETAPIR_SUVR"
  ),
  regex = "(suvr|standardized.*uptake.*ratio)",
  label = "amyloid SUVR",
  required = FALSE
)

amy_positive_col <- detect_first_amyloid_column(
  amyloid_df,
  preferred = c(
    "AMYLOID_POSITIVE", "AMYLOID_STATUS", "AMYLPOS",
    "PET_POSITIVE", "POSITIVE", "AB_STATUS", "AMYLOID_ELEVATED",
    "ELEVATED_AMYLOID", "AMYLOID_CLASS"
  ),
  regex = "(amyloid.*pos|amyloid.*status|amylpos|pet.*pos|ab.*status|amyloid.*elev|elev.*amyloid)",
  label = "amyloid positivity indicator",
  required = FALSE
)

amy_visual_col <- detect_first_amyloid_column(
  amyloid_df,
  preferred = c(
    "VISUAL_READ", "VISREAD", "VISUALREAD", "PET_VISUAL_READ",
    "AMYLOID_VISUAL_READ", "CENTRAL_VISUAL_READ", "READRESULT",
    "READ_RESULT", "QUALITATIVE_READ"
  ),
  regex = "(visual.*read|visread|read.*result|qualitative.*read)",
  label = "amyloid visual-read indicator",
  required = FALSE
)

# Always save the available amyloid columns for transparent review.
amyloid_column_audit <- tibble(
  role = c(
    "participant ID",
    "visit",
    "days from consent",
    "reported positivity",
    "visual read",
    "Centiloid",
    "SUVR"
  ),
  selected_column = c(
    amy_id_col,
    amy_visit_col,
    amy_days_col,
    amy_positive_col,
    amy_visual_col,
    amy_centiloid_col,
    amy_suvr_col
  )
)

write_tsv(
  amyloid_column_audit,
  file.path(
    out_dir,
    paste0(prefix, "_amyloid_column_audit.tsv")
  )
)

write_tsv(
  tibble(column_name = names(amyloid_df)),
  file.path(
    out_dir,
    paste0(prefix, "_all_amyloid_input_columns.tsv")
  )
)

if (
  is.na(amy_positive_col) &&
  is.na(amy_centiloid_col) &&
  is.na(amy_suvr_col)
) {
  stop(
    "No usable amyloid positivity, Centiloid, or SUVR column was identified. ",
    "The complete input-column list was written to: ",
    file.path(
      out_dir,
      paste0(prefix, "_all_amyloid_input_columns.tsv")
    )
  )
}

amyloid_long <- amyloid_df |>
  transmute(
    participant_id = as.character(.data[[amy_id_col]]),
    amyloid_visit_code = if (!is.na(amy_visit_col)) {
      normalize_visit(.data[[amy_visit_col]])
    } else {
      NA_character_
    },
    amyloid_year_from_consent = if (!is.na(amy_days_col)) {
      safe_numeric(.data[[amy_days_col]]) / 365.25
    } else {
      NA_real_
    },
    centiloid = if (!is.na(amy_centiloid_col)) {
      safe_numeric(.data[[amy_centiloid_col]])
    } else {
      NA_real_
    },
    suvr = if (!is.na(amy_suvr_col)) {
      safe_numeric(.data[[amy_suvr_col]])
    } else {
      NA_real_
    },
    amyloid_positive_reported = if (!is.na(amy_positive_col)) {
      normalize_binary_amyloid(.data[[amy_positive_col]])
    } else {
      NA
    },
    amyloid_visual_positive = if (!is.na(amy_visual_col)) {
      normalize_visual_read(.data[[amy_visual_col]])
    } else {
      NA
    }
  ) |>
  mutate(
    # Hierarchy:
    # 1. Use an explicit reported positivity field when available.
    # 2. Reproduce the A4 screening algorithm using SUVR plus visual read.
    # 3. Apply the equivalent Centiloid bands plus visual read.
    amyloid_positive = case_when(
      !is.na(amyloid_positive_reported) ~
        amyloid_positive_reported,
      
      is.finite(suvr) &
        suvr >= amyloid_positive_suvr_threshold ~
        TRUE,
      
      is.finite(suvr) &
        suvr >= amyloid_borderline_suvr_lower &
        suvr < amyloid_positive_suvr_threshold &
        amyloid_visual_positive %in% TRUE ~
        TRUE,
      
      is.finite(suvr) &
        suvr >= amyloid_borderline_suvr_lower &
        suvr < amyloid_positive_suvr_threshold &
        amyloid_visual_positive %in% FALSE ~
        FALSE,
      
      is.finite(suvr) &
        suvr < amyloid_borderline_suvr_lower ~
        FALSE,
      
      is.finite(centiloid) &
        centiloid >= amyloid_positive_centiloid_threshold ~
        TRUE,
      
      is.finite(centiloid) &
        centiloid >= amyloid_borderline_centiloid_lower &
        centiloid < amyloid_positive_centiloid_threshold &
        amyloid_visual_positive %in% TRUE ~
        TRUE,
      
      is.finite(centiloid) &
        centiloid >= amyloid_borderline_centiloid_lower &
        centiloid < amyloid_positive_centiloid_threshold &
        amyloid_visual_positive %in% FALSE ~
        FALSE,
      
      is.finite(centiloid) &
        centiloid < amyloid_borderline_centiloid_lower ~
        FALSE,
      
      TRUE ~ NA
    ),
    
    amyloid_classification_source = case_when(
      !is.na(amyloid_positive_reported) ~
        paste0("Reported field: ", amy_positive_col),
      
      is.finite(suvr) &
        suvr >= amyloid_positive_suvr_threshold ~
        paste0(
          "A4 quantitative rule: SUVR >= ",
          amyloid_positive_suvr_threshold
        ),
      
      is.finite(suvr) &
        suvr >= amyloid_borderline_suvr_lower &
        suvr < amyloid_positive_suvr_threshold &
        !is.na(amyloid_visual_positive) ~
        paste0(
          "A4 borderline rule: SUVR ",
          amyloid_borderline_suvr_lower,
          " to <",
          amyloid_positive_suvr_threshold,
          " plus visual read"
        ),
      
      is.finite(suvr) &
        suvr < amyloid_borderline_suvr_lower ~
        paste0(
          "A4 quantitative rule: SUVR < ",
          amyloid_borderline_suvr_lower
        ),
      
      is.finite(centiloid) &
        centiloid >= amyloid_positive_centiloid_threshold ~
        paste0(
          "SUVR-equivalent rule: Centiloid >= ",
          amyloid_positive_centiloid_threshold
        ),
      
      is.finite(centiloid) &
        centiloid >= amyloid_borderline_centiloid_lower &
        centiloid < amyloid_positive_centiloid_threshold &
        !is.na(amyloid_visual_positive) ~
        paste0(
          "SUVR-equivalent borderline rule: Centiloid ",
          amyloid_borderline_centiloid_lower,
          " to <",
          amyloid_positive_centiloid_threshold,
          " plus visual read"
        ),
      
      is.finite(centiloid) &
        centiloid < amyloid_borderline_centiloid_lower ~
        paste0(
          "SUVR-equivalent rule: Centiloid < ",
          amyloid_borderline_centiloid_lower
        ),
      
      is.finite(suvr) &
        suvr >= amyloid_borderline_suvr_lower &
        suvr < amyloid_positive_suvr_threshold &
        is.na(amyloid_visual_positive) ~
        "Borderline SUVR but visual read unavailable",
      
      is.finite(centiloid) &
        centiloid >= amyloid_borderline_centiloid_lower &
        centiloid < amyloid_positive_centiloid_threshold &
        is.na(amyloid_visual_positive) ~
        "Borderline Centiloid but visual read unavailable",
      
      TRUE ~ "Unavailable"
    )
  ) |>
  filter(
    !is.na(participant_id),
    !is.na(amyloid_positive)
  ) |>
  arrange(
    participant_id,
    amyloid_year_from_consent,
    amyloid_visit_code
  ) |>
  distinct(
    participant_id,
    amyloid_year_from_consent,
    amyloid_visit_code,
    .keep_all = TRUE
  )


amyloid_classification_summary <- amyloid_long |>
  count(
    amyloid_classification_source,
    amyloid_positive,
    name = "n_assessments"
  )

write_tsv(
  amyloid_classification_summary,
  file.path(
    out_dir,
    paste0(prefix, "_amyloid_classification_summary.tsv")
  )
)

# ------------------------------------------------------------------------------
# 9. Define participant-level amyloid trajectories
# ------------------------------------------------------------------------------

amyloid_subject_status <- amyloid_long |>
  group_by(participant_id) |>
  arrange(
    amyloid_year_from_consent,
    amyloid_visit_code,
    .by_group = TRUE
  ) |>
  summarise(
    n_amyloid_assessments = n(),
    earliest_amyloid_positive = first(amyloid_positive),
    earliest_amyloid_year_from_consent = first(amyloid_year_from_consent),
    latest_amyloid_year_from_consent = last(amyloid_year_from_consent),
    ever_amyloid_positive = any(amyloid_positive, na.rm = TRUE),
    conversion_amyloid_year_from_consent = if (
      !first(amyloid_positive) &&
      any(amyloid_positive, na.rm = TRUE)
    ) {
      min(
        amyloid_year_from_consent[amyloid_positive %in% TRUE],
        na.rm = TRUE
      )
    } else {
      NA_real_
    },
    amyloid_trajectory = case_when(
      first(amyloid_positive) ~ "Amyloid positive at baseline",
      !first(amyloid_positive) &
        any(amyloid_positive, na.rm = TRUE) ~
        "Converted to amyloid positive",
      !first(amyloid_positive) &
        !any(amyloid_positive, na.rm = TRUE) ~
        "Remained amyloid negative",
      TRUE ~ "Amyloid status unavailable"
    ),
    .groups = "drop"
  ) |>
  mutate(
    amyloid_trajectory = factor(
      amyloid_trajectory,
      levels = c(
        "Remained amyloid negative",
        "Converted to amyloid positive",
        "Amyloid positive at baseline",
        "Amyloid status unavailable"
      )
    )
  )

# ------------------------------------------------------------------------------
# 10. Assign nearest amyloid assessment to each EPOCH scan
# ------------------------------------------------------------------------------

amyloid_match_candidates <- epoch_scans |>
  select(
    prediction_source_row,
    participant_id,
    visit_code,
    years_from_consent
  ) |>
  inner_join(
    amyloid_long,
    by = "participant_id",
    relationship = "many-to-many"
  ) |>
  mutate(
    exact_visit_match = (
      !is.na(visit_code) &
        !is.na(amyloid_visit_code) &
        visit_code == amyloid_visit_code
    ),
    temporal_distance_years = abs(
      years_from_consent - amyloid_year_from_consent
    ),
    match_priority = case_when(
      exact_visit_match ~ 1,
      is.finite(temporal_distance_years) ~ 2,
      TRUE ~ 3
    ),
    match_distance = case_when(
      match_priority == 1 ~ 0,
      match_priority == 2 ~ temporal_distance_years,
      TRUE ~ Inf
    )
  ) |>
  arrange(
    prediction_source_row,
    match_priority,
    match_distance,
    amyloid_year_from_consent
  ) |>
  group_by(prediction_source_row) |>
  slice_head(n = 1) |>
  ungroup() |>
  mutate(
    amyloid_match_accepted = (
      exact_visit_match |
        (
          is.finite(temporal_distance_years) &
            temporal_distance_years <=
            maximum_amyloid_match_distance_years
        )
    )
  ) |>
  transmute(
    prediction_source_row,
    nearest_amyloid_positive = ifelse(
      amyloid_match_accepted,
      amyloid_positive,
      NA
    ),
    nearest_centiloid = ifelse(
      amyloid_match_accepted,
      centiloid,
      NA_real_
    ),
    nearest_suvr = ifelse(
      amyloid_match_accepted,
      suvr,
      NA_real_
    ),
    nearest_amyloid_year_from_consent = ifelse(
      amyloid_match_accepted,
      amyloid_year_from_consent,
      NA_real_
    ),
    amyloid_match_distance_years = ifelse(
      amyloid_match_accepted,
      temporal_distance_years,
      NA_real_
    ),
    amyloid_match_source = case_when(
      !amyloid_match_accepted ~ "No acceptable amyloid match",
      exact_visit_match ~ "Exact visit match",
      TRUE ~ "Nearest assessment by time"
    )
  )

analysis_scans <- epoch_scans |>
  left_join(
    amyloid_subject_status,
    by = "participant_id",
    relationship = "many-to-one"
  ) |>
  left_join(
    amyloid_match_candidates,
    by = "prediction_source_row",
    relationship = "one-to-one"
  ) |>
  mutate(
    # Derive amyloid status at each EPOCH scan using the participant-level
    # longitudinal amyloid trajectory. This avoids classifying most scans as
    # unavailable merely because no PET assessment occurred within the
    # nearest-match window.
    #
    # Rules:
    #   1. Remained negative: all EPOCH scans are negative.
    #   2. Positive at baseline: all EPOCH scans are positive.
    #   3. Converter: scans before conversion are negative and scans at/after
    #      the first positive amyloid assessment are positive.
    #   4. If conversion timing or MRI timing is unavailable, fall back to the
    #      nearest acceptable PET assessment; otherwise retain unavailable.
    scan_amyloid_status = case_when(
      amyloid_trajectory == "Remained amyloid negative" ~
        "Amyloid negative",
      
      amyloid_trajectory == "Amyloid positive at baseline" ~
        "Amyloid positive",
      
      amyloid_trajectory == "Converted to amyloid positive" &
        is.finite(years_from_consent) &
        is.finite(conversion_amyloid_year_from_consent) &
        years_from_consent < conversion_amyloid_year_from_consent ~
        "Amyloid negative",
      
      amyloid_trajectory == "Converted to amyloid positive" &
        is.finite(years_from_consent) &
        is.finite(conversion_amyloid_year_from_consent) &
        years_from_consent >= conversion_amyloid_year_from_consent ~
        "Amyloid positive",
      
      nearest_amyloid_positive %in% TRUE ~
        "Amyloid positive",
      
      nearest_amyloid_positive %in% FALSE ~
        "Amyloid negative",
      
      TRUE ~
        "Amyloid status unavailable"
    ),
    
    scan_amyloid_status_source = case_when(
      amyloid_trajectory == "Remained amyloid negative" ~
        "Participant remained amyloid negative",
      
      amyloid_trajectory == "Amyloid positive at baseline" ~
        "Participant amyloid positive at earliest assessment",
      
      amyloid_trajectory == "Converted to amyloid positive" &
        is.finite(years_from_consent) &
        is.finite(conversion_amyloid_year_from_consent) ~
        "Position relative to first positive amyloid assessment",
      
      !is.na(nearest_amyloid_positive) ~
        "Nearest acceptable amyloid assessment",
      
      TRUE ~
        "Unavailable"
    ),
    
    scan_amyloid_status = factor(
      scan_amyloid_status,
      levels = c(
        "Amyloid negative",
        "Amyloid positive",
        "Amyloid status unavailable"
      )
    ),
    
    treatment_group = factor(
      treatment_group,
      levels = c(
        "Placebo",
        "Drug",
        "Not randomized/Unknown"
      )
    )
  )

# ------------------------------------------------------------------------------
# 11. Restrict to participants with longitudinal EPOCH data
# ------------------------------------------------------------------------------

longitudinal_subjects <- analysis_scans |>
  group_by(participant_id) |>
  summarise(
    treatment_group = first(treatment_group),
    amyloid_trajectory = first(amyloid_trajectory),
    n_scans = n(),
    followup_span_years = max(
      years_since_epoch_baseline,
      na.rm = TRUE
    ) - min(
      years_since_epoch_baseline,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  filter(
    n_scans >= minimum_scans_per_subject,
    is.finite(followup_span_years),
    followup_span_years > minimum_followup_years
  )

traj_df <- analysis_scans |>
  semi_join(
    longitudinal_subjects,
    by = "participant_id"
  ) |>
  arrange(
    participant_id,
    years_since_epoch_baseline
  ) |>
  group_by(participant_id) |>
  mutate(
    scan_relation = ifelse(
      row_number() == 1,
      "Selected EPOCH baseline",
      "Follow-up"
    )
  ) |>
  ungroup()

if (nrow(traj_df) == 0) {
  stop(
    "No A4 participants had at least ",
    minimum_scans_per_subject,
    " qualified longitudinal EPOCH scans."
  )
}

conversion_df <- traj_df |>
  filter(
    amyloid_trajectory %in% c(
      "Remained amyloid negative",
      "Converted to amyloid positive"
    )
  ) |>
  droplevels()

if (
  restrict_conversion_analysis_to_baseline_negative &&
  nrow(conversion_df) == 0
) {
  warning(
    "No baseline-amyloid-negative participants with longitudinal EPOCH data ",
    "were available for the conversion analysis."
  )
}

randomized_df <- traj_df |>
  filter(
    treatment_group %in% c("Placebo", "Drug")
  ) |>
  droplevels()

# ------------------------------------------------------------------------------
# 12. Cohort summaries
# ------------------------------------------------------------------------------

cohort_flow <- tibble(
  stage = c(
    "Participants with scored A4 EPOCH scans",
    paste0(
      "Participants with >=",
      minimum_scans_per_subject,
      " longitudinal EPOCH scans"
    ),
    "Randomized placebo participants",
    "Randomized active-drug participants",
    "Baseline amyloid-negative participants who remained negative",
    "Baseline amyloid-negative participants who converted to positive",
    "Participants amyloid positive at earliest assessment"
  ),
  n_participants = c(
    n_distinct(epoch_scans$participant_id),
    n_distinct(traj_df$participant_id),
    n_distinct(
      traj_df$participant_id[
        traj_df$treatment_group == "Placebo"
      ]
    ),
    n_distinct(
      traj_df$participant_id[
        traj_df$treatment_group == "Drug"
      ]
    ),
    n_distinct(
      traj_df$participant_id[
        traj_df$amyloid_trajectory ==
          "Remained amyloid negative"
      ]
    ),
    n_distinct(
      traj_df$participant_id[
        traj_df$amyloid_trajectory ==
          "Converted to amyloid positive"
      ]
    ),
    n_distinct(
      traj_df$participant_id[
        traj_df$amyloid_trajectory ==
          "Amyloid positive at baseline"
      ]
    )
  )
)

treatment_amyloid_counts <- traj_df |>
  distinct(
    participant_id,
    treatment_group,
    amyloid_trajectory
  ) |>
  count(
    treatment_group,
    amyloid_trajectory,
    name = "n_participants",
    .drop = FALSE
  )

subject_summary <- traj_df |>
  group_by(participant_id) |>
  arrange(
    years_since_epoch_baseline,
    .by_group = TRUE
  ) |>
  summarise(
    treatment_group = first(treatment_group),
    treatment_original = first(treatment_original),
    amyloid_trajectory = first(amyloid_trajectory),
    n_amyloid_assessments = first(n_amyloid_assessments),
    conversion_amyloid_year_from_consent = first(
      conversion_amyloid_year_from_consent
    ),
    n_scans = n(),
    followup_span_years = max(years_since_epoch_baseline) -
      min(years_since_epoch_baseline),
    baseline_acceleration_years = first(acceleration_years),
    last_acceleration_years = last(acceleration_years),
    change_last_minus_baseline = last_acceleration_years -
      baseline_acceleration_years,
    .groups = "drop"
  )


scan_amyloid_status_summary <- traj_df |>
  count(
    amyloid_trajectory,
    scan_amyloid_status,
    scan_amyloid_status_source,
    name = "n_scans",
    .drop = FALSE
  )

# ------------------------------------------------------------------------------
# 13. Mixed-model trend estimates
# ------------------------------------------------------------------------------

treatment_trends <- randomized_df |>
  group_by(treatment_group) |>
  group_modify(~ fit_group_slope(.x)) |>
  ungroup() |>
  mutate(
    group_type = "Treatment arm",
    group = as.character(treatment_group),
    .before = 1
  )

amyloid_trends <- conversion_df |>
  group_by(amyloid_trajectory) |>
  group_modify(~ fit_group_slope(.x)) |>
  ungroup() |>
  mutate(
    group_type = "Amyloid trajectory",
    group = as.character(amyloid_trajectory),
    .before = 1
  )

interaction_models <- bind_rows(
  fit_interaction_model(
    randomized_df,
    "treatment_group",
    "Time-by-treatment interaction"
  ),
  fit_interaction_model(
    conversion_df,
    "amyloid_trajectory",
    "Time-by-amyloid-trajectory interaction"
  )
)

trend_annotations <- bind_rows(
  treatment_trends |>
    transmute(
      plot = "Treatment",
      group,
      n_subjects,
      n_scans,
      beta_slope_per_year,
      ci_low,
      ci_high,
      p_value,
      label = paste0(
        group,
        "\nN = ", n_subjects,
        "; scans = ", n_scans,
        "\nβ = ", format_beta(beta_slope_per_year), " y/y",
        "\n95% CI: ", format_beta(ci_low), " to ", format_beta(ci_high),
        "\n", format_p(p_value)
      )
    ),
  amyloid_trends |>
    transmute(
      plot = "Amyloid",
      group,
      n_subjects,
      n_scans,
      beta_slope_per_year,
      ci_low,
      ci_high,
      p_value,
      label = paste0(
        group,
        "\nN = ", n_subjects,
        "; scans = ", n_scans,
        "\nβ = ", format_beta(beta_slope_per_year), " y/y",
        "\n95% CI: ", format_beta(ci_low), " to ", format_beta(ci_high),
        "\n", format_p(p_value)
      )
    )
)

# ------------------------------------------------------------------------------
# 14. Plot 1: randomized treatment-arm trajectories
# ------------------------------------------------------------------------------

treatment_plot_ranges <- randomized_df |>
  summarise(
    x_min = min(years_since_epoch_baseline, na.rm = TRUE),
    x_max = max(years_since_epoch_baseline, na.rm = TRUE),
    y_min = min(acceleration_years, na.rm = TRUE),
    y_max = max(acceleration_years, na.rm = TRUE)
  ) |>
  mutate(
    x_range = pmax(x_max - x_min, 1),
    y_range = pmax(y_max - y_min, 1)
  )

treatment_annotation_positions <- trend_annotations |>
  filter(plot == "Treatment") |>
  arrange(group) |>
  mutate(
    annotation_index = row_number(),
    x = treatment_plot_ranges$x_min +
      0.03 * treatment_plot_ranges$x_range,
    y = treatment_plot_ranges$y_max -
      (
        0.03 +
          0.22 * (annotation_index - 1)
      ) * treatment_plot_ranges$y_range
  )

p_treatment <- ggplot(
  randomized_df,
  aes(
    x = years_since_epoch_baseline,
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
  geom_line(
    color = "grey65",
    alpha = 0.25,
    linewidth = 0.45
  ) +
  geom_point(
    aes(
      fill = treatment_group,
      shape = scan_relation
    ),
    color = "grey20",
    alpha = 0.78,
    size = 2.0,
    stroke = 0.35
  ) +
  geom_smooth(
    aes(
      color = treatment_group,
      group = treatment_group
    ),
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 1.25,
    alpha = 0.15
  ) +
  geom_label(
    data = treatment_annotation_positions,
    aes(
      x = x,
      y = y,
      label = label,
      color = group
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 2.85,
    lineheight = 0.98,
    label.size = 0.25,
    fill = "white",
    alpha = 0.92,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = treatment_palette[c("Placebo", "Drug")],
    name = "Treatment arm"
  ) +
  scale_fill_manual(
    values = treatment_palette[c("Placebo", "Drug")],
    name = "Treatment arm"
  ) +
  scale_shape_manual(
    values = c(
      "Selected EPOCH baseline" = 21,
      "Follow-up" = 24
    ),
    name = "Scan relation"
  ) +
  scale_x_continuous(
    name = "Years since first qualified A4 EPOCH scan",
    breaks = pretty_breaks(n = 6)
  ) +
  scale_y_continuous(
    name = "Raw-MUSE AD EPOCH acceleration (years)",
    labels = number_format(accuracy = 0.1),
    expand = expansion(mult = c(0.08, 0.25))
  ) +
  labs(
    title = "Longitudinal AD EPOCH trajectories by A4 treatment arm",
    subtitle = paste0(
      "Placebo and active-drug participants with at least ",
      minimum_scans_per_subject,
      " qualified EPOCH scans"
    ),
    caption = paste0(
      "Grey lines represent individual participants. Colored lines are ",
      "group-specific linear trends. β denotes the estimated annual change ",
      "from a participant-random-intercept mixed model when available."
    )
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10.3),
    plot.caption = element_text(size = 8.5, hjust = 0),
    legend.position = "top",
    legend.box = "vertical",
    legend.justification = "left",
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    panel.grid.major.y = element_line(
      color = "grey90",
      linewidth = 0.35
    ),
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------------------------------------
# 15. Plot 2: amyloid conversion trajectories among baseline-negative subjects
# ------------------------------------------------------------------------------

if (nrow(conversion_df) > 0) {
  amyloid_plot_ranges <- conversion_df |>
    summarise(
      x_min = min(years_since_epoch_baseline, na.rm = TRUE),
      x_max = max(years_since_epoch_baseline, na.rm = TRUE),
      y_min = min(acceleration_years, na.rm = TRUE),
      y_max = max(acceleration_years, na.rm = TRUE)
    ) |>
    mutate(
      x_range = pmax(x_max - x_min, 1),
      y_range = pmax(y_max - y_min, 1)
    )
  
  amyloid_annotation_positions <- trend_annotations |>
    filter(plot == "Amyloid") |>
    arrange(group) |>
    mutate(
      annotation_index = row_number(),
      x = amyloid_plot_ranges$x_min +
        0.03 * amyloid_plot_ranges$x_range,
      y = amyloid_plot_ranges$y_max -
        (
          0.03 +
            0.22 * (annotation_index - 1)
        ) * amyloid_plot_ranges$y_range
    )
  
  p_amyloid <- ggplot(
    conversion_df,
    aes(
      x = years_since_epoch_baseline,
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
    geom_line(
      color = "grey65",
      alpha = 0.25,
      linewidth = 0.45
    ) +
    geom_point(
      aes(
        fill = scan_amyloid_status,
        shape = scan_relation
      ),
      color = "grey20",
      alpha = 0.82,
      size = 2.0,
      stroke = 0.35
    ) +
    geom_smooth(
      aes(
        color = amyloid_trajectory,
        group = amyloid_trajectory
      ),
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      linewidth = 1.25,
      alpha = 0.15
    ) +
    geom_label(
      data = amyloid_annotation_positions,
      aes(
        x = x,
        y = y,
        label = label,
        color = group
      ),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 2.75,
      lineheight = 0.98,
      label.size = 0.25,
      fill = "white",
      alpha = 0.92,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = amyloid_palette[
        c(
          "Remained amyloid negative",
          "Converted to amyloid positive"
        )
      ],
      name = "Participant amyloid trajectory"
    ) +
    scale_fill_manual(
      values = scan_amyloid_palette,
      drop = FALSE,
      name = "Amyloid status at EPOCH scan"
    ) +
    scale_shape_manual(
      values = c(
        "Selected EPOCH baseline" = 21,
        "Follow-up" = 24
      ),
      name = "Scan relation"
    ) +
    scale_x_continuous(
      name = "Years since first qualified A4 EPOCH scan",
      breaks = pretty_breaks(n = 6)
    ) +
    scale_y_continuous(
      name = "Raw-MUSE AD EPOCH acceleration (years)",
      labels = number_format(accuracy = 0.1),
      expand = expansion(mult = c(0.08, 0.27))
    ) +
    labs(
      title = "Longitudinal AD EPOCH trajectories by amyloid conversion",
      subtitle = paste0(
        "Restricted to participants who were amyloid negative at their ",
        "earliest available amyloid assessment"
      ),
      caption = paste0(
        "Grey lines represent individual participants. Point fill indicates ",
        "amyloid status at the EPOCH scan, derived from the participant's ",
        "longitudinal amyloid trajectory and first positive assessment. ",
        "Colored fitted lines represent participant-level amyloid trajectories."
      )
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10.3),
      plot.caption = element_text(size = 8.5, hjust = 0),
      legend.position = "top",
      legend.box = "vertical",
      legend.justification = "left",
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      panel.grid.major.y = element_line(
        color = "grey90",
        linewidth = 0.35
      ),
      panel.grid.minor = element_blank()
    )
} else {
  p_amyloid <- ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = paste0(
        "No baseline-amyloid-negative participants with sufficient ",
        "longitudinal EPOCH data were available."
      ),
      size = 5
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void()
}

# ------------------------------------------------------------------------------
# 16. Plot 3: population means in time bins
# ------------------------------------------------------------------------------

population_treatment <- randomized_df |>
  mutate(
    time_bin = round(
      years_since_epoch_baseline / bin_width_years
    ) * bin_width_years
  ) |>
  group_by(treatment_group, time_bin) |>
  summarise(
    n_scans = n(),
    n_subjects = n_distinct(participant_id),
    mean_acceleration_years = mean(acceleration_years, na.rm = TRUE),
    sd_acceleration_years = sd(acceleration_years, na.rm = TRUE),
    se = sd_acceleration_years / sqrt(n_scans),
    ci_low = mean_acceleration_years - 1.96 * se,
    ci_high = mean_acceleration_years + 1.96 * se,
    .groups = "drop"
  ) |>
  filter(n_scans >= 2)

p_population_treatment <- ggplot(
  population_treatment,
  aes(
    x = time_bin,
    y = mean_acceleration_years,
    color = treatment_group,
    fill = treatment_group
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey55"
  ) +
  geom_ribbon(
    aes(ymin = ci_low, ymax = ci_high),
    alpha = 0.15,
    color = NA
  ) +
  geom_line(linewidth = 1.05) +
  geom_point(
    aes(size = n_subjects),
    alpha = 0.90
  ) +
  scale_color_manual(
    values = treatment_palette[c("Placebo", "Drug")]
  ) +
  scale_fill_manual(
    values = treatment_palette[c("Placebo", "Drug")]
  ) +
  scale_size_continuous(
    name = "N subjects",
    range = c(1.8, 4.4)
  ) +
  scale_x_continuous(
    name = "Years since first qualified A4 EPOCH scan",
    breaks = pretty_breaks(n = 6)
  ) +
  scale_y_continuous(
    name = "Mean raw-MUSE AD EPOCH acceleration (years)",
    labels = number_format(accuracy = 0.1)
  ) +
  labs(
    title = "Population AD EPOCH trajectories by treatment arm",
    subtitle = paste0(
      "Means and approximate 95% confidence intervals in ",
      bin_width_years,
      "-year bins"
    ),
    color = "Treatment arm",
    fill = "Treatment arm"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "top",
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    panel.grid.major.y = element_line(
      color = "grey90",
      linewidth = 0.35
    ),
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------------------------------------
# 17. Plot 4: joint treatment-by-amyloid descriptive trajectories
# ------------------------------------------------------------------------------

joint_df <- traj_df |>
  filter(
    treatment_group %in% c("Placebo", "Drug"),
    amyloid_trajectory %in% c(
      "Remained amyloid negative",
      "Converted to amyloid positive",
      "Amyloid positive at baseline"
    )
  ) |>
  droplevels()

joint_cell_counts <- joint_df |>
  distinct(
    participant_id,
    treatment_group,
    amyloid_trajectory
  ) |>
  count(
    treatment_group,
    amyloid_trajectory,
    name = "n_participants"
  )

if (
  nrow(joint_df) > 0 &&
  n_distinct(joint_df$treatment_group) >= 1 &&
  n_distinct(joint_df$amyloid_trajectory) >= 1
) {
  p_joint <- ggplot(
    joint_df,
    aes(
      x = years_since_epoch_baseline,
      y = acceleration_years,
      group = participant_id
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.40,
      color = "grey60"
    ) +
    geom_line(
      color = "grey70",
      alpha = 0.20,
      linewidth = 0.40
    ) +
    geom_point(
      aes(fill = scan_amyloid_status),
      shape = 21,
      color = "grey25",
      size = 1.65,
      alpha = 0.75,
      stroke = 0.30
    ) +
    geom_smooth(
      aes(
        color = amyloid_trajectory,
        group = interaction(
          treatment_group,
          amyloid_trajectory
        )
      ),
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      linewidth = 1.10
    ) +
    facet_wrap(~ treatment_group) +
    scale_color_manual(
      values = amyloid_palette,
      drop = FALSE,
      name = "Amyloid trajectory"
    ) +
    scale_fill_manual(
      values = scan_amyloid_palette,
      drop = FALSE,
      name = "Amyloid status at EPOCH scan"
    ) +
    scale_x_continuous(
      name = "Years since first qualified A4 EPOCH scan",
      breaks = pretty_breaks(n = 5)
    ) +
    scale_y_continuous(
      name = "Raw-MUSE AD EPOCH acceleration (years)",
      labels = number_format(accuracy = 0.1)
    ) +
    labs(
      title = "Joint treatment and amyloid-trajectory view",
      subtitle = paste0(
        "Descriptive visualization; interpret only cells with adequate ",
        "participant representation"
      )
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      legend.position = "top",
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      panel.grid.major.y = element_line(
        color = "grey90",
        linewidth = 0.35
      ),
      panel.grid.minor = element_blank()
    )
} else {
  p_joint <- ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = paste0(
        "No treatment-by-amyloid cells were available for a joint plot."
      ),
      size = 5
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void()
}

# ------------------------------------------------------------------------------
# 18. Subject-level slopes
# ------------------------------------------------------------------------------

subject_slopes <- traj_df |>
  group_by(participant_id) |>
  filter(n() >= minimum_scans_per_subject) |>
  summarise(
    treatment_group = first(treatment_group),
    amyloid_trajectory = first(amyloid_trajectory),
    n_scans = n(),
    followup_span_years = max(years_since_epoch_baseline) -
      min(years_since_epoch_baseline),
    slope_acceleration_years_per_year = {
      fit <- lm(
        acceleration_years ~ years_since_epoch_baseline
      )
      unname(coef(fit)[["years_since_epoch_baseline"]])
    },
    .groups = "drop"
  ) |>
  filter(is.finite(slope_acceleration_years_per_year))

p_slopes_treatment <- subject_slopes |>
  filter(treatment_group %in% c("Placebo", "Drug")) |>
  ggplot(
    aes(
      x = treatment_group,
      y = slope_acceleration_years_per_year,
      fill = treatment_group,
      color = treatment_group
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
  scale_fill_manual(
    values = treatment_palette,
    guide = "none"
  ) +
  scale_color_manual(
    values = treatment_palette,
    guide = "none"
  ) +
  labs(
    x = NULL,
    y = "Subject-level slope\n(acceleration years per year)",
    title = "Subject-level AD EPOCH slopes by treatment arm"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    panel.grid.major.y = element_line(
      color = "grey90",
      linewidth = 0.35
    ),
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------------------------------------
# 19. Save tables
# ------------------------------------------------------------------------------

write_tsv(
  participant_metadata,
  file.path(
    out_dir,
    paste0(prefix, "_participant_treatment_metadata.tsv")
  )
)

write_tsv(
  amyloid_long,
  file.path(
    out_dir,
    paste0(prefix, "_amyloid_longitudinal_classification.tsv")
  )
)

write_tsv(
  amyloid_subject_status,
  file.path(
    out_dir,
    paste0(prefix, "_participant_amyloid_trajectories.tsv")
  )
)

write_tsv(
  traj_df,
  file.path(
    out_dir,
    paste0(prefix, "_qualified_longitudinal_scans.tsv")
  )
)

write_tsv(
  conversion_df,
  file.path(
    out_dir,
    paste0(prefix, "_baseline_negative_conversion_analysis_scans.tsv")
  )
)

write_tsv(
  cohort_flow,
  file.path(
    out_dir,
    paste0(prefix, "_cohort_flow.tsv")
  )
)

write_tsv(
  treatment_amyloid_counts,
  file.path(
    out_dir,
    paste0(prefix, "_treatment_by_amyloid_counts.tsv")
  )
)

write_tsv(
  subject_summary,
  file.path(
    out_dir,
    paste0(prefix, "_subject_summary.tsv")
  )
)


write_tsv(
  scan_amyloid_status_summary,
  file.path(
    out_dir,
    paste0(prefix, "_scan_amyloid_status_summary.tsv")
  )
)

write_tsv(
  treatment_trends,
  file.path(
    out_dir,
    paste0(prefix, "_trend_models_by_treatment.tsv")
  )
)

write_tsv(
  amyloid_trends,
  file.path(
    out_dir,
    paste0(prefix, "_trend_models_by_amyloid_trajectory.tsv")
  )
)

write_tsv(
  interaction_models,
  file.path(
    out_dir,
    paste0(prefix, "_time_interaction_models.tsv")
  )
)

write_tsv(
  population_treatment,
  file.path(
    out_dir,
    paste0(prefix, "_population_timebins_by_treatment.tsv")
  )
)

write_tsv(
  joint_cell_counts,
  file.path(
    out_dir,
    paste0(prefix, "_joint_treatment_amyloid_cell_counts.tsv")
  )
)

write_tsv(
  subject_slopes,
  file.path(
    out_dir,
    paste0(prefix, "_subject_level_slopes.tsv")
  )
)

# ------------------------------------------------------------------------------
# 20. Save figures
# ------------------------------------------------------------------------------

ggsave(
  file.path(
    out_dir,
    paste0(prefix, "_individual_trajectories_by_treatment.pdf")
  ),
  p_treatment,
  width = 9.2,
  height = 6.0
)

ggsave(
  file.path(
    out_dir,
    paste0(prefix, "_individual_trajectories_by_treatment.png")
  ),
  p_treatment,
  width = 9.2,
  height = 6.0,
  dpi = 350,
  bg = "white"
)

ggsave(
  file.path(
    out_dir,
    paste0(prefix, "_individual_trajectories_by_amyloid_conversion.pdf")
  ),
  p_amyloid,
  width = 9.2,
  height = 6.0
)

ggsave(
  file.path(
    out_dir,
    paste0(prefix, "_individual_trajectories_by_amyloid_conversion.png")
  ),
  p_amyloid,
  width = 9.2,
  height = 6.0,
  dpi = 350,
  bg = "white"
)

ggsave(
  file.path(
    out_dir,
    paste0(prefix, "_population_trajectory_by_treatment.pdf")
  ),
  p_population_treatment,
  width = 8.7,
  height = 5.7
)

ggsave(
  file.path(
    out_dir,
    paste0(prefix, "_population_trajectory_by_treatment.png")
  ),
  p_population_treatment,
  width = 8.7,
  height = 5.7,
  dpi = 350,
  bg = "white"
)

ggsave(
  file.path(
    out_dir,
    paste0(prefix, "_joint_treatment_amyloid_trajectory.pdf")
  ),
  p_joint,
  width = 10.2,
  height = 5.9
)

ggsave(
  file.path(
    out_dir,
    paste0(prefix, "_joint_treatment_amyloid_trajectory.png")
  ),
  p_joint,
  width = 10.2,
  height = 5.9,
  dpi = 350,
  bg = "white"
)

ggsave(
  file.path(
    out_dir,
    paste0(prefix, "_subject_level_slopes_by_treatment.pdf")
  ),
  p_slopes_treatment,
  width = 7.2,
  height = 5.6
)

ggsave(
  file.path(
    out_dir,
    paste0(prefix, "_subject_level_slopes_by_treatment.png")
  ),
  p_slopes_treatment,
  width = 7.2,
  height = 5.6,
  dpi = 350,
  bg = "white"
)

# ------------------------------------------------------------------------------
# 21. Print summaries
# ------------------------------------------------------------------------------

message("============================================================")
message("A4 longitudinal AD EPOCH analysis complete.")
message("")
message("Cohort flow:")
print(cohort_flow)
message("")
message("Treatment-by-amyloid participant counts:")
print(treatment_amyloid_counts)
message("")
message("Scan-level amyloid status summary:")
print(scan_amyloid_status_summary)
message("")
message("Treatment-specific longitudinal trends:")
print(treatment_trends)
message("")
message("Amyloid-conversion longitudinal trends:")
print(amyloid_trends)
message("")
message("Time-by-group interaction models:")
print(interaction_models)
message("")
message("Amyloid classification used:")
message("  Reported positivity column: ", amy_positive_col)
message("  Visual-read column: ", amy_visual_col)
message("  Centiloid column: ", amy_centiloid_col)
message("  Centiloid automatic-positive threshold: ", amyloid_positive_centiloid_threshold)
message("  Centiloid borderline lower bound: ", amyloid_borderline_centiloid_lower)
message("  SUVR column: ", amy_suvr_col)
message("  SUVR automatic-positive threshold: ", amyloid_positive_suvr_threshold)
message("  SUVR borderline lower bound: ", amyloid_borderline_suvr_lower)
message("")
message("Outputs saved to: ", out_dir)
message("============================================================")