#!/usr/bin/env Rscript

# ==============================================================================
# BLSA AD EPOCH:
# Does an individualized longitudinal AD EPOCH rate of change predict subsequent
# CN -> MCI/AD conversion beyond baseline AD EPOCH, age, and sex?
# ==============================================================================
#
# RATIONALE
# ---------
# AD EPOCH was developed in ADNI and is evaluated here in the independent BLSA
# cohort. This script asks a distinct prognostic question:
#
#   Does a participant-specific longitudinal rate of change in AD EPOCH provide
#   information about subsequent CN -> MCI/AD conversion beyond the baseline
#   AD EPOCH value?
#
# PRIMARY DESIGN: FIXED LANDMARK / OBSERVATION WINDOW
# ---------------------------------------------------
# To ensure that the longitudinal slope is available BEFORE prediction begins:
#
#   1. Restrict to BLSA participants whose earliest mapped diagnosis is CN.
#   2. Define each participant's first qualified AD EPOCH MRI as time 0.
#   3. Use all qualified AD EPOCH scans within a fixed observation window
#      (default: 5 years after the first EPOCH MRI).
#   4. Require the participant to remain CN throughout that observation window.
#   5. Fit ONE longitudinal mixed-effects model across these repeated scans:
#
#        EPOCH_ij = beta0 + beta1*time_ij + b0_i + b1_i*time_ij + error_ij
#
#   6. Define each participant's individualized EPOCH rate of change as:
#
#        beta1 + b1_i
#
#      i.e., the population fixed slope plus the participant-specific random
#      slope (empirical-Bayes / BLUP estimate).
#
#   7. Begin conversion follow-up at the end of the fixed observation window.
#      Outcome:
#        event = first MCI or AD diagnosis after the landmark
#        censor = last mapped diagnostic follow-up while still CN
#
# PRIMARY COX COMPARISON
# ----------------------
# Reduced:
#   Surv(follow-up after landmark, conversion) ~
#       baseline EPOCH + baseline age + sex
#
# Full:
#   Surv(follow-up after landmark, conversion) ~
#       baseline EPOCH + individualized EPOCH slope + baseline age + sex
#
# The nested likelihood-ratio test is the primary answer to:
#
#   "Does longitudinal change in EPOCH provide prognostic information beyond
#    the baseline EPOCH value?"
#
# PREDICTIVE PERFORMANCE
# ----------------------
# Repeated stratified K-fold CV compares held-out Harrell C-indices for the
# reduced and full Cox models.
#
# Importantly, within each CV fold the longitudinal mixed model is re-fitted
# using TRAINING participants only. The fitted longitudinal model parameters
# are then used to obtain empirical-Bayes slopes for both training and held-out
# participants from their pre-landmark EPOCH measurements. This keeps the
# prediction evaluation fully out-of-fold.
#
# FIGURES
# -------
# A. Individualized mixed-model EPOCH slopes by eventual diagnosis group.
# B. Adjusted HRs for baseline EPOCH and individualized EPOCH slope.
# C. Held-out C-index: baseline model versus baseline + longitudinal slope.
#
# NOTES
# -----
# - Default observation window = 5 years.
# - Default minimum scans = 3, with at least 1 year of slope-estimation span.
# - If this leaves too few events, change `landmark_years` to 3 or 4 and/or
#   `minimum_scans_for_slope` to 2. Do not choose the window based on which one
#   gives the smallest P value; treat alternatives as sensitivity analyses.
# ==============================================================================


# ==============================================================================
# 1. User settings
# ==============================================================================

istaging_candidates <- c(
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/adni_lepoch/external_5_studies_istaging.tsv",
  "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch/external_5_studies_istaging.tsv"
)

prediction_candidates <- c(
  paste0(
    "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/adni_lepoch/",
    "results_external_longitudinal_ad_epoch_harmonized/",
    "external_5_studies_adni_brain_mri_ad_epoch_harmonized_scan_level_predictions.tsv"
  ),
  paste0(
    "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch/",
    "results_external_longitudinal_ad_epoch_harmonized/",
    "external_5_studies_adni_brain_mri_ad_epoch_harmonized_scan_level_predictions.tsv"
  )
)

out_dir_candidates <- c(
  paste0(
    "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/adni_lepoch/",
    "results_external_longitudinal_ad_epoch_comparison/",
    "BLSA_CN_conversion_mixed_model_individual_slopes"
  ),
  paste0(
    "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch/",
    "results_external_longitudinal_ad_epoch_comparison/",
    "BLSA_CN_conversion_mixed_model_individual_slopes"
  )
)

first_existing_file <- function(paths, label) {
  hit <- paths[file.exists(paths)]
  
  if (length(hit) == 0) {
    stop(
      "Could not find ", label, ". Checked:\n",
      paste(paths, collapse = "\n")
    )
  }
  
  hit[[1]]
}

istaging_file <- first_existing_file(
  istaging_candidates,
  "iSTAGING file"
)

harmonized_prediction_file <- first_existing_file(
  prediction_candidates,
  "harmonized AD EPOCH prediction file"
)

if (grepl("^/Users/", istaging_file)) {
  out_dir <- out_dir_candidates[[1]]
} else {
  out_dir <- out_dir_candidates[[2]]
}

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

prefix <- "BLSA_AD_EPOCH_mixed_slope_CN_conversion"

# Fixed observation window used to estimate each participant's longitudinal slope.
landmark_years <- 5

# Require enough scans and enough elapsed time to estimate a meaningful slope.
minimum_scans_for_slope <- 3L
minimum_slope_span_years <- 1.0

baseline_time_tolerance_years <- 0.05

# Repeated CV settings.
cv_folds <- 5L
cv_repeats <- 10L

# Paired participant-level bootstrap of out-of-fold Delta C-index.
bootstrap_replicates <- 1000L

random_seed <- 20260820L

group_levels <- c(
  "Remained CN",
  "Later MCI",
  "Later AD"
)

eventual_group_palette <- c(
  "Remained CN" = "#355C9A",
  "Later MCI"   = "#E3A018",
  "Later AD"    = "#B55239"
)

predictor_palette <- c(
  "Baseline AD EPOCH" = "#355C9A",
  "Individualized AD EPOCH slope" = "#B55239"
)

model_palette <- c(
  "Baseline EPOCH + age + sex" = "#355C9A",
  "Baseline + slope + age + sex" = "#B55239"
)

message("iSTAGING input: ", istaging_file)
message("Harmonized predictions: ", harmonized_prediction_file)
message("Output directory: ", out_dir)
message("Slope-estimation / landmark window: ", landmark_years, " years")
message("Minimum scans per participant: ", minimum_scans_for_slope)
message("Minimum slope-estimation span: ", minimum_slope_span_years, " years")


# ==============================================================================
# 2. Packages
# ==============================================================================

required_packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "ggplot2",
  "stringr",
  "scales",
  "tibble",
  "nlme",
  "survival",
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
  library(tibble)
  library(nlme)
  library(survival)
  library(patchwork)
})


# ==============================================================================
# 3. Helper functions
# ==============================================================================

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


normalize_study <- function(x) {
  x <- toupper(clean_optional_character(x))
  
  case_when(
    str_detect(x, "^BLSA") ~ "BLSA",
    TRUE ~ x
  )
}


normalize_diagnosis <- function(x) {
  upper <- toupper(clean_optional_character(x))
  
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


normalize_sex <- function(x) {
  upper <- toupper(clean_optional_character(x))
  
  case_when(
    upper %in% c(
      "F",
      "FEMALE",
      "WOMAN",
      "WOMEN"
    ) ~ "Female",
    
    upper %in% c(
      "M",
      "MALE",
      "MAN",
      "MEN"
    ) ~ "Male",
    
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
      output[indices[valid]] <- parsed[valid]
    }
  }
  
  output
}


first_nonmissing_character <- function(x) {
  x <- clean_optional_character(x)
  
  x <- x[!is.na(x)]
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  x[[1]]
}


max_date_or_na <- function(x) {
  if (all(is.na(x))) {
    return(as.Date(NA))
  }
  
  max(x, na.rm = TRUE)
}


max_num_or_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  max(x)
}


select_earliest_mapped_diagnosis <- function(df) {
  df |>
    mutate(
      date_missing = is.na(diagnosis_date)
    ) |>
    arrange(
      participant_id,
      date_missing,
      diagnosis_date,
      diagnosis_age
    ) |>
    group_by(participant_id) |>
    slice_head(n = 1) |>
    ungroup()
}


format_p <- function(p) {
  ifelse(
    is.na(p),
    "P = NA",
    ifelse(
      p < 0.001,
      paste0(
        "P = ",
        formatC(
          p,
          format = "e",
          digits = 2
        )
      ),
      paste0(
        "P = ",
        formatC(
          p,
          format = "f",
          digits = 3
        )
      )
    )
  )
}


print_table <- function(x, title) {
  cat("\n")
  cat(strrep("=", 104), "\n", sep = "")
  cat(title, "\n")
  cat(strrep("=", 104), "\n", sep = "")
  print(
    x,
    n = Inf,
    width = Inf
  )
}


varcorr_to_tibble <- function(fit) {
  vc <- nlme::VarCorr(fit)
  
  # VarCorr.lme is matrix-like but does not support as.data.frame() directly
  # in some nlme/R versions. Explicitly drop the S3 class and call the matrix
  # method so the code is portable across R 4.2.x installations.
  vc_df <- as.data.frame.matrix(
    unclass(vc),
    stringsAsFactors = FALSE
  )
  
  vc_df$term <- rownames(vc_df)
  rownames(vc_df) <- NULL
  
  vc_df <- vc_df |>
    tibble::as_tibble() |>
    dplyr::relocate(term)
  
  vc_df
}


standardization_parameters <- function(x) {
  x <- as.numeric(x)
  
  m <- mean(
    x,
    na.rm = TRUE
  )
  
  s <- sd(
    x,
    na.rm = TRUE
  )
  
  if (!is.finite(s) || s <= 0) {
    stop("Cannot standardize a variable with zero/non-finite SD.")
  }
  
  list(
    mean = m,
    sd = s
  )
}


apply_standardization <- function(
    x,
    pars
) {
  (
    as.numeric(x) -
      pars$mean
  ) /
    pars$sd
}


fit_random_slope_lme <- function(
    scan_data
) {
  correlated_fit <- tryCatch(
    nlme::lme(
      fixed =
        acceleration_years ~
        time_from_epoch_baseline,
      random =
        ~ time_from_epoch_baseline |
        participant_id,
      data = scan_data,
      method = "REML",
      na.action = na.omit,
      control = nlme::lmeControl(
        opt = "optim",
        msMaxIter = 500,
        msVerbose = FALSE,
        returnObject = TRUE
      )
    ),
    error = function(e) NULL
  )
  
  if (!is.null(correlated_fit)) {
    attr(
      correlated_fit,
      "random_structure_used"
    ) <- "correlated random intercept + random slope"
    
    return(correlated_fit)
  }
  
  message(
    "Correlated random-slope model failed; trying diagonal random-effects covariance."
  )
  
  diagonal_fit <- tryCatch(
    nlme::lme(
      fixed =
        acceleration_years ~
        time_from_epoch_baseline,
      random = list(
        participant_id =
          nlme::pdDiag(
            ~ time_from_epoch_baseline
          )
      ),
      data = scan_data,
      method = "REML",
      na.action = na.omit,
      control = nlme::lmeControl(
        opt = "optim",
        msMaxIter = 500,
        msVerbose = FALSE,
        returnObject = TRUE
      )
    ),
    error = function(e) NULL
  )
  
  if (is.null(diagonal_fit)) {
    stop(
      "Both correlated and diagonal random-slope mixed models failed."
    )
  }
  
  attr(
    diagonal_fit,
    "random_structure_used"
  ) <- "diagonal random intercept + random slope"
  
  diagonal_fit
}


extract_full_sample_blup_slopes <- function(
    fit
) {
  fixed <- nlme::fixef(fit)
  
  if (!("time_from_epoch_baseline" %in% names(fixed))) {
    stop(
      "Could not identify the fixed longitudinal time effect."
    )
  }
  
  re <- as.data.frame(
    nlme::ranef(fit)
  )
  
  re$participant_id <- rownames(re)
  rownames(re) <- NULL
  
  if (!("time_from_epoch_baseline" %in% names(re))) {
    stop(
      "Could not identify the random longitudinal slope."
    )
  }
  
  as_tibble(re) |>
    transmute(
      participant_id = as.character(
        participant_id
      ),
      fixed_population_slope =
        as.numeric(
          fixed[["time_from_epoch_baseline"]]
        ),
      random_slope =
        as.numeric(
          time_from_epoch_baseline
        ),
      individualized_slope_years_per_year =
        fixed_population_slope +
        random_slope
    )
}


# Empirical-Bayes random effects for arbitrary participants using a fitted lme.
# This allows held-out participants to receive individualized slopes from their
# own pre-landmark repeated EPOCH measurements without fitting the mixed model
# on their conversion outcomes.
estimate_eb_slopes_from_fitted_lme <- function(
    fit,
    scan_data
) {
  beta <- nlme::fixef(fit)
  
  required_beta <- c(
    "(Intercept)",
    "time_from_epoch_baseline"
  )
  
  if (!all(required_beta %in% names(beta))) {
    stop(
      "Unexpected fixed-effect names in longitudinal mixed model."
    )
  }
  
  D <- as.matrix(
    nlme::getVarCov(
      fit,
      type = "random.effects"
    )
  )
  
  if (!all(dim(D) == c(2, 2))) {
    stop(
      "Expected a 2 x 2 random-effects covariance matrix."
    )
  }
  
  sigma2 <- as.numeric(
    fit$sigma
  )^2
  
  scan_data |>
    group_by(participant_id) |>
    group_modify(
      function(.x, .y) {
        t_i <- as.numeric(
          .x$time_from_epoch_baseline
        )
        
        y_i <- as.numeric(
          .x$acceleration_years
        )
        
        X_i <- cbind(
          1,
          t_i
        )
        
        Z_i <- X_i
        
        beta_vec <- c(
          beta[["(Intercept)"]],
          beta[["time_from_epoch_baseline"]]
        )
        
        residual_i <-
          y_i -
          as.vector(
            X_i %*%
              beta_vec
          )
        
        V_i <-
          Z_i %*%
          D %*%
          t(Z_i) +
          diag(
            sigma2,
            nrow = length(t_i)
          )
        
        b_i <- as.vector(
          D %*%
            t(Z_i) %*%
            solve(
              V_i,
              residual_i
            )
        )
        
        tibble(
          n_scans_for_slope = length(t_i),
          slope_span_years =
            max(t_i) -
            min(t_i),
          fixed_population_slope =
            beta[[
              "time_from_epoch_baseline"
            ]],
          random_intercept = b_i[[1]],
          random_slope = b_i[[2]],
          individualized_slope_years_per_year =
            beta[[
              "time_from_epoch_baseline"
            ]] +
            b_i[[2]]
        )
      }
    ) |>
    ungroup()
}


cox_nested_lrt <- function(
    reduced_fit,
    full_fit,
    label
) {
  ll0 <- logLik(reduced_fit)
  ll1 <- logLik(full_fit)
  
  df0 <- attr(ll0, "df")
  df1 <- attr(ll1, "df")
  
  df_diff <- df1 - df0
  
  statistic <- 2 * (
    as.numeric(ll1) -
      as.numeric(ll0)
  )
  
  tibble(
    comparison = label,
    df_difference = df_diff,
    likelihood_ratio_chisq = statistic,
    p_value = pchisq(
      statistic,
      df = df_diff,
      lower.tail = FALSE
    ),
    AIC_reduced = AIC(reduced_fit),
    AIC_full = AIC(full_fit),
    delta_AIC_full_minus_reduced =
      AIC(full_fit) -
      AIC(reduced_fit),
    logLik_reduced = as.numeric(ll0),
    logLik_full = as.numeric(ll1)
  )
}


extract_cox_terms <- function(
    fit,
    term_labels
) {
  s <- summary(fit)
  
  cc <- as.data.frame(
    s$coefficients
  )
  
  ci <- as.data.frame(
    s$conf.int
  )
  
  cc$term <- rownames(cc)
  ci$term <- rownames(ci)
  
  rownames(cc) <- NULL
  rownames(ci) <- NULL
  
  cc |>
    as_tibble() |>
    select(
      term,
      coef,
      `se(coef)`,
      z,
      `Pr(>|z|)`
    ) |>
    left_join(
      ci |>
        as_tibble() |>
        select(
          term,
          `exp(coef)`,
          `lower .95`,
          `upper .95`
        ),
      by = "term"
    ) |>
    transmute(
      term,
      label = unname(
        term_labels[
          term
        ]
      ),
      beta = coef,
      se = `se(coef)`,
      z_value = z,
      HR = `exp(coef)`,
      ci_low = `lower .95`,
      ci_high = `upper .95`,
      p_value = `Pr(>|z|)`
    )
}


extract_mixed_fixed_effects <- function(
    fit
) {
  tt <- as.data.frame(
    summary(fit)$tTable
  )
  
  tt$term <- rownames(tt)
  rownames(tt) <- NULL
  
  as_tibble(tt) |>
    transmute(
      term,
      estimate = Value,
      se = Std.Error,
      df = DF,
      t_value = `t-value`,
      p_value = `p-value`
    )
}


harrell_c <- function(
    time,
    event,
    risk_score
) {
  keep <-
    is.finite(time) &
    !is.na(event) &
    is.finite(risk_score) &
    time > 0
  
  if (
    sum(
      event[keep] == 1
    ) == 0
  ) {
    return(NA_real_)
  }
  
  fit <- survival::concordance(
    survival::Surv(
      time[keep],
      event[keep]
    ) ~
      risk_score[keep],
    reverse = TRUE
  )
  
  as.numeric(
    fit$concordance
  )
}


make_stratified_folds <- function(
    event,
    k,
    seed
) {
  set.seed(seed)
  
  fold <- integer(
    length(event)
  )
  
  for (cls in sort(unique(event))) {
    ids <- which(
      event == cls
    )
    
    ids <- sample(
      ids,
      length(ids),
      replace = FALSE
    )
    
    fold[ids] <- rep(
      seq_len(k),
      length.out = length(ids)
    )
  }
  
  fold
}


prepare_train_test_features <- function(
    train_subjects,
    test_subjects,
    train_slopes,
    test_slopes
) {
  train <- train_subjects |>
    select(
      -any_of(
        c(
          "individualized_slope_years_per_year",
          "epoch_slope_z"
        )
      )
    ) |>
    left_join(
      train_slopes |>
        select(
          participant_id,
          individualized_slope_years_per_year
        ),
      by = "participant_id",
      relationship = "one-to-one"
    )
  
  test <- test_subjects |>
    select(
      -any_of(
        c(
          "individualized_slope_years_per_year",
          "epoch_slope_z"
        )
      )
    ) |>
    left_join(
      test_slopes |>
        select(
          participant_id,
          individualized_slope_years_per_year
        ),
      by = "participant_id",
      relationship = "one-to-one"
    )
  
  if (
    any(!is.finite(
      train$individualized_slope_years_per_year
    )) ||
    any(!is.finite(
      test$individualized_slope_years_per_year
    ))
  ) {
    stop(
      "Missing fold-specific individualized slopes after joining mixed-model results."
    )
  }
  
  pars_baseline <- standardization_parameters(
    train$baseline_epoch_years
  )
  
  pars_slope <- standardization_parameters(
    train$individualized_slope_years_per_year
  )
  
  pars_age <- standardization_parameters(
    train$baseline_scan_age
  )
  
  train <- train |>
    mutate(
      baseline_epoch_z =
        apply_standardization(
          baseline_epoch_years,
          pars_baseline
        ),
      epoch_slope_z =
        apply_standardization(
          individualized_slope_years_per_year,
          pars_slope
        ),
      baseline_age_z =
        apply_standardization(
          baseline_scan_age,
          pars_age
        ),
      sex = factor(
        sex,
        levels = c(
          "Female",
          "Male"
        )
      )
    )
  
  test <- test |>
    mutate(
      baseline_epoch_z =
        apply_standardization(
          baseline_epoch_years,
          pars_baseline
        ),
      epoch_slope_z =
        apply_standardization(
          individualized_slope_years_per_year,
          pars_slope
        ),
      baseline_age_z =
        apply_standardization(
          baseline_scan_age,
          pars_age
        ),
      sex = factor(
        sex,
        levels = c(
          "Female",
          "Male"
        )
      )
    )
  
  list(
    train = train,
    test = test
  )
}


repeated_cv_dynamic_cox <- function(
    subject_data,
    slope_scan_data,
    k = 5L,
    repeats = 10L,
    seed = 1L
) {
  n <- nrow(subject_data)
  
  event <- as.integer(
    subject_data$conversion_event
  )
  
  pred_reduced <- matrix(
    NA_real_,
    nrow = n,
    ncol = repeats
  )
  
  pred_full <- matrix(
    NA_real_,
    nrow = n,
    ncol = repeats
  )
  
  repeat_metrics <- vector(
    "list",
    repeats
  )
  
  reduced_formula <-
    survival::Surv(
      followup_years_after_landmark,
      conversion_event
    ) ~
    baseline_epoch_z +
    baseline_age_z +
    sex
  
  full_formula <-
    survival::Surv(
      followup_years_after_landmark,
      conversion_event
    ) ~
    baseline_epoch_z +
    epoch_slope_z +
    baseline_age_z +
    sex
  
  for (r in seq_len(repeats)) {
    message(
      "Cross-validation repeat ",
      r,
      " / ",
      repeats
    )
    
    folds <- make_stratified_folds(
      event = event,
      k = k,
      seed = seed + r - 1L
    )
    
    for (f in seq_len(k)) {
      test_idx <- which(
        folds == f
      )
      
      train_idx <- which(
        folds != f
      )
      
      train_subjects <- subject_data[
        train_idx,
        ,
        drop = FALSE
      ]
      
      test_subjects <- subject_data[
        test_idx,
        ,
        drop = FALSE
      ]
      
      train_ids <- as.character(
        train_subjects$participant_id
      )
      
      test_ids <- as.character(
        test_subjects$participant_id
      )
      
      train_scans <- slope_scan_data |>
        filter(
          participant_id %in%
            train_ids
        )
      
      test_scans <- slope_scan_data |>
        filter(
          participant_id %in%
            test_ids
        )
      
      longitudinal_fit <- fit_random_slope_lme(
        train_scans
      )
      
      train_slopes <- estimate_eb_slopes_from_fitted_lme(
        longitudinal_fit,
        train_scans
      )
      
      test_slopes <- estimate_eb_slopes_from_fitted_lme(
        longitudinal_fit,
        test_scans
      )
      
      fold_data <- prepare_train_test_features(
        train_subjects = train_subjects,
        test_subjects = test_subjects,
        train_slopes = train_slopes,
        test_slopes = test_slopes
      )
      
      train_df <- fold_data$train
      test_df <- fold_data$test
      
      reduced_fit <- survival::coxph(
        reduced_formula,
        data = train_df,
        ties = "efron",
        x = TRUE,
        model = TRUE
      )
      
      full_fit <- survival::coxph(
        full_formula,
        data = train_df,
        ties = "efron",
        x = TRUE,
        model = TRUE
      )
      
      pred_reduced[test_idx, r] <- as.numeric(
        predict(
          reduced_fit,
          newdata = test_df,
          type = "lp",
          reference = "zero"
        )
      )
      
      pred_full[test_idx, r] <- as.numeric(
        predict(
          full_fit,
          newdata = test_df,
          type = "lp",
          reference = "zero"
        )
      )
    }
    
    c0 <- harrell_c(
      time =
        subject_data$followup_years_after_landmark,
      event =
        subject_data$conversion_event,
      risk_score =
        pred_reduced[, r]
    )
    
    c1 <- harrell_c(
      time =
        subject_data$followup_years_after_landmark,
      event =
        subject_data$conversion_event,
      risk_score =
        pred_full[, r]
    )
    
    repeat_metrics[[r]] <- tibble(
      cv_repeat = r,
      C_index_baseline_model = c0,
      C_index_baseline_plus_slope = c1,
      delta_C_index = c1 - c0
    )
  }
  
  pred_reduced_mean <- rowMeans(
    pred_reduced,
    na.rm = TRUE
  )
  
  pred_full_mean <- rowMeans(
    pred_full,
    na.rm = TRUE
  )
  
  predictions <- tibble(
    participant_id =
      subject_data$participant_id,
    eventual_diagnosis_group =
      subject_data$eventual_diagnosis_group,
    conversion_event =
      subject_data$conversion_event,
    followup_years_after_landmark =
      subject_data$followup_years_after_landmark,
    OOF_risk_baseline_model =
      pred_reduced_mean,
    OOF_risk_baseline_plus_slope =
      pred_full_mean
  )
  
  overall_metrics <- tibble(
    model = c(
      "Baseline EPOCH + age + sex",
      "Baseline + slope + age + sex"
    ),
    repeated_OOF_C_index = c(
      harrell_c(
        subject_data$followup_years_after_landmark,
        subject_data$conversion_event,
        pred_reduced_mean
      ),
      harrell_c(
        subject_data$followup_years_after_landmark,
        subject_data$conversion_event,
        pred_full_mean
      )
    ),
    N = n,
    N_conversions = sum(
      subject_data$conversion_event == 1
    ),
    N_censored = sum(
      subject_data$conversion_event == 0
    ),
    cv_folds = k,
    cv_repeats = repeats
  )
  
  list(
    predictions = predictions,
    repeat_metrics = bind_rows(
      repeat_metrics
    ),
    overall_metrics = overall_metrics
  )
}


paired_bootstrap_delta_c <- function(
    time,
    event,
    risk_reduced,
    risk_full,
    B = 1000L,
    seed = 1L
) {
  observed_reduced <- harrell_c(
    time,
    event,
    risk_reduced
  )
  
  observed_full <- harrell_c(
    time,
    event,
    risk_full
  )
  
  observed_delta <-
    observed_full -
    observed_reduced
  
  set.seed(seed)
  
  n <- length(event)
  
  delta <- rep(
    NA_real_,
    B
  )
  
  for (b in seq_len(B)) {
    idx <- sample.int(
      n,
      size = n,
      replace = TRUE
    )
    
    if (
      sum(
        event[idx] == 1
      ) == 0
    ) {
      next
    }
    
    c0 <- harrell_c(
      time[idx],
      event[idx],
      risk_reduced[idx]
    )
    
    c1 <- harrell_c(
      time[idx],
      event[idx],
      risk_full[idx]
    )
    
    if (
      is.finite(c0) &&
      is.finite(c1)
    ) {
      delta[b] <- c1 - c0
    }
  }
  
  delta_ok <- delta[
    is.finite(delta)
  ]
  
  if (length(delta_ok) < 20) {
    stop(
      "Too few successful bootstrap replicates."
    )
  }
  
  ci <- quantile(
    delta_ok,
    probs = c(
      0.025,
      0.975
    ),
    na.rm = TRUE
  )
  
  p_left <- (
    sum(
      delta_ok <= 0
    ) + 1
  ) / (
    length(delta_ok) + 1
  )
  
  p_right <- (
    sum(
      delta_ok >= 0
    ) + 1
  ) / (
    length(delta_ok) + 1
  )
  
  p_two_sided <- min(
    1,
    2 * min(
      p_left,
      p_right
    )
  )
  
  tibble(
    C_index_reduced =
      observed_reduced,
    C_index_full =
      observed_full,
    delta_C_index_full_minus_reduced =
      observed_delta,
    bootstrap_CI_low =
      unname(ci[[1]]),
    bootstrap_CI_high =
      unname(ci[[2]]),
    bootstrap_p_two_sided =
      p_two_sided,
    successful_bootstrap_replicates =
      length(delta_ok)
  )
}


# ==============================================================================
# 4. Read iSTAGING diagnosis records
# ==============================================================================

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

sex_col_istaging <- detect_column(
  istaging_df,
  preferred = c(
    "Sex",
    "SEX",
    "Gender",
    "GENDER"
  ),
  regex = "(^|_)(sex|gender)$",
  label = "iSTAGING sex",
  required = FALSE
)

blsa_diagnosis_long <- istaging_df |>
  transmute(
    participant_id =
      as.character(
        .data[[id_col]]
      ),
    study =
      normalize_study(
        .data[[study_col]]
      ),
    diagnosis =
      normalize_diagnosis(
        .data[[dx_col]]
      ),
    diagnosis_date =
      if (!is.na(date_col)) {
        parse_date_flexibly(
          .data[[date_col]]
        )
      } else {
        as.Date(NA)
      },
    diagnosis_age =
      suppressWarnings(
        as.numeric(
          .data[[age_col]]
        )
      ),
    sex_from_istaging =
      if (!is.na(sex_col_istaging)) {
        normalize_sex(
          .data[[sex_col_istaging]]
        )
      } else {
        NA_character_
      }
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


# ==============================================================================
# 5. Define baseline-CN participants and diagnostic outcome summaries
# ==============================================================================

baseline_diagnosis <- blsa_diagnosis_long |>
  select_earliest_mapped_diagnosis()

baseline_cn_ids <- baseline_diagnosis |>
  filter(
    diagnosis == "CN"
  ) |>
  select(participant_id)

eventual_group <- blsa_diagnosis_long |>
  semi_join(
    baseline_cn_ids,
    by = "participant_id"
  ) |>
  group_by(participant_id) |>
  summarise(
    ever_mci =
      any(
        diagnosis == "MCI",
        na.rm = TRUE
      ),
    ever_ad =
      any(
        diagnosis == "AD",
        na.rm = TRUE
      ),
    eventual_diagnosis_group =
      case_when(
        ever_ad ~ "Later AD",
        ever_mci ~ "Later MCI",
        TRUE ~ "Remained CN"
      ),
    .groups = "drop"
  ) |>
  mutate(
    eventual_diagnosis_group =
      factor(
        eventual_diagnosis_group,
        levels = group_levels
      )
  )

first_conversion <- blsa_diagnosis_long |>
  semi_join(
    baseline_cn_ids,
    by = "participant_id"
  ) |>
  filter(
    diagnosis %in% c(
      "MCI",
      "AD"
    )
  ) |>
  mutate(
    date_missing =
      is.na(
        diagnosis_date
      )
  ) |>
  arrange(
    participant_id,
    date_missing,
    diagnosis_date,
    diagnosis_age
  ) |>
  group_by(participant_id) |>
  slice_head(n = 1) |>
  ungroup() |>
  transmute(
    participant_id,
    first_conversion_diagnosis =
      diagnosis,
    first_conversion_date =
      diagnosis_date,
    first_conversion_age =
      diagnosis_age
  )

last_diagnosis_followup <- blsa_diagnosis_long |>
  semi_join(
    baseline_cn_ids,
    by = "participant_id"
  ) |>
  group_by(participant_id) |>
  summarise(
    last_diagnosis_date =
      max_date_or_na(
        diagnosis_date
      ),
    last_diagnosis_age =
      max_num_or_na(
        diagnosis_age
      ),
    .groups = "drop"
  )

sex_from_istaging_tbl <- blsa_diagnosis_long |>
  semi_join(
    baseline_cn_ids,
    by = "participant_id"
  ) |>
  group_by(participant_id) |>
  summarise(
    sex_from_istaging =
      first_nonmissing_character(
        sex_from_istaging
      ),
    .groups = "drop"
  )


# ==============================================================================
# 6. Read harmonized BLSA AD EPOCH predictions
# ==============================================================================

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

sex_col_prediction <- detect_column(
  prediction_df,
  preferred = c(
    "Sex",
    "SEX",
    "Gender",
    "GENDER",
    "external_Sex",
    "external_SEX"
  ),
  regex = "(^|_)(sex|gender)$",
  label = "prediction sex",
  required = FALSE
)

blsa_epoch_scans <- prediction_df |>
  transmute(
    prediction_source_row =
      row_number(),
    participant_id =
      as.character(
        .data[[prediction_id_col]]
      ),
    study =
      normalize_study(
        .data[[prediction_study_col]]
      ),
    scan_date =
      if (!is.na(prediction_date_col)) {
        parse_date_flexibly(
          .data[[prediction_date_col]]
        )
      } else {
        as.Date(NA)
      },
    scan_age =
      suppressWarnings(
        as.numeric(
          .data[[prediction_age_col]]
        )
      ),
    acceleration_years =
      suppressWarnings(
        as.numeric(
          .data[[acceleration_years_col]]
        )
      ),
    sex_from_prediction =
      if (!is.na(sex_col_prediction)) {
        normalize_sex(
          .data[[sex_col_prediction]]
        )
      } else {
        NA_character_
      }
  ) |>
  filter(
    study == "BLSA",
    !is.na(participant_id),
    participant_id != "",
    is.finite(acceleration_years)
  )


# ==============================================================================
# 7. Resolve participant sex
# ==============================================================================

sex_from_prediction_tbl <- blsa_epoch_scans |>
  group_by(participant_id) |>
  summarise(
    sex_from_prediction =
      first_nonmissing_character(
        sex_from_prediction
      ),
    .groups = "drop"
  )

participant_sex <- full_join(
  sex_from_istaging_tbl,
  sex_from_prediction_tbl,
  by = "participant_id"
) |>
  mutate(
    sex =
      coalesce(
        sex_from_istaging,
        sex_from_prediction
      )
  ) |>
  select(
    participant_id,
    sex
  )

if (
  all(
    is.na(
      participant_sex$sex
    )
  )
) {
  stop(
    "Sex could not be identified from either iSTAGING or prediction data. ",
    "Please set the sex column manually in the script."
  )
}


# ==============================================================================
# 8. Identify each participant's first qualified EPOCH MRI
# ==============================================================================

first_epoch_scan <- blsa_epoch_scans |>
  semi_join(
    baseline_cn_ids,
    by = "participant_id"
  ) |>
  mutate(
    date_missing =
      is.na(
        scan_date
      )
  ) |>
  arrange(
    participant_id,
    date_missing,
    scan_date,
    scan_age,
    prediction_source_row
  ) |>
  group_by(participant_id) |>
  slice_head(n = 1) |>
  ungroup() |>
  transmute(
    participant_id,
    first_epoch_scan_date =
      scan_date,
    first_epoch_scan_age =
      scan_age
  )


# ==============================================================================
# 9. Express diagnosis and MRI times relative to first EPOCH MRI
# ==============================================================================

participant_timing <- first_epoch_scan |>
  left_join(
    first_conversion,
    by = "participant_id",
    relationship = "one-to-one"
  ) |>
  left_join(
    last_diagnosis_followup,
    by = "participant_id",
    relationship = "one-to-one"
  ) |>
  left_join(
    eventual_group,
    by = "participant_id",
    relationship = "one-to-one"
  ) |>
  left_join(
    participant_sex,
    by = "participant_id",
    relationship = "one-to-one"
  ) |>
  mutate(
    conversion_time_by_date =
      ifelse(
        !is.na(first_conversion_date) &
          !is.na(first_epoch_scan_date),
        as.numeric(
          first_conversion_date -
            first_epoch_scan_date
        ) / 365.25,
        NA_real_
      ),
    
    conversion_time_by_age =
      ifelse(
        is.finite(first_conversion_age) &
          is.finite(first_epoch_scan_age),
        first_conversion_age -
          first_epoch_scan_age,
        NA_real_
      ),
    
    conversion_time_from_epoch_baseline =
      coalesce(
        conversion_time_by_date,
        conversion_time_by_age
      ),
    
    last_followup_time_by_date =
      ifelse(
        !is.na(last_diagnosis_date) &
          !is.na(first_epoch_scan_date),
        as.numeric(
          last_diagnosis_date -
            first_epoch_scan_date
        ) / 365.25,
        NA_real_
      ),
    
    last_followup_time_by_age =
      ifelse(
        is.finite(last_diagnosis_age) &
          is.finite(first_epoch_scan_age),
        last_diagnosis_age -
          first_epoch_scan_age,
        NA_real_
      ),
    
    last_followup_time_from_epoch_baseline =
      coalesce(
        last_followup_time_by_date,
        last_followup_time_by_age
      )
  )

all_scans_aligned <- blsa_epoch_scans |>
  semi_join(
    baseline_cn_ids,
    by = "participant_id"
  ) |>
  inner_join(
    first_epoch_scan,
    by = "participant_id",
    relationship = "many-to-one"
  ) |>
  mutate(
    time_from_epoch_baseline_by_date =
      ifelse(
        !is.na(scan_date) &
          !is.na(first_epoch_scan_date),
        as.numeric(
          scan_date -
            first_epoch_scan_date
        ) / 365.25,
        NA_real_
      ),
    
    time_from_epoch_baseline_by_age =
      ifelse(
        is.finite(scan_age) &
          is.finite(first_epoch_scan_age),
        scan_age -
          first_epoch_scan_age,
        NA_real_
      ),
    
    time_from_epoch_baseline =
      coalesce(
        time_from_epoch_baseline_by_date,
        time_from_epoch_baseline_by_age
      )
  ) |>
  filter(
    is.finite(
      time_from_epoch_baseline
    ),
    time_from_epoch_baseline >=
      -baseline_time_tolerance_years
  ) |>
  mutate(
    time_from_epoch_baseline =
      ifelse(
        abs(
          time_from_epoch_baseline
        ) <=
          baseline_time_tolerance_years,
        0,
        time_from_epoch_baseline
      )
  ) |>
  arrange(
    participant_id,
    time_from_epoch_baseline,
    scan_date,
    scan_age,
    prediction_source_row
  ) |>
  group_by(
    participant_id,
    time_from_epoch_baseline
  ) |>
  slice_head(n = 1) |>
  ungroup()


# ==============================================================================
# 10. Fixed-landmark eligibility
# ==============================================================================

landmark_eligibility <- participant_timing |>
  mutate(
    conversion_before_or_at_landmark =
      is.finite(
        conversion_time_from_epoch_baseline
      ) &
      conversion_time_from_epoch_baseline <=
      landmark_years,
    
    known_at_risk_beyond_landmark =
      is.finite(
        last_followup_time_from_epoch_baseline
      ) &
      last_followup_time_from_epoch_baseline >
      landmark_years,
    
    eligible_at_landmark =
      !conversion_before_or_at_landmark &
      known_at_risk_beyond_landmark &
      !is.na(sex)
  )

# MRI scans used ONLY to estimate the slope. Every included participant is still
# CN through the landmark by construction.
candidate_slope_scans <- all_scans_aligned |>
  inner_join(
    landmark_eligibility |>
      filter(
        eligible_at_landmark
      ) |>
      select(
        participant_id,
        sex,
        eventual_diagnosis_group,
        conversion_time_from_epoch_baseline,
        last_followup_time_from_epoch_baseline
      ),
    by = "participant_id",
    relationship = "many-to-one"
  ) |>
  filter(
    time_from_epoch_baseline >= 0,
    time_from_epoch_baseline <=
      landmark_years
  ) |>
  arrange(
    participant_id,
    time_from_epoch_baseline
  )

slope_subject_qc <- candidate_slope_scans |>
  group_by(
    participant_id
  ) |>
  summarise(
    n_scans_for_slope = n(),
    slope_span_years =
      max(
        time_from_epoch_baseline
      ) -
      min(
        time_from_epoch_baseline
      ),
    .groups = "drop"
  ) |>
  mutate(
    passes_slope_qc =
      n_scans_for_slope >=
      minimum_scans_for_slope &
      slope_span_years >=
      minimum_slope_span_years
  )

analysis_ids <- slope_subject_qc |>
  filter(
    passes_slope_qc
  ) |>
  select(
    participant_id
  )

slope_scans <- candidate_slope_scans |>
  semi_join(
    analysis_ids,
    by = "participant_id"
  ) |>
  mutate(
    participant_id =
      factor(
        participant_id
      )
  ) |>
  arrange(
    participant_id,
    time_from_epoch_baseline
  )

if (
  n_distinct(
    slope_scans$participant_id
  ) < 20
) {
  stop(
    "Fewer than 20 participants remain after landmark/slope QC. ",
    "Consider a shorter landmark window or minimum_scans_for_slope = 2."
  )
}


# ==============================================================================
# 11. Subject-level outcome and baseline covariates
# ==============================================================================

subject_base <- slope_scans |>
  group_by(
    participant_id
  ) |>
  arrange(
    time_from_epoch_baseline,
    .by_group = TRUE
  ) |>
  summarise(
    sex =
      first(sex),
    
    eventual_diagnosis_group =
      first(
        eventual_diagnosis_group
      ),
    
    conversion_time_from_epoch_baseline =
      first(
        conversion_time_from_epoch_baseline
      ),
    
    last_followup_time_from_epoch_baseline =
      first(
        last_followup_time_from_epoch_baseline
      ),
    
    baseline_epoch_years =
      first(
        acceleration_years
      ),
    
    baseline_scan_age =
      first(
        scan_age
      ),
    
    n_scans_for_slope = n(),
    
    slope_span_years =
      max(
        time_from_epoch_baseline
      ) -
      min(
        time_from_epoch_baseline
      ),
    
    .groups = "drop"
  ) |>
  mutate(
    conversion_event =
      as.integer(
        is.finite(
          conversion_time_from_epoch_baseline
        ) &
          conversion_time_from_epoch_baseline >
          landmark_years
      ),
    
    followup_years_after_landmark =
      ifelse(
        conversion_event == 1,
        conversion_time_from_epoch_baseline -
          landmark_years,
        last_followup_time_from_epoch_baseline -
          landmark_years
      ),
    
    sex =
      factor(
        sex,
        levels = c(
          "Female",
          "Male"
        )
      ),
    
    eventual_diagnosis_group =
      factor(
        eventual_diagnosis_group,
        levels = group_levels
      )
  ) |>
  filter(
    is.finite(
      baseline_epoch_years
    ),
    is.finite(
      baseline_scan_age
    ),
    !is.na(sex),
    is.finite(
      followup_years_after_landmark
    ),
    followup_years_after_landmark > 0
  )

# Restrict slope scans to the exact subject-level analysis sample.
slope_scans <- slope_scans |>
  semi_join(
    subject_base |>
      select(
        participant_id
      ),
    by = "participant_id"
  )

if (
  sum(
    subject_base$conversion_event == 1
  ) < cv_folds
) {
  stop(
    "Too few post-landmark conversion events for the requested number of CV folds."
  )
}


# ==============================================================================
# 12. Fit full-sample longitudinal mixed model and derive BLUP slopes
# ==============================================================================

mixed_fit <- fit_random_slope_lme(
  slope_scans
)

message(
  "Longitudinal mixed-model random structure: ",
  attr(
    mixed_fit,
    "random_structure_used"
  )
)

mixed_fixed_effects <- extract_mixed_fixed_effects(
  mixed_fit
)

mixed_varcorr <- varcorr_to_tibble(mixed_fit)

full_sample_slopes <- extract_full_sample_blup_slopes(
  mixed_fit
)

analysis_df <- subject_base |>
  mutate(
    participant_id =
      as.character(
        participant_id
      )
  ) |>
  left_join(
    full_sample_slopes,
    by = "participant_id",
    relationship = "one-to-one"
  )

if (
  any(
    !is.finite(
      analysis_df$individualized_slope_years_per_year
    )
  )
) {
  stop(
    "Some participants are missing individualized mixed-model slopes."
  )
}

# Standardize continuous predictors in the full primary sample.
baseline_pars <- standardization_parameters(
  analysis_df$baseline_epoch_years
)

slope_pars <- standardization_parameters(
  analysis_df$individualized_slope_years_per_year
)

age_pars <- standardization_parameters(
  analysis_df$baseline_scan_age
)

analysis_df <- analysis_df |>
  mutate(
    baseline_epoch_z =
      apply_standardization(
        baseline_epoch_years,
        baseline_pars
      ),
    
    epoch_slope_z =
      apply_standardization(
        individualized_slope_years_per_year,
        slope_pars
      ),
    
    baseline_age_z =
      apply_standardization(
        baseline_scan_age,
        age_pars
      )
  )


# ==============================================================================
# 13. Descriptive tables
# ==============================================================================

cohort_summary <- analysis_df |>
  group_by(
    eventual_diagnosis_group
  ) |>
  summarise(
    n_subjects = n(),
    
    n_post_landmark_conversions =
      sum(
        conversion_event
      ),
    
    female_n =
      sum(
        sex == "Female"
      ),
    
    male_n =
      sum(
        sex == "Male"
      ),
    
    mean_baseline_age =
      mean(
        baseline_scan_age
      ),
    
    sd_baseline_age =
      sd(
        baseline_scan_age
      ),
    
    median_scans_for_slope =
      median(
        n_scans_for_slope
      ),
    
    mean_slope_span_years =
      mean(
        slope_span_years
      ),
    
    mean_baseline_EPOCH_years =
      mean(
        baseline_epoch_years
      ),
    
    sd_baseline_EPOCH_years =
      sd(
        baseline_epoch_years
      ),
    
    mean_individualized_slope =
      mean(
        individualized_slope_years_per_year
      ),
    
    sd_individualized_slope =
      sd(
        individualized_slope_years_per_year
      ),
    
    median_post_landmark_followup =
      median(
        followup_years_after_landmark
      ),
    
    .groups = "drop"
  )

sample_flow <- tibble(
  stage = c(
    "BLSA participants whose earliest mapped diagnosis is CN",
    paste0(
      "CN participants at risk and diagnostically observed beyond ",
      landmark_years,
      "-year landmark"
    ),
    paste0(
      "Participants with >= ",
      minimum_scans_for_slope,
      " EPOCH scans within landmark window"
    ),
    paste0(
      "Participants with slope span >= ",
      minimum_slope_span_years,
      " year(s)"
    ),
    "Final complete-case Cox analysis sample",
    "Post-landmark CN-to-MCI/AD conversions"
  ),
  
  N = c(
    nrow(
      baseline_cn_ids
    ),
    
    sum(
      landmark_eligibility$eligible_at_landmark,
      na.rm = TRUE
    ),
    
    sum(
      slope_subject_qc$n_scans_for_slope >=
        minimum_scans_for_slope
    ),
    
    sum(
      slope_subject_qc$passes_slope_qc
    ),
    
    nrow(
      analysis_df
    ),
    
    sum(
      analysis_df$conversion_event
    )
  )
)

print_table(
  sample_flow,
  "TABLE 1. Sample flow"
)

print_table(
  cohort_summary,
  "TABLE 2. Final landmark cohort by eventual diagnosis"
)

print_table(
  mixed_fixed_effects,
  "TABLE 3. Longitudinal mixed-model fixed effects"
)

cat("\nRandom-effects variance structure:\n")
print(
  nlme::VarCorr(
    mixed_fit
  )
)


# ==============================================================================
# 14. PRIMARY Cox models
# ==============================================================================

reduced_formula <-
  survival::Surv(
    followup_years_after_landmark,
    conversion_event
  ) ~
  baseline_epoch_z +
  baseline_age_z +
  sex

full_formula <-
  survival::Surv(
    followup_years_after_landmark,
    conversion_event
  ) ~
  baseline_epoch_z +
  epoch_slope_z +
  baseline_age_z +
  sex

reduced_fit <- survival::coxph(
  reduced_formula,
  data = analysis_df,
  ties = "efron",
  x = TRUE,
  model = TRUE
)

full_fit <- survival::coxph(
  full_formula,
  data = analysis_df,
  ties = "efron",
  x = TRUE,
  model = TRUE
)

primary_lrt <- cox_nested_lrt(
  reduced_fit,
  full_fit,
  label = paste0(
    "Baseline EPOCH + age + sex vs baseline EPOCH + ",
    "individualized longitudinal slope + age + sex"
  )
)

term_labels <- c(
  "baseline_epoch_z" =
    "Baseline AD EPOCH",
  
  "epoch_slope_z" =
    "Individualized AD EPOCH slope",
  
  "baseline_age_z" =
    "Baseline age",
  
  "sexMale" =
    "Male sex"
)

full_cox_terms <- extract_cox_terms(
  full_fit,
  term_labels = term_labels
)

ph_test <- as.data.frame(
  survival::cox.zph(
    full_fit,
    transform = "km"
  )$table
)

ph_test$term <- rownames(
  ph_test
)

rownames(
  ph_test
) <- NULL

ph_test <- as_tibble(
  ph_test
) |>
  rename(
    chisq = chisq,
    df = df,
    p_value = p
  )

print_table(
  primary_lrt,
  paste0(
    "TABLE 4. PRIMARY TEST: does individualized longitudinal EPOCH slope ",
    "add prognostic information beyond baseline EPOCH, age, and sex?"
  )
)

print_table(
  full_cox_terms,
  "TABLE 5. Full Cox model coefficients"
)

print_table(
  ph_test,
  "TABLE 6. Proportional-hazards diagnostics"
)


# ==============================================================================
# 15. Repeated out-of-fold prediction
# ==============================================================================

cv_results <- repeated_cv_dynamic_cox(
  subject_data = analysis_df,
  slope_scan_data = slope_scans |>
    mutate(
      participant_id =
        as.character(
          participant_id
        )
    ),
  k = cv_folds,
  repeats = cv_repeats,
  seed = random_seed
)

print_table(
  cv_results$overall_metrics,
  "TABLE 7. Repeated out-of-fold Harrell C-index"
)

print_table(
  cv_results$repeat_metrics,
  "TABLE 8. C-index by CV repeat"
)

bootstrap_delta_c <- paired_bootstrap_delta_c(
  time =
    cv_results$predictions$
    followup_years_after_landmark,
  
  event =
    cv_results$predictions$
    conversion_event,
  
  risk_reduced =
    cv_results$predictions$
    OOF_risk_baseline_model,
  
  risk_full =
    cv_results$predictions$
    OOF_risk_baseline_plus_slope,
  
  B = bootstrap_replicates,
  
  seed =
    random_seed +
    10000L
)

print_table(
  bootstrap_delta_c,
  paste0(
    "TABLE 9. Paired bootstrap test of incremental out-of-fold C-index ",
    "from individualized slope"
  )
)


# ==============================================================================
# 16. Plot A: individualized mixed-model slopes by eventual diagnosis
# ==============================================================================

group_n <- analysis_df |>
  count(
    eventual_diagnosis_group,
    name = "n"
  ) |>
  mutate(
    plot_label =
      paste0(
        as.character(
          eventual_diagnosis_group
        ),
        "\nN = ",
        n
      )
  )

group_labels <- setNames(
  group_n$plot_label,
  as.character(
    group_n$eventual_diagnosis_group
  )
)

p_slope <- ggplot(
  analysis_df,
  aes(
    x =
      eventual_diagnosis_group,
    y =
      individualized_slope_years_per_year,
    fill =
      eventual_diagnosis_group,
    color =
      eventual_diagnosis_group
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
    alpha = 0.23,
    trim = FALSE,
    linewidth = 0.45
  ) +
  geom_boxplot(
    width = 0.17,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.92,
    linewidth = 0.48
  ) +
  geom_jitter(
    width = 0.10,
    alpha = 0.50,
    size = 1.45
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 23,
    size = 3.0,
    fill = "#F7E8A4",
    color = "black"
  ) +
  scale_fill_manual(
    values =
      eventual_group_palette,
    guide = "none"
  ) +
  scale_color_manual(
    values =
      eventual_group_palette,
    guide = "none"
  ) +
  scale_x_discrete(
    name = NULL,
    labels = group_labels
  ) +
  scale_y_continuous(
    name = paste0(
      "Individualized AD EPOCH rate of change\n",
      "(acceleration years per year)"
    ),
    labels = number_format(
      accuracy = 0.05
    )
  ) +
  labs(
    title =
      "A  Individualized longitudinal AD EPOCH rates",
    
    subtitle =
      paste0(
        "BLUP slopes from repeated EPOCH scans during the first ",
        landmark_years,
        " years; all participants remained CN through the landmark"
      ),
    
    caption =
      paste0(
        "Each slope is the population mixed-model slope plus the ",
        "participant-specific random slope."
      )
  ) +
  theme_classic(
    base_size = 11.5
  ) +
  theme(
    plot.title =
      element_text(
        face = "bold",
        size = 13
      ),
    
    plot.subtitle =
      element_text(
        size = 9.8
      ),
    
    plot.caption =
      element_text(
        size = 8.3,
        hjust = 0
      ),
    
    axis.text =
      element_text(
        color = "black"
      ),
    
    axis.title =
      element_text(
        face = "bold"
      ),
    
    panel.grid.major.y =
      element_line(
        color = "grey90",
        linewidth = 0.35
      ),
    
    panel.grid.minor =
      element_blank()
  )


# ==============================================================================
# 17. Plot B: HR forest for baseline EPOCH and individualized slope
# ==============================================================================

forest_df <- full_cox_terms |>
  filter(
    term %in% c(
      "baseline_epoch_z",
      "epoch_slope_z"
    )
  ) |>
  mutate(
    label =
      factor(
        label,
        levels = rev(
          c(
            "Baseline AD EPOCH",
            "Individualized AD EPOCH slope"
          )
        )
      )
  )

primary_p <- primary_lrt$p_value[[1]]

p_forest <- ggplot(
  forest_df,
  aes(
    x = HR,
    y = label,
    color = label
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.50,
    color = "grey55"
  ) +
  geom_errorbarh(
    aes(
      xmin = ci_low,
      xmax = ci_high
    ),
    height = 0.13,
    linewidth = 0.75
  ) +
  geom_point(
    size = 3.2
  ) +
  scale_color_manual(
    values =
      predictor_palette,
    guide = "none"
  ) +
  scale_x_continuous(
    name =
      "Hazard ratio per 1-SD higher value",
    breaks =
      pretty_breaks(
        n = 5
      )
  ) +
  scale_y_discrete(
    name = NULL
  ) +
  labs(
    title =
      "B  Longitudinal change beyond baseline EPOCH",
    
    subtitle =
      paste0(
        "Nested Cox likelihood-ratio test for adding individualized slope: ",
        format_p(
          primary_p
        )
      ),
    
    caption =
      "Both effects are mutually adjusted; the full model also includes baseline age and sex."
  ) +
  theme_classic(
    base_size = 11.5
  ) +
  theme(
    plot.title =
      element_text(
        face = "bold",
        size = 13
      ),
    
    plot.subtitle =
      element_text(
        size = 9.8
      ),
    
    plot.caption =
      element_text(
        size = 8.3,
        hjust = 0
      ),
    
    axis.text =
      element_text(
        color = "black"
      ),
    
    axis.title =
      element_text(
        face = "bold"
      ),
    
    panel.grid.major.y =
      element_line(
        color = "grey90",
        linewidth = 0.35
      ),
    
    panel.grid.minor =
      element_blank()
  )


# ==============================================================================
# 18. Plot C: repeated out-of-fold C-index comparison
# ==============================================================================

cv_plot_df <- cv_results$repeat_metrics |>
  select(
    cv_repeat,
    C_index_baseline_model,
    C_index_baseline_plus_slope
  ) |>
  pivot_longer(
    cols = c(
      C_index_baseline_model,
      C_index_baseline_plus_slope
    ),
    names_to = "model_raw",
    values_to = "C_index"
  ) |>
  mutate(
    model =
      case_when(
        model_raw ==
          "C_index_baseline_model" ~
          "Baseline EPOCH + age + sex",
        
        model_raw ==
          "C_index_baseline_plus_slope" ~
          "Baseline + slope + age + sex",
        
        TRUE ~ model_raw
      ),
    
    model =
      factor(
        model,
        levels = c(
          "Baseline EPOCH + age + sex",
          "Baseline + slope + age + sex"
        )
      )
  )

delta_c <-
  bootstrap_delta_c$
  delta_C_index_full_minus_reduced[[1]]

delta_c_low <-
  bootstrap_delta_c$
  bootstrap_CI_low[[1]]

delta_c_high <-
  bootstrap_delta_c$
  bootstrap_CI_high[[1]]

delta_c_p <-
  bootstrap_delta_c$
  bootstrap_p_two_sided[[1]]

overall_c0 <-
  bootstrap_delta_c$
  C_index_reduced[[1]]

overall_c1 <-
  bootstrap_delta_c$
  C_index_full[[1]]

p_cindex <- ggplot(
  cv_plot_df,
  aes(
    x = model,
    y = C_index,
    group = cv_repeat
  )
) +
  geom_line(
    color = "grey75",
    alpha = 0.55,
    linewidth = 0.48
  ) +
  geom_point(
    aes(
      color = model
    ),
    alpha = 0.85,
    size = 2.35
  ) +
  stat_summary(
    aes(
      color = model
    ),
    fun = mean,
    geom = "point",
    shape = 18,
    size = 4.6
  ) +
  scale_color_manual(
    values = model_palette,
    guide = "none"
  ) +
  scale_x_discrete(
    name = NULL,
    labels = c(
      "Baseline EPOCH + age + sex" =
        paste0(
          "Baseline EPOCH + age + sex\nC = ",
          sprintf(
            "%.3f",
            overall_c0
          )
        ),
      
      "Baseline + slope + age + sex" =
        paste0(
          "Baseline + slope + age + sex\nC = ",
          sprintf(
            "%.3f",
            overall_c1
          )
        )
    )
  ) +
  scale_y_continuous(
    name =
      "Held-out Harrell C-index",
    labels =
      number_format(
        accuracy = 0.01
      ),
    breaks =
      pretty_breaks(
        n = 5
      )
  ) +
  labs(
    title =
      "C  Out-of-fold prediction of CN-to-MCI/AD conversion",
    
    subtitle =
      paste0(
        cv_folds,
        "-fold cross-validation x ",
        cv_repeats,
        " repeats; mixed model re-fitted within each training fold"
      ),
    
    caption =
      paste0(
        "\u0394C = ",
        sprintf(
          "%+.3f",
          delta_c
        ),
        " (95% CI ",
        sprintf(
          "%.3f",
          delta_c_low
        ),
        " to ",
        sprintf(
          "%.3f",
          delta_c_high
        ),
        "); ",
        format_p(
          delta_c_p
        ),
        ". Grey lines pair the two models within each CV repeat."
      )
  ) +
  theme_classic(
    base_size = 11.5
  ) +
  theme(
    plot.title =
      element_text(
        face = "bold",
        size = 13
      ),
    
    plot.subtitle =
      element_text(
        size = 9.5
      ),
    
    plot.caption =
      element_text(
        size = 8.3,
        hjust = 0
      ),
    
    axis.text =
      element_text(
        color = "black"
      ),
    
    axis.title =
      element_text(
        face = "bold"
      ),
    
    panel.grid.major.y =
      element_line(
        color = "grey90",
        linewidth = 0.35
      ),
    
    panel.grid.minor =
      element_blank()
  )


# ==============================================================================
# 19. Combined figure
# ==============================================================================

combined_plot <- (
  p_slope |
    p_forest |
    p_cindex
) +
  patchwork::plot_annotation(
    title =
      paste0(
        "Longitudinal AD EPOCH change predicts subsequent CN-to-MCI/AD conversion ",
        "beyond baseline EPOCH"
      ),
    
    subtitle =
      paste0(
        "Independent BLSA evaluation of the ADNI-trained AD EPOCH; ",
        landmark_years,
        "-year fixed landmark analysis"
      ),
    
    theme = theme(
      plot.title =
        element_text(
          face = "bold",
          size = 15
        ),
      
      plot.subtitle =
        element_text(
          size = 11
        )
    )
  )


# ==============================================================================
# 20. Save figures
# ==============================================================================

ggsave(
  file.path(
    out_dir,
    paste0(
      prefix,
      "_A_individualized_mixed_model_slopes.pdf"
    )
  ),
  p_slope,
  width = 7.4,
  height = 5.5
)

ggsave(
  file.path(
    out_dir,
    paste0(
      prefix,
      "_A_individualized_mixed_model_slopes.png"
    )
  ),
  p_slope,
  width = 7.4,
  height = 5.5,
  dpi = 350,
  bg = "white"
)

ggsave(
  file.path(
    out_dir,
    paste0(
      prefix,
      "_B_HR_baseline_and_slope.pdf"
    )
  ),
  p_forest,
  width = 7.3,
  height = 5.0
)

ggsave(
  file.path(
    out_dir,
    paste0(
      prefix,
      "_B_HR_baseline_and_slope.png"
    )
  ),
  p_forest,
  width = 7.3,
  height = 5.0,
  dpi = 350,
  bg = "white"
)

ggsave(
  file.path(
    out_dir,
    paste0(
      prefix,
      "_C_repeated_OOF_Cindex.pdf"
    )
  ),
  p_cindex,
  width = 7.3,
  height = 5.4
)

ggsave(
  file.path(
    out_dir,
    paste0(
      prefix,
      "_C_repeated_OOF_Cindex.png"
    )
  ),
  p_cindex,
  width = 7.3,
  height = 5.4,
  dpi = 350,
  bg = "white"
)

ggsave(
  file.path(
    out_dir,
    paste0(
      prefix,
      "_combined_ABC.pdf"
    )
  ),
  combined_plot,
  width = 19.8,
  height = 6.2
)

ggsave(
  file.path(
    out_dir,
    paste0(
      prefix,
      "_combined_ABC.png"
    )
  ),
  combined_plot,
  width = 19.8,
  height = 6.2,
  dpi = 350,
  bg = "white"
)


# ==============================================================================
# 21. Save analysis tables
# ==============================================================================

readr::write_tsv(
  sample_flow,
  file.path(
    out_dir,
    paste0(
      prefix,
      "_sample_flow.tsv"
    )
  )
)

readr::write_tsv(
  cohort_summary,
  file.path(
    out_dir,
    paste0(
      prefix,
      "_cohort_summary.tsv"
    )
  )
)

readr::write_tsv(
  slope_subject_qc,
  file.path(
    out_dir,
    paste0(
      prefix,
      "_slope_subject_QC.tsv"
    )
  )
)

readr::write_tsv(
  slope_scans |>
    mutate(
      participant_id =
        as.character(
          participant_id
        )
    ),
  file.path(
    out_dir,
    paste0(
      prefix,
      "_prelandmark_EPOCH_scans.tsv"
    )
  )
)

readr::write_tsv(
  analysis_df,
  file.path(
    out_dir,
    paste0(
      prefix,
      "_individualized_slopes_and_outcomes.tsv"
    )
  )
)

readr::write_tsv(
  mixed_fixed_effects,
  file.path(
    out_dir,
    paste0(
      prefix,
      "_mixed_model_fixed_effects.tsv"
    )
  )
)

# VarCorr contains mixed data types, so save after rowname conversion.
mixed_varcorr_out <- varcorr_to_tibble(mixed_fit)

readr::write_tsv(
  mixed_varcorr_out,
  file.path(
    out_dir,
    paste0(
      prefix,
      "_mixed_model_variance_components.tsv"
    )
  )
)

readr::write_tsv(
  primary_lrt,
  file.path(
    out_dir,
    paste0(
      prefix,
      "_PRIMARY_slope_beyond_baseline_LRT.tsv"
    )
  )
)

readr::write_tsv(
  full_cox_terms,
  file.path(
    out_dir,
    paste0(
      prefix,
      "_PRIMARY_full_Cox_coefficients.tsv"
    )
  )
)

readr::write_tsv(
  ph_test,
  file.path(
    out_dir,
    paste0(
      prefix,
      "_PRIMARY_Cox_PH_test.tsv"
    )
  )
)

readr::write_tsv(
  cv_results$overall_metrics,
  file.path(
    out_dir,
    paste0(
      prefix,
      "_repeated_OOF_Cindex.tsv"
    )
  )
)

readr::write_tsv(
  cv_results$repeat_metrics,
  file.path(
    out_dir,
    paste0(
      prefix,
      "_repeated_OOF_Cindex_by_repeat.tsv"
    )
  )
)

readr::write_tsv(
  cv_results$predictions,
  file.path(
    out_dir,
    paste0(
      prefix,
      "_repeated_OOF_predictions.tsv"
    )
  )
)

readr::write_tsv(
  bootstrap_delta_c,
  file.path(
    out_dir,
    paste0(
      prefix,
      "_bootstrap_delta_Cindex.tsv"
    )
  )
)


# ==============================================================================
# 22. Automated interpretation
# ==============================================================================

cat("\n")
cat(strrep("=", 104), "\n", sep = "")
cat("PRIMARY INTERPRETATION\n")
cat(strrep("=", 104), "\n", sep = "")

slope_row <- full_cox_terms |>
  filter(
    term == "epoch_slope_z"
  )

baseline_row <- full_cox_terms |>
  filter(
    term == "baseline_epoch_z"
  )

cat(
  paste0(
    "Landmark design: individualized slopes estimated from repeated AD EPOCH scans ",
    "during the first ",
    landmark_years,
    " years after the first EPOCH MRI; only participants who remained CN through ",
    "this landmark were included.\n"
  )
)

cat(
  paste0(
    "Final sample: N = ",
    nrow(
      analysis_df
    ),
    "; post-landmark CN-to-MCI/AD conversions = ",
    sum(
      analysis_df$conversion_event
    ),
    "; censored while CN = ",
    sum(
      analysis_df$conversion_event == 0
    ),
    ".\n"
  )
)

cat(
  paste0(
    "Baseline EPOCH: HR per 1 SD = ",
    sprintf(
      "%.3f",
      baseline_row$HR[[1]]
    ),
    " (95% CI ",
    sprintf(
      "%.3f",
      baseline_row$ci_low[[1]]
    ),
    " to ",
    sprintf(
      "%.3f",
      baseline_row$ci_high[[1]]
    ),
    "); ",
    format_p(
      baseline_row$p_value[[1]]
    ),
    ".\n"
  )
)

cat(
  paste0(
    "Individualized longitudinal slope: HR per 1 SD = ",
    sprintf(
      "%.3f",
      slope_row$HR[[1]]
    ),
    " (95% CI ",
    sprintf(
      "%.3f",
      slope_row$ci_low[[1]]
    ),
    " to ",
    sprintf(
      "%.3f",
      slope_row$ci_high[[1]]
    ),
    "); ",
    format_p(
      slope_row$p_value[[1]]
    ),
    ".\n"
  )
)

if (
  is.finite(
    primary_p
  ) &&
  primary_p < 0.05
) {
  cat(
    paste0(
      "YES: adding the individualized longitudinal slope significantly improved ",
      "the Cox model beyond baseline EPOCH, age, and sex (LRT ",
      format_p(
        primary_p
      ),
      ").\n"
    )
  )
} else {
  cat(
    paste0(
      "No statistically significant incremental model-fit evidence was detected ",
      "for the individualized longitudinal slope beyond baseline EPOCH, age, and sex ",
      "(LRT ",
      format_p(
        primary_p
      ),
      ").\n"
    )
  )
}

cat(
  paste0(
    "Out-of-fold discrimination: baseline model C = ",
    sprintf(
      "%.3f",
      overall_c0
    ),
    "; baseline + slope model C = ",
    sprintf(
      "%.3f",
      overall_c1
    ),
    "; Delta C = ",
    sprintf(
      "%+.3f",
      delta_c
    ),
    " (95% CI ",
    sprintf(
      "%.3f",
      delta_c_low
    ),
    " to ",
    sprintf(
      "%.3f",
      delta_c_high
    ),
    "); ",
    format_p(
      delta_c_p
    ),
    ".\n"
  )
)

cat(
  paste0(
    "These two results answer related but distinct questions: the LRT tests ",
    "incremental model information, whereas Delta C-index tests improvement in ",
    "held-out discrimination.\n"
  )
)

cat(
  "\nOutputs saved to: ",
  out_dir,
  "\n",
  sep = ""
)

cat(strrep("=", 104), "\n", sep = "")

# ========================================================================================================
# PRIMARY INTERPRETATION regarding output
# ========================================================================================================
# Landmark design: individualized slopes estimated from repeated AD EPOCH scans during the first 5 years after the first EPOCH MRI; only participants who remained CN through this landmark were included.
# Final sample: N = 330; post-landmark CN-to-MCI/AD conversions = 72; censored while CN = 258.
# Baseline EPOCH: HR per 1 SD = 0.989 (95% CI 0.734 to 1.334); P = 0.944.
# Individualized longitudinal slope: HR per 1 SD = 1.199 (95% CI 0.908 to 1.584); P = 0.201.
# No statistically significant incremental model-fit evidence was detected for the individualized longitudinal slope beyond baseline EPOCH, age, and sex (LRT P = 0.201).
# Out-of-fold discrimination: baseline model C = 0.769; baseline + slope model C = 0.769; Delta C = +0.001 (95% CI -0.010 to 0.010); P = 0.921.
# These two results answer related but distinct questions: the LRT tests incremental model information, whereas Delta C-index tests improvement in held-out discrimination.
# 
# Outputs saved to: /Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/adni_lepoch/results_external_longitudinal_ad_epoch_comparison/BLSA_CN_conversion_mixed_model_individual_slopes
# ========================================================================================================