#!/usr/bin/env Rscript

# ==============================================================================
# BLSA AD EPOCH longitudinal trajectory analysis:
# Formal test of whether longitudinal slopes differ among
#   1) Remained CN
#   2) Later MCI
#   3) Later AD
#
# This script reproduces the cohort construction used in the original BLSA
# trajectory script, but removes all plotting and adds a formal group-by-time
# interaction analysis.
#
# PRIMARY INFERENCE
# -----------------
# Outcome:
#   Harmonized AD EPOCH acceleration (years)
#
# Time:
#   Years since selected CN diagnosis baseline
#
# Group:
#   Eventual diagnosis group (reference = Remained CN)
#
# Primary mixed model:
#   acceleration_years ~ years_since_cn_baseline * eventual_diagnosis_group
#   random intercept for participant
#
# The global group-by-time interaction is tested by an ML likelihood-ratio test
# comparing nested mixed models with and without the interaction.
#
# Group-specific slopes and pairwise differences in slopes are estimated from
# the full interaction model refit with REML. Pairwise slope-difference P values
# are reported both unadjusted and Holm-adjusted.
#
# SENSITIVITY ANALYSIS
# --------------------
# A random-intercept + random-slope model is attempted. If it converges, the
# corresponding global interaction test and slope contrasts are printed.
# This is useful because the scientific question concerns longitudinal slopes,
# but the small Later AD group may make this model unstable.
#
# NOTE ON SMALL LATER-AD GROUP
# ----------------------------
# The formal test can still be run with a small Later AD group (e.g., N = 6),
# but estimates involving that group may have wide confidence intervals and low
# power. Interpret non-significant AD contrasts cautiously.
#
# OUTPUT
# ------
# Printed tables only. No figures and no output files are generated.
# ==============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(nlme)
  library(emmeans)
})

# ==============================================================================
# 1. User settings
# ==============================================================================

istaging_file <- paste0(
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/adni_lepoch/",
  "external_5_studies_istaging.tsv"
)

harmonized_prediction_file <- paste0(
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/adni_lepoch/",
  "results_external_longitudinal_ad_epoch_harmonized/",
  "external_5_studies_adni_brain_mri_ad_epoch_harmonized_scan_level_predictions.tsv"
)

minimum_scans_per_subject <- 2L
baseline_time_tolerance_years <- 0.05

# Reference group for all model coefficients and pairwise comparisons.
group_levels <- c(
  "Remained CN",
  "Later MCI",
  "Later AD"
)

# ==============================================================================
# 2. Helpers
# ==============================================================================

detect_column <- function(df, preferred, regex, label, required = TRUE) {
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
    
    str_detect(upper, "(^|[^A-Z])MCI([^A-Z]|$)") ~ "MCI",
    str_detect(upper, "ALZHEIMER|DEMENTIA") ~ "AD",
    str_detect(upper, "COGNITIVELY NORMAL|COGNITIVE NORMAL") ~ "CN",
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

select_earliest_mapped_diagnosis <- function(df) {
  df |>
    mutate(
      date_missing = is.na(diagnosis_date),
      baseline_selection_source = case_when(
        !date_missing ~ "Earliest mapped diagnosis Date",
        date_missing & !is.na(diagnosis_age) ~ "Youngest mapped diagnosis Age",
        TRUE ~ "Mapped diagnosis without Date or Age"
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
    ungroup()
}

fmt_num <- function(x, digits = 4) {
  ifelse(
    is.na(x),
    NA_character_,
    formatC(x, format = "f", digits = digits)
  )
}

fmt_p <- function(x) {
  ifelse(
    is.na(x),
    NA_character_,
    ifelse(
      x < 0.001,
      formatC(x, format = "e", digits = 3),
      formatC(x, format = "f", digits = 4)
    )
  )
}

print_table <- function(x, title) {
  cat("\n")
  cat(paste0(strrep("=", 88), "\n"))
  cat(title, "\n")
  cat(paste0(strrep("=", 88), "\n"))
  print(x, n = Inf, width = Inf)
}

safe_lme <- function(fixed, random, data, method) {
  tryCatch(
    nlme::lme(
      fixed = fixed,
      random = random,
      data = data,
      method = method,
      na.action = na.omit,
      control = nlme::lmeControl(
        opt = "optim",
        msMaxIter = 500,
        msVerbose = FALSE,
        returnObject = TRUE
      )
    ),
    error = function(e) {
      attr(e, "fit_error") <- TRUE
      e
    }
  )
}

is_successful_fit <- function(x) {
  inherits(x, "lme")
}

extract_lrt <- function(reduced_fit, full_fit, model_label) {
  cmp <- anova(reduced_fit, full_fit)
  
  lratio_col <- intersect(c("L.Ratio", "L.Ratio."), names(cmp))
  p_col <- intersect(c("p-value", "p.value"), names(cmp))
  
  lratio <- if (length(lratio_col) > 0) {
    as.numeric(cmp[[lratio_col[[1]]]][2])
  } else {
    NA_real_
  }
  
  p_value <- if (length(p_col) > 0) {
    as.numeric(cmp[[p_col[[1]]]][2])
  } else {
    NA_real_
  }
  
  df_diff <- as.numeric(cmp$df[2] - cmp$df[1])
  
  tibble(
    model = model_label,
    reduced_fixed_effects = "time + group",
    full_fixed_effects = "time * group",
    df_difference = df_diff,
    likelihood_ratio_chisq = lratio,
    p_value = p_value,
    AIC_reduced = AIC(reduced_fit),
    AIC_full = AIC(full_fit),
    BIC_reduced = BIC(reduced_fit),
    BIC_full = BIC(full_fit),
    logLik_reduced = as.numeric(logLik(reduced_fit)),
    logLik_full = as.numeric(logLik(full_fit))
  )
}

extract_fixed_effects <- function(fit) {
  tt <- as.data.frame(summary(fit)$tTable)
  tt$term <- rownames(tt)
  rownames(tt) <- NULL
  
  out <- as_tibble(tt) |>
    transmute(
      term,
      estimate = Value,
      se = Std.Error,
      df = DF,
      t_value = `t-value`,
      p_value = `p-value`,
      ci_low = estimate - qt(0.975, df = df) * se,
      ci_high = estimate + qt(0.975, df = df) * se
    )
  
  out
}

extract_group_slopes <- function(fit, analysis_data) {
  # IMPORTANT: pass the actual analysis dataset explicitly to emtrends().
  # This is required for models fitted inside safe_lme(), where the stored
  # nlme call contains `data = data` and emmeans cannot otherwise reconstruct
  # the original data after the helper function returns.
  tr <- emmeans::emtrends(
    fit,
    specs = "eventual_diagnosis_group",
    var = "years_since_cn_baseline",
    data = analysis_data
  )
  
  s <- as.data.frame(
    summary(
      tr,
      infer = c(TRUE, TRUE),
      adjust = "none"
    )
  )
  
  trend_col <- grep("\\.trend$", names(s), value = TRUE)
  
  if (length(trend_col) != 1) {
    stop("Could not uniquely identify the emtrends slope column.")
  }
  
  as_tibble(s) |>
    transmute(
      eventual_diagnosis_group,
      slope_acceleration_years_per_year = .data[[trend_col]],
      se = SE,
      df = df,
      ci_low = lower.CL,
      ci_high = upper.CL,
      t_ratio = t.ratio,
      p_value_vs_zero = p.value
    )
}

extract_pairwise_slope_differences <- function(fit, analysis_data) {
  # Explicit data avoids emmeans/ref_grid reconstruction failures for lme
  # models fitted through the safe_lme() wrapper.
  tr <- emmeans::emtrends(
    fit,
    specs = "eventual_diagnosis_group",
    var = "years_since_cn_baseline",
    data = analysis_data
  )
  
  pw <- as.data.frame(
    summary(
      pairs(tr, adjust = "none"),
      infer = c(TRUE, TRUE),
      adjust = "none"
    )
  )
  
  as_tibble(pw) |>
    transmute(
      contrast,
      slope_difference = estimate,
      se = SE,
      df = df,
      ci_low = lower.CL,
      ci_high = upper.CL,
      t_ratio = t.ratio,
      p_value_raw = p.value
    ) |>
    mutate(
      p_value_holm = p.adjust(
        p_value_raw,
        method = "holm"
      )
    )
}

extract_global_wald_interaction <- function(fit, model_label) {
  a <- as.data.frame(
    anova(
      fit,
      type = "marginal"
    )
  )
  
  a$term <- rownames(a)
  rownames(a) <- NULL
  
  interaction_row <- a |>
    filter(
      str_detect(
        term,
        "years_since_cn_baseline:eventual_diagnosis_group|eventual_diagnosis_group:years_since_cn_baseline"
      )
    )
  
  if (nrow(interaction_row) != 1) {
    return(
      tibble(
        model = model_label,
        term = "time x group",
        num_df = NA_real_,
        den_df = NA_real_,
        F_value = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  f_col <- intersect(c("F-value", "F.value"), names(interaction_row))
  p_col <- intersect(c("p-value", "p.value"), names(interaction_row))
  num_col <- intersect(c("numDF", "num.df"), names(interaction_row))
  den_col <- intersect(c("denDF", "den.df"), names(interaction_row))
  
  tibble(
    model = model_label,
    term = "time x group",
    num_df = if (length(num_col) > 0) as.numeric(interaction_row[[num_col[[1]]]]) else NA_real_,
    den_df = if (length(den_col) > 0) as.numeric(interaction_row[[den_col[[1]]]]) else NA_real_,
    F_value = if (length(f_col) > 0) as.numeric(interaction_row[[f_col[[1]]]]) else NA_real_,
    p_value = if (length(p_col) > 0) as.numeric(interaction_row[[p_col[[1]]]]) else NA_real_
  )
}

# ==============================================================================
# 3. Validate inputs
# ==============================================================================

if (!file.exists(istaging_file)) {
  stop("iSTAGING file does not exist: ", istaging_file)
}

if (!file.exists(harmonized_prediction_file)) {
  stop(
    "Harmonized prediction file does not exist: ",
    harmonized_prediction_file
  )
}

# ==============================================================================
# 4. Read BLSA longitudinal diagnosis records
# ==============================================================================

istaging_df <- readr::read_tsv(
  istaging_file,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "unique"
)

study_col <- detect_column(
  istaging_df,
  preferred = c("Study", "STUDY"),
  regex = "(^|_)study$",
  label = "iSTAGING Study"
)

id_col <- detect_column(
  istaging_df,
  preferred = c("PTID", "participant_id", "IID"),
  regex = "(^ptid$|participant.*id|^iid$)",
  label = "iSTAGING participant ID"
)

dx_col <- detect_column(
  istaging_df,
  preferred = c("DX_Binary", "Dx_binary", "dx_binary"),
  regex = "(^|_)dx[_\\.]*binary$",
  label = "iSTAGING DX_Binary"
)

date_col <- detect_column(
  istaging_df,
  preferred = c("Date", "scan_date", "MRI_Date"),
  regex = "(^|_)date$",
  label = "iSTAGING Date",
  required = FALSE
)

age_col <- detect_column(
  istaging_df,
  preferred = c("Age", "AGE"),
  regex = "^age$",
  label = "iSTAGING Age"
)

blsa_diagnosis_long <- istaging_df |>
  transmute(
    participant_id = as.character(.data[[id_col]]),
    study = normalize_study(.data[[study_col]]),
    diagnosis_original = clean_optional_character(.data[[dx_col]]),
    diagnosis = normalize_diagnosis(.data[[dx_col]]),
    diagnosis_date = if (!is.na(date_col)) {
      parse_date_flexibly(.data[[date_col]])
    } else {
      as.Date(NA)
    },
    diagnosis_age = suppressWarnings(
      as.numeric(.data[[age_col]])
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

# ==============================================================================
# 5. Define baseline-CN participants and eventual diagnosis groups
# ==============================================================================

baseline_diagnosis <- blsa_diagnosis_long |>
  select_earliest_mapped_diagnosis()

baseline_cn_ids <- baseline_diagnosis |>
  filter(diagnosis == "CN") |>
  select(participant_id)

cn_baseline <- baseline_diagnosis |>
  filter(diagnosis == "CN") |>
  transmute(
    participant_id,
    baseline_diagnosis = diagnosis,
    diagnosis_date,
    diagnosis_age,
    baseline_selection_source
  )

eventual_group <- blsa_diagnosis_long |>
  semi_join(
    baseline_cn_ids,
    by = "participant_id"
  ) |>
  group_by(participant_id) |>
  summarise(
    ever_mci = any(diagnosis == "MCI", na.rm = TRUE),
    ever_ad = any(diagnosis == "AD", na.rm = TRUE),
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
      levels = group_levels
    )
  )

# ==============================================================================
# 6. Read harmonized BLSA AD EPOCH scans
# ==============================================================================

prediction_df <- readr::read_tsv(
  harmonized_prediction_file,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "unique"
)

prediction_study_col <- detect_column(
  prediction_df,
  preferred = c("external_Study", "Study", "STUDY"),
  regex = "(^|_)study$",
  label = "prediction Study"
)

prediction_id_col <- detect_column(
  prediction_df,
  preferred = c("PTID", "participant_id", "IID"),
  regex = "(^ptid$|participant.*id|^iid$)",
  label = "prediction participant ID"
)

prediction_date_col <- detect_column(
  prediction_df,
  preferred = c("Date", "scan_date", "MRI_Date"),
  regex = "(^|_)date$",
  label = "prediction Date",
  required = FALSE
)

prediction_age_col <- detect_column(
  prediction_df,
  preferred = c("Age", "age_at_scan_used_for_model", "AGE"),
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

blsa_epoch_scans <- prediction_df |>
  transmute(
    prediction_source_row = row_number(),
    participant_id = as.character(.data[[prediction_id_col]]),
    study = normalize_study(.data[[prediction_study_col]]),
    scan_date = if (!is.na(prediction_date_col)) {
      parse_date_flexibly(.data[[prediction_date_col]])
    } else {
      as.Date(NA)
    },
    scan_age = suppressWarnings(
      as.numeric(.data[[prediction_age_col]])
    ),
    acceleration_years = suppressWarnings(
      as.numeric(.data[[acceleration_years_col]])
    )
  ) |>
  filter(
    study == "BLSA",
    !is.na(participant_id),
    participant_id != "",
    !is.na(acceleration_years)
  )

# ==============================================================================
# 7. Align scans to the selected CN baseline
# ==============================================================================

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
      !is.na(scan_date) & !is.na(diagnosis_date),
      as.numeric(scan_date - diagnosis_date) / 365.25,
      NA_real_
    ),
    years_since_cn_baseline_by_age = ifelse(
      !is.na(scan_age) & !is.na(diagnosis_age),
      scan_age - diagnosis_age,
      NA_real_
    ),
    years_since_cn_baseline = coalesce(
      years_since_cn_baseline_by_date,
      years_since_cn_baseline_by_age
    )
  ) |>
  filter(
    is.finite(years_since_cn_baseline),
    years_since_cn_baseline >= -baseline_time_tolerance_years
  ) |>
  mutate(
    years_since_cn_baseline = ifelse(
      abs(years_since_cn_baseline) <= baseline_time_tolerance_years,
      0,
      years_since_cn_baseline
    )
  )

# Match each scan to the closest available diagnosis record, preserving the
# same approach as the original plotting script.
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
      !is.na(scan_date) & !is.na(diagnosis_date),
      abs(as.numeric(scan_date - diagnosis_date)),
      NA_real_
    ),
    diagnosis_age_difference_years = ifelse(
      !is.na(scan_age) & !is.na(diagnosis_age),
      abs(scan_age - diagnosis_age),
      NA_real_
    ),
    diagnosis_match_priority = case_when(
      is.finite(diagnosis_date_difference_days) ~ 1,
      is.finite(diagnosis_age_difference_years) ~ 2,
      !is.na(diagnosis_date) ~ 3,
      !is.na(diagnosis_age) ~ 4,
      TRUE ~ 5
    ),
    diagnosis_match_distance = case_when(
      diagnosis_match_priority == 1 ~ diagnosis_date_difference_days,
      diagnosis_match_priority == 2 ~ diagnosis_age_difference_years,
      diagnosis_match_priority == 3 ~ as.numeric(diagnosis_date),
      diagnosis_match_priority == 4 ~ diagnosis_age,
      TRUE ~ Inf
    )
  ) |>
  arrange(
    prediction_source_row,
    diagnosis_match_priority,
    diagnosis_match_distance,
    diagnosis_date,
    diagnosis_age
  ) |>
  group_by(prediction_source_row) |>
  slice_head(n = 1) |>
  ungroup() |>
  transmute(
    prediction_source_row,
    scan_diagnosis = diagnosis
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
    eventual_diagnosis_group = factor(
      eventual_diagnosis_group,
      levels = group_levels
    )
  ) |>
  arrange(
    participant_id,
    years_since_cn_baseline,
    scan_date,
    scan_age
  )

# Remove duplicated participant-time observations, as in the original script.
aligned_scans <- aligned_scans |>
  group_by(
    participant_id,
    years_since_cn_baseline
  ) |>
  slice_head(n = 1) |>
  ungroup()

# ==============================================================================
# 8. Restrict to participants with longitudinal information
# ==============================================================================

longitudinal_subjects <- aligned_scans |>
  group_by(participant_id) |>
  summarise(
    eventual_diagnosis_group = first(eventual_diagnosis_group),
    n_scans = n(),
    followup_span_years =
      max(years_since_cn_baseline, na.rm = TRUE) -
      min(years_since_cn_baseline, na.rm = TRUE),
    .groups = "drop"
  ) |>
  filter(
    n_scans >= minimum_scans_per_subject,
    followup_span_years > 0,
    !is.na(eventual_diagnosis_group)
  )

traj_df <- aligned_scans |>
  semi_join(
    longitudinal_subjects,
    by = "participant_id"
  ) |>
  filter(
    !is.na(eventual_diagnosis_group),
    is.finite(years_since_cn_baseline),
    is.finite(acceleration_years)
  ) |>
  mutate(
    participant_id = factor(participant_id),
    eventual_diagnosis_group = factor(
      eventual_diagnosis_group,
      levels = group_levels
    )
  ) |>
  arrange(
    participant_id,
    years_since_cn_baseline
  )

if (nrow(traj_df) == 0) {
  stop("No qualified longitudinal BLSA observations remain for modeling.")
}

missing_groups <- setdiff(
  group_levels,
  unique(as.character(traj_df$eventual_diagnosis_group))
)

if (length(missing_groups) > 0) {
  stop(
    "The interaction cannot be estimated because the following group(s) are absent: ",
    paste(missing_groups, collapse = ", ")
  )
}

# ==============================================================================
# 9. Descriptive cohort table
# ==============================================================================

group_summary <- traj_df |>
  group_by(eventual_diagnosis_group) |>
  summarise(
    n_subjects = n_distinct(participant_id),
    n_scans = n(),
    scans_per_subject_median = median(
      as.numeric(table(participant_id)),
      na.rm = TRUE
    ),
    time_min_years = min(years_since_cn_baseline, na.rm = TRUE),
    time_max_years = max(years_since_cn_baseline, na.rm = TRUE),
    .groups = "drop"
  )

subject_followup <- traj_df |>
  group_by(
    participant_id,
    eventual_diagnosis_group
  ) |>
  summarise(
    n_scans = n(),
    followup_span_years =
      max(years_since_cn_baseline) -
      min(years_since_cn_baseline),
    .groups = "drop"
  )

group_followup_summary <- subject_followup |>
  group_by(eventual_diagnosis_group) |>
  summarise(
    n_subjects = n(),
    scans_per_subject_median = median(n_scans),
    scans_per_subject_q1 = quantile(n_scans, 0.25),
    scans_per_subject_q3 = quantile(n_scans, 0.75),
    followup_years_mean = mean(followup_span_years),
    followup_years_sd = sd(followup_span_years),
    followup_years_median = median(followup_span_years),
    followup_years_q1 = quantile(followup_span_years, 0.25),
    followup_years_q3 = quantile(followup_span_years, 0.75),
    .groups = "drop"
  )

print_table(
  group_followup_summary,
  "TABLE 1. Longitudinal sample size and follow-up by eventual diagnosis group"
)

later_ad_n <- group_followup_summary |>
  filter(eventual_diagnosis_group == "Later AD") |>
  pull(n_subjects)

if (length(later_ad_n) == 1 && later_ad_n < 10) {
  cat("\nIMPORTANT SMALL-SAMPLE NOTE:\n")
  cat(
    "The Later AD group contains only N = ",
    later_ad_n,
    " participants. Interaction and pairwise slope tests involving Later AD ",
    "are estimable, but power is limited and confidence intervals may be wide.\n",
    sep = ""
  )
}

# ==============================================================================
# 10. PRIMARY MODEL: random intercept for participant
# ==============================================================================

formula_no_interaction <-
  acceleration_years ~
  years_since_cn_baseline +
  eventual_diagnosis_group

formula_interaction <-
  acceleration_years ~
  years_since_cn_baseline *
  eventual_diagnosis_group

# ML fits for formal nested-model likelihood-ratio test.
primary_reduced_ml <- nlme::lme(
  fixed = formula_no_interaction,
  random = ~ 1 | participant_id,
  data = traj_df,
  method = "ML",
  na.action = na.omit,
  control = nlme::lmeControl(
    opt = "optim",
    msMaxIter = 500,
    returnObject = TRUE
  )
)

primary_full_ml <- nlme::lme(
  fixed = formula_interaction,
  random = ~ 1 | participant_id,
  data = traj_df,
  method = "ML",
  na.action = na.omit,
  control = nlme::lmeControl(
    opt = "optim",
    msMaxIter = 500,
    returnObject = TRUE
  )
)

primary_lrt <- extract_lrt(
  primary_reduced_ml,
  primary_full_ml,
  model_label = "Primary: random-intercept mixed model"
)

print_table(
  primary_lrt,
  paste0(
    "TABLE 2. PRIMARY GLOBAL TEST: Does adding the group-by-time interaction ",
    "improve model fit?"
  )
)

# REML full model for coefficient estimates, group-specific slopes, and slope
# contrasts.
primary_full_reml <- nlme::lme(
  fixed = formula_interaction,
  random = ~ 1 | participant_id,
  data = traj_df,
  method = "REML",
  na.action = na.omit,
  control = nlme::lmeControl(
    opt = "optim",
    msMaxIter = 500,
    returnObject = TRUE
  )
)

primary_wald <- extract_global_wald_interaction(
  primary_full_reml,
  model_label = "Primary: random-intercept mixed model"
)

print_table(
  primary_wald,
  "TABLE 3. PRIMARY GLOBAL WALD F TEST for the group-by-time interaction"
)

primary_fixed <- extract_fixed_effects(
  primary_full_reml
)

print_table(
  primary_fixed,
  paste0(
    "TABLE 4. Fixed-effect coefficients from the full primary interaction model ",
    "(REML)"
  )
)

primary_slopes <- extract_group_slopes(
  primary_full_reml,
  analysis_data = traj_df
)

print_table(
  primary_slopes,
  paste0(
    "TABLE 5. Estimated longitudinal AD EPOCH slope within each eventual ",
    "diagnosis group"
  )
)

primary_pairwise <- extract_pairwise_slope_differences(
  primary_full_reml,
  analysis_data = traj_df
)

print_table(
  primary_pairwise,
  paste0(
    "TABLE 6. PRIMARY PAIRWISE TESTS: Differences in longitudinal slopes ",
    "between groups"
  )
)

# ==============================================================================
# 11. SENSITIVITY MODEL: participant-specific random intercept + random slope
# ==============================================================================
# NOTE: The sensitivity fits are created inside safe_lme(). Their stored model
# call therefore contains `data = data`. All emmeans::emtrends() calls below
# explicitly receive traj_df via analysis_data, preventing ref_grid() from
# failing to reconstruct years_since_cn_baseline/eventual_diagnosis_group.

cat("\n")
cat(strrep("=", 88), "\n", sep = "")
cat("SENSITIVITY ANALYSIS: random intercept + random slope\n")
cat(strrep("=", 88), "\n", sep = "")

random_slope_structure <-
  ~ years_since_cn_baseline | participant_id

sensitivity_reduced_ml <- safe_lme(
  fixed = formula_no_interaction,
  random = random_slope_structure,
  data = traj_df,
  method = "ML"
)

sensitivity_full_ml <- safe_lme(
  fixed = formula_interaction,
  random = random_slope_structure,
  data = traj_df,
  method = "ML"
)

if (
  is_successful_fit(sensitivity_reduced_ml) &&
  is_successful_fit(sensitivity_full_ml)
) {
  
  sensitivity_lrt <- extract_lrt(
    sensitivity_reduced_ml,
    sensitivity_full_ml,
    model_label = "Sensitivity: random-intercept + random-slope mixed model"
  )
  
  print_table(
    sensitivity_lrt,
    paste0(
      "TABLE 7. SENSITIVITY GLOBAL TEST: group-by-time interaction with ",
      "participant-specific random slopes"
    )
  )
  
  sensitivity_full_reml <- safe_lme(
    fixed = formula_interaction,
    random = random_slope_structure,
    data = traj_df,
    method = "REML"
  )
  
  if (is_successful_fit(sensitivity_full_reml)) {
    
    sensitivity_slopes <- extract_group_slopes(
      sensitivity_full_reml,
      analysis_data = traj_df
    )
    
    print_table(
      sensitivity_slopes,
      paste0(
        "TABLE 8. SENSITIVITY estimated group-specific slopes from the ",
        "random-slope model"
      )
    )
    
    sensitivity_pairwise <- extract_pairwise_slope_differences(
      sensitivity_full_reml,
      analysis_data = traj_df
    )
    
    print_table(
      sensitivity_pairwise,
      paste0(
        "TABLE 9. SENSITIVITY pairwise differences in slopes from the ",
        "random-slope model"
      )
    )
    
    cat("\nRandom-effects variance structure from sensitivity model:\n")
    print(nlme::VarCorr(sensitivity_full_reml))
    
  } else {
    
    cat(
      "\nThe random-slope REML model did not converge. ",
      "Only the ML global sensitivity test is available.\n"
    )
  }
  
} else {
  
  cat(
    "The random-intercept + random-slope sensitivity model could not be ",
    "fit reliably. This may occur with sparse longitudinal data or a very ",
    "small Later AD group. The random-intercept model above remains the ",
    "pre-specified primary analysis.\n",
    sep = ""
  )
  
  if (!is_successful_fit(sensitivity_reduced_ml)) {
    cat(
      "Reduced random-slope model error: ",
      conditionMessage(sensitivity_reduced_ml),
      "\n",
      sep = ""
    )
  }
  
  if (!is_successful_fit(sensitivity_full_ml)) {
    cat(
      "Full random-slope model error: ",
      conditionMessage(sensitivity_full_ml),
      "\n",
      sep = ""
    )
  }
}

# ==============================================================================
# 12. Compact interpretation guide printed after the tables
# ==============================================================================

cat("\n")
cat(strrep("=", 88), "\n", sep = "")
cat("INTERPRETATION GUIDE\n")
cat(strrep("=", 88), "\n", sep = "")
cat(
  paste0(
    "1. The primary answer to the collaborator's question is TABLE 2: the ML ",
    "likelihood-ratio test of the global group-by-time interaction.\n",
    "2. If TABLE 2 is significant, the longitudinal slopes are not all equal ",
    "across Remained CN, Later MCI, and Later AD.\n",
    "3. TABLE 5 reports the estimated slope within each group. A within-group ",
    "P value tests whether that group's slope differs from zero; it does NOT ",
    "test whether slopes differ between groups.\n",
    "4. TABLE 6 directly tests the pairwise slope differences. Use the ",
    "Holm-adjusted P values for the three pairwise comparisons.\n",
    "5. Because Later AD is very small, emphasize its effect estimate and 95% ",
    "CI and avoid interpreting a non-significant contrast as evidence that ",
    "the slopes are equivalent.\n",
    "6. TABLES 7-9 are sensitivity analyses using participant-specific random ",
    "slopes when that model converges.\n"
  )
)

cat("\nAnalysis complete.\n")