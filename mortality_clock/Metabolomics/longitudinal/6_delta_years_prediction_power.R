# ============================================================
# Delta acceleration-years survival analysis for longitudinal
# metabolomics mortality clocks
#
# Goal:
#   1) Read instance 0 and instance 1 prediction outputs
#   2) Extract acceleration years:
#        {organ}_metabolomics_mortality_clock_acceleration_years
#   3) Compute:
#        delta_accel_years_1_minus_0 =
#          acceleration_years_instance1 - acceleration_years_instance0
#   4) Define landmark survival outcome from instance 1:
#        Cases     = death after instance 1 sample date and before censoring
#        Non-cases = alive/censored after instance 1
#   5) Test whether delta acceleration years add mortality-prediction power
#      beyond chronological age and baseline acceleration years:
#
#        Surv(time_from_instance1_years, event_after_instance1) ~
#          clock_acceleration_years_0_0 +
#          delta_accel_years_1_minus_0 +
#          chronological_age_0_0
#
#   6) Plot HR and 95% CI for delta acceleration years across clocks
#
# Organs:
#   Endocrine, Digestive, Hepatic, Immune
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(survival)
  library(broom)
  library(scales)
})

# -----------------------------
# 1. Paths and settings
# -----------------------------

baseline_root <- "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock"

longitudinal_root <- "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/mortality_clock/longitudinal/metabolomics"

id_match_csv <- "/Users/hao/cubic-home/Dataset/UKBB_UMelbourne/UKB_UMelbourne_vs_Penn_match_key.csv"

organ_labels <- c("Endocrine", "Digestive", "Hepatic", "Immune")
organ_clean_vec <- tolower(organ_labels)
names(organ_clean_vec) <- organ_labels

admin_censor_date_default <- as.Date("2022-11-30")

main_outdir <- file.path(
  longitudinal_root,
  "metabolomics_delta_acceleration_years_landmark_survival_analysis"
)
dir.create(main_outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 2. Helper functions
# -----------------------------

fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-300) return("<1e-300")
  if (p < 0.001) return(format(p, scientific = TRUE, digits = 2))
  formatC(p, format = "f", digits = 3)
}

fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

title_case <- function(x) {
  tools::toTitleCase(gsub("_", " ", x))
}

safe_as_date <- function(x) {
  suppressWarnings(as.Date(as.character(x)))
}

safe_min_date <- function(x) {
  x <- as.Date(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) return(as.Date(NA))
  min(x)
}

normalize_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("\\.0$", "", x)
  x[x %in% c("", "NA", "NaN", "NULL", "null")] <- NA_character_
  x
}

get_cindex <- function(fit) {
  cc <- summary(fit)$concordance
  if ("C" %in% names(cc)) return(unname(cc[["C"]]))
  unname(cc[1])
}

get_cindex_se <- function(fit) {
  cc <- summary(fit)$concordance
  if ("se(C)" %in% names(cc)) return(unname(cc[["se(C)"]]))
  if (length(cc) >= 2) return(unname(cc[2]))
  NA_real_
}

# -----------------------------
# 3. ID harmonization
# -----------------------------

load_id_match_key <- function(id_match_csv) {
  if (!file.exists(id_match_csv)) {
    warning("ID match file does not exist: ", id_match_csv)
    return(NULL)
  }

  m <- fread(id_match_csv)

  if (!all(c("id", "id_upenn") %in% names(m))) {
    warning(
      "ID match file must contain columns id and id_upenn. Found: ",
      paste(names(m), collapse = ", ")
    )
    return(NULL)
  }

  m <- m %>%
    transmute(
      id = normalize_id(id),
      id_upenn = normalize_id(id_upenn)
    ) %>%
    filter(!is.na(id), !is.na(id_upenn)) %>%
    distinct()

  m
}

make_id_candidates <- function(raw_ids, id_map) {
  raw_ids <- normalize_id(raw_ids)
  candidates <- list(raw = raw_ids)

  if (!is.null(id_map) && nrow(id_map) > 0) {
    map_id_to_upenn <- setNames(id_map$id_upenn, id_map$id)
    map_upenn_to_id <- setNames(id_map$id, id_map$id_upenn)

    x1 <- unname(map_id_to_upenn[raw_ids])
    x1 <- ifelse(is.na(x1) | x1 == "", raw_ids, x1)

    x2 <- unname(map_upenn_to_id[raw_ids])
    x2 <- ifelse(is.na(x2) | x2 == "", raw_ids, x2)

    candidates$id_to_id_upenn <- x1
    candidates$id_upenn_to_id <- x2
  }

  candidates
}

harmonize_pair_ids <- function(dat0, dat1, id_match_csv, outdir, organ_clean) {
  id_map <- load_id_match_key(id_match_csv)

  dat0$participant_id_original <- normalize_id(dat0$participant_id)
  dat1$participant_id_original <- normalize_id(dat1$participant_id)

  cand0 <- make_id_candidates(dat0$participant_id_original, id_map)
  cand1 <- make_id_candidates(dat1$participant_id_original, id_map)

  diag_tbl <- expand.grid(
    baseline_id_transform = names(cand0),
    instance1_id_transform = names(cand1),
    stringsAsFactors = FALSE
  ) %>%
    rowwise() %>%
    mutate(
      n_baseline = n_distinct(cand0[[baseline_id_transform]]),
      n_instance1 = n_distinct(cand1[[instance1_id_transform]]),
      n_overlap = length(intersect(
        unique(cand0[[baseline_id_transform]]),
        unique(cand1[[instance1_id_transform]])
      ))
    ) %>%
    ungroup() %>%
    arrange(desc(n_overlap))

  fwrite(
    diag_tbl,
    file.path(outdir, paste0(organ_clean, "_id_harmonization_overlap_diagnostics.tsv")),
    sep = "\t"
  )

  best <- diag_tbl[1, ]

  message(
    "Best ID harmonization for ", organ_clean, ": baseline = ",
    best$baseline_id_transform,
    ", instance1 = ",
    best$instance1_id_transform,
    ", overlap = ",
    best$n_overlap
  )

  dat0$participant_id <- cand0[[best$baseline_id_transform]]
  dat1$participant_id <- cand1[[best$instance1_id_transform]]

  list(
    dat0 = dat0,
    dat1 = dat1,
    diagnostics = diag_tbl,
    best = best
  )
}

# -----------------------------
# 4. Standardize and read prediction files
# -----------------------------

standardize_metabolomics_clock_columns <- function(x, organ_clean) {
  risk_col <- paste0(organ_clean, "_metabolomics_mortality_risk_score")
  z_col <- paste0(organ_clean, "_metabolomics_mortality_clock_acceleration_z")
  years_col <- paste0(organ_clean, "_metabolomics_mortality_clock_acceleration_years")
  clock_age_col <- paste0(organ_clean, "_metabolomics_mortality_clock_age_years")

  if (!years_col %in% names(x)) {
    candidate_years <- grep(
      paste0(
        organ_clean,
        ".*metabolomics.*mortality.*clock.*acceleration.*years$|acceleration_years$"
      ),
      names(x),
      value = TRUE
    )

    if (length(candidate_years) > 0) {
      message("Using candidate acceleration-years column for ", organ_clean, ": ", candidate_years[1])
      x[[years_col]] <- suppressWarnings(as.numeric(x[[candidate_years[1]]]))
    }
  }

  if (!years_col %in% names(x)) {
    stop("Could not find acceleration-years column for ", organ_clean, ": ", years_col)
  }

  if (!risk_col %in% names(x)) {
    candidate_risk <- grep(
      paste0(
        organ_clean,
        ".*metabolomics.*mortality.*risk_score$|metabolomics_mortality_risk_score$"
      ),
      names(x),
      value = TRUE
    )
    if (length(candidate_risk) > 0) {
      x[[risk_col]] <- suppressWarnings(as.numeric(x[[candidate_risk[1]]]))
    } else {
      x[[risk_col]] <- NA_real_
    }
  }

  if (!z_col %in% names(x)) {
    candidate_z <- grep(
      paste0(
        organ_clean,
        ".*metabolomics.*mortality.*clock.*acceleration.*z$|acceleration_z$"
      ),
      names(x),
      value = TRUE
    )
    if (length(candidate_z) > 0) {
      x[[z_col]] <- suppressWarnings(as.numeric(x[[candidate_z[1]]]))
    } else {
      x[[z_col]] <- NA_real_
    }
  }

  if (!clock_age_col %in% names(x)) {
    x[[clock_age_col]] <- NA_real_
  }

  x[[risk_col]] <- suppressWarnings(as.numeric(x[[risk_col]]))
  x[[z_col]] <- suppressWarnings(as.numeric(x[[z_col]]))
  x[[years_col]] <- suppressWarnings(as.numeric(x[[years_col]]))
  x[[clock_age_col]] <- suppressWarnings(as.numeric(x[[clock_age_col]]))

  x
}

read_metabolomics_clock_instance <- function(file, organ_clean, instance_id, instance_label, visit_order) {
  if (!file.exists(file)) {
    stop("Input file does not exist: ", file)
  }

  x <- fread(file, fill = TRUE, check.names = FALSE)
  x <- as.data.frame(x, check.names = FALSE)

  if (!"participant_id" %in% names(x)) {
    stop("participant_id column not found in: ", file)
  }

  x <- standardize_metabolomics_clock_columns(x, organ_clean)

  risk_col <- paste0(organ_clean, "_metabolomics_mortality_risk_score")
  z_col <- paste0(organ_clean, "_metabolomics_mortality_clock_acceleration_z")
  years_col <- paste0(organ_clean, "_metabolomics_mortality_clock_acceleration_years")
  clock_age_col <- paste0(organ_clean, "_metabolomics_mortality_clock_age_years")

  optional_cols <- c(
    "sample_date",
    "death_date",
    "admin_censor_date",
    "end_date",
    "age_at_baseline",
    "age_at_imaging",
    "sex",
    "split",
    risk_col,
    z_col,
    years_col,
    clock_age_col
  )

  for (cc in optional_cols) {
    if (!cc %in% names(x)) {
      x[[cc]] <- NA
    }
  }

  keep_cols <- c(
    "participant_id",
    "sample_date",
    "death_date",
    "admin_censor_date",
    "end_date",
    "age_at_baseline",
    "age_at_imaging",
    "sex",
    "split",
    risk_col,
    z_col,
    years_col,
    clock_age_col
  )

  x <- x[, keep_cols, drop = FALSE]

  x <- x %>%
    mutate(
      participant_id = normalize_id(participant_id),
      application_instance = instance_id,
      instance_label = instance_label,
      visit_order = visit_order,

      sample_date = safe_as_date(sample_date),
      death_date = safe_as_date(death_date),
      admin_censor_date = safe_as_date(admin_censor_date),
      end_date = safe_as_date(end_date),

      admin_censor_date = case_when(
        !is.na(admin_censor_date) ~ admin_censor_date,
        is.na(admin_censor_date) & !is.na(end_date) ~ end_date,
        TRUE ~ admin_censor_date_default
      ),

      chronological_age = case_when(
        !is.na(suppressWarnings(as.numeric(age_at_imaging))) ~ suppressWarnings(as.numeric(age_at_imaging)),
        !is.na(suppressWarnings(as.numeric(age_at_baseline))) ~ suppressWarnings(as.numeric(age_at_baseline)),
        TRUE ~ NA_real_
      ),

      clock_risk_score = suppressWarnings(as.numeric(.data[[risk_col]])),
      clock_acceleration_z = suppressWarnings(as.numeric(.data[[z_col]])),
      clock_acceleration_years = suppressWarnings(as.numeric(.data[[years_col]])),
      clock_age_years = suppressWarnings(as.numeric(.data[[clock_age_col]]))
    ) %>%
    select(
      participant_id,
      application_instance,
      instance_label,
      visit_order,
      sample_date,
      death_date,
      admin_censor_date,
      end_date,
      chronological_age,
      age_at_baseline,
      age_at_imaging,
      sex,
      split,
      clock_risk_score,
      clock_acceleration_z,
      clock_acceleration_years,
      clock_age_years
    )

  x
}

# -----------------------------
# 5. Analysis per organ
# -----------------------------

run_one_organ <- function(organ_label) {
  organ_clean <- organ_clean_vec[[organ_label]]

  message("\n============================================================")
  message("Running delta acceleration-years survival analysis for: ", organ_label)
  message("Clean organ name: ", organ_clean)
  message("============================================================")

  organ_outdir <- file.path(main_outdir, organ_label)
  dir.create(organ_outdir, recursive = TRUE, showWarnings = FALSE)

  baseline_file <- file.path(
    baseline_root,
    paste0(organ_label, "_metabolomics_mortality_clock"),
    paste0(organ_clean, "_metabolomics_mortality_clock_predictions.tsv")
  )

  instance1_file <- file.path(
    longitudinal_root,
    organ_label,
    paste0(organ_clean, "_metabolomics_mortality_clock_apply_instance_1_0_predictions.tsv")
  )

  dat0 <- read_metabolomics_clock_instance(
    file = baseline_file,
    organ_clean = organ_clean,
    instance_id = "0_0",
    instance_label = "Instance 0 baseline",
    visit_order = 0
  )

  dat1 <- read_metabolomics_clock_instance(
    file = instance1_file,
    organ_clean = organ_clean,
    instance_id = "1_0",
    instance_label = "Instance 1 follow-up",
    visit_order = 1
  )

  id_harm <- harmonize_pair_ids(
    dat0 = dat0,
    dat1 = dat1,
    id_match_csv = id_match_csv,
    outdir = organ_outdir,
    organ_clean = organ_clean
  )

  dat0 <- id_harm$dat0
  dat1 <- id_harm$dat1

  dat_all <- bind_rows(dat0, dat1)

  common_ids <- dat_all %>%
    filter(!is.na(clock_acceleration_years)) %>%
    distinct(participant_id, application_instance) %>%
    count(participant_id, name = "n_instances") %>%
    filter(n_instances == 2) %>%
    pull(participant_id)

  cat("Organ:", organ_label, "\n")
  cat("N in baseline instance 0:", n_distinct(dat0$participant_id), "\n")
  cat("N in follow-up instance 1:", n_distinct(dat1$participant_id), "\n")
  cat("N common with acceleration years at both instances:", length(common_ids), "\n")

  if (length(common_ids) == 0) {
    warning("No common participants with acceleration years for organ: ", organ_label)
    return(NULL)
  }

  dat_common <- dat_all %>%
    filter(participant_id %in% common_ids) %>%
    filter(!is.na(clock_acceleration_years)) %>%
    arrange(participant_id, application_instance, visit_order) %>%
    group_by(participant_id, application_instance) %>%
    slice(1) %>%
    ungroup()

  # Landmark survival outcome from instance 1
  landmark_tbl <- dat_common %>%
    filter(application_instance == "1_0") %>%
    group_by(participant_id) %>%
    summarise(
      sample_date_instance1 = safe_min_date(sample_date),
      death_date_instance1 = safe_min_date(death_date),
      admin_censor_date_instance1 = safe_min_date(admin_censor_date),
      .groups = "drop"
    ) %>%
    mutate(
      admin_censor_date_instance1 = if_else(
        is.na(admin_censor_date_instance1),
        admin_censor_date_default,
        admin_censor_date_instance1
      ),

      event_after_instance1 = case_when(
        !is.na(death_date_instance1) &
          !is.na(sample_date_instance1) &
          !is.na(admin_censor_date_instance1) ~
          death_date_instance1 > sample_date_instance1 &
          death_date_instance1 <= admin_censor_date_instance1,

        !is.na(death_date_instance1) &
          !is.na(sample_date_instance1) &
          is.na(admin_censor_date_instance1) ~
          death_date_instance1 > sample_date_instance1,

        TRUE ~ FALSE
      ),

      end_date_from_instance1 = case_when(
        event_after_instance1 ~ death_date_instance1,
        !is.na(admin_censor_date_instance1) ~ admin_censor_date_instance1,
        TRUE ~ as.Date(NA)
      ),

      time_from_instance1_days = as.numeric(end_date_from_instance1 - sample_date_instance1),
      time_from_instance1_years = time_from_instance1_days / 365.25,

      case_status = case_when(
        is.na(sample_date_instance1) ~ NA_character_,
        event_after_instance1 ~ "Cases",
        TRUE ~ "Non-cases"
      )
    )

  fwrite(
    landmark_tbl,
    file.path(organ_outdir, paste0(organ_clean, "_landmark_survival_from_instance1.tsv")),
    sep = "\t"
  )

  dat_common <- dat_common %>%
    left_join(landmark_tbl, by = "participant_id")

  fwrite(
    dat_common,
    file.path(organ_outdir, paste0(organ_clean, "_longitudinal_acceleration_years_common_long.tsv")),
    sep = "\t"
  )

  # Wide table with baseline, follow-up, and delta acceleration years
  dat_wide <- dat_common %>%
    select(
      participant_id,
      application_instance,
      clock_acceleration_years,
      clock_acceleration_z,
      chronological_age,
      clock_age_years,
      sample_date_instance1,
      death_date_instance1,
      admin_censor_date_instance1,
      end_date_from_instance1,
      time_from_instance1_years,
      event_after_instance1,
      case_status
    ) %>%
    pivot_wider(
      names_from = application_instance,
      values_from = c(
        clock_acceleration_years,
        clock_acceleration_z,
        chronological_age,
        clock_age_years
      )
    ) %>%
    mutate(
      delta_accel_years_1_minus_0 =
        clock_acceleration_years_1_0 - clock_acceleration_years_0_0,

      delta_accel_z_1_minus_0 =
        clock_acceleration_z_1_0 - clock_acceleration_z_0_0,

      delta_chrono_age_1_minus_0 =
        chronological_age_1_0 - chronological_age_0_0,

      delta_clock_age_1_minus_0 =
        clock_age_years_1_0 - clock_age_years_0_0,

      event_after_instance1 = as.integer(event_after_instance1)
    )

  fwrite(
    dat_wide,
    file.path(organ_outdir, paste0(organ_clean, "_wide_delta_acceleration_years.tsv")),
    sep = "\t"
  )

  # Cox model data
  cox_df <- dat_wide %>%
    filter(
      !is.na(time_from_instance1_years),
      time_from_instance1_years > 0,
      !is.na(event_after_instance1),
      !is.na(clock_acceleration_years_0_0),
      !is.na(delta_accel_years_1_minus_0),
      !is.na(chronological_age_0_0)
    )

  n_model <- nrow(cox_df)
  n_events <- sum(cox_df$event_after_instance1 == 1, na.rm = TRUE)

  cat("N in Cox model:", n_model, "\n")
  cat("N events after instance 1:", n_events, "\n")

  if (n_model < 50 || n_events < 10) {
    warning("Insufficient sample size or events for Cox model in organ: ", organ_label)

    return(list(
      organ_label = organ_label,
      organ_clean = organ_clean,
      dat_wide = dat_wide,
      delta_hr_tbl = data.frame(
        organ_label = organ_label,
        organ_clean = organ_clean,
        n = n_model,
        n_events = n_events,
        hr_delta = NA_real_,
        hr_delta_lower95 = NA_real_,
        hr_delta_upper95 = NA_real_,
        p_delta = NA_real_,
        lrt_chisq = NA_real_,
        lrt_df = NA_real_,
        lrt_p = NA_real_,
        cindex_baseline = NA_real_,
        cindex_delta = NA_real_,
        delta_cindex = NA_real_
      )
    ))
  }

  # Baseline-only model
  cox_baseline <- survival::coxph(
    survival::Surv(time_from_instance1_years, event_after_instance1) ~
      clock_acceleration_years_0_0 +
      chronological_age_0_0,
    data = cox_df
  )

  # Delta model
  cox_delta <- survival::coxph(
    survival::Surv(time_from_instance1_years, event_after_instance1) ~
      clock_acceleration_years_0_0 +
      delta_accel_years_1_minus_0 +
      chronological_age_0_0,
    data = cox_df
  )

  # Coefficient tables
  coef_baseline_tbl <- broom::tidy(
    cox_baseline,
    exponentiate = TRUE,
    conf.int = TRUE
  ) %>%
    mutate(
      organ_label = organ_label,
      organ_clean = organ_clean,
      model = "baseline_only"
    ) %>%
    select(organ_label, organ_clean, model, everything())

  coef_delta_tbl <- broom::tidy(
    cox_delta,
    exponentiate = TRUE,
    conf.int = TRUE
  ) %>%
    mutate(
      organ_label = organ_label,
      organ_clean = organ_clean,
      model = "baseline_plus_delta"
    ) %>%
    select(organ_label, organ_clean, model, everything())

  coef_all_tbl <- bind_rows(coef_baseline_tbl, coef_delta_tbl)

  fwrite(
    coef_all_tbl,
    file.path(organ_outdir, paste0(organ_clean, "_cox_model_coefficients.tsv")),
    sep = "\t"
  )

  # Save model summaries as text
  sink(file.path(organ_outdir, paste0(organ_clean, "_cox_baseline_only_summary.txt")))
  print(summary(cox_baseline))
  sink()

  sink(file.path(organ_outdir, paste0(organ_clean, "_cox_baseline_plus_delta_summary.txt")))
  print(summary(cox_delta))
  sink()

  # Likelihood-ratio test for added value of delta
  lrt_tbl <- as.data.frame(anova(cox_baseline, cox_delta, test = "LRT"))

  fwrite(
    lrt_tbl,
    file.path(organ_outdir, paste0(organ_clean, "_cox_delta_added_value_lrt.tsv")),
    sep = "\t",
    row.names = TRUE
  )

  lrt_chisq <- if ("Chisq" %in% names(lrt_tbl)) lrt_tbl$Chisq[2] else NA_real_
  lrt_df <- if ("Df" %in% names(lrt_tbl)) lrt_tbl$Df[2] else NA_real_
  lrt_p <- if ("Pr(>|Chi|)" %in% names(lrt_tbl)) lrt_tbl$`Pr(>|Chi|)`[2] else NA_real_

  delta_term_tbl <- coef_delta_tbl %>%
    filter(term == "delta_accel_years_1_minus_0")

  baseline_term_tbl <- coef_delta_tbl %>%
    filter(term == "clock_acceleration_years_0_0")

  age_term_tbl <- coef_delta_tbl %>%
    filter(term == "chronological_age_0_0")

  delta_hr_tbl <- data.frame(
    organ_label = organ_label,
    organ_clean = organ_clean,
    n = n_model,
    n_events = n_events,
    event_rate = n_events / n_model,

    hr_baseline_accel = baseline_term_tbl$estimate[1],
    hr_baseline_accel_lower95 = baseline_term_tbl$conf.low[1],
    hr_baseline_accel_upper95 = baseline_term_tbl$conf.high[1],
    p_baseline_accel = baseline_term_tbl$p.value[1],

    hr_delta = delta_term_tbl$estimate[1],
    hr_delta_lower95 = delta_term_tbl$conf.low[1],
    hr_delta_upper95 = delta_term_tbl$conf.high[1],
    p_delta = delta_term_tbl$p.value[1],

    hr_chronological_age = age_term_tbl$estimate[1],
    hr_chronological_age_lower95 = age_term_tbl$conf.low[1],
    hr_chronological_age_upper95 = age_term_tbl$conf.high[1],
    p_chronological_age = age_term_tbl$p.value[1],

    lrt_chisq = lrt_chisq,
    lrt_df = lrt_df,
    lrt_p = lrt_p,

    cindex_baseline = get_cindex(cox_baseline),
    cindex_baseline_se = get_cindex_se(cox_baseline),
    cindex_delta = get_cindex(cox_delta),
    cindex_delta_se = get_cindex_se(cox_delta),
    delta_cindex = get_cindex(cox_delta) - get_cindex(cox_baseline)
  )

  fwrite(
    delta_hr_tbl,
    file.path(organ_outdir, paste0(organ_clean, "_delta_acceleration_years_hr_summary.tsv")),
    sep = "\t"
  )

  invisible(list(
    organ_label = organ_label,
    organ_clean = organ_clean,
    dat_wide = dat_wide,
    cox_df = cox_df,
    cox_baseline = cox_baseline,
    cox_delta = cox_delta,
    coef_all_tbl = coef_all_tbl,
    delta_hr_tbl = delta_hr_tbl,
    lrt_tbl = lrt_tbl
  ))
}

# -----------------------------
# 6. Run all four clocks
# -----------------------------

all_results <- list()

for (organ_label in organ_labels) {
  all_results[[organ_label]] <- run_one_organ(organ_label)
}

# -----------------------------
# 7. Combine and save results
# -----------------------------

delta_hr_all <- bind_rows(
  lapply(all_results, function(x) {
    if (is.null(x)) return(NULL)
    x$delta_hr_tbl
  })
)

coef_all <- bind_rows(
  lapply(all_results, function(x) {
    if (is.null(x)) return(NULL)
    x$coef_all_tbl
  })
)

fwrite(
  delta_hr_all,
  file.path(main_outdir, "all_metabolomics_delta_acceleration_years_hr_summary.tsv"),
  sep = "\t"
)

fwrite(
  coef_all,
  file.path(main_outdir, "all_metabolomics_delta_acceleration_years_cox_coefficients.tsv"),
  sep = "\t"
)

print(delta_hr_all)

# -----------------------------
# 8. Forest plot for delta HR
# -----------------------------

plot_tbl <- delta_hr_all %>%
  filter(!is.na(hr_delta), !is.na(hr_delta_lower95), !is.na(hr_delta_upper95)) %>%
  mutate(
    organ_label = factor(organ_label, levels = rev(organ_labels)),
    hr_ci_label = paste0(
      "HR = ", fmt_num(hr_delta, 2),
      " (", fmt_num(hr_delta_lower95, 2),
      "-", fmt_num(hr_delta_upper95, 2), ")"
    ),
    p_label = paste0("P = ", vapply(p_delta, fmt_p, character(1))),
    lrt_label = paste0("LRT P = ", vapply(lrt_p, fmt_p, character(1))),
    cindex_label = paste0(
      "\u0394C = ",
      ifelse(is.na(delta_cindex), "NA", fmt_num(delta_cindex, 3))
    ),
    annotation_label = paste0(hr_ci_label, "\n", p_label, "; ", lrt_label)
  )

if (nrow(plot_tbl) > 0) {
  x_min_data <- min(plot_tbl$hr_delta_lower95, na.rm = TRUE)
  x_max_data <- max(plot_tbl$hr_delta_upper95, na.rm = TRUE)

  x_min <- min(0.90, x_min_data * 0.90)
  x_max <- max(1.20, x_max_data * 1.80)
  x_text <- x_max_data * 1.12

  p_forest <- ggplot(
    plot_tbl,
    aes(y = organ_label)
  ) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.7
    ) +
    geom_segment(
      aes(
        x = hr_delta_lower95,
        xend = hr_delta_upper95,
        yend = organ_label
      ),
      linewidth = 1.0,
      color = "black"
    ) +
    geom_point(
      aes(x = hr_delta),
      size = 3.8,
      shape = 21,
      fill = "white",
      color = "black",
      stroke = 0.9
    ) +
    geom_text(
      aes(
        x = x_text,
        label = annotation_label
      ),
      hjust = 0,
      size = 3.4,
      lineheight = 0.95
    ) +
    scale_x_log10(
      limits = c(x_min, x_max),
      breaks = c(0.90, 1.00, 1.05, 1.10, 1.20, 1.40, 1.60, 2.00),
      labels = number_format(accuracy = 0.01)
    ) +
    labs(
      x = "Hazard ratio per 1-year increase in \u0394 acceleration years",
      y = NULL,
      title = "Incremental mortality association of longitudinal metabolomics clock change",
      subtitle = paste0(
        "Landmark Cox models start at instance 1 and adjust for ",
        "chronological age at baseline and baseline acceleration years"
      )
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11),
      axis.text.y = element_text(face = "bold", size = 12),
      axis.text.x = element_text(size = 10),
      plot.margin = margin(10, 80, 10, 10)
    ) +
    coord_cartesian(clip = "off")

  ggsave(
    file.path(main_outdir, "all_metabolomics_delta_acceleration_years_hr_forest_plot.pdf"),
    p_forest,
    width = 10.5,
    height = 4.8
  )

  ggsave(
    file.path(main_outdir, "all_metabolomics_delta_acceleration_years_hr_forest_plot.png"),
    p_forest,
    width = 10.5,
    height = 4.8,
    dpi = 300
  )

  ggsave(
    file.path(main_outdir, "all_metabolomics_delta_acceleration_years_hr_forest_plot.svg"),
    p_forest,
    width = 10.5,
    height = 4.8
  )
}

# -----------------------------
# 9. Print final output locations
# -----------------------------

cat("\n============================================================\n")
cat("Finished delta acceleration-years survival analysis.\n")
cat("Main output directory:\n", main_outdir, "\n\n")
cat("Main combined outputs:\n")
cat("  all_metabolomics_delta_acceleration_years_hr_summary.tsv\n")
cat("  all_metabolomics_delta_acceleration_years_cox_coefficients.tsv\n")
cat("  all_metabolomics_delta_acceleration_years_hr_forest_plot.pdf/png/svg\n")
cat("============================================================\n")