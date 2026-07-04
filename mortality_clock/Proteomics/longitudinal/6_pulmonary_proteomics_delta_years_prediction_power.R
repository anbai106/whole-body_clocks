# ============================================================
# 6_pulmonary_proteomics_delta_years_prediction_power.R
#
# Landmark survival analysis for longitudinal change in the
# Pulmonary proteomics mortality-clock acceleration-years metric.
#
# Instances available for pulmonary proteomics:
#   0_0 = baseline/model instance
#   2_0 = follow-up proteomics instance 2
#   3_0 = follow-up proteomics instance 3
#
# This script tests whether delta acceleration-years adds
# mortality-prediction information beyond chronological age and
# the acceleration-years value at the start of the interval.
#
# Tested intervals:
#   2_0 - 0_0, landmark at instance 2_0
#   3_0 - 0_0, landmark at instance 3_0
#   3_0 - 2_0, landmark at instance 3_0
#
# Cox model for each interval:
#   Surv(time_from_endpoint_years, event_after_endpoint) ~
#     acceleration_years_start + delta_accel_years + chronological_age_start
#
# The forest plot displays HR and 95% CI for delta_accel_years.
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

longitudinal_root <- "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/mortality_clock/longitudinal/proteomics"

id_match_csv <- "/Users/hao/cubic-home/Dataset/UKBB_UMelbourne/UKB_UMelbourne_vs_Penn_match_key.csv"

organ_label <- "Pulmonary"
organ_clean <- "pulmonary"

admin_censor_date_default <- as.Date("2022-11-30")

main_outdir <- file.path(
  longitudinal_root,
  "pulmonary_proteomics_delta_acceleration_years_landmark_survival_analysis"
)
dir.create(main_outdir, recursive = TRUE, showWarnings = FALSE)

comparison_tbl <- tibble::tibble(
  comparison_id = c("2_minus_0", "3_minus_0", "3_minus_2"),
  comparison_label = c(
    "Instance 2 - baseline",
    "Instance 3 - baseline",
    "Instance 3 - instance 2"
  ),
  start_instance = c("0_0", "0_0", "2_0"),
  end_instance = c("2_0", "3_0", "3_0"),
  start_label = c("Instance 0", "Instance 0", "Instance 2"),
  end_label = c("Instance 2", "Instance 3", "Instance 3")
)

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

safe_as_date <- function(x) {
  suppressWarnings(as.Date(as.character(x)))
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
    warning("ID match file must contain columns id and id_upenn. Found: ", paste(names(m), collapse = ", "))
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

canonicalize_ids_to_upenn <- function(x, id_map) {
  x <- normalize_id(x)
  if (is.null(id_map) || nrow(id_map) == 0) return(x)

  map_id_to_upenn <- setNames(id_map$id_upenn, id_map$id)
  y <- unname(map_id_to_upenn[x])
  y <- ifelse(is.na(y) | y == "", x, y)
  normalize_id(y)
}

# -----------------------------
# 4. Read prediction files
# -----------------------------

standardize_pulmonary_clock_columns <- function(x) {
  risk_col <- paste0(organ_clean, "_proteomics_mortality_risk_score")
  z_col <- paste0(organ_clean, "_proteomics_mortality_clock_acceleration_z")
  years_col <- paste0(organ_clean, "_proteomics_mortality_clock_acceleration_years")
  age_col <- paste0(organ_clean, "_proteomics_mortality_clock_age_years")

  if (!risk_col %in% names(x)) {
    candidate_risk <- grep(
      "pulmonary.*proteomics.*mortality.*risk_score$|proteomics_mortality_risk_score$",
      names(x),
      value = TRUE
    )
    if (length(candidate_risk) > 0) {
      message("Using candidate risk-score column: ", candidate_risk[1])
      x[[risk_col]] <- suppressWarnings(as.numeric(x[[candidate_risk[1]]]))
    }
  }

  if (!years_col %in% names(x)) {
    candidate_years <- grep(
      "pulmonary.*proteomics.*mortality.*clock.*acceleration.*years$|acceleration_years$",
      names(x),
      value = TRUE
    )
    if (length(candidate_years) > 0) {
      message("Using candidate acceleration-years column: ", candidate_years[1])
      x[[years_col]] <- suppressWarnings(as.numeric(x[[candidate_years[1]]]))
    }
  }

  if (!z_col %in% names(x)) {
    candidate_z <- grep(
      "pulmonary.*proteomics.*mortality.*clock.*acceleration.*z$|acceleration_z$",
      names(x),
      value = TRUE
    )
    if (length(candidate_z) > 0) {
      x[[z_col]] <- suppressWarnings(as.numeric(x[[candidate_z[1]]]))
    } else {
      x[[z_col]] <- NA_real_
    }
  }

  if (!years_col %in% names(x)) {
    stop("Could not find acceleration-years column: ", years_col)
  }

  x[[risk_col]] <- if (risk_col %in% names(x)) suppressWarnings(as.numeric(x[[risk_col]])) else NA_real_
  x[[z_col]] <- suppressWarnings(as.numeric(x[[z_col]]))
  x[[years_col]] <- suppressWarnings(as.numeric(x[[years_col]]))
  x[[age_col]] <- if (age_col %in% names(x)) suppressWarnings(as.numeric(x[[age_col]])) else NA_real_

  x
}

read_pulmonary_clock_instance <- function(file, instance_id, instance_label, visit_order, visit_index) {
  if (!file.exists(file)) {
    stop("Input file does not exist: ", file)
  }

  x <- fread(file, fill = TRUE, check.names = FALSE)
  x <- as.data.frame(x, check.names = FALSE)

  if (!"participant_id" %in% names(x)) {
    stop("participant_id column not found in: ", file)
  }

  x <- standardize_pulmonary_clock_columns(x)

  risk_col <- paste0(organ_clean, "_proteomics_mortality_risk_score")
  z_col <- paste0(organ_clean, "_proteomics_mortality_clock_acceleration_z")
  years_col <- paste0(organ_clean, "_proteomics_mortality_clock_acceleration_years")
  age_col <- paste0(organ_clean, "_proteomics_mortality_clock_age_years")

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
    age_col
  )

  for (cc in optional_cols) {
    if (!cc %in% names(x)) x[[cc]] <- NA
  }

  x <- x[, c("participant_id", optional_cols), drop = FALSE]

  x <- x %>%
    mutate(
      participant_id = normalize_id(participant_id),
      participant_id_original = participant_id,
      application_instance = instance_id,
      instance_label = instance_label,
      visit_order = visit_order,
      visit_index = visit_index,
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
      clock_age_years = suppressWarnings(as.numeric(.data[[age_col]]))
    ) %>%
    select(
      participant_id,
      participant_id_original,
      application_instance,
      instance_label,
      visit_order,
      visit_index,
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
# 5. Load and wide-format data
# -----------------------------

baseline_file <- file.path(
  baseline_root,
  "Pulmonary_proteomics_mortality_clock",
  "pulmonary_proteomics_mortality_clock_predictions.tsv"
)

instance2_file <- file.path(
  longitudinal_root,
  "Pulmonary",
  "pulmonary_proteomics_mortality_clock_apply_instance_2_0_predictions.tsv"
)

instance3_file <- file.path(
  longitudinal_root,
  "Pulmonary",
  "pulmonary_proteomics_mortality_clock_apply_instance_3_0_predictions.tsv"
)

id_map <- load_id_match_key(id_match_csv)

dat0 <- read_pulmonary_clock_instance(baseline_file, "0_0", "Instance 0 baseline", 0, 0)
dat2 <- read_pulmonary_clock_instance(instance2_file, "2_0", "Instance 2 follow-up", 2, 1)
dat3 <- read_pulmonary_clock_instance(instance3_file, "3_0", "Instance 3 follow-up", 3, 2)

dat0$participant_id <- canonicalize_ids_to_upenn(dat0$participant_id, id_map)
dat2$participant_id <- canonicalize_ids_to_upenn(dat2$participant_id, id_map)
dat3$participant_id <- canonicalize_ids_to_upenn(dat3$participant_id, id_map)

dat_long <- bind_rows(dat0, dat2, dat3) %>%
  filter(!is.na(clock_acceleration_years)) %>%
  arrange(participant_id, visit_index) %>%
  group_by(participant_id, application_instance) %>%
  slice(1) %>%
  ungroup()

fwrite(
  dat_long,
  file.path(main_outdir, "pulmonary_proteomics_longitudinal_acceleration_years_long.tsv"),
  sep = "\t"
)

dat_wide <- dat_long %>%
  select(
    participant_id,
    application_instance,
    clock_acceleration_years,
    clock_acceleration_z,
    chronological_age,
    clock_age_years,
    sample_date,
    death_date,
    admin_censor_date
  ) %>%
  pivot_wider(
    names_from = application_instance,
    values_from = c(
      clock_acceleration_years,
      clock_acceleration_z,
      chronological_age,
      clock_age_years,
      sample_date,
      death_date,
      admin_censor_date
    )
  )

fwrite(
  dat_wide,
  file.path(main_outdir, "pulmonary_proteomics_longitudinal_acceleration_years_wide.tsv"),
  sep = "\t"
)

cat("N baseline:", n_distinct(dat0$participant_id), "\n")
cat("N instance 2:", n_distinct(dat2$participant_id), "\n")
cat("N instance 3:", n_distinct(dat3$participant_id), "\n")
cat("N wide:", nrow(dat_wide), "\n")

# -----------------------------
# 6. Interval-specific Cox models
# -----------------------------

make_interval_df <- function(dat_wide, comparison_id, comparison_label, start_instance, end_instance, start_label, end_label) {
  start_accel_col <- paste0("clock_acceleration_years_", start_instance)
  end_accel_col <- paste0("clock_acceleration_years_", end_instance)
  start_age_col <- paste0("chronological_age_", start_instance)
  end_age_col <- paste0("chronological_age_", end_instance)
  start_z_col <- paste0("clock_acceleration_z_", start_instance)
  end_z_col <- paste0("clock_acceleration_z_", end_instance)
  end_sample_date_col <- paste0("sample_date_", end_instance)
  end_death_date_col <- paste0("death_date_", end_instance)
  end_admin_censor_col <- paste0("admin_censor_date_", end_instance)

  required_cols <- c(start_accel_col, end_accel_col, start_age_col, end_sample_date_col, end_death_date_col, end_admin_censor_col)
  missing_cols <- setdiff(required_cols, names(dat_wide))
  if (length(missing_cols) > 0) {
    stop("Missing columns for ", comparison_id, ": ", paste(missing_cols, collapse = ", "))
  }

  out <- data.frame(
    participant_id = dat_wide$participant_id,
    comparison_id = comparison_id,
    comparison_label = comparison_label,
    start_instance = start_instance,
    end_instance = end_instance,
    start_label = start_label,
    end_label = end_label,
    acceleration_years_start = suppressWarnings(as.numeric(dat_wide[[start_accel_col]])),
    acceleration_years_end = suppressWarnings(as.numeric(dat_wide[[end_accel_col]])),
    acceleration_z_start = if (start_z_col %in% names(dat_wide)) suppressWarnings(as.numeric(dat_wide[[start_z_col]])) else NA_real_,
    acceleration_z_end = if (end_z_col %in% names(dat_wide)) suppressWarnings(as.numeric(dat_wide[[end_z_col]])) else NA_real_,
    chronological_age_start = suppressWarnings(as.numeric(dat_wide[[start_age_col]])),
    chronological_age_end = if (end_age_col %in% names(dat_wide)) suppressWarnings(as.numeric(dat_wide[[end_age_col]])) else NA_real_,
    sample_date_endpoint = safe_as_date(dat_wide[[end_sample_date_col]]),
    death_date_endpoint = safe_as_date(dat_wide[[end_death_date_col]]),
    admin_censor_date_endpoint = safe_as_date(dat_wide[[end_admin_censor_col]]),
    stringsAsFactors = FALSE
  )

  out <- out %>%
    mutate(
      admin_censor_date_endpoint = if_else(
        is.na(admin_censor_date_endpoint),
        admin_censor_date_default,
        admin_censor_date_endpoint
      ),
      delta_accel_years = acceleration_years_end - acceleration_years_start,
      delta_accel_z = acceleration_z_end - acceleration_z_start,
      delta_chrono_age = chronological_age_end - chronological_age_start,
      event_after_endpoint = case_when(
        !is.na(death_date_endpoint) &
          !is.na(sample_date_endpoint) &
          !is.na(admin_censor_date_endpoint) ~
          death_date_endpoint > sample_date_endpoint &
          death_date_endpoint <= admin_censor_date_endpoint,
        TRUE ~ FALSE
      ),
      end_date_endpoint = case_when(
        event_after_endpoint ~ death_date_endpoint,
        !is.na(admin_censor_date_endpoint) ~ admin_censor_date_endpoint,
        TRUE ~ as.Date(NA)
      ),
      time_from_endpoint_years = as.numeric(end_date_endpoint - sample_date_endpoint) / 365.25,
      event_after_endpoint = as.integer(event_after_endpoint)
    )

  out
}

run_interval_model <- function(interval_df) {
  comparison_id <- unique(interval_df$comparison_id)[1]
  comparison_label <- unique(interval_df$comparison_label)[1]
  start_instance <- unique(interval_df$start_instance)[1]
  end_instance <- unique(interval_df$end_instance)[1]

  cox_df <- interval_df %>%
    filter(
      !is.na(time_from_endpoint_years),
      time_from_endpoint_years > 0,
      !is.na(event_after_endpoint),
      !is.na(acceleration_years_start),
      !is.na(delta_accel_years),
      !is.na(chronological_age_start)
    )

  n_model <- nrow(cox_df)
  n_events <- sum(cox_df$event_after_endpoint == 1, na.rm = TRUE)

  cat("\nComparison:", comparison_label, "\n")
  cat("N model:", n_model, "\n")
  cat("N events:", n_events, "\n")

  if (n_model < 50 || n_events < 10) {
    warning("Insufficient sample size/events for ", comparison_label)
    return(list(
      comparison_id = comparison_id,
      comparison_label = comparison_label,
      cox_df = cox_df,
      coef_all_tbl = data.frame(),
      summary_tbl = data.frame(
        comparison_id = comparison_id,
        comparison_label = comparison_label,
        start_instance = start_instance,
        end_instance = end_instance,
        n = n_model,
        n_events = n_events,
        event_rate = ifelse(n_model > 0, n_events / n_model, NA_real_),
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

  cox_baseline <- survival::coxph(
    survival::Surv(time_from_endpoint_years, event_after_endpoint) ~
      acceleration_years_start +
      chronological_age_start,
    data = cox_df
  )

  cox_delta <- survival::coxph(
    survival::Surv(time_from_endpoint_years, event_after_endpoint) ~
      acceleration_years_start +
      delta_accel_years +
      chronological_age_start,
    data = cox_df
  )

  coef_baseline_tbl <- broom::tidy(cox_baseline, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(
      comparison_id = comparison_id,
      comparison_label = comparison_label,
      model = "baseline_only"
    ) %>%
    select(comparison_id, comparison_label, model, everything())

  coef_delta_tbl <- broom::tidy(cox_delta, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(
      comparison_id = comparison_id,
      comparison_label = comparison_label,
      model = "baseline_plus_delta"
    ) %>%
    select(comparison_id, comparison_label, model, everything())

  coef_all_tbl <- bind_rows(coef_baseline_tbl, coef_delta_tbl)

  sink(file.path(main_outdir, paste0("pulmonary_proteomics_", comparison_id, "_cox_baseline_only_summary.txt")))
  print(summary(cox_baseline))
  sink()

  sink(file.path(main_outdir, paste0("pulmonary_proteomics_", comparison_id, "_cox_baseline_plus_delta_summary.txt")))
  print(summary(cox_delta))
  sink()

  lrt_tbl <- as.data.frame(anova(cox_baseline, cox_delta, test = "LRT"))
  fwrite(
    lrt_tbl,
    file.path(main_outdir, paste0("pulmonary_proteomics_", comparison_id, "_cox_delta_added_value_lrt.tsv")),
    sep = "\t",
    row.names = TRUE
  )

  lrt_chisq <- if ("Chisq" %in% names(lrt_tbl)) lrt_tbl$Chisq[2] else NA_real_
  lrt_df <- if ("Df" %in% names(lrt_tbl)) lrt_tbl$Df[2] else NA_real_
  lrt_p <- if ("Pr(>|Chi|)" %in% names(lrt_tbl)) lrt_tbl$`Pr(>|Chi|)`[2] else NA_real_

  delta_term_tbl <- coef_delta_tbl %>% filter(term == "delta_accel_years")
  baseline_term_tbl <- coef_delta_tbl %>% filter(term == "acceleration_years_start")
  age_term_tbl <- coef_delta_tbl %>% filter(term == "chronological_age_start")

  summary_tbl <- data.frame(
    comparison_id = comparison_id,
    comparison_label = comparison_label,
    start_instance = start_instance,
    end_instance = end_instance,
    n = n_model,
    n_events = n_events,
    event_rate = n_events / n_model,

    mean_acceleration_years_start = mean(cox_df$acceleration_years_start, na.rm = TRUE),
    mean_delta_accel_years = mean(cox_df$delta_accel_years, na.rm = TRUE),
    sd_delta_accel_years = sd(cox_df$delta_accel_years, na.rm = TRUE),

    hr_start_accel = baseline_term_tbl$estimate[1],
    hr_start_accel_lower95 = baseline_term_tbl$conf.low[1],
    hr_start_accel_upper95 = baseline_term_tbl$conf.high[1],
    p_start_accel = baseline_term_tbl$p.value[1],

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

  list(
    comparison_id = comparison_id,
    comparison_label = comparison_label,
    cox_df = cox_df,
    cox_baseline = cox_baseline,
    cox_delta = cox_delta,
    coef_all_tbl = coef_all_tbl,
    summary_tbl = summary_tbl,
    lrt_tbl = lrt_tbl
  )
}

interval_dfs <- lapply(seq_len(nrow(comparison_tbl)), function(i) {
  row <- comparison_tbl[i, ]
  make_interval_df(
    dat_wide = dat_wide,
    comparison_id = row$comparison_id,
    comparison_label = row$comparison_label,
    start_instance = row$start_instance,
    end_instance = row$end_instance,
    start_label = row$start_label,
    end_label = row$end_label
  )
})

names(interval_dfs) <- comparison_tbl$comparison_id

interval_all <- bind_rows(interval_dfs)

fwrite(
  interval_all,
  file.path(main_outdir, "pulmonary_proteomics_delta_acceleration_years_interval_dataset.tsv"),
  sep = "\t"
)

all_results <- lapply(interval_dfs, run_interval_model)

summary_all <- bind_rows(lapply(all_results, function(x) x$summary_tbl))
coef_all <- bind_rows(lapply(all_results, function(x) x$coef_all_tbl))

fwrite(
  summary_all,
  file.path(main_outdir, "pulmonary_proteomics_delta_acceleration_years_hr_summary.tsv"),
  sep = "\t"
)

fwrite(
  coef_all,
  file.path(main_outdir, "pulmonary_proteomics_delta_acceleration_years_cox_coefficients.tsv"),
  sep = "\t"
)

print(summary_all)

# -----------------------------
# 7. Forest plot for delta HR
# -----------------------------

plot_tbl <- summary_all %>%
  filter(!is.na(hr_delta), !is.na(hr_delta_lower95), !is.na(hr_delta_upper95)) %>%
  mutate(
    comparison_label = factor(comparison_label, levels = rev(comparison_tbl$comparison_label)),
    hr_ci_label = paste0(
      "HR = ", fmt_num(hr_delta, 2),
      " (", fmt_num(hr_delta_lower95, 2),
      "-", fmt_num(hr_delta_upper95, 2), ")"
    ),
    p_label = paste0("P = ", vapply(p_delta, fmt_p, character(1))),
    lrt_label = paste0("LRT P = ", vapply(lrt_p, fmt_p, character(1))),
    cindex_label = paste0(
      "Delta C = ",
      ifelse(is.na(delta_cindex), "NA", fmt_num(delta_cindex, 3))
    ),
    n_label = paste0("N = ", comma(n), "; deaths = ", comma(n_events)),
    annotation_label = paste0(hr_ci_label, "\n", p_label, "; ", lrt_label, "\n", cindex_label, "; ", n_label)
  )

if (nrow(plot_tbl) > 0) {
  x_min_data <- min(plot_tbl$hr_delta_lower95, na.rm = TRUE)
  x_max_data <- max(plot_tbl$hr_delta_upper95, na.rm = TRUE)

  x_min <- min(0.80, x_min_data * 0.90)
  x_max <- max(1.25, x_max_data * 2.05)
  x_text <- x_max_data * 1.12

  p_forest <- ggplot(plot_tbl, aes(y = comparison_label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.7) +
    geom_segment(
      aes(x = hr_delta_lower95, xend = hr_delta_upper95, yend = comparison_label),
      linewidth = 1.0,
      color = "black"
    ) +
    geom_point(aes(x = hr_delta), size = 3.8, shape = 21, fill = "white", color = "black", stroke = 0.9) +
    geom_text(aes(x = x_text, label = annotation_label), hjust = 0, size = 3.2, lineheight = 0.95) +
    scale_x_log10(
      limits = c(x_min, x_max),
      breaks = c(0.80, 0.90, 1.00, 1.05, 1.10, 1.20, 1.40, 1.60, 2.00),
      labels = number_format(accuracy = 0.01)
    ) +
    labs(
      x = "Hazard ratio per 1-year increase in delta acceleration years",
      y = NULL,
      title = "Pulmonary proteomics mortality-clock longitudinal change and future mortality",
      subtitle = paste0(
        "Landmark Cox models start at the endpoint instance and adjust for ",
        "chronological age plus acceleration years at the start of each interval"
      )
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11),
      axis.text.y = element_text(face = "bold", size = 12),
      axis.text.x = element_text(size = 10),
      plot.margin = margin(10, 120, 10, 10)
    ) +
    coord_cartesian(clip = "off")

  ggsave(
    file.path(main_outdir, "pulmonary_proteomics_delta_acceleration_years_hr_forest_plot.pdf"),
    p_forest,
    width = 11.5,
    height = 4.8
  )

  ggsave(
    file.path(main_outdir, "pulmonary_proteomics_delta_acceleration_years_hr_forest_plot.png"),
    p_forest,
    width = 11.5,
    height = 4.8,
    dpi = 300
  )

  ggsave(
    file.path(main_outdir, "pulmonary_proteomics_delta_acceleration_years_hr_forest_plot.svg"),
    p_forest,
    width = 11.5,
    height = 4.8
  )
}

cat("\n============================================================\n")
cat("Finished Pulmonary proteomics delta acceleration-years survival analysis.\n")
cat("Main output directory:\n", main_outdir, "\n\n")
cat("Main outputs:\n")
cat("  pulmonary_proteomics_delta_acceleration_years_interval_dataset.tsv\n")
cat("  pulmonary_proteomics_delta_acceleration_years_hr_summary.tsv\n")
cat("  pulmonary_proteomics_delta_acceleration_years_cox_coefficients.tsv\n")
cat("  pulmonary_proteomics_delta_acceleration_years_hr_forest_plot.pdf/png/svg\n")
cat("============================================================\n")
