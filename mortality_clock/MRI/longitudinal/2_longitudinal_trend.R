# ============================================================
# Longitudinal MRI mortality-clock z-score (de)acceleration
# analysis for heart and pancreas MRI clocks
#
# Baseline/model MRI instance:
#   Instance 2 imaging clock predictions
#
# Longitudinal follow-up MRI instance:
#   Instance 3 repeat imaging predictions
#
# Main y-axis variable:
#   {organ}_mri_mortality_clock_acceleration_z
#
# Organs:
#   heart
#   pancreas
#
# Groups:
#   1) Whole sample
#   2) Cases: death event during administrative follow-up
#   3) Non-cases
#
# Notes:
#   - Uses clean longitudinal prediction files only.
#   - Does not require original MRI features.
#   - Handles the older baseline pancreas header issue where
#     risk_15 and acceleration_z may have been merged.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
  library(broom.mixed)
  library(scales)
})

# -----------------------------
# 1. Paths and settings
# -----------------------------

baseline_root <- "/Users/hao/cubic-home//Reproducibile_paper/WholeBodyClock"

longitudinal_root <- "/Users/hao/cubic-home//Reproducibile_paper/WholeBodyClock/mortality_clock/longitudinal/imaging"

organs <- c("heart") #### there is no overlap between instance 2 and instance 3 for pancreas

admin_censor_date_default <- as.Date("2022-11-30")

main_outdir <- file.path(
  longitudinal_root,
  "longitudinal_mri_mortality_clock_zscore_acceleration_analysis"
)
dir.create(main_outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 2. Van Gogh-style palette
# -----------------------------

instance_palette <- c(
  "Instance 2\nBaseline MRI"   = "#2F4B7C",
  "Instance 3\nRepeat imaging" = "#E0B43B"
)

instance_palette_fill <- alpha(instance_palette, 0.75)

group_palette <- c(
  "Whole sample" = "#2F4B7C",
  "Cases"        = "#B85C38",
  "Non-cases"    = "#3A7D6B"
)

# -----------------------------
# 3. Helper functions
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

safe_wilcox_p <- function(x, y) {
  idx <- complete.cases(x, y)
  if (sum(idx) < 2) return(NA_real_)
  
  tryCatch(
    wilcox.test(x[idx], y[idx], paired = TRUE)$p.value,
    error = function(e) NA_real_
  )
}

parse_event_column <- function(x) {
  x_chr <- as.character(x)
  
  case_when(
    x_chr %in% c("TRUE", "True", "true", "T", "1", "1.0") ~ TRUE,
    x_chr %in% c("FALSE", "False", "false", "F", "0", "0.0") ~ FALSE,
    TRUE ~ NA
  )
}

safe_min_date <- function(x) {
  x <- as.Date(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) return(as.Date(NA))
  min(x)
}

safe_as_date <- function(x) {
  suppressWarnings(as.Date(as.character(x)))
}

# -----------------------------
# 4. Repair/standardize MRI clock columns
# -----------------------------

standardize_mri_clock_columns <- function(x, organ) {
  z_col <- paste0(organ, "_mri_mortality_clock_acceleration_z")
  years_col <- paste0(organ, "_mri_mortality_clock_acceleration_years")
  age_col <- paste0(organ, "_mri_mortality_clock_age_years")
  
  # Some older baseline files may have a malformed merged header like:
  # risk_15pancreas_mri_mortality_clock_acceleration_z
  merged_z_col <- paste0("risk_15", z_col)
  
  if (!z_col %in% names(x)) {
    if (merged_z_col %in% names(x)) {
      message("Detected merged risk_15 + acceleration_z header for ", organ, ". Attempting repair.")
      
      # In the old malformed output, the real z-score values often appear
      # under the column named acceleration_years.
      if (years_col %in% names(x)) {
        x[[z_col]] <- suppressWarnings(as.numeric(x[[years_col]]))
      } else {
        x[[z_col]] <- suppressWarnings(as.numeric(x[[merged_z_col]]))
      }
    } else {
      candidate_z <- grep(
        paste0(organ, ".*mri.*mortality.*clock.*acceleration.*z$|acceleration_z$"),
        names(x),
        value = TRUE
      )
      
      if (length(candidate_z) > 0) {
        message("Using candidate z-score column for ", organ, ": ", candidate_z[1])
        x[[z_col]] <- suppressWarnings(as.numeric(x[[candidate_z[1]]]))
      }
    }
  }
  
  if (!z_col %in% names(x)) {
    stop("Could not find or repair z-score acceleration column for ", organ, ": ", z_col)
  }
  
  x[[z_col]] <- suppressWarnings(as.numeric(x[[z_col]]))
  
  if (years_col %in% names(x)) {
    x[[years_col]] <- suppressWarnings(as.numeric(x[[years_col]]))
  } else {
    x[[years_col]] <- NA_real_
  }
  
  if (age_col %in% names(x)) {
    x[[age_col]] <- suppressWarnings(as.numeric(x[[age_col]]))
  } else {
    x[[age_col]] <- NA_real_
  }
  
  x
}

# -----------------------------
# 5. Read and standardize each MRI clock instance
# -----------------------------

read_mri_clock_instance <- function(file, organ, instance_id, instance_label, visit_order, expected_instance_year) {
  if (!file.exists(file)) {
    stop("Input file does not exist: ", file)
  }
  
  x <- fread(file, fill = TRUE, check.names = FALSE)
  x <- as.data.frame(x, check.names = FALSE)
  
  if (!"participant_id" %in% names(x)) {
    stop("participant_id column not found in: ", file)
  }
  
  x <- standardize_mri_clock_columns(x, organ)
  
  z_col <- paste0(organ, "_mri_mortality_clock_acceleration_z")
  years_col <- paste0(organ, "_mri_mortality_clock_acceleration_years")
  age_col <- paste0(organ, "_mri_mortality_clock_age_years")
  
  # Baseline MRI prediction file often uses imaging_date instead of sample_date.
  if (!"sample_date" %in% names(x) && "imaging_date" %in% names(x)) {
    x$sample_date <- x$imaging_date
  }
  
  optional_cols <- c(
    "sample_date",
    "imaging_date",
    "death_date",
    "admin_censor_date",
    "end_date",
    "event",
    "age_at_baseline",
    "age_at_imaging",
    "sex",
    z_col,
    years_col,
    age_col
  )
  
  for (cc in optional_cols) {
    if (!cc %in% names(x)) {
      x[[cc]] <- NA
    }
  }
  
  keep_cols <- c(
    "participant_id",
    "sample_date",
    "imaging_date",
    "death_date",
    "admin_censor_date",
    "end_date",
    "event",
    "age_at_baseline",
    "age_at_imaging",
    "sex",
    z_col,
    years_col,
    age_col
  )
  
  x <- x[, keep_cols, drop = FALSE]
  
  x <- x %>%
    mutate(
      participant_id = as.character(participant_id),
      application_instance = instance_id,
      visit_order = visit_order,
      expected_instance_year = expected_instance_year,
      
      sample_date = safe_as_date(sample_date),
      imaging_date = safe_as_date(imaging_date),
      death_date = safe_as_date(death_date),
      admin_censor_date = safe_as_date(admin_censor_date),
      end_date = safe_as_date(end_date),
      
      # Use default admin censor date if missing.
      admin_censor_date = if_else(
        is.na(admin_censor_date),
        admin_censor_date_default,
        admin_censor_date
      ),
      
      sample_year = as.numeric(format(sample_date, "%Y")),
      
      chronological_age = case_when(
        !is.na(as.numeric(age_at_imaging)) ~ as.numeric(age_at_imaging),
        !is.na(as.numeric(age_at_baseline)) ~ as.numeric(age_at_baseline),
        TRUE ~ NA_real_
      ),
      
      clock_acceleration_z = as.numeric(.data[[z_col]]),
      clock_acceleration_years = as.numeric(.data[[years_col]]),
      clock_age_years = as.numeric(.data[[age_col]]),
      
      event_from_column = parse_event_column(event),
      event_from_dates = case_when(
        !is.na(death_date) & !is.na(sample_date) & !is.na(admin_censor_date) ~
          death_date > sample_date & death_date <= admin_censor_date,
        !is.na(death_date) & !is.na(sample_date) & is.na(admin_censor_date) ~
          death_date > sample_date,
        TRUE ~ NA
      ),
      event = case_when(
        !is.na(event_from_column) ~ event_from_column,
        !is.na(event_from_dates) ~ event_from_dates,
        TRUE ~ FALSE
      ),
      
      instance_label = instance_label
    ) %>%
    select(
      participant_id,
      application_instance,
      instance_label,
      visit_order,
      expected_instance_year,
      sample_date,
      imaging_date,
      sample_year,
      death_date,
      admin_censor_date,
      end_date,
      event,
      age_at_baseline,
      age_at_imaging,
      chronological_age,
      sex,
      clock_acceleration_z,
      clock_acceleration_years,
      clock_age_years
    )
  
  x
}

# -----------------------------
# 6. LMM helper functions
# -----------------------------

run_lmm_by_group <- function(df, group_name) {
  n_pid <- n_distinct(df$participant_id)
  n_obs <- nrow(df)
  n_visits <- n_distinct(df$visit_order)
  
  if (n_pid < 5 || n_obs < 10 || n_visits < 2) {
    return(data.frame(
      analysis_group = group_name,
      n_participants = n_pid,
      n_observations = n_obs,
      term = "visit_order",
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      df = NA_real_,
      p.value = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      label = paste0(group_name, ": LMM unavailable")
    ))
  }
  
  fit <- tryCatch(
    lmer(clock_acceleration_z ~ visit_order + (1 | participant_id), data = df),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(data.frame(
      analysis_group = group_name,
      n_participants = n_pid,
      n_observations = n_obs,
      term = "visit_order",
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      df = NA_real_,
      p.value = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      label = paste0(group_name, ": LMM failed")
    ))
  }
  
  tt <- broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE) %>%
    filter(term == "visit_order")
  
  data.frame(
    analysis_group = group_name,
    n_participants = n_pid,
    n_observations = n_obs,
    term = tt$term[1],
    estimate = tt$estimate[1],
    std.error = tt$std.error[1],
    statistic = tt$statistic[1],
    df = if ("df" %in% names(tt)) tt$df[1] else NA_real_,
    p.value = tt$p.value[1],
    conf.low = tt$conf.low[1],
    conf.high = tt$conf.high[1],
    label = paste0(
      group_name,
      ": LMM \u03B2 = ", fmt_num(tt$estimate[1], 2),
      " z/visit (95% CI ",
      fmt_num(tt$conf.low[1], 2), ", ",
      fmt_num(tt$conf.high[1], 2),
      "), P = ", fmt_p(tt$p.value[1])
    )
  )
}

run_lmm_year_by_group <- function(df, group_name) {
  df <- df %>% filter(!is.na(sample_year))
  
  n_pid <- n_distinct(df$participant_id)
  n_obs <- nrow(df)
  
  if (n_pid < 5 || n_obs < 10 || n_distinct(df$sample_year) < 2) {
    return(data.frame(
      analysis_group = group_name,
      n_participants = n_pid,
      n_observations = n_obs,
      term = "sample_year",
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      df = NA_real_,
      p.value = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_
    ))
  }
  
  fit <- tryCatch(
    lmer(clock_acceleration_z ~ sample_year + (1 | participant_id), data = df),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(data.frame(
      analysis_group = group_name,
      n_participants = n_pid,
      n_observations = n_obs,
      term = "sample_year",
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      df = NA_real_,
      p.value = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_
    ))
  }
  
  tt <- broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE) %>%
    filter(term == "sample_year")
  
  data.frame(
    analysis_group = group_name,
    n_participants = n_pid,
    n_observations = n_obs,
    term = tt$term[1],
    estimate = tt$estimate[1],
    std.error = tt$std.error[1],
    statistic = tt$statistic[1],
    df = if ("df" %in% names(tt)) tt$df[1] else NA_real_,
    p.value = tt$p.value[1],
    conf.low = tt$conf.low[1],
    conf.high = tt$conf.high[1]
  )
}

# -----------------------------
# 7. Main analysis function per organ
# -----------------------------

run_one_organ <- function(organ) {
  message("\n============================================================")
  message("Running longitudinal MRI z-score acceleration analysis for: ", organ)
  message("============================================================")
  
  baseline_file <- file.path(
    baseline_root,
    paste0(organ, "_mri_mortality_clock"),
    paste0(organ, "_mri_mortality_clock_predictions.tsv")
  )
  
  instance3_file <- file.path(
    longitudinal_root,
    organ,
    paste0(organ, "_mri_mortality_clock_apply_instance_3_0_predictions.tsv")
  )
  
  outdir <- file.path(main_outdir, organ)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  dat2 <- read_mri_clock_instance(
    file = baseline_file,
    organ = organ,
    instance_id = "2_0",
    instance_label = "Instance 2\nBaseline MRI",
    visit_order = 0,
    expected_instance_year = 2014
  )
  
  dat3 <- read_mri_clock_instance(
    file = instance3_file,
    organ = organ,
    instance_id = "3_0",
    instance_label = "Instance 3\nRepeat imaging",
    visit_order = 1,
    expected_instance_year = 2019
  )
  
  dat_all <- bind_rows(dat2, dat3) %>%
    mutate(
      application_instance = factor(application_instance, levels = c("2_0", "3_0")),
      instance_label = factor(
        instance_label,
        levels = c("Instance 2\nBaseline MRI", "Instance 3\nRepeat imaging")
      )
    )
  
  common_ids <- dat_all %>%
    filter(!is.na(clock_acceleration_z)) %>%
    distinct(participant_id, application_instance) %>%
    count(participant_id, name = "n_instances") %>%
    filter(n_instances == 2) %>%
    pull(participant_id)
  
  dat_common <- dat_all %>%
    filter(participant_id %in% common_ids) %>%
    filter(!is.na(clock_acceleration_z)) %>%
    arrange(participant_id, visit_order)
  
  cat("Organ:", organ, "\n")
  cat("N in baseline instance 2:", n_distinct(dat2$participant_id), "\n")
  cat("N in repeat instance 3:", n_distinct(dat3$participant_id), "\n")
  cat("N common across instance 2 and 3:", length(common_ids), "\n")
  
  # Case status from baseline instance 2.
  case_status_tbl <- dat_common %>%
    filter(application_instance == "2_0") %>%
    group_by(participant_id) %>%
    summarise(
      event_baseline = any(event %in% TRUE, na.rm = TRUE),
      death_date_baseline = safe_min_date(death_date),
      sample_date_baseline = safe_min_date(sample_date),
      admin_censor_date_baseline = safe_min_date(admin_censor_date),
      .groups = "drop"
    ) %>%
    mutate(
      event_from_dates_baseline = case_when(
        !is.na(death_date_baseline) & !is.na(sample_date_baseline) & !is.na(admin_censor_date_baseline) ~
          death_date_baseline > sample_date_baseline & death_date_baseline <= admin_censor_date_baseline,
        !is.na(death_date_baseline) & !is.na(sample_date_baseline) & is.na(admin_censor_date_baseline) ~
          death_date_baseline > sample_date_baseline,
        TRUE ~ FALSE
      ),
      event_any = event_baseline | event_from_dates_baseline,
      case_status = if_else(event_any, "Cases", "Non-cases")
    ) %>%
    select(participant_id, event_any, case_status)
  
  dat_common <- dat_common %>%
    left_join(case_status_tbl, by = "participant_id") %>%
    mutate(
      case_status = factor(case_status, levels = c("Cases", "Non-cases"))
    )
  
  cat("N cases among common population:", n_distinct(dat_common$participant_id[dat_common$case_status == "Cases"]), "\n")
  cat("N non-cases among common population:", n_distinct(dat_common$participant_id[dat_common$case_status == "Non-cases"]), "\n")
  
  dat_plot <- bind_rows(
    dat_common %>% mutate(analysis_group = "Whole sample"),
    dat_common %>% mutate(analysis_group = as.character(case_status))
  ) %>%
    mutate(
      analysis_group = factor(analysis_group, levels = c("Whole sample", "Cases", "Non-cases"))
    )
  
  fwrite(
    dat_common,
    file.path(outdir, paste0(organ, "_mri_mortality_z_common_population_long.tsv")),
    sep = "\t"
  )
  
  fwrite(
    dat_plot,
    file.path(outdir, paste0(organ, "_mri_mortality_z_plot_dataset_with_groups.tsv")),
    sep = "\t"
  )
  
  # -----------------------------
  # Summary statistics
  # -----------------------------
  
  summary_tbl <- dat_plot %>%
    group_by(analysis_group, application_instance, instance_label) %>%
    summarise(
      n_rows = n(),
      n_participants = n_distinct(participant_id),
      
      mean_z = mean(clock_acceleration_z, na.rm = TRUE),
      sd_z = sd(clock_acceleration_z, na.rm = TRUE),
      median_z = median(clock_acceleration_z, na.rm = TRUE),
      q25_z = quantile(clock_acceleration_z, 0.25, na.rm = TRUE),
      q75_z = quantile(clock_acceleration_z, 0.75, na.rm = TRUE),
      min_z = min(clock_acceleration_z, na.rm = TRUE),
      max_z = max(clock_acceleration_z, na.rm = TRUE),
      
      mean_chronological_age = mean(chronological_age, na.rm = TRUE),
      sd_chronological_age = sd(chronological_age, na.rm = TRUE),
      mean_clock_age = mean(clock_age_years, na.rm = TRUE),
      sd_clock_age = sd(clock_age_years, na.rm = TRUE),
      
      .groups = "drop"
    )
  
  fwrite(
    summary_tbl,
    file.path(outdir, paste0(organ, "_mri_mortality_z_instance_summary_by_group.tsv")),
    sep = "\t"
  )
  
  print(summary_tbl)
  
  # -----------------------------
  # Wide table and paired deltas
  # -----------------------------
  
  dat_wide <- dat_plot %>%
    select(
      participant_id,
      analysis_group,
      application_instance,
      clock_acceleration_z,
      chronological_age,
      clock_age_years
    ) %>%
    pivot_wider(
      names_from = application_instance,
      values_from = c(clock_acceleration_z, chronological_age, clock_age_years)
    ) %>%
    mutate(
      delta_z_3_minus_2 = clock_acceleration_z_3_0 - clock_acceleration_z_2_0,
      delta_chrono_age_3_minus_2 = chronological_age_3_0 - chronological_age_2_0,
      delta_clock_age_3_minus_2 = clock_age_years_3_0 - clock_age_years_2_0
    )
  
  fwrite(
    dat_wide,
    file.path(outdir, paste0(organ, "_mri_mortality_z_common_population_wide_deltas_by_group.tsv")),
    sep = "\t"
  )
  
  delta_summary_tbl <- dat_wide %>%
    group_by(analysis_group) %>%
    summarise(
      n = n(),
      
      mean_z_2 = mean(clock_acceleration_z_2_0, na.rm = TRUE),
      mean_z_3 = mean(clock_acceleration_z_3_0, na.rm = TRUE),
      
      mean_delta_z_3_minus_2 = mean(delta_z_3_minus_2, na.rm = TRUE),
      sd_delta_z_3_minus_2 = sd(delta_z_3_minus_2, na.rm = TRUE),
      median_delta_z_3_minus_2 = median(delta_z_3_minus_2, na.rm = TRUE),
      p_wilcox_z_3_vs_2 = safe_wilcox_p(clock_acceleration_z_3_0, clock_acceleration_z_2_0),
      
      mean_delta_chrono_age_3_minus_2 = mean(delta_chrono_age_3_minus_2, na.rm = TRUE),
      mean_delta_clock_age_3_minus_2 = mean(delta_clock_age_3_minus_2, na.rm = TRUE),
      
      .groups = "drop"
    )
  
  fwrite(
    delta_summary_tbl,
    file.path(outdir, paste0(organ, "_mri_mortality_z_delta_summary_by_group.tsv")),
    sep = "\t"
  )
  
  print(delta_summary_tbl)
  
  # -----------------------------
  # LMM trend tests
  # -----------------------------
  
  lmm_tbl <- bind_rows(
    run_lmm_by_group(dat_plot %>% filter(analysis_group == "Whole sample"), "Whole sample"),
    run_lmm_by_group(dat_plot %>% filter(analysis_group == "Cases"), "Cases"),
    run_lmm_by_group(dat_plot %>% filter(analysis_group == "Non-cases"), "Non-cases")
  ) %>%
    mutate(
      analysis_group = factor(analysis_group, levels = c("Whole sample", "Cases", "Non-cases"))
    )
  
  fwrite(
    lmm_tbl,
    file.path(outdir, paste0(organ, "_mri_mortality_z_lmm_visit_order_trend_by_group.tsv")),
    sep = "\t"
  )
  
  print(lmm_tbl)
  
  lmm_year_tbl <- bind_rows(
    run_lmm_year_by_group(dat_plot %>% filter(analysis_group == "Whole sample"), "Whole sample"),
    run_lmm_year_by_group(dat_plot %>% filter(analysis_group == "Cases"), "Cases"),
    run_lmm_year_by_group(dat_plot %>% filter(analysis_group == "Non-cases"), "Non-cases")
  ) %>%
    mutate(
      analysis_group = factor(analysis_group, levels = c("Whole sample", "Cases", "Non-cases"))
    )
  
  fwrite(
    lmm_year_tbl,
    file.path(outdir, paste0(organ, "_mri_mortality_z_lmm_sample_year_trend_by_group.tsv")),
    sep = "\t"
  )
  
  # Interaction test: do cases and non-cases differ in longitudinal change?
  interaction_tbl <- tryCatch({
    interaction_dat <- dat_common %>%
      filter(!is.na(case_status)) %>%
      mutate(case_status = relevel(case_status, ref = "Non-cases"))
    
    if (
      n_distinct(interaction_dat$participant_id) >= 10 &&
      n_distinct(interaction_dat$case_status) == 2
    ) {
      fit_int <- lmer(
        clock_acceleration_z ~ visit_order * case_status + (1 | participant_id),
        data = interaction_dat
      )
      
      broom.mixed::tidy(fit_int, effects = "fixed", conf.int = TRUE) %>%
        mutate(organ = organ) %>%
        select(organ, everything())
    } else {
      data.frame(
        organ = organ,
        term = "visit_order:case_statusCases",
        estimate = NA_real_,
        std.error = NA_real_,
        statistic = NA_real_,
        df = NA_real_,
        p.value = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_
      )
    }
  }, error = function(e) {
    data.frame(
      organ = organ,
      term = "visit_order:case_statusCases",
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      df = NA_real_,
      p.value = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_
    )
  })
  
  fwrite(
    interaction_tbl,
    file.path(outdir, paste0(organ, "_mri_mortality_z_lmm_case_status_interaction.tsv")),
    sep = "\t"
  )
  
  # -----------------------------
  # Annotation text
  # -----------------------------
  
  global_y_max <- max(dat_plot$clock_acceleration_z, na.rm = TRUE)
  global_y_min <- min(dat_plot$clock_acceleration_z, na.rm = TRUE)
  global_y_range <- global_y_max - global_y_min
  
  if (!is.finite(global_y_range) || global_y_range == 0) {
    global_y_range <- 1
  }
  
  ann_tbl <- lmm_tbl %>%
    mutate(
      x = 1,
      y = global_y_max + 0.16 * global_y_range,
      label_short = paste0(
        "LMM \u03B2 = ", fmt_num(estimate, 2),
        " z/visit\nP = ", vapply(p.value, fmt_p, character(1))
      )
    )
  
  delta_ann_tbl <- delta_summary_tbl %>%
    mutate(
      analysis_group = factor(analysis_group, levels = c("Whole sample", "Cases", "Non-cases")),
      x = 1,
      y = global_y_max + 0.08 * global_y_range,
      delta_label = paste0(
        "\u0394 z: 3-2 = ",
        fmt_num(mean_delta_z_3_minus_2, 2)
      )
    )
  
  caution_ann_tbl <- data.frame(
    analysis_group = factor(c("Whole sample", "Cases", "Non-cases"), levels = c("Whole sample", "Cases", "Non-cases")),
    x = 1,
    y = global_y_max + 0.01 * global_y_range,
    caution_label = "Caution: instance-3 application may include cross-instance calibration effects"
  )
  
  # -----------------------------
  # Figure 1: Distribution plot
  # -----------------------------
  
  p_dist <- ggplot(
    dat_plot,
    aes(x = instance_label, y = clock_acceleration_z, fill = instance_label, color = instance_label)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
    geom_violin(trim = FALSE, alpha = 0.28, linewidth = 0.5) +
    geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.70, linewidth = 0.5) +
    geom_jitter(width = 0.06, alpha = 0.12, size = 0.65, show.legend = FALSE) +
    stat_summary(fun = mean, geom = "point", size = 2.8, shape = 21, fill = "white", color = "black") +
    geom_text(
      data = ann_tbl,
      aes(x = x, y = y, label = label_short),
      inherit.aes = FALSE,
      size = 3.2,
      fontface = "bold"
    ) +
    geom_text(
      data = delta_ann_tbl,
      aes(x = x, y = y, label = delta_label),
      inherit.aes = FALSE,
      size = 3.0
    ) +
    geom_text(
      data = caution_ann_tbl,
      aes(x = x, y = y, label = caution_label),
      inherit.aes = FALSE,
      size = 2.7,
      fontface = "italic",
      color = "grey30"
    ) +
    facet_wrap(~ analysis_group, nrow = 1) +
    scale_fill_manual(values = instance_palette_fill, guide = "none") +
    scale_color_manual(values = instance_palette, guide = "none") +
    coord_cartesian(
      ylim = c(global_y_min, global_y_max + 0.25 * global_y_range),
      clip = "off"
    ) +
    labs(
      x = NULL,
      y = paste0(stringr::str_to_title(organ), " MRI mortality-clock (de)acceleration z-score"),
      title = paste0(stringr::str_to_title(organ), " MRI mortality-clock z-score (de)acceleration across UKB instances"),
      subtitle = paste0("Common participants across instance 2 and 3: N = ", length(common_ids))
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 11),
      strip.text = element_text(face = "bold", size = 12),
      axis.text.x = element_text(size = 10, face = "bold"),
      plot.margin = margin(10, 20, 10, 10)
    )
  
  ggsave(
    file.path(outdir, paste0(organ, "_mri_mortality_z_distribution_by_case_status.pdf")),
    p_dist,
    width = 13.5,
    height = 6.2
  )
  
  ggsave(
    file.path(outdir, paste0(organ, "_mri_mortality_z_distribution_by_case_status.png")),
    p_dist,
    width = 13.5,
    height = 6.2,
    dpi = 300
  )
  
  # -----------------------------
  # Figure 2: Mean +/- 95% CI trend
  # -----------------------------
  
  mean_se_tbl <- dat_plot %>%
    group_by(analysis_group, visit_order, instance_label) %>%
    summarise(
      n = n(),
      n_participants = n_distinct(participant_id),
      mean = mean(clock_acceleration_z, na.rm = TRUE),
      se = sd(clock_acceleration_z, na.rm = TRUE) / sqrt(n()),
      lower = mean - 1.96 * se,
      upper = mean + 1.96 * se,
      .groups = "drop"
    )
  
  fwrite(
    mean_se_tbl,
    file.path(outdir, paste0(organ, "_mri_mortality_z_mean_se_by_instance_and_group.tsv")),
    sep = "\t"
  )
  
  p_mean <- ggplot(
    mean_se_tbl,
    aes(x = visit_order, y = mean, group = 1)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
    geom_line(linewidth = 1.3, color = "#1C1C1C") +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = alpha("#6C8EBF", 0.20), color = NA) +
    geom_point(aes(fill = instance_label), shape = 21, size = 4.2, color = "black") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.05, linewidth = 0.6) +
    geom_text(
      data = ann_tbl,
      aes(x = x, y = y, label = label_short),
      inherit.aes = FALSE,
      size = 3.2,
      fontface = "bold"
    ) +
    geom_text(
      data = delta_ann_tbl,
      aes(x = x, y = y, label = delta_label),
      inherit.aes = FALSE,
      size = 3.0
    ) +
    geom_text(
      data = caution_ann_tbl,
      aes(x = x, y = y, label = caution_label),
      inherit.aes = FALSE,
      size = 2.7,
      fontface = "italic",
      color = "grey30"
    ) +
    facet_wrap(~ analysis_group, nrow = 1) +
    scale_fill_manual(values = instance_palette_fill, guide = "none") +
    scale_x_continuous(
      breaks = c(0, 1),
      labels = c("Instance 2\nBaseline MRI", "Instance 3\nRepeat imaging")
    ) +
    coord_cartesian(
      ylim = c(global_y_min, global_y_max + 0.25 * global_y_range),
      clip = "off"
    ) +
    labs(
      x = NULL,
      y = paste0("Mean ", stringr::str_to_title(organ), " MRI mortality-clock z-score"),
      title = paste0("Mean ", stringr::str_to_title(organ), " MRI mortality-clock z-score trajectory"),
      subtitle = paste0("Common participants across instance 2 and 3: N = ", length(common_ids))
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 11),
      strip.text = element_text(face = "bold", size = 12),
      axis.text.x = element_text(size = 10, face = "bold"),
      plot.margin = margin(10, 20, 10, 10)
    )
  
  ggsave(
    file.path(outdir, paste0(organ, "_mri_mortality_z_mean_trend_by_case_status.pdf")),
    p_mean,
    width = 13.5,
    height = 6.2
  )
  
  ggsave(
    file.path(outdir, paste0(organ, "_mri_mortality_z_mean_trend_by_case_status.png")),
    p_mean,
    width = 13.5,
    height = 6.2,
    dpi = 300
  )
  
  # -----------------------------
  # Figure 3: Spaghetti plot
  # -----------------------------
  
  set.seed(2026)
  
  plot_id_tbl <- dat_plot %>%
    distinct(analysis_group, participant_id) %>%
    group_by(analysis_group) %>%
    group_modify(~ {
      n_take <- min(250, nrow(.x))
      .x[sample(seq_len(nrow(.x)), n_take), , drop = FALSE]
    }) %>%
    ungroup()
  
  dat_spaghetti <- dat_plot %>%
    inner_join(plot_id_tbl, by = c("analysis_group", "participant_id"))
  
  mean_traj_tbl <- dat_plot %>%
    group_by(analysis_group, visit_order, instance_label) %>%
    summarise(
      mean = mean(clock_acceleration_z, na.rm = TRUE),
      .groups = "drop"
    )
  
  p_spaghetti <- ggplot(
    dat_spaghetti,
    aes(x = visit_order, y = clock_acceleration_z, group = participant_id)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
    geom_line(alpha = 0.10, color = "grey50") +
    geom_point(alpha = 0.15, size = 0.7, color = "grey50") +
    geom_line(
      data = mean_traj_tbl,
      aes(x = visit_order, y = mean),
      inherit.aes = FALSE,
      linewidth = 1.4,
      color = "#1C1C1C"
    ) +
    geom_point(
      data = mean_traj_tbl,
      aes(x = visit_order, y = mean, fill = instance_label),
      inherit.aes = FALSE,
      shape = 21,
      color = "black",
      size = 4
    ) +
    facet_wrap(~ analysis_group, nrow = 1) +
    scale_fill_manual(values = instance_palette_fill, guide = "none") +
    scale_x_continuous(
      breaks = c(0, 1),
      labels = c("Instance 2\nBaseline MRI", "Instance 3\nRepeat imaging")
    ) +
    labs(
      x = NULL,
      y = paste0(stringr::str_to_title(organ), " MRI mortality-clock z-score"),
      title = paste0("Within-person ", stringr::str_to_title(organ), " MRI mortality-clock z-score trajectories"),
      subtitle = "Thin lines show sampled participants; black line shows mean trajectory"
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 11),
      strip.text = element_text(face = "bold", size = 12),
      axis.text.x = element_text(size = 10, face = "bold")
    )
  
  ggsave(
    file.path(outdir, paste0(organ, "_mri_mortality_z_spaghetti_by_case_status.pdf")),
    p_spaghetti,
    width = 13.5,
    height = 6.2
  )
  
  ggsave(
    file.path(outdir, paste0(organ, "_mri_mortality_z_spaghetti_by_case_status.png")),
    p_spaghetti,
    width = 13.5,
    height = 6.2,
    dpi = 300
  )
  
  # -----------------------------
  # Figure 4: QC plot vs chronological age
  # -----------------------------
  
  p_vs_age <- ggplot(
    dat_plot,
    aes(x = chronological_age, y = clock_acceleration_z, color = instance_label)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.7) +
    geom_point(alpha = 0.25, size = 0.8) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
    facet_wrap(~ analysis_group, nrow = 1) +
    scale_color_manual(values = instance_palette, name = "Instance") +
    labs(
      x = "Chronological age at MRI assessment (years)",
      y = paste0(stringr::str_to_title(organ), " MRI mortality-clock z-score"),
      title = paste0(stringr::str_to_title(organ), " MRI mortality-clock z-score versus chronological age"),
      subtitle = "Dashed line indicates zero z-score acceleration"
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "bottom"
    )
  
  ggsave(
    file.path(outdir, paste0(organ, "_mri_mortality_z_vs_chronological_age_by_case_status.pdf")),
    p_vs_age,
    width = 13.5,
    height = 6.2
  )
  
  ggsave(
    file.path(outdir, paste0(organ, "_mri_mortality_z_vs_chronological_age_by_case_status.png")),
    p_vs_age,
    width = 13.5,
    height = 6.2,
    dpi = 300
  )
  
  # -----------------------------
  # Save annotation table
  # -----------------------------
  
  annotation_tbl <- lmm_tbl %>%
    select(
      analysis_group,
      n_participants,
      n_observations,
      estimate,
      conf.low,
      conf.high,
      p.value,
      label
    )
  
  fwrite(
    annotation_tbl,
    file.path(outdir, paste0(organ, "_mri_mortality_z_figure_annotation_text_by_group.tsv")),
    sep = "\t"
  )
  
  cat("\nFinished organ:", organ, "\n")
  cat("Output directory:\n", outdir, "\n\n")
  cat("Main outputs:\n")
  cat("  ", paste0(organ, "_mri_mortality_z_common_population_long.tsv"), "\n")
  cat("  ", paste0(organ, "_mri_mortality_z_plot_dataset_with_groups.tsv"), "\n")
  cat("  ", paste0(organ, "_mri_mortality_z_instance_summary_by_group.tsv"), "\n")
  cat("  ", paste0(organ, "_mri_mortality_z_common_population_wide_deltas_by_group.tsv"), "\n")
  cat("  ", paste0(organ, "_mri_mortality_z_delta_summary_by_group.tsv"), "\n")
  cat("  ", paste0(organ, "_mri_mortality_z_lmm_visit_order_trend_by_group.tsv"), "\n")
  cat("  ", paste0(organ, "_mri_mortality_z_lmm_case_status_interaction.tsv"), "\n")
  cat("  ", paste0(organ, "_mri_mortality_z_distribution_by_case_status.pdf/png"), "\n")
  cat("  ", paste0(organ, "_mri_mortality_z_mean_trend_by_case_status.pdf/png"), "\n")
  cat("  ", paste0(organ, "_mri_mortality_z_spaghetti_by_case_status.pdf/png"), "\n")
  cat("  ", paste0(organ, "_mri_mortality_z_vs_chronological_age_by_case_status.pdf/png"), "\n\n")
  
  invisible(list(
    dat_common = dat_common,
    dat_plot = dat_plot,
    summary_tbl = summary_tbl,
    delta_summary_tbl = delta_summary_tbl,
    lmm_tbl = lmm_tbl
  ))
}

# -----------------------------
# 8. Run both organs
# -----------------------------

all_results <- list()

for (organ in organs) {
  all_results[[organ]] <- run_one_organ(organ)
}

cat("\n============================================================\n")
cat("Finished longitudinal MRI mortality-clock z-score acceleration analysis for all organs.\n")
cat("Main output directory:\n", main_outdir, "\n")
cat("============================================================\n")