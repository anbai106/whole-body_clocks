# ============================================================
# 5_pulmonary_proteomics_longitudinal_trend_years.R
#
# Longitudinal Pulmonary proteomics mortality-clock
# acceleration-years analysis across UKB instances.
#
# Baseline/model instance:
#   0_0
#
# Longitudinal follow-up proteomics instances:
#   2_0 and 3_0
#
# Main y-axis variable:
#   pulmonary_proteomics_mortality_clock_acceleration_years
#
# Landmark case definition for descriptive plots:
#   Cases     = participants who died after their last available
#               follow-up proteomics visit among instance 2_0/3_0
#               and before administrative censoring.
#   Non-cases = participants alive/censored after their last available
#               follow-up proteomics visit.
#
# This is used only for descriptive longitudinal plots. Formal
# interval-specific survival analyses are implemented in script 6.
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
  library(broom)
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
  "longitudinal_pulmonary_proteomics_mortality_clock_acceleration_years_analysis"
)
dir.create(main_outdir, recursive = TRUE, showWarnings = FALSE)

instance_info <- tibble::tibble(
  application_instance = c("0_0", "2_0", "3_0"),
  instance_label = c("Instance 0\nBaseline", "Instance 2\nFollow-up", "Instance 3\nFollow-up"),
  visit_order = c(0, 2, 3),
  visit_index = c(0, 1, 2),
  expected_instance_year = c(2008, 2014, 2019)
)

# -----------------------------
# 2. Plot palette
# -----------------------------

instance_palette <- c(
  "Instance 0\nBaseline" = "#2F4B7C",
  "Instance 2\nFollow-up" = "#E0B43B",
  "Instance 3\nFollow-up" = "#7A5195"
)

instance_palette_fill <- alpha(instance_palette, 0.75)

group_palette <- c(
  "Whole sample" = "#2F4B7C",
  "Cases" = "#B85C38",
  "Non-cases" = "#3A7D6B"
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

title_case <- function(x) {
  tools::toTitleCase(gsub("_", " ", x))
}

safe_wilcox_p <- function(x, y) {
  idx <- complete.cases(x, y)
  if (sum(idx) < 2) return(NA_real_)
  tryCatch(
    wilcox.test(x[idx], y[idx], paired = TRUE)$p.value,
    error = function(e) NA_real_
  )
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

safe_max_date <- function(x) {
  x <- as.Date(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) return(as.Date(NA))
  max(x)
}

normalize_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("\\.0$", "", x)
  x[x %in% c("", "NA", "NaN", "NULL", "null")] <- NA_character_
  x
}

# -----------------------------
# 4. ID harmonization
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
# 5. Standardize pulmonary proteomics clock columns
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
      message("Using candidate z-score column: ", candidate_z[1])
      x[[z_col]] <- suppressWarnings(as.numeric(x[[candidate_z[1]]]))
    }
  }

  if (!years_col %in% names(x)) {
    stop("Could not find acceleration-years column: ", years_col)
  }

  x[[risk_col]] <- if (risk_col %in% names(x)) suppressWarnings(as.numeric(x[[risk_col]])) else NA_real_
  x[[z_col]] <- if (z_col %in% names(x)) suppressWarnings(as.numeric(x[[z_col]])) else NA_real_
  x[[years_col]] <- suppressWarnings(as.numeric(x[[years_col]]))
  x[[age_col]] <- if (age_col %in% names(x)) suppressWarnings(as.numeric(x[[age_col]])) else NA_real_

  x
}

# -----------------------------
# 6. Read one instance
# -----------------------------

read_pulmonary_clock_instance <- function(file, instance_id) {
  if (!file.exists(file)) {
    stop("Input file does not exist: ", file)
  }

  info <- instance_info %>% filter(application_instance == instance_id)
  if (nrow(info) != 1) stop("Unknown instance_id: ", instance_id)

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
    if (!cc %in% names(x)) {
      x[[cc]] <- NA
    }
  }

  x <- x[, c("participant_id", optional_cols), drop = FALSE]

  x <- x %>%
    mutate(
      participant_id = normalize_id(participant_id),
      participant_id_original = participant_id,
      application_instance = instance_id,
      instance_label = info$instance_label[1],
      visit_order = info$visit_order[1],
      visit_index = info$visit_index[1],
      expected_instance_year = info$expected_instance_year[1],

      sample_date = safe_as_date(sample_date),
      death_date = safe_as_date(death_date),
      admin_censor_date = safe_as_date(admin_censor_date),
      end_date = safe_as_date(end_date),

      admin_censor_date = case_when(
        !is.na(admin_censor_date) ~ admin_censor_date,
        is.na(admin_censor_date) & !is.na(end_date) ~ end_date,
        TRUE ~ admin_censor_date_default
      ),

      sample_year = as.numeric(format(sample_date, "%Y")),

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
      expected_instance_year,
      sample_date,
      sample_year,
      death_date,
      admin_censor_date,
      end_date,
      age_at_baseline,
      age_at_imaging,
      chronological_age,
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
# 7. LMM helpers
# -----------------------------

run_lmm_by_group <- function(df, group_name) {
  n_pid <- n_distinct(df$participant_id)
  n_obs <- nrow(df)
  n_visits <- n_distinct(df$visit_index)

  if (n_pid < 5 || n_obs < 10 || n_visits < 2) {
    return(data.frame(
      analysis_group = group_name,
      n_participants = n_pid,
      n_observations = n_obs,
      term = "visit_index",
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
    lmer(clock_acceleration_years ~ visit_index + (1 | participant_id), data = df),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(data.frame(
      analysis_group = group_name,
      n_participants = n_pid,
      n_observations = n_obs,
      term = "visit_index",
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
    filter(term == "visit_index")

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
      ": LMM beta = ", fmt_num(tt$estimate[1], 2),
      " acceleration-years/visit (95% CI ",
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
    lmer(clock_acceleration_years ~ sample_year + (1 | participant_id), data = df),
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
# 8. Main analysis
# -----------------------------

message("\n============================================================")
message("Running longitudinal Pulmonary proteomics acceleration-years analysis")
message("============================================================")

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

outdir <- main_outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

id_map <- load_id_match_key(id_match_csv)

dat0 <- read_pulmonary_clock_instance(baseline_file, "0_0")
dat2 <- read_pulmonary_clock_instance(instance2_file, "2_0")
dat3 <- read_pulmonary_clock_instance(instance3_file, "3_0")

# Canonicalize all IDs to id_upenn when an id_match file is available.
dat0$participant_id <- canonicalize_ids_to_upenn(dat0$participant_id, id_map)
dat2$participant_id <- canonicalize_ids_to_upenn(dat2$participant_id, id_map)
dat3$participant_id <- canonicalize_ids_to_upenn(dat3$participant_id, id_map)

id_diagnostics <- tibble::tibble(
  dataset = c("instance_0", "instance_2", "instance_3"),
  n_ids = c(n_distinct(dat0$participant_id), n_distinct(dat2$participant_id), n_distinct(dat3$participant_id)),
  n_overlap_with_baseline = c(
    n_distinct(dat0$participant_id),
    length(intersect(unique(dat0$participant_id), unique(dat2$participant_id))),
    length(intersect(unique(dat0$participant_id), unique(dat3$participant_id)))
  )
)

fwrite(
  id_diagnostics,
  file.path(outdir, "pulmonary_proteomics_id_harmonization_diagnostics.tsv"),
  sep = "\t"
)

print(id_diagnostics)

dat_all <- bind_rows(dat0, dat2, dat3) %>%
  mutate(
    application_instance = factor(application_instance, levels = c("0_0", "2_0", "3_0")),
    instance_label = factor(
      instance_label,
      levels = c("Instance 0\nBaseline", "Instance 2\nFollow-up", "Instance 3\nFollow-up")
    )
  )

# Keep participants with baseline and at least one follow-up instance.
eligible_ids <- dat_all %>%
  filter(!is.na(clock_acceleration_years)) %>%
  distinct(participant_id, application_instance) %>%
  mutate(
    has_baseline = application_instance == "0_0",
    has_followup = application_instance %in% c("2_0", "3_0")
  ) %>%
  group_by(participant_id) %>%
  summarise(
    has_baseline = any(has_baseline),
    has_followup = any(has_followup),
    n_instances = n_distinct(application_instance),
    .groups = "drop"
  ) %>%
  filter(has_baseline, has_followup) %>%
  pull(participant_id)

complete_triplet_ids <- dat_all %>%
  filter(!is.na(clock_acceleration_years)) %>%
  distinct(participant_id, application_instance) %>%
  count(participant_id, name = "n_instances") %>%
  filter(n_instances == 3) %>%
  pull(participant_id)

cat("N in baseline instance 0:", n_distinct(dat0$participant_id), "\n")
cat("N in follow-up instance 2:", n_distinct(dat2$participant_id), "\n")
cat("N in follow-up instance 3:", n_distinct(dat3$participant_id), "\n")
cat("N eligible baseline + at least one follow-up:", length(eligible_ids), "\n")
cat("N complete across 0, 2, and 3:", length(complete_triplet_ids), "\n")

dat_common <- dat_all %>%
  filter(participant_id %in% eligible_ids) %>%
  filter(!is.na(clock_acceleration_years)) %>%
  arrange(participant_id, visit_index)

# Landmark status from each participant's last available follow-up among 2_0/3_0.
case_status_tbl <- dat_common %>%
  filter(application_instance %in% c("2_0", "3_0")) %>%
  arrange(participant_id, visit_index) %>%
  group_by(participant_id) %>%
  slice_tail(n = 1) %>%
  ungroup() %>%
  transmute(
    participant_id,
    landmark_instance = as.character(application_instance),
    landmark_visit_index = visit_index,
    sample_date_landmark = sample_date,
    death_date_landmark = death_date,
    admin_censor_date_landmark = if_else(
      is.na(admin_censor_date),
      admin_censor_date_default,
      admin_censor_date
    ),
    event_after_landmark = case_when(
      !is.na(death_date_landmark) &
        !is.na(sample_date_landmark) &
        !is.na(admin_censor_date_landmark) ~
        death_date_landmark > sample_date_landmark &
        death_date_landmark <= admin_censor_date_landmark,
      TRUE ~ FALSE
    ),
    end_date_from_landmark = case_when(
      event_after_landmark ~ death_date_landmark,
      !is.na(admin_censor_date_landmark) ~ admin_censor_date_landmark,
      TRUE ~ as.Date(NA)
    ),
    time_from_landmark_years = as.numeric(end_date_from_landmark - sample_date_landmark) / 365.25,
    case_status = case_when(
      is.na(sample_date_landmark) ~ NA_character_,
      event_after_landmark ~ "Cases",
      TRUE ~ "Non-cases"
    )
  )

fwrite(
  case_status_tbl,
  file.path(outdir, "pulmonary_proteomics_landmark_case_status_from_last_followup.tsv"),
  sep = "\t"
)

dat_common <- dat_common %>%
  left_join(case_status_tbl, by = "participant_id") %>%
  mutate(case_status = factor(case_status, levels = c("Cases", "Non-cases")))

cat("N cases after last available follow-up:", n_distinct(dat_common$participant_id[dat_common$case_status == "Cases"]), "\n")
cat("N non-cases after last available follow-up:", n_distinct(dat_common$participant_id[dat_common$case_status == "Non-cases"]), "\n")

dat_plot <- bind_rows(
  dat_common %>% mutate(analysis_group = "Whole sample"),
  dat_common %>% mutate(analysis_group = as.character(case_status))
) %>%
  filter(!is.na(analysis_group)) %>%
  mutate(analysis_group = factor(analysis_group, levels = c("Whole sample", "Cases", "Non-cases")))

fwrite(
  dat_common,
  file.path(outdir, "pulmonary_proteomics_mortality_acceleration_years_common_population_long.tsv"),
  sep = "\t"
)

fwrite(
  dat_plot,
  file.path(outdir, "pulmonary_proteomics_mortality_acceleration_years_plot_dataset_with_groups.tsv"),
  sep = "\t"
)

# -----------------------------
# 9. Summary statistics and deltas
# -----------------------------

summary_tbl <- dat_plot %>%
  group_by(analysis_group, application_instance, instance_label) %>%
  summarise(
    n_rows = n(),
    n_participants = n_distinct(participant_id),
    mean_acceleration_years = mean(clock_acceleration_years, na.rm = TRUE),
    sd_acceleration_years = sd(clock_acceleration_years, na.rm = TRUE),
    median_acceleration_years = median(clock_acceleration_years, na.rm = TRUE),
    q25_acceleration_years = quantile(clock_acceleration_years, 0.25, na.rm = TRUE),
    q75_acceleration_years = quantile(clock_acceleration_years, 0.75, na.rm = TRUE),
    min_acceleration_years = min(clock_acceleration_years, na.rm = TRUE),
    max_acceleration_years = max(clock_acceleration_years, na.rm = TRUE),
    mean_z = mean(clock_acceleration_z, na.rm = TRUE),
    sd_z = sd(clock_acceleration_z, na.rm = TRUE),
    mean_chronological_age = mean(chronological_age, na.rm = TRUE),
    sd_chronological_age = sd(chronological_age, na.rm = TRUE),
    mean_clock_age = mean(clock_age_years, na.rm = TRUE),
    sd_clock_age = sd(clock_age_years, na.rm = TRUE),
    .groups = "drop"
  )

fwrite(
  summary_tbl,
  file.path(outdir, "pulmonary_proteomics_mortality_acceleration_years_instance_summary_by_group.tsv"),
  sep = "\t"
)

print(summary_tbl)

dat_wide <- dat_plot %>%
  select(
    participant_id,
    analysis_group,
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
  ) %>%
  mutate(
    delta_accel_years_2_minus_0 = clock_acceleration_years_2_0 - clock_acceleration_years_0_0,
    delta_accel_years_3_minus_0 = clock_acceleration_years_3_0 - clock_acceleration_years_0_0,
    delta_accel_years_3_minus_2 = clock_acceleration_years_3_0 - clock_acceleration_years_2_0,
    delta_z_2_minus_0 = clock_acceleration_z_2_0 - clock_acceleration_z_0_0,
    delta_z_3_minus_0 = clock_acceleration_z_3_0 - clock_acceleration_z_0_0,
    delta_z_3_minus_2 = clock_acceleration_z_3_0 - clock_acceleration_z_2_0,
    delta_chrono_age_2_minus_0 = chronological_age_2_0 - chronological_age_0_0,
    delta_chrono_age_3_minus_0 = chronological_age_3_0 - chronological_age_0_0,
    delta_chrono_age_3_minus_2 = chronological_age_3_0 - chronological_age_2_0,
    delta_clock_age_2_minus_0 = clock_age_years_2_0 - clock_age_years_0_0,
    delta_clock_age_3_minus_0 = clock_age_years_3_0 - clock_age_years_0_0,
    delta_clock_age_3_minus_2 = clock_age_years_3_0 - clock_age_years_2_0
  )

fwrite(
  dat_wide,
  file.path(outdir, "pulmonary_proteomics_mortality_acceleration_years_common_population_wide_deltas_by_group.tsv"),
  sep = "\t"
)

delta_summary_tbl <- dat_wide %>%
  group_by(analysis_group) %>%
  summarise(
    n = n(),
    mean_accel_0 = mean(clock_acceleration_years_0_0, na.rm = TRUE),
    mean_accel_2 = mean(clock_acceleration_years_2_0, na.rm = TRUE),
    mean_accel_3 = mean(clock_acceleration_years_3_0, na.rm = TRUE),
    mean_delta_accel_years_2_minus_0 = mean(delta_accel_years_2_minus_0, na.rm = TRUE),
    mean_delta_accel_years_3_minus_0 = mean(delta_accel_years_3_minus_0, na.rm = TRUE),
    mean_delta_accel_years_3_minus_2 = mean(delta_accel_years_3_minus_2, na.rm = TRUE),
    sd_delta_accel_years_2_minus_0 = sd(delta_accel_years_2_minus_0, na.rm = TRUE),
    sd_delta_accel_years_3_minus_0 = sd(delta_accel_years_3_minus_0, na.rm = TRUE),
    sd_delta_accel_years_3_minus_2 = sd(delta_accel_years_3_minus_2, na.rm = TRUE),
    p_wilcox_2_vs_0 = safe_wilcox_p(clock_acceleration_years_2_0, clock_acceleration_years_0_0),
    p_wilcox_3_vs_0 = safe_wilcox_p(clock_acceleration_years_3_0, clock_acceleration_years_0_0),
    p_wilcox_3_vs_2 = safe_wilcox_p(clock_acceleration_years_3_0, clock_acceleration_years_2_0),
    .groups = "drop"
  )

fwrite(
  delta_summary_tbl,
  file.path(outdir, "pulmonary_proteomics_mortality_acceleration_years_delta_summary_by_group.tsv"),
  sep = "\t"
)

print(delta_summary_tbl)

# -----------------------------
# 10. LMM trend tests
# -----------------------------

lmm_tbl <- bind_rows(
  run_lmm_by_group(dat_plot %>% filter(analysis_group == "Whole sample"), "Whole sample"),
  run_lmm_by_group(dat_plot %>% filter(analysis_group == "Cases"), "Cases"),
  run_lmm_by_group(dat_plot %>% filter(analysis_group == "Non-cases"), "Non-cases")
) %>%
  mutate(analysis_group = factor(analysis_group, levels = c("Whole sample", "Cases", "Non-cases")))

fwrite(
  lmm_tbl,
  file.path(outdir, "pulmonary_proteomics_mortality_acceleration_years_lmm_visit_index_trend_by_group.tsv"),
  sep = "\t"
)

print(lmm_tbl)

lmm_year_tbl <- bind_rows(
  run_lmm_year_by_group(dat_plot %>% filter(analysis_group == "Whole sample"), "Whole sample"),
  run_lmm_year_by_group(dat_plot %>% filter(analysis_group == "Cases"), "Cases"),
  run_lmm_year_by_group(dat_plot %>% filter(analysis_group == "Non-cases"), "Non-cases")
) %>%
  mutate(analysis_group = factor(analysis_group, levels = c("Whole sample", "Cases", "Non-cases")))

fwrite(
  lmm_year_tbl,
  file.path(outdir, "pulmonary_proteomics_mortality_acceleration_years_lmm_sample_year_trend_by_group.tsv"),
  sep = "\t"
)

interaction_tbl <- tryCatch({
  interaction_dat <- dat_common %>%
    filter(!is.na(case_status)) %>%
    mutate(case_status = relevel(case_status, ref = "Non-cases"))

  if (n_distinct(interaction_dat$participant_id) >= 10 && n_distinct(interaction_dat$case_status) == 2) {
    fit_int <- lmer(clock_acceleration_years ~ visit_index * case_status + (1 | participant_id), data = interaction_dat)
    broom.mixed::tidy(fit_int, effects = "fixed", conf.int = TRUE) %>%
      mutate(organ = organ_clean) %>%
      select(organ, everything())
  } else {
    data.frame(
      organ = organ_clean,
      term = "visit_index:case_statusCases",
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
    organ = organ_clean,
    term = "visit_index:case_statusCases",
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
  file.path(outdir, "pulmonary_proteomics_mortality_acceleration_years_lmm_case_status_interaction.tsv"),
  sep = "\t"
)

# -----------------------------
# 11. Annotation tables
# -----------------------------

global_y_max <- max(dat_plot$clock_acceleration_years, na.rm = TRUE)
global_y_min <- min(dat_plot$clock_acceleration_years, na.rm = TRUE)
global_y_range <- global_y_max - global_y_min
if (!is.finite(global_y_range) || global_y_range == 0) global_y_range <- 1

ann_tbl <- lmm_tbl %>%
  mutate(
    x = 1,
    y = global_y_max + 0.22 * global_y_range,
    label_short = paste0(
      "LMM beta = ", fmt_num(estimate, 2),
      " years/visit\nP = ", vapply(p.value, fmt_p, character(1))
    )
  )

delta_ann_tbl <- delta_summary_tbl %>%
  mutate(
    analysis_group = factor(analysis_group, levels = c("Whole sample", "Cases", "Non-cases")),
    x = 1,
    y = global_y_max + 0.11 * global_y_range,
    delta_label = paste0(
      "Delta 2-0 = ", fmt_num(mean_delta_accel_years_2_minus_0, 2),
      "; Delta 3-0 = ", fmt_num(mean_delta_accel_years_3_minus_0, 2)
    )
  )

caution_ann_tbl <- data.frame(
  analysis_group = factor(c("Whole sample", "Cases", "Non-cases"), levels = c("Whole sample", "Cases", "Non-cases")),
  x = 1,
  y = global_y_max + 0.03 * global_y_range,
  caution_label = "Cases defined by death after last available follow-up proteomics visit"
)

n_label_tbl <- dat_plot %>%
  group_by(analysis_group, instance_label) %>%
  summarise(n_participants = n_distinct(participant_id), .groups = "drop") %>%
  mutate(
    analysis_group = factor(analysis_group, levels = c("Whole sample", "Cases", "Non-cases")),
    y = global_y_min - 0.08 * global_y_range,
    n_label = paste0("N = ", comma(n_participants))
  )

group_n_tbl <- dat_plot %>%
  distinct(analysis_group, participant_id) %>%
  group_by(analysis_group) %>%
  summarise(n_participants = n_distinct(participant_id), .groups = "drop") %>%
  mutate(analysis_group = factor(analysis_group, levels = c("Whole sample", "Cases", "Non-cases")))

group_n_subtitle <- paste0(
  "Whole sample N = ", comma(group_n_tbl$n_participants[group_n_tbl$analysis_group == "Whole sample"]),
  "; Cases N = ", comma(group_n_tbl$n_participants[group_n_tbl$analysis_group == "Cases"]),
  "; Non-cases N = ", comma(group_n_tbl$n_participants[group_n_tbl$analysis_group == "Non-cases"])
)

# -----------------------------
# 12. Figure 1: Distribution
# -----------------------------

p_dist <- ggplot(
  dat_plot,
  aes(x = instance_label, y = clock_acceleration_years, fill = instance_label, color = instance_label)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  geom_violin(trim = FALSE, alpha = 0.32, linewidth = 0.5, scale = "width") +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.82, linewidth = 0.55, color = "black") +
  stat_summary(fun = mean, geom = "point", size = 2.8, shape = 21, fill = "white", color = "black", stroke = 0.7) +
  geom_text(data = n_label_tbl, aes(x = instance_label, y = y, label = n_label), inherit.aes = FALSE, size = 3.1, fontface = "bold", color = "black") +
  geom_text(data = ann_tbl, aes(x = x, y = y, label = label_short), inherit.aes = FALSE, size = 3.2, fontface = "bold") +
  geom_text(data = delta_ann_tbl, aes(x = x, y = y, label = delta_label), inherit.aes = FALSE, size = 2.9) +
  geom_text(data = caution_ann_tbl, aes(x = x, y = y, label = caution_label), inherit.aes = FALSE, size = 2.7, fontface = "italic", color = "grey30") +
  facet_wrap(~ analysis_group, nrow = 1, drop = TRUE) +
  scale_fill_manual(values = instance_palette_fill, guide = "none") +
  scale_color_manual(values = instance_palette, guide = "none") +
  coord_cartesian(ylim = c(global_y_min - 0.13 * global_y_range, global_y_max + 0.34 * global_y_range), clip = "off") +
  labs(
    x = NULL,
    y = "Pulmonary proteomics mortality-clock acceleration years",
    title = "Pulmonary proteomics mortality-clock acceleration years across longitudinal instances",
    subtitle = group_n_subtitle
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 10, face = "bold"),
    plot.margin = margin(10, 20, 25, 10)
  )

ggsave(file.path(outdir, "pulmonary_proteomics_mortality_acceleration_years_distribution_by_case_status.pdf"), p_dist, width = 14.5, height = 6.2)
ggsave(file.path(outdir, "pulmonary_proteomics_mortality_acceleration_years_distribution_by_case_status.png"), p_dist, width = 14.5, height = 6.2, dpi = 300)
ggsave(file.path(outdir, "pulmonary_proteomics_mortality_acceleration_years_distribution_by_case_status.svg"), p_dist, width = 14.5, height = 6.2)

# -----------------------------
# 13. Figure 2: Mean +/- 95% CI trend
# -----------------------------

mean_se_tbl <- dat_plot %>%
  group_by(analysis_group, visit_index, visit_order, instance_label) %>%
  summarise(
    n = n(),
    n_participants = n_distinct(participant_id),
    mean = mean(clock_acceleration_years, na.rm = TRUE),
    se = sd(clock_acceleration_years, na.rm = TRUE) / sqrt(n()),
    lower = mean - 1.96 * se,
    upper = mean + 1.96 * se,
    .groups = "drop"
  )

fwrite(
  mean_se_tbl,
  file.path(outdir, "pulmonary_proteomics_mortality_acceleration_years_mean_se_by_instance_and_group.tsv"),
  sep = "\t"
)

p_mean <- ggplot(mean_se_tbl, aes(x = visit_index, y = mean, group = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  geom_line(linewidth = 1.3, color = "#1C1C1C") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = alpha("#6C8EBF", 0.20), color = NA) +
  geom_point(aes(fill = instance_label), shape = 21, size = 4.2, color = "black") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.05, linewidth = 0.6) +
  geom_text(data = ann_tbl, aes(x = x, y = y, label = label_short), inherit.aes = FALSE, size = 3.2, fontface = "bold") +
  geom_text(data = delta_ann_tbl, aes(x = x, y = y, label = delta_label), inherit.aes = FALSE, size = 2.9) +
  geom_text(data = caution_ann_tbl, aes(x = x, y = y, label = caution_label), inherit.aes = FALSE, size = 2.7, fontface = "italic", color = "grey30") +
  facet_wrap(~ analysis_group, nrow = 1, drop = TRUE) +
  scale_fill_manual(values = instance_palette_fill, guide = "none") +
  scale_x_continuous(breaks = c(0, 1, 2), labels = c("Instance 0\nBaseline", "Instance 2\nFollow-up", "Instance 3\nFollow-up")) +
  coord_cartesian(ylim = c(global_y_min, global_y_max + 0.31 * global_y_range), clip = "off") +
  labs(
    x = NULL,
    y = "Mean pulmonary proteomics mortality-clock acceleration years",
    title = "Mean pulmonary proteomics mortality-clock acceleration-years trajectory",
    subtitle = group_n_subtitle
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 10, face = "bold"),
    plot.margin = margin(10, 20, 10, 10)
  )

ggsave(file.path(outdir, "pulmonary_proteomics_mortality_acceleration_years_mean_trend_by_case_status.pdf"), p_mean, width = 14.5, height = 6.2)
ggsave(file.path(outdir, "pulmonary_proteomics_mortality_acceleration_years_mean_trend_by_case_status.png"), p_mean, width = 14.5, height = 6.2, dpi = 300)
ggsave(file.path(outdir, "pulmonary_proteomics_mortality_acceleration_years_mean_trend_by_case_status.svg"), p_mean, width = 14.5, height = 6.2)

# -----------------------------
# 14. Figure 3: Spaghetti plot
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
  group_by(analysis_group, visit_index, instance_label) %>%
  summarise(mean = mean(clock_acceleration_years, na.rm = TRUE), .groups = "drop")

p_spaghetti <- ggplot(dat_spaghetti, aes(x = visit_index, y = clock_acceleration_years, group = participant_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  geom_line(alpha = 0.10, color = "grey50") +
  geom_point(alpha = 0.15, size = 0.7, color = "grey50") +
  geom_line(data = mean_traj_tbl, aes(x = visit_index, y = mean), inherit.aes = FALSE, linewidth = 1.4, color = "#1C1C1C") +
  geom_point(data = mean_traj_tbl, aes(x = visit_index, y = mean, fill = instance_label), inherit.aes = FALSE, shape = 21, color = "black", size = 4) +
  facet_wrap(~ analysis_group, nrow = 1, drop = TRUE) +
  scale_fill_manual(values = instance_palette_fill, guide = "none") +
  scale_x_continuous(breaks = c(0, 1, 2), labels = c("Instance 0\nBaseline", "Instance 2\nFollow-up", "Instance 3\nFollow-up")) +
  labs(
    x = NULL,
    y = "Pulmonary proteomics mortality-clock acceleration years",
    title = "Within-person pulmonary proteomics mortality-clock acceleration-years trajectories",
    subtitle = "Thin lines show sampled participants; black line shows mean trajectory"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 10, face = "bold")
  )

ggsave(file.path(outdir, "pulmonary_proteomics_mortality_acceleration_years_spaghetti_by_case_status.pdf"), p_spaghetti, width = 14.5, height = 6.2)
ggsave(file.path(outdir, "pulmonary_proteomics_mortality_acceleration_years_spaghetti_by_case_status.png"), p_spaghetti, width = 14.5, height = 6.2, dpi = 300)

# -----------------------------
# 15. Figure 4: acceleration years versus chronological age
# -----------------------------

p_vs_age <- ggplot(dat_plot, aes(x = chronological_age, y = clock_acceleration_years, color = instance_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.7) +
  geom_point(alpha = 0.25, size = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  facet_wrap(~ analysis_group, nrow = 1, drop = TRUE) +
  scale_color_manual(values = instance_palette, name = "Instance") +
  labs(
    x = "Chronological age at proteomics assessment (years)",
    y = "Pulmonary proteomics mortality-clock acceleration years",
    title = "Pulmonary proteomics mortality-clock acceleration years versus chronological age",
    subtitle = "Groups defined by death/censoring after the last available follow-up proteomics visit"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom"
  )

ggsave(file.path(outdir, "pulmonary_proteomics_mortality_acceleration_years_vs_chronological_age_by_case_status.pdf"), p_vs_age, width = 14.5, height = 6.2)
ggsave(file.path(outdir, "pulmonary_proteomics_mortality_acceleration_years_vs_chronological_age_by_case_status.png"), p_vs_age, width = 14.5, height = 6.2, dpi = 300)

annotation_tbl <- lmm_tbl %>%
  select(analysis_group, n_participants, n_observations, estimate, conf.low, conf.high, p.value, label)

fwrite(
  annotation_tbl,
  file.path(outdir, "pulmonary_proteomics_mortality_acceleration_years_figure_annotation_text_by_group.tsv"),
  sep = "\t"
)

cat("\n============================================================\n")
cat("Finished longitudinal Pulmonary proteomics acceleration-years analysis.\n")
cat("Main output directory:\n", main_outdir, "\n\n")
cat("Main outputs:\n")
cat("  pulmonary_proteomics_mortality_acceleration_years_common_population_long.tsv\n")
cat("  pulmonary_proteomics_mortality_acceleration_years_common_population_wide_deltas_by_group.tsv\n")
cat("  pulmonary_proteomics_mortality_acceleration_years_delta_summary_by_group.tsv\n")
cat("  pulmonary_proteomics_mortality_acceleration_years_distribution_by_case_status.pdf/png/svg\n")
cat("  pulmonary_proteomics_mortality_acceleration_years_mean_trend_by_case_status.pdf/png/svg\n")
cat("  pulmonary_proteomics_mortality_acceleration_years_spaghetti_by_case_status.pdf/png\n")
cat("  pulmonary_proteomics_mortality_acceleration_years_vs_chronological_age_by_case_status.pdf/png\n")
cat("============================================================\n")
