# Longitudinal metabolomics mortality-clock acceleration-years
# analysis across UKB instances
#
# Baseline/model metabolomics instance:
#   Instance 0_0
#
# Longitudinal follow-up metabolomics instance:
#   Instance 1_0
#
# Main y-axis variable:
#   {organ}_metabolomics_mortality_clock_acceleration_years
#
# IMPORTANT CASE DEFINITION FOR PAIRED LONGITUDINAL ANALYSIS:
#   Cases     = participants who died after instance 1 sample date
#               and before administrative censoring.
#   Non-cases = participants alive/censored after instance 1.
#
# This is a landmark definition at instance 1, because delta values
# are only observable among participants who survived and returned
# for instance 1.
#
# Organs:
#   Endocrine, Digestive, Hepatic, Immune
#
# Final plots/tables are generated for:
#   1) Whole sample
#   2) Cases after instance 1
#   3) Non-cases after instance 1
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
  library(survival)
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
  "longitudinal_metabolomics_mortality_clock_acceleration_years_analysis"
)
dir.create(main_outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 2. Plot palette
# -----------------------------

instance_palette <- c(
  "Instance 0\nBaseline" = "#2F4B7C",
  "Instance 1\nFollow-up" = "#E0B43B"
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

  print(diag_tbl)

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
# 5. Standardize metabolomics clock columns
# -----------------------------

standardize_metabolomics_clock_columns <- function(x, organ_clean) {
  risk_col <- paste0(organ_clean, "_metabolomics_mortality_risk_score")
  z_col <- paste0(organ_clean, "_metabolomics_mortality_clock_acceleration_z")
  years_col <- paste0(organ_clean, "_metabolomics_mortality_clock_acceleration_years")
  age_col <- paste0(organ_clean, "_metabolomics_mortality_clock_age_years")

  if (!risk_col %in% names(x)) {
    candidate_risk <- grep(
      paste0(organ_clean, ".*metabolomics.*mortality.*risk_score$|metabolomics_mortality_risk_score$"),
      names(x),
      value = TRUE
    )
    if (length(candidate_risk) > 0) {
      message("Using candidate risk-score column for ", organ_clean, ": ", candidate_risk[1])
      x[[risk_col]] <- suppressWarnings(as.numeric(x[[candidate_risk[1]]]))
    }
  }

  if (!z_col %in% names(x)) {
    merged_z_col <- paste0("risk_15", z_col)

    if (merged_z_col %in% names(x)) {
      message("Detected malformed merged risk_15 + acceleration_z header for ", organ_clean)
      x[[z_col]] <- suppressWarnings(as.numeric(x[[merged_z_col]]))
    } else {
      candidate_z <- grep(
        paste0(organ_clean, ".*metabolomics.*mortality.*clock.*acceleration.*z$|acceleration_z$"),
        names(x),
        value = TRUE
      )
      if (length(candidate_z) > 0) {
        message("Using candidate z-score column for ", organ_clean, ": ", candidate_z[1])
        x[[z_col]] <- suppressWarnings(as.numeric(x[[candidate_z[1]]]))
      }
    }
  }

  if (!z_col %in% names(x)) {
    stop("Could not find z-score acceleration column for ", organ_clean, ": ", z_col)
  }

  x[[risk_col]] <- if (risk_col %in% names(x)) suppressWarnings(as.numeric(x[[risk_col]])) else NA_real_
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
# 6. Read one instance
# -----------------------------

read_metabolomics_clock_instance <- function(file, organ_clean, instance_id, instance_label, visit_order, expected_instance_year) {
  if (!file.exists(file)) {
    stop("Input file does not exist: ", file)
  }

  x <- fread(file, fill = TRUE, check.names = FALSE)
  x <- as.data.frame(x, check.names = FALSE)

  if (!"participant_id" %in% names(x)) {
    stop("participant_id column not found in: ", file)
  }

  x <- standardize_metabolomics_clock_columns(x, organ_clean)

  z_col <- paste0(organ_clean, "_metabolomics_mortality_clock_acceleration_z")
  years_col <- paste0(organ_clean, "_metabolomics_mortality_clock_acceleration_years")
  age_col <- paste0(organ_clean, "_metabolomics_mortality_clock_age_years")
  risk_col <- paste0(organ_clean, "_metabolomics_mortality_risk_score")

  optional_cols <- c(
    "sample_date",
    "death_date",
    "admin_censor_date",
    "end_date",
    "event",
    "time_years",
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

  keep_cols <- c(
    "participant_id",
    "sample_date",
    "death_date",
    "admin_censor_date",
    "end_date",
    "event",
    "time_years",
    "age_at_baseline",
    "age_at_imaging",
    "sex",
    "split",
    risk_col,
    z_col,
    years_col,
    age_col
  )

  x <- x[, keep_cols, drop = FALSE]

  x <- x %>%
    mutate(
      participant_id = normalize_id(participant_id),
      application_instance = instance_id,
      visit_order = visit_order,
      expected_instance_year = expected_instance_year,

      sample_date = safe_as_date(sample_date),
      death_date = safe_as_date(death_date),
      admin_censor_date = safe_as_date(admin_censor_date),
      end_date = safe_as_date(end_date),

      admin_censor_date = if_else(
        is.na(admin_censor_date),
        admin_censor_date_default,
        admin_censor_date
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
      clock_age_years = suppressWarnings(as.numeric(.data[[age_col]])),

      # Keep this event variable for diagnostics only.
      # Formal longitudinal case/non-case definition below is recomputed
      # using instance 1 sample date as the landmark.
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
      sample_year,
      death_date,
      admin_censor_date,
      end_date,
      event,
      time_years,
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
    lmer(clock_acceleration_years ~ visit_order + (1 | participant_id), data = df),
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
# 8. Main analysis per organ
# -----------------------------

run_one_organ <- function(organ_label) {
  organ_clean <- organ_clean_vec[[organ_label]]

  message("\n============================================================")
  message("Running longitudinal metabolomics acceleration-years analysis for: ", organ_label)
  message("Clean organ name: ", organ_clean)
  message("============================================================")

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

  outdir <- file.path(main_outdir, organ_label)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  dat0 <- read_metabolomics_clock_instance(
    file = baseline_file,
    organ_clean = organ_clean,
    instance_id = "0_0",
    instance_label = "Instance 0\nBaseline",
    visit_order = 0,
    expected_instance_year = 2008
  )

  dat1 <- read_metabolomics_clock_instance(
    file = instance1_file,
    organ_clean = organ_clean,
    instance_id = "1_0",
    instance_label = "Instance 1\nFollow-up",
    visit_order = 1,
    expected_instance_year = 2012
  )

  id_harm <- harmonize_pair_ids(
    dat0 = dat0,
    dat1 = dat1,
    id_match_csv = id_match_csv,
    outdir = outdir,
    organ_clean = organ_clean
  )

  dat0 <- id_harm$dat0
  dat1 <- id_harm$dat1

  dat_all <- bind_rows(dat0, dat1) %>%
    mutate(
      application_instance = factor(application_instance, levels = c("0_0", "1_0")),
      instance_label = factor(
        instance_label,
        levels = c("Instance 0\nBaseline", "Instance 1\nFollow-up")
      )
    )

  common_ids <- dat_all %>%
    filter(!is.na(clock_acceleration_years)) %>%
    distinct(participant_id, application_instance) %>%
    count(participant_id, name = "n_instances") %>%
    filter(n_instances == 2) %>%
    pull(participant_id)

  cat("Organ:", organ_label, "\n")
  cat("N in baseline instance 0:", n_distinct(dat0$participant_id), "\n")
  cat("N in follow-up instance 1:", n_distinct(dat1$participant_id), "\n")
  cat("N common across instance 0 and 1:", length(common_ids), "\n")

  if (length(common_ids) == 0) {
    warning(
      "No common participants after ID harmonization for organ: ", organ_label,
      ". Skipping plots and writing diagnostics only."
    )

    fwrite(
      dat0 %>% distinct(participant_id, participant_id_original),
      file.path(outdir, paste0(organ_clean, "_baseline_ids_after_harmonization.tsv")),
      sep = "\t"
    )

    fwrite(
      dat1 %>% distinct(participant_id, participant_id_original),
      file.path(outdir, paste0(organ_clean, "_instance1_ids_after_harmonization.tsv")),
      sep = "\t"
    )

    return(invisible(NULL))
  }

  dat_common <- dat_all %>%
    filter(participant_id %in% common_ids) %>%
    filter(!is.na(clock_acceleration_years)) %>%
    arrange(participant_id, visit_order)

  # ------------------------------------------------------------
  # Correct longitudinal case definition:
  # Landmark status from instance 1.
  #
  # Cases:
  #   Died after instance 1 sample date and before admin censoring.
  #
  # Non-cases:
  #   Alive/censored after instance 1.
  #
  # This replaces the older baseline-based case definition.
  # ------------------------------------------------------------

  case_status_tbl <- dat_common %>%
    filter(application_instance == "1_0") %>%
    group_by(participant_id) %>%
    summarise(
      death_date_instance1 = safe_min_date(death_date),
      sample_date_instance1 = safe_min_date(sample_date),
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
    ) %>%
    select(
      participant_id,
      sample_date_instance1,
      death_date_instance1,
      admin_censor_date_instance1,
      end_date_from_instance1,
      time_from_instance1_years,
      event_after_instance1,
      case_status
    )

  fwrite(
    case_status_tbl,
    file.path(outdir, paste0(organ_clean, "_metabolomics_mortality_landmark_case_status_from_instance1.tsv")),
    sep = "\t"
  )

  dat_common <- dat_common %>%
    left_join(case_status_tbl, by = "participant_id") %>%
    mutate(
      case_status = factor(case_status, levels = c("Cases", "Non-cases"))
    )

  cat("N cases after instance 1 among common population:",
      n_distinct(dat_common$participant_id[dat_common$case_status == "Cases"]), "\n")
  cat("N non-cases after instance 1 among common population:",
      n_distinct(dat_common$participant_id[dat_common$case_status == "Non-cases"]), "\n")

  dat_plot <- bind_rows(
    dat_common %>% mutate(analysis_group = "Whole sample"),
    dat_common %>% mutate(analysis_group = as.character(case_status))
  ) %>%
    filter(!is.na(analysis_group)) %>%
    mutate(
      analysis_group = factor(analysis_group, levels = c("Whole sample", "Cases", "Non-cases"))
    )

  fwrite(
    dat_common,
    file.path(outdir, paste0(organ_clean, "_metabolomics_mortality_acceleration_years_common_population_long.tsv")),
    sep = "\t"
  )

  fwrite(
    dat_plot,
    file.path(outdir, paste0(organ_clean, "_metabolomics_mortality_acceleration_years_plot_dataset_with_groups.tsv")),
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

      mean_acceleration_years = mean(clock_acceleration_years, na.rm = TRUE),
      sd_acceleration_years = sd(clock_acceleration_years, na.rm = TRUE),

      mean_chronological_age = mean(chronological_age, na.rm = TRUE),
      sd_chronological_age = sd(chronological_age, na.rm = TRUE),
      mean_clock_age = mean(clock_age_years, na.rm = TRUE),
      sd_clock_age = sd(clock_age_years, na.rm = TRUE),

      .groups = "drop"
    )

  fwrite(
    summary_tbl,
    file.path(outdir, paste0(organ_clean, "_metabolomics_mortality_acceleration_years_instance_summary_by_group.tsv")),
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
      clock_acceleration_years,
      chronological_age,
      clock_age_years,
      sample_date_instance1,
      death_date_instance1,
      admin_censor_date_instance1,
      end_date_from_instance1,
      time_from_instance1_years,
      event_after_instance1
    ) %>%
    pivot_wider(
      names_from = application_instance,
      values_from = c(
        clock_acceleration_z,
        clock_acceleration_years,
        chronological_age,
        clock_age_years
      )
    )

  required_wide_cols <- c("clock_acceleration_years_0_0", "clock_acceleration_years_1_0")
  missing_wide_cols <- setdiff(required_wide_cols, names(dat_wide))

  if (length(missing_wide_cols) > 0) {
    warning(
      "Missing expected wide columns for ", organ_label, ": ",
      paste(missing_wide_cols, collapse = ", "),
      ". Skipping delta and plotting steps."
    )

    fwrite(
      dat_wide,
      file.path(outdir, paste0(organ_clean, "_metabolomics_mortality_acceleration_years_common_population_wide_INCOMPLETE.tsv")),
      sep = "\t"
    )

    return(invisible(NULL))
  }

  dat_wide <- dat_wide %>%
    mutate(
      delta_z_1_minus_0 = clock_acceleration_z_1_0 - clock_acceleration_z_0_0,
      delta_accel_years_1_minus_0 = clock_acceleration_years_1_0 - clock_acceleration_years_0_0,
      delta_chrono_age_1_minus_0 = chronological_age_1_0 - chronological_age_0_0,
      delta_clock_age_1_minus_0 = clock_age_years_1_0 - clock_age_years_0_0
    )

  fwrite(
    dat_wide,
    file.path(outdir, paste0(organ_clean, "_metabolomics_mortality_acceleration_years_common_population_wide_deltas_by_group.tsv")),
    sep = "\t"
  )

  delta_summary_tbl <- dat_wide %>%
    group_by(analysis_group) %>%
    summarise(
      n = n(),

      mean_acceleration_years_0 = mean(clock_acceleration_years_0_0, na.rm = TRUE),
      mean_acceleration_years_1 = mean(clock_acceleration_years_1_0, na.rm = TRUE),

      mean_delta_accel_years_1_minus_0 = mean(delta_accel_years_1_minus_0, na.rm = TRUE),
      sd_delta_accel_years_1_minus_0 = sd(delta_accel_years_1_minus_0, na.rm = TRUE),
      median_delta_accel_years_1_minus_0 = median(delta_accel_years_1_minus_0, na.rm = TRUE),
      p_wilcox_accel_years_1_vs_0 = safe_wilcox_p(clock_acceleration_years_1_0, clock_acceleration_years_0_0),

      mean_z_0 = mean(clock_acceleration_z_0_0, na.rm = TRUE),
      mean_z_1 = mean(clock_acceleration_z_1_0, na.rm = TRUE),
      mean_delta_z_1_minus_0 = mean(delta_z_1_minus_0, na.rm = TRUE),
      sd_delta_z_1_minus_0 = sd(delta_z_1_minus_0, na.rm = TRUE),
      median_delta_z_1_minus_0 = median(delta_z_1_minus_0, na.rm = TRUE),
      p_wilcox_z_1_vs_0 = safe_wilcox_p(clock_acceleration_z_1_0, clock_acceleration_z_0_0),

      mean_delta_chrono_age_1_minus_0 = mean(delta_chrono_age_1_minus_0, na.rm = TRUE),
      mean_delta_clock_age_1_minus_0 = mean(delta_clock_age_1_minus_0, na.rm = TRUE),

      .groups = "drop"
    )

  fwrite(
    delta_summary_tbl,
    file.path(outdir, paste0(organ_clean, "_metabolomics_mortality_acceleration_years_delta_summary_by_group.tsv")),
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
    file.path(outdir, paste0(organ_clean, "_metabolomics_mortality_acceleration_years_lmm_visit_order_trend_by_group.tsv")),
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
    file.path(outdir, paste0(organ_clean, "_metabolomics_mortality_acceleration_years_lmm_sample_year_trend_by_group.tsv")),
    sep = "\t"
  )

  # Case-status interaction using landmark case status from instance 1.
  interaction_tbl <- tryCatch({
    interaction_dat <- dat_common %>%
      filter(!is.na(case_status)) %>%
      mutate(case_status = relevel(case_status, ref = "Non-cases"))

    if (
      n_distinct(interaction_dat$participant_id) >= 10 &&
      n_distinct(interaction_dat$case_status) == 2
    ) {
      fit_int <- lmer(
        clock_acceleration_years ~ visit_order * case_status + (1 | participant_id),
        data = interaction_dat
      )

      broom.mixed::tidy(fit_int, effects = "fixed", conf.int = TRUE) %>%
        mutate(organ = organ_clean) %>%
        select(organ, everything())
    } else {
      data.frame(
        organ = organ_clean,
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
      organ = organ_clean,
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
    file.path(outdir, paste0(organ_clean, "_metabolomics_mortality_acceleration_years_lmm_case_status_interaction.tsv")),
    sep = "\t"
  )

  # -----------------------------
  # Annotation tables
  # -----------------------------

  global_y_max <- max(dat_plot$clock_acceleration_years, na.rm = TRUE)
  global_y_min <- min(dat_plot$clock_acceleration_years, na.rm = TRUE)
  global_y_range <- global_y_max - global_y_min

  if (!is.finite(global_y_range) || global_y_range == 0) {
    global_y_range <- 1
  }

  ann_tbl <- lmm_tbl %>%
    mutate(
      x = 1,
      y = global_y_max + 0.18 * global_y_range,
      label_short = paste0(
        "LMM \u03B2 = ", fmt_num(estimate, 2),
        " acceleration-years/visit\nP = ", vapply(p.value, fmt_p, character(1))
      )
    )

  delta_ann_tbl <- delta_summary_tbl %>%
    mutate(
      analysis_group = factor(analysis_group, levels = c("Whole sample", "Cases", "Non-cases")),
      x = 1,
      y = global_y_max + 0.10 * global_y_range,
      delta_label = paste0("\u0394 acceleration years: 1-0 = ", fmt_num(mean_delta_accel_years_1_minus_0, 2))
    )

  caution_ann_tbl <- data.frame(
    analysis_group = factor(
      c("Whole sample", "Cases", "Non-cases"),
      levels = c("Whole sample", "Cases", "Non-cases")
    ),
    x = 1,
    y = global_y_max + 0.02 * global_y_range,
    caution_label = "Landmark cases: death after instance 1; paired sample survived/returned for follow-up"
  )

  n_label_tbl <- dat_plot %>%
    group_by(analysis_group, instance_label) %>%
    summarise(
      n_participants = n_distinct(participant_id),
      .groups = "drop"
    ) %>%
    mutate(
      analysis_group = factor(analysis_group, levels = c("Whole sample", "Cases", "Non-cases")),
      y = global_y_min - 0.08 * global_y_range,
      n_label = paste0("N = ", comma(n_participants))
    )

  group_n_tbl <- dat_plot %>%
    distinct(analysis_group, participant_id) %>%
    group_by(analysis_group) %>%
    summarise(n_participants = n_distinct(participant_id), .groups = "drop") %>%
    mutate(
      analysis_group = factor(analysis_group, levels = c("Whole sample", "Cases", "Non-cases"))
    )

  group_n_subtitle <- paste0(
    "Whole sample N = ", comma(group_n_tbl$n_participants[group_n_tbl$analysis_group == "Whole sample"]),
    "; Cases after instance 1 N = ", comma(group_n_tbl$n_participants[group_n_tbl$analysis_group == "Cases"]),
    "; Non-cases after instance 1 N = ", comma(group_n_tbl$n_participants[group_n_tbl$analysis_group == "Non-cases"])
  )

  # -----------------------------
  # Figure 1: Distribution
  # Violin + embedded boxplot only; no participant-level dots.
  # -----------------------------

  p_dist <- ggplot(
    dat_plot,
    aes(x = instance_label, y = clock_acceleration_years, fill = instance_label, color = instance_label)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
    geom_violin(
      trim = FALSE,
      alpha = 0.32,
      linewidth = 0.5,
      scale = "width"
    ) +
    geom_boxplot(
      width = 0.18,
      outlier.shape = NA,
      alpha = 0.82,
      linewidth = 0.55,
      color = "black"
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      size = 2.8,
      shape = 21,
      fill = "white",
      color = "black",
      stroke = 0.7
    ) +
    geom_text(
      data = n_label_tbl,
      aes(x = instance_label, y = y, label = n_label),
      inherit.aes = FALSE,
      size = 3.1,
      fontface = "bold",
      color = "black"
    ) +
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
    facet_wrap(~ analysis_group, nrow = 1, drop = TRUE) +
    scale_fill_manual(values = instance_palette_fill, guide = "none") +
    scale_color_manual(values = instance_palette, guide = "none") +
    coord_cartesian(
      ylim = c(global_y_min - 0.13 * global_y_range, global_y_max + 0.28 * global_y_range),
      clip = "off"
    ) +
    labs(
      x = NULL,
      y = paste0(title_case(organ_clean), " metabolomics mortality-clock acceleration years"),
      title = paste0(title_case(organ_clean), " metabolomics mortality-clock acceleration years"),
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

  ggsave(
    file.path(outdir, paste0(organ_clean, "_metabolomics_mortality_acceleration_years_distribution_by_case_status.pdf")),
    p_dist,
    width = 13.5,
    height = 6.2
  )

  ggsave(
    file.path(outdir, paste0(organ_clean, "_metabolomics_mortality_acceleration_years_distribution_by_case_status.png")),
    p_dist,
    width = 13.5,
    height = 6.2,
    dpi = 300
  )

  ggsave(
    file.path(outdir, paste0(organ_clean, "_metabolomics_mortality_acceleration_years_distribution_by_case_status.svg")),
    p_dist,
    width = 13.5,
    height = 6.2
  )

  # -----------------------------
  # Figure 1b: Individual-level delta acceleration years
  # versus chronological age at baseline
  # -----------------------------

  delta_age_plot_tbl <- dat_wide %>%
    filter(
      !is.na(delta_accel_years_1_minus_0),
      !is.na(chronological_age_0_0)
    ) %>%
    mutate(
      analysis_group = factor(analysis_group, levels = c("Whole sample", "Cases", "Non-cases"))
    )

  fwrite(
    delta_age_plot_tbl,
    file.path(outdir, paste0(organ_clean, "_metabolomics_mortality_delta_acceleration_years_vs_baseline_age_plot_dataset.tsv")),
    sep = "\t"
  )

  delta_age_summary_tbl <- delta_age_plot_tbl %>%
    group_by(analysis_group) %>%
    summarise(
      n = n(),
      n_participants = n_distinct(participant_id),
      mean_baseline_age = mean(chronological_age_0_0, na.rm = TRUE),
      sd_baseline_age = sd(chronological_age_0_0, na.rm = TRUE),
      mean_delta_accel_years = mean(delta_accel_years_1_minus_0, na.rm = TRUE),
      sd_delta_accel_years = sd(delta_accel_years_1_minus_0, na.rm = TRUE),
      cor_pearson = suppressWarnings(cor(
        chronological_age_0_0,
        delta_accel_years_1_minus_0,
        use = "complete.obs",
        method = "pearson"
      )),
      cor_spearman = suppressWarnings(cor(
        chronological_age_0_0,
        delta_accel_years_1_minus_0,
        use = "complete.obs",
        method = "spearman"
      )),
      .groups = "drop"
    )

  fwrite(
    delta_age_summary_tbl,
    file.path(outdir, paste0(organ_clean, "_metabolomics_mortality_delta_acceleration_years_vs_baseline_age_summary.tsv")),
    sep = "\t"
  )

  delta_age_lm_tbl <- delta_age_plot_tbl %>%
    group_by(analysis_group) %>%
    group_modify(~ {
      if (nrow(.x) < 5) {
        return(data.frame(
          term = "chronological_age_0_0",
          estimate = NA_real_,
          std.error = NA_real_,
          statistic = NA_real_,
          p.value = NA_real_,
          conf.low = NA_real_,
          conf.high = NA_real_
        ))
      }

      fit <- lm(delta_accel_years_1_minus_0 ~ chronological_age_0_0, data = .x)

      broom::tidy(fit, conf.int = TRUE) %>%
        filter(term == "chronological_age_0_0")
    }) %>%
    ungroup() %>%
    mutate(
      label = paste0(
        "N = ", comma(delta_age_summary_tbl$n_participants[
          match(as.character(analysis_group), as.character(delta_age_summary_tbl$analysis_group))
        ]),
        "\n\u03B2 = ", fmt_num(estimate, 2),
        " years/year\nP = ", vapply(p.value, fmt_p, character(1))
      )
    )

  fwrite(
    delta_age_lm_tbl,
    file.path(outdir, paste0(organ_clean, "_metabolomics_mortality_delta_acceleration_years_vs_baseline_age_lm.tsv")),
    sep = "\t"
  )

  delta_y_max <- max(delta_age_plot_tbl$delta_accel_years_1_minus_0, na.rm = TRUE)
  delta_y_min <- min(delta_age_plot_tbl$delta_accel_years_1_minus_0, na.rm = TRUE)
  delta_y_range <- delta_y_max - delta_y_min

  if (!is.finite(delta_y_range) || delta_y_range == 0) {
    delta_y_range <- 1
  }

  delta_age_ann_tbl <- delta_age_lm_tbl %>%
    mutate(
      analysis_group = factor(analysis_group, levels = c("Whole sample", "Cases", "Non-cases")),
      x = min(delta_age_plot_tbl$chronological_age_0_0, na.rm = TRUE),
      y = delta_y_max + 0.08 * delta_y_range
    )

  p_delta_age <- ggplot(
    delta_age_plot_tbl,
    aes(x = chronological_age_0_0, y = delta_accel_years_1_minus_0)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.7) +
    geom_point(
      aes(color = analysis_group),
      alpha = 0.22,
      size = 0.75,
      stroke = 0
    ) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      linewidth = 0.9,
      color = "black"
    ) +
    geom_text(
      data = delta_age_ann_tbl,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      size = 3.1,
      fontface = "bold"
    ) +
    facet_wrap(~ analysis_group, nrow = 1, drop = TRUE) +
    scale_color_manual(values = group_palette, guide = "none") +
    coord_cartesian(
      ylim = c(delta_y_min, delta_y_max + 0.22 * delta_y_range),
      clip = "off"
    ) +
    labs(
      x = "Chronological age at baseline, instance 0 (years)",
      y = paste0(
        "\u0394 ",
        title_case(organ_clean),
        " metabolomics mortality-clock acceleration years\n(instance 1 - instance 0)"
      ),
      title = paste0(
        title_case(organ_clean),
        " longitudinal change in metabolomics mortality-clock acceleration years"
      ),
      subtitle = "Landmark groups defined by death/censoring after instance 1"
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 11),
      strip.text = element_text(face = "bold", size = 12),
      plot.margin = margin(10, 20, 10, 10)
    )

  ggsave(
    file.path(outdir, paste0(organ_clean, "_metabolomics_mortality_delta_acceleration_years_vs_baseline_age_by_case_status.pdf")),
    p_delta_age,
    width = 13.5,
    height = 6.2
  )

  ggsave(
    file.path(outdir, paste0(organ_clean, "_metabolomics_mortality_delta_acceleration_years_vs_baseline_age_by_case_status.png")),
    p_delta_age,
    width = 13.5,
    height = 6.2,
    dpi = 300
  )

  ggsave(
    file.path(outdir, paste0(organ_clean, "_metabolomics_mortality_delta_acceleration_years_vs_baseline_age_by_case_status.svg")),
    p_delta_age,
    width = 13.5,
    height = 6.2
  )

  # -----------------------------
  # Figure 1c: Individual-level delta acceleration years
  # versus baseline acceleration years
  #
  # y-axis:
  #   acceleration years at instance 1 - acceleration years at instance 0
  #
  # x-axis:
  #   acceleration years at baseline, instance 0
  #
  # Facets:
  #   Whole sample, Cases, Non-cases
  # -----------------------------

  delta_vs_baseline_accel_plot_tbl <- dat_wide %>%
    filter(
      !is.na(delta_accel_years_1_minus_0),
      !is.na(clock_acceleration_years_0_0)
    ) %>%
    mutate(
      analysis_group = factor(
        analysis_group,
        levels = c("Whole sample", "Cases", "Non-cases")
      )
    )

  fwrite(
    delta_vs_baseline_accel_plot_tbl,
    file.path(
      outdir,
      paste0(
        organ_clean,
        "_metabolomics_mortality_delta_acceleration_years_vs_baseline_acceleration_years_plot_dataset.tsv"
      )
    ),
    sep = "\t"
  )

  delta_vs_baseline_accel_summary_tbl <- delta_vs_baseline_accel_plot_tbl %>%
    group_by(analysis_group) %>%
    summarise(
      n = n(),
      n_participants = n_distinct(participant_id),
      mean_baseline_accel_years = mean(clock_acceleration_years_0_0, na.rm = TRUE),
      sd_baseline_accel_years = sd(clock_acceleration_years_0_0, na.rm = TRUE),
      mean_delta_accel_years = mean(delta_accel_years_1_minus_0, na.rm = TRUE),
      sd_delta_accel_years = sd(delta_accel_years_1_minus_0, na.rm = TRUE),
      cor_pearson = suppressWarnings(cor(
        clock_acceleration_years_0_0,
        delta_accel_years_1_minus_0,
        use = "complete.obs",
        method = "pearson"
      )),
      cor_spearman = suppressWarnings(cor(
        clock_acceleration_years_0_0,
        delta_accel_years_1_minus_0,
        use = "complete.obs",
        method = "spearman"
      )),
      .groups = "drop"
    )

  fwrite(
    delta_vs_baseline_accel_summary_tbl,
    file.path(
      outdir,
      paste0(
        organ_clean,
        "_metabolomics_mortality_delta_acceleration_years_vs_baseline_acceleration_years_summary.tsv"
      )
    ),
    sep = "\t"
  )

  delta_vs_baseline_accel_lm_tbl <- delta_vs_baseline_accel_plot_tbl %>%
    group_by(analysis_group) %>%
    group_modify(~ {
      if (nrow(.x) < 5) {
        return(data.frame(
          term = "clock_acceleration_years_0_0",
          estimate = NA_real_,
          std.error = NA_real_,
          statistic = NA_real_,
          p.value = NA_real_,
          conf.low = NA_real_,
          conf.high = NA_real_
        ))
      }

      fit <- lm(delta_accel_years_1_minus_0 ~ clock_acceleration_years_0_0, data = .x)

      broom::tidy(fit, conf.int = TRUE) %>%
        filter(term == "clock_acceleration_years_0_0")
    }) %>%
    ungroup() %>%
    mutate(
      label = paste0(
        "N = ", comma(delta_vs_baseline_accel_summary_tbl$n_participants[
          match(analysis_group, delta_vs_baseline_accel_summary_tbl$analysis_group)
        ]),
        "\n\u03B2 = ", fmt_num(estimate, 2),
        " years/year\nP = ", vapply(p.value, fmt_p, character(1))
      )
    )

  fwrite(
    delta_vs_baseline_accel_lm_tbl,
    file.path(
      outdir,
      paste0(
        organ_clean,
        "_metabolomics_mortality_delta_acceleration_years_vs_baseline_acceleration_years_lm.tsv"
      )
    ),
    sep = "\t"
  )

  delta_ba_y_max <- max(delta_vs_baseline_accel_plot_tbl$delta_accel_years_1_minus_0, na.rm = TRUE)
  delta_ba_y_min <- min(delta_vs_baseline_accel_plot_tbl$delta_accel_years_1_minus_0, na.rm = TRUE)
  delta_ba_y_range <- delta_ba_y_max - delta_ba_y_min

  if (!is.finite(delta_ba_y_range) || delta_ba_y_range == 0) {
    delta_ba_y_range <- 1
  }

  delta_ba_x_min <- min(delta_vs_baseline_accel_plot_tbl$clock_acceleration_years_0_0, na.rm = TRUE)

  delta_vs_baseline_accel_ann_tbl <- delta_vs_baseline_accel_lm_tbl %>%
    mutate(
      analysis_group = factor(
        analysis_group,
        levels = c("Whole sample", "Cases", "Non-cases")
      ),
      x = delta_ba_x_min,
      y = delta_ba_y_max + 0.08 * delta_ba_y_range
    )

  p_delta_vs_baseline_accel <- ggplot(
    delta_vs_baseline_accel_plot_tbl,
    aes(x = clock_acceleration_years_0_0, y = delta_accel_years_1_minus_0)
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.7
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dotted",
      color = "grey50",
      linewidth = 0.6
    ) +
    geom_point(
      aes(color = analysis_group),
      alpha = 0.22,
      size = 0.75,
      stroke = 0
    ) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      linewidth = 0.9,
      color = "black"
    ) +
    geom_text(
      data = delta_vs_baseline_accel_ann_tbl,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      size = 3.1,
      fontface = "bold"
    ) +
    facet_wrap(~ analysis_group, nrow = 1, drop = TRUE) +
    scale_color_manual(values = group_palette, guide = "none") +
    coord_cartesian(
      ylim = c(delta_ba_y_min, delta_ba_y_max + 0.22 * delta_ba_y_range),
      clip = "off"
    ) +
    labs(
      x = paste0(
        "Baseline ",
        title_case(organ_clean),
        " metabolomics mortality-clock acceleration years\n(instance 0)"
      ),
      y = paste0(
        "\u0394 ",
        title_case(organ_clean),
        " metabolomics mortality-clock acceleration years\n(instance 1 - instance 0)"
      ),
      title = paste0(
        title_case(organ_clean),
        " longitudinal change in metabolomics mortality-clock acceleration years"
      ),
      subtitle = "Individual-level change plotted against baseline acceleration years"
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 11),
      strip.text = element_text(face = "bold", size = 12),
      plot.margin = margin(10, 20, 10, 10)
    )

  ggsave(
    file.path(
      outdir,
      paste0(
        organ_clean,
        "_metabolomics_mortality_delta_acceleration_years_vs_baseline_acceleration_years_by_case_status.pdf"
      )
    ),
    p_delta_vs_baseline_accel,
    width = 13.5,
    height = 6.2
  )

  ggsave(
    file.path(
      outdir,
      paste0(
        organ_clean,
        "_metabolomics_mortality_delta_acceleration_years_vs_baseline_acceleration_years_by_case_status.png"
      )
    ),
    p_delta_vs_baseline_accel,
    width = 13.5,
    height = 6.2,
    dpi = 300
  )

  ggsave(
    file.path(
      outdir,
      paste0(
        organ_clean,
        "_metabolomics_mortality_delta_acceleration_years_vs_baseline_acceleration_years_by_case_status.svg"
      )
    ),
    p_delta_vs_baseline_accel,
    width = 13.5,
    height = 6.2
  )

  ### check whether the delta can provide additioinal power beyond baseline age and baseline clock acceleartion
  cox_delta <- coxph(
    Surv(time_from_instance1_years, event_after_instance1) ~
      clock_acceleration_years_0_0 +
      delta_accel_years_1_minus_0 +
      chronological_age_0_0,    data = dat_wide %>% filter(analysis_group == "Whole sample")
  )

  summary(cox_delta)

  # -----------------------------
  # Figure 2: Mean +/- 95% CI trend for acceleration years
  # -----------------------------

  mean_se_tbl <- dat_plot %>%
    group_by(analysis_group, visit_order, instance_label) %>%
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
    file.path(outdir, paste0(organ_clean, "_metabolomics_mortality_acceleration_years_mean_se_by_instance_and_group.tsv")),
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
    facet_wrap(~ analysis_group, nrow = 1, drop = TRUE) +
    scale_fill_manual(values = instance_palette_fill, guide = "none") +
    scale_x_continuous(
      breaks = c(0, 1),
      labels = c("Instance 0\nBaseline", "Instance 1\nFollow-up")
    ) +
    coord_cartesian(
      ylim = c(global_y_min, global_y_max + 0.25 * global_y_range),
      clip = "off"
    ) +
    labs(
      x = NULL,
      y = paste0("Mean ", title_case(organ_clean), " metabolomics mortality-clock acceleration years"),
      title = paste0("Mean ", title_case(organ_clean), " metabolomics mortality-clock acceleration-years trajectory"),
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

  ggsave(
    file.path(outdir, paste0(organ_clean, "_metabolomics_mortality_acceleration_years_mean_trend_by_case_status.pdf")),
    p_mean,
    width = 13.5,
    height = 6.2
  )

  ggsave(
    file.path(outdir, paste0(organ_clean, "_metabolomics_mortality_acceleration_years_mean_trend_by_case_status.png")),
    p_mean,
    width = 13.5,
    height = 6.2,
    dpi = 300
  )

  # -----------------------------
  # Figure 3: Spaghetti plot for acceleration years
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
      mean = mean(clock_acceleration_years, na.rm = TRUE),
      .groups = "drop"
    )

  p_spaghetti <- ggplot(
    dat_spaghetti,
    aes(x = visit_order, y = clock_acceleration_years, group = participant_id)
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
    facet_wrap(~ analysis_group, nrow = 1, drop = TRUE) +
    scale_fill_manual(values = instance_palette_fill, guide = "none") +
    scale_x_continuous(
      breaks = c(0, 1),
      labels = c("Instance 0\nBaseline", "Instance 1\nFollow-up")
    ) +
    labs(
      x = NULL,
      y = paste0(title_case(organ_clean), " metabolomics mortality-clock acceleration years"),
      title = paste0("Within-person ", title_case(organ_clean), " metabolomics mortality-clock acceleration-years trajectories"),
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
    file.path(outdir, paste0(organ_clean, "_metabolomics_mortality_acceleration_years_spaghetti_by_case_status.pdf")),
    p_spaghetti,
    width = 13.5,
    height = 6.2
  )

  ggsave(
    file.path(outdir, paste0(organ_clean, "_metabolomics_mortality_acceleration_years_spaghetti_by_case_status.png")),
    p_spaghetti,
    width = 13.5,
    height = 6.2,
    dpi = 300
  )

  # -----------------------------
  # Figure 4: acceleration years versus chronological age
  # -----------------------------

  p_vs_age <- ggplot(
    dat_plot,
    aes(x = chronological_age, y = clock_acceleration_years, color = instance_label)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.7) +
    geom_point(alpha = 0.25, size = 0.8) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
    facet_wrap(~ analysis_group, nrow = 1, drop = TRUE) +
    scale_color_manual(values = instance_palette, name = "Instance") +
    labs(
      x = "Chronological age at metabolomics assessment (years)",
      y = paste0(title_case(organ_clean), " metabolomics mortality-clock acceleration years"),
      title = paste0(title_case(organ_clean), " metabolomics mortality-clock acceleration years versus chronological age"),
      subtitle = "Landmark groups defined by death/censoring after instance 1"
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "bottom"
    )

  ggsave(
    file.path(outdir, paste0(organ_clean, "_metabolomics_mortality_acceleration_years_vs_chronological_age_by_case_status.pdf")),
    p_vs_age,
    width = 13.5,
    height = 6.2
  )

  ggsave(
    file.path(outdir, paste0(organ_clean, "_metabolomics_mortality_acceleration_years_vs_chronological_age_by_case_status.png")),
    p_vs_age,
    width = 13.5,
    height = 6.2,
    dpi = 300
  )

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
    file.path(outdir, paste0(organ_clean, "_metabolomics_mortality_acceleration_years_figure_annotation_text_by_group.tsv")),
    sep = "\t"
  )

  cat("\nFinished organ:", organ_label, "\n")
  cat("Output directory:\n", outdir, "\n\n")
  cat("Main outputs:\n")
  cat("  ", paste0(organ_clean, "_metabolomics_mortality_landmark_case_status_from_instance1.tsv"), "\n")
  cat("  ", paste0(organ_clean, "_metabolomics_mortality_acceleration_years_common_population_long.tsv"), "\n")
  cat("  ", paste0(organ_clean, "_metabolomics_mortality_acceleration_years_plot_dataset_with_groups.tsv"), "\n")
  cat("  ", paste0(organ_clean, "_metabolomics_mortality_acceleration_years_instance_summary_by_group.tsv"), "\n")
  cat("  ", paste0(organ_clean, "_metabolomics_mortality_acceleration_years_common_population_wide_deltas_by_group.tsv"), "\n")
  cat("  ", paste0(organ_clean, "_metabolomics_mortality_acceleration_years_delta_summary_by_group.tsv"), "\n")
  cat("  ", paste0(organ_clean, "_metabolomics_mortality_acceleration_years_lmm_visit_order_trend_by_group.tsv"), "\n")
  cat("  ", paste0(organ_clean, "_metabolomics_mortality_acceleration_years_lmm_sample_year_trend_by_group.tsv"), "\n")
  cat("  ", paste0(organ_clean, "_metabolomics_mortality_acceleration_years_lmm_case_status_interaction.tsv"), "\n")
  cat("  ", paste0(organ_clean, "_metabolomics_mortality_acceleration_years_distribution_by_case_status.pdf/png/svg"), "\n")
  cat("  ", paste0(organ_clean, "_metabolomics_mortality_delta_acceleration_years_vs_baseline_age_by_case_status.pdf/png/svg"), "\n")
  cat("  ", paste0(organ_clean, "_metabolomics_mortality_acceleration_years_mean_trend_by_case_status.pdf/png"), "\n")
  cat("  ", paste0(organ_clean, "_metabolomics_mortality_acceleration_years_spaghetti_by_case_status.pdf/png"), "\n")
  cat("  ", paste0(organ_clean, "_metabolomics_mortality_acceleration_years_vs_chronological_age_by_case_status.pdf/png"), "\n\n")

  invisible(list(
    dat_common = dat_common,
    dat_plot = dat_plot,
    dat_wide = dat_wide,
    summary_tbl = summary_tbl,
    delta_summary_tbl = delta_summary_tbl,
    lmm_tbl = lmm_tbl,
    case_status_tbl = case_status_tbl
  ))
}

# -----------------------------
# 9. Run all four metabolomics clocks
# -----------------------------

all_results <- list()

for (organ_label in organ_labels) {
  all_results[[organ_label]] <- run_one_organ(organ_label)
}

cat("\n============================================================\n")
cat("Finished longitudinal metabolomics mortality-clock acceleration-years analysis.\n")
cat("Main output directory:\n", main_outdir, "\n")
cat("============================================================\n")