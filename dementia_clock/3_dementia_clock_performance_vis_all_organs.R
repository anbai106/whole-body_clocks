# ============================================================
# Dementia L'EPOCH model performance visualization
# MRI + proteomics + metabolomics
# RStudio-ready direct-run script
#
# Adapted from all mortality-clock performance visualization.
#
# Expected per-clock files:
#   <prefix>_predictions.tsv
#   <prefix>_performance.json
#   <prefix>_model_comparison.tsv
#   <prefix>_incremental_value_delta_cindex.tsv
#
# Example:
#   brain_mri_dementia_clock/
#     brain_mri_dementia_clock_predictions.tsv
#     brain_mri_dementia_clock_performance.json
#     brain_mri_dementia_clock_model_comparison.tsv
#     brain_mri_dementia_clock_incremental_value_delta_cindex.tsv
# ============================================================
.libPaths('/gpfs/fs001/Users/hao/cubic-home/R/x86_64-pc-linux-gnu-library/4.3')
suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(jsonlite)
  library(patchwork)
  library(scales)
  library(glue)
})

# ============================================================
# 1. Settings
# ============================================================

possible_base_dirs <- c(
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock",
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock",
  getwd()
)

base_dir <- possible_base_dirs[file.exists(possible_base_dirs)][1]
if (is.na(base_dir) || is.null(base_dir)) stop("Please set base_dir manually.")

# Use a modest bootstrap number for quick iteration; increase to 1000 for final figures.
horizon_years_default <- 5
n_boot <- 500
skip_missing <- TRUE

# Main cross-clock output directory.
combined_outdir <- file.path(base_dir, "all_dementia_lepoch_model_performance")
dir.create(combined_outdir, recursive = TRUE, showWarnings = FALSE)

message("Base directory: ", base_dir)
message("Combined output directory: ", combined_outdir)

# ============================================================
# 2. Dementia L'EPOCH manifest
# ============================================================

make_manifest <- function() {
  mri_organs <- c(
    "brain",
    "heart",
    "adipose",
    "kidney",
    "liver",
    "pancreas",
    "spleen"
  )

  proteomics_organs <- c(
    "Reproductive_female",
    "Pulmonary",
    "Heart",
    "Brain",
    "Eye",
    "Hepatic",
    "Renal",
    "Reproductive_male",
    "Endocrine",
    "Immune",
    "Skin"
  )

  metabolomics_organs <- c(
    "Endocrine",
    "Digestive",
    "Hepatic",
    "Immune",
    "Metabolic"
  )

  bind_rows(
    tibble(
      modality = "MRI",
      modality_key = "mri",
      organ_folder_name = mri_organs,
      organ_key = stringr::str_to_lower(mri_organs),
      folder = paste0(stringr::str_to_lower(mri_organs), "_mri_dementia_clock"),
      prefix = paste0(stringr::str_to_lower(mri_organs), "_mri_dementia_clock")
    ),
    tibble(
      modality = "Proteomics",
      modality_key = "proteomics",
      organ_folder_name = proteomics_organs,
      organ_key = stringr::str_to_lower(proteomics_organs),
      folder = paste0(proteomics_organs, "_proteomics_dementia_clock"),
      prefix = paste0(stringr::str_to_lower(proteomics_organs), "_proteomics_dementia_clock")
    ),
    tibble(
      modality = "Metabolomics",
      modality_key = "metabolomics",
      organ_folder_name = metabolomics_organs,
      organ_key = stringr::str_to_lower(metabolomics_organs),
      folder = paste0(metabolomics_organs, "_metabolomics_dementia_clock"),
      prefix = paste0(stringr::str_to_lower(metabolomics_organs), "_metabolomics_dementia_clock")
    )
  ) %>%
    mutate(
      organ_label = organ_folder_name %>%
        stringr::str_replace_all("_", " ") %>%
        stringr::str_to_sentence(),
      clock_label = paste(organ_label, modality),
      clock_id = paste(organ_key, modality_key, sep = "__"),
      clock_dir = file.path(base_dir, folder),
      prediction_file = file.path(clock_dir, paste0(prefix, "_predictions.tsv")),
      performance_file = file.path(clock_dir, paste0(prefix, "_performance.json")),
      model_comparison_file = file.path(clock_dir, paste0(prefix, "_model_comparison.tsv")),
      delta_file = file.path(clock_dir, paste0(prefix, "_incremental_value_delta_cindex.tsv"))
    )
}

clock_manifest <- make_manifest()

message("Expected Dementia L'EPOCH clocks: ", nrow(clock_manifest))
print(clock_manifest %>% select(modality, organ_label, folder, prefix))

# ============================================================
# 3. General helpers
# ============================================================

theme_clock <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(color = "#17202A"),
      plot.title = element_text(
        face = "bold",
        size = base_size + 2,
        margin = margin(b = 6)
      ),
      plot.subtitle = element_text(
        size = base_size - 1,
        color = "#566573",
        margin = margin(b = 8)
      ),
      plot.caption = element_text(
        size = base_size - 3,
        color = "#7B7D7D"
      ),
      axis.title = element_text(face = "bold", size = base_size - 1),
      axis.text = element_text(color = "#2C3E50"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "#EAECEE", linewidth = 0.35),
      strip.text = element_text(face = "bold"),
      legend.title = element_blank(),
      legend.position = "bottom",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

split_cols <- c(
  "train" = "#2E86AB",
  "validation" = "#F18F01",
  "test" = "#6A994E"
)

modality_cols <- c(
  "MRI" = "#2E86AB",
  "Proteomics" = "#8E44AD",
  "Metabolomics" = "#D35400"
)

significance_cols <- c(
  "Significant" = "#2C7BB6",
  "Not significant" = "#9AA0A6",
  "Missing" = "#D0D3D4"
)

quartile_cols <- c(
  "Q1 lowest risk" = "#1B9E77",
  "Q2" = "#7570B3",
  "Q3" = "#E6AB02",
  "Q4 highest risk" = "#D95F02"
)

json_get_num <- function(x, field) {
  if (is.null(x[[field]])) return(NA_real_)
  out <- suppressWarnings(as.numeric(x[[field]]))
  if (length(out) == 0) return(NA_real_)
  out[[1]]
}

normalize_event_column <- function(x) {
  x_chr <- as.character(x)
  dplyr::case_when(
    x_chr %in% c("TRUE", "True", "true", "1", "1.0") ~ TRUE,
    x_chr %in% c("FALSE", "False", "false", "0", "0.0") ~ FALSE,
    TRUE ~ as.logical(x)
  )
}

first_existing_col <- function(df, candidates) {
  hits <- candidates[candidates %in% colnames(df)]
  if (length(hits) == 0) return(NA_character_)
  hits[[1]]
}

first_regex_col <- function(df, patterns, prefer_non_na = TRUE) {
  cols <- colnames(df)

  for (pat in patterns) {
    hits <- grep(pat, cols, value = TRUE, ignore.case = TRUE)
    if (length(hits) > 0) {
      if (prefer_non_na) {
        usable <- hits[purrr::map_lgl(hits, ~ any(is.finite(suppressWarnings(as.numeric(df[[.x]]))), na.rm = TRUE))]
        if (length(usable) > 0) return(usable[[1]])
      }
      return(hits[[1]])
    }
  }

  NA_character_
}

calc_cindex <- function(df) {
  if (nrow(df) < 5 || sum(df$event, na.rm = TRUE) < 2) return(NA_real_)

  # survival::concordance assumes larger predictor means longer survival.
  # Larger L'EPOCH risk score means higher dementia risk, so use -risk_score.
  fit <- survival::concordance(
    survival::Surv(time_years, event) ~ I(-risk_score),
    data = df
  )

  as.numeric(fit$concordance)
}

boot_cindex <- function(df, B = 500, seed = 2026) {
  set.seed(seed)
  observed <- calc_cindex(df)

  if (nrow(df) < 20 || sum(df$event, na.rm = TRUE) < 5) {
    return(
      tibble(
        cindex = observed,
        cindex_lower = NA_real_,
        cindex_upper = NA_real_,
        n_boot_ok = 0
      )
    )
  }

  vals <- replicate(B, {
    idx <- sample(seq_len(nrow(df)), size = nrow(df), replace = TRUE)
    d <- df[idx, , drop = FALSE]
    if (sum(d$event, na.rm = TRUE) < 2) return(NA_real_)
    calc_cindex(d)
  })

  vals <- vals[is.finite(vals)]

  tibble(
    cindex = observed,
    cindex_lower = as.numeric(quantile(vals, 0.025, na.rm = TRUE)),
    cindex_upper = as.numeric(quantile(vals, 0.975, na.rm = TRUE)),
    n_boot_ok = length(vals)
  )
}

km_risk_at_time <- function(df, tau) {
  if (nrow(df) < 5 || sum(df$event, na.rm = TRUE) < 1) {
    return(
      tibble(
        observed_risk = NA_real_,
        observed_lower = NA_real_,
        observed_upper = NA_real_
      )
    )
  }

  fit <- survival::survfit(
    survival::Surv(time_years, event) ~ 1,
    data = df,
    conf.type = "log-log"
  )

  ss <- summary(fit, times = tau, extend = TRUE)

  if (length(ss$surv) == 0) {
    return(
      tibble(
        observed_risk = NA_real_,
        observed_lower = NA_real_,
        observed_upper = NA_real_
      )
    )
  }

  tibble(
    observed_risk = 1 - ss$surv[[1]],
    observed_lower = 1 - ss$upper[[1]],
    observed_upper = 1 - ss$lower[[1]]
  )
}

tidy_survfit <- function(fit) {
  ss <- summary(fit)

  tibble(
    time = ss$time,
    survival = ss$surv,
    lower = ss$lower,
    upper = ss$upper,
    strata = if (is.null(ss$strata)) "All" else as.character(ss$strata)
  ) %>%
    mutate(strata = stringr::str_replace(strata, "^risk_quartile=", ""))
}

read_delta_from_files <- function(meta, perf) {
  delta <- json_get_num(perf, "delta_cindex_test_M3_vs_M1")
  lo <- json_get_num(perf, "delta_cindex_test_M3_vs_M1_ci_lower")
  hi <- json_get_num(perf, "delta_cindex_test_M3_vs_M1_ci_upper")
  p <- json_get_num(perf, "delta_cindex_test_M3_vs_M1_p_two_sided")

  if ((!is.finite(delta) || !is.finite(lo) || !is.finite(hi)) && file.exists(meta$delta_file)) {
    delta_tbl <- readr::read_tsv(
      meta$delta_file,
      show_col_types = FALSE,
      progress = FALSE
    )

    if (nrow(delta_tbl) > 0) {
      if ("delta_cindex" %in% colnames(delta_tbl)) {
        delta <- as.numeric(delta_tbl$delta_cindex[[1]])
      }
      if ("delta_cindex_ci_lower" %in% colnames(delta_tbl)) {
        lo <- as.numeric(delta_tbl$delta_cindex_ci_lower[[1]])
      }
      if ("delta_cindex_ci_upper" %in% colnames(delta_tbl)) {
        hi <- as.numeric(delta_tbl$delta_cindex_ci_upper[[1]])
      }
      if ("empirical_p_two_sided_delta_not_equal_0" %in% colnames(delta_tbl)) {
        p <- as.numeric(delta_tbl$empirical_p_two_sided_delta_not_equal_0[[1]])
      }
    }
  }

  tibble(
    delta_cindex_test_M3_vs_M1 = delta,
    delta_cindex_test_M3_vs_M1_ci_lower = lo,
    delta_cindex_test_M3_vs_M1_ci_upper = hi,
    delta_cindex_test_M3_vs_M1_p_two_sided = p
  )
}

# ============================================================
# 4. Dementia L'EPOCH column detection
# ============================================================

detect_dementia_columns <- function(pred, meta) {
  feature_token <- paste0(meta$organ_key, "_", meta$modality_key)

  expected_risk_score_col <- paste0(feature_token, "_dementia_risk_score")
  expected_accel_z_col <- paste0(feature_token, "_dementia_clock_acceleration_z")
  expected_accel_year_col <- paste0(feature_token, "_dementia_clock_acceleration_years")
  expected_clock_age_col <- paste0(feature_token, "_dementia_clock_age_years")

  # Risk score is required.
  risk_score_col <- first_existing_col(pred, c(
    expected_risk_score_col,
    paste0(feature_token, "_risk_score"),
    paste0(feature_token, "_lepoch_risk_score"),
    paste0(feature_token, "_dementia_lepoch_risk_score")
  ))

  if (is.na(risk_score_col)) {
    risk_score_col <- first_regex_col(
      pred,
      patterns = c(
        paste0("^", feature_token, ".*dementia.*risk.*score$"),
        paste0("^", feature_token, ".*lepoch.*risk.*score$"),
        paste0("^", feature_token, ".*risk.*score$"),
        "dementia.*risk.*score$",
        "lepoch.*risk.*score$",
        "risk_score$"
      )
    )
  }

  # Acceleration columns are useful but not mandatory.
  accel_z_col <- first_existing_col(pred, c(
    expected_accel_z_col,
    paste0(feature_token, "_clock_acceleration_z"),
    paste0(feature_token, "_dementia_acceleration_z"),
    paste0(feature_token, "_lepoch_acceleration_z")
  ))

  if (is.na(accel_z_col)) {
    accel_z_col <- first_regex_col(
      pred,
      patterns = c(
        paste0("^", feature_token, ".*dementia.*acceleration.*z$"),
        paste0("^", feature_token, ".*lepoch.*acceleration.*z$"),
        paste0("^", feature_token, ".*acceleration.*z$"),
        "dementia.*acceleration.*z$",
        "lepoch.*acceleration.*z$",
        "acceleration_z$"
      )
    )
  }

  accel_year_col <- first_existing_col(pred, c(
    expected_accel_year_col,
    paste0(feature_token, "_clock_acceleration_years"),
    paste0(feature_token, "_dementia_acceleration_years"),
    paste0(feature_token, "_lepoch_acceleration_years")
  ))

  if (is.na(accel_year_col)) {
    accel_year_col <- first_regex_col(
      pred,
      patterns = c(
        paste0("^", feature_token, ".*dementia.*acceleration.*years$"),
        paste0("^", feature_token, ".*lepoch.*acceleration.*years$"),
        paste0("^", feature_token, ".*acceleration.*years$"),
        "dementia.*acceleration.*years$",
        "lepoch.*acceleration.*years$",
        "acceleration_years$"
      )
    )
  }

  clock_age_col <- first_existing_col(pred, c(
    expected_clock_age_col,
    paste0(feature_token, "_clock_age_years"),
    paste0(feature_token, "_dementia_clock_age_years"),
    paste0(feature_token, "_lepoch_clock_age_years")
  ))

  if (is.na(clock_age_col)) {
    clock_age_col <- first_regex_col(
      pred,
      patterns = c(
        paste0("^", feature_token, ".*dementia.*clock.*age.*years$"),
        paste0("^", feature_token, ".*lepoch.*clock.*age.*years$"),
        paste0("^", feature_token, ".*clock.*age.*years$"),
        "dementia.*clock.*age.*years$",
        "lepoch.*clock.*age.*years$",
        "clock_age_years$"
      )
    )
  }

  list(
    risk_score_col = risk_score_col,
    accel_z_col = accel_z_col,
    accel_year_col = accel_year_col,
    clock_age_col = clock_age_col
  )
}

# ============================================================
# 5. Per-clock processing function
# ============================================================

plot_one_clock <- function(meta_row, horizon_years = 5, n_boot = 500, skip_missing = TRUE) {
  meta <- as.list(meta_row)

  if (!file.exists(meta$clock_dir) || !file.exists(meta$prediction_file) || !file.exists(meta$performance_file)) {
    msg <- paste0(
      "Missing input for ", meta$clock_label, " Dementia L'EPOCH:\n",
      "  clock_dir: ", meta$clock_dir, " exists=", file.exists(meta$clock_dir), "\n",
      "  predictions: ", meta$prediction_file, " exists=", file.exists(meta$prediction_file), "\n",
      "  performance: ", meta$performance_file, " exists=", file.exists(meta$performance_file)
    )

    if (skip_missing) {
      warning(msg)
      return(NULL)
    } else {
      stop(msg)
    }
  }

  message("\n============================================================")
  message("Processing ", meta$clock_label, " Dementia L'EPOCH")
  message("Directory: ", meta$clock_dir)
  message("============================================================")

  pred <- readr::read_tsv(
    meta$prediction_file,
    show_col_types = FALSE,
    progress = FALSE
  )

  perf <- jsonlite::fromJSON(meta$performance_file)

  detected <- detect_dementia_columns(pred, meta)

  risk_score_col <- detected$risk_score_col
  accel_z_col <- detected$accel_z_col
  accel_year_col <- detected$accel_year_col
  clock_age_col <- detected$clock_age_col

  if (is.na(risk_score_col) || !risk_score_col %in% colnames(pred)) {
    stop(
      "Could not identify Dementia L'EPOCH risk-score column for ",
      meta$clock_label,
      ". Available columns include:\n",
      paste(head(colnames(pred), 80), collapse = ", ")
    )
  }

  if (is.na(accel_z_col) || !accel_z_col %in% colnames(pred)) {
    warning(
      "Missing acceleration-z column for ",
      meta$clock_label,
      "; using standardized risk score as fallback."
    )
    accel_z_col <- ".fallback_clock_accel_z"
    pred[[accel_z_col]] <- as.numeric(scale(pred[[risk_score_col]]))
  }

  if (is.na(accel_year_col) || !accel_year_col %in% colnames(pred)) {
    accel_year_col <- ".missing_clock_accel_years"
    pred[[accel_year_col]] <- NA_real_
  }

  if (is.na(clock_age_col) || !clock_age_col %in% colnames(pred)) {
    clock_age_col <- ".missing_clock_age_years"
    pred[[clock_age_col]] <- NA_real_
  }

  age_col <- first_existing_col(
    pred,
    c(
      "age_at_imaging",
      "age_at_baseline",
      "age_at_assessment",
      "age",
      "Age",
      "Age_recruitment"
    )
  )

  if (is.na(age_col)) {
    age_col <- first_regex_col(
      pred,
      patterns = c(
        "^age_",
        "age.*baseline",
        "age.*imaging",
        "age"
      ),
      prefer_non_na = TRUE
    )
  }

  if (is.na(age_col)) {
    stop("Could not find age column for ", meta$clock_label)
  }

  required_cols <- c(
    "participant_id",
    "time_years",
    "event",
    "split",
    risk_score_col,
    accel_z_col
  )

  missing_cols <- setdiff(required_cols, colnames(pred))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns for ",
      meta$clock_label,
      ": ",
      paste(missing_cols, collapse = ", ")
    )
  }

  pred <- pred %>%
    mutate(
      event = normalize_event_column(event),
      split = factor(split, levels = c("train", "validation", "test")),
      age_at_clock = .data[[age_col]],
      risk_score = .data[[risk_score_col]],
      clock_accel_z = .data[[accel_z_col]],
      clock_accel_years = .data[[accel_year_col]],
      clock_age_years = .data[[clock_age_col]]
    ) %>%
    filter(
      !is.na(time_years),
      !is.na(event),
      !is.na(split),
      !is.na(risk_score),
      !is.na(clock_accel_z)
    )

  message("Detected risk-score column: ", risk_score_col)
  message("Detected acceleration-z column: ", accel_z_col)
  message("Detected acceleration-years column: ", accel_year_col)
  message("Detected clock-age-years column: ", clock_age_col)
  message("Detected age column: ", age_col)

  message("N total: ", nrow(pred))
  message("Dementia events total: ", sum(pred$event, na.rm = TRUE))
  message("Median follow-up: ", round(median(pred$time_years, na.rm = TRUE), 2), " years")

  # Detect usable absolute-risk column.
  local_horizon_years <- horizon_years
  risk_col_requested <- paste0("risk_", local_horizon_years, "y")
  risk_cols_available <- grep("^risk_[0-9.]+y$", colnames(pred), value = TRUE)
  risk_col <- NULL

  if (
    risk_col_requested %in% colnames(pred) &&
      any(is.finite(pred[[risk_col_requested]]), na.rm = TRUE)
  ) {
    risk_col <- risk_col_requested
  } else {
    usable_risk_cols <- risk_cols_available[
      purrr::map_lgl(
        risk_cols_available,
        ~ any(is.finite(pred[[.x]]), na.rm = TRUE)
      )
    ]

    if (length(usable_risk_cols) > 0) {
      risk_col <- usable_risk_cols[[1]]
      local_horizon_years <- as.numeric(
        stringr::str_match(risk_col, "^risk_([0-9.]+)y$")[, 2]
      )
      message("Requested absolute-risk column unavailable. Using: ", risk_col)
    } else {
      warning(
        "No usable absolute-risk column found for ",
        meta$clock_label,
        ". Calibration panel will be skipped."
      )
    }
  }

  # C-index and overfitting.
  message("Computing bootstrap C-index confidence intervals for ", meta$clock_label, "...")

  cindex_tbl <- pred %>%
    group_by(split) %>%
    group_modify(~ boot_cindex(.x, B = n_boot, seed = 2026)) %>%
    ungroup() %>%
    mutate(
      clock_id = meta$clock_id,
      modality = meta$modality,
      modality_key = meta$modality_key,
      organ_key = meta$organ_key,
      organ_label = meta$organ_label,
      clock_label = meta$clock_label,
      split = factor(split, levels = c("train", "validation", "test"))
    )

  json_cindex_tbl <- tibble(
    split = factor(
      c("train", "validation", "test"),
      levels = c("train", "validation", "test")
    ),
    cindex_json = c(
      json_get_num(perf, "cindex_train"),
      json_get_num(perf, "cindex_validation"),
      json_get_num(perf, "cindex_test")
    )
  )

  cindex_tbl <- cindex_tbl %>%
    left_join(json_cindex_tbl, by = "split")

  test_cindex <- cindex_tbl %>%
    filter(split == "test") %>%
    pull(cindex)

  overfit_tbl <- cindex_tbl %>%
    mutate(
      optimism_vs_test = cindex - test_cindex,
      n = purrr::map_int(
        as.character(split),
        ~ sum(pred$split == .x)
      ),
      events = purrr::map_int(
        as.character(split),
        ~ sum(pred$split == .x & pred$event)
      ),
      overfitting_flag = case_when(
        split == "train" & optimism_vs_test > 0.05 ~ "Possible overfitting",
        split == "train" & optimism_vs_test > 0.02 ~ "Mild optimism",
        split == "train" ~ "Low optimism",
        TRUE ~ ""
      )
    )

  # Preserve original per-clock folder, but use dementia-specific figure suffix.
  out_pdf <- file.path(meta$clock_dir, paste0(meta$prefix, "_lepoch_model_performance.pdf"))
  out_png <- file.path(meta$clock_dir, paste0(meta$prefix, "_lepoch_model_performance.png"))
  out_overfit <- file.path(meta$clock_dir, paste0(meta$prefix, "_lepoch_overfitting_summary.tsv"))
  out_quartile <- file.path(meta$clock_dir, paste0(meta$prefix, "_lepoch_risk_quartile_summary.tsv"))

  readr::write_tsv(overfit_tbl, out_overfit)

  # Risk quartiles using training-set cutoffs.
  train_risk <- pred %>%
    filter(split == "train") %>%
    pull(risk_score)

  risk_breaks <- quantile(
    train_risk,
    probs = seq(0, 1, by = 0.25),
    na.rm = TRUE,
    type = 8
  )

  risk_breaks[1] <- -Inf
  risk_breaks[length(risk_breaks)] <- Inf

  if (length(unique(risk_breaks)) < length(risk_breaks)) {
    warning(
      "Non-unique training risk quartile breaks for ",
      meta$clock_label,
      ". Using within-split quartiles."
    )

    pred <- pred %>%
      group_by(split) %>%
      mutate(
        risk_quartile = ntile(risk_score, 4),
        risk_quartile = factor(
          risk_quartile,
          levels = 1:4,
          labels = c(
            "Q1 lowest risk",
            "Q2",
            "Q3",
            "Q4 highest risk"
          )
        )
      ) %>%
      ungroup()
  } else {
    pred <- pred %>%
      mutate(
        risk_quartile = cut(
          risk_score,
          breaks = risk_breaks,
          include.lowest = TRUE,
          labels = c(
            "Q1 lowest risk",
            "Q2",
            "Q3",
            "Q4 highest risk"
          )
        )
      )
  }

  quartile_risk_tbl <- pred %>%
    filter(!is.na(risk_quartile)) %>%
    group_by(split, risk_quartile) %>%
    group_modify(~ {
      d <- .x
      km <- km_risk_at_time(d, local_horizon_years)

      tibble(
        n = nrow(d),
        events = sum(d$event, na.rm = TRUE),
        mean_risk_score = mean(d$risk_score, na.rm = TRUE),
        mean_clock_accel_z = mean(d$clock_accel_z, na.rm = TRUE)
      ) %>%
        bind_cols(km)
    }) %>%
    ungroup() %>%
    mutate(
      clock_id = meta$clock_id,
      modality = meta$modality,
      organ_key = meta$organ_key,
      organ_label = meta$organ_label,
      clock_label = meta$clock_label,
      horizon_years = local_horizon_years
    )

  readr::write_tsv(quartile_risk_tbl, out_quartile)

  # Calibration table.
  out_cal <- NULL

  if (!is.null(risk_col)) {
    cal_input <- pred %>%
      mutate(predicted_risk = .data[[risk_col]]) %>%
      filter(is.finite(predicted_risk))

    if (nrow(cal_input) > 0) {
      cal_tbl <- cal_input %>%
        group_by(split) %>%
        mutate(cal_bin = ntile(predicted_risk, 10)) %>%
        ungroup() %>%
        group_by(split, cal_bin) %>%
        group_modify(~ {
          d <- .x
          km <- km_risk_at_time(d, local_horizon_years)

          tibble(
            n = nrow(d),
            events = sum(d$event, na.rm = TRUE),
            predicted_risk = mean(d$predicted_risk, na.rm = TRUE)
          ) %>%
            bind_cols(km)
        }) %>%
        ungroup() %>%
        mutate(
          clock_id = meta$clock_id,
          modality = meta$modality,
          organ_key = meta$organ_key,
          organ_label = meta$organ_label,
          clock_label = meta$clock_label,
          horizon_years = local_horizon_years,
          risk_col = risk_col
        )

      out_cal <- file.path(
        meta$clock_dir,
        paste0(meta$prefix, "_lepoch_calibration_", local_horizon_years, "y.tsv")
      )

      readr::write_tsv(cal_tbl, out_cal)
    } else {
      cal_tbl <- tibble()
    }
  } else {
    cal_tbl <- tibble()
  }

  # Test-set Kaplan-Meier by risk quartile.
  test_for_km <- pred %>%
    filter(split == "test", !is.na(risk_quartile))

  if (sum(test_for_km$event, na.rm = TRUE) >= 2) {
    km_fit <- survival::survfit(
      survival::Surv(time_years, event) ~ risk_quartile,
      data = test_for_km
    )
    km_df <- tidy_survfit(km_fit)
  } else {
    km_df <- tibble()
    warning(
      "Too few test-set dementia events for Kaplan-Meier by quartile for ",
      meta$clock_label,
      "."
    )
  }

  time_origin_label <- ifelse(meta$modality == "MRI", "imaging", "baseline")

  # -----------------------------
  # Panels
  # -----------------------------

  p_cindex <- cindex_tbl %>%
    ggplot(aes(x = split, y = cindex, fill = split)) +
    geom_col(
      width = 0.62,
      alpha = 0.94,
      color = "white",
      linewidth = 0.4
    ) +
    geom_errorbar(
      aes(ymin = cindex_lower, ymax = cindex_upper),
      width = 0.14,
      linewidth = 0.65,
      color = "#1C2833"
    ) +
    geom_text(
      aes(label = sprintf("%.3f", cindex)),
      vjust = -0.55,
      size = 3.8,
      fontface = "bold",
      color = "#1C2833"
    ) +
    scale_fill_manual(values = split_cols) +
    coord_cartesian(
      ylim = c(
        0.45,
        max(cindex_tbl$cindex_upper, cindex_tbl$cindex, na.rm = TRUE) + 0.05
      )
    ) +
    labs(
      title = "A. Discrimination across splits",
      subtitle = glue("C-index with bootstrap 95% CI, B = {n_boot}"),
      x = NULL,
      y = "C-index"
    ) +
    theme_clock()

  gap_tbl <- overfit_tbl %>%
    filter(split %in% c("train", "validation")) %>%
    mutate(
      comparison = case_when(
        split == "train" ~ "Train - test",
        split == "validation" ~ "Validation - test",
        TRUE ~ as.character(split)
      )
    )

  p_gap <- gap_tbl %>%
    ggplot(aes(x = comparison, y = optimism_vs_test, fill = split)) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "#7B7D7D",
      linewidth = 0.5
    ) +
    geom_col(
      width = 0.55,
      alpha = 0.94,
      color = "white",
      linewidth = 0.4
    ) +
    geom_text(
      aes(label = sprintf("%+.3f", optimism_vs_test)),
      vjust = ifelse(gap_tbl$optimism_vs_test >= 0, -0.55, 1.25),
      size = 3.8,
      fontface = "bold"
    ) +
    scale_fill_manual(values = split_cols) +
    labs(
      title = "B. Optimism gap",
      subtitle = "Large positive train-test gap suggests overfitting",
      x = NULL,
      y = "C-index difference"
    ) +
    theme_clock()

  p_density <- pred %>%
    ggplot(aes(x = risk_score, fill = split, color = split)) +
    geom_density(alpha = 0.22, linewidth = 0.75) +
    scale_fill_manual(values = split_cols) +
    scale_color_manual(values = split_cols) +
    labs(
      title = "C. L'EPOCH risk-score distribution",
      subtitle = "Train, validation, and test distributions should broadly overlap",
      x = glue("{meta$clock_label} Dementia L'EPOCH risk score"),
      y = "Density"
    ) +
    theme_clock()

  if (nrow(km_df) > 0) {
    p_km <- km_df %>%
      ggplot(aes(x = time, y = 1 - survival, color = strata, fill = strata)) +
      geom_step(linewidth = 0.9) +
      geom_ribbon(
        aes(ymin = 1 - upper, ymax = 1 - lower),
        alpha = 0.11,
        color = NA
      ) +
      scale_color_manual(values = quartile_cols) +
      scale_fill_manual(values = quartile_cols) +
      scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
      labs(
        title = "D. Test-set dementia separation",
        subtitle = "Kaplan-Meier cumulative dementia incidence by training-defined risk quartile",
        x = glue("Years after {time_origin_label}"),
        y = "Cumulative dementia incidence"
      ) +
      theme_clock()
  } else {
    p_km <- ggplot() +
      annotate(
        "text",
        x = 0,
        y = 0,
        label = "Kaplan-Meier panel unavailable: too few test-set dementia events.",
        size = 4.5,
        fontface = "bold"
      ) +
      theme_void() +
      labs(title = "D. Test-set dementia separation")
  }

  p_quartile <- quartile_risk_tbl %>%
    ggplot(aes(x = risk_quartile, y = observed_risk, color = split, group = split)) +
    geom_pointrange(
      aes(ymin = observed_lower, ymax = observed_upper),
      position = position_dodge(width = 0.45),
      linewidth = 0.55,
      size = 0.45
    ) +
    geom_line(
      position = position_dodge(width = 0.45),
      linewidth = 0.65,
      alpha = 0.85
    ) +
    scale_color_manual(values = split_cols) +
    scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
    labs(
      title = glue("E. Observed {local_horizon_years}-year dementia risk by risk quartile"),
      subtitle = "Monotonic validation/test patterns support generalization",
      x = NULL,
      y = glue("Observed {local_horizon_years}-year dementia risk")
    ) +
    theme_clock() +
    theme(axis.text.x = element_text(angle = 18, hjust = 1))

  if (nrow(cal_tbl) > 0) {
    max_axis <- max(
      cal_tbl$predicted_risk,
      cal_tbl$observed_upper,
      cal_tbl$observed_risk,
      na.rm = TRUE
    )

    p_cal <- cal_tbl %>%
      ggplot(aes(x = predicted_risk, y = observed_risk, color = split)) +
      geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed",
        color = "#7B7D7D",
        linewidth = 0.55
      ) +
      geom_errorbar(
        aes(ymin = observed_lower, ymax = observed_upper),
        width = 0,
        alpha = 0.65,
        linewidth = 0.45
      ) +
      geom_point(size = 2.2, alpha = 0.9) +
      geom_smooth(
        method = "loess",
        se = FALSE,
        linewidth = 0.75,
        alpha = 0.8
      ) +
      scale_color_manual(values = split_cols) +
      scale_x_continuous(
        labels = percent_format(accuracy = 0.1),
        limits = c(0, max_axis * 1.08)
      ) +
      scale_y_continuous(
        labels = percent_format(accuracy = 0.1),
        limits = c(0, max_axis * 1.08)
      ) +
      labs(
        title = glue("F. Calibration at {local_horizon_years} years"),
        subtitle = "Observed Kaplan-Meier dementia risk versus predicted Cox risk by decile",
        x = glue("Mean predicted {local_horizon_years}-year dementia risk"),
        y = glue("Observed {local_horizon_years}-year dementia risk")
      ) +
      theme_clock()
  } else {
    p_cal <- ggplot() +
      annotate(
        "text",
        x = 0,
        y = 0,
        label = "Calibration skipped: no usable absolute-risk column.",
        size = 4.5,
        fontface = "bold"
      ) +
      xlim(-1, 1) +
      ylim(-1, 1) +
      labs(title = "F. Calibration unavailable", x = NULL, y = NULL) +
      theme_void()
  }

  admin_censor_date <- ifelse(
    is.null(perf$admin_censor_date),
    "not recorded",
    perf$admin_censor_date
  )

  subtitle_text <- glue(
    "N = {nrow(pred)}, dementia events = {sum(pred$event)}, ",
    "median follow-up = {round(median(pred$time_years, na.rm = TRUE), 2)} years; ",
    "administrative censor date = {admin_censor_date}"
  )

  final_fig <-
    (p_cindex | p_gap) /
    (p_density | p_km) /
    (p_quartile | p_cal) +
    plot_annotation(
      title = glue("{meta$clock_label} Dementia L'EPOCH performance"),
      subtitle = subtitle_text,
      caption = paste0(
        "Risk quartiles are defined from the training set and applied to validation/test. ",
        "Higher L'EPOCH risk score indicates greater dementia-event proximity."
      ),
      theme = theme(
        plot.title = element_text(face = "bold", size = 19, color = "#17202A"),
        plot.subtitle = element_text(size = 11, color = "#566573"),
        plot.caption = element_text(size = 9, color = "#7B7D7D")
      )
    )

  print(final_fig)

  if (capabilities("cairo")) {
    ggsave(
      filename = out_pdf,
      plot = final_fig,
      width = 15,
      height = 17,
      device = cairo_pdf
    )
  } else {
    ggsave(
      filename = out_pdf,
      plot = final_fig,
      width = 15,
      height = 17
    )
  }

  ggsave(
    filename = out_png,
    plot = final_fig,
    width = 15,
    height = 17,
    dpi = 320
  )

  # Incremental value.
  delta_tbl <- read_delta_from_files(meta, perf)

  m3_json_field_candidates <- c(
    paste0("cindex_test_M3_full_covariates_plus_", meta$organ_key, "_", meta$modality_key),
    paste0("cindex_test_M3_full_model_plus_", meta$organ_key, "_", meta$modality_key),
    "cindex_test_M3_full_model",
    "cindex_test_M3"
  )

  m3_cindex <- NA_real_

  for (field in m3_json_field_candidates) {
    val <- json_get_num(perf, field)
    if (is.finite(val)) {
      m3_cindex <- val
      break
    }
  }

  # If the dynamic M3 JSON field is unavailable, use the model comparison file.
  if (!is.finite(m3_cindex) && file.exists(meta$model_comparison_file)) {
    mc <- readr::read_tsv(
      meta$model_comparison_file,
      show_col_types = FALSE,
      progress = FALSE
    )

    m3_row <- mc %>%
      filter(split == "test", stringr::str_detect(model, "^M3_"))

    if (nrow(m3_row) > 0 && "cindex" %in% colnames(m3_row)) {
      m3_cindex <- as.numeric(m3_row$cindex[[1]])
    }
  }

  summary_row <- tibble(
    clock_id = meta$clock_id,
    modality = meta$modality,
    modality_key = meta$modality_key,
    organ_key = meta$organ_key,
    organ_label = meta$organ_label,
    clock_label = meta$clock_label,
    folder = meta$folder,
    prefix = meta$prefix,
    n_total = nrow(pred),
    n_events_total = sum(pred$event, na.rm = TRUE),
    median_followup_years = median(pred$time_years, na.rm = TRUE),

    cindex_train = cindex_tbl %>%
      filter(split == "train") %>%
      pull(cindex),

    cindex_validation = cindex_tbl %>%
      filter(split == "validation") %>%
      pull(cindex),

    cindex_test = cindex_tbl %>%
      filter(split == "test") %>%
      pull(cindex),

    cindex_test_lower = cindex_tbl %>%
      filter(split == "test") %>%
      pull(cindex_lower),

    cindex_test_upper = cindex_tbl %>%
      filter(split == "test") %>%
      pull(cindex_upper),

    train_minus_test = overfit_tbl %>%
      filter(split == "train") %>%
      pull(optimism_vs_test),

    validation_minus_test = overfit_tbl %>%
      filter(split == "validation") %>%
      pull(optimism_vs_test),

    cindex_test_M1_covariate_baseline = json_get_num(perf, "cindex_test_M1_covariate_baseline"),
    cindex_test_M3_full_model = m3_cindex,

    horizon_years_used = local_horizon_years,
    risk_col_used = ifelse(is.null(risk_col), NA_character_, risk_col),
    risk_score_col_used = risk_score_col,
    accel_z_col_used = accel_z_col,
    accel_year_col_used = accel_year_col,
    clock_age_col_used = clock_age_col,
    age_col_used = age_col,
    output_pdf = out_pdf,
    output_png = out_png
  ) %>%
    bind_cols(delta_tbl) %>%
    mutate(
      delta_significant = case_when(
        is.na(delta_cindex_test_M3_vs_M1) ~ "Missing",
        is.finite(delta_cindex_test_M3_vs_M1_ci_lower) &
          delta_cindex_test_M3_vs_M1_ci_lower > 0 ~ "Significant",
        is.finite(delta_cindex_test_M3_vs_M1_p_two_sided) &
          delta_cindex_test_M3_vs_M1_p_two_sided < 0.05 &
          delta_cindex_test_M3_vs_M1 > 0 ~ "Significant",
        TRUE ~ "Not significant"
      )
    )

  message("Saved figure:")
  message("  ", out_pdf)
  message("  ", out_png)
  message("Saved summary tables:")
  message("  ", out_overfit)
  message("  ", out_quartile)

  if (!is.null(out_cal)) {
    message("  ", out_cal)
  }

  list(
    meta = meta_row,
    summary = summary_row,
    overfit = overfit_tbl,
    quartile = quartile_risk_tbl,
    calibration = cal_tbl,
    cindex = cindex_tbl
  )
}

# ============================================================
# 6. Loop over all Dementia L'EPOCH clocks
# ============================================================

results <- purrr::map(
  seq_len(nrow(clock_manifest)),
  ~ plot_one_clock(
    meta_row = clock_manifest[.x, ],
    horizon_years = horizon_years_default,
    n_boot = n_boot,
    skip_missing = skip_missing
  )
)

results <- purrr::compact(results)

if (length(results) == 0) {
  stop("No Dementia L'EPOCH clocks were successfully processed.")
}

# ============================================================
# 7. Combined cross-clock outputs
# ============================================================

combined_summary <- purrr::map_dfr(results, "summary")
combined_overfit <- purrr::map_dfr(results, "overfit")
combined_quartile <- purrr::map_dfr(results, "quartile")
combined_cindex <- purrr::map_dfr(results, "cindex")
combined_calibration <- purrr::map_dfr(
  results,
  function(x) {
    if (nrow(x$calibration) == 0) tibble() else x$calibration
  }
)

out_combined_summary <- file.path(
  combined_outdir,
  "all_dementia_lepoch_summary.tsv"
)

out_combined_overfit <- file.path(
  combined_outdir,
  "all_dementia_lepoch_overfitting_summary.tsv"
)

out_combined_quartile <- file.path(
  combined_outdir,
  "all_dementia_lepoch_risk_quartile_summary.tsv"
)

out_combined_calibration <- file.path(
  combined_outdir,
  "all_dementia_lepoch_calibration_summary.tsv"
)

out_combined_cindex <- file.path(
  combined_outdir,
  "all_dementia_lepoch_cindex_summary.tsv"
)

out_combined_incremental <- file.path(
  combined_outdir,
  "all_dementia_lepoch_incremental_value_summary.tsv"
)

message("\n============================================================")
message("Writing combined Dementia L'EPOCH outputs")
message("============================================================")

readr::write_tsv(combined_summary, out_combined_summary)
readr::write_tsv(combined_overfit, out_combined_overfit)
readr::write_tsv(combined_quartile, out_combined_quartile)
readr::write_tsv(combined_cindex, out_combined_cindex)

readr::write_tsv(
  combined_summary %>%
    select(
      clock_id,
      modality,
      organ_label,
      clock_label,
      cindex_test_M1_covariate_baseline,
      cindex_test_M3_full_model,
      delta_cindex_test_M3_vs_M1,
      delta_cindex_test_M3_vs_M1_ci_lower,
      delta_cindex_test_M3_vs_M1_ci_upper,
      delta_cindex_test_M3_vs_M1_p_two_sided,
      delta_significant
    ),
  out_combined_incremental
)

if (nrow(combined_calibration) > 0) {
  readr::write_tsv(combined_calibration, out_combined_calibration)
}

message("  ", out_combined_summary)
message("  ", out_combined_overfit)
message("  ", out_combined_quartile)
message("  ", out_combined_cindex)
message("  ", out_combined_incremental)

if (nrow(combined_calibration) > 0) {
  message("  ", out_combined_calibration)
}

# ============================================================
# 8. Combined figure: discrimination, overfitting, and M3-M1
# ============================================================

clock_order_tbl <- combined_summary %>%
  mutate(
    modality = factor(
      modality,
      levels = c("MRI", "Proteomics", "Metabolomics")
    )
  ) %>%
  arrange(modality, desc(cindex_test)) %>%
  mutate(clock_plot = paste0(clock_label, "  "))

clock_levels <- rev(clock_order_tbl$clock_plot)

combined_cindex_plot_tbl <- combined_cindex %>%
  left_join(clock_order_tbl %>% select(clock_id, clock_plot), by = "clock_id") %>%
  mutate(
    clock_plot = factor(clock_plot, levels = clock_levels),
    split = factor(split, levels = c("train", "validation", "test")),
    modality = factor(
      modality,
      levels = c("MRI", "Proteomics", "Metabolomics")
    )
  )

combined_summary_plot_tbl <- combined_summary %>%
  left_join(clock_order_tbl %>% select(clock_id, clock_plot), by = "clock_id") %>%
  mutate(
    clock_plot = factor(clock_plot, levels = clock_levels),
    modality = factor(
      modality,
      levels = c("MRI", "Proteomics", "Metabolomics")
    ),
    delta_significant = factor(
      delta_significant,
      levels = c("Significant", "Not significant", "Missing")
    ),
    train_optimism_flag = case_when(
      train_minus_test > 0.05 ~ "High",
      train_minus_test > 0.02 ~ "Mild",
      TRUE ~ "Low/none"
    ),
    train_optimism_flag = factor(
      train_optimism_flag,
      levels = c("Low/none", "Mild", "High")
    )
  )

# A. C-index heatmap across train/validation/test.
p_cindex_heat <- combined_cindex_plot_tbl %>%
  ggplot(aes(x = split, y = clock_plot, fill = cindex)) +
  geom_tile(color = "white", linewidth = 0.55) +
  geom_text(
    aes(label = sprintf("%.3f", cindex)),
    size = 3.0,
    color = "#17202A"
  ) +
  scale_fill_gradient(
    low = "#EAF2F8",
    high = "#1F618D",
    limits = c(0.45, NA),
    oob = scales::squish
  ) +
  labs(
    title = "A. Discrimination across data splits",
    subtitle = "Transparent check for training, validation, and held-out test performance",
    x = NULL,
    y = NULL,
    fill = "C-index"
  ) +
  theme_clock(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 8.5),
    legend.position = "right"
  )

# B. Optimism gap.
optimism_tbl <- combined_summary_plot_tbl %>%
  select(
    clock_id,
    clock_plot,
    modality,
    train_minus_test,
    validation_minus_test
  ) %>%
  pivot_longer(
    cols = c(train_minus_test, validation_minus_test),
    names_to = "comparison",
    values_to = "optimism"
  ) %>%
  mutate(
    comparison = recode(
      comparison,
      train_minus_test = "Train - test",
      validation_minus_test = "Validation - test"
    ),
    comparison = factor(
      comparison,
      levels = c("Train - test", "Validation - test")
    )
  )

p_optimism <- optimism_tbl %>%
  ggplot(aes(x = optimism, y = clock_plot, color = comparison)) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "#7B7D7D",
    linewidth = 0.45
  ) +
  geom_vline(
    xintercept = 0.02,
    linetype = "dotted",
    color = "#B9770E",
    linewidth = 0.45
  ) +
  geom_vline(
    xintercept = 0.05,
    linetype = "dotted",
    color = "#922B21",
    linewidth = 0.45
  ) +
  geom_point(
    size = 2.0,
    alpha = 0.9,
    position = position_dodge(width = 0.55)
  ) +
  scale_color_manual(
    values = c(
      "Train - test" = "#2E86AB",
      "Validation - test" = "#F18F01"
    )
  ) +
  labs(
    title = "B. Optimism relative to test set",
    subtitle = "Dotted lines mark +0.02 and +0.05 C-index gaps",
    x = "C-index difference versus test",
    y = NULL
  ) +
  theme_clock(base_size = 11) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# C. Incremental value of omics/imaging features beyond covariates.
p_delta <- combined_summary_plot_tbl %>%
  ggplot(
    aes(
      x = delta_cindex_test_M3_vs_M1,
      y = clock_plot,
      color = delta_significant
    )
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "#7B7D7D",
    linewidth = 0.45
  ) +
  geom_errorbarh(
    aes(
      xmin = delta_cindex_test_M3_vs_M1_ci_lower,
      xmax = delta_cindex_test_M3_vs_M1_ci_upper
    ),
    height = 0.18,
    linewidth = 0.55,
    na.rm = TRUE
  ) +
  geom_point(
    aes(size = cindex_test),
    alpha = 0.95,
    na.rm = TRUE
  ) +
  scale_color_manual(values = significance_cols, drop = FALSE) +
  scale_size_continuous(
    range = c(1.8, 4.2),
    limits = range(combined_summary_plot_tbl$cindex_test, na.rm = TRUE)
  ) +
  labs(
    title = "C. Incremental value beyond covariates",
    subtitle = "Test-set ΔC-index = M3 full model - M1 covariate baseline; bars show bootstrap 95% CI",
    x = "ΔC-index on test set",
    y = NULL,
    size = "Test C-index"
  ) +
  theme_clock(base_size = 11) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# D. Practical quality map.
base_quality <- combined_summary_plot_tbl %>%
  ggplot(
    aes(
      x = train_minus_test,
      y = delta_cindex_test_M3_vs_M1,
      color = modality
    )
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "#BDC3C7",
    linewidth = 0.45
  ) +
  geom_vline(
    xintercept = 0.05,
    linetype = "dotted",
    color = "#922B21",
    linewidth = 0.45
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "#BDC3C7",
    linewidth = 0.45
  ) +
  geom_point(
    aes(size = cindex_test, shape = delta_significant),
    alpha = 0.9
  ) +
  scale_color_manual(values = modality_cols) +
  scale_shape_manual(
    values = c(
      "Significant" = 16,
      "Not significant" = 1,
      "Missing" = 4
    ),
    drop = FALSE
  ) +
  scale_size_continuous(range = c(2.2, 5.0)) +
  labs(
    title = "D. Practical Dementia L'EPOCH model-quality map",
    subtitle = "Ideal models appear near low optimism and positive ΔC-index",
    x = "Train - test C-index gap",
    y = "Test-set ΔC-index, M3 - M1",
    color = "Modality",
    shape = "Incremental value",
    size = "Test C-index"
  ) +
  theme_clock(base_size = 11)

if (requireNamespace("ggrepel", quietly = TRUE)) {
  p_quality <- base_quality +
    ggrepel::geom_text_repel(
      aes(label = organ_label),
      size = 2.7,
      max.overlaps = 30,
      min.segment.length = 0,
      box.padding = 0.25,
      show.legend = FALSE
    )
} else {
  warning("Package ggrepel is not installed. Panel D will use geom_text with check_overlap instead.")

  p_quality <- base_quality +
    geom_text(
      aes(label = organ_label),
      size = 2.4,
      check_overlap = TRUE,
      vjust = -0.7,
      show.legend = FALSE
    )
}

combined_fig <-
  (p_cindex_heat | p_optimism | p_delta) /
  p_quality +
  plot_layout(heights = c(2.3, 1.2)) +
  plot_annotation(
    title = "Performance summary of MRI, proteomics, and metabolomics Dementia L'EPOCH models",
    subtitle = glue(
      "Processed {nrow(combined_summary)} of {nrow(clock_manifest)} expected Dementia L'EPOCH models. ",
      "Panels show split-wise discrimination, overfitting risk, and incremental value beyond covariates."
    ),
    caption = "M1 = covariate baseline. M3 = covariates plus imaging/protein/metabolite features. Risk quartiles are training-defined and applied to validation/test.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 20, color = "#17202A"),
      plot.subtitle = element_text(size = 11, color = "#566573"),
      plot.caption = element_text(size = 9, color = "#7B7D7D")
    )
  )

out_combined_pdf <- file.path(
  combined_outdir,
  "all_dementia_lepoch_combined_model_performance.pdf"
)

out_combined_png <- file.path(
  combined_outdir,
  "all_dementia_lepoch_combined_model_performance.png"
)

if (capabilities("cairo")) {
  ggsave(
    filename = out_combined_pdf,
    plot = combined_fig,
    width = 18,
    height = 13,
    device = cairo_pdf
  )
} else {
  ggsave(
    filename = out_combined_pdf,
    plot = combined_fig,
    width = 18,
    height = 13
  )
}

ggsave(
  filename = out_combined_png,
  plot = combined_fig,
  width = 18,
  height = 13,
  dpi = 320
)

message("\n===== Combined Dementia L'EPOCH summary =====")

print(
  combined_summary %>%
    arrange(modality, desc(cindex_test)) %>%
    select(
      modality,
      organ_label,
      n_total,
      n_events_total,
      median_followup_years,
      cindex_train,
      cindex_validation,
      cindex_test,
      train_minus_test,
      validation_minus_test,
      delta_cindex_test_M3_vs_M1,
      delta_cindex_test_M3_vs_M1_ci_lower,
      delta_cindex_test_M3_vs_M1_ci_upper,
      delta_significant
    )
)

message("\nInterpretation guide:")
message("  - Similar train/validation/test C-index values suggest limited overfitting.")
message("  - Large positive train-test gaps, especially >0.05, suggest overfitting or split instability.")
message("  - Positive M3-M1 delta indicates that imaging/protein/metabolite features add value beyond covariates.")
message("  - A bootstrap CI for M3-M1 entirely above zero provides stronger evidence of incremental value.")
message("  - Per-clock PDF/PNG figures and summary TSVs are saved in each Dementia L'EPOCH folder.")
message("  - Combined figure saved to:")
message("    ", out_combined_pdf)
message("    ", out_combined_png)