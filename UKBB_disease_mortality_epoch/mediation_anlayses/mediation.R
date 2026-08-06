#!/usr/bin/env Rscript

# ==============================================================================
# Organ- and omics-specific statistical mediation analysis
# Baseline molecular disease EPOCH -> MRI mortality EPOCH -> observed mortality
#
# Primary analysis:
#   Exposure  : one baseline proteomics/metabolomics disease EPOCH acceleration-z
#   Mediator  : one MRI mortality EPOCH acceleration-z measured at imaging visit 2
#   Outcome   : 5-year all-cause mortality after the MRI visit (binary)
#   Estimator : lavaan WLSMV with a probit link for the ordered binary outcome
#
# Why a fixed 5-year outcome?
#   Standard lavaan SEM does not natively accommodate right-censored survival
#   outcomes. Defining mortality within a common fixed horizon preserves temporal
#   ordering and avoids treating participants censored before 5 years as known
#   survivors. Participants are eligible if they died within 5 years or had at
#   least 5 years of potential follow-up after MRI.
#
# Interpretation:
#   This is statistical mediation, not causal mediation. The MRI mortality EPOCH
#   was trained to forecast mortality, and unmeasured confounding may remain.
# ==============================================================================

.libPaths('/cbica/home/wenju/R/x86_64-pc-linux-gnu-library/4.3') ### this is to run on cluster, which should work now!

suppressPackageStartupMessages({
  library(lavaan)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(tibble)
})

# ------------------------------------------------------------------------------
# 1. User settings
# ------------------------------------------------------------------------------

root_dir <- "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock"

epoch_wide_file <- file.path(
  root_dir,
  "collected_significant_epoch_clocks",
  "significant_epoch_clocks_wide.tsv"
)

covariate_file <- paste0(
  "/cbica/home/wenju/Reproducibile_paper/PRS_UKBB/",
  "prediction/data/UKBB_fullsample_covariate.csv"
)

output_dir <- file.path(
  root_dir,
  "SEM_molecular_disease_EPOCH_to_MRI_mortality_EPOCH"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(covariate_file)) {
  stop("Covariate file not found: ", covariate_file)
}


# Fixed mortality horizon after MRI. Change to 10 for a 10-year sensitivity run.
mortality_horizon_years <- 5

# Minimum analysis size. Models below either threshold are skipped.
minimum_n <- 500
minimum_deaths <- 20

# Analysis mode:
#   "all_pairs"   = every molecular disease EPOCH x every MRI mortality EPOCH
#   "matched_only"= only mappings in organ_mapping below
analysis_mode <- "all_pairs"

# Full sample is requested. No train/test restriction is applied.
# Set to TRUE only for a sensitivity analysis restricted to mortality-clock test
# participants when the prediction file has a split column.
restrict_mortality_clock_test_set <- FALSE

# Complete-case SEM is used within each exposure-mediator pair.
# WLSMV handles the binary ordered mortality outcome.

# Covariates measured at or near imaging visit 2.
continuous_covariates <- c(
  "age_when_attended_assessment_centre_f21003_2_0",
  "body_mass_index_bmi_f23104_2_0",
  "genetic_principal_components_f22009_0_1",
  "genetic_principal_components_f22009_0_2",
  "genetic_principal_components_f22009_0_3",
  "genetic_principal_components_f22009_0_4",
  "genetic_principal_components_f22009_0_5",
  "genetic_principal_components_f22009_0_6",
  "genetic_principal_components_f22009_0_7",
  "genetic_principal_components_f22009_0_8",
  "genetic_principal_components_f22009_0_9",
  "genetic_principal_components_f22009_0_10"
)

categorical_covariates <- c(
  "sex_f31_0_0",
  "smoking_status_f20116_2_0",
  "uk_biobank_assessment_centre_f54_2_0"
)

# Optional organ mapping for analysis_mode == "matched_only".
# The mapping is intentionally explicit because some molecular systems do not
# have a one-to-one anatomical MRI counterpart.
organ_mapping <- tribble(
  ~exposure_organ,       ~mediator_organ,
  "brain",              "brain",
  "heart",              "heart",
  "hepatic",            "liver",
  "metabolic",          "pancreas",
  "endocrine",          "pancreas",
  "digestive",          "pancreas"
)

# ------------------------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------------------------

read_epoch_table <- function(path) {
  path <- path.expand(path)

  if (file.exists(path)) {
    if (grepl("\\.zip$", path, ignore.case = TRUE)) {
      stop("ZIP archives are not supported directly: ", path)
    }
    message("Reading uncompressed EPOCH TSV: ", path)
    return(read_tsv(
      path,
      show_col_types = FALSE,
      progress = FALSE,
      na = c("", "NA", "NaN", ".", "-9999")
    ))
  }

  gz_path <- paste0(path, ".gz")
  if (file.exists(gz_path)) {
    message("Uncompressed TSV not found; reading gzip-compressed EPOCH TSV: ", gz_path)
    return(read_tsv(
      gz_path,
      show_col_types = FALSE,
      progress = FALSE,
      na = c("", "NA", "NaN", ".", "-9999")
    ))
  }

  stop(
    "EPOCH table not found. Checked:\n  ", path,
    "\n  ", gz_path
  )
}

extract_epoch_metadata <- function(column_name) {
  # Expected standardized form from the collector:
  # epoch__brain__proteomics__dementia__acceleration_z
  pieces <- str_split(column_name, "__", simplify = TRUE)
  if (ncol(pieces) != 5 || pieces[1] != "epoch") {
    return(tibble(
      column = column_name,
      organ = NA_character_,
      modality = NA_character_,
      endpoint = NA_character_,
      measure = NA_character_
    ))
  }
  tibble(
    column = column_name,
    organ = pieces[2],
    modality = pieces[3],
    endpoint = pieces[4],
    measure = pieces[5]
  )
}

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

zscore_safe <- function(x) {
  x <- safe_numeric(x)
  sx <- sd(x, na.rm = TRUE)
  if (!is.finite(sx) || sx <= 0) return(rep(NA_real_, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / sx)
}

find_mortality_prediction_file <- function(root, mediator_column) {
  meta <- extract_epoch_metadata(mediator_column)
  organ <- meta$organ[[1]]

  directory_name <- paste0(organ, "_mri_mortality_clock")
  candidate_dir <- file.path(root, directory_name)

  if (!dir.exists(candidate_dir)) {
    dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
    normalized <- function(x) str_replace_all(tolower(basename(x)), "[^a-z0-9]+", "_")
    target <- str_replace_all(tolower(directory_name), "[^a-z0-9]+", "_")
    matched <- dirs[normalized(dirs) == target]
    if (length(matched) != 1) {
      stop("Could not uniquely resolve MRI mortality directory for ", mediator_column)
    }
    candidate_dir <- matched[[1]]
  }

  files <- list.files(
    candidate_dir,
    pattern = "_clock_predictions\\.tsv$",
    full.names = TRUE
  )

  if (length(files) != 1) {
    stop(
      "Expected exactly one mortality prediction TSV in ", candidate_dir,
      "; found ", length(files)
    )
  }

  files[[1]]
}

prepare_mortality_outcome <- function(prediction_file, horizon_years, restrict_test = FALSE) {
  dat <- read_tsv(
    prediction_file,
    show_col_types = FALSE,
    progress = FALSE,
    na = c("", "NA", "NaN", ".", "-9999")
  )

  # MRI mortality-clock files use imaging_date. Some older/non-MRI files use
  # sample_date. Accept either name, but always treat it as the landmark date
  # from which post-MRI mortality follow-up begins.
  landmark_candidates <- c("imaging_date", "sample_date")
  landmark_columns <- intersect(landmark_candidates, names(dat))

  if (length(landmark_columns) == 0) {
    stop(
      "No MRI landmark-date field found in ", prediction_file,
      ". Expected one of: ", paste(landmark_candidates, collapse = ", "),
      ". Available columns: ", paste(names(dat), collapse = ", ")
    )
  }

  # Prefer imaging_date for MRI files if both columns happen to be present.
  landmark_column <- landmark_columns[[1]]

  required <- c("participant_id", "death_date", "admin_censor_date")
  missing_required <- setdiff(required, names(dat))
  if (length(missing_required) > 0) {
    stop(
      "Missing required mortality fields in ", prediction_file, ": ",
      paste(missing_required, collapse = ", ")
    )
  }

  if (restrict_test) {
    if (!"split" %in% names(dat)) {
      warning(
        "restrict_test=TRUE, but no split column was found in ",
        prediction_file, "; using the full available sample."
      )
    } else {
      dat <- dat %>% filter(tolower(trimws(as.character(split))) == "test")
    }
  }

  # Rename the selected MRI date to a common internal name. Using .data[[...]]
  # avoids assuming that every MRI organ file uses the same original header.
  dat <- dat %>%
    transmute(
      participant_id = as.character(participant_id),
      landmark_date = as.Date(.data[[landmark_column]]),
      death_date = as.Date(death_date),
      admin_censor_date = as.Date(admin_censor_date),
      split = if ("split" %in% names(dat)) as.character(split) else NA_character_
    ) %>%
    filter(!is.na(participant_id), !is.na(landmark_date), !is.na(admin_censor_date)) %>%
    mutate(
      # A valid observed death must occur after MRI and no later than the
      # administrative censoring date.
      death_after_landmark = !is.na(death_date) & death_date > landmark_date,
      death_before_admin = death_after_landmark & death_date <= admin_censor_date,
      observed_end_date = if_else(
        death_before_admin,
        death_date,
        admin_censor_date
      ),
      followup_years = as.numeric(observed_end_date - landmark_date) / 365.25,
      death_time_years = if_else(
        death_before_admin,
        as.numeric(death_date - landmark_date) / 365.25,
        NA_real_
      ),
      died_within_horizon = death_before_admin &
        !is.na(death_time_years) &
        death_time_years <= horizon_years,
      observed_through_horizon = followup_years >= horizon_years,
      eligible_fixed_horizon = died_within_horizon | observed_through_horizon,
      mortality_horizon = case_when(
        died_within_horizon ~ 1L,
        observed_through_horizon ~ 0L,
        TRUE ~ NA_integer_
      )
    ) %>%
    filter(
      is.finite(followup_years),
      followup_years > 0,
      eligible_fixed_horizon,
      !is.na(mortality_horizon)
    ) %>%
    select(
      participant_id,
      landmark_date,
      mortality_horizon,
      followup_years,
      death_time_years,
      split
    )

  # One participant should contribute one MRI landmark record per organ file.
  # If duplicates exist, retain the earliest valid imaging landmark.
  if (anyDuplicated(dat$participant_id) > 0) {
    warning(
      "Duplicated participant IDs found in ", prediction_file,
      "; retaining the earliest landmark_date per participant."
    )
    dat <- dat %>%
      arrange(participant_id, landmark_date) %>%
      distinct(participant_id, .keep_all = TRUE)
  }

  message(
    "Mortality outcome source: ", prediction_file,
    " | landmark field: ", landmark_column,
    " | eligible N=", nrow(dat),
    " | deaths within ", horizon_years, " years=",
    sum(dat$mortality_horizon == 1L, na.rm = TRUE)
  )

  dat
}

make_dummy_covariates <- function(data, categorical_vars) {
  present <- intersect(categorical_vars, names(data))
  if (length(present) == 0) return(data)

  for (v in present) {
    data[[v]] <- factor(data[[v]])
  }

  mm <- model.matrix(
    as.formula(paste("~", paste(present, collapse = " + "))),
    data = data,
    na.action = na.pass
  )

  mm <- as.data.frame(mm, check.names = TRUE)
  mm$`(Intercept)` <- NULL

  data <- data %>% select(-all_of(present))
  bind_cols(data, mm)
}

extract_parameter <- function(pe, label) {
  row <- pe %>% filter(.data$label == !!label)
  if (nrow(row) == 0) {
    return(tibble(
      label = label,
      est = NA_real_, se = NA_real_, z = NA_real_, pvalue = NA_real_,
      ci.lower = NA_real_, ci.upper = NA_real_, std.all = NA_real_
    ))
  }
  row %>%
    slice(1) %>%
    transmute(
      label,
      est,
      se,
      z,
      pvalue,
      ci.lower,
      ci.upper,
      std.all
    )
}

fit_one_sem <- function(
    merged_data,
    exposure_column,
    mediator_column,
    continuous_covariates,
    categorical_covariates,
    minimum_n,
    minimum_deaths
) {
  model_data <- merged_data %>%
    transmute(
      participant_id,
      exposure_raw = .data[[exposure_column]],
      mediator_raw = .data[[mediator_column]],
      mortality_horizon = mortality_horizon,
      across(any_of(c(continuous_covariates, categorical_covariates)))
    )

  present_continuous <- intersect(continuous_covariates, names(model_data))
  present_categorical <- intersect(categorical_covariates, names(model_data))

  model_data <- model_data %>%
    mutate(
      exposure = zscore_safe(exposure_raw),
      mediator = zscore_safe(mediator_raw),
      across(all_of(present_continuous), safe_numeric)
    ) %>%
    select(-exposure_raw, -mediator_raw)

  model_data <- make_dummy_covariates(model_data, present_categorical)

  candidate_covariates <- setdiff(
    names(model_data),
    c("participant_id", "exposure", "mediator", "mortality_horizon")
  )

  # Remove covariates with no variation or all missing values.
  valid_covariates <- candidate_covariates[
    vapply(
      model_data[candidate_covariates],
      function(x) {
        x <- safe_numeric(x)
        sum(is.finite(x)) > 1 && sd(x, na.rm = TRUE) > 0
      },
      logical(1)
    )
  ]

  model_data <- model_data %>%
    select(
      participant_id,
      exposure,
      mediator,
      mortality_horizon,
      all_of(valid_covariates)
    ) %>%
    drop_na()

  n_total <- nrow(model_data)
  n_deaths <- sum(model_data$mortality_horizon == 1L)
  n_survivors <- sum(model_data$mortality_horizon == 0L)

  if (n_total < minimum_n || n_deaths < minimum_deaths) {
    return(list(
      status = "skipped_low_sample",
      n = n_total,
      deaths = n_deaths,
      survivors = n_survivors,
      fit = NULL,
      estimates = NULL,
      fit_measures = NULL,
      data = model_data,
      message = paste0("N=", n_total, ", deaths=", n_deaths)
    ))
  }

  # Treat mortality as an ordered binary outcome, giving a probit regression.
  model_data$mortality_horizon <- ordered(model_data$mortality_horizon)

  covariate_text <- if (length(valid_covariates) > 0) {
    paste(valid_covariates, collapse = " + ")
  } else {
    ""
  }

  mediator_rhs <- paste(c("a*exposure", valid_covariates), collapse = " + ")
  outcome_rhs <- paste(c("b*mediator", "cprime*exposure", valid_covariates), collapse = " + ")

  model_text <- paste0(
    "mediator ~ ", mediator_rhs, "\n",
    "mortality_horizon ~ ", outcome_rhs, "\n",
    "indirect := a*b\n",
    "direct := cprime\n",
    "total := cprime + (a*b)\n",
    "prop_mediated := indirect/total\n"
  )

  fit <- tryCatch(
    sem(
      model = model_text,
      data = model_data,
      ordered = "mortality_horizon",
      estimator = "WLSMV",
      parameterization = "theta",
      fixed.x = FALSE,
      missing = "listwise",
      meanstructure = TRUE
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(list(
      status = "model_error",
      n = n_total,
      deaths = n_deaths,
      survivors = n_survivors,
      fit = NULL,
      estimates = NULL,
      fit_measures = NULL,
      data = model_data,
      message = conditionMessage(fit)
    ))
  }

  converged <- lavInspect(fit, "converged")
  if (!isTRUE(converged)) {
    return(list(
      status = "not_converged",
      n = n_total,
      deaths = n_deaths,
      survivors = n_survivors,
      fit = fit,
      estimates = NULL,
      fit_measures = NULL,
      data = model_data,
      message = "lavaan did not converge"
    ))
  }

  pe <- parameterEstimates(
    fit,
    standardized = TRUE,
    ci = TRUE
  )

  requested_labels <- c(
    "a", "b", "cprime", "indirect", "direct", "total", "prop_mediated"
  )

  estimates <- map_dfr(requested_labels, ~ extract_parameter(pe, .x))

  fit_names <- c(
    "chisq.scaled", "df.scaled", "pvalue.scaled",
    "cfi.scaled", "tli.scaled", "rmsea.scaled", "srmr"
  )

  fm <- tryCatch(
    fitMeasures(fit, fit_names),
    error = function(e) setNames(rep(NA_real_, length(fit_names)), fit_names)
  )

  fit_measures <- as_tibble_row(as.list(fm))

  list(
    status = "ok",
    n = n_total,
    deaths = n_deaths,
    survivors = n_survivors,
    fit = fit,
    estimates = estimates,
    fit_measures = fit_measures,
    data = model_data,
    message = NA_character_
  )
}

# ------------------------------------------------------------------------------
# 3. Load and classify the 69 EPOCH variables
# ------------------------------------------------------------------------------

epoch <- read_epoch_table(epoch_wide_file) %>%
  mutate(participant_id = as.character(participant_id))

if (anyDuplicated(epoch$participant_id) > 0) {
  stop("EPOCH wide table contains duplicated participant_id values")
}

if (!"participant_id" %in% names(epoch)) {
  stop("The EPOCH wide table must contain participant_id")
}

epoch_columns <- names(epoch)[str_detect(names(epoch), "^epoch__")]
if (length(epoch_columns) == 0) {
  stop("No standardized epoch__... columns were found in the wide table")
}

metadata <- map_dfr(epoch_columns, extract_epoch_metadata) %>%
  filter(!is.na(organ), measure == "acceleration_z")

exposure_metadata <- metadata %>%
  filter(
    modality %in% c("proteomics", "metabolomics"),
    endpoint != "mortality"
  ) %>%
  arrange(organ, modality, endpoint)

mediator_metadata <- metadata %>%
  filter(
    modality == "mri",
    endpoint == "mortality"
  ) %>%
  arrange(organ)

if (nrow(exposure_metadata) == 0) {
  stop("No baseline proteomics/metabolomics disease EPOCH exposures found")
}
if (nrow(mediator_metadata) == 0) {
  stop("No MRI mortality EPOCH mediators found")
}

message("Molecular disease EPOCH exposures: ", nrow(exposure_metadata))
message("MRI mortality EPOCH mediators: ", nrow(mediator_metadata))
message("EPOCH table dimensions: ", nrow(epoch), " participants x ", ncol(epoch), " columns")
message("Exposure metadata columns: ", paste(names(exposure_metadata), collapse = ", "))
message("Mediator metadata columns: ", paste(names(mediator_metadata), collapse = ", "))

# ------------------------------------------------------------------------------
# 4. Load covariates
# ------------------------------------------------------------------------------

covariates <- read_csv(
  covariate_file,
  show_col_types = FALSE,
  na = c("", "NA", ".", "-9999")
) %>%
  rename(participant_id = eid) %>%
  mutate(participant_id = as.character(participant_id))

requested_covariates <- c(continuous_covariates, categorical_covariates)
missing_covariates <- setdiff(requested_covariates, names(covariates))
if (length(missing_covariates) > 0) {
  warning(
    "These requested covariates were not available and will be omitted: ",
    paste(missing_covariates, collapse = ", ")
  )
}

covariates <- covariates %>%
  select(participant_id, any_of(requested_covariates))

if (anyDuplicated(covariates$participant_id) > 0) {
  warning("Covariate file contains duplicated participant_id values; keeping the first row per participant")
  covariates <- covariates %>% distinct(participant_id, .keep_all = TRUE)
}

# ------------------------------------------------------------------------------
# 5. Build the exposure-mediator model grid
# ------------------------------------------------------------------------------

# Keep only uniquely named columns before constructing the Cartesian product.
# Both metadata tables contain a column called `measure`; passing them directly to
# tidyr::crossing() causes the duplicated-name error seen in the cluster log.
exposure_grid_metadata <- exposure_metadata %>%
  transmute(
    exposure_column = column,
    exposure_organ = organ,
    exposure_modality = modality,
    exposure_endpoint = endpoint,
    exposure_measure = measure
  )

mediator_grid_metadata <- mediator_metadata %>%
  transmute(
    mediator_column = column,
    mediator_organ = organ,
    mediator_modality = modality,
    mediator_endpoint = endpoint,
    mediator_measure = measure
  )

model_grid <- tidyr::crossing(
  exposure_grid_metadata,
  mediator_grid_metadata
)

if (analysis_mode == "matched_only") {
  model_grid <- model_grid %>%
    inner_join(
      organ_mapping,
      by = c("exposure_organ", "mediator_organ")
    )
} else if (analysis_mode != "all_pairs") {
  stop("analysis_mode must be either all_pairs or matched_only")
}

if (anyDuplicated(names(model_grid)) > 0) {
  stop(
    "Internal error: duplicated model_grid columns: ",
    paste(names(model_grid)[duplicated(names(model_grid))], collapse = ", ")
  )
}

if (nrow(model_grid) == 0) {
  stop("No exposure-mediator models remain after applying analysis_mode")
}

message("SEM models requested: ", nrow(model_grid))

# ------------------------------------------------------------------------------
# 6. Prepare one observed-mortality landmark dataset per MRI mediator organ
# ------------------------------------------------------------------------------

mortality_data_by_mediator <- list()

for (i in seq_len(nrow(mediator_metadata))) {
  mediator_col <- mediator_metadata$column[[i]]
  prediction_file <- find_mortality_prediction_file(root_dir, mediator_col)

  outcome_data <- prepare_mortality_outcome(
    prediction_file = prediction_file,
    horizon_years = mortality_horizon_years,
    restrict_test = restrict_mortality_clock_test_set
  )

  mortality_data_by_mediator[[mediator_col]] <- outcome_data

  message(
    "Prepared ", mortality_horizon_years, "-year mortality for ", mediator_col,
    ": N=", nrow(outcome_data),
    ", deaths=", sum(outcome_data$mortality_horizon == 1L)
  )
}

# ------------------------------------------------------------------------------
# 7. Run all organ- and omics-specific SEMs
# ------------------------------------------------------------------------------

all_estimates <- list()
all_fit_measures <- list()
all_model_qc <- list()

for (i in seq_len(nrow(model_grid))) {
  row <- model_grid[i, ]

  exposure_col <- row$exposure_column[[1]]
  mediator_col <- row$mediator_column[[1]]

  message(
    "[", i, "/", nrow(model_grid), "] ",
    row$exposure_organ[[1]], " ", row$exposure_modality[[1]], " ",
    row$exposure_endpoint[[1]], " -> ",
    row$mediator_organ[[1]], " MRI mortality EPOCH -> ",
    mortality_horizon_years, "-year mortality"
  )

  outcome_data <- mortality_data_by_mediator[[mediator_col]]

  merged <- epoch %>%
    select(
      participant_id,
      all_of(exposure_col),
      all_of(mediator_col)
    ) %>%
    inner_join(outcome_data, by = "participant_id") %>%
    left_join(covariates, by = "participant_id")

  result <- fit_one_sem(
    merged_data = merged,
    exposure_column = exposure_col,
    mediator_column = mediator_col,
    continuous_covariates = continuous_covariates,
    categorical_covariates = categorical_covariates,
    minimum_n = minimum_n,
    minimum_deaths = minimum_deaths
  )

  model_id <- paste(
    row$exposure_organ[[1]],
    row$exposure_modality[[1]],
    row$exposure_endpoint[[1]],
    "to",
    row$mediator_organ[[1]],
    "mri_mortality",
    sep = "__"
  )

  qc <- tibble(
    model_id = model_id,
    exposure_column = exposure_col,
    exposure_organ = row$exposure_organ[[1]],
    exposure_modality = row$exposure_modality[[1]],
    exposure_endpoint = row$exposure_endpoint[[1]],
    mediator_column = mediator_col,
    mediator_organ = row$mediator_organ[[1]],
    mortality_horizon_years = mortality_horizon_years,
    N = result$n,
    N_deaths = result$deaths,
    N_survivors = result$survivors,
    status = result$status,
    message = result$message
  )
  all_model_qc[[length(all_model_qc) + 1]] <- qc

  if (result$status == "ok") {
    estimates_out <- result$estimates %>%
      mutate(
        model_id = model_id,
        exposure_column = exposure_col,
        exposure_organ = row$exposure_organ[[1]],
        exposure_modality = row$exposure_modality[[1]],
        exposure_endpoint = row$exposure_endpoint[[1]],
        mediator_column = mediator_col,
        mediator_organ = row$mediator_organ[[1]],
        mortality_horizon_years = mortality_horizon_years,
        N = result$n,
        N_deaths = result$deaths,
        N_survivors = result$survivors,
        .before = 1
      )

    fit_out <- result$fit_measures %>%
      mutate(
        model_id = model_id,
        exposure_column = exposure_col,
        exposure_organ = row$exposure_organ[[1]],
        exposure_modality = row$exposure_modality[[1]],
        exposure_endpoint = row$exposure_endpoint[[1]],
        mediator_column = mediator_col,
        mediator_organ = row$mediator_organ[[1]],
        mortality_horizon_years = mortality_horizon_years,
        N = result$n,
        N_deaths = result$deaths,
        N_survivors = result$survivors,
        .before = 1
      )

    all_estimates[[length(all_estimates) + 1]] <- estimates_out
    all_fit_measures[[length(all_fit_measures) + 1]] <- fit_out
  }
}

# ------------------------------------------------------------------------------
# 8. Save outputs
# ------------------------------------------------------------------------------

estimates_table <- bind_rows(all_estimates)
fit_table <- bind_rows(all_fit_measures)
qc_table <- bind_rows(all_model_qc)

write_tsv(
  estimates_table,
  file.path(
    output_dir,
    paste0("SEM_path_estimates_", mortality_horizon_years, "y_mortality.tsv")
  )
)

write_tsv(
  fit_table,
  file.path(
    output_dir,
    paste0("SEM_fit_measures_", mortality_horizon_years, "y_mortality.tsv")
  )
)

write_tsv(
  qc_table,
  file.path(
    output_dir,
    paste0("SEM_model_QC_", mortality_horizon_years, "y_mortality.tsv")
  )
)

write_tsv(
  exposure_metadata,
  file.path(output_dir, "exposure_metadata.tsv")
)

write_tsv(
  mediator_metadata,
  file.path(output_dir, "mediator_metadata.tsv")
)

# Compact table with one row per model for the principal mediation quantities.
if (nrow(estimates_table) > 0) {
  compact <- estimates_table %>%
    filter(label %in% c("a", "b", "cprime", "indirect", "total", "prop_mediated")) %>%
    select(
      model_id,
      exposure_organ,
      exposure_modality,
      exposure_endpoint,
      mediator_organ,
      mortality_horizon_years,
      N,
      N_deaths,
      label,
      est,
      se,
      pvalue,
      ci.lower,
      ci.upper,
      std.all
    ) %>%
    pivot_wider(
      names_from = label,
      values_from = c(est, se, pvalue, ci.lower, ci.upper, std.all),
      names_glue = "{label}_{.value}"
    ) %>%
    mutate(
      indirect_FDR = p.adjust(indirect_pvalue, method = "BH"),
      a_FDR = p.adjust(a_pvalue, method = "BH"),
      b_FDR = p.adjust(b_pvalue, method = "BH")
    )

  write_tsv(
    compact,
    file.path(
      output_dir,
      paste0("SEM_compact_results_", mortality_horizon_years, "y_mortality.tsv")
    )
  )
}

message("\nCompleted.")
message("Output directory: ", output_dir)
message("Models requested: ", nrow(model_grid))
message("Models successfully fitted: ", sum(qc_table$status == "ok"))
message("Models skipped or failed: ", sum(qc_table$status != "ok"))