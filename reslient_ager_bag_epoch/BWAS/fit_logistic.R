#!/usr/bin/env Rscript

###############################################################################
# CRA vs comparator brain-wide association using logistic regression
#
# Primary intended analysis:
#   Candidate resilient agers (CRA) versus concordant unfavorable agers (CUA)
#   using brain-proteomics BAG-EPOCH discordance phenotypes.
#
# IMPORTANT FIX
# -------------
# The participant-level resilience TSV stores:
#
#   CRA_p20 / CUA_p20
#       logical indicator columns (TRUE/FALSE)
#
# and:
#
#   aging_phenotype_p20
#       "Candidate_resilient_ager"
#       "Concordant_unfavorable_ager"
#       "Concordant_favorable_ager"
#       "Latent_vulnerability_ager"
#       "Other"
#
# Therefore, the short labels "CRA" and "CUA" are NOT literal values of
# aging_phenotype_p20.
#
# This script now handles that automatically:
#
#   1) Preferred:
#      If phenotype_col = aging_phenotype_p20 and CRA_p20 / CUA_p20 exist,
#      it constructs the binary phenotype directly from those indicator columns:
#
#          CRA_p20 == TRUE -> case_binary = 1
#          CUA_p20 == TRUE -> case_binary = 0
#
#      Everyone else is excluded.
#
#   2) Fallback:
#      If the indicator columns are unavailable, it recognizes both abbreviated
#      and full phenotype labels, e.g.
#
#          CRA <-> Candidate_resilient_ager
#          CUA <-> Concordant_unfavorable_ager
#
# Model:
#   case_binary ~ IDP_z + covariates
#
# Coding:
#   CRA = 1
#   comparator (default CUA) = 0
#
# Therefore:
#   beta_IDP > 0 / OR > 1
#       = higher IDP is associated with greater odds of being CRA
#
#   beta_IDP < 0 / OR < 1
#       = higher IDP is associated with greater odds of being CUA
#
# This is an association analysis; positive beta should be interpreted as
# CRA-associated, not as proof that the imaging feature causally promotes CRA.
#
# Generalizability:
#   The script is imaging-modality agnostic. The IDP can come from DTI,
#   T1 gray-matter, resting/task fMRI, or another brain-imaging table.
#
# Covariates can come from:
#   1) phenotype TSV           -- phenotype_covariates
#   2) general covariate file  -- covariates
#   3) imaging/IDP TSV         -- idp_covariates
#
# This lets later T1 analyses add intracranial-volume/head-size covariates,
# and fMRI analyses add motion/quality covariates without modifying this R code.
###############################################################################

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

# -----------------------------------------------------------------------------
# Arguments
# -----------------------------------------------------------------------------
#
#  1 phenotype_tsv
#  2 idp_tsv
#  3 cov_tsv
#  4 output_dir
#  5 idp
#  6 phenotype_col
#  7 case_label
#  8 control_label
#  9 covariates_csv
# 10 factor_covariates_csv
# 11 phenotype_covariates_csv
# 12 idp_covariates_csv
# 13 outlier_sd
# 14 min_total_n
# 15 min_per_group
#
# Example:
# Rscript fit_CRA_vs_CUA_logistic.R \
#   phenotype.tsv imaging.tsv covariates.csv output_dir \
#   mean_fa_... aging_phenotype_p20 CRA CUA \
#   "" "sex" "age_at_imaging,sex,Brain_ProtBAG_z" "" \
#   4 100 30
# -----------------------------------------------------------------------------

if (length(args) < 5) {
  stop(
    paste0(
      "Usage:\n",
      "  Rscript fit_CRA_vs_CUA_logistic.R ",
      "<phenotype_tsv> <idp_tsv> <cov_tsv> <output_dir> <idp> ",
      "[phenotype_col] [case_label] [control_label] ",
      "[covariates_csv] [factor_covariates_csv] ",
      "[phenotype_covariates_csv] [idp_covariates_csv] ",
      "[outlier_sd] [min_total_n] [min_per_group]\n"
    )
  )
}

phenotype_tsv <- args[[1]]
idp_tsv       <- args[[2]]
cov_tsv       <- args[[3]]
output_dir    <- args[[4]]
idp           <- args[[5]]

phenotype_col <- if (length(args) >= 6 && nzchar(args[[6]])) {
  args[[6]]
} else {
  "aging_phenotype_p20"
}

case_label <- if (length(args) >= 7 && nzchar(args[[7]])) {
  args[[7]]
} else {
  "CRA"
}

control_label <- if (length(args) >= 8 && nzchar(args[[8]])) {
  args[[8]]
} else {
  "CUA"
}

covariates_csv <- if (length(args) >= 9) args[[9]] else ""
factor_covariates_csv <- if (length(args) >= 10) args[[10]] else ""

phenotype_covariates_csv <- if (length(args) >= 11) {
  args[[11]]
} else {
  "age_at_imaging,sex,Brain_ProtBAG_z"
}

idp_covariates_csv <- if (length(args) >= 12) args[[12]] else ""

outlier_sd <- if (length(args) >= 13 && nzchar(args[[13]])) {
  as.numeric(args[[13]])
} else {
  4
}

min_total_n <- if (length(args) >= 14 && nzchar(args[[14]])) {
  as.integer(args[[14]])
} else {
  100L
}

min_per_group <- if (length(args) >= 15 && nzchar(args[[15]])) {
  as.integer(args[[15]])
} else {
  30L
}

if (!is.finite(outlier_sd) || outlier_sd <= 0) {
  stop("outlier_sd must be a positive finite number.")
}

if (is.na(min_total_n) || min_total_n < 2) {
  stop("min_total_n must be >= 2.")
}

if (is.na(min_per_group) || min_per_group < 1) {
  stop("min_per_group must be >= 1.")
}

# -----------------------------------------------------------------------------
# Libraries
# -----------------------------------------------------------------------------

cluster_lib <- "/cbica/home/wenju/R/x86_64-pc-linux-gnu-library/4.2.2"

if (dir.exists(cluster_lib)) {
  .libPaths(c(cluster_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

split_csv <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(trimws(x))) {
    return(character(0))
  }

  out <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  out <- out[nzchar(out)]
  unique(out)
}

quote_name <- function(x) {
  paste0("`", gsub("`", "", x, fixed = TRUE), "`")
}

safe_file_token <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x
}

standardize_participant_id <- function(dt, source_name) {
  if ("participant_id" %in% names(dt)) {
    # already standardized
  } else if ("eid" %in% names(dt)) {
    setnames(dt, "eid", "participant_id")
  } else {
    stop(
      source_name,
      " must contain either 'participant_id' or 'eid'."
    )
  }

  dt[, participant_id := as.character(participant_id)]
  dt
}

check_unique_ids <- function(dt, source_name) {
  dup <- dt[
    duplicated(participant_id) |
      duplicated(participant_id, fromLast = TRUE)
  ]

  if (nrow(dup) > 0) {
    stop(
      source_name,
      " contains duplicated participant_id values. ",
      "Resolve duplicates before modeling. Example duplicated ID: ",
      dup$participant_id[[1]]
    )
  }

  invisible(TRUE)
}

check_columns <- function(dt, cols, source_name) {
  cols <- unique(cols)
  cols <- cols[nzchar(cols)]

  missing_cols <- setdiff(cols, names(dt))

  if (length(missing_cols) > 0) {
    stop(
      source_name,
      " is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
}

normalize_label <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- tolower(x)
  x <- gsub("[[:space:]-]+", "_", x)
  x <- gsub("_+", "_", x)
  x
}

parse_bool <- function(x) {
  if (is.logical(x)) {
    return(x)
  }

  if (is.numeric(x) || is.integer(x)) {
    out <- rep(NA, length(x))
    out[!is.na(x) & x == 1] <- TRUE
    out[!is.na(x) & x == 0] <- FALSE
    return(as.logical(out))
  }

  z <- tolower(trimws(as.character(x)))

  out <- rep(NA, length(z))

  out[z %in% c("true", "t", "1", "yes", "y")] <- TRUE
  out[z %in% c("false", "f", "0", "no", "n")] <- FALSE

  as.logical(out)
}

extract_threshold_tag <- function(phenotype_col) {
  hit <- regmatches(
    phenotype_col,
    regexpr("p[0-9]+$", phenotype_col, perl = TRUE)
  )

  if (length(hit) == 0 || !nzchar(hit)) {
    return(NA_character_)
  }

  hit
}

canonical_aliases <- function(label) {

  z <- normalize_label(label)

  if (z == "cra") {
    return(
      c(
        "cra",
        "candidate_resilient_ager",
        "candidate_resilient_agers"
      )
    )
  }

  if (z == "cua") {
    return(
      c(
        "cua",
        "concordant_unfavorable_ager",
        "concordant_unfavourable_ager",
        "concordant_unfavorable_agers",
        "concordant_unfavourable_agers"
      )
    )
  }

  if (z == "cfa") {
    return(
      c(
        "cfa",
        "concordant_favorable_ager",
        "concordant_favourable_ager",
        "concordant_favorable_agers",
        "concordant_favourable_agers"
      )
    )
  }

  if (z == "lva") {
    return(
      c(
        "lva",
        "latent_vulnerability_ager",
        "latent_vulnerability_agers"
      )
    )
  }

  unique(
    c(
      z,
      normalize_label(label)
    )
  )
}

capture_glm <- function(formula, data) {
  warning_messages <- character(0)

  fit <- withCallingHandlers(
    glm(
      formula = formula,
      data = data,
      family = binomial(link = "logit")
    ),
    warning = function(w) {
      warning_messages <<- c(
        warning_messages,
        conditionMessage(w)
      )
      invokeRestart("muffleWarning")
    }
  )

  list(
    fit = fit,
    warnings = unique(warning_messages)
  )
}

# -----------------------------------------------------------------------------
# Parse covariate specifications
# -----------------------------------------------------------------------------

covariates <- split_csv(covariates_csv)
factor_covariates <- split_csv(factor_covariates_csv)
phenotype_covariates <- split_csv(phenotype_covariates_csv)
idp_covariates <- split_csv(idp_covariates_csv)

if (length(intersect(covariates, phenotype_covariates)) > 0) {
  stop(
    "The same covariate appears in covariates_csv and ",
    "phenotype_covariates_csv: ",
    paste(intersect(covariates, phenotype_covariates), collapse = ", ")
  )
}

if (length(intersect(covariates, idp_covariates)) > 0) {
  stop(
    "The same covariate appears in covariates_csv and idp_covariates_csv: ",
    paste(intersect(covariates, idp_covariates), collapse = ", ")
  )
}

if (length(intersect(phenotype_covariates, idp_covariates)) > 0) {
  stop(
    "The same covariate appears in phenotype_covariates_csv and ",
    "idp_covariates_csv: ",
    paste(intersect(phenotype_covariates, idp_covariates), collapse = ", ")
  )
}

all_model_covariates <- unique(
  c(
    phenotype_covariates,
    covariates,
    idp_covariates
  )
)

unknown_factor_covariates <- setdiff(
  factor_covariates,
  all_model_covariates
)

if (length(unknown_factor_covariates) > 0) {
  stop(
    "factor_covariates_csv contains covariate(s) not included in the model: ",
    paste(unknown_factor_covariates, collapse = ", ")
  )
}

# -----------------------------------------------------------------------------
# File checks
# -----------------------------------------------------------------------------

for (path in c(phenotype_tsv, idp_tsv)) {
  if (!file.exists(path)) {
    stop("Input file not found: ", path)
  }
}

if (length(covariates) > 0 && !file.exists(cov_tsv)) {
  stop(
    "General covariates were requested but cov_tsv does not exist: ",
    cov_tsv
  )
}

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# -----------------------------------------------------------------------------
# Read phenotype data
# -----------------------------------------------------------------------------

message("Reading phenotype file:")
message("  ", phenotype_tsv)

phenotype <- fread(
  phenotype_tsv,
  data.table = TRUE,
  showProgress = FALSE
)

phenotype <- standardize_participant_id(
  phenotype,
  "phenotype_tsv"
)

check_unique_ids(
  phenotype,
  "phenotype_tsv"
)

check_columns(
  phenotype,
  c(
    "participant_id",
    phenotype_col,
    phenotype_covariates
  ),
  "phenotype_tsv"
)

# -----------------------------------------------------------------------------
# Construct CRA-vs-CUA (or other requested) binary phenotype
# -----------------------------------------------------------------------------

threshold_tag <- extract_threshold_tag(
  phenotype_col
)

case_indicator_col <- NA_character_
control_indicator_col <- NA_character_

if (!is.na(threshold_tag)) {

  candidate_case_indicator <- paste0(
    case_label,
    "_",
    threshold_tag
  )

  candidate_control_indicator <- paste0(
    control_label,
    "_",
    threshold_tag
  )

  if (
    candidate_case_indicator %in% names(phenotype) &&
      candidate_control_indicator %in% names(phenotype)
  ) {
    case_indicator_col <- candidate_case_indicator
    control_indicator_col <- candidate_control_indicator
  }
}

phenotype_creation_method <- NA_character_

if (
  !is.na(case_indicator_col) &&
    !is.na(control_indicator_col)
) {

  # ---------------------------------------------------------------------------
  # Preferred method: use explicit CRA_p20 / CUA_p20 indicator columns.
  # ---------------------------------------------------------------------------

  message("")
  message("Creating binary phenotype from indicator columns:")
  message("  Case indicator:    ", case_indicator_col)
  message("  Control indicator: ", control_indicator_col)

  phenotype[
    ,
    case_flag_internal := parse_bool(
      get(case_indicator_col)
    )
  ]

  phenotype[
    ,
    control_flag_internal := parse_bool(
      get(control_indicator_col)
    )
  ]

  # Participants must belong to exactly one of the two groups.
  phenotype[
    ,
    case_binary := fcase(
      case_flag_internal %in% TRUE &
        !(control_flag_internal %in% TRUE),
      1L,

      control_flag_internal %in% TRUE &
        !(case_flag_internal %in% TRUE),
      0L,

      default = NA_integer_
    )
  ]

  phenotype_creation_method <- paste0(
    "indicator_columns:",
    case_indicator_col,
    "_vs_",
    control_indicator_col
  )

} else {

  # ---------------------------------------------------------------------------
  # Fallback method: map full phenotype names in aging_phenotype_pXX.
  # ---------------------------------------------------------------------------

  message("")
  message(
    "Indicator columns were not found for the requested comparison; ",
    "falling back to phenotype-label mapping."
  )

  case_aliases <- canonical_aliases(
    case_label
  )

  control_aliases <- canonical_aliases(
    control_label
  )

  phenotype[
    ,
    phenotype_label_normalized := normalize_label(
      get(phenotype_col)
    )
  ]

  phenotype[
    ,
    case_binary := fcase(
      phenotype_label_normalized %in% case_aliases,
      1L,

      phenotype_label_normalized %in% control_aliases,
      0L,

      default = NA_integer_
    )
  ]

  phenotype_creation_method <- paste0(
    "label_mapping_from:",
    phenotype_col
  )
}

# Keep only requested case/control groups.
phenotype <- phenotype[
  !is.na(case_binary)
]

if (nrow(phenotype) == 0) {

  observed_labels <- sort(
    unique(
      as.character(
        fread(
          phenotype_tsv,
          select = phenotype_col,
          showProgress = FALSE
        )[[phenotype_col]]
      )
    )
  )

  stop(
    paste0(
      "No participants could be assigned to the requested comparison.\n\n",
      "Requested:\n",
      "  case = ", case_label, "\n",
      "  control = ", control_label, "\n",
      "  phenotype column = ", phenotype_col, "\n\n",
      "Observed values in ", phenotype_col, " include:\n  ",
      paste(
        head(observed_labels, 20),
        collapse = "\n  "
      ),
      "\n\n",
      "Expected for your resilience TSV:\n",
      "  CRA = Candidate_resilient_ager or CRA_", threshold_tag, "\n",
      "  CUA = Concordant_unfavorable_ager or CUA_", threshold_tag, "\n"
    )
  )
}

n_case_phenotype <- phenotype[
  case_binary == 1L,
  .N
]

n_control_phenotype <- phenotype[
  case_binary == 0L,
  .N
]

if (
  n_case_phenotype == 0 ||
    n_control_phenotype == 0
) {
  stop(
    "Binary phenotype construction produced an empty group: ",
    case_label,
    "=",
    n_case_phenotype,
    "; ",
    control_label,
    "=",
    n_control_phenotype
  )
}

message("")
message("Phenotype comparison created successfully:")
message("  Case    = ", case_label, " -> 1")
message("  Control = ", control_label, " -> 0")
message("  Method  = ", phenotype_creation_method)
message("  N case in phenotype file: ", n_case_phenotype)
message("  N control in phenotype file: ", n_control_phenotype)

# Retain only required phenotype columns after outcome creation.
phenotype_keep <- unique(
  c(
    "participant_id",
    "case_binary",
    phenotype_col,
    phenotype_covariates
  )
)

phenotype_keep <- phenotype_keep[
  phenotype_keep %in% names(phenotype)
]

phenotype <- phenotype[
  ,
  phenotype_keep,
  with = FALSE
]

# -----------------------------------------------------------------------------
# Read imaging IDP data
# -----------------------------------------------------------------------------

message("")
message("Reading imaging file:")
message("  ", idp_tsv)
message("IDP:")
message("  ", idp)

idp_data <- fread(
  idp_tsv,
  data.table = TRUE,
  showProgress = FALSE
)

idp_data <- standardize_participant_id(
  idp_data,
  "idp_tsv"
)

check_unique_ids(
  idp_data,
  "idp_tsv"
)

check_columns(
  idp_data,
  c(
    "participant_id",
    idp,
    idp_covariates
  ),
  "idp_tsv"
)

idp_data <- idp_data[
  ,
  c(
    "participant_id",
    idp,
    idp_covariates
  ),
  with = FALSE
]

setnames(
  idp_data,
  idp,
  "IDP_raw"
)

# -----------------------------------------------------------------------------
# Read optional general covariates
# -----------------------------------------------------------------------------

if (length(covariates) > 0) {

  message("")
  message("Reading general covariate file:")
  message("  ", cov_tsv)

  cov_data <- fread(
    cov_tsv,
    data.table = TRUE,
    showProgress = FALSE
  )

  cov_data <- standardize_participant_id(
    cov_data,
    "cov_tsv"
  )

  check_unique_ids(
    cov_data,
    "cov_tsv"
  )

  check_columns(
    cov_data,
    c(
      "participant_id",
      covariates
    ),
    "cov_tsv"
  )

  cov_data <- cov_data[
    ,
    c(
      "participant_id",
      covariates
    ),
    with = FALSE
  ]

} else {

  cov_data <- NULL
}

# -----------------------------------------------------------------------------
# Merge
# -----------------------------------------------------------------------------

df <- merge(
  phenotype,
  idp_data,
  by = "participant_id",
  all = FALSE,
  sort = FALSE
)

if (!is.null(cov_data)) {
  df <- merge(
    df,
    cov_data,
    by = "participant_id",
    all = FALSE,
    sort = FALSE
  )
}

n_after_merge <- nrow(df)

if (n_after_merge == 0) {
  stop(
    "No overlapping participants after merging phenotype and imaging data."
  )
}

# IDP must be numeric.
df[
  ,
  IDP_raw := suppressWarnings(
    as.numeric(IDP_raw)
  )
]

df[
  !is.finite(IDP_raw),
  IDP_raw := NA_real_
]

# Convert requested categorical covariates to factors.
for (v in factor_covariates) {
  df[[v]] <- factor(df[[v]])
}

model_variables <- unique(
  c(
    "case_binary",
    "IDP_raw",
    all_model_covariates
  )
)

check_columns(
  df,
  model_variables,
  "merged analysis table"
)

# Complete-case model sample.
complete_idx <- complete.cases(
  df[
    ,
    model_variables,
    with = FALSE
  ]
)

df_model <- df[
  complete_idx
]

n_complete_before_outlier <- nrow(
  df_model
)

if (n_complete_before_outlier < min_total_n) {
  stop(
    "Too few complete participants before IDP outlier filtering: ",
    n_complete_before_outlier,
    " < min_total_n=",
    min_total_n
  )
}

# -----------------------------------------------------------------------------
# IDP outlier filter
# -----------------------------------------------------------------------------

idp_mean_pre <- mean(
  df_model$IDP_raw,
  na.rm = TRUE
)

idp_sd_pre <- sd(
  df_model$IDP_raw,
  na.rm = TRUE
)

if (!is.finite(idp_sd_pre) || idp_sd_pre <= 0) {
  stop(
    "IDP has zero or invalid SD in the complete-case sample."
  )
}

lower_limit <- idp_mean_pre - outlier_sd * idp_sd_pre
upper_limit <- idp_mean_pre + outlier_sd * idp_sd_pre

n_before_outlier <- nrow(
  df_model
)

df_model <- df_model[
  IDP_raw >= lower_limit &
    IDP_raw <= upper_limit
]

n_after_outlier <- nrow(
  df_model
)

n_outliers_removed <- (
  n_before_outlier -
    n_after_outlier
)

if (n_after_outlier < min_total_n) {
  stop(
    "Too few participants after IDP outlier filtering: ",
    n_after_outlier,
    " < min_total_n=",
    min_total_n
  )
}

# Standardize IDP after outlier filtering.
idp_mean_final <- mean(
  df_model$IDP_raw,
  na.rm = TRUE
)

idp_sd_final <- sd(
  df_model$IDP_raw,
  na.rm = TRUE
)

if (!is.finite(idp_sd_final) || idp_sd_final <= 0) {
  stop(
    "IDP has zero or invalid SD after outlier filtering."
  )
}

df_model[
  ,
  IDP_z := (
    IDP_raw -
      idp_mean_final
  ) / idp_sd_final
]

n_case <- df_model[
  case_binary == 1L,
  .N
]

n_control <- df_model[
  case_binary == 0L,
  .N
]

if (
  n_case < min_per_group ||
    n_control < min_per_group
) {
  stop(
    "Too few participants in one comparison group after QC. ",
    case_label,
    "=",
    n_case,
    "; ",
    control_label,
    "=",
    n_control,
    "; min_per_group=",
    min_per_group
  )
}

# -----------------------------------------------------------------------------
# Logistic regression
# -----------------------------------------------------------------------------

rhs_terms <- c(
  "IDP_z",
  vapply(
    all_model_covariates,
    quote_name,
    FUN.VALUE = character(1)
  )
)

formula_text <- paste0(
  "case_binary ~ ",
  paste(
    rhs_terms,
    collapse = " + "
  )
)

model_formula <- as.formula(
  formula_text
)

message("")
message("Model formula:")
message("  ", formula_text)

message("")
message("Final analysis N:")
message("  ", case_label, ": ", n_case)
message("  ", control_label, ": ", n_control)
message("  Total: ", nrow(df_model))
message("  IDP outliers removed: ", n_outliers_removed)

fit_obj <- capture_glm(
  formula = model_formula,
  data = df_model
)

model <- fit_obj$fit
glm_warnings <- fit_obj$warnings

coef_table <- summary(model)$coefficients

if (!"IDP_z" %in% rownames(coef_table)) {
  stop(
    "IDP_z coefficient is absent from the fitted logistic model."
  )
}

idp_row <- coef_table[
  "IDP_z",
  ,
  drop = FALSE
]

beta <- unname(
  idp_row[
    1,
    "Estimate"
  ]
)

se <- unname(
  idp_row[
    1,
    "Std. Error"
  ]
)

z_value <- unname(
  idp_row[
    1,
    "z value"
  ]
)

p_value <- unname(
  idp_row[
    1,
    "Pr(>|z|)"
  ]
)

or <- exp(beta)
ci_lower <- exp(beta - 1.96 * se)
ci_upper <- exp(beta + 1.96 * se)

direction <- if (
  is.na(beta)
) {
  NA_character_
} else if (
  beta > 0
) {
  paste0(
    "Higher_IDP_associated_with_",
    case_label
  )
} else if (
  beta < 0
) {
  paste0(
    "Higher_IDP_associated_with_",
    control_label
  )
} else {
  "No_direction"
}

# Flag potentially extreme fitted probabilities / separation.
fitted_prob <- fitted(model)

possible_separation <- any(
  fitted_prob < 1e-6 |
    fitted_prob > 1 - 1e-6
)

warning_text <- if (length(glm_warnings) > 0) {
  paste(
    glm_warnings,
    collapse = " | "
  )
} else {
  ""
}

# -----------------------------------------------------------------------------
# Output
# -----------------------------------------------------------------------------

result <- data.table(
  comparison = paste0(
    case_label,
    "_vs_",
    control_label
  ),

  phenotype_column = phenotype_col,
  phenotype_creation_method = phenotype_creation_method,

  case_indicator_column = ifelse(
    is.na(case_indicator_col),
    "",
    case_indicator_col
  ),

  control_indicator_column = ifelse(
    is.na(control_indicator_col),
    "",
    control_indicator_col
  ),

  case_label = case_label,
  control_label = control_label,
  case_coding = 1L,
  control_coding = 0L,

  IDP = idp,

  Beta_log_odds_per_1SD_IDP = beta,
  SE = se,
  Z = z_value,
  P_value = p_value,

  OR_per_1SD_IDP = or,
  OR_CI95_lower = ci_lower,
  OR_CI95_upper = ci_upper,

  Direction = direction,

  N_case = n_case,
  N_control = n_control,
  N_total = nrow(df_model),

  N_case_in_phenotype_file = n_case_phenotype,
  N_control_in_phenotype_file = n_control_phenotype,
  N_after_merge = n_after_merge,
  N_complete_before_outlier_filter = n_complete_before_outlier,
  N_IDP_outliers_removed = n_outliers_removed,

  IDP_outlier_SD_threshold = outlier_sd,
  IDP_mean_before_outlier_filter = idp_mean_pre,
  IDP_SD_before_outlier_filter = idp_sd_pre,
  IDP_mean_analysis_sample = idp_mean_final,
  IDP_SD_analysis_sample = idp_sd_final,

  model_converged = isTRUE(
    model$converged
  ),

  possible_separation = possible_separation,
  AIC = AIC(model),

  phenotype_covariates = paste(
    phenotype_covariates,
    collapse = ","
  ),

  general_covariates = paste(
    covariates,
    collapse = ","
  ),

  idp_covariates = paste(
    idp_covariates,
    collapse = ","
  ),

  factor_covariates = paste(
    factor_covariates,
    collapse = ","
  ),

  model_formula = formula_text,
  glm_warnings = warning_text
)

comparison_token <- safe_file_token(
  paste0(
    case_label,
    "_vs_",
    control_label
  )
)

idp_token <- safe_file_token(
  idp
)

output_file <- file.path(
  output_dir,
  paste0(
    comparison_token,
    "_logistic_results_",
    idp_token,
    ".tsv"
  )
)

fwrite(
  result,
  output_file,
  sep = "\t",
  na = "NA",
  quote = FALSE
)

message("")
message("Result:")
message(
  "  beta = ",
  signif(beta, 5),
  "; OR = ",
  signif(or, 5),
  "; P = ",
  format.pval(
    p_value,
    digits = 5,
    eps = 1e-300
  )
)

message(
  "  Direction: ",
  direction
)

message("")
message("Wrote:")
message("  ", output_file)