#!/usr/bin/env Rscript

###############################################################################
# CRA vs comparator brain-wide association using logistic regression
#
# Primary intended analysis:
#   Candidate resilient agers (CRA) versus concordant unfavorable agers (CUA)
#   using the brain-proteomics BAG-EPOCH discordance phenotypes.
#
# Model:
#   CRA_case ~ IDP_z + covariates
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
#       = higher IDP is associated with greater odds of being the comparator
#
# IMPORTANT:
# This is an association model, not a causal model. "Positive beta" should be
# interpreted as CRA-associated rather than as proof that the imaging feature
# causally promotes resilience.
#
# Why CUA is the primary comparator:
#   CRA: high BAG + low EPOCH|BAG residual
#   CUA: high BAG + high EPOCH|BAG residual
#
# Both groups are in the high-BAG tail, so CRA vs CUA contrasts lower versus
# higher event-proximal vulnerability among people with similarly advanced
# generalized biological aging.
#
# Generalizability:
#   The script is imaging-modality agnostic. The IDP can come from DTI,
#   T1 gray-matter, resting/task fMRI, or another brain-imaging table.
#
#   Covariates can come from three places:
#     1) phenotype TSV      -- phenotype_covariates
#     2) general covariate file -- covariates
#     3) imaging/IDP TSV    -- idp_covariates
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
    # already correct
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
  dup <- dt[duplicated(participant_id) | duplicated(participant_id, fromLast = TRUE)]

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

# Prevent accidental duplication across sources.
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

phenotype <- phenotype[
  ,
  c(
    "participant_id",
    phenotype_col,
    phenotype_covariates
  ),
  with = FALSE
]

# Keep only the two groups being compared.
phenotype[, phenotype_label := as.character(get(phenotype_col))]

phenotype <- phenotype[
  phenotype_label %in% c(case_label, control_label)
]

if (nrow(phenotype) == 0) {
  stop(
    "No participants found with ",
    phenotype_col,
    " equal to ",
    case_label,
    " or ",
    control_label,
    "."
  )
}

phenotype[
  ,
  CRA_case := fifelse(
    phenotype_label == case_label,
    1L,
    0L
  )
]

n_case_phenotype <- phenotype[CRA_case == 1L, .N]
n_control_phenotype <- phenotype[CRA_case == 0L, .N]

message("")
message("Phenotype comparison:")
message("  Case    = ", case_label, " -> 1")
message("  Control = ", control_label, " -> 0")
message("  N case in phenotype file: ", n_case_phenotype)
message("  N control in phenotype file: ", n_control_phenotype)

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

# IDP must be numeric.
df[, IDP_raw := suppressWarnings(as.numeric(IDP_raw))]
df[!is.finite(IDP_raw), IDP_raw := NA_real_]

# Convert requested categorical covariates to factors.
for (v in factor_covariates) {
  df[[v]] <- factor(df[[v]])
}

model_variables <- unique(
  c(
    "CRA_case",
    "IDP_raw",
    all_model_covariates
  )
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

n_complete_before_outlier <- nrow(df_model)

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
  stop("IDP has zero or invalid SD in the complete-case sample.")
}

lower_limit <- idp_mean_pre - outlier_sd * idp_sd_pre
upper_limit <- idp_mean_pre + outlier_sd * idp_sd_pre

n_before_outlier <- nrow(df_model)

df_model <- df_model[
  IDP_raw >= lower_limit &
    IDP_raw <= upper_limit
]

n_after_outlier <- nrow(df_model)
n_outliers_removed <- n_before_outlier - n_after_outlier

if (n_after_outlier < min_total_n) {
  stop(
    "Too few participants after IDP outlier filtering: ",
    n_after_outlier,
    " < min_total_n=",
    min_total_n
  )
}

# Standardize after outlier filtering.
idp_mean_final <- mean(
  df_model$IDP_raw,
  na.rm = TRUE
)

idp_sd_final <- sd(
  df_model$IDP_raw,
  na.rm = TRUE
)

if (!is.finite(idp_sd_final) || idp_sd_final <= 0) {
  stop("IDP has zero or invalid SD after outlier filtering.")
}

df_model[
  ,
  IDP_z := (
    IDP_raw - idp_mean_final
  ) / idp_sd_final
]

n_case <- df_model[CRA_case == 1L, .N]
n_control <- df_model[CRA_case == 0L, .N]

if (n_case < min_per_group || n_control < min_per_group) {
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
  "CRA_case ~ ",
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

  model_converged = isTRUE(model$converged),
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
message("  Direction: ", direction)
message("")
message("Wrote:")
message("  ", output_file)