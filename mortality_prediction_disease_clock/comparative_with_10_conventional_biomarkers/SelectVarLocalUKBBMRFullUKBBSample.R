#!/usr/bin/env Rscript

# ==============================================================================
# Extract 10 conventional mortality biomarkers from the full UK Biobank file
# ==============================================================================
#
# This script follows the prefix-based variable-selection logic used in the
# existing UK Biobank extraction workflow, but it is optimized for the very
# large full-sample CSV. Rather than loading all columns into memory, it:
#
#   1) reads only the CSV header;
#   2) resolves all exact UKB columns beginning with each requested Var prefix;
#   3) loads only eid + those resolved columns using data.table::fread(select=);
#   4) writes the reduced full-sample CSV;
#   5) writes a manifest of resolved columns and a biomarker-level coverage TSV.
#
# The variable list contains 11 rows because grip strength is one conceptual
# biomarker represented by separate left- and right-hand UKB fields (46, 47).
#
# Run:
#   Rscript extract_10_conventional_mortality_biomarkers.R
#
# Optional positional arguments:
#   1. variable-list TSV
#   2. full UKB CSV
#   3. output directory
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

base_dir <- paste0(
  "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/",
  "comparative_with_10_conventional_biomarkers/1_prepare_data"
)

var_list_file <- if (length(args) >= 1) {
  args[1]
} else {
  file.path(
    base_dir,
    "phenotype_ukbb_variable_list_10_mortality_biomarkers.tsv"
  )
}

ukbb_full_file <- if (length(args) >= 2) {
  args[2]
} else {
  paste0(
    "/cbica/projects/ISTAGING/Pipelines/ClinicalDataConsolidation_201911/",
    "Data/External_Data/UKBiobank/ukb_230424_FullUKBSample.csv"
  )
}

out_dir <- if (length(args) >= 3) {
  args[3]
} else {
  base_dir
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_csv <- file.path(
  out_dir,
  "UKBB_fullsample_10_conventional_mortality_biomarkers.csv"
)

out_manifest <- file.path(
  out_dir,
  "UKBB_fullsample_10_conventional_mortality_biomarkers_resolved_columns.tsv"
)

out_coverage <- file.path(
  out_dir,
  "UKBB_fullsample_10_conventional_mortality_biomarkers_coverage.tsv"
)

# ------------------------------------------------------------------------------
# Checks
# ------------------------------------------------------------------------------

if (!file.exists(var_list_file)) {
  stop("Variable-list TSV not found: ", var_list_file)
}

if (!file.exists(ukbb_full_file)) {
  stop("Full UK Biobank CSV not found: ", ukbb_full_file)
}

cat("======================================================================\n")
cat("Extracting 10 conventional mortality biomarkers from UK Biobank\n")
cat("======================================================================\n")
cat("Variable list: ", var_list_file, "\n", sep = "")
cat("Full UKB CSV:  ", ukbb_full_file, "\n", sep = "")
cat("Output CSV:    ", out_csv, "\n", sep = "")
cat("======================================================================\n\n")

# ------------------------------------------------------------------------------
# Read requested variable prefixes
# ------------------------------------------------------------------------------

sel <- fread(
  var_list_file,
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  check.names = FALSE
)

required_varlist_cols <- c(
  "Type",
  "Category",
  "Sub-category",
  "Var",
  "Field ID"
)

missing_varlist_cols <- setdiff(required_varlist_cols, names(sel))

if (length(missing_varlist_cols) > 0) {
  stop(
    "Variable-list TSV is missing required columns: ",
    paste(missing_varlist_cols, collapse = ", ")
  )
}

sel$Var <- trimws(as.character(sel$Var))
sel <- sel[!is.na(sel$Var) & nzchar(sel$Var), , drop = FALSE]

if (nrow(sel) == 0) {
  stop("Variable-list TSV contains no usable Var prefixes.")
}

# ------------------------------------------------------------------------------
# Read only the header of the very large UKB file
# ------------------------------------------------------------------------------

cat("Reading UK Biobank header only...\n")

ukbb_header <- fread(
  ukbb_full_file,
  nrows = 0,
  data.table = FALSE,
  check.names = FALSE,
  showProgress = FALSE
)

ukbb_names <- names(ukbb_header)

if (!("eid" %in% ukbb_names)) {
  stop("The full UK Biobank file does not contain 'eid'.")
}

cat("Total columns in full UKB file: ", length(ukbb_names), "\n\n", sep = "")

# ------------------------------------------------------------------------------
# Resolve each prefix to all matching UKB columns
# ------------------------------------------------------------------------------

manifest_list <- vector("list", nrow(sel))
missing_prefixes <- character(0)

for (i in seq_len(nrow(sel))) {
  prefix <- sel$Var[i]
  matched <- ukbb_names[startsWith(ukbb_names, prefix)]

  if (length(matched) == 0) {
    missing_prefixes <- c(missing_prefixes, prefix)

    manifest_list[[i]] <- data.frame(
      Type = sel$Type[i],
      Category = sel$Category[i],
      Sub.category = sel[["Sub-category"]][i],
      Var = prefix,
      Field.ID = sel[["Field ID"]][i],
      Resolved_column = NA_character_,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else {
    manifest_list[[i]] <- data.frame(
      Type = rep(sel$Type[i], length(matched)),
      Category = rep(sel$Category[i], length(matched)),
      Sub.category = rep(sel[["Sub-category"]][i], length(matched)),
      Var = rep(prefix, length(matched)),
      Field.ID = rep(sel[["Field ID"]][i], length(matched)),
      Resolved_column = matched,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
}

manifest <- do.call(rbind, manifest_list)

# Restore original requested display column names in the manifest.
names(manifest)[names(manifest) == "Sub.category"] <- "Sub-category"
names(manifest)[names(manifest) == "Field.ID"] <- "Field ID"

fwrite(
  manifest,
  out_manifest,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

if (length(missing_prefixes) > 0) {
  cat("ERROR: the following requested prefixes were not found:\n")
  for (x in unique(missing_prefixes)) {
    cat("  - ", x, "\n", sep = "")
  }
  stop(
    "At least one requested biomarker prefix was not found. ",
    "See the resolved-column manifest: ",
    out_manifest
  )
}

resolved_cols <- unique(manifest$Resolved_column)
resolved_cols <- resolved_cols[!is.na(resolved_cols)]
selected_cols <- unique(c("eid", resolved_cols))

cat("Resolved requested prefixes: ", nrow(sel), "\n", sep = "")
cat("Resolved UKB columns:        ", length(resolved_cols), "\n", sep = "")
cat("Columns loaded incl. eid:    ", length(selected_cols), "\n\n", sep = "")

for (i in seq_len(nrow(sel))) {
  prefix <- sel$Var[i]
  matched <- resolved_cols[startsWith(resolved_cols, prefix)]
  cat(
    sprintf(
      "%-60s -> %d column(s)\n",
      prefix,
      length(matched)
    )
  )
}

# ------------------------------------------------------------------------------
# Read ONLY selected columns from the full 500k sample
# ------------------------------------------------------------------------------

cat("\nReading selected columns from full UK Biobank sample...\n")

ukbb_selected <- fread(
  ukbb_full_file,
  select = selected_cols,
  data.table = FALSE,
  check.names = FALSE,
  showProgress = TRUE
)

cat("Participants extracted: ", nrow(ukbb_selected), "\n", sep = "")
cat("Variables extracted:    ", ncol(ukbb_selected), "\n\n", sep = "")

# ------------------------------------------------------------------------------
# Biomarker-level coverage
# ------------------------------------------------------------------------------
# Coverage is defined as having at least one non-missing value among all exact
# UKB columns matching a requested Var prefix. For example, pulse-rate coverage
# uses all measurements beginning with pulse_rate_automated_reading_f102.
#
# Left and right grip strength are reported separately in this QC table because
# they are separate UKB fields, although they form one conceptual comparator.
# ------------------------------------------------------------------------------

coverage_list <- vector("list", nrow(sel))

for (i in seq_len(nrow(sel))) {
  prefix <- sel$Var[i]
  matched <- names(ukbb_selected)[startsWith(names(ukbb_selected), prefix)]

  if (length(matched) == 1) {
    has_value <- !is.na(ukbb_selected[[matched]])
  } else {
    has_value <- rowSums(!is.na(ukbb_selected[, matched, drop = FALSE])) > 0
  }

  n_available <- sum(has_value)
  n_total <- nrow(ukbb_selected)

  coverage_list[[i]] <- data.frame(
    Type = sel$Type[i],
    Category = sel$Category[i],
    `Sub-category` = sel[["Sub-category"]][i],
    Var = prefix,
    `Field ID` = sel[["Field ID"]][i],
    N_resolved_columns = length(matched),
    N_nonmissing_any_instance = n_available,
    N_total = n_total,
    Coverage_fraction = n_available / n_total,
    Coverage_percent = 100 * n_available / n_total,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

coverage <- do.call(rbind, coverage_list)

fwrite(
  coverage,
  out_coverage,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ------------------------------------------------------------------------------
# Write reduced full-sample file
# ------------------------------------------------------------------------------

cat("Writing reduced full-sample CSV...\n")

fwrite(
  ukbb_selected,
  out_csv,
  sep = ",",
  quote = TRUE,
  na = "NA"
)

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------

cat("\n======================================================================\n")
cat("Finished\n")
cat("======================================================================\n")
cat("Reduced full-sample CSV:\n  ", out_csv, "\n", sep = "")
cat("Resolved-column manifest:\n  ", out_manifest, "\n", sep = "")
cat("Coverage summary:\n  ", out_coverage, "\n", sep = "")
cat("\nCoverage by requested UKB field:\n")

for (i in seq_len(nrow(coverage))) {
  cat(
    sprintf(
      "  %-35s  %8d / %8d  (%6.2f%%)\n",
      coverage[["Sub-category"]][i],
      coverage$N_nonmissing_any_instance[i],
      coverage$N_total[i],
      coverage$Coverage_percent[i]
    )
  )
}
cat("======================================================================\n")