#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 6) {
  stop("Usage: MR_2_harmonization.R <exposure_gwas_tsv> <exposure_fuma_tsv> <exposure_key> <outcome_2sampleMR_tsv> <outcome_key> <output_dir_mr>")
}

exposure_gwas_tsv <- args[1]
exposure_fuma_tsv <- args[2]
exposure_key <- args[3]
output_2sampleMR_tsv <- args[4]
outcome_key <- args[5]
output_dir_mr <- args[6]

.libPaths('/gpfs/fs001/cbica/home/wenju/R/x86_64-pc-linux-gnu-library/4.3')

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(MRInstruments)
  library(stringr)
  library(dplyr)
  library(ieugwasr)
})

message("Exposure key: ", exposure_key)
message("Outcome key: ", outcome_key)
message("Exposure GWAS: ", exposure_gwas_tsv)
message("Exposure FUMA: ", exposure_fuma_tsv)
message("Outcome file: ", output_2sampleMR_tsv)
message("Output directory: ", output_dir_mr)

dir.create(output_dir_mr, recursive = TRUE, showWarnings = FALSE)

done_file <- file.path(output_dir_mr, paste0("DONE_", exposure_key, "_2_", outcome_key, ".tsv"))
harmonized_file <- file.path(output_dir_mr, paste0("harmonized_data_", exposure_key, "_2_", outcome_key, ".tsv"))
log_file <- file.path(output_dir_mr, paste0("harmonization_log_", exposure_key, "_2_", outcome_key, ".tsv"))

write_status <- function(status, n_fuma = NA_integer_, n_exp = NA_integer_, n_clumped = NA_integer_, n_outcome = NA_integer_, n_harmonized = NA_integer_) {
  status_df <- data.frame(
    exposure = exposure_key,
    outcome = outcome_key,
    status = status,
    n_fuma = n_fuma,
    n_exposure_read = n_exp,
    n_clumped = n_clumped,
    n_outcome_read = n_outcome,
    n_harmonized = n_harmonized,
    stringsAsFactors = FALSE
  )
  write.table(status_df, file = log_file, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  write.table(status, file = done_file, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
}

if (!file.exists(exposure_gwas_tsv)) {
  write_status("ERROR_MISSING_EXPOSURE_GWAS")
  stop("Missing exposure GWAS file: ", exposure_gwas_tsv)
}
if (!file.exists(exposure_fuma_tsv)) {
  write_status("ERROR_MISSING_FUMA_FILE")
  stop("Missing FUMA file: ", exposure_fuma_tsv)
}
if (!file.exists(output_2sampleMR_tsv)) {
  write_status("ERROR_MISSING_OUTCOME_FILE")
  stop("Missing outcome 2SampleMR file: ", output_2sampleMR_tsv)
}

# =============================
# Read FUMA independent significant SNPs
# =============================
fuma_df <- read.table(exposure_fuma_tsv, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, fill = TRUE, comment.char = "")

snp_col_candidates <- c("rsID", "SNP", "rsid", "ID")
snp_col <- snp_col_candidates[snp_col_candidates %in% names(fuma_df)][1]

if (is.na(snp_col)) {
  write_status("ERROR_NO_SNP_COLUMN_IN_FUMA", n_fuma = nrow(fuma_df))
  stop("No SNP column found in FUMA file. Expected one of: ", paste(snp_col_candidates, collapse = ", "))
}

snp_exp <- unique(na.omit(as.character(fuma_df[[snp_col]])))
snp_exp <- snp_exp[snp_exp != ""]

message("Number of FUMA SNPs before reading exposure GWAS: ", length(snp_exp))

if (length(snp_exp) <= 7) {
  write_status("DONE_TOO_FEW_FUMA_SNPS", n_fuma = length(snp_exp))
  quit(save = "no", status = 0)
}

# =============================
# Read exposure GWAS for FUMA SNPs
# =============================
exp_data <- read_outcome_data(
  snps = snp_exp,
  filename = exposure_gwas_tsv,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "AF1",
  pval_col = "P",
  samplesize_col = "N",
  pos_col = "POS"
)

if (nrow(exp_data) == 0) {
  write_status("DONE_NO_EXPOSURE_SNPS_READ", n_fuma = length(snp_exp), n_exp = 0)
  quit(save = "no", status = 0)
}

exp_data <- exp_data[, c(
  "pos.outcome", "SNP", "other_allele.outcome", "effect_allele.outcome",
  "samplesize.outcome", "beta.outcome", "se.outcome", "pval.outcome",
  "eaf.outcome", "outcome"
)]

exp_data$outcome <- exposure_key

exp_data <- format_data(
  exp_data,
  phenotype_col = "outcome",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  eaf_col = "eaf.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  pval_col = "pval.outcome",
  samplesize_col = "samplesize.outcome",
  pos_col = "pos.outcome",
  type = "exposure"
)

exp_data <- subset(exp_data, mr_keep.exposure == TRUE)

if (nrow(exp_data) <= 7) {
  write_status("DONE_TOO_FEW_EXPOSURE_SNPS_AFTER_FORMAT", n_fuma = length(snp_exp), n_exp = nrow(exp_data))
  quit(save = "no", status = 0)
}

# =============================
# Local LD clumping using 1000G EUR
# =============================
exp_for_clump <- exp_data
names(exp_for_clump)[names(exp_for_clump) == "SNP"] <- "rsid"
names(exp_for_clump)[names(exp_for_clump) == "pval.exposure"] <- "pval"
exp_for_clump$id <- exposure_key

exp_for_clump <- ld_clump(
  exp_for_clump,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p = 1,
  plink_bin = "/cbica/home/wenju/Software/plink_1.90_beta_linux_x86_64_20210606/plink",
  bfile = "/cbica/home/wenju/Dataset/2SampleMR_1000Genome_LD/EUR"
)

names(exp_for_clump)[names(exp_for_clump) == "rsid"] <- "SNP"
names(exp_for_clump)[names(exp_for_clump) == "pval"] <- "pval.exposure"

exp_data <- exp_for_clump

message("Number of clumped exposure IVs: ", nrow(exp_data))

if (nrow(exp_data) <= 7) {
  write_status("DONE_TOO_FEW_EXPOSURE_SNPS_AFTER_CLUMP", n_fuma = length(snp_exp), n_exp = length(unique(exp_data$SNP)), n_clumped = nrow(exp_data))
  quit(save = "no", status = 0)
}

# =============================
# Read outcome data
# =============================
outcome_data <- read_outcome_data(
  snps = exp_data$SNP,
  filename = output_2sampleMR_tsv,
  sep = "\t",
  snp_col = "ID",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "eaf_A1",
  pval_col = "P",
  samplesize_col = "N",
  phenotype_col = "pheno",
  pos_col = "POS"
)

outcome_data <- subset(outcome_data, mr_keep.outcome == TRUE)

message("Number of outcome SNPs read: ", nrow(outcome_data))

if (nrow(outcome_data) == 0) {
  write_status("DONE_NO_OUTCOME_SNPS", n_fuma = length(snp_exp), n_exp = nrow(exp_data), n_clumped = nrow(exp_data), n_outcome = 0)
  quit(save = "no", status = 0)
}

# =============================
# Harmonize exposure and outcome
# =============================
harmonized_data <- harmonise_data(
  exposure_dat = exp_data,
  outcome_dat = outcome_data
)

n_harmonized <- nrow(harmonized_data)
message("Number of harmonized IVs: ", n_harmonized)

if (n_harmonized > 7) {
  write.table(harmonized_data, file = harmonized_file, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  write_status("DONE_OK", n_fuma = length(snp_exp), n_exp = nrow(exp_data), n_clumped = nrow(exp_data), n_outcome = nrow(outcome_data), n_harmonized = n_harmonized)
} else {
  message("Too few IVs after harmonization.")
  write_status("DONE_TOO_FEW_HARMONIZED_IVS", n_fuma = length(snp_exp), n_exp = nrow(exp_data), n_clumped = nrow(exp_data), n_outcome = nrow(outcome_data), n_harmonized = n_harmonized)
}

message("Finished harmonization for ", exposure_key, " to ", outcome_key)
