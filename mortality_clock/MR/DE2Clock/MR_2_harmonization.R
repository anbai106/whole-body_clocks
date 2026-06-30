#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

exposure_gwas_tsv <- args[1]
plinkn_clumped_lead_snp_tsv <- args[2]
exposure_key <- args[3]
output_2sampleMR_tsv <- args[4]
outcome_key <- args[5]
output_dir <- args[6]

.libPaths('/gpfs/fs001/cbica/home/wenju/R/x86_64-pc-linux-gnu-library/4.3')
library(TwoSampleMR)
library(MRInstruments)
library(stringr)
library(dplyr)
library(ieugwasr)

### read exposure data from local GWAS
snp_exp <- read.table(plinkn_clumped_lead_snp_tsv, header=T)$SNP
exp_data <- read_outcome_data(
  exposure_gwas_tsv,
  snps = snp_exp,
  sep = "\t",
  snp_col = "ID",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "eaf_A1",
  pval_col = "P",
  samplesize_col = "N",
  phenotype_col = 'pheno',
  pos_col = "POS")
exp_data <- exp_data[,c("pos.outcome", "SNP", "other_allele.outcome", "effect_allele.outcome", "samplesize.outcome", "beta.outcome", "se.outcome", "pval.outcome",
                                          "eaf.outcome", "outcome")]
exp_data <- format_data(exp_data,
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
                                 type="exposure")
#### use cloud server
# exp_data <- clump_data(exp_data)

### run local
names(exp_data)[names(exp_data) == "SNP"] <- "rsid"
names(exp_data)[names(exp_data) == "pval.exposure"] <- "pval"
exp_data <- ld_clump(
  exp_data,
  plink_bin = "/cbica/home/wenju/Software/plink_1.90_beta_linux_x86_64_20210606/plink",
  bfile = "/cbica/home/wenju/Dataset/2SampleMR_1000Genome_LD/EUR")
names(exp_data)[names(exp_data) == "rsid"] <- "SNP"
names(exp_data)[names(exp_data) == "pval"] <- "pval.exposure"
exp_data <- exp_data[ ,  !names(exp_data) %in% c("id")]

table(exp_data$exposure)

### process the outcome variables from harmonized PSC
outcome_data <- read_outcome_data(
  output_2sampleMR_tsv,
  snps = exp_data$SNP,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "AF1",
  pval_col = "P",
  samplesize_col = "N",
  pos_col = "POS")

outcome_data<-subset(outcome_data, mr_keep.outcome=='TRUE')
exp_data<-subset(exp_data, mr_keep.exposure=='TRUE')

### harmonize the data
harmonized_data <- harmonise_data(
  exposure_dat = exp_data,
  outcome_dat = outcome_data
)

if(dim(harmonized_data)[1] > 7){
  ### save harmonized exposure and outcome into tsv
  print(paste0("After harmonizing data, there exist N IVs: ", dim(harmonized_data)[1]))
  write.table(harmonized_data,file=paste(output_dir, "/harmonized_data_", exposure_key, "_2_", outcome_key, ".tsv", sep=""),row.names=F,col.names=T,sep="\t",quote=F)
  write.table('DONE', file=paste(output_dir, "/DONE_", exposure_key, "_2_", outcome_key, ".tsv", sep=""), quote=FALSE, sep='\t', col.names = FALSE)
} else {
  write.table('DONE', file=paste(output_dir, "/DONE_", exposure_key, "_2_", outcome_key, ".tsv", sep=""), quote=FALSE, sep='\t', col.names = FALSE)
  message("There is no enought IVs from the exposure variables...")
}
print("Here...")