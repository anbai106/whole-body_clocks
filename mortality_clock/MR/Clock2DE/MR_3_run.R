#!/usr/bin/env Rscript

#### read the tsv contaiing the keywords for the search of exposure IVs
args = commandArgs(trailingOnly=TRUE)

pheno1 <- args[1]
pheno2 <- args[2]
output_dir_mr <- args[3]
harmonized_file <- args[4]

.libPaths('/gpfs/fs001/cbica/home/wenju/R/x86_64-pc-linux-gnu-library/4.3')
library(TwoSampleMR)
library(MRInstruments)
library(stringr)
library(dplyr)

print(paste0("Exposure variable: ", pheno1))
print(paste0("Outcome variable: ", pheno2))

hamonized_data <- read.table(harmonized_file, header=T, sep='\t', quote="")

### run MR
res_mr <- mr(hamonized_data)
write.table(res_mr,file=paste(output_dir_mr, "/MR_", pheno1, "_2_", pheno2, ".tsv", sep=""),row.names=F,col.names=T,sep="\t",quote=F)
res_mr_or <- generate_odds_ratios(res_mr)
write.table(res_mr_or,file=paste(output_dir_mr, "/MR_", pheno1, "_2_", pheno2, "_OR.tsv", sep=""),row.names=F,col.names=T,sep="\t",quote=F)

