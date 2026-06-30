#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

pheno1 <- args[1]
pheno2 <- args[2]
output_dir <- args[3]

.libPaths('/gpfs/fs001/cbica/home/wenju/R/x86_64-pc-linux-gnu-library/4.3')

library(TwoSampleMR)
library(MRInstruments)
library(stringr)
library(dplyr)
library(svglite)

## output dir
print(paste0("Exposure: ", pheno1))
print(paste0("Outcome: ", pheno2))
hamonized_tsv = paste(output_dir, pheno1,  paste('harmonized_data_', pheno1, '_2_', pheno2, '.tsv', sep=""), sep = '/')
hamonized_data <- read.table(hamonized_tsv, header=T, sep='\t', quote="")
mr_tsv = paste(output_dir, pheno1, paste('MR_', pheno1, '_2_', pheno2, '.tsv', sep=""), sep = '/')
res_mr <- read.table(mr_tsv, header=T, sep='\t', quote="")

#### Sensitivity analyses
# heterogeneity test
res_heterogeneity <- mr_heterogeneity(hamonized_data)
write.table(res_heterogeneity,file=paste(output_dir, '/', pheno1,  "/SC_", pheno1, "_2_", pheno2, "_heterogeneity.tsv", sep=""),row.names=F,col.names=T,sep="\t",quote=F)

# horizontal pleiotropy
res_pleiotropy <- mr_pleiotropy_test(hamonized_data)
write.table(res_pleiotropy,file=paste(output_dir, '/', pheno1,  "/SC_", pheno1, "_2_", pheno2, "_horizontal_pleiotropy.tsv", sep=""),row.names=F,col.names=T,sep="\t",quote=F)

# single SNP analyses
res_single <- mr_singlesnp(hamonized_data)
write.table(res_single,file=paste(output_dir, '/', pheno1,  "/SC_", pheno1, "_2_", pheno2, "_single_snp.tsv", sep=""),row.names=F,col.names=T,sep="\t",quote=F)

## LOO analysis
loo <- mr_leaveoneout(hamonized_data)
write.table(loo,file=paste(output_dir, '/', pheno1,  "/SC_", pheno1, "_2_", pheno2, "_LOO.tsv", sep=""),row.names=F,col.names=T,sep="\t",quote=F)

## forest plot
svglite(paste(output_dir, '/', pheno1,  "/SC_", pheno1, "_2_", pheno2, "_forest.svg", sep=""), width = 10, height = 10)
mr_forest_plot(res_single)
dev.off()

## loo plot
svglite(paste(output_dir, '/', pheno1,  "/SC_", pheno1, "_2_", pheno2, "_loo.svg", sep=""), width = 10, height = 10)
mr_leaveoneout_plot(loo)
dev.off()

## funnel plot
svglite(paste(output_dir, '/', pheno1,  "/SC_", pheno1, "_2_", pheno2, "_funnel.svg", sep=""), width = 10, height = 10)
mr_funnel_plot(res_single)
dev.off()

## scatter plot
svglite(paste(output_dir, '/', pheno1,  "/SC_", pheno1, "_2_", pheno2, "_scatter.svg", sep=""), width = 10, height = 10)
mr_scatter_plot(res_mr, hamonized_data)
dev.off()

print('Plots finished here...')
