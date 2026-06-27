#!/bin/bash

for organ in Endocrine Digestive Hepatic Immune Metabolic
do
    echo ${organ}
    rsync -avz wenju@cubic-login.uphs.upenn.edu:/cbica/home/wenju/Reproducibile_paper/UKBB_metabolomics/fastGWA_MetBAG/output/${organ}/organ_pheno_normalized_residualized.fastGWA.zip /Users/hao/${organ}_fuma.zip
done