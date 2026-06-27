#!/bin/bash

for organ in Endocrine Digestive Hepatic Immune Metabolic
do
  for sex in female male
  do
      in_file="/cbica/home/wenju/Reproducibile_paper/UKBB_metabolomics/fastGWA_MetBAG/output/${organ}_non_derived_${sex}/BAG_pheno_normalized_residualized.fastGWA"
      out_file="/cbica/home/wenju/Reproducibile_paper/FemaleAgingClock/fastGWA/output/${organ}_MetBAG_${sex}.fastGWA"
      # create/overwrite symlink: out_file -> in_file
      ln -sfn "$in_file" "$out_file"
  done
done