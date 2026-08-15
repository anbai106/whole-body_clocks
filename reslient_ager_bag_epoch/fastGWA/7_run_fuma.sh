#!/bin/bash

remote_user="wenju"
remote_host="cubic-login.uphs.upenn.edu"
remote_base="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/fastGWA/output/Brain_proteomics_mortality_clock/EPOCH_BAG_residual"
local_dir="/Users/hao/Downloads"

mkdir -p "${local_dir}"


clock=Brain_proteomics_mortality_clock
echo "Copying ${clock}..."

remote_file="${remote_base}/organ_pheno_normalized_residualized.fastGWA.zip"
local_file="${local_dir}/${clock}_fuma.zip"

rsync -avz \
    "${remote_user}@${remote_host}:${remote_file}" \
    "${local_file}"

echo "All zipped fastGWA files copied to:"
echo "${local_dir}"

echo "Number of copied zip files:"
ls "${local_dir}"/*_fuma.zip | wc -l