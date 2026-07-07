#!/bin/bash
set -euo pipefail

disease="copd"
remote_user="wenju"
remote_host="cubic-login.uphs.upenn.edu"
remote_base="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock"
local_dir="/Users/hao/Downloads/${disease}_clock_fastGWA_fuma"
mkdir -p ${local_dir}

mkdir -p "${local_dir}"

clock_list=(
    "Reproductive_female_proteomics_copd_clock"
    "heart_mri_copd_clock"
)

for clock in "${clock_list[@]}"; do
    echo "Copying ${clock}..."

    remote_file="${remote_base}/${clock}/fastGWA/output/organ_pheno_normalized_residualized.fastGWA.zip"
    local_file="${local_dir}/${clock}_fuma.zip"

    rsync -avz \
        "${remote_user}@${remote_host}:${remote_file}" \
        "${local_file}"
done

echo "All zipped fastGWA files copied to:"
echo "${local_dir}"

echo "Number of copied zip files:"
ls "${local_dir}"/*_fuma.zip | wc -l