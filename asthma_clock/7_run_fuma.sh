#!/bin/bash
set -euo pipefail

disease="asthma"
remote_user="wenju"
remote_host="cubic-login.uphs.upenn.edu"
remote_base="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock"
local_dir="/Users/hao/Downloads/${disease}_clock_fastGWA_fuma"
mkdir -p ${local_dir}

mkdir -p "${local_dir}"

clock_list=(
    "Digestive_metabolomics_${disease}_clock"
    "Hepatic_proteomics_${disease}_clock"
    "Reproductive_female_proteomics_${disease}_clock"
    "Endocrine_metabolomics_${disease}_clock"
    "Immune_metabolomics_${disease}_clock"
    "spleen_mri_${disease}_clock"
    "Endocrine_proteomics_${disease}_clock"
    "Immune_proteomics_${disease}_clock"
    "Hepatic_metabolomics_${disease}_clock"
    "Metabolic_metabolomics_${disease}_clock"
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