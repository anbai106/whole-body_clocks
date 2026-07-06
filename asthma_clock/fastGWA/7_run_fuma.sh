#!/bin/bash
set -euo pipefail

remote_user="wenju"
remote_host="cubic-login.uphs.upenn.edu"
remote_base="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/fastGWA/output"
local_dir="/Users/hao/Downloads/mortality_clock_fastGWA_fuma"

mkdir -p "${local_dir}"

clock_list=(
    adipose_mri_mortality_clock
    brain_mri_mortality_clock
    Brain_proteomics_mortality_clock
    Digestive_metabolomics_mortality_clock
    Endocrine_metabolomics_mortality_clock
    Endocrine_proteomics_mortality_clock
    Eye_proteomics_mortality_clock
    heart_mri_mortality_clock
    Heart_proteomics_mortality_clock
    Hepatic_metabolomics_mortality_clock
    Hepatic_proteomics_mortality_clock
    Immune_metabolomics_mortality_clock
    Immune_proteomics_mortality_clock
    kidney_mri_mortality_clock
    liver_mri_mortality_clock
    pancreas_mri_mortality_clock
    Pulmonary_proteomics_mortality_clock
    Renal_proteomics_mortality_clock
    Reproductive_female_proteomics_mortality_clock
    Reproductive_male_proteomics_mortality_clock
    Skin_proteomics_mortality_clock
    spleen_mri_mortality_clock
)

for clock in "${clock_list[@]}"; do
    echo "Copying ${clock}..."

    remote_file="${remote_base}/${clock}/organ_pheno_normalized_residualized.fastGWA.zip"
    local_file="${local_dir}/${clock}_fuma.zip"

    rsync -avz \
        "${remote_user}@${remote_host}:${remote_file}" \
        "${local_file}"
done

echo "All zipped fastGWA files copied to:"
echo "${local_dir}"

echo "Number of copied zip files:"
ls "${local_dir}"/*_fuma.zip | wc -l