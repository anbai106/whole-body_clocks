#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=zip_fastGWA
#SBATCH --array=0-21
#SBATCH --time=0-01:59:00
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/cbica/home/wenju/output/zip_fastGWA_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/zip_fastGWA_%A_%a.err

set -euo pipefail

base_dir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/fastGWA/output"

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

clock=${clock_list[$SLURM_ARRAY_TASK_ID]}

output_dir="${base_dir}/${clock}"
fastgwa="${output_dir}/organ_pheno_normalized_residualized.fastGWA"
zip_file="${output_dir}/organ_pheno_normalized_residualized.fastGWA.zip"

echo "SLURM array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Clock: ${clock}"
echo "Input fastGWA: ${fastgwa}"
echo "Output zip: ${zip_file}"

if [[ ! -d "${output_dir}" ]]; then
    echo "ERROR: output directory does not exist: ${output_dir}"
    exit 1
fi

if [[ ! -f "${fastgwa}" ]]; then
    echo "ERROR: fastGWA file does not exist: ${fastgwa}"
    exit 1
fi

# Remove incomplete old zip if present.
rm -f "${zip_file}"

# Compress the fastGWA file.
# -j stores only the file name inside the zip, not the full path.
# -9 uses maximum compression.
zip -9 -j "${zip_file}" "${fastgwa}"

# Test zip integrity.
zip -T "${zip_file}"

echo "Zip created successfully:"
ls -lh "${zip_file}"

echo "Done."
