#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=zip_fastGWA
#SBATCH --array=0-1
#SBATCH --time=0-03:59:00
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/cbica/home/wenju/output/zip_fastGWA_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/zip_fastGWA_%A_%a.err

set -euo pipefail
base_dir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock"

clock_list=(
    "Reproductive_female_proteomics_copd_clock"
    "heart_mri_copd_clock"
)

clock=${clock_list[$SLURM_ARRAY_TASK_ID]}

output_dir="${base_dir}/${clock}/fastGWA/output"
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
