#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=LDSC_mortality_clock
#SBATCH --array=0-21
#SBATCH --mem-per-cpu=12G
#SBATCH --time=0-01:59:00
#SBATCH --output=/cbica/home/wenju/output/LDSC_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/LDSC_%A_%a.err

set -euo pipefail

base_dir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/fastGWA/output"

# Collect all fastGWA files in a stable sorted order
mapfile -t source_files < <(
  find "${base_dir}" \
    -mindepth 2 \
    -maxdepth 2 \
    -type f \
    -name "organ_pheno_normalized_residualized.fastGWA" \
  | sort
)

# Sanity check
n_files=${#source_files[@]}
echo "Found ${n_files} fastGWA files."

if [[ "${SLURM_ARRAY_TASK_ID}" -ge "${n_files}" ]]; then
  echo "ERROR: SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} but only ${n_files} files found."
  exit 1
fi

source_file="${source_files[$SLURM_ARRAY_TASK_ID]}"

organ_dir=$(basename "$(dirname "${source_file}")")
organ="${organ_dir%_mortality_clock}"

echo "SLURM job ID: ${SLURM_JOB_ID}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Organ clock: ${organ}"
echo "Source file: ${source_file}"

bash /cbica/home/wenju/Project/whole-body_clocks/mortality_clock/LDSC/ldsc_prepare_data.sh "${source_file}"