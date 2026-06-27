#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=LDSC_munge_mortality_clock
#SBATCH --array=0-21
#SBATCH --mem-per-cpu=12G
#SBATCH --time=0-01:59:00
#SBATCH --output=/cbica/home/wenju/output/LDSC_munge_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/LDSC_munge_%A_%a.err

set -euo pipefail

base_dir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/fastGWA/output"

# Collect all LDSC-ready TSV files in a stable sorted order
mapfile -t files < <(
  find "${base_dir}" \
    -mindepth 2 \
    -maxdepth 2 \
    -type f \
    -name "organ_pheno_normalized_residualized.fastGWA.ldsc.tsv" \
  | sort
)

# Sanity check
n_files=${#files[@]}
echo "Found ${n_files} LDSC TSV files."

if [[ "${SLURM_ARRAY_TASK_ID}" -ge "${n_files}" ]]; then
  echo "ERROR: SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} but only ${n_files} files found."
  exit 1
fi

file="${files[$SLURM_ARRAY_TASK_ID]}"

organ_dir_path=$(dirname "${file}")
organ_dir=$(basename "${organ_dir_path}")
organ="${organ_dir%_mortality_clock}"

# Output prefix for ldsc_munge_stat.sh
# This will typically generate something like:
# organ_pheno_normalized_residualized.fastGWA.ldsc.sumstats.gz
output_file="${organ_dir_path}/organ_pheno_normalized_residualized.fastGWA.ldsc"

echo "SLURM job ID: ${SLURM_JOB_ID}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Organ clock: ${organ}"
echo "Input file: ${file}"
echo "Output prefix: ${output_file}"

bash /cbica/home/wenju/Project/whole-body_clocks/mortality_clock/LDSC/ldsc_munge_stat.sh "${file}" "${output_file}"
