#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=LDSC_main_mortality_clock
#SBATCH --array=0-21
#SBATCH --mem-per-cpu=12G
#SBATCH --time=0-05:59:00
#SBATCH --output=/cbica/home/wenju/output/LDSC_main_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/LDSC_main_%A_%a.err

set -euo pipefail

base_dir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/fastGWA/output"

# Collect all munged LDSC sumstats files in a stable sorted order
mapfile -t stat_files < <(
  find "${base_dir}" \
    -mindepth 2 \
    -maxdepth 2 \
    -type f \
    -name "organ_pheno_normalized_residualized.fastGWA.ldsc.sumstats.gz" \
  | sort
)

# Sanity check
n_files=${#stat_files[@]}
echo "Found ${n_files} LDSC sumstats files."

if [[ "${SLURM_ARRAY_TASK_ID}" -ge "${n_files}" ]]; then
  echo "ERROR: SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} but only ${n_files} files found."
  exit 1
fi

stat_gz1="${stat_files[$SLURM_ARRAY_TASK_ID]}"

organ_dir_path=$(dirname "${stat_gz1}")
organ_dir=$(basename "${organ_dir_path}")
organ="${organ_dir%_mortality_clock}"

# Create LDSC output subfolder inside the corresponding organ-clock folder
ldsc_out_dir="${organ_dir_path}/ldsc"
mkdir -p "${ldsc_out_dir}"

# Output prefix for LDSC main analysis
output_file="${ldsc_out_dir}/${organ}_vs"

echo "SLURM job ID: ${SLURM_JOB_ID}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Organ clock: ${organ}"
echo "Input sumstats: ${stat_gz1}"
echo "LDSC output directory: ${ldsc_out_dir}"
echo "Output prefix: ${output_file}"

bash /cbica/home/wenju/Project/UKBB_NMR_metabolomics/LDSC/ldsc_main_gc.sh "${stat_gz1}" "${output_file}"