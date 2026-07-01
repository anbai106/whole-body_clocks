#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=SA_clock_vs_BAG
#SBATCH --array=0-3572
#SBATCH --mem-per-cpu=12G
#SBATCH --time=0-05:59:00
#SBATCH --output=/cbica/home/wenju/output/SA_clock_vs_BAG_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/SA_clock_vs_BAG_%A_%a.err

set -euo pipefail

icd_list="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/SA/data/included_ICD_mortality_clock.tsv"
sa_data_dir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/SA/data"
out_dir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/SA/output_clock_vs_BAG/disease_free"
runner="/cbica/home/wenju/Project/whole-body_clocks/mortality_clock/survival_analysis/ICD10/disease_free/survival_analysis_clock_vs_bag.sh"

mkdir -p /cbica/home/wenju/output
mkdir -p "${out_dir}"

mapfile -t icd_array < <(awk 'NR > 1 {print $1}' "${icd_list}")

n_icd=${#icd_array[@]}
echo "Found ${n_icd} ICD endpoints."

if [[ "${SLURM_ARRAY_TASK_ID}" -ge "${n_icd}" ]]; then
  echo "ERROR: SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} but only ${n_icd} ICD endpoints found."
  exit 1
fi

icd="${icd_array[$SLURM_ARRAY_TASK_ID]}"
icd_tsv="${sa_data_dir}/${icd}_diagnosis_clock_disease_free.tsv"
output_tsv="${out_dir}/cox_compare_clock_vs_BAG_${icd}.tsv"

if [[ ! -f "${icd_tsv}" ]]; then
  echo "ERROR: missing ICD TSV: ${icd_tsv}"
  exit 1
fi

echo "SLURM job ID: ${SLURM_JOB_ID}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "ICD endpoint: ${icd}"
echo "Input TSV: ${icd_tsv}"
echo "Output TSV: ${output_tsv}"

if [[ ! -s "${output_tsv}" ]]; then
  bash "${runner}" "${icd_tsv}" "${output_tsv}"
else
  echo "Output already exists and is non-empty; skipping: ${output_tsv}"
fi
