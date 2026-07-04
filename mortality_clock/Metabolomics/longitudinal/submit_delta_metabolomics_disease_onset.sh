#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=delta_met_disease
#SBATCH --array=0-4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=0-06:00:00
#SBATCH --output=/cbica/home/wenju/output/delta_met_disease_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/delta_met_disease_%A_%a.err

set -euo pipefail

# User listed asthma twice; this array runs the five unique algorithmically-defined endpoints.
# UKB Category 42 fields used by the Python script:
#   all_cause_dementia: 42018
#   asthma: 42014
#   myocardial_infarction: 42000
#   copd: 42016
#   stroke: 42006
endpoints=(
  all_cause_dementia
  asthma
  myocardial_infarction
  copd
  stroke
)

out_dir=${OUT_DIR:-/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/SA/delta_metabolomics_algorithmic_disease_onset}
mkdir -p /cbica/home/wenju/output
mkdir -p "${out_dir}"

runner=/cbica/home/wenju/Project/whole-body_clocks/mortality_clock/Metabolomics/longitudinal/survival_analysis_delta_metabolomics_disease_onset.sh

n_endpoints=${#endpoints[@]}
if [[ "${SLURM_ARRAY_TASK_ID}" -ge "${n_endpoints}" ]]; then
  echo "ERROR: SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} but only ${n_endpoints} endpoints are configured."
  exit 1
fi

endpoint="${endpoints[$SLURM_ARRAY_TASK_ID]}"
output_tsv="${out_dir}/cox_delta_metabolomics_${endpoint}.tsv"

echo "SLURM job ID: ${SLURM_JOB_ID}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Endpoint: ${endpoint}"
echo "Output TSV: ${output_tsv}"
echo "Runner: ${runner}"
echo "OUT_DIR: ${out_dir}"

force_rerun=${FORCE_RERUN:-1}

if [[ "${force_rerun}" == "1" ]]; then
  echo "FORCE_RERUN=1; overwriting existing output if present: ${output_tsv}"
  rm -f "${output_tsv}"
  bash "${runner}" "${endpoint}" "${output_tsv}"
elif [[ ! -s "${output_tsv}" ]]; then
  bash "${runner}" "${endpoint}" "${output_tsv}"
else
  echo "Output already exists and is non-empty; skipping: ${output_tsv}"
fi
