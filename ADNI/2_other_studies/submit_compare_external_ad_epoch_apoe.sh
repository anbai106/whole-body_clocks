#!/usr/bin/env bash
#SBATCH --job-name=apoe_ad_epoch_external
#SBATCH --output=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch/apoe_external/logs/apoe_ad_epoch_%A_%a.out
#SBATCH --error=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch/apoe_external/logs/apoe_ad_epoch_%A_%a.err
#SBATCH --array=0-3
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G

# ==============================================================================
# APOE E4/E4 versus E2/E2 association with baseline harmonized AD EPOCH
#
# Array tasks:
#   0 = AIBL
#   1 = BLSA
#   2 = OASIS
#   3 = pooled AIBL + BLSA + OASIS
# ==============================================================================

PROJECT_ROOT="${PROJECT_ROOT:-/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch}"
OUTDIR="${OUTDIR:-${PROJECT_ROOT}/apoe_external}"
LOGDIR="${OUTDIR}/logs"
RESULTDIR="${OUTDIR}/results"

SAMPLE_FILE="${SAMPLE_FILE:-${PROJECT_ROOT}/external_5_studies_istaging.tsv}"
PREDICTION_FILE="${PREDICTION_FILE:-${PROJECT_ROOT}/results_external_longitudinal_ad_epoch_harmonized/external_5_studies_adni_brain_mri_ad_epoch_harmonized_scan_level_predictions.tsv}"

mkdir -p "${LOGDIR}" "${RESULTDIR}"

ANALYSIS_GROUPS=(
  "AIBL"
  "BLSA"
  "OASIS"
  "POOLED"
)

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"

if (( TASK_ID < 0 || TASK_ID >= ${#ANALYSIS_GROUPS[@]} )); then
  echo "ERROR: invalid SLURM_ARRAY_TASK_ID=${TASK_ID}" >&2
  exit 1
fi

ANALYSIS_GROUP="${ANALYSIS_GROUPS[${TASK_ID}]}"
TASK_RESULTDIR="${RESULTDIR}/${ANALYSIS_GROUP}"
mkdir -p "${TASK_RESULTDIR}"

for required_file in "${PY_SCRIPT}" "${SAMPLE_FILE}" "${PREDICTION_FILE}"; do
  if [[ ! -f "${required_file}" ]]; then
    echo "ERROR: required file does not exist: ${required_file}" >&2
    exit 1
  fi
done

echo "========================================================================"
echo "SLURM job ID:       ${SLURM_JOB_ID:-NA}"
echo "Array task ID:      ${TASK_ID}"
echo "Analysis group:     ${ANALYSIS_GROUP}"
echo "Sample file:        ${SAMPLE_FILE}"
echo "Prediction file:    ${PREDICTION_FILE}"
echo "Result directory:   ${TASK_RESULTDIR}"
echo "========================================================================"

source activate DNE

python "${PY_SCRIPT}" \
  --sample-file "${SAMPLE_FILE}" \
  --prediction-file "${PREDICTION_FILE}" \
  --outdir "${TASK_RESULTDIR}" \
  --analysis-group "${ANALYSIS_GROUP}" \
  --studies AIBL BLSA OASIS \
  --primary-outcome "ad_epoch_acceleration_years"

echo "Finished APOE analysis for ${ANALYSIS_GROUP}."

if command -v conda >/dev/null 2>&1; then
  conda deactivate || true
fi
