#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=brainProt_EPOCH_vs_10bio_disease
#SBATCH --array=0-1
#SBATCH --mem-per-cpu=24G
#SBATCH --time=0-08:00:00
#SBATCH --output=/cbica/home/wenju/output/brainProt_EPOCH_vs_10bio_disease_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/brainProt_EPOCH_vs_10bio_disease_%A_%a.err

# ==============================================================================
# Brain proteomics mortality EPOCH vs 10 conventional biomarkers
# for incident disease onset.
#
# Array tasks:
#   0 = G309 = Alzheimer's disease (G30.9)
#   1 = I500 = heart failure (I50.0)
#
# Primary analysis uses the FULL mortality-EPOCH prediction sample (train + validation + test) to maximize sample size.
# Every comparison within a disease uses ONE identical sample complete for:
#   Brain proteomics mortality EPOCH + all 10 biomarkers + valid disease follow-up.
#
# Models include:
#   1) baseline covariates
#   2) baseline + EPOCH
#   3) baseline + each conventional biomarker separately
#   4) baseline + all 10 biomarkers as a panel
#   5) baseline + all 10 biomarkers + EPOCH
#
# The last comparison directly tests whether EPOCH adds information beyond the
# complete 10-biomarker panel.
# ==============================================================================

BASE_DIR="${BASE_DIR:-/cbica/home/wenju/Reproducibile_paper/WholeBodyClock}"
EPOCH_DIR="${EPOCH_DIR:-${BASE_DIR}/Brain_proteomics_mortality_clock}"
DISEASE_DIR="${DISEASE_DIR:-${BASE_DIR}/mortality_clock/SA/data}"

# Primary/default: full prediction table, which contains train/validation/test rows
# plus the EPOCH acceleration-z column. For a held-out sensitivity analysis,
# keep this same full file and submit with ANALYSIS_SPLIT=test.
EPOCH_PREDICTIONS="${EPOCH_PREDICTIONS:-${EPOCH_DIR}/brain_proteomics_mortality_clock_predictions.tsv}"
ANALYSIS_SPLIT="${ANALYSIS_SPLIT:-all}"

# Exact EPOCH acceleration-z column present in the full prediction table.
# Set explicitly by default so the raw Cox mortality risk score is never selected by mistake.
EPOCH_COL="${EPOCH_COL:-brain_proteomics_mortality_clock_acceleration_z}"

DEATH_XLSX="${DEATH_XLSX:-/cbica/home/wenju/Dataset/UKBB_UMelbourne/Death_related_var_from_Ye.xlsx}"
ID_MATCH_CSV="${ID_MATCH_CSV:-/cbica/home/wenju/Dataset/UKBB_UMelbourne/UKB_UMelbourne_vs_Penn_match_key.csv}"
COVARIATE_CSV="${COVARIATE_CSV:-/cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv}"
BIOMARKER_CSV="${BIOMARKER_CSV:-${BASE_DIR}/comparative_with_10_conventional_biomarkers/1_prepare_data/UKBB_fullsample_10_conventional_mortality_biomarkers.csv}"

PY_SCRIPT="${PY_SCRIPT:-/cbica/home/wenju/Project/whole-body_clocks/mortality_clock/Proteomics/comparative_with_10_conventional_biomarkers/brain_proteomics_EPOCH_vs_10_biomarkers_disease_onset.py}"
OUTPUT_ROOT="${OUTPUT_ROOT:-${EPOCH_DIR}/comparative_with_10_conventional_biomarkers_disease_onset}"

ADMIN_CENSOR_DATE="${ADMIN_CENSOR_DATE:-2022-11-30}"
PENALIZER="${PENALIZER:-0.0}"
MIN_EVENTS="${MIN_EVENTS:-20}"
N_BOOTSTRAP="${N_BOOTSTRAP:-1000}"
RANDOM_STATE="${RANDOM_STATE:-2026}"
CONDA_ENV="${CONDA_ENV:-survival}"
OVERWRITE="${OVERWRITE:-0}"

# Optional exact-column overrides. Usually unnecessary.
FIELD53_0_COL="${FIELD53_0_COL:-}"
DEATH_DATE_COL="${DEATH_DATE_COL:-}"
DEATH_ID_COL="${DEATH_ID_COL:-}"
IDMATCH_SCORE_COL="${IDMATCH_SCORE_COL:-}"
IDMATCH_DEATH_COL="${IDMATCH_DEATH_COL:-}"
COVARIATE_ID_COL="${COVARIATE_ID_COL:-}"
BIOMARKER_ID_COL="${BIOMARKER_ID_COL:-}"
COVARIATE_COLS="${COVARIATE_COLS:-}"

codes=(G309 I500)
labels=("Alzheimer's disease" "Heart failure")
files=(
  "${DISEASE_DIR}/G309_diagnosis_clock_disease_free.tsv"
  "${DISEASE_DIR}/I500_diagnosis_clock_disease_free.tsv"
)

idx="${SLURM_ARRAY_TASK_ID}"
if [[ "${idx}" -lt 0 || "${idx}" -ge "${#codes[@]}" ]]; then
  echo "ERROR: invalid SLURM_ARRAY_TASK_ID=${idx}"
  exit 1
fi

DISEASE_CODE="${codes[$idx]}"
DISEASE_LABEL="${labels[$idx]}"
DISEASE_TSV="${files[$idx]}"
OUTPUT_DIR="${OUTPUT_ROOT}/${ANALYSIS_SPLIT}/${DISEASE_CODE}"
mkdir -p /cbica/home/wenju/output
mkdir -p "${OUTPUT_DIR}"

prefix="brain_proteomics_mortality_EPOCH_vs_10_biomarkers_${DISEASE_CODE,,}_${ANALYSIS_SPLIT}_disease_onset"
EXPECTED_SUMMARY="${OUTPUT_DIR}/${prefix}_individual_predictor_summary.tsv"

required_files=(
  "${EPOCH_PREDICTIONS}"
  "${DISEASE_TSV}"
  "${DEATH_XLSX}"
  "${ID_MATCH_CSV}"
  "${COVARIATE_CSV}"
  "${BIOMARKER_CSV}"
  "${PY_SCRIPT}"
)

for f in "${required_files[@]}"; do
  if [[ ! -f "${f}" ]]; then
    echo "ERROR: required file not found: ${f}"
    exit 1
  fi
done

if [[ -s "${EXPECTED_SUMMARY}" && "${OVERWRITE}" != "1" ]]; then
  echo "Output already exists: ${EXPECTED_SUMMARY}"
  echo "Set OVERWRITE=1 to rerun."
  exit 0
fi

# CUBIC environment used by the existing survival analyses.
module load python/anaconda/3
source activate "${CONDA_ENV}"

extra_args=()
[[ -n "${EPOCH_COL}" ]] && extra_args+=(--epoch-col "${EPOCH_COL}")
[[ -n "${FIELD53_0_COL}" ]] && extra_args+=(--field53-0-col "${FIELD53_0_COL}")
[[ -n "${DEATH_DATE_COL}" ]] && extra_args+=(--death-date-col "${DEATH_DATE_COL}")
[[ -n "${DEATH_ID_COL}" ]] && extra_args+=(--death-id-col "${DEATH_ID_COL}")
[[ -n "${IDMATCH_SCORE_COL}" ]] && extra_args+=(--idmatch-score-col "${IDMATCH_SCORE_COL}")
[[ -n "${IDMATCH_DEATH_COL}" ]] && extra_args+=(--idmatch-death-col "${IDMATCH_DEATH_COL}")
[[ -n "${COVARIATE_ID_COL}" ]] && extra_args+=(--covariate-id-col "${COVARIATE_ID_COL}")
[[ -n "${BIOMARKER_ID_COL}" ]] && extra_args+=(--biomarker-id-col "${BIOMARKER_ID_COL}")
[[ -n "${COVARIATE_COLS}" ]] && extra_args+=(--covariate-cols "${COVARIATE_COLS}")

echo "================================================================================"
echo "Brain proteomics mortality EPOCH vs 10 conventional biomarkers"
echo "Disease: ${DISEASE_LABEL} (${DISEASE_CODE})"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "EPOCH predictions: ${EPOCH_PREDICTIONS}"
echo "Analysis split: ${ANALYSIS_SPLIT}"
if [[ "${ANALYSIS_SPLIT}" == "all" ]]; then
  echo "Population note: FULL mortality-EPOCH sample (train + validation + test); not a strictly held-out validation."
else
  echo "Population note: held-out mortality-EPOCH test split sensitivity analysis."
fi
echo "EPOCH column: ${EPOCH_COL}"
echo "Disease TSV: ${DISEASE_TSV}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Started at: $(date)"
echo "================================================================================"

set +e
python "${PY_SCRIPT}" \
  --epoch-predictions "${EPOCH_PREDICTIONS}" \
  --analysis-split "${ANALYSIS_SPLIT}" \
  --disease-tsv "${DISEASE_TSV}" \
  --disease-code "${DISEASE_CODE}" \
  --disease-label "${DISEASE_LABEL}" \
  --death-xlsx "${DEATH_XLSX}" \
  --id-match-csv "${ID_MATCH_CSV}" \
  --covariate-csv "${COVARIATE_CSV}" \
  --biomarker-csv "${BIOMARKER_CSV}" \
  --output-dir "${OUTPUT_DIR}" \
  --admin-censor-date "${ADMIN_CENSOR_DATE}" \
  --penalizer "${PENALIZER}" \
  --min-events "${MIN_EVENTS}" \
  --n-bootstrap "${N_BOOTSTRAP}" \
  --random-state "$((RANDOM_STATE + idx * 10000))" \
  "${extra_args[@]}"
status=$?
set -e

if [[ "${status}" -ne 0 ]]; then
  echo "================================================================================"
  echo "ERROR: analysis failed for ${DISEASE_LABEL} (${DISEASE_CODE})"
  echo "Python exit code: ${status}"
  echo "Failed at: $(date)"
  echo "================================================================================"
  conda deactivate || true
  exit "${status}"
fi

echo "================================================================================"
echo "SUCCESS: ${DISEASE_LABEL} (${DISEASE_CODE})"
echo "Main summary: ${EXPECTED_SUMMARY}"
echo "Finished at: $(date)"
echo "================================================================================"

conda deactivate