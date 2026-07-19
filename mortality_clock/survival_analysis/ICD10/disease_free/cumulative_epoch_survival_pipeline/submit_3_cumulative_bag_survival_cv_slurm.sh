#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=SA_cum_BAG_PM
#SBATCH --array=0-3572
#SBATCH --mem-per-cpu=12G
#SBATCH --time=0-05:59:00
#SBATCH --output=/cbica/home/wenju/output/SA_cum_BAG_PM_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/SA_cum_BAG_PM_%A_%a.err

set -euo pipefail

# STEP 3B cumulative survival analysis for the matched BAGs corresponding to
# the 11 proteomics + 4 metabolomics mortality EPOCH clocks.
#
# Default behavior:
#   - Rank BAGs by their own individual BAG HR/P values from clock-vs-BAG results.
#   - Run Option 1 fixed-order 5-fold CV.
#   - Require complete cases for all 15 BAGs and all 15 corresponding EPOCH clocks
#     so the BAG analysis is paired with the EPOCH analysis sample whenever possible.
#
# Useful overrides:
#   RANK_SOURCE=epoch                  use EPOCH clock order instead of BAG order
#   PAIR_COMPLETE_CASE_WITH_EPOCH=0    use BAG-only complete cases
#   CV_FOLDS=0                         disable CV
#   OVERWRITE=1                        rerun existing endpoints

ICD_LIST="${ICD_LIST:-/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/SA/data/included_ICD_mortality_clock.tsv}"
SA_DATA_DIR="${SA_DATA_DIR:-/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/SA/data}"
INDIVIDUAL_RESULT_DIR="${INDIVIDUAL_RESULT_DIR:-/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/SA/output_clock_vs_BAG}"
OUT_DIR="${OUT_DIR:-/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/SA/output_cumulative_BAG_PM/disease_free_cv}"
RUNNER="${RUNNER:-/cbica/home/wenju/Project/whole-body_clocks/mortality_clock/survival_analysis/ICD10/disease_free/cumulative_epoch_survival_pipeline/run_3_cumulative_bag_survival_cv.sh}"
PYTHON_SCRIPT="${PYTHON_SCRIPT:-/cbica/home/wenju/Project/whole-body_clocks/mortality_clock/survival_analysis/ICD10/disease_free/cumulative_epoch_survival_pipeline/survival_analysis_cumulative_bag_cv.py}"
BAG_TSV="${BAG_TSV:-/cbica/home/wenju/Reproducibile_paper/SleepAging/data/MomoBAG.tsv}"

DATA_SUFFIX="${DATA_SUFFIX:-_diagnosis_clock_disease_free.tsv}"
INDIVIDUAL_PREFIX="${INDIVIDUAL_PREFIX:-cox_compare_clock_vs_BAG_}"

RANK_P_THRESHOLD="${RANK_P_THRESHOLD:-0.0033333333}"
RANK_SORT_MODE="${RANK_SORT_MODE:-hr_desc}"
NON_SIG_ORDER="${NON_SIG_ORDER:-random}"
RANK_SOURCE="${RANK_SOURCE:-bag}"
PAIR_COMPLETE_CASE_WITH_EPOCH="${PAIR_COMPLETE_CASE_WITH_EPOCH:-1}"
CV_FOLDS="${CV_FOLDS:-5}"

MIN_CASE="${MIN_CASE:-20}"
MIN_NONCASE="${MIN_NONCASE:-20}"
PENALIZER="${PENALIZER:-0.0}"
COVARIATE_SET="${COVARIATE_SET:-clinical}"
ADMIN_CENSOR_DATE="${ADMIN_CENSOR_DATE:-2022-11-30}"
OVERWRITE="${OVERWRITE:-0}"

mkdir -p /cbica/home/wenju/output
mkdir -p "${OUT_DIR}"

if [[ ! -r "${ICD_LIST}" ]]; then
    echo "ERROR: cannot read ICD list: ${ICD_LIST}" >&2
    exit 1
fi
if [[ ! -r "${RUNNER}" ]]; then
    echo "ERROR: cannot read runner: ${RUNNER}" >&2
    exit 1
fi
if [[ ! -r "${PYTHON_SCRIPT}" ]]; then
    echo "ERROR: cannot read Python script: ${PYTHON_SCRIPT}" >&2
    exit 1
fi
if [[ ! -r "${BAG_TSV}" ]]; then
    echo "ERROR: cannot read BAG TSV: ${BAG_TSV}" >&2
    exit 1
fi

mapfile -t ICD_ARRAY < <(awk -F '\t' 'NR > 1 && $1 != "" {gsub(/\r/, "", $1); print $1}' "${ICD_LIST}")

TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is not set}"
if (( TASK_ID < 0 || TASK_ID >= ${#ICD_ARRAY[@]} )); then
    echo "ERROR: array index ${TASK_ID} outside 0-$(( ${#ICD_ARRAY[@]} - 1 ))." >&2
    exit 1
fi

ICD="${ICD_ARRAY[$TASK_ID]}"
ICD_TSV="${SA_DATA_DIR}/${ICD}${DATA_SUFFIX}"
INDIVIDUAL_TSV="${INDIVIDUAL_RESULT_DIR}/${INDIVIDUAL_PREFIX}${ICD}.tsv"
OUTPUT_TSV="${OUT_DIR}/cox_cumulative_BAG_PM_${ICD}.tsv"
RANK_TSV="${OUT_DIR}/rank_order_BAG_PM_${ICD}.tsv"
AUDIT_TSV="${OUT_DIR}/audit_BAG_PM_${ICD}.tsv"

if [[ ! -s "${ICD_TSV}" ]]; then
    echo "ERROR: missing or empty ICD TSV: ${ICD_TSV}" >&2
    exit 1
fi
if [[ ! -s "${INDIVIDUAL_TSV}" ]]; then
    echo "ERROR: missing or empty individual clock-vs-BAG result TSV: ${INDIVIDUAL_TSV}" >&2
    exit 1
fi

if [[ "${OVERWRITE}" != "1" && -s "${OUTPUT_TSV}" && -s "${RANK_TSV}" ]]; then
    echo "Outputs already exist for ${ICD}; skipping. Set OVERWRITE=1 to rerun."
    exit 0
fi

echo "SLURM job ID: ${SLURM_JOB_ID}"
echo "Array task ID: ${TASK_ID}"
echo "ICD endpoint: ${ICD}"
echo "ICD TSV: ${ICD_TSV}"
echo "Individual result TSV: ${INDIVIDUAL_TSV}"
echo "BAG TSV: ${BAG_TSV}"
echo "Output TSV: ${OUTPUT_TSV}"
echo "Rank TSV: ${RANK_TSV}"
echo "Rank source: ${RANK_SOURCE}"
echo "Rank P threshold: ${RANK_P_THRESHOLD}"
echo "Rank sort mode: ${RANK_SORT_MODE}"
echo "Non-significant order: ${NON_SIG_ORDER}"
echo "CV folds: ${CV_FOLDS}"
echo "Paired complete case with EPOCH: ${PAIR_COMPLETE_CASE_WITH_EPOCH}"
echo "Covariate set: ${COVARIATE_SET}"
echo "Administrative censor date: ${ADMIN_CENSOR_DATE}"

export PYTHON_SCRIPT BAG_TSV RANK_P_THRESHOLD RANK_SORT_MODE NON_SIG_ORDER RANK_SOURCE PAIR_COMPLETE_CASE_WITH_EPOCH CV_FOLDS MIN_CASE MIN_NONCASE PENALIZER COVARIATE_SET ADMIN_CENSOR_DATE OVERWRITE

bash "${RUNNER}" "${ICD_TSV}" "${INDIVIDUAL_TSV}" "${OUTPUT_TSV}" "${RANK_TSV}" "${AUDIT_TSV}"

echo "Finished ${ICD}."
