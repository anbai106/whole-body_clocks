#!/bin/bash
set -euo pipefail

# Runner for one ICD endpoint: matched BAG cumulative CV analysis.
# Usage:
#   bash run_3_cumulative_bag_survival_cv.sh ICD_TSV INDIVIDUAL_RESULT_TSV OUTPUT_TSV RANK_TSV AUDIT_TSV

module load python/anaconda/3

ICD_TSV=$1
INDIVIDUAL_RESULT_TSV=$2
OUTPUT_TSV=$3
RANK_TSV=$4
AUDIT_TSV=$5

PYTHON_SCRIPT="${PYTHON_SCRIPT:-/cbica/home/wenju/Project/whole-body_clocks/mortality_clock/survival_analysis/ICD10/disease_free/cumulative_epoch_survival_pipeline/survival_analysis_cumulative_bag_cv.py}"
BAG_TSV="${BAG_TSV:-/cbica/home/wenju/Reproducibile_paper/SleepAging/data/MomoBAG.tsv}"

RANK_P_THRESHOLD="${RANK_P_THRESHOLD:-0.0033333333}"
RANK_SORT_MODE="${RANK_SORT_MODE:-hr_desc}"
NON_SIG_ORDER="${NON_SIG_ORDER:-random}"
MIN_CASE="${MIN_CASE:-20}"
MIN_NONCASE="${MIN_NONCASE:-20}"
PENALIZER="${PENALIZER:-0.0}"
COVARIATE_SET="${COVARIATE_SET:-clinical}"
ADMIN_CENSOR_DATE="${ADMIN_CENSOR_DATE:-2022-11-30}"
OVERWRITE="${OVERWRITE:-0}"

# CV / BAG-specific behavior is read by the Python script from environment.
CV_FOLDS="${CV_FOLDS:-5}"
RANK_SOURCE="${RANK_SOURCE:-bag}"
PAIR_COMPLETE_CASE_WITH_EPOCH="${PAIR_COMPLETE_CASE_WITH_EPOCH:-1}"

source activate survival

ARGS=(
    --icd_tsv "${ICD_TSV}"
    --individual_result_tsv "${INDIVIDUAL_RESULT_TSV}"
    --output_tsv "${OUTPUT_TSV}"
    --rank_output_tsv "${RANK_TSV}"
    --audit_tsv "${AUDIT_TSV}"
    --bag_tsv "${BAG_TSV}"
    --rank_p_threshold "${RANK_P_THRESHOLD}"
    --rank_sort_mode "${RANK_SORT_MODE}"
    --non_sig_order "${NON_SIG_ORDER}"
    --min_case "${MIN_CASE}"
    --min_noncase "${MIN_NONCASE}"
    --penalizer "${PENALIZER}"
    --covariate_set "${COVARIATE_SET}"
    --admin_censor_date "${ADMIN_CENSOR_DATE}"
)

if [[ "${OVERWRITE}" == "1" ]]; then
    ARGS+=(--overwrite)
fi

export CV_FOLDS RANK_SOURCE PAIR_COMPLETE_CASE_WITH_EPOCH BAG_TSV
python -u "${PYTHON_SCRIPT}" "${ARGS[@]}"

conda deactivate
