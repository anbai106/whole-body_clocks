#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=mi_immune_EPOCH_vs_10bio
#SBATCH --mem-per-cpu=24G
#SBATCH --time=05:00:00
#SBATCH --output=/cbica/home/wenju/output/mi_immune_EPOCH_vs_10bio_%j.out
#SBATCH --error=/cbica/home/wenju/output/mi_immune_EPOCH_vs_10bio_%j.err

# ==============================================================================
# Apple-to-apple landmark survival comparison for all-cause mortality:
#
#   MI immune-metabolomics disease-specific EPOCH
#       versus
#   10 conventional UK Biobank mortality biomarkers
#
# Every Cox model is fitted in ONE IDENTICAL population:
#   complete EPOCH + all 10 biomarkers + valid mortality follow-up.
#
# Metabolomics landmark:
#   UKB Field 53 instance 0_0
#
# Default common baseline covariates:
#   age at assessment
#   sex
#   genetic ethnic grouping
#   assessment centre
#   smoking status
#   BMI
#
# Systolic/diastolic BP are intentionally excluded from the default baseline
# because systolic BP is itself one of the 10 benchmark predictors.
#
# One separate model is fitted per predictor:
#   baseline covariates + predictor_z
#
# No EPOCH + biomarker joint model is fitted.
# ==============================================================================

# ------------------------------------------------------------------------------
# Core paths from the existing 47-clock mortality analysis
# ------------------------------------------------------------------------------

BASE_DIR="${BASE_DIR:-/gpfs/fs001/cbica/home/wenju/Reproducibile_paper/WholeBodyClock}"

ANALYSIS_DIR="${ANALYSIS_DIR:-${BASE_DIR}/all_disease_lepoch_incremental_value_scale_qc}"

SCORE_WIDE_TSV="${SCORE_WIDE_TSV:-${ANALYSIS_DIR}/stable_significant_disease_clock_acceleration_z_wide.tsv}"

SCORE_METADATA_TSV="${SCORE_METADATA_TSV:-${ANALYSIS_DIR}/stable_significant_disease_clock_acceleration_z_metadata.tsv}"

GOOD_CLOCK_TSV="${GOOD_CLOCK_TSV:-${ANALYSIS_DIR}/all_disease_lepoch_main_text_good_clocks.tsv}"

DEATH_XLSX="${DEATH_XLSX:-/cbica/home/wenju/Dataset/UKBB_UMelbourne/Death_related_var_from_Ye.xlsx}"

ID_MATCH_CSV="${ID_MATCH_CSV:-/cbica/home/wenju/Dataset/UKBB_UMelbourne/UKB_UMelbourne_vs_Penn_match_key.csv}"

COVARIATE_CSV="${COVARIATE_CSV:-/cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv}"

# ------------------------------------------------------------------------------
# 10 conventional mortality biomarkers prepared from the full UKB sample
# ------------------------------------------------------------------------------

BIOMARKER_CSV="${BIOMARKER_CSV:-${BASE_DIR}/comparative_with_10_conventional_biomarkers/1_prepare_data/UKBB_fullsample_10_conventional_mortality_biomarkers.csv}"

# ------------------------------------------------------------------------------
# Analysis code and output
# ------------------------------------------------------------------------------

PY_SCRIPT="${PY_SCRIPT:-/cbica/home/wenju/Project/whole-body_clocks/mortality_prediction_disease_clock/comparative_with_10_conventional_biomarkers/mortality_compare_mi_immune_metabolomics_vs_10_biomarkers.py}"

OUTPUT_DIR="${OUTPUT_DIR:-${BASE_DIR}/comparative_with_10_conventional_biomarkers/2_mortality_comparison_mi_immune_metabolomics}"

ADMIN_CENSOR_DATE="${ADMIN_CENSOR_DATE:-2022-11-30}"

PENALIZER="${PENALIZER:-0.01}"

MIN_EVENTS="${MIN_EVENTS:-20}"

# Use the environment that contains lifelines/scipy/pandas/openpyxl.
CONDA_ENV="${CONDA_ENV:-survival}"

# ------------------------------------------------------------------------------
# Target clock
# ------------------------------------------------------------------------------
#
# The Python script resolves the exact MI immune-metabolomics score column
# from the stable/significant 47-clock metadata. This avoids hard-coding a
# guessed column name.
#
# If you know the exact score column and want to force it, export:
#
#   SCORE_COL="exact_column_name"
#
# before submitting.
# ------------------------------------------------------------------------------

DISEASE_KEY="${DISEASE_KEY:-mi}"
ORGAN_KEY="${ORGAN_KEY:-immune}"
MODALITY_KEY="${MODALITY_KEY:-metabolomics}"
SCORE_COL="${SCORE_COL:-}"

# ------------------------------------------------------------------------------
# Optional exact ID/date/column overrides
# ------------------------------------------------------------------------------

FIELD53_0_COL="${FIELD53_0_COL:-}"
DEATH_DATE_COL="${DEATH_DATE_COL:-}"
DEATH_ID_COL="${DEATH_ID_COL:-}"
IDMATCH_SCORE_COL="${IDMATCH_SCORE_COL:-}"
IDMATCH_DEATH_COL="${IDMATCH_DEATH_COL:-}"
COVARIATE_ID_COL="${COVARIATE_ID_COL:-}"
BIOMARKER_ID_COL="${BIOMARKER_ID_COL:-}"

# ------------------------------------------------------------------------------
# Optional exact baseline covariate override
# ------------------------------------------------------------------------------
#
# Leave empty to use the recommended comparison baseline:
#   age + sex + ethnicity + assessment centre + smoking + BMI
#
# BP is omitted because systolic BP is one of the 10 comparison biomarkers.
#
# Example exact override:
# COVARIATE_COLS="age_when_attended_assessment_centre_f21003_0_0,sex_f31_0_0,..."
# ------------------------------------------------------------------------------

COVARIATE_COLS="${COVARIATE_COLS:-}"

# Set OVERWRITE=1 to rerun when the main summary already exists.
OVERWRITE="${OVERWRITE:-0}"

mkdir -p /cbica/home/wenju/output
mkdir -p "${OUTPUT_DIR}"

EXPECTED_SUMMARY="${OUTPUT_DIR}/mi_immune_metabolomics_EPOCH_vs_10_biomarkers_mortality_summary.tsv"

# ------------------------------------------------------------------------------
# Input checks
# ------------------------------------------------------------------------------

REQUIRED_FILES=(
  "${SCORE_WIDE_TSV}"
  "${SCORE_METADATA_TSV}"
  "${DEATH_XLSX}"
  "${ID_MATCH_CSV}"
  "${COVARIATE_CSV}"
  "${BIOMARKER_CSV}"
  "${PY_SCRIPT}"
)

for f in "${REQUIRED_FILES[@]}"; do
  if [[ ! -f "${f}" ]]; then
    echo "ERROR: required file not found:"
    echo "  ${f}"
    exit 1
  fi
done

if [[ -n "${GOOD_CLOCK_TSV}" && ! -f "${GOOD_CLOCK_TSV}" ]]; then
  echo "WARNING: GOOD_CLOCK_TSV not found; continuing without it:"
  echo "  ${GOOD_CLOCK_TSV}"
  GOOD_CLOCK_TSV=""
fi

if [[ -s "${EXPECTED_SUMMARY}" && "${OVERWRITE}" != "1" ]]; then
  echo "Output already exists:"
  echo "  ${EXPECTED_SUMMARY}"
  echo "Set OVERWRITE=1 to rerun."
  exit 0
fi

echo "================================================================================"
echo "MI immune-metabolomics EPOCH vs 10 conventional biomarkers"
echo "Apple-to-apple all-cause mortality landmark survival comparison"
echo "================================================================================"
echo "SCORE_WIDE_TSV:"
echo "  ${SCORE_WIDE_TSV}"
echo "SCORE_METADATA_TSV:"
echo "  ${SCORE_METADATA_TSV}"
echo "GOOD_CLOCK_TSV:"
echo "  ${GOOD_CLOCK_TSV:-<not used>}"
echo "BIOMARKER_CSV:"
echo "  ${BIOMARKER_CSV}"
echo "COVARIATE_CSV:"
echo "  ${COVARIATE_CSV}"
echo "DEATH_XLSX:"
echo "  ${DEATH_XLSX}"
echo "ID_MATCH_CSV:"
echo "  ${ID_MATCH_CSV}"
echo "OUTPUT_DIR:"
echo "  ${OUTPUT_DIR}"
echo "Landmark:"
echo "  UKB Field 53 instance 0_0"
echo "Administrative censor:"
echo "  ${ADMIN_CENSOR_DATE}"
echo "================================================================================"

# ------------------------------------------------------------------------------
# Environment
# ------------------------------------------------------------------------------

source activate "${CONDA_ENV}"

declare -a EXTRA_ARGS=()

if [[ -n "${GOOD_CLOCK_TSV}" ]]; then
  EXTRA_ARGS+=(--good-clock-tsv "${GOOD_CLOCK_TSV}")
fi

if [[ -n "${SCORE_COL}" ]]; then
  EXTRA_ARGS+=(--score-col "${SCORE_COL}")
fi

if [[ -n "${FIELD53_0_COL}" ]]; then
  EXTRA_ARGS+=(--field53-0-col "${FIELD53_0_COL}")
fi

if [[ -n "${DEATH_DATE_COL}" ]]; then
  EXTRA_ARGS+=(--death-date-col "${DEATH_DATE_COL}")
fi

if [[ -n "${DEATH_ID_COL}" ]]; then
  EXTRA_ARGS+=(--death-id-col "${DEATH_ID_COL}")
fi

if [[ -n "${IDMATCH_SCORE_COL}" ]]; then
  EXTRA_ARGS+=(--idmatch-score-col "${IDMATCH_SCORE_COL}")
fi

if [[ -n "${IDMATCH_DEATH_COL}" ]]; then
  EXTRA_ARGS+=(--idmatch-death-col "${IDMATCH_DEATH_COL}")
fi

if [[ -n "${COVARIATE_ID_COL}" ]]; then
  EXTRA_ARGS+=(--covariate-id-col "${COVARIATE_ID_COL}")
fi

if [[ -n "${BIOMARKER_ID_COL}" ]]; then
  EXTRA_ARGS+=(--biomarker-id-col "${BIOMARKER_ID_COL}")
fi

if [[ -n "${COVARIATE_COLS}" ]]; then
  EXTRA_ARGS+=(--covariate-cols "${COVARIATE_COLS}")
fi

python "${PY_SCRIPT}" \
  --score-wide-tsv "${SCORE_WIDE_TSV}" \
  --score-metadata-tsv "${SCORE_METADATA_TSV}" \
  --death-xlsx "${DEATH_XLSX}" \
  --id-match-csv "${ID_MATCH_CSV}" \
  --covariate-csv "${COVARIATE_CSV}" \
  --biomarker-csv "${BIOMARKER_CSV}" \
  --output-dir "${OUTPUT_DIR}" \
  --disease-key "${DISEASE_KEY}" \
  --organ-key "${ORGAN_KEY}" \
  --modality-key "${MODALITY_KEY}" \
  --admin-censor-date "${ADMIN_CENSOR_DATE}" \
  --penalizer "${PENALIZER}" \
  --min-events "${MIN_EVENTS}" \
  ${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"}

echo "================================================================================"
echo "Finished."
echo "Main comparison summary:"
echo "  ${EXPECTED_SUMMARY}"
echo "================================================================================"

conda deactivate
