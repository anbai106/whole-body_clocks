#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=EPOCH23_PH
#SBATCH --mem-per-cpu=16G
#SBATCH --time=08:00:00
#SBATCH --output=/cbica/home/wenju/output/EPOCH23_PH_%j.out
#SBATCH --error=/cbica/home/wenju/output/EPOCH23_PH_%j.err

# ==============================================================================
# Post-hoc proportional-hazards diagnostics for all 23 existing mortality EPOCH
# clocks.
#
# NO EPOCH model is retrained or retuned.
#
# Reads each existing:
#   *_mortality_clock_predictions.tsv
#
# Primary analysis:
#   held-out TEST split
#
# Default Cox model:
#   Surv(time_years, event) ~ EPOCH acceleration z + age + sex
#
# Age and sex are REQUIRED in the default age_sex analysis.
#
# IMPORTANT:
# "Reproductive female proteomics" and "Reproductive male proteomics" are clock
# labels / proteomic feature definitions, NOT female-only or male-only samples.
# The script NEVER filters participants based on clock name. Both males and
# females are retained and sex is included as a covariate for these clocks
# exactly as for the other mortality EPOCHs.
#
# Metrics:
#   - EPOCH HR + 95% CI + P
#   - EPOCH-specific scaled-Schoenfeld PH test
#   - BH-FDR across 23 EPOCH PH tests
#   - covariate-specific PH tests
#   - descriptive Schoenfeld residual/time rho
#   - Schoenfeld residual plots
#   - piecewise time-varying sensitivity at 5 years:
#         HR 0-5 years
#         HR >5 years
#         EPOCH x post-5y interaction P
#         1-df LRT P when both time-varying models fit unpenalized
#
# Master outputs:
#   mortality_EPOCH_23_PH_diagnostics.tsv
#   mortality_EPOCH_23_PH_covariate_tests.tsv
#   mortality_EPOCH_23_PH_manuscript_table.tsv
#   mortality_EPOCH_23_PH_run_manifest.tsv
#
# Dependency:
#   lifelines must be available in the survival_clock environment.
# ==============================================================================

BASE_DIR="${BASE_DIR:-/cbica/home/wenju/Reproducibile_paper/WholeBodyClock}"

PY_SCRIPT="${PY_SCRIPT:-/cbica/home/wenju/Project/whole-body_clocks/mortality_clock/PH_assumption_metrics/compute_23_mortality_EPOCH_PH_diagnostics.py}"

OUTPUT_DIR="${OUTPUT_DIR:-${BASE_DIR}/mortality_EPOCH_PH_diagnostics}"

CONDA_ENV="${CONDA_ENV:-survival_clock}"

TEST_SPLIT="${TEST_SPLIT:-test}"

# Recommended default for manuscript Cox HR diagnostics.
# Alternative: risk_score
PREDICTOR_MODE="${PREDICTOR_MODE:-acceleration}"

# Recommended and default: age_sex
# This mode requires age plus BOTH female and male participants in every clock.
# Alternative: none (unadjusted sensitivity analysis only).
COVARIATE_MODE="${COVARIATE_MODE:-age_sex}"

# "km" is closest to the usual cox.zph KM time transform.
TIME_TRANSFORM="${TIME_TRANSFORM:-km}"

# Prespecified time-varying sensitivity split.
TIME_SPLIT_YEARS="${TIME_SPLIT_YEARS:-5}"

MIN_EVENTS_PER_TIME_BAND="${MIN_EVENTS_PER_TIME_BAND:-20}"

mkdir -p /cbica/home/wenju/output
mkdir -p "${OUTPUT_DIR}"

if [[ ! -d "${BASE_DIR}" ]]; then
  echo "ERROR: BASE_DIR does not exist:"
  echo "  ${BASE_DIR}"
  exit 1
fi

if [[ ! -f "${PY_SCRIPT}" ]]; then
  echo "ERROR: Python script does not exist:"
  echo "  ${PY_SCRIPT}"
  exit 1
fi

echo "================================================================================"
echo "Post-hoc PH diagnostics for 23 mortality EPOCH clocks"
echo "================================================================================"
echo "BASE_DIR:"
echo "  ${BASE_DIR}"
echo "PY_SCRIPT:"
echo "  ${PY_SCRIPT}"
echo "OUTPUT_DIR:"
echo "  ${OUTPUT_DIR}"
echo "TEST_SPLIT:"
echo "  ${TEST_SPLIT}"
echo "PREDICTOR_MODE:"
echo "  ${PREDICTOR_MODE}"
echo "COVARIATE_MODE:"
echo "  ${COVARIATE_MODE}"
echo "TIME_TRANSFORM:"
echo "  ${TIME_TRANSFORM}"
echo "TIME_SPLIT_YEARS:"
echo "  ${TIME_SPLIT_YEARS}"
echo "MIN_EVENTS_PER_TIME_BAND:"
echo "  ${MIN_EVENTS_PER_TIME_BAND}"
echo "================================================================================"

source activate "${CONDA_ENV}"

python "${PY_SCRIPT}" \
  --base-dir "${BASE_DIR}" \
  --output-dir "${OUTPUT_DIR}" \
  --test-split "${TEST_SPLIT}" \
  --predictor-mode "${PREDICTOR_MODE}" \
  --covariate-mode "${COVARIATE_MODE}" \
  --time-transform "${TIME_TRANSFORM}" \
  --time-split-years "${TIME_SPLIT_YEARS}" \
  --min-events-per-time-band "${MIN_EVENTS_PER_TIME_BAND}"

conda deactivate

echo "================================================================================"
echo "Finished 23-clock PH diagnostic analysis."
echo "Master diagnostics:"
echo "  ${OUTPUT_DIR}/mortality_EPOCH_23_PH_diagnostics.tsv"
echo "Manuscript-oriented summary:"
echo "  ${OUTPUT_DIR}/mortality_EPOCH_23_PH_manuscript_table.tsv"
echo "Covariate-specific PH tests:"
echo "  ${OUTPUT_DIR}/mortality_EPOCH_23_PH_covariate_tests.tsv"
echo "Run manifest:"
echo "  ${OUTPUT_DIR}/mortality_EPOCH_23_PH_run_manifest.tsv"
echo "================================================================================"
