#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=EPOCH23_calibration
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --output=/cbica/home/wenju/output/EPOCH23_calibration_%j.out
#SBATCH --error=/cbica/home/wenju/output/EPOCH23_calibration_%j.err

# ==============================================================================
# Post-hoc calibration analysis for all 23 existing mortality EPOCH clocks.
#
# NO mortality EPOCH model is retrained.
#
# The Python script reads each existing:
#   *_mortality_clock_predictions.tsv
#
# and evaluates absolute-risk calibration in the held-out TEST split.
#
# Metrics:
#   - Cox calibration slope
#   - mean predicted mortality risk
#   - Kaplan-Meier observed mortality risk + 95% CI
#   - observed / expected risk ratio
#   - predicted - observed risk
#   - absolute calibration-in-the-large error
#   - IPCW Brier score
#   - null-KM Brier score and Brier skill score
#   - decile calibration error summaries
#   - held-out calibration plots
#
# Default horizon rule:
#   auto = use every risk_<time>y column that already exists in each clock's
#          saved prediction TSV. Thus this does NOT invent or reconstruct an
#          unsupported horizon.
#
# Master outputs:
#   mortality_EPOCH_23_calibration_by_horizon.tsv
#   mortality_EPOCH_23_calibration_slope.tsv
#   mortality_EPOCH_23_calibration_5y.tsv
#   mortality_EPOCH_23_calibration_manuscript_table.tsv
#   mortality_EPOCH_23_calibration_run_manifest.tsv
# ==============================================================================

set -euo pipefail

BASE_DIR="${BASE_DIR:-/cbica/home/wenju/Reproducibile_paper/WholeBodyClock}"

PY_SCRIPT="${PY_SCRIPT:-/cbica/home/wenju/Project/whole-body_clocks/mortality_clock/calibration_metrics/compute_23_mortality_EPOCH_calibration_metrics.py}"

OUTPUT_DIR="${OUTPUT_DIR:-${BASE_DIR}/mortality_EPOCH_calibration_metrics}"

CONDA_ENV="${CONDA_ENV:-survival_clock}"

# "auto" is recommended. It evaluates risk_5y/risk_10y/risk_15y only when the
# corresponding absolute-risk column already exists in a given prediction TSV.
HORIZONS="${HORIZONS:-auto}"

N_BINS="${N_BINS:-10}"

MINIMUM_BIN_N="${MINIMUM_BIN_N:-20}"

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
echo "Post-hoc calibration of 23 mortality EPOCH clocks"
echo "================================================================================"
echo "BASE_DIR:"
echo "  ${BASE_DIR}"
echo "PY_SCRIPT:"
echo "  ${PY_SCRIPT}"
echo "OUTPUT_DIR:"
echo "  ${OUTPUT_DIR}"
echo "HORIZONS:"
echo "  ${HORIZONS}"
echo "N_BINS:"
echo "  ${N_BINS}"
echo "MINIMUM_BIN_N:"
echo "  ${MINIMUM_BIN_N}"
echo "================================================================================"

source activate "${CONDA_ENV}"

python "${PY_SCRIPT}" \
  --base-dir "${BASE_DIR}" \
  --output-dir "${OUTPUT_DIR}" \
  --horizons "${HORIZONS}" \
  --n-bins "${N_BINS}" \
  --minimum-bin-n "${MINIMUM_BIN_N}"

conda deactivate

echo "================================================================================"
echo "Finished 23-clock calibration analysis."
echo "Master manuscript table:"
echo "  ${OUTPUT_DIR}/mortality_EPOCH_23_calibration_manuscript_table.tsv"
echo "Common 5-year summary:"
echo "  ${OUTPUT_DIR}/mortality_EPOCH_23_calibration_5y.tsv"
echo "Run manifest:"
echo "  ${OUTPUT_DIR}/mortality_EPOCH_23_calibration_run_manifest.tsv"
echo "================================================================================"
