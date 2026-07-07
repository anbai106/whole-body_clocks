#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=qmplot_mi
#SBATCH --array=0-11
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=1
#SBATCH --output=/cbica/home/wenju/output/qmplot_mi_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/qmplot_mi_%A_%a.err

# ============================================================
# Run qmplot for MI disease L'EPOCH fastGWA results
#
# One Slurm array task = one disease-clock fastGWA result.
#
# To generalize to another disease:
#   1. Change DISEASE below, e.g. DISEASE="dementia"
#   2. Change the #SBATCH job-name/output/error names manually
#   3. Change the #SBATCH --array range to match the number of clocks
#   4. Adjust CLOCK_PREFIXES if that disease has different available clocks
#
# For the MI list below:
#   N = 12 clocks
#   Therefore #SBATCH --array=0-11%10
# ============================================================

set -euo pipefail

# -----------------------------
# User settings
# -----------------------------

DISEASE="mi"

BASE_DIR="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock"

PROJECT_DIR="/cbica/home/wenju/Project/whole-body_clocks/${DISEASE}_clock/fastGWA"

PLOT_SCRIPT="${PROJECT_DIR}/5_qmplt.py"

SLURM_OUTPUT_DIR="/cbica/home/wenju/output"

CONDA_ENV="DNE"

mkdir -p "${SLURM_OUTPUT_DIR}"

# ============================================================
# Hard-coded clock list for MI
#
# Folder format:
#   <CLOCK_PREFIX>_<DISEASE>_clock
#
# Example:
#   Brain_proteomics_mi_clock
# ============================================================

CLOCK_PREFIXES=(
  "Brain_proteomics"
  "Digestive_metabolomics"
  "Endocrine_metabolomics"
  "Endocrine_proteomics"
  "Heart_proteomics"
  "Hepatic_metabolomics"
  "Hepatic_proteomics"
  "Immune_metabolomics"
  "Immune_proteomics"
  "Pulmonary_proteomics"
  "Reproductive_female_proteomics"
  "Reproductive_male_proteomics"
)

N_CLOCKS=${#CLOCK_PREFIXES[@]}

# ============================================================
# Basic checks
# ============================================================

echo "============================================================"
echo "Running qmplot Slurm task"
echo "============================================================"
echo "DISEASE: ${DISEASE}"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-NA}"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "Host: $(hostname)"
echo "Date: $(date)"
echo "Number of hard-coded clocks: ${N_CLOCKS}"
echo "============================================================"

if [[ ! -d "${BASE_DIR}" ]]; then
  echo "ERROR: BASE_DIR does not exist: ${BASE_DIR}" >&2
  exit 1
fi

if [[ ! -f "${PLOT_SCRIPT}" ]]; then
  echo "ERROR: PLOT_SCRIPT does not exist: ${PLOT_SCRIPT}" >&2
  exit 1
fi

if [[ "${SLURM_ARRAY_TASK_ID}" -ge "${N_CLOCKS}" ]]; then
  echo "ERROR: SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} but only ${N_CLOCKS} clocks are hard-coded." >&2
  echo "If you changed the clock list, update #SBATCH --array accordingly." >&2
  exit 1
fi

# ============================================================
# Select clock for this array task
# ============================================================

CLOCK_PREFIX="${CLOCK_PREFIXES[${SLURM_ARRAY_TASK_ID}]}"
CLOCK_NAME="${CLOCK_PREFIX}_${DISEASE}_clock"
CLOCK_DIR="${BASE_DIR}/${CLOCK_NAME}"

OUTPUT_RESULT="${CLOCK_DIR}/fastGWA/output/organ_pheno_normalized_residualized.fastGWA"
OUTPUT_DIR="${CLOCK_DIR}/fastGWA/qmplot"

echo "Selected clock:"
echo "  CLOCK_PREFIX:  ${CLOCK_PREFIX}"
echo "  CLOCK_NAME:    ${CLOCK_NAME}"
echo "  CLOCK_DIR:     ${CLOCK_DIR}"
echo "  OUTPUT_RESULT: ${OUTPUT_RESULT}"
echo "  OUTPUT_DIR:    ${OUTPUT_DIR}"

if [[ ! -d "${CLOCK_DIR}" ]]; then
  echo "ERROR: Clock directory does not exist: ${CLOCK_DIR}" >&2
  exit 1
fi

if [[ ! -f "${OUTPUT_RESULT}" ]]; then
  echo "ERROR: fastGWA result does not exist: ${OUTPUT_RESULT}" >&2
  exit 1
fi

mkdir -p "${OUTPUT_DIR}"

# ============================================================
# Run qmplot
# ============================================================

module load python/anaconda/3

source activate "${CONDA_ENV}"

echo "============================================================"
echo "Start qmplot"
echo "============================================================"

python -u "${PLOT_SCRIPT}" \
  --output_dir "${OUTPUT_DIR}" \
  --output_result "${OUTPUT_RESULT}" \
  --clock_name "${CLOCK_NAME}"

echo "============================================================"
echo "Finish qmplot"
echo "Date: $(date)"
echo "============================================================"

conda deactivate || true