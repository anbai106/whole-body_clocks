#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=qmplot_copd
#SBATCH --array=0-10%10
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=1
#SBATCH --output=/cbica/home/wenju/output/qmplot_copd_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/qmplot_copd_%A_%a.err

# ============================================================
# Run qmplot for copd disease L'EPOCH fastGWA results
#
# One Slurm array task = one disease-clock fastGWA result.
#
# Input pattern:
#   /cbica/home/wenju/Reproducibile_paper/WholeBodyClock/*_copd_clock/fastGWA/output/*fastGWA
#
# Output:
#   /cbica/home/wenju/Reproducibile_paper/WholeBodyClock/<clock_folder>/fastGWA/qmplot/
#
# Recommended:
#   bash submit_qmplot_copd.slurm
#
# Or directly:
#   sbatch submit_qmplot_copd.slurm
# ============================================================

set -euo pipefail

# -----------------------------
# User settings
# -----------------------------

DISEASE="${DISEASE:-copd}"

BASE_DIR="${BASE_DIR:-/cbica/home/wenju/Reproducibile_paper/WholeBodyClock}"

PROJECT_DIR="${PROJECT_DIR:-/cbica/home/wenju/Project/whole-body_clocks/copd_clock/fastGWA}"

PLOT_SCRIPT="${PLOT_SCRIPT:-${PROJECT_DIR}/5_qmplt.py}"

TASK_FILE="${TASK_FILE:-${PROJECT_DIR}/qmplot_${DISEASE}_fastgwa_tasks.tsv}"

SLURM_OUTPUT_DIR="${SLURM_OUTPUT_DIR:-/cbica/home/wenju/output}"

CONDA_ENV="${CONDA_ENV:-DNE}"

MAX_PARALLEL="${MAX_PARALLEL:-10}"

mkdir -p "${SLURM_OUTPUT_DIR}"
mkdir -p "$(dirname "${TASK_FILE}")"

# ============================================================
# Function: create task file
# ============================================================

make_task_file() {
  echo "Creating qmplot task file:"
  echo "  ${TASK_FILE}"

  if [[ ! -d "${BASE_DIR}" ]]; then
    echo "ERROR: BASE_DIR does not exist: ${BASE_DIR}" >&2
    exit 1
  fi

  if [[ ! -f "${PLOT_SCRIPT}" ]]; then
    echo "ERROR: PLOT_SCRIPT does not exist: ${PLOT_SCRIPT}" >&2
    exit 1
  fi

  shopt -s nullglob

  fastgwa_files=(
    "${BASE_DIR}"/*_"${DISEASE}"_clock/fastGWA/output/*fastGWA
  )

  shopt -u nullglob

  if [[ "${#fastgwa_files[@]}" -eq 0 ]]; then
    echo "ERROR: No fastGWA files found for disease=${DISEASE}" >&2
    echo "Expected pattern:" >&2
    echo "  ${BASE_DIR}/*_${DISEASE}_clock/fastGWA/output/*fastGWA" >&2
    exit 1
  fi

  {
    printf "array_id\tclock_name\tclock_dir\toutput_result\toutput_dir\n"

    i=0
    for output_result in "${fastgwa_files[@]}"; do
      output_result="$(readlink -f "${output_result}")"

      fastgwa_output_dir="$(dirname "${output_result}")"
      fastgwa_dir="$(dirname "${fastgwa_output_dir}")"
      clock_dir="$(dirname "${fastgwa_dir}")"
      clock_name="$(basename "${clock_dir}")"

      output_dir="${clock_dir}/fastGWA/qmplot"
      mkdir -p "${output_dir}"

      printf "%s\t%s\t%s\t%s\t%s\n" \
        "${i}" \
        "${clock_name}" \
        "${clock_dir}" \
        "${output_result}" \
        "${output_dir}"

      i=$((i + 1))
    done
  } > "${TASK_FILE}"

  echo "Task file created:"
  echo "  ${TASK_FILE}"
  echo "Number of tasks:"
  tail -n +2 "${TASK_FILE}" | wc -l

  echo "Preview:"
  column -t -s $'\t' "${TASK_FILE}" | head -n 20
}

# ============================================================
# Self-submit mode
# ============================================================

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "============================================================"
  echo "Preparing qmplot Slurm array for disease=${DISEASE}"
  echo "============================================================"
  echo "BASE_DIR:          ${BASE_DIR}"
  echo "PROJECT_DIR:       ${PROJECT_DIR}"
  echo "PLOT_SCRIPT:       ${PLOT_SCRIPT}"
  echo "TASK_FILE:         ${TASK_FILE}"
  echo "SLURM_OUTPUT_DIR:  ${SLURM_OUTPUT_DIR}"
  echo "CONDA_ENV:         ${CONDA_ENV}"
  echo "MAX_PARALLEL:      ${MAX_PARALLEL}"
  echo "============================================================"

  make_task_file

  N_TASKS=$(($(wc -l < "${TASK_FILE}") - 1))

  if [[ "${N_TASKS}" -lt 1 ]]; then
    echo "ERROR: No tasks created." >&2
    exit 1
  fi

  echo "Submitting ${N_TASKS} qmplot jobs with max parallel=${MAX_PARALLEL}"

  sbatch \
    --array=0-$((N_TASKS - 1))%${MAX_PARALLEL} \
    "$0"

  exit 0
fi

# ============================================================
# Slurm task mode
# ============================================================

echo "============================================================"
echo "Running qmplot Slurm task"
echo "============================================================"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "Host: $(hostname)"
echo "Date: $(date)"
echo "============================================================"

# If submitted directly with sbatch and task file is missing,
# task 0 creates it; other tasks wait.
if [[ ! -f "${TASK_FILE}" ]]; then
  if [[ "${SLURM_ARRAY_TASK_ID}" == "0" ]]; then
    echo "Task file missing. Task 0 will create it."
    make_task_file
  else
    echo "Task file missing. Waiting for task 0 to create it..."
    for i in $(seq 1 120); do
      if [[ -f "${TASK_FILE}" ]]; then
        break
      fi
      sleep 5
    done
  fi
fi

if [[ ! -f "${TASK_FILE}" ]]; then
  echo "ERROR: Task file still missing after waiting: ${TASK_FILE}" >&2
  exit 1
fi

mapfile -t TASK_LINES < <(tail -n +2 "${TASK_FILE}")

N_TASKS="${#TASK_LINES[@]}"

if [[ "${SLURM_ARRAY_TASK_ID}" -ge "${N_TASKS}" ]]; then
  echo "ERROR: SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} but only ${N_TASKS} tasks exist." >&2
  exit 1
fi

TASK_LINE="${TASK_LINES[${SLURM_ARRAY_TASK_ID}]}"

IFS=$'\t' read -r array_id clock_name clock_dir output_result output_dir <<< "${TASK_LINE}"

echo "Selected task:"
echo "  array_id:      ${array_id}"
echo "  clock_name:    ${clock_name}"
echo "  clock_dir:     ${clock_dir}"
echo "  output_result: ${output_result}"
echo "  output_dir:    ${output_dir}"

if [[ ! -f "${output_result}" ]]; then
  echo "ERROR: fastGWA result does not exist: ${output_result}" >&2
  exit 1
fi

mkdir -p "${output_dir}"

module load python/anaconda/3

# CUBIC-style conda activation
source activate "${CONDA_ENV}"

echo "============================================================"
echo "Start qmplot"
echo "============================================================"

python -u "${PLOT_SCRIPT}" \
  --output_dir "${output_dir}" \
  --output_result "${output_result}" \
  --clock_name "${clock_name}"

echo "============================================================"
echo "Finish qmplot"
echo "Date: $(date)"
echo "============================================================"

conda deactivate || true