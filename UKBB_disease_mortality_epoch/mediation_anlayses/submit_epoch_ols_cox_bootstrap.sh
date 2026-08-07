#!/bin/bash
#SBATCH --partition=bioinformatics
#SBATCH --job-name=EPOCH_MED_1MODEL
#SBATCH --array=0-314%315
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=0-06:00:00
#SBATCH --output=/cbica/home/wenju/output/EPOCH_MED_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/EPOCH_MED_%A_%a.err

set -euo pipefail

# ==============================================================================
# 315-way SLURM array: exactly ONE mediation model per array task
#
# Current model family:
#   45 molecular disease EPOCH exposures x 7 MRI mortality EPOCH mediators
#   = 315 independent mediation jobs.
#
# Each array task runs one model only:
#   SLURM_ARRAY_TASK_ID 0   -> model index 0
#   SLURM_ARRAY_TASK_ID 1   -> model index 1
#   ...
#   SLURM_ARRAY_TASK_ID 314 -> model index 314
#
# No aggregation is performed here. Each job writes one result TSV under:
#   <OUTPUT>/single_model_results/
# A separate collection script can be run after the array finishes.
#
# Default: full population.
# Test-only sensitivity:
#   ANALYSIS_MODE=test sbatch submit_epoch_ols_cox_bootstrap_315_single.sh
#
# NOTE: %315 permits up to 315 tasks to run concurrently. The scheduler/account
# limits may reduce actual concurrency automatically.
# ==============================================================================

SCRIPT=/cbica/home/wenju/Project/whole-body_clocks/UKBB_disease_mortality_epoch/mediation_anlayses/ols_cox_bootstrap_mediation.py
ROOT_OUTPUT=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock
EXPECTED_MODELS=315
BOOTSTRAP=${BOOTSTRAP:-1000}
MINIMUM_N=${MINIMUM_N:-500}
MINIMUM_DEATHS=${MINIMUM_DEATHS:-20}
N_PCS=${N_PCS:-10}
COX_PENALIZER=${COX_PENALIZER:-0.0}
ANALYSIS_MODE=${ANALYSIS_MODE:-full}
SAVE_BOOTSTRAP=${SAVE_BOOTSTRAP:-0}

if [[ ! -f "$SCRIPT" ]]; then
    echo "ERROR: Python script not found: $SCRIPT" >&2
    exit 1
fi

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: This launcher must be submitted as a SLURM array job." >&2
    exit 1
fi

if (( SLURM_ARRAY_TASK_ID < 0 || SLURM_ARRAY_TASK_ID >= EXPECTED_MODELS )); then
    echo "ERROR: SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} is outside 0-$((EXPECTED_MODELS - 1))." >&2
    exit 1
fi

mkdir -p /cbica/home/wenju/output

# Prevent numerical libraries from spawning extra threads inside each of the 315 jobs.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

# Activate survival-analysis environment.
if command -v conda >/dev/null 2>&1; then
    CONDA_BASE=$(conda info --base)
    # shellcheck disable=SC1091
    source "$CONDA_BASE/etc/profile.d/conda.sh"
    conda activate survival_clock
else
    # Compatibility with environments where legacy `source activate` is configured.
    set +u
    source activate survival_clock
    set -u
fi

# Lightweight dependency check for this task.
python - <<'PY'
import importlib.util
required = ["numpy", "pandas", "scipy", "statsmodels", "lifelines"]
missing = [p for p in required if importlib.util.find_spec(p) is None]
if missing:
    raise SystemExit("Missing required Python packages: " + ", ".join(missing))
PY

if [[ "$ANALYSIS_MODE" == "full" ]]; then
    SUBSET_ARGS=()
    OUTPUT="$ROOT_OUTPUT/mediation_OLS_Cox_bootstrap_full_single_models/job_${SLURM_ARRAY_JOB_ID}"
elif [[ "$ANALYSIS_MODE" == "test" ]]; then
    SUBSET_ARGS=(--test-only)
    OUTPUT="$ROOT_OUTPUT/mediation_OLS_Cox_bootstrap_test_single_models/job_${SLURM_ARRAY_JOB_ID}"
else
    echo "ERROR: ANALYSIS_MODE must be 'full' or 'test'; got '$ANALYSIS_MODE'." >&2
    exit 1
fi

SAVE_ARGS=()
if [[ "$SAVE_BOOTSTRAP" == "1" ]]; then
    SAVE_ARGS=(--save-bootstrap-samples)
fi

mkdir -p "$OUTPUT/single_model_results"

MODEL_INDEX="$SLURM_ARRAY_TASK_ID"
MODEL_NUMBER=$((MODEL_INDEX + 1))

echo "============================================================"
echo "Start: $(date)"
echo "Host: $(hostname)"
echo "Array job ID: ${SLURM_ARRAY_JOB_ID}"
echo "Array task ID / model index: ${MODEL_INDEX}"
echo "Model number: ${MODEL_NUMBER}/${EXPECTED_MODELS}"
echo "Analysis subset: ${ANALYSIS_MODE}"
echo "Bootstrap replicates: ${BOOTSTRAP}"
echo "Output directory: ${OUTPUT}"
echo "============================================================"

python -u "$SCRIPT" \
    "${SUBSET_ARGS[@]}" \
    "${SAVE_ARGS[@]}" \
    --model-index "$MODEL_INDEX" \
    --expected-model-count "$EXPECTED_MODELS" \
    --output-dir "$OUTPUT" \
    --bootstrap "$BOOTSTRAP" \
    --minimum-n "$MINIMUM_N" \
    --minimum-deaths "$MINIMUM_DEATHS" \
    --n-pcs "$N_PCS" \
    --cox-penalizer "$COX_PENALIZER"

echo "============================================================"
echo "Finish: $(date)"
echo "Completed model ${MODEL_NUMBER}/${EXPECTED_MODELS}"
echo "============================================================"