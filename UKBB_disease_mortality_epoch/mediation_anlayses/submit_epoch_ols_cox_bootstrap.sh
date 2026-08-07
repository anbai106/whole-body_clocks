#!/bin/bash
#SBATCH --partition=bioinformatics
#SBATCH --job-name=EPOCH_MED_1MODEL
#SBATCH --array=0-314%315
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=0-06:00:00
#SBATCH --output=/cbica/home/wenju/output/EPOCH_MED_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/EPOCH_MED_%A_%a.err

# ==============================================================================
# 315-way SLURM array: exactly ONE mediation model per array task
#
# Current model family:
#   45 molecular disease EPOCH exposures x 7 MRI mortality EPOCH mediators
#   = 315 independent mediation jobs.
#
# Each array task runs exactly one model:
#   SLURM_ARRAY_TASK_ID 0   -> model index 0
#   ...
#   SLURM_ARRAY_TASK_ID 314 -> model index 314
#
# IMPORTANT:
#   This version intentionally does NOT use:
#       set -euo pipefail
#   and it does NOT use empty Bash argument arrays such as SUBSET_ARGS=().
#   This avoids the "SUBSET_ARGS[@]: unbound variable" failure seen on CUBIC.
#
# No aggregation is performed here. Each task writes one TSV under:
#   <OUTPUT>/single_model_results/
#
# Default analysis: full population
# Test-only analysis:
#   ANALYSIS_MODE=test sbatch submit_epoch_ols_cox_bootstrap_315_single_no_strict.sh
#
# Optional environment overrides:
#   BOOTSTRAP=1000
#   MINIMUM_N=500
#   MINIMUM_DEATHS=20
#   COX_PENALIZER=0.0
#   SAVE_BOOTSTRAP=0
# ==============================================================================

SCRIPT=/cbica/home/wenju/Project/whole-body_clocks/UKBB_disease_mortality_epoch/mediation_anlayses/ols_cox_bootstrap_mediation.py
ROOT_OUTPUT=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock
EXPECTED_MODELS=315

BOOTSTRAP=${BOOTSTRAP:-1000}
MINIMUM_N=${MINIMUM_N:-500}
MINIMUM_DEATHS=${MINIMUM_DEATHS:-20}
COX_PENALIZER=${COX_PENALIZER:-0.0}
ANALYSIS_MODE=${ANALYSIS_MODE:-full}
SAVE_BOOTSTRAP=${SAVE_BOOTSTRAP:-0}

mkdir -p /cbica/home/wenju/output

if [[ ! -f "$SCRIPT" ]]; then
    echo "ERROR: Python script not found: $SCRIPT" >&2
    exit 1
fi

if [[ -z "$SLURM_ARRAY_TASK_ID" ]]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID is empty. Submit this script with sbatch as an array job." >&2
    exit 1
fi

if [[ "$SLURM_ARRAY_TASK_ID" -lt 0 || "$SLURM_ARRAY_TASK_ID" -ge "$EXPECTED_MODELS" ]]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID is outside 0-$((EXPECTED_MODELS - 1))." >&2
    exit 1
fi

# Prevent BLAS/numerical libraries from starting extra threads in each array task.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

# Activate the survival-analysis environment.
if command -v conda >/dev/null 2>&1; then
    CONDA_BASE=$(conda info --base 2>/dev/null)
    if [[ -n "$CONDA_BASE" && -f "$CONDA_BASE/etc/profile.d/conda.sh" ]]; then
        source "$CONDA_BASE/etc/profile.d/conda.sh"
        conda activate survival_clock
    else
        source activate survival_clock
    fi
else
    source activate survival_clock
fi

# Confirm that the required Python packages are available.
python - <<'PY'
import importlib.util
required = ["numpy", "pandas", "scipy", "statsmodels", "lifelines"]
missing = [p for p in required if importlib.util.find_spec(p) is None]
if missing:
    raise SystemExit("Missing required Python packages: " + ", ".join(missing))
print("Python dependency check passed.")
PY

if [[ $? -ne 0 ]]; then
    echo "ERROR: Python dependency check failed." >&2
    exit 1
fi

# Choose the analysis subset and output location.
# No Bash arrays are used here.
if [[ "$ANALYSIS_MODE" == "full" ]]; then
    OUTPUT="$ROOT_OUTPUT/mediation_OLS_Cox_bootstrap_full_single_models/job_${SLURM_ARRAY_JOB_ID}"
elif [[ "$ANALYSIS_MODE" == "test" ]]; then
    OUTPUT="$ROOT_OUTPUT/mediation_OLS_Cox_bootstrap_test_single_models/job_${SLURM_ARRAY_JOB_ID}"
else
    echo "ERROR: ANALYSIS_MODE must be 'full' or 'test'; got '$ANALYSIS_MODE'." >&2
    exit 1
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
echo "Save bootstrap samples: ${SAVE_BOOTSTRAP}"
echo "Output directory: ${OUTPUT}"
echo "============================================================"

# Run exactly one mediation model. We use explicit branches rather than optional
# Bash arrays so an empty optional argument can never trigger an unbound-variable
# error.
if [[ "$ANALYSIS_MODE" == "test" ]]; then
    if [[ "$SAVE_BOOTSTRAP" == "1" ]]; then
        python -u "$SCRIPT" \
            --test-only \
            --save-bootstrap-samples \
            --model-index "$MODEL_INDEX" \
            --expected-model-count "$EXPECTED_MODELS" \
            --output-dir "$OUTPUT" \
            --bootstrap "$BOOTSTRAP" \
            --minimum-n "$MINIMUM_N" \
            --minimum-deaths "$MINIMUM_DEATHS" \
            --cox-penalizer "$COX_PENALIZER"
    else
        python -u "$SCRIPT" \
            --test-only \
            --model-index "$MODEL_INDEX" \
            --expected-model-count "$EXPECTED_MODELS" \
            --output-dir "$OUTPUT" \
            --bootstrap "$BOOTSTRAP" \
            --minimum-n "$MINIMUM_N" \
            --minimum-deaths "$MINIMUM_DEATHS" \
            --cox-penalizer "$COX_PENALIZER"
    fi
else
    if [[ "$SAVE_BOOTSTRAP" == "1" ]]; then
        python -u "$SCRIPT" \
            --save-bootstrap-samples \
            --model-index "$MODEL_INDEX" \
            --expected-model-count "$EXPECTED_MODELS" \
            --output-dir "$OUTPUT" \
            --bootstrap "$BOOTSTRAP" \
            --minimum-n "$MINIMUM_N" \
            --minimum-deaths "$MINIMUM_DEATHS" \
            --cox-penalizer "$COX_PENALIZER"
    else
        python -u "$SCRIPT" \
            --model-index "$MODEL_INDEX" \
            --expected-model-count "$EXPECTED_MODELS" \
            --output-dir "$OUTPUT" \
            --bootstrap "$BOOTSTRAP" \
            --minimum-n "$MINIMUM_N" \
            --minimum-deaths "$MINIMUM_DEATHS" \
            --cox-penalizer "$COX_PENALIZER"
    fi
fi

PYTHON_EXIT=$?

if [[ "$PYTHON_EXIT" -ne 0 ]]; then
    echo "ERROR: mediation model ${MODEL_NUMBER}/${EXPECTED_MODELS} failed with exit code ${PYTHON_EXIT}." >&2
    exit "$PYTHON_EXIT"
fi

echo "============================================================"
echo "Finish: $(date)"
echo "Completed model ${MODEL_NUMBER}/${EXPECTED_MODELS}"
echo "============================================================"