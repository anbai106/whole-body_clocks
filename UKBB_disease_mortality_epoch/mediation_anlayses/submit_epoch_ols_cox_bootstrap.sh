#!/bin/bash
#SBATCH --partition=bioinformatics
#SBATCH --job-name=EPOCH_OLS_COX_MED
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=0-06:00:00
#SBATCH --output=/cbica/home/wenju/output/EPOCH_OLS_COX_MED_%j.out
#SBATCH --error=/cbica/home/wenju/output/EPOCH_OLS_COX_MED_%j.err

SCRIPT=/cbica/home/wenju/Project/whole-body_clocks/UKBB_disease_mortality_epoch/mediation_anlayses/ols_cox_bootstrap_mediation.py

if [[ ! -f "$SCRIPT" ]]; then
    echo "ERROR: Python script not found: $SCRIPT" >&2
    exit 1
fi

# Full-population analysis is the DEFAULT. No --test-only flag is passed.
OUTPUT=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mediation_OLS_Cox_bootstrap_full

# Number of participant-level bootstrap replicates per successfully fitted model.
# 1000 is a practical starting point. Bonferroni inference itself uses the
# continuous delta-method indirect p-value; bootstrap is used for percentile CIs.
BOOTSTRAP=1000

# Optional: if you use a conda/venv instead of a module, activate it here.
# source /path/to/venv/bin/activate
source activate survival_clock

python - <<'PY'
import importlib.util
required = ["numpy", "pandas", "scipy", "statsmodels", "lifelines"]
missing = [p for p in required if importlib.util.find_spec(p) is None]
if missing:
    raise SystemExit("Missing required Python packages: " + ", ".join(missing))
print("Python dependency check passed.")
PY

echo "Start: $(date)"
echo "Host: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID:-NA}"
echo "Analysis: full population"
echo "Bootstrap replicates/model: ${BOOTSTRAP}"

python -u "$SCRIPT" \
    --output-dir "$OUTPUT" \
    --bootstrap "$BOOTSTRAP" \
    --minimum-n 500 \
    --minimum-deaths 20 \
    --n-pcs 10 \
    --cox-penalizer 0.0

echo "Finish: $(date)"

# -----------------------------------------------------------------------------
# TEST-ONLY SENSITIVITY
# -----------------------------------------------------------------------------
# To run only the held-out mortality-clock test participants, either copy this
# SLURM file and replace the Python command above with the block below, or run it
# interactively. Full population remains the default behavior.
#
# TEST_OUTPUT=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mediation_OLS_Cox_bootstrap_test
# python -u "$SCRIPT" \
#     --test-only \
#     --output-dir "$TEST_OUTPUT" \
#     --bootstrap "$BOOTSTRAP" \
#     --minimum-n 500 \
#     --minimum-deaths 20 \
#     --n-pcs 10 \
#     --cox-penalizer 0.0
