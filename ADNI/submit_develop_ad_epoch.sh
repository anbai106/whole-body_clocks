#!/usr/bin/env bash
#SBATCH --job-name=adni_brain_ad_lepoch
#SBATCH --output=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch/logs/adni_brain_ad_lepoch_%j.out
#SBATCH --error=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch/logs/adni_brain_ad_lepoch_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

set -euo pipefail

# ============================================================
# Build ADNI brain MRI AD L'EPOCH
# Baseline CN MUSE GM ROI features -> time to MCI/AD conversion
#
# Important:
#   Time zero is defined as the first CN visit with usable MUSE GM
#   ROI features. Prefer Visit_Code == bl if it has MUSE, otherwise
#   use the earliest MUSE-available CN visit.
# ============================================================

WORKDIR="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch"
LOGDIR="${WORKDIR}/logs"
OUTDIR="${WORKDIR}/results_brain_mri_ad_lepoch"

mkdir -p "${WORKDIR}"
mkdir -p "${LOGDIR}"
mkdir -p "${OUTDIR}"

PY_SCRIPT="/cbica/home/wenju/Project/whole-body_clocks/ADNI/ad_epoch.py"

ADNI_FILE="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch/adni_istaging.tsv"

# ------------------------------------------------------------
# Python environment
# ------------------------------------------------------------

source activate survival_clock

python3 - <<'PY'
import importlib.util
missing = []
for pkg in ["pandas", "numpy", "sklearn", "sksurv", "joblib"]:
    if importlib.util.find_spec(pkg) is None:
        missing.append(pkg)
if missing:
    raise SystemExit(
        "Missing Python packages: " + ", ".join(missing) +
        "\nPlease activate your survival environment, e.g. source activate survival_clock"
    )
PY

# ------------------------------------------------------------
# Run ADNI brain MRI AD L'EPOCH
# ------------------------------------------------------------

python3 "${PY_SCRIPT}" \
  --input-file "${ADNI_FILE}" \
  --outdir "${OUTDIR}" \
  --prefix "adni_brain_mri_ad_lepoch" \
  --id-col "PTID" \
  --visit-col "Visit_Code" \
  --date-col "Date" \
  --dx-col "DX_Binary" \
  --baseline-dx "CN" \
  --event-dx "MCI,AD" \
  --covariates "Age,Sex,DLICV,SITE" \
  --test-size 0.20 \
  --validation-size 0.20 \
  --random-state 20260707 \
  --stratify-age-bins 5 \
  --max-feature-missing 0.30 \
  --min-baseline-roi-fraction 0.80 \
  --l1-ratios "0.1,0.25,0.5,0.75,1.0" \
  --n-alphas 120 \
  --alpha-min-ratio 0.001 \
  --min-nonzero-brain-features 5 \
  --min-followup-days 1 \
  --n-bootstrap-incremental 1000

echo "Finished ADNI brain MRI AD L'EPOCH."
echo "Results saved to: ${OUTDIR}"

conda deactivate