#!/usr/bin/env bash
#SBATCH --job-name=apply_adni_ad_lepoch_cn
#SBATCH --output=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch/logs/apply_adni_ad_lepoch_cn_%j.out
#SBATCH --error=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch/logs/apply_adni_ad_lepoch_cn_%j.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G

set -euo pipefail

# ============================================================
# Apply pre-trained ADNI brain MRI AD L'EPOCH to longitudinal CN scans only
#
# This is application-only:
#   - no Cox model refitting
#   - no train/validation/test split
#   - loads saved AD L'EPOCH model joblib
#   - applies saved preprocessor and Cox model
#   - computes risk score, absolute risk, acceleration_z,
#     acceleration_years, and clock_age_years
#
# Primary rule:
#   Only score CN-labeled MRI scans before MCI/AD conversion.
#
# Time zero per participant:
#   first CN visit with usable MUSE GM ROI features;
#   prefer bl if bl has MUSE, otherwise earliest usable CN MRI visit.
#
# Scored scans:
#   selected_baseline
#   pre_event_CN_followup
#   censored_CN_followup
#
# Not scored:
#   MCI scans
#   AD scans
#   any scan after first MCI/AD conversion
# ============================================================

WORKDIR="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch"
LOGDIR="${WORKDIR}/logs"
OUTDIR="${WORKDIR}/results_brain_mri_ad_lepoch_longitudinal_cn_only"

mkdir -p "${WORKDIR}"
mkdir -p "${LOGDIR}"
mkdir -p "${OUTDIR}"

PY_SCRIPT="/cbica/home/wenju/Project/whole-body_clocks/ADNI/ad_epoch_apply.py"

ADNI_FILE="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch/adni_istaging.tsv"

MODEL_JOBLIB="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch/results_brain_mri_ad_lepoch/adni_brain_mri_ad_lepoch_model.joblib"

source activate survival_clock

python "${PY_SCRIPT}" \
  --input-file "${ADNI_FILE}" \
  --model-joblib "${MODEL_JOBLIB}" \
  --outdir "${OUTDIR}" \
  --prefix "adni_brain_mri_ad_lepoch" \
  --id-col "PTID" \
  --visit-col "Visit_Code" \
  --date-col "Date" \
  --dx-col "DX_Binary" \
  --baseline-dx "CN" \
  --event-dx "MCI,AD" \
  --eligible-scan-dx "CN" \
  --min-baseline-roi-fraction 0.80 \
  --min-scan-roi-fraction 0.80 \
  --age-update-mode "from_baseline_date" \
  --risk-times "1,2,3,5" \
  --include-selected-baseline

echo "Finished applying ADNI brain MRI AD L'EPOCH to longitudinal CN scans only."
echo "Results saved to: ${OUTDIR}"

conda deactivate