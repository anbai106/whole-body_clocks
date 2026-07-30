#!/usr/bin/env bash
#SBATCH --job-name=ad_epoch_mci2ad_csf_epoch
#SBATCH --output=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch/logs/ad_epoch_mci2ad_csf_epoch_%j.out
#SBATCH --error=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch/logs/ad_epoch_mci2ad_csf_epoch_%j.err
#SBATCH --time=05:00:00
#SBATCH --mem=24G
#SBATCH --partition=bioinformatics

# ============================================================
# Cumulative MCI-to-AD CSF + AD EPOCH survival analysis
#
# Cohort:
#   Non-training baseline-MCI participants from the previously generated
#   MCI-to-AD EPOCH survival file. The application workflow should already
#   have excluded original AD EPOCH train/validation participants.
#
# Outcome:
#   First AD diagnosis after selected baseline-MCI MRI.
#
# Cumulative order:
#   M0: Age + Sex
#   M1: M0 + Abeta_CSF
#   M2: M1 + Tau_CSF + PTau_CSF
#   M3: M2 + AD EPOCH
#
# Relaxed CSF rule:
#   First chronologically available nonmissing measurement per subject and
#   marker, regardless of whether it occurs before or after baseline MRI.
#   Timing relative to baseline and event/censoring is saved in QC outputs.
#
# Primary inference:
#   All cumulative models are compared on one common complete-case cohort.
# ============================================================

WORKDIR="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch"
LOGDIR="${WORKDIR}/logs"
OUTDIR="${WORKDIR}/results_brain_mri_ad_lepoch_mci2ad_cumulative_csf_epoch"

PY_SCRIPT="/cbica/home/wenju/Project/whole-body_clocks/ADNI/ad_epoch_mci2ad_cumulative_csf_epoch.py"
SURVIVAL_FILE="${WORKDIR}/results_brain_mri_ad_lepoch_mci2ad_application/adni_brain_mri_ad_lepoch_mci2ad_epoch_survival.tsv"
ADNI_FILE="${WORKDIR}/adni_istaging.tsv"

mkdir -p "${WORKDIR}" "${LOGDIR}" "${OUTDIR}"

if command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)"
fi
conda activate survival_clock

for required_file in "${PY_SCRIPT}" "${SURVIVAL_FILE}" "${ADNI_FILE}"; do
  if [[ ! -s "${required_file}" ]]; then
    echo "ERROR: required file is missing or empty: ${required_file}" >&2
    exit 1
  fi
done

python3 "${PY_SCRIPT}" \
  --survival-file "${SURVIVAL_FILE}" \
  --adni-file "${ADNI_FILE}" \
  --outdir "${OUTDIR}" \
  --prefix "adni_mci2ad_cumulative_csf_epoch" \
  --id-col "PTID" \
  --visit-col "Visit_Code" \
  --date-col "Date" \
  --baseline-visit-col "selected_baseline_visit_code" \
  --baseline-date-col "selected_baseline_date" \
  --time-col "time_years" \
  --event-col "event" \
  --epoch-col "adni_brain_mri_ad_lepoch_risk_score" \
  --baseline-cols "Age,Sex" \
  --csf-cols "Abeta_CSF,Tau_CSF,PTau_CSF" \
  --cv-folds 5 \
  --random-state 20260730 \
  --cox-penalizer 0.01 \
  --bootstrap 1000 \
  --minimum-n 50 \
  --minimum-events 20

# Optional secondary analysis using each step's maximum available sample:
# add --allow-stage-specific-samples to the command above.

echo "Finished cumulative MCI-to-AD CSF + AD EPOCH survival analysis."
echo "Results saved to: ${OUTDIR}"
echo "Performance: ${OUTDIR}/adni_mci2ad_cumulative_csf_epoch_performance.tsv"
echo "Comparisons: ${OUTDIR}/adni_mci2ad_cumulative_csf_epoch_comparisons.tsv"

conda deactivate
