#!/usr/bin/env bash
#SBATCH --job-name=ad_epoch_cn_mciad_cumulative
#SBATCH --output=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch/logs/ad_epoch_cn_mciad_cumulative_%j.out
#SBATCH --error=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch/logs/ad_epoch_cn_mciad_cumulative_%j.err
#SBATCH --time=00:03:00
#SBATCH --mem=24G
#SBATCH --partition=bioinformatics


# ============================================================
# ADNI cumulative CN-to-MCI/AD survival analysis
#
# Common complete-case cohort and cumulative order:
#   M0: Age + Sex
#   M1: + Abeta_CSF
#   M2: + Tau_CSF + PTau_CSF
#   M3: + SPARE_BA
#   M4: + SPARE_AD
#   M5: + AD EPOCH risk score
#
# CSF rule:
#   first chronologically available nonmissing value per participant/marker.
#   Timing relative to the baseline MRI is written to QC output.
#
# The full original AD EPOCH prediction cohort is used, including train,
# validation, and test splits. No split is excluded.
# ============================================================

WORKDIR="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch"
LOGDIR="${WORKDIR}/logs"
OUTDIR="${WORKDIR}/results_brain_mri_ad_lepoch_cn_mciad_cumulative_biomarkers"

mkdir -p "${WORKDIR}" "${LOGDIR}" "${OUTDIR}"

PY_SCRIPT="/cbica/home/wenju/Project/whole-body_clocks/ADNI/ad_epoch_cn_mciad_cumulative_biomarkers.py"

# Original participant-level predictions from AD EPOCH development.
# This file contains the selected baseline-CN visit, CN-to-MCI/AD outcome,
# split labels, and frozen AD EPOCH score for all development participants.
EPOCH_PREDICTIONS="${WORKDIR}/results_brain_mri_ad_lepoch/adni_brain_mri_ad_lepoch_predictions.tsv"

# Full longitudinal ADNI/iSTAGING table used to retrieve baseline SPARE scores
# and each participant's first available CSF amyloid/tau measurements.
ADNI_FILE="${WORKDIR}/adni_istaging.tsv"

source activate survival_clock

python3 "${PY_SCRIPT}" \
  --epoch-predictions-file "${EPOCH_PREDICTIONS}" \
  --adni-file "${ADNI_FILE}" \
  --outdir "${OUTDIR}" \
  --prefix "adni_cn_mciad_cumulative_biomarkers" \
  --id-col "PTID" \
  --visit-col "Visit_Code" \
  --date-col "Date" \
  --baseline-visit-col "selected_baseline_visit_code" \
  --baseline-date-col "selected_baseline_date" \
  --time-col "time_years" \
  --event-col "event" \
  --epoch-col "adni_brain_mri_ad_lepoch_risk_score" \
  --age-col "Age" \
  --sex-col "Sex" \
  --amyloid-col "Abeta_CSF" \
  --tau-cols "Tau_CSF,PTau_CSF" \
  --spare-ba-col "SPARE_BA" \
  --spare-ad-col "SPARE_AD" \
  --cv-folds 5 \
  --random-state 20260730 \
  --cox-penalizer 0.01 \
  --bootstrap 1000 \
  --minimum-n 50 \
  --minimum-events 20

# Add --allow-stage-specific-samples above only for a secondary descriptive
# analysis. The primary cumulative result always uses a single common sample.

echo "Finished cumulative CN-to-MCI/AD biomarker analysis."
echo "Results saved to: ${OUTDIR}"

conda deactivate
