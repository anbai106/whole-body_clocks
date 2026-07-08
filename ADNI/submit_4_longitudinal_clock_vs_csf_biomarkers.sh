#!/usr/bin/env bash
#SBATCH --job-name=adni_lepoch_long_csf
#SBATCH --output=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch/logs/adni_lepoch_long_csf_%j.out
#SBATCH --error=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch/logs/adni_lepoch_long_csf_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G

set -euo pipefail

# ============================================================
# Test whether baseline CSF biomarkers predict longitudinal
# change in AD L'EPOCH clocks among CN-only follow-up scans.
#
# Longitudinal input:
#   ADNI CN-only longitudinal AD L'EPOCH predictions
#
# Baseline biomarker input:
#   ADNI baseline CSF analysis dataset with:
#     Abeta_CSF
#     Tau_CSF
#     PTau_CSF
#     Age, Sex, ICV, APOE4
#     conversion group
#
# Primary model:
#   Clock change ~ baseline CSF biomarker
#                  + baseline clock
#                  + conversion_group
#                  + Age + Sex + ICV + APOE4
#                  + follow-up span
#
# Main outputs:
#   subject-level clock-change metrics
#   combined full-sample association results
#   group-specific exploratory association results
# ============================================================

WORKDIR="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch"
LOGDIR="${WORKDIR}/logs"

OUTDIR="${WORKDIR}/longitudinal_change_vs_baseline_csf"

mkdir -p "${LOGDIR}" "${OUTDIR}"

PY_SCRIPT="/cbica/home/wenju//Users/hao/Project/whole-body_clocks/ADNI/ad_epoch_vs_CSF_longitudinal.py"

LONG_FILE="${WORKDIR}/results_brain_mri_ad_lepoch_longitudinal_cn_only/adni_brain_mri_ad_lepoch_longitudinal_cn_only_predictions.tsv"

BASELINE_CSF_FILE="${WORKDIR}/baseline_associations/csf_biomarkers/adni_brain_mri_ad_lepoch_baseline_csf_analysis_dataset.tsv"

source activate survival_clock

python3 "${PY_SCRIPT}" \
  --longitudinal-file "${LONG_FILE}" \
  --baseline-csf-file "${BASELINE_CSF_FILE}" \
  --outdir "${OUTDIR}" \
  --prefix "adni_brain_mri_ad_lepoch" \
  --id-col "PTID" \
  --time-col "years_since_selected_baseline" \
  --group-col "conversion_group_3level" \
  --clock-cols "adni_brain_mri_ad_lepoch_acceleration_z,adni_brain_mri_ad_lepoch_acceleration_years,adni_brain_mri_ad_lepoch_risk_score" \
  --csf-vars "Abeta_CSF,Tau_CSF,PTau_CSF" \
  --age-col "_age" \
  --sex-col "_sex_male" \
  --icv-col "_icv" \
  --apoe4-col "_apoe4" \
  --min-scans 2 \
  --min-followup-years 0.25 \
  --min-n 8

echo "Finished longitudinal AD L'EPOCH change vs baseline CSF analyses."
echo "Output directory: ${OUTDIR}"

conda deactivate