#!/usr/bin/env bash
#SBATCH --job-name=adni_lepoch_csf
#SBATCH --output=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch/logs/adni_lepoch_csf_%j.out
#SBATCH --error=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch/logs/adni_lepoch_csf_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G

set -euo pipefail

# ============================================================
# ADNI baseline AD L'EPOCH acceleration_z vs CSF biomarkers
#
# CSF outcomes:
#   Abeta_CSF
#   Tau_CSF
#   PTau_CSF
#
# Predictor:
#   adni_brain_mri_ad_lepoch_acceleration_z
#
# Groups:
#   Non-event & censored
#   CN-MCI
#   CN-AD
#
# Covariates:
#   Age
#   Sex
#   DLICV / ICV
#   APOE4
# ============================================================

WORKDIR="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch"
LOGDIR="${WORKDIR}/logs"
OUTDIR="${WORKDIR}/baseline_associations/csf_biomarkers"

mkdir -p "${LOGDIR}" "${OUTDIR}"

PY_SCRIPT="/cbica/home/wenju/Project/whole-body_clocks/ADNI/ad_epoch_vs_CSF.py"

PRED_FILE="${WORKDIR}/results_brain_mri_ad_lepoch/adni_brain_mri_ad_lepoch_predictions.tsv"
ADNI_FILE="${WORKDIR}/adni_istaging.tsv"

source activate survival_clock

python3 "${PY_SCRIPT}" \
  --predictions-file "${PRED_FILE}" \
  --adni-file "${ADNI_FILE}" \
  --outdir "${OUTDIR}" \
  --prefix "adni_brain_mri_ad_lepoch" \
  --id-col "PTID" \
  --visit-col "Visit_Code" \
  --date-col "Date" \
  --dx-col "DX_Binary" \
  --accel-col "adni_brain_mri_ad_lepoch_acceleration_z" \
  --age-col "Age" \
  --sex-col "Sex" \
  --icv-col "DLICV" \
  --apoe4-col "APOE4" \
  --csf-vars "Abeta_CSF,Tau_CSF,PTau_CSF" \
  --max-csf-baseline-distance-days 365 \
  --min-n 8

echo "Finished ADNI L'EPOCH baseline CSF biomarker association analyses."
echo "Output directory: ${OUTDIR}"

conda deactivate