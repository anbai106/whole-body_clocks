#!/usr/bin/env bash
#SBATCH --job-name=ad_epoch_mci2ad_biomarkers
#SBATCH --output=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch/logs/ad_epoch_mci2ad_biomarkers_%j.out
#SBATCH --error=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch/logs/ad_epoch_mci2ad_biomarkers_%j.err
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --partition=bioinformatics

# Compare the baseline AD EPOCH score with clinical, CSF, SPARE-AD, and
# SPARE-BA predictors of MCI-to-AD conversion.
#
# Clinical model (no baseline cognition):
#   Age + Sex + Education_Years + APOE_Genotype
#
# Biomarkers:
#   Abeta_CSF + Tau_CSF + PTau_CSF + SPARE_AD + SPARE_BA
#
# The MCI-to-AD survival file was generated after excluding original AD EPOCH
# train/validation participants and contains one selected baseline-MCI MRI row
# per participant.

WORKDIR="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch"
LOGDIR="${WORKDIR}/logs"
OUTDIR="${WORKDIR}/results_brain_mri_ad_lepoch_mci2ad_biomarker_survival"

PY_SCRIPT="/cbica/home/wenju/Project/whole-body_clocks/ADNI/ad_epoch_mci2ad_biomarker_survival.py"
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
  --prefix "adni_mci2ad_epoch_biomarkers" \
  --id-col "PTID" \
  --visit-col "Visit_Code" \
  --date-col "Date" \
  --baseline-visit-col "selected_baseline_visit_code" \
  --baseline-date-col "selected_baseline_date" \
  --time-col "time_years" \
  --event-col "event" \
  --epoch-col "adni_brain_mri_ad_lepoch_risk_score" \
  --clinical-cols "Age,Sex,Education_Years,APOE_Genotype" \
  --csf-cols "Abeta_CSF,Tau_CSF,PTau_CSF" \
  --spare-cols "SPARE_AD,SPARE_BA" \
  --csf-lookback-days 365 \
  --cv-folds 5 \
  --random-state 20260730 \
  --cox-penalizer 0.01 \
  --bootstrap 1000 \
  --minimum-n 50 \
  --minimum-events 20

echo "Finished MCI-to-AD AD EPOCH biomarker survival analyses."
echo "Results saved to: ${OUTDIR}"
echo "Performance: ${OUTDIR}/adni_mci2ad_epoch_biomarkers_performance.tsv"
echo "Comparisons: ${OUTDIR}/adni_mci2ad_epoch_biomarkers_comparisons.tsv"

conda deactivate
