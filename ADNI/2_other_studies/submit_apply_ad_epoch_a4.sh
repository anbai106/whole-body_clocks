#!/usr/bin/env bash
#SBATCH --job-name=apply_ad_epoch_a4
#SBATCH --partition=all
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --output=/cbica/home/wenju/output/apply_ad_epoch_a4_%j.out
#SBATCH --error=/cbica/home/wenju/output/apply_ad_epoch_a4_%j.err

# Apply the pretrained ADNI brain MRI AD EPOCH model to longitudinal A4 MUSE data.
# The model and saved preprocessing pipeline are never refitted.

WORKDIR="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch"
LOGDIR="/cbica/home/wenju/output"

PY_SCRIPT="/cbica/home/wenju/Project/whole-body_clocks/ADNI/2_other_studies/apply_ad_epoch_a4.py"
MODEL_JOBLIB="${WORKDIR}/results_brain_mri_ad_lepoch/adni_brain_mri_ad_lepoch_model.joblib"

MUSE_FILE="/cbica/projects/MULTI/processed/A4/tsv/MUSE.csv"
CLINICAL_ROOT="/cbica/projects/MULTI/download/A4/Clinical"
SUBJINFO_FILE="${CLINICAL_ROOT}/Derived_Data/SUBJINFO.csv"
MRI_VISITS_FILE="${CLINICAL_ROOT}/External_Data/imaging_volumetric_mri.csv"
SV_FILE="${CLINICAL_ROOT}/Derived_Data/SV.csv"

OUTDIR="${WORKDIR}/results_a4_longitudinal_ad_epoch"
PREFIX="a4_adni_brain_mri_ad_epoch"

MIN_ROI_FRACTION="${MIN_ROI_FRACTION:-0.80}"
RISK_TIMES="${RISK_TIMES:-1,2,3,5}"
SITE_LABEL="${SITE_LABEL:-A4}"

mkdir -p "${LOGDIR}" "${OUTDIR}"

source activate survival_clock

for file in \
    "${PY_SCRIPT}" \
    "${MODEL_JOBLIB}" \
    "${MUSE_FILE}" \
    "${SUBJINFO_FILE}" \
    "${MRI_VISITS_FILE}" \
    "${SV_FILE}"; do
    if [[ ! -s "${file}" ]]; then
        echo "ERROR: missing or empty file: ${file}" >&2
        exit 1
    fi
done

python3 "${PY_SCRIPT}" \
    --muse-file "${MUSE_FILE}" \
    --subjinfo-file "${SUBJINFO_FILE}" \
    --mri-visits-file "${MRI_VISITS_FILE}" \
    --sv-file "${SV_FILE}" \
    --model-joblib "${MODEL_JOBLIB}" \
    --outdir "${OUTDIR}" \
    --prefix "${PREFIX}" \
    --muse-id-col "ID" \
    --dlicv-label "702" \
    --site-label "${SITE_LABEL}" \
    --study-label "A4" \
    --min-roi-fraction "${MIN_ROI_FRACTION}" \
    --risk-times "${RISK_TIMES}"

echo "============================================================"
echo "A4 AD EPOCH application completed."
echo "Predictions: ${OUTDIR}/${PREFIX}_scan_level_predictions.tsv"
echo "Subject summary: ${OUTDIR}/${PREFIX}_subject_longitudinal_summary.tsv"
echo "Visit audit: ${OUTDIR}/${PREFIX}_visit_match_audit.tsv"
echo "ROI audit: ${OUTDIR}/${PREFIX}_model_roi_mapping_audit.tsv"
echo "Excluded scans: ${OUTDIR}/${PREFIX}_excluded_scans.tsv"
echo "============================================================"

conda deactivate
