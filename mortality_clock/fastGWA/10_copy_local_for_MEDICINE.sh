#!/bin/bash

set -euo pipefail

# ============================================================
# Copy Manhattan and QQ plots for 22 mortality EPOCH clocks
# from CUBIC to local Mac folder
#
# Run locally:
#   bash rsync_mortality_qmplot_from_cubic.sh
# ============================================================

REMOTE_USER="wenju"
REMOTE_HOST="cubic-login.uphs.upenn.edu"

REMOTE_BASE="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/fastGWA/output"

LOCAL_DIR="/Users/hao/Downloads/EPOCH_result/mortality_qmplot"

mkdir -p "${LOCAL_DIR}"

CLOCKS=(
  "adipose_mri_mortality_clock"
  "brain_mri_mortality_clock"
  "Brain_proteomics_mortality_clock"
  "Digestive_metabolomics_mortality_clock"
  "Endocrine_metabolomics_mortality_clock"
  "Endocrine_proteomics_mortality_clock"
  "Eye_proteomics_mortality_clock"
  "heart_mri_mortality_clock"
  "Heart_proteomics_mortality_clock"
  "Hepatic_metabolomics_mortality_clock"
  "Hepatic_proteomics_mortality_clock"
  "Immune_metabolomics_mortality_clock"
  "Immune_proteomics_mortality_clock"
  "kidney_mri_mortality_clock"
  "liver_mri_mortality_clock"
  "pancreas_mri_mortality_clock"
  "Pulmonary_proteomics_mortality_clock"
  "Renal_proteomics_mortality_clock"
  "Reproductive_female_proteomics_mortality_clock"
  "Reproductive_male_proteomics_mortality_clock"
  "Skin_proteomics_mortality_clock"
  "spleen_mri_mortality_clock"
)

echo "============================================================"
echo "Copying mortality EPOCH qmplot PNG files"
echo "Remote: ${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_BASE}"
echo "Local:  ${LOCAL_DIR}"
echo "Number of clocks: ${#CLOCKS[@]}"
echo "============================================================"

for clock in "${CLOCKS[@]}"; do
  echo "------------------------------------------------------------"
  echo "Clock: ${clock}"

  remote_qq="${REMOTE_BASE}/${clock}/QQ_plot.png"
  remote_manhattan="${REMOTE_BASE}/${clock}/manhattan_qmplot.png"

  local_qq="${LOCAL_DIR}/${clock}_QQ_plot.png"
  local_manhattan="${LOCAL_DIR}/${clock}_manhattan_qmplot.png"

  echo "Copying QQ plot..."
  rsync -avz \
    "${REMOTE_USER}@${REMOTE_HOST}:${remote_qq}" \
    "${local_qq}"

  echo "Copying Manhattan plot..."
  rsync -avz \
    "${REMOTE_USER}@${REMOTE_HOST}:${remote_manhattan}" \
    "${local_manhattan}"
done

echo "============================================================"
echo "Finished copying all mortality EPOCH qmplot PNG files."
echo "Local output folder:"
echo "  ${LOCAL_DIR}"
echo "Number of PNG files copied:"
find "${LOCAL_DIR}" -maxdepth 1 -name "*mortality_clock*png" | wc -l
echo "============================================================"