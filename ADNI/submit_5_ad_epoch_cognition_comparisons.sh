#!/bin/bash
#SBATCH --job-name=adni_epoch_cog_exact
#SBATCH --output=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch/baseline_associations/cognition_biomarker_comparison/logs/%x_%j.out
#SBATCH --error=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch/baseline_associations/cognition_biomarker_comparison/logs/%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1

set -euo pipefail

source activate survival_clock

PY_SCRIPT="/cbica/home/wenju/Project/whole-body_clocks/ADNI/ad_epoch_cognition.py"

ADNI_TSV="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch/adni_istaging.tsv"

BASELINE_EPOCH_DIR="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch/results_brain_mri_ad_lepoch"

LONGITUDINAL_EPOCH_DIR="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch/results_brain_mri_ad_lepoch_longitudinal_cn_only"

OUTDIR="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch/baseline_associations/cognition_biomarker_comparison"

mkdir -p "${OUTDIR}/logs"

python "${PY_SCRIPT}" \
  --adni_tsv "${ADNI_TSV}" \
  --baseline_epoch "${BASELINE_EPOCH_DIR}" \
  --longitudinal_epoch "${LONGITUDINAL_EPOCH_DIR}" \
  --outdir "${OUTDIR}" \
  --baseline_epoch_col "auto" \
  --longitudinal_epoch_col "auto" \
  --window_days 365 \
  --n_perm 10000 \
  --min_n 30 \
  --min_slope_scans 2 \
  --min_slope_followup_years 0.5 \
  --seed 2026

echo "Finished ADNI AD EPOCH cognition-biomarker comparison."
echo "Output directory: ${OUTDIR}"