#!/usr/bin/env bash
#SBATCH --job-name=apoe_dementia_clock
#SBATCH --output=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/apoe_status_ukbb/logs/apoe_dementia_clock_%A_%a.out
#SBATCH --error=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/apoe_status_ukbb/logs/apoe_dementia_clock_%A_%a.err
#SBATCH --array=0-10
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G

set -euo pipefail

# ============================================================
# Parallel APOE e4/e4 vs e2/e2 analysis for 11 dementia clocks
# ============================================================

OUTDIR="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/apoe_status_ukbb"
LOGDIR="${OUTDIR}/logs"
RESULTDIR="${OUTDIR}/apoe_e4e4_vs_e2e2_dementia_clock_results"

mkdir -p "${LOGDIR}"
mkdir -p "${RESULTDIR}"

PY_SCRIPT="/cbica/home/wenju/Project/whole-body_clocks/dementia_clock/APOE/2_apoe_clock_bag.py"

APOE_FILE="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/apoe_status_ukbb/ukbb_apoe_status.tsv"
LEPOCH_FILE="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/all_disease_lepoch_incremental_value_scale_qc/stable_significant_disease_clock_acceleration_z_wide.tsv"
BAG_FILE="/cbica/home/wenju/Reproducibile_paper/SleepAging/data/MomoBAG.tsv"
COV_FILE="/cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv"

# ============================================================
# The 11 dementia L'EPOCH clocks
# ============================================================

LEPOCH_CLOCKS=(
  "dementia_brain_proteomics_clock_acceleration_z"
  "dementia_hepatic_proteomics_clock_acceleration_z"
  "dementia_immune_proteomics_clock_acceleration_z"
  "dementia_reproductive_female_proteomics_clock_acceleration_z"
  "dementia_reproductive_male_proteomics_clock_acceleration_z"
  "dementia_endocrine_proteomics_clock_acceleration_z"
  "dementia_heart_proteomics_clock_acceleration_z"
  "dementia_metabolic_metabolomics_clock_acceleration_z"
  "dementia_digestive_metabolomics_clock_acceleration_z"
  "dementia_endocrine_metabolomics_clock_acceleration_z"
  "dementia_hepatic_metabolomics_clock_acceleration_z"
)

# ============================================================
# Matched BAGs for effect-size fold comparison
# ============================================================

MATCHED_BAGS=(
  "Brain_ProtBAG"
  "Hepatic_ProtBAG"
  "Immune_ProtBAG"
  "Reproductive_female_ProtBAG"
  "Reproductive_male_ProtBAG"
  "Endocrine_ProtBAG"
  "Heart_ProtBAG"
  "Metabolic_MetBAG"
  "Digestive_MetBAG"
  "Endocrine_MetBAG"
  "Hepatic_MetBAG"
)

TASK_ID="${SLURM_ARRAY_TASK_ID}"

LEPOCH_COL="${LEPOCH_CLOCKS[${TASK_ID}]}"
BAG_COL="${MATCHED_BAGS[${TASK_ID}]}"

echo "============================================================"
echo "SLURM job ID: ${SLURM_JOB_ID}"
echo "Array task ID: ${TASK_ID}"
echo "L'EPOCH clock: ${LEPOCH_COL}"
echo "Matched BAG: ${BAG_COL}"
echo "============================================================"

source activate DNE
python "${PY_SCRIPT}" \
  --apoe-file "${APOE_FILE}" \
  --lepoch-file "${LEPOCH_FILE}" \
  --bag-file "${BAG_FILE}" \
  --cov-file "${COV_FILE}" \
  --outdir "${RESULTDIR}" \
  --lepoch-col "${LEPOCH_COL}" \
  --bag-col "${BAG_COL}" \
  --apoe-id-col "FID" \
  --cov-id-col "eid" \
  --lepoch-id-col "participant_id" \
  --bag-id-col "participant_id"

echo "Finished task ${TASK_ID}: ${LEPOCH_COL}"
conda deactivate || true
