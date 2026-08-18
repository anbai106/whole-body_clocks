#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=BrainEPOCH_vs_9sub_FG
#SBATCH --mem=32G
#SBATCH --time=0-11:59:00
#SBATCH --output=/cbica/home/wenju/output/BrainEPOCH_vs_9sub_FG_%j.out
#SBATCH --error=/cbica/home/wenju/output/BrainEPOCH_vs_9sub_FG_%j.err

module load python/anaconda/3
source activate survival

mkdir -p /cbica/home/wenju/output

PY_SCRIPT="/cbica/home/wenju/Project/whole-body_clocks/mortality_clock/MRI/comparative_with_AI_biomarkers_brain_mri_clocks_to_disease/brain_mri_epoch_vs_9_subtypes_any_FG_survival.py"
EPOCH_TSV="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/brain_mri_mortality_clock/brain_mri_mortality_clock_predictions.tsv"
SUBTYPE_TSV="/cbica/projects/MULTI/processed/UKBB/derived_AI_biomakers_across_projects/UKBB_487894_participant_58_biomarker_matched_ID.tsv"
ICD10_CSV="/cbica/home/wenju/Reproducibile_paper/BrainEye/data/UKBB_fullsample_ICD10.csv"
UMEL_DEATH_XLSX="/cbica/home/wenju/Dataset/UKBB_UMelbourne/Death_related_var_from_Ye.xlsx"
UMEL_MATCH_CSV="/cbica/home/wenju/Dataset/UKBB_UMelbourne/UKB_UMelbourne_vs_Penn_match_key.csv"
COV_TSV="/cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv"
OUT_DIR="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/brain_mri_mortality_clock/comparative_with_9_disease_subtypes_any_FG"

mkdir -p "${OUT_DIR}"

echo "============================================================"
echo "Brain MRI mortality EPOCH vs 9 AI disease subtype scores"
echo "Endpoint: first incident inpatient ICD-10 F* or G* diagnosis"
echo "Primary analysis: held-out EPOCH test split"
echo "Fairness: one strict common complete-case sample for all predictors"
echo "============================================================"

python "${PY_SCRIPT}" \
  --epoch_tsv "${EPOCH_TSV}" \
  --subtype_tsv "${SUBTYPE_TSV}" \
  --subtype_id_col "id_upenn" \
  --icd10_csv "${ICD10_CSV}" \
  --umel_death_xlsx "${UMEL_DEATH_XLSX}" \
  --umel_match_csv "${UMEL_MATCH_CSV}" \
  --cov_tsv "${COV_TSV}" \
  --output_dir "${OUT_DIR}" \
  --split test \
  --admin_censor_date "2022-11-30" \
  --min_events 20 \
  --bootstrap 200 \
  --seed 20260818

echo "============================================================"
echo "Finished. Key outputs:"
echo "  ${OUT_DIR}/marker_models_common_sample.tsv"
echo "  ${OUT_DIR}/epoch_vs_subtype_pairwise.tsv"
echo "  ${OUT_DIR}/combined_subtypes_vs_epoch.tsv"
echo "  ${OUT_DIR}/bootstrap_cindex_common_sample.tsv"
echo "  ${OUT_DIR}/sample_flow.tsv"
echo "  ${OUT_DIR}/cohort_qc.tsv"
echo "============================================================"

conda deactivate
