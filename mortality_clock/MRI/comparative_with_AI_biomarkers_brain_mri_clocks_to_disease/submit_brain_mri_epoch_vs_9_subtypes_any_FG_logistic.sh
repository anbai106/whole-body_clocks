#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=BrainEPOCH_9sub_logistic
#SBATCH --mem=32G
#SBATCH --time=0-11:00:00
#SBATCH --output=/cbica/home/wenju/output/BrainEPOCH_9sub_logistic_%j.out
#SBATCH --error=/cbica/home/wenju/output/BrainEPOCH_9sub_logistic_%j.err

module load python/anaconda/3
source activate survival

mkdir -p /cbica/home/wenju/output

PY_SCRIPT="/cbica/home/wenju/Project/whole-body_clocks/mortality_clock/MRI/comparative_with_AI_biomarkers_brain_mri_clocks_to_disease/brain_mri_epoch_vs_9_subtypes_any_FG_logistic.py"

EPOCH_TSV="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/brain_mri_mortality_clock/brain_mri_mortality_clock_predictions.tsv"

SUBTYPE_TSV="/cbica/projects/MULTI/processed/UKBB/derived_AI_biomakers_across_projects/UKBB_487894_participant_58_biomarker_matched_ID.tsv"

ICD10_CSV="/cbica/home/wenju/Reproducibile_paper/BrainEye/data/UKBB_fullsample_ICD10.csv"

UMEL_DEATH_XLSX="/cbica/home/wenju/Dataset/UKBB_UMelbourne/Death_related_var_from_Ye.xlsx"

UMEL_MATCH_CSV="/cbica/home/wenju/Dataset/UKBB_UMelbourne/UKB_UMelbourne_vs_Penn_match_key.csv"

COV_TSV="/cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv"

BASE_OUT="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/brain_mri_mortality_clock/comparative_with_9_disease_subtypes_any_FG/logistic_regression"

OBSERVED_OUT="${BASE_OUT}/observed_followup"

FIXED3_OUT="${BASE_OUT}/fixed_3y"

mkdir -p "${OBSERVED_OUT}"
mkdir -p "${FIXED3_OUT}"

echo "======================================================================"
echo "Brain MRI mortality EPOCH vs 9 AI disease subtype scores"
echo "Binary logistic-regression sensitivity analysis"
echo "Primary EPOCH population: held-out test split"
echo "Fairness: exact same common complete-case sample within each endpoint"
echo "======================================================================"

echo ""
echo "======================================================================"
echo "ANALYSIS 1: observed-follow-up case/control logistic regression"
echo "Case = incident F/G diagnosis after MRI and before censoring"
echo "Control = no F/G diagnosis through observed censoring"
echo "======================================================================"

python "${PY_SCRIPT}" \
  --epoch_tsv "${EPOCH_TSV}" \
  --subtype_tsv "${SUBTYPE_TSV}" \
  --subtype_id_col "id_upenn" \
  --icd10_csv "${ICD10_CSV}" \
  --umel_death_xlsx "${UMEL_DEATH_XLSX}" \
  --umel_match_csv "${UMEL_MATCH_CSV}" \
  --cov_tsv "${COV_TSV}" \
  --output_dir "${OBSERVED_OUT}" \
  --split test \
  --endpoint_mode observed_followup \
  --admin_censor_date "2022-11-30" \
  --min_cases 20 \
  --min_controls 20 \
  --cv_folds 5 \
  --cv_repeats 5 \
  --bootstrap 500 \
  --seed 20260818

echo ""
echo "Observed-follow-up logistic analysis completed successfully."

echo ""
echo "======================================================================"
echo "ANALYSIS 2: fixed 3-year case/control logistic sensitivity analysis"
echo "Case = incident F/G diagnosis within 3 years after MRI"
echo "Control = F/G-disease-free and observed through the full 3-year horizon"
echo "This analysis equalizes follow-up opportunity across cases/controls."
echo "======================================================================"

python "${PY_SCRIPT}" \
  --epoch_tsv "${EPOCH_TSV}" \
  --subtype_tsv "${SUBTYPE_TSV}" \
  --subtype_id_col "id_upenn" \
  --icd10_csv "${ICD10_CSV}" \
  --umel_death_xlsx "${UMEL_DEATH_XLSX}" \
  --umel_match_csv "${UMEL_MATCH_CSV}" \
  --cov_tsv "${COV_TSV}" \
  --output_dir "${FIXED3_OUT}" \
  --split test \
  --endpoint_mode fixed_horizon \
  --horizon_years 3 \
  --admin_censor_date "2022-11-30" \
  --min_cases 20 \
  --min_controls 20 \
  --cv_folds 5 \
  --cv_repeats 5 \
  --bootstrap 500 \
  --seed 20260818

echo ""
echo "Fixed 3-year logistic analysis completed successfully."

echo ""
echo "======================================================================"
echo "ALL LOGISTIC ANALYSES FINISHED SUCCESSFULLY"
echo "======================================================================"
echo "Observed-follow-up results:"
echo "  ${OBSERVED_OUT}/single_marker_unadjusted.tsv"
echo "  ${OBSERVED_OUT}/single_marker_adjusted.tsv"
echo "  ${OBSERVED_OUT}/epoch_vs_subtype_adjusted_joint.tsv"
echo "  ${OBSERVED_OUT}/combined_adjusted_models.tsv"
echo "  ${OBSERVED_OUT}/cv_model_metrics.tsv"
echo "  ${OBSERVED_OUT}/bootstrap_discrimination_comparisons.tsv"
echo ""
echo "Fixed 3-year sensitivity results:"
echo "  ${FIXED3_OUT}/single_marker_unadjusted.tsv"
echo "  ${FIXED3_OUT}/single_marker_adjusted.tsv"
echo "  ${FIXED3_OUT}/epoch_vs_subtype_adjusted_joint.tsv"
echo "  ${FIXED3_OUT}/combined_adjusted_models.tsv"
echo "  ${FIXED3_OUT}/cv_model_metrics.tsv"
echo "  ${FIXED3_OUT}/bootstrap_discrimination_comparisons.tsv"
echo "======================================================================"

conda deactivate
