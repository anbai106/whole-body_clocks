#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=apply_metab_clock_i1
#SBATCH --mem-per-cpu=24G
#SBATCH --array=0-3
#SBATCH --output=/cbica/home/wenju/output/apply_metabolomics_clock_instance1_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/apply_metabolomics_clock_instance1_%A_%a.err

set -euo pipefail

source activate survival_clock

organs=(Endocrine Digestive Hepatic Immune)
organ_label=${organs[$SLURM_ARRAY_TASK_ID]}
organ_clean=$(echo "${organ_label}" | tr '[:upper:]' '[:lower:]')

script="/cbica/home/wenju/Project/whole-body_clocks/mortality_clock/Metabolomics/longitudinal/organ_metabolomics_clock.py"

input_tsv="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/longitudinal/metabolomics/data/metabolomics${organ_label}_1_0.tsv"

model_joblib="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/${organ_label}_metabolomics_mortality_clock/${organ_clean}_metabolomics_mortality_clock_model.joblib"

outdir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/longitudinal/metabolomics/${organ_label}"
mkdir -p "${outdir}"

exec > "${outdir}/apply_${organ_clean}_metabolomics_mortality_clock_instance1_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out" \
     2> "${outdir}/apply_${organ_clean}_metabolomics_mortality_clock_instance1_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"

echo "============================================================"
echo "Applying pretrained metabolomics mortality clock to instance 1"
echo "Organ label: ${organ_label}"
echo "Organ clean: ${organ_clean}"
echo "Model: ${model_joblib}"
echo "Input TSV: ${input_tsv}"
echo "Output directory: ${outdir}"
echo "Started at: $(date)"
echo "============================================================"

if [ ! -f "${script}" ]; then
  echo "ERROR: Python script not found: ${script}"
  exit 1
fi

if [ ! -f "${model_joblib}" ]; then
  echo "ERROR: pretrained model joblib not found: ${model_joblib}"
  exit 1
fi

if [ ! -f "${input_tsv}" ]; then
  echo "ERROR: input instance-1 TSV not found: ${input_tsv}"
  exit 1
fi

pref="${organ_clean}_metabolomics_mortality_clock"

echo "Removing stale output files before rerun..."
rm -f "${outdir}/${pref}_apply_instance_1_0_predictions.tsv"
rm -f "${outdir}/${pref}_apply_instance_1_0_summary.json"
rm -f "${outdir}/${pref}_apply_longitudinal_instances_combined_predictions.tsv"
rm -f "${outdir}/${pref}_apply_longitudinal_instances_combined_summary.json"

python -u "${script}" \
  --organ "${organ_label}" \
  --model-joblib "${model_joblib}" \
  --input-tsv "1_0:${input_tsv}" \
  --covariate-csv /cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv \
  --death-xlsx /cbica/home/wenju/Dataset/UKBB_UMelbourne/Death_related_var_from_Ye.xlsx \
  --id-match-csv /cbica/home/wenju/Dataset/UKBB_UMelbourne/UKB_UMelbourne_vs_Penn_match_key.csv \
  --admin-censor-date 2022-11-30 \
  --outdir "${outdir}" \
  --application-instance 1 \
  --model-instance 0 \
  --application-session-id none \
  --risk-times 5,10,15

prediction_file="${outdir}/${pref}_apply_instance_1_0_predictions.tsv"

echo "============================================================"
echo "Checking final clean output"
echo "Prediction file: ${prediction_file}"
echo "============================================================"

if [ ! -f "${prediction_file}" ]; then
  echo "ERROR: prediction file was not generated: ${prediction_file}"
  exit 1
fi

required_cols=(
  "participant_id"
  "application_instance"
  "application_source_file"
  "${organ_clean}_metabolomics_mortality_risk_score"
  "sample_date"
  "death_date"
  "age_at_baseline"
  "age_at_imaging"
  "sex"
  "bmi_at_baseline"
  "diastolic_bp_at_baseline"
  "systolic_bp_at_baseline"
  "smoking_status_at_baseline"
  "uk_biobank_assessment_centre_f54_0_0"
  "risk_5y"
  "risk_10y"
  "risk_15y"
  "${organ_clean}_metabolomics_mortality_clock_acceleration_z"
  "${organ_clean}_metabolomics_mortality_clock_acceleration_years"
  "${organ_clean}_metabolomics_mortality_clock_age_years"
  "n_model_metabolomics_features_expected"
  "n_model_metabolomics_features_present_in_input"
  "n_model_metabolomics_features_missing_from_input"
)

for col in "${required_cols[@]}"; do
  if ! head -n 1 "${prediction_file}" | tr '\t' '\n' | grep -Fxq "${col}"; then
    echo "ERROR: required column missing: ${col}"
    echo "Full header:"
    head -n 1 "${prediction_file}" | tr '\t' '\n'
    exit 1
  fi
done

bad_pattern="${organ_clean}_metabolomics_mortality_clock_acceleration_${organ_clean}_metabolomics_mortality_clock_acceleration"
if head -n 1 "${prediction_file}" | grep -q "${bad_pattern}"; then
  echo "ERROR: malformed duplicated acceleration column detected."
  head -n 1 "${prediction_file}" | tr '\t' '\n' | grep "${organ_clean}_metabolomics_mortality_clock_acceleration" || true
  exit 1
fi

if head -n 1 "${prediction_file}" | tr '\t' '\n' | grep -Fxq "session_id"; then
  echo "ERROR: session_id should not be in final prediction file."
  exit 1
fi

if head -n 1 "${prediction_file}" | tr '\t' '\n' | grep -Fxq "diagnosis"; then
  echo "ERROR: diagnosis should not be in final prediction file."
  exit 1
fi

echo "Checking that original input metabolomics features are NOT included..."
head -n 1 "${input_tsv}" | tr '\t' '\n' | grep -v -E '^(participant_id|session_id|diagnosis)$' | while IFS= read -r feature; do
  if head -n 1 "${prediction_file}" | tr '\t' '\n' | grep -Fxq "${feature}"; then
    echo "ERROR: input metabolomics feature is still present in final prediction file: ${feature}"
    exit 1
  fi
done

echo "Verified clean output. Final header:"
head -n 1 "${prediction_file}" | tr '\t' '\n'

echo "============================================================"
echo "Finished applying pretrained metabolomics mortality clock"
echo "Organ: ${organ_label}"
echo "Finished at: $(date)"
echo "============================================================"

conda deactivate