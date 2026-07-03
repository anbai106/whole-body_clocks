#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=apply_mri_clock_i3
#SBATCH --mem-per-cpu=24G
#SBATCH --array=0-1
#SBATCH --output=/cbica/home/wenju/output/apply_mri_clock_instance3_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/apply_mri_clock_instance3_%A_%a.err

set -euo pipefail

source activate survival_clock

organs=(heart pancreas)
organ=${organs[$SLURM_ARRAY_TASK_ID]}

script="/cbica/home/wenju/Project/whole-body_clocks/mortality_clock/MRI/longitudinal/organ_mri_clock.py"

model_joblib="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/${organ}_mri_mortality_clock/${organ}_mri_mortality_clock_model.joblib"

input_tsv="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/longitudinal/imaging/data/imaging_${organ}_3_0.tsv"

outdir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/longitudinal/imaging/${organ}"
mkdir -p "${outdir}"

exec > "${outdir}/apply_${organ}_mri_mortality_clock_instance3_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out" \
     2> "${outdir}/apply_${organ}_mri_mortality_clock_instance3_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"

echo "============================================================"
echo "Applying pretrained MRI mortality clock to instance 3"
echo "Organ: ${organ}"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
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
  echo "ERROR: input instance-3 TSV not found: ${input_tsv}"
  exit 1
fi

pref="${organ}_mri_mortality_clock"

echo "Removing stale output files before rerun..."
rm -f "${outdir}/${pref}_apply_instance_3_0_predictions.tsv"
rm -f "${outdir}/${pref}_apply_instance_3_0_summary.json"
rm -f "${outdir}/${pref}_apply_longitudinal_instances_combined_predictions.tsv"
rm -f "${outdir}/${pref}_apply_longitudinal_instances_combined_summary.json"

python -u "${script}" \
  --organ "${organ}" \
  --model-joblib "${model_joblib}" \
  --input-tsv "3_0:${input_tsv}" \
  --covariate-csv /cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv \
  --death-xlsx /cbica/home/wenju/Dataset/UKBB_UMelbourne/Death_related_var_from_Ye.xlsx \
  --id-match-csv /cbica/home/wenju/Dataset/UKBB_UMelbourne/UKB_UMelbourne_vs_Penn_match_key.csv \
  --admin-censor-date 2022-11-30 \
  --outdir "${outdir}" \
  --application-session-id ses-M3 \
  --imaging-instance 3 \
  --model-instance 2 \
  --feature-start-column diagnosis \
  --risk-times 5,10,15 \
  --complete-case-organ-features

prediction_file="${outdir}/${pref}_apply_instance_3_0_predictions.tsv"

echo "============================================================"
echo "Checking final clean output"
echo "Prediction file: ${prediction_file}"
echo "============================================================"

if [ ! -f "${prediction_file}" ]; then
  echo "ERROR: prediction file was not generated: ${prediction_file}"
  exit 1
fi

required_cols=(
  "${organ}_mri_mortality_risk_score"
  "${organ}_mri_mortality_clock_acceleration_z"
  "${organ}_mri_mortality_clock_acceleration_years"
  "${organ}_mri_mortality_clock_age_years"
  "n_model_mri_features_expected"
  "n_model_mri_features_present_in_input"
  "n_model_mri_features_missing_from_input"
)

for col in "${required_cols[@]}"; do
  if ! head -n 1 "${prediction_file}" | tr '\t' '\n' | grep -Fxq "${col}"; then
    echo "ERROR: required column missing: ${col}"
    echo "Header columns containing ${organ}_mri_mortality:"
    head -n 1 "${prediction_file}" | tr '\t' '\n' | grep "${organ}_mri_mortality" || true
    exit 1
  fi
done

bad_pattern="${organ}_mri_mortality_clock_acceleration_${organ}_mri_mortality_clock_acceleration"
if head -n 1 "${prediction_file}" | grep -q "${bad_pattern}"; then
  echo "ERROR: malformed duplicated acceleration column detected."
  head -n 1 "${prediction_file}" | tr '\t' '\n' | grep "${organ}_mri_mortality_clock_acceleration" || true
  exit 1
fi

echo "Checking that original MRI feature columns are NOT included..."

if [ "${organ}" = "heart" ]; then
  if head -n 1 "${prediction_file}" | tr '\t' '\n' | grep -E "^(lv_|rv_|la_|ra_|ascending_aorta_|descending_aorta_)"; then
    echo "ERROR: heart MRI feature columns are still present in final prediction file."
    exit 1
  fi
fi

if [ "${organ}" = "pancreas" ]; then
  if head -n 1 "${prediction_file}" | tr '\t' '\n' | grep -E "^Pancreas_"; then
    echo "ERROR: pancreas MRI feature columns are still present in final prediction file."
    exit 1
  fi
fi

if head -n 1 "${prediction_file}" | tr '\t' '\n' | grep -Fxq "session_id"; then
  echo "ERROR: session_id should not be in final clean prediction file."
  exit 1
fi

if head -n 1 "${prediction_file}" | tr '\t' '\n' | grep -Fxq "diagnosis"; then
  echo "ERROR: diagnosis should not be in final clean prediction file."
  exit 1
fi

echo "Verified required columns:"
for col in "${required_cols[@]}"; do
  echo "  ${col}"
done

echo "Final clean header:"
head -n 1 "${prediction_file}" | tr '\t' '\n'

echo "============================================================"
echo "Finished applying pretrained MRI mortality clock"
echo "Organ: ${organ}"
echo "Finished at: $(date)"
echo "============================================================"

conda deactivate