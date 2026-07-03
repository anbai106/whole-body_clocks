#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=apply_mri_clock_i3
#SBATCH --mem-per-cpu=24G
#SBATCH --array=0-1
#SBATCH --output=/cbica/home/wenju/output/apply_mri_clock_instance3_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/apply_mri_clock_instance3_%A_%a.err

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

python "${script}" \
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
  --complete-case-organ-features \
  --include-features-in-output

echo "============================================================"
echo "Finished applying pretrained MRI mortality clock"
echo "Organ: ${organ}"
echo "Finished at: $(date)"
echo "============================================================"

conda deactivate