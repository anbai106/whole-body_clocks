#!/usr/bin/env bash
#SBATCH --job-name=apply_ad_epoch_mci2ad
#SBATCH --output=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch/logs/apply_ad_epoch_mci2ad_%j.out
#SBATCH --error=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch/logs/apply_ad_epoch_mci2ad_%j.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G

# Apply the frozen ADNI brain MRI AD L'EPOCH model to participants who:
#   1. were not used to fit the original model;
#   2. have a first qualifying MCI diagnosis with usable MUSE MRI data; and
#   3. have follow-up for AD conversion or censoring.
#
# Final model fitting in ad_epoch.py used the train + validation splits.
# Therefore, this application excludes split == train or validation from the
# original training prediction file. Test participants are retained because
# they were not used to fit the final Cox model.

WORKDIR="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch"
LOGDIR="${WORKDIR}/logs"
OUTDIR="${WORKDIR}/results_brain_mri_ad_lepoch_mci2ad_application"

mkdir -p "${WORKDIR}" "${LOGDIR}" "${OUTDIR}"

PY_SCRIPT="/cbica/home/wenju/Project/whole-body_clocks/ADNI/ad_epoch_apply_mci2ad.py"
ADNI_FILE="${WORKDIR}/adni_istaging.tsv"
MODEL_JOBLIB="${WORKDIR}/results_brain_mri_ad_lepoch/adni_brain_mri_ad_lepoch_model.joblib"
TRAINING_PARTICIPANTS_FILE="${WORKDIR}/results_brain_mri_ad_lepoch/adni_brain_mri_ad_lepoch_predictions.tsv"

# Activate the environment robustly in a non-interactive SLURM shell.
if command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)"
fi
conda activate survival_clock

for required_file in "${PY_SCRIPT}" "${ADNI_FILE}" "${MODEL_JOBLIB}" "${TRAINING_PARTICIPANTS_FILE}"; do
  if [[ ! -s "${required_file}" ]]; then
    echo "ERROR: required file is missing or empty: ${required_file}" >&2
    exit 1
  fi
done

python3 "${PY_SCRIPT}" \
  --input-file "${ADNI_FILE}" \
  --model-joblib "${MODEL_JOBLIB}" \
  --training-participants-file "${TRAINING_PARTICIPANTS_FILE}" \
  --outdir "${OUTDIR}" \
  --prefix "adni_brain_mri_ad_lepoch" \
  --id-col "PTID" \
  --visit-col "Visit_Code" \
  --date-col "Date" \
  --dx-col "DX_Binary" \
  --training-id-col "PTID" \
  --training-split-col "split" \
  --exclude-splits "train,validation" \
  --baseline-dx "MCI" \
  --event-dx "AD" \
  --min-baseline-roi-fraction 0.80 \
  --min-followup-days 1 \
  --risk-times "1,2,3,5"

echo "Finished applying AD L'EPOCH to non-training baseline-MCI participants."
echo "Primary output: ${OUTDIR}/adni_brain_mri_ad_lepoch_mci2ad_epoch_survival.tsv"

conda deactivate
