#!/usr/bin/env bash
#SBATCH --job-name=apply_ad_epoch_external
#SBATCH --partition=all
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=36G
#SBATCH --output=/cbica/home/wenju/output/apply_ad_epoch_external_%j.out
#SBATCH --error=/cbica/home/wenju/output/apply_ad_epoch_external_%j.err

set -euo pipefail

# Apply the pretrained ADNI brain MRI AD EPOCH to every eligible longitudinal
# scan in the external iSTAGING dataset. This does not refit the model.

WORKDIR="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch"
OUTDIR="${WORKDIR}/results_external_longitudinal_ad_epoch"
LOGDIR="/cbica/home/wenju/output"

PY_SCRIPT="/cbica/home/wenju/Project/whole-body_clocks/ADNI/2_other_studies/ad_epoch_apply_external_longitudinal.py"

# Update this filename to the actual external five-study iSTAGING TSV.
EXTERNAL_FILE="${WORKDIR}/external_5_studies_istaging.tsv"

MODEL_JOBLIB="${WORKDIR}/results_brain_mri_ad_lepoch/adni_brain_mri_ad_lepoch_model.joblib"

mkdir -p "${OUTDIR}" "${LOGDIR}"

source activate survival_clock

python3 - <<'PY'
import importlib.util
required = ["pandas", "numpy", "sklearn", "sksurv", "joblib"]
missing = [x for x in required if importlib.util.find_spec(x) is None]
if missing:
    raise SystemExit("Missing packages: " + ", ".join(missing))
PY

for f in "${EXTERNAL_FILE}" "${MODEL_JOBLIB}" "${PY_SCRIPT}"; do
    if [[ ! -s "${f}" ]]; then
        echo "ERROR: missing or empty file: ${f}" >&2
        exit 1
    fi
done

python3 "${PY_SCRIPT}" \
  --input-file "${EXTERNAL_FILE}" \
  --model-joblib "${MODEL_JOBLIB}" \
  --outdir "${OUTDIR}" \
  --prefix "external_5_studies_adni_brain_mri_ad_epoch" \
  --id-col "PTID" \
  --visit-col "Visit_Code" \
  --date-col "Date" \
  --age-col "Age" \
  --sex-col "Sex" \
  --dlicv-col "DLICV" \
  --site-col "SITE" \
  --study-col "Study" \
  --dx-col "DX_Binary" \
  --delta-baseline-col "Delta_Baseline" \
  --min-roi-fraction 0.80 \
  --age-mode "row" \
  --risk-times "1,2,3,5"

echo "Finished external longitudinal AD EPOCH application."
echo "Outputs: ${OUTDIR}"

conda deactivate
