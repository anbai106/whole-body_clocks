#!/usr/bin/env bash
#SBATCH --job-name=apply_ad_epoch_ext_dual
#SBATCH --partition=all
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --output=/cbica/home/wenju/output/apply_ad_epoch_external_dual_%j.out
#SBATCH --error=/cbica/home/wenju/output/apply_ad_epoch_external_dual_%j.err

set -euo pipefail

# Apply the same pretrained ADNI brain MRI AD EPOCH model twice:
#   1. raw MUSE_Volume_* features, matching the original training feature type;
#   2. harmonized H_MUSE_Volume_* features, as a distribution-shift sensitivity analysis.
#
# Both analyses are run by default and are saved to separate output directories
# with distinct filename prefixes. The model is never refitted.

WORKDIR="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch"
LOGDIR="/cbica/home/wenju/output"

PY_SCRIPT="/cbica/home/wenju/Project/whole-body_clocks/ADNI/2_other_studies/ad_epoch_apply_external_longitudinal_dual.py"
EXTERNAL_FILE="${WORKDIR}/external_5_studies_istaging.tsv"
MODEL_JOBLIB="${WORKDIR}/results_brain_mri_ad_lepoch/adni_brain_mri_ad_lepoch_model.joblib"

RAW_OUTDIR="${WORKDIR}/results_external_longitudinal_ad_epoch_raw"
HARMONIZED_OUTDIR="${WORKDIR}/results_external_longitudinal_ad_epoch_harmonized"

MIN_ROI_FRACTION="${MIN_ROI_FRACTION:-0.80}"
AGE_MODE="${AGE_MODE:-row}"
RISK_TIMES="${RISK_TIMES:-1,2,3,5}"

# Default behavior is to run both modes. Override with, for example:
#   sbatch --export=ALL,RUN_MODES=raw submit_apply_ad_epoch_external_longitudinal_dual.sh
#   sbatch --export=ALL,RUN_MODES=harmonized submit_apply_ad_epoch_external_longitudinal_dual.sh
RUN_MODES="${RUN_MODES:-raw,harmonized}"

# Strict harmonized mode is the default: no raw fallback is used. Set to 1 to
# fill missing harmonized ROI values from raw MUSE values within the harmonized run.
HARMONIZED_ALLOW_RAW_FALLBACK="${HARMONIZED_ALLOW_RAW_FALLBACK:-0}"

mkdir -p "${LOGDIR}" "${RAW_OUTDIR}" "${HARMONIZED_OUTDIR}"

source activate survival_clock

python3 - <<'PY'
import importlib.util
required = ["pandas", "numpy", "sklearn", "sksurv", "joblib"]
missing = [name for name in required if importlib.util.find_spec(name) is None]
if missing:
    raise SystemExit("Missing packages: " + ", ".join(missing))
PY

for file in "${EXTERNAL_FILE}" "${MODEL_JOBLIB}" "${PY_SCRIPT}"; do
    if [[ ! -s "${file}" ]]; then
        echo "ERROR: missing or empty file: ${file}" >&2
        exit 1
    fi
done

run_application() {
    local roi_source="$1"
    local outdir="$2"
    local prefix="$3"
    local allow_raw_fallback="${4:-0}"

    local -a cmd=(
        python3 "${PY_SCRIPT}"
        --input-file "${EXTERNAL_FILE}"
        --model-joblib "${MODEL_JOBLIB}"
        --outdir "${outdir}"
        --prefix "${prefix}"
        --id-col "PTID"
        --visit-col "Visit_Code"
        --date-col "Date"
        --age-col "Age"
        --sex-col "Sex"
        --dlicv-col "DLICV"
        --site-col "SITE"
        --study-col "Study"
        --dx-col "DX_Binary"
        --delta-baseline-col "Delta_Baseline"
        --roi-source "${roi_source}"
        --min-roi-fraction "${MIN_ROI_FRACTION}"
        --age-mode "${AGE_MODE}"
        --risk-times "${RISK_TIMES}"
    )

    if [[ "${allow_raw_fallback}" == "1" ]]; then
        cmd+=(--allow-raw-roi-fallback)
    fi

    "${cmd[@]}"
}

IFS=',' read -r -a MODES <<< "${RUN_MODES}"

for mode in "${MODES[@]}"; do
    mode="$(echo "${mode}" | tr -d '[:space:]')"

    case "${mode}" in
        raw)
            run_application \
                "raw" \
                "${RAW_OUTDIR}" \
                "external_5_studies_adni_brain_mri_ad_epoch_raw"
            ;;

        harmonized)
            harmonized_extra=()
            if [[ "${HARMONIZED_ALLOW_RAW_FALLBACK}" == "1" ]]; then
                harmonized_extra+=(--allow-raw-roi-fallback)
            fi

            run_application \
                "harmonized" \
                "${HARMONIZED_OUTDIR}" \
                "external_5_studies_adni_brain_mri_ad_epoch_harmonized" \
                "${harmonized_extra[@]}"
            ;;

        "")
            ;;

        *)
            echo "ERROR: unsupported RUN_MODES entry: ${mode}" >&2
            echo "Allowed values: raw,harmonized" >&2
            exit 1
            ;;
    esac
done

echo "============================================================"
echo "All requested external AD EPOCH applications completed."
echo "RUN_MODES: ${RUN_MODES}"
echo "Raw outputs: ${RAW_OUTDIR}"
echo "Harmonized outputs: ${HARMONIZED_OUTDIR}"
echo "Harmonized raw fallback: ${HARMONIZED_ALLOW_RAW_FALLBACK}"
echo "============================================================"

conda deactivate
