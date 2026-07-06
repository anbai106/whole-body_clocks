#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=disease_clock_mortality
#SBATCH --mem-per-cpu=24G
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --array=0-46
#SBATCH --output=/cbica/home/wenju/output/disease_clock_mortality_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/disease_clock_mortality_%A_%a.err

# ============================================================
# Submit/run mortality survival analysis for 47 stable and
# significant disease clocks.
#
# Important:
#   - MRI disease clocks use UKB Field 53 instance 2_0.
#   - Proteomics/metabolomics disease clocks use UKB Field 53 instance 0_0.
#   - Field 53 is read from the Melbourne death-related file.
#   - Covariates are read from the full-sample covariate file.
#   - Output goes to:
#       <BASE_DIR>/<clock_folder>/survival_analysis_mortality/
#
# Recommended:
#   bash submit_47_disease_clock_mortality_survival.slurm
# ============================================================

set -euo pipefail

BASE_DIR="${BASE_DIR:-/gpfs/fs001/cbica/home/wenju/Reproducibile_paper/WholeBodyClock}"

ANALYSIS_DIR="${ANALYSIS_DIR:-${BASE_DIR}/all_disease_lepoch_incremental_value_scale_qc}"

SCORE_WIDE_TSV="${SCORE_WIDE_TSV:-${ANALYSIS_DIR}/stable_significant_disease_clock_acceleration_z_wide.tsv}"

SCORE_METADATA_TSV="${SCORE_METADATA_TSV:-${ANALYSIS_DIR}/stable_significant_disease_clock_acceleration_z_metadata.tsv}"

GOOD_CLOCK_TSV="${GOOD_CLOCK_TSV:-${ANALYSIS_DIR}/all_disease_lepoch_main_text_good_clocks.tsv}"

TASKS_TSV="${TASKS_TSV:-${ANALYSIS_DIR}/mortality_survival_47_stable_significant_clock_tasks.tsv}"

PY_SCRIPT=/cbica/home/wenju/Project/whole-body_clocks/mortality_prediction_disease_clock/mortality_prediction.py

DEATH_XLSX="${DEATH_XLSX:-/cbica/home/wenju/Dataset/UKBB_UMelbourne/Death_related_var_from_Ye.xlsx}"

ID_MATCH_CSV="${ID_MATCH_CSV:-/cbica/home/wenju/Dataset/UKBB_UMelbourne/UKB_UMelbourne_vs_Penn_match_key.csv}"

COVARIATE_CSV="${COVARIATE_CSV:-/cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv}"

ADMIN_CENSOR_DATE="${ADMIN_CENSOR_DATE:-2022-11-30}"

SLURM_OUTPUT_DIR="${SLURM_OUTPUT_DIR:-/cbica/home/wenju/output}"

CONDA_ENV="${CONDA_ENV:-survival_clock}"

# Optional exact column overrides. Leave empty for auto-detection.
FIELD53_0_COL="${FIELD53_0_COL:-}"
FIELD53_2_COL="${FIELD53_2_COL:-}"
DEATH_DATE_COL="${DEATH_DATE_COL:-}"
DEATH_ID_COL="${DEATH_ID_COL:-}"
IDMATCH_SCORE_COL="${IDMATCH_SCORE_COL:-}"
IDMATCH_DEATH_COL="${IDMATCH_DEATH_COL:-}"
COVARIATE_ID_COL="${COVARIATE_ID_COL:-}"

# Optional exact covariates.
# Example:
# COVARIATE_COLS="age,sex,assessment_center,genotype_array,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
COVARIATE_COLS="${COVARIATE_COLS:-}"

PENALIZER="${PENALIZER:-0.01}"

mkdir -p "${SLURM_OUTPUT_DIR}"
mkdir -p "${ANALYSIS_DIR}"

make_tasks_file() {
  python - "${SCORE_METADATA_TSV}" "${GOOD_CLOCK_TSV}" "${SCORE_WIDE_TSV}" "${TASKS_TSV}" "${BASE_DIR}" <<'PY'
import sys
import re
from pathlib import Path
import pandas as pd

metadata_tsv = Path(sys.argv[1])
good_tsv = Path(sys.argv[2])
wide_tsv = Path(sys.argv[3])
tasks_tsv = Path(sys.argv[4])
base_dir = Path(sys.argv[5])

def clean_col(x):
    x = str(x).strip().lower()
    x = re.sub(r"[^a-z0-9]+", "_", x)
    x = re.sub(r"_+", "_", x)
    x = x.strip("_")
    return x

def infer_modality(folder, modality):
    text = "{} {}".format(folder, modality).lower()
    if "mri" in text:
        return "MRI"
    if "proteomics" in text:
        return "Proteomics"
    if "metabolomics" in text:
        return "Metabolomics"
    return str(modality)

def make_score_col(disease, clock_label):
    return "{}_{}_clock_acceleration_z".format(
        clean_col(disease),
        clean_col(clock_label)
    )

if not metadata_tsv.exists():
    raise FileNotFoundError("Missing metadata TSV: {}".format(metadata_tsv))

if not wide_tsv.exists():
    raise FileNotFoundError("Missing wide score TSV: {}".format(wide_tsv))

meta = pd.read_csv(metadata_tsv, sep="\t")
wide_cols = list(pd.read_csv(wide_tsv, sep="\t", nrows=0).columns)

if "participant_id" not in wide_cols:
    raise ValueError("participant_id is not in the wide score TSV.")

if "status" in meta.columns:
    meta = meta[meta["status"] == "collected"].copy()

if good_tsv.exists():
    good = pd.read_csv(good_tsv, sep="\t")
    needed = {"disease", "folder"}
    if needed.issubset(set(good.columns)):
        good_key = good[["disease", "folder"]].drop_duplicates()
        good_key["disease"] = good_key["disease"].astype(str).str.lower()
        good_key["folder"] = good_key["folder"].astype(str)

        meta["disease"] = meta["disease"].astype(str).str.lower()
        meta["folder"] = meta["folder"].astype(str)

        meta = meta.merge(good_key, on=["disease", "folder"], how="inner")

required = ["disease", "folder"]
for col in required:
    if col not in meta.columns:
        raise ValueError("Required column missing from metadata: {}".format(col))

rows = []

for _, r in meta.drop_duplicates(subset=["disease", "folder"]).iterrows():
    disease = str(r["disease"]).lower()
    folder = str(r["folder"])

    clock_label = str(r["clock_label"]) if "clock_label" in meta.columns else folder
    clock_id = str(r["clock_id"]) if "clock_id" in meta.columns else ""
    modality = str(r["modality"]) if "modality" in meta.columns else infer_modality(folder, "")
    organ_label = str(r["organ_label"]) if "organ_label" in meta.columns else ""

    score_col = None

    if "score_col_wide" in meta.columns and pd.notna(r["score_col_wide"]):
        candidate = str(r["score_col_wide"])
        candidate_single = candidate.replace("__", "_")

        if candidate in wide_cols:
            score_col = candidate
        elif candidate_single in wide_cols:
            score_col = candidate_single

    if score_col is None:
        candidate = make_score_col(disease, clock_label)
        if candidate in wide_cols:
            score_col = candidate

    if score_col is None:
        print(
            "WARNING: could not find score column for disease={} folder={} clock_label={}".format(
                disease,
                folder,
                clock_label
            ),
            file=sys.stderr,
        )
        continue

    outdir = base_dir / folder / "survival_analysis_mortality"
    outdir.mkdir(parents=True, exist_ok=True)

    rows.append({
        "array_id": len(rows),
        "disease": disease,
        "folder": folder,
        "clock_label": clock_label,
        "clock_id": clock_id,
        "modality": infer_modality(folder, modality),
        "organ_label": organ_label,
        "score_col_wide": score_col,
        "output_dir": str(outdir),
    })

tasks = pd.DataFrame(rows)

if tasks.empty:
    raise ValueError("No stable/significant clock tasks were created.")

tasks.to_csv(tasks_tsv, sep="\t", index=False)

print("Wrote task file:", tasks_tsv)
print("Number of tasks:", tasks.shape[0])

if tasks.shape[0] != 47:
    print(
        "WARNING: Expected 47 stable/significant disease clocks, but task file has {}.".format(tasks.shape[0]),
        file=sys.stderr,
    )

print(tasks[["array_id", "disease", "modality", "folder", "score_col_wide"]].head(10).to_string(index=False))
PY
}

# ============================================================
# Self-submit mode
# ============================================================

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "============================================================"
  echo "Preparing and submitting 47-clock mortality Slurm array"
  echo "============================================================"
  echo "BASE_DIR:           ${BASE_DIR}"
  echo "ANALYSIS_DIR:       ${ANALYSIS_DIR}"
  echo "SCORE_WIDE_TSV:     ${SCORE_WIDE_TSV}"
  echo "SCORE_METADATA_TSV: ${SCORE_METADATA_TSV}"
  echo "GOOD_CLOCK_TSV:     ${GOOD_CLOCK_TSV}"
  echo "TASKS_TSV:          ${TASKS_TSV}"
  echo "PY_SCRIPT:          ${PY_SCRIPT}"
  echo "DEATH_XLSX:         ${DEATH_XLSX}"
  echo "ID_MATCH_CSV:       ${ID_MATCH_CSV}"
  echo "COVARIATE_CSV:      ${COVARIATE_CSV}"
  echo "ADMIN_CENSOR_DATE:  ${ADMIN_CENSOR_DATE}"
  echo "SLURM_OUTPUT_DIR:   ${SLURM_OUTPUT_DIR}"
  echo "============================================================"

  if [[ ! -f "${SCORE_WIDE_TSV}" ]]; then
    echo "ERROR: Missing SCORE_WIDE_TSV: ${SCORE_WIDE_TSV}" >&2
    exit 1
  fi

  if [[ ! -f "${SCORE_METADATA_TSV}" ]]; then
    echo "ERROR: Missing SCORE_METADATA_TSV: ${SCORE_METADATA_TSV}" >&2
    exit 1
  fi

  if [[ ! -f "${PY_SCRIPT}" ]]; then
    echo "ERROR: Missing PY_SCRIPT: ${PY_SCRIPT}" >&2
    exit 1
  fi

  make_tasks_file

  N_TASKS=$(($(wc -l < "${TASKS_TSV}") - 1))

  if [[ "${N_TASKS}" -lt 1 ]]; then
    echo "ERROR: No tasks created." >&2
    exit 1
  fi

  echo "Submitting ${N_TASKS} array jobs."

  sbatch --array=0-$((N_TASKS - 1))%12 "$0"

  exit 0
fi

# ============================================================
# Slurm task mode
# ============================================================

echo "============================================================"
echo "Running disease-clock mortality survival task"
echo "============================================================"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "Host: $(hostname)"
echo "Date: $(date)"
echo "============================================================"

if [[ ! -f "${TASKS_TSV}" ]]; then
  if [[ "${SLURM_ARRAY_TASK_ID}" == "0" ]]; then
    echo "Task file missing. Task 0 will create it."
    make_tasks_file
  else
    echo "Waiting for task file to be created by task 0..."
    for i in $(seq 1 120); do
      if [[ -f "${TASKS_TSV}" ]]; then
        break
      fi
      sleep 5
    done
  fi
fi

if [[ ! -f "${TASKS_TSV}" ]]; then
  echo "ERROR: Task file still missing: ${TASKS_TSV}" >&2
  exit 1
fi

source ~/.bashrc || true
conda activate "${CONDA_ENV}" || true

# ------------------------------------------------------------
# Critical fix:
# Declare the array, then use safe expansion later:
#   ${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"}
# This avoids the Bash set -u error:
#   EXTRA_ARGS[@]: unbound variable
# ------------------------------------------------------------
declare -a EXTRA_ARGS=()

if [[ -n "${FIELD53_0_COL}" ]]; then
  EXTRA_ARGS+=(--field53-0-col "${FIELD53_0_COL}")
fi

if [[ -n "${FIELD53_2_COL}" ]]; then
  EXTRA_ARGS+=(--field53-2-col "${FIELD53_2_COL}")
fi

if [[ -n "${DEATH_DATE_COL}" ]]; then
  EXTRA_ARGS+=(--death-date-col "${DEATH_DATE_COL}")
fi

if [[ -n "${DEATH_ID_COL}" ]]; then
  EXTRA_ARGS+=(--death-id-col "${DEATH_ID_COL}")
fi

if [[ -n "${IDMATCH_SCORE_COL}" ]]; then
  EXTRA_ARGS+=(--idmatch-score-col "${IDMATCH_SCORE_COL}")
fi

if [[ -n "${IDMATCH_DEATH_COL}" ]]; then
  EXTRA_ARGS+=(--idmatch-death-col "${IDMATCH_DEATH_COL}")
fi

if [[ -n "${COVARIATE_ID_COL}" ]]; then
  EXTRA_ARGS+=(--covariate-id-col "${COVARIATE_ID_COL}")
fi

if [[ -n "${COVARIATE_COLS}" ]]; then
  EXTRA_ARGS+=(--covariate-cols "${COVARIATE_COLS}")
fi

python "${PY_SCRIPT}" \
  --base-dir "${BASE_DIR}" \
  --score-wide-tsv "${SCORE_WIDE_TSV}" \
  --tasks-tsv "${TASKS_TSV}" \
  --task-index "${SLURM_ARRAY_TASK_ID}" \
  --death-xlsx "${DEATH_XLSX}" \
  --id-match-csv "${ID_MATCH_CSV}" \
  --covariate-csv "${COVARIATE_CSV}" \
  --admin-censor-date "${ADMIN_CENSOR_DATE}" \
  --penalizer "${PENALIZER}" \
  ${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"}

echo "============================================================"
echo "Finished task ${SLURM_ARRAY_TASK_ID}"
echo "Date: $(date)"
echo "============================================================"