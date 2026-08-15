#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=G309_EPOCH_vs_proteins
#SBATCH --mem-per-cpu=16G
#SBATCH --time=0-05:59:00
#SBATCH --output=/cbica/home/wenju/output/G309_EPOCH_vs_proteins_%j.out
#SBATCH --error=/cbica/home/wenju/output/G309_EPOCH_vs_proteins_%j.err

###############################################################################
# G309 incident-onset survival analysis:
# brain-proteomics mortality EPOCH vs each underlying brain-enriched protein.
#
# FULL BRAIN-PROTEOMICS SAMPLE:
#   1) training/training_4589.tsv
#   2) PT/patient_pop.tsv
#   3) test/ind_test_500.tsv
#
# One independent Cox model per predictor:
#   covariates + EPOCH_z
#   covariates + PMCH_z
#   covariates + MOG_z
#   ...
#
# There is NO joint EPOCH + protein model.
#
# Primary effect:
#   HR per 1-SD higher predictor.
#
# Primary comparison:
#   common complete-case sample across EPOCH + every protein + covariates.
###############################################################################

module load python/anaconda/3
source activate survival

mkdir -p /cbica/home/wenju/output

# -----------------------------------------------------------------------------
# Endpoint
# -----------------------------------------------------------------------------

ICD="G309"

ICD_TSV="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/SA/data/${ICD}_diagnosis_clock_disease_free.tsv"

# -----------------------------------------------------------------------------
# Brain-proteomics EPOCH
# -----------------------------------------------------------------------------

ORGAN="Brain"

EPOCH_TSV="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/Brain_proteomics_mortality_clock/brain_proteomics_mortality_clock_predictions.tsv"

EPOCH_COL="brain_proteomics_mortality_clock_acceleration_z"

# -----------------------------------------------------------------------------
# FULL underlying brain-enriched protein sample
# -----------------------------------------------------------------------------
#
# These three files are combined inside the Python script.
#
# If participant IDs were ever duplicated across source files, the Python
# script preferentially matches the organ table to the EPOCH prediction using:
#
#   participant_id + organ_source_file
#
# because brain_proteomics_mortality_clock_predictions.tsv contains
# organ_source_file.
# -----------------------------------------------------------------------------

ORGAN_TSV="/cbica/home/wenju/Reproducibile_paper/UKBB_Proteomics/MLNI/data/${ORGAN}/training/training_4589.tsv,/cbica/home/wenju/Reproducibile_paper/UKBB_Proteomics/MLNI/data/${ORGAN}/PT/patient_pop.tsv,/cbica/home/wenju/Reproducibile_paper/UKBB_Proteomics/MLNI/data/${ORGAN}/test/ind_test_500.tsv"

# -----------------------------------------------------------------------------
# Analysis settings
# -----------------------------------------------------------------------------

# Recommended for the forest plot:
# EPOCH and every protein use exactly the same complete-case participants.
SAMPLE_MODE="common"

# Use all EPOCH predictions because ORGAN_TSV contains training + PT + test.
EPOCH_SPLIT="all"

SESSION_ID="ses-M0"

MIN_CASE="20"
MIN_NONCASE="20"

# -----------------------------------------------------------------------------
# Output
# -----------------------------------------------------------------------------

OUT_DIR="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/SA/output_EPOCH_vs_underlying_proteins/G309"

OUTPUT_TSV="${OUT_DIR}/brain_proteomics_EPOCH_vs_underlying_proteins_G309.tsv"

mkdir -p "${OUT_DIR}"

# -----------------------------------------------------------------------------
# Python script
# -----------------------------------------------------------------------------

PY_SCRIPT="/cbica/home/wenju/Project/whole-body_clocks/mortality_clock/survival_analysis/ICD10/disease_free/comparative_analysis_with_underlying_features/survival_analysis_brain_proteomics_epoch_vs_proteins_G309.py"

# -----------------------------------------------------------------------------
# Input checks
# -----------------------------------------------------------------------------

if [[ ! -f "${ICD_TSV}" ]]; then
  echo "ERROR: missing G309 disease file:"
  echo "  ${ICD_TSV}"
  conda deactivate || true
  exit 1
fi

if [[ ! -f "${EPOCH_TSV}" ]]; then
  echo "ERROR: missing EPOCH prediction file:"
  echo "  ${EPOCH_TSV}"
  conda deactivate || true
  exit 1
fi

if [[ ! -f "${PY_SCRIPT}" ]]; then
  echo "ERROR: missing Python script:"
  echo "  ${PY_SCRIPT}"
  conda deactivate || true
  exit 1
fi

IFS=',' read -r -a ORGAN_FILES <<< "${ORGAN_TSV}"

for organ_file in "${ORGAN_FILES[@]}"; do
  if [[ ! -f "${organ_file}" ]]; then
    echo "ERROR: missing organ feature file:"
    echo "  ${organ_file}"
    conda deactivate || true
    exit 1
  fi
done

echo "======================================================================"
echo "G309: brain-proteomics EPOCH vs underlying proteins"
echo "======================================================================"
echo "ICD TSV:"
echo "  ${ICD_TSV}"
echo
echo "EPOCH TSV:"
echo "  ${EPOCH_TSV}"
echo
echo "EPOCH column:"
echo "  ${EPOCH_COL}"
echo
echo "Full organ feature sample:"
for organ_file in "${ORGAN_FILES[@]}"; do
  echo "  ${organ_file}"
done
echo
echo "Sample mode: ${SAMPLE_MODE}"
echo "EPOCH split: ${EPOCH_SPLIT}"
echo "Output:"
echo "  ${OUTPUT_TSV}"
echo "======================================================================"

if [[ -s "${OUTPUT_TSV}" ]]; then
  echo "Output already exists and is non-empty; skipping:"
  echo "  ${OUTPUT_TSV}"
  conda deactivate
  exit 0
fi

python "${PY_SCRIPT}" \
  --icd-tsv "${ICD_TSV}" \
  --epoch-tsv "${EPOCH_TSV}" \
  --organ-tsv "${ORGAN_TSV}" \
  --output-tsv "${OUTPUT_TSV}" \
  --epoch-col "${EPOCH_COL}" \
  --sample-mode "${SAMPLE_MODE}" \
  --split "${EPOCH_SPLIT}" \
  --session-id "${SESSION_ID}" \
  --min-case "${MIN_CASE}" \
  --min-noncase "${MIN_NONCASE}"

echo "======================================================================"
echo "Finished G309 survival analysis."
echo "Output:"
echo "  ${OUTPUT_TSV}"
echo "======================================================================"

conda deactivate