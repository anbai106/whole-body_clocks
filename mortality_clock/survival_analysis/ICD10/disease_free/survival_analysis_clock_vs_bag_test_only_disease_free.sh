#!/bin/bash
set -euo pipefail

module load python/anaconda/3

if [[ "$#" -ne 2 ]]; then
  echo "Usage: $0 <disease_free_icd_tsv> <output_tsv>" >&2
  exit 2
fi

icd_tsv=$1
output_tsv=$2

source activate survival

PY_SCRIPT="/cbica/home/wenju/Project/whole-body_clocks/mortality_clock/survival_analysis/ICD10/disease_free/survival_analysis_clock_vs_bag_test_only_disease_free.py"
TEST_ROOT="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock"

echo "Start held-out test-only clock-vs-BAG survival comparison with disease-free CN controls"
python "${PY_SCRIPT}" \
  --icd_tsv "${icd_tsv}" \
  --output_tsv "${output_tsv}" \
  --test_predictions_root "${TEST_ROOT}" \
  --min_case 20 \
  --min_noncase 20

echo "Finish!"
conda deactivate
