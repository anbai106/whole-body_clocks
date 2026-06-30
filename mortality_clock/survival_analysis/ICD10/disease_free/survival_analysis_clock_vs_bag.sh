#!/bin/bash
set -euo pipefail

module load python/anaconda/3

icd_tsv=$1
output_tsv=$2

source activate survival

echo "Start clock-vs-BAG disease survival comparison"
python /cbica/home/wenju/Project/whole-body_clocks/mortality_clock/survival_analysis/ICD10/disease_free/survival_analysis_clock_vs_bag.py \
  --icd_tsv "${icd_tsv}" \
  --output_tsv "${output_tsv}"
echo "Finish!"

conda deactivate
