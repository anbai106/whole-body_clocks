#!/bin/bash
set -euo pipefail

module load python/anaconda/3

endpoint=$1
output_tsv=$2

delta_root=${DELTA_ROOT:-/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/longitudinal/metabolomics/metabolomics_delta_acceleration_years_landmark_survival_analysis}
umel_death_xlsx=${UMEL_DEATH_XLSX:-/cbica/home/wenju/Dataset/UKBB_UMelbourne/Death_related_var_from_Ye.xlsx}
umel_match_csv=${UMEL_MATCH_CSV:-/cbica/home/wenju/Dataset/UKBB_UMelbourne/UKB_UMelbourne_vs_Penn_match_key.csv}
cov_tsv=${COV_TSV:-/cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv}

# Default follows the column specified by the user.
# To use pure acceleration-year change, run with:
#   DELTA_COLUMN=delta_accel_years_1_minus_0 bash this_script.sh <endpoint> <output.tsv>
delta_column=${DELTA_COLUMN:-delta_clock_age_1_minus_0}

python_script=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/longitudinal/metabolomics/metabolomics_delta_acceleration_years_landmark_survival_analysis/survival_analysis_delta_metabolomics_disease_onset.py

source activate survival

echo "============================================================"
echo "Delta metabolomics disease-onset survival analysis"
echo "Endpoint: ${endpoint}"
echo "Delta root: ${delta_root}"
echo "Delta column: ${delta_column}"
echo "Output TSV: ${output_tsv}"
echo "Python script: ${python_script}"
echo "Started at: $(date)"
echo "============================================================"

python "${python_script}" \
  --endpoint "${endpoint}" \
  --output_tsv "${output_tsv}" \
  --delta_root "${delta_root}" \
  --delta_column "${delta_column}" \
  --umel_death_xlsx "${umel_death_xlsx}" \
  --umel_match_csv "${umel_match_csv}" \
  --cov_tsv "${cov_tsv}" \
  --include_bp

echo "============================================================"
echo "Finished endpoint: ${endpoint}"
echo "Finished at: $(date)"
echo "============================================================"

conda deactivate
