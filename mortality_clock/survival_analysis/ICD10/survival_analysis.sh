#!/bin/bash
module load python/anaconda/3

icd_tsv=$1
output_tsv=$2

source activate survival
echo "Start training"
python /cbica/home/wenju/Project/FemaleAgingClock/survival_analysis/ICD10/full_covariate/female/survival_analysis.py --icd_tsv ${icd_tsv} --output_tsv ${output_tsv}
echo "Finish!"
conda deactivate
