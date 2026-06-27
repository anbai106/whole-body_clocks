#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=ArrayJob
#SBATCH --array=0-1998
#SBATCH --mem-per-cpu=12G
#SBATCH --time=0-05:59:00
#SBATCH --output=/cbica/home/wenju/output/SA_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/SA_%A_%a.err

############################## END OF DEFAULT EMBEDDED SGE COMMANDS #######################

### MetBAG
icd_array=( $( awk '{print $1}' /cbica/home/wenju/Reproducibile_paper/FemaleAgingClock/survival_analysis/data/included_ICD_disease_femalebag.tsv | grep -v "value" | xargs) )
icd=${icd_array[$SLURM_ARRAY_TASK_ID]}

echo "Run SA to derive HR for ICD: ${icd}..."
icd_tsv=/cbica/home/wenju/Reproducibile_paper/FemaleAgingClock/survival_analysis/data/${icd}_diagnosis_female.tsv
mkdir -p /cbica/home/wenju/Reproducibile_paper/FemaleAgingClock/survival_analysis/output/full_cov
output_tsv=/cbica/home/wenju/Reproducibile_paper/FemaleAgingClock/survival_analysis/output/full_cov/cox_hr_${icd}_female.tsv
if [[ ! -f ${output_tsv} ]]; then
  bash /cbica/home/wenju/Project/FemaleAgingClock/survival_analysis/ICD10/full_covariate/female/survival_analysis.sh ${icd_tsv} ${output_tsv}
fi


