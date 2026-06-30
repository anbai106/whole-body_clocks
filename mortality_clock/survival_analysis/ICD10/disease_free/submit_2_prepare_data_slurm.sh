#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=ArrayJob
#SBATCH --array=0-3572
#SBATCH --mem-per-cpu=24G
#SBATCH --time=0-05:59:00
#SBATCH --output=/cbica/home/wenju/output/SA_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/SA_%A_%a.err

### MetBAG
icd_array=( $( awk '{print $1}' /cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/SA/data/included_ICD_mortality_clock.tsv | grep -v "value" | xargs) )
disease=${icd_array[$SLURM_ARRAY_TASK_ID]}

source activate survival
python /cbica/home/wenju/Project/whole-body_clocks/mortality_clock/survival_analysis/ICD10/disease_free/2_prepare_population.py --disease ${disease}
conda deactivate