#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=SAnalysis
#SBATCH --mem-per-cpu=24G
#SBATCH --output=/cbica/home/wenju/output/ICD_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/ICD_%A_%a.err

############################## END OF DEFAULT EMBEDDED SGE COMMANDS #######################
### brain
source activate survival
python /cbica/home/wenju/Project/whole-body_clocks/mortality_clock/survival_analysis/ICD10/1_prepare_population.py
conda deactivate