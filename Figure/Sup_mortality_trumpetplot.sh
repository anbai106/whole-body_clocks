#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=collect
#SBATCH --mem-per-cpu=24G
#SBATCH --output=/cbica/home/wenju/output/collect%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/collect%A_%a.err

############################## END OF DEFAULT EMBEDDED SGE COMMANDS #######################

source activate MLNI_012
echo "Start applying..."
python /cbica/home/wenju/Project/whole-body_clocks/Figure/Sup_mortality_trumpetplot.py
echo "Finish!"
conda deactivate
