#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=SEM
#SBATCH --mem-per-cpu=24G
#SBATCH --time=0-04:59:00
#SBATCH --output=/cbica/home/wenju/output/SEM_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/SEM_%A_%a.err

module load R/4.3
echo "Start training"
Rscript /cbica/home/wenju/Project/whole-body_clocks/UKBB_disease_mortality_epoch/mediation_anlayses/mediation.R
echo "Finish!"
