#!/bin/bash
#SBATCH --partition=bioinformatics
#SBATCH --job-name=SEM_EPOCH
#SBATCH --mem=24G
#SBATCH --time=0-04:59:00
#SBATCH --output=/cbica/home/wenju/output/SEM_%j.out
#SBATCH --error=/cbica/home/wenju/output/SEM_%j.err

module load R/4.3

echo "Start: $(date)"
echo "Host: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID:-NA}"

Rscript /cbica/home/wenju/Project/whole-body_clocks/UKBB_disease_mortality_epoch/mediation_anlayses/mediation.R

echo "Finish: $(date)"
