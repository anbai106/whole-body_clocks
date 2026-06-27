#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=ArrayJob
#SBATCH --array=0-4
#SBATCH --time=1-01:59:00
#SBATCH --mem-per-cpu=24G
#SBATCH --output=/cbica/home/wenju/output/MLNI_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/MLNI_%A_%a.err

numbers=(Endocrine Digestive Hepatic Immune Metabolic)
organ=${numbers[$SLURM_ARRAY_TASK_ID]}

echo ${organ}
zip /cbica/home/wenju/Reproducibile_paper/UKBB_metabolomics/fastGWA_MetBAG/output/${organ}/organ_pheno_normalized_residualized.fastGWA.zip /cbica/home/wenju/Reproducibile_paper/UKBB_metabolomics/fastGWA_MetBAG/output/${organ}/organ_pheno_normalized_residualized.fastGWA