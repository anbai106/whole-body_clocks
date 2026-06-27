#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=ArrayJob
#SBATCH --array=0-4
#SBATCH --mem-per-cpu=12G
#SBATCH --time=0-01:59:00
#SBATCH --output=/cbica/home/wenju/output/LDSC_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/LDSC_%A_%a.err

numbers=(Endocrine Digestive Hepatic Immune Metabolic)
organ=${numbers[$SLURM_ARRAY_TASK_ID]}

file=/cbica/home/wenju/Reproducibile_paper/UKBB_metabolomics/fastGWA_MetBAG/output/${organ}/organ_pheno_normalized_residualized.fastGWA.ldsc.tsv
output_file=/cbica/home/wenju/Reproducibile_paper/UKBB_metabolomics/gc/${organ}/${organ}
mkdir -p /cbica/home/wenju/Reproducibile_paper/UKBB_metabolomics/gc/${organ}
bash /cbica/home/wenju/Project/UKBB_NMR_metabolomics/LDSC/ldsc_munge_stat.sh ${file} ${output_file}
