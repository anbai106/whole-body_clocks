#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=ArrayJob
#SBATCH --array=0-4
#SBATCH --mem-per-cpu=12G
#SBATCH --time=0-05:59:00
#SBATCH --output=/cbica/home/wenju/output/LDSC_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/LDSC_%A_%a.err

numbers=(Endocrine Digestive Hepatic Immune Metabolic)
organ=${numbers[$SLURM_ARRAY_TASK_ID]}

stat_gz1=/cbica/home/wenju/Reproducibile_paper/UKBB_metabolomics/gc/${organ}/${organ}.sumstats.gz
output_file=/cbica/home/wenju/Reproducibile_paper/UKBB_metabolomics/gc/${organ}/${organ}_vs
bash /cbica/home/wenju/Project/UKBB_NMR_metabolomics/LDSC/ldsc_main_gc.sh ${stat_gz1} ${output_file}

