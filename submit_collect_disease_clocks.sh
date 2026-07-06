#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=QC
#SBATCH --mem-per-cpu=12G
#SBATCH --time=0-5:59:00
#SBATCH --output=/cbica/home/wenju/output/collect_disease_clocks_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/collect_disease_clocks_%A_%a.err

module load python/anaconda/3

source activate DNE

python /cbica/home/wenju/Project/whole-body_clocks/collect_disease_clocks_stable_significant.py

conda deactivate