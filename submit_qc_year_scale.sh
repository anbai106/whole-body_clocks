#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=QC
#SBATCH --mem-per-cpu=12G
#SBATCH --time=0-01:59:00
#SBATCH --output=/cbica/home/wenju/output/QC_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/QC_%A_%a.err

module load python/anaconda/3

source activate DNE

python /cbica/home/wenju/Project/Project/whole-body_clocks/qc_year_scale.py

conda deactivate