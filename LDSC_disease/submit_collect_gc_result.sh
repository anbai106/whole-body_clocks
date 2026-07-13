#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=LDSC_collect
#SBATCH --mem-per-cpu=12G
#SBATCH --time=0-09:59:00
#SBATCH --output=/cbica/home/wenju/output/LDSC_collect_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/LDSC_collect_%A_%a.err

module load python/anaconda/3
source activate DNE
bash /cbica/home/wenju/Project/whole-body_clocks/LDSC_disease/collect_gc_results.py
conda deactivate