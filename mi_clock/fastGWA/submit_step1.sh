#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=prepare_data
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --output=/cbica/home/wenju/output/prepare_data_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/prepare_data_%A_%a.err

source activate DNE
python /cbica/home/wenju/Project/whole-body_clocks/mi_clock/fastGWA/1_prepare_data.py

echo "Finish!"
conda deactivate
