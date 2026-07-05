#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=med_step1
#SBATCH --mem-per-cpu=16G
#SBATCH --time=0-06:00:00
#SBATCH --output=/cbica/home/wenju/output/step1_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/step1_%A_%a.err

set -euo pipefail

source activate survival_clock

python /cbica/home/wenju/Project/whole-body_clocks/mortality_clock/Metabolomics/longitudinal/medication/1_prepare_medication_cluster_baseline.py

conda deactivate