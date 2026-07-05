#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=model_performance
#SBATCH --mem-per-cpu=12G
#SBATCH --time=0-9:59:00
#SBATCH --output=/cbica/home/wenju/output/model_performance_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/model_performance_%A_%a.err

module load R/4.3

Rscript /cbica/home/wenju/Project/whole-body_clocks/mortality_clock/3_mortality_clock_performance_vis_all_organs.R
