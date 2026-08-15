#!/bin/bash
#SBATCH --partition=bioinformatics
#SBATCH --job-name=select
#SBATCH --mem-per-cpu=24G
#SBATCH --time=12:00:00
#SBATCH --output=/cbica/home/wenju/output/select_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/select_%A_%a.err


module load python/anaconda/3
source activate r_3.6.3

output_dir=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/comparative_with_10_conventional_biomarkers/1_prepare_data
mkdir -p ${output_dir}

Rscript /cbica/home/wenju/Project/whole-body_clocks/mortality_prediction_disease_clock/comparative_with_10_conventional_biomarkers/SelectVarLocalUKBBMRFullUKBBSample.R
conda deactivate