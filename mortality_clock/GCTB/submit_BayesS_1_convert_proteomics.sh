#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=ArrayJob
#SBATCH --array=0-10
#SBATCH --mem-per-cpu=12G
#SBATCH --time=0-01:59:00
#SBATCH --output=/cbica/home/wenju/output/SBayesS_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/SBayesS_%A_%a.err

numbers=(Reproductive_female Pulmonary Heart Brain Eye Hepatic Renal Reproductive_male Endocrine Immune Skin)
organ=${numbers[$SLURM_ARRAY_TASK_ID]}

output_dir=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/GCTB/${organ}_proteomics_mortality_clock
mkdir -p $output_dir

source_file=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/fastGWA/output/${organ}_proteomics_mortality_clock/organ_pheno_normalized_residualized.fastGWA
target_file="${output_dir}/${organ}_gctb.tsv"
if [ ! -f ${target_file} ]; then
  echo "Generate GCTB .../ for ${organ}..."
  bash /cbica/home/wenju/Project/whole-body_clocks/mortality_clock/GCTB/BayesS_1_convert.sh $source_file $target_file
else
  :
fi






