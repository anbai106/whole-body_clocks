#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=Arrayjob
#SBATCH --array=0-1
#SBATCH --time=0-02:59:00
#SBATCH --mem-per-cpu=24G
#SBATCH --output=/cbica/home/wenju/output/qmplot_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/qmplot_%A_%a.err

organ=Metabolic
sexes=(female male)
sex=${sexes[$SLURM_ARRAY_TASK_ID]}

output_dir="/cbica/home/wenju/Reproducibile_paper/UKBB_metabolomics/fastGWA_MetBAG/output/${organ}/${sex}"
output_result="${output_dir}/BAG_pheno_normalized_residualized.fastGWA"
if [[ -f ${output_result} ]]; then
  result_png="${output_dir}/QQ_plot.png"
  if [[ ! -f ${result_png} ]]; then
    echo "qmplot for ${organ}..."
    bash /cbica/home/wenju/Project/UKBB_NMR_metabolomics/fastGWA/BAG/split_sex/qmplot.sh $output_dir $output_result
  fi
fi
