#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=BAG2DE
#SBATCH --array=0-520
#SBATCH --mem-per-cpu=12G
#SBATCH --time=1-01:59:00
#SBATCH --output=/cbica/home/wenju/output/MLNI_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/MLNI_%A_%a.err

module load R/4.3

output_dir=/cbica/home/wenju/Reproducibile_paper/AbdoImaging/MR/BAG2DE
de_array=( $( awk '{print $1}' /cbica/home/wenju/Dataset/FinnGen/GWAS_summary_stats/summary_stats_R9_manifest_5000_cases.tsv | grep -v "phenocode" | xargs) )
finngen=${de_array[$SLURM_ARRAY_TASK_ID]}

for organ in brain adipose heart kidney liver pancreas spleen
do
  output_dir_mr=${output_dir}/${organ}
  harmonized_file=$output_dir_mr/harmonized_data_${organ}_2_${finngen}.tsv
  if [ -f "${harmonized_file}" ]; then
    output_file=${output_dir_mr}/MR_${organ}_2_${finngen}_OR.tsv
    if [ ! -f "${output_file}" ]; then
      echo "Run 2SampleMR from ${organ} to ${finngen}}..."
      Rscript /cbica/home/wenju/Project/AbdoImaging/MR/BAG2DE/MR_3_run.R ${organ} ${finngen} ${output_dir_mr} ${harmonized_file}
    fi
  fi
done