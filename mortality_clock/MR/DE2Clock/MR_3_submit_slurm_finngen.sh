#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=DE2BAG
#SBATCH --array=0-520
#SBATCH --mem-per-cpu=12G
#SBATCH --time=1-00:59:00
#SBATCH --output=/cbica/home/wenju/output/array_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/array_%A_%a.err

############################## END OF DEFAULT EMBEDDED SGE COMMANDS #######################
### define input var

module load R/4.3

output_dir=/cbica/home/wenju/Dataset/FinnGen/2SampleMR/FINNGEN2MRIBAG
mkdir -p $output_dir
de_array=( $( awk '{print $1}' /cbica/home/wenju/Dataset/FinnGen/GWAS_summary_stats/summary_stats_R9_manifest_5000_cases.tsv | grep -v "phenocode" | xargs) )
de=${de_array[$SLURM_ARRAY_TASK_ID]}
for organ in brain adipose heart kidney liver pancreas spleen
do
  harmonized_file=${output_dir}/harmonized_data_${de}_2_${organ}.tsv
  if [ -f "${harmonized_file}" ]; then
    output_file=${output_dir}/MR_${de}_2_${organ}_OR.tsv
    if [ ! -f "${output_file}" ]; then
      echo "Run 2SampleMR from ${de} to ${organ}..."
      Rscript /cbica/home/wenju/Project/AbdoImaging/MR/DE2BAG/MR_3.R ${de} ${organ} ${output_dir} ${harmonized_file}
    fi
  fi
done