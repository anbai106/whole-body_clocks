#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=DE2BAG
#SBATCH --array=0-3
#SBATCH --mem-per-cpu=12G
#SBATCH --output=/cbica/home/wenju/output/array_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/array_%A_%a.err

############################## END OF DEFAULT EMBEDDED SGE COMMANDS #######################
### define input var

module load R/4.3

output_dir=/cbica/projects/MULTI/processed/PGC/GWAS_summary_stats/Result/2SampleMR/PGC2MRIBAG
mkdir -p $output_dir
de_array=( $( awk '{print $2}' /cbica/projects/MULTI/processed/PGC/PGC_MUTATE_MR.tsv | grep -v "Phenotype" | xargs) )
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