#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=DE2BAG
#SBATCH --array=0-3
#SBATCH --mem-per-cpu=12G
#SBATCH --output=/cbica/home/wenju/output/DE2Pheno_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/DE2Pheno_%A_%a.err

############################## END OF DEFAULT EMBEDDED SGE COMMANDS #######################
### define input var

module load R/4.3

output_dir=/cbica/projects/MULTI/processed/PGC/GWAS_summary_stats/Result/2SampleMR/PGC2MRIBAG
mkdir -p $output_dir
de_array=( $( awk '{print $2}' /cbica/projects/MULTI/processed/PGC/PGC_MUTATE_MR.tsv | grep -v "Phenotype" | xargs) )
de=${de_array[$SLURM_ARRAY_TASK_ID]}
for organ in brain adipose heart kidney liver pancreas spleen
do
    mr_tsv=${output_dir}/MR_${de}_2_${organ}.tsv
    if [ -f "${mr_tsv}" ]; then
    p_value_thres=0.00023364485981308412 ### 0.05/214 DEs

    # Debugging: Check p-values
    echo "Checking p-values in ${mr_tsv}..."
    awk -F'\t' 'NR > 1 && $9 != "" {print $9}' "${mr_tsv}" | head
    # Extract all p-values in the 9th column and check if any are below the threshold
    significant_pval=$(awk -F'\t' 'NR > 1 && $9 != "" {print $9}' "${mr_tsv}" | \
                       awk '{if ($1+0 < '"$p_value_thres"') print $1}' | head -n 1)
    if [ -n "$significant_pval" ]; then
      echo "Run 2SampleMR QC from ${de} to ${organ}..."
      Rscript /cbica/home/wenju/Project/AbdoImaging/MR/DE2BAG/MR_5_qc.R ${de} ${organ} ${output_dir}
    else
      echo "Not significant for ${de}..."
    fi
  else
    echo "File ${mr_tsv} not found."
  fi
done