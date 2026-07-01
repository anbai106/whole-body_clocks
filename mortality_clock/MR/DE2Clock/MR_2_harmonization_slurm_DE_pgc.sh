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

#### PGC
output_dir=/cbica/projects/MULTI/processed/PGC/GWAS_summary_stats/Result/2SampleMR/PGC2MRIBAG
mkdir -p $output_dir
de_array=( $( awk '{print $2}' /cbica/projects/MULTI/processed/PGC/PGC_MUTATE_MR.tsv | grep -v "Phenotype" | xargs) )
de=${de_array[$SLURM_ARRAY_TASK_ID]}
mkdir -p $output_dir
exposure_gwas_tsv=/cbica/projects/MULTI/processed/PGC/GWAS_summary_stats/${de}/2sampleMR.tsv
plinkn_clumped_lead_snp_tsv=/cbica/projects/MULTI/processed/PGC/GWAS_summary_stats/Result/2SampleMR/plink_clumped/${de}.clumped
num_rows=$(wc -l < ${plinkn_clumped_lead_snp_tsv})
if [ -f ${plinkn_clumped_lead_snp_tsv} ] && [ ${num_rows} -gt 7 ]; then
  for organ in brain adipose heart kidney liver pancreas spleen
  do
    output_2sampleMR_tsv=/cbica/home/wenju/Reproducibile_paper/AbdoImaging/fastGWA/${organ}/organ_pheno_normalized_residualized.fastGWA
    if_done_file=${output_dir}/DONE_${de}_2_${organ}.tsv
    if [[ ! -f ${if_done_file} ]]; then
      Rscript /cbica/home/wenju/Project/AbdoImaging/MR/DE2BAG/MR_2_harmonization.R ${exposure_gwas_tsv} ${plinkn_clumped_lead_snp_tsv} ${de} ${output_2sampleMR_tsv} ${organ} ${output_dir}
    fi
  done
fi


