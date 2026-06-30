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

#### FinnGen
output_dir=/cbica/home/wenju/Dataset/FinnGen/2SampleMR/FINNGEN2MRIBAG
mkdir -p $output_dir
de_array=( $( awk '{print $1}' /cbica/home/wenju/Dataset/FinnGen/GWAS_summary_stats/summary_stats_R9_manifest_5000_cases.tsv | grep -v "phenocode" | xargs) )
de=${de_array[$SLURM_ARRAY_TASK_ID]}
mkdir -p $output_dir
exposure_gwas_tsv=/cbica/home/wenju/Dataset/FinnGen/2SampleMR/finngen_R9_${de}_2SampleMR.tsv
plinkn_clumped_lead_snp_tsv=/cbica/home/wenju/Dataset/FinnGen/2SampleMR/plink_clumped/finngen_R9_${de}.clumped
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


