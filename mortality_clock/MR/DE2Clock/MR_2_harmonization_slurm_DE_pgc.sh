#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=DE2Clock
#SBATCH --array=0-3
#SBATCH --mem-per-cpu=12G
#SBATCH --output=/cbica/home/wenju/output/DE2Clock_MR_2_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/DE2Clock_MR_2_%A_%a.err

############################## END OF DEFAULT EMBEDDED SGE COMMANDS #######################
### define input var

module load R/4.3

#### PGC
output_dir=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/MR/DE2Clock/PGC
mkdir -p $output_dir
de_array=( $( awk '{print $2}' /cbica/projects/MULTI/processed/PGC/PGC_MUTATE_MR.tsv | grep -v "Phenotype" | xargs) )
pgc=${de_array[$SLURM_ARRAY_TASK_ID]}
mkdir -p $output_dir

exposure_gwas_tsv=/cbica/projects/MULTI/processed/PGC/GWAS_summary_stats/${pgc}/2sampleMR.tsv
plink_clumped_lead_snp_tsv=/cbica/projects/MULTI/processed/PGC/GWAS_summary_stats/Result/2SampleMR/plink_clumped/${pgc}.clumped
num_rows=$(wc -l < ${plink_clumped_lead_snp_tsv})
if [ -f ${plink_clumped_lead_snp_tsv} ] && [ ${num_rows} -gt 7 ]; then
  ### 6 mri
  for organ in brain adipose heart kidney liver pancreas spleen
  do
    output_dir_mr=${output_dir}/${organ}_mri_mortality_clock/harmonization
    output_2sampleMR_tsv=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/fastGWA/output/${organ}_mri_mortality_clock/organ_pheno_normalized_residualized.fastGWA
    if_done_file=${output_dir_mr}/DONE_${pgc}_2_${organ}.tsv
    if [[ ! -f ${if_done_file} ]]; then
      Rscript /cbica/home/wenju/Project/whole-body_clocks/mortality_clock/MR/DE2Clock/MR_2_harmonization.R ${exposure_gwas_tsv} ${plink_clumped_lead_snp_tsv} ${pgc} ${output_2sampleMR_tsv} ${organ} ${output_dir_mr}
    fi
  done

  ### 11 proteomics
  for organ in Reproductive_female Pulmonary Heart Brain Eye Hepatic Renal Reproductive_male Endocrine Immune Skin
  do
    output_dir_mr=${output_dir}/${organ}_proteomics_mortality_clock/harmonization
    output_2sampleMR_tsv=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/fastGWA/output/${organ}_proteomics_mortality_clock/organ_pheno_normalized_residualized.fastGWA
    if_done_file=${output_dir_mr}/DONE_${pgc}_2_${organ}.tsv
    if [[ ! -f ${if_done_file} ]]; then
      Rscript /cbica/home/wenju/Project/whole-body_clocks/mortality_clock/MR/DE2Clock/MR_2_harmonization.R ${exposure_gwas_tsv} ${plink_clumped_lead_snp_tsv} ${pgc} ${output_2sampleMR_tsv} ${organ} ${output_dir_mr}
    fi
  done

  ### 4 metabolomics
  for organ in Endocrine Digestive Hepatic Immune
  do
    output_dir_mr=${output_dir}/${organ}_metabolomics_mortality_clock/harmonization
    output_2sampleMR_tsv=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/fastGWA/output/${organ}_metabolomics_mortality_clock/organ_pheno_normalized_residualized.fastGWA
    if_done_file=${output_dir_mr}/DONE_${pgc}_2_${organ}.tsv
    if [[ ! -f ${if_done_file} ]]; then
      Rscript /cbica/home/wenju/Project/whole-body_clocks/mortality_clock/MR/DE2Clock/MR_2_harmonization.R ${exposure_gwas_tsv} ${plink_clumped_lead_snp_tsv} ${pgc} ${output_2sampleMR_tsv} ${organ} ${output_dir_mr}
    fi
  done
fi


