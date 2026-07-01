#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=DE2Clock
#SBATCH --array=0-520
#SBATCH --mem-per-cpu=12G
#SBATCH --time=0-09:59:00
#SBATCH --output=/cbica/home/wenju/output/DE2Clock_MR_2_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/DE2Clock_MR_2_%A_%a.err

############################## END OF DEFAULT EMBEDDED SGE COMMANDS #######################
### define input var

module load R/4.3

#### FinnGen
output_dir=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/MR/DE2Clock/FinnGen
mkdir -p $output_dir
de_array=( $( awk '{print $1}' /cbica/projects/MULTI/processed/FinnGen/GWAS_summary_stats/summary_stats_R9_manifest_5000_cases.tsv | grep -v "phenocode" | xargs) )
finngen=${de_array[$SLURM_ARRAY_TASK_ID]}
mkdir -p $output_dir

exposure_gwas_tsv=/cbica/projects/MULTI/processed/FinnGen/2SampleMR/finngen_R9_${finngen}_2SampleMR.tsv
plink_clumped_lead_snp_tsv=/cbica/projects/MULTI/processed/FinnGen/2SampleMR/plink_clumped/finngen_R9_${finngen}.clumped
num_rows=$(wc -l < ${plink_clumped_lead_snp_tsv})
if [ -f ${plink_clumped_lead_snp_tsv} ] && [ ${num_rows} -gt 7 ]; then

  ### 7 mri
  for organ in brain adipose heart kidney liver pancreas spleen
  do
    output_dir_mr=${output_dir}/${organ}_mri_mortality_clock/harmonization
    output_2sampleMR_tsv=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/fastGWA/output/${organ}_mri_mortality_clock/organ_pheno_normalized_residualized.fastGWA
    if_done_file=${output_dir_mr}/DONE_${finngen}_2_${organ}.tsv
    if [[ ! -f ${if_done_file} ]]; then
      Rscript /cbica/home/wenju/Project/whole-body_clocks/mortality_clock/MR/DE2Clock/MR_2_harmonization.R ${exposure_gwas_tsv} ${plink_clumped_lead_snp_tsv} ${finngen} ${output_2sampleMR_tsv} ${organ} ${output_dir_mr}
    fi
  done

  ### 11 proteomics
  for organ in Reproductive_female Pulmonary Heart Brain Eye Hepatic Renal Reproductive_male Endocrine Immune Skin
  do
    output_dir_mr=${output_dir}/${organ}_proteomics_mortality_clock/harmonization
    output_2sampleMR_tsv=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/fastGWA/output/${organ}_proteomics_mortality_clock/organ_pheno_normalized_residualized.fastGWA
    if_done_file=${output_dir_mr}/DONE_${finngen}_2_${organ}.tsv
    if [[ ! -f ${if_done_file} ]]; then
      Rscript /cbica/home/wenju/Project/whole-body_clocks/mortality_clock/MR/DE2Clock/MR_2_harmonization.R ${exposure_gwas_tsv} ${plink_clumped_lead_snp_tsv} ${finngen} ${output_2sampleMR_tsv} ${organ} ${output_dir_mr}
    fi
  done

  ### 4 metabolomics
  for organ in Endocrine Digestive Hepatic Immune
  do
    output_dir_mr=${output_dir}/${organ}_metabolomics_mortality_clock/harmonization
    output_2sampleMR_tsv=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/fastGWA/output/${organ}_metabolomics_mortality_clock/organ_pheno_normalized_residualized.fastGWA
    if_done_file=${output_dir_mr}/DONE_${finngen}_2_${organ}.tsv
    if [[ ! -f ${if_done_file} ]]; then
      Rscript /cbica/home/wenju/Project/whole-body_clocks/mortality_clock/MR/DE2Clock/MR_2_harmonization.R ${exposure_gwas_tsv} ${plink_clumped_lead_snp_tsv} ${finngen} ${output_2sampleMR_tsv} ${organ} ${output_dir_mr}
    fi
  done
fi


