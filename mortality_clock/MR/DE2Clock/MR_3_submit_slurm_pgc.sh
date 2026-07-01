#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=DE2Clock
#SBATCH --array=0-21
#SBATCH --mem-per-cpu=12G
#SBATCH --time=0-04:59:00
#SBATCH --output=/cbica/home/wenju/output/DE2Clock_MR_3_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/DE2Clock_MR_3_%A_%a.err
############################## END OF DEFAULT EMBEDDED SGE COMMANDS #######################
### define input var

module load R/4.3

numbers=(adipose_mri_mortality_clock Brain_proteomics_mortality_clock Endocrine_metabolomics_mortality_clock Eye_proteomics_mortality_clock Heart_proteomics_mortality_clock Hepatic_proteomics_mortality_clock Immune_proteomics_mortality_clock liver_mri_mortality_clock Pulmonary_proteomics_mortality_clock Reproductive_female_proteomics_mortality_clock Skin_proteomics_mortality_clock brain_mri_mortality_clock Digestive_metabolomics_mortality_clock Endocrine_proteomics_mortality_clock heart_mri_mortality_clock Hepatic_metabolomics_mortality_clock Immune_metabolomics_mortality_clock kidney_mri_mortality_clock pancreas_mri_mortality_clock Renal_proteomics_mortality_clock Reproductive_male_proteomics_mortality_clock spleen_mri_mortality_clock)
organ=${numbers[$SLURM_ARRAY_TASK_ID]}
organ_first="${organ%%_*}"

output_dir=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/MR/DE2Clock/PGC
output_dir_mr=${output_dir}/${organ}/MR
output_dir_har=${output_dir}/${organ}/harmonization
mkdir -p $output_dir_mr
mkdir -p $output_dir_har

### pgc
for pgc in AD ADHD BIP SCZ
do
  harmonized_file=${output_dir_har}/harmonized_data_${pgc}_2_${organ_first}.tsv
  if [ -f "${harmonized_file}" ]; then
    output_file=${output_dir_mr}/MR_${pgc}_2_${organ_first}_OR.tsv
    if [ ! -f "${output_file}" ]; then
      echo "Run 2SampleMR from ${pgc} to ${organ_first}..."
      Rscript /cbica/home/wenju/Project/whole-body_clocks/mortality_clock/MR/DE2Clock/MR_3.R ${pgc} ${organ_first} ${output_dir_mr} ${harmonized_file}
    fi
  fi
done