#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=Clock2DE
#SBATCH --array=0-21
#SBATCH --mem-per-cpu=12G
#SBATCH --time=0-04:59:00
#SBATCH --output=/cbica/home/wenju/output/MLNI_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/MLNI_%A_%a.err

numbers=(adipose_mri_mortality_clock Brain_proteomics_mortality_clock Endocrine_metabolomics_mortality_clock Eye_proteomics_mortality_clock Heart_proteomics_mortality_clock Hepatic_proteomics_mortality_clock Immune_proteomics_mortality_clock liver_mri_mortality_clock Pulmonary_proteomics_mortality_clock Reproductive_female_proteomics_mortality_clock Skin_proteomics_mortality_clock brain_mri_mortality_clock Digestive_metabolomics_mortality_clock Endocrine_proteomics_mortality_clock heart_mri_mortality_clock Hepatic_metabolomics_mortality_clock Immune_metabolomics_mortality_clock kidney_mri_mortality_clock pancreas_mri_mortality_clock Renal_proteomics_mortality_clock Reproductive_male_proteomics_mortality_clock spleen_mri_mortality_clock)
organ=${numbers[$SLURM_ARRAY_TASK_ID]}
organ_first="${organ%%_*}"
module load R/4.3

output_dir=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/MR/Clock2DE/PGC
output_dir_mr=${output_dir}/${organ}

### pgc
for pgc in AD ADHD BIP SCZ
do
  harmonized_file=$output_dir_mr/harmonized_data_${organ}_2_${pgc}.tsv
  if [ -f "${harmonized_file}" ]; then
    output_file=${output_dir_mr}/MR_${organ_first}_2_${pgc}_OR.tsv
    if [ ! -f "${output_file}" ]; then
      echo "Run 2SampleMR from ${organ} to ${pgc}}..."
      Rscript /cbica/home/wenju/Project/whole-body_clocks/mortality_clock/MR/Clock2DE/MR_3_run.R ${organ_first} ${pgc} ${output_dir_mr} ${harmonized_file}
    fi
  fi
done