#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=Clock2DE
#SBATCH --array=0-520
#SBATCH --mem-per-cpu=12G
#SBATCH --time=0-04:59:00
#SBATCH --output=/cbica/home/wenju/output/Clock2DE_3_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/Clock2DE_3_%A_%a.err

module load R/4.3

output_dir=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/MR/Clock2DE/FinnGen
de_array=( $( awk '{print $1}' /cbica/projects/MULTI/processed/FinnGen/GWAS_summary_stats/summary_stats_R9_manifest_5000_cases.tsv | grep -v "phenocode" | xargs) )
finngen=${de_array[$SLURM_ARRAY_TASK_ID]}

### 7 mri
for organ in brain adipose heart kidney liver pancreas spleen
do
  output_dir_mr=${output_dir}/${organ}_mri_mortality_clock/MR
  output_dir_har=${output_dir}/${organ}_mri_mortality_clock/harmonization
  harmonized_file=$output_dir_har/harmonized_data_${organ}_mri_mortality_clock_2_${finngen}.tsv
  if [ -f "${harmonized_file}" ]; then
    output_file=${output_dir_mr}/MR_${organ}_2_${finngen}_OR.tsv
    if [ ! -f "${output_file}" ]; then
      echo "Run 2SampleMR from ${organ} to ${finngen}}..."
      Rscript /cbica/home/wenju/Project/whole-body_clocks/mortality_clock/MR/Clock2DE/MR_3_run.R ${organ} ${finngen} ${output_dir_mr} ${harmonized_file}
    fi
  fi
done

### 11 proteomics
for organ in Reproductive_female Pulmonary Heart Brain Eye Hepatic Renal Reproductive_male Endocrine Immune Skin
do
  output_dir_mr=${output_dir}/${organ}_proteomics_mortality_clock/MR
  output_dir_har=${output_dir}/${organ}_proteomics_mortality_clock/harmonization
  harmonized_file=$output_dir_har/harmonized_data_${organ}_proteomics_mortality_clock_2_${finngen}.tsv
  if [ -f "${harmonized_file}" ]; then
    output_file=${output_dir_mr}/MR_${organ}_2_${finngen}_OR.tsv
    if [ ! -f "${output_file}" ]; then
      echo "Run 2SampleMR from ${organ} to ${finngen}}..."
      Rscript /cbica/home/wenju/Project/whole-body_clocks/mortality_clock/MR/Clock2DE/MR_3_run.R ${organ} ${finngen} ${output_dir_mr} ${harmonized_file}
    fi
  fi
done

### 4 proteomics
for organ in Endocrine Digestive Hepatic Immune
do
  output_dir_mr=${output_dir}/${organ}_metabolomics_mortality_clock/MR
  output_dir_har=${output_dir}/${organ}_metabolomics_mortality_clock/harmonization
  harmonized_file=$output_dir_har/harmonized_data_${organ}_metabolomics_mortality_clock_2_${finngen}.tsv
  if [ -f "${harmonized_file}" ]; then
    output_file=${output_dir_mr}/MR_${organ}_2_${finngen}_OR.tsv
    if [ ! -f "${output_file}" ]; then
      echo "Run 2SampleMR from ${organ} to ${finngen}}..."
      Rscript /cbica/home/wenju/Project/whole-body_clocks/mortality_clock/MR/Clock2DE/MR_3_run.R ${organ} ${finngen} ${output_dir_mr} ${harmonized_file}
    fi
  fi
done