#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=Clock2DE
#SBATCH --array=0-21
#SBATCH --mem-per-cpu=12G
#SBATCH --time=1-00:59:00
#SBATCH --output=/cbica/home/wenju/output/Clock2DE%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/Clock2DE%A_%a.err

############################## END OF DEFAULT EMBEDDED SGE COMMANDS #######################
output_dir=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/MR/Clock2DE/PGC
numbers=(adipose_mri_mortality_clock Brain_proteomics_mortality_clock Endocrine_metabolomics_mortality_clock Eye_proteomics_mortality_clock Heart_proteomics_mortality_clock Hepatic_proteomics_mortality_clock Immune_proteomics_mortality_clock liver_mri_mortality_clock Pulmonary_proteomics_mortality_clock Reproductive_female_proteomics_mortality_clock Skin_proteomics_mortality_clock brain_mri_mortality_clock Digestive_metabolomics_mortality_clock Endocrine_proteomics_mortality_clock heart_mri_mortality_clock Hepatic_metabolomics_mortality_clock Immune_metabolomics_mortality_clock kidney_mri_mortality_clock pancreas_mri_mortality_clock Renal_proteomics_mortality_clock Reproductive_male_proteomics_mortality_clock spleen_mri_mortality_clock)
organ=${numbers[$SLURM_ARRAY_TASK_ID]}
organ_first="${organ%%_*}"

module load R/4.3

### PGC
for pgc in AD ADHD BIP SCZ
do
  output_dir_mr=${output_dir}/${organ}/MR
  mr_tsv=${output_dir_mr}/MR_${organ_first}_2_${pgc}.tsv
  output_dir_qc=${output_dir}/${organ}/QC
  mkdir -p ${output_dir_qc}
  if [ -f "${mr_tsv}" ]; then
    p_value_thres=0.00009523809523809524 # 0.05/525

    # Debugging: Check p-values
    echo "Checking p-values in ${mr_tsv}..."
    awk -F'\t' 'NR > 1 && $9 != "" {print $9}' "${mr_tsv}" | head

    # Extract all p-values in the 9th column and check if any are below the threshold
    significant_pval=$(awk -F'\t' 'NR > 1 && $9 != "" {print $9}' "${mr_tsv}" | \
                       awk '{if ($1+0 < '"$p_value_thres"') print $1}' | head -n 1)

    if [ -n "$significant_pval" ]; then
      echo "Run 2SampleMR QC from ${organ_first} to ${pgc}..."
      Rscript /cbica/home/wenju/Project/whole-body_clocks/mortality_clock/MR/Clock2DE/MR_5_qc.R ${organ_first} ${pgc} ${output_dir_qc}
    else
      echo "Not significant for ${organ_first}..."
    fi
  else
    echo "File ${mr_tsv} not found."
  fi
done