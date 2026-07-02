#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=Clock2DE
#SBATCH --array=0-520
#SBATCH --mem-per-cpu=12G
#SBATCH --time=1-00:59:00
#SBATCH --output=/cbica/home/wenju/output/Clock2DE%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/Clock2DE%A_%a.err

############################## END OF DEFAULT EMBEDDED SGE COMMANDS #######################
output_dir=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/MR/Clock2DE/FinnGen
de_array=( $( awk '{print $1}' /cbica/projects/MULTI/processed/FinnGen/GWAS_summary_stats/summary_stats_R9_manifest_5000_cases.tsv | grep -v "phenocode" | xargs) )
finngen=${de_array[$SLURM_ARRAY_TASK_ID]}

module load R/4.3

### 7 mri
for organ in brain adipose heart kidney liver pancreas spleen
do
  output_dir_mr=${output_dir}/${organ}_mri_mortality_clock/MR
  output_dir_qc=${output_dir}/${organ}_mri_mortality_clock/QC
  mkdir -p ${output_dir_qc}
  mr_tsv=${output_dir_mr}/MR_${organ}_2_${finngen}_OR.tsv
  if [ -f "${mr_tsv}" ]; then
    p_value_thres=0.00009523809523809524 # 0.05/525

    # Debugging: Check p-values
    echo "Checking p-values in ${mr_tsv}..."
    awk -F'\t' 'NR > 1 && $9 != "" {print $9}' "${mr_tsv}" | head

    # Extract all p-values in the 9th column and check if any are below the threshold
    significant_pval=$(awk -F'\t' 'NR > 1 && $9 != "" {print $9}' "${mr_tsv}" | \
                       awk '{if ($1+0 < '"$p_value_thres"') print $1}' | head -n 1)

    if [ -n "$significant_pval" ]; then
      echo "Run 2SampleMR QC from ${organ} to ${finngen}..."
      Rscript /cbica/home/wenju/Project/whole-body_clocks/mortality_clock/MR/Clock2DE/MR_5_qc.R ${organ} ${finngen} ${output_dir_qc}
    else
      echo "Not significant for ${organ}..."
    fi
  else
    echo "File ${mr_tsv} not found."
  fi
done

### 11 proteomics
for organ in Reproductive_female Pulmonary Heart Brain Eye Hepatic Renal Reproductive_male Endocrine Immune Skin
do
  output_dir_mr=${output_dir}/${organ}_proteomics_mortality_clock/MR
  output_dir_qc=${output_dir}/${organ}_proteomics_mortality_clock/QC
  mkdir -p ${output_dir_qc}
  mr_tsv=${output_dir_mr}/MR_${organ}_2_${finngen}_OR.tsv
  if [ -f "${mr_tsv}" ]; then
    p_value_thres=0.00009523809523809524 # 0.05/525

    # Debugging: Check p-values
    echo "Checking p-values in ${mr_tsv}..."
    awk -F'\t' 'NR > 1 && $9 != "" {print $9}' "${mr_tsv}" | head

    # Extract all p-values in the 9th column and check if any are below the threshold
    significant_pval=$(awk -F'\t' 'NR > 1 && $9 != "" {print $9}' "${mr_tsv}" | \
                       awk '{if ($1+0 < '"$p_value_thres"') print $1}' | head -n 1)

    if [ -n "$significant_pval" ]; then
      echo "Run 2SampleMR QC from ${organ} to ${finngen}..."
      Rscript /cbica/home/wenju/Project/whole-body_clocks/mortality_clock/MR/Clock2DE/MR_5_qc.R ${organ} ${finngen} ${output_dir_qc}
    else
      echo "Not significant for ${organ}..."
    fi
  else
    echo "File ${mr_tsv} not found."
  fi
done

### 4 proteomics
for organ in Endocrine Digestive Hepatic Immune
do
  output_dir_mr=${output_dir}/${organ}_metabolomics_mortality_clock/MR
  output_dir_qc=${output_dir}/${organ}_metabolomics_mortality_clock/QC
  mkdir -p ${output_dir_qc}
  mr_tsv=${output_dir_mr}/MR_${organ}_2_${finngen}_OR.tsv
  if [ -f "${mr_tsv}" ]; then
    p_value_thres=0.00009523809523809524 # 0.05/525

    # Debugging: Check p-values
    echo "Checking p-values in ${mr_tsv}..."
    awk -F'\t' 'NR > 1 && $9 != "" {print $9}' "${mr_tsv}" | head

    # Extract all p-values in the 9th column and check if any are below the threshold
    significant_pval=$(awk -F'\t' 'NR > 1 && $9 != "" {print $9}' "${mr_tsv}" | \
                       awk '{if ($1+0 < '"$p_value_thres"') print $1}' | head -n 1)

    if [ -n "$significant_pval" ]; then
      echo "Run 2SampleMR QC from ${organ} to ${finngen}..."
      Rscript /cbica/home/wenju/Project/whole-body_clocks/mortality_clock/MR/Clock2DE/MR_5_qc.R ${organ} ${finngen} ${output_dir_qc}
    else
      echo "Not significant for ${organ}..."
    fi
  else
    echo "File ${mr_tsv} not found."
  fi
done