#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=DE2BAG
#SBATCH --array=0-520
#SBATCH --mem-per-cpu=12G
#SBATCH --time=1-00:59:00
#SBATCH --output=/cbica/home/wenju/output/DE2Pheno_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/DE2Pheno_%A_%a.err

############################## END OF DEFAULT EMBEDDED SGE COMMANDS #######################
### finngenfine input var

module load R/4.3

#### FinnGen
output_dir=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/MR/DE2Clock/FinnGen
mkdir -p $output_dir
finngen_array=( $( awk '{print $1}' /cbica/projects/MULTI/processed/FinnGen/GWAS_summary_stats/summary_stats_R9_manifest_5000_cases.tsv | grep -v "phenocofinngen" | xargs) )
finngen=${finngen_array[$SLURM_ARRAY_TASK_ID]}

### 7 MRI
for organ in brain adipose heart kidney liver pancreas spleen
do
    output_dir_har=${output_dir}/${organ}_mri_mortality_clock/harmonization
    output_dir_mr=${output_dir}/${organ}_mri_mortality_clock/MR
    output_dir_qc=${output_dir}/${organ}_mri_mortality_clock/QC
    mkdir -p ${output_dir_qc}
    mr_tsv=${output_dir_mr}/MR_${finngen}_2_${organ}.tsv
    if [ -f "${mr_tsv}" ]; then
    p_value_thres=0.000234742 ### 0.05/N DEs

    # Debugging: Check p-values
    echo "Checking p-values in ${mr_tsv}..."
    awk -F'\t' 'NR > 1 && $9 != "" {print $9}' "${mr_tsv}" | head
    # Extract all p-values in the 9th column and check if any are below the threshold
    significant_pval=$(awk -F'\t' 'NR > 1 && $9 != "" {print $9}' "${mr_tsv}" | \
                       awk '{if ($1+0 < '"$p_value_thres"') print $1}' | head -n 1)
    if [ -n "$significant_pval" ]; then
      echo "Run 2SampleMR QC from ${finngen} to ${organ}..."
      Rscript /cbica/home/wenju/Project/whole-body_clocks/mortality_clock/MR/DE2Clock/MR_5_qc.R ${finngen} ${organ} ${output_dir_qc} ${output_dir_har} ${output_dir_mr}
    else
      echo "Not significant for ${finngen}..."
    fi
  else
    echo "File ${mr_tsv} not found."
  fi
done

### 11 proteomics
for organ in Reproductive_female Pulmonary Heart Brain Eye Hepatic Renal Reproductive_male Endocrine Immune Skin
do
    output_dir_har=${output_dir}/${organ}_proteomics_mortality_clock/harmonization
    output_dir_mr=${output_dir}/${organ}_proteomics_mortality_clock/MR
    output_dir_qc=${output_dir}/${organ}_proteomics_mortality_clock/QC
    mkdir -p ${output_dir_qc}
    mr_tsv=${output_dir_mr}/MR_${finngen}_2_${organ}.tsv
    if [ -f "${mr_tsv}" ]; then
    p_value_thres=0.000234742 ### 0.05/214 DEs

    # Debugging: Check p-values
    echo "Checking p-values in ${mr_tsv}..."
    awk -F'\t' 'NR > 1 && $9 != "" {print $9}' "${mr_tsv}" | head
    # Extract all p-values in the 9th column and check if any are below the threshold
    significant_pval=$(awk -F'\t' 'NR > 1 && $9 != "" {print $9}' "${mr_tsv}" | \
                       awk '{if ($1+0 < '"$p_value_thres"') print $1}' | head -n 1)
    if [ -n "$significant_pval" ]; then
      echo "Run 2SampleMR QC from ${finngen} to ${organ}..."
      Rscript /cbica/home/wenju/Project/whole-body_clocks/mortality_clock/MR/DE2Clock/MR_5_qc.R ${finngen} ${organ} ${output_dir_qc} ${output_dir_har} ${output_dir_mr}
    else
      echo "Not significant for ${finngen}..."
    fi
  else
    echo "File ${mr_tsv} not found."
  fi
done

### 4 proteomics
for organ in Endocrine Digestive Hepatic Immune
do
    output_dir_har=${output_dir}/${organ}_metabolomics_mortality_clock/harmonization
    output_dir_mr=${output_dir}/${organ}_metabolomics_mortality_clock/MR
    output_dir_qc=${output_dir}/${organ}_metabolomics_mortality_clock/QC
    mkdir -p ${output_dir_qc}
    mr_tsv=${output_dir_mr}/MR_${finngen}_2_${organ}.tsv
    if [ -f "${mr_tsv}" ]; then
      p_value_thres=0.000234742 ### 0.05/214 DEs

      # Debugging: Check p-values
      echo "Checking p-values in ${mr_tsv}..."
      awk -F'\t' 'NR > 1 && $9 != "" {print $9}' "${mr_tsv}" | head
      # Extract all p-values in the 9th column and check if any are below the threshold
      significant_pval=$(awk -F'\t' 'NR > 1 && $9 != "" {print $9}' "${mr_tsv}" | \
                         awk '{if ($1+0 < '"$p_value_thres"') print $1}' | head -n 1)
      if [ -n "$significant_pval" ]; then
        echo "Run 2SampleMR QC from ${finngen} to ${organ}..."
        Rscript /cbica/home/wenju/Project/whole-body_clocks/mortality_clock/MR/DE2Clock/MR_5_qc.R ${finngen} ${organ} ${output_dir_qc} ${output_dir_har} ${output_dir_mr}
      else
        echo "Not significant for ${finngen}..."
      fi
    else
    echo "File ${mr_tsv} not found."
    fi
done