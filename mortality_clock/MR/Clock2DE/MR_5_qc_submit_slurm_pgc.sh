#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=Clock2DE
#SBATCH --array=0-6
#SBATCH --mem-per-cpu=12G
#SBATCH --time=1-00:59:00
#SBATCH --output=/cbica/home/wenju/output/Clock2DE%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/Clock2DE%A_%a.err

############################## END OF DEFAULT EMBEDDED SGE COMMANDS #######################
output_dir=/cbica/home/wenju/Reproducibile_paper/AbdoImaging/MR/Clock2DE
numbers=(brain adipose heart kidney liver pancreas spleen)
organ=${numbers[$SLURM_ARRAY_TASK_ID]}

module load R/4.3

### PGC
for pgc in AD ADHD BIP SCZ
do
  mr_tsv=${output_dir}/${organ}/MR_${organ}_2_${pgc}.tsv
  if [ -f "${mr_tsv}" ]; then
    p_value_thres=0.00009523809523809524 # 0.05/525

    # Debugging: Check p-values
    echo "Checking p-values in ${mr_tsv}..."
    awk -F'\t' 'NR > 1 && $9 != "" {print $9}' "${mr_tsv}" | head

    # Extract all p-values in the 9th column and check if any are below the threshold
    significant_pval=$(awk -F'\t' 'NR > 1 && $9 != "" {print $9}' "${mr_tsv}" | \
                       awk '{if ($1+0 < '"$p_value_thres"') print $1}' | head -n 1)

    if [ -n "$significant_pval" ]; then
      echo "Run 2SampleMR QC from ${organ} to ${pgc}..."
      Rscript /cbica/home/wenju/Project/AbdoImaging/MR/Clock2DE/MR_5_qc.R ${organ} ${pgc} ${output_dir}
    else
      echo "Not significant for ${organ}..."
    fi
  else
    echo "File ${mr_tsv} not found."
  fi
done