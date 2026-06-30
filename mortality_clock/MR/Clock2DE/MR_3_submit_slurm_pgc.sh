#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=BAG2DE
#SBATCH --array=0-6
#SBATCH --mem-per-cpu=12G
#SBATCH --time=1-01:59:00
#SBATCH --output=/cbica/home/wenju/output/MLNI_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/MLNI_%A_%a.err

numbers=(brain adipose heart kidney liver pancreas spleen)
organ=${numbers[$SLURM_ARRAY_TASK_ID]}

module load R/4.3

output_dir=/cbica/home/wenju/Reproducibile_paper/AbdoImaging/MR/BAG2DE
output_dir_mr=${output_dir}/${organ}

### pgc
for pgc in AD ADHD BIP SCZ
do
  harmonized_file=$output_dir_mr/harmonized_data_${organ}_2_${pgc}.tsv
  if [ -f "${harmonized_file}" ]; then
    output_file=${output_dir_mr}/MR_${organ}_2_${pgc}_OR.tsv
    if [ ! -f "${output_file}" ]; then
      echo "Run 2SampleMR from ${organ} to ${pgc}}..."
      Rscript /cbica/home/wenju/Project/AbdoImaging/MR/BAG2DE/MR_3_run.R ${organ} ${pgc} ${output_dir_mr} ${harmonized_file}
    fi
  fi
done