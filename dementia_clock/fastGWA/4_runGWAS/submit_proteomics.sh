#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=ArrayJob
#SBATCH --array=0-3
#SBATCH --time=0-10:59:00
#SBATCH --mem-per-cpu=24G
#SBATCH --output=/cbica/home/wenju/output/fastGWA_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/fastGWA_%A_%a.err

numbers=(Reproductive_female Hepatic Endocrine Immune)
organ=${numbers[$SLURM_ARRAY_TASK_ID]}
disease=asthma
modality=proteomics
echo "Run fastGWA for organ: ${organ}"
input_dir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/${organ}_${modality}_${disease}_clock/fastGWA/data"
output_dir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/${organ}_${modality}_${disease}_clock/fastGWA/output"
mkdir -p ${output_dir}
cd $output_dir

## A step to derive the bed binary file and also the replicate participants.
bfile="/cbica/home/wenju/Reproducibile_paper/UKBB_Proteomics/fastGWA/protein_keep_for_BAG_with_related_indi"
pheno="${input_dir}/EPOCH_pheno_normalized_residualized_with_related_indi.phen"
output_file="${output_dir}/organ_pheno_normalized_residualized"
if [[ -f ${output_file}.fastGWA ]]; then
    echo "GWAS has been run..."
    :
else
  echo "GWAS for ${organ}..."
  sparse_grm=/cbica/home/wenju/Reproducibile_paper/UKBB_Proteomics/fastGWA/GRM_GCTA/fastGWA_grm_sparse_0.05
  bash /cbica/home/wenju/Project/whole-body_clocks/asthma_clock/fastGWA/4_runGWAS/run.sh $output_dir $bfile $sparse_grm $pheno $output_file
fi