#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=sex_stratified
#SBATCH --time=0-09:59:00
#SBATCH --array=0-1
#SBATCH --mem-per-cpu=24G
#SBATCH --output=/cbica/home/wenju/output/sex_stratified_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/sex_stratified_%A_%a.err

numbers=(female male)
sex=${numbers[$SLURM_ARRAY_TASK_ID]}
organ=Metabolic

echo "Run fastGWA for organ: ${organ} on ${sex} data..."
input_dir="/cbica/home/wenju/Reproducibile_paper/UKBB_metabolomics/fastGWA_MetBAG/data/R1_sex_pooled/"
output_dir="/cbica/home/wenju/Reproducibile_paper/UKBB_metabolomics/fastGWA_MetBAG/output/R1_sex_pooled/${organ}_non_derived_${sex}"
mkdir -p ${output_dir}
cd $output_dir

## A step to derive the bed binary file and also the replicate participants.
bfile="/cbica/home/wenju/Reproducibile_paper/UKBB_metabolomics/data/metabolite_binary"
pheno="${input_dir}/${organ}_non_derived_${sex}/BAG_pheno_normalized_residualized_with_related_indi.phen"
output_file="${output_dir}/BAG_pheno_normalized_residualized"
if [[ -f ${output_file}.fastGWA ]]; then
    echo "GWAS has been run..."
    :
else
  echo "GWAS for ${organ}..."
  sparse_grm=/cbica/home/wenju/Reproducibile_paper/UKBB_metabolomics/data/GRM_GCTA/fastGWA_grm_sparse_0.05
  bash /cbica/home/wenju/Project/UKBB_NMR_metabolomics/fastGWA/BAG/sex_pooled_metabolic_MetBAG_R1/4_runGWAS/run.sh $output_dir $bfile $sparse_grm $pheno $output_file
fi