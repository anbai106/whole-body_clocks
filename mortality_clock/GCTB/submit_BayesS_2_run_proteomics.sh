#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=ArrayJob
#SBATCH --array=0-10
#SBATCH --mem-per-cpu=12G
#SBATCH --time=0-01:59:00
#SBATCH --output=/cbica/home/wenju/output/SBayesS_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/SBayesS_%A_%a.err

numbers=(Reproductive_female Pulmonary Heart Brain Eye Hepatic Renal Reproductive_male Endocrine Immune Skin)
organ=${numbers[$SLURM_ARRAY_TASK_ID]}

output_dir=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/GCTB/${organ}_proteomics_mortality_clock
ldm='/cbica/home/wenju/Software/gctb_2.05beta_Linux/ref_UKBB_sparse_matrix_Zeng2021/ukbEURu_imp_v3_HM3_n50k.chisq10.ldm.sparse'

gwas_summary="${output_dir}/${organ}_gctb.tsv"
output_name="${output_dir}/${organ}"
output_file="${output_name}.parRes"
if [ ! -f ${output_file} ]; then
  echo "Rrun GCTB for ${organ}..."
  bash /cbica/home/wenju/Project/whole-body_clocks/mortality_clock/GCTB/BayesS_2_run.sh $gwas_summary $ldm $output_name
else
  :
fi



