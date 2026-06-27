######START OF EMBEDDED SGE COMMANDS ##########################
#$ -S /bin/bash  ### Specifies the interpreting shell for the job
#$ -N pyOPNMF #### Name of the job
#$ -M anbai106@hotmail.com #### email to nofity with following options/scenarios
#$ -m a #### send mail in case the job is aborted
#$ -m b #### send mail when job begins
#$ -m e #### send mail when job ends
#$ -m s #### send mail when job is suspende
#$ -l h_vmem=40G  ### required cpu memory
#$ -o /cbica/home/wenju/output ### output path
#$ -e /cbica/home/wenju/output ### error path

############################## END OF DEFAULT EMBEDDED SGE COMMANDS #######################

module load python/anaconda/3

stat_gz=$1
output_file=$2
ref_ld_chr="/cbica/home/wenju/Project/ldsc/pre-computed_Eu_GWAS/eur_w_ld_chr/"  ### the / at the end is very important, corresponding to the issues here on github: https://github.com/bulik/ldsc/issues/250
w_ld_chr="/cbica/home/wenju/Project/ldsc/pre-computed_Eu_GWAS/eur_w_ld_chr/"

source activate ldsc
echo "Start training"
python /cbica/home/wenju/Project/ldsc/ldsc.py --h2 ${stat_gz} --ref-ld-chr ${ref_ld_chr} --w-ld-chr ${w_ld_chr} --out ${output_file}

echo "Finish!"
conda deactivate