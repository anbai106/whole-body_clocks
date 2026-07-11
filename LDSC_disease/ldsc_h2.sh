#!/bin/bash

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