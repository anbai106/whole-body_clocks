#!/bin/bash
module load python/anaconda/3

gwas_summary=$1
ldm=$2
output_name=$3

echo "Start training"
/cbica/home/wenju/Software/gctb_2.05beta_Linux/gctb --sbayes S --ldm ${ldm} --gwas-summary ${gwas_summary} --out ${output_name}
echo "Finish!"
