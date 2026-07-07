#!/bin/bash

output_dir=$1
bfile=$2
sparse_grm=$3
pheno=$4
output_file=$5

cd ${output_dir}

# To run a fastGWA analysis based on the sparse GRM generated above
/cbica/home/wenju/Software/gcta_1.93.2beta/gcta64 --bfile ${bfile} --grm-sparse ${sparse_grm} --fastGWA-mlm --pheno ${pheno} --thread-num 8 --out ${output_file}