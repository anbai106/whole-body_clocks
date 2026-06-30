#!/bin/bash

module load plink/2.20210701

output_dir=$1
k=$2
out_name=$3
pheno=$4
keep=$5
grm=$6

cd $output_dir

/cbica/home/wenju/Software/gcta_1.93.2beta/gcta64 --HEreg \
    --grm "${grm}" \
    --pheno "${pheno}" \
    --keep "${keep}" \
    --mpheno "${k}" \
    --thread-num 4 \
    --out "${out_name}"