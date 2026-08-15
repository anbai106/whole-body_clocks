#!/bin/bash

sleep_dir=$1
idp_tsv=$2
output_dir=$3
idp=$4
cov_tsv=$5

module load R/4.2.2
echo "Start training"
Rscript /cbica/home/wenju/Project/SleepAging/NSS/BWAS/fit_logistic.R ${sleep_dir} ${idp_tsv} ${output_dir} ${idp} ${cov_tsv}
echo "Finish!"
