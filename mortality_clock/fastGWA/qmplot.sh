#!/bin/bash
module load python/anaconda/3

output_dir=$1
output_result=$2

source activate DNE
echo "Start training"
python /cbica/home/wenju/Project/UKBB_NMR_metabolomics/fastGWA/BAG/5_qmplt_manhatton.py --output_dir ${output_dir} --output_result ${output_result}

echo "Finish!"
conda deactivate
