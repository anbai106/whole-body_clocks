#!/bin/bash

module load python/anaconda/3

source_file=$1
target_file=$2

source activate DNE
echo "Start training"
python  /cbica/home/wenju/Project/AbdoImaging/GCTB/convert.py --source_file ${source_file} --target_file ${target_file}

echo "Finish!"
conda deactivate
