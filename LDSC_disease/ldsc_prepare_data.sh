module load python/anaconda/3

source_file=$1


source activate sopnmf_atlas
echo "Start training"
python /cbica/home/wenju/Project/whole-body_clocks/LDSC_disease/ldsc_prepare_data.py --source_file ${source_file}

echo "Finish!"
conda deactivate
