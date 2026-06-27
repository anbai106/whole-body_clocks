module load python/anaconda/3

file=$1
output_file=$2

source activate ldsc
echo "Start training"

python /cbica/home/wenju/Project/ldsc/munge_sumstats.py --sumstats ${file}  --out ${output_file} --merge-alleles /cbica/home/wenju/Project/ldsc/pre-computed_Eu_GWAS/w_hm3.snplist --chunksize 50000

echo "Finish!"
conda deactivate