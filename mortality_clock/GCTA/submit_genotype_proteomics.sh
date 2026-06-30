#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=ArrayJob
#SBATCH --array=0-10
#SBATCH --mem-per-cpu=96G
#SBATCH --time=0-10:59:00
#SBATCH --output=/cbica/home/wenju/output/GCTA_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/GCTA_%A_%a.err

numbers=(Reproductive_female Pulmonary Heart Brain Eye Hepatic Renal Reproductive_male Endocrine Immune Skin)
organ=${numbers[$SLURM_ARRAY_TASK_ID]}

echo "Run GCTA for brain phenotype: ${organ}"
input_dir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/fastGWA/data"
output_dir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/GCTA_h2/${organ}_proteomics_mortality_clock"
mkdir -p ${output_dir}
cd $output_dir

## A step to derive the bed binary file and also the replicate participants.
keep="${input_dir}/${organ}_proteomics_mortality_clock/EPOCH_keep_for_fastgwa.txt"
pheno="${input_dir}/${organ}_proteomics_mortality_clock/EPOCH_pheno_normalized_residualized_with_related_indi.phen"
out_name=${organ}
grm=/cbica/home/wenju/Reproducibile_paper/UKBB_Proteomics/h2/GRM_GCTA/gcta_h2_grm ## need rerun GRM for h2
if [[ -f ${output_dir}/${organ}_proteomics_mortality_clock/${organ}.hsq ]]; then
    echo "GCTA has been run..."
    :
else
    bash /cbica/home/wenju/Project/whole-body_clocks/mortality_clock/GCTA/run_genotype.sh $output_dir 1 $out_name $pheno $keep ${grm}
fi