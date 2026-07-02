#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=organ_mri_clock
#SBATCH --mem-per-cpu=24G
#SBATCH --array=0-5
#SBATCH --output=/cbica/home/wenju/output/organ_mri_mi_clock_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/organ_mri_mi_clock_%A_%a.err

source activate survival_clock

numbers=(adipose heart kidney liver pancreas spleen)
organ=${numbers[$SLURM_ARRAY_TASK_ID]}

outdir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/${organ}_mri_mi_clock"
mkdir -p "${outdir}"

# Redirect main job log dynamically into each organ-specific output folder
exec > "${outdir}/${organ}_mri_mi_clock_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out" \
     2> "${outdir}/${organ}_mri_mi_clock_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"

echo "============================================================"
echo "Starting organ MRI mi clock"
echo "Organ: ${organ}"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "Output directory: ${outdir}"
echo "Started at: $(date)"
echo "============================================================"

python /cbica/home/wenju/Project/whole-body_clocks/mi_clock/MRI/organ_mri_clock.py \
  --mi-xlsx /cbica/home/wenju/Dataset/UKBB_UMelbourne/Death_related_var_from_Ye.xlsx \
  --id-match-csv /cbica/home/wenju/Dataset/UKBB_UMelbourne/UKB_UMelbourne_vs_Penn_match_key.csv \
  --organ-tsv "/cbica/home/wenju/Reproducibile_paper/AbdoImaging/Data_split_by_Huizi/MLNI_dataset_splits/${organ}_train_val_test.tsv,/cbica/home/wenju/Reproducibile_paper/AbdoImaging/Data_split_by_Huizi/MLNI_dataset_splits/${organ}_independent_test.tsv,/cbica/home/wenju/Reproducibile_paper/AbdoImaging/Data_split_by_Huizi/raw_data_MLNI/${organ}_patient_mlni.tsv" \
  --covariate-csv /cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv \
  --admin-censor-date 2022-11-30 \
  --outdir "${outdir}" \
  --organ "${organ}"

echo "============================================================"
echo "Finished organ MRI mi clock"
echo "Organ: ${organ}"
echo "Finished at: $(date)"
echo "============================================================"

conda deactivate