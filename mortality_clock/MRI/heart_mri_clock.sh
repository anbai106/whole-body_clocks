#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=heart_mri_clock
#SBATCH --mem-per-cpu=24G
#SBATCH --output=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/heart_mri_mortality_clock_%A_%a.out
#SBATCH --error=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/heart_mri_mortality_clock_%A_%a.err

############################## END OF DEFAULT EMBEDDED SGE COMMANDS #######################
source activate survival_clock

mkdir -p /cbica/home/wenju/Reproducibile_paper/WholeBodyClock/heart_mri_mortality_clock

python /cbica/home/wenju/Project/whole-body_clocks/mortality_clock/MRI/heart_mri_clock.py \
  --death-xlsx /cbica/home/wenju/Dataset/UKBB_UMelbourne/Death_related_var_from_Ye.xlsx \
  --id-match-csv /cbica/home/wenju/Dataset/UKBB_UMelbourne/UKB_UMelbourne_vs_Penn_match_key.csv \
  --heart-tsv "/cbica/home/wenju/Reproducibile_paper/AbdoImaging/Data_split_by_Huizi/MLNI_dataset_splits/heart_train_val_test.tsv,/cbica/home/wenju/Reproducibile_paper/AbdoImaging/Data_split_by_Huizi/MLNI_dataset_splits/heart_independent_test.tsv,/cbica/home/wenju/Reproducibile_paper/AbdoImaging/Data_split_by_Huizi/raw_data_MLNI/heart_patient_mlni.tsv" \
  --covariate-csv /cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv \
  --admin-censor-date 2022-11-30 \
  --outdir /cbica/home/wenju/Reproducibile_paper/WholeBodyClock/heart_mri_mortality_clock

conda deactivate
