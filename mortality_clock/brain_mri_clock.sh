#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=SAnalysis
#SBATCH --mem-per-cpu=24G
#SBATCH --output=/cbica/home/wenju/output/SAnalysis_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/SAnalysis_%A_%a.err

############################## END OF DEFAULT EMBEDDED SGE COMMANDS #######################
source activate survival_clock
python /cbica/home/wenju/Project/whole-body_clocks/mortality_clock/brain_mri_clock.py \
  --death-xlsx /cbica/home/wenju/Dataset/UKBB_UMelbourne/Death_related_var_from_Ye.xlsx \
  --id-match-csv /cbica/home/wenju/Dataset/UKBB_UMelbourne/UKB_UMelbourne_vs_Penn_match_key.csv \
  --brain-tsv /cbica/home/wenju/Reproducibile_paper/BrainAge/data/imaging/T1_MUSE_GM.tsv \
  --covariate-csv /cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv \
  --admin-censor-date 2022-11-30 \
  --outdir /cbica/home/wenju/Reproducibile_paper/WholeBodyClock/brain_mri_mortality_clock
conda deactivate