#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=brain_mri_asthma_clock
#SBATCH --mem-per-cpu=24G
#SBATCH --time=0-13:59:00
#SBATCH --output=/cbica/home/wenju/output/brain_mri_asthma_clock_%A.out
#SBATCH --error=/cbica/home/wenju/output/brain_mri_asthma_clock_%A.err

set -euo pipefail

module load python/anaconda/3
source activate survival_clock

outdir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/brain_mri_asthma_clock"
mkdir -p "${outdir}"
mkdir -p /cbica/home/wenju/output

# Keep an organ-specific copy of the log inside the output folder.
exec > "${outdir}/brain_mri_asthma_clock_${SLURM_JOB_ID}.out" \
     2> "${outdir}/brain_mri_asthma_clock_${SLURM_JOB_ID}.err"

echo "============================================================"
echo "Starting brain MRI asthma L'EPOCH clock"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}"
echo "Output directory: ${outdir}"
echo "Started at: $(date)"
echo "============================================================"

python /cbica/home/wenju/Project/whole-body_clocks/asthma_clock/MRI/brain_mri_asthma_clock.py \
  --asthma-xlsx /cbica/home/wenju/Dataset/UKBB_UMelbourne/Death_related_var_from_Ye.xlsx \
  --id-match-csv /cbica/home/wenju/Dataset/UKBB_UMelbourne/UKB_UMelbourne_vs_Penn_match_key.csv \
  --brain-tsv /cbica/home/wenju/Reproducibile_paper/BrainAge/data/imaging/T1_MUSE_GM.tsv \
  --covariate-csv /cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv \
  --admin-censor-date 2022-11-30 \
  --outdir "${outdir}" \
  --imaging-session-id 1

echo "============================================================"
echo "Finished brain MRI asthma L'EPOCH clock"
echo "Finished at: $(date)"
echo "============================================================"

conda deactivate
