#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=organ_proteomics_stroke_clock
#SBATCH --mem-per-cpu=24G
#SBATCH --array=0-10
#SBATCH --output=/cbica/home/wenju/output/organ_proteomics_stroke_clock_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/organ_proteomics_stroke_clock_%A_%a.err

set -euo pipefail

module load python/anaconda/3
source activate survival_clock

numbers=(Reproductive_female Pulmonary Heart Brain Eye Hepatic Renal Reproductive_male Endocrine Immune Skin)

if [[ "${SLURM_ARRAY_TASK_ID}" -ge "${#numbers[@]}" ]]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} but only ${#numbers[@]} organs are defined."
    exit 1
fi

organ=${numbers[$SLURM_ARRAY_TASK_ID]}

outdir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/${organ}_proteomics_stroke_clock"
mkdir -p "${outdir}"
mkdir -p /cbica/home/wenju/output

# Redirect main job log dynamically into each organ-specific output folder
exec > "${outdir}/${organ}_proteomics_stroke_clock_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out" \
     2> "${outdir}/${organ}_proteomics_stroke_clock_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"

echo "============================================================"
echo "Starting organ proteomics stroke clock"
echo "Organ: ${organ}"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "Output directory: ${outdir}"
echo "Started at: $(date)"
echo "============================================================"

python /cbica/home/wenju/Project/whole-body_clocks/stroke_clock/Proteomics/organ_proteomics_clock.py \
  --stroke-xlsx /cbica/home/wenju/Dataset/UKBB_UMelbourne/Death_related_var_from_Ye.xlsx \
  --id-match-csv /cbica/home/wenju/Dataset/UKBB_UMelbourne/UKB_UMelbourne_vs_Penn_match_key.csv \
  --organ-tsv "/cbica/home/wenju/Reproducibile_paper/UKBB_Proteomics/MLNI/data/${organ}/training/training_4589.tsv,/cbica/home/wenju/Reproducibile_paper/UKBB_Proteomics/MLNI/data/${organ}/PT/patient_pop.tsv,/cbica/home/wenju/Reproducibile_paper/UKBB_Proteomics/MLNI/data/${organ}/test/ind_test_500.tsv" \
  --covariate-csv /cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv \
  --admin-censor-date 2022-11-30 \
  --outdir "${outdir}" \
  --organ "${organ}" \
  --omics-session-id ses-M0

echo "============================================================"
echo "Finished organ proteomics stroke clock"
echo "Organ: ${organ}"
echo "Finished at: $(date)"
echo "============================================================"

conda deactivate
