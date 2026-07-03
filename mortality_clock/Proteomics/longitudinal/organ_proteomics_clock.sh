#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=apply_pulmonary_proteomics_mortality_clock_longitudinal
#SBATCH --mem-per-cpu=24G
#SBATCH --time=0-04:00:00
#SBATCH --output=/cbica/home/wenju/output/apply_pulmonary_proteomics_mortality_clock_longitudinal_%A.out
#SBATCH --error=/cbica/home/wenju/output/apply_pulmonary_proteomics_mortality_clock_longitudinal_%A.err

set -euo pipefail

module load python/anaconda/3
source activate survival_clock

outdir="Reproducibile_paper/WholeBodyClock/mortality_clock/longitudinal/proteomics/Pulmonary"
mkdir -p "${outdir}"
mkdir -p /cbica/home/wenju/output

exec > "${outdir}/apply_pulmonary_proteomics_mortality_clock_longitudinal_${SLURM_JOB_ID}.out" \
     2> "${outdir}/apply_pulmonary_proteomics_mortality_clock_longitudinal_${SLURM_JOB_ID}.err"

echo "============================================================"
echo "Applying pre-trained Pulmonary proteomics mortality clock to longitudinal instances"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}"
echo "Output directory: ${outdir}"
echo "Started at: $(date)"
echo "============================================================"

python /cbica/home/wenju/Project/whole-body_clocks/mortality_clock/Proteomics/longitudinal/organ_proteomics_clock.py \
  --model-joblib /cbica/home/wenju/Reproducibile_paper/WholeBodyClock/Pulmonary_proteomics_mortality_clock/pulmonary_proteomics_mortality_clock_model.joblib \
  --input-tsv 2_0:/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/longitudinal/proteomics/data/proteimics_Pulmonary_2_0.tsv \
  --input-tsv 3_0:/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/longitudinal/proteomics/data/proteimics_Pulmonary_3_0.tsv \
  --covariate-csv /cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv \
  --assessment-xlsx /cbica/home/wenju/Dataset/UKBB_UMelbourne/Death_related_var_from_Ye.xlsx \
  --id-match-csv /cbica/home/wenju/Dataset/UKBB_UMelbourne/UKB_UMelbourne_vs_Penn_match_key.csv \
  --outdir "${outdir}" \
  --organ pulmonary

echo "============================================================"
echo "Finished applying Pulmonary proteomics mortality clock"
echo "Finished at: $(date)"
echo "============================================================"

conda deactivate
