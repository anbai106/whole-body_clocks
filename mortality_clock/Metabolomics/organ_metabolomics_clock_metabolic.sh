#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=organ_metabolomics_clock
#SBATCH --mem-per-cpu=24G
#SBATCH --array=0-0
#SBATCH --output=/cbica/home/wenju/output/organ_metabolomics_mortality_clock_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/organ_metabolomics_mortality_clock_%A_%a.err

set -euo pipefail

source activate survival_clock

# Update this list if your metabolomics systems differ.
# The array index must match: #SBATCH --array=0-(number_of_organs - 1)
numbers=(Metabolic)
organ=${numbers[$SLURM_ARRAY_TASK_ID]}

# Main paths
project_root="/cbica/home/wenju/Project/whole-body_clocks/mortality_clock/Metabolomics"
script_path="${project_root}/organ_metabolomics_clock_metabolic.py"
wholebody_root="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock"
outdir="${wholebody_root}/${organ}_metabolomics_mortality_clock"
mkdir -p "${outdir}"
mkdir -p /cbica/home/wenju/output

# Redirect the real job log dynamically into each organ-specific output folder.
# The #SBATCH output/error paths above are only temporary bootstrap logs because
# Slurm parses #SBATCH lines before the bash variable ${organ} is defined.
exec > "${outdir}/${organ}_metabolomics_mortality_clock_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out" \
     2> "${outdir}/${organ}_metabolomics_mortality_clock_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"

echo "============================================================"
echo "Starting organ metabolomics mortality clock"
echo "Organ: ${organ}"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "Output directory: ${outdir}"
echo "Started at: $(date)"
echo "============================================================"

# Numerical-stability settings for Coxnet metabolomics models.
# These arguments correspond to the revised Python script that fixes:
#   ArithmeticError: Numerical error, because weights are too large.
alpha_min_ratio="0.01"
covariate_penalty_factor="0.05"
l1_ratios="0.05,0.1,0.25,0.5"

# The broad Metabolic clock uses all metabolites and is more collinear, so use
# a slightly stronger minimum alpha to avoid unstable Coxnet paths.
if [[ "${organ}" == "Metabolic" ]]; then
  alpha_min_ratio="0.05"
  l1_ratios="0.05,0.1,0.25"
fi

echo "Coxnet settings:"
echo "  alpha_min_ratio=${alpha_min_ratio}"
echo "  covariate_penalty_factor=${covariate_penalty_factor}"
echo "  l1_ratios=${l1_ratios}"

# Edit this root path if your metabolomics feature files are stored elsewhere.
data_root="/cbica/home/wenju/Reproducibile_paper/UKBB_metabolomics/MLNI/data"
organ_tsv="${data_root}/${organ}/training/training_4589.tsv,${data_root}/${organ}/PT/patient_pop.tsv,${data_root}/${organ}/test/ind_test_500.tsv"

python "${script_path}" \
  --death-xlsx /cbica/home/wenju/Dataset/UKBB_UMelbourne/Death_related_var_from_Ye.xlsx \
  --id-match-csv /cbica/home/wenju/Dataset/UKBB_UMelbourne/UKB_UMelbourne_vs_Penn_match_key.csv \
  --organ-tsv "${organ_tsv}" \
  --covariate-csv /cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv \
  --admin-censor-date 2022-11-30 \
  --outdir "${outdir}" \
  --organ "${organ}" \
  --omics-session-id none \
  --alpha-min-ratio "${alpha_min_ratio}" \
  --covariate-penalty-factor "${covariate_penalty_factor}" \
  --l1-ratios "${l1_ratios}"

echo "============================================================"
echo "Finished organ metabolomics mortality clock"
echo "Organ: ${organ}"
echo "Finished at: $(date)"
echo "============================================================"

conda deactivate
