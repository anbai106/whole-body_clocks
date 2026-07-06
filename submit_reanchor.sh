#!/bin/bash
#SBATCH --job-name=reanchor_lepoch_years
#SBATCH --output=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/lepoch_year_scale_reanchoring/reanchor_lepoch_years_%j.out
#SBATCH --error=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/lepoch_year_scale_reanchoring/reanchor_lepoch_years_%j.err
#SBATCH --partition=all
#SBATCH --time=08:00:00
#SBATCH --mem=32G

set -euo pipefail

echo "============================================================"
echo "Re-anchor disease L'EPOCH year-scale metrics and age coefficient"
echo "Job ID: ${SLURM_JOB_ID:-NA}"
echo "Host: $(hostname)"
echo "Date: $(date)"
echo "============================================================"

# -----------------------------
# 1. Environment
# -----------------------------

source ~/.bashrc || true
conda activate survival_clock

# -----------------------------
# 2. Paths
# -----------------------------

base_dir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock"

script_py="/cbica/home/wenju/Project/whole-body_clocks/reanchor_age_coefficient_4_problematic_disease_clock.py"

out_dir="${base_dir}/lepoch_year_scale_reanchoring"

mkdir -p "${out_dir}"

# -----------------------------
# 3. Run mode
# -----------------------------
# APPLY mode overwrites files after creating backups.
# For dry run, set:
#   apply_flag=""
# -----------------------------

apply_flag="--apply"

# -----------------------------
# 4. Diseases to process
# -----------------------------

diseases=("asthma" "dementia" "copd" "mi" "stroke")

# -----------------------------
# 5. Run
# -----------------------------
# This will update:
#   *_clock_predictions.tsv
#   *_clock_test_predictions.tsv
#   *_clock_coefficients.tsv
#   *_clock_nonzero_coefficients.tsv if age row exists
#
# Backups are saved in the same original clock folders with:
#   .pre_reanchor_year_scale.<timestamp>.bak
# -----------------------------

python "${script_py}" \
  --base_dir "${base_dir}" \
  --diseases "${diseases[@]}" \
  --out_dir "${out_dir}" \
  ${apply_flag} \
  --also_test_predictions \
  --update_coefficients \
  --update_nonzero_coefficients \
  --qc_file_kind problematic \
  --backup_suffix pre_reanchor_year_scale \
  --min_events 20 \
  --age_beta_min_abs 0.005 \
  --winsorize_score_quantiles 0.001 0.999

echo "============================================================"
echo "Finished re-anchoring."
echo "Date: $(date)"
echo "Output directory:"
echo "  ${out_dir}"
echo "Check latest summary:"
echo "  ${out_dir}/lepoch_year_scale_reanchoring_applied_latest.tsv"
echo "============================================================"