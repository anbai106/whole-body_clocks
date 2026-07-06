#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=qmplot_mortality
#SBATCH --array=0-21
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --output=/cbica/home/wenju/output/qmplot_mortality_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/qmplot_mortality_%A_%a.err

set -euo pipefail

base_dir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/fastGWA/output"
runner="/cbica/home/wenju/Project/whole-body_clocks/mortality_clock/fastGWA/qmplot_mortality_clock.sh"

# These are the 22 mortality-clock fastGWA output folders reported by ls output.
clock_dirs=(
  "adipose_mri_mortality_clock"
  "brain_mri_mortality_clock"
  "Brain_proteomics_mortality_clock"
  "Digestive_metabolomics_mortality_clock"
  "Endocrine_metabolomics_mortality_clock"
  "Endocrine_proteomics_mortality_clock"
  "Eye_proteomics_mortality_clock"
  "heart_mri_mortality_clock"
  "Heart_proteomics_mortality_clock"
  "Hepatic_metabolomics_mortality_clock"
  "Hepatic_proteomics_mortality_clock"
  "Immune_metabolomics_mortality_clock"
  "Immune_proteomics_mortality_clock"
  "kidney_mri_mortality_clock"
  "liver_mri_mortality_clock"
  "pancreas_mri_mortality_clock"
  "Pulmonary_proteomics_mortality_clock"
  "Renal_proteomics_mortality_clock"
  "Reproductive_female_proteomics_mortality_clock"
  "Reproductive_male_proteomics_mortality_clock"
  "Skin_proteomics_mortality_clock"
  "spleen_mri_mortality_clock"
)

clock_dir=${clock_dirs[$SLURM_ARRAY_TASK_ID]}
output_dir="${base_dir}/${clock_dir}"
output_result="${output_dir}/organ_pheno_normalized_residualized.fastGWA"

manhattan_png="${output_dir}/manhattan_qmplot.png"
qq_png="${output_dir}/QQ_plot.png"

printf 'SLURM_ARRAY_TASK_ID: %s\n' "${SLURM_ARRAY_TASK_ID}"
printf 'Clock directory: %s\n' "${clock_dir}"
printf 'fastGWA file: %s\n' "${output_result}"

if [[ ! -f "${output_result}" ]]; then
  echo "ERROR: fastGWA file does not exist: ${output_result}"
  exit 1
fi

# Skip only when both plots already exist and are non-empty.
if [[ -s "${manhattan_png}" && -s "${qq_png}" ]]; then
  echo "Plots already exist for ${clock_dir}; skipping."
  exit 0
fi

echo "Running qmplot for ${clock_dir}"
bash "${runner}" "${output_dir}" "${output_result}" "${clock_dir}"

echo "Completed ${clock_dir}"
