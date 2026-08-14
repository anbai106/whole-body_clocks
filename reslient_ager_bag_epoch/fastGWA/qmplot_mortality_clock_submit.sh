#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=qmplot_mortality
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --output=/cbica/home/wenju/output/qmplot_mortality_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/qmplot_mortality_%A_%a.err

set -euo pipefail

output_dir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/fastGWA/output/Brain_proteomics_mortality_clock/EPOCH_BAG_residual"
clock_name=Brain_proteomics_residual
output_result="${output_dir}/organ_pheno_normalized_residualized.fastGWA"

manhattan_png="${output_dir}/manhattan_qmplot.png"
qq_png="${output_dir}/QQ_plot.png"

printf 'Clock directory: %s\n' "${output_dir}"
printf 'fastGWA file: %s\n' "${output_result}"

if [[ ! -f "${output_result}" ]]; then
  echo "ERROR: fastGWA file does not exist: ${output_result}"
  exit 1
fi

# Skip only when both plots already exist and are non-empty.
if [[ -s "${manhattan_png}" && -s "${qq_png}" ]]; then
  echo "Plots already exist; skipping."
  exit 0
fi

# Path on CUBIC where you should place the revised Python script.
plot_script="/cbica/home/wenju/Project/whole-body_clocks/reslient_ager_bag_epoch/fastGWA/5_qmplt_manhatton_mortality_clock.py"

source activate DNE

echo "Start qmplot"
echo "output_dir: ${output_dir}"
echo "output_result: ${output_result}"
echo "clock_name: ${clock_name}"

python -u "${plot_script}" \
  --output_dir "${output_dir}" \
  --output_result "${output_result}" \
  --clock_name "${clock_name}"

echo "Finish qmplot"
conda deactivate
