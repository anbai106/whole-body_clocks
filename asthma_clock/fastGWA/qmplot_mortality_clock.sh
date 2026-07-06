#!/bin/bash
set -euo pipefail

module load python/anaconda/3

output_dir=$1
output_result=$2
clock_name=${3:-$(basename "${output_dir}")}

# Path on CUBIC where you should place the revised Python script.
plot_script="/cbica/home/wenju/Project/whole-body_clocks/mortality_clock/fastGWA/5_qmplt_manhatton_mortality_clock.py"

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
