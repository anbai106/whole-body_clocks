#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=Clock2DE
#SBATCH --array=0-520
#SBATCH --mem-per-cpu=12G
#SBATCH --time=0-09:59:00
#SBATCH --output=/cbica/home/wenju/output/MR_2_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/MR_2_%A_%a.err

set -euo pipefail

module load R/4.3

# =============================
# User-defined paths
# =============================
fastgwa_base="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/fastGWA/output"
fuma_base="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/fuma"
mr_output_base="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/MR/Clock2DE/FinnGen"
manifest="/cbica/projects/MULTI/processed/FinnGen/GWAS_summary_stats/summary_stats_R9_manifest_5000_cases.tsv"
r_script="/cbica/home/wenju/Project/whole-body_clocks/mortality_clock/MR/Clock2DE/MR_2_harmonization.R"

# =============================
# Select FinnGen endpoint by SLURM array index
# =============================
mapfile -t de_array < <(awk 'NR > 1 {print $1}' "${manifest}")

if [[ "${SLURM_ARRAY_TASK_ID}" -ge "${#de_array[@]}" ]]; then
  echo "ERROR: SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} exceeds number of FinnGen endpoints=${#de_array[@]}"
  exit 1
fi

finngen="${de_array[$SLURM_ARRAY_TASK_ID]}"

echo "SLURM task ID: ${SLURM_ARRAY_TASK_ID}"
echo "FinnGen endpoint: ${finngen}"

# =============================
# Build the 22-clock list dynamically from fastGWA folders
# Each folder should contain organ_pheno_normalized_residualized.fastGWA
# =============================
mapfile -t fastgwa_files < <(find "${fastgwa_base}" -mindepth 2 -maxdepth 2 -type f -name "*.fastGWA" | sort)

if [[ "${#fastgwa_files[@]}" -eq 0 ]]; then
  echo "ERROR: No fastGWA files found under ${fastgwa_base}"
  exit 1
fi

echo "Number of exposure fastGWA files found: ${#fastgwa_files[@]}"

# =============================
# Loop across all mortality-clock exposures
# =============================
for exposure_gwas_tsv in "${fastgwa_files[@]}"; do
  clock="$(basename "$(dirname "${exposure_gwas_tsv}")")"

  exposure_fuma_tsv="${fuma_base}/${clock}/IndSigSNPs.txt"
  output_dir_mr="${mr_output_base}/${clock}"
  mkdir -p "${output_dir_mr}"

  # This outcome file is expected to have been created by the FinnGen preprocessing/2SampleMR step.
  output_2sampleMR_tsv="/cbica/projects/MULTI/processed/FinnGen/2SampleMR/finngen_R9_${finngen}_2SampleMR.tsv"

  harmonized_file="${output_dir_mr}/harmonized_data_${clock}_2_${finngen}.tsv"
  done_file="${output_dir_mr}/DONE_${clock}_2_${finngen}.tsv"

  echo "----------------------------------------"
  echo "Clock: ${clock}"
  echo "Exposure GWAS: ${exposure_gwas_tsv}"
  echo "FUMA IV file: ${exposure_fuma_tsv}"
  echo "Outcome file: ${output_2sampleMR_tsv}"
  echo "Output dir: ${output_dir_mr}"

  if [[ -f "${harmonized_file}" && -f "${done_file}" ]]; then
    echo "Already done: ${clock} to ${finngen}"
    continue
  fi

  if [[ ! -f "${exposure_gwas_tsv}" ]]; then
    echo "WARNING: missing exposure GWAS file: ${exposure_gwas_tsv}; skipping."
    continue
  fi

  if [[ ! -f "${exposure_fuma_tsv}" ]]; then
    echo "WARNING: missing FUMA IndSigSNPs file: ${exposure_fuma_tsv}; skipping."
    continue
  fi

  if [[ ! -f "${output_2sampleMR_tsv}" ]]; then
    echo "WARNING: missing outcome 2SampleMR file: ${output_2sampleMR_tsv}; skipping."
    continue
  fi

  # Require at least 8 total rows in IndSigSNPs.txt, following your previous logic.
  # This is conservative because the file includes a header.
  num_rows=$(wc -l < "${exposure_fuma_tsv}")
  if [[ "${num_rows}" -le 7 ]]; then
    echo "Too few rows in FUMA IV file (${num_rows}); skipping ${clock}."
    echo "DONE_TOO_FEW_FUMA_IV_ROWS" > "${done_file}"
    continue
  fi

  echo "Harmonizing ${clock} to ${finngen}..."

  Rscript "${r_script}" \
    "${exposure_gwas_tsv}" \
    "${exposure_fuma_tsv}" \
    "${clock}" \
    "${output_2sampleMR_tsv}" \
    "${finngen}" \
    "${output_dir_mr}"

  echo "Finished ${clock} to ${finngen}."
done

echo "All clocks processed for FinnGen endpoint ${finngen}."
