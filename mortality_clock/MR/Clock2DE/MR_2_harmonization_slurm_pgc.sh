#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=Clock2PGC
#SBATCH --array=0-21
#SBATCH --mem-per-cpu=12G
#SBATCH --time=0-04:59:00
#SBATCH --output=/cbica/home/wenju/output/MR_PGC_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/MR_PGC_%A_%a.err

set -euo pipefail

module load R/4.3

# =========================
# Paths
# =========================

fastgwa_base="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/fastGWA/output"
fuma_base="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/fuma"
output_dir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/MR/Clock2DE/PGC"

r_script="/cbica/home/wenju/Project/whole-body_clocks/mortality_clock/MR/Clock2DE/MR_2_harmonization.R"

# PGC outcome files
pgc_list=(AD ADHD BIP SCZ)

# =========================
# Get all 22 mortality clocks
# =========================

mapfile -t fastgwa_files < <(
    find "${fastgwa_base}" -mindepth 2 -maxdepth 2 -name "organ_pheno_normalized_residualized.fastGWA" | sort
)

n_fastgwa=${#fastgwa_files[@]}

echo "Number of fastGWA files found: ${n_fastgwa}"

if [[ "${n_fastgwa}" -ne 22 ]]; then
    echo "ERROR: Expected 22 fastGWA files, but found ${n_fastgwa}."
    printf '%s\n' "${fastgwa_files[@]}"
    exit 1
fi

exposure_gwas_tsv="${fastgwa_files[$SLURM_ARRAY_TASK_ID]}"
clock=$(basename "$(dirname "${exposure_gwas_tsv}")")

exposure_fuma_tsv="${fuma_base}/${clock}/IndSigSNPs.txt"
output_dir_mr="${output_dir}/${clock}/harmonization"

mkdir -p "${output_dir_mr}"

echo "SLURM task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Clock: ${clock}"
echo "Exposure GWAS: ${exposure_gwas_tsv}"
echo "Exposure FUMA SNP file: ${exposure_fuma_tsv}"
echo "Output directory: ${output_dir_mr}"

# =========================
# Check exposure files
# =========================

if [[ ! -f "${exposure_gwas_tsv}" ]]; then
    echo "ERROR: Missing exposure GWAS file: ${exposure_gwas_tsv}"
    exit 1
fi

if [[ ! -f "${exposure_fuma_tsv}" ]]; then
    echo "WARNING: Missing FUMA IndSigSNPs file: ${exposure_fuma_tsv}"
    echo "Skipping ${clock}."
    exit 0
fi

num_rows=$(wc -l < "${exposure_fuma_tsv}")

if [[ "${num_rows}" -le 7 ]]; then
    echo "WARNING: Too few rows in ${exposure_fuma_tsv}: ${num_rows}"
    echo "Skipping ${clock}."
    exit 0
fi

# =========================
# Harmonize exposure with PGC outcomes
# =========================

for pgc in "${pgc_list[@]}"; do

    output_2sampleMR_tsv="/cbica/projects/MULTI/processed/PGC/GWAS_summary_stats/${pgc}/2sampleMR.tsv"

    harmonized_file="${output_dir_mr}/harmonized_data_${clock}_2_${pgc}.tsv"
    done_file="${output_dir_mr}/DONE_${clock}_2_${pgc}.tsv"

    echo "----------------------------------------"
    echo "PGC outcome: ${pgc}"
    echo "Outcome 2SampleMR file: ${output_2sampleMR_tsv}"
    echo "Harmonized output: ${harmonized_file}"
    echo "Done file: ${done_file}"

    if [[ ! -f "${output_2sampleMR_tsv}" ]]; then
        echo "WARNING: Missing PGC outcome file: ${output_2sampleMR_tsv}"
        echo "Skipping ${clock} -> ${pgc}."
        continue
    fi

    if [[ -f "${harmonized_file}" && -f "${done_file}" ]]; then
        echo "Harmonization already done for ${clock} -> ${pgc}."
        continue
    fi

    echo "Harmonizing ${clock} -> ${pgc}..."

    Rscript "${r_script}" \
        "${exposure_gwas_tsv}" \
        "${exposure_fuma_tsv}" \
        "${clock}" \
        "${output_2sampleMR_tsv}" \
        "${pgc}" \
        "${output_dir_mr}"

    echo "Finished harmonizing ${clock} -> ${pgc}."

done

echo "All PGC harmonization tasks finished for ${clock}."