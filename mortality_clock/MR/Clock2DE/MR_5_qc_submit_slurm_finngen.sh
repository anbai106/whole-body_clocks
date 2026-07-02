#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=Clock2DE_QC
#SBATCH --array=0-520
#SBATCH --mem-per-cpu=12G
#SBATCH --time=0-05:59:00
#SBATCH --output=/cbica/home/wenju/output/Clock2DE_QC_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/Clock2DE_QC_%A_%a.err

set -euo pipefail

module load R/4.3

base_dir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/MR/Clock2DE/FinnGen"
manifest="/cbica/projects/MULTI/processed/FinnGen/GWAS_summary_stats/summary_stats_R9_manifest_5000_cases.tsv"
qc_script="/cbica/home/wenju/Project/whole-body_clocks/mortality_clock/MR/Clock2DE/MR_5_qc.R"

mkdir -p /cbica/home/wenju/output

# Bonferroni threshold for 525 FinnGen endpoints
p_value_thres="0.00009523809523809524"

# Read FinnGen phenocodes
mapfile -t de_array < <(awk 'NR > 1 {print $1}' "${manifest}")

n_de=${#de_array[@]}
echo "Found ${n_de} FinnGen endpoints."

if [[ "${SLURM_ARRAY_TASK_ID}" -ge "${n_de}" ]]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}, but only ${n_de} FinnGen endpoints found."
    exit 1
fi

finngen="${de_array[$SLURM_ARRAY_TASK_ID]}"

echo "SLURM job ID: ${SLURM_JOB_ID}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "FinnGen endpoint: ${finngen}"
echo "P-value threshold: ${p_value_thres}"

# ------------------------------------------------------------
# Define clock folders and MR file prefixes.
#
# Format:
#   organ_for_Rscript|clock_folder|mr_file_prefix
#
# Example:
#   brain|brain_mri_mortality_clock|brain
# means:
#   MR file: ${base_dir}/brain_mri_mortality_clock/MR/MR_brain_2_${finngen}_OR.tsv
#   Rscript argument: brain
# ------------------------------------------------------------

clock_specs=(
    # 7 MRI mortality clocks
    "brain|brain_mri_mortality_clock|brain"
    "adipose|adipose_mri_mortality_clock|adipose"
    "heart|heart_mri_mortality_clock|heart"
    "kidney|kidney_mri_mortality_clock|kidney"
    "liver|liver_mri_mortality_clock|liver"
    "pancreas|pancreas_mri_mortality_clock|pancreas"
    "spleen|spleen_mri_mortality_clock|spleen"

    # 11 proteomics mortality clocks
    "Reproductive_female|Reproductive_female_proteomics_mortality_clock|Reproductive_female"
    "Pulmonary|Pulmonary_proteomics_mortality_clock|Pulmonary"
    "Heart|Heart_proteomics_mortality_clock|Heart"
    "Brain|Brain_proteomics_mortality_clock|Brain"
    "Eye|Eye_proteomics_mortality_clock|Eye"
    "Hepatic|Hepatic_proteomics_mortality_clock|Hepatic"
    "Renal|Renal_proteomics_mortality_clock|Renal"
    "Reproductive_male|Reproductive_male_proteomics_mortality_clock|Reproductive_male"
    "Endocrine|Endocrine_proteomics_mortality_clock|Endocrine"
    "Immune|Immune_proteomics_mortality_clock|Immune"
    "Skin|Skin_proteomics_mortality_clock|Skin"

    # 4 metabolomics mortality clocks
    "Endocrine|Endocrine_metabolomics_mortality_clock|Endocrine"
    "Digestive|Digestive_metabolomics_mortality_clock|Digestive"
    "Hepatic|Hepatic_metabolomics_mortality_clock|Hepatic"
    "Immune|Immune_metabolomics_mortality_clock|Immune"
)

for spec in "${clock_specs[@]}"; do

    IFS="|" read -r organ clock_folder mr_prefix <<< "${spec}"

    output_dir_mr="${base_dir}/${clock_folder}/MR"
    output_dir_qc="${base_dir}/${clock_folder}/QC"

    mkdir -p "${output_dir_qc}"

    if [[ ! -d "${output_dir_mr}" ]]; then
        echo "MR folder not found: ${output_dir_mr}. Skipping."
        continue
    fi

    if [[ -z "$(find "${output_dir_mr}" -maxdepth 1 -type f -name "*.tsv" -print -quit)" ]]; then
        echo "MR folder is empty: ${output_dir_mr}. Skipping."
        continue
    fi

    mr_tsv="${output_dir_mr}/MR_${mr_prefix}_2_${finngen}_OR.tsv"

    if [[ ! -s "${mr_tsv}" ]]; then
        echo "MR result not found or empty: ${mr_tsv}. Skipping."
        continue
    fi

    echo "Checking all MR estimator P values in: ${mr_tsv}"

    # ------------------------------------------------------------
    # Check whether ANY p-value from the 5 MR estimators is below
    # the Bonferroni threshold.
    #
    # Expected columns:
    # method = column 5
    # pval   = column 9
    #
    # This checks all rows after the header, including:
    #   MR Egger
    #   Weighted median
    #   Inverse variance weighted
    #   Simple mode
    #   Weighted mode
    # ------------------------------------------------------------

    significant_rows=$(
        awk -F'\t' -v thres="${p_value_thres}" '
            NR > 1 && $9 != "" && $9 != "NA" && $9 != "NaN" {
                p = $9 + 0
                if (p < thres) {
                    print $5 "\t" $9
                }
            }
        ' "${mr_tsv}"
    )

    if [[ -n "${significant_rows}" ]]; then
        echo "At least one MR estimator is significant for ${clock_folder} -> ${finngen}."
        echo "Significant method(s) and P value(s):"
        echo "${significant_rows}"

        echo "Running 2SampleMR QC for ${clock_folder} -> ${finngen}..."
        Rscript "${qc_script}" "${organ}" "${finngen}" "${output_dir_qc}"
    else
        echo "No MR estimator P value below threshold for ${clock_folder} -> ${finngen}. Skipping QC."
    fi

done

echo "Finished QC screening for FinnGen endpoint: ${finngen}"