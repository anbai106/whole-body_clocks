#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=Clock2DE_PGC_QC
#SBATCH --array=0-21
#SBATCH --mem-per-cpu=12G
#SBATCH --time=0-05:59:00
#SBATCH --output=/cbica/home/wenju/output/Clock2DE_PGC_QC_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/Clock2DE_PGC_QC_%A_%a.err

set -euo pipefail

module load R/4.3

base_dir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/MR/Clock2DE/PGC"
qc_script="/cbica/home/wenju/Project/whole-body_clocks/mortality_clock/MR/Clock2DE/MR_5_qc.R"

mkdir -p /cbica/home/wenju/output

# Keep your original threshold.
# You can change this to 0.05/(22*4) or another PGC-specific correction later if desired.
p_value_thres="0.00009523809523809524"

clock_specs=(
    # 7 MRI mortality clocks
    "adipose_mri_mortality_clock|adipose"
    "brain_mri_mortality_clock|brain"
    "heart_mri_mortality_clock|heart"
    "kidney_mri_mortality_clock|kidney"
    "liver_mri_mortality_clock|liver"
    "pancreas_mri_mortality_clock|pancreas"
    "spleen_mri_mortality_clock|spleen"

    # 11 proteomics mortality clocks
    "Brain_proteomics_mortality_clock|Brain"
    "Endocrine_proteomics_mortality_clock|Endocrine"
    "Eye_proteomics_mortality_clock|Eye"
    "Heart_proteomics_mortality_clock|Heart"
    "Hepatic_proteomics_mortality_clock|Hepatic"
    "Immune_proteomics_mortality_clock|Immune"
    "Pulmonary_proteomics_mortality_clock|Pulmonary"
    "Renal_proteomics_mortality_clock|Renal"
    "Reproductive_female_proteomics_mortality_clock|Reproductive_female"
    "Reproductive_male_proteomics_mortality_clock|Reproductive_male"
    "Skin_proteomics_mortality_clock|Skin"

    # 4 metabolomics mortality clocks
    "Digestive_metabolomics_mortality_clock|Digestive"
    "Endocrine_metabolomics_mortality_clock|Endocrine"
    "Hepatic_metabolomics_mortality_clock|Hepatic"
    "Immune_metabolomics_mortality_clock|Immune"
)

n_clocks=${#clock_specs[@]}

if [[ "${SLURM_ARRAY_TASK_ID}" -ge "${n_clocks}" ]]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}, but only ${n_clocks} clocks found."
    exit 1
fi

IFS="|" read -r organ_folder organ_prefix <<< "${clock_specs[$SLURM_ARRAY_TASK_ID]}"

echo "SLURM job ID: ${SLURM_JOB_ID}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Clock folder: ${organ_folder}"
echo "MR file prefix / Rscript organ argument: ${organ_prefix}"
echo "P-value threshold: ${p_value_thres}"

output_dir_mr="${base_dir}/${organ_folder}/MR"
output_dir_qc="${base_dir}/${organ_folder}/QC"

mkdir -p "${output_dir_qc}"

if [[ ! -d "${output_dir_mr}" ]]; then
    echo "MR folder not found: ${output_dir_mr}. Exiting."
    exit 0
fi

if [[ -z "$(find "${output_dir_mr}" -maxdepth 1 -type f -name "*.tsv" -print -quit)" ]]; then
    echo "MR folder is empty: ${output_dir_mr}. Exiting."
    exit 0
fi

for pgc in AD ADHD BIP SCZ; do

    mr_tsv="${output_dir_mr}/MR_${organ_prefix}_2_${pgc}.tsv"

    if [[ ! -s "${mr_tsv}" ]]; then
        echo "MR result not found or empty: ${mr_tsv}. Skipping."
        continue
    fi

    echo "Checking all MR estimator P values in: ${mr_tsv}"

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
        echo "At least one MR estimator is significant for ${organ_folder} -> ${pgc}."
        echo "Significant method(s) and P value(s):"
        echo "${significant_rows}"

        echo "Running 2SampleMR QC for ${organ_prefix} -> ${pgc}..."
        Rscript "${qc_script}" "${organ_prefix}" "${pgc}" "${output_dir_qc}"
    else
        echo "No MR estimator P value below threshold for ${organ_folder} -> ${pgc}. Skipping QC."
    fi

done

echo "Finished PGC QC screening for clock: ${organ_folder}"