#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=LDSC_disease_clock
#SBATCH --array=0-46
#SBATCH --mem-per-cpu=12G
#SBATCH --time=0-01:59:00
#SBATCH --output=/cbica/home/wenju/output/LDSC_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/LDSC_%A_%a.err

set -euo pipefail

base_dir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock"

# Collect all fastGWA files in a stable sorted order
source_files=(Brain_proteomics_dementia_clock        Heart_proteomics_dementia_clock      Immune_proteomics_stroke_clock
Brain_proteomics_mi_clock              Heart_proteomics_mi_clock            Metabolic_metabolomics_asthma_clock
Brain_proteomics_stroke_clock          Hepatic_metabolomics_asthma_clock    Metabolic_metabolomics_dementia_clock
Digestive_metabolomics_asthma_clock    Hepatic_metabolomics_dementia_clock  Metabolic_metabolomics_stroke_clock
Digestive_metabolomics_dementia_clock  Hepatic_metabolomics_mi_clock        Pulmonary_proteomics_mi_clock
Digestive_metabolomics_mi_clock        Hepatic_metabolomics_stroke_clock    Pulmonary_proteomics_stroke_clock
Digestive_metabolomics_stroke_clock    Hepatic_proteomics_asthma_clock      Reproductive_female_proteomics_asthma_clock
Endocrine_metabolomics_asthma_clock    Hepatic_proteomics_dementia_clock    Reproductive_female_proteomics_copd_clock
Endocrine_metabolomics_dementia_clock  Hepatic_proteomics_mi_clock          Reproductive_female_proteomics_dementia_clock
Endocrine_metabolomics_mi_clock        Hepatic_proteomics_stroke_clock      Reproductive_female_proteomics_mi_clock
Endocrine_metabolomics_stroke_clock    Immune_metabolomics_asthma_clock     Reproductive_female_proteomics_stroke_clock
Endocrine_proteomics_asthma_clock      Immune_metabolomics_mi_clock         Reproductive_male_proteomics_dementia_clock
Endocrine_proteomics_dementia_clock    Immune_metabolomics_stroke_clock     Reproductive_male_proteomics_mi_clock
Endocrine_proteomics_mi_clock          Immune_proteomics_asthma_clock       Reproductive_male_proteomics_stroke_clock
Endocrine_proteomics_stroke_clock      Immune_proteomics_dementia_clock     spleen_mri_asthma_clock
heart_mri_copd_clock                   Immune_proteomics_mi_clock)

organ_dir="${source_files[$SLURM_ARRAY_TASK_ID]}"

source_file=${base_dir}/${organ_dir}/fastGWA/output/organ_pheno_normalized_residualized.fastGWA

echo "SLURM job ID: ${SLURM_JOB_ID}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Source file: ${source_file}"

bash /cbica/home/wenju/Project/whole-body_clocks/LDSC_disease/ldsc_prepare_data.sh "${source_file}"