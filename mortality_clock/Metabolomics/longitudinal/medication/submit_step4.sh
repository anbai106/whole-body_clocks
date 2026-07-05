#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=resp_med_immune_delta
#SBATCH --mem-per-cpu=24G
#SBATCH --time=0-06:00:00
#SBATCH --output=/cbica/home/wenju/output/resp_med_immune_delta_%j.out
#SBATCH --error=/cbica/home/wenju/output/resp_med_immune_delta_%j.err

set -euo pipefail

module load python/anaconda/3

mkdir -p /cbica/home/wenju/output

# ------------------------------------------------------------
# Input files from your previous medication-cluster pipeline
# ------------------------------------------------------------

input_long_tsv=${INPUT_LONG_TSV:-/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/SA/delta_metabolomics_algorithmic_disease_onset/medication_cluster_delta_clock_inputs/metabolomics_delta_clock_medication_cluster_requested5_long.tsv}

med_long_tsv=${MED_LONG_TSV:-/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/SA/delta_metabolomics_algorithmic_disease_onset/medication_cluster_delta_clock_inputs/medication_instance0_long_classified.tsv}

participant_cluster_tsv=${PARTICIPANT_CLUSTER_TSV:-/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/SA/delta_metabolomics_algorithmic_disease_onset/medication_cluster_delta_clock_inputs/medication_participant_clusters.tsv}

umel_death_xlsx=${UMEL_DEATH_XLSX:-/cbica/home/wenju/Dataset/UKBB_UMelbourne/Death_related_var_from_Ye.xlsx}

umel_match_csv=${UMEL_MATCH_CSV:-/cbica/home/wenju/Dataset/UKBB_UMelbourne/UKB_UMelbourne_vs_Penn_match_key.csv}

# ------------------------------------------------------------
# Output and script paths
# ------------------------------------------------------------

out_dir=${OUT_DIR:-/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/SA/delta_metabolomics_algorithmic_disease_onset/respiratory_medication_subtype_immune_delta_results}

python_script=/cbica/home/wenju/Project/whole-body_clocks/mortality_clock/Metabolomics/longitudinal/medication/4_immune_respirotary_sensitivity.py

outcome=${OUTCOME:-delta_clock_age_years}

echo "============================================================"
echo "Respiratory medication subtype analysis for immune delta clock"
echo "Input long TSV: ${input_long_tsv}"
echo "Medication long TSV: ${med_long_tsv}"
echo "Participant cluster TSV: ${participant_cluster_tsv}"
echo "UMel disease/date Excel: ${umel_death_xlsx}"
echo "UMel/Penn match CSV: ${umel_match_csv}"
echo "Output directory: ${out_dir}"
echo "Python script: ${python_script}"
echo "Outcome: ${outcome}"
echo "Started at: $(date)"
echo "============================================================"

source activate survival_clock

python "${python_script}" \
  --input_long_tsv "${input_long_tsv}" \
  --med_long_tsv "${med_long_tsv}" \
  --participant_cluster_tsv "${participant_cluster_tsv}" \
  --umel_death_xlsx "${umel_death_xlsx}" \
  --umel_match_csv "${umel_match_csv}" \
  --out_dir "${out_dir}" \
  --outcome "${outcome}" \
  --robust_cov HC3 \
  --min_n_exposure 20

conda deactivate

echo "============================================================"
echo "Finished at: $(date)"
echo "============================================================"