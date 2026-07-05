#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=med_delta_lm
#SBATCH --mem-per-cpu=16G
#SBATCH --time=0-04:00:00
#SBATCH --output=/cbica/home/wenju/output/med_delta_lm_%j.out
#SBATCH --error=/cbica/home/wenju/output/med_delta_lm_%j.err

set -euo pipefail

module load python/anaconda/3

mkdir -p /cbica/home/wenju/output

input_tsv=${INPUT_TSV:-/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/SA/delta_metabolomics_algorithmic_disease_onset/medication_cluster_delta_clock_inputs/metabolomics_delta_clock_medication_cluster_requested5_long.tsv}

out_dir=${OUT_DIR:-/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/SA/delta_metabolomics_algorithmic_disease_onset/medication_cluster_delta_clock_lm_results}

python_script=/cbica/home/wenju/Project/whole-body_clocks/mortality_clock/Metabolomics/longitudinal/medication/2_run_lm_across_clusters_and_delta_clocks.py

outcome=${OUTCOME:-delta_clock_age_years}

echo "============================================================"
echo "Medication-cluster delta-clock linear model"
echo "Input TSV: ${input_tsv}"
echo "Output directory: ${out_dir}"
echo "Python script: ${python_script}"
echo "Outcome: ${outcome}"
echo "Started at: $(date)"
echo "============================================================"

source activate DNE

python "${python_script}" \
  --input_tsv "${input_tsv}" \
  --out_dir "${out_dir}" \
  --outcome "${outcome}" \
  --robust_cov HC3 \
  --min_n_cluster 20 \
  --save_complete_cases

conda deactivate

echo "============================================================"
echo "Finished at: $(date)"
echo "============================================================"