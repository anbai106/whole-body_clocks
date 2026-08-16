#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=endocrine_horizon_disease
#SBATCH --array=0-4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=0-08:00:00
#SBATCH --output=/cbica/home/wenju/output/endocrine_horizon_disease_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/endocrine_horizon_disease_%A_%a.err

module load python/anaconda/3
source activate survival

# The five unique algorithmically-defined disease endpoints used in the
# previous longitudinal metabolomics analysis.
endpoints=(
  all_cause_dementia
  asthma
  myocardial_infarction
  copd
  stroke
)

n_endpoints=${#endpoints[@]}
if [[ "${SLURM_ARRAY_TASK_ID}" -ge "${n_endpoints}" ]]; then
  echo "ERROR: SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} but only ${n_endpoints} endpoints are configured."
  exit 1
fi

endpoint="${endpoints[$SLURM_ARRAY_TASK_ID]}"

python_script=${PYTHON_SCRIPT:-/cbica/home/wenju/Project/whole-body_clocks/mortality_clock/Metabolomics/influence_of_time_horizon/survival_analysis_endocrine_horizon_clocks_disease_onset.py}

horizon_predictions=${HORIZON_PREDICTIONS:-/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/Endocrine_metabolomics_mortality_horizon_clocks/endocrine_metabolomics_mortality_horizon_clocks_predictions.tsv}

umel_death_xlsx=${UMEL_DEATH_XLSX:-/cbica/home/wenju/Dataset/UKBB_UMelbourne/Death_related_var_from_Ye.xlsx}
umel_match_csv=${UMEL_MATCH_CSV:-/cbica/home/wenju/Dataset/UKBB_UMelbourne/UKB_UMelbourne_vs_Penn_match_key.csv}
cov_tsv=${COV_TSV:-/cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv}

# Primary default: held-out test participants from mortality-EPOCH development.
# For a larger sensitivity analysis, submit with:
#   ANALYSIS_SPLIT=all sbatch submit_endocrine_horizon_clocks_disease_onset.sh
analysis_split=${ANALYSIS_SPLIT:-test}

n_bootstrap=${N_BOOTSTRAP:-1000}
force_rerun=${FORCE_RERUN:-1}

out_dir=${OUT_DIR:-/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/Endocrine_metabolomics_mortality_horizon_clocks/disease_onset_${analysis_split}}
mkdir -p /cbica/home/wenju/output
mkdir -p "${out_dir}"

output_tsv="${out_dir}/cox_endocrine_horizon_clocks_${endpoint}_${analysis_split}.tsv"

# Redirect detailed task log into the analysis folder as well as preserving the
# top-level SBATCH stdout/stderr destination until this point.
exec > "${out_dir}/endocrine_horizon_disease_${endpoint}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out" \
     2> "${out_dir}/endocrine_horizon_disease_${endpoint}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"

echo "============================================================"
echo "Endocrine mortality-EPOCH horizon -> disease onset analysis"
echo "Endpoint: ${endpoint}"
echo "Analysis split: ${analysis_split}"
echo "Horizon predictions: ${horizon_predictions}"
echo "Output TSV: ${output_tsv}"
echo "Bootstrap replicates: ${n_bootstrap}"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "Started at: $(date)"
echo "============================================================"

if [[ "${force_rerun}" != "1" && -s "${output_tsv}" ]]; then
  echo "Existing non-empty output found and FORCE_RERUN=${force_rerun}; skipping."
  conda deactivate || true
  exit 0
fi

if [[ "${force_rerun}" == "1" ]]; then
  prefix="${output_tsv%.tsv}"
  rm -f "${output_tsv}" \
        "${prefix}_paired_delta_cindex.tsv" \
        "${prefix}_clock_correlations.tsv" \
        "${prefix}_analysis_dataset.tsv" \
        "${prefix}_common_complete_cases.tsv" \
        "${prefix}_run_summary.json" \
        "${prefix}_horizon_clock_forest_plot.png" \
        "${prefix}_horizon_clock_forest_plot.pdf" \
        "${prefix}_horizon_clock_forest_plot.svg"
fi

set +e
python "${python_script}" \
  --endpoint "${endpoint}" \
  --output_tsv "${output_tsv}" \
  --horizon_predictions "${horizon_predictions}" \
  --analysis_split "${analysis_split}" \
  --umel_death_xlsx "${umel_death_xlsx}" \
  --umel_match_csv "${umel_match_csv}" \
  --cov_tsv "${cov_tsv}" \
  --admin_censor_date 2022-11-30 \
  --include_bp \
  --n_bootstrap "${n_bootstrap}"
status=$?
set -e

if [[ "${status}" -ne 0 ]]; then
  echo "============================================================"
  echo "ERROR: disease-onset analysis failed"
  echo "Endpoint: ${endpoint}"
  echo "Python exit code: ${status}"
  echo "Failed at: $(date)"
  echo "============================================================"
  conda deactivate || true
  exit "${status}"
fi

echo "============================================================"
echo "SUCCESS: finished Endocrine horizon-clock disease-onset analysis"
echo "Endpoint: ${endpoint}"
echo "Primary result: ${output_tsv}"
echo "Finished at: $(date)"
echo "============================================================"

conda deactivate || true
