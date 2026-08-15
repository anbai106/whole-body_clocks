#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=endocrine_epoch_horizons
#SBATCH --mem-per-cpu=24G
#SBATCH --time=11:59:00
#SBATCH --output=/cbica/home/wenju/output/endocrine_epoch_horizons_%j.out
#SBATCH --error=/cbica/home/wenju/output/endocrine_epoch_horizons_%j.err

source activate survival_clock

organ="Endocrine"
outdir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/${organ}_metabolomics_mortality_horizon_clocks"
mkdir -p "${outdir}"

# Keep the detailed run log with the analysis outputs as well as the initial
# scheduler log defined above.
exec > "${outdir}/${organ}_metabolomics_mortality_horizon_clocks_${SLURM_JOB_ID}.out" \
     2> "${outdir}/${organ}_metabolomics_mortality_horizon_clocks_${SLURM_JOB_ID}.err"

echo "============================================================"
echo "Starting Endocrine metabolomics EPOCH mortality-horizon experiment"
echo "Organ: ${organ}"
echo "Training horizons: 5 years, 10 years, full available follow-up"
echo "Common mortality evaluation times: 5 and 10 years"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}"
echo "Output directory: ${outdir}"
echo "Started at: $(date)"
echo "============================================================"

python /cbica/home/wenju/Project/whole-body_clocks/mortality_clock/Metabolomics/influence_of_time_horizon/endocrine_metabolomics_horizon_clocks.py \
  --death-xlsx /cbica/home/wenju/Dataset/UKBB_UMelbourne/Death_related_var_from_Ye.xlsx \
  --id-match-csv /cbica/home/wenju/Dataset/UKBB_UMelbourne/UKB_UMelbourne_vs_Penn_match_key.csv \
  --organ-tsv "/cbica/home/wenju/Reproducibile_paper/UKBB_metabolomics/MLNI/data/${organ}/PT/patient_pop_non_derived.tsv,/cbica/home/wenju/Reproducibile_paper/UKBB_metabolomics/MLNI/data/${organ}/test/ind_test_5000_non_derived.tsv,/cbica/home/wenju/Reproducibile_paper/UKBB_metabolomics/MLNI/data/${organ}/training/training_28142_non_derived.tsv" \
  --covariate-csv /cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv \
  --admin-censor-date 2022-11-30 \
  --outdir "${outdir}" \
  --organ "${organ}" \
  --feature-start-column diagnosis \
  --horizons 5,10,full \
  --evaluation-times 5,10 \
  --split-stratify-horizon 5 \
  --test-size 0.20 \
  --validation-size 0.20 \
  --random-state 2026 \
  --stratify-age-bins 5 \
  --max-feature-missing 0.20 \
  --l1-ratios 0.1,0.25,0.5,0.75,1.0 \
  --n-alphas 100 \
  --min-followup-days 1 \
  --final-alpha-backoff-multipliers 1,2,5,10 \
  --n-bootstrap-comparison 1000 \
  --n-calibration-groups 10 \
  --ibs-grid-points 30 \
  --ibs-start-years 0.5 || {
    status=$?
    echo "============================================================"
    echo "ERROR: Endocrine metabolomics EPOCH analysis failed"
    echo "Python exit code: ${status}"
    echo "Failed at: $(date)"
    echo "============================================================"
    conda deactivate || true
    exit "${status}"
  }

echo "============================================================"
echo "SUCCESS: Finished Endocrine metabolomics EPOCH mortality-horizon experiment"
echo "Finished at: $(date)"
echo "============================================================"

conda deactivate