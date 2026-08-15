#!/bin/bash
#SBATCH --partition=bioinformatics
#SBATCH --job-name=CRA_CUA_DTI
#SBATCH --array=0-47
#SBATCH --mem-per-cpu=12G
#SBATCH --time=0-02:59:00
#SBATCH --output=/cbica/home/wenju/output/CRA_CUA_DTI_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/CRA_CUA_DTI_%A_%a.err

###############################################################################
# Brain-wide association:
#   Brain-proteomics EPOCH-BAG candidate resilient agers (CRA)
#   versus concordant unfavorable agers (CUA)
#
# Outcome coding:
#   CRA = 1
#   CUA = 0
#
# Model in fit_CRA_vs_CUA_logistic.R:
#   CRA_case ~ standardized imaging IDP + covariates
#
# Therefore:
#   positive beta / OR > 1 = higher imaging feature associated with CRA
#   negative beta / OR < 1 = higher imaging feature associated with CUA
#
# Primary categorical threshold:
#   p20, using aging_phenotype_p20
#
# This DTI submit script contains 192 IDPs, so the correct array is 0-191.
# The uploaded older script used 0-192, which had one out-of-range task.
#
# The R script itself is modality-agnostic. For T1 or fMRI, keep the same R
# script and change only the imaging TSV, IDP list, output directory, and
# modality-specific covariates in a new submit script.
###############################################################################

# -----------------------------------------------------------------------------
# DTI IDPs from the uploaded script
# -----------------------------------------------------------------------------

idp_array=(
  "mean_fa_in_middle_cerebellar_peduncle_on_fa_skeleton_f25056_2_0"
  "mean_fa_in_pontine_crossing_tract_on_fa_skeleton_f25057_2_0"
  "mean_fa_in_genu_of_corpus_callosum_on_fa_skeleton_f25058_2_0"
  "mean_fa_in_body_of_corpus_callosum_on_fa_skeleton_f25059_2_0"
  "mean_fa_in_splenium_of_corpus_callosum_on_fa_skeleton_f25060_2_0"
  "mean_fa_in_fornix_on_fa_skeleton_f25061_2_0"
  "mean_fa_in_corticospinal_tract_on_fa_skeleton_right_f25062_2_0"
  "mean_fa_in_corticospinal_tract_on_fa_skeleton_left_f25063_2_0"
  "mean_fa_in_medial_lemniscus_on_fa_skeleton_right_f25064_2_0"
  "mean_fa_in_medial_lemniscus_on_fa_skeleton_left_f25065_2_0"
  "mean_fa_in_inferior_cerebellar_peduncle_on_fa_skeleton_right_f25066_2_0"
  "mean_fa_in_inferior_cerebellar_peduncle_on_fa_skeleton_left_f25067_2_0"
  "mean_fa_in_superior_cerebellar_peduncle_on_fa_skeleton_right_f25068_2_0"
  "mean_fa_in_superior_cerebellar_peduncle_on_fa_skeleton_left_f25069_2_0"
  "mean_fa_in_cerebral_peduncle_on_fa_skeleton_right_f25070_2_0"
  "mean_fa_in_cerebral_peduncle_on_fa_skeleton_left_f25071_2_0"
  "mean_fa_in_anterior_limb_of_internal_capsule_on_fa_skeleton_right_f25072_2_0"
  "mean_fa_in_anterior_limb_of_internal_capsule_on_fa_skeleton_left_f25073_2_0"
  "mean_fa_in_posterior_limb_of_internal_capsule_on_fa_skeleton_right_f25074_2_0"
  "mean_fa_in_posterior_limb_of_internal_capsule_on_fa_skeleton_left_f25075_2_0"
  "mean_fa_in_retrolenticular_part_of_internal_capsule_on_fa_skeleton_right_f25076_2_0"
  "mean_fa_in_retrolenticular_part_of_internal_capsule_on_fa_skeleton_left_f25077_2_0"
  "mean_fa_in_anterior_corona_radiata_on_fa_skeleton_right_f25078_2_0"
  "mean_fa_in_anterior_corona_radiata_on_fa_skeleton_left_f25079_2_0"
  "mean_fa_in_superior_corona_radiata_on_fa_skeleton_right_f25080_2_0"
  "mean_fa_in_superior_corona_radiata_on_fa_skeleton_left_f25081_2_0"
  "mean_fa_in_posterior_corona_radiata_on_fa_skeleton_right_f25082_2_0"
  "mean_fa_in_posterior_corona_radiata_on_fa_skeleton_left_f25083_2_0"
  "mean_fa_in_posterior_thalamic_radiation_on_fa_skeleton_right_f25084_2_0"
  "mean_fa_in_posterior_thalamic_radiation_on_fa_skeleton_left_f25085_2_0"
  "mean_fa_in_sagittal_stratum_on_fa_skeleton_right_f25086_2_0"
  "mean_fa_in_sagittal_stratum_on_fa_skeleton_left_f25087_2_0"
  "mean_fa_in_external_capsule_on_fa_skeleton_right_f25088_2_0"
  "mean_fa_in_external_capsule_on_fa_skeleton_left_f25089_2_0"
  "mean_fa_in_cingulum_cingulate_gyrus_on_fa_skeleton_right_f25090_2_0"
  "mean_fa_in_cingulum_cingulate_gyrus_on_fa_skeleton_left_f25091_2_0"
  "mean_fa_in_cingulum_hippocampus_on_fa_skeleton_right_f25092_2_0"
  "mean_fa_in_cingulum_hippocampus_on_fa_skeleton_left_f25093_2_0"
  "mean_fa_in_fornix_cresstria_terminalis_on_fa_skeleton_right_f25094_2_0"
  "mean_fa_in_fornix_cresstria_terminalis_on_fa_skeleton_left_f25095_2_0"
  "mean_fa_in_superior_longitudinal_fasciculus_on_fa_skeleton_right_f25096_2_0"
  "mean_fa_in_superior_longitudinal_fasciculus_on_fa_skeleton_left_f25097_2_0"
  "mean_fa_in_superior_frontooccipital_fasciculus_on_fa_skeleton_right_f25098_2_0"
  "mean_fa_in_superior_frontooccipital_fasciculus_on_fa_skeleton_left_f25099_2_0"
  "mean_fa_in_uncinate_fasciculus_on_fa_skeleton_right_f25100_2_0"
  "mean_fa_in_uncinate_fasciculus_on_fa_skeleton_left_f25101_2_0"
  "mean_fa_in_tapetum_on_fa_skeleton_right_f25102_2_0"
  "mean_fa_in_tapetum_on_fa_skeleton_left_f25103_2_0"
)

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "ERROR: SLURM_ARRAY_TASK_ID is not set." >&2
  exit 1
fi

if (( SLURM_ARRAY_TASK_ID < 0 || SLURM_ARRAY_TASK_ID >= ${#idp_array[@]} )); then
  echo "ERROR: SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} is outside 0-$((${#idp_array[@]} - 1))." >&2
  exit 1
fi

idp="${idp_array[$SLURM_ARRAY_TASK_ID]}"

# -----------------------------------------------------------------------------
# Input/output paths
# -----------------------------------------------------------------------------

PHENOTYPE_TSV="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/Brain_proteomics_mortality_clock/brain_proteomics_BAG_EPOCH_discordance_resilience.tsv"

IDP_TSV="/cbica/home/wenju/Reproducibile_paper/BrainAge/data/imaging/dmri_fullmetric.tsv"

# Retained for optional clinical/sensitivity covariates.
COV_TSV="/cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv"

OUTPUT_DIR="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/BWAS/Brain_proteomics_EPOCH_BAG_resilience/CRA_vs_CUA/DTI"

R_SCRIPT="/cbica/home/wenju/Project/whole-body_clocks/reslient_ager_bag_epoch/BWAS/fit_logistic.R"

mkdir -p "${OUTPUT_DIR}"
mkdir -p "/cbica/home/wenju/output"

# -----------------------------------------------------------------------------
# Phenotype definition
# -----------------------------------------------------------------------------

PHENOTYPE_COL="aging_phenotype_p20"
CASE_LABEL="CRA"
CONTROL_LABEL="CUA"

# -----------------------------------------------------------------------------
# Covariates
# -----------------------------------------------------------------------------
#
# PRIMARY model:
#   age_at_imaging + sex + Brain_ProtBAG_z
#
# age_at_imaging and sex are retained in the brain-proteomics EPOCH phenotype
# table, and Brain_ProtBAG_z further controls residual within-high-BAG
# differences between CRA and CUA.
#
# We intentionally do NOT include BMI, blood pressure, smoking, lifestyle,
# etc. in the primary model because these can themselves be mechanisms or
# correlates of resilience. They can be added in a sensitivity analysis.
#
# If a true MRI scanning-site variable is available, add it as a FACTOR
# covariate. Do not substitute a baseline assessment-centre field unless it
# actually represents the imaging site.
#
# For later T1 gray-matter analyses:
#   add intracranial volume / head-size to IDP_COVARIATES if available in the
#   imaging TSV.
#
# For later fMRI analyses:
#   add relevant motion / imaging-QC variables to IDP_COVARIATES.
# -----------------------------------------------------------------------------

PHENOTYPE_COVARIATES="age_at_imaging,sex,Brain_ProtBAG_z"

# Covariates read from COV_TSV. Empty for the primary model.
COVARIATES=""

# Modality-specific covariates read from IDP_TSV. Empty for DTI primary model.
IDP_COVARIATES=""

# Any categorical covariates included above should also be named here.
FACTOR_COVARIATES="sex"

# -----------------------------------------------------------------------------
# QC
# -----------------------------------------------------------------------------

OUTLIER_SD="4"
MIN_TOTAL_N="30"
MIN_PER_GROUP="10"

# -----------------------------------------------------------------------------
# R environment
# -----------------------------------------------------------------------------

module load R/4.2.2

# -----------------------------------------------------------------------------
# Run
# -----------------------------------------------------------------------------

comparison="${CASE_LABEL}_vs_${CONTROL_LABEL}"
output_file="${OUTPUT_DIR}/${comparison}_logistic_results_${idp}.tsv"

echo "============================================================"
echo "CRA/CUA DTI association"
echo "============================================================"
echo "SLURM task:       ${SLURM_ARRAY_TASK_ID} / $((${#idp_array[@]} - 1))"
echo "IDP:              ${idp}"
echo "Phenotype column: ${PHENOTYPE_COL}"
echo "Case:             ${CASE_LABEL} = 1"
echo "Control:          ${CONTROL_LABEL} = 0"
echo "Output:           ${output_file}"
echo "============================================================"

if [[ -s "${output_file}" ]]; then
  echo "Result already exists; skipping:"
  echo "  ${output_file}"
  exit 0
fi

Rscript "${R_SCRIPT}" \
  "${PHENOTYPE_TSV}" \
  "${IDP_TSV}" \
  "${COV_TSV}" \
  "${OUTPUT_DIR}" \
  "${idp}" \
  "${PHENOTYPE_COL}" \
  "${CASE_LABEL}" \
  "${CONTROL_LABEL}" \
  "${COVARIATES}" \
  "${FACTOR_COVARIATES}" \
  "${PHENOTYPE_COVARIATES}" \
  "${IDP_COVARIATES}" \
  "${OUTLIER_SD}" \
  "${MIN_TOTAL_N}" \
  "${MIN_PER_GROUP}"

echo "Finished:"
echo "  ${idp}"