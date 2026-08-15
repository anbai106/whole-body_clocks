#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=CRA_CUA_T1
#SBATCH --array=0-118
#SBATCH --mem-per-cpu=12G
#SBATCH --time=0-02:59:00
#SBATCH --output=/cbica/home/wenju/output/CRA_CUA_T1_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/CRA_CUA_T1_%A_%a.err

###############################################################################
# T1 gray-matter brain-wide association:
#   Brain-proteomics EPOCH-BAG candidate resilient agers (CRA)
#   versus concordant unfavorable agers (CUA)
#
# Outcome coding:
#   CRA = 1
#   CUA = 0
#
# Positive IDP beta / OR > 1:
#   higher regional gray-matter volume is associated with greater odds of CRA.
#
# Negative IDP beta / OR < 1:
#   higher regional gray-matter volume is associated with greater odds of CUA.
#
# Primary categorical phenotype:
#   p20, using CRA_p20 versus CUA_p20. The R script constructs this binary
#   phenotype automatically from the resilience TSV.
###############################################################################

# -----------------------------------------------------------------------------
# T1 MUSE gray-matter IDPs
# -----------------------------------------------------------------------------

idp_array=(
  "MUSE_Volume_23"
  "MUSE_Volume_30"
  "MUSE_Volume_31"
  "MUSE_Volume_32"
  "MUSE_Volume_36"
  "MUSE_Volume_37"
  "MUSE_Volume_38"
  "MUSE_Volume_39"
  "MUSE_Volume_47"
  "MUSE_Volume_48"
  "MUSE_Volume_55"
  "MUSE_Volume_56"
  "MUSE_Volume_57"
  "MUSE_Volume_58"
  "MUSE_Volume_59"
  "MUSE_Volume_60"
  "MUSE_Volume_71"
  "MUSE_Volume_72"
  "MUSE_Volume_73"
  "MUSE_Volume_75"
  "MUSE_Volume_76"
  "MUSE_Volume_100"
  "MUSE_Volume_101"
  "MUSE_Volume_102"
  "MUSE_Volume_103"
  "MUSE_Volume_104"
  "MUSE_Volume_105"
  "MUSE_Volume_106"
  "MUSE_Volume_107"
  "MUSE_Volume_108"
  "MUSE_Volume_109"
  "MUSE_Volume_112"
  "MUSE_Volume_113"
  "MUSE_Volume_114"
  "MUSE_Volume_115"
  "MUSE_Volume_116"
  "MUSE_Volume_117"
  "MUSE_Volume_118"
  "MUSE_Volume_119"
  "MUSE_Volume_120"
  "MUSE_Volume_121"
  "MUSE_Volume_122"
  "MUSE_Volume_123"
  "MUSE_Volume_124"
  "MUSE_Volume_125"
  "MUSE_Volume_128"
  "MUSE_Volume_129"
  "MUSE_Volume_132"
  "MUSE_Volume_133"
  "MUSE_Volume_134"
  "MUSE_Volume_135"
  "MUSE_Volume_136"
  "MUSE_Volume_137"
  "MUSE_Volume_138"
  "MUSE_Volume_139"
  "MUSE_Volume_140"
  "MUSE_Volume_141"
  "MUSE_Volume_142"
  "MUSE_Volume_143"
  "MUSE_Volume_144"
  "MUSE_Volume_145"
  "MUSE_Volume_146"
  "MUSE_Volume_147"
  "MUSE_Volume_148"
  "MUSE_Volume_149"
  "MUSE_Volume_150"
  "MUSE_Volume_151"
  "MUSE_Volume_152"
  "MUSE_Volume_153"
  "MUSE_Volume_154"
  "MUSE_Volume_155"
  "MUSE_Volume_156"
  "MUSE_Volume_157"
  "MUSE_Volume_160"
  "MUSE_Volume_161"
  "MUSE_Volume_162"
  "MUSE_Volume_163"
  "MUSE_Volume_164"
  "MUSE_Volume_165"
  "MUSE_Volume_166"
  "MUSE_Volume_167"
  "MUSE_Volume_168"
  "MUSE_Volume_169"
  "MUSE_Volume_170"
  "MUSE_Volume_171"
  "MUSE_Volume_172"
  "MUSE_Volume_173"
  "MUSE_Volume_174"
  "MUSE_Volume_175"
  "MUSE_Volume_176"
  "MUSE_Volume_177"
  "MUSE_Volume_178"
  "MUSE_Volume_179"
  "MUSE_Volume_180"
  "MUSE_Volume_181"
  "MUSE_Volume_182"
  "MUSE_Volume_183"
  "MUSE_Volume_184"
  "MUSE_Volume_185"
  "MUSE_Volume_186"
  "MUSE_Volume_187"
  "MUSE_Volume_190"
  "MUSE_Volume_191"
  "MUSE_Volume_192"
  "MUSE_Volume_193"
  "MUSE_Volume_194"
  "MUSE_Volume_195"
  "MUSE_Volume_196"
  "MUSE_Volume_197"
  "MUSE_Volume_198"
  "MUSE_Volume_199"
  "MUSE_Volume_200"
  "MUSE_Volume_201"
  "MUSE_Volume_202"
  "MUSE_Volume_203"
  "MUSE_Volume_204"
  "MUSE_Volume_205"
  "MUSE_Volume_206"
  "MUSE_Volume_207"
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

IDP_TSV="/cbica/home/wenju/Reproducibile_paper/BrainAge/data/imaging/T1_MUSE_GM.tsv"

COV_TSV="/cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv"

OUTPUT_DIR="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/BWAS/Brain_proteomics_EPOCH_BAG_resilience/CRA_vs_CUA/T1"

# Keep the R script name exactly as requested.
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
# Primary model:
#   CRA vs CUA ~ T1_IDP_z + age_at_imaging + sex + Brain_ProtBAG_z
#
# age_at_imaging, sex, and Brain_ProtBAG_z are read directly from the
# resilience phenotype TSV.
#
# IMPORTANT FOR T1:
# If T1_MUSE_GM.tsv contains a total intracranial volume / head-size variable,
# add its exact column name to IDP_COVARIATES below, e.g.
#
#   IDP_COVARIATES="your_exact_ICV_column"
#
# Do not invent an ICV column name; leave this empty unless the exact variable
# exists in T1_MUSE_GM.tsv.
# -----------------------------------------------------------------------------

PHENOTYPE_COVARIATES="age_at_imaging,sex,Brain_ProtBAG_z"

# Optional general covariates read from COV_TSV.
COVARIATES=""

# Optional T1-specific covariates read from IDP_TSV.
IDP_COVARIATES=""

# Categorical covariates included in the model.
FACTOR_COVARIATES="sex"

# -----------------------------------------------------------------------------
# QC
# -----------------------------------------------------------------------------

OUTLIER_SD="4"
MIN_TOTAL_N="100"
MIN_PER_GROUP="30"

# -----------------------------------------------------------------------------
# R environment
# -----------------------------------------------------------------------------

module load R/4.2.2

# -----------------------------------------------------------------------------
# Run one T1 IDP
# -----------------------------------------------------------------------------

comparison="${CASE_LABEL}_vs_${CONTROL_LABEL}"
output_file="${OUTPUT_DIR}/${comparison}_logistic_results_${idp}.tsv"

echo "============================================================"
echo "CRA/CUA T1 gray-matter association"
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
