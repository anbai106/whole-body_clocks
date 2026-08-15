#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=CRA_CUA_FC
#SBATCH --array=0-209
#SBATCH --mem-per-cpu=12G
#SBATCH --time=0-02:59:00
#SBATCH --output=/cbica/home/wenju/output/CRA_CUA_FC_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/CRA_CUA_FC_%A_%a.err

###############################################################################
# Functional MRI brain-wide association:
#   Brain-proteomics EPOCH-BAG candidate resilient agers (CRA)
#   versus concordant unfavorable agers (CUA)
#
# Functional-connectivity features:
#   f_1 ... f_210 from fmri_25_component.tsv
#
# Outcome coding:
#   CRA = 1
#   CUA = 0
#
# Positive IDP beta / OR > 1:
#   higher functional-connectivity feature value is associated with greater
#   odds of CRA.
#
# Negative IDP beta / OR < 1:
#   higher functional-connectivity feature value is associated with greater
#   odds of CUA.
#
# Primary categorical phenotype:
#   p25, using CRA_p25 versus CUA_p25. The R script constructs this binary
#   phenotype automatically from the resilience TSV.
#
# There are 210 FC features, so valid SLURM indices are 0-209.
###############################################################################

# -----------------------------------------------------------------------------
# Functional MRI features
# -----------------------------------------------------------------------------

idp_array=(
  "f_1"
  "f_2"
  "f_3"
  "f_4"
  "f_5"
  "f_6"
  "f_7"
  "f_8"
  "f_9"
  "f_10"
  "f_11"
  "f_12"
  "f_13"
  "f_14"
  "f_15"
  "f_16"
  "f_17"
  "f_18"
  "f_19"
  "f_20"
  "f_21"
  "f_22"
  "f_23"
  "f_24"
  "f_25"
  "f_26"
  "f_27"
  "f_28"
  "f_29"
  "f_30"
  "f_31"
  "f_32"
  "f_33"
  "f_34"
  "f_35"
  "f_36"
  "f_37"
  "f_38"
  "f_39"
  "f_40"
  "f_41"
  "f_42"
  "f_43"
  "f_44"
  "f_45"
  "f_46"
  "f_47"
  "f_48"
  "f_49"
  "f_50"
  "f_51"
  "f_52"
  "f_53"
  "f_54"
  "f_55"
  "f_56"
  "f_57"
  "f_58"
  "f_59"
  "f_60"
  "f_61"
  "f_62"
  "f_63"
  "f_64"
  "f_65"
  "f_66"
  "f_67"
  "f_68"
  "f_69"
  "f_70"
  "f_71"
  "f_72"
  "f_73"
  "f_74"
  "f_75"
  "f_76"
  "f_77"
  "f_78"
  "f_79"
  "f_80"
  "f_81"
  "f_82"
  "f_83"
  "f_84"
  "f_85"
  "f_86"
  "f_87"
  "f_88"
  "f_89"
  "f_90"
  "f_91"
  "f_92"
  "f_93"
  "f_94"
  "f_95"
  "f_96"
  "f_97"
  "f_98"
  "f_99"
  "f_100"
  "f_101"
  "f_102"
  "f_103"
  "f_104"
  "f_105"
  "f_106"
  "f_107"
  "f_108"
  "f_109"
  "f_110"
  "f_111"
  "f_112"
  "f_113"
  "f_114"
  "f_115"
  "f_116"
  "f_117"
  "f_118"
  "f_119"
  "f_120"
  "f_121"
  "f_122"
  "f_123"
  "f_124"
  "f_125"
  "f_126"
  "f_127"
  "f_128"
  "f_129"
  "f_130"
  "f_131"
  "f_132"
  "f_133"
  "f_134"
  "f_135"
  "f_136"
  "f_137"
  "f_138"
  "f_139"
  "f_140"
  "f_141"
  "f_142"
  "f_143"
  "f_144"
  "f_145"
  "f_146"
  "f_147"
  "f_148"
  "f_149"
  "f_150"
  "f_151"
  "f_152"
  "f_153"
  "f_154"
  "f_155"
  "f_156"
  "f_157"
  "f_158"
  "f_159"
  "f_160"
  "f_161"
  "f_162"
  "f_163"
  "f_164"
  "f_165"
  "f_166"
  "f_167"
  "f_168"
  "f_169"
  "f_170"
  "f_171"
  "f_172"
  "f_173"
  "f_174"
  "f_175"
  "f_176"
  "f_177"
  "f_178"
  "f_179"
  "f_180"
  "f_181"
  "f_182"
  "f_183"
  "f_184"
  "f_185"
  "f_186"
  "f_187"
  "f_188"
  "f_189"
  "f_190"
  "f_191"
  "f_192"
  "f_193"
  "f_194"
  "f_195"
  "f_196"
  "f_197"
  "f_198"
  "f_199"
  "f_200"
  "f_201"
  "f_202"
  "f_203"
  "f_204"
  "f_205"
  "f_206"
  "f_207"
  "f_208"
  "f_209"
  "f_210"
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

IDP_TSV="/cbica/home/wenju/Reproducibile_paper/BrainAge/data/imaging/fmri_25_component.tsv"

COV_TSV="/cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv"

OUTPUT_DIR="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/BWAS/Brain_proteomics_EPOCH_BAG_resilience/CRA_vs_CUA/FC"

# Keep the same generalized R script.
R_SCRIPT="/cbica/home/wenju/Project/whole-body_clocks/reslient_ager_bag_epoch/BWAS/fit_logistic.R"

mkdir -p "${OUTPUT_DIR}"
mkdir -p "/cbica/home/wenju/output"

# -----------------------------------------------------------------------------
# Phenotype definition
# -----------------------------------------------------------------------------

PHENOTYPE_COL="aging_phenotype_p25"
CASE_LABEL="CRA"
CONTROL_LABEL="CUA"

# -----------------------------------------------------------------------------
# Covariates
# -----------------------------------------------------------------------------
#
# Primary model:
#   CRA vs CUA ~ FC_IDP_z + age_at_imaging + sex + Brain_ProtBAG_z
#
# age_at_imaging, sex, and Brain_ProtBAG_z are read from the resilience TSV.
#
# IMPORTANT FOR fMRI:
# Motion and imaging-quality covariates can materially affect functional
# connectivity. If fmri_25_component.tsv contains appropriate motion/QC
# variables, add their exact column names to IDP_COVARIATES, separated by
# commas, for example:
#
#   IDP_COVARIATES="mean_framewise_displacement,another_qc_variable"
#
# Do not invent variable names. Leave this empty unless those exact columns
# are present in the imaging TSV.
# -----------------------------------------------------------------------------

PHENOTYPE_COVARIATES="age_at_imaging,sex,Brain_ProtBAG_z"

# Optional general covariates read from COV_TSV.
COVARIATES=""

# Optional fMRI-specific motion / imaging-QC covariates read from IDP_TSV.
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
# Run one functional-connectivity feature
# -----------------------------------------------------------------------------

comparison="${CASE_LABEL}_vs_${CONTROL_LABEL}"
output_file="${OUTPUT_DIR}/${comparison}_logistic_results_${idp}.tsv"

echo "============================================================"
echo "CRA/CUA functional MRI association"
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