#!/usr/bin/env bash

# ============================================================
# Copy disease-clock Manhattan and QQ plots from CUBIC to local Mac
#
# Remote input:
#   /cbica/home/wenju/Reproducibile_paper/WholeBodyClock/
#     *_disease_clock/fastGWA/qmplot/manhattan_qmplot.png
#     *_disease_clock/fastGWA/qmplot/QQ_plot.png
#
# Local output:
#   /Users/hao/Dropbox/2026_EPOCH/QMPLOT/disease_clock
#
# Output filenames:
#   <organ>_<modality>_<disease>_manhattan_plot.png
#   <organ>_<modality>_<disease>_qq_plot.png
# ============================================================

set -euo pipefail

# ------------------------------------------------------------
# 1. Settings
# ------------------------------------------------------------

REMOTE_HOST="wenju@cubic-login.uphs.upenn.edu"
REMOTE_BASE="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock"
LOCAL_DIR="/Users/hao/Dropbox/2026_EPOCH/QMPLOT/disease_clock"

mkdir -p "${LOCAL_DIR}"

echo "Remote host: ${REMOTE_HOST}"
echo "Remote base: ${REMOTE_BASE}"
echo "Local output: ${LOCAL_DIR}"

# ------------------------------------------------------------
# 2. Collect remote plot paths
# ------------------------------------------------------------

TMP_LIST="$(mktemp)"
TMP_MANIFEST="$(mktemp)"

ssh "${REMOTE_HOST}" "REMOTE_BASE='${REMOTE_BASE}' bash -s" > "${TMP_LIST}" <<'REMOTE_SCRIPT'
set -euo pipefail
shopt -s nullglob

diseases=(asthma copd dementia mi stroke)

for disease in "${diseases[@]}"; do
  for f in "${REMOTE_BASE}"/*_"${disease}"_clock/fastGWA/qmplot/manhattan_qmplot.png; do
    [[ -f "${f}" ]] && printf "%s\n" "${f}"
  done

  for f in "${REMOTE_BASE}"/*_"${disease}"_clock/fastGWA/qmplot/QQ_plot.png; do
    [[ -f "${f}" ]] && printf "%s\n" "${f}"
  done
done
REMOTE_SCRIPT

N_REMOTE="$(wc -l < "${TMP_LIST}" | tr -d ' ')"
echo "Detected remote PNG files: ${N_REMOTE}"

if [[ "${N_REMOTE}" -eq 0 ]]; then
  echo "ERROR: No remote plot files found. Please check REMOTE_HOST and REMOTE_BASE."
  rm -f "${TMP_LIST}" "${TMP_MANIFEST}"
  exit 1
fi

# ------------------------------------------------------------
# 3. Copy and rename files
# ------------------------------------------------------------

echo -e "remote_file\tlocal_file\tplot_type\tclock_name" > "${TMP_MANIFEST}"

while IFS= read -r REMOTE_PNG; do
  [[ -z "${REMOTE_PNG}" ]] && continue

  # Remote path structure:
  #   <clock_folder>/fastGWA/qmplot/<plot_file>
  PLOT_FILE="$(basename "${REMOTE_PNG}")"
  CLOCK_FOLDER="$(basename "$(dirname "$(dirname "$(dirname "${REMOTE_PNG}")")")")"

  # Remove final "_clock"
  # Example:
  #   Brain_proteomics_dementia_clock -> Brain_proteomics_dementia
  CLOCK_NAME="${CLOCK_FOLDER%_clock}"

  case "${PLOT_FILE}" in
    "manhattan_qmplot.png")
      OUT_FILE="${CLOCK_NAME}_manhattan_plot.png"
      PLOT_TYPE="manhattan"
      ;;
    "QQ_plot.png")
      OUT_FILE="${CLOCK_NAME}_qq_plot.png"
      PLOT_TYPE="qq"
      ;;
    *)
      echo "Skipping unrecognized plot file: ${REMOTE_PNG}"
      continue
      ;;
  esac

  LOCAL_FILE="${LOCAL_DIR}/${OUT_FILE}"

  echo "Copying:"
  echo "  ${REMOTE_PNG}"
  echo "  -> ${LOCAL_FILE}"

  scp -p "${REMOTE_HOST}:${REMOTE_PNG}" "${LOCAL_FILE}"

  echo -e "${REMOTE_PNG}\t${LOCAL_FILE}\t${PLOT_TYPE}\t${CLOCK_NAME}" >> "${TMP_MANIFEST}"

done < "${TMP_LIST}"

# ------------------------------------------------------------
# 4. Save copy manifest and print QC
# ------------------------------------------------------------

MANIFEST_OUT="${LOCAL_DIR}/copy_manifest_disease_clock_qmplots.tsv"
cp "${TMP_MANIFEST}" "${MANIFEST_OUT}"

N_MANHATTAN="$(find "${LOCAL_DIR}" -maxdepth 1 -type f -name '*_manhattan_plot.png' | wc -l | tr -d ' ')"
N_QQ="$(find "${LOCAL_DIR}" -maxdepth 1 -type f -name '*_qq_plot.png' | wc -l | tr -d ' ')"

echo ""
echo "Done."
echo "Local output folder:"
echo "  ${LOCAL_DIR}"
echo "Copy manifest:"
echo "  ${MANIFEST_OUT}"
echo ""
echo "QC counts:"
echo "  Manhattan plots: ${N_MANHATTAN}"
echo "  QQ plots:        ${N_QQ}"

rm -f "${TMP_LIST}" "${TMP_MANIFEST}"