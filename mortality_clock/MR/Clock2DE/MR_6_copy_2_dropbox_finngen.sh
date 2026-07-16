#!/usr/bin/env bash
set -euo pipefail

REMOTE_USER="wenju"
REMOTE_HOST="cubic-login.uphs.upenn.edu"
REMOTE_ROOT="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/MR/Clock2DE/FinnGen"

LOCAL_ROOT="/Users/hao/Dropbox/2026_EPOCH/Supplementary_Dataset_MR/EPOCH2DE/FinnGen"

SCRIPT_NAME="$(basename "$0")"

mkdir -p "${LOCAL_ROOT}"

echo "Remote source:"
echo "  ${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_ROOT}"
echo
echo "Local destination:"
echo "  ${LOCAL_ROOT}"
echo

DRY_RUN="${DRY_RUN:-1}"

RSYNC_FLAGS="-avhm --progress --prune-empty-dirs"

if [[ "${DRY_RUN}" == "1" ]]; then
  RSYNC_FLAGS="${RSYNC_FLAGS} --dry-run"
  echo "Running in DRY-RUN mode. No files will be copied."
  echo "To actually copy, run:"
  echo "  DRY_RUN=0 bash ${SCRIPT_NAME}"
  echo
else
  echo "Running in COPY mode."
  echo
fi

rsync ${RSYNC_FLAGS} \
  --include='*/' \
  --include='*/QC/SC_*' \
  --exclude='*' \
  "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_ROOT}/" \
  "${LOCAL_ROOT}/"

echo
echo "Done."

echo
echo "Number of local SC_* files:"
find "${LOCAL_ROOT}" -path "*/QC/SC_*" -type f | wc -l

echo
echo "Subfolders copied:"
find "${LOCAL_ROOT}" -mindepth 1 -maxdepth 1 -type d | sort