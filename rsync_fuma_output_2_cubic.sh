#!/bin/bash

set -euo pipefail

# ============================================================
# Local and remote settings
# ============================================================

LOCAL_FUMA_BASE="/Users/hao/Downloads/FUMA"

REMOTE_HOST="wenju@cubic-login.uphs.upenn.edu"
REMOTE_BASE="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock"

LOG_FILE="${LOCAL_FUMA_BASE}/rsync_fuma_to_cubic_manifest.tsv"

# Set DRY_RUN=1 to test without copying/unzipping.
DRY_RUN="${DRY_RUN:-0}"

# ============================================================
# Initialize manifest
# ============================================================

echo -e "disease\tlocal_zip\tzip_file\tclock_name\tremote_clock_dir\tremote_fuma_dir\tstatus" > "${LOG_FILE}"

# ============================================================
# Check local folder
# ============================================================

if [[ ! -d "${LOCAL_FUMA_BASE}" ]]; then
  echo "ERROR: Cannot find local FUMA folder: ${LOCAL_FUMA_BASE}" >&2
  exit 1
fi

N_LOCAL_ZIPS=$(find "${LOCAL_FUMA_BASE}" -mindepth 2 -maxdepth 2 -type f -name "*.zip" | wc -l | tr -d ' ')
echo "Found ${N_LOCAL_ZIPS} local FUMA zip files under ${LOCAL_FUMA_BASE}"

if [[ "${N_LOCAL_ZIPS}" -eq 0 ]]; then
  echo "ERROR: No zip files found." >&2
  exit 1
fi

# ============================================================
# Main loop
# Important:
# Use process substitution instead of a pipe, and use ssh -n,
# otherwise ssh may consume stdin and terminate the loop after 1 file.
# ============================================================

N_OK=0
N_SKIPPED=0
N_FAILED=0

while IFS= read -r LOCAL_ZIP; do

  DISEASE="$(basename "$(dirname "${LOCAL_ZIP}")")"
  ZIP_FILE="$(basename "${LOCAL_ZIP}")"

  # Remove .zip suffix.
  CLOCK_NAME="${ZIP_FILE%.zip}"

  # Handle files like Reproductive_male_proteomics_mi_clock_fuma.zip
  # so that they map to Reproductive_male_proteomics_mi_clock.
  CLOCK_NAME="${CLOCK_NAME%_fuma}"

  REMOTE_CLOCK_DIR="${REMOTE_BASE}/${CLOCK_NAME}"
  REMOTE_FUMA_DIR="${REMOTE_CLOCK_DIR}/fuma"
  REMOTE_ZIP="${REMOTE_FUMA_DIR}/${ZIP_FILE}"

  echo "============================================================"
  echo "Disease:        ${DISEASE}"
  echo "Local zip:      ${LOCAL_ZIP}"
  echo "Zip file:       ${ZIP_FILE}"
  echo "Clock name:     ${CLOCK_NAME}"
  echo "Remote clock:   ${REMOTE_CLOCK_DIR}"
  echo "Remote fuma:    ${REMOTE_FUMA_DIR}"

  if [[ "${DRY_RUN}" == "1" ]]; then
    echo "[DRY RUN] Would check remote clock folder, create fuma folder, rsync zip, and unzip."
    echo -e "${DISEASE}\t${LOCAL_ZIP}\t${ZIP_FILE}\t${CLOCK_NAME}\t${REMOTE_CLOCK_DIR}\t${REMOTE_FUMA_DIR}\tdry_run" >> "${LOG_FILE}"
    continue
  fi

  # Check that the corresponding clock folder exists on CUBIC.
  if ! ssh -n "${REMOTE_HOST}" "test -d '${REMOTE_CLOCK_DIR}'"; then
    echo "WARNING: Remote clock folder does not exist, skipping: ${REMOTE_CLOCK_DIR}" >&2
    echo -e "${DISEASE}\t${LOCAL_ZIP}\t${ZIP_FILE}\t${CLOCK_NAME}\t${REMOTE_CLOCK_DIR}\t${REMOTE_FUMA_DIR}\tmissing_remote_clock_dir" >> "${LOG_FILE}"
    N_SKIPPED=$((N_SKIPPED + 1))
    continue
  fi

  # Create fuma folder.
  ssh -n "${REMOTE_HOST}" "mkdir -p '${REMOTE_FUMA_DIR}'"

  # Copy zip file.
  if ! rsync -avh --progress "${LOCAL_ZIP}" "${REMOTE_HOST}:${REMOTE_ZIP}"; then
    echo "ERROR: rsync failed for ${LOCAL_ZIP}" >&2
    echo -e "${DISEASE}\t${LOCAL_ZIP}\t${ZIP_FILE}\t${CLOCK_NAME}\t${REMOTE_CLOCK_DIR}\t${REMOTE_FUMA_DIR}\trsync_failed" >> "${LOG_FILE}"
    N_FAILED=$((N_FAILED + 1))
    continue
  fi

  # Unzip on CUBIC inside the corresponding fuma folder.
  if ! ssh -n "${REMOTE_HOST}" "cd '${REMOTE_FUMA_DIR}' && unzip -oq '${ZIP_FILE}'"; then
    echo "ERROR: unzip failed for ${REMOTE_ZIP}" >&2
    echo -e "${DISEASE}\t${LOCAL_ZIP}\t${ZIP_FILE}\t${CLOCK_NAME}\t${REMOTE_CLOCK_DIR}\t${REMOTE_FUMA_DIR}\tunzip_failed" >> "${LOG_FILE}"
    N_FAILED=$((N_FAILED + 1))
    continue
  fi

  echo "DONE: ${CLOCK_NAME}"
  echo -e "${DISEASE}\t${LOCAL_ZIP}\t${ZIP_FILE}\t${CLOCK_NAME}\t${REMOTE_CLOCK_DIR}\t${REMOTE_FUMA_DIR}\tok" >> "${LOG_FILE}"
  N_OK=$((N_OK + 1))

done < <(find "${LOCAL_FUMA_BASE}" -mindepth 2 -maxdepth 2 -type f -name "*.zip" | sort)

echo "============================================================"
echo "All done."
echo "Local zip files found: ${N_LOCAL_ZIPS}"
echo "Copied and unzipped OK: ${N_OK}"
echo "Skipped missing remote clock dirs: ${N_SKIPPED}"
echo "Failed: ${N_FAILED}"
echo "Manifest written to:"
echo "${LOG_FILE}"