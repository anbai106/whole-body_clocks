#!/usr/bin/env bash

set -euo pipefail

REMOTE_USER_HOST="wenju@cubic-login.uphs.upenn.edu"
REMOTE_BASE="Reproducibile_paper/WholeBodyClock/mortality_clock/fuma"
SOURCE_FILE="EPOCH_open_circular_R.svg"

LOCAL_DEST="/Users/hao/Dropbox/2026_EPOCH/Fig/orig/circular_manhattan_mortality_epoch"

mkdir -p "${LOCAL_DEST}"

echo "Searching for SVG files on ${REMOTE_USER_HOST}..."

mapfile_output="$(
    ssh "${REMOTE_USER_HOST}" \
        "find '${REMOTE_BASE}' -mindepth 2 -maxdepth 2 \
        -type f \
        -path '*mortality_clock/${SOURCE_FILE}' \
        -print | sort"
)"

if [[ -z "${mapfile_output}" ]]; then
    echo "ERROR: No matching SVG files were found."
    exit 1
fi

file_count="$(printf '%s\n' "${mapfile_output}" | sed '/^$/d' | wc -l | tr -d ' ')"

echo "Found ${file_count} SVG files."

copied=0
failed=0

while IFS= read -r remote_file; do
    [[ -z "${remote_file}" ]] && continue

    # Extract the parent directory, which is the clock name.
    clock_name="$(basename "$(dirname "${remote_file}")")"

    local_file="${LOCAL_DEST}/${clock_name}.svg"

    echo "Copying:"
    echo "  ${remote_file}"
    echo "  -> ${local_file}"

    if scp "${REMOTE_USER_HOST}:${remote_file}" "${local_file}"; then
        copied=$((copied + 1))
    else
        echo "WARNING: Failed to copy ${remote_file}" >&2
        failed=$((failed + 1))
    fi
done <<< "${mapfile_output}"

echo
echo "Copy completed."
echo "Successfully copied: ${copied}"
echo "Failed:              ${failed}"
echo "Destination:         ${LOCAL_DEST}"

actual_count="$(
    find "${LOCAL_DEST}" -maxdepth 1 -type f -name '*mortality_clock.svg' |
        wc -l |
        tr -d ' '
)"

echo "Matching SVG files now in destination: ${actual_count}"

if [[ "${failed}" -gt 0 ]]; then
    exit 1
fi

if [[ "${copied}" -ne "${file_count}" ]]; then
    echo "WARNING: The number copied does not match the number found." >&2
    exit 1
fi