#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Download / verify / re-download NHANES 1999-2018 data
# plus 2019 public-use linked mortality files.
#
# Output:
#   /Users/hao/Dropbox/NHANES/
#
# Normal run:
#   bash /Users/hao/Dropbox/NHANES/1_download_data.sh
#
# If curl gives SSL certificate errors, run:
#   ALLOW_INSECURE_SSL=1 bash /Users/hao/Dropbox/NHANES/1_download_data.sh
#
# This script:
#   1. Fetches expected XPT links from CDC component pages
#   2. Verifies existing files by size when remote size is available
#   3. Re-downloads missing, empty, or incomplete files
#   4. Downloads NHANES linked mortality files
#   5. Writes manifest and verification reports
# ============================================================

OUTDIR="/Users/hao/Dropbox/NHANES"
mkdir -p "${OUTDIR}"

BASE_NHANES_PAGE="https://wwwn.cdc.gov/nchs/nhanes/search/datapage.aspx"
BASE_MORT_URL="https://ftp.cdc.gov/pub/health_statistics/NCHS/datalinkage/linked_mortality"

EXPECTED_MANIFEST="${OUTDIR}/expected_download_manifest.tsv"
DOWNLOADED_REPORT="${OUTDIR}/downloaded_or_redownloaded_files.tsv"
FAILED_REPORT="${OUTDIR}/failed_downloads.tsv"
FINAL_MISSING_REPORT="${OUTDIR}/final_missing_or_incomplete_files.tsv"
SUMMARY_REPORT="${OUTDIR}/verification_summary.tsv"

echo -e "cycle\tcomponent\tfilename\turl\tlocal_path\tremote_size" > "${EXPECTED_MANIFEST}"
echo -e "cycle\tcomponent\tfilename\turl\tlocal_path\treason\told_local_size\tremote_size" > "${DOWNLOADED_REPORT}"
echo -e "cycle\tcomponent\tfilename\turl\tlocal_path\terror" > "${FAILED_REPORT}"
echo -e "cycle\tcomponent\tfilename\turl\tlocal_path\treason\tlocal_size\tremote_size" > "${FINAL_MISSING_REPORT}"

cycles=(
  "1999-2000|1999|1999_2000"
  "2001-2002|2001|2001_2002"
  "2003-2004|2003|2003_2004"
  "2005-2006|2005|2005_2006"
  "2007-2008|2007|2007_2008"
  "2009-2010|2009|2009_2010"
  "2011-2012|2011|2011_2012"
  "2013-2014|2013|2013_2014"
  "2015-2016|2015|2015_2016"
  "2017-2018|2017|2017_2018"
)

components=(
  "Demographics"
  "Questionnaire"
  "Examination"
  "Laboratory"
  "Dietary"
)

if [[ "${ALLOW_INSECURE_SSL:-0}" == "1" ]]; then
  echo "WARNING: ALLOW_INSECURE_SSL=1 is enabled."
  echo "curl will use -k and bypass SSL certificate verification."
  echo
fi

curl_common() {
  if [[ "${ALLOW_INSECURE_SSL:-0}" == "1" ]]; then
    curl -k -L --fail --retry 5 --retry-delay 3 --connect-timeout 30 "$@"
  else
    curl -L --fail --retry 5 --retry-delay 3 --connect-timeout 30 "$@"
  fi
}

curl_head_common() {
  if [[ "${ALLOW_INSECURE_SSL:-0}" == "1" ]]; then
    curl -k -sI -L --retry 3 --connect-timeout 30 "$@"
  else
    curl -sI -L --retry 3 --connect-timeout 30 "$@"
  fi
}

get_remote_size() {
  local url="$1"
  local size

  size="$(
    curl_head_common "${url}" \
      | awk 'BEGIN{IGNORECASE=1} /^content-length:/ {gsub("\r","",$2); s=$2} END{print s}'
  )"

  if [[ "${size}" =~ ^[0-9]+$ ]]; then
    echo "${size}"
  else
    echo "NA"
  fi
}

get_local_size() {
  local file="$1"

  if [[ -f "${file}" ]]; then
    wc -c < "${file}" | tr -d ' '
  else
    echo "0"
  fi
}

needs_download() {
  local outfile="$1"
  local remote_size="$2"
  local local_size

  if [[ ! -s "${outfile}" ]]; then
    return 0
  fi

  local_size="$(get_local_size "${outfile}")"

  if [[ "${remote_size}" != "NA" && "${local_size}" != "${remote_size}" ]]; then
    return 0
  fi

  return 1
}

download_one_file() {
  local cycle_label="$1"
  local component="$2"
  local url="$3"
  local outfile="$4"

  local filename
  local remote_size
  local local_size
  local tmpfile
  local reason
  local new_size

  filename="$(basename "${outfile}")"
  remote_size="$(get_remote_size "${url}")"
  local_size="$(get_local_size "${outfile}")"

  echo -e "${cycle_label}\t${component}\t${filename}\t${url}\t${outfile}\t${remote_size}" >> "${EXPECTED_MANIFEST}"

  if needs_download "${outfile}" "${remote_size}"; then
    if [[ ! -s "${outfile}" ]]; then
      reason="missing_or_empty"
    else
      reason="size_mismatch"
    fi

    echo -e "${cycle_label}\t${component}\t${filename}\t${url}\t${outfile}\t${reason}\t${local_size}\t${remote_size}" >> "${DOWNLOADED_REPORT}"

    echo "Downloading: ${filename}"
    echo "  ${url}"

    mkdir -p "$(dirname "${outfile}")"
    tmpfile="${outfile}.part"
    rm -f "${tmpfile}"

    if curl_common -o "${tmpfile}" "${url}"; then
      if [[ ! -s "${tmpfile}" ]]; then
        rm -f "${tmpfile}"
        echo -e "${cycle_label}\t${component}\t${filename}\t${url}\t${outfile}\tdownloaded_file_empty" >> "${FAILED_REPORT}"
        return 1
      fi

      new_size="$(get_local_size "${tmpfile}")"

      if [[ "${remote_size}" != "NA" && "${new_size}" != "${remote_size}" ]]; then
        rm -f "${tmpfile}"
        echo -e "${cycle_label}\t${component}\t${filename}\t${url}\t${outfile}\tdownloaded_size_mismatch_${new_size}_vs_${remote_size}" >> "${FAILED_REPORT}"
        return 1
      fi

      mv "${tmpfile}" "${outfile}"
      echo "  Saved: ${outfile}"
    else
      rm -f "${tmpfile}"
      echo -e "${cycle_label}\t${component}\t${filename}\t${url}\t${outfile}\tcurl_failed" >> "${FAILED_REPORT}"
      return 1
    fi
  else
    echo "Verified existing file: ${filename}"
  fi
}

extract_xpt_links() {
  local page_url="$1"
  local html_file="$2"

  rm -f "${html_file}"

  curl_common -o "${html_file}" "${page_url}"

  python3 - "${page_url}" "${html_file}" <<'PY'
import sys
import re
from urllib.parse import urljoin

page_url = sys.argv[1]
html_file = sys.argv[2]

with open(html_file, "r", encoding="utf-8", errors="ignore") as f:
    html = f.read()

links = re.findall(r'href=["\']([^"\']+\.XPT)["\']', html, flags=re.IGNORECASE)

seen = set()
for link in links:
    full = urljoin(page_url, link)
    if full not in seen:
        seen.add(full)
        print(full)
PY
}

echo "============================================================"
echo "NHANES verify/redownload started"
echo "Output directory:"
echo "  ${OUTDIR}"
echo "============================================================"
echo

# ------------------------------------------------------------
# NHANES component XPT files
# ------------------------------------------------------------
for item in "${cycles[@]}"; do
  IFS="|" read -r cycle_label begin_year mort_cycle <<< "${item}"

  echo
  echo "------------------------------------------------------------"
  echo "Cycle: ${cycle_label}"
  echo "------------------------------------------------------------"

  cycle_dir="${OUTDIR}/${cycle_label}"
  mkdir -p "${cycle_dir}"

  for component in "${components[@]}"; do
    component_dir="${cycle_dir}/${component}"
    mkdir -p "${component_dir}"

    page_url="${BASE_NHANES_PAGE}?Component=${component}&CycleBeginYear=${begin_year}"

    echo
    echo "Component: ${component}"
    echo "Page: ${page_url}"

    html_file="${component_dir}/datapage.html"
    link_file="${component_dir}/download_links.txt"

    if extract_xpt_links "${page_url}" "${html_file}" > "${link_file}"; then
      n_links="$(wc -l < "${link_file}" | tr -d ' ')"
      echo "Expected XPT files from CDC page: ${n_links}"
    else
      echo "WARNING: could not fetch or parse ${page_url}"
      echo -e "${cycle_label}\t${component}\tNA\t${page_url}\t${link_file}\tcould_not_fetch_component_page" >> "${FAILED_REPORT}"
      continue
    fi

    if [[ "${n_links}" -eq 0 ]]; then
      echo "WARNING: no XPT links found for ${cycle_label} ${component}"
      echo -e "${cycle_label}\t${component}\tNA\t${page_url}\t${link_file}\tno_xpt_links_found" >> "${FAILED_REPORT}"
      continue
    fi

    while IFS= read -r url; do
      [[ -z "${url}" ]] && continue

      filename="$(basename "${url}")"
      outfile="${component_dir}/${filename}"

      download_one_file "${cycle_label}" "${component}" "${url}" "${outfile}" || true
    done < "${link_file}"
  done
done

# ------------------------------------------------------------
# Linked mortality files
# ------------------------------------------------------------
mort_dir="${OUTDIR}/linked_mortality_2019_public"
mkdir -p "${mort_dir}"

echo
echo "============================================================"
echo "Verifying/downloading linked mortality files"
echo "============================================================"

for item in "${cycles[@]}"; do
  IFS="|" read -r cycle_label begin_year mort_cycle <<< "${item}"

  mort_file="NHANES_${mort_cycle}_MORT_2019_PUBLIC.dat"
  mort_url="${BASE_MORT_URL}/${mort_file}"
  mort_outfile="${mort_dir}/${mort_file}"

  download_one_file "${cycle_label}" "Mortality" "${mort_url}" "${mort_outfile}" || true
done

# ------------------------------------------------------------
# Documentation and read-in programs
# ------------------------------------------------------------
echo
echo "============================================================"
echo "Verifying/downloading mortality documentation"
echo "============================================================"

download_one_file "all" "Mortality_docs" \
  "${BASE_MORT_URL}/R_ReadInProgramAllSurveys.R" \
  "${mort_dir}/R_ReadInProgramAllSurveys.R" || true

download_one_file "all" "Mortality_docs" \
  "${BASE_MORT_URL}/SAS_ReadInProgramAllSurveys.sas" \
  "${mort_dir}/SAS_ReadInProgramAllSurveys.sas" || true

download_one_file "all" "Mortality_docs" \
  "${BASE_MORT_URL}/Stata_ReadInProgramAllSurveys.do" \
  "${mort_dir}/Stata_ReadInProgramAllSurveys.do" || true

download_one_file "all" "Mortality_docs" \
  "https://www.cdc.gov/nchs/data/datalinkage/public-use-linked-mortality-files-data-dictionary.pdf" \
  "${mort_dir}/public-use-linked-mortality-files-data-dictionary.pdf" || true

# ------------------------------------------------------------
# Final completeness check against expected manifest
# ------------------------------------------------------------
echo
echo "============================================================"
echo "Final completeness check"
echo "============================================================"

tail -n +2 "${EXPECTED_MANIFEST}" | while IFS=$'\t' read -r cycle_label component filename url local_path remote_size; do
  local_size="$(get_local_size "${local_path}")"

  if [[ ! -s "${local_path}" ]]; then
    echo -e "${cycle_label}\t${component}\t${filename}\t${url}\t${local_path}\tmissing_or_empty\t${local_size}\t${remote_size}" >> "${FINAL_MISSING_REPORT}"
  elif [[ "${remote_size}" != "NA" && "${local_size}" != "${remote_size}" ]]; then
    echo -e "${cycle_label}\t${component}\t${filename}\t${url}\t${local_path}\tsize_mismatch\t${local_size}\t${remote_size}" >> "${FINAL_MISSING_REPORT}"
  fi
done

expected_xpt_files="$(
  tail -n +2 "${EXPECTED_MANIFEST}" \
    | awk -F'\t' 'tolower($3) ~ /\.xpt$/ {n++} END{print n+0}'
)"

local_xpt_files="$(
  find "${OUTDIR}" -type f -iname "*.xpt" | wc -l | tr -d ' '
)"

mortality_dat_files="$(
  find "${mort_dir}" -type f -name "NHANES_*_MORT_2019_PUBLIC.dat" | wc -l | tr -d ' '
)"

expected_mortality_dat_files=10

downloaded_or_redownloaded_rows="$(( $(wc -l < "${DOWNLOADED_REPORT}" | tr -d ' ') - 1 ))"
failed_download_rows="$(( $(wc -l < "${FAILED_REPORT}" | tr -d ' ') - 1 ))"
final_missing_or_incomplete_rows="$(( $(wc -l < "${FINAL_MISSING_REPORT}" | tr -d ' ') - 1 ))"

{
  echo -e "metric\tvalue"
  echo -e "expected_xpt_files_from_cdc_pages\t${expected_xpt_files}"
  echo -e "local_xpt_files\t${local_xpt_files}"
  echo -e "expected_mortality_dat_files\t${expected_mortality_dat_files}"
  echo -e "local_mortality_dat_files\t${mortality_dat_files}"
  echo -e "downloaded_or_redownloaded_rows\t${downloaded_or_redownloaded_rows}"
  echo -e "failed_download_rows\t${failed_download_rows}"
  echo -e "final_missing_or_incomplete_rows\t${final_missing_or_incomplete_rows}"
  echo -e "folder_size\t$(du -sh "${OUTDIR}" | awk '{print $1}')"
} > "${SUMMARY_REPORT}"

cat "${SUMMARY_REPORT}"

echo
echo "Reports:"
echo "  Expected manifest: ${EXPECTED_MANIFEST}"
echo "  Downloaded or re-downloaded: ${DOWNLOADED_REPORT}"
echo "  Failed downloads: ${FAILED_REPORT}"
echo "  Final missing/incomplete: ${FINAL_MISSING_REPORT}"
echo "  Verification summary: ${SUMMARY_REPORT}"

echo
echo "Top-level folders:"
find "${OUTDIR}" -mindepth 1 -maxdepth 1 -type d | sort

echo
echo "Done."