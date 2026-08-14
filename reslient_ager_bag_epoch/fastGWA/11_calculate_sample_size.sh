### need to run this on CUBIC

cd ~/Reproducibile_paper/WholeBodyClock

OUT="Result/mortality_epoch_fastGWA_pheno_sample_sizes.tsv"
mkdir -p Result

printf "Clock\tN\n" > "${OUT}"

find mortality_clock/fastGWA/data \
  -mindepth 2 -maxdepth 2 \
  -type f \
  -name "EPOCH_pheno_normalized_residualized_with_related_indi.phen" \
  | sort \
  | while read -r f; do
      clock="$(basename "$(dirname "${f}")")"
      n="$(awk 'NF>=3 && $3!="NA" && $3!="NaN" && $3!="nan" {n++} END{print n+0}' "${f}")"
      printf "%s\t%s\n" "${clock}" "${n}" >> "${OUT}"
    done

column -t "${OUT}"

MIN_N="$(awk 'NR>1 {if (NR==2 || $2<min) min=$2} END{print min}' "${OUT}")"
MAX_N="$(awk 'NR>1 {if (NR==2 || $2>max) max=$2} END{print max}' "${OUT}")"
N_CLOCKS="$(awk 'NR>1 {n++} END{print n+0}' "${OUT}")"

echo
echo "Number of mortality EPOCH clocks with phenotype files: ${N_CLOCKS}"
echo "Sample-size range: ${MIN_N} <= N <= ${MAX_N}"
echo
echo "Sentence:"
echo "We conducted GWAS (Method 3a) for the 22 mortality EPOCH clocks (${MIN_N} <= N <= ${MAX_N} participants with European ancestries) and identified 365 genomic locus-EPOCH pairs at the Bonferroni-corrected genome-wide significance threshold (P < 5 x 10^-8 / 22 = 2.27 x 10^-9)."