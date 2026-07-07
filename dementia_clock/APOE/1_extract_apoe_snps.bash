#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=apoe_status
#SBATCH --time=0-04:59:00
#SBATCH --mem-per-cpu=24G
#SBATCH --output=/cbica/home/wenju/output/apoe_status_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/apoe_status_%A_%a.err

module load plink/2.20210701

set -euo pipefail

# ============================================================
# Define APOE genotype/status in UK Biobank using PLINK2
# SNPs:
#   rs429358
#   rs7412
#
# APOE allele definitions:
#   epsilon2: rs429358 T + rs7412 T
#   epsilon3: rs429358 T + rs7412 C
#   epsilon4: rs429358 C + rs7412 C
#
# Carrier group definitions:
#   APOE3_ref        = e3/e3
#   APOE4_carrier    = e3/e4 or e4/e4
#   APOE2_carrier    = e2/e2 or e2/e3
#   APOE2_APOE4_mixed = e2/e4
# ============================================================

THREADS=4

BFILE="/cbica/projects/MULTI/processed/UKBB/UKBB_Pe/genetics/S3_apply_all/chr_all_AllUKBBPeople"
APOE_SNPS="/cbica/home/wenju/Project/Surreal_GAN_genetic_paper/ADNI/GWAS/APOE/apoe_snps.txt"

OUTDIR="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/apoe_status_ukbb"
mkdir -p "${OUTDIR}"

OUTPREFIX="${OUTDIR}/ukbb_apoe_rs429358_rs7412"

# ------------------------------------------------------------
# 1. Tell PLINK2 which allele to count.
#
# We want:
#   rs429358_C dosage = number of epsilon4-defining C alleles
#   rs7412_T dosage   = number of epsilon2-defining T alleles
# ------------------------------------------------------------

COUNT_ALLELES="${OUTDIR}/apoe_count_alleles.txt"

cat > "${COUNT_ALLELES}" <<EOF
rs429358 C
rs7412 T
EOF

# ------------------------------------------------------------
# 2. Extract APOE SNPs and export additive allele counts.
#
# Output:
#   ${OUTPREFIX}.raw
#
# The .raw file should contain one row per individual and columns
# corresponding to allele counts for rs429358 and rs7412.
# ------------------------------------------------------------

plink2 \
  --threads "${THREADS}" \
  --bfile "${BFILE}" \
  --extract "${APOE_SNPS}" \
  --export A \
  --export-allele "${COUNT_ALLELES}" \
  --out "${OUTPREFIX}"

# ------------------------------------------------------------
# 3. Map SNP allele counts to APOE genotype/status.
#
# Expected dosage combinations:
#
# rs429358_C  rs7412_T  APOE genotype
#     0          2      e2/e2
#     0          1      e2/e3
#     0          0      e3/e3
#     1          1      e2/e4
#     1          0      e3/e4
#     2          0      e4/e4
#
# Other combinations are biologically rare/ambiguous or may indicate
# allele-coding problems/missingness and are marked NA/ambiguous.
# ------------------------------------------------------------

python3 <<EOF
import pandas as pd
import re
from pathlib import Path

raw_file = Path("${OUTPREFIX}.raw")
out_file = Path("${OUTDIR}/ukbb_apoe_status.tsv")
summary_file = Path("${OUTDIR}/ukbb_apoe_status_counts.tsv")

df = pd.read_csv(raw_file, sep=r"\\s+", engine="python")

# PLINK2 .raw usually has ID columns plus genotype columns.
# Depending on PLINK2 version/settings, columns may look like:
#   rs429358_C, rs7412_T
# or may contain additional suffixes. We detect them robustly.

def find_snp_col(df, rsid):
    matches = [c for c in df.columns if c == rsid or c.startswith(rsid + "_")]
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one column for {rsid}, found {len(matches)}: {matches}\\n"
            f"Available columns: {list(df.columns)}"
        )
    return matches[0]

col_429358 = find_snp_col(df, "rs429358")
col_7412 = find_snp_col(df, "rs7412")

# Rename to explicit allele-count names.
df = df.rename(columns={
    col_429358: "rs429358_C_count",
    col_7412: "rs7412_T_count"
})

# Convert to numeric. Missing genotypes remain NA.
df["rs429358_C_count"] = pd.to_numeric(df["rs429358_C_count"], errors="coerce")
df["rs7412_T_count"] = pd.to_numeric(df["rs7412_T_count"], errors="coerce")

def assign_apoe(row):
    c429 = row["rs429358_C_count"]
    t7412 = row["rs7412_T_count"]

    if pd.isna(c429) or pd.isna(t7412):
        return "missing"

    c429 = int(c429)
    t7412 = int(t7412)

    if c429 == 0 and t7412 == 2:
        return "e2/e2"
    elif c429 == 0 and t7412 == 1:
        return "e2/e3"
    elif c429 == 0 and t7412 == 0:
        return "e3/e3"
    elif c429 == 1 and t7412 == 1:
        return "e2/e4"
    elif c429 == 1 and t7412 == 0:
        return "e3/e4"
    elif c429 == 2 and t7412 == 0:
        return "e4/e4"
    else:
        return "ambiguous"

def assign_group(gt):
    if gt == "e3/e3":
        return "APOE3_ref"
    elif gt in ["e3/e4", "e4/e4"]:
        return "APOE4_carrier"
    elif gt in ["e2/e2", "e2/e3"]:
        return "APOE2_carrier"
    elif gt == "e2/e4":
        return "APOE2_APOE4_mixed"
    elif gt == "missing":
        return "missing"
    else:
        return "ambiguous"

def e2_count(gt):
    return {"e2/e2": 2, "e2/e3": 1, "e2/e4": 1,
            "e3/e3": 0, "e3/e4": 0, "e4/e4": 0}.get(gt, pd.NA)

def e4_count(gt):
    return {"e4/e4": 2, "e3/e4": 1, "e2/e4": 1,
            "e3/e3": 0, "e2/e3": 0, "e2/e2": 0}.get(gt, pd.NA)

df["APOE_genotype"] = df.apply(assign_apoe, axis=1)
df["APOE_group"] = df["APOE_genotype"].apply(assign_group)
df["APOE_e2_count"] = df["APOE_genotype"].apply(e2_count)
df["APOE_e4_count"] = df["APOE_genotype"].apply(e4_count)

# Keep standard ID columns if present.
id_cols = [c for c in ["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"] if c in df.columns]

keep_cols = (
    id_cols
    + [
        "rs429358_C_count",
        "rs7412_T_count",
        "APOE_genotype",
        "APOE_group",
        "APOE_e2_count",
        "APOE_e4_count",
    ]
)

df_out = df[keep_cols]
df_out.to_csv(out_file, sep="\\t", index=False)

summary = (
    df_out
    .groupby(["APOE_genotype", "APOE_group"], dropna=False)
    .size()
    .reset_index(name="N")
    .sort_values(["APOE_group", "APOE_genotype"])
)

summary.to_csv(summary_file, sep="\\t", index=False)

print(f"Wrote individual APOE status to: {out_file}")
print(f"Wrote APOE status counts to:     {summary_file}")
print()
print(summary.to_string(index=False))
EOF

echo "Done."