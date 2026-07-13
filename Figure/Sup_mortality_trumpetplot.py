#!/usr/bin/env python3

import argparse
import os
import re
from glob import glob
from typing import List, Optional, Dict

import numpy as np
import pandas as pd


# ============================================================
# 1. Default paths
# ============================================================

DEFAULT_BASE_DIR = "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock"

DEFAULT_FASTGWA_DIR = os.path.join(
    DEFAULT_BASE_DIR,
    "mortality_clock",
    "fastGWA",
    "output",
)

DEFAULT_FUMA_DIR = os.path.join(
    DEFAULT_BASE_DIR,
    "mortality_clock",
    "fuma",
)

DEFAULT_OUTDIR = os.path.join(
    DEFAULT_BASE_DIR,
    "Result",
    "TrumpetPlots_mortality_epoch",
)


# ============================================================
# 2. Expected mortality EPOCH clocks
# ============================================================

EXPECTED_CLOCKS = [
    "adipose_mri_mortality_clock",
    "brain_mri_mortality_clock",
    "Brain_proteomics_mortality_clock",
    "Digestive_metabolomics_mortality_clock",
    "Endocrine_metabolomics_mortality_clock",
    "Endocrine_proteomics_mortality_clock",
    "Eye_proteomics_mortality_clock",
    "heart_mri_mortality_clock",
    "Heart_proteomics_mortality_clock",
    "Hepatic_metabolomics_mortality_clock",
    "Hepatic_proteomics_mortality_clock",
    "Immune_metabolomics_mortality_clock",
    "Immune_proteomics_mortality_clock",
    "kidney_mri_mortality_clock",
    "liver_mri_mortality_clock",
    "pancreas_mri_mortality_clock",
    "Pulmonary_proteomics_mortality_clock",
    "Renal_proteomics_mortality_clock",
    "Reproductive_female_proteomics_mortality_clock",
    "Reproductive_male_proteomics_mortality_clock",
    "Skin_proteomics_mortality_clock",
    "spleen_mri_mortality_clock",
]


# ============================================================
# 3. Helpers
# ============================================================

def find_first_existing_col(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    exact = {c: c for c in df.columns}
    lower = {c.lower(): c for c in df.columns}

    for c in candidates:
        if c in exact:
            return exact[c]

    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]

    return None


def clean_token(x: str) -> str:
    x = str(x)
    x = x.replace("-", "_").replace(" ", "_")
    x = re.sub(r"_+", "_", x)
    return x


def infer_modality(clock_folder: str) -> str:
    x = clock_folder.replace("_mortality_clock", "").lower()

    if x.endswith("_mri"):
        return "MRI"
    if x.endswith("_proteomics"):
        return "Proteomics"
    if x.endswith("_metabolomics"):
        return "Metabolomics"

    return "Unknown"


def infer_organ_raw(clock_folder: str) -> str:
    x = clock_folder.replace("_mortality_clock", "")

    for suffix in ["_metabolomics", "_proteomics", "_mri"]:
        if x.lower().endswith(suffix):
            return re.sub(f"{suffix}$", "", x, flags=re.IGNORECASE)

    return x


def format_organ_label(organ_raw: str) -> str:
    x = clean_token(organ_raw)

    words = []
    for token in x.split("_"):
        t = token.lower()

        if t == "female":
            words.append("female")
        elif t == "male":
            words.append("male")
        else:
            words.append(t.capitalize())

    return " ".join(words)


def safe_read_tsv(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", low_memory=False)


def ensure_required_cols(df: pd.DataFrame, required: Dict[str, List[str]], path: str) -> Dict[str, str]:
    col_map = {}

    for standard_name, candidates in required.items():
        col = find_first_existing_col(df, candidates)

        if col is None:
            raise ValueError(
                f"Cannot find required column '{standard_name}' in {path}. "
                f"Tried candidates: {candidates}. "
                f"Available columns: {list(df.columns)}"
            )

        col_map[standard_name] = col

    return col_map


# ============================================================
# 4. Input discovery
# ============================================================

def build_clock_file_table(
    fastgwa_dir: str,
    fuma_dir: str,
    use_expected_only: bool = True,
) -> pd.DataFrame:

    rows = []

    if use_expected_only:
        clock_folders = EXPECTED_CLOCKS
    else:
        fastgwa_files = sorted(
            glob(os.path.join(fastgwa_dir, "*", "organ_pheno_normalized_residualized.fastGWA"))
        )
        clock_folders = sorted({os.path.basename(os.path.dirname(x)) for x in fastgwa_files})

    for clock_folder in clock_folders:
        fastgwa_file = os.path.join(
            fastgwa_dir,
            clock_folder,
            "organ_pheno_normalized_residualized.fastGWA",
        )

        iss_file = os.path.join(
            fuma_dir,
            clock_folder,
            "IndSigSNPs.txt",
        )

        organ_raw = infer_organ_raw(clock_folder)
        organ_label = format_organ_label(organ_raw)
        modality = infer_modality(clock_folder)

        rows.append(
            {
                "clock_folder": clock_folder,
                "organ_raw": organ_raw,
                "organ_label": organ_label,
                "modality": modality,
                "fastgwa_file": fastgwa_file,
                "iss_file": iss_file,
                "fastgwa_exists": os.path.exists(fastgwa_file),
                "iss_exists": os.path.exists(iss_file),
            }
        )

    return pd.DataFrame(rows)


# ============================================================
# 5. Core processing
# ============================================================

def read_independent_significant_snps(iss_file: str) -> pd.DataFrame:
    df_iss = safe_read_tsv(iss_file)

    rsid_col = find_first_existing_col(
        df_iss,
        ["rsID", "rsid", "SNP", "snp", "leadSNP", "LeadSNP", "rsIDuniq"],
    )

    if rsid_col is None:
        raise ValueError(
            f"Cannot find rsID/SNP column in {iss_file}. "
            f"Available columns: {list(df_iss.columns)}"
        )

    out = df_iss[[rsid_col]].copy()
    out.columns = ["rsID"]

    out["rsID"] = out["rsID"].astype(str)
    out = out.dropna(subset=["rsID"])
    out = out.loc[out["rsID"].str.lower() != "nan"].copy()
    out = out.drop_duplicates(subset=["rsID"])

    return out


def read_fastgwa_beta_afreq_n(fastgwa_file: str) -> pd.DataFrame:
    df_gwas = safe_read_tsv(fastgwa_file)

    col_map = ensure_required_cols(
        df_gwas,
        required={
            "SNP": ["SNP", "rsID", "rsid"],
            "freq": ["AF1", "A1FREQ", "A1_FREQ", "freq", "Freq", "MAF"],
            "A1_beta": ["BETA", "beta", "Effect", "effect"],
            "N": ["N", "n", "N_eff", "Neff", "N_total"],
        },
        path=fastgwa_file,
    )

    keep_cols = [
        col_map["SNP"],
        col_map["freq"],
        col_map["A1_beta"],
        col_map["N"],
    ]

    optional_cols = {}

    a1_col = find_first_existing_col(df_gwas, ["A1", "EA", "effect_allele"])
    a2_col = find_first_existing_col(df_gwas, ["A2", "NEA", "non_effect_allele"])
    se_col = find_first_existing_col(df_gwas, ["SE", "se", "BETA_SE"])
    p_col = find_first_existing_col(df_gwas, ["P", "p", "PVAL", "pval"])

    if a1_col is not None:
        keep_cols.append(a1_col)
        optional_cols[a1_col] = "A1"

    if a2_col is not None:
        keep_cols.append(a2_col)
        optional_cols[a2_col] = "A2"

    if se_col is not None:
        keep_cols.append(se_col)
        optional_cols[se_col] = "SE"

    if p_col is not None:
        keep_cols.append(p_col)
        optional_cols[p_col] = "P"

    out = df_gwas[keep_cols].copy()

    rename_map = {
        col_map["SNP"]: "rsID",
        col_map["freq"]: "freq",
        col_map["A1_beta"]: "A1_beta",
        col_map["N"]: "N",
    }

    rename_map.update(optional_cols)
    out = out.rename(columns=rename_map)

    out["rsID"] = out["rsID"].astype(str)
    out["freq"] = pd.to_numeric(out["freq"], errors="coerce")
    out["A1_beta"] = pd.to_numeric(out["A1_beta"], errors="coerce")
    out["N"] = pd.to_numeric(out["N"], errors="coerce")

    # Keep valid trumpet-plot inputs.
    out = out.dropna(subset=["rsID", "freq", "A1_beta", "N"]).copy()
    out = out.loc[(out["freq"] > 0) & (out["freq"] < 1)].copy()

    return out


def prepare_trumpet_input_for_clock(
    clock_folder: str,
    organ_raw: str,
    organ_label: str,
    modality: str,
    fastgwa_file: str,
    iss_file: str,
    outdir: str,
) -> pd.DataFrame:

    df_iss = read_independent_significant_snps(iss_file)
    df_gwas = read_fastgwa_beta_afreq_n(fastgwa_file)

    df_final = df_gwas.merge(
        df_iss,
        how="inner",
        on="rsID",
    )

    df_final["Gene"] = "NaN"
    df_final["Analysis"] = "GWAS"

    df_final["clock_folder"] = clock_folder
    df_final["organ_raw"] = organ_raw
    df_final["organ_label"] = organ_label
    df_final["modality"] = modality

    # TrumpetPlots required columns first.
    required_order = ["rsID", "freq", "A1_beta", "N", "Gene", "Analysis"]

    optional_order = [
        "A1",
        "A2",
        "SE",
        "P",
        "clock_folder",
        "organ_raw",
        "organ_label",
        "modality",
    ]

    col_order = required_order + [c for c in optional_order if c in df_final.columns]
    df_final = df_final[col_order].copy()

    out_file = os.path.join(
        outdir,
        f"TrumpetPlots_input_{clock_folder}.tsv",
    )

    df_final.to_csv(out_file, index=False, sep="\t", encoding="utf-8")

    return df_final


# ============================================================
# 6. Main
# ============================================================

def run(
    fastgwa_dir: str,
    fuma_dir: str,
    outdir: str,
    use_expected_only: bool,
    expected_n_clocks: int,
):

    os.makedirs(outdir, exist_ok=True)

    print("============================================================")
    print("Preparing TrumpetPlots input for 22 mortality EPOCH clocks")
    print("============================================================")
    print(f"fastGWA directory:     {fastgwa_dir}")
    print(f"FUMA directory:        {fuma_dir}")
    print(f"Output directory:      {outdir}")
    print(f"Use expected list:     {use_expected_only}")
    print(f"Expected clocks:       {expected_n_clocks}")

    df_clock = build_clock_file_table(
        fastgwa_dir=fastgwa_dir,
        fuma_dir=fuma_dir,
        use_expected_only=use_expected_only,
    )

    manifest_out = os.path.join(
        outdir,
        "TrumpetPlots_mortality_epoch_input_manifest.tsv",
    )

    df_clock.to_csv(manifest_out, index=False, sep="\t", encoding="utf-8")

    if df_clock.shape[0] != expected_n_clocks:
        print(
            f"WARNING: Found/discovered {df_clock.shape[0]} clock folders, "
            f"but expected {expected_n_clocks}."
        )

    missing = df_clock.loc[
        (~df_clock["fastgwa_exists"]) | (~df_clock["iss_exists"])
    ].copy()

    missing_out = os.path.join(
        outdir,
        "TrumpetPlots_mortality_epoch_missing_inputs.tsv",
    )

    missing.to_csv(missing_out, index=False, sep="\t", encoding="utf-8")

    if not missing.empty:
        print("WARNING: Some input files are missing.")
        print(missing[["clock_folder", "fastgwa_exists", "iss_exists"]].to_string(index=False))
        print(f"Missing-input table written to: {missing_out}")

    all_rows = []
    summary_rows = []

    for _, row in df_clock.iterrows():
        clock_folder = row["clock_folder"]

        print("------------------------------------------------------------")
        print(f"Clock: {clock_folder}")

        if not row["fastgwa_exists"]:
            print(f"  SKIP: missing fastGWA file: {row['fastgwa_file']}")
            summary_rows.append({
                "clock_folder": clock_folder,
                "status": "missing_fastgwa",
                "n_independent_significant_snps": np.nan,
                "n_rows_after_merge": 0,
                "fastgwa_file": row["fastgwa_file"],
                "iss_file": row["iss_file"],
                "output_file": "",
            })
            continue

        if not row["iss_exists"]:
            print(f"  SKIP: missing IndSigSNPs file: {row['iss_file']}")
            summary_rows.append({
                "clock_folder": clock_folder,
                "status": "missing_IndSigSNPs",
                "n_independent_significant_snps": np.nan,
                "n_rows_after_merge": 0,
                "fastgwa_file": row["fastgwa_file"],
                "iss_file": row["iss_file"],
                "output_file": "",
            })
            continue

        try:
            df_iss = read_independent_significant_snps(row["iss_file"])
            n_iss = df_iss.shape[0]

            df_final = prepare_trumpet_input_for_clock(
                clock_folder=clock_folder,
                organ_raw=row["organ_raw"],
                organ_label=row["organ_label"],
                modality=row["modality"],
                fastgwa_file=row["fastgwa_file"],
                iss_file=row["iss_file"],
                outdir=outdir,
            )

            out_file = os.path.join(
                outdir,
                f"TrumpetPlots_input_{clock_folder}.tsv",
            )

            print(f"  Independent significant SNPs: {n_iss}")
            print(f"  Rows after merge:             {df_final.shape[0]}")
            print(f"  Output:                       {out_file}")

            all_rows.append(df_final)

            summary_rows.append({
                "clock_folder": clock_folder,
                "status": "ok",
                "organ_raw": row["organ_raw"],
                "organ_label": row["organ_label"],
                "modality": row["modality"],
                "n_independent_significant_snps": n_iss,
                "n_rows_after_merge": df_final.shape[0],
                "fastgwa_file": row["fastgwa_file"],
                "iss_file": row["iss_file"],
                "output_file": out_file,
            })

        except Exception as e:
            print(f"  ERROR: {e}")

            summary_rows.append({
                "clock_folder": clock_folder,
                "status": "error",
                "organ_raw": row["organ_raw"],
                "organ_label": row["organ_label"],
                "modality": row["modality"],
                "n_independent_significant_snps": np.nan,
                "n_rows_after_merge": 0,
                "fastgwa_file": row["fastgwa_file"],
                "iss_file": row["iss_file"],
                "output_file": "",
                "error": str(e),
            })

    summary = pd.DataFrame(summary_rows)

    summary_out = os.path.join(
        outdir,
        "TrumpetPlots_mortality_epoch_input_summary.tsv",
    )

    summary.to_csv(summary_out, index=False, sep="\t", encoding="utf-8")

    if len(all_rows) > 0:
        df_combined = pd.concat(all_rows, ignore_index=True)

        combined_out = os.path.join(
            outdir,
            "TrumpetPlots_input_22_mortality_epoch_combined.tsv",
        )

        df_combined.to_csv(combined_out, index=False, sep="\t", encoding="utf-8")

    else:
        combined_out = ""

    print("============================================================")
    print("Finished.")
    print(f"Manifest:       {manifest_out}")
    print(f"Missing inputs: {missing_out}")
    print(f"Summary:        {summary_out}")

    if combined_out != "":
        print(f"Combined input: {combined_out}")

    print("Per-clock files:")
    print(f"  {outdir}/TrumpetPlots_input_<clock_folder>.tsv")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Prepare TrumpetPlots input files for 22 mortality EPOCH clocks "
            "by merging FUMA IndSigSNPs.txt with fastGWA AF/BETA/N."
        )
    )

    parser.add_argument(
        "--fastgwa_dir",
        default=DEFAULT_FASTGWA_DIR,
        type=str,
        help="Directory containing mortality_clock/fastGWA/output/<clock>/organ_pheno_normalized_residualized.fastGWA.",
    )

    parser.add_argument(
        "--fuma_dir",
        default=DEFAULT_FUMA_DIR,
        type=str,
        help="Directory containing mortality_clock/fuma/<clock>/IndSigSNPs.txt.",
    )

    parser.add_argument(
        "--outdir",
        default=DEFAULT_OUTDIR,
        type=str,
        help="Output directory for TrumpetPlots input files.",
    )

    parser.add_argument(
        "--discover_all",
        action="store_true",
        help=(
            "If set, discover all clock folders from fastGWA output. "
            "Default uses the hard-coded expected 22 mortality EPOCH clocks."
        ),
    )

    parser.add_argument(
        "--expected_n_clocks",
        default=22,
        type=int,
        help="Expected number of mortality EPOCH clocks.",
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    run(
        fastgwa_dir=args.fastgwa_dir,
        fuma_dir=args.fuma_dir,
        outdir=args.outdir,
        use_expected_only=not args.discover_all,
        expected_n_clocks=args.expected_n_clocks,
    )