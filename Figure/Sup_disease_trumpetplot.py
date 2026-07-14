#!/usr/bin/env python3

import argparse
import os
import re
from glob import glob
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


# ============================================================
# 1. Default paths
# ============================================================

DEFAULT_BASE_DIR = "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock"

DEFAULT_OUTDIR = os.path.join(
    DEFAULT_BASE_DIR,
    "Result",
    "TrumpetPlots_47_disease_epoch",
)

FASTGWA_BASENAME = "organ_pheno_normalized_residualized.fastGWA"
FUMA_BASENAME = "IndSigSNPs.txt"

DISEASE_ORDER = ["asthma", "copd", "dementia", "mi", "stroke"]
MODALITY_ORDER = ["mri", "proteomics", "metabolomics"]

DISEASE_LABEL = {
    "asthma": "Asthma",
    "copd": "COPD",
    "dementia": "Dementia",
    "mi": "MI",
    "stroke": "Stroke",
}

MODALITY_LABEL = {
    "mri": "MRI",
    "proteomics": "Proteomics",
    "metabolomics": "Metabolomics",
}


# ============================================================
# 2. Helper functions
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
    return x.strip("_")


def canonical_organ_key(organ_raw: str) -> str:
    x = clean_token(organ_raw).lower()

    if x in ["liver", "hepatic"]:
        return "hepatic"
    if x in ["kidney", "renal"]:
        return "renal"
    if x in ["reproductive_female", "female_reproductive"]:
        return "reproductive_female"
    if x in ["reproductive_male", "male_reproductive"]:
        return "reproductive_male"

    return x


def format_organ_label(organ_raw: str) -> str:
    key = canonical_organ_key(organ_raw)

    label_map = {
        "brain": "Brain",
        "heart": "Heart",
        "eye": "Eye",
        "hepatic": "Hepatic",
        "renal": "Renal",
        "pulmonary": "Pulmonary",
        "endocrine": "Endocrine",
        "immune": "Immune",
        "skin": "Skin",
        "digestive": "Digestive",
        "metabolic": "Metabolic",
        "adipose": "Adipose",
        "pancreas": "Pancreas",
        "spleen": "Spleen",
        "reproductive_female": "Reproductive female",
        "reproductive_male": "Reproductive male",
    }

    if key in label_map:
        return label_map[key]

    words = []
    for token in key.split("_"):
        if token in ["female", "male"]:
            words.append(token)
        else:
            words.append(token.capitalize())

    return " ".join(words)


def parse_disease_epoch_clock_folder(clock_folder: str) -> Dict[str, object]:
    """
    Parse folders such as:
      Brain_proteomics_dementia_clock
      heart_mri_copd_clock
      spleen_mri_asthma_clock
      Reproductive_female_proteomics_mi_clock
    """

    pattern = re.compile(
        r"^(?P<organ>.+)_(?P<modality>mri|proteomics|metabolomics)_(?P<disease>asthma|copd|dementia|mi|stroke)_clock$",
        flags=re.IGNORECASE,
    )

    m = pattern.match(clock_folder)

    if m is None:
        return {
            "parse_ok": False,
            "clock_folder": clock_folder,
            "organ_raw": "",
            "organ_key": "",
            "organ_label": "",
            "modality": "",
            "modality_label": "",
            "disease": "",
            "disease_label": "",
            "parse_error": "folder_name_does_not_match_expected_disease_epoch_pattern",
        }

    organ_raw = m.group("organ")
    modality = m.group("modality").lower()
    disease = m.group("disease").lower()

    organ_key = canonical_organ_key(organ_raw)

    return {
        "parse_ok": True,
        "clock_folder": clock_folder,
        "organ_raw": organ_raw,
        "organ_key": organ_key,
        "organ_label": format_organ_label(organ_raw),
        "modality": MODALITY_LABEL.get(modality, modality),
        "modality_key": modality,
        "disease": disease,
        "disease_label": DISEASE_LABEL.get(disease, disease.capitalize()),
        "parse_error": "",
    }


def safe_read_tsv(path: str, usecols: Optional[List[str]] = None) -> pd.DataFrame:
    return pd.read_csv(
        path,
        sep="\t",
        low_memory=False,
        usecols=usecols,
    )


def read_header(path: str) -> List[str]:
    return list(pd.read_csv(path, sep="\t", nrows=0).columns)


def ensure_required_cols_from_header(
    columns: List[str],
    required: Dict[str, List[str]],
    path: str,
) -> Dict[str, str]:

    exact = {c: c for c in columns}
    lower = {c.lower(): c for c in columns}
    col_map = {}

    for standard_name, candidates in required.items():
        found = None

        for c in candidates:
            if c in exact:
                found = exact[c]
                break

        if found is None:
            for c in candidates:
                if c.lower() in lower:
                    found = lower[c.lower()]
                    break

        if found is None:
            raise ValueError(
                f"Cannot find required column '{standard_name}' in {path}. "
                f"Tried candidates: {candidates}. "
                f"Available columns: {columns}"
            )

        col_map[standard_name] = found

    return col_map


# ============================================================
# 3. Input discovery
# ============================================================

def discover_disease_epoch_fastgwa_files(
    base_dir: str,
    allowed_diseases: List[str],
) -> pd.DataFrame:

    fastgwa_pattern = os.path.join(
        base_dir,
        "*_clock",
        "fastGWA",
        "output",
        FASTGWA_BASENAME,
    )

    fastgwa_files = sorted(glob(fastgwa_pattern))

    rows = []

    for fastgwa_file in fastgwa_files:
        rel_parts = os.path.relpath(fastgwa_file, base_dir).split(os.sep)

        if len(rel_parts) < 5:
            continue

        clock_folder = rel_parts[0]
        parsed = parse_disease_epoch_clock_folder(clock_folder)

        if not parsed["parse_ok"]:
            rows.append({
                **parsed,
                "fastgwa_file": fastgwa_file,
                "iss_file": "",
                "fastgwa_exists": os.path.exists(fastgwa_file),
                "iss_exists": False,
            })
            continue

        if parsed["disease"] not in allowed_diseases:
            continue

        iss_file = os.path.join(
            base_dir,
            clock_folder,
            "fuma",
            FUMA_BASENAME,
        )

        rows.append({
            **parsed,
            "fastgwa_file": fastgwa_file,
            "iss_file": iss_file,
            "fastgwa_exists": os.path.exists(fastgwa_file),
            "iss_exists": os.path.exists(iss_file),
        })

    df = pd.DataFrame(rows)

    if df.empty:
        return df

    disease_rank = {d: i for i, d in enumerate(DISEASE_ORDER)}
    modality_rank = {m: i for i, m in enumerate(["MRI", "Proteomics", "Metabolomics"])}

    df["disease_order"] = df["disease"].map(disease_rank).fillna(999).astype(int)
    df["modality_order"] = df["modality"].map(modality_rank).fillna(999).astype(int)

    df = df.sort_values(
        ["disease_order", "modality_order", "organ_label", "clock_folder"],
        kind="mergesort",
    ).reset_index(drop=True)

    return df


# ============================================================
# 4. Core processing
# ============================================================

def read_independent_significant_snps(iss_file: str) -> pd.DataFrame:
    try:
        df_iss = safe_read_tsv(iss_file)
    except pd.errors.EmptyDataError:
        return pd.DataFrame({"rsID": []})

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
    out = out.loc[out["rsID"] != ""].copy()
    out = out.drop_duplicates(subset=["rsID"])

    return out


def read_fastgwa_beta_afreq_n(fastgwa_file: str) -> pd.DataFrame:
    columns = read_header(fastgwa_file)

    col_map = ensure_required_cols_from_header(
        columns,
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

    for candidates, standard_name in [
        (["A1", "EA", "effect_allele"], "A1"),
        (["A2", "NEA", "non_effect_allele"], "A2"),
        (["SE", "se", "BETA_SE"], "SE"),
        (["P", "p", "PVAL", "pval"], "P"),
    ]:
        col = None
        lower = {c.lower(): c for c in columns}
        exact = {c: c for c in columns}

        for c in candidates:
            if c in exact:
                col = exact[c]
                break

        if col is None:
            for c in candidates:
                if c.lower() in lower:
                    col = lower[c.lower()]
                    break

        if col is not None:
            keep_cols.append(col)
            optional_cols[col] = standard_name

    keep_cols = list(dict.fromkeys(keep_cols))

    df_gwas = safe_read_tsv(fastgwa_file, usecols=keep_cols)

    rename_map = {
        col_map["SNP"]: "rsID",
        col_map["freq"]: "freq",
        col_map["A1_beta"]: "A1_beta",
        col_map["N"]: "N",
    }

    rename_map.update(optional_cols)

    out = df_gwas.rename(columns=rename_map).copy()

    out["rsID"] = out["rsID"].astype(str)
    out["freq"] = pd.to_numeric(out["freq"], errors="coerce")
    out["A1_beta"] = pd.to_numeric(out["A1_beta"], errors="coerce")
    out["N"] = pd.to_numeric(out["N"], errors="coerce")

    out = out.dropna(subset=["rsID", "freq", "A1_beta", "N"]).copy()
    out = out.loc[(out["freq"] > 0) & (out["freq"] < 1)].copy()
    out = out.loc[out["rsID"].str.lower() != "nan"].copy()
    out = out.drop_duplicates(subset=["rsID"])

    return out


def prepare_trumpet_input_for_clock(
    clock_row: pd.Series,
    outdir: str,
) -> Tuple[pd.DataFrame, int, str]:

    clock_folder = clock_row["clock_folder"]

    df_iss = read_independent_significant_snps(clock_row["iss_file"])
    n_iss = df_iss.shape[0]

    if df_iss.empty:
        df_final = pd.DataFrame(
            columns=[
                "rsID",
                "freq",
                "A1_beta",
                "N",
                "Gene",
                "Analysis",
                "A1",
                "A2",
                "SE",
                "P",
                "clock_folder",
                "organ_raw",
                "organ_key",
                "organ_label",
                "modality",
                "disease",
                "disease_label",
            ]
        )
    else:
        df_gwas = read_fastgwa_beta_afreq_n(clock_row["fastgwa_file"])

        df_final = df_gwas.merge(
            df_iss,
            how="inner",
            on="rsID",
        )

        df_final["Gene"] = "NaN"
        df_final["Analysis"] = "GWAS"

        df_final["clock_folder"] = clock_folder
        df_final["organ_raw"] = clock_row["organ_raw"]
        df_final["organ_key"] = clock_row["organ_key"]
        df_final["organ_label"] = clock_row["organ_label"]
        df_final["modality"] = clock_row["modality"]
        df_final["disease"] = clock_row["disease"]
        df_final["disease_label"] = clock_row["disease_label"]

        required_order = ["rsID", "freq", "A1_beta", "N", "Gene", "Analysis"]

        optional_order = [
            "A1",
            "A2",
            "SE",
            "P",
            "clock_folder",
            "organ_raw",
            "organ_key",
            "organ_label",
            "modality",
            "disease",
            "disease_label",
        ]

        col_order = required_order + [c for c in optional_order if c in df_final.columns]
        df_final = df_final[col_order].copy()

    out_file = os.path.join(
        outdir,
        f"TrumpetPlots_input_{clock_folder}.tsv",
    )

    df_final.to_csv(out_file, index=False, sep="\t", encoding="utf-8")

    return df_final, n_iss, out_file


# ============================================================
# 5. Main
# ============================================================

def run(
    base_dir: str,
    outdir: str,
    diseases: List[str],
    expected_n_clocks: int,
):

    os.makedirs(outdir, exist_ok=True)

    print("============================================================")
    print("Preparing TrumpetPlots input for 47 disease EPOCH clocks")
    print("============================================================")
    print(f"Base directory:        {base_dir}")
    print(f"Output directory:      {outdir}")
    print(f"Diseases:              {','.join(diseases)}")
    print(f"Expected clocks:       {expected_n_clocks}")

    df_clock = discover_disease_epoch_fastgwa_files(
        base_dir=base_dir,
        allowed_diseases=diseases,
    )

    manifest_out = os.path.join(
        outdir,
        "TrumpetPlots_47_disease_epoch_input_manifest.tsv",
    )

    df_clock.to_csv(manifest_out, index=False, sep="\t", encoding="utf-8")

    if df_clock.empty:
        raise RuntimeError(
            f"No disease EPOCH fastGWA files found under: {base_dir}"
        )

    n_parse_ok = int(df_clock["parse_ok"].sum())
    n_found = int(df_clock.shape[0])

    print(f"Discovered clock folders: {n_found}")
    print(f"Parse-ok clock folders:   {n_parse_ok}")

    if n_parse_ok != expected_n_clocks:
        print(
            f"WARNING: Found {n_parse_ok} parse-ok disease EPOCH clock folders, "
            f"but expected {expected_n_clocks}."
        )

    parse_failed = df_clock.loc[~df_clock["parse_ok"]].copy()

    parse_failed_out = os.path.join(
        outdir,
        "TrumpetPlots_47_disease_epoch_parse_failed.tsv",
    )

    parse_failed.to_csv(parse_failed_out, index=False, sep="\t", encoding="utf-8")

    df_clock_ok = df_clock.loc[df_clock["parse_ok"]].copy()

    missing = df_clock_ok.loc[
        (~df_clock_ok["fastgwa_exists"]) | (~df_clock_ok["iss_exists"])
    ].copy()

    missing_out = os.path.join(
        outdir,
        "TrumpetPlots_47_disease_epoch_missing_inputs.tsv",
    )

    missing.to_csv(missing_out, index=False, sep="\t", encoding="utf-8")

    if not missing.empty:
        print("WARNING: Some input files are missing.")
        print(
            missing[
                [
                    "clock_folder",
                    "disease_label",
                    "organ_label",
                    "modality",
                    "fastgwa_exists",
                    "iss_exists",
                ]
            ].to_string(index=False)
        )
        print(f"Missing-input table written to: {missing_out}")

    all_rows = []
    summary_rows = []

    for _, row in df_clock_ok.iterrows():
        clock_folder = row["clock_folder"]

        print("------------------------------------------------------------")
        print(f"Clock:    {clock_folder}")
        print(f"Disease:  {row['disease_label']}")
        print(f"Organ:    {row['organ_label']}")
        print(f"Modality: {row['modality']}")

        if not row["fastgwa_exists"]:
            print(f"  SKIP: missing fastGWA file: {row['fastgwa_file']}")

            summary_rows.append({
                "clock_folder": clock_folder,
                "status": "missing_fastgwa",
                "disease": row["disease"],
                "disease_label": row["disease_label"],
                "organ_raw": row["organ_raw"],
                "organ_key": row["organ_key"],
                "organ_label": row["organ_label"],
                "modality": row["modality"],
                "n_independent_significant_snps": np.nan,
                "n_rows_after_merge": 0,
                "fastgwa_file": row["fastgwa_file"],
                "iss_file": row["iss_file"],
                "output_file": "",
                "error": "",
            })

            continue

        if not row["iss_exists"]:
            print(f"  SKIP: missing IndSigSNPs file: {row['iss_file']}")

            summary_rows.append({
                "clock_folder": clock_folder,
                "status": "missing_IndSigSNPs",
                "disease": row["disease"],
                "disease_label": row["disease_label"],
                "organ_raw": row["organ_raw"],
                "organ_key": row["organ_key"],
                "organ_label": row["organ_label"],
                "modality": row["modality"],
                "n_independent_significant_snps": np.nan,
                "n_rows_after_merge": 0,
                "fastgwa_file": row["fastgwa_file"],
                "iss_file": row["iss_file"],
                "output_file": "",
                "error": "",
            })

            continue

        try:
            df_final, n_iss, out_file = prepare_trumpet_input_for_clock(
                clock_row=row,
                outdir=outdir,
            )

            n_merge = int(df_final.shape[0])

            print(f"  Independent significant SNPs: {n_iss}")
            print(f"  Rows after merge:             {n_merge}")
            print(f"  Output:                       {out_file}")

            if n_merge > 0:
                all_rows.append(df_final)

            status = "ok" if n_merge > 0 else "ok_empty_after_merge"

            summary_rows.append({
                "clock_folder": clock_folder,
                "status": status,
                "disease": row["disease"],
                "disease_label": row["disease_label"],
                "organ_raw": row["organ_raw"],
                "organ_key": row["organ_key"],
                "organ_label": row["organ_label"],
                "modality": row["modality"],
                "n_independent_significant_snps": n_iss,
                "n_rows_after_merge": n_merge,
                "fastgwa_file": row["fastgwa_file"],
                "iss_file": row["iss_file"],
                "output_file": out_file,
                "error": "",
            })

        except Exception as e:
            print(f"  ERROR: {e}")

            summary_rows.append({
                "clock_folder": clock_folder,
                "status": "error",
                "disease": row["disease"],
                "disease_label": row["disease_label"],
                "organ_raw": row["organ_raw"],
                "organ_key": row["organ_key"],
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
        "TrumpetPlots_47_disease_epoch_input_summary.tsv",
    )

    summary.to_csv(summary_out, index=False, sep="\t", encoding="utf-8")

    if len(all_rows) > 0:
        df_combined = pd.concat(all_rows, ignore_index=True)

        combined_out = os.path.join(
            outdir,
            "TrumpetPlots_input_47_disease_epoch_combined.tsv",
        )

        df_combined.to_csv(combined_out, index=False, sep="\t", encoding="utf-8")
    else:
        combined_out = ""

    by_disease_out = os.path.join(
        outdir,
        "TrumpetPlots_47_disease_epoch_summary_by_disease.tsv",
    )

    by_modality_out = os.path.join(
        outdir,
        "TrumpetPlots_47_disease_epoch_summary_by_modality.tsv",
    )

    if not summary.empty:
        summary.groupby("disease_label", dropna=False).agg(
            n_clocks=("clock_folder", "nunique"),
            n_ok=("status", lambda x: int((x == "ok").sum())),
            n_empty=("status", lambda x: int((x == "ok_empty_after_merge").sum())),
            n_error=("status", lambda x: int((x == "error").sum())),
            total_independent_significant_snps=("n_independent_significant_snps", "sum"),
            total_rows_after_merge=("n_rows_after_merge", "sum"),
        ).reset_index().to_csv(by_disease_out, index=False, sep="\t", encoding="utf-8")

        summary.groupby("modality", dropna=False).agg(
            n_clocks=("clock_folder", "nunique"),
            n_ok=("status", lambda x: int((x == "ok").sum())),
            n_empty=("status", lambda x: int((x == "ok_empty_after_merge").sum())),
            n_error=("status", lambda x: int((x == "error").sum())),
            total_independent_significant_snps=("n_independent_significant_snps", "sum"),
            total_rows_after_merge=("n_rows_after_merge", "sum"),
        ).reset_index().to_csv(by_modality_out, index=False, sep="\t", encoding="utf-8")

    print("============================================================")
    print("Finished.")
    print(f"Manifest:           {manifest_out}")
    print(f"Parse failures:     {parse_failed_out}")
    print(f"Missing inputs:     {missing_out}")
    print(f"Summary:            {summary_out}")
    print(f"By disease summary: {by_disease_out}")
    print(f"By modality summary:{by_modality_out}")

    if combined_out != "":
        print(f"Combined input:     {combined_out}")

    print("Per-clock files:")
    print(f"  {outdir}/TrumpetPlots_input_<clock_folder>.tsv")
    print("============================================================")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Prepare TrumpetPlots input files for 47 disease EPOCH clocks "
            "by merging FUMA IndSigSNPs.txt with fastGWA AF/BETA/N."
        )
    )

    parser.add_argument(
        "--base_dir",
        default=DEFAULT_BASE_DIR,
        type=str,
        help=(
            "WholeBodyClock base directory containing <clock_folder>/fastGWA/output/"
            "organ_pheno_normalized_residualized.fastGWA and <clock_folder>/fuma/IndSigSNPs.txt."
        ),
    )

    parser.add_argument(
        "--outdir",
        default=DEFAULT_OUTDIR,
        type=str,
        help="Output directory for TrumpetPlots input files.",
    )

    parser.add_argument(
        "--diseases",
        default="asthma,copd,dementia,mi,stroke",
        type=str,
        help="Comma-separated diseases to include.",
    )

    parser.add_argument(
        "--expected_n_clocks",
        default=47,
        type=int,
        help="Expected number of disease EPOCH clocks.",
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    diseases = [
        x.strip().lower()
        for x in args.diseases.split(",")
        if x.strip() != ""
    ]

    run(
        base_dir=args.base_dir,
        outdir=args.outdir,
        diseases=diseases,
        expected_n_clocks=args.expected_n_clocks,
    )