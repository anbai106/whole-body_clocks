#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Collect GCTA/HEreg SNP heritability estimates for 23 BAGs.

Sources:
  1) MRIBAGs, 7 organs:
     /Users/hao/cubic-home/Reproducibile_paper/AbdoImaging/Result/GCTA_h2_results_genotype_array.tsv

  2) ProtBAGs, 11 organs:
     /Users/hao/cubic-home/Reproducibile_paper/UKBB_Proteomics/h2/*/*.hsq

  3) MetBAGs, 5 systems:
     /Users/hao/cubic-home/Reproducibile_paper/UKBB_metabolomics/h2/*/*.HEreg

Outputs:
  /Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/Result/BAG_23_GCTA_h2_summary.tsv
  /Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/Result/BAG_23_GCTA_h2_summary_wide.tsv
  /Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/Result/BAG_23_GCTA_h2_parse_report.tsv

Notes:
  - For .hsq files, this script extracts the row Source == "V(G)/Vp".
  - For .HEreg files, this script uses HE-CP by default.
  - HE-SD is also parsed and saved in the parse report, but HE-CP is used as the main h2 estimate.
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd


# ============================================================
# 1. Default paths
# ============================================================

DEFAULT_MRI_TSV = Path(
    "/Users/hao/cubic-home/Reproducibile_paper/AbdoImaging/Result/"
    "GCTA_h2_results_genotype_array.tsv"
)

DEFAULT_PROT_DIR = Path(
    "/Users/hao/cubic-home/Reproducibile_paper/UKBB_Proteomics/h2"
)

DEFAULT_MET_DIR = Path(
    "/Users/hao/cubic-home/Reproducibile_paper/UKBB_metabolomics/h2"
)

DEFAULT_OUT_DIR = Path(
    "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/Result"
)


# ============================================================
# 2. Expected BAGs
# ============================================================

EXPECTED_MRI = [
    "brain",
    "adipose",
    "heart",
    "kidney",
    "liver",
    "pancreas",
    "spleen",
]

EXPECTED_PROTEOMICS = [
    "Brain",
    "Endocrine",
    "Eye",
    "Heart",
    "Hepatic",
    "Immune",
    "Pulmonary",
    "Renal",
    "Reproductive_female",
    "Reproductive_male",
    "Skin",
]

EXPECTED_METABOLOMICS = [
    "Digestive",
    "Endocrine",
    "Hepatic",
    "Immune",
    "Metabolic",
]


# ============================================================
# 3. Helpers
# ============================================================

def safe_float(x) -> float:
    """Convert a value to float; return NaN if conversion fails."""
    if x is None:
        return math.nan
    try:
        s = str(x).strip()
        if s == "" or s.lower() in {"na", "nan", "none", "null"}:
            return math.nan
        return float(s)
    except Exception:
        return math.nan


def clean_space(x: str) -> str:
    return re.sub(r"\s+", " ", str(x).strip())


def normalize_organ_for_label(x: str) -> str:
    """Normalize organ/system names for consistent output."""
    x = str(x).strip()
    x_lower = x.lower()

    mapping = {
        "brain": "Brain",
        "adipose": "Adipose",
        "heart": "Heart",
        "kidney": "Kidney",
        "renal": "Renal",
        "liver": "Liver",
        "hepatic": "Hepatic",
        "pancreas": "Pancreas",
        "spleen": "Spleen",
        "endocrine": "Endocrine",
        "eye": "Eye",
        "immune": "Immune",
        "pulmonary": "Pulmonary",
        "skin": "Skin",
        "digestive": "Digestive",
        "metabolic": "Metabolic",
        "reproductive_female": "Reproductive_female",
        "reproductive_male": "Reproductive_male",
    }

    return mapping.get(x_lower, x)


def make_bag_id(organ: str, modality: str) -> str:
    organ_clean = str(organ).strip().replace(" ", "_")
    modality_clean = str(modality).strip().lower()
    return f"{organ_clean}_{modality_clean}_BAG"


def parse_pvalue_zero_ok(x) -> float:
    """
    Parse p-value. The text '0' from GCTA/HEreg is a valid very-small p-value.
    Keep it as 0.0 rather than NA.
    """
    return safe_float(x)


def read_text_file(path: Path) -> List[str]:
    with path.open("r", encoding="utf-8", errors="replace") as f:
        return f.read().splitlines()


# ============================================================
# 4. Parse MRIBAG TSV
# ============================================================

def collect_mri_h2(mri_tsv: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    report = []

    if not mri_tsv.exists():
        report.append({
            "source_type": "MRI",
            "source_file": str(mri_tsv),
            "parse_status": "missing_file",
            "message": "MRI GCTA h2 TSV not found",
        })
        return pd.DataFrame(rows), pd.DataFrame(report)

    df = pd.read_csv(mri_tsv, sep="\t")

    required = {"mribag", "Heritability", "SE", "Pvalue"}
    missing = required - set(df.columns)

    if missing:
        report.append({
            "source_type": "MRI",
            "source_file": str(mri_tsv),
            "parse_status": "missing_required_columns",
            "message": f"Missing columns: {', '.join(sorted(missing))}",
        })
        return pd.DataFrame(rows), pd.DataFrame(report)

    for _, r in df.iterrows():
        organ_raw = str(r["mribag"]).strip()
        organ = normalize_organ_for_label(organ_raw)

        h2 = safe_float(r["Heritability"])
        se = safe_float(r["SE"])
        pval = parse_pvalue_zero_ok(r["Pvalue"])

        status = "ok" if math.isfinite(h2) else "missing_h2"

        rows.append({
            "bag_id": make_bag_id(organ, "MRI"),
            "organ": organ,
            "modality": "MRI",
            "method": "GCTA_REML_genotype_array",
            "h2": h2,
            "se": se,
            "pvalue": pval,
            "n": math.nan,
            "source_type": "MRI_GCTA_summary_TSV",
            "source_file": str(mri_tsv),
            "source_row_id": organ_raw,
            "parse_status": status,
            "he_block": "",
        })

        report.append({
            "source_type": "MRI",
            "source_file": str(mri_tsv),
            "source_row_id": organ_raw,
            "parse_status": status,
            "message": "",
        })

    return pd.DataFrame(rows), pd.DataFrame(report)


# ============================================================
# 5. Parse .hsq files
# ============================================================

def parse_hsq_file(hsq_file: Path, organ: str, modality: str) -> Tuple[dict, List[dict]]:
    report_rows = []

    try:
        df = pd.read_csv(hsq_file, sep="\t")
    except Exception as e:
        row = {
            "bag_id": make_bag_id(organ, modality),
            "organ": organ,
            "modality": modality,
            "method": "GCTA_REML_hsq",
            "h2": math.nan,
            "se": math.nan,
            "pvalue": math.nan,
            "n": math.nan,
            "source_type": "GCTA_hsq",
            "source_file": str(hsq_file),
            "source_row_id": "V(G)/Vp",
            "parse_status": "read_error",
            "he_block": "",
        }
        report_rows.append({
            "source_type": modality,
            "source_file": str(hsq_file),
            "source_row_id": organ,
            "parse_status": "read_error",
            "message": str(e),
        })
        return row, report_rows

    required = {"Source", "Variance", "SE"}
    if not required.issubset(set(df.columns)):
        row = {
            "bag_id": make_bag_id(organ, modality),
            "organ": organ,
            "modality": modality,
            "method": "GCTA_REML_hsq",
            "h2": math.nan,
            "se": math.nan,
            "pvalue": math.nan,
            "n": math.nan,
            "source_type": "GCTA_hsq",
            "source_file": str(hsq_file),
            "source_row_id": "V(G)/Vp",
            "parse_status": "missing_required_columns",
            "he_block": "",
        }
        report_rows.append({
            "source_type": modality,
            "source_file": str(hsq_file),
            "source_row_id": organ,
            "parse_status": "missing_required_columns",
            "message": f"Columns found: {', '.join(df.columns)}",
        })
        return row, report_rows

    df["Source"] = df["Source"].astype(str).str.strip()

    h2_rows = df.loc[df["Source"] == "V(G)/Vp"]
    p_rows = df.loc[df["Source"] == "Pval"]
    n_rows = df.loc[df["Source"] == "n"]

    if h2_rows.empty:
        h2 = se = math.nan
        status = "missing_VG_over_Vp_row"
    else:
        h2 = safe_float(h2_rows.iloc[0]["Variance"])
        se = safe_float(h2_rows.iloc[0]["SE"])
        status = "ok" if math.isfinite(h2) else "missing_h2_value"

    pval = parse_pvalue_zero_ok(p_rows.iloc[0]["Variance"]) if not p_rows.empty else math.nan
    n = safe_float(n_rows.iloc[0]["Variance"]) if not n_rows.empty else math.nan

    row = {
        "bag_id": make_bag_id(organ, modality),
        "organ": organ,
        "modality": modality,
        "method": "GCTA_REML_hsq",
        "h2": h2,
        "se": se,
        "pvalue": pval,
        "n": n,
        "source_type": "GCTA_hsq",
        "source_file": str(hsq_file),
        "source_row_id": "V(G)/Vp",
        "parse_status": status,
        "he_block": "",
    }

    report_rows.append({
        "source_type": modality,
        "source_file": str(hsq_file),
        "source_row_id": organ,
        "parse_status": status,
        "message": "",
    })

    return row, report_rows


def collect_proteomics_h2(prot_dir: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    report = []

    if not prot_dir.exists():
        report.append({
            "source_type": "Proteomics",
            "source_file": str(prot_dir),
            "source_row_id": "",
            "parse_status": "missing_directory",
            "message": "Proteomics h2 directory not found",
        })
        return pd.DataFrame(rows), pd.DataFrame(report)

    for organ in EXPECTED_PROTEOMICS:
        hsq_file = prot_dir / organ / f"{organ}.hsq"

        if not hsq_file.exists():
            report.append({
                "source_type": "Proteomics",
                "source_file": str(hsq_file),
                "source_row_id": organ,
                "parse_status": "missing_file",
                "message": "Expected .hsq file not found",
            })
            rows.append({
                "bag_id": make_bag_id(organ, "Proteomics"),
                "organ": organ,
                "modality": "Proteomics",
                "method": "GCTA_REML_hsq",
                "h2": math.nan,
                "se": math.nan,
                "pvalue": math.nan,
                "n": math.nan,
                "source_type": "GCTA_hsq",
                "source_file": str(hsq_file),
                "source_row_id": "V(G)/Vp",
                "parse_status": "missing_file",
                "he_block": "",
            })
            continue

        row, rep = parse_hsq_file(
            hsq_file=hsq_file,
            organ=organ,
            modality="Proteomics",
        )

        rows.append(row)
        report.extend(rep)

    return pd.DataFrame(rows), pd.DataFrame(report)


# ============================================================
# 6. Parse .HEreg files
# ============================================================

def parse_hereg_blocks(hereg_file: Path) -> Dict[str, pd.DataFrame]:
    """
    Parse HEreg files with blocks such as:

    HE-CP
    Coefficient     Estimate        SE_OLS          SE_Jackknife    P_OLS           P_Jackknife
    Intercept       ...
    V(G)/Vp         ...

    HE-SD
    ...
    """
    lines = read_text_file(hereg_file)
    blocks: Dict[str, pd.DataFrame] = {}

    block_positions = []
    for i, line in enumerate(lines):
        s = line.strip()
        if s in {"HE-CP", "HE-SD"}:
            block_positions.append((s, i))

    for b_idx, (block_name, start_i) in enumerate(block_positions):
        header_i = start_i + 1
        data_start = start_i + 2

        if header_i >= len(lines):
            continue

        if b_idx + 1 < len(block_positions):
            data_end = block_positions[b_idx + 1][1]
        else:
            data_end = len(lines)

        header = clean_space(lines[header_i]).split(" ")
        data_lines = [
            clean_space(x)
            for x in lines[data_start:data_end]
            if clean_space(x) != ""
        ]

        parsed_rows = []
        for dl in data_lines:
            parts = dl.split(" ")
            if len(parts) < len(header):
                continue

            # Join any extra columns onto the first field defensively,
            # although current files appear clean.
            if len(parts) > len(header):
                extra = len(parts) - len(header)
                first = "_".join(parts[: extra + 1])
                parts = [first] + parts[extra + 1:]

            parsed_rows.append(dict(zip(header, parts)))

        if parsed_rows:
            blocks[block_name] = pd.DataFrame(parsed_rows)

    return blocks


def parse_hereg_file(
    hereg_file: Path,
    organ: str,
    modality: str,
    preferred_block: str = "HE-CP",
) -> Tuple[dict, List[dict], pd.DataFrame]:
    """
    Return:
      - main row using preferred HE block
      - parse report rows
      - all block estimates table
    """
    report_rows = []
    all_block_rows = []

    try:
        blocks = parse_hereg_blocks(hereg_file)
    except Exception as e:
        main = {
            "bag_id": make_bag_id(organ, modality),
            "organ": organ,
            "modality": modality,
            "method": "HEreg",
            "h2": math.nan,
            "se": math.nan,
            "pvalue": math.nan,
            "n": math.nan,
            "source_type": "HEreg",
            "source_file": str(hereg_file),
            "source_row_id": "V(G)/Vp",
            "parse_status": "read_error",
            "he_block": preferred_block,
        }
        report_rows.append({
            "source_type": modality,
            "source_file": str(hereg_file),
            "source_row_id": organ,
            "parse_status": "read_error",
            "message": str(e),
        })
        return main, report_rows, pd.DataFrame(all_block_rows)

    if not blocks:
        main = {
            "bag_id": make_bag_id(organ, modality),
            "organ": organ,
            "modality": modality,
            "method": "HEreg",
            "h2": math.nan,
            "se": math.nan,
            "pvalue": math.nan,
            "n": math.nan,
            "source_type": "HEreg",
            "source_file": str(hereg_file),
            "source_row_id": "V(G)/Vp",
            "parse_status": "missing_HE_blocks",
            "he_block": preferred_block,
        }
        report_rows.append({
            "source_type": modality,
            "source_file": str(hereg_file),
            "source_row_id": organ,
            "parse_status": "missing_HE_blocks",
            "message": "No HE-CP or HE-SD blocks found",
        })
        return main, report_rows, pd.DataFrame(all_block_rows)

    for block_name, df in blocks.items():
        df_cols = {c.lower(): c for c in df.columns}
        coef_col = df_cols.get("coefficient")
        est_col = df_cols.get("estimate")
        se_jack_col = df_cols.get("se_jackknife")
        se_ols_col = df_cols.get("se_ols")
        p_jack_col = df_cols.get("p_jackknife")
        p_ols_col = df_cols.get("p_ols")

        if coef_col is None or est_col is None:
            continue

        df[coef_col] = df[coef_col].astype(str).str.strip()
        h2_rows = df.loc[df[coef_col] == "V(G)/Vp"]

        if h2_rows.empty:
            all_block_rows.append({
                "bag_id": make_bag_id(organ, modality),
                "organ": organ,
                "modality": modality,
                "he_block": block_name,
                "h2": math.nan,
                "se_jackknife": math.nan,
                "se_ols": math.nan,
                "p_jackknife": math.nan,
                "p_ols": math.nan,
                "source_file": str(hereg_file),
                "parse_status": "missing_VG_over_Vp_row",
            })
            continue

        r = h2_rows.iloc[0]
        all_block_rows.append({
            "bag_id": make_bag_id(organ, modality),
            "organ": organ,
            "modality": modality,
            "he_block": block_name,
            "h2": safe_float(r[est_col]),
            "se_jackknife": safe_float(r[se_jack_col]) if se_jack_col else math.nan,
            "se_ols": safe_float(r[se_ols_col]) if se_ols_col else math.nan,
            "p_jackknife": parse_pvalue_zero_ok(r[p_jack_col]) if p_jack_col else math.nan,
            "p_ols": parse_pvalue_zero_ok(r[p_ols_col]) if p_ols_col else math.nan,
            "source_file": str(hereg_file),
            "parse_status": "ok",
        })

    all_blocks_df = pd.DataFrame(all_block_rows)

    if all_blocks_df.empty:
        main_status = "missing_VG_over_Vp_row"
        h2 = se = pval = math.nan
        used_block = preferred_block
    else:
        available_blocks = set(all_blocks_df["he_block"].dropna().astype(str))
        if preferred_block in available_blocks:
            used_block = preferred_block
        elif "HE-CP" in available_blocks:
            used_block = "HE-CP"
        elif "HE-SD" in available_blocks:
            used_block = "HE-SD"
        else:
            used_block = str(all_blocks_df["he_block"].iloc[0])

        main_block = all_blocks_df.loc[
            (all_blocks_df["he_block"] == used_block)
            & (all_blocks_df["parse_status"] == "ok")
        ]

        if main_block.empty:
            main_status = "no_ok_preferred_block"
            h2 = se = pval = math.nan
        else:
            mb = main_block.iloc[0]
            h2 = safe_float(mb["h2"])
            se = safe_float(mb["se_jackknife"])
            if not math.isfinite(se):
                se = safe_float(mb["se_ols"])

            pval = parse_pvalue_zero_ok(mb["p_jackknife"])
            if not math.isfinite(pval):
                pval = parse_pvalue_zero_ok(mb["p_ols"])

            main_status = "ok" if math.isfinite(h2) else "missing_h2_value"

    main = {
        "bag_id": make_bag_id(organ, modality),
        "organ": organ,
        "modality": modality,
        "method": "HEreg",
        "h2": h2,
        "se": se,
        "pvalue": pval,
        "n": math.nan,
        "source_type": "HEreg",
        "source_file": str(hereg_file),
        "source_row_id": "V(G)/Vp",
        "parse_status": main_status,
        "he_block": used_block,
    }

    report_rows.append({
        "source_type": modality,
        "source_file": str(hereg_file),
        "source_row_id": organ,
        "parse_status": main_status,
        "message": f"Used block: {used_block}",
    })

    return main, report_rows, all_blocks_df


def collect_metabolomics_h2(
    met_dir: Path,
    preferred_block: str = "HE-CP",
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rows = []
    report = []
    all_blocks = []

    if not met_dir.exists():
        report.append({
            "source_type": "Metabolomics",
            "source_file": str(met_dir),
            "source_row_id": "",
            "parse_status": "missing_directory",
            "message": "Metabolomics h2 directory not found",
        })
        return pd.DataFrame(rows), pd.DataFrame(report), pd.DataFrame(all_blocks)

    for organ in EXPECTED_METABOLOMICS:
        hereg_file = met_dir / organ / f"{organ}.HEreg"

        if not hereg_file.exists():
            report.append({
                "source_type": "Metabolomics",
                "source_file": str(hereg_file),
                "source_row_id": organ,
                "parse_status": "missing_file",
                "message": "Expected .HEreg file not found",
            })
            rows.append({
                "bag_id": make_bag_id(organ, "Metabolomics"),
                "organ": organ,
                "modality": "Metabolomics",
                "method": "HEreg",
                "h2": math.nan,
                "se": math.nan,
                "pvalue": math.nan,
                "n": math.nan,
                "source_type": "HEreg",
                "source_file": str(hereg_file),
                "source_row_id": "V(G)/Vp",
                "parse_status": "missing_file",
                "he_block": preferred_block,
            })
            continue

        row, rep, block_df = parse_hereg_file(
            hereg_file=hereg_file,
            organ=organ,
            modality="Metabolomics",
            preferred_block=preferred_block,
        )

        rows.append(row)
        report.extend(rep)

        if not block_df.empty:
            all_blocks.append(block_df)

    all_blocks_df = pd.concat(all_blocks, ignore_index=True) if all_blocks else pd.DataFrame()

    return pd.DataFrame(rows), pd.DataFrame(report), all_blocks_df


# ============================================================
# 7. Main
# ============================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Collect 23 BAG GCTA/HEreg heritability estimates."
    )

    parser.add_argument(
        "--mri-tsv",
        default=str(DEFAULT_MRI_TSV),
        help="MRI BAG GCTA summary TSV.",
    )

    parser.add_argument(
        "--prot-dir",
        default=str(DEFAULT_PROT_DIR),
        help="Proteomics BAG h2 directory containing organ/*.hsq.",
    )

    parser.add_argument(
        "--met-dir",
        default=str(DEFAULT_MET_DIR),
        help="Metabolomics BAG h2 directory containing organ/*.HEreg.",
    )

    parser.add_argument(
        "--out-dir",
        default=str(DEFAULT_OUT_DIR),
        help="Output directory.",
    )

    parser.add_argument(
        "--met-hereg-block",
        default="HE-CP",
        choices=["HE-CP", "HE-SD"],
        help="Which HEreg block to use for metabolomics main estimate.",
    )

    args = parser.parse_args()

    mri_tsv = Path(args.mri_tsv).expanduser()
    prot_dir = Path(args.prot_dir).expanduser()
    met_dir = Path(args.met_dir).expanduser()
    out_dir = Path(args.out_dir).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("Collecting 23 BAG GCTA/HEreg heritability estimates")
    print("=" * 80)
    print(f"MRI TSV        : {mri_tsv}")
    print(f"Proteomics dir : {prot_dir}")
    print(f"Metabolomics dir: {met_dir}")
    print(f"Met HE block   : {args.met_hereg_block}")
    print(f"Output dir     : {out_dir}")
    print("=" * 80)

    mri_df, mri_report = collect_mri_h2(mri_tsv)
    prot_df, prot_report = collect_proteomics_h2(prot_dir)
    met_df, met_report, met_blocks_df = collect_metabolomics_h2(
        met_dir,
        preferred_block=args.met_hereg_block,
    )

    summary = pd.concat(
        [mri_df, prot_df, met_df],
        ignore_index=True,
        sort=False,
    )

    report = pd.concat(
        [mri_report, prot_report, met_report],
        ignore_index=True,
        sort=False,
    )

    # Standardize final columns and ordering.
    final_cols = [
        "bag_id",
        "organ",
        "modality",
        "method",
        "h2",
        "se",
        "pvalue",
        "n",
        "source_type",
        "source_file",
        "source_row_id",
        "he_block",
        "parse_status",
    ]

    for c in final_cols:
        if c not in summary.columns:
            summary[c] = math.nan if c in {"h2", "se", "pvalue", "n"} else ""

    summary = summary[final_cols].copy()

    modality_order = {
        "MRI": 1,
        "Proteomics": 2,
        "Metabolomics": 3,
    }

    summary["modality_order"] = summary["modality"].map(modality_order).fillna(99)
    summary = summary.sort_values(["modality_order", "organ", "bag_id"]).drop(columns=["modality_order"])

    # Add expected flag.
    expected_bags = (
        [make_bag_id(normalize_organ_for_label(x), "MRI") for x in EXPECTED_MRI]
        + [make_bag_id(x, "Proteomics") for x in EXPECTED_PROTEOMICS]
        + [make_bag_id(x, "Metabolomics") for x in EXPECTED_METABOLOMICS]
    )

    summary["expected_bag"] = summary["bag_id"].isin(expected_bags)

    # Save outputs.
    out_summary = out_dir / "BAG_23_GCTA_h2_summary.tsv"
    out_wide = out_dir / "BAG_23_GCTA_h2_summary_wide.tsv"
    out_report = out_dir / "BAG_23_GCTA_h2_parse_report.tsv"
    out_met_blocks = out_dir / "BAG_23_metabolomics_HEreg_all_blocks.tsv"

    summary.to_csv(out_summary, sep="\t", index=False)

    wide = summary.pivot_table(
        index=["bag_id", "organ", "modality"],
        values=["h2", "se", "pvalue", "n"],
        aggfunc="first",
    ).reset_index()

    wide.to_csv(out_wide, sep="\t", index=False)

    report.to_csv(out_report, sep="\t", index=False)

    if not met_blocks_df.empty:
        met_blocks_df.to_csv(out_met_blocks, sep="\t", index=False)

    # Validation summary.
    n_total = len(summary)
    n_ok = int((summary["parse_status"] == "ok").sum())
    n_expected = len(expected_bags)
    n_unique_bag = summary["bag_id"].nunique()

    missing_expected = sorted(set(expected_bags) - set(summary["bag_id"]))
    bad_rows = summary.loc[summary["parse_status"] != "ok"].copy()

    print("\nSummary:")
    print(f"  Rows collected             : {n_total}")
    print(f"  Unique BAG IDs             : {n_unique_bag}")
    print(f"  Expected BAG count         : {n_expected}")
    print(f"  Rows with parse_status ok  : {n_ok}")
    print("\nCounts by modality:")
    print(summary.groupby(["modality", "parse_status"]).size().reset_index(name="n").to_string(index=False))

    if missing_expected:
        print("\nMissing expected BAG IDs:")
        for x in missing_expected:
            print(f"  - {x}")

    if not bad_rows.empty:
        print("\nRows with non-ok parse_status:")
        print(
            bad_rows[
                ["bag_id", "organ", "modality", "parse_status", "source_file"]
            ].to_string(index=False)
        )

    print("\nSaved:")
    print(f"  {out_summary}")
    print(f"  {out_wide}")
    print(f"  {out_report}")
    if not met_blocks_df.empty:
        print(f"  {out_met_blocks}")

    if n_unique_bag != n_expected:
        print(
            f"\nWARNING: Expected {n_expected} unique BAGs but collected {n_unique_bag}. "
            "Please check the parse report."
        )
    elif n_ok != n_expected:
        print(
            f"\nWARNING: Collected {n_unique_bag} BAGs, but only {n_ok} parsed successfully. "
            "Please check the parse report."
        )
    else:
        print("\nAll 23 expected BAG h2 estimates were collected successfully.")


if __name__ == "__main__":
    main()