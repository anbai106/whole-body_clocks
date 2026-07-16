#!/usr/bin/env python3

import argparse
import os
import re
from pathlib import Path
import pandas as pd
import numpy as np

CYCLES = [
    "1999-2000", "2001-2002", "2003-2004", "2005-2006", "2007-2008",
    "2009-2010", "2011-2012", "2013-2014", "2015-2016", "2017-2018"
]

COMPONENTS = ["Demographics", "Questionnaire", "Examination", "Laboratory", "Dietary"]

DIAGNOSIS_KEYWORDS = [
    "asthma", "copd", "emphysema", "bronchitis", "respiratory",
    "diabetes", "insulin", "glucose", "glycohemoglobin", "hba1c",
    "hypertension", "blood pressure", "cholesterol",
    "heart", "angina", "coronary", "myocardial", "stroke",
    "kidney", "renal", "albumin", "creatinine",
    "liver", "hepatitis", "cancer", "tumor",
    "arthritis", "osteoporosis",
    "depression", "sleep", "smoking", "alcohol",
    "medication", "prescription", "hospital", "health status"
]

def read_xpt(path: Path):
    """
    Read a SAS transport file. Prefer pyreadstat if installed because it can
    expose variable labels. Fallback to pandas.read_sas.
    """
    try:
        import pyreadstat
        df, meta = pyreadstat.read_xport(str(path), metadataonly=False)
        labels = dict(zip(meta.column_names, meta.column_labels))
        return df, labels
    except Exception:
        df = pd.read_sas(path, format="xport", encoding="latin1")
        labels = {c: "" for c in df.columns}
        return df, labels

def safe_nunique(s):
    try:
        return s.nunique(dropna=True)
    except Exception:
        return np.nan

def parse_cycle_from_path(path: Path):
    for c in CYCLES:
        if c in path.parts:
            return c
    return "UNKNOWN"

def parse_component_from_path(path: Path):
    for comp in COMPONENTS:
        if comp in path.parts:
            return comp
    return "UNKNOWN"

def read_mortality_file(path: Path):
    """
    Public-use linked mortality files are fixed-width.
    These column positions follow the NCHS public-use LMF layout.
    """
    colspecs = [
        (0, 6),    # SEQN
        (14, 15),  # ELIGSTAT
        (15, 16),  # MORTSTAT
        (16, 19),  # UCOD_LEADING
        (19, 20),  # DIABETES
        (20, 21),  # HYPERTEN
        (42, 45),  # PERMTH_INT
        (45, 48),  # PERMTH_EXM
    ]
    names = [
        "SEQN", "ELIGSTAT", "MORTSTAT", "UCOD_LEADING",
        "DIABETES", "HYPERTEN", "PERMTH_INT", "PERMTH_EXM"
    ]
    df = pd.read_fwf(path, colspecs=colspecs, names=names, dtype=str)
    for c in ["SEQN", "ELIGSTAT", "MORTSTAT", "DIABETES", "HYPERTEN", "PERMTH_INT", "PERMTH_EXM"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df

def infer_mortality_cycle(filename: str):
    m = re.search(r"NHANES_(\d{4})_(\d{4})_MORT_2019_PUBLIC\.dat", filename)
    if m:
        return f"{m.group(1)}-{m.group(2)}"
    return "UNKNOWN"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nhanes_root", default="/Users/hao/Dropbox/NHANES")
    parser.add_argument("--outdir", default="/Users/hao/Dropbox/NHANES/output_dir/meta_data")
    args = parser.parse_args()

    nhanes_root = Path(args.nhanes_root)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    file_rows = []
    var_rows = []
    sample_rows = []

    xpt_files = sorted(nhanes_root.glob("*/*/*.xpt")) + sorted(nhanes_root.glob("*/*/*.XPT"))

    print(f"Found XPT files: {len(xpt_files)}")

    for i, path in enumerate(xpt_files, start=1):
        cycle = parse_cycle_from_path(path)
        component = parse_component_from_path(path)
        file_name = path.name

        print(f"[{i}/{len(xpt_files)}] Reading {cycle} {component} {file_name}")

        try:
            df, labels = read_xpt(path)
        except Exception as e:
            file_rows.append({
                "cycle": cycle,
                "component": component,
                "file": file_name,
                "path": str(path),
                "n_rows": np.nan,
                "n_cols": np.nan,
                "read_status": "failed",
                "error": str(e)
            })
            continue

        n_rows, n_cols = df.shape

        file_rows.append({
            "cycle": cycle,
            "component": component,
            "file": file_name,
            "path": str(path),
            "n_rows": n_rows,
            "n_cols": n_cols,
            "read_status": "ok",
            "error": ""
        })

        sample_rows.append({
            "cycle": cycle,
            "component": component,
            "file": file_name,
            "n_rows": n_rows,
            "n_cols": n_cols
        })

        for col in df.columns:
            s = df[col]
            label = labels.get(col, "")
            nonmissing = int(s.notna().sum())
            missing = int(s.isna().sum())
            missing_rate = missing / n_rows if n_rows > 0 else np.nan

            text_for_keyword = f"{col} {label}".lower()
            matched_keywords = [k for k in DIAGNOSIS_KEYWORDS if k in text_for_keyword]

            var_rows.append({
                "cycle": cycle,
                "component": component,
                "file": file_name,
                "variable": col,
                "label": label,
                "dtype": str(s.dtype),
                "n_rows": n_rows,
                "nonmissing": nonmissing,
                "missing": missing,
                "missing_rate": missing_rate,
                "n_unique": safe_nunique(s),
                "diagnosis_keyword_hit": int(len(matched_keywords) > 0),
                "matched_keywords": ";".join(matched_keywords)
            })

    file_inv = pd.DataFrame(file_rows)
    var_inv = pd.DataFrame(var_rows)
    sample_summary = pd.DataFrame(sample_rows)

    file_inv.to_csv(outdir / "nhanes_file_inventory.tsv", sep="\t", index=False)
    var_inv.to_csv(outdir / "nhanes_variable_inventory.tsv", sep="\t", index=False)
    sample_summary.to_csv(outdir / "nhanes_sample_size_by_file.tsv", sep="\t", index=False)

    if not var_inv.empty:
        matrix = (
            var_inv
            .assign(present=1)
            .pivot_table(
                index=["component", "file", "variable", "label"],
                columns="cycle",
                values="present",
                aggfunc="max",
                fill_value=0
            )
            .reset_index()
        )
        matrix["n_cycles_present"] = matrix[CYCLES].sum(axis=1)
        matrix.to_csv(outdir / "nhanes_variable_availability_matrix.tsv", sep="\t", index=False)

        diagnosis_candidates = (
            var_inv[var_inv["diagnosis_keyword_hit"] == 1]
            .sort_values(["component", "file", "variable", "cycle"])
        )
        diagnosis_candidates.to_csv(outdir / "nhanes_diagnosis_variable_candidates.tsv", sep="\t", index=False)

        broad_consistent = (
            matrix[matrix["n_cycles_present"] >= 7]
            .sort_values(["component", "file", "n_cycles_present", "variable"], ascending=[True, True, False, True])
        )
        broad_consistent.to_csv(outdir / "nhanes_variables_present_in_at_least_7_cycles.tsv", sep="\t", index=False)

    # Mortality summary
    mort_dir = nhanes_root / "linked_mortality_2019_public"
    mort_files = sorted(mort_dir.glob("NHANES_*_MORT_2019_PUBLIC.dat"))

    mort_summary_rows = []

    for mf in mort_files:
        cycle = infer_mortality_cycle(mf.name)
        try:
            mdf = read_mortality_file(mf)
        except Exception as e:
            mort_summary_rows.append({
                "cycle": cycle,
                "file": mf.name,
                "n": np.nan,
                "eligible_n": np.nan,
                "death_n": np.nan,
                "death_rate": np.nan,
                "median_followup_months_int": np.nan,
                "median_followup_months_exm": np.nan,
                "read_status": "failed",
                "error": str(e)
            })
            continue

        eligible = mdf["ELIGSTAT"].eq(1)
        death = mdf["MORTSTAT"].eq(1)

        mort_summary_rows.append({
            "cycle": cycle,
            "file": mf.name,
            "n": len(mdf),
            "eligible_n": int(eligible.sum()),
            "death_n": int((eligible & death).sum()),
            "death_rate_among_eligible": float((eligible & death).sum() / eligible.sum()) if eligible.sum() > 0 else np.nan,
            "median_followup_months_int": float(mdf.loc[eligible, "PERMTH_INT"].median()),
            "median_followup_months_exm": float(mdf.loc[eligible, "PERMTH_EXM"].median()),
            "read_status": "ok",
            "error": ""
        })

    pd.DataFrame(mort_summary_rows).to_csv(outdir / "nhanes_mortality_summary_by_cycle.tsv", sep="\t", index=False)

    print("\nDone.")
    print(f"Outputs written to: {outdir}")

if __name__ == "__main__":
    main()