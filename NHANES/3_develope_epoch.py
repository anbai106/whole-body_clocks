#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
STEP II: NHANES Model 2 non-disease-input mortality EPOCH

Goal
----
Train a mortality EPOCH model in NHANES using baseline variables that exclude
explicit disease-status modules. The model predicts future all-cause mortality
from baseline survey/exam/lab data.

Design
------
- Outcome: all-cause mortality
- Time: PERMTH_EXM / 12, because this model uses MEC exam/lab variables
- Event: MORTSTAT == 1
- Eligible mortality linkage: ELIGSTAT == 1
- Primary age range: age >= 40
- Explicit disease-status modules excluded:
    MCQ, DIQ, BPQ, CDQ, CKQ, KIQ, HEQ, OSQ
- Main downstream score:
    mortality_epoch_acceleration_z

Main outputs
------------
nhanes_model2_analysis_table.tsv.gz
nhanes_model2_epoch_scores.tsv
nhanes_model2_performance.tsv
nhanes_model2_feature_summary.tsv
nhanes_model2_coefficients.tsv
nhanes_model2_loaded_file_sources.tsv
nhanes_model2_config.json
nhanes_model2_model_bundle.joblib
"""

import argparse
import json
import re
import sys
import warnings
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

try:
    from sksurv.linear_model import CoxnetSurvivalAnalysis
    from sksurv.metrics import concordance_index_censored
    from sksurv.util import Surv
except ImportError as exc:
    raise ImportError(
        "This script requires scikit-survival. Install with:\n"
        "  conda install -c conda-forge scikit-survival\n"
    ) from exc


# ============================================================
# Global configuration
# ============================================================

CYCLES = [
    {
        "cycle": "1999-2000",
        "begin_year": 1999,
        "mort_cycle": "1999_2000",
        "split": "train",
    },
    {
        "cycle": "2001-2002",
        "begin_year": 2001,
        "mort_cycle": "2001_2002",
        "split": "train",
    },
    {
        "cycle": "2003-2004",
        "begin_year": 2003,
        "mort_cycle": "2003_2004",
        "split": "train",
    },
    {
        "cycle": "2005-2006",
        "begin_year": 2005,
        "mort_cycle": "2005_2006",
        "split": "train",
    },
    {
        "cycle": "2007-2008",
        "begin_year": 2007,
        "mort_cycle": "2007_2008",
        "split": "train",
    },
    {
        "cycle": "2009-2010",
        "begin_year": 2009,
        "mort_cycle": "2009_2010",
        "split": "train",
    },
    {
        "cycle": "2011-2012",
        "begin_year": 2011,
        "mort_cycle": "2011_2012",
        "split": "validation",
    },
    {
        "cycle": "2013-2014",
        "begin_year": 2013,
        "mort_cycle": "2013_2014",
        "split": "validation",
    },
    {
        "cycle": "2015-2016",
        "begin_year": 2015,
        "mort_cycle": "2015_2016",
        "split": "test",
    },
    {
        "cycle": "2017-2018",
        "begin_year": 2017,
        "mort_cycle": "2017_2018",
        "split": "test",
    },
]

# Explicit disease-status modules to exclude.
EXCLUDED_DISEASE_MODULE_PATTERNS = [
    r"^MCQ",   # medical conditions: cancer, heart disease, stroke, asthma, etc.
    r"^DIQ",   # diabetes questionnaire
    r"^BPQ",   # blood pressure/cholesterol diagnosis questionnaire
    r"^CDQ",   # cardiovascular questionnaire
    r"^CKQ",   # kidney conditions
    r"^KIQ",   # kidney/urologic conditions
    r"^HEQ",   # hepatitis
    r"^OSQ",   # osteoporosis
]

# Non-disease questionnaire modules allowed.
ALLOWED_QUESTIONNAIRE_MODULE_PATTERNS = [
    r"^ALQ",       # alcohol
    r"^DBQ",       # diet behavior
    r"^DPQ",       # depressive symptoms
    r"^FSQ",       # food security
    r"^HIQ",       # health insurance
    r"^HSQ",       # health status / recent health
    r"^HUQ",       # healthcare utilization/access
    r"^INQ",       # income
    r"^OCQ",       # occupation
    r"^PAQ",       # physical activity
    r"^PFQ",       # physical functioning
    r"^SLQ",       # sleep
    r"^SMQ",       # smoking
    r"^WHQ",       # weight history
    r"^RXQ_RX",    # prescription medications; used as medication count
    r"^RXQASA",    # aspirin use, later cycles
]

DEMO_VARS = [
    "SEQN",
    "SDDSRVYR",
    "RIDSTATR",
    "RIAGENDR",
    "RIDAGEYR",
    "RIDRETH1",
    "DMDEDUC2",
    "INDFMPIR",
    "WTINT2YR",
    "WTMEC2YR",
    "SDMVPSU",
    "SDMVSTRA",
]

EXAM_VARS = [
    "SEQN",
    "BMXWT",
    "BMXHT",
    "BMXBMI",
    "BMXWAIST",
    "BPXSY1",
    "BPXSY2",
    "BPXSY3",
    "BPXSY4",
    "BPXDI1",
    "BPXDI2",
    "BPXDI3",
    "BPXDI4",
    "BPXPLS",
]

LAB_VARS = [
    "SEQN",
    # glycemia/metabolism
    "LBXGH",
    "LBXGLU",
    "LBXIN",
    # lipids
    "LBXTC",
    "LBDHDD",
    "LBDHDL",
    "LBDLDL",
    "LBXTR",
    "LBDTRSI",
    # renal / biochemistry
    "LBXSCR",
    "LBDSCRSI",
    "LBXSBU",
    "LBXSAL",
    "LBXSUA",
    # liver
    "LBXSATSI",
    "LBXSASSI",
    "LBXSGTSI",
    "LBXSAPSI",
    "LBXSTB",
    # hematology
    "LBXHGB",
    "LBXWBCSI",
    "LBXPLTSI",
    "LBXRBCSI",
    "LBXHCT",
    # inflammation
    "LBXCRP",
    "LBDCRP",
    "LBXHSCRP",
    # smoking biomarker
    "LBXCOT",
    # urine kidney markers
    "URXUMA",
    "URXUCR",
]

ID_OUTCOME_DESIGN_VARS = [
    "SEQN",
    "cycle",
    "cycle_begin_year",
    "split",
    "ELIGSTAT",
    "MORTSTAT",
    "death",
    "UCOD_LEADING",
    "DIABETES_MORT_FLAG",
    "HYPERTEN_MORT_FLAG",
    "PERMTH_INT",
    "PERMTH_EXM",
    "followup_years_int",
    "followup_years_exm",
    "eligible_mortality",
    "RIDSTATR",
    "WTINT2YR",
    "WTMEC2YR",
    "SDMVPSU",
    "SDMVSTRA",
    "SDDSRVYR",
]

CORE_COVARS = [
    "RIDAGEYR",
    "RIAGENDR",
    "RIDRETH1",
    "DMDEDUC2",
    "INDFMPIR",
]


# ============================================================
# Helper functions
# ============================================================

def make_one_hot_encoder():
    """Support both newer and older scikit-learn versions."""
    try:
        return OneHotEncoder(handle_unknown="ignore", sparse_output=False)
    except TypeError:
        return OneHotEncoder(handle_unknown="ignore", sparse=False)


def matches_any_pattern(text, patterns):
    text = str(text).upper()
    return any(re.search(pattern, text, flags=re.IGNORECASE) for pattern in patterns)


def module_name(path):
    return Path(path).stem.upper()


def decode_object_columns(df):
    object_cols = df.select_dtypes(include=["object"]).columns
    for col in object_cols:
        df[col] = df[col].map(
            lambda x: x.decode("utf-8", errors="ignore") if isinstance(x, bytes) else x
        )
    return df


def read_xpt_safe(path):
    path = Path(path)
    try:
        try:
            import pyreadstat
            df, _ = pyreadstat.read_xport(str(path))
        except Exception:
            df = pd.read_sas(path, format="xport", encoding="latin1")
    except Exception as exc:
        warnings.warn(f"Could not read XPT file: {path}\nError: {exc}")
        return None

    df = decode_object_columns(df)
    df.columns = [str(c).upper() for c in df.columns]
    return df


def collapse_by_seqn(df):
    if df is None or "SEQN" not in df.columns:
        return None
    if not df["SEQN"].duplicated().any():
        return df
    return df.groupby("SEQN", as_index=False).first()


def safe_merge(left, right):
    if right is None:
        return left
    if left is None:
        return right

    new_cols = [c for c in right.columns if c == "SEQN" or c not in left.columns]
    if len(new_cols) <= 1:
        return left

    return left.merge(right[new_cols], on="SEQN", how="outer")


def to_int_or_nan(x):
    try:
        x = str(x).strip()
        if x == "":
            return np.nan
        return int(x)
    except Exception:
        return np.nan


def read_mortality_file(path):
    """
    Read NHANES public-use fixed-width linked mortality file.

    Column positions:
    SEQN:          1-6
    ELIGSTAT:      15
    MORTSTAT:      16
    UCOD_LEADING:  17-19
    DIABETES:      20
    HYPERTEN:      21
    PERMTH_INT:    43-45
    PERMTH_EXM:    46-48
    """
    records = []
    with open(path, "r", encoding="latin1", errors="ignore") as f:
        for line in f:
            records.append(
                {
                    "SEQN": to_int_or_nan(line[0:6]),
                    "ELIGSTAT": to_int_or_nan(line[14:15]),
                    "MORTSTAT": to_int_or_nan(line[15:16]),
                    "UCOD_LEADING": to_int_or_nan(line[16:19]),
                    "DIABETES_MORT_FLAG": to_int_or_nan(line[19:20]),
                    "HYPERTEN_MORT_FLAG": to_int_or_nan(line[20:21]),
                    "PERMTH_INT": to_int_or_nan(line[42:45]),
                    "PERMTH_EXM": to_int_or_nan(line[45:48]),
                }
            )

    return pd.DataFrame.from_records(records)


def normalize_dtypes(df):
    """
    Convert columns that are mostly numeric to numeric while preserving true
    string/object columns when present.
    """
    for col in df.columns:
        if col in ["cycle", "split"]:
            continue
        if pd.api.types.is_numeric_dtype(df[col]):
            continue

        nonmissing = df[col].notna().sum()
        if nonmissing == 0:
            continue

        numeric = pd.to_numeric(df[col], errors="coerce")
        numeric_nonmissing = numeric.notna().sum()

        if numeric_nonmissing / max(nonmissing, 1) >= 0.90:
            df[col] = numeric

    return df


def clean_special_missing_codes(df):
    """
    Replace common NHANES special missing codes for categorical-like variables.
    This is intentionally conservative: only numeric variables with relatively
    low cardinality are modified.
    """
    special_codes = {
        7,
        9,
        77,
        99,
        777,
        999,
        7777,
        9999,
        77777,
        99999,
    }

    protected = {
        "SEQN",
        "RIDAGEYR",
        "PERMTH_INT",
        "PERMTH_EXM",
        "followup_years_int",
        "followup_years_exm",
    }

    for col in df.columns:
        if col in protected:
            continue
        if not pd.api.types.is_numeric_dtype(df[col]):
            continue

        vals = df[col].dropna()
        if vals.empty:
            continue

        if vals.nunique() <= 50:
            df.loc[df[col].isin(special_codes), col] = np.nan

    return df


def row_mean(df, cols):
    cols = [c for c in cols if c in df.columns]
    if not cols:
        return pd.Series(np.nan, index=df.index)
    return df[cols].mean(axis=1, skipna=True)


def derive_variables(df):
    """Derive interpretable non-disease clinical variables."""
    sys_cols = ["BPXSY1", "BPXSY2", "BPXSY3", "BPXSY4"]
    dia_cols = ["BPXDI1", "BPXDI2", "BPXDI3", "BPXDI4"]

    if any(c in df.columns for c in sys_cols):
        df["mean_systolic_bp"] = row_mean(df, sys_cols)

    if any(c in df.columns for c in dia_cols):
        df["mean_diastolic_bp"] = row_mean(df, dia_cols)

    if {"mean_systolic_bp", "mean_diastolic_bp"}.issubset(df.columns):
        df["pulse_pressure"] = df["mean_systolic_bp"] - df["mean_diastolic_bp"]

    if "HDL_cholesterol" not in df.columns:
        if "LBDHDD" in df.columns:
            df["HDL_cholesterol"] = df["LBDHDD"]
        elif "LBDHDL" in df.columns:
            df["HDL_cholesterol"] = df["LBDHDL"]

    if "creatinine_umol_l" not in df.columns:
        if "LBDSCRSI" in df.columns:
            df["creatinine_umol_l"] = df["LBDSCRSI"]
        elif "LBXSCR" in df.columns:
            df["creatinine_umol_l"] = df["LBXSCR"] * 88.4

    if {"LBXTC", "HDL_cholesterol"}.issubset(df.columns):
        df["non_hdl_cholesterol"] = df["LBXTC"] - df["HDL_cholesterol"]

    if {"URXUMA", "URXUCR"}.issubset(df.columns):
        df["urine_albumin_creatinine_ratio"] = np.where(
            df["URXUCR"] > 0,
            df["URXUMA"] / df["URXUCR"],
            np.nan,
        )

    if {"LBXSASSI", "LBXSATSI"}.issubset(df.columns):
        df["ast_alt_ratio"] = np.where(
            df["LBXSATSI"] > 0,
            df["LBXSASSI"] / df["LBXSATSI"],
            np.nan,
        )

    # CKD-EPI 2009 eGFR, derived from creatinine. This is a lab-derived renal
    # function biomarker, not self-reported disease status.
    if {"creatinine_umol_l", "RIDAGEYR", "RIAGENDR", "RIDRETH1"}.issubset(df.columns):
        scr_mg_dl = df["creatinine_umol_l"] / 88.4
        female = df["RIAGENDR"] == 2
        black = df["RIDRETH1"] == 4

        kappa = np.where(female, 0.7, 0.9)
        alpha = np.where(female, -0.329, -0.411)
        scr_k = scr_mg_dl / kappa

        egfr = (
            141
            * np.power(np.minimum(scr_k, 1), alpha)
            * np.power(np.maximum(scr_k, 1), -1.209)
            * np.power(0.993, df["RIDAGEYR"])
            * np.where(female, 1.018, 1.0)
            * np.where(black, 1.159, 1.0)
        )

        egfr = np.where(np.isfinite(egfr), egfr, np.nan)
        df["egfr_ckdepi_2009"] = egfr

    return df


def read_cycle_data(nhanes_root, cycle_info):
    cycle = cycle_info["cycle"]
    begin_year = cycle_info["begin_year"]
    mort_cycle = cycle_info["mort_cycle"]
    split = cycle_info["split"]

    cycle_dir = Path(nhanes_root) / cycle
    mort_file = (
        Path(nhanes_root)
        / "linked_mortality_2019_public"
        / f"NHANES_{mort_cycle}_MORT_2019_PUBLIC.dat"
    )

    if not cycle_dir.exists():
        raise FileNotFoundError(f"Missing NHANES cycle directory: {cycle_dir}")

    if not mort_file.exists():
        raise FileNotFoundError(f"Missing mortality file: {mort_file}")

    print(f"\nReading cycle: {cycle}")

    xpt_files = sorted([p for p in cycle_dir.rglob("*") if p.suffix.lower() == ".xpt"])

    merged = None
    source_rows = []

    for path in xpt_files:
        comp = path.parent.name
        mod = module_name(path)

        read_this = False
        keep_vars = None
        rx_count_only = False

        if comp == "Demographics":
            read_this = True
            keep_vars = DEMO_VARS

        elif comp == "Examination":
            read_this = True
            keep_vars = EXAM_VARS

        elif comp == "Laboratory":
            read_this = True
            keep_vars = LAB_VARS

        elif comp == "Questionnaire":
            if matches_any_pattern(mod, EXCLUDED_DISEASE_MODULE_PATTERNS):
                read_this = False
            elif re.search(r"^RXQ_RX", mod, flags=re.IGNORECASE):
                read_this = True
                rx_count_only = True
            elif matches_any_pattern(mod, ALLOWED_QUESTIONNAIRE_MODULE_PATTERNS):
                read_this = True
                keep_vars = None

        if not read_this:
            continue

        dt = read_xpt_safe(path)
        if dt is None or "SEQN" not in dt.columns:
            continue

        if rx_count_only:
            dt2 = (
                dt.groupby("SEQN", as_index=False)
                .size()
                .rename(columns={"size": "RX_MED_COUNT"})
            )
        elif keep_vars is not None:
            present = [v for v in keep_vars if v in dt.columns]
            if len(present) <= 1:
                continue
            dt2 = dt[present].copy()
        else:
            dt2 = dt.copy()

        dt2 = collapse_by_seqn(dt2)
        if dt2 is None:
            continue

        source_rows.append(
            {
                "cycle": cycle,
                "component": comp,
                "file": path.name,
                "module": mod,
                "n_variables_loaded": len([c for c in dt2.columns if c != "SEQN"]),
                "variables_loaded": ";".join([c for c in dt2.columns if c != "SEQN"]),
            }
        )

        merged = safe_merge(merged, dt2)

    if merged is None:
        raise RuntimeError(f"No usable data loaded for cycle: {cycle}")

    mort = read_mortality_file(mort_file)

    merged = merged.merge(mort, on="SEQN", how="left")
    merged["cycle"] = cycle
    merged["cycle_begin_year"] = begin_year
    merged["split"] = split
    merged["eligible_mortality"] = merged["ELIGSTAT"] == 1
    merged["death"] = (merged["MORTSTAT"] == 1).astype(float)
    merged["followup_years_int"] = merged["PERMTH_INT"] / 12.0
    merged["followup_years_exm"] = merged["PERMTH_EXM"] / 12.0

    if "RX_MED_COUNT" in merged.columns:
        merged["RX_MED_COUNT"] = merged["RX_MED_COUNT"].fillna(0)

    merged = normalize_dtypes(merged)
    merged = clean_special_missing_codes(merged)
    merged = derive_variables(merged)

    return merged, pd.DataFrame(source_rows)


def make_feature_summary(df, candidate_vars, train_mask):
    train = df.loc[train_mask].copy()
    rows = []

    for var in candidate_vars:
        x = train[var]
        rows.append(
            {
                "variable": var,
                "dtype": str(x.dtype),
                "n_train": len(x),
                "nonmissing_train": int(x.notna().sum()),
                "missing_train": int(x.isna().sum()),
                "missing_rate_train": float(x.isna().mean()),
                "n_unique_train": int(x.nunique(dropna=True)),
            }
        )

    return pd.DataFrame(rows)


def source_variable_from_feature_name(feature_name, model_vars):
    """
    Map transformed feature name back to original variable name.
    Examples:
      num__RIDAGEYR -> RIDAGEYR
      cat__RIAGENDR_1.0 -> RIAGENDR
    """
    clean = feature_name
    if "__" in clean:
        clean = clean.split("__", 1)[1]

    for var in sorted(model_vars, key=len, reverse=True):
        if clean == var or clean.startswith(f"{var}_"):
            return var

    return clean


def compute_cindex(time, event, risk):
    mask = np.isfinite(time) & np.isfinite(event) & np.isfinite(risk)
    if mask.sum() < 20:
        return np.nan
    if np.sum(event[mask] == 1) < 5:
        return np.nan

    try:
        return float(
            concordance_index_censored(
                event_indicator=event[mask].astype(bool),
                event_time=time[mask].astype(float),
                estimate=risk[mask].astype(float),
            )[0]
        )
    except Exception:
        return np.nan


def prepare_model_dataframe(df, model_vars, categorical_vars, numeric_vars):
    out = df[model_vars].copy()

    for col in numeric_vars:
        out[col] = pd.to_numeric(out[col], errors="coerce")

    for col in categorical_vars:
        out[col] = out[col].astype("object")
        out[col] = out[col].where(out[col].notna(), "Missing")
        out[col] = out[col].astype(str)

    return out


def build_residualizer_inputs(score_df, categorical_vars, numeric_vars):
    tmp = score_df[categorical_vars + numeric_vars].copy()

    for col in numeric_vars:
        tmp[col] = pd.to_numeric(tmp[col], errors="coerce")

    for col in categorical_vars:
        tmp[col] = tmp[col].astype("object")
        tmp[col] = tmp[col].where(tmp[col].notna(), "Missing")
        tmp[col] = tmp[col].astype(str)

    return tmp


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--nhanes_root", default="/Users/hao/Dropbox/NHANES",
        help="Root directory, e.g. /Users/hao/Dropbox/NHANES",
    )
    parser.add_argument(
        "--outdir",
            default="/Users/hao/Dropbox/NHANES/output_dir/model2_nondisease_mortality_epoch",
        help="Output directory for Model 2 NHANES EPOCH",
    )
    parser.add_argument("--min_age", type=float, default=40.0)
    parser.add_argument("--max_missing_rate_train", type=float, default=0.60)
    parser.add_argument("--min_nonmissing_train", type=int, default=200)
    parser.add_argument("--l1_ratio", type=float, default=0.50)
    parser.add_argument("--n_alphas", type=int, default=100)
    parser.add_argument("--alpha_min_ratio", type=float, default=0.01)
    parser.add_argument("--max_iter", type=int, default=100000)
    parser.add_argument(
        "--validation_tie_margin",
        type=float,
        default=0.002,
        help="Choose the sparsest alpha within this C-index margin of the best validation C-index.",
    )

    args = parser.parse_args()

    nhanes_root = Path(args.nhanes_root)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("STEP II: NHANES Model 2 non-disease-input mortality EPOCH")
    print("=" * 80)
    print(f"NHANES root: {nhanes_root}")
    print(f"Output dir:  {outdir}")
    print(f"Minimum age: {args.min_age}")
    print("=" * 80)

    # --------------------------------------------------------
    # Read all cycles
    # --------------------------------------------------------
    cycle_tables = []
    source_tables = []

    for ci in CYCLES:
        dt, src = read_cycle_data(nhanes_root, ci)
        cycle_tables.append(dt)
        source_tables.append(src)

    nhanes = pd.concat(cycle_tables, axis=0, ignore_index=True, sort=False)
    sources = pd.concat(source_tables, axis=0, ignore_index=True, sort=False)

    sources.to_csv(
        outdir / "nhanes_model2_loaded_file_sources.tsv",
        sep="\t",
        index=False,
    )

    print("\nRaw merged NHANES table")
    print(f"  Rows:    {nhanes.shape[0]}")
    print(f"  Columns: {nhanes.shape[1]}")

    # --------------------------------------------------------
    # Analytic sample
    # --------------------------------------------------------
    keep = (
        (nhanes["eligible_mortality"] == True)
        & nhanes["death"].notna()
        & nhanes["followup_years_exm"].notna()
        & (nhanes["followup_years_exm"] > 0)
        & nhanes["RIDAGEYR"].notna()
        & (nhanes["RIDAGEYR"] >= args.min_age)
    )

    if "RIDSTATR" in nhanes.columns:
        keep = keep & (nhanes["RIDSTATR"].isna() | (nhanes["RIDSTATR"] == 2))

    nhanes = nhanes.loc[keep].copy()
    nhanes = normalize_dtypes(nhanes)

    print("\nAnalytic sample")
    print(f"  Rows:   {nhanes.shape[0]}")
    print(f"  Deaths: {int((nhanes['death'] == 1).sum())}")

    nhanes.to_csv(
        outdir / "nhanes_model2_analysis_table.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    # --------------------------------------------------------
    # Candidate predictor selection
    # --------------------------------------------------------
    present_core_covars = [c for c in CORE_COVARS if c in nhanes.columns]
    missing_core_covars = [c for c in CORE_COVARS if c not in nhanes.columns]

    if missing_core_covars:
        print(f"WARNING: missing core covariates: {missing_core_covars}")

    drop_vars = set(ID_OUTCOME_DESIGN_VARS + present_core_covars)

    candidate_vars = [
        c for c in nhanes.columns
        if c not in drop_vars
        and not matches_any_pattern(c, EXCLUDED_DISEASE_MODULE_PATTERNS)
        and not c.startswith("UCOD")
        and c not in ["DIABETES_MORT_FLAG", "HYPERTEN_MORT_FLAG"]
    ]

    train_mask = nhanes["split"] == "train"
    validation_mask = nhanes["split"] == "validation"
    test_mask = nhanes["split"] == "test"

    feature_summary = make_feature_summary(nhanes, candidate_vars, train_mask)

    feature_summary["selected"] = (
        (feature_summary["nonmissing_train"] >= args.min_nonmissing_train)
        & (feature_summary["missing_rate_train"] <= args.max_missing_rate_train)
        & (feature_summary["n_unique_train"] >= 2)
    )

    selected_features = feature_summary.loc[
        feature_summary["selected"], "variable"
    ].tolist()

    feature_summary.to_csv(
        outdir / "nhanes_model2_feature_summary.tsv",
        sep="\t",
        index=False,
    )

    print("\nFeature selection")
    print(f"  Candidate variables: {len(candidate_vars)}")
    print(f"  Selected variables:  {len(selected_features)}")

    # --------------------------------------------------------
    # Build model matrix
    # --------------------------------------------------------
    model_vars = present_core_covars + selected_features
    model_vars = list(dict.fromkeys(model_vars))

    # Core categorical variables. Age and INDFMPIR are kept numeric.
    core_categorical = [c for c in ["RIAGENDR", "RIDRETH1", "DMDEDUC2"] if c in model_vars]

    # Non-core object variables become categorical. Most NHANES variables are numeric.
    object_categorical = [
        c for c in selected_features
        if c in nhanes.columns and not pd.api.types.is_numeric_dtype(nhanes[c])
    ]

    categorical_vars = list(dict.fromkeys(core_categorical + object_categorical))
    numeric_vars = [c for c in model_vars if c not in categorical_vars]

    print("\nModel variable types")
    print(f"  Numeric variables:     {len(numeric_vars)}")
    print(f"  Categorical variables: {len(categorical_vars)}")

    model_input = prepare_model_dataframe(
        nhanes,
        model_vars=model_vars,
        categorical_vars=categorical_vars,
        numeric_vars=numeric_vars,
    )

    numeric_pipeline = Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
        ]
    )

    categorical_pipeline = Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="constant", fill_value="Missing")),
            ("onehot", make_one_hot_encoder()),
        ]
    )

    transformers = []
    if numeric_vars:
        transformers.append(("num", numeric_pipeline, numeric_vars))
    if categorical_vars:
        transformers.append(("cat", categorical_pipeline, categorical_vars))

    preprocessor = ColumnTransformer(
        transformers=transformers,
        remainder="drop",
        verbose_feature_names_out=True,
    )

    X_train_input = model_input.loc[train_mask, :]
    preprocessor.fit(X_train_input)

    X_all = preprocessor.transform(model_input)
    X_all = np.asarray(X_all, dtype=np.float64)

    try:
        feature_names = preprocessor.get_feature_names_out()
    except Exception:
        feature_names = np.array([f"x{i}" for i in range(X_all.shape[1])])

    X_train = X_all[train_mask.to_numpy(), :]

    # Remove transformed columns with zero variance in training.
    train_sd = np.nanstd(X_train, axis=0)
    keep_cols = train_sd > 1e-12
    X_all = X_all[:, keep_cols]
    X_train = X_train[:, keep_cols]
    feature_names = np.asarray(feature_names)[keep_cols]

    print("\nTransformed design matrix")
    print(f"  Rows:    {X_all.shape[0]}")
    print(f"  Columns: {X_all.shape[1]}")

    # Penalty factor: core covariates unpenalized, EPOCH features penalized.
    source_vars = np.array(
        [source_variable_from_feature_name(f, model_vars) for f in feature_names]
    )
    penalty_factor = np.where(np.isin(source_vars, present_core_covars), 0.0, 1.0)

    print("\nPenalty structure")
    print(f"  Unpenalized transformed core covariate columns: {int((penalty_factor == 0).sum())}")
    print(f"  Penalized transformed EPOCH feature columns:    {int((penalty_factor == 1).sum())}")

    y_train = Surv.from_arrays(
        event=nhanes.loc[train_mask, "death"].astype(bool).to_numpy(),
        time=nhanes.loc[train_mask, "followup_years_exm"].astype(float).to_numpy(),
    )

    # --------------------------------------------------------
    # Train Cox elastic-net
    # --------------------------------------------------------
    print("\nTraining Cox elastic-net model...")

    coxnet = CoxnetSurvivalAnalysis(
        l1_ratio=args.l1_ratio,
        n_alphas=args.n_alphas,
        alpha_min_ratio=args.alpha_min_ratio,
        penalty_factor=penalty_factor,
        normalize=False,
        max_iter=args.max_iter,
        fit_baseline_model=False,
    )

    coxnet.fit(X_train, y_train)

    coef_path = np.asarray(coxnet.coef_)
    if coef_path.ndim == 1:
        coef_path = coef_path.reshape((-1, 1))

    alphas = np.asarray(coxnet.alphas_)

    # --------------------------------------------------------
    # Select alpha using validation C-index
    # --------------------------------------------------------
    X_val = X_all[validation_mask.to_numpy(), :]
    val_time = nhanes.loc[validation_mask, "followup_years_exm"].to_numpy(dtype=float)
    val_event = nhanes.loc[validation_mask, "death"].to_numpy(dtype=float)

    val_cindex = []
    for j in range(coef_path.shape[1]):
        risk = X_val @ coef_path[:, j]
        val_cindex.append(compute_cindex(val_time, val_event, risk))

    val_cindex = np.asarray(val_cindex, dtype=float)

    if np.all(~np.isfinite(val_cindex)):
        warnings.warn("Validation C-index is unavailable. Falling back to alpha with largest deviance ratio.")
        if hasattr(coxnet, "deviance_ratio_"):
            best_idx = int(np.nanargmax(coxnet.deviance_ratio_))
        else:
            best_idx = coef_path.shape[1] // 2
    else:
        best_c = np.nanmax(val_cindex)
        eligible = np.where(val_cindex >= best_c - args.validation_tie_margin)[0]
        # Coxnet alphas usually decrease from strong to weak regularization.
        # Pick the largest alpha among near-best models for parsimony.
        best_idx = eligible[np.argmax(alphas[eligible])]

    alpha_use = float(alphas[best_idx])
    beta = coef_path[:, best_idx]

    print("\nSelected alpha")
    print(f"  alpha index: {best_idx}")
    print(f"  alpha:       {alpha_use}")
    print(f"  validation C-index: {val_cindex[best_idx]:.4f}")

    # --------------------------------------------------------
    # Compute EPOCH scores
    # --------------------------------------------------------
    lp_total = X_all @ beta

    feature_mask = penalty_factor == 1
    lp_feature_only = X_all[:, feature_mask] @ beta[feature_mask]

    # Convert feature-only log-risk to mortality-risk-equivalent years.
    beta_age_raw = np.nan
    beta_age_scaled = np.nan

    if "RIDAGEYR" in numeric_vars:
        age_feature_name = "num__RIDAGEYR"
        age_idx_transformed = np.where(feature_names == age_feature_name)[0]

        if len(age_idx_transformed) == 1:
            beta_age_scaled = float(beta[age_idx_transformed[0]])

            try:
                num_pipe = preprocessor.named_transformers_["num"]
                scaler = num_pipe.named_steps["scaler"]
                age_idx_raw = numeric_vars.index("RIDAGEYR")
                age_sd_raw = float(scaler.scale_[age_idx_raw])
                beta_age_raw = beta_age_scaled / age_sd_raw
            except Exception as exc:
                warnings.warn(f"Could not derive raw age coefficient: {exc}")

    if np.isfinite(beta_age_raw) and abs(beta_age_raw) > 1e-10:
        epoch_year_equivalent = lp_feature_only / beta_age_raw
        residual_target_name = "mortality_epoch_year_equivalent"
    else:
        warnings.warn(
            "Age coefficient unavailable or near zero. "
            "Using feature-only log-risk for residualized acceleration."
        )
        epoch_year_equivalent = np.full_like(lp_feature_only, np.nan, dtype=float)
        residual_target_name = "mortality_epoch_lp_feature_only"

    score_cols = [
        "SEQN",
        "cycle",
        "cycle_begin_year",
        "split",
        "death",
        "followup_years_exm",
        "RIDAGEYR",
        "RIAGENDR",
        "RIDRETH1",
        "DMDEDUC2",
        "INDFMPIR",
        "WTMEC2YR",
        "SDMVPSU",
        "SDMVSTRA",
    ]

    score_cols = [c for c in score_cols if c in nhanes.columns]
    scores = nhanes[score_cols].copy()

    scores["mortality_epoch_lp_total"] = lp_total
    scores["mortality_epoch_lp_feature_only"] = lp_feature_only
    scores["mortality_epoch_year_equivalent"] = epoch_year_equivalent

    # --------------------------------------------------------
    # Residualize EPOCH score for downstream disease analyses
    # --------------------------------------------------------
    residual_categorical = [
        c for c in ["RIAGENDR", "RIDRETH1", "DMDEDUC2", "cycle"]
        if c in scores.columns
    ]
    residual_numeric = [
        c for c in ["RIDAGEYR", "INDFMPIR"]
        if c in scores.columns
    ]

    residual_target = (
        scores["mortality_epoch_year_equivalent"].copy()
        if residual_target_name == "mortality_epoch_year_equivalent"
        else scores["mortality_epoch_lp_feature_only"].copy()
    )

    resid_input = build_residualizer_inputs(
        scores,
        categorical_vars=residual_categorical,
        numeric_vars=residual_numeric,
    )

    residual_numeric_pipeline = Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
        ]
    )

    residual_categorical_pipeline = Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="constant", fill_value="Missing")),
            ("onehot", make_one_hot_encoder()),
        ]
    )

    residual_transformers = []
    if residual_numeric:
        residual_transformers.append(("num", residual_numeric_pipeline, residual_numeric))
    if residual_categorical:
        residual_transformers.append(("cat", residual_categorical_pipeline, residual_categorical))

    residual_preprocessor = ColumnTransformer(
        transformers=residual_transformers,
        remainder="drop",
        verbose_feature_names_out=True,
    )

    valid_resid = np.isfinite(residual_target.to_numpy(dtype=float))

    # Residualization uses all analysis rows and no mortality outcome.
    # This removes age/demographic/cycle structure from the score for downstream
    # disease-linking analyses.
    residual_preprocessor.fit(resid_input.loc[valid_resid, :])
    R_all = residual_preprocessor.transform(resid_input)
    R_all = np.asarray(R_all, dtype=float)

    residual_model = LinearRegression()
    residual_model.fit(R_all[valid_resid, :], residual_target.loc[valid_resid].to_numpy(dtype=float))

    pred_resid = residual_model.predict(R_all)
    acceleration = residual_target.to_numpy(dtype=float) - pred_resid

    scores["mortality_epoch_acceleration_years"] = acceleration

    train_acc = scores.loc[scores["split"] == "train", "mortality_epoch_acceleration_years"]
    train_acc_mean = float(np.nanmean(train_acc))
    train_acc_sd = float(np.nanstd(train_acc, ddof=1))

    if np.isfinite(train_acc_sd) and train_acc_sd > 0:
        scores["mortality_epoch_acceleration_z"] = (
            scores["mortality_epoch_acceleration_years"] - train_acc_mean
        ) / train_acc_sd
    else:
        scores["mortality_epoch_acceleration_z"] = np.nan

    # --------------------------------------------------------
    # Performance
    # --------------------------------------------------------
    perf_rows = []

    for split_name in ["train", "validation", "test", "all"]:
        if split_name == "all":
            idx = np.ones(len(scores), dtype=bool)
        else:
            idx = scores["split"].to_numpy() == split_name

        d = scores.loc[idx]

        perf_rows.append(
            {
                "split": split_name,
                "n": int(d.shape[0]),
                "deaths": int((d["death"] == 1).sum()),
                "median_followup_years": float(np.nanmedian(d["followup_years_exm"])),
                "cindex_lp_total": compute_cindex(
                    d["followup_years_exm"].to_numpy(dtype=float),
                    d["death"].to_numpy(dtype=float),
                    d["mortality_epoch_lp_total"].to_numpy(dtype=float),
                ),
                "cindex_lp_feature_only": compute_cindex(
                    d["followup_years_exm"].to_numpy(dtype=float),
                    d["death"].to_numpy(dtype=float),
                    d["mortality_epoch_lp_feature_only"].to_numpy(dtype=float),
                ),
                "cindex_acceleration_z": compute_cindex(
                    d["followup_years_exm"].to_numpy(dtype=float),
                    d["death"].to_numpy(dtype=float),
                    d["mortality_epoch_acceleration_z"].to_numpy(dtype=float),
                ),
            }
        )

    performance = pd.DataFrame(perf_rows)

    # --------------------------------------------------------
    # Coefficient table
    # --------------------------------------------------------
    coef_df = pd.DataFrame(
        {
            "transformed_feature": feature_names,
            "source_variable": source_vars,
            "beta": beta,
            "penalty_factor": penalty_factor,
            "is_core_covariate": penalty_factor == 0,
            "is_epoch_feature": penalty_factor == 1,
            "nonzero": beta != 0,
            "alpha_used": alpha_use,
        }
    )

    coef_df = coef_df.sort_values(
        ["nonzero", "is_core_covariate", "source_variable", "transformed_feature"],
        ascending=[False, False, True, True],
    )

    # --------------------------------------------------------
    # Save outputs
    # --------------------------------------------------------
    scores.to_csv(outdir / "nhanes_model2_epoch_scores.tsv", sep="\t", index=False)
    performance.to_csv(outdir / "nhanes_model2_performance.tsv", sep="\t", index=False)
    coef_df.to_csv(outdir / "nhanes_model2_coefficients.tsv", sep="\t", index=False)

    alpha_df = pd.DataFrame(
        {
            "alpha_index": np.arange(len(alphas)),
            "alpha": alphas,
            "validation_cindex": val_cindex,
            "selected": np.arange(len(alphas)) == best_idx,
        }
    )
    alpha_df.to_csv(outdir / "nhanes_model2_alpha_path_validation.tsv", sep="\t", index=False)

    config = {
        "model_name": "NHANES_Model2_non_disease_input_mortality_EPOCH",
        "nhanes_root": str(nhanes_root),
        "outdir": str(outdir),
        "min_age": args.min_age,
        "time_origin": "MEC exam date",
        "time_variable": "PERMTH_EXM / 12",
        "event_variable": "MORTSTAT == 1",
        "eligibility_variable": "ELIGSTAT == 1",
        "excluded_disease_module_patterns": EXCLUDED_DISEASE_MODULE_PATTERNS,
        "allowed_questionnaire_module_patterns": ALLOWED_QUESTIONNAIRE_MODULE_PATTERNS,
        "core_covariates": present_core_covars,
        "n_selected_features": len(selected_features),
        "l1_ratio": args.l1_ratio,
        "n_alphas": args.n_alphas,
        "alpha_min_ratio": args.alpha_min_ratio,
        "alpha_used": alpha_use,
        "alpha_index_used": int(best_idx),
        "validation_cindex_at_selected_alpha": (
            None if not np.isfinite(val_cindex[best_idx]) else float(val_cindex[best_idx])
        ),
        "beta_age_scaled": beta_age_scaled,
        "beta_age_raw_per_year": beta_age_raw,
        "residual_target": residual_target_name,
        "residualization_covariates": residual_numeric + residual_categorical,
        "residualization_fit_scope": "all_analysis_rows_no_mortality_outcome_used",
        "train_acceleration_mean": train_acc_mean,
        "train_acceleration_sd": train_acc_sd,
        "split_definition": CYCLES,
    }

    with open(outdir / "nhanes_model2_config.json", "w") as f:
        json.dump(config, f, indent=2)

    joblib.dump(
        {
            "preprocessor": preprocessor,
            "coxnet_model": coxnet,
            "alpha_used": alpha_use,
            "alpha_index_used": best_idx,
            "beta_used": beta,
            "feature_names": feature_names,
            "source_variables": source_vars,
            "penalty_factor": penalty_factor,
            "model_vars": model_vars,
            "numeric_vars": numeric_vars,
            "categorical_vars": categorical_vars,
            "selected_features": selected_features,
            "core_covariates": present_core_covars,
            "residual_preprocessor": residual_preprocessor,
            "residual_model": residual_model,
            "config": config,
        },
        outdir / "nhanes_model2_model_bundle.joblib",
    )

    print("\nPerformance:")
    print(performance.to_string(index=False))

    print("\nSaved outputs:")
    print(f"  {outdir / 'nhanes_model2_epoch_scores.tsv'}")
    print(f"  {outdir / 'nhanes_model2_performance.tsv'}")
    print(f"  {outdir / 'nhanes_model2_coefficients.tsv'}")
    print(f"  {outdir / 'nhanes_model2_model_bundle.joblib'}")

    print("\nRecommended downstream score:")
    print("  mortality_epoch_acceleration_z")


if __name__ == "__main__":
    main()