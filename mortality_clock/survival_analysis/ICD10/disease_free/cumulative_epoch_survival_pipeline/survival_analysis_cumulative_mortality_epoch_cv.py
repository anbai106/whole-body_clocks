#!/usr/bin/env python3
"""
STEP 3: cumulative Cox survival analysis for proteomics and metabolomics
mortality EPOCH clocks predicting incident disease onset.

This script is designed to run after the individual clock-vs-BAG survival
analysis has already been completed for each ICD endpoint. For each disease:

1. Read the disease-specific ICD clock TSV.
2. Read the existing individual clock survival result file.
3. Keep the 11 proteomics mortality EPOCH clocks and 4 metabolomics mortality
   EPOCH clocks.
4. Rank clocks as follows:
     a) significant clocks first, using the individual clock P value and HR;
     b) non-significant clocks after the significant clocks in a reproducible
        random order.
5. Fit cumulative Cox models using the same complete-case analysis population
   for all cumulative steps.
6. Report apparent C-index, delta C-index, likelihood-ratio tests, and added-clock HRs.
7. Perform Option 1 five-fold cross-validation: keep the full-sample clock
   ranking fixed, fit each cumulative Cox model inside each training fold, and
   evaluate out-of-fold Harrell C-index in the held-out fold.

The default significance threshold for ranking is P < 0.05. To use Bonferroni
within the 15 clocks, run with --rank-p-threshold 0.0033333333.

Cross-validation uses five folds by default and does not add any command-line
arguments. To disable or change folds without changing the Slurm/bash interface,
set the environment variable CV_FOLDS, for example CV_FOLDS=0 or CV_FOLDS=10.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import os
import sys
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from scipy.stats import chi2

warnings.filterwarnings("ignore")

# -----------------------------------------------------------------------------
# Clock list: 11 proteomics mortality EPOCH clocks + 4 metabolomics EPOCH clocks.
# These names match the columns used in the disease clock TSVs and the existing
# clock-vs-BAG outputs.
# -----------------------------------------------------------------------------
CLOCKS: List[Dict[str, str]] = [
    {"pair_id": "Reproductive_female_proteomics", "organ": "Reproductive_female", "modality": "Proteomics", "clock": "Reproductive_female_proteomics"},
    {"pair_id": "Pulmonary_proteomics", "organ": "Pulmonary", "modality": "Proteomics", "clock": "Pulmonary_proteomics"},
    {"pair_id": "Heart_proteomics", "organ": "Heart", "modality": "Proteomics", "clock": "Heart_proteomics"},
    {"pair_id": "Brain_proteomics", "organ": "Brain", "modality": "Proteomics", "clock": "Brain_proteomics"},
    {"pair_id": "Eye_proteomics", "organ": "Eye", "modality": "Proteomics", "clock": "Eye_proteomics"},
    {"pair_id": "Hepatic_proteomics", "organ": "Hepatic", "modality": "Proteomics", "clock": "Hepatic_proteomics"},
    {"pair_id": "Renal_proteomics", "organ": "Renal", "modality": "Proteomics", "clock": "Renal_proteomics"},
    {"pair_id": "Reproductive_male_proteomics", "organ": "Reproductive_male", "modality": "Proteomics", "clock": "Reproductive_male_proteomics"},
    {"pair_id": "Endocrine_proteomics", "organ": "Endocrine", "modality": "Proteomics", "clock": "Endocrine_proteomics"},
    {"pair_id": "Immune_proteomics", "organ": "Immune", "modality": "Proteomics", "clock": "Immune_proteomics"},
    {"pair_id": "Skin_proteomics", "organ": "Skin", "modality": "Proteomics", "clock": "Skin_proteomics"},
    {"pair_id": "Endocrine_metabolomics", "organ": "Endocrine", "modality": "Metabolomics", "clock": "Endocrine_metabolomics"},
    {"pair_id": "Digestive_metabolomics", "organ": "Digestive", "modality": "Metabolomics", "clock": "Digestive_metabolomics"},
    {"pair_id": "Hepatic_metabolomics", "organ": "Hepatic", "modality": "Metabolomics", "clock": "Hepatic_metabolomics"},
    {"pair_id": "Immune_metabolomics", "organ": "Immune", "modality": "Metabolomics", "clock": "Immune_metabolomics"},
]

CLOCK_COLS = [x["clock"] for x in CLOCKS]
CLOCK_META_BY_CLOCK = {x["clock"]: x for x in CLOCKS}
CLOCK_META_BY_PAIR = {x["pair_id"]: x for x in CLOCKS}

BASELINE_DATE_COL = "date_of_attending_assessment_centre_f53_0_0"
IMAGING_DATE_COL = "date_of_attending_assessment_centre_f53_2_0"
DEATH_DATE_COL = "death_date_f40000_0_0"

UMEL_BASELINE_DATE_RAW = "53-0.0"
UMEL_IMAGING_DATE_RAW = "53-2.0"
UMEL_DEATH_DATE_RAW = "40000-0.0"

AGE_RECRUIT_COL = "age_at_recruitment_f21022_0_0"
SMOKING_COL = "smoking_status_f20116_0_0"
BMI_COL = "body_mass_index_bmi_f23104_0_0"
DIASTOLIC_COL = "diastolic_blood_pressure_automated_reading_f4079_0_0"
SYSTOLIC_COL = "systolic_blood_pressure_automated_reading_f4080_0_0"

TIME_COL = "time_baseline"
EVENT_COL = "case"

# Option 1 cumulative cross-validation settings. This keeps the full-sample
# clock ranking fixed, then estimates discrimination using out-of-fold
# predictions from Cox models fit within each training fold.
# No command-line arguments are added; override only through the environment.
CV_FOLDS = int(os.getenv("CV_FOLDS", "5"))
CV_SEED_OFFSET = 7919


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Cumulative Cox survival analysis for the 11 proteomics and 4 "
            "metabolomics mortality EPOCH clocks."
        )
    )
    p.add_argument("--icd_tsv", required=True, help="Disease-specific ICD clock TSV.")
    p.add_argument("--individual_result_tsv", required=True, help="Existing clock-vs-BAG individual survival result TSV for this disease.")
    p.add_argument("--output_tsv", required=True, help="Output cumulative Cox TSV.")
    p.add_argument("--rank_output_tsv", default=None, help="Optional output TSV with the clock ranking order.")
    p.add_argument("--audit_tsv", default=None, help="Optional output TSV with merge/QC audit counts.")

    p.add_argument(
        "--cov_tsv",
        default="/cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv",
        help="UKBB covariate CSV containing eid and common covariates.",
    )
    p.add_argument(
        "--date_tsv",
        default="/cbica/home/wenju/Reproducibile_paper/Multiorgan_Subtype/data/PWAS/UKBB_fullsample_death_variables.csv",
        help="Fallback date file containing UKBB assessment/death dates.",
    )
    p.add_argument(
        "--umel_death_xlsx",
        default="/cbica/home/wenju/Dataset/UKBB_UMelbourne/Death_related_var_from_Ye.xlsx",
        help="UMelbourne UKBB Excel file containing field 53 and death dates.",
    )
    p.add_argument(
        "--umel_match_csv",
        default="/cbica/home/wenju/Dataset/UKBB_UMelbourne/UKB_UMelbourne_vs_Penn_match_key.csv",
        help="Mapping key from UMelbourne ID to Penn/UPenn participant ID.",
    )
    p.add_argument("--no_umelbourne_f53", action="store_true", help="Disable UMelbourne field-53 import.")
    p.add_argument("--admin_censor_date", default="2022-11-30", help="Administrative censor date, YYYY-MM-DD.")

    p.add_argument("--min_case", type=int, default=20, help="Minimum incident cases required.")
    p.add_argument("--min_noncase", type=int, default=20, help="Minimum noncases/censored participants required.")
    p.add_argument("--penalizer", type=float, default=0.0, help="Initial lifelines Cox penalizer.")
    p.add_argument("--l1_ratio", type=float, default=0.0, help="L1 ratio used with penalizer.")
    p.add_argument("--epv_warning_threshold", type=float, default=10.0, help="Flag models below this events-per-parameter threshold.")

    p.add_argument(
        "--covariate_set",
        choices=["minimal", "clinical"],
        default="clinical",
        help=(
            "minimal = age+sex; clinical = age+sex+smoking+BMI+diastolic+systolic "
            "when available."
        ),
    )
    p.add_argument(
        "--rank_p_threshold",
        type=float,
        default=0.05,
        help="P-value threshold for defining significant individual clocks used in ranking.",
    )
    p.add_argument(
        "--rank_sort_mode",
        choices=["hr_desc", "abs_log_hr"],
        default="hr_desc",
        help="How to sort significant clocks. hr_desc sorts by clock_hr decreasing.",
    )
    p.add_argument(
        "--non_sig_order",
        choices=["random", "input"],
        default="random",
        help="How to append non-significant clocks after significant clocks.",
    )
    p.add_argument("--random_seed", type=int, default=20260718, help="Seed for reproducible non-significant clock order.")
    p.add_argument("--overwrite", action="store_true", help="Overwrite existing nonempty output.")
    return p.parse_args()


def disease_id_from_path(path: str) -> str:
    base = os.path.basename(path)
    suffixes = [
        "_diagnosis_clock_disease_free.tsv",
        "_diagnosis_clock.tsv",
        "_diagnosis.tsv",
        ".tsv",
    ]
    for suffix in suffixes:
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return os.path.splitext(base)[0]


def normalize_participant_id(df: pd.DataFrame, col: str = "participant_id") -> pd.DataFrame:
    if col not in df.columns:
        raise ValueError(f"Missing ID column: {col}")
    out = df.copy()
    out[col] = pd.to_numeric(out[col], errors="coerce").astype("Int64")
    out = out[out[col].notna()].copy()
    return out


def clean_event_dates(series: pd.Series) -> pd.Series:
    x = series.copy()
    x = x.replace([0, 0.0, "0", "0.0", "", "NA", "NaN", "nan", "None", "-1", -1], np.nan)
    return pd.to_datetime(x, errors="coerce")


def parse_ukb_date(series: pd.Series) -> pd.Series:
    x = series.copy()
    x = x.replace([0, 0.0, "0", "0.0", "", "NA", "NaN", "nan", "None", "-1", -1], np.nan)
    parsed = pd.to_datetime(x, errors="coerce")
    numeric = pd.to_numeric(x, errors="coerce")
    excel_mask = numeric.between(20000, 60000)
    if excel_mask.any():
        excel_dates = pd.to_datetime(numeric, unit="D", origin="1899-12-30", errors="coerce")
        parsed = parsed.where(~excel_mask, excel_dates)
    return parsed


def coalesce_series(primary: pd.Series, fallback: pd.Series) -> pd.Series:
    out = primary.copy()
    out = out.where(out.notna(), fallback)
    return out


def read_covariates(path: str) -> pd.DataFrame:
    cov_all = pd.read_csv(path, low_memory=False)
    if "eid" not in cov_all.columns:
        raise ValueError(f"Cannot find eid in covariate file: {path}")

    sex_candidates = ["sex_f31_0_0", "genetic_sex_f22001_0_0", "Sex", "sex"]
    sex_col = next((c for c in sex_candidates if c in cov_all.columns), None)

    keep = ["eid"]
    for col in [AGE_RECRUIT_COL, SMOKING_COL, BMI_COL, DIASTOLIC_COL, SYSTOLIC_COL, sex_col]:
        if col is not None and col in cov_all.columns and col not in keep:
            keep.append(col)

    cov = cov_all[keep].copy()
    rename = {
        "eid": "participant_id",
        AGE_RECRUIT_COL: "Age_baseline",
        SMOKING_COL: "Smoking",
        BMI_COL: "BMI",
        DIASTOLIC_COL: "Diastolic",
        SYSTOLIC_COL: "Systolic",
    }
    if sex_col is not None:
        rename[sex_col] = "Sex"
    cov = cov.rename(columns=rename)
    cov = normalize_participant_id(cov, "participant_id")
    for c in cov.columns:
        if c != "participant_id":
            cov[c] = pd.to_numeric(cov[c], errors="coerce")
    return cov


def read_fallback_assessment_dates(path: str) -> pd.DataFrame:
    d = pd.read_csv(path, low_memory=False)
    if "eid" in d.columns:
        d = d.rename(columns={"eid": "participant_id"})
    d = normalize_participant_id(d, "participant_id")
    keep = ["participant_id"]
    for c in [BASELINE_DATE_COL, IMAGING_DATE_COL, DEATH_DATE_COL]:
        if c in d.columns:
            keep.append(c)
    d = d[keep].copy()
    for c in [BASELINE_DATE_COL, IMAGING_DATE_COL, DEATH_DATE_COL]:
        if c not in d.columns:
            d[c] = pd.NaT
        d[c] = parse_ukb_date(d[c])
    return d


def read_umelbourne_assessment_dates(death_xlsx: str, match_csv: str) -> pd.DataFrame:
    if not os.path.exists(death_xlsx):
        raise FileNotFoundError(f"UMelbourne death/date Excel not found: {death_xlsx}")
    if not os.path.exists(match_csv):
        raise FileNotFoundError(f"UMelbourne-to-Penn match key not found: {match_csv}")

    df_ukb_death = pd.read_excel(death_xlsx, engine="openpyxl")
    df_id_match = pd.read_csv(match_csv, low_memory=False)

    required_death_cols = ["eid", UMEL_BASELINE_DATE_RAW]
    missing_death = [c for c in required_death_cols if c not in df_ukb_death.columns]
    if missing_death:
        raise ValueError(f"Missing columns in UMelbourne date Excel: {missing_death}")

    required_match_cols = ["id", "id_upenn"]
    missing_match = [c for c in required_match_cols if c not in df_id_match.columns]
    if missing_match:
        raise ValueError(f"Missing columns in UMelbourne/Penn match key: {missing_match}")

    df_ukb_death = df_ukb_death.rename(columns={"eid": "participant_id_umel"})
    df_id_match = df_id_match.rename(columns={"id": "participant_id_umel", "id_upenn": "participant_id"})

    df_ukb_death = normalize_participant_id(df_ukb_death, "participant_id_umel")
    df_id_match = normalize_participant_id(df_id_match, "participant_id_umel")
    df_id_match = normalize_participant_id(df_id_match, "participant_id")

    raw_cols = ["participant_id_umel", UMEL_BASELINE_DATE_RAW]
    if UMEL_IMAGING_DATE_RAW in df_ukb_death.columns:
        raw_cols.append(UMEL_IMAGING_DATE_RAW)
    if UMEL_DEATH_DATE_RAW in df_ukb_death.columns:
        raw_cols.append(UMEL_DEATH_DATE_RAW)

    merged = df_id_match[["participant_id_umel", "participant_id"]].merge(
        df_ukb_death[raw_cols], on="participant_id_umel", how="inner"
    )
    rename = {
        UMEL_BASELINE_DATE_RAW: BASELINE_DATE_COL,
        UMEL_IMAGING_DATE_RAW: IMAGING_DATE_COL,
        UMEL_DEATH_DATE_RAW: DEATH_DATE_COL,
    }
    merged = merged.rename(columns={k: v for k, v in rename.items() if k in merged.columns})

    keep = ["participant_id", BASELINE_DATE_COL]
    for c in [IMAGING_DATE_COL, DEATH_DATE_COL]:
        if c in merged.columns:
            keep.append(c)
    merged = merged[keep].copy()
    for c in [BASELINE_DATE_COL, IMAGING_DATE_COL, DEATH_DATE_COL]:
        if c not in merged.columns:
            merged[c] = pd.NaT
        merged[c] = parse_ukb_date(merged[c])

    merged = merged.sort_values("participant_id").groupby("participant_id", as_index=False).first()
    return merged


def read_assessment_dates(args: argparse.Namespace) -> pd.DataFrame:
    fallback = read_fallback_assessment_dates(args.date_tsv)
    if args.no_umelbourne_f53:
        dates = fallback.copy()
    else:
        umel = read_umelbourne_assessment_dates(args.umel_death_xlsx, args.umel_match_csv)
        dates = umel.merge(fallback, on="participant_id", how="outer", suffixes=("_umel", "_fallback"))
        for c in [BASELINE_DATE_COL, IMAGING_DATE_COL, DEATH_DATE_COL]:
            umel_c = f"{c}_umel"
            fallback_c = f"{c}_fallback"
            if umel_c in dates.columns and fallback_c in dates.columns:
                dates[c] = coalesce_series(dates[umel_c], dates[fallback_c])
            elif umel_c in dates.columns:
                dates[c] = dates[umel_c]
            elif fallback_c in dates.columns:
                dates[c] = dates[fallback_c]
            else:
                dates[c] = pd.NaT
        dates = dates[["participant_id", BASELINE_DATE_COL, IMAGING_DATE_COL, DEATH_DATE_COL]].copy()

    for c in [BASELINE_DATE_COL, IMAGING_DATE_COL, DEATH_DATE_COL]:
        if c not in dates.columns:
            dates[c] = pd.NaT
        dates[c] = parse_ukb_date(dates[c])

    print("Assessment-date coverage after coalescing sources:", flush=True)
    print(f"  {BASELINE_DATE_COL}: {int(dates[BASELINE_DATE_COL].notna().sum()):,}", flush=True)
    print(f"  {IMAGING_DATE_COL}: {int(dates[IMAGING_DATE_COL].notna().sum()):,}", flush=True)
    print(f"  {DEATH_DATE_COL}: {int(dates[DEATH_DATE_COL].notna().sum()):,}", flush=True)
    return dates


def covariate_columns(data: pd.DataFrame, args: argparse.Namespace) -> List[str]:
    if args.covariate_set == "minimal":
        candidates = ["Age_baseline", "Sex"]
    else:
        candidates = ["Age_baseline", "Sex", "Smoking", "BMI", "Diastolic", "Systolic"]
    out: List[str] = []
    for c in candidates:
        if c in data.columns:
            vals = pd.to_numeric(data[c], errors="coerce")
            if vals.notna().sum() > 0 and vals.nunique(dropna=True) > 1:
                data[c] = vals
                out.append(c)
    if "Age_baseline" not in out or "Sex" not in out:
        raise ValueError(f"Age_baseline and Sex must be available. Found covariates: {out}")
    return out


def construct_survival_data(args: argparse.Namespace) -> Tuple[pd.DataFrame, str, pd.DataFrame]:
    disease_id = disease_id_from_path(args.icd_tsv)

    df_clock_all = pd.read_csv(args.icd_tsv, sep="\t", low_memory=False)
    required = ["participant_id", "case", "date"] + CLOCK_COLS
    missing = [c for c in required if c not in df_clock_all.columns]
    if missing:
        raise ValueError(f"Missing columns in ICD clock TSV {args.icd_tsv}: {missing}")
    df_clock = df_clock_all[required].copy()
    df_clock = normalize_participant_id(df_clock, "participant_id")
    n_clock = len(df_clock)

    cov = read_covariates(args.cov_tsv)
    dates = read_assessment_dates(args)

    data = df_clock.merge(cov, on="participant_id", how="left")
    n_after_cov = len(data)
    data = data.merge(dates, on="participant_id", how="left")
    n_after_dates = len(data)

    data["case"] = pd.to_numeric(data["case"], errors="coerce").fillna(0).astype(int)
    data["case"] = (data["case"] == 1).astype(int)
    data["event_date"] = clean_event_dates(data["date"])
    data.loc[data["case"] == 0, "event_date"] = pd.NaT

    admin_censor_date = pd.Timestamp(args.admin_censor_date).normalize()
    data["admin_censor_date"] = admin_censor_date
    if DEATH_DATE_COL in data.columns:
        data[DEATH_DATE_COL] = parse_ukb_date(data[DEATH_DATE_COL])
        data["censor_date"] = data[DEATH_DATE_COL].where(
            data[DEATH_DATE_COL].notna() & (data[DEATH_DATE_COL] < admin_censor_date),
            data["admin_censor_date"],
        )
    else:
        data["censor_date"] = data["admin_censor_date"]

    data[TIME_COL] = np.where(
        data["case"] == 1,
        (data["event_date"] - data[BASELINE_DATE_COL]).dt.days,
        (data["censor_date"] - data[BASELINE_DATE_COL]).dt.days,
    )

    data["disease_id"] = disease_id
    covars = covariate_columns(data, args)

    audit_rows = []
    def add_audit(step: str, df: pd.DataFrame) -> None:
        audit_rows.append(
            {
                "disease_id": disease_id,
                "step": step,
                "N": int(len(df)),
                "N_case": int(pd.to_numeric(df.get("case", pd.Series(dtype=float)), errors="coerce").fillna(0).sum()) if "case" in df.columns else np.nan,
                "N_noncase": int((pd.to_numeric(df.get("case", pd.Series(dtype=float)), errors="coerce").fillna(0) == 0).sum()) if "case" in df.columns else np.nan,
            }
        )

    add_audit("clock_input", df_clock)
    add_audit("after_covariate_merge", data.iloc[:n_after_cov])
    add_audit("after_date_merge", data.iloc[:n_after_dates])

    needed = [TIME_COL, EVENT_COL, *covars, *CLOCK_COLS]
    analysis = data[needed + ["participant_id", "disease_id", "event_date", BASELINE_DATE_COL, "censor_date"]].copy()
    for c in [TIME_COL, EVENT_COL, *covars, *CLOCK_COLS]:
        analysis[c] = pd.to_numeric(analysis[c], errors="coerce")
    analysis = analysis.replace([np.inf, -np.inf], np.nan)
    add_audit("before_time_and_complete_case_filter", analysis)

    analysis = analysis[analysis[TIME_COL] > 0].copy()
    add_audit("after_time_gt_0", analysis)

    analysis = analysis.dropna(subset=[TIME_COL, EVENT_COL, *covars, *CLOCK_COLS]).copy()
    add_audit("after_complete_case_all_15_clocks", analysis)

    n_case = int(analysis[EVENT_COL].sum())
    n_noncase = int((analysis[EVENT_COL] == 0).sum())
    if n_case < args.min_case or n_noncase < args.min_noncase:
        raise ValueError(
            f"Insufficient events after complete-case filtering: cases={n_case}, noncases={n_noncase}. "
            f"Require min_case={args.min_case}, min_noncase={args.min_noncase}."
        )

    for clock in CLOCK_COLS:
        vals = pd.to_numeric(analysis[clock], errors="coerce")
        sd = vals.std(ddof=1)
        if not np.isfinite(sd) or sd <= 0:
            raise ValueError(f"{clock} has zero or undefined variance in the complete-case sample.")
        analysis[f"z__{clock}"] = (vals - vals.mean()) / sd

    audit = pd.DataFrame(audit_rows)
    print(
        f"Analysis population for {disease_id}: N={len(analysis):,}, "
        f"cases={n_case:,}, noncases={n_noncase:,}, covariates={','.join(covars)}",
        flush=True,
    )
    return analysis, disease_id, audit


def choose_p_value(row: pd.Series) -> float:
    for col in ["clock_p", "clock_p_value", "p_clock", "lrt_p_clock_vs_base"]:
        if col in row.index:
            try:
                return float(row[col])
            except Exception:
                pass
    return float("nan")


def choose_hr(row: pd.Series) -> float:
    for col in ["clock_hr", "hazard_ratio_per_1SD", "HR_per_1SD"]:
        if col in row.index:
            try:
                return float(row[col])
            except Exception:
                pass
    return float("nan")


def stable_seed(base_seed: int, disease_id: str) -> int:
    digest = hashlib.md5(disease_id.encode("utf-8")).hexdigest()
    disease_int = int(digest[:8], 16)
    return int((base_seed + disease_int) % (2**32 - 1))


def build_clock_ranking(args: argparse.Namespace, disease_id: str) -> pd.DataFrame:
    ind = pd.read_csv(args.individual_result_tsv, sep="\t", low_memory=False)
    rows: List[Dict[str, object]] = []
    for meta in CLOCKS:
        pair_id = meta["pair_id"]
        clock = meta["clock"]
        hit = pd.DataFrame()
        if "pair_id" in ind.columns:
            hit = ind[ind["pair_id"].astype(str) == pair_id]
        if hit.empty and "clock" in ind.columns:
            hit = ind[ind["clock"].astype(str) == clock]
        if hit.empty:
            rows.append(
                {
                    **meta,
                    "individual_status": "missing_in_individual_result",
                    "individual_clock_hr": np.nan,
                    "individual_clock_p": np.nan,
                    "individual_lrt_p_clock_vs_base": np.nan,
                    "rank_significant": False,
                    "rank_sort_score": -np.inf,
                }
            )
            continue
        row = hit.iloc[0]
        status = str(row.get("status", "ok"))
        hr = choose_hr(row)
        p_value = choose_p_value(row)
        lrt_p = float(row.get("lrt_p_clock_vs_base", np.nan)) if "lrt_p_clock_vs_base" in row.index else np.nan
        significant = (status == "ok") and np.isfinite(p_value) and (p_value < args.rank_p_threshold)
        if args.rank_sort_mode == "abs_log_hr":
            sort_score = abs(math.log(hr)) if np.isfinite(hr) and hr > 0 else -np.inf
        else:
            sort_score = hr if np.isfinite(hr) else -np.inf
        rows.append(
            {
                **meta,
                "individual_status": status,
                "individual_clock_hr": hr,
                "individual_clock_p": p_value,
                "individual_lrt_p_clock_vs_base": lrt_p,
                "individual_clock_cindex": float(row.get("clock_cindex", np.nan)) if "clock_cindex" in row.index else np.nan,
                "individual_delta_cindex_clock_minus_base": float(row.get("delta_cindex_clock_minus_base", np.nan)) if "delta_cindex_clock_minus_base" in row.index else np.nan,
                "rank_significant": bool(significant),
                "rank_sort_score": sort_score,
            }
        )

    rank = pd.DataFrame(rows)
    sig = rank[rank["rank_significant"]].copy()
    nonsig = rank[~rank["rank_significant"]].copy()

    sig = sig.sort_values(["rank_sort_score", "individual_clock_p", "pair_id"], ascending=[False, True, True])
    if args.non_sig_order == "random" and len(nonsig) > 1:
        rng = np.random.default_rng(stable_seed(args.random_seed, disease_id))
        nonsig = nonsig.iloc[rng.permutation(len(nonsig))].copy()
    else:
        nonsig = nonsig.sort_values("pair_id")

    ranked = pd.concat([sig, nonsig], axis=0, ignore_index=True)
    ranked.insert(0, "disease_id", disease_id)
    ranked.insert(1, "cumulative_rank", np.arange(1, len(ranked) + 1))
    ranked["rank_p_threshold"] = args.rank_p_threshold
    ranked["rank_sort_mode"] = args.rank_sort_mode
    ranked["non_sig_order"] = args.non_sig_order
    ranked["random_seed_used"] = stable_seed(args.random_seed, disease_id)
    return ranked


def fit_cox(
    data: pd.DataFrame,
    predictors: Sequence[str],
    penalizer: float,
    l1_ratio: float,
) -> Tuple[CoxPHFitter, pd.DataFrame, float, str]:
    cols = [TIME_COL, EVENT_COL, *predictors]
    fit_df = data.loc[:, cols].replace([np.inf, -np.inf], np.nan).dropna().copy()
    for c in cols:
        fit_df[c] = pd.to_numeric(fit_df[c], errors="coerce")
    fit_df = fit_df.dropna()

    if int(fit_df[EVENT_COL].sum()) < 1:
        raise ValueError("No events available for Cox model.")

    last_error = ""
    for pen in [penalizer, max(penalizer, 0.001), max(penalizer, 0.01), max(penalizer, 0.1)]:
        try:
            cph = CoxPHFitter(penalizer=pen, l1_ratio=l1_ratio)
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                cph.fit(fit_df, duration_col=TIME_COL, event_col=EVENT_COL, show_progress=False)
            warning_text = " | ".join(dict.fromkeys(str(w.message) for w in caught))
            return cph, fit_df, pen, warning_text
        except Exception as exc:
            last_error = f"{type(exc).__name__}: {exc}"
    raise RuntimeError(last_error)


def get_cindex(cph: CoxPHFitter, fit_df: pd.DataFrame, predictors: Sequence[str]) -> float:
    risk = np.log(cph.predict_partial_hazard(fit_df.loc[:, list(predictors)]).to_numpy(dtype=float).reshape(-1))
    return float(concordance_index(fit_df[TIME_COL].to_numpy(dtype=float), -risk, fit_df[EVENT_COL].to_numpy(dtype=int)))


def extract_hr(cph: CoxPHFitter, predictor: str) -> Dict[str, float]:
    if predictor not in cph.summary.index:
        return {"beta": np.nan, "se": np.nan, "hr": np.nan, "ci_lo": np.nan, "ci_hi": np.nan, "p": np.nan}
    row = cph.summary.loc[predictor]
    beta = float(row.get("coef", np.nan))
    se = float(row.get("se(coef)", np.nan))
    p = float(row.get("p", np.nan))
    return {
        "beta": beta,
        "se": se,
        "hr": float(row.get("exp(coef)", math.exp(beta) if np.isfinite(beta) else np.nan)),
        "ci_lo": float(row.get("exp(coef) lower 95%", np.nan)),
        "ci_hi": float(row.get("exp(coef) upper 95%", np.nan)),
        "p": p,
    }


def lrt_pvalue(full: CoxPHFitter, reduced: CoxPHFitter, df: int, full_pen: float, reduced_pen: float) -> Tuple[float, float]:
    if full_pen > 0 or reduced_pen > 0:
        return np.nan, np.nan
    stat = max(0.0, 2.0 * (float(full.log_likelihood_) - float(reduced.log_likelihood_)))
    return stat, float(chi2.sf(stat, df))


def make_stratified_fold_ids(events: pd.Series, folds: int, seed: int) -> np.ndarray:
    """Create deterministic stratified fold assignments for event/non-event status."""
    if folds < 2:
        raise ValueError("CV folds must be at least 2.")

    event_values = pd.to_numeric(events, errors="coerce").fillna(0).astype(int).to_numpy()
    n_events = int(event_values.sum())
    n_nonevents = int(len(event_values) - n_events)
    if n_events < folds or n_nonevents < folds:
        raise ValueError(
            f"Cannot create {folds} stratified folds with "
            f"{n_events} events and {n_nonevents} non-events."
        )

    rng = np.random.default_rng(seed)
    fold_ids = np.empty(len(event_values), dtype=int)
    for class_value in (0, 1):
        idx = np.flatnonzero(event_values == class_value)
        rng.shuffle(idx)
        for fold, fold_idx in enumerate(np.array_split(idx, folds)):
            fold_ids[fold_idx] = fold
    return fold_ids


def add_fold_standardized_clock_columns(
    train: pd.DataFrame,
    test: pd.DataFrame,
    selected_clocks: Sequence[str],
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, str]]:
    """Standardize clock columns using training-fold mean/SD and apply to test."""
    train_cv = train.copy()
    test_cv = test.copy()
    z_by_clock: Dict[str, str] = {}

    for clock in selected_clocks:
        z_col = f"z__{clock}"
        train_values = pd.to_numeric(train_cv[clock], errors="coerce")
        mean = float(train_values.mean())
        sd = float(train_values.std(ddof=1))
        if not np.isfinite(sd) or sd <= 0:
            raise ValueError(f"{clock} has zero or undefined training-fold SD.")

        train_cv[z_col] = (pd.to_numeric(train_cv[clock], errors="coerce") - mean) / sd
        test_cv[z_col] = (pd.to_numeric(test_cv[clock], errors="coerce") - mean) / sd
        z_by_clock[clock] = z_col

    return train_cv, test_cv, z_by_clock


def cross_validated_cindex_fixed_order(
    args: argparse.Namespace,
    data: pd.DataFrame,
    covars: Sequence[str],
    selected_clocks: Sequence[str],
    fold_ids: np.ndarray,
) -> Tuple[float, str, str, str]:
    """
    Compute out-of-fold Harrell C-index for one cumulative model.

    This implements Option 1: the clock order is fixed from the full-sample
    individual analyses, while each cumulative Cox model is fit only in the
    training folds and evaluated in the held-out fold. Clock standardization is
    also learned in the training fold and then applied to the held-out fold.
    """
    if len(fold_ids) != len(data):
        raise ValueError("fold_ids length does not match data length.")

    out_of_fold_log_risk = np.full(len(data), np.nan, dtype=float)
    warnings_by_fold: List[str] = []
    penalties_by_fold: List[str] = []

    for fold in sorted(np.unique(fold_ids)):
        train_idx = np.flatnonzero(fold_ids != fold)
        test_idx = np.flatnonzero(fold_ids == fold)
        train = data.iloc[train_idx].copy()
        test = data.iloc[test_idx].copy()

        train_cv, test_cv, z_by_clock = add_fold_standardized_clock_columns(
            train, test, selected_clocks
        )
        predictors = [*covars, *[z_by_clock[x] for x in selected_clocks]]

        cph, df_train, pen_used, warn = fit_cox(
            train_cv, predictors, args.penalizer, args.l1_ratio
        )
        if warn:
            warnings_by_fold.append(f"fold {fold}: {warn}")
        penalties_by_fold.append(f"fold {fold}: {pen_used}")

        test_model = test_cv.loc[:, [TIME_COL, EVENT_COL, *predictors]].replace(
            [np.inf, -np.inf], np.nan
        ).dropna().copy()
        if len(test_model) != len(test):
            raise ValueError(
                f"Fold {fold} lost rows during CV test prediction after complete-case filtering."
            )

        partial_hazard = cph.predict_partial_hazard(test_model.loc[:, list(predictors)])
        out_of_fold_log_risk[test_idx] = np.log(
            partial_hazard.to_numpy(dtype=float).reshape(-1)
        )

    if np.isnan(out_of_fold_log_risk).any():
        raise ValueError("Cross-validation did not generate predictions for every participant.")

    cv_cindex = concordance_index(
        data[TIME_COL].to_numpy(dtype=float),
        -out_of_fold_log_risk,
        data[EVENT_COL].to_numpy(dtype=int),
    )
    warning_text = " | ".join(dict.fromkeys(warnings_by_fold))
    penalty_text = " | ".join(dict.fromkeys(penalties_by_fold))
    return float(cv_cindex), "ok", warning_text, penalty_text


def build_cv_results_by_step(
    args: argparse.Namespace,
    data: pd.DataFrame,
    disease_id: str,
    covars: Sequence[str],
    rank: pd.DataFrame,
) -> Dict[int, Dict[str, object]]:
    """Compute five-fold out-of-fold C-index for baseline and each cumulative step."""
    if CV_FOLDS < 2:
        print("CV_FOLDS < 2; skipping cross-validated C-index.", flush=True)
        return {}

    cv_seed = stable_seed(args.random_seed + CV_SEED_OFFSET, disease_id)
    fold_ids = make_stratified_fold_ids(data[EVENT_COL], CV_FOLDS, cv_seed)
    print(
        f"Running Option 1 fixed-order {CV_FOLDS}-fold CV "
        f"with fold seed {cv_seed}...",
        flush=True,
    )

    results: Dict[int, Dict[str, object]] = {}

    # Step 0: covariate-only baseline model.
    try:
        cv_c, cv_status, cv_warning, cv_penalties = cross_validated_cindex_fixed_order(
            args=args,
            data=data,
            covars=covars,
            selected_clocks=[],
            fold_ids=fold_ids,
        )
    except Exception as exc:
        cv_c = np.nan
        cv_status = "cv_failed"
        cv_warning = ""
        cv_penalties = ""
        cv_error = f"{type(exc).__name__}: {exc}"
    else:
        cv_error = ""

    results[0] = {
        "cv_c_index": cv_c,
        "cv_status": cv_status,
        "cv_warning": cv_warning,
        "cv_penalizer_used_by_fold": cv_penalties,
        "cv_error": cv_error,
    }

    selected: List[str] = []
    for _, rank_row in rank.iterrows():
        step = int(rank_row["cumulative_rank"])
        selected.append(str(rank_row["clock"]))
        print(
            f"CV cumulative model {step}/{len(rank)}: {','.join(selected)}",
            flush=True,
        )
        try:
            cv_c, cv_status, cv_warning, cv_penalties = cross_validated_cindex_fixed_order(
                args=args,
                data=data,
                covars=covars,
                selected_clocks=selected,
                fold_ids=fold_ids,
            )
        except Exception as exc:
            cv_c = np.nan
            cv_status = "cv_failed"
            cv_warning = ""
            cv_penalties = ""
            cv_error = f"{type(exc).__name__}: {exc}"
        else:
            cv_error = ""

        results[step] = {
            "cv_c_index": cv_c,
            "cv_status": cv_status,
            "cv_warning": cv_warning,
            "cv_penalizer_used_by_fold": cv_penalties,
            "cv_error": cv_error,
        }

    base_cv = results.get(0, {}).get("cv_c_index", np.nan)
    previous_cv = base_cv
    for step in sorted(results):
        current_cv = results[step].get("cv_c_index", np.nan)
        results[step]["cv_folds"] = CV_FOLDS
        results[step]["cv_seed"] = cv_seed
        results[step]["cv_base_c_index"] = base_cv
        results[step]["delta_cv_c_index_vs_base"] = (
            current_cv - base_cv if np.isfinite(current_cv) and np.isfinite(base_cv) else np.nan
        )
        results[step]["delta_cv_c_index_vs_previous"] = (
            current_cv - previous_cv
            if step != 0 and np.isfinite(current_cv) and np.isfinite(previous_cv)
            else np.nan
        )
        if np.isfinite(current_cv):
            previous_cv = current_cv

    return results


def empty_cv_result() -> Dict[str, object]:
    return {
        "cv_folds": CV_FOLDS,
        "cv_seed": np.nan,
        "cv_c_index": np.nan,
        "cv_base_c_index": np.nan,
        "delta_cv_c_index_vs_base": np.nan,
        "delta_cv_c_index_vs_previous": np.nan,
        "cv_status": "not_run",
        "cv_warning": "",
        "cv_penalizer_used_by_fold": "",
        "cv_error": "",
    }


def common_result_columns(
    args: argparse.Namespace,
    disease_id: str,
    data: pd.DataFrame,
    covars: List[str],
) -> Dict[str, object]:
    followup_years = pd.to_numeric(data[TIME_COL], errors="coerce") / 365.25
    event_followup_years = followup_years[data[EVENT_COL] == 1]
    return {
        "disease_id": disease_id,
        "N": int(len(data)),
        "N_case": int(data[EVENT_COL].sum()),
        "N_noncase": int((data[EVENT_COL] == 0).sum()),
        "event_rate": float(data[EVENT_COL].mean()),
        "median_followup_years": float(followup_years.median()),
        "followup_years_min": float(followup_years.min()),
        "followup_years_max": float(followup_years.max()),
        "event_followup_years_min": float(event_followup_years.min()) if len(event_followup_years) else np.nan,
        "event_followup_years_max": float(event_followup_years.max()) if len(event_followup_years) else np.nan,
        "covariates": ",".join(covars),
        "covariate_set": args.covariate_set,
        "admin_censor_date": args.admin_censor_date,
        "rank_p_threshold": args.rank_p_threshold,
        "rank_sort_mode": args.rank_sort_mode,
        "non_sig_order": args.non_sig_order,
        "analysis_clock_set": "11_proteomics_plus_4_metabolomics",
        "complete_case_mode": "all_15_clocks_for_all_steps",
    }


def run_cumulative(args: argparse.Namespace, data: pd.DataFrame, disease_id: str, rank: pd.DataFrame) -> pd.DataFrame:
    covars = covariate_columns(data, args)
    z_by_clock = {clock: f"z__{clock}" for clock in CLOCK_COLS}
    base_predictors = covars

    # Option 1 cross-validation is computed once for the fixed full-sample ranking.
    # These values are appended as additional columns while keeping the original
    # apparent/in-sample output columns unchanged.
    try:
        cv_by_step = build_cv_results_by_step(args, data, disease_id, covars, rank)
    except Exception as exc:
        print(
            f"WARNING: fixed-order cross-validation failed for {disease_id}: "
            f"{type(exc).__name__}: {exc}",
            file=sys.stderr,
            flush=True,
        )
        cv_by_step = {}

    print("Fitting base Cox model...", flush=True)
    cph_base, df_base, pen_base, warn_base = fit_cox(data, base_predictors, args.penalizer, args.l1_ratio)
    base_cindex = get_cindex(cph_base, df_base, base_predictors)
    base_ll = float(cph_base.log_likelihood_)
    base_aic = float(cph_base.AIC_partial_)

    rows: List[Dict[str, object]] = []
    common = common_result_columns(args, disease_id, data, covars)
    cv0 = {**empty_cv_result(), **cv_by_step.get(0, {})}
    rows.append(
        {
            **common,
            "cumulative_step": 0,
            "added_pair_id": "BASE",
            "added_clock": "BASE",
            "added_organ": "BASE",
            "added_modality": "BASE",
            "clocks_in_model": "",
            "n_clocks": 0,
            "n_model_parameters": len(base_predictors),
            "events_per_parameter": common["N_case"] / max(1, len(base_predictors)),
            "low_EPV_flag": common["N_case"] / max(1, len(base_predictors)) < args.epv_warning_threshold,
            "rank_significant": np.nan,
            "individual_clock_hr": np.nan,
            "individual_clock_p": np.nan,
            "individual_lrt_p_clock_vs_base": np.nan,
            "added_clock_hr_in_cumulative_model": np.nan,
            "added_clock_ci_lo": np.nan,
            "added_clock_ci_hi": np.nan,
            "added_clock_p": np.nan,
            "c_index": base_cindex,
            "base_c_index": base_cindex,
            "delta_c_index_vs_base": 0.0,
            "delta_c_index_vs_previous": np.nan,
            **cv0,
            "log_likelihood": base_ll,
            "partial_AIC": base_aic,
            "lr_chi2_vs_base": 0.0,
            "lr_p_vs_base": np.nan,
            "sequential_lr_chi2_vs_previous": np.nan,
            "sequential_lr_p_vs_previous": np.nan,
            "penalizer_used": pen_base,
            "fit_warning": warn_base,
            "status": "ok",
            "error": "",
        }
    )

    selected: List[str] = []
    previous_cindex = base_cindex
    previous_ll = base_ll
    previous_cph = cph_base
    previous_pen = pen_base

    for _, rank_row in rank.iterrows():
        step = int(rank_row["cumulative_rank"])
        clock = str(rank_row["clock"])
        pair_id = str(rank_row["pair_id"])
        z_col = z_by_clock[clock]
        selected.append(clock)
        z_predictors = [z_by_clock[x] for x in selected]
        predictors = [*base_predictors, *z_predictors]
        parameter_count = len(predictors)
        epv = common["N_case"] / max(1, parameter_count)
        print(f"Fitting cumulative Cox model {step}/{len(rank)}: {','.join(selected)}", flush=True)

        try:
            cph, df_fit, pen_used, warn = fit_cox(data, predictors, args.penalizer, args.l1_ratio)
            cindex = get_cindex(cph, df_fit, predictors)
            ll = float(cph.log_likelihood_)
            aic = float(cph.AIC_partial_)
            hr_stats = extract_hr(cph, z_col)
            lr_stat_base, lr_p_base = lrt_pvalue(cph, cph_base, step, pen_used, pen_base)
            lr_stat_prev, lr_p_prev = lrt_pvalue(cph, previous_cph, 1, pen_used, previous_pen)
            status = "ok"
            error = ""
        except Exception as exc:
            cindex = np.nan
            ll = np.nan
            aic = np.nan
            hr_stats = {"hr": np.nan, "ci_lo": np.nan, "ci_hi": np.nan, "p": np.nan}
            lr_stat_base = lr_p_base = lr_stat_prev = lr_p_prev = np.nan
            pen_used = np.nan
            warn = ""
            status = "cox_fit_failed"
            error = f"{type(exc).__name__}: {exc}"

        cv_step = {**empty_cv_result(), **cv_by_step.get(step, {})}
        rows.append(
            {
                **common,
                "cumulative_step": step,
                "added_pair_id": pair_id,
                "added_clock": clock,
                "added_organ": rank_row.get("organ", ""),
                "added_modality": rank_row.get("modality", ""),
                "clocks_in_model": ",".join(selected),
                "n_clocks": len(selected),
                "n_model_parameters": parameter_count,
                "events_per_parameter": epv,
                "low_EPV_flag": epv < args.epv_warning_threshold,
                "rank_significant": bool(rank_row.get("rank_significant", False)),
                "individual_clock_hr": float(rank_row.get("individual_clock_hr", np.nan)),
                "individual_clock_p": float(rank_row.get("individual_clock_p", np.nan)),
                "individual_lrt_p_clock_vs_base": float(rank_row.get("individual_lrt_p_clock_vs_base", np.nan)),
                "added_clock_hr_in_cumulative_model": hr_stats["hr"],
                "added_clock_ci_lo": hr_stats["ci_lo"],
                "added_clock_ci_hi": hr_stats["ci_hi"],
                "added_clock_p": hr_stats["p"],
                "c_index": cindex,
                "base_c_index": base_cindex,
                "delta_c_index_vs_base": cindex - base_cindex if np.isfinite(cindex) else np.nan,
                "delta_c_index_vs_previous": cindex - previous_cindex if np.isfinite(cindex) and np.isfinite(previous_cindex) else np.nan,
                **cv_step,
                "log_likelihood": ll,
                "partial_AIC": aic,
                "lr_chi2_vs_base": lr_stat_base,
                "lr_p_vs_base": lr_p_base,
                "sequential_lr_chi2_vs_previous": lr_stat_prev,
                "sequential_lr_p_vs_previous": lr_p_prev,
                "penalizer_used": pen_used,
                "fit_warning": warn,
                "status": status,
                "error": error,
            }
        )
        if status != "ok":
            print(f"Stopping after failed cumulative model at step {step}: {error}", file=sys.stderr, flush=True)
            break
        previous_cindex = cindex
        previous_ll = ll
        previous_cph = cph
        previous_pen = pen_used

    return pd.DataFrame(rows)

def write_tsv(df: pd.DataFrame, path: Optional[str]) -> None:
    if not path:
        return
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(f".{out.name}.{os.getpid()}.tmp")
    df.to_csv(tmp, sep="\t", index=False, na_rep="NA")
    os.replace(tmp, out)


def main() -> int:
    args = parse_args()
    out = Path(args.output_tsv)
    if out.exists() and out.stat().st_size > 0 and not args.overwrite:
        print(f"Output exists; skipping. Use --overwrite to rerun: {out}", flush=True)
        return 0

    disease_id = disease_id_from_path(args.icd_tsv)
    print(f"Disease endpoint: {disease_id}", flush=True)
    print(f"ICD TSV: {args.icd_tsv}", flush=True)
    print(f"Individual result TSV: {args.individual_result_tsv}", flush=True)

    data, disease_id, audit = construct_survival_data(args)
    rank = build_clock_ranking(args, disease_id)
    print(
        "Clock ranking: " + " > ".join(rank["pair_id"].astype(str).tolist()),
        flush=True,
    )

    res = run_cumulative(args, data, disease_id, rank)
    write_tsv(res, args.output_tsv)
    write_tsv(rank, args.rank_output_tsv)
    write_tsv(audit, args.audit_tsv)

    print(f"Wrote cumulative results: {args.output_tsv}", flush=True)
    if args.rank_output_tsv:
        print(f"Wrote rank order: {args.rank_output_tsv}", flush=True)
    if args.audit_tsv:
        print(f"Wrote audit table: {args.audit_tsv}", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {type(error).__name__}: {error}", file=sys.stderr, flush=True)
        raise
