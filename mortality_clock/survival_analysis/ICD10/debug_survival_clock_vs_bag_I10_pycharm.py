#!/usr/bin/env python3
"""
PyCharm-friendly debugger for the mortality-clock versus BAG survival comparison.

Purpose
-------
This script reproduces the inputs that the SLURM/bash pipeline would use for one
ICD endpoint, with I10 as the default example. It is designed for interactive
line-by-line debugging in PyCharm and for identifying exactly where sample size
is lost for each mortality-clock/BAG pair.

Default paths are taken from the submitted bash scripts:
  - included ICD list
  - SA data directory
  - output directory
  - default BAG/covariate/date files used by survival_analysis_clock_vs_bag.py

Typical PyCharm use
-------------------
1. Open this file in PyCharm.
2. Set breakpoints inside main(), load_inputs(), audit_pair(), and optionally
   analyze_pair().
3. Run without arguments to debug I10.
4. Inspect the printed sample-size audit and the TSV audit file.

Command-line examples
---------------------
python debug_survival_clock_vs_bag_I10_pycharm.py
python debug_survival_clock_vs_bag_I10_pycharm.py --icd I10 --run_cox
python debug_survival_clock_vs_bag_I10_pycharm.py --icd I10 --pair_id Endocrine_metabolomics
"""

import argparse
import os
import warnings
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from scipy.stats import chi2, norm, wilcoxon

warnings.filterwarnings("ignore")
pd.set_option("display.max_columns", 200)
pd.set_option("display.width", 220)


# -----------------------------------------------------------------------------
# Paths copied from the bash scripts, with I10 as a PyCharm-friendly default.
# -----------------------------------------------------------------------------
DEFAULT_ICD = "I10"
DEFAULT_ICD_LIST = "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/SA/data/included_ICD_mortality_clock.tsv"
DEFAULT_SA_DATA_DIR = "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/SA/data"
DEFAULT_OUT_DIR = "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/SA/output_clock_vs_BAG"

DEFAULT_BAG_TSV = "/cbica/home/wenju/Reproducibile_paper/SleepAging/data/MomoBAG.tsv"
DEFAULT_COV_TSV = "/cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv"
DEFAULT_DATE_TSV = "/cbica/home/wenju/Reproducibile_paper/Multiorgan_Subtype/data/PWAS/UKBB_fullsample_death_variables.csv"

GLOBAL_END_DATE = pd.Timestamp("2022-11-30")


CLOCK_BAG_PAIRS: List[Dict[str, str]] = [
    # MRI clocks, measured at imaging visit
    {"pair_id": "Brain_mri", "organ": "Brain", "modality": "MRI", "time_origin": "imaging", "clock": "Brain_mri", "bag": "Brain_MRIBAG"},
    {"pair_id": "Adipose_mri", "organ": "Adipose", "modality": "MRI", "time_origin": "imaging", "clock": "Adipose_mri", "bag": "Adipose_MRIBAG"},
    {"pair_id": "Heart_mri", "organ": "Heart", "modality": "MRI", "time_origin": "imaging", "clock": "Heart_mri", "bag": "Heart_MRIBAG"},
    {"pair_id": "Kidney_mri", "organ": "Kidney", "modality": "MRI", "time_origin": "imaging", "clock": "Kidney_mri", "bag": "Kidney_MRIBAG"},
    {"pair_id": "Liver_mri", "organ": "Liver", "modality": "MRI", "time_origin": "imaging", "clock": "Liver_mri", "bag": "Liver_MRIBAG"},
    {"pair_id": "Pancreas_mri", "organ": "Pancreas", "modality": "MRI", "time_origin": "imaging", "clock": "Pancreas_mri", "bag": "Pancreas_MRIBAG"},
    {"pair_id": "Spleen_mri", "organ": "Spleen", "modality": "MRI", "time_origin": "imaging", "clock": "Spleen_mri", "bag": "Spleen_MRIBAG"},

    # Proteomic clocks, assumed measured at baseline/blood assessment
    {"pair_id": "Reproductive_female_proteomics", "organ": "Reproductive_female", "modality": "Proteomics", "time_origin": "baseline", "clock": "Reproductive_female_proteomics", "bag": "Reproductive_female_ProtBAG"},
    {"pair_id": "Pulmonary_proteomics", "organ": "Pulmonary", "modality": "Proteomics", "time_origin": "baseline", "clock": "Pulmonary_proteomics", "bag": "Pulmonary_ProtBAG"},
    {"pair_id": "Heart_proteomics", "organ": "Heart", "modality": "Proteomics", "time_origin": "baseline", "clock": "Heart_proteomics", "bag": "Heart_ProtBAG"},
    {"pair_id": "Brain_proteomics", "organ": "Brain", "modality": "Proteomics", "time_origin": "baseline", "clock": "Brain_proteomics", "bag": "Brain_ProtBAG"},
    {"pair_id": "Eye_proteomics", "organ": "Eye", "modality": "Proteomics", "time_origin": "baseline", "clock": "Eye_proteomics", "bag": "Eye_ProtBAG"},
    {"pair_id": "Hepatic_proteomics", "organ": "Hepatic", "modality": "Proteomics", "time_origin": "baseline", "clock": "Hepatic_proteomics", "bag": "Hepatic_ProtBAG"},
    {"pair_id": "Renal_proteomics", "organ": "Renal", "modality": "Proteomics", "time_origin": "baseline", "clock": "Renal_proteomics", "bag": "Renal_ProtBAG"},
    {"pair_id": "Reproductive_male_proteomics", "organ": "Reproductive_male", "modality": "Proteomics", "time_origin": "baseline", "clock": "Reproductive_male_proteomics", "bag": "Reproductive_male_ProtBAG"},
    {"pair_id": "Endocrine_proteomics", "organ": "Endocrine", "modality": "Proteomics", "time_origin": "baseline", "clock": "Endocrine_proteomics", "bag": "Endocrine_ProtBAG"},
    {"pair_id": "Immune_proteomics", "organ": "Immune", "modality": "Proteomics", "time_origin": "baseline", "clock": "Immune_proteomics", "bag": "Immune_ProtBAG"},
    {"pair_id": "Skin_proteomics", "organ": "Skin", "modality": "Proteomics", "time_origin": "baseline", "clock": "Skin_proteomics", "bag": "Skin_ProtBAG"},

    # Metabolomic clocks, assumed measured at baseline/blood assessment
    {"pair_id": "Endocrine_metabolomics", "organ": "Endocrine", "modality": "Metabolomics", "time_origin": "baseline", "clock": "Endocrine_metabolomics", "bag": "Endocrine_MetBAG"},
    {"pair_id": "Digestive_metabolomics", "organ": "Digestive", "modality": "Metabolomics", "time_origin": "baseline", "clock": "Digestive_metabolomics", "bag": "Digestive_MetBAG"},
    {"pair_id": "Hepatic_metabolomics", "organ": "Hepatic", "modality": "Metabolomics", "time_origin": "baseline", "clock": "Hepatic_metabolomics", "bag": "Hepatic_MetBAG"},
    {"pair_id": "Immune_metabolomics", "organ": "Immune", "modality": "Metabolomics", "time_origin": "baseline", "clock": "Immune_metabolomics", "bag": "Immune_MetBAG"},
]

MORTALITY_CLOCK_COLS = [p["clock"] for p in CLOCK_BAG_PAIRS]
BAG_COLS = [p["bag"] for p in CLOCK_BAG_PAIRS]

BASELINE_DATE_COL = "date_of_attending_assessment_centre_f53_0_0"
IMAGING_DATE_COL = "date_of_attending_assessment_centre_f53_2_0"
AGE_RECRUIT_COL = "age_at_recruitment_f21022_0_0"
SMOKING_COL = "smoking_status_f20116_0_0"
BMI_COL = "body_mass_index_bmi_f23104_0_0"


# -----------------------------------------------------------------------------
# Argument parsing with direct PyCharm defaults.
# -----------------------------------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser(
        description="Debug I10 clock-vs-BAG survival-analysis sample-size losses."
    )
    parser.add_argument("--icd", default=DEFAULT_ICD, type=str, help="ICD endpoint to debug; default: I10")
    parser.add_argument("--icd_list", default=DEFAULT_ICD_LIST, type=str)
    parser.add_argument("--sa_data_dir", default=DEFAULT_SA_DATA_DIR, type=str)
    parser.add_argument("--out_dir", default=DEFAULT_OUT_DIR, type=str)
    parser.add_argument("--bag_tsv", default=DEFAULT_BAG_TSV, type=str)
    parser.add_argument("--cov_tsv", default=DEFAULT_COV_TSV, type=str)
    parser.add_argument("--date_tsv", default=DEFAULT_DATE_TSV, type=str)
    parser.add_argument("--min_case", default=20, type=int)
    parser.add_argument("--min_noncase", default=20, type=int)
    parser.add_argument("--penalizer", default=0.0, type=float)
    parser.add_argument(
        "--pair_id",
        default="ALL",
        type=str,
        help="Use ALL or one pair_id, e.g. Brain_mri, Renal_proteomics, Endocrine_metabolomics.",
    )
    parser.add_argument(
        "--run_cox",
        action="store_true",
        help="If set, also run the Cox comparison after the sample-size audit.",
    )
    parser.add_argument(
        "--max_print_rows",
        default=30,
        type=int,
        help="Number of audit rows printed to console.",
    )
    return parser.parse_args()


# -----------------------------------------------------------------------------
# Utility functions.
# -----------------------------------------------------------------------------
def build_paths(args) -> Tuple[str, str, str]:
    icd_tsv = os.path.join(args.sa_data_dir, f"{args.icd}_diagnosis_clock.tsv")
    output_tsv = os.path.join(args.out_dir, f"cox_compare_clock_vs_BAG_{args.icd}.tsv")
    audit_tsv = os.path.join(args.out_dir, f"debug_audit_clock_vs_BAG_{args.icd}.tsv")
    return icd_tsv, output_tsv, audit_tsv


def require_file(path: str, label: str):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing {label}: {path}")


def normalize_participant_id(df: pd.DataFrame, id_col: str = "participant_id") -> pd.DataFrame:
    """Make participant IDs merge-safe while preserving original numeric values when possible."""
    if id_col not in df.columns:
        raise ValueError(f"Missing ID column: {id_col}")
    out = df.copy()
    out[id_col] = pd.to_numeric(out[id_col], errors="coerce").astype("Int64")
    return out


def clean_event_dates(series: pd.Series) -> pd.Series:
    x = series.copy()
    x = x.replace([0, 0.0, "0", "0.0", "", "NA", "NaN", "nan", "None", "-1", -1], np.nan)
    return pd.to_datetime(x, errors="coerce")


def safe_nunique(series: pd.Series) -> int:
    try:
        return int(series.nunique(dropna=True))
    except Exception:
        return 0


def n_nonmissing(df: pd.DataFrame, col: str) -> int:
    if col not in df.columns:
        return 0
    return int(df[col].notna().sum())


def id_set(df: pd.DataFrame, value_col: Optional[str] = None) -> set:
    if value_col is None:
        sub = df[["participant_id"]].dropna()
    else:
        sub = df.loc[df[value_col].notna(), ["participant_id"]].dropna()
    return set(sub["participant_id"].astype(int).tolist())


def read_covariates(path: str) -> pd.DataFrame:
    cov_all = pd.read_csv(path)
    if "eid" not in cov_all.columns:
        raise ValueError(f"Cannot find eid in covariate file: {path}")

    sex_candidates = ["sex_f31_0_0", "genetic_sex_f22001_0_0", "Sex", "sex"]
    sex_col = next((c for c in sex_candidates if c in cov_all.columns), None)

    keep = ["eid"]
    for col in [AGE_RECRUIT_COL, SMOKING_COL, BMI_COL, sex_col]:
        if col is not None and col in cov_all.columns and col not in keep:
            keep.append(col)

    cov = cov_all[keep].copy().rename(columns={"eid": "participant_id"})
    cov = normalize_participant_id(cov, "participant_id")

    rename = {
        AGE_RECRUIT_COL: "Age_baseline",
        SMOKING_COL: "Smoking",
        BMI_COL: "BMI",
    }
    if sex_col is not None:
        rename[sex_col] = "Sex"
    cov = cov.rename(columns=rename)

    for c in cov.columns:
        if c != "participant_id":
            cov[c] = pd.to_numeric(cov[c], errors="coerce")
    return cov


def read_assessment_dates(path: str) -> pd.DataFrame:
    d = pd.read_csv(path)
    if "eid" in d.columns:
        d = d.rename(columns={"eid": "participant_id"})
    d = normalize_participant_id(d, "participant_id")

    keep = ["participant_id", BASELINE_DATE_COL, IMAGING_DATE_COL]
    missing = [c for c in keep if c not in d.columns]
    if missing:
        raise ValueError(f"Missing date columns in {path}: {missing}")

    d = d[keep].copy()
    d[BASELINE_DATE_COL] = pd.to_datetime(d[BASELINE_DATE_COL], errors="coerce")
    d[IMAGING_DATE_COL] = pd.to_datetime(d[IMAGING_DATE_COL], errors="coerce")
    return d


def read_bag(path: str) -> pd.DataFrame:
    df_bag_all = pd.read_csv(path, sep="\t")
    df_bag_all = normalize_participant_id(df_bag_all, "participant_id")

    # MomoBAG stores the brain MRI BAG as Brain_PhenoBAG. Create the analysis alias.
    if "Brain_MRIBAG" not in df_bag_all.columns and "Brain_PhenoBAG" in df_bag_all.columns:
        df_bag_all = df_bag_all.rename(columns={"Brain_PhenoBAG": "Brain_MRIBAG"})

    required = ["participant_id"] + BAG_COLS
    missing = [c for c in required if c not in df_bag_all.columns]
    if missing:
        raise ValueError(f"Missing columns in BAG TSV {path}: {missing}")
    return df_bag_all[required].copy()


def read_icd_clock(path: str) -> pd.DataFrame:
    df_clock_all = pd.read_csv(path, sep="\t")
    df_clock_all = normalize_participant_id(df_clock_all, "participant_id")

    required = ["participant_id", "case", "date"] + MORTALITY_CLOCK_COLS
    missing = [c for c in required if c not in df_clock_all.columns]
    if missing:
        raise ValueError(f"Missing columns in ICD clock TSV {path}: {missing}")
    return df_clock_all[required].copy()


def base_covariates(df: pd.DataFrame, age_col: str) -> List[str]:
    candidates = [age_col, "Sex", "Smoking", "BMI"]
    covars = []
    for c in candidates:
        if c in df.columns:
            vals = pd.to_numeric(df[c], errors="coerce")
            if vals.notna().sum() > 0 and vals.nunique(dropna=True) > 1:
                df[c] = vals
                covars.append(c)
    return covars


def standardize_inplace(df: pd.DataFrame, src: str, dst: str) -> bool:
    vals = pd.to_numeric(df[src], errors="coerce")
    sd = vals.std()
    if not np.isfinite(sd) or sd == 0:
        return False
    df[dst] = (vals - vals.mean()) / sd
    return True


# -----------------------------------------------------------------------------
# Loading and construction.
# -----------------------------------------------------------------------------
def load_inputs(args) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, str, str, str]:
    icd_tsv, output_tsv, audit_tsv = build_paths(args)

    print("\n=== Input paths derived from the bash scripts ===")
    print(f"ICD endpoint: {args.icd}")
    print(f"ICD TSV:     {icd_tsv}")
    print(f"Output TSV:  {output_tsv}")
    print(f"Audit TSV:   {audit_tsv}")
    print(f"BAG TSV:     {args.bag_tsv}")
    print(f"Cov TSV:     {args.cov_tsv}")
    print(f"Date TSV:    {args.date_tsv}")

    for label, path in [
        ("ICD TSV", icd_tsv),
        ("BAG TSV", args.bag_tsv),
        ("covariate TSV/CSV", args.cov_tsv),
        ("date TSV/CSV", args.date_tsv),
    ]:
        require_file(path, label)

    df_clock = read_icd_clock(icd_tsv)
    df_bag = read_bag(args.bag_tsv)
    cov = read_covariates(args.cov_tsv)
    dates = read_assessment_dates(args.date_tsv)

    print("\n=== File-level row counts ===")
    print(f"ICD clock rows: {len(df_clock):,}; unique IDs: {df_clock['participant_id'].nunique():,}")
    print(f"BAG rows:       {len(df_bag):,}; unique IDs: {df_bag['participant_id'].nunique():,}")
    print(f"Cov rows:       {len(cov):,}; unique IDs: {cov['participant_id'].nunique():,}")
    print(f"Date rows:      {len(dates):,}; unique IDs: {dates['participant_id'].nunique():,}")
    print(f"Date baseline nonmissing {BASELINE_DATE_COL}: {dates[BASELINE_DATE_COL].notna().sum():,}")
    print(f"Date imaging nonmissing  {IMAGING_DATE_COL}: {dates[IMAGING_DATE_COL].notna().sum():,}")
    for c in ["Age_baseline", "Sex", "Smoking", "BMI"]:
        if c in cov.columns:
            print(f"Covariate nonmissing {c}: {cov[c].notna().sum():,}")

    return df_clock, df_bag, cov, dates, icd_tsv, output_tsv, audit_tsv


def construct_survival_data(
    df_clock: pd.DataFrame, df_bag: pd.DataFrame, cov: pd.DataFrame, dates: pd.DataFrame, icd: str
) -> pd.DataFrame:
    data = df_clock.merge(df_bag, on="participant_id", how="left", validate="one_to_one")
    data = data.merge(cov, on="participant_id", how="left", validate="one_to_one")
    data = data.merge(dates, on="participant_id", how="left", validate="one_to_one")

    data["case"] = pd.to_numeric(data["case"], errors="coerce").fillna(0).astype(int)
    data["case"] = (data["case"] == 1).astype(int)

    data["event_date"] = clean_event_dates(data["date"])
    data.loc[data["case"] == 0, "event_date"] = pd.NaT

    data["time_baseline"] = np.where(
        data["case"] == 1,
        (data["event_date"] - data[BASELINE_DATE_COL]).dt.days,
        (GLOBAL_END_DATE - data[BASELINE_DATE_COL]).dt.days,
    )

    data["time_imaging"] = np.where(
        data["case"] == 1,
        (data["event_date"] - data[IMAGING_DATE_COL]).dt.days,
        (GLOBAL_END_DATE - data[IMAGING_DATE_COL]).dt.days,
    )

    # Use recruitment/baseline age plus elapsed time to imaging. This mirrors the production script.
    if "Age_baseline" in data.columns:
        data["Age_imaging"] = data["Age_baseline"] + (
            (data[IMAGING_DATE_COL] - data[BASELINE_DATE_COL]).dt.days / 365.25
        )
    else:
        data["Age_imaging"] = np.nan

    data["disease_id"] = icd
    return data


# -----------------------------------------------------------------------------
# Sample-size audit.
# -----------------------------------------------------------------------------
def audit_pair(
    pair: Dict[str, str],
    df_clock: pd.DataFrame,
    df_bag: pd.DataFrame,
    cov: pd.DataFrame,
    dates: pd.DataFrame,
    data: pd.DataFrame,
) -> Dict[str, object]:
    clock = pair["clock"]
    bag = pair["bag"]
    time_col = "time_imaging" if pair["time_origin"] == "imaging" else "time_baseline"
    age_col = "Age_imaging" if pair["time_origin"] == "imaging" else "Age_baseline"
    date_col = IMAGING_DATE_COL if pair["time_origin"] == "imaging" else BASELINE_DATE_COL

    # ID-set overlaps before full merge.
    clock_ids = id_set(df_clock, clock)
    bag_ids = id_set(df_bag, bag)
    cov_ids = id_set(cov)
    date_ids = id_set(dates, date_col)

    ids_clock_bag = clock_ids & bag_ids
    ids_clock_bag_cov = ids_clock_bag & cov_ids
    ids_clock_bag_date = ids_clock_bag & date_ids
    ids_clock_bag_cov_date = ids_clock_bag & cov_ids & date_ids

    covars_initial = base_covariates(data.copy(), age_col)
    needed_pre = [clock, bag, time_col, "case"] + covars_initial

    # Stage-by-stage loss inside merged data.
    stage = data[["participant_id", clock, bag, time_col, "case", date_col] + [c for c in covars_initial if c in data.columns]].copy()
    for c in [clock, bag, time_col, "case"] + covars_initial:
        if c in stage.columns:
            stage[c] = pd.to_numeric(stage[c], errors="coerce")

    n_clock_nonmissing_after_merge = n_nonmissing(stage, clock)
    n_bag_nonmissing_after_merge = n_nonmissing(stage, bag)
    n_date_nonmissing_after_merge = n_nonmissing(stage, date_col)
    n_time_nonmissing_after_merge = n_nonmissing(stage, time_col)

    n_clock_bag_after_merge = int(stage[[clock, bag]].dropna().shape[0])
    n_clock_bag_time_after_merge = int(stage[[clock, bag, time_col]].dropna().shape[0])

    n_cov_complete = int(stage[[c for c in covars_initial if c in stage.columns]].dropna().shape[0]) if covars_initial else len(stage)
    n_needed_complete_before_time_filter = int(stage[needed_pre].replace([np.inf, -np.inf], np.nan).dropna().shape[0])

    df_complete = stage[needed_pre].replace([np.inf, -np.inf], np.nan).dropna().copy()
    df_incident = df_complete[df_complete[time_col] > 0].copy()

    n_case_after_incident = int(df_incident["case"].sum()) if len(df_incident) else 0
    n_noncase_after_incident = int((df_incident["case"] == 0).sum()) if len(df_incident) else 0

    followup_years = pd.to_numeric(df_incident[time_col], errors="coerce") / 365.25 if len(df_incident) else pd.Series(dtype=float)
    event_followup_years = followup_years[df_incident["case"] == 1] if len(df_incident) else pd.Series(dtype=float)

    # Case/event-date sanity among clock+BAG complete participants.
    tmp_cb = data.loc[data[[clock, bag]].notna().all(axis=1), ["case", "event_date", time_col]].copy()
    n_cases_clock_bag = int(tmp_cb["case"].sum()) if len(tmp_cb) else 0
    n_cases_clock_bag_event_date_nonmissing = int(tmp_cb.loc[tmp_cb["case"] == 1, "event_date"].notna().sum()) if len(tmp_cb) else 0
    n_cases_clock_bag_time_positive = int(((tmp_cb["case"] == 1) & (pd.to_numeric(tmp_cb[time_col], errors="coerce") > 0)).sum()) if len(tmp_cb) else 0

    return {
        "pair_id": pair["pair_id"],
        "organ": pair["organ"],
        "modality": pair["modality"],
        "time_origin": pair["time_origin"],
        "clock": clock,
        "bag": bag,
        "time_col": time_col,
        "date_col": date_col,
        "age_col": age_col,
        "covars_used_initial": ",".join(covars_initial),

        "file_rows_icd": len(df_clock),
        "file_rows_bag": len(df_bag),
        "file_rows_cov": len(cov),
        "file_rows_dates": len(dates),

        "clock_nonmissing_in_icd": len(clock_ids),
        "bag_nonmissing_in_bag": len(bag_ids),
        "date_nonmissing_in_date_file": len(date_ids),
        "ids_clock_and_bag_before_merge": len(ids_clock_bag),
        "ids_clock_bag_cov_before_merge": len(ids_clock_bag_cov),
        "ids_clock_bag_date_before_merge": len(ids_clock_bag_date),
        "ids_clock_bag_cov_date_before_merge": len(ids_clock_bag_cov_date),

        "clock_nonmissing_after_merge": n_clock_nonmissing_after_merge,
        "bag_nonmissing_after_merge": n_bag_nonmissing_after_merge,
        "date_nonmissing_after_merge": n_date_nonmissing_after_merge,
        "time_nonmissing_after_merge": n_time_nonmissing_after_merge,
        "clock_bag_complete_after_merge": n_clock_bag_after_merge,
        "clock_bag_time_complete_after_merge": n_clock_bag_time_after_merge,
        "covariate_complete_after_merge": n_cov_complete,
        "needed_complete_before_time_filter": n_needed_complete_before_time_filter,
        "final_N_after_time_gt0": len(df_incident),
        "final_N_case": n_case_after_incident,
        "final_N_noncase": n_noncase_after_incident,
        "event_rate_final": float(n_case_after_incident / len(df_incident)) if len(df_incident) else np.nan,

        "followup_years_min": float(followup_years.min()) if followup_years.notna().any() else np.nan,
        "followup_years_max": float(followup_years.max()) if followup_years.notna().any() else np.nan,
        "event_followup_years_min": float(event_followup_years.min()) if event_followup_years.notna().any() else np.nan,
        "event_followup_years_max": float(event_followup_years.max()) if event_followup_years.notna().any() else np.nan,

        "cases_among_clock_bag_complete": n_cases_clock_bag,
        "cases_clock_bag_event_date_nonmissing": n_cases_clock_bag_event_date_nonmissing,
        "cases_clock_bag_time_positive": n_cases_clock_bag_time_positive,
    }


def run_audit(args, df_clock, df_bag, cov, dates, data) -> pd.DataFrame:
    pairs = CLOCK_BAG_PAIRS
    if args.pair_id != "ALL":
        pairs = [p for p in CLOCK_BAG_PAIRS if p["pair_id"] == args.pair_id]
        if not pairs:
            known = ", ".join([p["pair_id"] for p in CLOCK_BAG_PAIRS])
            raise ValueError(f"Unknown pair_id={args.pair_id}. Known pair_id values: {known}")

    rows = [audit_pair(pair, df_clock, df_bag, cov, dates, data) for pair in pairs]
    audit = pd.DataFrame(rows)

    print("\n=== Sample-size audit by pair ===")
    cols_to_print = [
        "pair_id", "modality", "time_origin",
        "clock_nonmissing_in_icd", "bag_nonmissing_in_bag",
        "ids_clock_and_bag_before_merge", "ids_clock_bag_cov_date_before_merge",
        "needed_complete_before_time_filter", "final_N_after_time_gt0",
        "final_N_case", "final_N_noncase", "covars_used_initial",
    ]
    print(audit[cols_to_print].head(args.max_print_rows).to_string(index=False))

    return audit


# -----------------------------------------------------------------------------
# Optional Cox analysis copied from the production script for one-file debugging.
# -----------------------------------------------------------------------------
def fit_cox(df: pd.DataFrame, time_col: str, covars: List[str], predictors: List[str], penalizer: float = 0.0):
    cols = [time_col, "case"] + covars + predictors
    fit_df = df[cols].copy()
    for c in cols:
        fit_df[c] = pd.to_numeric(fit_df[c], errors="coerce")
    fit_df = fit_df.replace([np.inf, -np.inf], np.nan).dropna().copy()
    fit_df = fit_df[fit_df[time_col] > 0].copy()

    usable_covars = []
    for c in covars:
        if fit_df[c].nunique(dropna=True) > 1:
            usable_covars.append(c)
    model_cols = [time_col, "case"] + usable_covars + predictors
    fit_df = fit_df[model_cols].copy()

    if fit_df["case"].sum() == 0 or (fit_df["case"] == 0).sum() == 0:
        raise ValueError("No cases or no non-cases in fit_df.")

    last_error = None
    for pen in [penalizer, 0.001, 0.01, 0.1]:
        try:
            cph = CoxPHFitter(penalizer=pen)
            cph.fit(fit_df, duration_col=time_col, event_col="case", show_progress=False)
            return cph, fit_df, usable_covars, pen
        except Exception as e:
            last_error = e
    raise last_error


def extract_hr(cph: CoxPHFitter, var: str) -> Dict[str, float]:
    beta = float(cph.params_.loc[var])
    se = float(cph.standard_errors_.loc[var])
    return {
        "beta": beta,
        "se": se,
        "hr": float(np.exp(beta)),
        "ci_lo": float(np.exp(beta - 1.96 * se)),
        "ci_hi": float(np.exp(beta + 1.96 * se)),
        "p": float(cph.summary.loc[var, "p"]),
    }


def get_cindex(cph: CoxPHFitter, fit_df: pd.DataFrame, time_col: str) -> float:
    risk = cph.predict_partial_hazard(fit_df).values.ravel()
    return float(concordance_index(fit_df[time_col], -risk, fit_df["case"]))


def lrt_pvalue(cph_full: CoxPHFitter, cph_reduced: CoxPHFitter, df_diff: int) -> float:
    stat = 2.0 * (float(cph_full.log_likelihood_) - float(cph_reduced.log_likelihood_))
    if not np.isfinite(stat) or stat < 0:
        return np.nan
    return float(chi2.sf(stat, df_diff))


def empty_result(icd: str, pair: Dict[str, str], status: str, error: str = "") -> Dict[str, object]:
    out = {
        "disease_id": icd,
        "pair_id": pair["pair_id"],
        "organ": pair["organ"],
        "modality": pair["modality"],
        "time_origin": pair["time_origin"],
        "mortality_clock": pair["clock"],
        "bag": pair["bag"],
        "status": status,
        "error": error,
    }
    numeric_cols = [
        "N", "N_case", "N_noncase", "event_rate",
        "followup_years_min", "followup_years_max",
        "event_followup_years_min", "event_followup_years_max",
        "clock_bag_pearson",
        "clock_beta", "clock_se", "clock_hr", "clock_ci_lo", "clock_ci_hi", "clock_p",
        "bag_beta", "bag_se", "bag_hr", "bag_ci_lo", "bag_ci_hi", "bag_p",
        "clock_joint_beta", "clock_joint_se", "clock_joint_hr", "clock_joint_p",
        "bag_joint_beta", "bag_joint_se", "bag_joint_hr", "bag_joint_p",
        "joint_beta_diff_clock_minus_bag", "joint_se_diff", "joint_z_diff", "joint_p_diff",
        "base_cindex", "clock_cindex", "bag_cindex", "both_cindex",
        "delta_cindex_clock_minus_bag", "delta_cindex_clock_minus_base",
        "delta_cindex_bag_minus_base", "delta_cindex_both_minus_base",
        "lrt_p_clock_vs_base", "lrt_p_bag_vs_base", "lrt_p_both_vs_base",
        "penalizer_base", "penalizer_clock", "penalizer_bag", "penalizer_both",
        "disease_level_n_success_pairs", "disease_level_mean_delta_abs_beta_clock_minus_bag",
        "disease_level_wilcoxon_abs_beta_p", "disease_level_mean_delta_cindex_clock_minus_bag",
        "disease_level_wilcoxon_cindex_p",
    ]
    for c in numeric_cols:
        out[c] = np.nan
    return out


def analyze_pair(data: pd.DataFrame, icd: str, pair: Dict[str, str], args) -> Dict[str, object]:
    clock = pair["clock"]
    bag = pair["bag"]
    time_col = "time_imaging" if pair["time_origin"] == "imaging" else "time_baseline"
    age_col = "Age_imaging" if pair["time_origin"] == "imaging" else "Age_baseline"

    if clock not in data.columns or bag not in data.columns:
        return empty_result(icd, pair, "missing_columns", f"Missing {clock} or {bag}")

    covars_initial = base_covariates(data, age_col)
    needed = [clock, bag, time_col, "case"] + covars_initial
    df = data[needed].copy()
    for c in needed:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.replace([np.inf, -np.inf], np.nan).dropna().copy()
    df = df[df[time_col] > 0].copy()

    followup_years = pd.to_numeric(df[time_col], errors="coerce") / 365.25
    event_followup_years = followup_years[df["case"] == 1]
    followup_summary = {
        "followup_years_min": float(followup_years.min()) if followup_years.notna().any() else np.nan,
        "followup_years_max": float(followup_years.max()) if followup_years.notna().any() else np.nan,
        "event_followup_years_min": float(event_followup_years.min()) if event_followup_years.notna().any() else np.nan,
        "event_followup_years_max": float(event_followup_years.max()) if event_followup_years.notna().any() else np.nan,
    }

    n_case = int(df["case"].sum())
    n_noncase = int((df["case"] == 0).sum())
    if n_case < args.min_case or n_noncase < args.min_noncase:
        out = empty_result(icd, pair, "insufficient_events")
        out.update({"N": len(df), "N_case": n_case, "N_noncase": n_noncase})
        out.update(followup_summary)
        return out

    if not standardize_inplace(df, clock, "clock_z"):
        return empty_result(icd, pair, "zero_variance_clock")
    if not standardize_inplace(df, bag, "bag_z"):
        return empty_result(icd, pair, "zero_variance_bag")

    covars = base_covariates(df, age_col)

    try:
        cph_base, df_base, covars_base, pen_base = fit_cox(df, time_col, covars, [], args.penalizer)
        cph_clock, df_clock, _, pen_clock = fit_cox(df, time_col, covars_base, ["clock_z"], args.penalizer)
        cph_bag, df_bag, _, pen_bag = fit_cox(df, time_col, covars_base, ["bag_z"], args.penalizer)
        cph_both, df_both, _, pen_both = fit_cox(df, time_col, covars_base, ["clock_z", "bag_z"], args.penalizer)
    except Exception as e:
        out = empty_result(icd, pair, "cox_fit_failed", str(e))
        out.update({"N": len(df), "N_case": n_case, "N_noncase": n_noncase})
        out.update(followup_summary)
        return out

    clock_stats = extract_hr(cph_clock, "clock_z")
    bag_stats = extract_hr(cph_bag, "bag_z")
    clock_joint = extract_hr(cph_both, "clock_z")
    bag_joint = extract_hr(cph_both, "bag_z")

    try:
        var_mat = cph_both.variance_matrix_
        var_diff = (
            float(var_mat.loc["clock_z", "clock_z"])
            + float(var_mat.loc["bag_z", "bag_z"])
            - 2.0 * float(var_mat.loc["clock_z", "bag_z"])
        )
        se_diff = np.sqrt(var_diff) if var_diff > 0 else np.nan
        beta_diff = clock_joint["beta"] - bag_joint["beta"]
        z_diff = beta_diff / se_diff if np.isfinite(se_diff) and se_diff > 0 else np.nan
        p_diff = 2.0 * norm.sf(abs(z_diff)) if np.isfinite(z_diff) else np.nan
    except Exception:
        beta_diff, se_diff, z_diff, p_diff = np.nan, np.nan, np.nan, np.nan

    base_cindex = get_cindex(cph_base, df_base, time_col)
    clock_cindex = get_cindex(cph_clock, df_clock, time_col)
    bag_cindex = get_cindex(cph_bag, df_bag, time_col)
    both_cindex = get_cindex(cph_both, df_both, time_col)

    out = empty_result(icd, pair, "ok")
    out.update({
        "N": int(len(df)),
        "N_case": n_case,
        "N_noncase": n_noncase,
        "event_rate": float(n_case / len(df)) if len(df) else np.nan,
        **followup_summary,
        "clock_bag_pearson": float(df[["clock_z", "bag_z"]].corr().iloc[0, 1]),
        "clock_beta": clock_stats["beta"],
        "clock_se": clock_stats["se"],
        "clock_hr": clock_stats["hr"],
        "clock_ci_lo": clock_stats["ci_lo"],
        "clock_ci_hi": clock_stats["ci_hi"],
        "clock_p": clock_stats["p"],
        "bag_beta": bag_stats["beta"],
        "bag_se": bag_stats["se"],
        "bag_hr": bag_stats["hr"],
        "bag_ci_lo": bag_stats["ci_lo"],
        "bag_ci_hi": bag_stats["ci_hi"],
        "bag_p": bag_stats["p"],
        "clock_joint_beta": clock_joint["beta"],
        "clock_joint_se": clock_joint["se"],
        "clock_joint_hr": clock_joint["hr"],
        "clock_joint_p": clock_joint["p"],
        "bag_joint_beta": bag_joint["beta"],
        "bag_joint_se": bag_joint["se"],
        "bag_joint_hr": bag_joint["hr"],
        "bag_joint_p": bag_joint["p"],
        "joint_beta_diff_clock_minus_bag": beta_diff,
        "joint_se_diff": se_diff,
        "joint_z_diff": z_diff,
        "joint_p_diff": p_diff,
        "base_cindex": base_cindex,
        "clock_cindex": clock_cindex,
        "bag_cindex": bag_cindex,
        "both_cindex": both_cindex,
        "delta_cindex_clock_minus_bag": clock_cindex - bag_cindex,
        "delta_cindex_clock_minus_base": clock_cindex - base_cindex,
        "delta_cindex_bag_minus_base": bag_cindex - base_cindex,
        "delta_cindex_both_minus_base": both_cindex - base_cindex,
        "lrt_p_clock_vs_base": lrt_pvalue(cph_clock, cph_base, 1),
        "lrt_p_bag_vs_base": lrt_pvalue(cph_bag, cph_base, 1),
        "lrt_p_both_vs_base": lrt_pvalue(cph_both, cph_base, 2),
        "penalizer_base": pen_base,
        "penalizer_clock": pen_clock,
        "penalizer_bag": pen_bag,
        "penalizer_both": pen_both,
    })
    return out


def add_disease_level_tests(res: pd.DataFrame) -> pd.DataFrame:
    ok = res[res["status"] == "ok"].copy()
    n_ok = len(ok)
    mean_delta_abs_beta = np.nan
    p_abs_beta = np.nan
    mean_delta_cindex = np.nan
    p_cindex = np.nan

    if n_ok >= 5:
        delta_abs_beta = ok["clock_beta"].abs().values - ok["bag_beta"].abs().values
        mean_delta_abs_beta = float(np.nanmean(delta_abs_beta))
        try:
            p_abs_beta = float(wilcoxon(delta_abs_beta, zero_method="wilcox", alternative="two-sided").pvalue)
        except Exception:
            p_abs_beta = np.nan

        delta_cindex = ok["delta_cindex_clock_minus_bag"].values
        mean_delta_cindex = float(np.nanmean(delta_cindex))
        try:
            p_cindex = float(wilcoxon(delta_cindex, zero_method="wilcox", alternative="two-sided").pvalue)
        except Exception:
            p_cindex = np.nan

    res["disease_level_n_success_pairs"] = n_ok
    res["disease_level_mean_delta_abs_beta_clock_minus_bag"] = mean_delta_abs_beta
    res["disease_level_wilcoxon_abs_beta_p"] = p_abs_beta
    res["disease_level_mean_delta_cindex_clock_minus_bag"] = mean_delta_cindex
    res["disease_level_wilcoxon_cindex_p"] = p_cindex
    return res


def run_cox_if_requested(args, data: pd.DataFrame, output_tsv: str):
    if not args.run_cox:
        print("\nSkipping Cox analysis because --run_cox was not set.")
        return

    pairs = CLOCK_BAG_PAIRS
    if args.pair_id != "ALL":
        pairs = [p for p in CLOCK_BAG_PAIRS if p["pair_id"] == args.pair_id]

    print("\n=== Running Cox comparison ===")
    rows = [analyze_pair(data, args.icd, pair, args) for pair in pairs]
    res = pd.DataFrame(rows)
    res = add_disease_level_tests(res)
    os.makedirs(os.path.dirname(output_tsv), exist_ok=True)
    res.to_csv(output_tsv, sep="\t", index=False)
    print(f"Saved Cox output: {output_tsv}")
    print(res[["pair_id", "status", "N", "N_case", "N_noncase", "clock_hr", "bag_hr", "delta_cindex_clock_minus_bag"]].to_string(index=False))


# -----------------------------------------------------------------------------
# Main.
# -----------------------------------------------------------------------------
def main():
    args = parse_args()

    df_clock, df_bag, cov, dates, icd_tsv, output_tsv, audit_tsv = load_inputs(args)
    data = construct_survival_data(df_clock, df_bag, cov, dates, args.icd)

    print("\n=== Merged analysis table ===")
    print(f"Rows after merging ICD + BAG + covariates + dates: {len(data):,}")
    print(f"I10 cases in merged data: {int(data['case'].sum()):,}")
    print(f"I10 cases with usable event date: {int(data.loc[data['case'] == 1, 'event_date'].notna().sum()):,}")

    audit = run_audit(args, df_clock, df_bag, cov, dates, data)
    os.makedirs(os.path.dirname(audit_tsv), exist_ok=True)
    audit.to_csv(audit_tsv, sep="\t", index=False)
    print(f"\nSaved audit output: {audit_tsv}")

    print("\n=== How to diagnose the bottleneck ===")
    print("If clock_nonmissing_in_icd is already small, the ICD diagnosis-clock TSV is upstream-restricted.")
    print("If bag_nonmissing_in_bag is small, the BAG TSV is the bottleneck.")
    print("If ids_clock_and_bag_before_merge is large but ids_clock_bag_cov_date_before_merge is small, cov/date merge is the bottleneck.")
    print("If needed_complete_before_time_filter is large but final_N_after_time_gt0 is small, incident-date/time-origin filtering is the bottleneck.")

    run_cox_if_requested(args, data, output_tsv)


if __name__ == "__main__":
    main()
