#!/usr/bin/env python3
"""
Test longitudinal metabolomics mortality-clock delta biomarkers for incident algorithmically-defined diseases.

This script is adapted from the baseline mortality-clock disease-onset survival pipeline, but replaces
baseline mortality clocks as the main predictor with longitudinal delta biomarkers from four metabolomics
mortality clocks:

  Endocrine, Digestive, Hepatic, Immune

Input delta files are expected to be produced by the longitudinal metabolomics delta script:

  <delta_root>/Endocrine/endocrine_wide_delta_acceleration_years.tsv
  <delta_root>/Digestive/digestive_wide_delta_acceleration_years.tsv
  <delta_root>/Hepatic/hepatic_wide_delta_acceleration_years.tsv
  <delta_root>/Immune/immune_wide_delta_acceleration_years.tsv

Primary model for each clock and endpoint:

  Cox(time from instance 1 to incident disease/censoring) ~
      common covariates + baseline mortality-clock acceleration + delta biomarker

The main test is whether the delta biomarker adds predictive power beyond the common covariates and
baseline mortality-clock acceleration, tested using both:
  1) Wald P-value for the delta biomarker in the full model
  2) likelihood-ratio test comparing reduced model vs full model

By default, --delta_column is delta_clock_age_1_minus_0 because the user specifically requested that
column. To use the pure acceleration-year change instead, set:
  --delta_column delta_accel_years_1_minus_0
"""

import argparse
import json
import os
import re
import warnings
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from scipy.stats import chi2

warnings.filterwarnings("ignore")


ORGANS: Dict[str, str] = {
    "Endocrine": "endocrine",
    "Digestive": "digestive",
    "Hepatic": "hepatic",
    "Immune": "immune",
}

# UK Biobank algorithmically-defined outcome date fields under Category 42.
# Asthma was listed twice in the user message, so this script includes the unique endpoints.
ENDPOINTS: Dict[str, Dict[str, str]] = {
    "all_cause_dementia": {
        "label": "All-cause dementia",
        "date_field": "42018",
        "source_field": "42019",
    },
    "asthma": {
        "label": "Asthma",
        "date_field": "42014",
        "source_field": "42015",
    },
    "myocardial_infarction": {
        "label": "Myocardial infarction",
        "date_field": "42000",
        "source_field": "42001",
    },
    "copd": {
        "label": "COPD",
        "date_field": "42016",
        "source_field": "42017",
    },
    "stroke": {
        "label": "Stroke",
        "date_field": "42006",
        "source_field": "42007",
    },
}

DEATH_FIELD = "40000"
ADMIN_CENSOR_DATE_DEFAULT = pd.Timestamp("2022-11-30")

AGE_RECRUIT_COL = "age_at_recruitment_f21022_0_0"
SMOKING_COL = "smoking_status_f20116_0_0"
BMI_COL = "body_mass_index_bmi_f23104_0_0"
DIASTOLIC_COL = "diastolic_blood_pressure_automated_reading_f4079_0_0"
SYSTOLIC_COL = "systolic_blood_pressure_automated_reading_f4080_0_0"


def parse_args():
    p = argparse.ArgumentParser(
        description="Cox survival analysis for metabolomics mortality-clock delta biomarkers and algorithmic disease outcomes."
    )
    p.add_argument("--endpoint", required=True, choices=sorted(ENDPOINTS.keys()))
    p.add_argument("--output_tsv", required=True)
    p.add_argument(
        "--delta_root",
        default=(
            "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/"
            "mortality_clock/longitudinal/metabolomics/"
            "metabolomics_delta_acceleration_years_landmark_survival_analysis"
        ),
        help="Directory containing organ-specific *_wide_delta_acceleration_years.tsv files.",
    )
    p.add_argument(
        "--delta_column",
        default="delta_clock_age_1_minus_0",
        choices=["delta_clock_age_1_minus_0", "delta_accel_years_1_minus_0"],
        help=(
            "Column used as the longitudinal delta biomarker. Default follows the user-provided column. "
            "Use delta_accel_years_1_minus_0 for pure acceleration-year change."
        ),
    )
    p.add_argument(
        "--umel_death_xlsx",
        default="/cbica/home/wenju/Dataset/UKBB_UMelbourne/Death_related_var_from_Ye.xlsx",
        help="UMelbourne UKBB Excel file containing algorithmically-defined disease outcome dates.",
    )
    p.add_argument(
        "--umel_match_csv",
        default="/cbica/home/wenju/Dataset/UKBB_UMelbourne/UKB_UMelbourne_vs_Penn_match_key.csv",
        help="Mapping key from UMelbourne participant ID to Penn/UPenn participant ID.",
    )
    p.add_argument(
        "--cov_tsv",
        default="/cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv",
        help="UKB covariate CSV containing eid and common covariates.",
    )
    p.add_argument("--admin_censor_date", default="2022-11-30")
    p.add_argument("--min_case", default=20, type=int)
    p.add_argument("--min_noncase", default=20, type=int)
    p.add_argument("--penalizer", default=0.0, type=float)
    p.add_argument(
        "--include_bp",
        action="store_true",
        help="Include baseline systolic and diastolic BP as additional common covariates when available.",
    )
    p.add_argument(
        "--no_plot",
        action="store_true",
        help="Disable forest-plot generation.",
    )
    return p.parse_args()


def normalize_participant_id(df: pd.DataFrame, col: str = "participant_id") -> pd.DataFrame:
    if col not in df.columns:
        raise ValueError(f"Missing ID column: {col}")
    out = df.copy()
    out[col] = pd.to_numeric(out[col], errors="coerce").astype("Int64")
    out = out[out[col].notna()].copy()
    out[col] = out[col].astype(int)
    return out


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


def coalesce_columns(df: pd.DataFrame, cols: List[str]) -> pd.Series:
    available = [c for c in cols if c in df.columns]
    if not available:
        return pd.Series(pd.NaT, index=df.index)
    out = df[available[0]].copy()
    for c in available[1:]:
        out = out.where(out.notna(), df[c])
    return out


def find_ukb_field_column(df: pd.DataFrame, field_id: str, required: bool = True) -> str:
    """Find a UKB field column with robust support for names like 42018-0.0 or ...f42018_0_0."""
    field_id = str(field_id)
    col_str = {c: str(c) for c in df.columns}

    exact_candidates = [
        field_id,
        f"{field_id}-0.0",
        f"{field_id}_0_0",
        f"f{field_id}",
        f"f{field_id}_0_0",
    ]
    for wanted in exact_candidates:
        for c, s in col_str.items():
            if s == wanted:
                return c

    patterns = [
        rf"(^|[^0-9]){re.escape(field_id)}([^0-9]|$)",
        rf"f{re.escape(field_id)}(_|$)",
    ]
    matches = []
    for c, s in col_str.items():
        for pat in patterns:
            if re.search(pat, s):
                matches.append(c)
                break

    # Prefer instance-0 array if multiple are present.
    if matches:
        for c in matches:
            s = str(c)
            if s.endswith("-0.0") or s.endswith("_0_0"):
                return c
        return matches[0]

    if required:
        preview = ", ".join(list(map(str, df.columns[:40])))
        raise ValueError(f"Could not find UKB field {field_id}. First columns: {preview}")
    return ""


def read_algorithmic_outcome_dates(args) -> pd.DataFrame:
    endpoint_info = ENDPOINTS[args.endpoint]
    if not os.path.exists(args.umel_death_xlsx):
        raise FileNotFoundError(args.umel_death_xlsx)
    if not os.path.exists(args.umel_match_csv):
        raise FileNotFoundError(args.umel_match_csv)

    raw = pd.read_excel(args.umel_death_xlsx, engine="openpyxl")
    match = pd.read_csv(args.umel_match_csv)

    if "eid" not in raw.columns:
        raise ValueError("UMelbourne Excel file must contain eid.")
    if not all(c in match.columns for c in ["id", "id_upenn"]):
        raise ValueError("ID match CSV must contain id and id_upenn.")

    disease_col = find_ukb_field_column(raw, endpoint_info["date_field"], required=True)
    source_col = find_ukb_field_column(raw, endpoint_info["source_field"], required=False)
    death_col = find_ukb_field_column(raw, DEATH_FIELD, required=False)

    raw = raw.rename(columns={"eid": "participant_id_umel"})
    match = match.rename(columns={"id": "participant_id_umel", "id_upenn": "participant_id"})

    raw = normalize_participant_id(raw, "participant_id_umel")
    match = normalize_participant_id(match, "participant_id_umel")
    match = normalize_participant_id(match, "participant_id")

    keep_raw = ["participant_id_umel", disease_col]
    if source_col:
        keep_raw.append(source_col)
    if death_col:
        keep_raw.append(death_col)

    d = match[["participant_id_umel", "participant_id"]].merge(
        raw[keep_raw], on="participant_id_umel", how="inner"
    )

    out = pd.DataFrame({"participant_id": d["participant_id"].values})
    out["event_date_raw"] = d[disease_col].values
    out["event_date"] = parse_ukb_date(d[disease_col])
    out["event_source"] = d[source_col].values if source_col else np.nan
    out["death_date_algorithmic_file"] = parse_ukb_date(d[death_col]) if death_col else pd.NaT
    out["endpoint"] = args.endpoint
    out["endpoint_label"] = endpoint_info["label"]
    out["endpoint_date_field"] = endpoint_info["date_field"]
    out["endpoint_source_field"] = endpoint_info["source_field"]

    # Collapse duplicated mappings, keeping first non-null values.
    out = (
        out.sort_values("participant_id")
        .groupby("participant_id", as_index=False)
        .first()
    )
    return out


def read_covariates(path: str, include_bp: bool) -> pd.DataFrame:
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    cov_all = pd.read_csv(path)
    if "eid" not in cov_all.columns and "participant_id" not in cov_all.columns:
        raise ValueError("Covariate file must contain eid or participant_id.")
    if "eid" in cov_all.columns:
        cov_all = cov_all.rename(columns={"eid": "participant_id"})

    sex_candidates = ["sex_f31_0_0", "genetic_sex_f22001_0_0", "Sex", "sex"]
    sex_col = next((c for c in sex_candidates if c in cov_all.columns), None)

    keep = ["participant_id"]
    rename = {}

    for src, dst in [
        (AGE_RECRUIT_COL, "Age_recruitment"),
        (SMOKING_COL, "Smoking"),
        (BMI_COL, "BMI"),
        (sex_col, "Sex"),
    ]:
        if src is not None and src in cov_all.columns:
            keep.append(src)
            rename[src] = dst

    if include_bp:
        for src, dst in [(DIASTOLIC_COL, "Diastolic"), (SYSTOLIC_COL, "Systolic")]:
            if src in cov_all.columns:
                keep.append(src)
                rename[src] = dst

    cov = cov_all[keep].copy().rename(columns=rename)
    cov = normalize_participant_id(cov, "participant_id")

    for c in cov.columns:
        if c == "participant_id":
            continue
        cov[c] = pd.to_numeric(cov[c], errors="coerce")

    return cov.drop_duplicates("participant_id", keep="first")


def read_one_delta_file(delta_root: str, organ_label: str, organ_clean: str, delta_column: str) -> pd.DataFrame:
    path = os.path.join(delta_root, organ_label, f"{organ_clean}_wide_delta_acceleration_years.tsv")
    if not os.path.exists(path):
        raise FileNotFoundError(path)

    df = pd.read_csv(path, sep="\t")
    df = normalize_participant_id(df, "participant_id")

    required = [
        "participant_id",
        "sample_date_instance1",
        "clock_acceleration_years_0_0",
        "clock_acceleration_years_1_0",
        "chronological_age_0_0",
        "chronological_age_1_0",
        "delta_chrono_age_1_minus_0",
        delta_column,
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in {path}: {missing}")

    keep = [
        "participant_id",
        "sample_date_instance1",
        "death_date_instance1",
        "admin_censor_date_instance1",
        "end_date_from_instance1",
        "time_from_instance1_years",
        "event_after_instance1",
        "case_status",
        "clock_acceleration_years_0_0",
        "clock_acceleration_years_1_0",
        "chronological_age_0_0",
        "chronological_age_1_0",
        "delta_chrono_age_1_minus_0",
        "delta_accel_years_1_minus_0",
        "delta_clock_age_1_minus_0",
    ]
    keep = [c for c in keep if c in df.columns]
    df = df[keep].copy()

    for c in ["sample_date_instance1", "death_date_instance1", "admin_censor_date_instance1", "end_date_from_instance1"]:
        if c in df.columns:
            df[c] = parse_ukb_date(df[c])

    rename = {
        "clock_acceleration_years_0_0": f"{organ_clean}_baseline_accel_years",
        "clock_acceleration_years_1_0": f"{organ_clean}_followup_accel_years",
        delta_column: f"{organ_clean}_delta_biomarker",
    }
    if "delta_accel_years_1_minus_0" in df.columns:
        rename["delta_accel_years_1_minus_0"] = f"{organ_clean}_delta_accel_years"
    if "delta_clock_age_1_minus_0" in df.columns:
        rename["delta_clock_age_1_minus_0"] = f"{organ_clean}_delta_clock_age_years"

    # Shared columns get organ-specific temporary names, then coalesced later.
    rename.update({
        "sample_date_instance1": f"{organ_clean}_sample_date_instance1",
        "death_date_instance1": f"{organ_clean}_death_date_instance1",
        "admin_censor_date_instance1": f"{organ_clean}_admin_censor_date_instance1",
        "end_date_from_instance1": f"{organ_clean}_end_date_from_instance1",
        "time_from_instance1_years": f"{organ_clean}_mortality_time_from_instance1_years",
        "event_after_instance1": f"{organ_clean}_mortality_event_after_instance1",
        "case_status": f"{organ_clean}_mortality_case_status",
        "chronological_age_0_0": f"{organ_clean}_chronological_age_0_0",
        "chronological_age_1_0": f"{organ_clean}_chronological_age_1_0",
        "delta_chrono_age_1_minus_0": f"{organ_clean}_delta_chrono_age_1_minus_0",
    })

    df = df.rename(columns=rename)
    return df.drop_duplicates("participant_id", keep="first")


def read_delta_biomarkers(delta_root: str, delta_column: str) -> pd.DataFrame:
    merged = None
    for organ_label, organ_clean in ORGANS.items():
        d = read_one_delta_file(delta_root, organ_label, organ_clean, delta_column)
        if merged is None:
            merged = d
        else:
            merged = merged.merge(d, on="participant_id", how="outer")

    assert merged is not None

    sample_cols = [f"{o}_sample_date_instance1" for o in ORGANS.values()]
    admin_cols = [f"{o}_admin_censor_date_instance1" for o in ORGANS.values()]
    death_cols = [f"{o}_death_date_instance1" for o in ORGANS.values()]
    age0_cols = [f"{o}_chronological_age_0_0" for o in ORGANS.values()]
    age1_cols = [f"{o}_chronological_age_1_0" for o in ORGANS.values()]
    delta_age_cols = [f"{o}_delta_chrono_age_1_minus_0" for o in ORGANS.values()]

    merged["sample_date_instance1"] = coalesce_columns(merged, sample_cols)
    merged["admin_censor_date_instance1"] = coalesce_columns(merged, admin_cols)
    merged["death_date_instance1"] = coalesce_columns(merged, death_cols)

    for out_col, cols in [
        ("chronological_age_0_0", age0_cols),
        ("chronological_age_1_0", age1_cols),
        ("delta_chrono_age_1_minus_0", delta_age_cols),
    ]:
        vals = [pd.to_numeric(merged[c], errors="coerce") for c in cols if c in merged.columns]
        if vals:
            merged[out_col] = pd.concat(vals, axis=1).bfill(axis=1).iloc[:, 0]
        else:
            merged[out_col] = np.nan

    return merged


def standardize(values: pd.Series) -> Tuple[pd.Series, float, float]:
    x = pd.to_numeric(values, errors="coerce")
    mean = float(x.mean())
    sd = float(x.std())
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.nan, index=values.index), mean, sd
    return (x - mean) / sd, mean, sd


def fit_cox(df: pd.DataFrame, time_col: str, event_col: str, covars: List[str], predictors: List[str], penalizer: float):
    cols = [time_col, event_col] + covars + predictors
    fit_df = df[cols].copy()
    for c in cols:
        fit_df[c] = pd.to_numeric(fit_df[c], errors="coerce")
    fit_df = fit_df.replace([np.inf, -np.inf], np.nan).dropna().copy()
    fit_df = fit_df[fit_df[time_col] > 0].copy()

    usable_covars = []
    for c in covars:
        if fit_df[c].nunique(dropna=True) > 1:
            usable_covars.append(c)

    model_cols = [time_col, event_col] + usable_covars + predictors
    fit_df = fit_df[model_cols].copy()

    if int(fit_df[event_col].sum()) == 0 or int((fit_df[event_col] == 0).sum()) == 0:
        raise ValueError("No cases or no non-cases in complete-case fit_df.")

    last_error = None
    for pen in [penalizer, 0.001, 0.01, 0.1]:
        try:
            cph = CoxPHFitter(penalizer=pen)
            cph.fit(fit_df, duration_col=time_col, event_col=event_col, show_progress=False)
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


def get_cindex(cph: CoxPHFitter, fit_df: pd.DataFrame, time_col: str, event_col: str) -> float:
    risk = cph.predict_partial_hazard(fit_df).values.ravel()
    return float(concordance_index(fit_df[time_col], -risk, fit_df[event_col]))


def lrt_stats(cph_full: CoxPHFitter, cph_reduced: CoxPHFitter, df_diff: int = 1) -> Tuple[float, float]:
    stat = 2.0 * (float(cph_full.log_likelihood_) - float(cph_reduced.log_likelihood_))
    if not np.isfinite(stat) or stat < 0:
        return np.nan, np.nan
    return stat, float(chi2.sf(stat, df_diff))


def prepare_survival_dataset(args) -> pd.DataFrame:
    admin_censor_date = pd.Timestamp(args.admin_censor_date)

    delta = read_delta_biomarkers(args.delta_root, args.delta_column)
    outcome = read_algorithmic_outcome_dates(args)
    cov = read_covariates(args.cov_tsv, args.include_bp)

    data = delta.merge(outcome, on="participant_id", how="left")
    data = data.merge(cov, on="participant_id", how="left")

    data["sample_date_instance1"] = parse_ukb_date(data["sample_date_instance1"])
    data["admin_censor_date_instance1"] = parse_ukb_date(data["admin_censor_date_instance1"])
    data["death_date_algorithmic_file"] = parse_ukb_date(data["death_date_algorithmic_file"])
    data["event_date"] = parse_ukb_date(data["event_date"])

    data["admin_censor_date"] = admin_censor_date
    data["admin_censor_date"] = data["admin_censor_date_instance1"].where(
        data["admin_censor_date_instance1"].notna(), data["admin_censor_date"]
    )

    data["death_censor_date"] = data["death_date_algorithmic_file"].where(
        data["death_date_algorithmic_file"].notna(), data["death_date_instance1"]
    )

    data["censor_date"] = data["admin_censor_date"]
    death_before_admin = data["death_censor_date"].notna() & (data["death_censor_date"] < data["admin_censor_date"])
    data.loc[death_before_admin, "censor_date"] = data.loc[death_before_admin, "death_censor_date"]

    data["prevalent_or_prior_to_landmark"] = (
        data["event_date"].notna()
        & data["sample_date_instance1"].notna()
        & (data["event_date"] <= data["sample_date_instance1"])
    ).astype(int)

    data["incident_event_after_instance1"] = (
        data["event_date"].notna()
        & data["sample_date_instance1"].notna()
        & (data["event_date"] > data["sample_date_instance1"])
        & (data["event_date"] <= data["censor_date"])
    ).astype(int)

    data["end_date_disease"] = data["censor_date"]
    event_mask = data["incident_event_after_instance1"] == 1
    data.loc[event_mask, "end_date_disease"] = data.loc[event_mask, "event_date"]

    data["time_from_instance1_to_disease_or_censor_days"] = (
        data["end_date_disease"] - data["sample_date_instance1"]
    ).dt.days
    data["time_from_instance1_to_disease_or_censor_years"] = (
        data["time_from_instance1_to_disease_or_censor_days"] / 365.25
    )

    data["disease_case"] = data["incident_event_after_instance1"].astype(int)
    data["endpoint"] = args.endpoint
    data["endpoint_label"] = ENDPOINTS[args.endpoint]["label"]
    data["delta_column_used"] = args.delta_column

    return data


def base_covariates(df: pd.DataFrame, include_bp: bool) -> List[str]:
    # Landmark model: age at instance 1 is preferred because follow-up starts at instance 1.
    candidates = ["chronological_age_1_0", "Sex", "Smoking", "BMI"]
    if include_bp:
        candidates.extend(["Diastolic", "Systolic"])

    covars = []
    for c in candidates:
        if c not in df.columns:
            continue
        vals = pd.to_numeric(df[c], errors="coerce")
        if vals.notna().sum() > 0 and vals.nunique(dropna=True) > 1:
            df[c] = vals
            covars.append(c)
    return covars


def empty_result(args, organ_label: str, organ_clean: str, status: str, error: str = "") -> Dict[str, object]:
    endpoint_info = ENDPOINTS[args.endpoint]
    out = {
        "endpoint": args.endpoint,
        "endpoint_label": endpoint_info["label"],
        "endpoint_date_field": endpoint_info["date_field"],
        "organ_label": organ_label,
        "organ_clean": organ_clean,
        "delta_column_used": args.delta_column,
        "status": status,
        "error": error,
    }
    numeric = [
        "N", "N_case", "N_noncase", "N_prevalent_excluded", "event_rate",
        "followup_years_min", "followup_years_max", "event_followup_years_min", "event_followup_years_max",
        "delta_raw_mean", "delta_raw_sd", "baseline_accel_mean", "baseline_accel_sd",
        "baseline_beta", "baseline_se", "baseline_hr", "baseline_ci_lo", "baseline_ci_hi", "baseline_p",
        "delta_beta", "delta_se", "delta_hr", "delta_ci_lo", "delta_ci_hi", "delta_p",
        "reduced_cindex", "full_cindex", "delta_cindex_full_minus_reduced",
        "lrt_chisq_delta_vs_reduced", "lrt_p_delta_vs_reduced",
        "penalizer_reduced", "penalizer_full",
    ]
    for c in numeric:
        out[c] = np.nan
    return out


def analyze_one_organ(data_all: pd.DataFrame, args, organ_label: str, organ_clean: str) -> Tuple[Dict[str, object], pd.DataFrame]:
    delta_col = f"{organ_clean}_delta_biomarker"
    baseline_col = f"{organ_clean}_baseline_accel_years"

    if delta_col not in data_all.columns or baseline_col not in data_all.columns:
        return empty_result(args, organ_label, organ_clean, "missing_delta_or_baseline_columns"), pd.DataFrame()

    data = data_all.copy()
    data = data[data["prevalent_or_prior_to_landmark"] != 1].copy()
    data = data[data["time_from_instance1_to_disease_or_censor_days"] > 0].copy()

    data["delta_z"], delta_mean, delta_sd = standardize(data[delta_col])
    data["baseline_clock_z"], baseline_mean, baseline_sd = standardize(data[baseline_col])

    if not np.isfinite(delta_sd) or delta_sd == 0:
        return empty_result(args, organ_label, organ_clean, "zero_variance_delta"), pd.DataFrame()
    if not np.isfinite(baseline_sd) or baseline_sd == 0:
        return empty_result(args, organ_label, organ_clean, "zero_variance_baseline_clock"), pd.DataFrame()

    covars = base_covariates(data, args.include_bp)
    time_col = "time_from_instance1_to_disease_or_censor_days"
    event_col = "disease_case"

    needed = [time_col, event_col, delta_col, baseline_col, "delta_z", "baseline_clock_z"] + covars
    fit_pool = data[needed + ["participant_id"]].copy()
    for c in needed:
        fit_pool[c] = pd.to_numeric(fit_pool[c], errors="coerce")
    fit_pool = fit_pool.replace([np.inf, -np.inf], np.nan).dropna().copy()
    fit_pool = fit_pool[fit_pool[time_col] > 0].copy()

    n_case = int(fit_pool[event_col].sum())
    n_noncase = int((fit_pool[event_col] == 0).sum())
    n_prev = int(data_all["prevalent_or_prior_to_landmark"].sum())

    if n_case < args.min_case or n_noncase < args.min_noncase:
        out = empty_result(args, organ_label, organ_clean, "insufficient_events")
        out.update({"N": len(fit_pool), "N_case": n_case, "N_noncase": n_noncase, "N_prevalent_excluded": n_prev})
        return out, fit_pool

    try:
        cph_reduced, df_reduced, covars_used, pen_red = fit_cox(
            fit_pool,
            time_col=time_col,
            event_col=event_col,
            covars=covars,
            predictors=["baseline_clock_z"],
            penalizer=args.penalizer,
        )
        cph_full, df_full, _, pen_full = fit_cox(
            fit_pool,
            time_col=time_col,
            event_col=event_col,
            covars=covars_used,
            predictors=["baseline_clock_z", "delta_z"],
            penalizer=args.penalizer,
        )
    except Exception as e:
        out = empty_result(args, organ_label, organ_clean, "cox_fit_failed", str(e))
        out.update({"N": len(fit_pool), "N_case": n_case, "N_noncase": n_noncase, "N_prevalent_excluded": n_prev})
        return out, fit_pool

    baseline_stats = extract_hr(cph_full, "baseline_clock_z")
    delta_stats = extract_hr(cph_full, "delta_z")
    red_c = get_cindex(cph_reduced, df_reduced, time_col, event_col)
    full_c = get_cindex(cph_full, df_full, time_col, event_col)
    lrt_chisq, lrt_p = lrt_stats(cph_full, cph_reduced, df_diff=1)

    followup_years = pd.to_numeric(fit_pool[time_col], errors="coerce") / 365.25
    event_followup_years = followup_years[fit_pool[event_col] == 1]

    out = empty_result(args, organ_label, organ_clean, "ok")
    out.update({
        "N": int(len(fit_pool)),
        "N_case": n_case,
        "N_noncase": n_noncase,
        "N_prevalent_excluded": n_prev,
        "event_rate": float(n_case / len(fit_pool)) if len(fit_pool) > 0 else np.nan,
        "followup_years_min": float(followup_years.min()) if followup_years.notna().any() else np.nan,
        "followup_years_max": float(followup_years.max()) if followup_years.notna().any() else np.nan,
        "event_followup_years_min": float(event_followup_years.min()) if event_followup_years.notna().any() else np.nan,
        "event_followup_years_max": float(event_followup_years.max()) if event_followup_years.notna().any() else np.nan,
        "delta_raw_mean": delta_mean,
        "delta_raw_sd": delta_sd,
        "baseline_accel_mean": baseline_mean,
        "baseline_accel_sd": baseline_sd,
        "baseline_beta": baseline_stats["beta"],
        "baseline_se": baseline_stats["se"],
        "baseline_hr": baseline_stats["hr"],
        "baseline_ci_lo": baseline_stats["ci_lo"],
        "baseline_ci_hi": baseline_stats["ci_hi"],
        "baseline_p": baseline_stats["p"],
        "delta_beta": delta_stats["beta"],
        "delta_se": delta_stats["se"],
        "delta_hr": delta_stats["hr"],
        "delta_ci_lo": delta_stats["ci_lo"],
        "delta_ci_hi": delta_stats["ci_hi"],
        "delta_p": delta_stats["p"],
        "reduced_cindex": red_c,
        "full_cindex": full_c,
        "delta_cindex_full_minus_reduced": full_c - red_c,
        "lrt_chisq_delta_vs_reduced": lrt_chisq,
        "lrt_p_delta_vs_reduced": lrt_p,
        "penalizer_reduced": pen_red,
        "penalizer_full": pen_full,
        "covariates_used": ",".join(covars_used),
    })

    return out, fit_pool


def fmt_p(p: float) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    if p < 1e-300:
        return "<1e-300"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3f}"


def make_forest_plot(res: pd.DataFrame, args):
    ok = res[(res["status"] == "ok") & res["delta_hr"].notna()].copy()
    if ok.empty:
        return

    try:
        import matplotlib.pyplot as plt
    except Exception as e:
        warnings.warn(f"matplotlib is not available; skipping plot: {e}")
        return

    order = ["Endocrine", "Digestive", "Hepatic", "Immune"]
    ok["organ_order"] = ok["organ_label"].map({x: i for i, x in enumerate(order)})
    ok = ok.sort_values("organ_order", ascending=False).copy()

    y = np.arange(ok.shape[0])
    hr = ok["delta_hr"].astype(float).values
    lo = ok["delta_ci_lo"].astype(float).values
    hi = ok["delta_ci_hi"].astype(float).values
    labels = ok["organ_label"].values

    fig, ax = plt.subplots(figsize=(8.8, 4.4))
    ax.axvline(1.0, linestyle="--", linewidth=1.0)
    ax.errorbar(hr, y, xerr=[hr - lo, hi - hr], fmt="o", capsize=3, linewidth=1.2)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xscale("log")
    ax.set_xlabel("Hazard ratio per 1 SD higher delta biomarker")
    ax.set_title(f"{ENDPOINTS[args.endpoint]['label']}: delta metabolomics mortality-clock biomarkers")

    xmax = max(np.nanmax(hi), 1.05)
    xmin = min(np.nanmin(lo), 0.95)
    ax.set_xlim(max(0.2, xmin * 0.85), xmax * 1.8)

    text_x = xmax * 1.08
    for i, row in enumerate(ok.itertuples(index=False)):
        txt = (
            f"HR {row.delta_hr:.2f} ({row.delta_ci_lo:.2f}-{row.delta_ci_hi:.2f}); "
            f"P={fmt_p(row.delta_p)}; LRT P={fmt_p(row.lrt_p_delta_vs_reduced)}"
        )
        ax.text(text_x, i, txt, va="center", fontsize=8.5)

    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()

    prefix = os.path.splitext(args.output_tsv)[0]
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(f"{prefix}_delta_biomarker_forest_plot.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()
    os.makedirs(os.path.dirname(args.output_tsv), exist_ok=True)

    data = prepare_survival_dataset(args)

    dataset_tsv = os.path.splitext(args.output_tsv)[0] + "_analysis_dataset.tsv"
    data.to_csv(dataset_tsv, sep="\t", index=False)

    rows = []
    complete_case_dir = os.path.splitext(args.output_tsv)[0] + "_complete_cases"
    os.makedirs(complete_case_dir, exist_ok=True)

    for organ_label, organ_clean in ORGANS.items():
        row, fit_pool = analyze_one_organ(data, args, organ_label, organ_clean)
        rows.append(row)
        if fit_pool is not None and not fit_pool.empty:
            fit_pool.to_csv(
                os.path.join(complete_case_dir, f"{args.endpoint}_{organ_clean}_complete_cases.tsv"),
                sep="\t",
                index=False,
            )

    res = pd.DataFrame(rows)
    res.to_csv(args.output_tsv, sep="\t", index=False)

    summary_json = os.path.splitext(args.output_tsv)[0] + "_run_summary.json"
    with open(summary_json, "w") as f:
        json.dump(
            {
                "endpoint": args.endpoint,
                "endpoint_info": ENDPOINTS[args.endpoint],
                "delta_root": args.delta_root,
                "delta_column_used": args.delta_column,
                "output_tsv": args.output_tsv,
                "analysis_dataset_tsv": dataset_tsv,
                "n_rows_analysis_dataset": int(data.shape[0]),
                "n_with_instance1_sample_date": int(data["sample_date_instance1"].notna().sum()),
                "n_prevalent_or_prior_to_landmark": int(data["prevalent_or_prior_to_landmark"].sum()),
                "n_incident_events_after_instance1": int(data["disease_case"].sum()),
            },
            f,
            indent=2,
            default=str,
        )

    if not args.no_plot:
        make_forest_plot(res, args)

    print("Finished delta metabolomics disease-onset survival analysis")
    print(f"Endpoint: {args.endpoint}")
    print(f"Output TSV: {args.output_tsv}")
    print(f"Analysis dataset: {dataset_tsv}")


if __name__ == "__main__":
    main()
