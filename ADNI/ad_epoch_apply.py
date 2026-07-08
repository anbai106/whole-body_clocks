#!/usr/bin/env python3
"""
Apply pre-trained ADNI brain MRI AD L'EPOCH to longitudinal CN MRI scans only.

This script does not refit the model. It loads the saved model bundle from
the ADNI brain MRI AD L'EPOCH training script and applies it to qualified
MUSE MRI scans.

Primary analysis rule:
  Score only CN-labeled MRI scans before MCI/AD conversion.

Time zero:
  first CN visit with usable MUSE GM ROI features;
  prefer Visit_Code == bl if it has usable MUSE features,
  otherwise earliest usable CN MRI visit.

Scored scans:
  selected_baseline
  pre_event_CN_followup
  censored_CN_followup

Excluded scans:
  MCI-labeled scans
  AD-labeled scans
  any scans after first MCI/AD conversion

Primary output:
  adni_brain_mri_ad_lepoch_longitudinal_cn_only_predictions.tsv
"""

import argparse
import json
import re
import warnings
from pathlib import Path

import joblib
import numpy as np
import pandas as pd


RISK_COL = "adni_brain_mri_ad_lepoch_risk_score"
ACCEL_Z_COL = "adni_brain_mri_ad_lepoch_acceleration_z"
ACCEL_YEARS_COL = "adni_brain_mri_ad_lepoch_acceleration_years"
CLOCK_AGE_COL = "adni_brain_mri_ad_lepoch_clock_age_years"


# ============================================================
# Arguments
# ============================================================

def parse_args():
    p = argparse.ArgumentParser()

    p.add_argument("--input-file", required=True)
    p.add_argument("--model-joblib", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--prefix", default="adni_brain_mri_ad_lepoch")

    p.add_argument("--id-col", default="PTID")
    p.add_argument("--visit-col", default="Visit_Code")
    p.add_argument("--date-col", default="Date")
    p.add_argument("--dx-col", default="DX_Binary")

    p.add_argument("--baseline-dx", default="CN")
    p.add_argument("--event-dx", default="MCI,AD")
    p.add_argument("--eligible-scan-dx", default="CN")

    p.add_argument("--min-baseline-roi-fraction", type=float, default=0.80)
    p.add_argument("--min-scan-roi-fraction", type=float, default=0.80)

    p.add_argument(
        "--age-update-mode",
        choices=["from_baseline_date", "row", "none"],
        default="from_baseline_date",
        help=(
            "How to define Age for longitudinal scan scoring. "
            "'from_baseline_date' uses selected baseline Age plus years since baseline; "
            "'row' uses each scan row's Age; 'none' leaves Age unchanged."
        )
    )

    p.add_argument("--risk-times", default="1,2,3,5")

    p.add_argument(
        "--include-selected-baseline",
        action="store_true",
        help="Include the selected baseline scan in the longitudinal output."
    )

    p.add_argument(
        "--complete-case-model-features",
        action="store_true",
        help="Require zero missing model MUSE ROI features per scan."
    )

    return p.parse_args()


# ============================================================
# Helpers
# ============================================================

def log(msg):
    print(msg, flush=True)


def read_table(path):
    path = Path(path)
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path, low_memory=False)
    return pd.read_csv(path, sep="\t", low_memory=False)


def parse_list_arg(x):
    if x is None or str(x).strip() == "":
        return []
    return [v.strip() for v in str(x).split(",") if v.strip()]


def parse_risk_times(x):
    vals = []
    for s in str(x).split(","):
        s = s.strip()
        if s:
            vals.append(float(s))
    return vals


def normalize_dx(x):
    if pd.isna(x):
        return np.nan

    x = str(x).strip().upper()

    dx_map = {
        "CN": "CN",
        "NL": "CN",
        "NORMAL": "CN",
        "CONTROL": "CN",
        "HC": "CN",
        "MCI": "MCI",
        "EMCI": "MCI",
        "LMCI": "MCI",
        "AD": "AD",
        "DEMENTIA": "AD"
    }

    return dx_map.get(x, x)


def parse_date_series(s):
    return pd.to_datetime(s, errors="coerce")


def visit_code_to_month(x):
    if pd.isna(x):
        return np.nan

    x = str(x).strip().lower()

    if x in ["bl", "base", "baseline", "m00", "m0", "screen", "screening", "sc"]:
        return 0.0

    m = re.search(r"m(?:onth)?\s*0*([0-9]+)", x)
    if m:
        return float(m.group(1))

    m = re.search(r"([0-9]+)", x)
    if m:
        return float(m.group(1))

    return np.nan


def clean_numeric_series(s):
    if pd.api.types.is_numeric_dtype(s):
        return pd.to_numeric(s, errors="coerce")

    s2 = (
        s.astype("object")
         .where(s.notna(), np.nan)
         .astype(str)
         .str.strip()
         .str.replace(",", "", regex=False)
    )

    s2 = s2.replace({
        "": np.nan,
        "nan": np.nan,
        "NaN": np.nan,
        "NA": np.nan,
        "N/A": np.nan,
        "None": np.nan,
        "null": np.nan
    })

    return pd.to_numeric(s2, errors="coerce")


def normalize_sex_series(s):
    return (
        s.astype(str)
         .str.strip()
         .replace({
             "0": "Female",
             "0.0": "Female",
             "1": "Male",
             "1.0": "Male",
             "F": "Female",
             "M": "Male",
             "female": "Female",
             "male": "Male",
             "Female": "Female",
             "Male": "Male"
         })
    )


def compute_row_roi_coverage(df, roi_cols):
    roi_numeric = df[roi_cols].apply(clean_numeric_series, axis=0)
    n_nonmissing = roi_numeric.notna().sum(axis=1)
    frac_nonmissing = n_nonmissing / float(len(roi_cols))
    return n_nonmissing, frac_nonmissing


def get_model_roi_features(bundle):
    for key in ["selected_muse_gm_rois", "available_muse_gm_rois", "hardcoded_muse_gm_rois"]:
        if key in bundle and bundle[key]:
            return list(bundle[key])

    numeric_cols = list(bundle.get("numeric_cols", []))
    roi_cols = [c for c in numeric_cols if str(c).startswith("MUSE_Volume_")]

    if len(roi_cols) == 0:
        raise ValueError(
            "Could not identify MUSE ROI features from model bundle. "
            "Expected selected_muse_gm_rois or numeric_cols containing MUSE_Volume_*."
        )

    return roi_cols


def ensure_expected_columns(df, expected_cols):
    df = df.copy()
    for c in expected_cols:
        if c not in df.columns:
            warnings.warn(f"Model expected column {c} is missing in application data. Creating as NA.")
            df[c] = np.nan
    return df


# ============================================================
# Longitudinal CN-only scan selection
# ============================================================

def prepare_longitudinal_cn_table(
    df_raw,
    args,
    model_roi_cols,
    baseline_dx,
    event_dx_set,
    eligible_scan_dx_set
):
    d = df_raw.copy()

    id_col = args.id_col
    visit_col = args.visit_col
    date_col = args.date_col
    dx_col = args.dx_col

    missing_roi = [c for c in model_roi_cols if c not in d.columns]
    present_roi = [c for c in model_roi_cols if c in d.columns]

    if len(present_roi) == 0:
        raise ValueError("None of the model MUSE GM ROI columns are present in the input ADNI file.")

    if missing_roi:
        warnings.warn(
            f"{len(missing_roi)} model MUSE GM ROI columns are missing in the input file. "
            f"First missing: {missing_roi[:10]}"
        )

    d["_dx_norm"] = d[dx_col].apply(normalize_dx)
    d["_date"] = parse_date_series(d[date_col])
    d["_visit_month"] = d[visit_col].apply(visit_code_to_month)

    d["_date_sort"] = d["_date"].map(lambda x: x.toordinal() if pd.notna(x) else np.nan)
    d["_visit_sort"] = d["_visit_month"] * 30.4375
    d["_sort_key"] = d["_date_sort"].fillna(d["_visit_sort"])

    d = d.loc[d[id_col].notna()].copy()
    d = d.loc[d["_sort_key"].notna()].copy()

    d["_muse_gm_roi_nonmissing_n"], d["_muse_gm_roi_nonmissing_fraction"] = compute_row_roi_coverage(
        d,
        present_roi
    )

    d["_has_usable_muse_gm_for_baseline"] = (
        d["_muse_gm_roi_nonmissing_fraction"] >= args.min_baseline_roi_fraction
    )

    d["_has_usable_muse_gm_for_scan"] = (
        d["_muse_gm_roi_nonmissing_fraction"] >= args.min_scan_roi_fraction
    )

    d = d.sort_values([id_col, "_sort_key"], kind="mergesort")

    scored_rows = []
    subject_rows = []
    skipped_subjects = []

    for pid, g in d.groupby(id_col, sort=False):
        g = g.copy().sort_values("_sort_key", kind="mergesort")

        baseline_candidates = g.loc[
            g["_has_usable_muse_gm_for_baseline"] &
            (g["_dx_norm"] == baseline_dx)
        ].copy()

        if baseline_candidates.empty:
            skipped_subjects.append({
                id_col: pid,
                "skip_reason": "no_CN_scan_with_sufficient_MUSE_GM_ROI_coverage",
                "n_rows": int(g.shape[0]),
                "max_muse_gm_roi_fraction": float(g["_muse_gm_roi_nonmissing_fraction"].max())
            })
            continue

        bl_mask = baseline_candidates[visit_col].astype(str).str.lower().isin(
            ["bl", "base", "baseline", "m00", "m0"]
        )

        if bl_mask.any():
            baseline = baseline_candidates.loc[bl_mask].iloc[0]
        else:
            baseline = baseline_candidates.iloc[0]

        baseline_sort = baseline["_sort_key"]
        baseline_date = baseline["_date"]
        baseline_month = baseline["_visit_month"]

        prior_event = g.loc[
            (g["_sort_key"] < baseline_sort) &
            g["_dx_norm"].isin(event_dx_set)
        ]

        if not prior_event.empty:
            skipped_subjects.append({
                id_col: pid,
                "skip_reason": "MCI_or_AD_before_selected_CN_MUSE_baseline",
                "n_rows": int(g.shape[0]),
                "max_muse_gm_roi_fraction": float(g["_muse_gm_roi_nonmissing_fraction"].max())
            })
            continue

        after_baseline = g.loc[g["_sort_key"] > baseline_sort].copy()
        event_rows = after_baseline.loc[after_baseline["_dx_norm"].isin(event_dx_set)].copy()

        if not event_rows.empty:
            event_row = event_rows.iloc[0]
            has_event = True
            event_sort = event_row["_sort_key"]
            event_dx = event_row["_dx_norm"]
            event_date = event_row["_date"]
            event_visit = event_row[visit_col]
        else:
            event_row = None
            has_event = False
            event_sort = np.nan
            event_dx = np.nan
            event_date = pd.NaT
            event_visit = np.nan

        if has_event:
            end_row = event_row
        elif after_baseline.empty:
            end_row = baseline
        else:
            end_row = after_baseline.iloc[-1]

        if pd.notna(baseline_date) and pd.notna(end_row["_date"]):
            time_to_event_or_censor_years = (end_row["_date"] - baseline_date).days / 365.25
        elif pd.notna(baseline_month) and pd.notna(end_row["_visit_month"]):
            time_to_event_or_censor_years = (end_row["_visit_month"] - baseline_month) / 12.0
        else:
            time_to_event_or_censor_years = np.nan

        conversion_group = "Censored / non-event"
        if has_event and event_dx == "MCI":
            conversion_group = "CN to MCI"
        elif has_event and event_dx == "AD":
            conversion_group = "CN to AD"
        elif has_event:
            conversion_group = f"CN to {event_dx}"

        # Key rule:
        #   score only CN scans from selected baseline onward and before first MCI/AD event.
        scan_candidates = g.loc[
            (g["_sort_key"] >= baseline_sort) &
            g["_has_usable_muse_gm_for_scan"] &
            g["_dx_norm"].isin(eligible_scan_dx_set)
        ].copy()

        if not args.include_selected_baseline:
            scan_candidates = scan_candidates.loc[scan_candidates["_sort_key"] > baseline_sort].copy()

        if has_event:
            scan_candidates = scan_candidates.loc[scan_candidates["_sort_key"] < event_sort].copy()

        if scan_candidates.empty:
            skipped_subjects.append({
                id_col: pid,
                "skip_reason": "selected_baseline_exists_but_no_CN_scans_to_score_before_event",
                "n_rows": int(g.shape[0]),
                "max_muse_gm_roi_fraction": float(g["_muse_gm_roi_nonmissing_fraction"].max())
            })
            continue

        baseline_age = np.nan
        if "Age" in baseline.index:
            baseline_age = clean_numeric_series(pd.Series([baseline["Age"]])).iloc[0]

        for _, scan in scan_candidates.iterrows():
            scan_sort = scan["_sort_key"]
            scan_date = scan["_date"]
            scan_month = scan["_visit_month"]

            if pd.notna(baseline_date) and pd.notna(scan_date):
                years_since_baseline = (scan_date - baseline_date).days / 365.25
            elif pd.notna(baseline_month) and pd.notna(scan_month):
                years_since_baseline = (scan_month - baseline_month) / 12.0
            else:
                years_since_baseline = np.nan

            if has_event and pd.notna(event_date) and pd.notna(scan_date):
                years_to_event = (event_date - scan_date).days / 365.25
            elif has_event and pd.notna(event_row["_visit_month"]) and pd.notna(scan_month):
                years_to_event = (event_row["_visit_month"] - scan_month) / 12.0
            else:
                years_to_event = np.nan

            if scan_sort == baseline_sort:
                scan_relation = "selected_baseline"
            elif has_event and scan_sort < event_sort:
                scan_relation = "pre_event_CN_followup"
            else:
                scan_relation = "censored_CN_followup"

            row = scan.copy()

            row["selected_baseline_visit_code"] = baseline[visit_col]
            row["selected_baseline_date"] = baseline_date
            row["selected_baseline_dx"] = baseline_dx
            row["years_since_selected_baseline"] = float(years_since_baseline) if np.isfinite(years_since_baseline) else np.nan

            row["event_from_selected_baseline"] = bool(has_event)
            row["event_dx"] = event_dx
            row["event_visit_code"] = event_visit
            row["event_date"] = event_date
            row["time_to_event_or_censor_years"] = time_to_event_or_censor_years
            row["years_to_event"] = float(years_to_event) if np.isfinite(years_to_event) else np.nan

            row["conversion_group"] = conversion_group
            row["scan_relation_to_event"] = scan_relation
            row["scan_dx"] = scan["_dx_norm"]
            row["scan_date"] = scan_date
            row["scan_visit_month"] = scan_month

            row["n_model_muse_gm_rois_expected"] = len(model_roi_cols)
            row["n_model_muse_gm_rois_present_in_file"] = len(present_roi)
            row["n_model_muse_gm_rois_missing_from_file"] = len(missing_roi)
            row["n_model_muse_gm_rois_nonmissing_in_scan"] = int(scan["_muse_gm_roi_nonmissing_n"])
            row["fraction_model_muse_gm_rois_nonmissing_in_scan"] = float(scan["_muse_gm_roi_nonmissing_fraction"])

            row["baseline_age_raw"] = baseline_age

            scored_rows.append(row)

        subject_rows.append({
            id_col: pid,
            "selected_baseline_visit_code": baseline[visit_col],
            "selected_baseline_date": baseline_date,
            "selected_baseline_dx": baseline_dx,
            "event_from_selected_baseline": bool(has_event),
            "event_dx": event_dx,
            "event_visit_code": event_visit,
            "event_date": event_date,
            "conversion_group": conversion_group,
            "time_to_event_or_censor_years": time_to_event_or_censor_years,
            "n_CN_scans_scored_before_event": int(scan_candidates.shape[0]),
            "baseline_muse_gm_roi_nonmissing_fraction": float(baseline["_muse_gm_roi_nonmissing_fraction"])
        })

    scored_df = pd.DataFrame(scored_rows)
    subject_df = pd.DataFrame(subject_rows)
    skipped_df = pd.DataFrame(skipped_subjects)

    return scored_df, subject_df, skipped_df, present_roi, missing_roi


# ============================================================
# Prediction and acceleration reconstruction
# ============================================================

def categorical_match(series, category_value):
    s_obj = series.astype("object")
    s_str = s_obj.astype(str).str.strip()
    cat_str = str(category_value).strip()

    mask = s_str == cat_str

    try:
        cat_float = float(cat_str)
        s_num = pd.to_numeric(s_obj, errors="coerce")
        mask = mask | np.isclose(s_num.astype(float), cat_float, equal_nan=False)
    except Exception:
        pass

    return mask.astype(float)


def parse_categorical_coef_name(term, categorical_covs):
    if not term.startswith("cat__"):
        return None, None

    stem = term[len("cat__"):]

    for cov in sorted(categorical_covs, key=len, reverse=True):
        prefix = f"{cov}_"
        if stem.startswith(prefix):
            category = stem[len(prefix):]
            return cov, category

    return None, None


def numeric_fill_for_residualization(series, name):
    x = clean_numeric_series(series)
    if x.isna().all():
        warnings.warn(f"Residualization covariate {name} is all missing. Using 0.")
        return pd.Series(np.repeat(0.0, len(x)), index=x.index)

    med = float(np.nanmedian(x))
    return x.fillna(med)


def compute_clock_transforms_from_saved_info(df, risk, clock_info):
    out = df.copy()
    out[RISK_COL] = risk

    if not clock_info:
        warnings.warn("clock_transform_info missing from model bundle. Acceleration columns set to NA.")
        out[ACCEL_Z_COL] = np.nan
        out[ACCEL_YEARS_COL] = np.nan
        out[CLOCK_AGE_COL] = np.nan
        return out

    intercept = float(clock_info.get("risk_score_covariate_model_intercept", 0.0))
    coef_dict = clock_info.get("risk_score_covariate_model_coef", {})

    covariates = list(clock_info.get("residualization_covariates", []))
    numeric_covs = list(clock_info.get("numeric_residualization_covariates", []))
    categorical_covs = list(clock_info.get("categorical_residualization_covariates", []))

    if not numeric_covs and not categorical_covs:
        for c in covariates:
            if c in out.columns and c not in ["Sex", "SITE"]:
                numeric_covs.append(c)
            else:
                categorical_covs.append(c)

    expected = np.repeat(intercept, out.shape[0]).astype(float)

    for term, beta in coef_dict.items():
        beta = float(beta)

        if term.startswith("num__"):
            cov = term[len("num__"):]
            if cov not in out.columns:
                warnings.warn(f"Missing numeric residualization covariate: {cov}. Using 0.")
                vals = np.repeat(0.0, out.shape[0])
            else:
                vals = numeric_fill_for_residualization(out[cov], cov).values.astype(float)

            expected += beta * vals

        elif term.startswith("cat__"):
            cov, category = parse_categorical_coef_name(term, categorical_covs)

            if cov is None:
                warnings.warn(f"Could not parse categorical residualization term: {term}. Ignoring.")
                continue

            if cov not in out.columns:
                warnings.warn(f"Missing categorical residualization covariate: {cov}. Term set to 0.")
                vals = np.repeat(0.0, out.shape[0])
            else:
                vals = categorical_match(out[cov], category).values.astype(float)

            expected += beta * vals

    mean_train = float(clock_info.get("risk_score_residual_mean_train", 0.0))
    sd_train = float(clock_info.get("risk_score_residual_sd_train", np.nan))
    beta_age = clock_info.get("adjusted_age_coefficient_risk_score_per_year", None)

    resid_raw = risk - expected
    resid = resid_raw - mean_train

    if np.isfinite(sd_train) and sd_train > 0:
        out[ACCEL_Z_COL] = resid / sd_train
    else:
        out[ACCEL_Z_COL] = np.nan

    if beta_age is not None and np.isfinite(float(beta_age)) and abs(float(beta_age)) > 1e-8:
        beta_age = float(beta_age)
        out[ACCEL_YEARS_COL] = resid / beta_age

        if "Age" in out.columns:
            out[CLOCK_AGE_COL] = clean_numeric_series(out["Age"]) + out[ACCEL_YEARS_COL]
        else:
            out[CLOCK_AGE_COL] = np.nan
    else:
        out[ACCEL_YEARS_COL] = np.nan
        out[CLOCK_AGE_COL] = np.nan

    return out


def predict_absolute_risk(model, X, times_years):
    out = {}

    try:
        surv_funcs = model.predict_survival_function(X)
    except Exception as e:
        warnings.warn(f"Could not compute absolute risks: {e}")
        for t in times_years:
            out[f"risk_{t:g}y"] = np.repeat(np.nan, X.shape[0])
        return pd.DataFrame(out)

    for t in times_years:
        vals = []
        for sf in surv_funcs:
            try:
                vals.append(1.0 - float(sf(t)))
            except Exception:
                vals.append(np.nan)

        out[f"risk_{t:g}y"] = vals

    return pd.DataFrame(out)


def prepare_model_input(scored_df, bundle, args):
    df = scored_df.copy()

    numeric_cols = list(bundle.get("numeric_cols", []))
    categorical_cols = list(bundle.get("categorical_cols", []))

    if len(numeric_cols) + len(categorical_cols) == 0:
        raise ValueError("Model bundle does not contain numeric_cols/categorical_cols.")

    df = ensure_expected_columns(df, numeric_cols + categorical_cols)

    if "Age" in numeric_cols:
        if args.age_update_mode == "from_baseline_date":
            baseline_age = clean_numeric_series(df["baseline_age_raw"]) if "baseline_age_raw" in df.columns else clean_numeric_series(df["Age"])
            years_since = clean_numeric_series(df["years_since_selected_baseline"])
            derived_age = baseline_age + years_since

            row_age = clean_numeric_series(df["Age"])
            df["Age_raw_from_file"] = row_age
            df["Age"] = derived_age.where(derived_age.notna(), row_age)
            df["age_at_scan_used_for_model"] = df["Age"]

        elif args.age_update_mode == "row":
            df["Age"] = clean_numeric_series(df["Age"])
            df["age_at_scan_used_for_model"] = df["Age"]

        elif args.age_update_mode == "none":
            df["age_at_scan_used_for_model"] = clean_numeric_series(df["Age"])

    if "DLICV" in numeric_cols and "DLICV" in df.columns:
        df["DLICV"] = clean_numeric_series(df["DLICV"])

    if "Sex" in categorical_cols and "Sex" in df.columns:
        df["Sex"] = normalize_sex_series(df["Sex"])

    if "SITE" in categorical_cols and "SITE" in df.columns:
        df["SITE"] = df["SITE"].astype(str).str.strip()

    for c in numeric_cols:
        df[c] = clean_numeric_series(df[c])

    for c in categorical_cols:
        df[c] = df[c].astype("object")

    return df, numeric_cols, categorical_cols


# ============================================================
# Main
# ============================================================

def main():
    args = parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    risk_times = parse_risk_times(args.risk_times)

    baseline_dx = normalize_dx(args.baseline_dx)
    event_dx_set = set(normalize_dx(x) for x in parse_list_arg(args.event_dx))
    eligible_scan_dx_set = set(normalize_dx(x) for x in parse_list_arg(args.eligible_scan_dx))

    log("============================================================")
    log("Apply ADNI brain MRI AD L'EPOCH longitudinally: CN scans only")
    log("============================================================")
    log(f"Input file: {args.input_file}")
    log(f"Model joblib: {args.model_joblib}")
    log(f"Output directory: {outdir}")
    log(f"Baseline DX: {baseline_dx}")
    log(f"Event DX: {sorted(event_dx_set)}")
    log(f"Eligible scan DX for scoring: {sorted(eligible_scan_dx_set)}")
    log(f"Age update mode: {args.age_update_mode}")
    log("Rule: score only CN scans before first MCI/AD conversion.")
    log("============================================================")

    df_raw = read_table(args.input_file)

    required = [args.id_col, args.visit_col, args.date_col, args.dx_col]
    missing_required = [c for c in required if c not in df_raw.columns]
    if missing_required:
        raise ValueError(f"Missing required input columns: {missing_required}")

    bundle = joblib.load(args.model_joblib)

    preprocessor = bundle["preprocessor"]
    model = bundle["model"]
    clock_info = bundle.get("clock_transform_info", None)

    model_roi_cols = get_model_roi_features(bundle)

    scored_df, subject_df, skipped_df, roi_present, roi_missing = prepare_longitudinal_cn_table(
        df_raw=df_raw,
        args=args,
        model_roi_cols=model_roi_cols,
        baseline_dx=baseline_dx,
        event_dx_set=event_dx_set,
        eligible_scan_dx_set=eligible_scan_dx_set
    )

    if scored_df.empty:
        raise ValueError("No qualified longitudinal CN scans were found for application.")

    if args.complete_case_model_features:
        before = scored_df.shape[0]
        scored_df = scored_df.loc[
            scored_df["fraction_model_muse_gm_rois_nonmissing_in_scan"] >= 1.0
        ].copy()
        log(f"Complete-case model ROI filtering: {before} -> {scored_df.shape[0]} scans")

    if scored_df.empty:
        raise ValueError("No scans remain after complete-case feature filtering.")

    model_df, numeric_cols, categorical_cols = prepare_model_input(scored_df, bundle, args)

    X_raw = model_df[numeric_cols + categorical_cols].copy()
    X = preprocessor.transform(X_raw)

    risk = np.asarray(model.predict(X)).reshape(-1)

    pred = compute_clock_transforms_from_saved_info(
        df=model_df,
        risk=risk,
        clock_info=clock_info
    )

    abs_risk = predict_absolute_risk(model, X, risk_times)
    pred = pd.concat([pred.reset_index(drop=True), abs_risk.reset_index(drop=True)], axis=1)

    id_cols = [
        args.id_col,
        args.visit_col,
        args.date_col,
        args.dx_col,
        "scan_dx",
        "scan_date",
        "scan_visit_month",
        "selected_baseline_visit_code",
        "selected_baseline_date",
        "selected_baseline_dx",
        "years_since_selected_baseline",
        "event_from_selected_baseline",
        "event_dx",
        "event_visit_code",
        "event_date",
        "time_to_event_or_censor_years",
        "years_to_event",
        "conversion_group",
        "scan_relation_to_event",
        "Age",
        "age_at_scan_used_for_model",
        "Sex",
        "DLICV",
        "SITE"
    ]

    score_cols = [
        RISK_COL,
        ACCEL_Z_COL,
        ACCEL_YEARS_COL,
        CLOCK_AGE_COL
    ] + [f"risk_{t:g}y" for t in risk_times]

    qc_cols = [
        "n_model_muse_gm_rois_expected",
        "n_model_muse_gm_rois_present_in_file",
        "n_model_muse_gm_rois_missing_from_file",
        "n_model_muse_gm_rois_nonmissing_in_scan",
        "fraction_model_muse_gm_rois_nonmissing_in_scan"
    ]

    output_cols = []
    for c in id_cols + score_cols + qc_cols:
        if c in pred.columns and c not in output_cols:
            output_cols.append(c)

    pred_out = pred[output_cols].copy()

    pred_out = pred_out.rename(columns={
        args.id_col: "PTID",
        args.visit_col: "Visit_Code",
        args.date_col: "Date",
        args.dx_col: "DX_Binary"
    })

    subject_out = subject_df.rename(columns={args.id_col: "PTID"})
    skipped_out = skipped_df.rename(columns={args.id_col: "PTID"})

    out_pred = outdir / f"{args.prefix}_longitudinal_cn_only_predictions.tsv"
    out_subject = outdir / f"{args.prefix}_longitudinal_cn_only_subject_event_summary.tsv"
    out_skipped = outdir / f"{args.prefix}_longitudinal_cn_only_skipped_subjects.tsv"

    pred_out.to_csv(out_pred, sep="\t", index=False)
    subject_out.to_csv(out_subject, sep="\t", index=False)
    skipped_out.to_csv(out_skipped, sep="\t", index=False)

    scan_summary = pred_out.groupby(
        ["conversion_group", "scan_relation_to_event", "scan_dx"],
        dropna=False
    ).size().reset_index(name="n_scans")

    scan_summary.to_csv(
        outdir / f"{args.prefix}_longitudinal_cn_only_scan_summary.tsv",
        sep="\t",
        index=False
    )

    subject_summary = subject_out.groupby(["conversion_group"], dropna=False).agg(
        n_subjects=("PTID", "nunique"),
        n_events=("event_from_selected_baseline", "sum"),
        median_time_to_event_or_censor_years=("time_to_event_or_censor_years", "median"),
        median_n_CN_scans_scored_before_event=("n_CN_scans_scored_before_event", "median")
    ).reset_index()

    subject_summary.to_csv(
        outdir / f"{args.prefix}_longitudinal_cn_only_subject_summary.tsv",
        sep="\t",
        index=False
    )

    summary = {
        "input_file": str(args.input_file),
        "model_joblib": str(args.model_joblib),
        "output_prediction_file": str(out_pred),
        "n_subjects_with_scored_CN_scans": int(pred_out["PTID"].nunique()),
        "n_scored_CN_scans": int(pred_out.shape[0]),
        "n_model_roi_features_expected": int(len(model_roi_cols)),
        "n_model_roi_features_present_in_input": int(len(roi_present)),
        "n_model_roi_features_missing_in_input": int(len(roi_missing)),
        "missing_model_roi_features": roi_missing,
        "baseline_dx": baseline_dx,
        "event_dx": sorted(event_dx_set),
        "eligible_scan_dx": sorted(eligible_scan_dx_set),
        "min_baseline_roi_fraction": float(args.min_baseline_roi_fraction),
        "min_scan_roi_fraction": float(args.min_scan_roi_fraction),
        "age_update_mode": args.age_update_mode,
        "include_selected_baseline": bool(args.include_selected_baseline),
        "risk_times": risk_times,
        "primary_score": RISK_COL,
        "acceleration_z": ACCEL_Z_COL,
        "acceleration_years": ACCEL_YEARS_COL,
        "note": (
            "Scores are generated without refitting. Only CN-labeled MRI scans before first MCI/AD conversion "
            "are scored. MCI, AD, event, and post-event scans are excluded."
        )
    }

    with open(outdir / f"{args.prefix}_longitudinal_cn_only_application_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    pd.DataFrame([summary]).to_csv(
        outdir / f"{args.prefix}_longitudinal_cn_only_application_summary.tsv",
        sep="\t",
        index=False
    )

    log("============================================================")
    log("Longitudinal CN-only AD L'EPOCH application complete.")
    log(f"Predictions: {out_pred}")
    log(f"Subject event summary: {out_subject}")
    log(f"Skipped subjects: {out_skipped}")
    log("Scan summary:")
    log(scan_summary.to_string(index=False))
    log("Subject summary:")
    log(subject_summary.to_string(index=False))
    log("============================================================")


if __name__ == "__main__":
    main()