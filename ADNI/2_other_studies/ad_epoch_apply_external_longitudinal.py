#!/usr/bin/env python3
"""Apply a pretrained ADNI brain MRI AD EPOCH to external longitudinal scans.

The saved ADNI preprocessing pipeline and Cox model are applied without
refitting. Every scan passing ROI coverage QC is scored, regardless of CN,
MCI, AD/CD, or other diagnosis. Repeated scans are retained.

The Python script:

loads the pretrained ADNI .joblib model without refitting;
scores every eligible longitudinal scan, regardless of CN, MCI, AD/CD, or other diagnosis;
retains repeated scans and calculates scan order and years since each participant’s earliest eligible scan;
uses scan-level age by default;
uses the raw MUSE_Volume_* variables expected by the pretrained model, not the harmonized H_MUSE_Volume_* variables;
requires at least 80% of the model-selected ROIs to be nonmissing;
safely accepts all external SITE labels because the saved ADNI encoder uses handle_unknown="ignore";
treats missing SITE values as an explicit unseen external category rather than imputing them to an ADNI site;
generates the Cox risk score, acceleration z score, year-scale acceleration, clock age, and 1-, 2-, 3-, and 5-year absolute risks.

"""

from __future__ import annotations

import argparse
import json
import re
import warnings
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

RISK_COL = "adni_brain_mri_ad_epoch_risk_score"
ACCEL_Z_COL = "adni_brain_mri_ad_epoch_acceleration_z"
ACCEL_YEARS_COL = "adni_brain_mri_ad_epoch_acceleration_years"
CLOCK_AGE_COL = "adni_brain_mri_ad_epoch_clock_age_years"


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input-file", required=True)
    p.add_argument("--model-joblib", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--prefix", default="external_adni_brain_mri_ad_epoch")
    p.add_argument("--id-col", default="PTID")
    p.add_argument("--visit-col", default="Visit_Code")
    p.add_argument("--date-col", default="Date")
    p.add_argument("--age-col", default="Age")
    p.add_argument("--sex-col", default="Sex")
    p.add_argument("--dlicv-col", default="DLICV")
    p.add_argument("--site-col", default="SITE")
    p.add_argument("--study-col", default="Study")
    p.add_argument("--dx-col", default="DX_Binary")
    p.add_argument("--delta-baseline-col", default="Delta_Baseline")
    p.add_argument("--study-values", default="")
    p.add_argument("--min-roi-fraction", type=float, default=0.80)
    p.add_argument("--complete-case-model-rois", action="store_true")
    p.add_argument("--age-mode", choices=["row", "baseline_plus_time"], default="row")
    p.add_argument("--risk-times", default="1,2,3,5")
    return p.parse_args()


def log(x):
    print(x, flush=True)


def read_table(path):
    path = Path(path)
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path, low_memory=False)
    return pd.read_csv(path, sep="\t", low_memory=False)


def parse_list(x):
    return [v.strip() for v in str(x).split(",") if v.strip()]


def parse_risk_times(x):
    return [float(v) for v in parse_list(x)]


def clean_numeric(s):
    if pd.api.types.is_numeric_dtype(s):
        return pd.to_numeric(s, errors="coerce")
    x = (
        s.astype("object").where(s.notna(), np.nan).astype(str).str.strip()
        .str.replace(",", "", regex=False)
        .replace({"": np.nan, "nan": np.nan, "NaN": np.nan, "NA": np.nan,
                  "N/A": np.nan, "None": np.nan, "null": np.nan, "-1": np.nan})
    )
    return pd.to_numeric(x, errors="coerce")


def clean_string(s, missing=None):
    x = s.astype("object").where(s.notna(), np.nan)
    x = x.apply(lambda v: str(v).strip() if pd.notna(v) else np.nan)
    x = x.replace({"": np.nan, "nan": np.nan, "NaN": np.nan, "NA": np.nan,
                   "N/A": np.nan, "None": np.nan, "null": np.nan})
    return x.fillna(missing) if missing is not None else x


def normalize_sex(s):
    return clean_string(s).replace({
        "0": "Female", "0.0": "Female", "1": "Male", "1.0": "Male",
        "F": "Female", "M": "Male", "female": "Female", "male": "Male",
        "FEMALE": "Female", "MALE": "Male",
    })


def normalize_dx(x):
    if pd.isna(x):
        return np.nan
    y = str(x).strip().upper()
    return {"NL": "CN", "NORMAL": "CN", "CONTROL": "CN", "HC": "CN",
            "EMCI": "MCI", "LMCI": "MCI", "DEMENTIA": "AD", "CD": "AD"}.get(y, y)


def visit_month(x):
    if pd.isna(x):
        return np.nan
    y = str(x).strip().lower()
    if y in {"bl", "base", "baseline", "m00", "m0", "screen", "screening", "sc"}:
        return 0.0
    m = re.search(r"m(?:onth)?\s*0*([0-9]+)", y)
    if m:
        return float(m.group(1))
    m = re.search(r"([0-9]+)", y)
    return float(m.group(1)) if m else np.nan


def model_rois(bundle):
    for key in ["selected_muse_gm_rois", "available_muse_gm_rois", "hardcoded_muse_gm_rois"]:
        if bundle.get(key):
            return list(bundle[key])
    cols = [c for c in bundle.get("numeric_cols", []) if str(c).startswith("MUSE_Volume_")]
    if not cols:
        raise ValueError("No MUSE ROI list found in model bundle.")
    return cols


def categorical_match(series, category):
    s = clean_string(series)
    mask = s.astype(str).eq(str(category).strip()).to_numpy(dtype=bool)
    try:
        num = pd.to_numeric(series, errors="coerce").to_numpy(dtype=float)
        mask |= np.isclose(num, float(category), equal_nan=False)
    except Exception:
        pass
    return mask.astype(float)


def parse_cat_term(term, categorical_covs):
    if not term.startswith("cat__"):
        return None, None
    stem = term[5:]
    for cov in sorted(categorical_covs, key=len, reverse=True):
        prefix = cov + "_"
        if stem.startswith(prefix):
            return cov, stem[len(prefix):]
    return None, None


def clock_transforms(df, risk, info):
    out = df.copy()
    out[RISK_COL] = risk
    if not info:
        out[ACCEL_Z_COL] = np.nan
        out[ACCEL_YEARS_COL] = np.nan
        out[CLOCK_AGE_COL] = np.nan
        return out

    expected = np.repeat(float(info.get("risk_score_covariate_model_intercept", 0.0)), len(out))
    cat_covs = list(info.get("categorical_residualization_covariates", []))
    for term, beta in dict(info.get("risk_score_covariate_model_coef", {})).items():
        beta = float(beta)
        if term.startswith("num__"):
            cov = term[5:]
            if cov not in out.columns:
                vals = np.zeros(len(out))
            else:
                vals = clean_numeric(out[cov])
                vals = vals.fillna(float(np.nanmedian(vals)) if vals.notna().any() else 0.0).to_numpy()
            expected += beta * vals
        elif term.startswith("cat__"):
            cov, cat = parse_cat_term(term, cat_covs)
            if cov is not None and cov in out.columns:
                expected += beta * categorical_match(out[cov], cat)

    resid_raw = risk - expected
    resid = resid_raw - float(info.get("risk_score_residual_mean_train", 0.0))
    resid_sd = float(info.get("risk_score_residual_sd_train", np.nan))
    out["adni_expected_risk_from_saved_covariates"] = expected
    out["adni_risk_residual_raw"] = resid_raw
    out[ACCEL_Z_COL] = resid / resid_sd if np.isfinite(resid_sd) and resid_sd > 0 else np.nan

    beta_age = info.get("adjusted_age_coefficient_risk_score_per_year")
    beta_age = float(beta_age) if beta_age is not None else np.nan
    if np.isfinite(beta_age) and abs(beta_age) > 1e-8:
        out[ACCEL_YEARS_COL] = resid / beta_age
        out[CLOCK_AGE_COL] = clean_numeric(out["Age"]) + out[ACCEL_YEARS_COL]
    else:
        out[ACCEL_YEARS_COL] = np.nan
        out[CLOCK_AGE_COL] = np.nan
    return out


def absolute_risk(model, X, times):
    data = {}
    try:
        funcs = model.predict_survival_function(X)
        for t in times:
            vals = []
            for f in funcs:
                try:
                    vals.append(1.0 - float(f(t)))
                except Exception:
                    vals.append(np.nan)
            data[f"risk_{t:g}y"] = vals
    except Exception as exc:
        warnings.warn(f"Absolute-risk prediction failed: {exc}")
        for t in times:
            data[f"risk_{t:g}y"] = np.repeat(np.nan, X.shape[0])
    return pd.DataFrame(data)


def longitudinal_metadata(df, args):
    d = df.copy()
    d["_scan_date"] = pd.to_datetime(d[args.date_col], errors="coerce")
    d["_visit_month"] = d[args.visit_col].apply(visit_month)
    delta = clean_numeric(d[args.delta_baseline_col]) if args.delta_baseline_col in d.columns else pd.Series(np.nan, index=d.index)
    d["_delta_years"] = delta
    d["_sort"] = d["_scan_date"].map(lambda x: x.toordinal() if pd.notna(x) else np.nan)
    d["_sort"] = d["_sort"].fillna(d["_delta_years"] * 365.25).fillna(d["_visit_month"] * 30.4375)
    d = d.sort_values([args.id_col, "_sort"], kind="mergesort", na_position="last")

    rows = []
    for pid, g in d.groupby(args.id_col, sort=False):
        ordered = g[g["_sort"].notna()]
        b = ordered.iloc[0] if not ordered.empty else g.iloc[0]
        rows.append({args.id_col: pid, "_base_date": b["_scan_date"],
                     "_base_visit_month": b["_visit_month"],
                     "_base_age": clean_numeric(pd.Series([b[args.age_col]])).iloc[0],
                     "external_selected_baseline_visit": b[args.visit_col],
                     "external_selected_baseline_date": b["_scan_date"]})
    d = d.merge(pd.DataFrame(rows), on=args.id_col, how="left")

    yrs = pd.Series(np.nan, index=d.index, dtype=float)
    m = d["_scan_date"].notna() & d["_base_date"].notna()
    yrs.loc[m] = (d.loc[m, "_scan_date"] - d.loc[m, "_base_date"]).dt.days / 365.25
    m = yrs.isna() & d["_delta_years"].notna()
    yrs.loc[m] = d.loc[m, "_delta_years"]
    m = yrs.isna() & d["_visit_month"].notna() & d["_base_visit_month"].notna()
    yrs.loc[m] = (d.loc[m, "_visit_month"] - d.loc[m, "_base_visit_month"]) / 12.0
    d["years_since_external_baseline"] = yrs
    d["longitudinal_scan_number"] = d.groupby(args.id_col, sort=False).cumcount() + 1
    d["is_external_baseline_scan"] = d["longitudinal_scan_number"].eq(1)

    row_age = clean_numeric(d[args.age_col])
    if args.age_mode == "baseline_plus_time":
        derived = clean_numeric(d["_base_age"]) + d["years_since_external_baseline"]
        d["Age_original_external"] = row_age
        d["Age"] = derived.where(derived.notna(), row_age)
    else:
        d["Age"] = row_age
    d["age_at_scan_used_for_model"] = d["Age"]
    return d


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    raw = read_table(args.input_file)
    required = [args.id_col, args.visit_col, args.date_col, args.age_col]
    missing = [c for c in required if c not in raw.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    raw = raw[raw[args.id_col].notna()].copy()
    raw["_source_row_number"] = np.arange(1, len(raw) + 1)

    if args.study_values:
        studies = set(parse_list(args.study_values))
        raw = raw[clean_string(raw[args.study_col]).isin(studies)].copy()
    if raw.empty:
        raise ValueError("No scans remain after study filtering.")

    bundle = joblib.load(args.model_joblib)
    preprocessor = bundle["preprocessor"]
    model = bundle["model"]
    rois = model_rois(bundle)

    present = [c for c in rois if c in raw.columns]
    missing_rois = [c for c in rois if c not in raw.columns]
    if not present:
        raise ValueError("No raw MUSE_Volume_* model features are present.")

    roi_data = pd.DataFrame(index=raw.index)
    for roi in rois:
        roi_data[roi] = clean_numeric(raw[roi]) if roi in raw.columns else np.nan
    raw["n_model_rois_expected"] = len(rois)
    raw["n_model_rois_present_in_file"] = len(present)
    raw["n_model_rois_missing_from_file"] = len(missing_rois)
    raw["n_model_rois_nonmissing_in_scan"] = roi_data.notna().sum(axis=1)
    raw["fraction_model_rois_nonmissing_in_scan"] = raw["n_model_rois_nonmissing_in_scan"] / float(len(rois))

    threshold = 1.0 if args.complete_case_model_rois else args.min_roi_fraction
    eligible = raw["fraction_model_rois_nonmissing_in_scan"].ge(threshold)
    excluded = raw[~eligible].copy()
    excluded["exclusion_reason"] = "insufficient_model_MUSE_ROI_coverage"
    scored = raw[eligible].copy()
    if scored.empty:
        raise ValueError("No scans passed ROI coverage QC.")

    scored = longitudinal_metadata(scored, args)
    scored["diagnosis_normalized"] = scored[args.dx_col].apply(normalize_dx) if args.dx_col in scored.columns else np.nan

    numeric_cols = list(bundle.get("numeric_cols", []))
    categorical_cols = list(bundle.get("categorical_cols", []))
    if not numeric_cols and not categorical_cols:
        raise ValueError("Model bundle lacks numeric_cols/categorical_cols.")
    created_missing = []
    for c in numeric_cols + categorical_cols:
        if c not in scored.columns:
            scored[c] = np.nan
            created_missing.append(c)

    if "DLICV" in numeric_cols:
        scored["DLICV"] = clean_numeric(scored[args.dlicv_col]) if args.dlicv_col in scored.columns else np.nan
    if "Sex" in categorical_cols:
        scored["Sex"] = normalize_sex(scored[args.sex_col]) if args.sex_col in scored.columns else np.nan
    if "SITE" in categorical_cols:
        # All external SITE labels are unseen by the ADNI encoder and therefore
        # receive zero across the ADNI SITE dummy variables.
        scored["SITE"] = clean_string(scored[args.site_col], "__EXTERNAL_UNKNOWN__") if args.site_col in scored.columns else "__EXTERNAL_UNKNOWN__"
    for c in numeric_cols:
        scored[c] = clean_numeric(scored[c])
    for c in categorical_cols:
        scored[c] = scored[c].astype("object")

    X = preprocessor.transform(scored[numeric_cols + categorical_cols])
    risk = np.asarray(model.predict(X)).reshape(-1)
    pred = clock_transforms(scored, risk, bundle.get("clock_transform_info"))
    pred = pd.concat([pred.reset_index(drop=True), absolute_risk(model, X, parse_risk_times(args.risk_times))], axis=1)

    keep = ["_source_row_number", args.id_col, args.visit_col, args.date_col,
            args.study_col, args.site_col, args.dx_col, "diagnosis_normalized",
            "longitudinal_scan_number", "is_external_baseline_scan",
            "external_selected_baseline_visit", "external_selected_baseline_date",
            "years_since_external_baseline", "Age", "age_at_scan_used_for_model",
            args.sex_col, args.dlicv_col, "SITE", "MRI_ID", "MRID",
            "MRI_Magnetic_Field_Strength", "MRI_Manufacturer", "MRI_Scanner_Model",
            "MRI_Protocol", RISK_COL, ACCEL_Z_COL, ACCEL_YEARS_COL, CLOCK_AGE_COL,
            "adni_expected_risk_from_saved_covariates", "adni_risk_residual_raw"]
    keep += [f"risk_{t:g}y" for t in parse_risk_times(args.risk_times)]
    keep += ["n_model_rois_expected", "n_model_rois_present_in_file",
             "n_model_rois_missing_from_file", "n_model_rois_nonmissing_in_scan",
             "fraction_model_rois_nonmissing_in_scan"]
    keep = [c for c in keep if c in pred.columns]
    scan = pred[keep].copy()

    scan_file = outdir / f"{args.prefix}_scan_level_predictions.tsv"
    subject_file = outdir / f"{args.prefix}_subject_longitudinal_summary.tsv"
    site_file = outdir / f"{args.prefix}_study_site_summary.tsv"
    excluded_file = outdir / f"{args.prefix}_excluded_scans.tsv"
    roi_file = outdir / f"{args.prefix}_model_roi_compatibility.tsv"
    scan.to_csv(scan_file, sep="\t", index=False)
    excluded.to_csv(excluded_file, sep="\t", index=False)
    pd.DataFrame({"model_roi": rois, "present_in_external_file": [r in present for r in rois]}).to_csv(roi_file, sep="\t", index=False)

    subject = scan.groupby(args.id_col, dropna=False).agg(
        n_scans=(args.id_col, "size"),
        baseline_risk=(RISK_COL, "first"), last_risk=(RISK_COL, "last"),
        baseline_acceleration_z=(ACCEL_Z_COL, "first"), last_acceleration_z=(ACCEL_Z_COL, "last"),
        max_followup_years=("years_since_external_baseline", "max")
    ).reset_index()
    subject["change_risk_first_to_last"] = subject["last_risk"] - subject["baseline_risk"]
    subject["change_acceleration_z_first_to_last"] = subject["last_acceleration_z"] - subject["baseline_acceleration_z"]
    if args.study_col in scan.columns:
        first_study = scan.groupby(args.id_col, dropna=False)[args.study_col].first().reset_index()
        subject = first_study.merge(subject, on=args.id_col, how="right")
    subject.to_csv(subject_file, sep="\t", index=False)

    groups = [c for c in [args.study_col, args.site_col] if c in scan.columns]
    site_summary = scan.groupby(groups, dropna=False).agg(
        n_subjects=(args.id_col, "nunique"), n_scans=(args.id_col, "size"),
        median_age=("Age", "median"), median_roi_coverage=("fraction_model_rois_nonmissing_in_scan", "median"),
        mean_risk=(RISK_COL, "mean"), sd_risk=(RISK_COL, "std"),
        mean_acceleration_z=(ACCEL_Z_COL, "mean"), sd_acceleration_z=(ACCEL_Z_COL, "std")
    ).reset_index() if groups else pd.DataFrame({"n_subjects": [scan[args.id_col].nunique()], "n_scans": [len(scan)]})
    site_summary.to_csv(site_file, sep="\t", index=False)

    summary = {
        "input_file": args.input_file, "model_joblib": args.model_joblib,
        "n_input_scans": int(len(raw)), "n_scored_scans": int(len(scan)),
        "n_excluded_scans": int(len(excluded)), "n_scored_subjects": int(scan[args.id_col].nunique()),
        "n_model_rois_expected": len(rois), "n_model_rois_present": len(present),
        "missing_model_rois": missing_rois, "created_missing_expected_columns": created_missing,
        "minimum_roi_fraction": threshold, "model_refit": False,
        "diagnosis_rule": "All eligible scans are scored regardless of diagnosis.",
        "site_rule": "External SITE values are unseen and are ignored by the saved ADNI one-hot encoder.",
        "feature_rule": "Raw MUSE_Volume_* features are used; H_MUSE_Volume_* fields are not substituted."
    }
    with open(outdir / f"{args.prefix}_application_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    pd.DataFrame([summary]).to_csv(outdir / f"{args.prefix}_application_summary.tsv", sep="\t", index=False)

    log(f"Scored {len(scan):,} scans from {scan[args.id_col].nunique():,} participants.")
    log(f"Predictions: {scan_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
