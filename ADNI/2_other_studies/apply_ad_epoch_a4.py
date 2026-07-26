#!/usr/bin/env python3
"""
Apply the pretrained ADNI brain MRI AD EPOCH model to longitudinal A4 MUSE data.

A4-specific input structure
---------------------------
MUSE.csv
    ID = <BID>_<visit>, for example B10081264_009
    ROI columns are numeric MUSE labels, for example 702, 701, 601, ...
    DLICV is MUSE label 702.

SUBJINFO.csv
    BID, AGEYR, SEX
    A4 derived SEX coding: 1 = Female, 2 = Male.

imaging_volumetric_mri.csv
    BID, VISCODE, Date_DAYS_CONSENT, Date_DAYS_T0

SV.csv
    BID, VISITCD, SVUSEDTC_DAYS_CONSENT, SVUSEDTC_DAYS_T0
    Used only as a visit-timing fallback when the MRI table lacks a matching row.

Scientific behavior
-------------------
1. The saved ADNI preprocessor and Cox model are applied without refitting.
2. A4 numeric ROI columns are renamed internally to the MUSE_Volume_<label>
   names expected by the saved model bundle.
3. Age at scan is AGEYR + days-from-consent / 365.25.
4. DLICV is taken from A4 MUSE column 702, preserving the same MUSE-derived
   definition and scale used for model application.
5. SITE is set to a constant external label (default: A4). Because the saved
   ADNI OneHotEncoder uses handle_unknown='ignore', the unseen A4 site is encoded
   as zero across the learned ADNI site indicators.
6. Rows that cannot be matched to participant demographics, visit timing, or
   sufficient model ROI coverage are written to explicit audit files.
"""

from __future__ import annotations

import argparse
import json
import re
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import joblib
import numpy as np
import pandas as pd

RISK_COL = "adni_brain_mri_ad_epoch_risk_score"
ACCEL_Z_COL = "adni_brain_mri_ad_epoch_acceleration_z"
ACCEL_YEARS_COL = "adni_brain_mri_ad_epoch_acceleration_years"
CLOCK_AGE_COL = "adni_brain_mri_ad_epoch_clock_age_years"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--muse-file", required=True)
    p.add_argument("--subjinfo-file", required=True)
    p.add_argument("--mri-visits-file", required=True)
    p.add_argument("--sv-file", default="")
    p.add_argument("--model-joblib", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--prefix", default="a4_adni_brain_mri_ad_epoch")

    p.add_argument("--muse-id-col", default="ID")
    p.add_argument("--dlicv-label", default="702")
    p.add_argument("--site-label", default="A4")
    p.add_argument("--study-label", default="A4")
    p.add_argument("--min-roi-fraction", type=float, default=0.80)
    p.add_argument("--complete-case-model-rois", action="store_true")
    p.add_argument("--risk-times", default="1,2,3,5")
    return p.parse_args()


def log(message: str) -> None:
    print(message, flush=True)


def read_csv(path: str) -> pd.DataFrame:
    return pd.read_csv(path, low_memory=False)


def clean_numeric(series: pd.Series) -> pd.Series:
    if pd.api.types.is_numeric_dtype(series):
        return pd.to_numeric(series, errors="coerce")
    cleaned = (
        series.astype("object")
        .where(series.notna(), np.nan)
        .astype(str)
        .str.strip()
        .str.replace(",", "", regex=False)
        .replace(
            {
                "": np.nan,
                "nan": np.nan,
                "NaN": np.nan,
                "NA": np.nan,
                "N/A": np.nan,
                "None": np.nan,
                "null": np.nan,
                "-1": np.nan,
            }
        )
    )
    return pd.to_numeric(cleaned, errors="coerce")


def clean_string(series: pd.Series, missing: Optional[str] = None) -> pd.Series:
    cleaned = series.astype("object").where(series.notna(), np.nan)
    cleaned = cleaned.apply(lambda x: str(x).strip() if pd.notna(x) else np.nan)
    cleaned = cleaned.replace(
        {
            "": np.nan,
            "nan": np.nan,
            "NaN": np.nan,
            "NA": np.nan,
            "N/A": np.nan,
            "None": np.nan,
            "null": np.nan,
        }
    )
    return cleaned.fillna(missing) if missing is not None else cleaned


def normalize_a4_sex(series: pd.Series) -> pd.Series:
    """Normalize A4 SUBJINFO SEX coding: 1=Female, 2=Male."""
    s = clean_string(series)
    return s.replace(
        {
            "1": "Female",
            "1.0": "Female",
            "2": "Male",
            "2.0": "Male",
            "F": "Female",
            "M": "Male",
            "female": "Female",
            "male": "Male",
            "FEMALE": "Female",
            "MALE": "Male",
        }
    )


def parse_risk_times(value: str) -> List[float]:
    return [float(x.strip()) for x in str(value).split(",") if x.strip()]


def dedupe_preserve_order(values: Sequence[str]) -> List[str]:
    seen = set()
    output: List[str] = []
    for value in values:
        if value not in seen:
            seen.add(value)
            output.append(value)
    return output


def get_model_rois(bundle: Dict[str, object]) -> List[str]:
    for key in (
        "selected_muse_gm_rois",
        "available_muse_gm_rois",
        "hardcoded_muse_gm_rois",
    ):
        values = bundle.get(key)
        if values:
            return list(values)
    rois = [
        str(x)
        for x in bundle.get("numeric_cols", [])
        if str(x).startswith("MUSE_Volume_")
    ]
    if not rois:
        raise ValueError("No MUSE ROI list was found in the pretrained model bundle.")
    return rois


def muse_label_from_model_roi(model_roi: str) -> str:
    match = re.fullmatch(r"MUSE_Volume_(.+)", str(model_roi))
    if not match:
        raise ValueError(f"Unexpected model ROI name: {model_roi}")
    return match.group(1)


def parse_a4_id(series: pd.Series) -> pd.DataFrame:
    parsed = series.astype(str).str.extract(r"^(?P<BID>.+)_(?P<MUSE_VISIT_RAW>[^_]+)$")
    parsed["MUSE_VISIT_NUM"] = pd.to_numeric(parsed["MUSE_VISIT_RAW"], errors="coerce")
    return parsed


def prepare_visit_timing(mri: pd.DataFrame, sv: Optional[pd.DataFrame]) -> pd.DataFrame:
    required_mri = {"BID", "VISCODE"}
    missing = required_mri.difference(mri.columns)
    if missing:
        raise ValueError(f"MRI visit file is missing columns: {sorted(missing)}")

    keep = [
        c
        for c in ["BID", "VISCODE", "Date_DAYS_CONSENT", "Date_DAYS_T0"]
        if c in mri.columns
    ]
    m = mri.loc[:, keep].copy()
    m["BID"] = clean_string(m["BID"])
    m["VISIT_NUM"] = clean_numeric(m["VISCODE"])
    m["days_from_consent_mri"] = (
        clean_numeric(m["Date_DAYS_CONSENT"])
        if "Date_DAYS_CONSENT" in m.columns
        else np.nan
    )
    m["days_from_t0_mri"] = (
        clean_numeric(m["Date_DAYS_T0"])
        if "Date_DAYS_T0" in m.columns
        else np.nan
    )
    m = m.sort_values(["BID", "VISIT_NUM"]).drop_duplicates(
        ["BID", "VISIT_NUM"], keep="first"
    )

    if sv is None:
        m["days_from_consent_sv"] = np.nan
        m["days_from_t0_sv"] = np.nan
        return m

    required_sv = {"BID", "VISITCD"}
    missing = required_sv.difference(sv.columns)
    if missing:
        raise ValueError(f"SV file is missing columns: {sorted(missing)}")

    keep_sv = [
        c
        for c in [
            "BID",
            "VISITCD",
            "SVUSEDTC_DAYS_CONSENT",
            "SVSTDTC_DAYS_CONSENT",
            "SVUSEDTC_DAYS_T0",
            "SVSTDTC_DAYS_T0",
        ]
        if c in sv.columns
    ]
    s = sv.loc[:, keep_sv].copy()
    s["BID"] = clean_string(s["BID"])
    s["VISIT_NUM"] = clean_numeric(s["VISITCD"])

    consent_candidates = [
        c for c in ["SVUSEDTC_DAYS_CONSENT", "SVSTDTC_DAYS_CONSENT"] if c in s.columns
    ]
    t0_candidates = [c for c in ["SVUSEDTC_DAYS_T0", "SVSTDTC_DAYS_T0"] if c in s.columns]

    s["days_from_consent_sv"] = np.nan
    for c in consent_candidates:
        s["days_from_consent_sv"] = s["days_from_consent_sv"].fillna(clean_numeric(s[c]))
    s["days_from_t0_sv"] = np.nan
    for c in t0_candidates:
        s["days_from_t0_sv"] = s["days_from_t0_sv"].fillna(clean_numeric(s[c]))

    s = s.sort_values(["BID", "VISIT_NUM"]).drop_duplicates(
        ["BID", "VISIT_NUM"], keep="first"
    )
    return m.merge(
        s[["BID", "VISIT_NUM", "days_from_consent_sv", "days_from_t0_sv"]],
        on=["BID", "VISIT_NUM"],
        how="outer",
    )


def categorical_match(series: pd.Series, category: object) -> np.ndarray:
    strings = clean_string(series)
    mask = strings.astype(str).eq(str(category).strip()).to_numpy(dtype=bool)
    try:
        numeric = pd.to_numeric(series, errors="coerce").to_numpy(dtype=float)
        mask |= np.isclose(numeric, float(category), equal_nan=False)
    except Exception:
        pass
    return mask.astype(float)


def parse_cat_term(term: str, categorical_covariates: Sequence[str]) -> Tuple[Optional[str], Optional[str]]:
    if not term.startswith("cat__"):
        return None, None
    stem = term[len("cat__") :]
    for covariate in sorted(categorical_covariates, key=len, reverse=True):
        prefix = covariate + "_"
        if stem.startswith(prefix):
            return covariate, stem[len(prefix) :]
    return None, None


def compute_clock_transforms(
    df: pd.DataFrame,
    risk: np.ndarray,
    info: Optional[Dict[str, object]],
) -> pd.DataFrame:
    out = df.copy()
    out[RISK_COL] = risk
    if not info:
        warnings.warn("clock_transform_info is missing; acceleration outputs set to NA.")
        out[ACCEL_Z_COL] = np.nan
        out[ACCEL_YEARS_COL] = np.nan
        out[CLOCK_AGE_COL] = np.nan
        return out

    expected = np.repeat(float(info.get("risk_score_covariate_model_intercept", 0.0)), len(out))
    categorical_covariates = list(info.get("categorical_residualization_covariates", []))

    for term, beta in dict(info.get("risk_score_covariate_model_coef", {})).items():
        beta = float(beta)
        if term.startswith("num__"):
            covariate = term[len("num__") :]
            if covariate not in out.columns:
                values = np.zeros(len(out), dtype=float)
            else:
                values = clean_numeric(out[covariate])
                median = float(np.nanmedian(values)) if values.notna().any() else 0.0
                values = values.fillna(median).to_numpy(dtype=float)
            expected += beta * values
        elif term.startswith("cat__"):
            covariate, category = parse_cat_term(term, categorical_covariates)
            if covariate is not None and covariate in out.columns:
                expected += beta * categorical_match(out[covariate], category)

    residual_raw = risk - expected
    residual_centered = residual_raw - float(info.get("risk_score_residual_mean_train", 0.0))
    residual_sd = float(info.get("risk_score_residual_sd_train", np.nan))

    out["adni_expected_risk_from_saved_covariates"] = expected
    out["adni_risk_residual_raw"] = residual_raw
    out[ACCEL_Z_COL] = (
        residual_centered / residual_sd
        if np.isfinite(residual_sd) and residual_sd > 0
        else np.nan
    )

    beta_age = info.get("adjusted_age_coefficient_risk_score_per_year")
    beta_age = float(beta_age) if beta_age is not None else np.nan
    if np.isfinite(beta_age) and abs(beta_age) > 1e-8:
        out[ACCEL_YEARS_COL] = residual_centered / beta_age
        out[CLOCK_AGE_COL] = clean_numeric(out["Age"]) + out[ACCEL_YEARS_COL]
    else:
        out[ACCEL_YEARS_COL] = np.nan
        out[CLOCK_AGE_COL] = np.nan
    return out


def predict_absolute_risk(model: object, X: np.ndarray, times: Sequence[float]) -> pd.DataFrame:
    output: Dict[str, object] = {}
    try:
        functions = model.predict_survival_function(X)
        for time in times:
            values = []
            for function in functions:
                try:
                    values.append(1.0 - float(function(time)))
                except Exception:
                    values.append(np.nan)
            output[f"risk_{time:g}y"] = values
    except Exception as exc:
        warnings.warn(f"Absolute-risk prediction failed: {exc}")
        for time in times:
            output[f"risk_{time:g}y"] = np.repeat(np.nan, X.shape[0])
    return pd.DataFrame(output)


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    muse = read_csv(args.muse_file)
    subj = read_csv(args.subjinfo_file)
    mri = read_csv(args.mri_visits_file)
    sv = read_csv(args.sv_file) if args.sv_file else None

    if args.muse_id_col not in muse.columns:
        raise ValueError(f"MUSE file lacks ID column: {args.muse_id_col}")
    for col in ["BID", "AGEYR", "SEX"]:
        if col not in subj.columns:
            raise ValueError(f"SUBJINFO file lacks required column: {col}")

    parsed = parse_a4_id(muse[args.muse_id_col])
    muse = pd.concat([muse.copy(), parsed], axis=1)
    muse["_source_row_number"] = np.arange(1, len(muse) + 1)
    muse["BID"] = clean_string(muse["BID"])

    timing = prepare_visit_timing(mri, sv)
    timing["BID"] = clean_string(timing["BID"])

    subj_small = subj[["BID", "AGEYR", "SEX"]].copy()
    subj_small["BID"] = clean_string(subj_small["BID"])
    subj_small["AGEYR"] = clean_numeric(subj_small["AGEYR"])
    subj_small["Sex"] = normalize_a4_sex(subj_small["SEX"])
    subj_small = subj_small.drop_duplicates("BID", keep="first")

    data = muse.merge(subj_small, on="BID", how="left", validate="many_to_one")
    data = data.merge(
        timing,
        left_on=["BID", "MUSE_VISIT_NUM"],
        right_on=["BID", "VISIT_NUM"],
        how="left",
        validate="many_to_one",
    )

    data["days_from_consent_used"] = clean_numeric(data.get("days_from_consent_mri", pd.Series(np.nan, index=data.index)))
    data["timing_source"] = np.where(data["days_from_consent_used"].notna(), "imaging_volumetric_mri", "missing")
    if "days_from_consent_sv" in data.columns:
        sv_mask = data["days_from_consent_used"].isna() & clean_numeric(data["days_from_consent_sv"]).notna()
        data.loc[sv_mask, "days_from_consent_used"] = clean_numeric(data.loc[sv_mask, "days_from_consent_sv"])
        data.loc[sv_mask, "timing_source"] = "SV_fallback"

    data["Age"] = data["AGEYR"] + data["days_from_consent_used"] / 365.25
    data["DLICV"] = clean_numeric(data[args.dlicv_label]) if args.dlicv_label in data.columns else np.nan
    data["SITE"] = args.site_label
    data["Study"] = args.study_label
    data["Visit_Code"] = data["MUSE_VISIT_RAW"]
    data["Date_DAYS_CONSENT"] = data["days_from_consent_used"]
    data["years_since_consent"] = data["days_from_consent_used"] / 365.25

    bundle = joblib.load(args.model_joblib)
    preprocessor = bundle["preprocessor"]
    model = bundle["model"]
    model_rois = get_model_rois(bundle)

    roi_audit_rows: List[Dict[str, object]] = []
    for model_roi in model_rois:
        label = muse_label_from_model_roi(model_roi)
        present = label in data.columns
        data[model_roi] = clean_numeric(data[label]) if present else np.nan
        roi_audit_rows.append(
            {
                "model_roi": model_roi,
                "a4_muse_label": label,
                "a4_column_present": bool(present),
                "n_nonmissing": int(data[model_roi].notna().sum()),
                "n_missing": int(data[model_roi].isna().sum()),
            }
        )
    roi_audit = pd.DataFrame(roi_audit_rows)

    data["n_model_rois_expected"] = len(model_rois)
    data["n_model_rois_nonmissing_in_scan"] = data[model_rois].notna().sum(axis=1)
    data["fraction_model_rois_nonmissing_in_scan"] = (
        data["n_model_rois_nonmissing_in_scan"] / float(len(model_rois))
    )

    threshold = 1.0 if args.complete_case_model_rois else args.min_roi_fraction

    exclusion_reasons = []
    for _, row in data.iterrows():
        reasons: List[str] = []
        if pd.isna(row["BID"]):
            reasons.append("invalid_ID_parse")
        if pd.isna(row["AGEYR"]):
            reasons.append("missing_SUBJINFO_AGEYR")
        if pd.isna(row["Sex"]):
            reasons.append("missing_or_unrecognized_SUBJINFO_SEX")
        if pd.isna(row["days_from_consent_used"]):
            reasons.append("unmatched_visit_timing")
        if pd.isna(row["Age"]):
            reasons.append("missing_age_at_scan")
        if pd.isna(row["DLICV"]):
            reasons.append("missing_MUSE_702_DLICV")
        if row["fraction_model_rois_nonmissing_in_scan"] < threshold:
            reasons.append("insufficient_model_ROI_coverage")
        exclusion_reasons.append(";".join(reasons))

    data["exclusion_reason"] = exclusion_reasons
    eligible = data["exclusion_reason"].eq("")
    excluded = data.loc[~eligible].copy()
    scored = data.loc[eligible].copy()

    if scored.empty:
        raise ValueError(
            "No A4 scans passed QC. Inspect the excluded-scan and visit-match audit outputs."
        )

    scored = scored.sort_values(["BID", "days_from_consent_used", "_source_row_number"], kind="mergesort")
    scored["longitudinal_scan_number"] = scored.groupby("BID", sort=False).cumcount() + 1
    scored["is_external_baseline_scan"] = scored["longitudinal_scan_number"].eq(1)
    baseline_days = scored.groupby("BID", sort=False)["days_from_consent_used"].transform("min")
    scored["years_since_external_baseline"] = (
        scored["days_from_consent_used"] - baseline_days
    ) / 365.25

    numeric_cols = list(bundle.get("numeric_cols", []))
    categorical_cols = list(bundle.get("categorical_cols", []))
    if not numeric_cols and not categorical_cols:
        raise ValueError("Model bundle lacks numeric_cols/categorical_cols.")

    created_missing: List[str] = []
    for column in numeric_cols + categorical_cols:
        if column not in scored.columns:
            scored[column] = np.nan
            created_missing.append(column)

    for column in numeric_cols:
        scored[column] = clean_numeric(scored[column])
    for column in categorical_cols:
        scored[column] = scored[column].astype("object")

    X_raw = scored[numeric_cols + categorical_cols].copy()
    X = preprocessor.transform(X_raw)
    risk = np.asarray(model.predict(X)).reshape(-1)

    prediction = compute_clock_transforms(
        scored,
        risk,
        bundle.get("clock_transform_info"),
    )
    prediction = pd.concat(
        [
            prediction.reset_index(drop=True),
            predict_absolute_risk(model, X, parse_risk_times(args.risk_times)),
        ],
        axis=1,
    )

    requested_output_columns = [
        "_source_row_number",
        args.muse_id_col,
        "BID",
        "Visit_Code",
        "MUSE_VISIT_NUM",
        "VISIT_NUM",
        "timing_source",
        "Date_DAYS_CONSENT",
        "years_since_consent",
        "years_since_external_baseline",
        "longitudinal_scan_number",
        "is_external_baseline_scan",
        "Study",
        "SITE",
        "AGEYR",
        "Age",
        "SEX",
        "Sex",
        args.dlicv_label,
        "DLICV",
        RISK_COL,
        ACCEL_Z_COL,
        ACCEL_YEARS_COL,
        CLOCK_AGE_COL,
        "adni_expected_risk_from_saved_covariates",
        "adni_risk_residual_raw",
    ]
    requested_output_columns += [
        f"risk_{time:g}y" for time in parse_risk_times(args.risk_times)
    ]
    requested_output_columns += [
        "n_model_rois_expected",
        "n_model_rois_nonmissing_in_scan",
        "fraction_model_rois_nonmissing_in_scan",
    ]

    output_columns = dedupe_preserve_order(
        [x for x in requested_output_columns if x in prediction.columns]
    )
    scan = prediction.loc[:, output_columns].copy()

    scan_file = outdir / f"{args.prefix}_scan_level_predictions.tsv"
    subject_file = outdir / f"{args.prefix}_subject_longitudinal_summary.tsv"
    excluded_file = outdir / f"{args.prefix}_excluded_scans.tsv"
    roi_file = outdir / f"{args.prefix}_model_roi_mapping_audit.tsv"
    visit_file = outdir / f"{args.prefix}_visit_match_audit.tsv"

    scan.to_csv(scan_file, sep="\t", index=False)
    excluded.to_csv(excluded_file, sep="\t", index=False)
    roi_audit.to_csv(roi_file, sep="\t", index=False)

    visit_audit_columns = [
        c
        for c in [
            args.muse_id_col,
            "BID",
            "MUSE_VISIT_RAW",
            "MUSE_VISIT_NUM",
            "VISIT_NUM",
            "days_from_consent_mri",
            "days_from_consent_sv",
            "days_from_consent_used",
            "timing_source",
            "exclusion_reason",
        ]
        if c in data.columns
    ]
    data[visit_audit_columns].to_csv(visit_file, sep="\t", index=False)

    subject = (
        scan.groupby("BID", dropna=False, sort=False)
        .agg(
            n_scans=("BID", "size"),
            baseline_age=("Age", "first"),
            last_age=("Age", "last"),
            baseline_risk=(RISK_COL, "first"),
            last_risk=(RISK_COL, "last"),
            baseline_acceleration_z=(ACCEL_Z_COL, "first"),
            last_acceleration_z=(ACCEL_Z_COL, "last"),
            max_followup_years=("years_since_external_baseline", "max"),
        )
        .reset_index()
    )
    subject["change_risk_first_to_last"] = subject["last_risk"] - subject["baseline_risk"]
    subject["change_acceleration_z_first_to_last"] = (
        subject["last_acceleration_z"] - subject["baseline_acceleration_z"]
    )
    subject.to_csv(subject_file, sep="\t", index=False)

    summary = {
        "muse_file": args.muse_file,
        "subjinfo_file": args.subjinfo_file,
        "mri_visits_file": args.mri_visits_file,
        "sv_file": args.sv_file,
        "model_joblib": args.model_joblib,
        "n_input_scans": int(len(data)),
        "n_scored_scans": int(len(scan)),
        "n_excluded_scans": int(len(excluded)),
        "n_scored_subjects": int(scan["BID"].nunique()),
        "n_model_rois_expected": int(len(model_rois)),
        "minimum_roi_fraction": float(threshold),
        "n_mri_timing_matches": int((data["timing_source"] == "imaging_volumetric_mri").sum()),
        "n_sv_timing_fallbacks": int((data["timing_source"] == "SV_fallback").sum()),
        "n_unmatched_timing": int((data["timing_source"] == "missing").sum()),
        "created_missing_expected_columns": created_missing,
        "model_refit": False,
        "age_rule": "AGEYR + days_from_consent_used / 365.25",
        "sex_rule": "A4 SUBJINFO SEX: 1=Female, 2=Male",
        "dlicv_rule": f"MUSE label {args.dlicv_label}",
        "site_rule": f"Constant external label: {args.site_label}",
        "visit_rule": (
            "MUSE ID suffix is converted to numeric and matched to imaging VISCODE; "
            "SV VISITCD is used only as a timing fallback."
        ),
    }

    with open(outdir / f"{args.prefix}_application_summary.json", "w") as handle:
        json.dump(summary, handle, indent=2)
    pd.DataFrame([summary]).to_csv(
        outdir / f"{args.prefix}_application_summary.tsv",
        sep="\t",
        index=False,
    )

    log(f"Scored {len(scan):,} scans from {scan['BID'].nunique():,} A4 participants.")
    log(f"Excluded scans: {len(excluded):,}")
    log(f"Predictions: {scan_file}")
    log(f"Visit-match audit: {visit_file}")
    log(f"ROI mapping audit: {roi_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
