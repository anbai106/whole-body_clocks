#!/usr/bin/env python3
"""
Apply a pretrained ADNI brain MRI AD EPOCH model to external longitudinal MRI scans.

Key behavior
------------
1. Application only: the saved ADNI preprocessor and Cox model are not refit.
2. All eligible longitudinal scans are scored regardless of CN/MCI/AD/CD status.
3. Harmonized H_MUSE_Volume_* features can be mapped to the raw MUSE_Volume_*
   names expected by the pretrained model.
4. External SITE labels are preserved for summaries but are passed to the saved
   ADNI encoder as unseen categories. The training encoder uses
   OneHotEncoder(handle_unknown='ignore'), so unseen sites are represented by
   zeros across the learned ADNI SITE indicators.
5. Output column names are explicitly deduplicated, preventing the pandas error
   "Grouper for 'SITE' not 1-dimensional".

Important scientific note
-------------------------
The pretrained ADNI model was fitted using raw MUSE_Volume_* values. Harmonized
features may reduce scanner/study distribution shift, but they are not the exact
feature distribution used to fit the saved scaler and Cox coefficients. This
script therefore records the source used for every model ROI and supports three
modes:
  harmonized : require/use H_MUSE_Volume_* whenever available; optional raw fallback
  raw        : use only MUSE_Volume_*
  auto       : prefer H_MUSE_Volume_* and otherwise use MUSE_Volume_*

This application script supports raw and harmonized feature modes. The companion Slurm script runs both modes by default and saves them separately. Use --roi-source raw for strict feature replication and --roi-source harmonized for the harmonized sensitivity analysis.
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

    p.add_argument(
        "--roi-source",
        choices=["harmonized", "raw", "auto"],
        default="harmonized",
        help=(
            "harmonized: prefer H_MUSE_Volume_*; raw: use MUSE_Volume_*; "
            "auto: prefer harmonized and otherwise raw."
        ),
    )
    p.add_argument(
        "--allow-raw-roi-fallback",
        action="store_true",
        help=(
            "When --roi-source harmonized, allow raw MUSE_Volume_* values for "
            "model ROIs whose H_MUSE counterpart is unavailable or missing."
        ),
    )
    p.add_argument("--min-roi-fraction", type=float, default=0.80)
    p.add_argument("--complete-case-model-rois", action="store_true")
    p.add_argument(
        "--age-mode",
        choices=["row", "baseline_plus_time"],
        default="row",
    )
    p.add_argument("--risk-times", default="1,2,3,5")
    return p.parse_args()


def log(message: str) -> None:
    print(message, flush=True)


def read_table(path: str) -> pd.DataFrame:
    p = Path(path)
    if p.suffix.lower() == ".csv":
        return pd.read_csv(p, low_memory=False)
    return pd.read_csv(p, sep="\t", low_memory=False)


def parse_list(value: str) -> List[str]:
    return [x.strip() for x in str(value).split(",") if x.strip()]


def parse_risk_times(value: str) -> List[float]:
    return [float(x) for x in parse_list(value)]


def dedupe_preserve_order(values: Sequence[str]) -> List[str]:
    seen = set()
    output: List[str] = []
    for value in values:
        if value not in seen:
            seen.add(value)
            output.append(value)
    return output


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


def normalize_sex(series: pd.Series) -> pd.Series:
    return clean_string(series).replace(
        {
            "0": "Female",
            "0.0": "Female",
            "1": "Male",
            "1.0": "Male",
            "F": "Female",
            "M": "Male",
            "female": "Female",
            "male": "Male",
            "FEMALE": "Female",
            "MALE": "Male",
        }
    )


def normalize_dx(value: object) -> object:
    if pd.isna(value):
        return np.nan
    x = str(value).strip().upper()
    return {
        "NL": "CN",
        "NORMAL": "CN",
        "CONTROL": "CN",
        "HC": "CN",
        "EMCI": "MCI",
        "LMCI": "MCI",
        "DEMENTIA": "AD",
        "CD": "AD",
    }.get(x, x)


def visit_month(value: object) -> float:
    if pd.isna(value):
        return np.nan
    x = str(value).strip().lower()
    if x in {"bl", "base", "baseline", "m00", "m0", "screen", "screening", "sc"}:
        return 0.0
    match = re.search(r"m(?:onth)?\s*0*([0-9]+)", x)
    if match:
        return float(match.group(1))
    match = re.search(r"([0-9]+)", x)
    return float(match.group(1)) if match else np.nan


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
        x
        for x in bundle.get("numeric_cols", [])
        if str(x).startswith("MUSE_Volume_")
    ]
    if not rois:
        raise ValueError("No MUSE ROI list found in the pretrained model bundle.")
    return rois


def harmonized_name(raw_roi_name: str) -> str:
    if not raw_roi_name.startswith("MUSE_Volume_"):
        raise ValueError(f"Unexpected model ROI name: {raw_roi_name}")
    return "H_" + raw_roi_name


def construct_model_roi_matrix(
    raw: pd.DataFrame,
    model_rois: Sequence[str],
    roi_source: str,
    allow_raw_fallback: bool,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Create columns named exactly as the model expects and a source audit table."""
    model_values = pd.DataFrame(index=raw.index)
    audit_rows: List[Dict[str, object]] = []

    for model_roi in model_rois:
        raw_col = model_roi
        harmonized_col = harmonized_name(model_roi)

        raw_exists = raw_col in raw.columns
        harmonized_exists = harmonized_col in raw.columns
        raw_values = clean_numeric(raw[raw_col]) if raw_exists else pd.Series(np.nan, index=raw.index)
        harmonized_values = (
            clean_numeric(raw[harmonized_col])
            if harmonized_exists
            else pd.Series(np.nan, index=raw.index)
        )

        if roi_source == "raw":
            selected = raw_values
            source_label = np.where(selected.notna(), "raw", "missing")
        elif roi_source == "auto":
            selected = harmonized_values.where(harmonized_values.notna(), raw_values)
            source_label = np.where(
                harmonized_values.notna(),
                "harmonized",
                np.where(raw_values.notna(), "raw_fallback", "missing"),
            )
        else:  # harmonized
            if allow_raw_fallback:
                selected = harmonized_values.where(harmonized_values.notna(), raw_values)
                source_label = np.where(
                    harmonized_values.notna(),
                    "harmonized",
                    np.where(raw_values.notna(), "raw_fallback", "missing"),
                )
            else:
                selected = harmonized_values
                source_label = np.where(selected.notna(), "harmonized", "missing")

        model_values[model_roi] = selected
        counts = pd.Series(source_label).value_counts(dropna=False).to_dict()
        audit_rows.append(
            {
                "model_roi": model_roi,
                "raw_column": raw_col,
                "harmonized_column": harmonized_col,
                "raw_column_present": raw_exists,
                "harmonized_column_present": harmonized_exists,
                "n_harmonized_values_used": int(counts.get("harmonized", 0)),
                "n_raw_values_used": int(counts.get("raw", 0)),
                "n_raw_fallback_values_used": int(counts.get("raw_fallback", 0)),
                "n_missing_values": int(counts.get("missing", 0)),
                "requested_roi_source": roi_source,
                "allow_raw_fallback": bool(allow_raw_fallback),
            }
        )

    return model_values, pd.DataFrame(audit_rows)


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

    expected = np.repeat(
        float(info.get("risk_score_covariate_model_intercept", 0.0)),
        len(out),
    )
    categorical_covariates = list(
        info.get("categorical_residualization_covariates", [])
    )

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
    residual_centered = residual_raw - float(
        info.get("risk_score_residual_mean_train", 0.0)
    )
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


def add_longitudinal_metadata(df: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    d = df.copy()
    d["_scan_date"] = pd.to_datetime(d[args.date_col], errors="coerce")
    d["_visit_month"] = d[args.visit_col].apply(visit_month)
    d["_delta_years"] = (
        clean_numeric(d[args.delta_baseline_col])
        if args.delta_baseline_col in d.columns
        else pd.Series(np.nan, index=d.index)
    )

    d["_sort"] = d["_scan_date"].map(
        lambda x: x.toordinal() if pd.notna(x) else np.nan
    )
    d["_sort"] = (
        d["_sort"]
        .fillna(d["_delta_years"] * 365.25)
        .fillna(d["_visit_month"] * 30.4375)
    )
    d = d.sort_values(
        [args.id_col, "_sort"],
        kind="mergesort",
        na_position="last",
    )

    baseline_rows: List[Dict[str, object]] = []
    for participant_id, group in d.groupby(args.id_col, sort=False):
        ordered = group.loc[group["_sort"].notna()]
        baseline = ordered.iloc[0] if not ordered.empty else group.iloc[0]
        baseline_rows.append(
            {
                args.id_col: participant_id,
                "_base_date": baseline["_scan_date"],
                "_base_visit_month": baseline["_visit_month"],
                "_base_age": clean_numeric(pd.Series([baseline[args.age_col]])).iloc[0],
                "external_selected_baseline_visit": baseline[args.visit_col],
                "external_selected_baseline_date": baseline["_scan_date"],
            }
        )
    d = d.merge(pd.DataFrame(baseline_rows), on=args.id_col, how="left")

    years = pd.Series(np.nan, index=d.index, dtype=float)
    date_mask = d["_scan_date"].notna() & d["_base_date"].notna()
    years.loc[date_mask] = (
        d.loc[date_mask, "_scan_date"] - d.loc[date_mask, "_base_date"]
    ).dt.days / 365.25

    delta_mask = years.isna() & d["_delta_years"].notna()
    years.loc[delta_mask] = d.loc[delta_mask, "_delta_years"]

    visit_mask = (
        years.isna()
        & d["_visit_month"].notna()
        & d["_base_visit_month"].notna()
    )
    years.loc[visit_mask] = (
        d.loc[visit_mask, "_visit_month"]
        - d.loc[visit_mask, "_base_visit_month"]
    ) / 12.0

    d["years_since_external_baseline"] = years
    d["longitudinal_scan_number"] = d.groupby(args.id_col, sort=False).cumcount() + 1
    d["is_external_baseline_scan"] = d["longitudinal_scan_number"].eq(1)

    row_age = clean_numeric(d[args.age_col])
    if args.age_mode == "baseline_plus_time":
        derived_age = clean_numeric(d["_base_age"]) + d["years_since_external_baseline"]
        d["Age_original_external"] = row_age
        d["Age"] = derived_age.where(derived_age.notna(), row_age)
    else:
        d["Age"] = row_age
    d["age_at_scan_used_for_model"] = d["Age"]
    return d


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    raw = read_table(args.input_file)
    required = [args.id_col, args.visit_col, args.date_col, args.age_col]
    missing_required = [x for x in required if x not in raw.columns]
    if missing_required:
        raise ValueError(f"Missing required input columns: {missing_required}")

    raw = raw.loc[raw[args.id_col].notna()].copy()
    raw["_source_row_number"] = np.arange(1, len(raw) + 1)

    # Preserve external descriptors under unique names. This avoids duplicate
    # SITE columns when the model-ready SITE field is added later.
    raw["external_SITE"] = (
        clean_string(raw[args.site_col], "__EXTERNAL_UNKNOWN__")
        if args.site_col in raw.columns
        else "__EXTERNAL_UNKNOWN__"
    )
    raw["external_Study"] = (
        clean_string(raw[args.study_col], "__EXTERNAL_UNKNOWN_STUDY__")
        if args.study_col in raw.columns
        else "__EXTERNAL_UNKNOWN_STUDY__"
    )

    if args.study_values:
        requested_studies = set(parse_list(args.study_values))
        raw = raw.loc[raw["external_Study"].isin(requested_studies)].copy()
    if raw.empty:
        raise ValueError("No scans remain after study filtering.")

    bundle = joblib.load(args.model_joblib)
    preprocessor = bundle["preprocessor"]
    model = bundle["model"]
    model_rois = get_model_rois(bundle)

    roi_matrix, roi_audit = construct_model_roi_matrix(
        raw=raw,
        model_rois=model_rois,
        roi_source=args.roi_source,
        allow_raw_fallback=args.allow_raw_roi_fallback,
    )

    for roi in model_rois:
        raw[roi] = roi_matrix[roi]

    raw["model_roi_source_mode"] = args.roi_source
    raw["model_roi_raw_fallback_allowed"] = bool(args.allow_raw_roi_fallback)
    raw["n_model_rois_expected"] = len(model_rois)
    raw["n_model_rois_nonmissing_in_scan"] = roi_matrix.notna().sum(axis=1)
    raw["fraction_model_rois_nonmissing_in_scan"] = (
        raw["n_model_rois_nonmissing_in_scan"] / float(len(model_rois))
    )

    threshold = 1.0 if args.complete_case_model_rois else args.min_roi_fraction
    eligible = raw["fraction_model_rois_nonmissing_in_scan"].ge(threshold)
    excluded = raw.loc[~eligible].copy()
    excluded["exclusion_reason"] = "insufficient_selected_model_ROI_coverage"
    scored = raw.loc[eligible].copy()
    if scored.empty:
        raise ValueError(
            "No scans passed ROI coverage QC. Inspect the ROI source audit and "
            "consider --roi-source auto or --allow-raw-roi-fallback."
        )

    scored = add_longitudinal_metadata(scored, args)
    scored["diagnosis_normalized"] = (
        scored[args.dx_col].apply(normalize_dx)
        if args.dx_col in scored.columns
        else np.nan
    )

    numeric_cols = list(bundle.get("numeric_cols", []))
    categorical_cols = list(bundle.get("categorical_cols", []))
    if not numeric_cols and not categorical_cols:
        raise ValueError("Model bundle lacks numeric_cols/categorical_cols.")

    created_missing: List[str] = []
    for column in numeric_cols + categorical_cols:
        if column not in scored.columns:
            scored[column] = np.nan
            created_missing.append(column)

    if "DLICV" in numeric_cols:
        scored["DLICV"] = (
            clean_numeric(scored[args.dlicv_col])
            if args.dlicv_col in scored.columns
            else np.nan
        )
    if "Sex" in categorical_cols:
        scored["Sex"] = (
            normalize_sex(scored[args.sex_col])
            if args.sex_col in scored.columns
            else np.nan
        )
    if "SITE" in categorical_cols:
        # Separate model input from external_SITE. Every external label is unseen
        # relative to ADNI and is ignored by the saved one-hot encoder.
        scored["SITE"] = scored["external_SITE"].astype("object")

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
        args.id_col,
        args.visit_col,
        args.date_col,
        "external_Study",
        "external_SITE",
        args.dx_col,
        "diagnosis_normalized",
        "longitudinal_scan_number",
        "is_external_baseline_scan",
        "external_selected_baseline_visit",
        "external_selected_baseline_date",
        "years_since_external_baseline",
        "Age",
        "age_at_scan_used_for_model",
        args.sex_col,
        args.dlicv_col,
        "SITE",
        "MRI_ID",
        "MRID",
        "MRI_Magnetic_Field_Strength",
        "MRI_Manufacturer",
        "MRI_Scanner_Model",
        "MRI_Protocol",
        "model_roi_source_mode",
        "model_roi_raw_fallback_allowed",
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

    # Defensive assertion: duplicate names would recreate the SITE grouper bug.
    duplicate_output_names = scan.columns[scan.columns.duplicated()].tolist()
    if duplicate_output_names:
        raise RuntimeError(
            f"Duplicate output columns remain unexpectedly: {duplicate_output_names}"
        )

    scan_file = outdir / f"{args.prefix}_scan_level_predictions.tsv"
    subject_file = outdir / f"{args.prefix}_subject_longitudinal_summary.tsv"
    site_file = outdir / f"{args.prefix}_study_site_summary.tsv"
    excluded_file = outdir / f"{args.prefix}_excluded_scans.tsv"
    roi_file = outdir / f"{args.prefix}_model_roi_source_audit.tsv"

    scan.to_csv(scan_file, sep="\t", index=False)
    excluded.to_csv(excluded_file, sep="\t", index=False)
    roi_audit.to_csv(roi_file, sep="\t", index=False)

    subject = (
        scan.groupby(args.id_col, dropna=False, sort=False)
        .agg(
            n_scans=(args.id_col, "size"),
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

    descriptors = (
        scan.groupby(args.id_col, dropna=False, sort=False)[
            ["external_Study", "external_SITE"]
        ]
        .first()
        .reset_index()
    )
    subject = descriptors.merge(subject, on=args.id_col, how="right")
    subject.to_csv(subject_file, sep="\t", index=False)

    # Group exclusively by unique external descriptor columns. Never group by
    # the model-ready SITE field, which could otherwise collide with site-col.
    grouping_columns = ["external_Study", "external_SITE"]
    site_summary = (
        scan.groupby(grouping_columns, dropna=False, sort=False)
        .agg(
            n_subjects=(args.id_col, "nunique"),
            n_scans=(args.id_col, "size"),
            median_age=("Age", "median"),
            median_roi_coverage=("fraction_model_rois_nonmissing_in_scan", "median"),
            mean_risk=(RISK_COL, "mean"),
            sd_risk=(RISK_COL, "std"),
            mean_acceleration_z=(ACCEL_Z_COL, "mean"),
            sd_acceleration_z=(ACCEL_Z_COL, "std"),
        )
        .reset_index()
    )
    site_summary.to_csv(site_file, sep="\t", index=False)

    summary = {
        "input_file": args.input_file,
        "model_joblib": args.model_joblib,
        "n_input_scans": int(len(raw)),
        "n_scored_scans": int(len(scan)),
        "n_excluded_scans": int(len(excluded)),
        "n_scored_subjects": int(scan[args.id_col].nunique()),
        "n_model_rois_expected": int(len(model_rois)),
        "requested_roi_source": args.roi_source,
        "allow_raw_roi_fallback": bool(args.allow_raw_roi_fallback),
        "minimum_roi_fraction": float(threshold),
        "created_missing_expected_columns": created_missing,
        "model_refit": False,
        "diagnosis_rule": "All eligible scans are scored regardless of diagnosis.",
        "site_rule": (
            "External site is preserved as external_SITE. The model-ready SITE "
            "contains the same external label and is treated as unseen by the "
            "saved ADNI OneHotEncoder(handle_unknown='ignore')."
        ),
        "feature_rule": (
            "Model inputs retain the raw MUSE_Volume_* names expected by the "
            "saved ADNI preprocessor, while values are sourced according to "
            "--roi-source and documented in the ROI source audit."
        ),
    }
    with open(outdir / f"{args.prefix}_application_summary.json", "w") as handle:
        json.dump(summary, handle, indent=2)
    pd.DataFrame([summary]).to_csv(
        outdir / f"{args.prefix}_application_summary.tsv",
        sep="\t",
        index=False,
    )

    log(f"Scored {len(scan):,} scans from {scan[args.id_col].nunique():,} participants.")
    log(f"ROI source mode: {args.roi_source}")
    log(f"Predictions: {scan_file}")
    log(f"Study/site summary: {site_file}")
    log(f"ROI source audit: {roi_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
