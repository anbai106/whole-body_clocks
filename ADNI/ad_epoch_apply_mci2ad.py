#!/usr/bin/env python3
"""
Apply a pretrained ADNI brain MRI AD L'EPOCH model to baseline-MCI participants
who were not used to fit the original model, and construct an MCI-to-AD
conversion survival-analysis dataset.

Primary rules
-------------
1. Load the frozen preprocessor, Cox model, and clock-transform information.
2. Exclude participants used to fit the original model. By default, IDs with
   split == train or validation in the training prediction file are excluded.
   If the exclusion file has no split column, every listed ID is excluded.
3. Among remaining participants, define time zero as the first MCI-labeled visit
   with sufficient MUSE GM ROI coverage.
4. Exclude participants with AD at or before the selected MCI baseline.
5. Event = first AD diagnosis after baseline MCI.
6. Censoring = last available follow-up visit without AD.
7. Score only the selected baseline MCI scan, yielding one row per participant.

Primary output
--------------
<prefix>_mci2ad_epoch_survival.tsv
"""

from __future__ import annotations

import argparse
import json
import re
import warnings
from pathlib import Path
from typing import Iterable

import joblib
import numpy as np
import pandas as pd


RISK_COL = "adni_brain_mri_ad_lepoch_risk_score"
ACCEL_Z_COL = "adni_brain_mri_ad_lepoch_acceleration_z"
ACCEL_YEARS_COL = "adni_brain_mri_ad_lepoch_acceleration_years"
CLOCK_AGE_COL = "adni_brain_mri_ad_lepoch_clock_age_years"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Apply a frozen ADNI AD L'EPOCH model to non-training baseline-MCI "
            "participants and build an MCI-to-AD survival dataset."
        )
    )
    parser.add_argument("--input-file", required=True)
    parser.add_argument("--model-joblib", required=True)
    parser.add_argument("--training-participants-file", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--prefix", default="adni_brain_mri_ad_lepoch")

    parser.add_argument("--id-col", default="PTID")
    parser.add_argument("--visit-col", default="Visit_Code")
    parser.add_argument("--date-col", default="Date")
    parser.add_argument("--dx-col", default="DX_Binary")
    parser.add_argument("--training-id-col", default="PTID")
    parser.add_argument("--training-split-col", default="split")
    parser.add_argument(
        "--exclude-splits",
        default="train,validation",
        help=(
            "Comma-separated model-development splits to exclude. If the split "
            "column is absent, all IDs in the training-participants file are excluded."
        ),
    )

    parser.add_argument("--baseline-dx", default="MCI")
    parser.add_argument("--event-dx", default="AD")
    parser.add_argument("--min-baseline-roi-fraction", type=float, default=0.80)
    parser.add_argument("--min-followup-days", type=int, default=1)
    parser.add_argument("--risk-times", default="1,2,3,5")
    parser.add_argument(
        "--complete-case-model-features",
        action="store_true",
        help="Require all model MUSE ROI features to be observed at baseline.",
    )
    return parser.parse_args()


def log(message: str) -> None:
    print(message, flush=True)


def read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path, low_memory=False)
    return pd.read_csv(path, sep="\t", low_memory=False)


def parse_list_arg(value: str | None) -> list[str]:
    if value is None or str(value).strip() == "":
        return []
    return [item.strip() for item in str(value).split(",") if item.strip()]


def parse_risk_times(value: str) -> list[float]:
    return [float(x.strip()) for x in value.split(",") if x.strip()]


def normalize_dx(value: object) -> str | float:
    if pd.isna(value):
        return np.nan
    text = str(value).strip().upper()
    mapping = {
        "CN": "CN", "NL": "CN", "NORMAL": "CN", "CONTROL": "CN", "HC": "CN",
        "MCI": "MCI", "EMCI": "MCI", "LMCI": "MCI",
        "AD": "AD", "DEMENTIA": "AD",
    }
    return mapping.get(text, text)


def parse_date_series(series: pd.Series) -> pd.Series:
    return pd.to_datetime(series, errors="coerce")


def visit_code_to_month(value: object) -> float:
    if pd.isna(value):
        return np.nan
    text = str(value).strip().lower()
    if text in {"bl", "base", "baseline", "m00", "m0", "screen", "screening", "sc"}:
        return 0.0
    match = re.search(r"m(?:onth)?\s*0*([0-9]+)", text)
    if match:
        return float(match.group(1))
    match = re.search(r"([0-9]+)", text)
    return float(match.group(1)) if match else np.nan


def clean_numeric_series(series: pd.Series) -> pd.Series:
    if pd.api.types.is_numeric_dtype(series):
        return pd.to_numeric(series, errors="coerce")
    cleaned = (
        series.astype("object")
        .where(series.notna(), np.nan)
        .astype(str)
        .str.strip()
        .str.replace(",", "", regex=False)
        .replace({
            "": np.nan, "nan": np.nan, "NaN": np.nan, "NA": np.nan,
            "N/A": np.nan, "None": np.nan, "null": np.nan,
        })
    )
    return pd.to_numeric(cleaned, errors="coerce")


def normalize_sex_series(series: pd.Series) -> pd.Series:
    return (
        series.astype(str).str.strip().replace({
            "0": "Female", "0.0": "Female", "1": "Male", "1.0": "Male",
            "F": "Female", "M": "Male", "female": "Female", "male": "Male",
            "Female": "Female", "Male": "Male",
        })
    )


def get_model_roi_features(bundle: dict) -> list[str]:
    for key in ("selected_muse_gm_rois", "available_muse_gm_rois", "hardcoded_muse_gm_rois"):
        values = bundle.get(key)
        if values:
            return list(values)
    numeric_cols = list(bundle.get("numeric_cols", []))
    roi_cols = [col for col in numeric_cols if str(col).startswith("MUSE_Volume_")]
    if not roi_cols:
        raise ValueError("Could not identify MUSE ROI features from the model bundle.")
    return roi_cols


def ensure_expected_columns(df: pd.DataFrame, expected_cols: Iterable[str]) -> pd.DataFrame:
    out = df.copy()
    for col in expected_cols:
        if col not in out.columns:
            warnings.warn(f"Model expected column {col} is missing; creating it as NA.")
            out[col] = np.nan
    return out


def compute_row_roi_coverage(df: pd.DataFrame, roi_cols: list[str]) -> tuple[pd.Series, pd.Series]:
    numeric = df[roi_cols].apply(clean_numeric_series, axis=0)
    n_nonmissing = numeric.notna().sum(axis=1)
    fraction = n_nonmissing / float(len(roi_cols))
    return n_nonmissing, fraction


def load_excluded_training_ids(
    path: str | Path,
    id_col: str,
    split_col: str,
    exclude_splits: list[str],
) -> tuple[set[str], pd.DataFrame, str]:
    table = read_table(path)
    if id_col not in table.columns:
        raise ValueError(
            f"Training participant file lacks ID column {id_col}. "
            f"Observed columns: {list(table.columns)}"
        )

    table = table.loc[table[id_col].notna()].copy()
    table[id_col] = table[id_col].astype(str).str.strip()

    if split_col in table.columns and exclude_splits:
        normalized_splits = {x.strip().lower() for x in exclude_splits}
        split_values = table[split_col].astype(str).str.strip().str.lower()
        excluded_rows = table.loc[split_values.isin(normalized_splits)].copy()
        rule = f"Excluded IDs with {split_col} in {sorted(normalized_splits)}"
    else:
        excluded_rows = table.copy()
        rule = "Split column unavailable or no splits requested; excluded every listed ID"

    excluded_ids = set(excluded_rows[id_col].dropna().astype(str).str.strip())
    return excluded_ids, excluded_rows, rule


def construct_mci2ad_baseline_dataset(
    df_raw: pd.DataFrame,
    args: argparse.Namespace,
    model_roi_cols: list[str],
    excluded_training_ids: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[str], list[str]]:
    d = df_raw.copy()
    required = [args.id_col, args.visit_col, args.date_col, args.dx_col]
    missing_required = [col for col in required if col not in d.columns]
    if missing_required:
        raise ValueError(f"Missing required input columns: {missing_required}")

    d[args.id_col] = d[args.id_col].astype("object")
    d = d.loc[d[args.id_col].notna()].copy()
    d[args.id_col] = d[args.id_col].astype(str).str.strip()

    excluded_rows = d.loc[d[args.id_col].isin(excluded_training_ids)].copy()
    d = d.loc[~d[args.id_col].isin(excluded_training_ids)].copy()

    missing_roi = [col for col in model_roi_cols if col not in d.columns]
    present_roi = [col for col in model_roi_cols if col in d.columns]
    if not present_roi:
        raise ValueError("None of the model MUSE ROI columns are present in the ADNI input file.")
    if missing_roi:
        warnings.warn(
            f"{len(missing_roi)} model ROI columns are absent from the application file. "
            f"First missing: {missing_roi[:10]}"
        )

    d["_dx_norm"] = d[args.dx_col].apply(normalize_dx)
    d["_date"] = parse_date_series(d[args.date_col])
    d["_visit_month"] = d[args.visit_col].apply(visit_code_to_month)
    d["_date_sort"] = d["_date"].map(lambda x: x.toordinal() if pd.notna(x) else np.nan)
    d["_visit_sort"] = d["_visit_month"] * 30.4375
    d["_sort_key"] = d["_date_sort"].fillna(d["_visit_sort"])
    d = d.loc[d["_sort_key"].notna()].copy()

    d["_muse_roi_nonmissing_n"], d["_muse_roi_nonmissing_fraction"] = compute_row_roi_coverage(
        d, present_roi
    )
    d["_has_usable_muse"] = (
        d["_muse_roi_nonmissing_fraction"] >= args.min_baseline_roi_fraction
    )
    d = d.sort_values([args.id_col, "_sort_key"], kind="mergesort")

    baseline_dx = normalize_dx(args.baseline_dx)
    event_dx = normalize_dx(args.event_dx)
    baseline_rows: list[pd.Series] = []
    skipped_rows: list[dict] = []

    for pid, group in d.groupby(args.id_col, sort=False):
        g = group.sort_values("_sort_key", kind="mergesort").copy()
        candidates = g.loc[
            g["_has_usable_muse"] & (g["_dx_norm"] == baseline_dx)
        ].copy()

        if candidates.empty:
            skipped_rows.append({
                args.id_col: pid,
                "skip_reason": "no_MCI_visit_with_sufficient_MUSE_ROI_coverage",
                "n_rows": int(g.shape[0]),
                "max_muse_roi_nonmissing_fraction": float(g["_muse_roi_nonmissing_fraction"].max()),
            })
            continue

        # First available MCI diagnosis with usable MUSE data; no preference for Visit_Code == bl.
        baseline = candidates.iloc[0]
        baseline_sort = baseline["_sort_key"]
        baseline_date = baseline["_date"]
        baseline_month = baseline["_visit_month"]

        ad_at_or_before_baseline = g.loc[
            (g["_sort_key"] <= baseline_sort) & (g["_dx_norm"] == event_dx)
        ]
        if not ad_at_or_before_baseline.empty:
            skipped_rows.append({
                args.id_col: pid,
                "skip_reason": "AD_at_or_before_selected_MCI_baseline",
                "n_rows": int(g.shape[0]),
                "max_muse_roi_nonmissing_fraction": float(g["_muse_roi_nonmissing_fraction"].max()),
            })
            continue

        follow = g.loc[g["_sort_key"] > baseline_sort].copy()
        if follow.empty:
            skipped_rows.append({
                args.id_col: pid,
                "skip_reason": "no_followup_after_selected_MCI_baseline",
                "n_rows": int(g.shape[0]),
                "max_muse_roi_nonmissing_fraction": float(g["_muse_roi_nonmissing_fraction"].max()),
            })
            continue

        event_rows = follow.loc[follow["_dx_norm"] == event_dx].copy()
        if not event_rows.empty:
            end_row = event_rows.iloc[0]
            event = True
        else:
            end_row = follow.iloc[-1]
            event = False

        if pd.notna(baseline_date) and pd.notna(end_row["_date"]):
            time_days = float((end_row["_date"] - baseline_date).days)
        elif pd.notna(baseline_month) and pd.notna(end_row["_visit_month"]):
            time_days = float((end_row["_visit_month"] - baseline_month) * 30.4375)
        else:
            time_days = np.nan

        if not np.isfinite(time_days) or time_days < args.min_followup_days:
            skipped_rows.append({
                args.id_col: pid,
                "skip_reason": "invalid_or_too_short_followup_after_selected_MCI_baseline",
                "n_rows": int(g.shape[0]),
                "max_muse_roi_nonmissing_fraction": float(g["_muse_roi_nonmissing_fraction"].max()),
            })
            continue

        row = baseline.copy()
        row["selected_baseline_visit_code"] = baseline[args.visit_col]
        row["selected_baseline_date"] = baseline_date
        row["selected_baseline_dx"] = baseline_dx
        row["baseline_muse_roi_nonmissing_n"] = int(baseline["_muse_roi_nonmissing_n"])
        row["baseline_muse_roi_nonmissing_fraction"] = float(
            baseline["_muse_roi_nonmissing_fraction"]
        )
        row["event"] = bool(event)
        row["event_dx"] = event_dx if event else np.nan
        row["event_or_censor_dx"] = end_row["_dx_norm"]
        row["event_or_censor_visit_code"] = end_row[args.visit_col]
        row["event_or_censor_date"] = end_row["_date"]
        row["time_days"] = time_days
        row["time_years"] = time_days / 365.25
        row["conversion_group"] = "MCI to AD" if event else "MCI censored without AD"
        row["excluded_from_original_model_fit"] = False
        row["n_model_muse_rois_expected"] = len(model_roi_cols)
        row["n_model_muse_rois_present_in_file"] = len(present_roi)
        row["n_model_muse_rois_missing_from_file"] = len(missing_roi)
        baseline_rows.append(row)

    baseline_df = pd.DataFrame(baseline_rows)
    skipped_df = pd.DataFrame(skipped_rows)
    excluded_subject_df = (
        excluded_rows[[args.id_col]].drop_duplicates().assign(
            exclusion_reason="participant_used_in_original_model_fit"
        )
        if not excluded_rows.empty
        else pd.DataFrame(columns=[args.id_col, "exclusion_reason"])
    )
    return baseline_df, skipped_df, excluded_subject_df, present_roi, missing_roi


def categorical_match(series: pd.Series, category_value: object) -> pd.Series:
    obj = series.astype("object")
    text = obj.astype(str).str.strip()
    category_text = str(category_value).strip()
    mask = text == category_text
    try:
        category_float = float(category_text)
        numeric = pd.to_numeric(obj, errors="coerce")
        mask = mask | np.isclose(numeric.astype(float), category_float, equal_nan=False)
    except Exception:
        pass
    return mask.astype(float)


def parse_categorical_coef_name(term: str, categorical_covs: list[str]) -> tuple[str | None, str | None]:
    if not term.startswith("cat__"):
        return None, None
    stem = term[len("cat__"):]
    for covariate in sorted(categorical_covs, key=len, reverse=True):
        prefix = f"{covariate}_"
        if stem.startswith(prefix):
            return covariate, stem[len(prefix):]
    return None, None


def numeric_fill_for_residualization(series: pd.Series, name: str) -> pd.Series:
    numeric = clean_numeric_series(series)
    if numeric.isna().all():
        warnings.warn(f"Residualization covariate {name} is all missing; using zero.")
        return pd.Series(np.zeros(len(numeric)), index=numeric.index)
    return numeric.fillna(float(np.nanmedian(numeric)))


def compute_clock_transforms_from_saved_info(
    df: pd.DataFrame,
    risk: np.ndarray,
    clock_info: dict | None,
) -> pd.DataFrame:
    out = df.copy()
    out[RISK_COL] = risk
    if not clock_info:
        warnings.warn("clock_transform_info is missing; acceleration columns set to NA.")
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
        for col in covariates:
            if col in out.columns and col not in {"Sex", "SITE"}:
                numeric_covs.append(col)
            else:
                categorical_covs.append(col)

    expected = np.repeat(intercept, out.shape[0]).astype(float)
    for term, beta_value in coef_dict.items():
        beta = float(beta_value)
        if term.startswith("num__"):
            covariate = term[len("num__"):]
            values = (
                numeric_fill_for_residualization(out[covariate], covariate).values.astype(float)
                if covariate in out.columns
                else np.zeros(out.shape[0])
            )
            expected += beta * values
        elif term.startswith("cat__"):
            covariate, category = parse_categorical_coef_name(term, categorical_covs)
            if covariate is None or covariate not in out.columns:
                continue
            expected += beta * categorical_match(out[covariate], category).values.astype(float)

    residual_mean = float(clock_info.get("risk_score_residual_mean_train", 0.0))
    residual_sd = float(clock_info.get("risk_score_residual_sd_train", np.nan))
    beta_age = clock_info.get("adjusted_age_coefficient_risk_score_per_year")
    residual = risk - expected - residual_mean

    out[ACCEL_Z_COL] = residual / residual_sd if np.isfinite(residual_sd) and residual_sd > 0 else np.nan
    if beta_age is not None and np.isfinite(float(beta_age)) and abs(float(beta_age)) > 1e-8:
        beta_age = float(beta_age)
        out[ACCEL_YEARS_COL] = residual / beta_age
        out[CLOCK_AGE_COL] = (
            clean_numeric_series(out["Age"]) + out[ACCEL_YEARS_COL]
            if "Age" in out.columns
            else np.nan
        )
    else:
        out[ACCEL_YEARS_COL] = np.nan
        out[CLOCK_AGE_COL] = np.nan
    return out


def predict_absolute_risk(model, X, times_years: list[float]) -> pd.DataFrame:
    output: dict[str, list[float] | np.ndarray] = {}
    try:
        survival_functions = model.predict_survival_function(X)
    except Exception as exc:
        warnings.warn(f"Could not compute absolute risks: {exc}")
        for time in times_years:
            output[f"risk_{time:g}y"] = np.repeat(np.nan, X.shape[0])
        return pd.DataFrame(output)

    for time in times_years:
        values = []
        for function in survival_functions:
            try:
                values.append(1.0 - float(function(time)))
            except Exception:
                values.append(np.nan)
        output[f"risk_{time:g}y"] = values
    return pd.DataFrame(output)


def prepare_model_input(
    baseline_df: pd.DataFrame,
    bundle: dict,
) -> tuple[pd.DataFrame, list[str], list[str]]:
    df = baseline_df.copy()
    numeric_cols = list(bundle.get("numeric_cols", []))
    categorical_cols = list(bundle.get("categorical_cols", []))
    if not numeric_cols and not categorical_cols:
        raise ValueError("Model bundle lacks numeric_cols and categorical_cols.")

    df = ensure_expected_columns(df, numeric_cols + categorical_cols)
    for col in numeric_cols:
        df[col] = clean_numeric_series(df[col])
    if "Sex" in categorical_cols and "Sex" in df.columns:
        df["Sex"] = normalize_sex_series(df["Sex"])
    if "SITE" in categorical_cols and "SITE" in df.columns:
        df["SITE"] = df["SITE"].astype(str).str.strip()
    for col in categorical_cols:
        df[col] = df[col].astype("object")
    return df, numeric_cols, categorical_cols


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    risk_times = parse_risk_times(args.risk_times)

    log("============================================================")
    log("Apply frozen AD L'EPOCH to non-training baseline-MCI participants")
    log("Outcome: MCI-to-AD conversion")
    log("============================================================")

    bundle = joblib.load(args.model_joblib)
    preprocessor = bundle["preprocessor"]
    model = bundle["model"]
    clock_info = bundle.get("clock_transform_info")
    model_roi_cols = get_model_roi_features(bundle)

    excluded_ids, excluded_training_rows, exclusion_rule = load_excluded_training_ids(
        args.training_participants_file,
        args.training_id_col,
        args.training_split_col,
        parse_list_arg(args.exclude_splits),
    )
    log(f"Training-exclusion rule: {exclusion_rule}")
    log(f"Unique participant IDs excluded before application: {len(excluded_ids):,}")

    df_raw = read_table(args.input_file)
    baseline_df, skipped_df, excluded_subject_df, roi_present, roi_missing = (
        construct_mci2ad_baseline_dataset(
            df_raw=df_raw,
            args=args,
            model_roi_cols=model_roi_cols,
            excluded_training_ids=excluded_ids,
        )
    )

    excluded_training_rows.to_csv(
        outdir / f"{args.prefix}_mci2ad_training_rows_used_for_exclusion.tsv",
        sep="\t", index=False,
    )
    excluded_subject_df.rename(columns={args.id_col: "PTID"}).to_csv(
        outdir / f"{args.prefix}_mci2ad_excluded_training_participants.tsv",
        sep="\t", index=False,
    )
    skipped_df.rename(columns={args.id_col: "PTID"}).to_csv(
        outdir / f"{args.prefix}_mci2ad_skipped_participants.tsv",
        sep="\t", index=False,
    )

    if baseline_df.empty:
        raise ValueError("No eligible non-training baseline-MCI participants remained.")

    if args.complete_case_model_features:
        before = baseline_df.shape[0]
        baseline_df = baseline_df.loc[
            baseline_df["baseline_muse_roi_nonmissing_fraction"] >= 1.0
        ].copy()
        log(f"Complete-case ROI filtering: {before} -> {baseline_df.shape[0]}")
    if baseline_df.empty:
        raise ValueError("No participants remain after complete-case ROI filtering.")

    model_df, numeric_cols, categorical_cols = prepare_model_input(baseline_df, bundle)
    X_raw = model_df[numeric_cols + categorical_cols].copy()
    X = preprocessor.transform(X_raw)
    risk = np.asarray(model.predict(X)).reshape(-1)

    predictions = compute_clock_transforms_from_saved_info(model_df, risk, clock_info)
    absolute_risk = predict_absolute_risk(model, X, risk_times)
    predictions = pd.concat(
        [predictions.reset_index(drop=True), absolute_risk.reset_index(drop=True)],
        axis=1,
    )

    preferred_cols = [
        args.id_col, args.visit_col, args.date_col, args.dx_col,
        "selected_baseline_visit_code", "selected_baseline_date", "selected_baseline_dx",
        "baseline_muse_roi_nonmissing_n", "baseline_muse_roi_nonmissing_fraction",
        "event", "event_dx", "event_or_censor_dx", "event_or_censor_visit_code",
        "event_or_censor_date", "time_days", "time_years", "conversion_group",
        "Age", "Sex", "DLICV", "SITE",
        RISK_COL, ACCEL_Z_COL, ACCEL_YEARS_COL, CLOCK_AGE_COL,
    ] + [f"risk_{time:g}y" for time in risk_times] + [
        "n_model_muse_rois_expected", "n_model_muse_rois_present_in_file",
        "n_model_muse_rois_missing_from_file",
    ]
    output_cols = [col for col in preferred_cols if col in predictions.columns]
    output = predictions[output_cols].copy().rename(columns={
        args.id_col: "PTID",
        args.visit_col: "Visit_Code",
        args.date_col: "Date",
        args.dx_col: "DX_Binary",
    })

    if output["PTID"].duplicated().any():
        duplicates = output.loc[output["PTID"].duplicated(), "PTID"].unique().tolist()
        raise RuntimeError(f"Output is not one row per participant. Duplicate IDs: {duplicates[:10]}")

    output_path = outdir / f"{args.prefix}_mci2ad_epoch_survival.tsv"
    output.to_csv(output_path, sep="\t", index=False)

    summary_table = output.groupby("conversion_group", dropna=False).agg(
        n_participants=("PTID", "nunique"),
        n_events=("event", "sum"),
        median_followup_years=("time_years", "median"),
        median_epoch_risk=(RISK_COL, "median"),
    ).reset_index()
    summary_table.to_csv(
        outdir / f"{args.prefix}_mci2ad_epoch_summary.tsv",
        sep="\t", index=False,
    )

    metadata = {
        "input_file": str(args.input_file),
        "model_joblib": str(args.model_joblib),
        "training_participants_file": str(args.training_participants_file),
        "training_exclusion_rule": exclusion_rule,
        "exclude_splits": parse_list_arg(args.exclude_splits),
        "n_unique_training_ids_excluded": int(len(excluded_ids)),
        "n_eligible_mci_baseline_participants": int(output.shape[0]),
        "n_mci_to_ad_events": int(output["event"].sum()),
        "n_censored_without_ad": int((~output["event"].astype(bool)).sum()),
        "median_followup_years": float(output["time_years"].median()),
        "baseline_definition": "First MCI diagnosis with sufficient MUSE data after excluding model-fit participants",
        "event_definition": "First AD diagnosis after selected MCI baseline",
        "censoring_definition": "Last available follow-up visit without AD",
        "n_model_roi_features_expected": int(len(model_roi_cols)),
        "n_model_roi_features_present": int(len(roi_present)),
        "n_model_roi_features_missing": int(len(roi_missing)),
        "missing_model_roi_features": roi_missing,
        "primary_output": str(output_path),
        "primary_score": RISK_COL,
        "recommended_adjusted_score": ACCEL_Z_COL,
    }
    with open(outdir / f"{args.prefix}_mci2ad_application_summary.json", "w") as handle:
        json.dump(metadata, handle, indent=2)
    pd.DataFrame([metadata]).to_csv(
        outdir / f"{args.prefix}_mci2ad_application_summary.tsv",
        sep="\t", index=False,
    )

    log("============================================================")
    log("MCI-to-AD AD L'EPOCH application complete")
    log(f"Primary survival TSV: {output_path}")
    log(summary_table.to_string(index=False))
    log("============================================================")


if __name__ == "__main__":
    main()
