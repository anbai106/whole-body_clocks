#!/usr/bin/env python3
"""
Compare baseline AD EPOCH with clinical, CSF, and MRI-derived biomarkers for
MCI-to-AD conversion in ADNI.

Primary cohort
--------------
One row per participant from the previously generated MCI-to-AD EPOCH survival
file. Participants used to fit the original AD EPOCH model should already have
been excluded by the application script.

Time-to-event outcome
---------------------
Time zero: first qualifying MCI MRI with usable MUSE data.
Event: first subsequent AD diagnosis.
Censoring: last follow-up without AD.

Predictors
----------
Clinical covariates (no baseline cognition):
    Age, Sex, Education_Years, APOE_Genotype

Primary AD EPOCH predictor:
    adni_brain_mri_ad_lepoch_risk_score

Comparison biomarkers:
    Abeta_CSF, Tau_CSF, PTau_CSF, SPARE_AD, SPARE_BA

Baseline matching
-----------------
* Clinical covariates and SPARE scores are taken from the exact selected
  baseline-MCI MRI row whenever possible.
* CSF values are first taken from the exact baseline row. If absent, the latest
  measurement on or before the baseline date within --csf-lookback-days is used.
  Post-baseline CSF is not used, preventing future-information leakage.

Analysis
--------
For each biomarker comparison set, all candidate models are evaluated on the
same participant subset. The script reports:
* full-sample Cox hazard ratios, 95% confidence intervals, and P values;
* stratified K-fold out-of-fold Harrell C-index;
* paired bootstrap confidence intervals for delta C-index;
* likelihood-ratio tests for adding AD EPOCH to nested biomarker models.

Required packages
-----------------
pandas, numpy, scipy, scikit-learn, lifelines
"""

from __future__ import annotations

import argparse
import json
import re
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from scipy.stats import chi2
from sklearn.model_selection import StratifiedKFold


DEFAULT_EPOCH_COL = "adni_brain_mri_ad_lepoch_risk_score"
CLINICAL_COLS = ["Age", "Sex", "Education_Years", "APOE_Genotype"]
CSF_COLS = ["Abeta_CSF", "Tau_CSF", "PTau_CSF"]
SPARE_COLS = ["SPARE_AD", "SPARE_BA"]


@dataclass(frozen=True)
class ModelSpec:
    name: str
    predictors: tuple[str, ...]


@dataclass(frozen=True)
class AnalysisGroup:
    name: str
    biomarkers: tuple[str, ...]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare baseline AD EPOCH with clinical, CSF, SPARE-AD, and "
            "SPARE-BA predictors of MCI-to-AD conversion."
        )
    )
    parser.add_argument("--survival-file", required=True)
    parser.add_argument("--adni-file", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--prefix", default="adni_mci2ad_epoch_biomarkers")

    parser.add_argument("--id-col", default="PTID")
    parser.add_argument("--visit-col", default="Visit_Code")
    parser.add_argument("--date-col", default="Date")
    parser.add_argument("--baseline-visit-col", default="selected_baseline_visit_code")
    parser.add_argument("--baseline-date-col", default="selected_baseline_date")
    parser.add_argument("--time-col", default="time_years")
    parser.add_argument("--event-col", default="event")
    parser.add_argument("--epoch-col", default=DEFAULT_EPOCH_COL)

    parser.add_argument("--clinical-cols", default=",".join(CLINICAL_COLS))
    parser.add_argument("--csf-cols", default=",".join(CSF_COLS))
    parser.add_argument("--spare-cols", default=",".join(SPARE_COLS))
    parser.add_argument("--csf-lookback-days", type=int, default=365)

    parser.add_argument("--cv-folds", type=int, default=5)
    parser.add_argument("--random-state", type=int, default=20260730)
    parser.add_argument("--cox-penalizer", type=float, default=0.01)
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--minimum-n", type=int, default=50)
    parser.add_argument("--minimum-events", type=int, default=20)
    return parser.parse_args()


def log(message: str) -> None:
    print(message, flush=True)


def parse_list(value: str | None) -> list[str]:
    if value is None or str(value).strip() == "":
        return []
    return [item.strip() for item in str(value).split(",") if item.strip()]


def read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty input file: {path}")
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path, low_memory=False)
    return pd.read_csv(path, sep="\t", low_memory=False)


def clean_id(series: pd.Series) -> pd.Series:
    return series.astype("object").where(series.notna(), np.nan).astype(str).str.strip()


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
                "<200": np.nan,
                ">1700": np.nan,
            }
        )
    )
    # Extract a numeric value from strings that contain units or inequality text.
    extracted = cleaned.astype(str).str.extract(r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", expand=False)
    extracted = extracted.where(cleaned.notna(), np.nan)
    return pd.to_numeric(extracted, errors="coerce")


def normalize_sex(series: pd.Series) -> pd.Series:
    return (
        series.astype("object")
        .where(series.notna(), np.nan)
        .astype(str)
        .str.strip()
        .replace(
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
    )


def normalize_apoe(series: pd.Series) -> pd.Series:
    def one(value: object) -> object:
        if pd.isna(value):
            return np.nan
        text = str(value).strip().upper()
        text = text.replace("APOE", "").replace("Ε", "E").replace("EPSILON", "")
        text = re.sub(r"[^234]", "", text)
        if len(text) >= 2:
            alleles = sorted([text[0], text[1]])
            return f"e{alleles[0]}/e{alleles[1]}"
        original = str(value).strip()
        return original if original else np.nan

    return series.map(one)


def parse_bool_event(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.astype(int)
    numeric = pd.to_numeric(series, errors="coerce")
    result = numeric.copy()
    missing = result.isna()
    text = series.astype(str).str.strip().str.lower()
    mapping = {
        "true": 1,
        "t": 1,
        "yes": 1,
        "y": 1,
        "event": 1,
        "false": 0,
        "f": 0,
        "no": 0,
        "n": 0,
        "censored": 0,
    }
    result.loc[missing] = text.loc[missing].map(mapping)
    return result.astype("Int64")


def first_nonmissing(values: pd.Series) -> object:
    observed = values.dropna()
    return observed.iloc[0] if not observed.empty else np.nan


def select_exact_baseline_rows(
    survival: pd.DataFrame,
    adni: pd.DataFrame,
    args: argparse.Namespace,
    requested_cols: Sequence[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return one exact or nearest same-day baseline MRI row per participant."""
    s = survival[[args.id_col, args.baseline_visit_col, args.baseline_date_col]].copy()
    s[args.id_col] = clean_id(s[args.id_col])
    s["_baseline_date"] = pd.to_datetime(s[args.baseline_date_col], errors="coerce")
    s["_baseline_visit"] = s[args.baseline_visit_col].astype(str).str.strip().str.lower()

    a = adni.copy()
    a[args.id_col] = clean_id(a[args.id_col])
    a["_row_date"] = pd.to_datetime(a[args.date_col], errors="coerce")
    a["_row_visit"] = a[args.visit_col].astype(str).str.strip().str.lower()

    available = [col for col in requested_cols if col in a.columns]
    missing = sorted(set(requested_cols) - set(available))
    if missing:
        warnings.warn(f"Requested baseline columns absent from ADNI file: {missing}")

    rows: list[dict] = []
    qc: list[dict] = []
    grouped = {pid: group.copy() for pid, group in a.groupby(args.id_col, sort=False)}

    for rec in s.itertuples(index=False):
        pid = str(getattr(rec, args.id_col)).strip()
        baseline_visit = str(getattr(rec, args.baseline_visit_col)).strip().lower()
        baseline_date = getattr(rec, "_baseline_date")
        g = grouped.get(pid)

        selected = None
        match_method = "unmatched"
        date_difference = np.nan

        if g is not None and not g.empty:
            by_visit = g.loc[g["_row_visit"] == baseline_visit]
            if not by_visit.empty:
                if pd.notna(baseline_date) and by_visit["_row_date"].notna().any():
                    differences = (by_visit["_row_date"] - baseline_date).abs().dt.days
                    selected = by_visit.loc[differences.idxmin()]
                    date_difference = float(differences.loc[selected.name])
                else:
                    selected = by_visit.iloc[0]
                match_method = "PTID+Visit_Code"
            elif pd.notna(baseline_date) and g["_row_date"].notna().any():
                differences = (g["_row_date"] - baseline_date).abs().dt.days
                minimum_index = differences.idxmin()
                if float(differences.loc[minimum_index]) == 0:
                    selected = g.loc[minimum_index]
                    date_difference = 0.0
                    match_method = "PTID+exact_Date"

        output = {args.id_col: pid}
        for col in requested_cols:
            output[col] = selected[col] if selected is not None and col in selected.index else np.nan
        rows.append(output)
        qc.append(
            {
                args.id_col: pid,
                "baseline_match_method": match_method,
                "baseline_match_date_difference_days": date_difference,
                "baseline_row_found": selected is not None,
            }
        )

    return pd.DataFrame(rows), pd.DataFrame(qc)


def supplement_subject_level_clinical(
    merged: pd.DataFrame,
    adni: pd.DataFrame,
    args: argparse.Namespace,
    clinical_cols: Sequence[str],
) -> pd.DataFrame:
    """Fill stable clinical fields from any participant row when baseline is missing."""
    out = merged.copy()
    stable = [col for col in clinical_cols if col in adni.columns]
    if not stable:
        return out

    source = adni[[args.id_col] + stable].copy()
    source[args.id_col] = clean_id(source[args.id_col])
    source = source.groupby(args.id_col, as_index=False).agg({col: first_nonmissing for col in stable})
    source = source.rename(columns={col: f"_subject_{col}" for col in stable})
    out = out.merge(source, on=args.id_col, how="left", validate="one_to_one")

    for col in stable:
        backup = f"_subject_{col}"
        if col not in out.columns:
            out[col] = out[backup]
        else:
            out[col] = out[col].where(out[col].notna(), out[backup])
        out = out.drop(columns=backup)
    return out


def attach_prebaseline_csf(
    merged: pd.DataFrame,
    adni: pd.DataFrame,
    args: argparse.Namespace,
    csf_cols: Sequence[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fill missing CSF from latest observation at/before baseline within window."""
    out = merged.copy()
    a = adni.copy()
    a[args.id_col] = clean_id(a[args.id_col])
    a["_row_date"] = pd.to_datetime(a[args.date_col], errors="coerce")

    csf_available = [col for col in csf_cols if col in a.columns]
    grouped = {pid: g.sort_values("_row_date") for pid, g in a.groupby(args.id_col, sort=False)}
    qc_rows: list[dict] = []

    for index, row in out.iterrows():
        pid = str(row[args.id_col]).strip()
        baseline_date = pd.to_datetime(row[args.baseline_date_col], errors="coerce")
        g = grouped.get(pid)
        record = {args.id_col: pid}

        for col in csf_cols:
            source = "baseline_exact" if col in out.columns and pd.notna(row.get(col, np.nan)) else "missing"
            day_difference = 0.0 if source == "baseline_exact" else np.nan

            if source == "missing" and col in csf_available and g is not None and pd.notna(baseline_date):
                candidates = g.loc[
                    g[col].notna()
                    & g["_row_date"].notna()
                    & (g["_row_date"] <= baseline_date)
                    & ((baseline_date - g["_row_date"]).dt.days <= args.csf_lookback_days)
                ]
                if not candidates.empty:
                    selected = candidates.iloc[-1]
                    out.at[index, col] = selected[col]
                    source = "latest_prebaseline"
                    day_difference = float((baseline_date - selected["_row_date"]).days)

            record[f"{col}_source"] = source
            record[f"{col}_days_before_baseline"] = day_difference
        qc_rows.append(record)

    return out, pd.DataFrame(qc_rows)


def prepare_analysis_dataset(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame]:
    survival = read_table(args.survival_file)
    adni = read_table(args.adni_file)

    required_survival = [
        args.id_col,
        args.baseline_visit_col,
        args.baseline_date_col,
        args.time_col,
        args.event_col,
        args.epoch_col,
    ]
    missing = [col for col in required_survival if col not in survival.columns]
    if missing:
        raise ValueError(f"Survival file is missing required columns: {missing}")

    clinical_cols = parse_list(args.clinical_cols)
    csf_cols = parse_list(args.csf_cols)
    spare_cols = parse_list(args.spare_cols)
    requested = list(dict.fromkeys(clinical_cols + csf_cols + spare_cols))

    survival[args.id_col] = clean_id(survival[args.id_col])
    if survival[args.id_col].duplicated().any():
        duplicated = survival.loc[survival[args.id_col].duplicated(), args.id_col].head(10).tolist()
        raise ValueError(f"Survival file must have one row per participant. Duplicates include: {duplicated}")

    exact, match_qc = select_exact_baseline_rows(survival, adni, args, requested)
    # Avoid duplicate suffixes: retain outcome/EPOCH columns from survival, source
    # requested baseline fields from ADNI, and fill survival fields only when needed.
    drop_from_survival = [col for col in requested if col in survival.columns]
    base = survival.drop(columns=drop_from_survival, errors="ignore")
    merged = base.merge(exact, on=args.id_col, how="left", validate="one_to_one")
    merged = supplement_subject_level_clinical(merged, adni, args, clinical_cols)
    merged, csf_qc = attach_prebaseline_csf(merged, adni, args, csf_cols)
    match_qc = match_qc.merge(csf_qc, on=args.id_col, how="left", validate="one_to_one")

    merged[args.time_col] = clean_numeric(merged[args.time_col])
    merged[args.event_col] = parse_bool_event(merged[args.event_col])
    merged[args.epoch_col] = clean_numeric(merged[args.epoch_col])

    for col in ["Age", "Education_Years"] + csf_cols + spare_cols:
        if col in merged.columns:
            merged[col] = clean_numeric(merged[col])
    if "Sex" in merged.columns:
        merged["Sex"] = normalize_sex(merged["Sex"])
    if "APOE_Genotype" in merged.columns:
        merged["APOE_Genotype"] = normalize_apoe(merged["APOE_Genotype"])

    valid = (
        merged[args.id_col].notna()
        & merged[args.time_col].notna()
        & (merged[args.time_col] > 0)
        & merged[args.event_col].notna()
        & merged[args.epoch_col].notna()
    )
    removed = int((~valid).sum())
    if removed:
        warnings.warn(f"Removed {removed} rows with invalid outcome, time, ID, or AD EPOCH score.")
    merged = merged.loc[valid].copy()
    merged[args.event_col] = merged[args.event_col].astype(int)

    return merged, match_qc


def predictor_types(predictors: Sequence[str]) -> tuple[list[str], list[str]]:
    categorical = [col for col in predictors if col in {"Sex", "APOE_Genotype"}]
    numeric = [col for col in predictors if col not in categorical]
    return numeric, categorical


def fit_preprocessor(
    train: pd.DataFrame,
    test: pd.DataFrame,
    predictors: Sequence[str],
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    numeric, categorical = predictor_types(predictors)
    train_out = pd.DataFrame(index=train.index)
    test_out = pd.DataFrame(index=test.index)
    metadata: dict = {"numeric": {}, "categorical": {}}

    for col in numeric:
        train_values = clean_numeric(train[col])
        test_values = clean_numeric(test[col])
        median = float(train_values.median()) if train_values.notna().any() else 0.0
        train_values = train_values.fillna(median)
        test_values = test_values.fillna(median)
        mean = float(train_values.mean())
        sd = float(train_values.std(ddof=0))
        if not np.isfinite(sd) or sd <= 0:
            sd = 1.0
        train_out[col] = (train_values - mean) / sd
        test_out[col] = (test_values - mean) / sd
        metadata["numeric"][col] = {"median": median, "mean": mean, "sd": sd}

    for col in categorical:
        train_values = train[col].astype("object")
        test_values = test[col].astype("object")
        mode = train_values.dropna().mode()
        fill = str(mode.iloc[0]) if not mode.empty else "Missing"
        train_values = train_values.where(train_values.notna(), fill).astype(str)
        test_values = test_values.where(test_values.notna(), fill).astype(str)
        levels = sorted(train_values.unique().tolist())
        reference = levels[0] if levels else fill
        metadata["categorical"][col] = {"fill": fill, "levels": levels, "reference": reference}
        for level in levels[1:]:
            safe = re.sub(r"[^A-Za-z0-9]+", "_", level).strip("_") or "level"
            name = f"{col}__{safe}"
            train_out[name] = (train_values == level).astype(float)
            test_out[name] = (test_values == level).astype(float)

    # Remove zero-variance columns in the training fold.
    keep = [col for col in train_out.columns if train_out[col].nunique(dropna=False) > 1]
    return train_out[keep], test_out.reindex(columns=keep, fill_value=0.0), metadata


def fit_cox_model(
    train: pd.DataFrame,
    test: pd.DataFrame,
    predictors: Sequence[str],
    time_col: str,
    event_col: str,
    penalizer: float,
) -> tuple[CoxPHFitter, pd.DataFrame, pd.DataFrame, dict]:
    x_train, x_test, metadata = fit_preprocessor(train, test, predictors)
    if x_train.shape[1] == 0:
        raise ValueError("No nonconstant transformed predictors remain.")

    fit_data = x_train.copy()
    fit_data[time_col] = train[time_col].astype(float).values
    fit_data[event_col] = train[event_col].astype(int).values

    last_error: Exception | None = None
    for candidate_penalizer in [penalizer, max(penalizer, 0.05), 0.1, 0.5, 1.0]:
        try:
            model = CoxPHFitter(penalizer=candidate_penalizer)
            model.fit(fit_data, duration_col=time_col, event_col=event_col, show_progress=False)
            metadata["penalizer_used"] = candidate_penalizer
            return model, x_train, x_test, metadata
        except Exception as exc:
            last_error = exc
    raise RuntimeError(f"Cox fitting failed after penalizer fallbacks: {last_error}")


def harrell_cindex(time: np.ndarray, event: np.ndarray, risk: np.ndarray) -> float:
    # lifelines concordance_index expects larger prediction = longer survival;
    # Cox partial hazard is larger for shorter survival, hence the negative sign.
    return float(concordance_index(time, -np.asarray(risk), event))


def cross_validated_predictions(
    data: pd.DataFrame,
    spec: ModelSpec,
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, float, int]:
    event_counts = data[args.event_col].value_counts()
    if len(event_counts) < 2:
        raise ValueError("Both events and censored observations are required.")
    folds = min(args.cv_folds, int(event_counts.min()))
    if folds < 2:
        raise ValueError("Insufficient event/censor counts for cross-validation.")

    splitter = StratifiedKFold(n_splits=folds, shuffle=True, random_state=args.random_state)
    predictions = np.repeat(np.nan, data.shape[0])
    fold_ids = np.repeat(-1, data.shape[0])

    for fold, (train_index, test_index) in enumerate(
        splitter.split(np.zeros(data.shape[0]), data[args.event_col].values), start=1
    ):
        train = data.iloc[train_index].copy()
        test = data.iloc[test_index].copy()
        model, _, x_test, _ = fit_cox_model(
            train,
            test,
            spec.predictors,
            args.time_col,
            args.event_col,
            args.cox_penalizer,
        )
        predictions[test_index] = model.predict_partial_hazard(x_test).values.reshape(-1)
        fold_ids[test_index] = fold

    cindex = harrell_cindex(
        data[args.time_col].to_numpy(float),
        data[args.event_col].to_numpy(int),
        predictions,
    )
    output = data[[args.id_col, args.time_col, args.event_col]].copy()
    output["analysis_model"] = spec.name
    output["fold"] = fold_ids
    output["oof_risk"] = predictions
    return output, cindex, folds


def fit_full_model(
    data: pd.DataFrame,
    spec: ModelSpec,
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, float, int, float]:
    model, x, _, metadata = fit_cox_model(
        data,
        data,
        spec.predictors,
        args.time_col,
        args.event_col,
        args.cox_penalizer,
    )
    summary = model.summary.reset_index().rename(columns={"covariate": "term", "index": "term"})
    rename = {
        "exp(coef)": "hazard_ratio",
        "exp(coef) lower 95%": "hazard_ratio_ci_lower",
        "exp(coef) upper 95%": "hazard_ratio_ci_upper",
        "p": "p_value",
        "coef": "log_hazard_coefficient",
        "se(coef)": "standard_error",
    }
    summary = summary.rename(columns=rename)
    summary["analysis_model"] = spec.name
    summary["n"] = data.shape[0]
    summary["n_events"] = int(data[args.event_col].sum())
    summary["penalizer"] = metadata["penalizer_used"]
    wanted = [
        "analysis_model",
        "term",
        "log_hazard_coefficient",
        "hazard_ratio",
        "hazard_ratio_ci_lower",
        "hazard_ratio_ci_upper",
        "standard_error",
        "p_value",
        "n",
        "n_events",
        "penalizer",
    ]
    summary = summary[[col for col in wanted if col in summary.columns]]
    risk = model.predict_partial_hazard(x).values.reshape(-1)
    apparent_cindex = harrell_cindex(
        data[args.time_col].to_numpy(float),
        data[args.event_col].to_numpy(int),
        risk,
    )
    return summary, float(model.log_likelihood_), x.shape[1], apparent_cindex


def paired_bootstrap_delta(
    data: pd.DataFrame,
    prediction_a: np.ndarray,
    prediction_b: np.ndarray,
    n_bootstrap: int,
    random_state: int,
) -> dict:
    rng = np.random.default_rng(random_state)
    time = data["_analysis_time"].to_numpy(float)
    event = data["_analysis_event"].to_numpy(int)
    prediction_a = np.asarray(prediction_a, dtype=float)
    prediction_b = np.asarray(prediction_b, dtype=float)

    observed_a = harrell_cindex(time, event, prediction_a)
    observed_b = harrell_cindex(time, event, prediction_b)
    observed_delta = observed_a - observed_b
    values: list[float] = []

    for _ in range(n_bootstrap):
        indices = rng.integers(0, len(data), size=len(data))
        if np.unique(event[indices]).size < 2:
            continue
        try:
            delta = harrell_cindex(time[indices], event[indices], prediction_a[indices]) - harrell_cindex(
                time[indices], event[indices], prediction_b[indices]
            )
            if np.isfinite(delta):
                values.append(delta)
        except Exception:
            continue

    array = np.asarray(values, dtype=float)
    if array.size:
        lower, upper = np.quantile(array, [0.025, 0.975])
        p_le_zero = float(np.mean(array <= 0))
        p_ge_zero = float(np.mean(array >= 0))
        p_two = min(1.0, 2.0 * min(p_le_zero, p_ge_zero))
    else:
        lower = upper = p_two = np.nan

    return {
        "cindex_model_a": observed_a,
        "cindex_model_b": observed_b,
        "delta_cindex_a_minus_b": observed_delta,
        "delta_ci_lower": lower,
        "delta_ci_upper": upper,
        "bootstrap_p_two_sided": p_two,
        "n_bootstrap_successful": int(array.size),
    }


def make_groups(csf_cols: Sequence[str], spare_cols: Sequence[str]) -> list[AnalysisGroup]:
    groups = [AnalysisGroup("epoch_all", tuple())]
    for col in list(csf_cols) + list(spare_cols):
        groups.append(AnalysisGroup(col.lower(), (col,)))
    if csf_cols:
        groups.append(AnalysisGroup("csf_panel", tuple(csf_cols)))
    if spare_cols:
        groups.append(AnalysisGroup("spare_panel", tuple(spare_cols)))
    if csf_cols or spare_cols:
        groups.append(AnalysisGroup("all_biomarkers", tuple(list(csf_cols) + list(spare_cols))))
    return groups


def make_model_specs(
    clinical_cols: Sequence[str],
    biomarkers: Sequence[str],
    epoch_col: str,
) -> list[ModelSpec]:
    clinical = tuple(clinical_cols)
    biomarkers = tuple(biomarkers)
    specs = [
        ModelSpec("clinical", clinical),
        ModelSpec("clinical_plus_epoch", clinical + (epoch_col,)),
    ]
    if biomarkers:
        specs.extend(
            [
                ModelSpec("clinical_plus_biomarkers", clinical + biomarkers),
                ModelSpec("clinical_plus_biomarkers_plus_epoch", clinical + biomarkers + (epoch_col,)),
            ]
        )
    return specs


def run_analyses(data: pd.DataFrame, args: argparse.Namespace) -> dict[str, pd.DataFrame]:
    clinical_cols = [col for col in parse_list(args.clinical_cols) if col in data.columns]
    csf_cols = [col for col in parse_list(args.csf_cols) if col in data.columns]
    spare_cols = [col for col in parse_list(args.spare_cols) if col in data.columns]

    missing_clinical = sorted(set(parse_list(args.clinical_cols)) - set(clinical_cols))
    if missing_clinical:
        warnings.warn(f"Clinical covariates absent and omitted: {missing_clinical}")

    performance_rows: list[dict] = []
    coefficient_tables: list[pd.DataFrame] = []
    prediction_tables: list[pd.DataFrame] = []
    comparison_rows: list[dict] = []
    status_rows: list[dict] = []

    for group in make_groups(csf_cols, spare_cols):
        required_complete = [args.epoch_col] + list(group.biomarkers)
        subset = data.dropna(subset=required_complete + [args.time_col, args.event_col]).copy()
        n = subset.shape[0]
        events = int(subset[args.event_col].sum())
        censored = int(n - events)

        if n < args.minimum_n or events < args.minimum_events or censored < args.minimum_events:
            status_rows.append(
                {
                    "analysis_group": group.name,
                    "status": "skipped_insufficient_sample",
                    "n": n,
                    "n_events": events,
                    "n_censored": censored,
                    "required_biomarkers": ",".join(group.biomarkers),
                }
            )
            continue

        status_rows.append(
            {
                "analysis_group": group.name,
                "status": "analyzed",
                "n": n,
                "n_events": events,
                "n_censored": censored,
                "required_biomarkers": ",".join(group.biomarkers),
            }
        )

        specs = make_model_specs(clinical_cols, group.biomarkers, args.epoch_col)
        model_results: dict[str, dict] = {}

        for spec in specs:
            try:
                oof, cv_cindex, folds = cross_validated_predictions(subset, spec, args)
                coefficients, log_likelihood, n_parameters, apparent_cindex = fit_full_model(subset, spec, args)
                oof["analysis_group"] = group.name
                coefficient_tables.append(coefficients.assign(analysis_group=group.name))
                prediction_tables.append(oof)
                performance_rows.append(
                    {
                        "analysis_group": group.name,
                        "model": spec.name,
                        "predictors": ",".join(spec.predictors),
                        "n": n,
                        "n_events": events,
                        "n_censored": censored,
                        "cv_folds": folds,
                        "cv_cindex": cv_cindex,
                        "apparent_cindex": apparent_cindex,
                        "log_likelihood": log_likelihood,
                        "n_parameters": n_parameters,
                        "status": "success",
                    }
                )
                model_results[spec.name] = {
                    "oof": oof,
                    "cv_cindex": cv_cindex,
                    "log_likelihood": log_likelihood,
                    "n_parameters": n_parameters,
                }
            except Exception as exc:
                performance_rows.append(
                    {
                        "analysis_group": group.name,
                        "model": spec.name,
                        "predictors": ",".join(spec.predictors),
                        "n": n,
                        "n_events": events,
                        "n_censored": censored,
                        "status": "failed",
                        "message": str(exc),
                    }
                )

        pairs = [("clinical_plus_epoch", "clinical")]
        if group.biomarkers:
            pairs.extend(
                [
                    ("clinical_plus_biomarkers_plus_epoch", "clinical_plus_biomarkers"),
                    ("clinical_plus_biomarkers_plus_epoch", "clinical_plus_epoch"),
                    ("clinical_plus_biomarkers", "clinical"),
                ]
            )

        for model_a, model_b in pairs:
            if model_a not in model_results or model_b not in model_results:
                continue
            pred_a = model_results[model_a]["oof"]["oof_risk"].to_numpy(float)
            pred_b = model_results[model_b]["oof"]["oof_risk"].to_numpy(float)
            bootstrap_data = subset[[args.time_col, args.event_col]].copy().rename(
                columns={args.time_col: "_analysis_time", args.event_col: "_analysis_event"}
            )
            delta = paired_bootstrap_delta(
                bootstrap_data,
                pred_a,
                pred_b,
                args.bootstrap,
                args.random_state,
            )
            ll_a = model_results[model_a]["log_likelihood"]
            ll_b = model_results[model_b]["log_likelihood"]
            df_a = model_results[model_a]["n_parameters"]
            df_b = model_results[model_b]["n_parameters"]
            lr_stat = 2.0 * (ll_a - ll_b)
            lr_df = df_a - df_b
            lr_p = float(chi2.sf(max(lr_stat, 0.0), lr_df)) if lr_df > 0 else np.nan
            comparison_rows.append(
                {
                    "analysis_group": group.name,
                    "model_a": model_a,
                    "model_b": model_b,
                    "n": n,
                    "n_events": events,
                    **delta,
                    "likelihood_ratio_chi2": lr_stat,
                    "likelihood_ratio_df": lr_df,
                    "likelihood_ratio_p": lr_p,
                }
            )

    return {
        "performance": pd.DataFrame(performance_rows),
        "coefficients": pd.concat(coefficient_tables, ignore_index=True) if coefficient_tables else pd.DataFrame(),
        "oof_predictions": pd.concat(prediction_tables, ignore_index=True) if prediction_tables else pd.DataFrame(),
        "comparisons": pd.DataFrame(comparison_rows),
        "status": pd.DataFrame(status_rows),
    }


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    log("============================================================")
    log("ADNI MCI-to-AD: AD EPOCH versus clinical, CSF, and SPARE biomarkers")
    log("Clinical covariates exclude baseline cognition")
    log("============================================================")

    data, match_qc = prepare_analysis_dataset(args)
    data_path = outdir / f"{args.prefix}_baseline_analysis_dataset.tsv"
    qc_path = outdir / f"{args.prefix}_baseline_matching_qc.tsv"
    data.to_csv(data_path, sep="\t", index=False)
    match_qc.to_csv(qc_path, sep="\t", index=False)

    variables = parse_list(args.clinical_cols) + [args.epoch_col] + parse_list(args.csf_cols) + parse_list(args.spare_cols)
    missingness = pd.DataFrame(
        [
            {
                "variable": col,
                "present": col in data.columns,
                "n_nonmissing": int(data[col].notna().sum()) if col in data.columns else 0,
                "n_missing": int(data[col].isna().sum()) if col in data.columns else data.shape[0],
                "fraction_nonmissing": float(data[col].notna().mean()) if col in data.columns else 0.0,
            }
            for col in variables
        ]
    )
    missingness.to_csv(outdir / f"{args.prefix}_missingness.tsv", sep="\t", index=False)

    results = run_analyses(data, args)
    for name, table in results.items():
        table.to_csv(outdir / f"{args.prefix}_{name}.tsv", sep="\t", index=False)

    summary = {
        "survival_file": str(args.survival_file),
        "adni_file": str(args.adni_file),
        "n_participants": int(data.shape[0]),
        "n_events": int(data[args.event_col].sum()),
        "n_censored": int(data.shape[0] - data[args.event_col].sum()),
        "clinical_covariates": parse_list(args.clinical_cols),
        "baseline_cognition_included": False,
        "epoch_predictor": args.epoch_col,
        "csf_biomarkers": parse_list(args.csf_cols),
        "spare_biomarkers": parse_list(args.spare_cols),
        "csf_rule": (
            "Exact selected baseline row first; otherwise latest measurement on or before "
            f"baseline within {args.csf_lookback_days} days; no post-baseline CSF used."
        ),
        "cv_folds_requested": args.cv_folds,
        "bootstrap_requested": args.bootstrap,
        "cox_penalizer_requested": args.cox_penalizer,
        "primary_output_performance": str(outdir / f"{args.prefix}_performance.tsv"),
        "primary_output_comparisons": str(outdir / f"{args.prefix}_comparisons.tsv"),
    }
    with open(outdir / f"{args.prefix}_summary.json", "w") as handle:
        json.dump(summary, handle, indent=2)

    log(f"Participants analyzed before biomarker-specific filtering: {data.shape[0]}")
    log(f"Events: {int(data[args.event_col].sum())}")
    log(f"Censored: {int(data.shape[0] - data[args.event_col].sum())}")
    log("Analysis-group status:")
    log(results["status"].to_string(index=False))
    log("Main model performance:")
    if not results["performance"].empty:
        cols = [col for col in ["analysis_group", "model", "n", "n_events", "cv_cindex", "status"] if col in results["performance"].columns]
        log(results["performance"][cols].to_string(index=False))
    log("============================================================")
    log(f"Baseline dataset: {data_path}")
    log(f"Model performance: {outdir / f'{args.prefix}_performance.tsv'}")
    log(f"Nested/delta comparisons: {outdir / f'{args.prefix}_comparisons.tsv'}")
    log("============================================================")


if __name__ == "__main__":
    main()
