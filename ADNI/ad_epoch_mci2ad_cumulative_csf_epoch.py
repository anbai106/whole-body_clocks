#!/usr/bin/env python3
"""
Cumulative survival analysis for MCI-to-AD conversion in ADNI using demographic, CSF, and AD EPOCH predictors.

Cohort
------
One row per participant from the previously generated MCI-to-AD EPOCH
survival file. That file should already exclude participants used to fit the
original AD EPOCH model (for example, original train and validation IDs).

Outcome
-------
Time zero: first qualifying MCI MRI with usable MUSE data.
Event: first subsequent AD diagnosis.
Censoring: last follow-up without AD.

Cumulative model order
----------------------
M0: Age + Sex
M1: M0 + Abeta_CSF
M2: M1 + Tau_CSF + PTau_CSF
M3: M2 + AD EPOCH

Relaxed CSF definition
----------------------
For each participant and each CSF marker separately, use the first
chronologically available nonmissing measurement anywhere in ADNI. The
measurement may occur before or after the selected baseline-MCI MRI. Timing
relative to baseline and event/censoring is retained in QC outputs so the
analysis can be described as a cumulative/incremental-information analysis
rather than a strictly prospective baseline-biomarker analysis.

Primary comparison
------------------
All cumulative models are fit on one common complete-case cohort containing
all variables required through M3. This guarantees paired comparisons across
steps. A sample-flow table also reports the available N at each stage.

Outputs
-------
* baseline_analysis_dataset.tsv
* common_complete_case_dataset.tsv
* biomarker_matching_qc.tsv
* first_csf_timing_summary.tsv
* missingness.tsv
* sample_flow.tsv
* performance.tsv
* comparisons.tsv
* coefficients.tsv
* oof_predictions.tsv
* summary.json

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
from typing import Sequence

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from scipy.stats import chi2
from sklearn.model_selection import StratifiedKFold


DEFAULT_EPOCH_COL = "adni_brain_mri_ad_lepoch_risk_score"
DEFAULT_CSF = ["Abeta_CSF", "Tau_CSF", "PTau_CSF"]


@dataclass(frozen=True)
class ModelSpec:
    step: int
    name: str
    added_block: str
    predictors: tuple[str, ...]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Cumulative MCI-to-AD survival analysis with demographic covariates, CSF biomarkers, and AD EPOCH."
    )
    p.add_argument("--survival-file", required=True)
    p.add_argument("--adni-file", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--prefix", default="adni_mci2ad_cumulative_biomarkers")

    p.add_argument("--id-col", default="PTID")
    p.add_argument("--visit-col", default="Visit_Code")
    p.add_argument("--date-col", default="Date")
    p.add_argument("--baseline-visit-col", default="selected_baseline_visit_code")
    p.add_argument("--baseline-date-col", default="selected_baseline_date")
    p.add_argument("--time-col", default="time_years")
    p.add_argument("--event-col", default="event")
    p.add_argument("--epoch-col", default=DEFAULT_EPOCH_COL)

    p.add_argument("--baseline-cols", default="Age,Sex")
    p.add_argument("--csf-cols", default=",".join(DEFAULT_CSF))

    p.add_argument("--cv-folds", type=int, default=5)
    p.add_argument("--random-state", type=int, default=20260730)
    p.add_argument("--cox-penalizer", type=float, default=0.01)
    p.add_argument("--bootstrap", type=int, default=1000)
    p.add_argument("--minimum-n", type=int, default=50)
    p.add_argument("--minimum-events", type=int, default=20)
    p.add_argument(
        "--allow-stage-specific-samples",
        action="store_true",
        help=(
            "Also run a secondary analysis in which each cumulative step uses its own "
            "maximum available sample. Primary paired comparisons always use the common sample."
        ),
    )
    return p.parse_args()


def log(msg: str) -> None:
    print(msg, flush=True)


def parse_list(x: str | None) -> list[str]:
    if x is None or str(x).strip() == "":
        return []
    return [v.strip() for v in str(x).split(",") if v.strip()]


def read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty input file: {path}")
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path, low_memory=False)
    return pd.read_csv(path, sep="\t", low_memory=False)


def clean_id(s: pd.Series) -> pd.Series:
    return s.astype("object").where(s.notna(), np.nan).astype(str).str.strip()


def clean_numeric(s: pd.Series) -> pd.Series:
    if pd.api.types.is_numeric_dtype(s):
        return pd.to_numeric(s, errors="coerce")
    cleaned = (
        s.astype("object")
        .where(s.notna(), np.nan)
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
            }
        )
    )
    extracted = cleaned.astype(str).str.extract(
        r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", expand=False
    )
    extracted = extracted.where(cleaned.notna(), np.nan)
    return pd.to_numeric(extracted, errors="coerce")


def normalize_sex(s: pd.Series) -> pd.Series:
    return (
        s.astype("object")
        .where(s.notna(), np.nan)
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


def parse_event(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.astype(int)
    out = pd.to_numeric(s, errors="coerce")
    missing = out.isna()
    text = s.astype(str).str.strip().str.lower()
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
    out.loc[missing] = text.loc[missing].map(mapping)
    return out.astype("Int64")


def select_baseline_rows(
    survival: pd.DataFrame,
    adni: pd.DataFrame,
    args: argparse.Namespace,
    requested_cols: Sequence[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Select the ADNI row corresponding to each participant's selected MCI baseline."""
    s = survival[
        [args.id_col, args.baseline_visit_col, args.baseline_date_col]
    ].copy()
    s[args.id_col] = clean_id(s[args.id_col])
    s["_baseline_visit"] = (
        s[args.baseline_visit_col].astype(str).str.strip().str.lower()
    )
    s["_baseline_date"] = pd.to_datetime(
        s[args.baseline_date_col], errors="coerce"
    )

    a = adni.copy()
    a[args.id_col] = clean_id(a[args.id_col])
    a["_row_visit"] = a[args.visit_col].astype(str).str.strip().str.lower()
    a["_row_date"] = pd.to_datetime(a[args.date_col], errors="coerce")

    missing = [c for c in requested_cols if c not in a.columns]
    if missing:
        warnings.warn(f"Requested baseline columns absent from ADNI file: {missing}")

    grouped = {pid: g.copy() for pid, g in a.groupby(args.id_col, sort=False)}
    rows: list[dict] = []
    qc: list[dict] = []

    for _, rec in s.iterrows():
        pid = str(rec[args.id_col]).strip()
        baseline_visit = str(rec["_baseline_visit"]).strip().lower()
        baseline_date = rec["_baseline_date"]
        g = grouped.get(pid)

        selected = None
        method = "unmatched"
        date_diff = np.nan

        if g is not None and not g.empty:
            by_visit = g.loc[g["_row_visit"] == baseline_visit].copy()
            if not by_visit.empty:
                if pd.notna(baseline_date) and by_visit["_row_date"].notna().any():
                    diff = (by_visit["_row_date"] - baseline_date).abs().dt.days
                    idx = diff.idxmin()
                    selected = by_visit.loc[idx]
                    date_diff = float(diff.loc[idx])
                else:
                    selected = by_visit.iloc[0]
                method = "PTID+Visit_Code"
            elif pd.notna(baseline_date) and g["_row_date"].notna().any():
                diff = (g["_row_date"] - baseline_date).abs().dt.days
                idx = diff.idxmin()
                if float(diff.loc[idx]) == 0.0:
                    selected = g.loc[idx]
                    date_diff = 0.0
                    method = "PTID+exact_Date"

        out = {args.id_col: pid}
        for col in requested_cols:
            out[col] = (
                selected[col]
                if selected is not None and col in selected.index
                else np.nan
            )
        rows.append(out)
        qc.append(
            {
                args.id_col: pid,
                "selected_baseline_visit_code": rec[args.baseline_visit_col],
                "selected_baseline_date": baseline_date,
                "baseline_match_method": method,
                "baseline_match_date_difference_days": date_diff,
                "baseline_row_found": selected is not None,
            }
        )

    return pd.DataFrame(rows), pd.DataFrame(qc)


def first_available_measurements(
    adni: pd.DataFrame,
    args: argparse.Namespace,
    marker_cols: Sequence[str],
) -> pd.DataFrame:
    """Return first chronologically available nonmissing measurement per marker and subject."""
    a = adni.copy()
    a[args.id_col] = clean_id(a[args.id_col])
    a["_measure_date"] = pd.to_datetime(a[args.date_col], errors="coerce")
    a["_visit_order"] = np.arange(a.shape[0], dtype=int)

    rows: list[dict] = []
    for pid, g in a.groupby(args.id_col, sort=False):
        record: dict = {args.id_col: pid}
        for marker in marker_cols:
            if marker not in g.columns:
                record[marker] = np.nan
                record[f"{marker}_first_visit_code"] = np.nan
                record[f"{marker}_first_date"] = pd.NaT
                continue

            candidates = g.loc[g[marker].notna()].copy()
            if candidates.empty:
                record[marker] = np.nan
                record[f"{marker}_first_visit_code"] = np.nan
                record[f"{marker}_first_date"] = pd.NaT
                continue

            candidates = candidates.sort_values(
                ["_measure_date", "_visit_order"], na_position="last", kind="mergesort"
            )
            selected = candidates.iloc[0]
            record[marker] = selected[marker]
            record[f"{marker}_first_visit_code"] = selected[args.visit_col]
            record[f"{marker}_first_date"] = selected["_measure_date"]
        rows.append(record)

    return pd.DataFrame(rows)


def prepare_dataset(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    survival = read_table(args.survival_file)
    adni = read_table(args.adni_file)

    required = [
        args.id_col,
        args.baseline_visit_col,
        args.baseline_date_col,
        args.time_col,
        args.event_col,
        args.epoch_col,
    ]
    missing = [c for c in required if c not in survival.columns]
    if missing:
        raise ValueError(f"Survival file is missing required columns: {missing}")

    survival[args.id_col] = clean_id(survival[args.id_col])
    if survival[args.id_col].duplicated().any():
        dup = survival.loc[
            survival[args.id_col].duplicated(), args.id_col
        ].head(10).tolist()
        raise ValueError(f"Survival file must have one row per participant. Duplicates: {dup}")

    baseline_cols = parse_list(args.baseline_cols)
    csf_cols = parse_list(args.csf_cols)
    baseline_requested = baseline_cols

    baseline_data, baseline_qc = select_baseline_rows(
        survival, adni, args, baseline_requested
    )
    first_csf = first_available_measurements(adni, args, csf_cols)

    drop_cols = [c for c in baseline_requested + csf_cols if c in survival.columns]
    base = survival.drop(columns=drop_cols, errors="ignore").copy()
    data = base.merge(baseline_data, on=args.id_col, how="left", validate="one_to_one")
    data = data.merge(first_csf, on=args.id_col, how="left", validate="one_to_one")

    data[args.time_col] = clean_numeric(data[args.time_col])
    data[args.event_col] = parse_event(data[args.event_col])
    data[args.epoch_col] = clean_numeric(data[args.epoch_col])

    for col in ["Age"] + csf_cols:
        if col in data.columns:
            data[col] = clean_numeric(data[col])
    if "Sex" in data.columns:
        data["Sex"] = normalize_sex(data["Sex"])

    data[args.baseline_date_col] = pd.to_datetime(
        data[args.baseline_date_col], errors="coerce"
    )
    if "event_or_censor_date" in data.columns:
        data["event_or_censor_date"] = pd.to_datetime(
            data["event_or_censor_date"], errors="coerce"
        )

    timing_rows: list[dict] = []
    for _, row in data.iterrows():
        pid = row[args.id_col]
        baseline_date = row[args.baseline_date_col]
        end_date = row.get("event_or_censor_date", pd.NaT)
        record: dict = {args.id_col: pid}
        for marker in csf_cols:
            date_col = f"{marker}_first_date"
            measure_date = pd.to_datetime(row.get(date_col, pd.NaT), errors="coerce")
            if pd.notna(measure_date) and pd.notna(baseline_date):
                days_from_baseline = float((measure_date - baseline_date).days)
                if days_from_baseline < 0:
                    relation = "prebaseline"
                elif days_from_baseline == 0:
                    relation = "same_day"
                else:
                    relation = "postbaseline"
            else:
                days_from_baseline = np.nan
                relation = "missing_date_or_marker"

            if pd.notna(measure_date) and pd.notna(end_date):
                after_end = bool(measure_date > end_date)
                days_from_end = float((measure_date - end_date).days)
            else:
                after_end = np.nan
                days_from_end = np.nan

            record[f"{marker}_days_from_baseline"] = days_from_baseline
            record[f"{marker}_timing_relation"] = relation
            record[f"{marker}_after_event_or_censor"] = after_end
            record[f"{marker}_days_from_event_or_censor"] = days_from_end
        timing_rows.append(record)

    timing_qc = pd.DataFrame(timing_rows)
    baseline_qc = baseline_qc.merge(timing_qc, on=args.id_col, how="left", validate="one_to_one")

    valid = (
        data[args.id_col].notna()
        & data[args.time_col].notna()
        & (data[args.time_col] > 0)
        & data[args.event_col].notna()
        & data[args.epoch_col].notna()
    )
    if (~valid).any():
        warnings.warn(
            f"Removed {int((~valid).sum())} rows with invalid ID, time, event, or EPOCH score."
        )
    data = data.loc[valid].copy()
    data[args.event_col] = data[args.event_col].astype(int)

    timing_summary_rows: list[dict] = []
    for marker in csf_cols:
        rel_col = f"{marker}_timing_relation"
        if rel_col in baseline_qc.columns:
            counts = baseline_qc[rel_col].value_counts(dropna=False)
            for relation, n in counts.items():
                timing_summary_rows.append(
                    {"marker": marker, "timing_relation": relation, "n": int(n)}
                )
        after_col = f"{marker}_after_event_or_censor"
        if after_col in baseline_qc.columns:
            counts = baseline_qc[after_col].value_counts(dropna=False)
            for status, n in counts.items():
                timing_summary_rows.append(
                    {
                        "marker": marker,
                        "timing_relation": f"after_event_or_censor={status}",
                        "n": int(n),
                    }
                )

    return data, baseline_qc, pd.DataFrame(timing_summary_rows)


def make_specs(args: argparse.Namespace) -> list[ModelSpec]:
    baseline = tuple(parse_list(args.baseline_cols))
    csf = parse_list(args.csf_cols)
    abeta = tuple(csf[:1])
    tau = tuple(csf[1:])

    m0 = baseline
    m1 = m0 + abeta
    m2 = m1 + tau
    m3 = m2 + (args.epoch_col,)

    return [
        ModelSpec(0, "M0_age_sex", "baseline", m0),
        ModelSpec(1, "M1_plus_amyloid", "Abeta_CSF", m1),
        ModelSpec(2, "M2_plus_tau", "Tau_CSF+PTau_CSF", m2),
        ModelSpec(3, "M3_plus_AD_EPOCH", args.epoch_col, m3),
    ]


def predictor_types(predictors: Sequence[str]) -> tuple[list[str], list[str]]:
    categorical = [c for c in predictors if c == "Sex"]
    numeric = [c for c in predictors if c not in categorical]
    return numeric, categorical


def preprocess(
    train: pd.DataFrame,
    test: pd.DataFrame,
    predictors: Sequence[str],
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    numeric, categorical = predictor_types(predictors)
    x_train = pd.DataFrame(index=train.index)
    x_test = pd.DataFrame(index=test.index)
    meta: dict = {"numeric": {}, "categorical": {}}

    for col in numeric:
        tr = clean_numeric(train[col])
        te = clean_numeric(test[col])
        median = float(tr.median()) if tr.notna().any() else 0.0
        tr = tr.fillna(median)
        te = te.fillna(median)
        mean = float(tr.mean())
        sd = float(tr.std(ddof=0))
        if not np.isfinite(sd) or sd <= 0:
            sd = 1.0
        x_train[col] = (tr - mean) / sd
        x_test[col] = (te - mean) / sd
        meta["numeric"][col] = {"median": median, "mean": mean, "sd": sd}

    for col in categorical:
        tr = train[col].astype("object")
        te = test[col].astype("object")
        mode = tr.dropna().mode()
        fill = str(mode.iloc[0]) if not mode.empty else "Missing"
        tr = tr.where(tr.notna(), fill).astype(str)
        te = te.where(te.notna(), fill).astype(str)
        levels = sorted(tr.unique().tolist())
        reference = levels[0] if levels else fill
        meta["categorical"][col] = {
            "fill": fill,
            "levels": levels,
            "reference": reference,
        }
        for level in levels[1:]:
            safe = re.sub(r"[^A-Za-z0-9]+", "_", level).strip("_") or "level"
            name = f"{col}__{safe}"
            x_train[name] = (tr == level).astype(float)
            x_test[name] = (te == level).astype(float)

    keep = [c for c in x_train.columns if x_train[c].nunique(dropna=False) > 1]
    return x_train[keep], x_test.reindex(columns=keep, fill_value=0.0), meta


def fit_cox(
    train: pd.DataFrame,
    test: pd.DataFrame,
    predictors: Sequence[str],
    args: argparse.Namespace,
) -> tuple[CoxPHFitter, pd.DataFrame, pd.DataFrame, dict]:
    x_train, x_test, meta = preprocess(train, test, predictors)
    if x_train.shape[1] == 0:
        raise ValueError("No nonconstant transformed predictors remain.")

    fit_data = x_train.copy()
    fit_data[args.time_col] = train[args.time_col].astype(float).values
    fit_data[args.event_col] = train[args.event_col].astype(int).values

    last_error: Exception | None = None
    candidates = []
    for value in [args.cox_penalizer, max(args.cox_penalizer, 0.05), 0.1, 0.5, 1.0]:
        if value not in candidates:
            candidates.append(value)

    for penalizer in candidates:
        try:
            model = CoxPHFitter(penalizer=penalizer)
            model.fit(
                fit_data,
                duration_col=args.time_col,
                event_col=args.event_col,
                show_progress=False,
            )
            meta["penalizer_used"] = penalizer
            return model, x_train, x_test, meta
        except Exception as exc:
            last_error = exc

    raise RuntimeError(f"Cox fitting failed after penalizer fallbacks: {last_error}")


def cindex(time: np.ndarray, event: np.ndarray, risk: np.ndarray) -> float:
    return float(concordance_index(time, -np.asarray(risk), event))


def cv_predictions(
    data: pd.DataFrame,
    spec: ModelSpec,
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, float, int]:
    counts = data[args.event_col].value_counts()
    if len(counts) < 2:
        raise ValueError("Both event and censored observations are required.")
    folds = min(args.cv_folds, int(counts.min()))
    if folds < 2:
        raise ValueError("Insufficient event/censor counts for cross-validation.")

    splitter = StratifiedKFold(
        n_splits=folds, shuffle=True, random_state=args.random_state
    )
    risk = np.full(data.shape[0], np.nan)
    fold_id = np.full(data.shape[0], -1, dtype=int)

    for fold, (tr_idx, te_idx) in enumerate(
        splitter.split(np.zeros(data.shape[0]), data[args.event_col].values), start=1
    ):
        tr = data.iloc[tr_idx].copy()
        te = data.iloc[te_idx].copy()
        model, _, x_test, _ = fit_cox(tr, te, spec.predictors, args)
        risk[te_idx] = model.predict_partial_hazard(x_test).values.reshape(-1)
        fold_id[te_idx] = fold

    score = cindex(
        data[args.time_col].to_numpy(float),
        data[args.event_col].to_numpy(int),
        risk,
    )
    out = data[[args.id_col, args.time_col, args.event_col]].copy()
    out["model_step"] = spec.step
    out["model"] = spec.name
    out["fold"] = fold_id
    out["oof_risk"] = risk
    return out, score, folds


def full_model(
    data: pd.DataFrame,
    spec: ModelSpec,
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, float, int, float, float]:
    model, x, _, meta = fit_cox(data, data, spec.predictors, args)
    summary = model.summary.reset_index().rename(
        columns={"covariate": "term", "index": "term"}
    )
    summary = summary.rename(
        columns={
            "coef": "log_hazard_coefficient",
            "exp(coef)": "hazard_ratio",
            "exp(coef) lower 95%": "hazard_ratio_ci_lower",
            "exp(coef) upper 95%": "hazard_ratio_ci_upper",
            "se(coef)": "standard_error",
            "p": "p_value",
        }
    )
    summary["model_step"] = spec.step
    summary["model"] = spec.name
    summary["added_block"] = spec.added_block
    summary["n"] = data.shape[0]
    summary["n_events"] = int(data[args.event_col].sum())
    summary["penalizer"] = meta["penalizer_used"]
    wanted = [
        "model_step",
        "model",
        "added_block",
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
    summary = summary[[c for c in wanted if c in summary.columns]]

    risk = model.predict_partial_hazard(x).values.reshape(-1)
    apparent = cindex(
        data[args.time_col].to_numpy(float),
        data[args.event_col].to_numpy(int),
        risk,
    )
    return (
        summary,
        float(model.log_likelihood_),
        x.shape[1],
        apparent,
        float(meta["penalizer_used"]),
    )


def bootstrap_delta(
    data: pd.DataFrame,
    risk_new: np.ndarray,
    risk_old: np.ndarray,
    args: argparse.Namespace,
) -> dict:
    rng = np.random.default_rng(args.random_state)
    time = data[args.time_col].to_numpy(float)
    event = data[args.event_col].to_numpy(int)
    risk_new = np.asarray(risk_new, dtype=float)
    risk_old = np.asarray(risk_old, dtype=float)

    c_new = cindex(time, event, risk_new)
    c_old = cindex(time, event, risk_old)
    delta_obs = c_new - c_old
    values: list[float] = []

    for _ in range(args.bootstrap):
        idx = rng.integers(0, len(data), size=len(data))
        if np.unique(event[idx]).size < 2:
            continue
        try:
            d = cindex(time[idx], event[idx], risk_new[idx]) - cindex(
                time[idx], event[idx], risk_old[idx]
            )
            if np.isfinite(d):
                values.append(d)
        except Exception:
            continue

    arr = np.asarray(values, dtype=float)
    if arr.size:
        lo, hi = np.quantile(arr, [0.025, 0.975])
        p_le = float(np.mean(arr <= 0))
        p_ge = float(np.mean(arr >= 0))
        p_two = min(1.0, 2.0 * min(p_le, p_ge))
    else:
        lo = hi = p_two = np.nan

    return {
        "cindex_new": c_new,
        "cindex_previous": c_old,
        "delta_cindex": delta_obs,
        "delta_ci_lower": lo,
        "delta_ci_upper": hi,
        "bootstrap_p_two_sided": p_two,
        "n_bootstrap_successful": int(arr.size),
    }


def sufficient(data: pd.DataFrame, args: argparse.Namespace) -> tuple[bool, int, int, int]:
    n = data.shape[0]
    events = int(data[args.event_col].sum())
    censored = int(n - events)
    ok = (
        n >= args.minimum_n
        and events >= args.minimum_events
        and censored >= args.minimum_events
    )
    return ok, n, events, censored


def run_model_sequence(
    data: pd.DataFrame,
    specs: Sequence[ModelSpec],
    args: argparse.Namespace,
    analysis_type: str,
) -> dict[str, pd.DataFrame]:
    performance_rows: list[dict] = []
    coefficient_tables: list[pd.DataFrame] = []
    prediction_tables: list[pd.DataFrame] = []
    comparison_rows: list[dict] = []
    status_rows: list[dict] = []
    results: dict[str, dict] = {}

    for spec in specs:
        if analysis_type == "common_complete_case":
            subset = data
        else:
            subset = data.dropna(
                subset=list(spec.predictors)
                + [args.time_col, args.event_col, args.epoch_col]
            ).copy()

        ok, n, events, censored = sufficient(subset, args)
        if not ok:
            status_rows.append(
                {
                    "analysis_type": analysis_type,
                    "model_step": spec.step,
                    "model": spec.name,
                    "status": "skipped_insufficient_sample",
                    "n": n,
                    "n_events": events,
                    "n_censored": censored,
                }
            )
            continue

        try:
            oof, cv_c, folds = cv_predictions(subset, spec, args)
            coef, ll, n_par, apparent, penalizer = full_model(subset, spec, args)
            oof["analysis_type"] = analysis_type
            coef["analysis_type"] = analysis_type
            prediction_tables.append(oof)
            coefficient_tables.append(coef)
            performance_rows.append(
                {
                    "analysis_type": analysis_type,
                    "model_step": spec.step,
                    "model": spec.name,
                    "added_block": spec.added_block,
                    "predictors": ",".join(spec.predictors),
                    "n": n,
                    "n_events": events,
                    "n_censored": censored,
                    "cv_folds": folds,
                    "cv_cindex": cv_c,
                    "apparent_cindex": apparent,
                    "log_likelihood": ll,
                    "n_parameters": n_par,
                    "penalizer": penalizer,
                    "status": "success",
                }
            )
            status_rows.append(
                {
                    "analysis_type": analysis_type,
                    "model_step": spec.step,
                    "model": spec.name,
                    "status": "analyzed",
                    "n": n,
                    "n_events": events,
                    "n_censored": censored,
                }
            )
            results[spec.name] = {
                "data": subset,
                "oof": oof,
                "ll": ll,
                "n_parameters": n_par,
            }
        except Exception as exc:
            performance_rows.append(
                {
                    "analysis_type": analysis_type,
                    "model_step": spec.step,
                    "model": spec.name,
                    "added_block": spec.added_block,
                    "predictors": ",".join(spec.predictors),
                    "n": n,
                    "n_events": events,
                    "n_censored": censored,
                    "status": "failed",
                    "message": str(exc),
                }
            )
            status_rows.append(
                {
                    "analysis_type": analysis_type,
                    "model_step": spec.step,
                    "model": spec.name,
                    "status": "failed",
                    "n": n,
                    "n_events": events,
                    "n_censored": censored,
                    "message": str(exc),
                }
            )

    if analysis_type == "common_complete_case":
        for previous, current in zip(specs[:-1], specs[1:]):
            if previous.name not in results or current.name not in results:
                continue
            old = results[previous.name]
            new = results[current.name]
            same_ids = old["oof"][args.id_col].tolist() == new["oof"][args.id_col].tolist()
            if not same_ids:
                raise RuntimeError(
                    f"Participant order mismatch between {previous.name} and {current.name}."
                )
            delta = bootstrap_delta(
                data,
                new["oof"]["oof_risk"].to_numpy(float),
                old["oof"]["oof_risk"].to_numpy(float),
                args,
            )
            lr_stat = 2.0 * (new["ll"] - old["ll"])
            lr_df = new["n_parameters"] - old["n_parameters"]
            lr_p = float(chi2.sf(max(lr_stat, 0.0), lr_df)) if lr_df > 0 else np.nan
            comparison_rows.append(
                {
                    "analysis_type": analysis_type,
                    "previous_step": previous.step,
                    "previous_model": previous.name,
                    "new_step": current.step,
                    "new_model": current.name,
                    "added_block": current.added_block,
                    "n": data.shape[0],
                    "n_events": int(data[args.event_col].sum()),
                    **delta,
                    "likelihood_ratio_chi2": lr_stat,
                    "likelihood_ratio_df": lr_df,
                    "likelihood_ratio_p": lr_p,
                }
            )

    return {
        "performance": pd.DataFrame(performance_rows),
        "coefficients": pd.concat(coefficient_tables, ignore_index=True)
        if coefficient_tables
        else pd.DataFrame(),
        "oof_predictions": pd.concat(prediction_tables, ignore_index=True)
        if prediction_tables
        else pd.DataFrame(),
        "comparisons": pd.DataFrame(comparison_rows),
        "status": pd.DataFrame(status_rows),
    }


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    log("============================================================")
    log("ADNI MCI-to-AD cumulative CSF + AD EPOCH survival analysis")
    log("Cohort: non-training baseline-MCI participants")
    log("CSF rule: first available nonmissing measurement per subject/marker")
    log("============================================================")

    data, match_qc, timing_summary = prepare_dataset(args)
    specs = make_specs(args)

    all_required = list(dict.fromkeys(specs[-1].predictors))
    common = data.dropna(
        subset=all_required + [args.time_col, args.event_col]
    ).copy()

    sample_flow_rows: list[dict] = []
    for spec in specs:
        stage = data.dropna(
            subset=list(spec.predictors) + [args.time_col, args.event_col]
        ).copy()
        sample_flow_rows.append(
            {
                "model_step": spec.step,
                "model": spec.name,
                "added_block": spec.added_block,
                "predictors": ",".join(spec.predictors),
                "n_stage_available": int(stage.shape[0]),
                "n_events_stage_available": int(stage[args.event_col].sum()),
                "n_censored_stage_available": int(stage.shape[0] - stage[args.event_col].sum()),
                "n_common_complete_case": int(common.shape[0]),
                "n_events_common_complete_case": int(common[args.event_col].sum()),
            }
        )

    data.to_csv(
        outdir / f"{args.prefix}_baseline_analysis_dataset.tsv", sep="\t", index=False
    )
    common.to_csv(
        outdir / f"{args.prefix}_common_complete_case_dataset.tsv", sep="\t", index=False
    )
    match_qc.to_csv(
        outdir / f"{args.prefix}_biomarker_matching_qc.tsv", sep="\t", index=False
    )
    timing_summary.to_csv(
        outdir / f"{args.prefix}_first_csf_timing_summary.tsv", sep="\t", index=False
    )
    pd.DataFrame(sample_flow_rows).to_csv(
        outdir / f"{args.prefix}_sample_flow.tsv", sep="\t", index=False
    )

    variables = list(dict.fromkeys(all_required))
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
    missingness.to_csv(
        outdir / f"{args.prefix}_missingness.tsv", sep="\t", index=False
    )

    ok, n, events, censored = sufficient(common, args)
    if not ok:
        raise ValueError(
            "Common complete-case cohort is too small for the requested thresholds: "
            f"N={n}, events={events}, censored={censored}. "
            "Relax thresholds or inspect sample_flow/missingness outputs."
        )

    primary = run_model_sequence(common, specs, args, "common_complete_case")
    result_sets = [primary]
    if args.allow_stage_specific_samples:
        secondary = run_model_sequence(data, specs, args, "stage_specific")
        result_sets.append(secondary)

    for key in ["performance", "coefficients", "oof_predictions", "comparisons", "status"]:
        tables = [r[key] for r in result_sets if not r[key].empty]
        combined = pd.concat(tables, ignore_index=True) if tables else pd.DataFrame()
        combined.to_csv(outdir / f"{args.prefix}_{key}.tsv", sep="\t", index=False)

    summary = {
        "survival_file": str(args.survival_file),
        "adni_file": str(args.adni_file),
        "outcome": "First AD diagnosis after selected baseline-MCI MRI",
        "full_mci2ad_sample_n": int(data.shape[0]),
        "full_mci2ad_sample_events": int(data[args.event_col].sum()),
        "common_complete_case_n": int(common.shape[0]),
        "common_complete_case_events": int(common[args.event_col].sum()),
        "common_complete_case_censored": int(common.shape[0] - common[args.event_col].sum()),
        "baseline_model": parse_list(args.baseline_cols),
        "cumulative_order": [
            parse_list(args.csf_cols)[0] if parse_list(args.csf_cols) else None,
            "+".join(parse_list(args.csf_cols)[1:]),
            args.epoch_col,
        ],
        "csf_selection_rule": (
            "First chronologically available nonmissing measurement per participant and marker; "
            "measurements may occur before or after the selected baseline-MCI MRI; timing and "
            "post-event/censor status are retained in QC outputs."
        ),
        "primary_comparison_rule": (
            "All nested models fit on one common complete-case cohort for paired incremental comparisons."
        ),
        "epoch_training_overlap": (
            "The input MCI-to-AD survival file is expected to exclude original AD EPOCH train/validation participants."
        ),
        "cv_folds_requested": args.cv_folds,
        "bootstrap_requested": args.bootstrap,
        "cox_penalizer_requested": args.cox_penalizer,
        "stage_specific_secondary_enabled": bool(args.allow_stage_specific_samples),
        "primary_performance_file": str(outdir / f"{args.prefix}_performance.tsv"),
        "primary_comparison_file": str(outdir / f"{args.prefix}_comparisons.tsv"),
    }
    with open(outdir / f"{args.prefix}_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    perf = primary["performance"]
    comp = primary["comparisons"]
    log(f"Full MCI-to-AD sample: N={data.shape[0]}, events={int(data[args.event_col].sum())}")
    log(
        f"Common complete-case sample: N={common.shape[0]}, "
        f"events={int(common[args.event_col].sum())}, "
        f"censored={int(common.shape[0] - common[args.event_col].sum())}"
    )
    if not perf.empty:
        show = [c for c in ["model_step", "model", "n", "n_events", "cv_cindex", "status"] if c in perf.columns]
        log("Cumulative model performance:")
        log(perf[show].to_string(index=False))
    if not comp.empty:
        show = [c for c in ["new_step", "new_model", "added_block", "delta_cindex", "delta_ci_lower", "delta_ci_upper", "bootstrap_p_two_sided", "likelihood_ratio_p"] if c in comp.columns]
        log("Sequential comparisons:")
        log(comp[show].to_string(index=False))
    log("============================================================")
    log(f"Results saved to: {outdir}")
    log(f"Performance: {outdir / f'{args.prefix}_performance.tsv'}")
    log(f"Comparisons: {outdir / f'{args.prefix}_comparisons.tsv'}")
    log("============================================================")


if __name__ == "__main__":
    main()
