#!/usr/bin/env python3
# ============================================================
# Re-anchor disease L'EPOCH year scale for problematic clocks
#
# Revised version:
#   In addition to correcting prediction-file metrics:
#     *_clock_acceleration_years
#     *_clock_acceleration_z
#     *_clock_age_years
#
#   this script also overwrites the age coefficient in:
#     *_clock_coefficients.tsv
#
#   and, if present, in:
#     *_clock_nonzero_coefficients.tsv
#
# This is useful because the downstream QC script reads age_beta from
# *_coefficients.tsv to evaluate whether the year-scale transform is stable.
#
# Safety:
#   Before any overwrite, the original TSV is backed up in the same folder.
#
# Core idea:
#   1. Keep the original *_risk_score unchanged.
#   2. Fit a stable unpenalized age-only Cox model in the training split:
#        Surv(time_years, event) ~ age
#   3. Use this stable age beta to re-anchor the year scale:
#        acceleration_years = centered_risk_score / beta_age_ref
#   4. Recompute:
#        *_clock_acceleration_years
#        *_clock_acceleration_z
#        *_clock_age_years
#   5. Overwrite age coefficient in *_coefficients.tsv so future QC
#      uses the stable age anchor.
# ============================================================

import argparse
import datetime
import glob
import json
import os
import re
import shutil
import sys
import warnings
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")


# ============================================================
# 1. Arguments
# ============================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Re-anchor disease L'EPOCH year-scale metrics and overwrite age coefficient."
    )

    parser.add_argument(
        "--base_dir",
        required=True,
        help="Base WholeBodyClock directory containing disease clock folders and QC directories."
    )

    parser.add_argument(
        "--diseases",
        nargs="+",
        default=["asthma", "dementia", "copd", "mi", "stroke"],
        help="Disease labels to process."
    )

    parser.add_argument(
        "--out_dir",
        default=None,
        help="Output directory for re-anchoring summary. Default: <base_dir>/lepoch_year_scale_reanchoring"
    )

    parser.add_argument(
        "--apply",
        action="store_true",
        help="Actually overwrite prediction and coefficient files. Without this flag, runs dry-run only."
    )

    parser.add_argument(
        "--also_test_predictions",
        action="store_true",
        help="Also update *_test_predictions.tsv using the same re-anchoring parameters."
    )

    parser.add_argument(
        "--update_coefficients",
        action="store_true",
        help="Overwrite age coefficient in *_coefficients.tsv."
    )

    parser.add_argument(
        "--update_nonzero_coefficients",
        action="store_true",
        help="Also overwrite age coefficient in *_nonzero_coefficients.tsv if the age row exists."
    )

    parser.add_argument(
        "--backup_suffix",
        default="pre_reanchor_year_scale",
        help="Suffix used for backup files."
    )

    parser.add_argument(
        "--min_events",
        type=int,
        default=20,
        help="Minimum number of events required to fit the reference age Cox model."
    )

    parser.add_argument(
        "--age_beta_min_abs",
        type=float,
        default=0.005,
        help="Minimum acceptable absolute age beta."
    )

    parser.add_argument(
        "--winsorize_score_quantiles",
        nargs=2,
        type=float,
        default=[0.001, 0.999],
        help="Quantiles for winsorizing centered risk score before converting to years."
    )

    parser.add_argument(
        "--no_winsorize",
        action="store_true",
        help="Disable winsorization of centered risk score."
    )

    parser.add_argument(
        "--use_score_calibration_slope",
        action="store_true",
        help="Fit Cox Surv(time,event) ~ centered risk score and multiply score by calibration slope."
    )

    parser.add_argument(
        "--qc_file_kind",
        default="problematic",
        choices=["problematic", "summary"],
        help="Use *_problematic_year_scale_clocks.tsv or *_year_scale_qc_summary.tsv."
    )

    parser.add_argument(
        "--force_process_all_summary_rows",
        action="store_true",
        help="Only relevant with --qc_file_kind summary. Process all rows rather than only FAIL/WARN."
    )

    return parser.parse_args()


# ============================================================
# 2. Generic helpers
# ============================================================

def safe_read_tsv(path: str) -> Optional[pd.DataFrame]:
    if path is None or not os.path.exists(path):
        return None

    try:
        return pd.read_csv(path, sep="\t", low_memory=False)
    except Exception:
        try:
            return pd.read_csv(path, sep=None, engine="python", low_memory=False)
        except Exception as e:
            print(f"WARNING: failed to read {path}: {e}")
            return None


def write_tsv_atomic(df: pd.DataFrame, path: str):
    tmp = path + ".tmp"
    df.to_csv(tmp, sep="\t", index=False)
    os.replace(tmp, path)


def as_numeric(x) -> pd.Series:
    return pd.to_numeric(x, errors="coerce")


def normalize_event(x) -> pd.Series:
    s = pd.Series(x).astype(str).str.strip().str.lower()
    return s.isin(["1", "1.0", "true", "t", "yes", "y"]).astype(int)


def finite_series(x) -> pd.Series:
    y = as_numeric(x)
    y = y.replace([np.inf, -np.inf], np.nan)
    return y


def summarize_numeric(x, prefix: str) -> Dict:
    y = finite_series(x)
    y = y[np.isfinite(y)]

    if len(y) == 0:
        return {
            f"{prefix}_n": 0,
            f"{prefix}_mean": np.nan,
            f"{prefix}_sd": np.nan,
            f"{prefix}_min": np.nan,
            f"{prefix}_p001": np.nan,
            f"{prefix}_p01": np.nan,
            f"{prefix}_p05": np.nan,
            f"{prefix}_median": np.nan,
            f"{prefix}_p95": np.nan,
            f"{prefix}_p99": np.nan,
            f"{prefix}_p999": np.nan,
            f"{prefix}_max": np.nan,
            f"{prefix}_iqr": np.nan,
            f"{prefix}_p01_p99_range": np.nan,
        }

    q = y.quantile([0, 0.001, 0.01, 0.05, 0.5, 0.95, 0.99, 0.999, 1.0])

    return {
        f"{prefix}_n": int(len(y)),
        f"{prefix}_mean": float(y.mean()),
        f"{prefix}_sd": float(y.std(ddof=1)) if len(y) > 1 else np.nan,
        f"{prefix}_min": float(q.loc[0.0]),
        f"{prefix}_p001": float(q.loc[0.001]),
        f"{prefix}_p01": float(q.loc[0.01]),
        f"{prefix}_p05": float(q.loc[0.05]),
        f"{prefix}_median": float(q.loc[0.5]),
        f"{prefix}_p95": float(q.loc[0.95]),
        f"{prefix}_p99": float(q.loc[0.99]),
        f"{prefix}_p999": float(q.loc[0.999]),
        f"{prefix}_max": float(q.loc[1.0]),
        f"{prefix}_iqr": float(y.quantile(0.75) - y.quantile(0.25)),
        f"{prefix}_p01_p99_range": float(q.loc[0.99] - q.loc[0.01]),
    }


def first_existing_col(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def find_col_by_regex(df: pd.DataFrame, patterns: List[str]) -> Optional[str]:
    cols = list(df.columns)
    for pat in patterns:
        hits = [c for c in cols if re.search(pat, str(c), flags=re.IGNORECASE)]
        if len(hits) > 0:
            return hits[0]
    return None


def find_unique_file(clock_dir: str, pattern: str) -> Optional[str]:
    files = sorted(glob.glob(os.path.join(clock_dir, pattern)))
    if len(files) == 0:
        return None
    nonempty = [f for f in files if os.path.getsize(f) > 0]
    return nonempty[0] if len(nonempty) > 0 else files[0]


def backup_file(path: str, backup_suffix: str) -> str:
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_path = f"{path}.{backup_suffix}.{timestamp}.bak"
    shutil.copy2(path, backup_path)
    return backup_path


# ============================================================
# 3. Column detection
# ============================================================

def detect_age_col(df: pd.DataFrame) -> Optional[str]:
    candidates = [
        "age_at_baseline",
        "age_at_imaging",
        "age_at_clock",
        "chronological_age",
        "Age_recruitment",
        "age",
        "Age",
    ]

    hit = first_existing_col(df, candidates)
    if hit is not None:
        return hit

    return find_col_by_regex(
        df,
        patterns=[
            r"^age$",
            r"age_at_baseline",
            r"age_at_imaging",
            r"age_at_clock",
            r"chronological_age",
            r"age.*recruit",
        ],
    )


def detect_risk_score_col(df: pd.DataFrame, disease: str, qc_value: str = "") -> Optional[str]:
    if qc_value and qc_value in df.columns:
        return qc_value

    cols = list(df.columns)

    candidates = [
        c for c in cols
        if re.search(r"risk_score$", str(c), flags=re.IGNORECASE)
    ]

    disease_candidates = [
        c for c in candidates
        if disease.lower() in str(c).lower()
    ]

    if len(disease_candidates) > 0:
        return disease_candidates[0]

    if len(candidates) > 0:
        return candidates[0]

    return find_col_by_regex(
        df,
        patterns=[
            r"risk.*score",
            r"linear.*predictor",
            r"lp$",
        ],
    )


def infer_metric_columns(
    df: pd.DataFrame,
    risk_score_col: str,
    qc_accel_year_col: str = "",
) -> Tuple[str, str, str]:
    if qc_accel_year_col and qc_accel_year_col in df.columns:
        accel_year_col = qc_accel_year_col
    else:
        if risk_score_col.endswith("_risk_score"):
            prefix = risk_score_col.replace("_risk_score", "")
            accel_year_col = prefix + "_clock_acceleration_years"
        else:
            accel_year_col = risk_score_col + "_clock_acceleration_years"

    if "clock_acceleration_years" in accel_year_col:
        accel_z_col = accel_year_col.replace(
            "clock_acceleration_years",
            "clock_acceleration_z"
        )
        clock_age_col = accel_year_col.replace(
            "clock_acceleration_years",
            "clock_age_years"
        )
    elif "acceleration_years" in accel_year_col:
        accel_z_col = accel_year_col.replace(
            "acceleration_years",
            "acceleration_z"
        )
        clock_age_col = accel_year_col.replace(
            "acceleration_years",
            "clock_age_years"
        )
    else:
        accel_z_col = accel_year_col + "_z"
        clock_age_col = accel_year_col + "_clock_age_years"

    return accel_year_col, accel_z_col, clock_age_col


# ============================================================
# 4. Coefficient table utilities
# ============================================================

def detect_coefficient_columns(coef_df: pd.DataFrame) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    term_candidates = [
        "feature",
        "term",
        "variable",
        "covariate",
        "name",
        "Feature",
        "Variable",
        "parameter",
        "coef_name",
    ]

    beta_candidates = [
        "coefficient",
        "coef",
        "beta",
        "Coefficient",
        "Coef",
        "Beta",
        "estimate",
        "value",
    ]

    abs_candidates = [
        "abs_coefficient",
        "abs_coef",
        "abs_beta",
        "absolute_coefficient",
    ]

    term_col = first_existing_col(coef_df, term_candidates)
    beta_col = first_existing_col(coef_df, beta_candidates)
    abs_col = first_existing_col(coef_df, abs_candidates)

    if term_col is None:
        term_col = find_col_by_regex(
            coef_df,
            patterns=[r"feature", r"term", r"variable", r"covariate", r"parameter", r"name"],
        )

    if beta_col is None:
        beta_col = find_col_by_regex(
            coef_df,
            patterns=[r"coef", r"beta", r"estimate", r"value"],
        )

    if abs_col is None:
        abs_col = find_col_by_regex(
            coef_df,
            patterns=[r"abs.*coef", r"abs.*beta", r"absolute.*coef"],
        )

    return term_col, beta_col, abs_col


def age_term_mask(feature_series: pd.Series) -> pd.Series:
    s = feature_series.astype(str)

    patterns = [
        r"^age$",
        r"(^|[^a-zA-Z])age([^a-zA-Z]|$)",
        r"age_at_",
        r"chronological_age",
        r"Age_recruitment",
        r"age_at_baseline",
        r"age_at_imaging",
        r"age_at_clock",
    ]

    mask = pd.Series(False, index=s.index)
    for pat in patterns:
        mask = mask | s.str.contains(pat, case=False, regex=True, na=False)

    return mask


def choose_age_row_index(coef_df: pd.DataFrame, term_col: str) -> Optional[int]:
    mask = age_term_mask(coef_df[term_col])
    candidates = coef_df[mask].copy()

    if candidates.empty:
        return None

    def priority(term: str) -> int:
        t = str(term).lower()

        ordered = [
            r"^num__age_at_baseline$",
            r"^age_at_baseline$",
            r"^num__age_at_imaging$",
            r"^age_at_imaging$",
            r"^num__age_at_clock$",
            r"^age_at_clock$",
            r"^num__chronological_age$",
            r"^chronological_age$",
            r"^num__age$",
            r"^age$",
            r"age_recruitment",
            r"age_at_recruitment",
        ]

        for i, pat in enumerate(ordered):
            if re.search(pat, t):
                return i

        return 100

    candidates["_priority"] = candidates[term_col].apply(priority)
    candidates = candidates.sort_values(["_priority", term_col])

    return int(candidates.index[0])


def update_age_coefficient_file(
    coef_path: Optional[str],
    beta_age_ref: float,
    backup_suffix: str,
    apply_changes: bool,
) -> Dict:
    out = {
        "coefficient_file": coef_path or "",
        "coefficient_file_updated": False,
        "coefficient_file_backup": "",
        "coefficient_age_term": "",
        "coefficient_age_beta_before": np.nan,
        "coefficient_age_beta_after": beta_age_ref,
        "coefficient_update_error": "",
    }

    if coef_path is None or not os.path.exists(coef_path):
        out["coefficient_update_error"] = "coefficient file missing"
        return out

    coef_df = safe_read_tsv(coef_path)
    if coef_df is None or coef_df.empty:
        out["coefficient_update_error"] = "coefficient file unreadable or empty"
        return out

    term_col, beta_col, abs_col = detect_coefficient_columns(coef_df)

    if term_col is None or beta_col is None:
        out["coefficient_update_error"] = "could not detect feature/coefficient columns"
        return out

    idx = choose_age_row_index(coef_df, term_col)

    if idx is None:
        out["coefficient_update_error"] = "no age row found"
        return out

    before = pd.to_numeric(pd.Series([coef_df.loc[idx, beta_col]]), errors="coerce").iloc[0]

    out["coefficient_age_term"] = str(coef_df.loc[idx, term_col])
    out["coefficient_age_beta_before"] = float(before) if np.isfinite(before) else np.nan

    if apply_changes:
        out["coefficient_file_backup"] = backup_file(coef_path, backup_suffix)

        coef_df.loc[idx, beta_col] = beta_age_ref

        if abs_col is not None and abs_col in coef_df.columns:
            coef_df.loc[idx, abs_col] = abs(beta_age_ref)

        if "is_nonzero" in coef_df.columns:
            coef_df.loc[idx, "is_nonzero"] = bool(beta_age_ref != 0)

        write_tsv_atomic(coef_df, coef_path)
        out["coefficient_file_updated"] = True

    return out


# ============================================================
# 5. Cox fitting for stable age beta
# ============================================================

def fit_unpenalized_cox_beta(
    df: pd.DataFrame,
    covariate_col: str,
    min_events: int = 20,
    prefer_train: bool = True,
) -> Tuple[float, str, int, int]:
    try:
        from statsmodels.duration.hazard_regression import PHReg
    except Exception as e:
        raise RuntimeError(
            "statsmodels PHReg is required for unpenalized Cox fitting. "
            f"Import failed: {e}"
        )

    required = ["time_years", "event", covariate_col]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise RuntimeError(f"Missing required columns for Cox fitting: {missing}")

    d = df.copy()
    d[".time"] = finite_series(d["time_years"])
    d[".event"] = normalize_event(d["event"])
    d[".x"] = finite_series(d[covariate_col])

    d = d[
        np.isfinite(d[".time"])
        & np.isfinite(d[".x"])
        & (d[".time"] > 0)
        & d[".event"].isin([0, 1])
    ].copy()

    if d.empty:
        raise RuntimeError("No usable rows for Cox fitting.")

    candidate_sets = []

    if prefer_train and "split" in d.columns:
        train = d[d["split"].astype(str).str.lower() == "train"].copy()
        candidate_sets.append(("train", train))

    candidate_sets.append(("all", d))

    last_error = None

    for source_name, sub in candidate_sets:
        n_events = int(sub[".event"].sum())
        n_used = int(sub.shape[0])

        if n_used < 50 or n_events < min_events:
            last_error = (
                f"{source_name}: insufficient sample/events "
                f"(N={n_used}, events={n_events})"
            )
            continue

        try:
            x = sub[[".x"]].astype(float)
            y = sub[".time"].astype(float)
            status = sub[".event"].astype(int)

            model = PHReg(y, x, status=status)
            fit = model.fit(disp=False)

            beta = float(fit.params[0])

            if not np.isfinite(beta):
                last_error = f"{source_name}: fitted beta is non-finite"
                continue

            return beta, source_name, n_used, n_events

        except Exception as e:
            last_error = f"{source_name}: {e}"
            continue

    raise RuntimeError(f"Could not fit Cox model. Last error: {last_error}")


# ============================================================
# 6. Re-anchoring parameters and metric updates
# ============================================================

def estimate_reanchor_params(
    pred: pd.DataFrame,
    risk_score_col: str,
    age_col: str,
    min_events: int,
    age_beta_min_abs: float,
    winsorize_score_quantiles: Optional[Tuple[float, float]],
    use_score_calibration_slope: bool,
) -> Dict:
    beta_age_ref, beta_source, beta_n, beta_events = fit_unpenalized_cox_beta(
        pred,
        covariate_col=age_col,
        min_events=min_events,
        prefer_train=True,
    )

    if not np.isfinite(beta_age_ref):
        raise RuntimeError("Stable reference age beta is non-finite.")

    if beta_age_ref <= 0:
        raise RuntimeError(f"Stable reference age beta is non-positive: {beta_age_ref}")

    if abs(beta_age_ref) < age_beta_min_abs:
        raise RuntimeError(f"Stable reference age beta is near zero: {beta_age_ref}")

    score = finite_series(pred[risk_score_col])

    if "split" in pred.columns:
        train_mask = pred["split"].astype(str).str.lower() == "train"
        if train_mask.sum() == 0:
            train_mask = pd.Series(True, index=pred.index)
    else:
        train_mask = pd.Series(True, index=pred.index)

    score_center = score[train_mask].median(skipna=True)
    if not np.isfinite(score_center):
        score_center = score.median(skipna=True)

    centered_score = score - score_center

    score_slope = 1.0
    score_slope_source = "not_used"

    if use_score_calibration_slope:
        tmp = pred.copy()
        tmp[".centered_score"] = centered_score

        try:
            score_slope, score_slope_source, _, _ = fit_unpenalized_cox_beta(
                tmp,
                covariate_col=".centered_score",
                min_events=min_events,
                prefer_train=True,
            )

            if not np.isfinite(score_slope) or score_slope <= 0:
                score_slope = 1.0
                score_slope_source = "invalid_or_nonpositive_slope_fallback_to_1"

        except Exception as e:
            score_slope = 1.0
            score_slope_source = f"fit_failed_fallback_to_1: {e}"

    clip_low = np.nan
    clip_high = np.nan

    if winsorize_score_quantiles is not None:
        q_low, q_high = winsorize_score_quantiles
        train_centered = centered_score[train_mask]
        clip_low = float(train_centered.quantile(q_low))
        clip_high = float(train_centered.quantile(q_high))

    score_for_years = centered_score.copy()

    if np.isfinite(clip_low) and np.isfinite(clip_high) and clip_low < clip_high:
        score_for_years = score_for_years.clip(lower=clip_low, upper=clip_high)

    accel_years = (score_for_years * score_slope) / beta_age_ref

    train_accel = accel_years[train_mask]
    z_mean = float(train_accel.mean(skipna=True))
    z_sd = float(train_accel.std(skipna=True, ddof=1))

    if not np.isfinite(z_sd) or z_sd <= 0:
        raise RuntimeError(f"Invalid training SD for acceleration years: {z_sd}")

    return {
        "beta_age_ref": float(beta_age_ref),
        "beta_age_ref_source": beta_source,
        "beta_age_ref_n": int(beta_n),
        "beta_age_ref_events": int(beta_events),
        "risk_score_center": float(score_center),
        "score_calibration_slope": float(score_slope),
        "score_calibration_slope_source": score_slope_source,
        "winsorize_clip_low": clip_low,
        "winsorize_clip_high": clip_high,
        "z_mean_train_accel_years": z_mean,
        "z_sd_train_accel_years": z_sd,
        "age_col": age_col,
        "risk_score_col": risk_score_col,
    }


def compute_reanchored_metrics(
    df: pd.DataFrame,
    params: Dict,
    age_col: str,
    risk_score_col: str,
) -> Tuple[pd.Series, pd.Series, pd.Series]:
    score = finite_series(df[risk_score_col])
    age = finite_series(df[age_col])

    centered_score = score - params["risk_score_center"]

    clip_low = params.get("winsorize_clip_low", np.nan)
    clip_high = params.get("winsorize_clip_high", np.nan)

    score_for_years = centered_score.copy()

    if np.isfinite(clip_low) and np.isfinite(clip_high) and clip_low < clip_high:
        score_for_years = score_for_years.clip(lower=clip_low, upper=clip_high)

    accel_years = (
        score_for_years
        * params["score_calibration_slope"]
        / params["beta_age_ref"]
    )

    accel_z = (
        accel_years - params["z_mean_train_accel_years"]
    ) / params["z_sd_train_accel_years"]

    clock_age_years = age + accel_years

    return accel_years, accel_z, clock_age_years


def update_prediction_file(
    pred_path: str,
    disease: str,
    params: Dict,
    qc_risk_score_col: str,
    qc_accel_year_col: str,
    backup_suffix: str,
    apply_changes: bool,
) -> Dict:
    if pred_path is None or not os.path.exists(pred_path):
        raise RuntimeError(f"Prediction file missing: {pred_path}")

    pred = safe_read_tsv(pred_path)
    if pred is None or pred.empty:
        raise RuntimeError(f"Prediction file is empty/unreadable: {pred_path}")

    risk_score_col = detect_risk_score_col(
        pred,
        disease=disease,
        qc_value=qc_risk_score_col,
    )

    if risk_score_col is None or risk_score_col not in pred.columns:
        raise RuntimeError(f"Could not detect risk score column in {pred_path}")

    age_col = detect_age_col(pred)
    if age_col is None or age_col not in pred.columns:
        raise RuntimeError(f"Could not detect age column in {pred_path}")

    accel_year_col, accel_z_col, clock_age_col = infer_metric_columns(
        pred,
        risk_score_col=risk_score_col,
        qc_accel_year_col=qc_accel_year_col,
    )

    before = {}
    if accel_year_col in pred.columns:
        before.update(summarize_numeric(pred[accel_year_col], "before_accel_years"))
    else:
        before.update(summarize_numeric(pd.Series([], dtype=float), "before_accel_years"))

    accel_years, accel_z, clock_age_years = compute_reanchored_metrics(
        pred,
        params=params,
        age_col=age_col,
        risk_score_col=risk_score_col,
    )

    after = summarize_numeric(accel_years, "after_accel_years")

    backup_path = ""

    if apply_changes:
        backup_path = backup_file(pred_path, backup_suffix)

        pred[accel_year_col] = accel_years
        pred[accel_z_col] = accel_z
        pred[clock_age_col] = clock_age_years

        write_tsv_atomic(pred, pred_path)

    return {
        "prediction_file": pred_path,
        "prediction_file_updated": bool(apply_changes),
        "prediction_file_backup": backup_path,
        "risk_score_col_used": risk_score_col,
        "age_col_used": age_col,
        "accel_year_col_updated": accel_year_col,
        "accel_z_col_updated": accel_z_col,
        "clock_age_col_updated": clock_age_col,
        **before,
        **after,
    }


# ============================================================
# 7. QC file loading
# ============================================================

def get_qc_file_path(base_dir: str, disease: str, qc_file_kind: str) -> str:
    qc_dir = os.path.join(base_dir, f"{disease}_lepoch_year_scale_qc")

    if qc_file_kind == "problematic":
        return os.path.join(qc_dir, f"{disease}_problematic_year_scale_clocks.tsv")

    return os.path.join(qc_dir, f"{disease}_year_scale_qc_summary.tsv")


def load_qc_rows(
    base_dir: str,
    disease: str,
    qc_file_kind: str,
    force_process_all_summary_rows: bool,
) -> pd.DataFrame:
    qc_file = get_qc_file_path(base_dir, disease, qc_file_kind)

    if not os.path.exists(qc_file):
        print(f"WARNING: QC file does not exist for {disease}: {qc_file}")
        return pd.DataFrame()

    qc = safe_read_tsv(qc_file)

    if qc is None or qc.empty:
        print(f"WARNING: QC file is empty for {disease}: {qc_file}")
        return pd.DataFrame()

    qc["disease"] = disease
    qc["qc_file"] = qc_file

    if qc_file_kind == "summary" and not force_process_all_summary_rows:
        if "final_year_scale_qc_status" in qc.columns:
            status = qc["final_year_scale_qc_status"].astype(str)
            qc = qc[
                status.str.startswith("FAIL")
                | status.str.startswith("WARN")
            ].copy()

    return qc


def resolve_clock_dir(base_dir: str, row: pd.Series) -> str:
    folder = str(row.get("clock_folder", "")).strip()
    if folder:
        candidate = os.path.join(base_dir, folder)
        if os.path.isdir(candidate):
            return candidate

    qc_clock_dir = str(row.get("clock_dir", "")).strip()
    if qc_clock_dir and os.path.isdir(qc_clock_dir):
        return qc_clock_dir

    return qc_clock_dir


# ============================================================
# 8. Process one clock
# ============================================================

def process_one_clock(row: pd.Series, args) -> Dict:
    disease = str(row.get("disease", row.get("analysis_label", ""))).strip().lower()

    clock_dir = resolve_clock_dir(args.base_dir, row)
    clock_folder = os.path.basename(clock_dir.rstrip("/"))

    if not os.path.isdir(clock_dir):
        raise RuntimeError(f"Clock directory does not exist: {clock_dir}")

    pred_file = find_unique_file(clock_dir, "*_predictions.tsv")
    test_pred_file = find_unique_file(clock_dir, "*_test_predictions.tsv")
    coeff_file = find_unique_file(clock_dir, "*_coefficients.tsv")
    nonzero_file = find_unique_file(clock_dir, "*_nonzero_coefficients.tsv")

    if pred_file is None:
        raise RuntimeError(f"No *_predictions.tsv found in {clock_dir}")

    qc_risk_score_col = str(row.get("risk_score_col_used", "")).strip()
    qc_accel_year_col = str(row.get("existing_accel_year_col_used", "")).strip()

    pred = safe_read_tsv(pred_file)
    if pred is None or pred.empty:
        raise RuntimeError(f"Main prediction file is empty/unreadable: {pred_file}")

    risk_score_col = detect_risk_score_col(
        pred,
        disease=disease,
        qc_value=qc_risk_score_col,
    )

    if risk_score_col is None or risk_score_col not in pred.columns:
        raise RuntimeError(f"Could not identify risk-score column in {pred_file}")

    age_col = detect_age_col(pred)
    if age_col is None or age_col not in pred.columns:
        raise RuntimeError(f"Could not identify age column in {pred_file}")

    if args.no_winsorize:
        winsorize_quantiles = None
    else:
        winsorize_quantiles = tuple(args.winsorize_score_quantiles)

    params = estimate_reanchor_params(
        pred=pred,
        risk_score_col=risk_score_col,
        age_col=age_col,
        min_events=args.min_events,
        age_beta_min_abs=args.age_beta_min_abs,
        winsorize_score_quantiles=winsorize_quantiles,
        use_score_calibration_slope=args.use_score_calibration_slope,
    )

    main_update = update_prediction_file(
        pred_path=pred_file,
        disease=disease,
        params=params,
        qc_risk_score_col=risk_score_col,
        qc_accel_year_col=qc_accel_year_col,
        backup_suffix=args.backup_suffix,
        apply_changes=args.apply,
    )

    test_update = {}
    if args.also_test_predictions and test_pred_file is not None:
        try:
            test_update_raw = update_prediction_file(
                pred_path=test_pred_file,
                disease=disease,
                params=params,
                qc_risk_score_col=risk_score_col,
                qc_accel_year_col=qc_accel_year_col,
                backup_suffix=args.backup_suffix,
                apply_changes=args.apply,
            )
            test_update = {
                "test_prediction_file": test_update_raw["prediction_file"],
                "test_prediction_file_updated": test_update_raw["prediction_file_updated"],
                "test_prediction_file_backup": test_update_raw["prediction_file_backup"],
            }
        except Exception as e:
            test_update = {
                "test_prediction_file": test_pred_file,
                "test_prediction_file_updated": False,
                "test_prediction_file_backup": "",
                "test_prediction_file_error": str(e),
            }

    coefficient_update = {}
    if args.update_coefficients:
        coefficient_update = update_age_coefficient_file(
            coef_path=coeff_file,
            beta_age_ref=params["beta_age_ref"],
            backup_suffix=args.backup_suffix,
            apply_changes=args.apply,
        )

    nonzero_update = {}
    if args.update_nonzero_coefficients:
        nz_update = update_age_coefficient_file(
            coef_path=nonzero_file,
            beta_age_ref=params["beta_age_ref"],
            backup_suffix=args.backup_suffix,
            apply_changes=args.apply,
        )
        nonzero_update = {
            "nonzero_coefficient_file": nz_update.get("coefficient_file", ""),
            "nonzero_coefficient_file_updated": nz_update.get("coefficient_file_updated", False),
            "nonzero_coefficient_file_backup": nz_update.get("coefficient_file_backup", ""),
            "nonzero_coefficient_age_term": nz_update.get("coefficient_age_term", ""),
            "nonzero_coefficient_age_beta_before": nz_update.get("coefficient_age_beta_before", np.nan),
            "nonzero_coefficient_age_beta_after": nz_update.get("coefficient_age_beta_after", np.nan),
            "nonzero_coefficient_update_error": nz_update.get("coefficient_update_error", ""),
        }

    result = {
        "disease": disease,
        "clock_folder": clock_folder,
        "clock_dir": clock_dir,
        "dry_run": not args.apply,
        "qc_status_original": row.get("final_year_scale_qc_status", ""),
        "qc_reason_original": row.get("final_year_scale_qc_reason", ""),
        "original_age_beta_from_qc": row.get("age_beta", np.nan),
        "delta_cindex_from_qc": row.get("delta_cindex", np.nan),
        "delta_cindex_significant_positive_from_qc": row.get(
            "delta_cindex_significant_positive",
            ""
        ),
        **params,
        **main_update,
        **test_update,
        **coefficient_update,
        **nonzero_update,
        "status": "updated" if args.apply else "dry_run_ok",
        "error": "",
    }

    return result


# ============================================================
# 9. Main
# ============================================================

def main():
    args = parse_args()
    args.base_dir = os.path.abspath(args.base_dir)

    if args.out_dir is None:
        args.out_dir = os.path.join(args.base_dir, "lepoch_year_scale_reanchoring")
    else:
        args.out_dir = os.path.abspath(args.out_dir)

    os.makedirs(args.out_dir, exist_ok=True)

    print("============================================================")
    print("Disease L'EPOCH year-scale re-anchoring")
    print("============================================================")
    print("Base directory:", args.base_dir)
    print("Diseases:", ", ".join(args.diseases))
    print("Output directory:", args.out_dir)
    print("Mode:", "APPLY / overwrite with backup" if args.apply else "DRY RUN / no files overwritten")
    print("Also update test predictions:", args.also_test_predictions)
    print("Update coefficients:", args.update_coefficients)
    print("Update nonzero coefficients:", args.update_nonzero_coefficients)
    print("Winsorization:", "disabled" if args.no_winsorize else args.winsorize_score_quantiles)
    print("Use score calibration slope:", args.use_score_calibration_slope)
    print("============================================================")

    all_qc_rows = []

    for disease in args.diseases:
        qc = load_qc_rows(
            base_dir=args.base_dir,
            disease=disease,
            qc_file_kind=args.qc_file_kind,
            force_process_all_summary_rows=args.force_process_all_summary_rows,
        )

        if not qc.empty:
            all_qc_rows.append(qc)

    if len(all_qc_rows) == 0:
        raise RuntimeError("No QC rows found to process.")

    qc_all = pd.concat(all_qc_rows, ignore_index=True)

    print(f"Total QC rows to process: {qc_all.shape[0]}")

    records = []

    for i, row in qc_all.iterrows():
        disease = str(row.get("disease", row.get("analysis_label", "")))
        clock_folder = str(row.get("clock_folder", ""))

        print(f"\n[{i + 1}/{qc_all.shape[0]}] {disease}: {clock_folder}")

        try:
            rec = process_one_clock(row, args)
            print("  status:", rec["status"])
            print("  beta_age_ref:", rec["beta_age_ref"])
            print("  prediction:", rec["prediction_file"])
            if args.update_coefficients:
                print("  coefficient age term:", rec.get("coefficient_age_term", ""))
                print("  coefficient beta before:", rec.get("coefficient_age_beta_before", ""))
                print("  coefficient beta after:", rec.get("coefficient_age_beta_after", ""))
            if args.apply:
                print("  prediction backup:", rec["prediction_file_backup"])
                if args.update_coefficients:
                    print("  coefficient backup:", rec.get("coefficient_file_backup", ""))

        except Exception as e:
            rec = {
                "disease": disease,
                "clock_folder": clock_folder,
                "clock_dir": row.get("clock_dir", ""),
                "dry_run": not args.apply,
                "status": "failed",
                "error": str(e),
            }
            print("  ERROR:", e)

        records.append(rec)

    summary = pd.DataFrame(records)

    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    mode = "applied" if args.apply else "dry_run"

    summary_path = os.path.join(
        args.out_dir,
        f"lepoch_year_scale_reanchoring_{mode}_{timestamp}.tsv"
    )

    latest_path = os.path.join(
        args.out_dir,
        f"lepoch_year_scale_reanchoring_{mode}_latest.tsv"
    )

    metadata_path = os.path.join(
        args.out_dir,
        f"lepoch_year_scale_reanchoring_{mode}_{timestamp}.json"
    )

    summary.to_csv(summary_path, sep="\t", index=False)
    summary.to_csv(latest_path, sep="\t", index=False)

    metadata = {
        "base_dir": args.base_dir,
        "diseases": args.diseases,
        "out_dir": args.out_dir,
        "apply": args.apply,
        "also_test_predictions": args.also_test_predictions,
        "update_coefficients": args.update_coefficients,
        "update_nonzero_coefficients": args.update_nonzero_coefficients,
        "backup_suffix": args.backup_suffix,
        "min_events": args.min_events,
        "age_beta_min_abs": args.age_beta_min_abs,
        "winsorize_score_quantiles": None if args.no_winsorize else args.winsorize_score_quantiles,
        "use_score_calibration_slope": args.use_score_calibration_slope,
        "qc_file_kind": args.qc_file_kind,
        "n_processed": int(summary.shape[0]),
        "status_counts": summary["status"].value_counts(dropna=False).to_dict(),
        "summary_path": summary_path,
        "latest_path": latest_path,
    }

    with open(metadata_path, "w") as f:
        json.dump(metadata, f, indent=2, default=str)

    print("\n============================================================")
    print("Finished.")
    print("Summary:", summary_path)
    print("Latest summary:", latest_path)
    print("Metadata:", metadata_path)
    print("Status counts:")
    print(summary["status"].value_counts(dropna=False).to_string())
    print("============================================================")

    failed = summary[summary["status"] == "failed"]
    if not failed.empty:
        print("\nFailed clocks:")
        cols = ["disease", "clock_folder", "error"]
        print(failed[cols].to_string(index=False))
        sys.exit(1)


if __name__ == "__main__":
    main()