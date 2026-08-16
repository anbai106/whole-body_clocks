#!/usr/bin/env python3
"""
Post-hoc calibration metrics for the 23 existing mortality EPOCH clocks.

IMPORTANT
---------
This script DOES NOT retrain or refit the original 23 EPOCH models.
It uses the already-saved participant-level *_mortality_clock_predictions.tsv
files and evaluates calibration in the held-out TEST split.

For each mortality EPOCH it reports:

1. Global Cox calibration slope using the original EPOCH linear predictor
   - ideal approximately 1
   - <1: predictions/risk contrasts tend to be too extreme
   - >1: predictions/risk contrasts tend to be too weak

2. At every absolute-risk horizon already present in the prediction file
   (for example risk_5y, risk_10y, risk_15y):
   - mean predicted mortality risk
   - Kaplan-Meier observed mortality risk
   - pointwise 95% CI for observed KM risk when supported by installed sksurv
   - observed / expected (O/E) risk ratio; ideal approximately 1
   - predicted minus observed risk; ideal approximately 0
   - absolute calibration-in-the-large error
   - IPCW Brier score; lower is better
   - null Kaplan-Meier Brier score
   - Brier skill score versus null
   - decile calibration table
   - weighted mean absolute decile calibration error
   - maximum absolute decile calibration error
   - test-set calibration plot

3. Master summary tables across all 23 clocks.

The IPCW Brier score uses train+validation participants to estimate the
censoring distribution. Test participants are evaluated only in the held-out
test split.

The script is deliberately based on saved predictions, so no raw MRI,
proteomics, metabolomics, covariate preprocessing, hyperparameter tuning,
or original Coxnet fitting is rerun.
"""

from __future__ import print_function

import argparse
import json
import math
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from sksurv.linear_model import CoxPHSurvivalAnalysis
    from sksurv.metrics import brier_score, concordance_index_censored
    from sksurv.nonparametric import kaplan_meier_estimator
    from sksurv.util import Surv
except ImportError as exc:
    raise ImportError(
        "This script requires scikit-survival. Run it in the same "
        "'survival_clock' environment used for the mortality EPOCHs."
    ) from exc

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError as exc:
    raise ImportError("This script requires matplotlib.") from exc


# =============================================================================
# 1. Exact 23-clock manifest from the existing WholeBodyClock output folders
# =============================================================================

CLOCKS = [
    {"folder": "adipose_mri_mortality_clock", "clock": "Adipose MRI", "modality": "MRI"},
    {"folder": "brain_mri_mortality_clock", "clock": "Brain MRI", "modality": "MRI"},
    {"folder": "Brain_proteomics_mortality_clock", "clock": "Brain proteomics", "modality": "Proteomics"},
    {"folder": "Digestive_metabolomics_mortality_clock", "clock": "Digestive metabolomics", "modality": "Metabolomics"},
    {"folder": "Endocrine_metabolomics_mortality_clock", "clock": "Endocrine metabolomics", "modality": "Metabolomics"},
    {"folder": "Endocrine_proteomics_mortality_clock", "clock": "Endocrine proteomics", "modality": "Proteomics"},
    {"folder": "Eye_proteomics_mortality_clock", "clock": "Eye proteomics", "modality": "Proteomics"},
    {"folder": "heart_mri_mortality_clock", "clock": "Heart MRI", "modality": "MRI"},
    {"folder": "Heart_proteomics_mortality_clock", "clock": "Heart proteomics", "modality": "Proteomics"},
    {"folder": "Hepatic_metabolomics_mortality_clock", "clock": "Hepatic metabolomics", "modality": "Metabolomics"},
    {"folder": "Hepatic_proteomics_mortality_clock", "clock": "Hepatic proteomics", "modality": "Proteomics"},
    {"folder": "Immune_metabolomics_mortality_clock", "clock": "Immune metabolomics", "modality": "Metabolomics"},
    {"folder": "Immune_proteomics_mortality_clock", "clock": "Immune proteomics", "modality": "Proteomics"},
    {"folder": "kidney_mri_mortality_clock", "clock": "Kidney MRI", "modality": "MRI"},
    {"folder": "liver_mri_mortality_clock", "clock": "Liver MRI", "modality": "MRI"},
    {"folder": "Metabolic_metabolomics_mortality_clock", "clock": "Metabolic metabolomics", "modality": "Metabolomics"},
    {"folder": "pancreas_mri_mortality_clock", "clock": "Pancreas MRI", "modality": "MRI"},
    {"folder": "Pulmonary_proteomics_mortality_clock", "clock": "Pulmonary proteomics", "modality": "Proteomics"},
    {"folder": "Renal_proteomics_mortality_clock", "clock": "Renal proteomics", "modality": "Proteomics"},
    {"folder": "Reproductive_female_proteomics_mortality_clock", "clock": "Reproductive female proteomics", "modality": "Proteomics"},
    {"folder": "Reproductive_male_proteomics_mortality_clock", "clock": "Reproductive male proteomics", "modality": "Proteomics"},
    {"folder": "Skin_proteomics_mortality_clock", "clock": "Skin proteomics", "modality": "Proteomics"},
    {"folder": "spleen_mri_mortality_clock", "clock": "Spleen MRI", "modality": "MRI"},
]


# =============================================================================
# 2. CLI and helpers
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Compute post-hoc held-out calibration metrics for the 23 "
            "mortality EPOCH clocks from existing prediction TSVs."
        )
    )
    p.add_argument(
        "--base-dir",
        default="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock",
    )
    p.add_argument(
        "--output-dir",
        default=None,
        help=(
            "Master output directory. Default: "
            "<base-dir>/mortality_EPOCH_calibration_metrics"
        ),
    )
    p.add_argument(
        "--horizons",
        default="auto",
        help=(
            "'auto' uses all risk_<time>y columns found in each predictions "
            "file. Alternatively provide comma-separated horizons, e.g. 5 "
            "or 5,10,15. Missing horizons are skipped, not reconstructed."
        ),
    )
    p.add_argument("--n-bins", type=int, default=10)
    p.add_argument("--test-split", default="test")
    p.add_argument("--reference-splits", default="train,validation")
    p.add_argument(
        "--minimum-bin-n",
        type=int,
        default=20,
        help="Minimum target size per calibration bin.",
    )
    p.add_argument(
        "--fail-on-missing-clock",
        action="store_true",
        help="Stop instead of recording a failed manifest row.",
    )
    return p.parse_args()


def info(msg):
    print(msg, flush=True)


def warn(msg):
    print("WARNING: {}".format(msg), file=sys.stderr, flush=True)


def sanitize_prefix(x):
    x = re.sub(r"_predictions\.tsv$", "", str(x), flags=re.IGNORECASE)
    x = re.sub(r"[^A-Za-z0-9._-]+", "_", x)
    return x.strip("_")


def parse_event(series):
    if pd.api.types.is_bool_dtype(series):
        return series.astype(bool)

    numeric = pd.to_numeric(series, errors="coerce")
    if numeric.notna().mean() >= 0.90:
        return numeric.fillna(0).astype(int).astype(bool)

    s = series.astype(str).str.strip().str.lower()
    true_vals = {"1", "true", "t", "yes", "y", "event", "dead", "death"}
    false_vals = {"0", "false", "f", "no", "n", "censored", "alive"}

    recognized = s.isin(true_vals | false_vals)
    if recognized.mean() < 0.90:
        bad = sorted(s[~recognized].dropna().unique())[:10]
        raise ValueError(
            "Could not reliably parse event values. Examples: {}".format(bad)
        )

    return s.isin(true_vals)


def make_surv(df):
    return Surv.from_arrays(
        event=df["event_bool"].astype(bool).values,
        time=df["time_years"].astype(float).values,
    )


def resolve_prediction_file(clock_dir):
    hits = sorted(clock_dir.glob("*_mortality_clock_predictions.tsv"))

    # Exact ending excludes the *.before_* backup files shown in some folders.
    hits = [
        x for x in hits
        if ".before_" not in x.name
        and x.name.endswith("_mortality_clock_predictions.tsv")
    ]

    if len(hits) == 1:
        return hits[0]

    if len(hits) == 0:
        raise FileNotFoundError(
            "No primary *_mortality_clock_predictions.tsv found in {}".format(
                clock_dir
            )
        )

    raise RuntimeError(
        "Multiple primary prediction files found in {}: {}".format(
            clock_dir,
            [x.name for x in hits],
        )
    )


def detect_horizons(df, requested):
    found = {}

    for col in df.columns:
        m = re.fullmatch(r"risk_([0-9]+(?:\.[0-9]+)?)y", str(col))
        if m:
            found[float(m.group(1))] = col

    if str(requested).strip().lower() == "auto":
        return dict(sorted(found.items()))

    wanted = [
        float(x.strip())
        for x in str(requested).split(",")
        if x.strip()
    ]

    out = {}
    for h in wanted:
        if h in found:
            out[h] = found[h]
        else:
            warn(
                "Requested horizon {}y is not present as an absolute-risk "
                "column in this clock and will be skipped.".format(h)
            )
    return out


def detect_primary_risk_score(df):
    candidates = [
        str(c)
        for c in df.columns
        if str(c).endswith("_mortality_risk_score")
        and not str(c).startswith("risk_score_")
    ]

    if len(candidates) == 1:
        return candidates[0]

    # Exclude any comparison-model aliases.
    candidates = [
        c for c in candidates
        if "M0_" not in c
        and "M1_" not in c
        and "M2_" not in c
        and "M3_" not in c
    ]

    if len(candidates) == 1:
        return candidates[0]

    if len(candidates) == 0:
        return None

    # Prefer shortest name as a conservative fallback.
    candidates = sorted(candidates, key=lambda x: (len(x), x))
    warn(
        "Multiple primary risk-score candidates found; using {} from {}".format(
            candidates[0], candidates
        )
    )
    return candidates[0]


# =============================================================================
# 3. Calibration estimators
# =============================================================================

def km_risk_with_ci(event, time, horizon):
    event = np.asarray(event, dtype=bool)
    time = np.asarray(time, dtype=float)

    ok = np.isfinite(time) & (time > 0)
    event = event[ok]
    time = time[ok]

    if time.size == 0:
        return np.nan, np.nan, np.nan

    # Newer scikit-survival versions support log-log pointwise CI.
    try:
        km_t, km_s, ci = kaplan_meier_estimator(
            event,
            time,
            conf_type="log-log",
            conf_level=0.95,
        )
        idx = np.searchsorted(km_t, float(horizon), side="right") - 1

        if idx < 0:
            surv = 1.0
            surv_lo = 1.0
            surv_hi = 1.0
        else:
            surv = float(km_s[idx])
            surv_lo = float(ci[0, idx])
            surv_hi = float(ci[1, idx])

        # Risk = 1 - survival, so CI endpoints reverse.
        risk = 1.0 - surv
        risk_lo = 1.0 - surv_hi
        risk_hi = 1.0 - surv_lo

        return float(risk), float(risk_lo), float(risk_hi)

    except TypeError:
        # Compatibility fallback for older versions.
        km_t, km_s = kaplan_meier_estimator(event, time)
        idx = np.searchsorted(km_t, float(horizon), side="right") - 1

        surv = 1.0 if idx < 0 else float(km_s[idx])
        return float(1.0 - surv), np.nan, np.nan


def calibration_slope(test_df, risk_score_col):
    if risk_score_col is None:
        return np.nan

    d = test_df[
        ["time_years", "event_bool", risk_score_col]
    ].copy()

    d[risk_score_col] = pd.to_numeric(
        d[risk_score_col],
        errors="coerce",
    )

    d = d.dropna()
    d = d[
        np.isfinite(d["time_years"])
        & (d["time_years"] > 0)
        & np.isfinite(d[risk_score_col])
    ]

    if d.shape[0] < 20 or int(d["event_bool"].sum()) < 2:
        return np.nan

    X = d[[risk_score_col]].values.astype(float)
    y = make_surv(d)

    try:
        model = CoxPHSurvivalAnalysis(alpha=0.0)
        model.fit(X, y)
    except Exception:
        model = CoxPHSurvivalAnalysis(alpha=1e-8)
        model.fit(X, y)

    return float(np.asarray(model.coef_).reshape(-1)[0])


def test_cindex(test_df, risk_score_col):
    if risk_score_col is None:
        return np.nan

    d = test_df[
        ["time_years", "event_bool", risk_score_col]
    ].copy()

    d[risk_score_col] = pd.to_numeric(
        d[risk_score_col],
        errors="coerce",
    )
    d = d.dropna()
    d = d[
        np.isfinite(d["time_years"])
        & np.isfinite(d[risk_score_col])
    ]

    if d.empty:
        return np.nan

    return float(
        concordance_index_censored(
            d["event_bool"].astype(bool).values,
            d["time_years"].astype(float).values,
            d[risk_score_col].astype(float).values,
        )[0]
    )


def prepare_test_for_brier(reference_df, test_df, horizon):
    """
    scikit-survival requires test follow-up to lie within the follow-up range
    used to estimate censoring. For a fixed earlier horizon, participants whose
    follow-up exceeds the maximum reference follow-up can be administratively
    censored just below that maximum without changing their status at the
    requested horizon.
    """
    ref = reference_df[
        ["time_years", "event_bool"]
    ].dropna().copy()

    tst = test_df[
        ["time_years", "event_bool"]
    ].dropna().copy()

    max_ref = float(ref["time_years"].max())
    max_test = float(tst["time_years"].max())

    if not np.isfinite(max_ref) or max_ref <= 0:
        raise ValueError("Invalid reference maximum follow-up.")

    # Need evaluation horizon strictly inside censoring support.
    eps = max(1e-6, abs(max_ref) * 1e-10)
    cap = max_ref - eps

    if float(horizon) >= cap:
        raise ValueError(
            "Horizon {}y is outside train+validation censoring support "
            "(max reference follow-up {:.4f}y).".format(
                horizon, max_ref
            )
        )

    # Administrative cap only affects observations after the evaluation window.
    over = tst["time_years"] > cap
    if over.any():
        tst.loc[over, "time_years"] = cap
        tst.loc[over, "event_bool"] = False

    return ref, tst, max_ref, max_test


def ipcw_brier(reference_df, test_df, predicted_risk, horizon):
    d = test_df[
        ["time_years", "event_bool"]
    ].copy()
    d["predicted_risk"] = np.asarray(predicted_risk, dtype=float)

    d = d.dropna()
    d = d[
        np.isfinite(d["time_years"])
        & np.isfinite(d["predicted_risk"])
    ].copy()

    if d.empty:
        return np.nan, ""

    try:
        ref, tst, _, _ = prepare_test_for_brier(
            reference_df,
            d,
            horizon,
        )

        y_ref = make_surv(ref)
        y_test = make_surv(tst)

        pred_survival = np.clip(
            1.0 - d["predicted_risk"].values.astype(float),
            0.0,
            1.0,
        )

        _, bs = brier_score(
            y_ref,
            y_test,
            pred_survival.reshape(-1, 1),
            np.asarray([float(horizon)]),
        )

        return float(bs[0]), ""

    except Exception as exc:
        return np.nan, str(exc)


def calibration_bins(test_df, risk_col, horizon, n_bins, minimum_bin_n):
    d = test_df[
        ["participant_id", "time_years", "event_bool", risk_col]
    ].copy()

    d[risk_col] = pd.to_numeric(d[risk_col], errors="coerce")

    d = d.dropna()
    d = d[
        np.isfinite(d["time_years"])
        & (d["time_years"] > 0)
        & np.isfinite(d[risk_col])
    ].copy()

    if d.empty:
        return pd.DataFrame()

    # Prevent excessive binning in smaller subcohorts.
    max_bins_by_n = max(2, int(d.shape[0] // max(1, minimum_bin_n)))
    q = min(
        int(n_bins),
        max_bins_by_n,
        int(d[risk_col].nunique()),
    )

    if q < 2:
        return pd.DataFrame()

    try:
        d["calibration_bin"] = pd.qcut(
            d[risk_col],
            q=q,
            labels=False,
            duplicates="drop",
        )
    except Exception:
        return pd.DataFrame()

    rows = []

    for bin_id, g in d.groupby("calibration_bin", observed=True):
        mean_pred = float(g[risk_col].mean())

        obs, obs_lo, obs_hi = km_risk_with_ci(
            g["event_bool"].values,
            g["time_years"].values,
            horizon,
        )

        rows.append(
            {
                "horizon_years": float(horizon),
                "calibration_bin": int(bin_id) + 1,
                "n": int(g.shape[0]),
                "mean_predicted_risk": mean_pred,
                "observed_risk_km": obs,
                "observed_risk_km_ci_lower": obs_lo,
                "observed_risk_km_ci_upper": obs_hi,
                "predicted_minus_observed": (
                    mean_pred - obs
                    if np.isfinite(obs)
                    else np.nan
                ),
                "absolute_calibration_error": (
                    abs(mean_pred - obs)
                    if np.isfinite(obs)
                    else np.nan
                ),
                "deaths_by_horizon": int(
                    (
                        g["event_bool"]
                        & (g["time_years"] <= float(horizon))
                    ).sum()
                ),
                "n_followed_at_least_to_horizon": int(
                    (g["time_years"] >= float(horizon)).sum()
                ),
            }
        )

    return pd.DataFrame(rows)


def decile_error_metrics(bins):
    if bins.empty:
        return np.nan, np.nan

    d = bins[
        np.isfinite(bins["absolute_calibration_error"])
    ].copy()

    if d.empty:
        return np.nan, np.nan

    weighted_mae = float(
        np.average(
            d["absolute_calibration_error"].values,
            weights=d["n"].values,
        )
    )

    max_abs = float(d["absolute_calibration_error"].max())

    return weighted_mae, max_abs


def save_calibration_plot(bins, clock_name, horizon, png, pdf):
    if bins.empty:
        return

    d = bins[
        np.isfinite(bins["mean_predicted_risk"])
        & np.isfinite(bins["observed_risk_km"])
    ].copy()

    if d.empty:
        return

    max_axis = max(
        float(d["mean_predicted_risk"].max()),
        float(d["observed_risk_km"].max()),
        0.01,
    )
    max_axis = min(1.0, max_axis * 1.12)

    fig, ax = plt.subplots(figsize=(5.4, 5.2))

    ax.plot(
        [0, max_axis],
        [0, max_axis],
        linestyle="--",
        linewidth=1.0,
        color="0.45",
        label="Perfect calibration",
    )

    yerr = None
    if (
        d["observed_risk_km_ci_lower"].notna().all()
        and d["observed_risk_km_ci_upper"].notna().all()
    ):
        lower = (
            d["observed_risk_km"]
            - d["observed_risk_km_ci_lower"]
        ).clip(lower=0)
        upper = (
            d["observed_risk_km_ci_upper"]
            - d["observed_risk_km"]
        ).clip(lower=0)
        yerr = np.vstack([lower.values, upper.values])

    ax.errorbar(
        d["mean_predicted_risk"].values,
        d["observed_risk_km"].values,
        yerr=yerr,
        marker="o",
        linewidth=1.4,
        capsize=2.5 if yerr is not None else 0,
        label="Held-out test set",
    )

    for _, r in d.iterrows():
        ax.annotate(
            str(int(r["calibration_bin"])),
            (
                float(r["mean_predicted_risk"]),
                float(r["observed_risk_km"]),
            ),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=8,
        )

    ax.set_xlim(0, max_axis)
    ax.set_ylim(0, max_axis)
    ax.set_xlabel("Mean predicted mortality risk")
    ax.set_ylabel("Kaplan-Meier observed mortality risk")
    ax.set_title(
        "{} mortality EPOCH: {:g}-year calibration".format(
            clock_name,
            float(horizon),
        )
    )
    ax.legend(frameon=False)
    fig.tight_layout()

    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# 4. One-clock analysis
# =============================================================================

def analyze_clock(clock_info, base_dir, master_dir, args):
    clock_dir = base_dir / clock_info["folder"]

    if not clock_dir.exists():
        raise FileNotFoundError(
            "Clock directory does not exist: {}".format(clock_dir)
        )

    pred_file = resolve_prediction_file(clock_dir)
    prefix = sanitize_prefix(pred_file.name)

    per_clock_dir = clock_dir / "calibration_metrics"
    per_clock_dir.mkdir(parents=True, exist_ok=True)

    info("")
    info("=" * 92)
    info("Clock: {}".format(clock_info["clock"]))
    info("Modality: {}".format(clock_info["modality"]))
    info("Predictions: {}".format(pred_file))
    info("=" * 92)

    df = pd.read_csv(
        pred_file,
        sep="\t",
        low_memory=False,
    )

    required = {
        "participant_id",
        "split",
        "time_years",
        "event",
    }
    missing = sorted(required - set(df.columns))

    if missing:
        raise ValueError(
            "Missing required prediction columns: {}".format(missing)
        )

    df["participant_id"] = df["participant_id"].astype(str)
    df["split"] = df["split"].astype(str).str.lower()
    df["time_years"] = pd.to_numeric(
        df["time_years"],
        errors="coerce",
    )
    df["event_bool"] = parse_event(df["event"])

    test_name = str(args.test_split).lower()
    reference_names = [
        x.strip().lower()
        for x in str(args.reference_splits).split(",")
        if x.strip()
    ]

    test = df[
        df["split"] == test_name
    ].copy()

    reference = df[
        df["split"].isin(reference_names)
    ].copy()

    test = test[
        np.isfinite(test["time_years"])
        & (test["time_years"] > 0)
    ].copy()

    reference = reference[
        np.isfinite(reference["time_years"])
        & (reference["time_years"] > 0)
    ].copy()

    if test.empty:
        raise ValueError("No held-out test rows found.")

    if reference.empty:
        raise ValueError(
            "No train/validation rows found. Use the full predictions.tsv "
            "rather than only test_predictions.tsv."
        )

    risk_score_col = detect_primary_risk_score(df)
    horizon_cols = detect_horizons(df, args.horizons)

    if not horizon_cols:
        raise ValueError(
            "No usable absolute-risk columns such as risk_5y were found."
        )

    cindex_value = test_cindex(test, risk_score_col)
    slope_value = calibration_slope(test, risk_score_col)

    slope_row = {
        "clock": clock_info["clock"],
        "modality": clock_info["modality"],
        "folder": clock_info["folder"],
        "prediction_file": str(pred_file),
        "risk_score_col": risk_score_col,
        "n_test": int(test.shape[0]),
        "n_events_test": int(test["event_bool"].sum()),
        "median_followup_test_years": float(
            test["time_years"].median()
        ),
        "max_followup_test_years": float(
            test["time_years"].max()
        ),
        "test_cindex_from_saved_risk_score": cindex_value,
        "cox_calibration_slope": slope_value,
        "ideal_calibration_slope": 1.0,
        "calibration_slope_interpretation": (
            "Ideal approximately 1; values below 1 suggest overly extreme "
            "risk contrasts, whereas values above 1 suggest risk contrasts "
            "that are too weak."
        ),
    }

    pd.DataFrame([slope_row]).to_csv(
        per_clock_dir / "{}_calibration_slope.tsv".format(prefix),
        sep="\t",
        index=False,
        na_rep="NA",
    )

    horizon_rows = []
    all_bins = []

    for horizon, risk_col in horizon_cols.items():
        d = test[
            [
                "participant_id",
                "time_years",
                "event_bool",
                risk_col,
            ]
        ].copy()

        d[risk_col] = pd.to_numeric(
            d[risk_col],
            errors="coerce",
        )

        d = d[
            np.isfinite(d["time_years"])
            & np.isfinite(d[risk_col])
        ].copy()

        if d.empty:
            warn(
                "{}: no usable predictions at {}y".format(
                    clock_info["clock"],
                    horizon,
                )
            )
            continue

        pred = d[risk_col].astype(float).values

        mean_pred = float(np.mean(pred))

        obs, obs_lo, obs_hi = km_risk_with_ci(
            d["event_bool"].values,
            d["time_years"].values,
            horizon,
        )

        oe = (
            float(obs / mean_pred)
            if np.isfinite(obs)
            and np.isfinite(mean_pred)
            and mean_pred > 0
            else np.nan
        )

        pred_minus_obs = (
            float(mean_pred - obs)
            if np.isfinite(obs)
            else np.nan
        )

        abs_cil_error = (
            float(abs(mean_pred - obs))
            if np.isfinite(obs)
            else np.nan
        )

        bs_model, bs_model_error = ipcw_brier(
            reference,
            d,
            pred,
            horizon,
        )

        # Null model: everyone receives the train+validation KM mortality risk.
        null_risk, _, _ = km_risk_with_ci(
            reference["event_bool"].values,
            reference["time_years"].values,
            horizon,
        )

        null_pred = np.repeat(
            null_risk,
            d.shape[0],
        )

        bs_null, bs_null_error = ipcw_brier(
            reference,
            d,
            null_pred,
            horizon,
        )

        brier_skill = (
            float(1.0 - bs_model / bs_null)
            if np.isfinite(bs_model)
            and np.isfinite(bs_null)
            and bs_null > 0
            else np.nan
        )

        bins = calibration_bins(
            d,
            risk_col,
            horizon,
            n_bins=args.n_bins,
            minimum_bin_n=args.minimum_bin_n,
        )

        weighted_mae, max_abs_error = decile_error_metrics(bins)

        label = "{:g}".format(float(horizon))

        if not bins.empty:
            bins.insert(0, "clock", clock_info["clock"])
            bins.insert(1, "modality", clock_info["modality"])
            bins.insert(2, "folder", clock_info["folder"])
            bins.insert(3, "risk_column", risk_col)

            bins.to_csv(
                per_clock_dir
                / "{}_calibration_bins_{}y.tsv".format(
                    prefix,
                    label,
                ),
                sep="\t",
                index=False,
                na_rep="NA",
            )

            all_bins.append(bins)

            save_calibration_plot(
                bins,
                clock_name=clock_info["clock"],
                horizon=horizon,
                png=(
                    per_clock_dir
                    / "{}_calibration_plot_{}y.png".format(
                        prefix,
                        label,
                    )
                ),
                pdf=(
                    per_clock_dir
                    / "{}_calibration_plot_{}y.pdf".format(
                        prefix,
                        label,
                    )
                ),
            )

        row = {
            "clock": clock_info["clock"],
            "modality": clock_info["modality"],
            "folder": clock_info["folder"],
            "prediction_file": str(pred_file),
            "risk_score_col": risk_score_col,
            "risk_column": risk_col,
            "horizon_years": float(horizon),
            "n_test": int(test.shape[0]),
            "n_test_with_absolute_risk": int(d.shape[0]),
            "absolute_risk_prediction_coverage": float(
                d.shape[0] / test.shape[0]
            ),
            "n_events_test_total": int(test["event_bool"].sum()),
            "deaths_by_horizon": int(
                (
                    d["event_bool"]
                    & (d["time_years"] <= float(horizon))
                ).sum()
            ),
            "n_followed_at_least_to_horizon": int(
                (d["time_years"] >= float(horizon)).sum()
            ),
            "fraction_followed_at_least_to_horizon": float(
                (d["time_years"] >= float(horizon)).mean()
            ),
            "median_followup_test_years": float(
                d["time_years"].median()
            ),
            "max_followup_test_years": float(
                d["time_years"].max()
            ),
            "mean_predicted_mortality_risk": mean_pred,
            "observed_mortality_risk_km": obs,
            "observed_mortality_risk_km_ci_lower": obs_lo,
            "observed_mortality_risk_km_ci_upper": obs_hi,
            "observed_to_expected_ratio": oe,
            "ideal_observed_to_expected_ratio": 1.0,
            "predicted_minus_observed_risk": pred_minus_obs,
            "absolute_calibration_in_the_large_error": abs_cil_error,
            "ideal_predicted_minus_observed_risk": 0.0,
            "ipcw_brier_score": bs_model,
            "ipcw_brier_score_null_km": bs_null,
            "brier_skill_score_vs_null": brier_skill,
            "weighted_mean_absolute_decile_calibration_error": weighted_mae,
            "maximum_absolute_decile_calibration_error": max_abs_error,
            "n_calibration_bins": int(
                bins["calibration_bin"].nunique()
            ) if not bins.empty else 0,
            "test_cindex_from_saved_risk_score": cindex_value,
            "cox_calibration_slope": slope_value,
            "ideal_calibration_slope": 1.0,
            "brier_model_error_if_any": bs_model_error,
            "brier_null_error_if_any": bs_null_error,
        }

        horizon_rows.append(row)

        info(
            "{:>4s}y | N={:>6,d} | pred={:.4f} | KM obs={:.4f} | "
            "O/E={:.3f} | slope={} | Brier={}".format(
                label,
                d.shape[0],
                mean_pred,
                obs if np.isfinite(obs) else np.nan,
                oe if np.isfinite(oe) else np.nan,
                "{:.3f}".format(slope_value)
                if np.isfinite(slope_value)
                else "NA",
                "{:.5f}".format(bs_model)
                if np.isfinite(bs_model)
                else "NA",
            )
        )

    if not horizon_rows:
        raise RuntimeError(
            "No calibration horizons were successfully evaluated."
        )

    horizon_df = pd.DataFrame(horizon_rows)

    horizon_df.to_csv(
        per_clock_dir / "{}_calibration_metrics.tsv".format(prefix),
        sep="\t",
        index=False,
        na_rep="NA",
    )

    if all_bins:
        pd.concat(
            all_bins,
            ignore_index=True,
        ).to_csv(
            per_clock_dir
            / "{}_calibration_bins_all_horizons.tsv".format(prefix),
            sep="\t",
            index=False,
            na_rep="NA",
        )

    compact = {
        "clock": clock_info["clock"],
        "modality": clock_info["modality"],
        "folder": clock_info["folder"],
        "prediction_file": str(pred_file),
        "risk_score_col": risk_score_col,
        "n_test": int(test.shape[0]),
        "n_events_test": int(test["event_bool"].sum()),
        "test_cindex_from_saved_risk_score": (
            cindex_value
            if np.isfinite(cindex_value)
            else None
        ),
        "cox_calibration_slope": (
            slope_value
            if np.isfinite(slope_value)
            else None
        ),
        "calibration_by_horizon": json.loads(
            horizon_df.to_json(orient="records")
        ),
    }

    with open(
        per_clock_dir / "{}_calibration_metrics.json".format(prefix),
        "w",
    ) as f:
        json.dump(compact, f, indent=2)

    return horizon_df, pd.DataFrame([slope_row]), {
        "clock": clock_info["clock"],
        "modality": clock_info["modality"],
        "folder": clock_info["folder"],
        "status": "success",
        "prediction_file": str(pred_file),
        "n_test": int(test.shape[0]),
        "n_events_test": int(test["event_bool"].sum()),
        "horizons_evaluated": ",".join(
            "{:g}".format(float(x))
            for x in sorted(horizon_df["horizon_years"].unique())
        ),
        "output_dir": str(per_clock_dir),
        "error": "",
    }


# =============================================================================
# 5. Batch all 23 clocks and aggregate
# =============================================================================

def main():
    args = parse_args()

    base_dir = Path(args.base_dir).resolve()

    if args.output_dir:
        master_dir = Path(args.output_dir).resolve()
    else:
        master_dir = base_dir / "mortality_EPOCH_calibration_metrics"

    master_dir.mkdir(parents=True, exist_ok=True)

    if len(CLOCKS) != 23:
        raise RuntimeError(
            "Internal manifest does not contain exactly 23 clocks."
        )

    all_horizon = []
    all_slope = []
    manifest = []

    info("=" * 92)
    info("POST-HOC CALIBRATION OF 23 MORTALITY EPOCH CLOCKS")
    info("=" * 92)
    info("Base directory: {}".format(base_dir))
    info("Master output: {}".format(master_dir))
    info("Horizon rule: {}".format(args.horizons))
    info("Primary evaluation split: {}".format(args.test_split))
    info("IPCW reference splits: {}".format(args.reference_splits))
    info("=" * 92)

    for i, clock_info in enumerate(CLOCKS, start=1):
        info("")
        info("[{:02d}/23] {}".format(i, clock_info["clock"]))

        try:
            h_df, s_df, m = analyze_clock(
                clock_info,
                base_dir,
                master_dir,
                args,
            )
            all_horizon.append(h_df)
            all_slope.append(s_df)
            manifest.append(m)

        except Exception as exc:
            msg = "{}".format(exc)
            warn(
                "{} FAILED: {}".format(
                    clock_info["clock"],
                    msg,
                )
            )
            manifest.append(
                {
                    "clock": clock_info["clock"],
                    "modality": clock_info["modality"],
                    "folder": clock_info["folder"],
                    "status": "failed",
                    "prediction_file": "",
                    "n_test": np.nan,
                    "n_events_test": np.nan,
                    "horizons_evaluated": "",
                    "output_dir": "",
                    "error": msg,
                }
            )

            if args.fail_on_missing_clock:
                raise

    manifest_df = pd.DataFrame(manifest)

    manifest_df.to_csv(
        master_dir / "mortality_EPOCH_23_calibration_run_manifest.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )

    if not all_horizon:
        raise RuntimeError(
            "Calibration failed for all clocks. See run manifest."
        )

    horizon_master = pd.concat(
        all_horizon,
        ignore_index=True,
    )

    slope_master = pd.concat(
        all_slope,
        ignore_index=True,
    )

    horizon_master.to_csv(
        master_dir / "mortality_EPOCH_23_calibration_by_horizon.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )

    slope_master.to_csv(
        master_dir / "mortality_EPOCH_23_calibration_slope.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )

    # Common 5-year table for manuscript-wide comparison.
    five_year = horizon_master[
        np.isclose(
            horizon_master["horizon_years"].astype(float),
            5.0,
        )
    ].copy()

    if not five_year.empty:
        five_year.to_csv(
            master_dir / "mortality_EPOCH_23_calibration_5y.tsv",
            sep="\t",
            index=False,
            na_rep="NA",
        )

    # Compact manuscript-oriented table.
    report_cols = [
        "clock",
        "modality",
        "horizon_years",
        "n_test",
        "n_events_test_total",
        "mean_predicted_mortality_risk",
        "observed_mortality_risk_km",
        "observed_mortality_risk_km_ci_lower",
        "observed_mortality_risk_km_ci_upper",
        "observed_to_expected_ratio",
        "predicted_minus_observed_risk",
        "absolute_calibration_in_the_large_error",
        "cox_calibration_slope",
        "ipcw_brier_score",
        "ipcw_brier_score_null_km",
        "brier_skill_score_vs_null",
        "weighted_mean_absolute_decile_calibration_error",
        "maximum_absolute_decile_calibration_error",
        "fraction_followed_at_least_to_horizon",
        "test_cindex_from_saved_risk_score",
    ]

    report = horizon_master[
        [c for c in report_cols if c in horizon_master.columns]
    ].copy()

    report.to_csv(
        master_dir / "mortality_EPOCH_23_calibration_manuscript_table.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )

    n_success = int((manifest_df["status"] == "success").sum())
    n_failed = int((manifest_df["status"] == "failed").sum())

    info("")
    info("=" * 92)
    info("CALIBRATION ANALYSIS FINISHED")
    info("=" * 92)
    info("Successful clocks: {}/23".format(n_success))
    info("Failed clocks: {}".format(n_failed))
    info("Master files:")
    info(
        "  {}".format(
            master_dir
            / "mortality_EPOCH_23_calibration_by_horizon.tsv"
        )
    )
    info(
        "  {}".format(
            master_dir
            / "mortality_EPOCH_23_calibration_slope.tsv"
        )
    )
    if not five_year.empty:
        info(
            "  {}".format(
                master_dir
                / "mortality_EPOCH_23_calibration_5y.tsv"
            )
        )
    info(
        "  {}".format(
            master_dir
            / "mortality_EPOCH_23_calibration_manuscript_table.tsv"
        )
    )
    info(
        "  {}".format(
            master_dir
            / "mortality_EPOCH_23_calibration_run_manifest.tsv"
        )
    )
    info("=" * 92)


if __name__ == "__main__":
    main()
