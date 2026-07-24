#!/usr/bin/env python3
"""
STEP 3: Validate and summarize the MHAS phenotype-based mortality EPOCH clock

Inputs from Step 2:
  /Users/hao/Dropbox/MHAS/step2_mortality_epoch_model/
    mhas_mortality_epoch_predictions.tsv
    mhas_mortality_epoch_performance.tsv

Optional input from Step 1 for adjusted Cox models:
  /Users/hao/Dropbox/MHAS/step1_2001_mortality_epoch/
    mhas_2001_step1_model_input_primary_nondisease.tsv

Main outputs:
  /Users/hao/Dropbox/MHAS/step3_mortality_epoch_validation/
    mhas_step3_cindex_bootstrap_ci.tsv
    mhas_step3_cox_hr_per_sd.tsv
    mhas_step3_cox_hr_by_quartile.tsv
    mhas_step3_km_mortality_by_quartile.tsv
    mhas_step3_calibration_by_quartile_horizon.tsv
    mhas_step3_validation_summary.tsv
    mhas_step3_fig_km_test.pdf/png
    mhas_step3_fig_calibration_test_10yr.pdf/png
    mhas_step3_fig_cindex_summary.pdf/png
    mhas_step3_audit.txt

Purpose:
  Step 1 builds the analytic cohort.
  Step 2 trains the mortality EPOCH clock.
  Step 3 validates the trained score and creates manuscript-ready summary tables:
    1) bootstrap C-index confidence intervals,
    2) Cox HRs per 1-SD higher EPOCH acceleration,
    3) Cox HRs across EPOCH quartiles,
    4) Kaplan-Meier mortality curves by quartile,
    5) observed calibration by quartile at fixed time horizons.

Primary reporting:
  The test split should be treated as the main internal hold-out evaluation.
"""

import argparse
import json
import math
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

try:
    from lifelines import CoxPHFitter, KaplanMeierFitter
    from lifelines.statistics import logrank_test, multivariate_logrank_test
    from lifelines.utils import concordance_index
except Exception as e:
    raise ImportError(
        "This script requires lifelines. Install it with:\n"
        "  pip install lifelines\n"
        "or:\n"
        "  conda install -c conda-forge lifelines"
    ) from e


# -----------------------------
# Utility functions
# -----------------------------
def log(msg: str) -> None:
    print(msg, flush=True)


def safe_mkdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def clean_colname(x: str) -> str:
    return (
        str(x)
        .replace(" ", "_")
        .replace("/", "_")
        .replace(":", "_")
        .replace("(", "")
        .replace(")", "")
        .replace("[", "")
        .replace("]", "")
        .replace(",", "")
    )


def harrell_cindex(time, event, score):
    """
    Larger score = higher risk. lifelines concordance_index expects larger predicted
    value = longer survival, so pass -score.
    """
    time = pd.Series(time).astype(float)
    event = pd.Series(event).astype(int)
    score = pd.Series(score).astype(float)
    mask = time.notna() & event.notna() & score.notna()
    if int(mask.sum()) < 10 or event[mask].nunique() < 2:
        return np.nan
    return float(concordance_index(time[mask], -score[mask], event[mask]))


def bootstrap_cindex(df, score_col, split_name, n_boot=1000, seed=20260721):
    """
    Bootstrap C-index CI within a split.
    Resamples participants with replacement. Skips bootstrap samples without
    at least one event and one non-event.
    """
    rng = np.random.default_rng(seed)
    sub = df.copy()
    sub = sub[
        sub["followup_years"].notna()
        & sub["event_death"].notna()
        & sub[score_col].notna()
    ].copy()

    if len(sub) < 20 or sub["event_death"].nunique() < 2:
        return {
            "split": split_name,
            "score": score_col,
            "n": len(sub),
            "deaths": int(sub["event_death"].sum()) if len(sub) else 0,
            "cindex": np.nan,
            "ci_lower": np.nan,
            "ci_upper": np.nan,
            "n_boot_success": 0,
        }

    c_obs = harrell_cindex(sub["followup_years"], sub["event_death"], sub[score_col])

    vals = []
    n = len(sub)
    arr_idx = np.arange(n)
    for _ in range(n_boot):
        idx = rng.choice(arr_idx, size=n, replace=True)
        b = sub.iloc[idx]
        if b["event_death"].nunique() < 2:
            continue
        val = harrell_cindex(b["followup_years"], b["event_death"], b[score_col])
        if pd.notna(val):
            vals.append(val)

    if len(vals) < max(20, n_boot * 0.1):
        ci_l = np.nan
        ci_u = np.nan
    else:
        ci_l, ci_u = np.percentile(vals, [2.5, 97.5])

    return {
        "split": split_name,
        "score": score_col,
        "n": int(n),
        "deaths": int(sub["event_death"].sum()),
        "cindex": c_obs,
        "ci_lower": float(ci_l) if pd.notna(ci_l) else np.nan,
        "ci_upper": float(ci_u) if pd.notna(ci_u) else np.nan,
        "n_boot_success": int(len(vals)),
    }


def get_adjustment_covariates(step1_df):
    """
    Use age and sex dummy columns if present in Step 1 model input.
    These match the residualization used in Step 2 and provide adjusted HRs.
    """
    covars = []
    if "age_2001" in step1_df.columns:
        covars.append("age_2001")
    sex_cols = [c for c in step1_df.columns if c.startswith("sex_")]
    covars += sex_cols
    # Avoid all-zero or constant columns later when fitting.
    return covars


def prepare_cox_df(df, cols):
    out = df[cols].copy()
    for c in cols:
        if c not in ["mortality_epoch_quartile"]:
            out[c] = pd.to_numeric(out[c], errors="coerce")
    out = out.replace([np.inf, -np.inf], np.nan).dropna()
    return out


def fit_cox_per_sd(df, score_col, split_name, adjustment_covars=None):
    """
    Cox HR per 1-SD score increase. Uses the supplied score as is if already z-scored.
    For non-z scores, standardizes within the analyzed data.
    """
    adjustment_covars = adjustment_covars or []
    cols = ["followup_years", "event_death", score_col] + adjustment_covars
    sub = prepare_cox_df(df, cols)

    if len(sub) < 20 or sub["event_death"].nunique() < 2:
        return {
            "split": split_name,
            "score": score_col,
            "adjusted": bool(adjustment_covars),
            "n": len(sub),
            "deaths": int(sub["event_death"].sum()) if len(sub) else 0,
            "hr": np.nan,
            "ci_lower": np.nan,
            "ci_upper": np.nan,
            "p": np.nan,
            "status": "too_few_events_or_rows",
        }

    score_sd = sub[score_col].std(ddof=0)
    if score_sd == 0 or pd.isna(score_sd):
        return {
            "split": split_name,
            "score": score_col,
            "adjusted": bool(adjustment_covars),
            "n": len(sub),
            "deaths": int(sub["event_death"].sum()),
            "hr": np.nan,
            "ci_lower": np.nan,
            "ci_upper": np.nan,
            "p": np.nan,
            "status": "zero_score_sd",
        }

    model_score = f"{score_col}__per_sd"
    sub[model_score] = (sub[score_col] - sub[score_col].mean()) / score_sd

    model_cols = ["followup_years", "event_death", model_score] + adjustment_covars
    # Drop constant covariates.
    keep_covars = []
    for c in adjustment_covars:
        if sub[c].nunique(dropna=True) > 1:
            keep_covars.append(c)
    model_cols = ["followup_years", "event_death", model_score] + keep_covars

    try:
        cph = CoxPHFitter(penalizer=0.001)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cph.fit(
                sub[model_cols],
                duration_col="followup_years",
                event_col="event_death",
                show_progress=False,
            )

        s = cph.summary.loc[model_score]
        return {
            "split": split_name,
            "score": score_col,
            "adjusted": bool(keep_covars),
            "adjustment_covariates": ",".join(keep_covars),
            "n": int(len(sub)),
            "deaths": int(sub["event_death"].sum()),
            "hr": float(np.exp(s["coef"])),
            "ci_lower": float(np.exp(s["coef lower 95%"])),
            "ci_upper": float(np.exp(s["coef upper 95%"])),
            "p": float(s["p"]),
            "status": "ok",
        }
    except Exception as e:
        return {
            "split": split_name,
            "score": score_col,
            "adjusted": bool(adjustment_covars),
            "adjustment_covariates": ",".join(adjustment_covars),
            "n": int(len(sub)),
            "deaths": int(sub["event_death"].sum()),
            "hr": np.nan,
            "ci_lower": np.nan,
            "ci_upper": np.nan,
            "p": np.nan,
            "status": f"error: {str(e)[:300]}",
        }


def fit_cox_quartiles(df, split_name, adjustment_covars=None):
    """
    Cox model with EPOCH quartiles as categorical predictors, Q1 as reference.
    """
    adjustment_covars = adjustment_covars or []
    needed = ["followup_years", "event_death", "mortality_epoch_quartile"] + adjustment_covars
    sub = df[needed].copy()
    sub["followup_years"] = pd.to_numeric(sub["followup_years"], errors="coerce")
    sub["event_death"] = pd.to_numeric(sub["event_death"], errors="coerce")
    sub = sub.dropna(subset=["followup_years", "event_death", "mortality_epoch_quartile"])
    sub["event_death"] = sub["event_death"].astype(int)

    if len(sub) < 20 or sub["event_death"].nunique() < 2:
        return pd.DataFrame([{
            "split": split_name,
            "comparison": "quartiles",
            "n": len(sub),
            "deaths": int(sub["event_death"].sum()) if len(sub) else 0,
            "hr": np.nan,
            "ci_lower": np.nan,
            "ci_upper": np.nan,
            "p": np.nan,
            "status": "too_few_events_or_rows",
        }])

    # Make dummy variables with Q1_lowest as reference.
    q = sub["mortality_epoch_quartile"].astype(str)
    dummies = pd.get_dummies(q, prefix="quartile", dtype=int)
    for col in ["quartile_Q2", "quartile_Q3", "quartile_Q4_highest"]:
        if col not in dummies.columns:
            dummies[col] = 0
    dummies = dummies[["quartile_Q2", "quartile_Q3", "quartile_Q4_highest"]]

    model_df = pd.concat(
        [sub[["followup_years", "event_death"]].reset_index(drop=True), dummies.reset_index(drop=True)],
        axis=1,
    )

    keep_covars = []
    for c in adjustment_covars:
        x = pd.to_numeric(sub[c], errors="coerce")
        if x.notna().sum() > 0:
            x = x.fillna(x.median())
            if x.nunique(dropna=True) > 1:
                model_df[c] = x.values
                keep_covars.append(c)

    rows = []
    try:
        cph = CoxPHFitter(penalizer=0.001)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cph.fit(model_df, duration_col="followup_years", event_col="event_death", show_progress=False)

        for col, label in [
            ("quartile_Q2", "Q2 vs Q1_lowest"),
            ("quartile_Q3", "Q3 vs Q1_lowest"),
            ("quartile_Q4_highest", "Q4_highest vs Q1_lowest"),
        ]:
            if col in cph.summary.index:
                s = cph.summary.loc[col]
                rows.append({
                    "split": split_name,
                    "comparison": label,
                    "adjusted": bool(keep_covars),
                    "adjustment_covariates": ",".join(keep_covars),
                    "n": int(len(model_df)),
                    "deaths": int(model_df["event_death"].sum()),
                    "hr": float(np.exp(s["coef"])),
                    "ci_lower": float(np.exp(s["coef lower 95%"])),
                    "ci_upper": float(np.exp(s["coef upper 95%"])),
                    "p": float(s["p"]),
                    "status": "ok",
                })

    except Exception as e:
        rows.append({
            "split": split_name,
            "comparison": "quartiles",
            "adjusted": bool(keep_covars),
            "adjustment_covariates": ",".join(keep_covars),
            "n": int(len(model_df)),
            "deaths": int(model_df["event_death"].sum()),
            "hr": np.nan,
            "ci_lower": np.nan,
            "ci_upper": np.nan,
            "p": np.nan,
            "status": f"error: {str(e)[:300]}",
        })

    return pd.DataFrame(rows)


def observed_km_risk_at_horizon(time, event, horizon):
    """
    Observed cumulative mortality by Kaplan-Meier at a fixed time horizon.
    Returns 1 - S(horizon). If no observations, returns NaN.
    """
    time = pd.Series(time).astype(float)
    event = pd.Series(event).astype(int)
    mask = time.notna() & event.notna()
    if mask.sum() < 5:
        return np.nan
    kmf = KaplanMeierFitter()
    kmf.fit(time[mask], event_observed=event[mask])
    try:
        surv = float(kmf.survival_function_at_times([horizon]).iloc[0])
        return 1.0 - surv
    except Exception:
        return np.nan


def make_km_table(df, split_name, horizons):
    sub = df[df["mortality_epoch_quartile"].notna()].copy()
    rows = []

    for quartile, g in sub.groupby("mortality_epoch_quartile"):
        row = {
            "split": split_name,
            "quartile": quartile,
            "n": int(len(g)),
            "deaths": int(g["event_death"].sum()),
            "event_rate_raw": float(g["event_death"].mean()) if len(g) else np.nan,
            "median_followup_years": float(g["followup_years"].median()) if len(g) else np.nan,
        }
        for h in horizons:
            row[f"observed_mortality_{h:g}yr"] = observed_km_risk_at_horizon(
                g["followup_years"], g["event_death"], h
            )
        rows.append(row)

    return pd.DataFrame(rows)


def logrank_by_quartile(df, split_name):
    sub = df[df["mortality_epoch_quartile"].notna()].copy()
    if len(sub) < 20 or sub["event_death"].nunique() < 2 or sub["mortality_epoch_quartile"].nunique() < 2:
        return {
            "split": split_name,
            "test": "multivariate_logrank_by_quartile",
            "p": np.nan,
            "test_statistic": np.nan,
            "status": "too_few_groups_or_events",
        }
    try:
        res = multivariate_logrank_test(
            event_durations=sub["followup_years"],
            groups=sub["mortality_epoch_quartile"],
            event_observed=sub["event_death"],
        )
        return {
            "split": split_name,
            "test": "multivariate_logrank_by_quartile",
            "p": float(res.p_value),
            "test_statistic": float(res.test_statistic),
            "status": "ok",
        }
    except Exception as e:
        return {
            "split": split_name,
            "test": "multivariate_logrank_by_quartile",
            "p": np.nan,
            "test_statistic": np.nan,
            "status": f"error: {str(e)[:300]}",
        }


def plot_km_by_quartile(df, out_prefix, title):
    sub = df[df["mortality_epoch_quartile"].notna()].copy()
    if sub.empty:
        return

    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    kmf = KaplanMeierFitter()

    order = ["Q1_lowest", "Q2", "Q3", "Q4_highest"]
    for q in order:
        g = sub[sub["mortality_epoch_quartile"].astype(str) == q]
        if len(g) < 5:
            continue
        kmf.fit(g["followup_years"], event_observed=g["event_death"], label=f"{q} (n={len(g)})")
        # Plot cumulative mortality = 1 - survival.
        surv = kmf.survival_function_
        ax.step(surv.index, 1.0 - surv.iloc[:, 0], where="post", label=f"{q} (n={len(g)})")

    ax.set_xlabel("Years since 2001 baseline")
    ax.set_ylabel("Cumulative mortality")
    ax.set_title(title)
    ax.legend(frameon=False)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(f"{out_prefix}.pdf")
    fig.savefig(f"{out_prefix}.png", dpi=300)
    plt.close(fig)


def plot_calibration_bar(cal_df, split_name, horizon, out_prefix):
    col = f"observed_mortality_{horizon:g}yr"
    sub = cal_df[(cal_df["split"] == split_name) & cal_df[col].notna()].copy()
    if sub.empty:
        return
    order = ["Q1_lowest", "Q2", "Q3", "Q4_highest"]
    sub["quartile"] = pd.Categorical(sub["quartile"], categories=order, ordered=True)
    sub = sub.sort_values("quartile")

    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    ax.bar(sub["quartile"].astype(str), sub[col])
    ax.set_xlabel("Mortality EPOCH acceleration quartile")
    ax.set_ylabel(f"Observed {horizon:g}-year mortality")
    ax.set_title(f"{split_name}: observed mortality by EPOCH quartile")
    ax.set_ylim(0, min(1, max(0.05, float(sub[col].max()) * 1.25)))
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(f"{out_prefix}.pdf")
    fig.savefig(f"{out_prefix}.png", dpi=300)
    plt.close(fig)


def plot_cindex_summary(cindex_df, out_prefix):
    sub = cindex_df[cindex_df["score"].isin(["lp_total", "mortality_epoch_acceleration_z", "clinical_baseline_lp"])].copy()
    sub = sub[sub["split"].isin(["train", "validation", "test"])]
    if sub.empty:
        return

    score_order = ["clinical_baseline_lp", "lp_total", "mortality_epoch_acceleration_z"]
    labels = {
        "clinical_baseline_lp": "Clinical baseline",
        "lp_total": "EPOCH LP",
        "mortality_epoch_acceleration_z": "EPOCH acceleration",
    }
    split_order = ["train", "validation", "test"]

    rows = []
    for split in split_order:
        for score in score_order:
            r = sub[(sub["split"] == split) & (sub["score"] == score)]
            if len(r):
                rr = r.iloc[0]
                rows.append({
                    "label": f"{split}\n{labels.get(score, score)}",
                    "cindex": rr["cindex"],
                    "ci_lower": rr["ci_lower"],
                    "ci_upper": rr["ci_upper"],
                })
    p = pd.DataFrame(rows)
    if p.empty:
        return

    x = np.arange(len(p))
    y = p["cindex"].astype(float).values
    yerr_low = y - p["ci_lower"].astype(float).values
    yerr_high = p["ci_upper"].astype(float).values - y
    yerr = np.vstack([yerr_low, yerr_high])

    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    ax.errorbar(x, y, yerr=yerr, fmt="o", capsize=3)
    ax.set_xticks(x)
    ax.set_xticklabels(p["label"], rotation=45, ha="right")
    ax.set_ylabel("Harrell C-index")
    ax.set_title("Mortality EPOCH discrimination with bootstrap 95% CI")
    ax.set_ylim(max(0.45, np.nanmin(p["ci_lower"]) - 0.05), min(1.0, np.nanmax(p["ci_upper"]) + 0.05))
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(f"{out_prefix}.pdf")
    fig.savefig(f"{out_prefix}.png", dpi=300)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--predictions",
        default="/Users/hao/Dropbox/MHAS/step2_mortality_epoch_model/mhas_mortality_epoch_predictions.tsv",
        help="Step 2 prediction file."
    )
    parser.add_argument(
        "--step1-input",
        default="/Users/hao/Dropbox/MHAS/step1_2001_mortality_epoch/mhas_2001_step1_model_input_primary_nondisease.tsv",
        help="Step 1 model-input file. Used only for age/sex adjusted Cox models."
    )
    parser.add_argument(
        "--out-dir",
        default="/Users/hao/Dropbox/MHAS/step3_mortality_epoch_validation",
        help="Output directory."
    )
    parser.add_argument("--n-boot", type=int, default=1000, help="Bootstrap iterations for C-index CI.")
    parser.add_argument("--seed", type=int, default=20260721)
    parser.add_argument(
        "--horizons",
        default="5,10,15,20",
        help="Comma-separated time horizons in years for observed mortality calibration."
    )
    parser.add_argument(
        "--primary-split",
        default="test",
        choices=["train", "validation", "test", "all"],
        help="Primary split for manuscript-style plots."
    )
    parser.add_argument(
        "--skip-plots",
        action="store_true",
        help="Skip PDF/PNG plots."
    )
    args = parser.parse_args()

    pred_path = Path(args.predictions)
    step1_path = Path(args.step1_input)
    out_dir = Path(args.out_dir)
    safe_mkdir(out_dir)

    if not pred_path.exists():
        raise FileNotFoundError(f"Prediction file not found: {pred_path}")

    log(f"Reading Step 2 predictions: {pred_path}")
    pred = pd.read_csv(pred_path, sep="\t", low_memory=False)

    required = [
        "participant_id", "followup_years", "event_death", "split",
        "lp_total", "mortality_epoch_acceleration_z", "mortality_epoch_quartile"
    ]
    missing = [c for c in required if c not in pred.columns]
    if missing:
        raise ValueError(f"Prediction file is missing required columns: {missing}")

    pred["followup_years"] = pd.to_numeric(pred["followup_years"], errors="coerce")
    pred["event_death"] = pd.to_numeric(pred["event_death"], errors="coerce")
    pred = pred[pred["followup_years"].notna() & pred["event_death"].notna()].copy()
    pred["event_death"] = pred["event_death"].astype(int)

    # Optional Step 1 covariates for adjusted HR.
    adjustment_covars = []
    if step1_path.exists():
        log(f"Reading Step 1 model input for covariates: {step1_path}")
        step1 = pd.read_csv(step1_path, sep="\t", low_memory=False)
        covars = get_adjustment_covariates(step1)
        keep_cols = ["participant_id"] + covars
        covars = [c for c in covars if c in step1.columns]
        if covars:
            step1_small = step1[["participant_id"] + covars].copy()
            pred = pred.merge(step1_small, on="participant_id", how="left")
            adjustment_covars = covars
    else:
        log("Step 1 input not found; adjusted Cox models will be skipped.")

    log(f"Adjustment covariates: {adjustment_covars if adjustment_covars else 'none'}")

    horizons = [float(x) for x in args.horizons.split(",") if x.strip()]

    # Score columns to evaluate.
    score_cols = ["lp_total", "mortality_epoch_acceleration_z"]
    if "clinical_baseline_lp" in pred.columns:
        score_cols.append("clinical_baseline_lp")

    # C-index with bootstrap CIs by split.
    cindex_rows = []
    for split_name in ["train", "validation", "test", "all"]:
        sub = pred.copy() if split_name == "all" else pred[pred["split"] == split_name].copy()
        for score in score_cols:
            log(f"Bootstrap C-index: split={split_name}, score={score}")
            cindex_rows.append(
                bootstrap_cindex(
                    sub,
                    score_col=score,
                    split_name=split_name,
                    n_boot=args.n_boot,
                    seed=args.seed + len(cindex_rows),
                )
            )

    cindex_df = pd.DataFrame(cindex_rows)

    # Cox HR per 1-SD score increase.
    hr_rows = []
    for split_name in ["train", "validation", "test", "all"]:
        sub = pred.copy() if split_name == "all" else pred[pred["split"] == split_name].copy()
        for score in ["mortality_epoch_acceleration_z", "lp_total"]:
            hr_rows.append(fit_cox_per_sd(sub, score, split_name, adjustment_covars=[]))
            if adjustment_covars:
                hr_rows.append(fit_cox_per_sd(sub, score, split_name, adjustment_covars=adjustment_covars))
    hr_df = pd.DataFrame(hr_rows)

    # Quartile HRs.
    qhr_tables = []
    for split_name in ["train", "validation", "test", "all"]:
        sub = pred.copy() if split_name == "all" else pred[pred["split"] == split_name].copy()
        qhr_tables.append(fit_cox_quartiles(sub, split_name, adjustment_covars=[]))
        if adjustment_covars:
            qhr_tables.append(fit_cox_quartiles(sub, split_name, adjustment_covars=adjustment_covars))
    qhr_df = pd.concat(qhr_tables, ignore_index=True)

    # KM observed risk and log-rank tests.
    km_tables = []
    logrank_rows = []
    for split_name in ["train", "validation", "test", "all"]:
        sub = pred.copy() if split_name == "all" else pred[pred["split"] == split_name].copy()
        km_tables.append(make_km_table(sub, split_name, horizons))
        logrank_rows.append(logrank_by_quartile(sub, split_name))
    km_df = pd.concat(km_tables, ignore_index=True)
    logrank_df = pd.DataFrame(logrank_rows)

    # Validation summary.
    summary_rows = []
    for split_name in ["train", "validation", "test", "all"]:
        sub = pred.copy() if split_name == "all" else pred[pred["split"] == split_name].copy()
        row = {
            "split": split_name,
            "n": int(len(sub)),
            "deaths": int(sub["event_death"].sum()),
            "event_rate": float(sub["event_death"].mean()) if len(sub) else np.nan,
            "median_followup_years": float(sub["followup_years"].median()) if len(sub) else np.nan,
        }
        for score in score_cols:
            r = cindex_df[(cindex_df["split"] == split_name) & (cindex_df["score"] == score)]
            if len(r):
                rr = r.iloc[0]
                prefix = clean_colname(score)
                row[f"{prefix}_cindex"] = rr["cindex"]
                row[f"{prefix}_cindex_ci_lower"] = rr["ci_lower"]
                row[f"{prefix}_cindex_ci_upper"] = rr["ci_upper"]

        h = hr_df[
            (hr_df["split"] == split_name)
            & (hr_df["score"] == "mortality_epoch_acceleration_z")
            & (hr_df["adjusted"] == False)
        ]
        if len(h):
            hh = h.iloc[0]
            row["epoch_acceleration_hr_per_sd"] = hh["hr"]
            row["epoch_acceleration_hr_ci_lower"] = hh["ci_lower"]
            row["epoch_acceleration_hr_ci_upper"] = hh["ci_upper"]
            row["epoch_acceleration_hr_p"] = hh["p"]

        lr = logrank_df[logrank_df["split"] == split_name]
        if len(lr):
            row["quartile_logrank_p"] = lr.iloc[0]["p"]
        summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows)

    # Save tables.
    cindex_out = out_dir / "mhas_step3_cindex_bootstrap_ci.tsv"
    hr_out = out_dir / "mhas_step3_cox_hr_per_sd.tsv"
    qhr_out = out_dir / "mhas_step3_cox_hr_by_quartile.tsv"
    km_out = out_dir / "mhas_step3_km_mortality_by_quartile.tsv"
    logrank_out = out_dir / "mhas_step3_logrank_by_quartile.tsv"
    summary_out = out_dir / "mhas_step3_validation_summary.tsv"
    audit_out = out_dir / "mhas_step3_audit.txt"

    cindex_df.to_csv(cindex_out, sep="\t", index=False)
    hr_df.to_csv(hr_out, sep="\t", index=False)
    qhr_df.to_csv(qhr_out, sep="\t", index=False)
    km_df.to_csv(km_out, sep="\t", index=False)
    logrank_df.to_csv(logrank_out, sep="\t", index=False)
    summary_df.to_csv(summary_out, sep="\t", index=False)

    # Plots.
    if not args.skip_plots:
        primary = pred.copy() if args.primary_split == "all" else pred[pred["split"] == args.primary_split].copy()

        plot_km_by_quartile(
            primary,
            out_prefix=str(out_dir / f"mhas_step3_fig_km_{args.primary_split}"),
            title=f"{args.primary_split}: cumulative mortality by EPOCH acceleration quartile",
        )

        # Prefer 10-year calibration if available; otherwise use first horizon.
        cal_horizon = 10.0 if 10.0 in horizons else horizons[0]
        plot_calibration_bar(
            km_df,
            split_name=args.primary_split,
            horizon=cal_horizon,
            out_prefix=str(out_dir / f"mhas_step3_fig_calibration_{args.primary_split}_{cal_horizon:g}yr"),
        )

        plot_cindex_summary(
            cindex_df,
            out_prefix=str(out_dir / "mhas_step3_fig_cindex_summary"),
        )

    audit = f"""MHAS STEP 3: mortality EPOCH validation

Inputs
------
Predictions:
{pred_path}

Step 1 model input:
{step1_path if step1_path.exists() else "not found"}

Output directory
----------------
{out_dir}

Data summary
------------
N total: {len(pred):,}
Deaths total: {int(pred["event_death"].sum()):,}
Median follow-up years: {float(pred["followup_years"].median()):.2f}

Scores evaluated
----------------
{", ".join(score_cols)}

Adjustment covariates for adjusted Cox models
---------------------------------------------
{", ".join(adjustment_covars) if adjustment_covars else "None"}

Bootstrap settings
------------------
n_boot: {args.n_boot}
seed: {args.seed}

Time horizons for observed mortality calibration
------------------------------------------------
{", ".join([str(h) for h in horizons])} years

Primary split for figures
-------------------------
{args.primary_split}

Validation summary
------------------
{summary_df.to_string(index=False)}

Output files
------------
C-index bootstrap CI:
{cindex_out}

Cox HR per 1-SD:
{hr_out}

Cox HR by quartile:
{qhr_out}

Observed KM mortality by quartile/horizon:
{km_out}

Log-rank tests:
{logrank_out}

Validation summary:
{summary_out}

Figures:
{out_dir / f"mhas_step3_fig_km_{args.primary_split}.pdf"}
{out_dir / f"mhas_step3_fig_calibration_{args.primary_split}_{(10.0 if 10.0 in horizons else horizons[0]):g}yr.pdf"}
{out_dir / "mhas_step3_fig_cindex_summary.pdf"}

Interpretation note
-------------------
Use the test split as the primary internal hold-out estimate. The all-sample results
are useful for descriptive summaries but are not independent validation.
"""
    audit_out.write_text(audit)
    log("\n" + audit)
    log("STEP 3 finished successfully.")


if __name__ == "__main__":
    main()
