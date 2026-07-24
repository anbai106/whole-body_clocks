#!/usr/bin/env python3
"""
STEP 4 UPDATED: Baseline prevalent disease logistic regression for MHAS mortality EPOCH

This replaces the earlier incident/future-disease Step 4 analysis.

Why this update?
----------------
The earlier Step 4 used *_incident_after_w1 outcomes. Because detailed dates of
diagnosis are not available and missing later disease assessments can be strongly
related to mortality/survival, the prospective "future disease" logistic analysis
can be biased.

This updated script matches the NHANES-style disease-linking analysis:
  baseline prevalent disease status ~ mortality EPOCH score + covariates

Model:
  baseline_disease_w1 ~ mortality_epoch_acceleration_z + age_2001 + sex

Primary outcome columns:
  *_baseline_w1

Inputs:
  Step 2 predictions:
    /Users/hao/Dropbox/MHAS/step2_mortality_epoch_model/mhas_mortality_epoch_predictions.tsv

  Step 1 disease labels:
    /Users/hao/Dropbox/MHAS/step1_2001_mortality_epoch/mhas_2001_step1_downstream_disease_labels.tsv

  Step 1 model input, optional but recommended for age/sex adjustment:
    /Users/hao/Dropbox/MHAS/step1_2001_mortality_epoch/mhas_2001_step1_model_input_primary_nondisease.tsv

Outputs, compatible with existing Step 5/6 filenames:
  /Users/hao/Dropbox/MHAS/step4_epoch_disease_associations/
    mhas_step4_disease_or_per_sd.tsv
    mhas_step4_disease_or_by_quartile.tsv
    mhas_step4_disease_event_rate_by_quartile.tsv
    mhas_step4_disease_association_summary.tsv
    mhas_step4_fig_disease_forest_all.pdf/png
    mhas_step4_fig_top_disease_quartile_rate_all.pdf/png
    mhas_step4_audit.txt

Run:
  cd /Users/hao/Project/whole-body_clocks/MHAS
  /Users/hao/opt/anaconda3/envs/DNE/bin/python 4_step4_mhas_epoch_baseline_prevalent_disease.py

Optional:
  --primary-split all      # default, recommended for baseline prevalent disease
  --primary-split test     # if you want to keep strict hold-out partitioning
"""

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    import statsmodels.api as sm
    HAVE_STATSMODELS = True
except Exception:
    HAVE_STATSMODELS = False

try:
    from sklearn.metrics import roc_auc_score
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    HAVE_SKLEARN = True
except Exception:
    HAVE_SKLEARN = False


def log(msg: str) -> None:
    print(msg, flush=True)


def safe_mkdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def bh_fdr(pvals):
    p = pd.Series(pvals, dtype="float64")
    q = pd.Series(np.nan, index=p.index, dtype="float64")
    valid = p.notna()
    if valid.sum() == 0:
        return q
    pv = p[valid].values
    order = np.argsort(pv)
    ranked = pv[order]
    m = len(ranked)
    qvals = ranked * m / np.arange(1, m + 1)
    qvals = np.minimum.accumulate(qvals[::-1])[::-1]
    qvals = np.minimum(qvals, 1.0)
    out = np.empty(m)
    out[order] = qvals
    q.loc[valid] = out
    return q


def parse_endpoint_name(col):
    suffix = "_baseline_w1"
    if col.endswith(suffix):
        return col[: -len(suffix)]
    return col


def get_adjustment_covariates(step1_df, mode="age_sex"):
    covars = []
    if "age_2001" in step1_df.columns:
        covars.append("age_2001")
    covars += [c for c in step1_df.columns if c.startswith("sex_")]

    if mode == "expanded":
        prefixes = [
            "current_smoking_2001_",
            "ever_smoked_2001_",
            "bmi_self_report_2001",
            "bmi_measured_2001",
            "self_rated_health_2001_",
        ]
        for c in step1_df.columns:
            if any(c == p or c.startswith(p) for p in prefixes):
                covars.append(c)

    seen = set()
    return [c for c in covars if not (c in seen or seen.add(c))]


def prepare_xy(df, outcome_col, score_col, covars=None):
    covars = covars or []
    cols = ["participant_id", outcome_col, score_col] + covars
    cols = [c for c in cols if c in df.columns]

    sub = df[cols].copy()
    sub[outcome_col] = pd.to_numeric(sub[outcome_col], errors="coerce")
    sub[score_col] = pd.to_numeric(sub[score_col], errors="coerce")

    for c in covars:
        if c in sub.columns:
            sub[c] = pd.to_numeric(sub[c], errors="coerce")

    sub = sub.replace([np.inf, -np.inf], np.nan).dropna(subset=[outcome_col, score_col])
    sub[outcome_col] = sub[outcome_col].astype(int)

    for c in covars:
        if c in sub.columns:
            med = sub[c].median()
            if pd.isna(med):
                sub = sub.drop(columns=[c])
            else:
                sub[c] = sub[c].fillna(med)

    keep_covars = []
    for c in covars:
        if c in sub.columns and sub[c].nunique(dropna=True) > 1:
            keep_covars.append(c)

    return sub, keep_covars


def fit_logistic_or(df, outcome_col, score_col, split_name, endpoint, covars=None,
                    min_cases=20, min_controls=20):
    covars = covars or []
    sub, keep_covars = prepare_xy(df, outcome_col, score_col, covars=covars)

    n = len(sub)
    cases = int(sub[outcome_col].sum()) if n else 0
    controls = n - cases

    base = {
        "analysis_type": "baseline_prevalent",
        "endpoint": endpoint,
        "outcome_col": outcome_col,
        "score": score_col,
        "split": split_name,
        "adjusted": bool(keep_covars),
        "adjustment_covariates": ",".join(keep_covars),
        "exclude_prevalent": False,
        "n": n,
        "cases": cases,
        "controls": controls,
        "event_rate": cases / n if n else np.nan,
        "auc_score_only": np.nan,
        "or_per_sd": np.nan,
        "ci_lower": np.nan,
        "ci_upper": np.nan,
        "p": np.nan,
        "status": "not_run",
    }

    if cases < min_cases or controls < min_controls:
        base["status"] = "too_few_cases_or_controls"
        return base

    try:
        base["auc_score_only"] = float(roc_auc_score(sub[outcome_col], sub[score_col]))
    except Exception:
        pass

    score_sd = sub[score_col].std(ddof=0)
    if score_sd == 0 or pd.isna(score_sd):
        base["status"] = "zero_score_sd"
        return base

    sub["score_per_sd"] = (sub[score_col] - sub[score_col].mean()) / score_sd
    model_cols = ["score_per_sd"] + keep_covars
    X = sub[model_cols].copy()
    y = sub[outcome_col].astype(int)

    if HAVE_STATSMODELS:
        try:
            X_sm = sm.add_constant(X, has_constant="add")
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                fit = sm.Logit(y, X_sm).fit(disp=False, maxiter=300)

            coef = fit.params["score_per_sd"]
            conf = fit.conf_int().loc["score_per_sd"]
            pval = fit.pvalues["score_per_sd"]

            base.update({
                "or_per_sd": float(np.exp(coef)),
                "ci_lower": float(np.exp(conf[0])),
                "ci_upper": float(np.exp(conf[1])),
                "p": float(pval),
                "status": "ok",
            })
            return base
        except Exception as e:
            base["status"] = f"statsmodels_error: {str(e)[:250]}"

    if HAVE_SKLEARN:
        try:
            scaler = StandardScaler()
            X_scaled = scaler.fit_transform(X)
            clf = LogisticRegression(max_iter=1000, penalty="l2", solver="lbfgs")
            clf.fit(X_scaled, y)
            score_coef = clf.coef_[0][0]
            base.update({
                "or_per_sd": float(np.exp(score_coef)),
                "status": base["status"] + "; sklearn_fallback_no_ci",
            })
        except Exception as e:
            base["status"] = base["status"] + f"; sklearn_error: {str(e)[:250]}"

    return base


def fit_quartile_or(df, outcome_col, endpoint, split_name, covars=None,
                    min_cases=20, min_controls=20):
    covars = covars or []
    qcol = "mortality_epoch_quartile"
    needed = ["participant_id", outcome_col, qcol] + covars
    needed = [c for c in needed if c in df.columns]

    sub = df[needed].copy()
    sub[outcome_col] = pd.to_numeric(sub[outcome_col], errors="coerce")
    sub = sub.dropna(subset=[outcome_col, qcol])
    sub[outcome_col] = sub[outcome_col].astype(int)

    n = len(sub)
    cases = int(sub[outcome_col].sum()) if n else 0
    controls = n - cases

    if cases < min_cases or controls < min_controls:
        return pd.DataFrame([{
            "analysis_type": "baseline_prevalent",
            "endpoint": endpoint,
            "outcome_col": outcome_col,
            "split": split_name,
            "comparison": "quartiles",
            "adjusted": bool(covars),
            "adjustment_covariates": ",".join(covars),
            "n": n,
            "cases": cases,
            "controls": controls,
            "or": np.nan,
            "ci_lower": np.nan,
            "ci_upper": np.nan,
            "p": np.nan,
            "status": "too_few_cases_or_controls",
        }])

    dummies = pd.get_dummies(sub[qcol].astype(str), prefix="quartile", dtype=int)
    for col in ["quartile_Q2", "quartile_Q3", "quartile_Q4_highest"]:
        if col not in dummies.columns:
            dummies[col] = 0
    dummies = dummies[["quartile_Q2", "quartile_Q3", "quartile_Q4_highest"]]

    X = dummies.copy()
    keep_covars = []
    for c in covars:
        if c in sub.columns:
            x = pd.to_numeric(sub[c], errors="coerce")
            med = x.median()
            if not pd.isna(med):
                x = x.fillna(med)
                if x.nunique(dropna=True) > 1:
                    X[c] = x.values
                    keep_covars.append(c)

    if HAVE_STATSMODELS:
        try:
            X_sm = sm.add_constant(X, has_constant="add")
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                fit = sm.Logit(sub[outcome_col], X_sm).fit(disp=False, maxiter=300)

            rows = []
            for col, label in [
                ("quartile_Q2", "Q2 vs Q1_lowest"),
                ("quartile_Q3", "Q3 vs Q1_lowest"),
                ("quartile_Q4_highest", "Q4_highest vs Q1_lowest"),
            ]:
                coef = fit.params[col]
                conf = fit.conf_int().loc[col]
                rows.append({
                    "analysis_type": "baseline_prevalent",
                    "endpoint": endpoint,
                    "outcome_col": outcome_col,
                    "split": split_name,
                    "comparison": label,
                    "adjusted": bool(keep_covars),
                    "adjustment_covariates": ",".join(keep_covars),
                    "n": n,
                    "cases": cases,
                    "controls": controls,
                    "or": float(np.exp(coef)),
                    "ci_lower": float(np.exp(conf[0])),
                    "ci_upper": float(np.exp(conf[1])),
                    "p": float(fit.pvalues[col]),
                    "status": "ok",
                })
            return pd.DataFrame(rows)
        except Exception as e:
            return pd.DataFrame([{
                "analysis_type": "baseline_prevalent",
                "endpoint": endpoint,
                "outcome_col": outcome_col,
                "split": split_name,
                "comparison": "quartiles",
                "adjusted": bool(keep_covars),
                "adjustment_covariates": ",".join(keep_covars),
                "n": n,
                "cases": cases,
                "controls": controls,
                "or": np.nan,
                "ci_lower": np.nan,
                "ci_upper": np.nan,
                "p": np.nan,
                "status": f"error: {str(e)[:250]}",
            }])

    return pd.DataFrame([{
        "analysis_type": "baseline_prevalent",
        "endpoint": endpoint,
        "outcome_col": outcome_col,
        "split": split_name,
        "comparison": "quartiles",
        "adjusted": bool(keep_covars),
        "adjustment_covariates": ",".join(keep_covars),
        "n": n,
        "cases": cases,
        "controls": controls,
        "or": np.nan,
        "ci_lower": np.nan,
        "ci_upper": np.nan,
        "p": np.nan,
        "status": "statsmodels_not_available",
    }])


def event_rate_by_quartile(df, outcome_col, endpoint, split_name):
    qcol = "mortality_epoch_quartile"
    sub = df[["participant_id", outcome_col, qcol]].copy()
    sub[outcome_col] = pd.to_numeric(sub[outcome_col], errors="coerce")
    sub = sub.dropna(subset=[outcome_col, qcol])
    sub[outcome_col] = sub[outcome_col].astype(int)

    rows = []
    order = ["Q1_lowest", "Q2", "Q3", "Q4_highest"]
    for q in order:
        g = sub[sub[qcol].astype(str) == q]
        rows.append({
            "analysis_type": "baseline_prevalent",
            "endpoint": endpoint,
            "outcome_col": outcome_col,
            "split": split_name,
            "quartile": q,
            "n": int(len(g)),
            "cases": int(g[outcome_col].sum()) if len(g) else 0,
            "event_rate": float(g[outcome_col].mean()) if len(g) else np.nan,
        })
    return pd.DataFrame(rows)


def plot_disease_forest(results, out_prefix, split_name="all"):
    sub = results[
        (results["split"] == split_name)
        & (results["score"] == "mortality_epoch_acceleration_z")
        & (results["status"].astype(str).str.startswith("ok"))
    ].copy()

    if "adjusted" in sub.columns and sub["adjusted"].any():
        sub = sub[sub["adjusted"] == True].copy()
    else:
        sub = sub[sub["adjusted"] == False].copy()

    sub = sub.dropna(subset=["or_per_sd", "ci_lower", "ci_upper"])
    if sub.empty:
        return

    sub = sub.sort_values("or_per_sd", ascending=True)
    fig_h = max(4.5, 0.45 * len(sub) + 1.5)
    fig, ax = plt.subplots(figsize=(7.8, fig_h))
    y = np.arange(len(sub))
    x = sub["or_per_sd"].astype(float).values
    lo = sub["ci_lower"].astype(float).values
    hi = sub["ci_upper"].astype(float).values

    ax.errorbar(x, y, xerr=[x - lo, hi - x], fmt="o", capsize=3)
    ax.axvline(1.0, linestyle="--", linewidth=1)
    ax.set_yticks(y)
    ax.set_yticklabels(sub["endpoint"])
    ax.set_xlabel("Odds ratio per 1-SD higher mortality EPOCH acceleration")
    ax.set_title(f"{split_name}: baseline prevalent disease associations")
    ax.grid(True, axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(f"{out_prefix}.pdf")
    fig.savefig(f"{out_prefix}.png", dpi=300)
    plt.close(fig)


def plot_top_disease_quartile_rates(rate_df, association_df, out_prefix, split_name="all"):
    sub_assoc = association_df[
        (association_df["split"] == split_name)
        & (association_df["score"] == "mortality_epoch_acceleration_z")
        & association_df["p"].notna()
    ].copy()
    if sub_assoc.empty:
        return
    if sub_assoc["adjusted"].any():
        sub_assoc = sub_assoc[sub_assoc["adjusted"] == True].copy()

    sub_assoc = sub_assoc.sort_values("p")
    endpoint = sub_assoc.iloc[0]["endpoint"]
    sub = rate_df[(rate_df["split"] == split_name) & (rate_df["endpoint"] == endpoint)].copy()
    if sub.empty:
        return

    order = ["Q1_lowest", "Q2", "Q3", "Q4_highest"]
    sub["quartile"] = pd.Categorical(sub["quartile"], categories=order, ordered=True)
    sub = sub.sort_values("quartile")

    fig, ax = plt.subplots(figsize=(6.5, 4.6))
    ax.bar(sub["quartile"].astype(str), sub["event_rate"])
    ax.set_xlabel("Mortality EPOCH acceleration quartile")
    ax.set_ylabel("Baseline disease prevalence")
    ax.set_title(f"{split_name}: baseline {endpoint} prevalence by EPOCH quartile")
    ax.set_ylim(0, min(1, max(0.05, float(sub["event_rate"].max()) * 1.25)))
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
        "--disease-labels",
        default="/Users/hao/Dropbox/MHAS/step1_2001_mortality_epoch/mhas_2001_step1_downstream_disease_labels.tsv",
        help="Step 1 disease-label file."
    )
    parser.add_argument(
        "--step1-input",
        default="/Users/hao/Dropbox/MHAS/step1_2001_mortality_epoch/mhas_2001_step1_model_input_primary_nondisease.tsv",
        help="Step 1 model input for covariates."
    )
    parser.add_argument(
        "--out-dir",
        default="/Users/hao/Dropbox/MHAS/step4_epoch_disease_associations",
        help="Output directory. Default intentionally matches Step 5/6 expectations."
    )
    parser.add_argument(
        "--score-cols",
        default="mortality_epoch_acceleration_z,lp_total",
        help="Comma-separated score columns to test."
    )
    parser.add_argument(
        "--covariate-mode",
        default="age_sex",
        choices=["none", "age_sex", "expanded"],
        help="Covariate adjustment mode."
    )
    parser.add_argument("--min-cases", type=int, default=20)
    parser.add_argument("--min-controls", type=int, default=20)
    parser.add_argument(
        "--primary-split",
        default="all",
        choices=["train", "validation", "test", "all"],
        help="Primary split for summary/plots. Default is all for baseline prevalent disease."
    )
    parser.add_argument("--skip-plots", action="store_true")
    args = parser.parse_args()

    pred_path = Path(args.predictions)
    disease_path = Path(args.disease_labels)
    step1_path = Path(args.step1_input)
    out_dir = Path(args.out_dir)
    safe_mkdir(out_dir)

    if not pred_path.exists():
        raise FileNotFoundError(f"Prediction file not found: {pred_path}")
    if not disease_path.exists():
        raise FileNotFoundError(f"Disease-label file not found: {disease_path}")

    log(f"Reading predictions: {pred_path}")
    pred = pd.read_csv(pred_path, sep="\t", low_memory=False)
    log(f"Reading disease labels: {disease_path}")
    disease = pd.read_csv(disease_path, sep="\t", low_memory=False)

    df = pred.merge(disease, on="participant_id", how="inner")
    log(f"Merged N: {len(df):,}")

    covars = []
    if args.covariate_mode != "none" and step1_path.exists():
        step1 = pd.read_csv(step1_path, sep="\t", low_memory=False)
        covars = get_adjustment_covariates(step1, mode=args.covariate_mode)
        covars = [c for c in covars if c in step1.columns]
        if covars:
            df = df.merge(step1[["participant_id"] + covars], on="participant_id", how="left")
    log(f"Covariates: {covars if covars else 'none'}")

    score_cols = [x.strip() for x in args.score_cols.split(",") if x.strip()]
    score_cols = [s for s in score_cols if s in df.columns]
    if not score_cols:
        raise ValueError("No requested score columns are present in predictions file.")

    outcome_cols = [c for c in df.columns if c.endswith("_baseline_w1")]
    if not outcome_cols:
        raise ValueError("No *_baseline_w1 columns found in disease-label file.")

    valid_outcomes = []
    for c in outcome_cols:
        x = pd.to_numeric(df[c], errors="coerce")
        if x.notna().sum() > 0 and x.dropna().nunique() >= 2:
            valid_outcomes.append(c)
        else:
            log(f"Skipping {c}: no variation or all missing.")
    outcome_cols = valid_outcomes
    if not outcome_cols:
        raise ValueError("No valid baseline disease outcomes with variation.")

    splits = ["train", "validation", "test", "all"]
    per_sd_rows = []
    quartile_tables = []
    rate_tables = []

    for outcome_col in outcome_cols:
        endpoint = parse_endpoint_name(outcome_col)

        for split_name in splits:
            sub = df.copy() if split_name == "all" else df[df["split"] == split_name].copy()

            for score in score_cols:
                per_sd_rows.append(
                    fit_logistic_or(
                        sub,
                        outcome_col=outcome_col,
                        score_col=score,
                        split_name=split_name,
                        endpoint=endpoint,
                        covars=[],
                        min_cases=args.min_cases,
                        min_controls=args.min_controls,
                    )
                )
                if covars:
                    per_sd_rows.append(
                        fit_logistic_or(
                            sub,
                            outcome_col=outcome_col,
                            score_col=score,
                            split_name=split_name,
                            endpoint=endpoint,
                            covars=covars,
                            min_cases=args.min_cases,
                            min_controls=args.min_controls,
                        )
                    )

            quartile_tables.append(
                fit_quartile_or(
                    sub,
                    outcome_col=outcome_col,
                    endpoint=endpoint,
                    split_name=split_name,
                    covars=[],
                    min_cases=args.min_cases,
                    min_controls=args.min_controls,
                )
            )
            if covars:
                quartile_tables.append(
                    fit_quartile_or(
                        sub,
                        outcome_col=outcome_col,
                        endpoint=endpoint,
                        split_name=split_name,
                        covars=covars,
                        min_cases=args.min_cases,
                        min_controls=args.min_controls,
                    )
                )

            rate_tables.append(
                event_rate_by_quartile(
                    sub,
                    outcome_col=outcome_col,
                    endpoint=endpoint,
                    split_name=split_name,
                )
            )

    per_sd = pd.DataFrame(per_sd_rows)
    quartile = pd.concat(quartile_tables, ignore_index=True) if quartile_tables else pd.DataFrame()
    rates = pd.concat(rate_tables, ignore_index=True) if rate_tables else pd.DataFrame()

    per_sd["fdr_bh"] = np.nan
    for _, idx in per_sd.groupby(["split", "score", "adjusted"]).groups.items():
        per_sd.loc[idx, "fdr_bh"] = bh_fdr(per_sd.loc[idx, "p"]).values

    if not quartile.empty:
        quartile["fdr_bh"] = np.nan
        for _, idx in quartile.groupby(["split", "comparison", "adjusted"]).groups.items():
            quartile.loc[idx, "fdr_bh"] = bh_fdr(quartile.loc[idx, "p"]).values

    primary = per_sd[
        (per_sd["split"] == args.primary_split)
        & (per_sd["score"] == "mortality_epoch_acceleration_z")
    ].copy()
    if primary["adjusted"].any():
        primary = primary[primary["adjusted"] == True].copy()
    else:
        primary = primary[primary["adjusted"] == False].copy()
    summary = primary.sort_values(["p", "endpoint"], na_position="last").copy()

    per_sd_out = out_dir / "mhas_step4_disease_or_per_sd.tsv"
    quartile_out = out_dir / "mhas_step4_disease_or_by_quartile.tsv"
    rates_out = out_dir / "mhas_step4_disease_event_rate_by_quartile.tsv"
    summary_out = out_dir / "mhas_step4_disease_association_summary.tsv"
    audit_out = out_dir / "mhas_step4_audit.txt"

    per_sd.to_csv(per_sd_out, sep="\t", index=False)
    quartile.to_csv(quartile_out, sep="\t", index=False)
    rates.to_csv(rates_out, sep="\t", index=False)
    summary.to_csv(summary_out, sep="\t", index=False)

    if not args.skip_plots:
        plot_disease_forest(
            per_sd,
            out_prefix=str(out_dir / f"mhas_step4_fig_disease_forest_{args.primary_split}"),
            split_name=args.primary_split,
        )
        plot_top_disease_quartile_rates(
            rates,
            per_sd,
            out_prefix=str(out_dir / f"mhas_step4_fig_top_disease_quartile_rate_{args.primary_split}"),
            split_name=args.primary_split,
        )

    top_lines = summary.head(20).to_string(index=False) if len(summary) else "No valid endpoint results."

    audit = f"""MHAS STEP 4 UPDATED: baseline prevalent disease logistic regression

Inputs
------
Predictions:
{pred_path}

Disease labels:
{disease_path}

Step 1 covariates:
{step1_path if step1_path.exists() else "not found"}

Output directory
----------------
{out_dir}

Merged sample
-------------
N merged: {len(df):,}

Design
------
Analysis type: baseline prevalent disease logistic regression
Outcome columns: {len(outcome_cols)}
Outcomes: {", ".join([parse_endpoint_name(c) for c in outcome_cols])}

Model
-----
baseline disease at Wave 1 ~ mortality EPOCH score + covariates

Score columns tested:
{", ".join(score_cols)}

Covariate mode:
{args.covariate_mode}

Covariates:
{", ".join(covars) if covars else "None"}

Minimum cases/controls:
{args.min_cases}/{args.min_controls}

Primary split:
{args.primary_split}

Top primary-split associations
------------------------------
{top_lines}

Output files
------------
Per-SD ORs:
{per_sd_out}

Quartile ORs:
{quartile_out}

Disease prevalence by EPOCH quartile:
{rates_out}

Primary summary:
{summary_out}

Figures:
{out_dir / f"mhas_step4_fig_disease_forest_{args.primary_split}.pdf"}
{out_dir / f"mhas_step4_fig_top_disease_quartile_rate_{args.primary_split}.pdf"}

Notes
-----
This updated Step 4 intentionally uses baseline prevalent disease status rather
than incident/future disease labels. This matches the NHANES-style analysis where
mortality EPOCH was linked to baseline disease status. Because exact disease
diagnosis dates are unavailable in harmonized MHAS, baseline prevalent disease
logistic regression is safer than treating missing future disease assessments as
disease-free controls.
"""
    audit_out.write_text(audit)
    log("\n" + audit)
    log("STEP 4 baseline prevalent disease analysis finished successfully.")


if __name__ == "__main__":
    main()
