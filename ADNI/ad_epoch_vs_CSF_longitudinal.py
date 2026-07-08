#!/usr/bin/env python3
# ============================================================
# ADNI longitudinal AD L'EPOCH change versus baseline CSF biomarkers
#
# Purpose
# -------
# Test whether baseline CSF biomarkers are associated with subsequent
# longitudinal change in AD L'EPOCH clocks among CN-labeled follow-up
# MRI scans before MCI/AD conversion.
#
# Input 1: longitudinal AD L'EPOCH prediction file
# ------------------------------------------------
# Expected variables:
#   PTID
#   years_since_selected_baseline
#   conversion_group or conversion_group_3level
#   scan_dx
#   scan_relation_to_event
#
# Default longitudinal AD clock variables:
#   adni_brain_mri_ad_lepoch_acceleration_z
#   adni_brain_mri_ad_lepoch_acceleration_years
#   adni_brain_mri_ad_lepoch_risk_score
#
# Input 2: baseline CSF analysis dataset
# --------------------------------------
# Expected subject/covariate variables:
#   PTID
#   conversion_group_3level
#   _age
#   _sex_male
#   _icv
#   _apoe4
#
# Default baseline CSF biomarker variables:
#   Abeta_CSF
#   Tau_CSF
#   PTau_CSF
#
# Longitudinal clock-change metrics
# ---------------------------------
# For each participant and each clock:
#   baseline_clock_value
#   last_clock_value
#   delta_last_minus_baseline
#   annualized_delta
#   slope_per_year from lm(clock ~ years_since_selected_baseline)
#   n_scans
#   followup_span_years
#
# Primary combined full-sample model
# ----------------------------------
# For each clock-change metric and each CSF biomarker:
#
#   clock_change_metric ~ baseline_CSF_biomarker
#                         + baseline_clock_value
#                         + conversion_group
#                         + Age + Sex + ICV + APOE4
#                         + followup_span_years
#
# Reference group:
#   Non-event & censored
#
# Secondary exploratory group-specific model
# ------------------------------------------
# Within each conversion group:
#
#   clock_change_metric ~ baseline_CSF_biomarker
#                         + baseline_clock_value
#                         + Age + Sex + ICV + APOE4
#                         + followup_span_years
#
# Major output files
# ------------------
#   {prefix}_longitudinal_clock_change_subject_metrics.tsv
#   {prefix}_longitudinal_change_vs_baseline_csf_analysis_dataset.tsv
#   {prefix}_longitudinal_change_vs_baseline_csf_combined_associations.tsv
#   {prefix}_longitudinal_change_vs_baseline_csf_group_specific_associations.tsv
#   {prefix}_longitudinal_change_vs_baseline_csf_all_associations.tsv
#   {prefix}_longitudinal_change_vs_baseline_csf_subject_summary.tsv
#   {prefix}_longitudinal_change_vs_baseline_csf_analysis_summary.json
#
# Notes
# -----
# - Raw P-values are saved.
# - BH-adjusted P-values across successful tests are saved.
# - For Abeta_CSF, lower values are more pathological, so an additional
#   pathology-direction effect is reported:
#       pathology_direction_beta = -beta for Abeta_CSF
#       pathology_direction_beta =  beta for Tau_CSF and PTau_CSF
#   Therefore, positive pathology-direction beta means worse baseline CSF
#   pathology is associated with greater subsequent AD L'EPOCH increase.
# ============================================================

import argparse
import json
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


# -----------------------------
# 1. Arguments
# -----------------------------

def parse_args():
    p = argparse.ArgumentParser()

    p.add_argument("--longitudinal-file", required=True)
    p.add_argument("--baseline-csf-file", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--prefix", default="adni_brain_mri_ad_lepoch")

    p.add_argument("--id-col", default="PTID")
    p.add_argument("--time-col", default="years_since_selected_baseline")
    p.add_argument("--group-col", default="conversion_group_3level")

    p.add_argument(
        "--clock-cols",
        default=(
            "adni_brain_mri_ad_lepoch_acceleration_z,"
            "adni_brain_mri_ad_lepoch_acceleration_years,"
            "adni_brain_mri_ad_lepoch_risk_score"
        ),
        help="Comma-separated longitudinal AD clock columns."
    )

    p.add_argument(
        "--csf-vars",
        default="Abeta_CSF,Tau_CSF,PTau_CSF",
        help="Comma-separated baseline CSF variables."
    )

    p.add_argument("--age-col", default="_age")
    p.add_argument("--sex-col", default="_sex_male")
    p.add_argument("--icv-col", default="_icv")
    p.add_argument("--apoe4-col", default="_apoe4")

    p.add_argument("--min-scans", type=int, default=2)
    p.add_argument("--min-followup-years", type=float, default=0.25)
    p.add_argument("--min-n", type=int, default=8)

    return p.parse_args()


# -----------------------------
# 2. Helpers
# -----------------------------

def log(msg):
    print(msg, flush=True)


def read_table(path):
    path = Path(path)
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path, low_memory=False)
    return pd.read_csv(path, sep="\t", low_memory=False)


def parse_list_arg(x):
    if x is None:
        return []
    x = str(x).strip()
    if x == "" or x.lower() == "none":
        return []
    return [v.strip() for v in x.split(",") if v.strip()]


def sanitize_name(x):
    return re.sub(r"[^a-z0-9]+", "_", str(x).lower()).strip("_")


def clean_numeric_series(s):
    if pd.api.types.is_numeric_dtype(s):
        return pd.to_numeric(s, errors="coerce")

    s2 = (
        s.astype("object")
         .where(s.notna(), np.nan)
         .astype(str)
         .str.strip()
         .str.replace(",", "", regex=False)
         .str.replace(">", "", regex=False)
         .str.replace("<", "", regex=False)
    )

    s2 = s2.replace({
        "": np.nan,
        "nan": np.nan,
        "NaN": np.nan,
        "NA": np.nan,
        "N/A": np.nan,
        "None": np.nan,
        "null": np.nan,
        ".": np.nan
    })

    return pd.to_numeric(s2, errors="coerce")


def normalize_id_series(s):
    return s.astype(str).str.strip()


def normalize_group(x):
    if pd.isna(x):
        return np.nan

    x0 = str(x).strip()

    xl = x0.lower()

    if "mci" in xl:
        return "CN-MCI"

    if "ad" in xl and "censored" not in xl:
        return "CN-AD"

    if "censored" in xl or "non" in xl:
        return "Non-event & censored"

    return x0


def adjust_p_bh(pvals):
    pvals = np.asarray(pvals, dtype=float)
    out = np.full_like(pvals, np.nan, dtype=float)

    mask = np.isfinite(pvals)
    if mask.sum() == 0:
        return out

    pv = pvals[mask]
    order = np.argsort(pv)
    ranked = np.empty_like(pv, dtype=float)
    m = len(pv)

    prev = 1.0
    for i in range(m - 1, -1, -1):
        rank = i + 1
        val = min(prev, pv[order[i]] * m / rank)
        ranked[order[i]] = val
        prev = val

    out[mask] = ranked
    return out


def pathology_direction_multiplier(csf_var):
    v = str(csf_var).lower()

    # Lower CSF Abeta is more pathological.
    if "abeta" in v or "amyloid" in v:
        return -1.0

    # Higher tau/p-tau is more pathological.
    return 1.0


# -----------------------------
# 3. Longitudinal clock-change metrics
# -----------------------------

def summarize_one_subject_clock(g, id_col, time_col, clock_col, min_scans, min_followup_years):
    pid = g[id_col].iloc[0]

    d = g[[id_col, time_col, clock_col]].copy()
    d[time_col] = clean_numeric_series(d[time_col])
    d[clock_col] = clean_numeric_series(d[clock_col])
    d = d.replace([np.inf, -np.inf], np.nan).dropna()
    d = d.sort_values(time_col, kind="mergesort")

    out = {
        id_col: pid,
        "clock_col": clock_col,
        "n_scans": int(d.shape[0]),
        "min_year": np.nan,
        "max_year": np.nan,
        "followup_span_years": np.nan,
        "baseline_clock_value": np.nan,
        "last_clock_value": np.nan,
        "mean_clock_value": np.nan,
        "sd_clock_value": np.nan,
        "delta_last_minus_baseline": np.nan,
        "annualized_delta": np.nan,
        "slope_per_year": np.nan,
        "slope_intercept": np.nan,
        "metric_status": "not_enough_scans"
    }

    if d.shape[0] < min_scans:
        return out

    t = d[time_col].values.astype(float)
    y = d[clock_col].values.astype(float)

    span = float(np.nanmax(t) - np.nanmin(t))

    out["min_year"] = float(np.nanmin(t))
    out["max_year"] = float(np.nanmax(t))
    out["followup_span_years"] = span
    out["mean_clock_value"] = float(np.nanmean(y))
    out["sd_clock_value"] = float(np.nanstd(y, ddof=1)) if len(y) > 1 else np.nan

    if span < min_followup_years:
        out["metric_status"] = "followup_span_too_short"
        return out

    baseline_idx = int(np.argmin(np.abs(t)))
    last_idx = int(np.argmax(t))

    baseline_value = float(y[baseline_idx])
    last_value = float(y[last_idx])
    delta = last_value - baseline_value

    out["baseline_clock_value"] = baseline_value
    out["last_clock_value"] = last_value
    out["delta_last_minus_baseline"] = delta
    out["annualized_delta"] = delta / span if span > 0 else np.nan

    if len(np.unique(t)) < 2:
        out["metric_status"] = "time_has_no_variation"
        return out

    X = np.column_stack([np.ones(len(t)), t])

    try:
        beta, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
        out["slope_intercept"] = float(beta[0])
        out["slope_per_year"] = float(beta[1])
        out["metric_status"] = "ok"
    except Exception:
        out["metric_status"] = "slope_fit_failed"

    return out


def build_subject_clock_metrics(long_df, args, clock_cols):
    rows = []

    for clock_col in clock_cols:
        if clock_col not in long_df.columns:
            warnings.warn(f"Clock column missing from longitudinal file and skipped: {clock_col}")
            continue

        for _, g in long_df.groupby(args.id_col, sort=False):
            rows.append(
                summarize_one_subject_clock(
                    g=g,
                    id_col=args.id_col,
                    time_col=args.time_col,
                    clock_col=clock_col,
                    min_scans=args.min_scans,
                    min_followup_years=args.min_followup_years
                )
            )

    return pd.DataFrame(rows)


# -----------------------------
# 4. Prepare merged analysis dataset
# -----------------------------

def prepare_longitudinal_input(long_df, args):
    d = long_df.copy()

    if args.id_col not in d.columns:
        raise ValueError(f"{args.id_col} is missing from longitudinal file.")

    if args.time_col not in d.columns:
        raise ValueError(f"{args.time_col} is missing from longitudinal file.")

    d[args.id_col] = normalize_id_series(d[args.id_col])
    d[args.time_col] = clean_numeric_series(d[args.time_col])

    # Safety filter: if scan_dx exists, keep only CN rows.
    if "scan_dx" in d.columns:
        d = d.loc[d["scan_dx"].astype(str).str.upper() == "CN"].copy()

    # Safety filter: if scan_relation_to_event exists, remove event/post-event labels.
    if "scan_relation_to_event" in d.columns:
        allowed = {
            "selected_baseline",
            "pre_event_CN_followup",
            "censored_CN_followup",
            "pre_event_followup",
            "censored_followup"
        }
        d = d.loc[d["scan_relation_to_event"].astype(str).isin(allowed)].copy()

    d = d.loc[d[args.time_col].notna()].copy()

    return d


def prepare_baseline_csf_input(csf_df, args, csf_vars):
    d = csf_df.copy()

    if args.id_col not in d.columns:
        raise ValueError(f"{args.id_col} is missing from baseline CSF file.")

    required = [
        args.age_col,
        args.sex_col,
        args.icv_col,
        args.apoe4_col
    ] + csf_vars

    missing = [c for c in required if c not in d.columns]
    if missing:
        raise ValueError(f"Missing required columns in baseline CSF file: {missing}")

    d[args.id_col] = normalize_id_series(d[args.id_col])

    group_col = None

    for c in [args.group_col, "conversion_group_3level", "conversion_group"]:
        if c in d.columns:
            group_col = c
            break

    if group_col is None:
        raise ValueError(
            "Could not find conversion group column in baseline CSF file. "
            "Expected conversion_group_3level or conversion_group."
        )

    d["conversion_group_3level"] = d[group_col].apply(normalize_group)

    valid_groups = ["Non-event & censored", "CN-MCI", "CN-AD"]
    d = d.loc[d["conversion_group_3level"].isin(valid_groups)].copy()

    d["_age_model"] = clean_numeric_series(d[args.age_col])
    d["_sex_male_model"] = clean_numeric_series(d[args.sex_col])
    d["_icv_model"] = clean_numeric_series(d[args.icv_col])
    d["_apoe4_model"] = clean_numeric_series(d[args.apoe4_col])

    for c in csf_vars:
        d[c] = clean_numeric_series(d[c])

    keep_cols = [
        args.id_col,
        "conversion_group_3level",
        "_age_model",
        "_sex_male_model",
        "_icv_model",
        "_apoe4_model"
    ] + csf_vars

    if "csf_match_method" in d.columns:
        keep_cols.append("csf_match_method")
    if "csf_days_from_selected_baseline" in d.columns:
        keep_cols.append("csf_days_from_selected_baseline")

    d = d[keep_cols].drop_duplicates(args.id_col, keep="first").copy()

    return d


def build_analysis_dataset(subject_metrics, baseline_csf, args, csf_vars):
    d = subject_metrics.copy()
    d[args.id_col] = normalize_id_series(d[args.id_col])

    merged = d.merge(
        baseline_csf,
        on=args.id_col,
        how="inner"
    )

    merged["group_CN_MCI_vs_non_event"] = (
        merged["conversion_group_3level"].astype(str) == "CN-MCI"
    ).astype(float)

    merged["group_CN_AD_vs_non_event"] = (
        merged["conversion_group_3level"].astype(str) == "CN-AD"
    ).astype(float)

    # Long format by biomarker.
    rows = []

    for csf_var in csf_vars:
        tmp = merged.copy()
        tmp["csf_var"] = csf_var
        tmp["csf_value"] = clean_numeric_series(tmp[csf_var])
        tmp["csf_pathology_multiplier"] = pathology_direction_multiplier(csf_var)
        tmp["csf_pathology_value"] = tmp["csf_value"] * tmp["csf_pathology_multiplier"]
        rows.append(tmp)

    long = pd.concat(rows, axis=0, ignore_index=True)

    return merged, long


# -----------------------------
# 5. Regression
# -----------------------------

def fit_ols_from_dataframe(df, outcome_col, predictor_col, covariate_cols, min_n):
    needed = [outcome_col, predictor_col] + covariate_cols

    sub = df[needed].copy()

    for c in needed:
        sub[c] = clean_numeric_series(sub[c])

    sub = sub.replace([np.inf, -np.inf], np.nan).dropna()

    n = sub.shape[0]

    if n < min_n:
        return {
            "n": n,
            "status": "skipped_too_few_complete_cases"
        }

    y = sub[outcome_col].values.astype(float)
    predictor = sub[predictor_col].values.astype(float)

    kept_terms = ["intercept", predictor_col]
    X_parts = [np.ones(n), predictor]
    dropped_terms = []

    if np.nanstd(predictor) <= 1e-12:
        return {
            "n": n,
            "status": "skipped_predictor_has_zero_variance",
            "terms_kept": "intercept",
            "terms_dropped_zero_variance": predictor_col
        }

    for c in covariate_cols:
        x = sub[c].values.astype(float)

        if np.nanstd(x) <= 1e-12:
            dropped_terms.append(c)
            continue

        kept_terms.append(c)
        X_parts.append(x)

    X = np.column_stack(X_parts)

    if n <= X.shape[1] + 1:
        return {
            "n": n,
            "status": "skipped_insufficient_residual_df",
            "terms_kept": ",".join(kept_terms),
            "terms_dropped_zero_variance": ",".join(dropped_terms)
        }

    if np.linalg.matrix_rank(X) < X.shape[1]:
        return {
            "n": n,
            "status": "skipped_rank_deficient",
            "terms_kept": ",".join(kept_terms),
            "terms_dropped_zero_variance": ",".join(dropped_terms)
        }

    beta, _, _, _ = np.linalg.lstsq(X, y, rcond=None)

    yhat = X @ beta
    resid = y - yhat

    df_resid = n - X.shape[1]
    sse = float(np.sum(resid ** 2))
    sst = float(np.sum((y - np.mean(y)) ** 2))

    mse = sse / df_resid
    cov_beta = mse * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(cov_beta))

    idx = kept_terms.index(predictor_col)

    t_stat = beta[idx] / se[idx]
    p_value = 2.0 * stats.t.sf(np.abs(t_stat), df=df_resid)

    ci_low = beta[idx] - stats.t.ppf(0.975, df_resid) * se[idx]
    ci_high = beta[idx] + stats.t.ppf(0.975, df_resid) * se[idx]

    r2 = 1.0 - sse / sst if sst > 0 else np.nan
    partial_r2 = (t_stat ** 2) / ((t_stat ** 2) + df_resid)
    partial_r = np.sign(beta[idx]) * np.sqrt(partial_r2)

    sd_y = np.std(y, ddof=1)
    sd_x = np.std(predictor, ddof=1)

    std_beta = beta[idx] * sd_x / sd_y if sd_y > 0 and sd_x > 0 else np.nan

    out = {
        "n": n,
        "status": "ok",
        "beta_predictor": float(beta[idx]),
        "se_predictor": float(se[idx]),
        "ci_low_predictor": float(ci_low),
        "ci_high_predictor": float(ci_high),
        "t_predictor": float(t_stat),
        "p_raw_predictor": float(p_value),
        "std_beta_predictor": float(std_beta),
        "partial_r_predictor": float(partial_r),
        "partial_r2_predictor": float(partial_r2),
        "model_r2": float(r2),
        "df_resid": int(df_resid),
        "mean_outcome": float(np.mean(y)),
        "sd_outcome": float(sd_y),
        "mean_predictor": float(np.mean(predictor)),
        "sd_predictor": float(sd_x),
        "terms_kept": ",".join(kept_terms),
        "terms_dropped_zero_variance": ",".join(dropped_terms)
    }

    for term_name, term_beta, term_se in zip(kept_terms, beta, se):
        safe = sanitize_name(term_name)
        out[f"coef_{safe}"] = float(term_beta)
        out[f"se_{safe}"] = float(term_se)

    return out


def make_result_with_pathology_direction(res, csf_var):
    mult = pathology_direction_multiplier(csf_var)

    if res.get("status") == "ok":
        res["beta_pathology_direction"] = res["beta_predictor"] * mult
        res["std_beta_pathology_direction"] = res["std_beta_predictor"] * mult
        res["partial_r_pathology_direction"] = res["partial_r_predictor"] * mult
    else:
        res["beta_pathology_direction"] = np.nan
        res["std_beta_pathology_direction"] = np.nan
        res["partial_r_pathology_direction"] = np.nan

    res["pathology_direction_note"] = (
        "positive means worse baseline CSF pathology predicts greater clock increase; "
        "Abeta_CSF sign is reversed because lower Abeta is pathological"
    )

    return res


def run_combined_associations(analysis_long, min_n):
    rows = []

    change_metrics = [
        "slope_per_year",
        "annualized_delta",
        "delta_last_minus_baseline"
    ]

    covariates = [
        "baseline_clock_value",
        "group_CN_MCI_vs_non_event",
        "group_CN_AD_vs_non_event",
        "_age_model",
        "_sex_male_model",
        "_icv_model",
        "_apoe4_model",
        "followup_span_years"
    ]

    for clock_col in sorted(analysis_long["clock_col"].dropna().unique()):
        for csf_var in sorted(analysis_long["csf_var"].dropna().unique()):
            sub = analysis_long.loc[
                (analysis_long["clock_col"] == clock_col) &
                (analysis_long["csf_var"] == csf_var) &
                (analysis_long["metric_status"] == "ok")
            ].copy()

            for metric in change_metrics:
                res = fit_ols_from_dataframe(
                    df=sub,
                    outcome_col=metric,
                    predictor_col="csf_value",
                    covariate_cols=covariates,
                    min_n=min_n
                )

                res.update({
                    "analysis_type": "combined_group_adjusted",
                    "group": "All groups combined",
                    "group_reference": "Non-event & censored",
                    "clock_col": clock_col,
                    "change_metric": metric,
                    "csf_var": csf_var,
                    "predictor": "baseline_csf_value",
                    "model_formula": (
                        f"{metric} ~ baseline_CSF + baseline_clock_value + "
                        "conversion_group + Age + Sex + ICV + APOE4 + followup_span_years"
                    ),
                    "n_non_event_censored": int(
                        sub.loc[
                            sub["conversion_group_3level"] == "Non-event & censored"
                        ][metric].notna().sum()
                    ),
                    "n_cn_mci": int(
                        sub.loc[
                            sub["conversion_group_3level"] == "CN-MCI"
                        ][metric].notna().sum()
                    ),
                    "n_cn_ad": int(
                        sub.loc[
                            sub["conversion_group_3level"] == "CN-AD"
                        ][metric].notna().sum()
                    )
                })

                res = make_result_with_pathology_direction(res, csf_var)
                rows.append(res)

    return pd.DataFrame(rows)


def run_group_specific_associations(analysis_long, min_n):
    rows = []

    change_metrics = [
        "slope_per_year",
        "annualized_delta",
        "delta_last_minus_baseline"
    ]

    groups = [
        "Non-event & censored",
        "CN-MCI",
        "CN-AD"
    ]

    covariates = [
        "baseline_clock_value",
        "_age_model",
        "_sex_male_model",
        "_icv_model",
        "_apoe4_model",
        "followup_span_years"
    ]

    for clock_col in sorted(analysis_long["clock_col"].dropna().unique()):
        for csf_var in sorted(analysis_long["csf_var"].dropna().unique()):
            for group_name in groups:
                sub = analysis_long.loc[
                    (analysis_long["clock_col"] == clock_col) &
                    (analysis_long["csf_var"] == csf_var) &
                    (analysis_long["conversion_group_3level"] == group_name) &
                    (analysis_long["metric_status"] == "ok")
                ].copy()

                for metric in change_metrics:
                    res = fit_ols_from_dataframe(
                        df=sub,
                        outcome_col=metric,
                        predictor_col="csf_value",
                        covariate_cols=covariates,
                        min_n=min_n
                    )

                    res.update({
                        "analysis_type": "group_specific",
                        "group": group_name,
                        "group_reference": np.nan,
                        "clock_col": clock_col,
                        "change_metric": metric,
                        "csf_var": csf_var,
                        "predictor": "baseline_csf_value",
                        "model_formula": (
                            f"{metric} ~ baseline_CSF + baseline_clock_value + "
                            "Age + Sex + ICV + APOE4 + followup_span_years"
                        )
                    })

                    res = make_result_with_pathology_direction(res, csf_var)
                    rows.append(res)

    return pd.DataFrame(rows)


def add_bh_adjustment(df):
    df = df.copy()

    if "p_raw_predictor" not in df.columns:
        df["p_bh_all_tests"] = np.nan
        return df

    mask = df["status"].eq("ok") & df["p_raw_predictor"].notna()
    df["p_bh_all_tests"] = np.nan

    if mask.sum() > 0:
        df.loc[mask, "p_bh_all_tests"] = adjust_p_bh(
            df.loc[mask, "p_raw_predictor"].values
        )

    return df


# -----------------------------
# 6. Main
# -----------------------------

def main():
    args = parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    clock_cols = parse_list_arg(args.clock_cols)
    csf_vars = parse_list_arg(args.csf_vars)

    log("============================================================")
    log("ADNI longitudinal AD L'EPOCH change vs baseline CSF")
    log("============================================================")
    log(f"Longitudinal file: {args.longitudinal_file}")
    log(f"Baseline CSF file: {args.baseline_csf_file}")
    log(f"Output directory: {outdir}")
    log(f"Clock columns: {clock_cols}")
    log(f"CSF variables: {csf_vars}")
    log("============================================================")

    long_raw = read_table(args.longitudinal_file)
    csf_raw = read_table(args.baseline_csf_file)

    log(f"Longitudinal rows: {long_raw.shape[0]}, columns: {long_raw.shape[1]}")
    log(f"Baseline CSF rows: {csf_raw.shape[0]}, columns: {csf_raw.shape[1]}")

    long_df = prepare_longitudinal_input(long_raw, args)
    baseline_csf = prepare_baseline_csf_input(csf_raw, args, csf_vars)

    subject_metrics = build_subject_clock_metrics(
        long_df=long_df,
        args=args,
        clock_cols=clock_cols
    )

    subject_metrics.to_csv(
        outdir / f"{args.prefix}_longitudinal_clock_change_subject_metrics.tsv",
        sep="\t",
        index=False
    )

    merged_wide, analysis_long = build_analysis_dataset(
        subject_metrics=subject_metrics,
        baseline_csf=baseline_csf,
        args=args,
        csf_vars=csf_vars
    )

    merged_wide.to_csv(
        outdir / f"{args.prefix}_longitudinal_change_vs_baseline_csf_subject_dataset_wide.tsv",
        sep="\t",
        index=False
    )

    analysis_long.to_csv(
        outdir / f"{args.prefix}_longitudinal_change_vs_baseline_csf_analysis_dataset.tsv",
        sep="\t",
        index=False
    )

    subject_summary = analysis_long.groupby(
        ["clock_col", "metric_status", "conversion_group_3level"],
        dropna=False
    ).agg(
        n_rows=(args.id_col, "size"),
        n_subjects=(args.id_col, "nunique"),
        median_n_scans=("n_scans", "median"),
        median_followup_span_years=("followup_span_years", "median")
    ).reset_index()

    subject_summary.to_csv(
        outdir / f"{args.prefix}_longitudinal_change_vs_baseline_csf_subject_summary.tsv",
        sep="\t",
        index=False
    )

    combined_assoc = run_combined_associations(
        analysis_long=analysis_long,
        min_n=args.min_n
    )

    group_assoc = run_group_specific_associations(
        analysis_long=analysis_long,
        min_n=args.min_n
    )

    combined_assoc = add_bh_adjustment(combined_assoc)
    group_assoc = add_bh_adjustment(group_assoc)

    all_assoc = pd.concat(
        [combined_assoc, group_assoc],
        axis=0,
        ignore_index=True,
        sort=False
    )
    all_assoc = add_bh_adjustment(all_assoc)

    combined_assoc.to_csv(
        outdir / f"{args.prefix}_longitudinal_change_vs_baseline_csf_combined_associations.tsv",
        sep="\t",
        index=False
    )

    group_assoc.to_csv(
        outdir / f"{args.prefix}_longitudinal_change_vs_baseline_csf_group_specific_associations.tsv",
        sep="\t",
        index=False
    )

    all_assoc.to_csv(
        outdir / f"{args.prefix}_longitudinal_change_vs_baseline_csf_all_associations.tsv",
        sep="\t",
        index=False
    )

    summary = {
        "longitudinal_file": args.longitudinal_file,
        "baseline_csf_file": args.baseline_csf_file,
        "n_longitudinal_rows_after_safety_filters": int(long_df.shape[0]),
        "n_baseline_csf_subjects": int(baseline_csf[args.id_col].nunique()),
        "n_subject_clock_metric_rows": int(subject_metrics.shape[0]),
        "n_merged_subject_clock_rows": int(merged_wide.shape[0]),
        "n_analysis_long_rows": int(analysis_long.shape[0]),
        "clock_cols": clock_cols,
        "csf_vars": csf_vars,
        "change_metrics": [
            "slope_per_year",
            "annualized_delta",
            "delta_last_minus_baseline"
        ],
        "primary_combined_model": (
            "clock_change_metric ~ baseline_CSF + baseline_clock_value + "
            "conversion_group + Age + Sex + ICV + APOE4 + followup_span_years"
        ),
        "secondary_group_specific_model": (
            "clock_change_metric ~ baseline_CSF + baseline_clock_value + "
            "Age + Sex + ICV + APOE4 + followup_span_years"
        ),
        "group_reference": "Non-event & censored",
        "min_scans": args.min_scans,
        "min_followup_years": args.min_followup_years,
        "min_n": args.min_n,
        "outputs": {
            "subject_clock_metrics": str(
                outdir / f"{args.prefix}_longitudinal_clock_change_subject_metrics.tsv"
            ),
            "analysis_dataset_long": str(
                outdir / f"{args.prefix}_longitudinal_change_vs_baseline_csf_analysis_dataset.tsv"
            ),
            "combined_associations": str(
                outdir / f"{args.prefix}_longitudinal_change_vs_baseline_csf_combined_associations.tsv"
            ),
            "group_specific_associations": str(
                outdir / f"{args.prefix}_longitudinal_change_vs_baseline_csf_group_specific_associations.tsv"
            ),
            "all_associations": str(
                outdir / f"{args.prefix}_longitudinal_change_vs_baseline_csf_all_associations.tsv"
            ),
            "subject_summary": str(
                outdir / f"{args.prefix}_longitudinal_change_vs_baseline_csf_subject_summary.tsv"
            )
        }
    }

    with open(
        outdir / f"{args.prefix}_longitudinal_change_vs_baseline_csf_analysis_summary.json",
        "w"
    ) as f:
        json.dump(summary, f, indent=2)

    log("============================================================")
    log("Done.")
    log(f"Subject-level metrics: {summary['outputs']['subject_clock_metrics']}")
    log(f"Combined associations: {summary['outputs']['combined_associations']}")
    log(f"Group-specific associations: {summary['outputs']['group_specific_associations']}")
    log(f"All associations: {summary['outputs']['all_associations']}")
    log("Subject summary:")
    log(subject_summary.to_string(index=False))
    log("============================================================")


if __name__ == "__main__":
    main()