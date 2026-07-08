#!/usr/bin/env python3
# ============================================================
# ADNI baseline AD L'EPOCH acceleration_z association with CSF biomarkers
#
# Variables used
# --------------
# Required subject / visit variables:
#   PTID
#   Visit_Code
#   Date
#   DX_Binary
#
# Main predictor:
#   adni_brain_mri_ad_lepoch_acceleration_z
#
# Outcome groups:
#   Non-event & censored
#   CN-MCI
#   CN-AD
#
# Covariates:
#   Age
#   Sex
#   DLICV or ICV
#   APOE4
#
# CSF biomarker outcomes:
#   Abeta_CSF
#   Tau_CSF
#   PTau_CSF
#
# Pipeline workflow
# -----------------
# 1. Load the baseline AD L'EPOCH prediction file.
# 2. Load the full ADNI iSTAGING table.
# 3. Construct the 3 conversion groups from event/event_or_censor_dx:
#      Non-event & censored
#      CN-MCI
#      CN-AD
# 4. Attach baseline CSF biomarkers to each selected baseline MRI row.
#      Priority:
#        a) same PTID + selected baseline visit
#        b) same PTID + ADNI baseline visit
#        c) nearest CSF visit to selected MRI baseline date within the specified window
# 5. Derive APOE4 from APOE4, APOE genotype, or APGEN1/APGEN2 when available.
# 6. For each CSF biomarker and each group, fit:
#      CSF biomarker ~ acceleration_z + Age + Sex + ICV + APOE4
# 7. Save analysis dataset, variable map, group summaries, and association results.
#
# Major output files
# ------------------
#   {prefix}_baseline_csf_analysis_dataset.tsv
#   {prefix}_baseline_csf_variable_map.tsv
#   {prefix}_baseline_csf_group_summary.tsv
#   {prefix}_baseline_csf_associations.tsv
#   {prefix}_baseline_csf_csf_match_summary.tsv
#   {prefix}_baseline_csf_analysis_summary.json
#
# Notes
# -----
# - Raw P-values are always saved.
# - BH-adjusted P-values across all CSF tests are also saved for reference.
# - If a covariate has zero variance within a group, it is dropped only for that
#   group-specific regression to avoid rank deficiency; dropped covariates are recorded.
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

    p.add_argument("--predictions-file", required=True)
    p.add_argument("--adni-file", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--prefix", default="adni_brain_mri_ad_lepoch")

    p.add_argument("--id-col", default="PTID")
    p.add_argument("--visit-col", default="Visit_Code")
    p.add_argument("--date-col", default="Date")
    p.add_argument("--dx-col", default="DX_Binary")

    p.add_argument("--accel-col", default="adni_brain_mri_ad_lepoch_acceleration_z")
    p.add_argument("--age-col", default="Age")
    p.add_argument("--sex-col", default="Sex")
    p.add_argument("--icv-col", default="DLICV")
    p.add_argument("--apoe4-col", default="APOE4")

    p.add_argument(
        "--csf-vars",
        default="Abeta_CSF,Tau_CSF,PTau_CSF",
        help="Comma-separated CSF biomarker variables."
    )

    p.add_argument(
        "--max-csf-baseline-distance-days",
        type=float,
        default=365.0,
        help="Maximum allowed days between selected MRI baseline and fallback nearest CSF visit."
    )

    p.add_argument("--min-n", type=int, default=8)

    return p.parse_args()


# -----------------------------
# 2. Generic helpers
# -----------------------------

def log(msg):
    print(msg, flush=True)


def read_table(path):
    path = Path(path)
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path, low_memory=False)
    return pd.read_csv(path, sep="\t", low_memory=False)


def sanitize_name(x):
    return re.sub(r"[^a-z0-9]+", "", str(x).lower())


def parse_list_arg(x):
    if x is None:
        return []
    x = str(x).strip()
    if x == "" or x.lower() == "none":
        return []
    return [v.strip() for v in x.split(",") if v.strip()]


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


def normalize_visit_series(s):
    return s.astype(str).str.strip().str.lower()


def normalize_dx(x):
    if pd.isna(x):
        return np.nan

    x = str(x).strip().upper()

    dx_map = {
        "CN": "CN",
        "NL": "CN",
        "NORMAL": "CN",
        "CONTROL": "CN",
        "HC": "CN",
        "MCI": "MCI",
        "EMCI": "MCI",
        "LMCI": "MCI",
        "AD": "AD",
        "DEMENTIA": "AD"
    }

    return dx_map.get(x, x)


def coerce_event_bool(s):
    return s.astype(str).str.lower().isin(["true", "1", "yes", "y"])


def resolve_column(df, requested, aliases=None, required=False):
    aliases = aliases or []
    candidates = [requested] + aliases

    cols = list(df.columns)
    lower_map = {c.lower(): c for c in cols}
    sanit_map = {sanitize_name(c): c for c in cols}

    for c in candidates:
        if c in df.columns:
            return c
        if str(c).lower() in lower_map:
            return lower_map[str(c).lower()]
        key = sanitize_name(c)
        if key in sanit_map:
            return sanit_map[key]

    if required:
        raise ValueError(f"Could not find required column: {requested}")

    return None


def parse_date_series(s):
    return pd.to_datetime(s, errors="coerce")


def first_nonmissing(x):
    x = x.dropna()
    if len(x) == 0:
        return np.nan
    return x.iloc[0]


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


# -----------------------------
# 3. Conversion group
# -----------------------------

def add_conversion_group(df):
    df = df.copy()

    if "conversion_group" in df.columns:
        existing = df["conversion_group"].astype(str)

        df["conversion_group_3level"] = np.select(
            [
                existing.str.contains("MCI", case=False, na=False),
                existing.str.contains("AD", case=False, na=False),
                existing.str.contains("censored|non", case=False, na=False)
            ],
            [
                "CN-MCI",
                "CN-AD",
                "Non-event & censored"
            ],
            default=existing
        )

        return df

    if "event" not in df.columns:
        raise ValueError("Neither conversion_group nor event column is present in prediction file.")

    event = coerce_event_bool(df["event"])

    if "event_or_censor_dx" in df.columns:
        end_dx = df["event_or_censor_dx"].apply(normalize_dx)
    elif "event_dx" in df.columns:
        end_dx = df["event_dx"].apply(normalize_dx)
    else:
        end_dx = pd.Series(np.nan, index=df.index)

    df["conversion_group_3level"] = np.where(
        ~event,
        "Non-event & censored",
        np.where(
            end_dx == "MCI",
            "CN-MCI",
            np.where(end_dx == "AD", "CN-AD", "Other event")
        )
    )

    return df


# -----------------------------
# 4. CSF matching
# -----------------------------

def is_baseline_visit(v):
    if pd.isna(v):
        return False
    x = str(v).strip().lower()
    return x in {"bl", "base", "baseline", "m00", "m0"}


def attach_csf_to_predictions(pred, raw, args, csf_vars):
    pred = pred.copy()
    raw = raw.copy()

    id_col = args.id_col
    visit_col = args.visit_col
    date_col = args.date_col

    pred["_row_id"] = np.arange(pred.shape[0])
    pred["_pid"] = normalize_id_series(pred[id_col])

    raw["_pid"] = normalize_id_series(raw[id_col])

    if visit_col in raw.columns:
        raw["_visit_norm"] = normalize_visit_series(raw[visit_col])
    else:
        raw["_visit_norm"] = np.nan

    if date_col in raw.columns:
        raw["_date"] = parse_date_series(raw[date_col])
    else:
        raw["_date"] = pd.NaT

    for c in csf_vars:
        raw[c] = clean_numeric_series(raw[c])

    raw["_n_csf_nonmissing"] = raw[csf_vars].notna().sum(axis=1)
    raw_csf = raw.loc[raw["_n_csf_nonmissing"] > 0].copy()

    if raw_csf.empty:
        warnings.warn("No rows in raw ADNI file have non-missing CSF biomarkers.")

    selected_rows = []

    baseline_visit_col = None
    for c in ["selected_baseline_visit_code", visit_col]:
        if c in pred.columns:
            baseline_visit_col = c
            break

    baseline_date_col = None
    for c in ["selected_baseline_date", date_col]:
        if c in pred.columns:
            baseline_date_col = c
            break

    for _, row in pred.iterrows():
        pid = row["_pid"]
        row_id = row["_row_id"]

        baseline_visit = np.nan
        if baseline_visit_col is not None:
            baseline_visit = row[baseline_visit_col]

        baseline_visit_norm = str(baseline_visit).strip().lower() if pd.notna(baseline_visit) else ""

        baseline_date = pd.NaT
        if baseline_date_col is not None:
            baseline_date = pd.to_datetime(row[baseline_date_col], errors="coerce")

        candidates = raw_csf.loc[raw_csf["_pid"] == pid].copy()

        out = {
            "_row_id": row_id,
            "csf_match_method": "no_csf_available",
            "csf_match_visit": np.nan,
            "csf_match_date": pd.NaT,
            "csf_days_from_selected_baseline": np.nan,
            "n_csf_biomarkers_nonmissing": 0
        }

        for c in csf_vars:
            out[c] = np.nan

        if candidates.empty:
            selected_rows.append(out)
            continue

        chosen = None
        method = None

        # Priority 1: same selected baseline visit.
        if baseline_visit_norm:
            same_visit = candidates.loc[candidates["_visit_norm"] == baseline_visit_norm].copy()
            if not same_visit.empty:
                chosen = same_visit.iloc[0]
                method = "same_selected_baseline_visit"

        # Priority 2: ADNI baseline visit.
        if chosen is None:
            baseline_rows = candidates.loc[candidates["_visit_norm"].apply(is_baseline_visit)].copy()
            if not baseline_rows.empty:
                if baseline_rows["_date"].notna().any():
                    baseline_rows = baseline_rows.sort_values("_date", kind="mergesort")
                chosen = baseline_rows.iloc[0]
                method = "adni_baseline_visit"

        # Priority 3: nearest CSF visit to selected MRI baseline date.
        if chosen is None and pd.notna(baseline_date) and candidates["_date"].notna().any():
            candidates["_abs_days_from_selected_baseline"] = (
                candidates["_date"] - baseline_date
            ).abs().dt.days

            nearest = candidates.loc[
                candidates["_abs_days_from_selected_baseline"] <= args.max_csf_baseline_distance_days
            ].copy()

            if not nearest.empty:
                nearest = nearest.sort_values("_abs_days_from_selected_baseline", kind="mergesort")
                chosen = nearest.iloc[0]
                method = "nearest_csf_visit_within_window"

        # Priority 4: first available CSF visit, but flagged as not baseline-proximal.
        if chosen is None:
            tmp = candidates.copy()
            if tmp["_date"].notna().any():
                tmp = tmp.sort_values("_date", kind="mergesort")
            chosen = tmp.iloc[0]
            method = "first_available_csf_visit_no_window_match"

        out["csf_match_method"] = method
        out["csf_match_visit"] = chosen[visit_col] if visit_col in chosen.index else np.nan
        out["csf_match_date"] = chosen["_date"] if "_date" in chosen.index else pd.NaT

        if pd.notna(baseline_date) and pd.notna(out["csf_match_date"]):
            out["csf_days_from_selected_baseline"] = (
                out["csf_match_date"] - baseline_date
            ).days

        out["n_csf_biomarkers_nonmissing"] = int(chosen[csf_vars].notna().sum())

        for c in csf_vars:
            out[c] = chosen[c]

        selected_rows.append(out)

    csf_attach = pd.DataFrame(selected_rows)

    merged = pred.merge(csf_attach, on="_row_id", how="left")
    merged = merged.drop(columns=["_row_id", "_pid"], errors="ignore")

    return merged, csf_attach


# -----------------------------
# 5. APOE4 handling
# -----------------------------

def derive_apoe4_by_subject(raw, id_col, requested_col):
    raw = raw.copy()
    raw["_pid"] = normalize_id_series(raw[id_col])

    explicit_aliases = [
        "APOE4",
        "APOE4_bl",
        "APOE4_STATUS",
        "APOE4STATUS",
        "APOE_e4",
        "APOE_E4",
        "APOE4_COUNT",
        "APOE4_count"
    ]

    col = resolve_column(raw, requested_col, aliases=explicit_aliases, required=False)

    if col is not None:
        tmp = raw[["_pid", col]].copy()
        tmp["_apoe4"] = clean_numeric_series(tmp[col])

        apoe = tmp.groupby("_pid", as_index=False).agg(
            _apoe4=("_apoe4", first_nonmissing)
        )
        apoe["_apoe4_source_column"] = col
        apoe["_apoe4_method"] = "numeric_APOE4_column"
        return apoe

    genotype_aliases = [
        "APOE",
        "APOEGEN",
        "APOE_GENOTYPE",
        "APOE Genotype",
        "APOE_Genotype",
        "apoe_genotype"
    ]

    geno_col = resolve_column(raw, "APOEGEN", aliases=genotype_aliases, required=False)

    if geno_col is not None:
        tmp = raw[["_pid", geno_col]].copy()
        tmp["_geno"] = tmp[geno_col].astype(str)
        tmp["_apoe4"] = tmp["_geno"].apply(lambda z: len(re.findall(r"4", z))).astype(float)

        apoe = tmp.groupby("_pid", as_index=False).agg(
            _apoe4=("_apoe4", first_nonmissing)
        )
        apoe["_apoe4_source_column"] = geno_col
        apoe["_apoe4_method"] = "counted_allele_4_from_genotype_string"
        return apoe

    allele1 = resolve_column(raw, "APGEN1", aliases=["APOE_A1", "APOE1"], required=False)
    allele2 = resolve_column(raw, "APGEN2", aliases=["APOE_A2", "APOE2"], required=False)

    if allele1 is not None and allele2 is not None:
        tmp = raw[["_pid", allele1, allele2]].copy()
        a1 = clean_numeric_series(tmp[allele1])
        a2 = clean_numeric_series(tmp[allele2])
        tmp["_apoe4"] = (a1 == 4).astype(float) + (a2 == 4).astype(float)
        tmp.loc[a1.isna() | a2.isna(), "_apoe4"] = np.nan

        apoe = tmp.groupby("_pid", as_index=False).agg(
            _apoe4=("_apoe4", first_nonmissing)
        )
        apoe["_apoe4_source_column"] = f"{allele1}+{allele2}"
        apoe["_apoe4_method"] = "counted_allele_4_from_two_allele_columns"
        return apoe

    raise ValueError(
        "Could not identify APOE4. Expected APOE4, APOE genotype, or APGEN1/APGEN2 columns."
    )


# -----------------------------
# 6. Core analysis dataset
# -----------------------------

def prepare_analysis_dataset(pred_csf, raw, args, csf_vars):
    d = pred_csf.copy()
    d = add_conversion_group(d)

    d["_pid"] = normalize_id_series(d[args.id_col])

    apoe = derive_apoe4_by_subject(
        raw=raw,
        id_col=args.id_col,
        requested_col=args.apoe4_col
    )

    d = d.merge(apoe, on="_pid", how="left")

    accel_col = resolve_column(d, args.accel_col, required=True)
    age_col = resolve_column(
        d,
        args.age_col,
        aliases=["AGE", "age_at_scan_used_for_model", "baseline_age_raw"],
        required=True
    )
    sex_col = resolve_column(d, args.sex_col, aliases=["SEX"], required=True)
    icv_col = resolve_column(
        d,
        args.icv_col,
        aliases=["ICV", "DLICV_baseline", "IntracranialVolume", "ICV_baseline"],
        required=True
    )

    d["_acceleration_z"] = clean_numeric_series(d[accel_col])
    d["_age"] = clean_numeric_series(d[age_col])
    d["_icv"] = clean_numeric_series(d[icv_col])
    d["_apoe4"] = clean_numeric_series(d["_apoe4"])

    sex = d[sex_col].astype(str).str.strip()
    sex_norm = sex.replace({
        "0": "Female",
        "0.0": "Female",
        "1": "Male",
        "1.0": "Male",
        "F": "Female",
        "M": "Male",
        "female": "Female",
        "male": "Male",
        "Female": "Female",
        "Male": "Male"
    })

    d["_sex"] = sex_norm
    d["_sex_male"] = np.where(
        sex_norm == "Male",
        1.0,
        np.where(sex_norm == "Female", 0.0, np.nan)
    )

    variable_map = pd.DataFrame([
        {
            "role": "main_predictor",
            "requested": args.accel_col,
            "resolved_column": accel_col,
            "method": "exact_or_alias"
        },
        {
            "role": "covariate_age",
            "requested": args.age_col,
            "resolved_column": age_col,
            "method": "exact_or_alias"
        },
        {
            "role": "covariate_sex",
            "requested": args.sex_col,
            "resolved_column": sex_col,
            "method": "exact_or_alias"
        },
        {
            "role": "covariate_icv",
            "requested": args.icv_col,
            "resolved_column": icv_col,
            "method": "exact_or_alias"
        },
        {
            "role": "covariate_apoe4",
            "requested": args.apoe4_col,
            "resolved_column": str(d["_apoe4_source_column"].dropna().iloc[0]) if d["_apoe4_source_column"].notna().any() else "NA",
            "method": str(d["_apoe4_method"].dropna().iloc[0]) if d["_apoe4_method"].notna().any() else "NA"
        }
    ])

    for c in csf_vars:
        variable_map = pd.concat([
            variable_map,
            pd.DataFrame([{
                "role": "csf_outcome",
                "requested": c,
                "resolved_column": c,
                "method": "exact",
                "n_nonmissing": int(clean_numeric_series(d[c]).notna().sum())
            }])
        ], axis=0, ignore_index=True)

    return d, variable_map


# -----------------------------
# 7. Regression
# -----------------------------

def fit_ols_one_group(d, outcome, group_name, min_n):
    sub = d.loc[d["conversion_group_3level"] == group_name].copy()

    sub["_outcome"] = clean_numeric_series(sub[outcome])

    model_cols = [
        "_outcome",
        "_acceleration_z",
        "_age",
        "_sex_male",
        "_icv",
        "_apoe4"
    ]

    sub = sub[model_cols].replace([np.inf, -np.inf], np.nan).dropna()

    n = sub.shape[0]

    if n < min_n:
        return {
            "group": group_name,
            "outcome": outcome,
            "n": n,
            "status": "skipped_too_few_complete_cases"
        }

    y = sub["_outcome"].values.astype(float)

    predictors = [
        ("acceleration_z", sub["_acceleration_z"].values.astype(float), True),
        ("age", sub["_age"].values.astype(float), False),
        ("sex_male", sub["_sex_male"].values.astype(float), False),
        ("icv", sub["_icv"].values.astype(float), False),
        ("apoe4", sub["_apoe4"].values.astype(float), False)
    ]

    kept_terms = ["intercept"]
    X_parts = [np.ones(n)]
    dropped_terms = []

    for name, values, required_predictor in predictors:
        if np.nanstd(values) <= 1e-12:
            if required_predictor:
                return {
                    "group": group_name,
                    "outcome": outcome,
                    "n": n,
                    "status": "skipped_acceleration_z_has_zero_variance"
                }
            dropped_terms.append(name)
            continue

        kept_terms.append(name)
        X_parts.append(values)

    X = np.column_stack(X_parts)

    if n <= X.shape[1] + 1:
        return {
            "group": group_name,
            "outcome": outcome,
            "n": n,
            "status": "skipped_insufficient_residual_df",
            "terms_kept": ",".join(kept_terms),
            "terms_dropped_zero_variance": ",".join(dropped_terms)
        }

    if np.linalg.matrix_rank(X) < X.shape[1]:
        return {
            "group": group_name,
            "outcome": outcome,
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

    idx = kept_terms.index("acceleration_z")
    t_stat = beta[idx] / se[idx]
    p_value = 2.0 * stats.t.sf(np.abs(t_stat), df=df_resid)

    ci_low = beta[idx] - stats.t.ppf(0.975, df_resid) * se[idx]
    ci_high = beta[idx] + stats.t.ppf(0.975, df_resid) * se[idx]

    r2 = 1.0 - sse / sst if sst > 0 else np.nan
    partial_r2 = (t_stat ** 2) / ((t_stat ** 2) + df_resid)
    partial_r = np.sign(beta[idx]) * np.sqrt(partial_r2)

    x_accel = sub["_acceleration_z"].values.astype(float)
    sd_y = np.std(y, ddof=1)
    sd_x = np.std(x_accel, ddof=1)

    std_beta = beta[idx] * sd_x / sd_y if sd_y > 0 and sd_x > 0 else np.nan

    return {
        "group": group_name,
        "outcome": outcome,
        "n": n,
        "status": "ok",
        "beta_acceleration_z": beta[idx],
        "se_acceleration_z": se[idx],
        "ci_low_acceleration_z": ci_low,
        "ci_high_acceleration_z": ci_high,
        "t_acceleration_z": t_stat,
        "p_raw_acceleration_z": p_value,
        "std_beta_acceleration_z": std_beta,
        "partial_r_acceleration_z": partial_r,
        "partial_r2_acceleration_z": partial_r2,
        "model_r2": r2,
        "df_resid": df_resid,
        "mean_outcome": float(np.mean(y)),
        "sd_outcome": float(sd_y),
        "mean_acceleration_z": float(np.mean(x_accel)),
        "sd_acceleration_z": float(sd_x),
        "terms_kept": ",".join(kept_terms),
        "terms_dropped_zero_variance": ",".join(dropped_terms),
        "model_formula": "CSF ~ acceleration_z + Age + Sex + ICV + APOE4"
    }


def run_associations(d, csf_vars, min_n):
    groups = ["Non-event & censored", "CN-MCI", "CN-AD"]
    rows = []

    for outcome in csf_vars:
        for group_name in groups:
            rows.append(
                fit_ols_one_group(
                    d=d,
                    outcome=outcome,
                    group_name=group_name,
                    min_n=min_n
                )
            )

    out = pd.DataFrame(rows)

    if "p_raw_acceleration_z" in out.columns:
        mask = out["status"].eq("ok") & out["p_raw_acceleration_z"].notna()
        out["p_bh_within_all_csf_tests"] = np.nan
        if mask.sum() > 0:
            out.loc[mask, "p_bh_within_all_csf_tests"] = adjust_p_bh(
                out.loc[mask, "p_raw_acceleration_z"].values
            )

    return out


# -----------------------------
# 8. Main
# -----------------------------

def main():
    args = parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    csf_vars = parse_list_arg(args.csf_vars)

    log("============================================================")
    log("ADNI baseline L'EPOCH acceleration_z vs CSF biomarkers")
    log("============================================================")
    log(f"Predictions file: {args.predictions_file}")
    log(f"ADNI file: {args.adni_file}")
    log(f"Output directory: {outdir}")
    log(f"CSF variables: {csf_vars}")
    log("============================================================")

    pred = read_table(args.predictions_file)
    raw = read_table(args.adni_file)

    log(f"Prediction rows: {pred.shape[0]}, columns: {pred.shape[1]}")
    log(f"Raw ADNI rows: {raw.shape[0]}, columns: {raw.shape[1]}")

    missing_pred_cols = [
        c for c in [args.id_col, args.accel_col, args.age_col, args.sex_col, args.icv_col]
        if c not in pred.columns
    ]
    if missing_pred_cols:
        raise ValueError(f"Missing required columns in prediction file: {missing_pred_cols}")

    missing_raw_cols = [c for c in [args.id_col] + csf_vars if c not in raw.columns]
    if missing_raw_cols:
        raise ValueError(f"Missing required columns in ADNI file: {missing_raw_cols}")

    pred_csf, csf_attach = attach_csf_to_predictions(
        pred=pred,
        raw=raw,
        args=args,
        csf_vars=csf_vars
    )

    analysis_df, variable_map = prepare_analysis_dataset(
        pred_csf=pred_csf,
        raw=raw,
        args=args,
        csf_vars=csf_vars
    )

    keep_cols = [
        args.id_col,
        "split" if "split" in analysis_df.columns else None,
        args.visit_col if args.visit_col in analysis_df.columns else None,
        args.date_col if args.date_col in analysis_df.columns else None,
        args.dx_col if args.dx_col in analysis_df.columns else None,
        "selected_baseline_visit_code" if "selected_baseline_visit_code" in analysis_df.columns else None,
        "selected_baseline_date" if "selected_baseline_date" in analysis_df.columns else None,
        "event" if "event" in analysis_df.columns else None,
        "event_or_censor_dx" if "event_or_censor_dx" in analysis_df.columns else None,
        "conversion_group_3level",
        "csf_match_method",
        "csf_match_visit",
        "csf_match_date",
        "csf_days_from_selected_baseline",
        "n_csf_biomarkers_nonmissing",
        "_acceleration_z",
        "_age",
        "_sex",
        "_sex_male",
        "_icv",
        "_apoe4"
    ] + csf_vars

    keep_cols = [c for c in keep_cols if c is not None and c in analysis_df.columns]

    analysis_out = analysis_df[keep_cols].copy()

    analysis_out.to_csv(
        outdir / f"{args.prefix}_baseline_csf_analysis_dataset.tsv",
        sep="\t",
        index=False
    )

    variable_map.to_csv(
        outdir / f"{args.prefix}_baseline_csf_variable_map.tsv",
        sep="\t",
        index=False
    )

    csf_match_summary = analysis_df.groupby(
        ["conversion_group_3level", "csf_match_method"],
        dropna=False
    ).agg(
        n_rows=(args.id_col, "size"),
        n_subjects=(args.id_col, "nunique"),
        median_abs_days_from_selected_baseline=(
            "csf_days_from_selected_baseline",
            lambda x: np.nanmedian(np.abs(pd.to_numeric(x, errors="coerce")))
        )
    ).reset_index()

    csf_match_summary.to_csv(
        outdir / f"{args.prefix}_baseline_csf_csf_match_summary.tsv",
        sep="\t",
        index=False
    )

    group_summary = analysis_df.groupby(
        "conversion_group_3level",
        dropna=False
    ).agg(
        n_subjects=(args.id_col, "nunique"),
        n_rows=(args.id_col, "size"),
        n_accel_nonmissing=("_acceleration_z", lambda x: int(pd.notna(x).sum())),
        mean_acceleration_z=("_acceleration_z", "mean"),
        sd_acceleration_z=("_acceleration_z", "std"),
        n_age_nonmissing=("_age", lambda x: int(pd.notna(x).sum())),
        n_icv_nonmissing=("_icv", lambda x: int(pd.notna(x).sum())),
        n_apoe4_nonmissing=("_apoe4", lambda x: int(pd.notna(x).sum())),
        n_abeta_csf_nonmissing=("Abeta_CSF", lambda x: int(pd.notna(x).sum())) if "Abeta_CSF" in analysis_df.columns else (args.id_col, "size"),
        n_tau_csf_nonmissing=("Tau_CSF", lambda x: int(pd.notna(x).sum())) if "Tau_CSF" in analysis_df.columns else (args.id_col, "size"),
        n_ptau_csf_nonmissing=("PTau_CSF", lambda x: int(pd.notna(x).sum())) if "PTau_CSF" in analysis_df.columns else (args.id_col, "size"),
    ).reset_index()

    group_summary.to_csv(
        outdir / f"{args.prefix}_baseline_csf_group_summary.tsv",
        sep="\t",
        index=False
    )

    assoc = run_associations(
        d=analysis_df,
        csf_vars=csf_vars,
        min_n=args.min_n
    )

    assoc.to_csv(
        outdir / f"{args.prefix}_baseline_csf_associations.tsv",
        sep="\t",
        index=False
    )

    summary = {
        "predictions_file": args.predictions_file,
        "adni_file": args.adni_file,
        "n_rows_prediction_file": int(pred.shape[0]),
        "n_rows_analysis_dataset": int(analysis_out.shape[0]),
        "csf_variables": csf_vars,
        "groups": ["Non-event & censored", "CN-MCI", "CN-AD"],
        "main_predictor": args.accel_col,
        "covariates": ["Age", "Sex", "ICV", "APOE4"],
        "model": "CSF biomarker ~ acceleration_z + Age + Sex + ICV + APOE4",
        "csf_matching_priority": [
            "same_selected_baseline_visit",
            "adni_baseline_visit",
            "nearest_csf_visit_within_window",
            "first_available_csf_visit_no_window_match"
        ],
        "max_csf_baseline_distance_days": args.max_csf_baseline_distance_days,
        "min_n": args.min_n,
        "major_outputs": {
            "analysis_dataset": str(outdir / f"{args.prefix}_baseline_csf_analysis_dataset.tsv"),
            "variable_map": str(outdir / f"{args.prefix}_baseline_csf_variable_map.tsv"),
            "group_summary": str(outdir / f"{args.prefix}_baseline_csf_group_summary.tsv"),
            "csf_match_summary": str(outdir / f"{args.prefix}_baseline_csf_csf_match_summary.tsv"),
            "associations": str(outdir / f"{args.prefix}_baseline_csf_associations.tsv")
        }
    }

    with open(outdir / f"{args.prefix}_baseline_csf_analysis_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    log("============================================================")
    log("Done.")
    log(f"Analysis dataset: {outdir / f'{args.prefix}_baseline_csf_analysis_dataset.tsv'}")
    log(f"Association table: {outdir / f'{args.prefix}_baseline_csf_associations.tsv'}")
    log("CSF match summary:")
    log(csf_match_summary.to_string(index=False))
    log("Group summary:")
    log(group_summary.to_string(index=False))
    log("============================================================")


if __name__ == "__main__":
    main()