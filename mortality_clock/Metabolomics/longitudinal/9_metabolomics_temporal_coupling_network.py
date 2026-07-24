#!/usr/bin/env python3
"""
9_metabolomics_temporal_coupling_network_revised.py

Temporal coupling network among four metabolomics mortality EPOCH clocks.

This revised script reads:
  Instance 0 / baseline from:
    <base_dir>/<Organ>_metabolomics_mortality_clock/<organ>_metabolomics_mortality_clock_predictions.tsv

  Instance 1 / follow-up from:
    <base_dir>/mortality_clock/longitudinal/metabolomics/<Organ>/<organ>_metabolomics_mortality_clock_apply_instance_1_0_predictions.tsv

This fixes the previous issue where the longitudinal instance-1 file contained only
application_instance == 1_0 and therefore could not provide instance 0 rows.

Primary simultaneous-baseline-predictor model:
  annualized_delta_EPOCH_j ~ baseline_Endocrine_EPOCH_z
                           + baseline_Digestive_EPOCH_z
                           + baseline_Hepatic_EPOCH_z
                           + baseline_Immune_EPOCH_z
                           + covariates

Default covariates included when present:
  Continuous: age_at_baseline, bmi_at_baseline, diastolic_bp_at_baseline,
              systolic_bp_at_baseline
  Categorical: sex
  followup_years is added only if raw delta is used and follow-up interval varies.

Sensitivity model:
  EPOCH_j_instance1 ~ all four baseline EPOCH clocks + covariates

Example on local Mac with CUBIC mounted:
  cd /Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock
  python /Users/hao/Project/whole-body_clocks/mortality_clock/Metabolomics/longitudinal/9_metabolomics_temporal_coupling_network_revised.py

Example on CUBIC:
  cd /gpfs/fs001/cbica/home/wenju/Reproducibile_paper/WholeBodyClock
  python 9_metabolomics_temporal_coupling_network_revised.py --base-dir /gpfs/fs001/cbica/home/wenju/Reproducibile_paper/WholeBodyClock
"""

import argparse
import json
import re
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import statsmodels.api as sm
except Exception as exc:
    raise ImportError("This script requires statsmodels: pip install statsmodels") from exc


def log(msg: str) -> None:
    print(msg, flush=True)


def safe_mkdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def read_table(path: Path) -> pd.DataFrame:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    name = path.name.lower()
    if name.endswith(".tsv.gz") or name.endswith(".txt.gz"):
        return pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)
    if name.endswith(".csv.gz"):
        return pd.read_csv(path, compression="gzip", low_memory=False)
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path, low_memory=False)
    return pd.read_csv(path, sep="\t", low_memory=False)


def clean_col(x: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(x).lower())


def first_existing_col(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    colmap = {clean_col(c): c for c in df.columns}
    for c in candidates:
        if c in df.columns:
            return c
        cc = clean_col(c)
        if cc in colmap:
            return colmap[cc]
    return None


def zscore(x: pd.Series) -> pd.Series:
    x = pd.to_numeric(x, errors="coerce")
    mu = x.mean()
    sd = x.std(ddof=1)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.nan, index=x.index)
    return (x - mu) / sd


def winsorize_series(x: pd.Series, lower_q: float, upper_q: float) -> pd.Series:
    x = pd.to_numeric(x, errors="coerce")
    if x.notna().sum() < 5:
        return x
    lo = x.quantile(lower_q)
    hi = x.quantile(upper_q)
    return x.clip(lo, hi)


def bh_fdr(pvals) -> pd.Series:
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


def p_to_stars(p: float) -> str:
    if not np.isfinite(p):
        return ""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ""


def organ_lower(organ: str) -> str:
    return organ.lower()


def find_baseline_prediction_file(base_dir: Path, organ: str) -> Path:
    ol = organ_lower(organ)
    organ_dir = base_dir / f"{organ}_metabolomics_mortality_clock"
    preferred = organ_dir / f"{ol}_metabolomics_mortality_clock_predictions.tsv"
    if preferred.exists():
        return preferred
    matches = sorted(organ_dir.glob("*metabolomics_mortality_clock_predictions.tsv"))
    if matches:
        return matches[0]
    raise FileNotFoundError(f"Could not find baseline instance-0 prediction file for {organ} in {organ_dir}")


def find_instance1_prediction_file(base_dir: Path, longitudinal_root: str, organ: str) -> Path:
    ol = organ_lower(organ)
    organ_dir = base_dir / longitudinal_root / organ
    preferred = organ_dir / f"{ol}_metabolomics_mortality_clock_apply_instance_1_0_predictions.tsv"
    if preferred.exists():
        return preferred
    matches = sorted(organ_dir.glob("*apply_instance_1_0_predictions.tsv"))
    if matches:
        return matches[0]
    matches = sorted(organ_dir.glob("*apply_longitudinal_instances_combined_predictions.tsv"))
    if matches:
        return matches[0]
    raise FileNotFoundError(f"Could not find instance-1 prediction file for {organ} in {organ_dir}")


def detect_id_col(df: pd.DataFrame, user_col: Optional[str] = None) -> str:
    if user_col and user_col in df.columns:
        return user_col
    col = first_existing_col(df, ["participant_id", "eid", "EID", "FID", "IID", "id", "ID", "sample_id", "individual_id"])
    if col:
        return col
    raise ValueError("Could not detect participant ID column. Pass --id-col.")


def detect_date_col(df: pd.DataFrame, user_col: Optional[str] = None) -> Optional[str]:
    if user_col and user_col in df.columns:
        return user_col
    col = first_existing_col(
        df,
        ["sample_date", "assessment_date", "date", "Date", "visit_date", "blood_date", "metabolomics_date",
         "baseline_date", "instance_date", "attending_assessment_centre_date", "date_of_attending_assessment_centre"],
    )
    if col:
        return col
    for c in df.columns:
        lc = str(c).lower()
        if "date" in lc and "death" not in lc and "censor" not in lc and "update" not in lc:
            return c
    return None


def detect_score_col(df: pd.DataFrame, organ: str, score_suffix: str, user_col: Optional[str] = None) -> str:
    if user_col and user_col in df.columns:
        return user_col
    ol = organ_lower(organ)
    suffix_tokens = [t for t in score_suffix.lower().split("_") if t]
    strong = [
        f"{ol}_metabolomics_mortality_clock_{score_suffix}",
        f"{ol}_metabolomics_mortality_clock_acceleration_years",
        f"{ol}_metabolomics_mortality_clock_acceleration_z",
        f"{ol}_metabolomics_mortality_risk_score",
        f"{ol}_metabolomics_mortality_clock_risk_score",
    ]
    for cand in strong:
        col = first_existing_col(df, [cand])
        if col and all(t in str(col).lower() for t in suffix_tokens):
            return col
    fuzzy = []
    for c in df.columns:
        lc = str(c).lower()
        if ol in lc and "metabolomics" in lc and "mortality" in lc:
            if all(t in lc for t in suffix_tokens) and "before" not in lc:
                fuzzy.append(c)
    if len(fuzzy) == 1:
        return fuzzy[0]
    if len(fuzzy) > 1:
        return sorted(fuzzy, key=lambda x: len(str(x)))[0]
    fallback = []
    for c in df.columns:
        lc = str(c).lower()
        if all(t in lc for t in suffix_tokens) and "before" not in lc:
            fallback.append(c)
    if len(fallback) == 1:
        return fallback[0]
    raise ValueError(f"Could not detect score column for {organ} with suffix {score_suffix}. Available columns include: {list(df.columns)[:80]}")


def sort_and_deduplicate_by_date(df: pd.DataFrame, id_col: str, date_col: Optional[str], keep: str) -> pd.DataFrame:
    d = df.copy()
    if date_col and date_col in d.columns:
        d["_date_sort"] = pd.to_datetime(d[date_col], errors="coerce")
        ascending = keep == "first"
        d = d.sort_values([id_col, "_date_sort"], ascending=[True, ascending])
        return d.drop_duplicates(id_col, keep="first").drop(columns=["_date_sort"])
    return d.drop_duplicates(id_col, keep="first")

def build_wide_for_organ(
    base_dir: Path,
    longitudinal_root: str,
    organ: str,
    score_suffix: str,
    id_col_arg: Optional[str],
    baseline_date_col_arg: Optional[str],
    instance1_date_col_arg: Optional[str],
    covariate_candidates: List[str],
) -> Tuple[pd.DataFrame, Dict]:
    baseline_path = find_baseline_prediction_file(base_dir, organ)
    instance1_path = find_instance1_prediction_file(base_dir, longitudinal_root, organ)

    t0 = read_table(baseline_path)
    t1 = read_table(instance1_path)

    id_col_t0 = detect_id_col(t0, id_col_arg)
    id_col_t1 = detect_id_col(t1, id_col_arg)

    date_col_t0 = detect_date_col(t0, baseline_date_col_arg)
    date_col_t1 = detect_date_col(t1, instance1_date_col_arg)

    score_col_t0 = detect_score_col(t0, organ, score_suffix)
    score_col_t1 = detect_score_col(t1, organ, score_suffix)

    t0 = sort_and_deduplicate_by_date(t0, id_col_t0, date_col_t0, keep="first")
    t1 = sort_and_deduplicate_by_date(t1, id_col_t1, date_col_t1, keep="first")

    keep0 = [id_col_t0, score_col_t0]
    keep1 = [id_col_t1, score_col_t1]

    if date_col_t0 and date_col_t0 in t0.columns:
        keep0.append(date_col_t0)
    if date_col_t1 and date_col_t1 in t1.columns:
        keep1.append(date_col_t1)

    covars_found = []
    for c in covariate_candidates:
        col = first_existing_col(t0, [c])
        if col and col not in keep0:
            keep0.append(col)
            covars_found.append(col)

    t0s = t0[keep0].copy().rename(columns={id_col_t0: "participant_id", score_col_t0: f"{organ}_t0"})
    t1s = t1[keep1].copy().rename(columns={id_col_t1: "participant_id", score_col_t1: f"{organ}_t1"})

    if date_col_t0 and date_col_t0 in t0s.columns:
        t0s = t0s.rename(columns={date_col_t0: f"{organ}_date_t0"})
    if date_col_t1 and date_col_t1 in t1s.columns:
        t1s = t1s.rename(columns={date_col_t1: f"{organ}_date_t1"})

    wide = t0s.merge(t1s, on="participant_id", how="inner")
    wide[f"{organ}_t0"] = pd.to_numeric(wide[f"{organ}_t0"], errors="coerce")
    wide[f"{organ}_t1"] = pd.to_numeric(wide[f"{organ}_t1"], errors="coerce")
    wide[f"{organ}_delta"] = wide[f"{organ}_t1"] - wide[f"{organ}_t0"]

    info = {
        "organ": organ,
        "baseline_instance0_file": str(baseline_path),
        "instance1_file": str(instance1_path),
        "baseline_n_rows_input": int(t0.shape[0]),
        "instance1_n_rows_input": int(t1.shape[0]),
        "id_col_t0": id_col_t0,
        "id_col_t1": id_col_t1,
        "date_col_t0": date_col_t0,
        "date_col_t1": date_col_t1,
        "score_col_t0": score_col_t0,
        "score_col_t1": score_col_t1,
        "n_paired_organ": int(wide.shape[0]),
        "covariates_found_in_baseline": covars_found,
    }
    return wide, info


def infer_followup_years(df: pd.DataFrame, organs: List[str], fallback_years: Optional[float]) -> pd.Series:
    all_years = []
    for organ in organs:
        c0 = f"{organ}_date_t0"
        c1 = f"{organ}_date_t1"
        if c0 in df.columns and c1 in df.columns:
            d0 = pd.to_datetime(df[c0], errors="coerce")
            d1 = pd.to_datetime(df[c1], errors="coerce")
            years = (d1 - d0).dt.days / 365.25
            if years.notna().sum() > 0 and (years > 0).sum() > 0:
                all_years.append(years.rename(organ))
    if all_years:
        years_df = pd.concat(all_years, axis=1)
        return years_df.mean(axis=1, skipna=True)
    if fallback_years is not None and np.isfinite(fallback_years) and fallback_years > 0:
        return pd.Series(fallback_years, index=df.index, dtype="float64")
    return pd.Series(np.nan, index=df.index, dtype="float64")


def build_analysis_dataset(args) -> Tuple[pd.DataFrame, List[Dict]]:
    base_dir = Path(args.base_dir)
    organs = [x.strip() for x in args.organs.split(",") if x.strip()]

    covariate_candidates = [
        "age_at_baseline", "sex", "bmi_at_baseline",
        "diastolic_bp_at_baseline", "systolic_bp_at_baseline"
    ]

    wide_all = None
    infos = []

    for organ in organs:
        log(f"Reading {organ}: instance 0 from baseline predictions, instance 1 from longitudinal apply file")
        wide, info = build_wide_for_organ(
            base_dir=base_dir,
            longitudinal_root=args.longitudinal_root,
            organ=organ,
            score_suffix=args.score_suffix,
            id_col_arg=args.id_col,
            baseline_date_col_arg=args.baseline_date_col,
            instance1_date_col_arg=args.instance1_date_col,
            covariate_candidates=covariate_candidates,
        )
        infos.append(info)

        if wide_all is None:
            wide_all = wide
        else:
            duplicate_cols = [c for c in wide.columns if c in wide_all.columns and c != "participant_id"]
            wide = wide.drop(columns=duplicate_cols, errors="ignore")
            wide_all = wide_all.merge(wide, on="participant_id", how="inner")

    if wide_all is None or wide_all.empty:
        raise RuntimeError("No analysis dataset could be constructed.")

    required_clock_cols = []
    for organ in organs:
        for suffix in ["t0", "t1", "delta"]:
            c = f"{organ}_{suffix}"
            wide_all[c] = pd.to_numeric(wide_all[c], errors="coerce")
            required_clock_cols.append(c)

    wide_all = wide_all.replace([np.inf, -np.inf], np.nan)
    n_before = wide_all.shape[0]
    wide_all = wide_all.dropna(subset=required_clock_cols).copy()
    n_after = wide_all.shape[0]

    if args.winsorize > 0:
        q = float(args.winsorize)
        for c in required_clock_cols:
            wide_all[c] = winsorize_series(wide_all[c], q, 1.0 - q)

    fallback_years = args.default_followup_years if args.default_followup_years > 0 else None
    wide_all["followup_years"] = infer_followup_years(wide_all, organs, fallback_years)

    for organ in organs:
        wide_all[f"{organ}_t0_z"] = zscore(wide_all[f"{organ}_t0"])
        wide_all[f"{organ}_t1_z"] = zscore(wide_all[f"{organ}_t1"])
        wide_all[f"{organ}_delta_z"] = zscore(wide_all[f"{organ}_delta"])
        if wide_all["followup_years"].notna().sum() > 0:
            wide_all[f"{organ}_delta_per_year"] = wide_all[f"{organ}_delta"] / wide_all["followup_years"]
            wide_all[f"{organ}_delta_per_year_z"] = zscore(wide_all[f"{organ}_delta_per_year"])
        else:
            wide_all[f"{organ}_delta_per_year"] = np.nan
            wide_all[f"{organ}_delta_per_year_z"] = np.nan

    infos.append({
        "analysis_dataset_n_before_complete_case": int(n_before),
        "analysis_dataset_n_after_complete_case": int(n_after),
        "required_clock_cols": required_clock_cols,
    })
    return wide_all, infos


def make_binary_sex(series: pd.Series) -> pd.Series:
    txt = series.astype(str).str.strip().str.lower()
    out = pd.Series(np.nan, index=series.index, dtype=float)
    out[txt.isin(["female", "f", "0", "0.0"])] = 1.0
    out[txt.isin(["male", "m", "1", "1.0"])] = 0.0
    x = pd.to_numeric(series, errors="coerce")
    out[(out.isna()) & (x == 2)] = 1.0
    return out


def detect_covariates(df: pd.DataFrame, args, use_raw_delta: bool) -> Tuple[pd.DataFrame, List[str], Dict[str, str]]:
    """
    Build covariate design matrix.

    Default covariates:
      Continuous:
        age_at_baseline
        bmi_at_baseline
        diastolic_bp_at_baseline
        systolic_bp_at_baseline
        followup_years, only if raw delta is primary and it varies

      Categorical:
        sex
    """
    if args.covariates.lower() == "none":
        return pd.DataFrame(index=df.index), [], {}

    if args.covariates.lower() != "auto":
        requested = [x.strip() for x in args.covariates.split(",") if x.strip()]
    else:
        requested = [
            "age_at_baseline", "sex", "bmi_at_baseline",
            "diastolic_bp_at_baseline", "systolic_bp_at_baseline"
        ]
        if use_raw_delta and "followup_years" in df.columns and df["followup_years"].nunique(dropna=True) > 1:
            requested.append("followup_years")

    continuous_like = {
        "age_at_baseline", "bmi_at_baseline",
        "diastolic_bp_at_baseline", "systolic_bp_at_baseline", "followup_years",
    }
    categorical_like = {
        "sex"
    }

    parts = []
    used = []
    treatment = {}

    for req in requested:
        col = first_existing_col(df, [req])
        if col is None:
            warnings.warn(f"Covariate not found and skipped: {req}")
            continue

        col_key = clean_col(col)

        if req in continuous_like or col in continuous_like or col_key in {clean_col(x) for x in continuous_like}:
            x = pd.to_numeric(df[col], errors="coerce")
            if x.notna().sum() >= 5 and x.nunique(dropna=True) > 1:
                x = x.fillna(x.median())
                parts.append(pd.DataFrame({f"cov__{req}": zscore(x)}, index=df.index))
                used.append(col)
                treatment[col] = "continuous_z"
            continue

        if req == "sex" or col_key == "sex":
            sex = make_binary_sex(df[col])
            if sex.notna().sum() >= 5 and sex.nunique(dropna=True) > 1:
                sex = sex.fillna(sex.mode(dropna=True).iloc[0])
                parts.append(pd.DataFrame({"cov__female": sex.astype(float)}, index=df.index))
                used.append(col)
                treatment[col] = "binary_female_1"
            continue

        if req in categorical_like or col in categorical_like or col_key in {clean_col(x) for x in categorical_like}:
            s = df[col].astype("object").where(df[col].notna(), "Missing").astype(str)
            if s.nunique(dropna=True) > 1:
                dummies = pd.get_dummies(s, prefix=f"cov__{req}", drop_first=True, dtype=float)
                dummies.index = df.index
                parts.append(dummies)
                used.append(col)
                treatment[col] = "categorical_dummy"
            continue

        x = pd.to_numeric(df[col], errors="coerce")
        if x.notna().mean() > 0.80 and x.nunique(dropna=True) > 1:
            x = x.fillna(x.median())
            parts.append(pd.DataFrame({f"cov__{req}": zscore(x)}, index=df.index))
            treatment[col] = "continuous_z_fallback"
        else:
            s = df[col].astype("object").where(df[col].notna(), "Missing").astype(str)
            if s.nunique(dropna=True) > 1:
                dummies = pd.get_dummies(s, prefix=f"cov__{req}", drop_first=True, dtype=float)
                dummies.index = df.index
                parts.append(dummies)
                treatment[col] = "categorical_dummy_fallback"
        used.append(col)

    if not parts:
        return pd.DataFrame(index=df.index), [], treatment

    Xcov = pd.concat(parts, axis=1)
    Xcov = Xcov.loc[:, Xcov.nunique(dropna=True) > 1]
    return Xcov, used, treatment

def fit_ols_hc3(y: pd.Series, X: pd.DataFrame) -> Optional[object]:
    data = pd.concat([y.rename("y"), X], axis=1).replace([np.inf, -np.inf], np.nan).dropna()
    if data.shape[0] < 20:
        return None
    y2 = data["y"].astype(float)
    X2 = data.drop(columns=["y"]).astype(float)
    X2 = sm.add_constant(X2, has_constant="add")
    return sm.OLS(y2, X2).fit(cov_type="HC3")


def run_coupling_models(df: pd.DataFrame, organs: List[str], args) -> Tuple[pd.DataFrame, pd.DataFrame, Dict]:
    has_followup = df["followup_years"].notna().sum() > 0 and (df["followup_years"] > 0).sum() > 0
    use_annualized = args.change_scale == "annualized" and has_followup
    change_suffix = "delta_per_year" if use_annualized else "delta"
    change_label = "annualized_change" if use_annualized else "raw_change"

    Xcov, used_covariates, covariate_treatment = detect_covariates(
        df=df, args=args, use_raw_delta=(not use_annualized)
    )

    baseline_cols = [f"{organ}_t0_z" for organ in organs]
    Xbase = df[baseline_cols].copy()
    Xbase.columns = [f"baseline_{organ}_z" for organ in organs]
    X = pd.concat([Xbase, Xcov], axis=1)

    rows_change = []
    rows_followup = []

    for outcome_organ in organs:
        y_change = df[f"{outcome_organ}_{change_suffix}"].copy()
        y_followup = df[f"{outcome_organ}_t1"].copy()
        fit_change = fit_ols_hc3(y_change, X)
        fit_followup = fit_ols_hc3(y_followup, X)

        for predictor_organ in organs:
            term = f"baseline_{predictor_organ}_z"
            edge_type = "within_system" if predictor_organ == outcome_organ else "cross_system"

            if fit_change is not None and term in fit_change.params.index:
                beta = float(fit_change.params[term])
                se = float(fit_change.bse[term])
                p = float(fit_change.pvalues[term])
                lo = float(beta - 1.96 * se)
                hi = float(beta + 1.96 * se)
                n = int(fit_change.nobs)
                r2 = float(fit_change.rsquared)
            else:
                beta = se = p = lo = hi = r2 = np.nan
                n = 0

            rows_change.append({
                "model": "simultaneous_baseline_predictors",
                "outcome_type": change_label,
                "outcome_system": outcome_organ,
                "predictor_system": predictor_organ,
                "edge_type": edge_type,
                "outcome": f"{outcome_organ}_{change_suffix}",
                "predictor": term,
                "beta": beta,
                "se_HC3": se,
                "ci_lower": lo,
                "ci_upper": hi,
                "p": p,
                "n": n,
                "r2": r2,
                "covariates_used": ",".join(used_covariates),
                "interpretation": f"{change_label} in {outcome_organ} EPOCH per 1-SD higher baseline {predictor_organ} EPOCH, mutually adjusted for all four baseline clocks.",
            })

            if fit_followup is not None and term in fit_followup.params.index:
                beta = float(fit_followup.params[term])
                se = float(fit_followup.bse[term])
                p = float(fit_followup.pvalues[term])
                lo = float(beta - 1.96 * se)
                hi = float(beta + 1.96 * se)
                n = int(fit_followup.nobs)
                r2 = float(fit_followup.rsquared)
            else:
                beta = se = p = lo = hi = r2 = np.nan
                n = 0

            rows_followup.append({
                "model": "lagged_followup_sensitivity",
                "outcome_type": "instance1_level",
                "outcome_system": outcome_organ,
                "predictor_system": predictor_organ,
                "edge_type": edge_type,
                "outcome": f"{outcome_organ}_t1",
                "predictor": term,
                "beta": beta,
                "se_HC3": se,
                "ci_lower": lo,
                "ci_upper": hi,
                "p": p,
                "n": n,
                "r2": r2,
                "covariates_used": ",".join(used_covariates),
                "interpretation": f"Instance-1 {outcome_organ} EPOCH per 1-SD higher baseline {predictor_organ} EPOCH, mutually adjusted for all four baseline clocks.",
            })

    change_df = pd.DataFrame(rows_change)
    followup_df = pd.DataFrame(rows_followup)

    for out in [change_df, followup_df]:
        out["q_fdr_bh"] = bh_fdr(out["p"]).values
        out["bonferroni_p_threshold"] = 0.05 / max(1, out["p"].notna().sum())
        out["fdr_significant"] = out["q_fdr_bh"] < 0.05
        out["bonferroni_significant"] = out["p"] < out["bonferroni_p_threshold"]

    meta = {
        "used_covariates": used_covariates,
        "covariate_treatment": covariate_treatment,
        "design_matrix_columns": list(X.columns),
        "primary_change_suffix": change_suffix,
        "primary_change_label": change_label,
        "change_scale_requested": args.change_scale,
        "annualized_used": bool(use_annualized),
    }
    return change_df, followup_df, meta


def matrix_from_results(res: pd.DataFrame, value_col: str, organs: List[str]) -> pd.DataFrame:
    mat = res.pivot(index="outcome_system", columns="predictor_system", values=value_col)
    return mat.reindex(index=organs, columns=organs)


def save_matrices(res: pd.DataFrame, organs: List[str], out_dir: Path, prefix: str) -> None:
    for value_col, name in [
        ("beta", "beta"),
        ("p", "p"),
        ("q_fdr_bh", "q"),
        ("fdr_significant", "fdr_significant"),
        ("bonferroni_significant", "bonferroni_significant"),
    ]:
        mat = matrix_from_results(res, value_col, organs)
        mat.to_csv(out_dir / f"{prefix}_matrix_{name}.tsv", sep="\t")


def plot_heatmap(res: pd.DataFrame, organs: List[str], out_dir: Path, prefix: str, title: str) -> None:
    beta = matrix_from_results(res, "beta", organs)
    pmat = matrix_from_results(res, "p", organs)
    qmat = matrix_from_results(res, "q_fdr_bh", organs)

    vals = beta.values.astype(float)
    vmax = np.nanmax(np.abs(vals))
    if not np.isfinite(vmax) or vmax == 0:
        vmax = 1.0

    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    im = ax.imshow(vals, cmap="RdBu_r", vmin=-vmax, vmax=vmax)

    ax.set_xticks(np.arange(len(organs)))
    ax.set_xticklabels(organs, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(organs)))
    yprefix = "Delta" if "change" in prefix else "Follow-up"
    ax.set_yticklabels([f"{yprefix} {x}" for x in organs])

    ax.set_xlabel("Baseline predictor clock")
    ax.set_ylabel("Outcome clock")
    ax.set_title(title, fontsize=12, fontweight="bold")

    for i in range(len(organs)):
        for j in range(len(organs)):
            b = beta.iloc[i, j]
            p = pmat.iloc[i, j]
            q = qmat.iloc[i, j]
            if pd.isna(b):
                label = ""
            else:
                label = f"{b:.2f}{p_to_stars(p)}"
                if np.isfinite(q) and q < 0.05:
                    label += "\nFDR"
            ax.text(j, i, label, ha="center", va="center", fontsize=8, color="black")

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Beta per 1-SD baseline EPOCH", rotation=90)
    fig.tight_layout()

    for ext in ["pdf", "png", "svg"]:
        out = out_dir / f"{prefix}_heatmap_beta.{ext}"
        if ext == "png":
            fig.savefig(out, dpi=300)
        else:
            fig.savefig(out)
    plt.close(fig)


def make_edges_table(res: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "predictor_system", "outcome_system", "edge_type", "beta", "se_HC3",
        "ci_lower", "ci_upper", "p", "q_fdr_bh", "fdr_significant",
        "bonferroni_significant", "n", "r2", "covariates_used",
    ]
    edge = res[cols].copy()
    edge = edge.sort_values(["q_fdr_bh", "p", "edge_type"], na_position="last")
    return edge.rename(columns={"predictor_system": "from_baseline_system", "outcome_system": "to_change_system"})


def write_audit(out_dir: Path, args, organs: List[str], infos: List[Dict], model_meta: Dict,
                analysis_n: int, change_df: pd.DataFrame, followup_df: pd.DataFrame) -> None:
    top = change_df.sort_values(["q_fdr_bh", "p"], na_position="last").head(20)

    audit_json = {
        "analysis": "metabolomics temporal coupling network",
        "base_dir": args.base_dir,
        "longitudinal_root": args.longitudinal_root,
        "organs": organs,
        "score_suffix": args.score_suffix,
        "analysis_n_complete_case": int(analysis_n),
        "model_meta": model_meta,
        "input_file_info": infos,
        "multiple_testing": {
            "primary_tests": int(change_df["p"].notna().sum()),
            "fdr_method": "Benjamini-Hochberg across the 4x4 primary change matrix",
            "bonferroni_threshold": float(change_df["bonferroni_p_threshold"].dropna().iloc[0]) if change_df["bonferroni_p_threshold"].notna().any() else None,
        },
        "important_assumptions": [
            "Longitudinal associations are not causal effects.",
            "Diagonal change-score associations can reflect regression to the mean.",
            "All four baseline clocks are included simultaneously.",
            "Complete-case analysis requires all four clocks at both instances.",
            "Instance-1 returner/survivor bias may be present.",
            "Baseline instance 0 and follow-up instance 1 are read from different files.",
        ],
        "top_primary_change_results": top.to_dict(orient="records"),
    }
    with open(out_dir / "metabolomics_temporal_coupling_audit.json", "w") as f:
        json.dump(audit_json, f, indent=2)

    lines = [
        "Metabolomics temporal coupling network audit",
        "===========================================",
        "",
        f"Base directory: {args.base_dir}",
        f"Longitudinal root: {args.longitudinal_root}",
        f"Organs: {', '.join(organs)}",
        f"Score suffix: {args.score_suffix}",
        f"Complete-case N: {analysis_n}",
        "",
        "Input design:",
        "  Instance 0 is read from <Organ>_metabolomics_mortality_clock/*_predictions.tsv",
        "  Instance 1 is read from mortality_clock/longitudinal/metabolomics/<Organ>/*_apply_instance_1_0_predictions.tsv",
        "",
        "Primary model:",
        "  Delta EPOCH_j ~ baseline Endocrine_z + baseline Digestive_z + baseline Hepatic_z + baseline Immune_z + covariates",
        "",
        "Sensitivity model:",
        "  Instance1 EPOCH_j ~ baseline Endocrine_z + baseline Digestive_z + baseline Hepatic_z + baseline Immune_z + covariates",
        "",
        "Covariates used:",
        ", ".join(model_meta.get("used_covariates", [])) if model_meta.get("used_covariates") else "None detected/used",
        "",
        "Covariate treatment:",
        json.dumps(model_meta.get("covariate_treatment", {}), indent=2),
        "",
        "Input files:",
    ]
    for info in infos:
        lines.append(json.dumps(info, indent=2))
    lines.extend([
        "",
        "Top primary change results:",
        top.to_string(index=False) if not top.empty else "No valid results.",
        "",
        "Assumptions and cautions:",
        "- These are longitudinal association models, not causal effect estimates.",
        "- Diagonal change-score coefficients may reflect regression to the mean.",
        "- Use the lagged-follow-up sensitivity model to assess robustness.",
        "- Cross-system edges may reflect shared systemic processes rather than direct organ-to-organ propagation.",
        "- Participants with instance 1 data may be healthier than those without repeat metabolomics.",
        "- Technical/batch covariates are included only if detected in input files.",
    ])
    (out_dir / "metabolomics_temporal_coupling_audit.txt").write_text("\n".join(lines))


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-dir", default="/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock")
    parser.add_argument("--longitudinal-root", default="mortality_clock/longitudinal/metabolomics")
    parser.add_argument("--out-dir", default=None)
    parser.add_argument("--organs", default="Endocrine,Digestive,Hepatic,Immune")
    parser.add_argument("--score-suffix", default="acceleration_years", help="e.g., acceleration_years or acceleration_z")
    parser.add_argument("--id-col", default=None)
    parser.add_argument("--baseline-date-col", default=None)
    parser.add_argument("--instance1-date-col", default=None)
    parser.add_argument("--covariates", default="auto", help="auto, none, or comma-separated covariates.")
    parser.add_argument("--change-scale", choices=["annualized", "raw"], default="annualized")
    parser.add_argument("--default-followup-years", type=float, default=-1.0)
    parser.add_argument("--winsorize", type=float, default=0.0)
    parser.add_argument("--min-n", type=int, default=50)
    parser.add_argument("--skip-plots", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    base_dir = Path(args.base_dir)
    organs = [x.strip() for x in args.organs.split(",") if x.strip()]
    if args.out_dir:
        out_dir = Path(args.out_dir)
    else:
        out_dir = base_dir / "mortality_clock/longitudinal/metabolomics/metabolomics_temporal_coupling_network"
    safe_mkdir(out_dir)

    log("Building complete-case longitudinal metabolomics EPOCH dataset...")
    df, infos = build_analysis_dataset(args)
    if df.shape[0] < args.min_n:
        raise RuntimeError(f"Complete-case N={df.shape[0]} is smaller than --min-n={args.min_n}")

    dataset_out = out_dir / "metabolomics_temporal_coupling_analysis_dataset.tsv"
    df.to_csv(dataset_out, sep="\t", index=False)
    log(f"Saved analysis dataset: {dataset_out}")
    log(f"Complete-case N: {df.shape[0]:,}")

    log("Fitting simultaneous baseline-predictor models...")
    change_df, followup_df, model_meta = run_coupling_models(df, organs, args)

    change_out = out_dir / "metabolomics_temporal_coupling_coefficients_change.tsv"
    followup_out = out_dir / "metabolomics_temporal_coupling_coefficients_followup.tsv"
    edge_out = out_dir / "metabolomics_temporal_coupling_edges_change.tsv"

    change_df.to_csv(change_out, sep="\t", index=False)
    followup_df.to_csv(followup_out, sep="\t", index=False)
    edge_df = make_edges_table(change_df)
    edge_df.to_csv(edge_out, sep="\t", index=False)

    save_matrices(change_df, organs, out_dir, "metabolomics_temporal_coupling_change")
    save_matrices(followup_df, organs, out_dir, "metabolomics_temporal_coupling_followup")

    if not args.skip_plots:
        plot_heatmap(change_df, organs, out_dir, "metabolomics_temporal_coupling_change",
                     "Temporal coupling of metabolomics mortality EPOCH change")
        plot_heatmap(followup_df, organs, out_dir, "metabolomics_temporal_coupling_followup",
                     "Lagged follow-up sensitivity model")

    write_audit(out_dir, args, organs, infos, model_meta, df.shape[0], change_df, followup_df)

    log("Done.")
    log("Primary output files:")
    for p in [
        dataset_out,
        change_out,
        followup_out,
        edge_out,
        out_dir / "metabolomics_temporal_coupling_change_matrix_beta.tsv",
        out_dir / "metabolomics_temporal_coupling_change_matrix_p.tsv",
        out_dir / "metabolomics_temporal_coupling_change_matrix_q.tsv",
        out_dir / "metabolomics_temporal_coupling_audit.txt",
    ]:
        log(f"  {p}")


if __name__ == "__main__":
    main()
