#!/usr/bin/env python3
"""
Endocrine metabolomics mortality EPOCH horizon experiment.

This script adapts the original organ_metabolomics_clock.py pipeline to answer
three linked questions while holding the cohort, predictors, preprocessing, and
train/validation/test split fixed:

Step 1. Train otherwise identical Endocrine metabolomics mortality clocks using
        5-year, 10-year, and full available mortality follow-up.
Step 2. Quantify how the resulting clocks and selected/weighted features change
        across mortality-training horizons.
Step 3. Evaluate all clocks on the same held-out test participants at common
        mortality evaluation horizons (5 and 10 years by default).

Important design choices
------------------------
* The analytic cohort is created once from baseline Endocrine metabolomics.
* The participant split is created once and reused for every mortality horizon.
* Missingness filtering and preprocessing are learned once from the common
  training split and reused for every horizon.
* Only the survival outcome definition and the horizon-specific Coxnet tuning
  are allowed to change across 5y / 10y / full models.
* Discrimination can be compared at 5y and 10y for every clock because Cox risk
  scores are time-independent rankings.
* Absolute-risk calibration/Brier metrics are only computed at times supported
  by the model's training horizon. Thus a 5-year-trained Cox model is NOT
  extrapolated to 10-year absolute risk.

Outputs are prefixed with:
  endocrine_metabolomics_mortality_horizon_clocks
unless --organ is changed.
"""

import argparse
import glob
import json
import os
import re
import warnings
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

try:
    from sksurv.linear_model import CoxnetSurvivalAnalysis, CoxPHSurvivalAnalysis
    from sksurv.metrics import (
        brier_score,
        concordance_index_censored,
        concordance_index_ipcw,
        cumulative_dynamic_auc,
    )
    from sksurv.nonparametric import kaplan_meier_estimator
    from sksurv.util import Surv
except ImportError as e:
    raise ImportError(
        "This script requires scikit-survival with IPCW/AUC/Brier metrics. Install with:\n"
        "  conda install -c conda-forge scikit-survival"
    ) from e


# -----------------------------------------------------------------------------
# Arguments and general utilities
# -----------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--organ", default="Endocrine", help="Expected use: Endocrine.")
    p.add_argument(
        "--death-xlsx",
        default="/cbica/home/wenju/Dataset/UKBB_UMelbourne/Death_related_var_from_Ye.xlsx",
    )
    p.add_argument(
        "--id-match-csv",
        default="/cbica/home/wenju/Dataset/UKBB_UMelbourne/UKB_UMelbourne_vs_Penn_match_key.csv",
    )
    p.add_argument("--organ-tsv", required=True, help="One TSV or comma-separated list/globs of TSVs.")
    p.add_argument(
        "--covariate-csv",
        default="/cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv",
    )
    p.add_argument("--admin-censor-date", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument(
        "--omics-session-id",
        "--imaging-session-id",
        dest="imaging_session_id",
        default="none",
        help="Baseline omics session_id to keep if present; default disables filtering.",
    )
    p.add_argument(
        "--feature-start-column",
        default="diagnosis",
        help="All columns after this column are organ metabolomics features.",
    )

    # Horizon experiment.
    p.add_argument(
        "--horizons",
        default="5,10,full",
        help="Comma-separated mortality training horizons in years plus optional 'full'.",
    )
    p.add_argument(
        "--evaluation-times",
        default="5,10",
        help="Common held-out mortality evaluation times in years.",
    )
    p.add_argument(
        "--split-stratify-horizon",
        default="5",
        help="Outcome horizon used only to stratify the ONE common split; default=5. Use 'full' if desired.",
    )

    # Original model settings.
    p.add_argument("--test-size", type=float, default=0.20)
    p.add_argument("--validation-size", type=float, default=0.20)
    p.add_argument("--random-state", type=int, default=2026)
    p.add_argument("--stratify-age-bins", type=int, default=5)
    p.add_argument("--max-feature-missing", type=float, default=0.20)
    p.add_argument("--l1-ratios", default="0.1,0.25,0.5,0.75,1.0")
    p.add_argument("--n-alphas", type=int, default=100)
    p.add_argument("--min-followup-days", type=int, default=1)
    p.add_argument(
        "--final-alpha-backoff-multipliers",
        default="1,2,5,10",
        help=(
            "Deterministic multipliers applied to the validation-selected alpha only if "
            "the train+validation Coxnet refit is numerically unstable. The first stable "
            "alpha is used; test-set performance is never used to choose it. The default "
            "stops at 10x; if that is still unstable, the script fails rather than silently "
            "moving far from the validation-selected model."
        ),
    )

    # Step 3 evaluation.
    p.add_argument("--n-bootstrap-comparison", type=int, default=1000)
    p.add_argument("--n-calibration-groups", type=int, default=10)
    p.add_argument("--ibs-grid-points", type=int, default=30)
    p.add_argument("--ibs-start-years", type=float, default=0.5)
    return p.parse_args()


def clean_name(x):
    x = re.sub(r"[^A-Za-z0-9]+", "_", str(x).strip().lower())
    x = re.sub(r"_+", "_", x).strip("_")
    if not x:
        raise ValueError("--organ is empty after sanitization.")
    return x


def output_prefix(organ):
    return f"{organ}_metabolomics_mortality_horizon_clocks"


def make_onehot_encoder():
    try:
        return OneHotEncoder(handle_unknown="ignore", sparse_output=False)
    except TypeError:
        return OneHotEncoder(handle_unknown="ignore", sparse=False)


def parse_horizons(text):
    out = []
    seen = set()
    for piece in [x.strip() for x in str(text).split(",") if x.strip()]:
        low = piece.lower()
        if low in {"full", "max", "maximum", "all"}:
            label, years = "full", None
        else:
            years = float(piece)
            if years <= 0:
                raise ValueError(f"Horizon must be >0 years: {piece}")
            label = f"{years:g}y"
        if label not in seen:
            out.append((label, years))
            seen.add(label)
    if not out:
        raise ValueError("No valid --horizons supplied.")
    return out


def parse_float_list(text, name):
    vals = [float(x.strip()) for x in str(text).split(",") if x.strip()]
    if not vals or any(v <= 0 for v in vals):
        raise ValueError(f"{name} must contain positive numbers.")
    return sorted(set(vals))


def json_safe(x):
    if isinstance(x, dict):
        return {str(k): json_safe(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [json_safe(v) for v in x]
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating,)):
        return None if not np.isfinite(x) else float(x)
    if isinstance(x, np.ndarray):
        return [json_safe(v) for v in x.tolist()]
    if isinstance(x, (pd.Timestamp,)):
        return str(x)
    return x


# -----------------------------------------------------------------------------
# Data loading: retained from the original mortality-clock pipeline
# -----------------------------------------------------------------------------

def load_death_data(death_xlsx, id_match_csv):
    d = pd.read_excel(death_xlsx)
    m = pd.read_csv(id_match_csv)
    d = d.rename(columns={"eid": "participant_id_umel"})
    m = m.rename(columns={"id": "participant_id_umel", "id_upenn": "participant_id"})
    d = m.merge(d, on="participant_id_umel", how="inner")
    need = ["participant_id", "53-0.0", "40000-0.0"]
    miss = [c for c in need if c not in d.columns]
    if miss:
        raise ValueError(f"Death file is missing columns: {miss}")
    d = d[need].copy()
    d["baseline_date"] = pd.to_datetime(d["53-0.0"], errors="coerce")
    d["sample_date"] = d["baseline_date"]
    d["death_date"] = pd.to_datetime(d["40000-0.0"], errors="coerce")
    return d


def expand_paths(arg):
    pieces = [x.strip() for x in str(arg).split(",") if x.strip()]
    if not pieces:
        raise ValueError("--organ-tsv cannot be empty.")
    paths = []
    for piece in pieces:
        expanded = sorted(glob.glob(piece)) if any(ch in piece for ch in ["*", "?", "["]) else [piece]
        if not expanded:
            raise FileNotFoundError(f"No files matched: {piece}")
        paths.extend(expanded)
    return list(dict.fromkeys(paths))


def load_organ_data(organ_tsv, organ, imaging_session_id):
    frames = []
    for i, path in enumerate(expand_paths(organ_tsv)):
        if not os.path.exists(path):
            raise FileNotFoundError(f"Organ TSV file not found: {path}")
        part = pd.read_csv(path, sep="\t")
        if "participant_id" not in part.columns:
            raise ValueError(f"participant_id is missing from {path}")
        part["organ_source_file"] = Path(path).name
        part["organ_source_order"] = i
        part["organ_source_row"] = np.arange(part.shape[0])
        frames.append(part)
        print(f"Loaded {organ} TSV: {path}; rows={part.shape[0]}, cols={part.shape[1]}")

    df = pd.concat(frames, axis=0, ignore_index=True, sort=False)
    print(f"Concatenated {organ} TSVs: rows={df.shape[0]}, cols={df.shape[1]}")

    if imaging_session_id and str(imaging_session_id).lower() not in {"none", "null", ""}:
        if "session_id" in df.columns:
            before = df.shape[0]
            df = df.loc[df["session_id"].astype(str) == str(imaging_session_id)].copy()
            print(f"Filtered {organ} to session_id={imaging_session_id}: {before} -> {df.shape[0]} rows")
        else:
            warnings.warn("--omics-session-id was provided, but session_id is missing.")

    df = df.sort_values(
        ["participant_id", "organ_source_order", "organ_source_row"], kind="mergesort"
    )
    dup = int(df["participant_id"].duplicated().sum())
    if dup > 0:
        warnings.warn(f"Found {dup} duplicated participant_id rows. Keeping first row by input order.")
        df = df.drop_duplicates("participant_id", keep="first")
    return df


def infer_feature_columns(df_organ, feature_start_column):
    if feature_start_column not in df_organ.columns:
        raise ValueError(
            f"Feature start column '{feature_start_column}' not found. "
            f"First columns: {list(df_organ.columns[:20])}"
        )
    start = list(df_organ.columns).index(feature_start_column) + 1
    excluded = {
        "organ_source_file",
        "organ_source_order",
        "organ_source_row",
        "participant_id",
        "session_id",
        feature_start_column,
    }
    features = [c for c in list(df_organ.columns[start:]) if c not in excluded]
    if not features:
        raise ValueError(f"No metabolomics features found after {feature_start_column}.")
    print(f"Feature rule: all columns after '{feature_start_column}'. N={len(features)}")
    print(f"First feature: {features[0]}")
    print(f"Last feature:  {features[-1]}")
    return features


def load_covariates(path):
    if path is None or str(path).lower() in {"none", ""}:
        return None
    if not os.path.exists(path):
        warnings.warn(f"Covariate file not found: {path}. Continuing without it.")
        return None
    cov = pd.read_csv(path)
    if "eid" not in cov.columns:
        warnings.warn("Covariate file does not contain eid. Continuing without it.")
        return None
    return cov.rename(columns={"eid": "participant_id"})


def mean_existing_numeric_columns(df, cols, out_col):
    present = [c for c in cols if c in df.columns]
    if not present:
        return df
    for c in present:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df[out_col] = df[present].mean(axis=1, skipna=True)
    return df


def add_basic_covariates(df):
    df = df.copy()
    age_col = "age_when_attended_assessment_centre_f21003_0_0"
    if age_col in df.columns:
        df["age_at_baseline"] = pd.to_numeric(df[age_col], errors="coerce")
    elif "age_at_baseline" in df.columns:
        df["age_at_baseline"] = pd.to_numeric(df["age_at_baseline"], errors="coerce")
    elif "diagnosis" in df.columns:
        warnings.warn("Using numeric diagnosis as age_at_baseline fallback. Please verify this is correct.")
        df["age_at_baseline"] = pd.to_numeric(df["diagnosis"], errors="coerce")
    else:
        raise ValueError("Could not infer age_at_baseline.")

    if "sex_f31_0_0" in df.columns:
        df["sex"] = df["sex_f31_0_0"].astype(str).str.strip()
    elif "sex" in df.columns:
        df["sex"] = df["sex"].astype(str).str.strip()
    elif "Sex" in df.columns:
        df["sex"] = df["Sex"].astype(str).str.strip()
    else:
        raise ValueError("Could not infer sex. Expected sex_f31_0_0 in covariates.")

    df["sex"] = df["sex"].replace(
        {
            "0": "Female",
            "0.0": "Female",
            "1": "Male",
            "1.0": "Male",
            "F": "Female",
            "M": "Male",
            "female": "Female",
            "male": "Male",
        }
    )

    if "body_mass_index_bmi_f23104_0_0" in df.columns:
        df["bmi_at_baseline"] = pd.to_numeric(
            df["body_mass_index_bmi_f23104_0_0"], errors="coerce"
        )
    df = mean_existing_numeric_columns(
        df,
        [
            "diastolic_blood_pressure_automated_reading_f4079_0_0",
            "diastolic_blood_pressure_automated_reading_f4079_0_1",
        ],
        "diastolic_bp_at_baseline",
    )
    df = mean_existing_numeric_columns(
        df,
        [
            "systolic_blood_pressure_automated_reading_f4080_0_0",
            "systolic_blood_pressure_automated_reading_f4080_0_1",
        ],
        "systolic_bp_at_baseline",
    )
    if "smoking_status_f20116_0_0" in df.columns:
        df["smoking_status_at_baseline"] = df["smoking_status_f20116_0_0"].astype("category")
    if "uk_biobank_assessment_centre_f54_0_0" in df.columns:
        df["uk_biobank_assessment_centre_f54_0_0"] = df[
            "uk_biobank_assessment_centre_f54_0_0"
        ].astype("category")

    # Backward-compatible alias retained from the original script.
    df["age_at_imaging"] = df["age_at_baseline"]
    return df


def build_design_matrix(df, organ_feature_cols):
    numeric_covariates = ["age_at_baseline"]
    for c in ["bmi_at_baseline", "diastolic_bp_at_baseline", "systolic_bp_at_baseline"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
            numeric_covariates.append(c)

    categorical_covariates = ["sex"]
    if "smoking_status_at_baseline" in df.columns:
        df["smoking_status_at_baseline"] = df["smoking_status_at_baseline"].astype("category")
        categorical_covariates.append("smoking_status_at_baseline")
    if "uk_biobank_assessment_centre_f54_0_0" in df.columns:
        df["uk_biobank_assessment_centre_f54_0_0"] = df[
            "uk_biobank_assessment_centre_f54_0_0"
        ].astype("category")
        categorical_covariates.append("uk_biobank_assessment_centre_f54_0_0")

    for c in organ_feature_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    return (
        df,
        numeric_covariates + organ_feature_cols,
        categorical_covariates,
        numeric_covariates,
        organ_feature_cols,
    )


# -----------------------------------------------------------------------------
# Step 1: one common cohort, then horizon-specific survival outcomes
# -----------------------------------------------------------------------------

def construct_full_survival_dataset(df):
    """Construct prospective baseline-to-death/censoring outcome once."""
    df = df.copy()
    df["death_before_or_on_sample"] = (
        df["death_date"].notna()
        & df["sample_date"].notna()
        & (df["death_date"] <= df["sample_date"])
    )
    n_pre = int(df["death_before_or_on_sample"].sum())
    if n_pre:
        warnings.warn(f"Excluding {n_pre} participants with death before/on baseline sample date.")

    df = df.loc[df["sample_date"].notna()].copy()
    df = df.loc[~df["death_before_or_on_sample"]].copy()
    df = df.loc[df["sample_date"] <= df["admin_censor_date"]].copy()

    df["event_full"] = (
        df["death_date"].notna()
        & (df["death_date"] > df["sample_date"])
        & (df["death_date"] <= df["admin_censor_date"])
    )
    df["end_date_full"] = df["admin_censor_date"]
    df.loc[df["event_full"], "end_date_full"] = df.loc[df["event_full"], "death_date"]
    df["time_days_full"] = (df["end_date_full"] - df["sample_date"]).dt.days
    df["time_years_full"] = df["time_days_full"] / 365.25
    return df


def calendar_horizon_date(series, years):
    """Calendar-year horizon date; optimized for the intended 5y/10y integer horizons."""
    if float(years).is_integer():
        off = pd.DateOffset(years=int(round(years)))
        return series.apply(lambda x: x + off if pd.notna(x) else pd.NaT)
    days = int(round(float(years) * 365.25))
    return series + pd.to_timedelta(days, unit="D")


def add_horizon_outcome(df, label, years):
    """Add event_<label>, time_years_<label>, and censor/end-date columns."""
    df = df.copy()
    if years is None:
        # Full follow-up columns already exist.
        return df

    horizon_date = calendar_horizon_date(df["sample_date"], years)
    horizon_censor_date = pd.concat(
        [horizon_date.rename("horizon"), df["admin_censor_date"].rename("admin")], axis=1
    ).min(axis=1)

    event_col = f"event_{label}"
    end_col = f"end_date_{label}"
    days_col = f"time_days_{label}"
    time_col = f"time_years_{label}"
    censor_col = f"censor_date_{label}"

    df[censor_col] = horizon_censor_date
    df[event_col] = (
        df["death_date"].notna()
        & (df["death_date"] > df["sample_date"])
        & (df["death_date"] <= horizon_censor_date)
    )
    df[end_col] = horizon_censor_date
    df.loc[df[event_col], end_col] = df.loc[df[event_col], "death_date"]
    df[days_col] = (df[end_col] - df["sample_date"]).dt.days

    # Express time in years, but cap at the nominal horizon. This avoids small
    # leap-year differences making a nominal 5y/10y endpoint slightly >5 or >10.
    raw_years = df[days_col] / 365.25
    df[time_col] = np.minimum(raw_years.astype(float), float(years))
    return df


def outcome_columns(label):
    return f"event_{label}", f"time_years_{label}"


def survival_from_df(df, label):
    event_col, time_col = outcome_columns(label)
    return Surv.from_arrays(
        event=df[event_col].astype(bool).values,
        time=df[time_col].astype(float).values,
    )


def truncate_survival(y, tau):
    event = np.asarray(y["event"]).astype(bool)
    time = np.asarray(y["time"]).astype(float)
    event_tau = event & (time <= float(tau))
    time_tau = np.minimum(time, float(tau))
    return Surv.from_arrays(event=event_tau, time=time_tau)


def make_stratify_vector(df, event_col, age_bins=5):
    if event_col not in df.columns:
        raise ValueError(f"Stratification event column not found: {event_col}")
    counts = df[event_col].value_counts()
    if len(counts) < 2 or counts.min() < 2:
        return None
    if age_bins and age_bins > 1 and "age_at_baseline" in df.columns:
        try:
            age_bin = pd.qcut(
                df["age_at_baseline"].rank(method="first"),
                q=age_bins,
                labels=False,
                duplicates="drop",
            )
            lab = df[event_col].astype(int).astype(str) + "_age" + age_bin.astype(str)
            if lab.value_counts().min() >= 2:
                return lab
        except Exception:
            pass
    return df[event_col]


def make_preprocessor(numeric_cols, categorical_cols):
    num = Pipeline(
        [("imputer", SimpleImputer(strategy="median")), ("scaler", StandardScaler())]
    )
    cat = Pipeline(
        [
            ("imputer", SimpleImputer(strategy="most_frequent")),
            ("onehot", make_onehot_encoder()),
        ]
    )
    return ColumnTransformer(
        [("num", num, numeric_cols), ("cat", cat, categorical_cols)], remainder="drop"
    )


def get_feature_names(preprocessor):
    names = [f"num__{c}" for c in preprocessor.named_transformers_["num"].feature_names_in_]
    cat_pipe = preprocessor.named_transformers_["cat"]
    ohe = cat_pipe.named_steps["onehot"]
    cat_features = preprocessor.transformers_[1][2]
    try:
        cat_names = ohe.get_feature_names_out(cat_features)
    except AttributeError:
        cat_names = ohe.get_feature_names(cat_features)
    names.extend([f"cat__{c}" for c in cat_names])
    return np.array(names)


def drop_high_missing_features(
    X_train_raw, other_raw_list, numeric_cols, categorical_cols, max_missing
):
    keep_numeric, dropped = [], []
    for c in numeric_cols:
        miss = X_train_raw[c].isna().mean()
        if miss <= max_missing:
            keep_numeric.append(c)
        else:
            dropped.append((c, float(miss)))
    if dropped:
        print(f"Dropped {len(dropped)} numeric columns with missingness > {max_missing}.")
        for c, miss in dropped[:20]:
            print(f"  dropped: {c}, missing={miss:.3f}")
        if len(dropped) > 20:
            print("  ...")
    cols = keep_numeric + categorical_cols
    return (
        X_train_raw[cols].copy(),
        [x[cols].copy() for x in other_raw_list],
        keep_numeric,
        categorical_cols,
        dropped,
    )


def organ_feature_mask(feature_names, organ_feature_cols):
    organ_set = {f"num__{c}" for c in organ_feature_cols}
    return np.array([str(f) in organ_set for f in feature_names], dtype=bool)


# -----------------------------------------------------------------------------
# Coxnet model fitting: same model class/tuning logic as original script
# -----------------------------------------------------------------------------

def compute_harrell_c(y, risk):
    return float(
        concordance_index_censored(
            y["event"], y["time"], np.asarray(risk).reshape(-1)
        )[0]
    )


def fit_and_select_coxnet(
    X_train,
    y_train,
    X_val,
    y_val,
    feature_names,
    organ_feature_cols,
    organ,
    l1_ratios,
    n_alphas,
):
    penalty_factor = np.ones(len(feature_names), dtype=float)
    is_organ = organ_feature_mask(feature_names, organ_feature_cols)
    penalty_factor[~is_organ] = 0.0

    best = {
        "cindex": -np.inf,
        "l1_ratio": None,
        "alpha": None,
        "coef": None,
        "used_penalty_factor": True,
    }

    for l1_ratio in l1_ratios:
        print(f"Fitting Coxnet path for l1_ratio={l1_ratio}")
        try:
            model = CoxnetSurvivalAnalysis(
                l1_ratio=l1_ratio,
                n_alphas=n_alphas,
                alpha_min_ratio="auto",
                penalty_factor=penalty_factor,
                fit_baseline_model=False,
                max_iter=100000,
            )
        except TypeError:
            warnings.warn(
                "Installed scikit-survival does not support penalty_factor. Covariates will be penalized."
            )
            model = CoxnetSurvivalAnalysis(
                l1_ratio=l1_ratio,
                n_alphas=n_alphas,
                alpha_min_ratio="auto",
                fit_baseline_model=False,
                max_iter=100000,
            )
            best["used_penalty_factor"] = False

        model.fit(X_train, y_train)
        coefs = model.coef_
        if coefs.ndim == 1:
            coefs = coefs[:, None]

        for j, alpha in enumerate(model.alphas_):
            risk_val = np.dot(X_val, coefs[:, j])
            cindex = compute_harrell_c(y_val, risk_val)
            if np.isfinite(cindex) and cindex > best["cindex"]:
                best.update(
                    {
                        "cindex": float(cindex),
                        "l1_ratio": float(l1_ratio),
                        "alpha": float(alpha),
                        "coef": coefs[:, j].copy(),
                        "alpha_index": int(j),
                        "n_alphas_fitted": int(len(model.alphas_)),
                        "alpha_path_max": float(np.max(model.alphas_)),
                        "alpha_path_min": float(np.min(model.alphas_)),
                        "selected_at_smallest_alpha": bool(
                            np.isclose(float(alpha), float(np.min(model.alphas_)))
                        ),
                    }
                )
        print(
            f"  best so far: C-index={best['cindex']:.4f}, "
            f"l1_ratio={best['l1_ratio']}, alpha={best['alpha']}"
        )

    if best["alpha"] is None:
        raise RuntimeError("Failed to select a Coxnet model.")
    return best, penalty_factor


def fit_final_model(
    X_trainval,
    y_trainval,
    best,
    penalty_factor,
    alpha_backoff_multipliers,
):
    """Refit the validation-selected Coxnet model on train+validation robustly.

    Coxnet can become numerically unstable when a very small alpha that was
    reachable along a warm-start regularization path is refit as a single-alpha
    model on the larger train+validation set. We therefore:

    1. Preserve the validation-selected l1_ratio and target alpha.
    2. Refit using a short descending alpha path ending at the target alpha,
       which supplies warm starts similar to the original tuning fit.
    3. If scikit-survival raises the specific numerical ArithmeticError, increase
       alpha deterministically and retry. The first numerically stable alpha is
       used. No validation or test performance is consulted during this backoff.

    The actual alpha used is returned and saved in the output metadata.
    """
    selected_alpha = float(best["alpha"])
    l1_ratio = float(best["l1_ratio"])
    warm_start_multipliers = (100.0, 30.0, 10.0, 3.0, 1.0)
    attempts = []

    multipliers = sorted({float(x) for x in alpha_backoff_multipliers if float(x) >= 1.0})
    if not multipliers or multipliers[0] != 1.0:
        multipliers = [1.0] + multipliers

    last_exc = None
    for backoff in multipliers:
        target_alpha = selected_alpha * backoff
        alpha_path = [target_alpha * m for m in warm_start_multipliers]
        # Ensure strict descending order and remove accidental duplicates.
        alpha_path = list(dict.fromkeys(alpha_path))
        alpha_path = sorted(alpha_path, reverse=True)

        print(
            f"  Final refit attempt: selected alpha={selected_alpha:.6g}, "
            f"multiplier={backoff:g}, target alpha={target_alpha:.6g}"
        )
        try:
            try:
                model = CoxnetSurvivalAnalysis(
                    l1_ratio=l1_ratio,
                    alphas=alpha_path,
                    penalty_factor=penalty_factor,
                    fit_baseline_model=True,
                    max_iter=100000,
                )
            except TypeError:
                warnings.warn(
                    "Installed scikit-survival does not support penalty_factor in the "
                    "final refit. Covariates will be penalized."
                )
                model = CoxnetSurvivalAnalysis(
                    l1_ratio=l1_ratio,
                    alphas=alpha_path,
                    fit_baseline_model=True,
                    max_iter=100000,
                )

            model.fit(X_trainval, y_trainval)

            actual_alpha = float(np.asarray(model.alphas_).reshape(-1)[-1])
            risk = np.asarray(model.predict(X_trainval, alpha=actual_alpha)).reshape(-1)
            if not np.all(np.isfinite(risk)):
                raise FloatingPointError(
                    "Final Coxnet refit produced non-finite train+validation risk scores."
                )

            coef_arr = np.asarray(model.coef_)
            coef_last = coef_arr if coef_arr.ndim == 1 else coef_arr[:, -1]
            if not np.all(np.isfinite(coef_last)):
                raise FloatingPointError(
                    "Final Coxnet refit produced non-finite coefficients."
                )

            attempts.append(
                {
                    "alpha_multiplier": float(backoff),
                    "target_alpha": float(target_alpha),
                    "status": "success",
                    "error": None,
                }
            )
            info = {
                "selected_alpha_from_validation": selected_alpha,
                "final_alpha_used": actual_alpha,
                "final_alpha_multiplier": float(actual_alpha / selected_alpha),
                "alpha_was_increased_for_numerical_stability": bool(
                    actual_alpha > selected_alpha * (1.0 + 1e-10)
                ),
                "warm_start_alpha_multipliers": list(warm_start_multipliers),
                "max_abs_trainval_linear_predictor": float(np.max(np.abs(risk))),
                "attempts": attempts,
            }
            if info["alpha_was_increased_for_numerical_stability"]:
                warnings.warn(
                    "Final train+validation refit required a larger alpha for numerical "
                    f"stability: {selected_alpha:.6g} -> {actual_alpha:.6g} "
                    f"({info['final_alpha_multiplier']:.3g}x). This adjustment was based "
                    "only on numerical convergence, not on validation/test performance."
                )
            return model, info

        except (ArithmeticError, FloatingPointError) as exc:
            last_exc = exc
            attempts.append(
                {
                    "alpha_multiplier": float(backoff),
                    "target_alpha": float(target_alpha),
                    "status": "failed_numerically",
                    "error": str(exc),
                }
            )
            warnings.warn(
                f"Final Coxnet refit failed numerically at alpha={target_alpha:.6g}: {exc}"
            )
            continue

    raise RuntimeError(
        "Final Coxnet refit remained numerically unstable after all alpha backoff "
        f"attempts. Selected alpha={selected_alpha:.6g}. Last error: {last_exc}"
    )


def predict_risk_score(model, X):
    return np.asarray(model.predict(X)).reshape(-1)


# -----------------------------------------------------------------------------
# Clock transforms retained from the original EPOCH pipeline
# -----------------------------------------------------------------------------

def add_clock_age_and_acceleration(
    pred_df, organ, horizon_label, risk_col, covariate_cols
):
    df = pred_df.copy()
    z_col = f"{organ}_metabolomics_mortality_clock_acceleration_z_{horizon_label}"
    yrs_col = f"{organ}_metabolomics_mortality_clock_acceleration_years_{horizon_label}"
    age_col = f"{organ}_metabolomics_mortality_clock_age_years_{horizon_label}"

    covariate_cols = [c for c in covariate_cols if c in df.columns]
    if "age_at_baseline" not in covariate_cols and "age_at_baseline" in df.columns:
        covariate_cols = ["age_at_baseline"] + covariate_cols

    train = df.loc[df["split"] == "train"].copy()
    if train.shape[0] < 10:
        df[z_col] = df[yrs_col] = df[age_col] = np.nan
        return df, None

    numeric_covs, categorical_covs = [], []
    for c in covariate_cols:
        if c == "sex" or not pd.api.types.is_numeric_dtype(df[c]):
            categorical_covs.append(c)
        else:
            numeric_covs.append(c)

    transformers = []
    if numeric_covs:
        transformers.append(
            ("num", Pipeline([("imputer", SimpleImputer(strategy="median"))]), numeric_covs)
        )
    if categorical_covs:
        transformers.append(
            (
                "cat",
                Pipeline(
                    [
                        ("imputer", SimpleImputer(strategy="most_frequent")),
                        ("onehot", make_onehot_encoder()),
                    ]
                ),
                categorical_covs,
            )
        )

    prep = ColumnTransformer(transformers=transformers, remainder="drop")
    X_train_raw = train[covariate_cols].copy()
    X_all_raw = df[covariate_cols].copy()

    for c in numeric_covs:
        X_train_raw[c] = pd.to_numeric(X_train_raw[c], errors="coerce")
        X_all_raw[c] = pd.to_numeric(X_all_raw[c], errors="coerce")
    for c in categorical_covs:
        X_train_raw[c] = X_train_raw[c].astype("object")
        X_all_raw[c] = X_all_raw[c].astype("object")

    Xtr = prep.fit_transform(X_train_raw)
    Xall = prep.transform(X_all_raw)
    lr = LinearRegression().fit(Xtr, train[risk_col].values)
    expected = lr.predict(Xall)
    resid_raw = df[risk_col].values - expected

    train_index = df["split"].values == "train"
    mean_train = float(np.nanmean(resid_raw[train_index]))
    sd_train = float(np.nanstd(resid_raw[train_index]))
    resid = resid_raw - mean_train
    df[z_col] = resid / sd_train if sd_train > 0 else np.nan

    feat_names = []
    feat_names.extend([f"num__{c}" for c in numeric_covs])
    if categorical_covs:
        ohe = prep.named_transformers_["cat"].named_steps["onehot"]
        try:
            feat_names.extend([f"cat__{c}" for c in ohe.get_feature_names_out(categorical_covs)])
        except AttributeError:
            feat_names.extend([f"cat__{c}" for c in ohe.get_feature_names(categorical_covs)])

    beta_age = np.nan
    if "num__age_at_baseline" in feat_names:
        beta_age = float(lr.coef_[feat_names.index("num__age_at_baseline")])

    if np.isfinite(beta_age) and abs(beta_age) > 1e-8:
        df[yrs_col] = resid / beta_age
        df[age_col] = df["age_at_baseline"] + df[yrs_col]
    else:
        warnings.warn(
            f"{horizon_label}: adjusted age coefficient is near zero/unavailable; "
            "year-scale acceleration set to missing."
        )
        df[yrs_col] = np.nan
        df[age_col] = np.nan

    info = {
        "horizon_label": horizon_label,
        "risk_col": risk_col,
        "z_col": z_col,
        "years_col": yrs_col,
        "clock_age_col": age_col,
        "residualization_covariates": covariate_cols,
        "numeric_residualization_covariates": numeric_covs,
        "categorical_residualization_covariates": categorical_covs,
        "risk_score_covariate_model_intercept": float(lr.intercept_),
        "risk_score_covariate_model_coef": {
            k: float(v) for k, v in zip(feat_names, lr.coef_)
        },
        "adjusted_age_coefficient_risk_score_per_year": (
            float(beta_age) if np.isfinite(beta_age) else None
        ),
        "risk_score_residual_mean_train": mean_train,
        "risk_score_residual_sd_train": sd_train,
        "note": (
            "Clock acceleration is the residual of the Cox risk score after adjustment "
            "for retained non-organ covariates; transform is learned on the shared train split."
        ),
    }
    return df, info


# -----------------------------------------------------------------------------
# Step 2: cross-horizon clock and coefficient comparisons
# -----------------------------------------------------------------------------

def pairwise_labels(labels):
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            yield labels[i], labels[j]


def safe_corr(a, b, method="pearson"):
    a = pd.Series(np.asarray(a, dtype=float))
    b = pd.Series(np.asarray(b, dtype=float))
    ok = a.notna() & b.notna()
    if int(ok.sum()) < 3:
        return np.nan
    return float(a[ok].corr(b[ok], method=method))


def build_score_correlation_tables(pred_df, organ, horizon_labels):
    rows = []
    for score_type in ["risk_score", "acceleration_z"]:
        if score_type == "risk_score":
            cols = {
                h: f"{organ}_metabolomics_mortality_risk_score_{h}" for h in horizon_labels
            }
        else:
            cols = {
                h: f"{organ}_metabolomics_mortality_clock_acceleration_z_{h}"
                for h in horizon_labels
            }

        for subset_name, sub in [
            ("test", pred_df.loc[pred_df["split"] == "test"]),
            ("all", pred_df),
        ]:
            for a, b in pairwise_labels(horizon_labels):
                rows.append(
                    {
                        "score_type": score_type,
                        "subset": subset_name,
                        "horizon_a": a,
                        "horizon_b": b,
                        "n_pair": int(sub[[cols[a], cols[b]]].dropna().shape[0]),
                        "pearson_r": safe_corr(sub[cols[a]], sub[cols[b]], "pearson"),
                        "spearman_rho": safe_corr(sub[cols[a]], sub[cols[b]], "spearman"),
                    }
                )
    return pd.DataFrame(rows)


def coefficient_comparison(coefficients, feature_names, organ_feature_cols, horizon_labels):
    rows = []
    is_organ = organ_feature_mask(feature_names, organ_feature_cols)
    masks = {
        "all_features": np.ones(len(feature_names), dtype=bool),
        "endocrine_metabolomics_only": is_organ,
    }

    for a, b in pairwise_labels(horizon_labels):
        ca = np.asarray(coefficients[a], dtype=float)
        cb = np.asarray(coefficients[b], dtype=float)
        for scope, mask in masks.items():
            xa, xb = ca[mask], cb[mask]
            nz_a, nz_b = xa != 0, xb != 0
            union = int(np.sum(nz_a | nz_b))
            inter = int(np.sum(nz_a & nz_b))
            rows.append(
                {
                    "horizon_a": a,
                    "horizon_b": b,
                    "scope": scope,
                    "n_features": int(np.sum(mask)),
                    "n_nonzero_a": int(np.sum(nz_a)),
                    "n_nonzero_b": int(np.sum(nz_b)),
                    "n_nonzero_intersection": inter,
                    "n_nonzero_union": union,
                    "nonzero_jaccard": float(inter / union) if union > 0 else np.nan,
                    "coefficient_pearson_r": safe_corr(xa, xb, "pearson"),
                    "coefficient_spearman_rho": safe_corr(xa, xb, "spearman"),
                }
            )
    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# Step 3: held-out mortality evaluation at common horizons
# -----------------------------------------------------------------------------

def safe_uno_c(y_train_full, y_test_full, risk, tau):
    try:
        return float(
            concordance_index_ipcw(
                y_train_full,
                y_test_full,
                np.asarray(risk).reshape(-1),
                tau=float(tau),
            )[0]
        )
    except Exception as exc:
        warnings.warn(f"Uno C-index failed at tau={tau}: {exc}")
        return np.nan


def safe_td_auc(y_train_full, y_test_full, risk, tau):
    try:
        auc, _ = cumulative_dynamic_auc(
            y_train_full,
            y_test_full,
            np.asarray(risk).reshape(-1),
            np.asarray([float(tau)]),
        )
        return float(np.asarray(auc).reshape(-1)[0])
    except Exception as exc:
        warnings.warn(f"Time-dependent AUC failed at t={tau}: {exc}")
        return np.nan


def predict_survival_probability(model, X, tau):
    try:
        surv_funcs = model.predict_survival_function(X)
    except Exception as exc:
        warnings.warn(f"Could not generate survival functions: {exc}")
        return np.full(X.shape[0], np.nan, dtype=float)

    out = np.full(X.shape[0], np.nan, dtype=float)
    for i, sf in enumerate(surv_funcs):
        try:
            out[i] = float(sf(float(tau)))
        except Exception:
            out[i] = np.nan
    return out


def safe_brier(y_train_full, y_test_full, surv_prob, tau):
    ok = np.isfinite(surv_prob)
    if int(np.sum(ok)) != len(surv_prob):
        return np.nan
    try:
        _, scores = brier_score(
            y_train_full,
            y_test_full,
            np.asarray(surv_prob, dtype=float).reshape(-1, 1),
            np.asarray([float(tau)]),
        )
        return float(np.asarray(scores).reshape(-1)[0])
    except Exception as exc:
        warnings.warn(f"Brier score failed at t={tau}: {exc}")
        return np.nan


def safe_ibs(
    y_train_full,
    y_test_full,
    model,
    X_test,
    tau,
    start_years=0.5,
    grid_points=30,
):
    tau = float(tau)
    start = float(start_years)
    if tau <= start:
        return np.nan
    times = np.linspace(start, tau, int(max(3, grid_points)))

    try:
        surv_funcs = model.predict_survival_function(X_test)
    except Exception as exc:
        warnings.warn(f"IBS survival prediction failed: {exc}")
        return np.nan

    probs = np.empty((X_test.shape[0], len(times)), dtype=float)
    probs[:] = np.nan
    for i, sf in enumerate(surv_funcs):
        for j, t in enumerate(times):
            try:
                probs[i, j] = float(sf(float(t)))
            except Exception:
                return np.nan

    if not np.isfinite(probs).all():
        return np.nan

    try:
        _, scores = brier_score(y_train_full, y_test_full, probs, times)
        scores = np.asarray(scores, dtype=float)
        return float(np.trapz(scores, times) / (times[-1] - times[0]))
    except Exception as exc:
        warnings.warn(f"IBS calculation failed through t={tau}: {exc}")
        return np.nan


def km_observed_risk(y, tau):
    event = np.asarray(y["event"]).astype(bool)
    time = np.asarray(y["time"]).astype(float)
    if len(time) == 0:
        return np.nan
    try:
        t, s = kaplan_meier_estimator(event, time)
        idx = np.where(t <= float(tau))[0]
        surv = 1.0 if len(idx) == 0 else float(s[idx[-1]])
        return 1.0 - surv
    except Exception as exc:
        warnings.warn(f"Kaplan-Meier risk failed at t={tau}: {exc}")
        return np.nan


def safe_calibration_slope(y_test_full, risk, tau):
    """Cox calibration slope of the original linear predictor on test data, truncated at tau."""
    y_tau = truncate_survival(y_test_full, tau)
    x = np.asarray(risk, dtype=float).reshape(-1, 1)
    try:
        try:
            model = CoxPHSurvivalAnalysis(alpha=0.0, ties="breslow")
        except TypeError:
            model = CoxPHSurvivalAnalysis(alpha=0.0)
        model.fit(x, y_tau)
        return float(np.asarray(model.coef_).reshape(-1)[0])
    except Exception as exc:
        warnings.warn(f"Calibration slope failed at t={tau}: {exc}")
        return np.nan


def calibration_groups(y_test_full, pred_risk, model_label, tau, n_groups=10):
    pred_risk = np.asarray(pred_risk, dtype=float)
    ok = np.isfinite(pred_risk)
    if int(np.sum(ok)) < max(20, n_groups * 3):
        return pd.DataFrame()

    y_event = np.asarray(y_test_full["event"])[ok]
    y_time = np.asarray(y_test_full["time"])[ok]
    risk = pred_risk[ok]
    rank = pd.Series(risk).rank(method="first")
    q = int(min(n_groups, len(risk)))
    groups = pd.qcut(rank, q=q, labels=False, duplicates="drop").to_numpy()

    rows = []
    for g in sorted(pd.unique(groups)):
        ind = groups == g
        yg = Surv.from_arrays(event=y_event[ind], time=y_time[ind])
        rows.append(
            {
                "model_horizon": model_label,
                "evaluation_horizon_years": float(tau),
                "calibration_group": int(g) + 1,
                "n": int(np.sum(ind)),
                "mean_predicted_risk": float(np.mean(risk[ind])),
                "observed_km_risk": km_observed_risk(yg, tau),
            }
        )
    return pd.DataFrame(rows)


def evaluate_models_common_horizons(
    models,
    model_train_horizon_years,
    risks_test,
    X_test,
    y_trainval_full,
    y_test_full,
    evaluation_times,
    n_calibration_groups,
    ibs_start_years,
    ibs_grid_points,
):
    rows = []
    cal_frames = []

    full_test_harrell = {
        label: compute_harrell_c(y_test_full, risks_test[label]) for label in models
    }

    for label, model in models.items():
        train_h = model_train_horizon_years[label]
        for tau in evaluation_times:
            uno = safe_uno_c(y_trainval_full, y_test_full, risks_test[label], tau)
            auc = safe_td_auc(y_trainval_full, y_test_full, risks_test[label], tau)
            cal_slope = safe_calibration_slope(y_test_full, risks_test[label], tau)

            # A training horizon shorter than tau does not identify the baseline
            # survival beyond that horizon. Do not extrapolate absolute risk.
            supports_absolute_risk = train_h is None or float(tau) <= float(train_h) + 1e-8

            if supports_absolute_risk:
                surv_prob = predict_survival_probability(model, X_test, tau)
                pred_abs_risk = 1.0 - surv_prob
                brier = safe_brier(y_trainval_full, y_test_full, surv_prob, tau)
                ibs = safe_ibs(
                    y_trainval_full,
                    y_test_full,
                    model,
                    X_test,
                    tau,
                    start_years=ibs_start_years,
                    grid_points=ibs_grid_points,
                )
                obs = km_observed_risk(y_test_full, tau)
                mean_pred = (
                    float(np.nanmean(pred_abs_risk))
                    if np.isfinite(pred_abs_risk).any()
                    else np.nan
                )
                cal_in_large = mean_pred - obs if np.isfinite(mean_pred) and np.isfinite(obs) else np.nan
                cal = calibration_groups(
                    y_test_full,
                    pred_abs_risk,
                    label,
                    tau,
                    n_groups=n_calibration_groups,
                )
                if not cal.empty:
                    cal_frames.append(cal)
            else:
                brier = ibs = obs = mean_pred = cal_in_large = np.nan

            event = np.asarray(y_test_full["event"]).astype(bool)
            time = np.asarray(y_test_full["time"]).astype(float)
            rows.append(
                {
                    "model_horizon": label,
                    "evaluation_horizon_years": float(tau),
                    "n_test": int(len(y_test_full)),
                    "n_events_by_evaluation_horizon": int(np.sum(event & (time <= float(tau)))),
                    "harrell_c_full_followup": float(full_test_harrell[label]),
                    "uno_c": uno,
                    "time_dependent_auc": auc,
                    "absolute_risk_supported": bool(supports_absolute_risk),
                    "brier_score": brier,
                    "ibs_from_start_to_horizon": ibs,
                    "ibs_start_years": float(ibs_start_years),
                    "mean_predicted_absolute_risk": mean_pred,
                    "observed_km_risk": obs,
                    "calibration_in_large_pred_minus_observed": cal_in_large,
                    "cox_calibration_slope": cal_slope,
                }
            )

    cal_df = pd.concat(cal_frames, ignore_index=True) if cal_frames else pd.DataFrame()
    return pd.DataFrame(rows), cal_df


def paired_bootstrap_delta_uno(
    y_trainval_full,
    y_test_full,
    risk_a,
    risk_b,
    label_a,
    label_b,
    tau,
    n_boot,
    random_state,
):
    risk_a = np.asarray(risk_a, dtype=float).reshape(-1)
    risk_b = np.asarray(risk_b, dtype=float).reshape(-1)
    point_a = safe_uno_c(y_trainval_full, y_test_full, risk_a, tau)
    point_b = safe_uno_c(y_trainval_full, y_test_full, risk_b, tau)
    delta = point_a - point_b if np.isfinite(point_a) and np.isfinite(point_b) else np.nan

    rng = np.random.default_rng(random_state)
    event = np.asarray(y_test_full["event"]).astype(bool)
    time = np.asarray(y_test_full["time"]).astype(float)
    boots = []
    n = len(time)

    for _ in range(int(n_boot)):
        idx = rng.integers(0, n, size=n)
        if np.sum(event[idx] & (time[idx] <= float(tau))) < 2:
            continue
        yb = Surv.from_arrays(event=event[idx], time=time[idx])
        try:
            ca = concordance_index_ipcw(
                y_trainval_full, yb, risk_a[idx], tau=float(tau)
            )[0]
            cb = concordance_index_ipcw(
                y_trainval_full, yb, risk_b[idx], tau=float(tau)
            )[0]
            d = float(ca - cb)
            if np.isfinite(d):
                boots.append(d)
        except Exception:
            continue

    boots = np.asarray(boots, dtype=float)
    if boots.size:
        lo, hi = np.quantile(boots, [0.025, 0.975])
        p_le0 = float(np.mean(boots <= 0.0))
        p_ge0 = float(np.mean(boots >= 0.0))
        p2 = float(min(1.0, 2.0 * min(p_le0, p_ge0)))
    else:
        lo = hi = p2 = p_le0 = np.nan

    return {
        "comparison": f"{label_a}_minus_{label_b}",
        "model_a": label_a,
        "model_b": label_b,
        "evaluation_horizon_years": float(tau),
        "uno_c_a": point_a,
        "uno_c_b": point_b,
        "delta_uno_c_a_minus_b": delta,
        "delta_uno_c_ci_lower": float(lo) if np.isfinite(lo) else np.nan,
        "delta_uno_c_ci_upper": float(hi) if np.isfinite(hi) else np.nan,
        "n_bootstrap_requested": int(n_boot),
        "n_bootstrap_successful": int(boots.size),
        "empirical_p_two_sided_delta_not_equal_0": p2,
        "empirical_p_one_sided_delta_le_0": p_le0,
    }


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    args = parse_args()
    organ = clean_name(args.organ)
    pref = output_prefix(organ)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    admin_censor_date = pd.to_datetime(args.admin_censor_date)
    horizons = parse_horizons(args.horizons)
    horizon_labels = [h[0] for h in horizons]
    horizon_years = {label: years for label, years in horizons}
    evaluation_times = parse_float_list(args.evaluation_times, "--evaluation-times")
    l1_ratios = tuple(float(x) for x in args.l1_ratios.split(","))
    final_alpha_backoff_multipliers = parse_float_list(
        args.final_alpha_backoff_multipliers, "--final-alpha-backoff-multipliers"
    )

    print("============================================================")
    print("Endocrine metabolomics mortality EPOCH horizon experiment")
    print(f"Organ: {organ}")
    print(f"Training horizons: {horizons}")
    print(f"Common evaluation times: {evaluation_times}")
    print("============================================================")

    # ----- Load and merge once -----
    print("Loading death/assessment data...")
    death = load_death_data(args.death_xlsx, args.id_match_csv)
    death["admin_censor_date"] = admin_censor_date

    print(f"Loading {organ} metabolomics data...")
    organ_df = load_organ_data(args.organ_tsv, organ, args.imaging_session_id)
    organ_feature_cols = infer_feature_columns(organ_df, args.feature_start_column)

    print("Loading optional covariates...")
    cov = load_covariates(args.covariate_csv)

    print("Merging data...")
    df = organ_df.merge(death, on="participant_id", how="inner")
    if cov is not None:
        df = df.merge(cov, on="participant_id", how="left", suffixes=("", "_cov"))

    print("Constructing ONE prospective baseline cohort and full-follow-up outcome...")
    df = construct_full_survival_dataset(df)
    df = df.loc[df["time_days_full"] >= args.min_followup_days].copy()

    # Add every requested horizon outcome to this same cohort.
    for label, years in horizons:
        df = add_horizon_outcome(df, label, years)

    print("Adding generic covariates...")
    df = add_basic_covariates(df)
    df, numeric_cols, categorical_cols, _, organ_feature_cols = build_design_matrix(
        df, organ_feature_cols
    )
    df = df.dropna(
        subset=["participant_id", "time_years_full", "event_full", "age_at_baseline", "sex"]
    ).copy()

    # Verify all requested outcome columns are valid.
    for label in horizon_labels:
        event_col, time_col = outcome_columns(label)
        if event_col not in df.columns or time_col not in df.columns:
            raise RuntimeError(f"Missing outcome columns for horizon {label}")
        df = df.loc[df[time_col] > 0].copy()

    print("Common analytic cohort:")
    print(f"  N = {df.shape[0]}")
    print(f"  Full-follow-up deaths = {int(df['event_full'].sum())}")
    print(f"  Median full follow-up years = {df['time_years_full'].median():.2f}")
    print(f"  Maximum full follow-up years = {df['time_years_full'].max():.2f}")
    for label in horizon_labels:
        e, t = outcome_columns(label)
        print(
            f"  {label}: deaths={int(df[e].sum())}, "
            f"median follow-up={df[t].median():.2f}y, max={df[t].max():.2f}y"
        )

    # ----- ONE common split -----
    split_key = str(args.split_stratify_horizon).strip().lower()
    if split_key in {"full", "max", "maximum", "all"}:
        strat_event_col = "event_full"
        split_key = "full"
    else:
        split_year = float(split_key)
        split_key = f"{split_year:g}y"
        strat_event_col = f"event_{split_key}"
        if strat_event_col not in df.columns:
            warnings.warn(
                f"Requested split stratification horizon {split_key} was not trained; falling back to full."
            )
            strat_event_col = "event_full"
            split_key = "full"

    print(f"Creating ONE train/validation/test split, stratified on {strat_event_col} + age bins...")
    strat_all = make_stratify_vector(df, strat_event_col, args.stratify_age_bins)
    try:
        df_trainval, df_test = train_test_split(
            df,
            test_size=args.test_size,
            random_state=args.random_state,
            stratify=strat_all,
        )
    except ValueError as exc:
        warnings.warn(f"Age+event stratification failed ({exc}); retrying with event only.")
        df_trainval, df_test = train_test_split(
            df,
            test_size=args.test_size,
            random_state=args.random_state,
            stratify=df[strat_event_col] if df[strat_event_col].nunique() > 1 else None,
        )

    strat_trainval = make_stratify_vector(
        df_trainval, strat_event_col, args.stratify_age_bins
    )
    try:
        df_train, df_val = train_test_split(
            df_trainval,
            test_size=args.validation_size,
            random_state=args.random_state,
            stratify=strat_trainval,
        )
    except ValueError as exc:
        warnings.warn(f"Train/validation stratification failed ({exc}); retrying with event only.")
        df_train, df_val = train_test_split(
            df_trainval,
            test_size=args.validation_size,
            random_state=args.random_state,
            stratify=(
                df_trainval[strat_event_col]
                if df_trainval[strat_event_col].nunique() > 1
                else None
            ),
        )

    print(f"  Train N={df_train.shape[0]}")
    print(f"  Val   N={df_val.shape[0]}")
    print(f"  Test  N={df_test.shape[0]}")
    for label in horizon_labels:
        event_col, _ = outcome_columns(label)
        print(
            f"  {label} events: train={int(df_train[event_col].sum())}, "
            f"val={int(df_val[event_col].sum())}, test={int(df_test[event_col].sum())}"
        )

    # Save shared split assignments early.
    split_assign = pd.DataFrame(
        {
            "participant_id": df["participant_id"].values,
            "split": "",
        },
        index=df.index,
    )
    split_assign.loc[df_train.index, "split"] = "train"
    split_assign.loc[df_val.index, "split"] = "validation"
    split_assign.loc[df_test.index, "split"] = "test"
    split_assign.reset_index(drop=True).to_csv(
        outdir / f"{pref}_split_assignments.tsv", sep="\t", index=False
    )

    # ----- ONE missingness screen and ONE preprocessor -----
    all_cols = numeric_cols + categorical_cols
    X_train_raw = df_train[all_cols].copy()
    X_val_raw = df_val[all_cols].copy()
    X_test_raw = df_test[all_cols].copy()
    X_trainval_raw = df_trainval[all_cols].copy()

    X_train_raw, other, numeric_cols_kept, categorical_cols_kept, dropped_numeric = (
        drop_high_missing_features(
            X_train_raw,
            [X_val_raw, X_test_raw, X_trainval_raw],
            numeric_cols,
            categorical_cols,
            args.max_feature_missing,
        )
    )
    X_val_raw, X_test_raw, X_trainval_raw = other

    organ_feature_cols_kept = [c for c in organ_feature_cols if c in numeric_cols_kept]
    residualization_covariates = [
        c
        for c in (numeric_cols_kept + categorical_cols_kept)
        if c not in organ_feature_cols
    ]

    print("Shared retained non-organ covariates:")
    for c in residualization_covariates:
        print(f"  {c}")
    print(
        f"Retained metabolomics features after training-set missingness filter: "
        f"{len(organ_feature_cols_kept)} / {len(organ_feature_cols)}"
    )

    preprocessor = make_preprocessor(numeric_cols_kept, categorical_cols_kept)
    X_train = preprocessor.fit_transform(X_train_raw)
    X_val = preprocessor.transform(X_val_raw)
    X_test = preprocessor.transform(X_test_raw)
    X_trainval = preprocessor.transform(X_trainval_raw)
    feature_names = get_feature_names(preprocessor)

    # Full-follow-up survival arrays are the common reference for Step 3.
    y_trainval_full = survival_from_df(df_trainval, "full")
    y_test_full = survival_from_df(df_test, "full")

    # Prepare prediction frame once.
    base_cols = [
        "participant_id",
        "baseline_date",
        "sample_date",
        "death_date",
        "admin_censor_date",
        "age_at_baseline",
        "age_at_imaging",
        "sex",
    ]
    if "organ_source_file" in df.columns:
        base_cols.append("organ_source_file")
    extra_covs = [c for c in residualization_covariates if c in df.columns and c not in base_cols]
    survival_cols = []
    for label in horizon_labels:
        e, t = outcome_columns(label)
        survival_cols.extend([e, t])
    for c in ["event_full", "time_years_full"]:
        if c not in survival_cols:
            survival_cols.append(c)

    pred_all = df[base_cols + extra_covs + survival_cols].copy()
    pred_all["split"] = ""
    pred_all.loc[df_train.index, "split"] = "train"
    pred_all.loc[df_val.index, "split"] = "validation"
    pred_all.loc[df_test.index, "split"] = "test"

    # Containers shared across horizons.
    models = {}
    best_by_horizon = {}
    penalty_factor_by_horizon = {}
    risks = {"train": {}, "validation": {}, "test": {}, "trainval": {}}
    coefficients = {}
    clock_transform_info = {}
    final_fit_info_by_horizon = {}
    fit_summary_rows = []

    # ----- Train each horizon-specific clock -----
    for h_index, (label, years) in enumerate(horizons):
        print("\n" + "=" * 72)
        print(f"STEP 1: training mortality EPOCH horizon = {label}")
        print("=" * 72)

        y_train = survival_from_df(df_train, label)
        y_val = survival_from_df(df_val, label)
        y_test = survival_from_df(df_test, label)
        y_trainval = survival_from_df(df_trainval, label)

        if np.sum(y_train["event"]) < 20:
            warnings.warn(f"{label}: very few training mortality events; model may be unstable.")

        print("Tuning elastic-net Cox model...")
        best, penalty_factor = fit_and_select_coxnet(
            X_train,
            y_train,
            X_val,
            y_val,
            feature_names,
            organ_feature_cols,
            organ,
            l1_ratios,
            args.n_alphas,
        )
        print("Best validation model:")
        print(json.dumps({k: v for k, v in best.items() if k != "coef"}, indent=2))

        print("Refitting final model on shared train+validation participants...")
        model, final_fit_info = fit_final_model(
            X_trainval,
            y_trainval,
            best,
            penalty_factor,
            final_alpha_backoff_multipliers,
        )

        models[label] = model
        best_by_horizon[label] = {k: v for k, v in best.items() if k != "coef"}
        final_fit_info_by_horizon[label] = final_fit_info
        penalty_factor_by_horizon[label] = penalty_factor

        for split_name, X in [
            ("train", X_train),
            ("validation", X_val),
            ("test", X_test),
            ("trainval", X_trainval),
        ]:
            risks[split_name][label] = predict_risk_score(model, X)

        coef_arr = np.asarray(model.coef_)
        coef = coef_arr.reshape(-1) if coef_arr.ndim == 1 else coef_arr[:, -1]
        coefficients[label] = coef

        split_y = {
            "train": y_train,
            "validation": y_val,
            "test": y_test,
            "trainval": y_trainval,
        }
        for split_name in ["train", "validation", "test", "trainval"]:
            yy = split_y[split_name]
            fit_summary_rows.append(
                {
                    "model_horizon": label,
                    "training_horizon_years": years,
                    "split": split_name,
                    "n": int(len(yy)),
                    "n_events": int(np.sum(yy["event"])),
                    "harrell_c_on_own_training_horizon_outcome": compute_harrell_c(
                        yy, risks[split_name][label]
                    ),
                    "best_l1_ratio": float(best["l1_ratio"]),
                    "selected_alpha_from_validation": float(best["alpha"]),
                    "final_alpha_used": float(final_fit_info["final_alpha_used"]),
                    "final_alpha_multiplier": float(final_fit_info["final_alpha_multiplier"]),
                    "alpha_was_increased_for_numerical_stability": bool(
                        final_fit_info["alpha_was_increased_for_numerical_stability"]
                    ),
                    "best_validation_cindex_during_tuning": float(best["cindex"]),
                    "n_nonzero_coefficients": int(np.sum(coef != 0)),
                }
            )

        # Put risk scores into the common all-participant prediction table.
        risk_col = f"{organ}_metabolomics_mortality_risk_score_{label}"
        pred_all[risk_col] = np.nan
        pred_all.loc[df_train.index, risk_col] = risks["train"][label]
        pred_all.loc[df_val.index, risk_col] = risks["validation"][label]
        pred_all.loc[df_test.index, risk_col] = risks["test"][label]

        # Also save available absolute mortality risks. No extrapolation beyond
        # the model's training horizon.
        for tau in evaluation_times:
            supports = years is None or tau <= years + 1e-8
            abs_col = f"{organ}_metabolomics_mortality_absolute_risk_{tau:g}y_from_{label}_model"
            pred_all[abs_col] = np.nan
            if supports:
                for idx, X in [
                    (df_train.index, X_train),
                    (df_val.index, X_val),
                    (df_test.index, X_test),
                ]:
                    pred_all.loc[idx, abs_col] = 1.0 - predict_survival_probability(
                        model, X, tau
                    )

        # Preserve the original EPOCH clock acceleration transform separately
        # for each mortality-training horizon.
        pred_all, transform_info = add_clock_age_and_acceleration(
            pred_all,
            organ,
            label,
            risk_col,
            residualization_covariates,
        )
        clock_transform_info[label] = transform_info

        # Save one coefficient file per horizon for direct inspection.
        is_org = organ_feature_mask(feature_names, organ_feature_cols)
        coef_df = pd.DataFrame(
            {
                "feature": feature_names,
                "coefficient": coef,
                "abs_coefficient": np.abs(coef),
                "penalty_factor": penalty_factor,
                "is_nonzero": coef != 0,
                f"is_{organ}_metabolomics_feature": is_org,
                "model_horizon": label,
            }
        ).sort_values("abs_coefficient", ascending=False)
        coef_df.to_csv(
            outdir / f"{pref}_coefficients_{label}.tsv", sep="\t", index=False
        )
        coef_df.loc[coef_df["is_nonzero"]].to_csv(
            outdir / f"{pref}_nonzero_coefficients_{label}.tsv", sep="\t", index=False
        )

    # ----- Step 2: compare scores and coefficients -----
    print("\n" + "=" * 72)
    print("STEP 2: comparing EPOCH scores and model coefficients across horizons")
    print("=" * 72)

    fit_summary_df = pd.DataFrame(fit_summary_rows)
    fit_summary_df.to_csv(outdir / f"{pref}_fit_summary.tsv", sep="\t", index=False)

    score_corr_df = build_score_correlation_tables(pred_all, organ, horizon_labels)
    score_corr_df.to_csv(
        outdir / f"{pref}_cross_horizon_score_correlations.tsv", sep="\t", index=False
    )

    coef_compare_df = coefficient_comparison(
        coefficients, feature_names, organ_feature_cols, horizon_labels
    )
    coef_compare_df.to_csv(
        outdir / f"{pref}_cross_horizon_coefficient_comparison.tsv", sep="\t", index=False
    )

    # Wide coefficient table is convenient for scatterplots later.
    coef_wide = pd.DataFrame({"feature": feature_names})
    coef_wide[f"is_{organ}_metabolomics_feature"] = organ_feature_mask(
        feature_names, organ_feature_cols
    )
    for label in horizon_labels:
        coef_wide[f"coefficient_{label}"] = coefficients[label]
    coef_wide.to_csv(
        outdir / f"{pref}_coefficients_across_horizons.tsv", sep="\t", index=False
    )

    # ----- Step 3: common held-out mortality evaluation -----
    print("\n" + "=" * 72)
    print("STEP 3: evaluating all clocks on the SAME test participants")
    print("=" * 72)

    mortality_eval_df, calibration_df = evaluate_models_common_horizons(
        models=models,
        model_train_horizon_years=horizon_years,
        risks_test=risks["test"],
        X_test=X_test,
        y_trainval_full=y_trainval_full,
        y_test_full=y_test_full,
        evaluation_times=evaluation_times,
        n_calibration_groups=args.n_calibration_groups,
        ibs_start_years=args.ibs_start_years,
        ibs_grid_points=args.ibs_grid_points,
    )
    mortality_eval_df.to_csv(
        outdir / f"{pref}_common_test_mortality_evaluation.tsv", sep="\t", index=False
    )
    if not calibration_df.empty:
        calibration_df.to_csv(
            outdir / f"{pref}_calibration_groups.tsv", sep="\t", index=False
        )

    # Primary paired comparison: Uno C differences on the same test participants.
    delta_rows = []
    for tau in evaluation_times:
        for i, (a, b) in enumerate(pairwise_labels(horizon_labels)):
            delta_rows.append(
                paired_bootstrap_delta_uno(
                    y_trainval_full=y_trainval_full,
                    y_test_full=y_test_full,
                    risk_a=risks["test"][a],
                    risk_b=risks["test"][b],
                    label_a=a,
                    label_b=b,
                    tau=tau,
                    n_boot=args.n_bootstrap_comparison,
                    random_state=args.random_state + int(round(tau * 100)) + i,
                )
            )
    delta_uno_df = pd.DataFrame(delta_rows)
    delta_uno_df.to_csv(
        outdir / f"{pref}_paired_delta_uno_c.tsv", sep="\t", index=False
    )

    # Save outcome dataset and final predictions only after every horizon is added.
    survival_keep = [
        "participant_id",
        "baseline_date",
        "sample_date",
        "death_date",
        "admin_censor_date",
        "age_at_baseline",
        "sex",
    ]
    for label in horizon_labels:
        e, t = outcome_columns(label)
        survival_keep.extend([e, t])
    survival_keep = list(dict.fromkeys([c for c in survival_keep if c in df.columns]))
    df[survival_keep].to_csv(
        outdir / f"{pref}_survival_dataset.tsv", sep="\t", index=False
    )

    pred_all.to_csv(outdir / f"{pref}_predictions.tsv", sep="\t", index=False)
    pred_all.loc[pred_all["split"] == "test"].to_csv(
        outdir / f"{pref}_test_predictions.tsv", sep="\t", index=False
    )

    # Horizon-level cohort/event summary.
    event_summary_rows = []
    for label, years in horizons:
        e, t = outcome_columns(label)
        for split_name, part in [
            ("all", df),
            ("train", df_train),
            ("validation", df_val),
            ("test", df_test),
            ("trainval", df_trainval),
        ]:
            event_summary_rows.append(
                {
                    "model_horizon": label,
                    "training_horizon_years": years,
                    "split": split_name,
                    "n": int(part.shape[0]),
                    "n_events": int(part[e].sum()),
                    "n_censored": int((~part[e].astype(bool)).sum()),
                    "median_followup_years": float(part[t].median()),
                    "max_followup_years": float(part[t].max()),
                }
            )
    event_summary_df = pd.DataFrame(event_summary_rows)
    event_summary_df.to_csv(
        outdir / f"{pref}_event_summary.tsv", sep="\t", index=False
    )

    # Save a single bundle containing the shared preprocessor plus all horizon models.
    model_bundle = {
        "organ": organ,
        "out_prefix": pref,
        "study_design": (
            "Fixed baseline cohort, fixed participant split, fixed missingness filtering and "
            "preprocessing; mortality follow-up horizon changes across models."
        ),
        "horizons": horizon_years,
        "evaluation_times": evaluation_times,
        "split_stratify_horizon": split_key,
        "split_stratify_event_col": strat_event_col,
        "preprocessor": preprocessor,
        "models": models,
        "feature_names": feature_names,
        "numeric_cols_kept": numeric_cols_kept,
        "categorical_cols_kept": categorical_cols_kept,
        "organ_feature_cols": organ_feature_cols,
        "organ_feature_cols_kept": organ_feature_cols_kept,
        "dropped_numeric": dropped_numeric,
        "residualization_covariates": residualization_covariates,
        "best_by_horizon": best_by_horizon,
        "final_fit_info_by_horizon": final_fit_info_by_horizon,
        "penalty_factor_by_horizon": penalty_factor_by_horizon,
        "clock_transform_info": clock_transform_info,
        "organ_tsv_input": args.organ_tsv,
        "feature_start_column": args.feature_start_column,
        "admin_censor_date": str(admin_censor_date.date()),
        "random_state": args.random_state,
    }
    joblib.dump(model_bundle, outdir / f"{pref}_models.joblib")

    performance = {
        "organ": organ,
        "n_total": int(df.shape[0]),
        "n_train": int(df_train.shape[0]),
        "n_validation": int(df_val.shape[0]),
        "n_test": int(df_test.shape[0]),
        "n_trainval": int(df_trainval.shape[0]),
        "full_followup_deaths": int(df["event_full"].sum()),
        "median_full_followup_years": float(df["time_years_full"].median()),
        "max_full_followup_years": float(df["time_years_full"].max()),
        "horizons": horizon_years,
        "evaluation_times": evaluation_times,
        "common_split": True,
        "common_preprocessor": True,
        "split_stratify_horizon": split_key,
        "split_stratify_event_col": strat_event_col,
        "best_by_horizon": best_by_horizon,
        "final_fit_info_by_horizon": final_fit_info_by_horizon,
        "fit_summary": json.loads(fit_summary_df.to_json(orient="records")),
        "common_test_mortality_evaluation": json.loads(
            mortality_eval_df.to_json(orient="records")
        ),
        "paired_delta_uno_c": json.loads(delta_uno_df.to_json(orient="records")),
        "cross_horizon_score_correlations": json.loads(
            score_corr_df.to_json(orient="records")
        ),
        "cross_horizon_coefficient_comparison": json.loads(
            coef_compare_df.to_json(orient="records")
        ),
        "n_original_organ_features": int(len(organ_feature_cols)),
        "n_retained_organ_features": int(len(organ_feature_cols_kept)),
        "n_numeric_cols_kept": int(len(numeric_cols_kept)),
        "n_categorical_cols_kept": int(len(categorical_cols_kept)),
        "residualization_covariates": residualization_covariates,
        "organ_tsv_input": args.organ_tsv,
        "feature_start_column": args.feature_start_column,
        "admin_censor_date": str(admin_censor_date.date()),
        "time_zero": "UKB baseline assessment date, field 53-0.0",
        "event_date": "UKB death date, field 40000-0.0",
        "important_note": (
            "For cross-horizon comparison, all models are evaluated on the same held-out test "
            "participants. Absolute-risk metrics are not extrapolated beyond a model's training horizon."
        ),
    }
    with open(outdir / f"{pref}_performance.json", "w") as f:
        json.dump(json_safe(performance), f, indent=2)

    print("\nDone.")
    print(f"Outputs written to: {outdir}")
    print("Main output files:")
    for name in [
        "survival_dataset.tsv",
        "split_assignments.tsv",
        "predictions.tsv",
        "test_predictions.tsv",
        "event_summary.tsv",
        "fit_summary.tsv",
        "cross_horizon_score_correlations.tsv",
        "cross_horizon_coefficient_comparison.tsv",
        "coefficients_across_horizons.tsv",
        "common_test_mortality_evaluation.tsv",
        "calibration_groups.tsv",
        "paired_delta_uno_c.tsv",
        "models.joblib",
        "performance.json",
    ]:
        print(f"  {outdir / f'{pref}_{name}'}")


if __name__ == "__main__":
    main()