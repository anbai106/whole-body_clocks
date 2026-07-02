#!/usr/bin/env python3
"""
General organ proteomics mi clock using elastic-net Cox survival modeling.

The organ feature TSV can be one file or a comma-separated list of files. The
script concatenates them, optionally filters to session_id if requested, and treats all
columns after --feature-start-column (default: diagnosis) as organ proteomics features.
Baseline proteomics/metabolomics variables use UKB instance 0_0, and survival time zero
is the baseline assessment date, UKB field 53-0.0.

Output naming is controlled by --organ, e.g. --organ heart creates:
  heart_proteomics_mi_clock_predictions.tsv
  heart_proteomics_mi_clock_performance.json
  etc.
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
    from sksurv.metrics import concordance_index_censored
    from sksurv.util import Surv
except ImportError as e:
    raise ImportError(
        "This script requires scikit-survival. Install with:\n"
        "  conda install -c conda-forge scikit-survival"
    ) from e


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--organ", required=True, help="Organ name for dynamic output naming, e.g. heart, liver, kidney.")
    p.add_argument("--mi-xlsx", default="/cbica/home/wenju/Dataset/UKBB_UMelbourne/Death_related_var_from_Ye.xlsx")
    p.add_argument("--id-match-csv", default="/cbica/home/wenju/Dataset/UKBB_UMelbourne/UKB_UMelbourne_vs_Penn_match_key.csv")
    p.add_argument("--organ-tsv", required=True, help="One TSV or comma-separated list/globs of TSVs.")
    p.add_argument("--covariate-csv", default="/cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv")
    p.add_argument("--admin-censor-date", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--omics-session-id", "--imaging-session-id", dest="imaging_session_id", default="none", help="Baseline omics session_id to keep if present; default disables filtering.")
    p.add_argument("--feature-start-column", default="diagnosis", help="All columns after this column are organ proteomics features.")
    p.add_argument("--test-size", type=float, default=0.20)
    p.add_argument("--validation-size", type=float, default=0.20)
    p.add_argument("--random-state", type=int, default=2026)
    p.add_argument("--stratify-age-bins", type=int, default=5)
    p.add_argument("--max-feature-missing", type=float, default=0.20)
    p.add_argument("--l1-ratios", default="0.1,0.25,0.5,0.75,1.0")
    p.add_argument("--n-alphas", type=int, default=100)
    p.add_argument("--min-followup-days", type=int, default=1)
    p.add_argument("--n-bootstrap-incremental", type=int, default=1000)
    return p.parse_args()


def clean_name(x):
    x = re.sub(r"[^A-Za-z0-9]+", "_", str(x).strip().lower())
    x = re.sub(r"_+", "_", x).strip("_")
    if not x:
        raise ValueError("--organ is empty after sanitization.")
    return x


def make_onehot_encoder():
    try:
        return OneHotEncoder(handle_unknown="ignore", sparse_output=False)
    except TypeError:
        return OneHotEncoder(handle_unknown="ignore", sparse=False)


def output_prefix(organ):
    return f"{organ}_proteomics_mi_clock"


def load_mi_data(mi_xlsx, id_match_csv):
    """
    Load UKB assessment dates, mi date, and death date.

    Proteomics is treated as baseline omics, so time zero is the baseline
    assessment date from UMelbourne field 53-0.0. mi event date is
    UKB field 42000-0.0. Death date, if available, is used as a competing
    censoring time for mi-free participants.
    """
    d = pd.read_excel(mi_xlsx)
    m = pd.read_csv(id_match_csv)

    d = d.rename(columns={"eid": "participant_id_umel"})
    m = m.rename(columns={"id": "participant_id_umel", "id_upenn": "participant_id"})
    d = m.merge(d, on="participant_id_umel", how="inner")

    required = ["participant_id", "53-0.0", "42000-0.0"]
    missing = [c for c in required if c not in d.columns]
    if missing:
        raise ValueError(f"mi/assessment file is missing required columns: {missing}")

    keep = required.copy()

    # UKB death date is commonly field 40000-0.0.
    # Keep it if available so mi-free participants can be censored at death.
    if "40000-0.0" in d.columns:
        keep.append("40000-0.0")
    else:
        warnings.warn(
            "Death date field 40000-0.0 was not found. "
            "mi-free participants will be censored only at the administrative censor date."
        )

    d = d[keep].copy()
    d["baseline_date"] = pd.to_datetime(d["53-0.0"], errors="coerce")
    d["sample_date"] = d["baseline_date"]
    d["mi_date"] = pd.to_datetime(d["42000-0.0"], errors="coerce")

    if "40000-0.0" in d.columns:
        d["death_date"] = pd.to_datetime(d["40000-0.0"], errors="coerce")
    else:
        d["death_date"] = pd.NaT

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

    df = df.sort_values(["participant_id", "organ_source_order", "organ_source_row"], kind="mergesort")
    dup = int(df["participant_id"].duplicated().sum())
    if dup > 0:
        warnings.warn(f"Found {dup} duplicated participant_id rows. Keeping first row by input order.")
        df = df.drop_duplicates("participant_id", keep="first")
    return df


def infer_feature_columns(df_organ, feature_start_column):
    if feature_start_column not in df_organ.columns:
        raise ValueError(f"Feature start column '{feature_start_column}' not found. First columns: {list(df_organ.columns[:20])}")
    start = list(df_organ.columns).index(feature_start_column) + 1
    excluded = {"organ_source_file", "organ_source_order", "organ_source_row", "participant_id", "session_id", feature_start_column}
    features = [c for c in list(df_organ.columns[start:]) if c not in excluded]
    if not features:
        raise ValueError(f"No proteomics features found after {feature_start_column}.")
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


def construct_survival_dataset(df):
    """
    Construct prospective mi survival outcome using baseline sample date as time zero.

    Event:
        Incident mi after baseline sample date and before censoring.

    Exclusions:
        - Missing baseline/sample date
        - mi before or on baseline sample date
        - Baseline/sample date after censoring

    Censoring:
        - Administrative censor date
        - Death date, if available and before administrative censor date
    """
    df = df.copy()

    df["mi_before_or_on_sample"] = (
        df["mi_date"].notna()
        & df["sample_date"].notna()
        & (df["mi_date"] <= df["sample_date"])
    )

    n_pre = int(df["mi_before_or_on_sample"].sum())
    if n_pre:
        warnings.warn(f"Excluding {n_pre} participants with mi before/on baseline sample date.")

    df = df.loc[df["sample_date"].notna()].copy()
    df = df.loc[~df["mi_before_or_on_sample"]].copy()

    # Censor at the earlier of admin censor date and death date.
    df["censor_date"] = df["admin_censor_date"]

    if "death_date" in df.columns:
        has_death_before_admin = (
            df["death_date"].notna()
            & (df["death_date"] < df["admin_censor_date"])
        )
        df.loc[has_death_before_admin, "censor_date"] = df.loc[has_death_before_admin, "death_date"]

    # Remove participants whose baseline/sample date is after censoring.
    df = df.loc[df["sample_date"] <= df["censor_date"]].copy()

    # Incident mi must occur after sample date and on/before censor date.
    df["event"] = (
        df["mi_date"].notna()
        & (df["mi_date"] > df["sample_date"])
        & (df["mi_date"] <= df["censor_date"])
    )

    df["end_date"] = df["censor_date"]
    df.loc[df["event"], "end_date"] = df.loc[df["event"], "mi_date"]

    df["time_days"] = (df["end_date"] - df["sample_date"]).dt.days
    df["time_years"] = df["time_days"] / 365.25

    df = df.loc[df["time_days"] > 0].copy()

    return df


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
    df["sex"] = df["sex"].replace({"0": "Female", "0.0": "Female", "1": "Male", "1.0": "Male", "F": "Female", "M": "Male", "female": "Female", "male": "Male"})

    if "body_mass_index_bmi_f23104_0_0" in df.columns:
        df["bmi_at_baseline"] = pd.to_numeric(df["body_mass_index_bmi_f23104_0_0"], errors="coerce")
    df = mean_existing_numeric_columns(df, ["diastolic_blood_pressure_automated_reading_f4079_0_0", "diastolic_blood_pressure_automated_reading_f4079_0_1"], "diastolic_bp_at_baseline")
    df = mean_existing_numeric_columns(df, ["systolic_blood_pressure_automated_reading_f4080_0_0", "systolic_blood_pressure_automated_reading_f4080_0_1"], "systolic_bp_at_baseline")
    if "smoking_status_f20116_0_0" in df.columns:
        df["smoking_status_at_baseline"] = df["smoking_status_f20116_0_0"].astype("category")

    if "uk_biobank_assessment_centre_f54_0_0" in df.columns:
        df["uk_biobank_assessment_centre_f54_0_0"] = df["uk_biobank_assessment_centre_f54_0_0"].astype("category")

    # Backward-compatible alias for generic plotting scripts that expect age_at_imaging.
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
        df["uk_biobank_assessment_centre_f54_0_0"] = df["uk_biobank_assessment_centre_f54_0_0"].astype("category")
        categorical_covariates.append("uk_biobank_assessment_centre_f54_0_0")
    for c in organ_feature_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df, numeric_covariates + organ_feature_cols, categorical_covariates, numeric_covariates, organ_feature_cols


def make_stratify_vector(df, age_bins=5):
    counts = df["event"].value_counts()
    if len(counts) < 2 or counts.min() < 2:
        return None
    if age_bins and age_bins > 1 and "age_at_baseline" in df.columns:
        try:
            age_bin = pd.qcut(df["age_at_baseline"].rank(method="first"), q=age_bins, labels=False, duplicates="drop")
            lab = df["event"].astype(int).astype(str) + "_age" + age_bin.astype(str)
            if lab.value_counts().min() >= 2:
                return lab
        except Exception:
            pass
    return df["event"]


def make_preprocessor(numeric_cols, categorical_cols):
    num = Pipeline([("imputer", SimpleImputer(strategy="median")), ("scaler", StandardScaler())])
    cat = Pipeline([("imputer", SimpleImputer(strategy="most_frequent")), ("onehot", make_onehot_encoder())])
    return ColumnTransformer([("num", num, numeric_cols), ("cat", cat, categorical_cols)], remainder="drop")


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


def drop_high_missing_features(X_train_raw, other_raw_list, numeric_cols, categorical_cols, max_missing):
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
    return X_train_raw[cols].copy(), [x[cols].copy() for x in other_raw_list], keep_numeric, categorical_cols, dropped


def organ_feature_mask(feature_names, organ_feature_cols):
    organ_set = {f"num__{c}" for c in organ_feature_cols}
    return np.array([str(f) in organ_set for f in feature_names], dtype=bool)


def compute_cindex(y, risk):
    return float(concordance_index_censored(y["event"], y["time"], np.asarray(risk).reshape(-1))[0])


def fit_and_select_coxnet(X_train, y_train, X_val, y_val, feature_names, organ_feature_cols, organ, l1_ratios, n_alphas):
    penalty_factor = np.ones(len(feature_names), dtype=float)
    is_organ = organ_feature_mask(feature_names, organ_feature_cols)
    penalty_factor[~is_organ] = 0.0
    best = {"cindex": -np.inf, "l1_ratio": None, "alpha": None, "coef": None, "used_penalty_factor": True}
    for l1_ratio in l1_ratios:
        print(f"Fitting Coxnet path for l1_ratio={l1_ratio}")
        try:
            model = CoxnetSurvivalAnalysis(l1_ratio=l1_ratio, n_alphas=n_alphas, alpha_min_ratio="auto", penalty_factor=penalty_factor, fit_baseline_model=False, max_iter=100000)
        except TypeError:
            warnings.warn("Installed scikit-survival does not support penalty_factor. Covariates will be penalized.")
            model = CoxnetSurvivalAnalysis(l1_ratio=l1_ratio, n_alphas=n_alphas, alpha_min_ratio="auto", fit_baseline_model=False, max_iter=100000)
            best["used_penalty_factor"] = False
        model.fit(X_train, y_train)
        coefs = model.coef_
        if coefs.ndim == 1:
            coefs = coefs[:, None]
        for j, alpha in enumerate(model.alphas_):
            risk_val = np.dot(X_val, coefs[:, j])
            cindex = compute_cindex(y_val, risk_val)
            if np.isfinite(cindex) and cindex > best["cindex"]:
                best.update({"cindex": float(cindex), "l1_ratio": float(l1_ratio), "alpha": float(alpha), "coef": coefs[:, j].copy()})
        print(f"  best so far: C-index={best['cindex']:.4f}, l1_ratio={best['l1_ratio']}, alpha={best['alpha']}")
    if best["alpha"] is None:
        raise RuntimeError("Failed to select a Coxnet model.")
    return best, penalty_factor


def fit_final_model(X_trainval, y_trainval, best, penalty_factor):
    try:
        model = CoxnetSurvivalAnalysis(l1_ratio=best["l1_ratio"], alphas=[best["alpha"]], penalty_factor=penalty_factor, fit_baseline_model=True, max_iter=100000)
    except TypeError:
        model = CoxnetSurvivalAnalysis(l1_ratio=best["l1_ratio"], alphas=[best["alpha"]], fit_baseline_model=True, max_iter=100000)
    model.fit(X_trainval, y_trainval)
    return model


def predict_risk_score(model, X):
    return np.asarray(model.predict(X)).reshape(-1)


def fit_coxph_or_ridge_fallback(X, y, model_name="cox_model"):
    try:
        try:
            model = CoxPHSurvivalAnalysis(alpha=0.0, ties="breslow")
        except TypeError:
            model = CoxPHSurvivalAnalysis(alpha=0.0)
        model.fit(X, y)
        return model, "CoxPHSurvivalAnalysis(alpha=0.0)"
    except Exception as exc:
        warnings.warn(f"{model_name}: unpenalized CoxPH failed ({exc}). Falling back to weakly penalized Coxnet.")
    last = None
    for alpha in [1e-6, 1e-5, 1e-4, 1e-3, 1e-2]:
        try:
            model = CoxnetSurvivalAnalysis(l1_ratio=0.01, alphas=[alpha], fit_baseline_model=True, max_iter=100000)
            model.fit(X, y)
            return model, f"CoxnetSurvivalAnalysis(l1_ratio=0.01, alpha={alpha})"
        except Exception as exc:
            last = exc
    raise RuntimeError(f"{model_name}: all Cox fitting attempts failed. Last error: {last}")


def paired_bootstrap_delta_cindex(y, risk_full, risk_base, organ, n_boot, random_state):
    rng = np.random.default_rng(random_state)
    event = np.asarray(y["event"]).astype(bool)
    time = np.asarray(y["time"]).astype(float)
    risk_full = np.asarray(risk_full).reshape(-1)
    risk_base = np.asarray(risk_base).reshape(-1)
    c_full = compute_cindex(y, risk_full)
    c_base = compute_cindex(y, risk_base)
    delta = c_full - c_base
    boots = []
    n = len(time)
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        if np.sum(event[idx]) < 2:
            continue
        yb = Surv.from_arrays(event=event[idx], time=time[idx])
        try:
            d = compute_cindex(yb, risk_full[idx]) - compute_cindex(yb, risk_base[idx])
            if np.isfinite(d):
                boots.append(d)
        except Exception:
            continue
    boots = np.asarray(boots, dtype=float)
    comparison = f"M3_full_covariates_plus_{organ}_proteomics_vs_M1_covariate_baseline"
    if boots.size == 0:
        return {"comparison": comparison, "cindex_full": c_full, "cindex_baseline": c_base, "delta_cindex": delta, "delta_cindex_ci_lower": np.nan, "delta_cindex_ci_upper": np.nan, "n_bootstrap_requested": int(n_boot), "n_bootstrap_successful": 0, "empirical_p_two_sided_delta_not_equal_0": np.nan, "empirical_p_one_sided_delta_le_0": np.nan, "interpretation": "Bootstrap failed."}
    lo, hi = np.quantile(boots, [0.025, 0.975])
    p_le0 = float(np.mean(boots <= 0.0))
    p_ge0 = float(np.mean(boots >= 0.0))
    p2 = float(min(1.0, 2.0 * min(p_le0, p_ge0)))
    if delta > 0 and lo > 0:
        interp = f"{organ} proteomics improves test-set C-index beyond the covariate baseline."
    elif delta > 0:
        interp = f"{organ} proteomics has positive delta C-index, but the bootstrap CI includes or approaches zero."
    else:
        interp = f"No evidence that {organ} proteomics improves test-set C-index beyond the covariate baseline."
    return {"comparison": comparison, "cindex_full": float(c_full), "cindex_baseline": float(c_base), "delta_cindex": float(delta), "delta_cindex_ci_lower": float(lo), "delta_cindex_ci_upper": float(hi), "n_bootstrap_requested": int(n_boot), "n_bootstrap_successful": int(boots.size), "empirical_p_two_sided_delta_not_equal_0": p2, "empirical_p_one_sided_delta_le_0": p_le0, "interpretation": interp}


def run_incremental_value_analysis(X_train, X_val, X_test, X_trainval, y_train, y_val, y_test, y_trainval, feature_names, organ_feature_cols, organ, final_model, n_bootstrap, random_state):
    is_org = organ_feature_mask(feature_names, organ_feature_cols)
    is_age = feature_names == "num__age_at_baseline"
    is_sex = np.char.startswith(feature_names.astype(str), "cat__sex")
    m2 = f"M2_{organ}_proteomics_only"
    m3 = f"M3_full_covariates_plus_{organ}_proteomics"
    idx = {
        "M0_age_sex": np.where(is_age | is_sex)[0],
        "M1_covariate_baseline": np.where(~is_org)[0],
        m2: np.where(is_org)[0],
        m3: np.arange(len(feature_names)),
    }
    split_data = {"train": (X_train, y_train), "validation": (X_val, y_val), "test": (X_test, y_test), "trainval": (X_trainval, y_trainval)}
    fitting_methods = {}
    risk_predictions = {s: {} for s in split_data}
    fitted_models = {}
    for model_name in ["M0_age_sex", "M1_covariate_baseline", m2]:
        ind = idx[model_name]
        print(f"Fitting incremental model: {model_name} ({len(ind)} features)")
        model, method = fit_coxph_or_ridge_fallback(X_trainval[:, ind], y_trainval, model_name=model_name)
        fitted_models[model_name] = model
        fitting_methods[model_name] = method
        for split, (X, _) in split_data.items():
            risk_predictions[split][model_name] = predict_risk_score(model, X[:, ind])
    fitting_methods[m3] = "Selected elastic-net Cox model from main pipeline"
    for split, (X, _) in split_data.items():
        risk_predictions[split][m3] = predict_risk_score(final_model, X)
    labels = {"M0_age_sex": "Age + sex", "M1_covariate_baseline": "Covariate baseline", m2: f"{organ} proteomics only", m3: f"Covariates + {organ} proteomics"}
    rows = []
    for model_name in ["M0_age_sex", "M1_covariate_baseline", m2, m3]:
        for split, (_, y) in split_data.items():
            rows.append({"model": model_name, "model_label": labels[model_name], "split": split, "n_features": int(len(idx[model_name])), "training_data": "train+validation", "fitting_method": fitting_methods[model_name], "cindex": compute_cindex(y, risk_predictions[split][model_name]), "n": int(len(y)), "n_events": int(np.sum(y["event"]))})
    comp = pd.DataFrame(rows)
    delta = pd.DataFrame([paired_bootstrap_delta_cindex(y_test, risk_predictions["test"][m3], risk_predictions["test"]["M1_covariate_baseline"], organ, n_bootstrap, random_state)])
    return {"model_comparison_df": comp, "delta_cindex_df": delta, "fitted_models": fitted_models, "fitting_methods": fitting_methods, "risk_predictions": risk_predictions, "feature_indices": idx, "model_m2": m2, "model_m3": m3}


def predict_absolute_risk(model, X, times_years):
    out = {}
    try:
        surv_funcs = model.predict_survival_function(X)
    except Exception as e:
        warnings.warn(f"Could not compute absolute risks: {e}")
        for t in times_years:
            out[f"risk_{t:g}y"] = np.repeat(np.nan, X.shape[0])
        return pd.DataFrame(out)
    for t in times_years:
        vals = []
        for sf in surv_funcs:
            try:
                vals.append(1.0 - float(sf(t)))
            except Exception:
                vals.append(np.nan)
        out[f"risk_{t:g}y"] = vals
    return pd.DataFrame(out)


def add_clock_age_and_acceleration(pred_df, organ, covariate_cols):
    df = pred_df.copy()
    risk_col = f"{organ}_proteomics_mi_risk_score"
    z_col = f"{organ}_proteomics_mi_clock_acceleration_z"
    yrs_col = f"{organ}_proteomics_mi_clock_acceleration_years"
    age_col = f"{organ}_proteomics_mi_clock_age_years"
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
        transformers.append(("num", Pipeline([("imputer", SimpleImputer(strategy="median"))]), numeric_covs))
    if categorical_covs:
        transformers.append(("cat", Pipeline([("imputer", SimpleImputer(strategy="most_frequent")), ("onehot", make_onehot_encoder())]), categorical_covs))
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
        warnings.warn("Adjusted age coefficient is near zero or unavailable. Year-scale acceleration set to missing.")
        df[yrs_col] = np.nan
        df[age_col] = np.nan
    info = {"residualization_covariates": covariate_cols, "numeric_residualization_covariates": numeric_covs, "categorical_residualization_covariates": categorical_covs, "risk_score_covariate_model_intercept": float(lr.intercept_), "risk_score_covariate_model_coef": {k: float(v) for k, v in zip(feat_names, lr.coef_)}, "adjusted_age_coefficient_risk_score_per_year": float(beta_age) if np.isfinite(beta_age) else None, "risk_score_residual_mean_train": mean_train, "risk_score_residual_sd_train": sd_train, "note": "Clock acceleration is the residual of the Cox risk score after adjustment for retained non-organ covariates."}
    return df, info


def main():
    args = parse_args()
    organ = clean_name(args.organ)
    pref = output_prefix(organ)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    admin_censor_date = pd.to_datetime(args.admin_censor_date)
    l1_ratios = tuple(float(x) for x in args.l1_ratios.split(","))
    print(f"Building {organ} proteomics mi clock")

    print("Loading mi/assessment data...")
    mi = load_mi_data(args.mi_xlsx, args.id_match_csv)
    mi["admin_censor_date"] = admin_censor_date

    print(f"Loading {organ} proteomics data...")
    organ_df = load_organ_data(args.organ_tsv, organ, args.imaging_session_id)
    organ_feature_cols = infer_feature_columns(organ_df, args.feature_start_column)

    print("Loading optional covariates...")
    cov = load_covariates(args.covariate_csv)

    print("Merging data...")
    df = organ_df.merge(mi, on="participant_id", how="inner")
    if cov is not None:
        df = df.merge(cov, on="participant_id", how="left", suffixes=("", "_cov"))

    print("Constructing survival outcome...")
    df = construct_survival_dataset(df)
    df = df.loc[df["time_days"] >= args.min_followup_days].copy()

    print("Adding generic covariates...")
    df = add_basic_covariates(df)
    df, numeric_cols, categorical_cols, _, organ_feature_cols = build_design_matrix(df, organ_feature_cols)
    df = df.dropna(subset=["participant_id", "time_years", "event", "age_at_baseline", "sex"]).copy()

    print("Final prospective baseline-omics mi dataset:")
    print(f"  N = {df.shape[0]}")
    print(f"  Incident mi events after baseline = {int(df['event'].sum())}")
    print(f"  Censored = {int((~df['event']).sum())}")
    print(f"  Median follow-up years = {df['time_years'].median():.2f}")
    if df["event"].sum() < 20:
        warnings.warn("Very few mi events. The fitted clock may be unstable.")

    dataset_cols = ["participant_id", "baseline_date", "sample_date", "mi_date", "death_date", "admin_censor_date", "censor_date", "end_date", "time_days", "time_years", "event", "age_at_baseline", "age_at_imaging", "sex"]
    if "organ_source_file" in df.columns:
        dataset_cols.append("organ_source_file")
    df[dataset_cols].to_csv(outdir / f"{pref}_survival_dataset.tsv", sep="\t", index=False)

    print("Splitting into train/validation/test...")
    df_trainval, df_test = train_test_split(df, test_size=args.test_size, random_state=args.random_state, stratify=make_stratify_vector(df, args.stratify_age_bins))
    df_train, df_val = train_test_split(df_trainval, test_size=args.validation_size, random_state=args.random_state, stratify=make_stratify_vector(df_trainval, args.stratify_age_bins))
    print(f"  Train N={df_train.shape[0]}, events={int(df_train['event'].sum())}")
    print(f"  Val   N={df_val.shape[0]}, events={int(df_val['event'].sum())}")
    print(f"  Test  N={df_test.shape[0]}, events={int(df_test['event'].sum())}")

    all_cols = numeric_cols + categorical_cols
    X_train_raw = df_train[all_cols].copy()
    X_val_raw = df_val[all_cols].copy()
    X_test_raw = df_test[all_cols].copy()
    X_trainval_raw = df_trainval[all_cols].copy()
    X_train_raw, other, numeric_cols_kept, categorical_cols_kept, dropped_numeric = drop_high_missing_features(X_train_raw, [X_val_raw, X_test_raw, X_trainval_raw], numeric_cols, categorical_cols, args.max_feature_missing)
    X_val_raw, X_test_raw, X_trainval_raw = other

    residualization_covariates = [c for c in (numeric_cols_kept + categorical_cols_kept) if c not in organ_feature_cols]
    print("Residualizing clock acceleration on retained non-organ covariates:")
    for c in residualization_covariates:
        print(f"  {c}")

    preprocessor = make_preprocessor(numeric_cols_kept, categorical_cols_kept)
    X_train = preprocessor.fit_transform(X_train_raw)
    X_val = preprocessor.transform(X_val_raw)
    X_test = preprocessor.transform(X_test_raw)
    X_trainval = preprocessor.transform(X_trainval_raw)
    feature_names = get_feature_names(preprocessor)

    y_train = Surv.from_arrays(event=df_train["event"].astype(bool).values, time=df_train["time_years"].astype(float).values)
    y_val = Surv.from_arrays(event=df_val["event"].astype(bool).values, time=df_val["time_years"].astype(float).values)
    y_test = Surv.from_arrays(event=df_test["event"].astype(bool).values, time=df_test["time_years"].astype(float).values)
    y_trainval = Surv.from_arrays(event=df_trainval["event"].astype(bool).values, time=df_trainval["time_years"].astype(float).values)

    print("Tuning elastic-net Cox model...")
    best, penalty_factor = fit_and_select_coxnet(X_train, y_train, X_val, y_val, feature_names, organ_feature_cols, organ, l1_ratios, args.n_alphas)
    print("Best validation model:")
    print(json.dumps(best | {"coef": "omitted"}, indent=2))
    print("Refitting final model on train+validation...")
    final_model = fit_final_model(X_trainval, y_trainval, best, penalty_factor)

    print("Generating predictions...")
    risk_train = predict_risk_score(final_model, X_train)
    risk_val = predict_risk_score(final_model, X_val)
    risk_test = predict_risk_score(final_model, X_test)
    risk_trainval = predict_risk_score(final_model, X_trainval)
    cindex_train = compute_cindex(y_train, risk_train)
    cindex_val = compute_cindex(y_val, risk_val)
    cindex_test = compute_cindex(y_test, risk_test)
    cindex_trainval = compute_cindex(y_trainval, risk_trainval)
    print(f"Train C-index:     {cindex_train:.4f}")
    print(f"Validation C-index:{cindex_val:.4f}")
    print(f"Train+Val C-index: {cindex_trainval:.4f}")
    print(f"Test C-index:      {cindex_test:.4f}")

    print(f"Running incremental-value analysis: does {organ} proteomics add value beyond covariates?")
    inc = run_incremental_value_analysis(X_train, X_val, X_test, X_trainval, y_train, y_val, y_test, y_trainval, feature_names, organ_feature_cols, organ, final_model, args.n_bootstrap_incremental, args.random_state)
    model_comparison_df = inc["model_comparison_df"]
    delta_cindex_df = inc["delta_cindex_df"]
    model_comparison_records = json.loads(model_comparison_df.to_json(orient="records"))
    delta_cindex_records = json.loads(delta_cindex_df.to_json(orient="records"))
    model_comparison_df.to_csv(outdir / f"{pref}_model_comparison.tsv", sep="\t", index=False)
    delta_cindex_df.to_csv(outdir / f"{pref}_incremental_value_delta_cindex.tsv", sep="\t", index=False)
    print(delta_cindex_df.to_string(index=False))

    risk_col = f"{organ}_proteomics_mi_risk_score"
    def make_pred_frame(part, split, risk):
        base = ["participant_id", "sample_date", "mi_date", "death_date", "admin_censor_date", "censor_date", "end_date", "time_years", "event", "age_at_baseline", "age_at_imaging", "sex"]
        if "organ_source_file" in part.columns:
            base.append("organ_source_file")
        extra = [c for c in residualization_covariates if c in part.columns and c not in base]
        out = part[base + extra].copy()
        out["split"] = split
        out[risk_col] = risk
        return out

    pred_train = make_pred_frame(df_train, "train", risk_train)
    pred_val = make_pred_frame(df_val, "validation", risk_val)
    pred_test = make_pred_frame(df_test, "test", risk_test)
    m2, m3 = inc["model_m2"], inc["model_m3"]
    comparison_cols = {"M0_age_sex": "risk_score_M0_age_sex", "M1_covariate_baseline": "risk_score_M1_covariate_baseline", m2: f"risk_score_M2_{organ}_proteomics_only", m3: f"risk_score_M3_full_covariates_plus_{organ}_proteomics"}
    for model_name, col in comparison_cols.items():
        pred_train[col] = inc["risk_predictions"]["train"][model_name]
        pred_val[col] = inc["risk_predictions"]["validation"][model_name]
        pred_test[col] = inc["risk_predictions"]["test"][model_name]

    risk_times = [5.0, 10.0, 15.0]
    pred_train = pd.concat([pred_train.reset_index(drop=True), predict_absolute_risk(final_model, X_train, risk_times)], axis=1)
    pred_val = pd.concat([pred_val.reset_index(drop=True), predict_absolute_risk(final_model, X_val, risk_times)], axis=1)
    pred_test = pd.concat([pred_test.reset_index(drop=True), predict_absolute_risk(final_model, X_test, risk_times)], axis=1)
    pred_all = pd.concat([pred_train, pred_val, pred_test], ignore_index=True)
    print("Adding approximate mi-clock age and acceleration...")
    pred_all, clock_transform_info = add_clock_age_and_acceleration(pred_all, organ, residualization_covariates)
    pred_all.to_csv(outdir / f"{pref}_predictions.tsv", sep="\t", index=False)
    pred_test.to_csv(outdir / f"{pref}_test_predictions.tsv", sep="\t", index=False)

    print("Saving coefficients...")
    coef = np.asarray(final_model.coef_).reshape(-1)
    is_org = organ_feature_mask(feature_names, organ_feature_cols)
    coef_df = pd.DataFrame({"feature": feature_names, "coefficient": coef, "abs_coefficient": np.abs(coef), "penalty_factor": penalty_factor, "is_nonzero": coef != 0, f"is_{organ}_proteomics_feature": is_org}).sort_values("abs_coefficient", ascending=False)
    coef_df.to_csv(outdir / f"{pref}_coefficients.tsv", sep="\t", index=False)
    nonzero_coef_df = coef_df.loc[coef_df["is_nonzero"]].copy()
    nonzero_coef_df.to_csv(outdir / f"{pref}_nonzero_coefficients.tsv", sep="\t", index=False)

    model_bundle = {"organ": organ, "out_prefix": pref, "preprocessor": preprocessor, "model": final_model, "feature_names": feature_names, "numeric_cols_kept": numeric_cols_kept, "categorical_cols_kept": categorical_cols_kept, "organ_feature_cols": organ_feature_cols, "dropped_numeric": dropped_numeric, "residualization_covariates": residualization_covariates, "best": best, "penalty_factor": penalty_factor, "clock_transform_info": clock_transform_info, "incremental_value_model_comparison": model_comparison_records, "incremental_value_delta_cindex": delta_cindex_records, "incremental_value_fitting_methods": inc["fitting_methods"], "organ_tsv_input": args.organ_tsv, "feature_start_column": args.feature_start_column, "admin_censor_date": str(admin_censor_date.date())}
    joblib.dump(model_bundle, outdir / f"{pref}_model.joblib")

    delta_row = delta_cindex_df.iloc[0]
    performance = {"organ": organ, "n_total": int(df.shape[0]), "n_events_total": int(df["event"].sum()), "n_censored_total": int((~df["event"]).sum()), "median_followup_years": float(df["time_years"].median()), "n_train": int(df_train.shape[0]), "n_events_train": int(df_train["event"].sum()), "n_validation": int(df_val.shape[0]), "n_events_validation": int(df_val["event"].sum()), "n_test": int(df_test.shape[0]), "n_events_test": int(df_test["event"].sum()), "cindex_train": float(cindex_train), "cindex_validation": float(cindex_val), "cindex_trainval": float(cindex_trainval), "cindex_test": float(cindex_test), "incremental_value_model_comparison": model_comparison_records, "incremental_value_delta_cindex": delta_cindex_records, "cindex_test_M1_covariate_baseline": float(model_comparison_df.loc[(model_comparison_df["model"] == "M1_covariate_baseline") & (model_comparison_df["split"] == "test"), "cindex"].iloc[0]), f"cindex_test_M3_full_covariates_plus_{organ}_proteomics": float(model_comparison_df.loc[(model_comparison_df["model"] == m3) & (model_comparison_df["split"] == "test"), "cindex"].iloc[0]), "delta_cindex_test_M3_vs_M1": float(delta_row["delta_cindex"]), "delta_cindex_test_M3_vs_M1_ci_lower": float(delta_row["delta_cindex_ci_lower"]), "delta_cindex_test_M3_vs_M1_ci_upper": float(delta_row["delta_cindex_ci_upper"]), "delta_cindex_test_M3_vs_M1_p_two_sided": float(delta_row["empirical_p_two_sided_delta_not_equal_0"]), "best_l1_ratio": float(best["l1_ratio"]), "best_alpha": float(best["alpha"]), "best_validation_cindex_during_tuning": float(best["cindex"]), "used_penalty_factor": bool(best["used_penalty_factor"]), "n_original_organ_features": int(len(organ_feature_cols)), f"n_original_{organ}_features": int(len(organ_feature_cols)), "n_numeric_cols_kept": int(len(numeric_cols_kept)), "n_categorical_cols_kept": int(len(categorical_cols_kept)), "n_nonzero_coefficients": int(nonzero_coef_df.shape[0]), "n_residualization_covariates": int(len(residualization_covariates)), "residualization_covariates": residualization_covariates, "organ_tsv_input": args.organ_tsv, "feature_start_column": args.feature_start_column, "admin_censor_date": str(admin_censor_date.date()), "time_zero": "UKB baseline assessment date / proteomics sample date, field 53-0.0", "event_date": "UKB mi date, field 42000-0.0", "note": f"Primary score is {organ}_proteomics_mi_risk_score from elastic-net Cox. Clock age/acceleration are post-hoc residualized transforms adjusted for retained non-organ covariates."}
    with open(outdir / f"{pref}_performance.json", "w") as f:
        json.dump(performance, f, indent=2)

    print("Done.")
    print(f"Outputs written to: {outdir}")
    print("Main output files:")
    for name in ["predictions.tsv", "test_predictions.tsv", "coefficients.tsv", "nonzero_coefficients.tsv", "model_comparison.tsv", "incremental_value_delta_cindex.tsv", "model.joblib", "performance.json"]:
        print(f"  {outdir / f'{pref}_{name}'}")


if __name__ == "__main__":
    main()
