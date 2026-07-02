#!/usr/bin/env python3
"""
Build a brain MRI dementia L'EPOCH clock in UK Biobank using elastic-net Cox survival modeling.

Time zero:
  UKB imaging assessment date, field 53-2.0

Event:
  incident dementia date, field 42018-0.0, after imaging and before censoring

Censoring:
  earlier of death date, field 40000-0.0 if available, and endpoint-specific administrative censor date

Brain-specific covariates:
  age at imaging, sex, DLICV/intracranial-volume adjustment, and optional UKB assessment center.

Primary output:
  brain_mri_dementia_risk_score from the survival model.

Recommended downstream phenotype:
  brain_mri_dementia_clock_acceleration_z, the residualized risk score after adjusting for retained non-brain covariates in the training set.
"""

import argparse
import json
import os
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
        "  conda install -c conda-forge scikit-survival\n"
        "or use your existing survival_clock conda environment."
    ) from e

PREFIX = "brain_mri_dementia_clock"
RISK_COL = "brain_mri_dementia_risk_score"
ACCEL_Z_COL = "brain_mri_dementia_clock_acceleration_z"
ACCEL_YEARS_COL = "brain_mri_dementia_clock_acceleration_years"
CLOCK_AGE_COL = "brain_mri_dementia_clock_age_years"


def parse_args():
    parser = argparse.ArgumentParser(description="Build brain MRI dementia L'EPOCH clock.")
    parser.add_argument(
        "--dementia-xlsx",
        default="/cbica/home/wenju/Dataset/UKBB_UMelbourne/Death_related_var_from_Ye.xlsx",
        help="Excel file containing UKB assessment dates, dementia date, and optionally death date.",
    )
    parser.add_argument(
        "--id-match-csv",
        default="/cbica/home/wenju/Dataset/UKBB_UMelbourne/UKB_UMelbourne_vs_Penn_match_key.csv",
        help="CSV matching UMelbourne IDs to Penn participant IDs.",
    )
    parser.add_argument(
        "--brain-tsv",
        default="/cbica/home/wenju/Reproducibile_paper/BrainAge/data/imaging/T1_MUSE_GM.tsv",
        help="Brain MRI MUSE/WMLS gray-matter volume TSV.",
    )
    parser.add_argument(
        "--covariate-csv",
        default="/cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv",
        help="Optional UKB covariate CSV.",
    )
    parser.add_argument(
        "--admin-censor-date",
        required=True,
        help="Endpoint-specific administrative censor date for dementia ascertainment, e.g. 2022-11-30.",
    )
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument(
        "--imaging-session-id",
        type=int,
        default=1,
        help="Imaging session_id in brain TSV. Default 1 for UKB imaging visit in the MUSE file.",
    )
    parser.add_argument("--test-size", type=float, default=0.20)
    parser.add_argument("--validation-size", type=float, default=0.20)
    parser.add_argument("--random-state", type=int, default=2026)
    parser.add_argument("--stratify-age-bins", type=int, default=5)
    parser.add_argument("--max-feature-missing", type=float, default=0.20)
    parser.add_argument("--l1-ratios", default="0.1,0.25,0.5,0.75,1.0")
    parser.add_argument("--n-alphas", type=int, default=100)
    parser.add_argument("--min-followup-days", type=int, default=1)
    parser.add_argument("--n-bootstrap-incremental", type=int, default=1000)
    return parser.parse_args()


def make_onehot_encoder():
    try:
        return OneHotEncoder(handle_unknown="ignore", sparse_output=False)
    except TypeError:
        return OneHotEncoder(handle_unknown="ignore", sparse=False)


def load_dementia_data(dementia_xlsx, id_match_csv):
    d = pd.read_excel(dementia_xlsx)
    m = pd.read_csv(id_match_csv)

    d = d.rename(columns={"eid": "participant_id_umel"})
    m = m.rename(columns={"id": "participant_id_umel", "id_upenn": "participant_id"})
    d = m.merge(d, on="participant_id_umel", how="inner")

    required = ["participant_id", "53-0.0", "53-2.0", "42018-0.0"]
    missing = [c for c in required if c not in d.columns]
    if missing:
        raise ValueError(f"Dementia/assessment file is missing required columns: {missing}")

    keep = required.copy()
    if "40000-0.0" in d.columns:
        keep.append("40000-0.0")
    else:
        warnings.warn(
            "Death date field 40000-0.0 was not found. Dementia-free participants "
            "will be censored only at the administrative censor date."
        )

    d = d[keep].copy()
    d["baseline_date"] = pd.to_datetime(d["53-0.0"], errors="coerce")
    d["imaging_date"] = pd.to_datetime(d["53-2.0"], errors="coerce")
    d["dementia_date"] = pd.to_datetime(d["42018-0.0"], errors="coerce")
    d["death_date"] = pd.to_datetime(d["40000-0.0"], errors="coerce") if "40000-0.0" in d.columns else pd.NaT
    return d


def load_brain_data(brain_tsv, imaging_session_id):
    b = pd.read_csv(brain_tsv, sep="\t")
    if "participant_id" not in b.columns:
        raise ValueError("Brain TSV must contain participant_id.")

    if "session_id" in b.columns and imaging_session_id is not None:
        before = b.shape[0]
        b = b.loc[b["session_id"] == imaging_session_id].copy()
        print(f"Filtered brain file to session_id={imaging_session_id}: {before} -> {b.shape[0]} rows")

    sort_cols = ["participant_id"] + (["session_id"] if "session_id" in b.columns else [])
    b = b.sort_values(sort_cols, kind="mergesort")
    n_dup = int(b["participant_id"].duplicated().sum())
    if n_dup > 0:
        warnings.warn(f"Found {n_dup} duplicated participant_id rows in brain file. Keeping first row.")
        b = b.drop_duplicates("participant_id", keep="first")
    return b


def load_covariates(covariate_csv):
    if covariate_csv is None or str(covariate_csv).lower() in {"none", ""}:
        return None
    if not os.path.exists(covariate_csv):
        warnings.warn(f"Covariate file not found: {covariate_csv}. Continuing without it.")
        return None
    cov = pd.read_csv(covariate_csv)
    if "eid" not in cov.columns:
        warnings.warn("Covariate file does not contain eid. Continuing without it.")
        return None
    return cov.rename(columns={"eid": "participant_id"})


def construct_survival_dataset(df):
    """
    Construct prospective dementia survival outcome using imaging date as time zero.

    Exclude prevalent dementia before/on imaging.
    Event is incident dementia after imaging and before censoring.
    Censor at the earlier of death date, if available, and administrative censor date.
    """
    df = df.copy()

    df["dementia_before_or_on_imaging"] = (
        df["dementia_date"].notna()
        & df["imaging_date"].notna()
        & (df["dementia_date"] <= df["imaging_date"])
    )
    n_pre = int(df["dementia_before_or_on_imaging"].sum())
    if n_pre > 0:
        warnings.warn(f"Excluding {n_pre} participants with dementia before/on imaging date.")

    df = df.loc[df["imaging_date"].notna()].copy()
    df = df.loc[~df["dementia_before_or_on_imaging"]].copy()

    df["censor_date"] = df["admin_censor_date"]
    if "death_date" in df.columns:
        died_before_admin = df["death_date"].notna() & (df["death_date"] < df["admin_censor_date"])
        df.loc[died_before_admin, "censor_date"] = df.loc[died_before_admin, "death_date"]

    df = df.loc[df["imaging_date"] <= df["censor_date"]].copy()

    df["event"] = (
        df["dementia_date"].notna()
        & (df["dementia_date"] > df["imaging_date"])
        & (df["dementia_date"] <= df["censor_date"])
    )
    df["end_date"] = df["censor_date"]
    df.loc[df["event"], "end_date"] = df.loc[df["event"], "dementia_date"]
    df["time_days"] = (df["end_date"] - df["imaging_date"]).dt.days
    df["time_years"] = df["time_days"] / 365.25
    return df


def infer_brain_feature_columns(df):
    features = [c for c in df.columns if c.startswith("MUSE_Volume_") or c.startswith("WMLS_Volume_")]
    if not features:
        raise ValueError("No brain MRI features found. Expected columns starting with MUSE_Volume_ or WMLS_Volume_.")
    return features


def add_basic_covariates(df):
    df = df.copy()

    if "Age" not in df.columns:
        raise ValueError("Brain file must contain Age column for age at imaging.")
    if "Sex" not in df.columns:
        raise ValueError("Brain file must contain Sex column.")

    df["age_at_imaging"] = pd.to_numeric(df["Age"], errors="coerce")
    df["sex"] = df["Sex"].astype(str).str.strip().replace(
        {"0": "Female", "0.0": "Female", "1": "Male", "1.0": "Male", "F": "Female", "M": "Male", "female": "Female", "male": "Male"}
    )

    # DLICV is the key brain-specific covariate for volume measures.
    if "DLICV" in df.columns:
        df["DLICV"] = pd.to_numeric(df["DLICV"], errors="coerce")
    else:
        warnings.warn("DLICV not found in brain file. Brain volumes will not be adjusted for intracranial volume in the model.")

    return df


def build_design_matrix(df, brain_feature_cols):
    df = df.copy()
    numeric_covariates = ["age_at_imaging"]
    if "DLICV" in df.columns:
        numeric_covariates.append("DLICV")

    # Add optional UKB imaging assessment center covariate from the covariate file if available.
    categorical_covariates = ["sex"]
    for c in ["uk_biobank_assessment_centre_f54_2_0"]:
        if c in df.columns:
            df[c] = df[c].astype("category")
            categorical_covariates.append(c)

    for c in brain_feature_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    numeric_cols = numeric_covariates + brain_feature_cols
    categorical_cols = categorical_covariates
    return df, numeric_cols, categorical_cols, numeric_covariates, brain_feature_cols


def make_stratify_vector(df, age_bins=5):
    counts = df["event"].value_counts()
    if len(counts) < 2 or counts.min() < 2:
        return None
    if age_bins and age_bins > 1 and "age_at_imaging" in df.columns:
        try:
            age_bin = pd.qcut(df["age_at_imaging"].rank(method="first"), q=age_bins, labels=False, duplicates="drop")
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
    keep_numeric = []
    dropped = []
    for c in numeric_cols:
        miss = X_train_raw[c].isna().mean()
        if miss <= max_missing:
            keep_numeric.append(c)
        else:
            dropped.append((c, float(miss)))
    if dropped:
        print(f"Dropped {len(dropped)} numeric columns with missingness > {max_missing}.")
        for c, miss in dropped[:30]:
            print(f"  dropped: {c}, missing={miss:.3f}")
        if len(dropped) > 30:
            print("  ...")
    cols = keep_numeric + categorical_cols
    return X_train_raw[cols].copy(), [x[cols].copy() for x in other_raw_list], keep_numeric, categorical_cols, dropped


def brain_feature_mask(feature_names):
    feature_names = np.asarray(feature_names).astype(str)
    return np.char.startswith(feature_names, "num__MUSE_Volume_") | np.char.startswith(feature_names, "num__WMLS_Volume_")


def compute_cindex(y, risk):
    return float(concordance_index_censored(y["event"], y["time"], np.asarray(risk).reshape(-1))[0])


def fit_and_select_coxnet(X_train, y_train, X_val, y_val, feature_names, l1_ratios, n_alphas):
    is_brain = brain_feature_mask(feature_names)
    penalty_factor = np.ones(len(feature_names), dtype=float)
    penalty_factor[~is_brain] = 0.0
    best = {"cindex": -np.inf, "l1_ratio": None, "alpha": None, "coef": None, "used_penalty_factor": True}

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
            warnings.warn("Installed scikit-survival does not support penalty_factor. Covariates will be penalized.")
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
            cindex = compute_cindex(y_val, risk_val)
            if np.isfinite(cindex) and cindex > best["cindex"]:
                best.update({"cindex": float(cindex), "l1_ratio": float(l1_ratio), "alpha": float(alpha), "coef": coefs[:, j].copy()})
        print(f"  best so far: C-index={best['cindex']:.4f}, l1_ratio={best['l1_ratio']}, alpha={best['alpha']}")

    if best["alpha"] is None:
        raise RuntimeError("Failed to select a Coxnet model. Check event counts and input data.")
    return best, penalty_factor


def fit_final_model(X_trainval, y_trainval, best, penalty_factor):
    try:
        model = CoxnetSurvivalAnalysis(
            l1_ratio=best["l1_ratio"],
            alphas=[best["alpha"]],
            penalty_factor=penalty_factor,
            fit_baseline_model=True,
            max_iter=100000,
        )
    except TypeError:
        model = CoxnetSurvivalAnalysis(
            l1_ratio=best["l1_ratio"],
            alphas=[best["alpha"]],
            fit_baseline_model=True,
            max_iter=100000,
        )
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


def get_incremental_model_feature_indices(feature_names):
    feature_names = np.asarray(feature_names).astype(str)
    is_brain = brain_feature_mask(feature_names)
    is_age = feature_names == "num__age_at_imaging"
    is_sex = np.char.startswith(feature_names, "cat__sex")
    return {
        "M0_age_sex": np.where(is_age | is_sex)[0],
        "M1_covariate_baseline": np.where(~is_brain)[0],
        "M2_brain_mri_only": np.where(is_brain)[0],
        "M3_full_covariates_plus_brain_mri": np.arange(len(feature_names)),
    }


def paired_bootstrap_delta_cindex(y, risk_full, risk_baseline, n_boot, random_state):
    rng = np.random.default_rng(random_state)
    event = np.asarray(y["event"]).astype(bool)
    time = np.asarray(y["time"]).astype(float)
    risk_full = np.asarray(risk_full).reshape(-1)
    risk_baseline = np.asarray(risk_baseline).reshape(-1)
    observed_full = compute_cindex(y, risk_full)
    observed_baseline = compute_cindex(y, risk_baseline)
    observed_delta = observed_full - observed_baseline

    boot = []
    n = len(time)
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        if np.sum(event[idx]) < 2:
            continue
        yb = Surv.from_arrays(event=event[idx], time=time[idx])
        try:
            d = compute_cindex(yb, risk_full[idx]) - compute_cindex(yb, risk_baseline[idx])
            if np.isfinite(d):
                boot.append(d)
        except Exception:
            continue
    boot = np.asarray(boot, dtype=float)
    out = {
        "comparison": "M3_full_covariates_plus_brain_mri_vs_M1_covariate_baseline",
        "cindex_full": float(observed_full),
        "cindex_baseline": float(observed_baseline),
        "delta_cindex": float(observed_delta),
        "n_bootstrap_requested": int(n_boot),
        "n_bootstrap_successful": int(boot.size),
    }
    if boot.size == 0:
        out.update({"delta_cindex_ci_lower": np.nan, "delta_cindex_ci_upper": np.nan, "empirical_p_two_sided_delta_not_equal_0": np.nan, "empirical_p_one_sided_delta_le_0": np.nan, "interpretation": "Bootstrap failed."})
        return out
    lo, hi = np.quantile(boot, [0.025, 0.975])
    p_le0 = float(np.mean(boot <= 0.0))
    p_ge0 = float(np.mean(boot >= 0.0))
    p2 = float(min(1.0, 2.0 * min(p_le0, p_ge0)))
    if observed_delta > 0 and lo > 0:
        interp = "Brain MRI improves test-set C-index beyond the covariate baseline."
    elif observed_delta > 0:
        interp = "Brain MRI has positive delta C-index, but the bootstrap CI includes or approaches zero."
    else:
        interp = "No evidence that brain MRI improves test-set C-index beyond the covariate baseline."
    out.update({"delta_cindex_ci_lower": float(lo), "delta_cindex_ci_upper": float(hi), "empirical_p_two_sided_delta_not_equal_0": p2, "empirical_p_one_sided_delta_le_0": p_le0, "interpretation": interp})
    return out


def run_incremental_value_analysis(X_train, X_val, X_test, X_trainval, y_train, y_val, y_test, y_trainval, feature_names, final_model, n_bootstrap, random_state):
    idx = get_incremental_model_feature_indices(feature_names)
    split_data = {"train": (X_train, y_train), "validation": (X_val, y_val), "test": (X_test, y_test), "trainval": (X_trainval, y_trainval)}
    model_labels = {
        "M0_age_sex": "Age + sex",
        "M1_covariate_baseline": "Covariate baseline",
        "M2_brain_mri_only": "Brain MRI only",
        "M3_full_covariates_plus_brain_mri": "Covariates + brain MRI",
    }
    fitted_models = {}
    fitting_methods = {}
    risk_predictions = {s: {} for s in split_data}

    for model_name in ["M0_age_sex", "M1_covariate_baseline", "M2_brain_mri_only"]:
        ind = idx[model_name]
        print(f"Fitting incremental model: {model_name} ({len(ind)} features)")
        model, method = fit_coxph_or_ridge_fallback(X_trainval[:, ind], y_trainval, model_name=model_name)
        fitted_models[model_name] = model
        fitting_methods[model_name] = method
        for split, (X, _) in split_data.items():
            risk_predictions[split][model_name] = predict_risk_score(model, X[:, ind])

    fitting_methods["M3_full_covariates_plus_brain_mri"] = "Selected elastic-net Cox model from main pipeline"
    for split, (X, _) in split_data.items():
        risk_predictions[split]["M3_full_covariates_plus_brain_mri"] = predict_risk_score(final_model, X)

    rows = []
    for model_name in ["M0_age_sex", "M1_covariate_baseline", "M2_brain_mri_only", "M3_full_covariates_plus_brain_mri"]:
        for split, (_, y) in split_data.items():
            rows.append({
                "model": model_name,
                "model_label": model_labels[model_name],
                "split": split,
                "n_features": int(len(idx[model_name])),
                "training_data": "train+validation",
                "fitting_method": fitting_methods[model_name],
                "cindex": compute_cindex(y, risk_predictions[split][model_name]),
                "n": int(len(y)),
                "n_events": int(np.sum(y["event"])),
            })
    model_comparison_df = pd.DataFrame(rows)
    delta_df = pd.DataFrame([paired_bootstrap_delta_cindex(y_test, risk_predictions["test"]["M3_full_covariates_plus_brain_mri"], risk_predictions["test"]["M1_covariate_baseline"], n_bootstrap, random_state)])
    return {"model_comparison_df": model_comparison_df, "delta_cindex_df": delta_df, "fitted_models": fitted_models, "fitting_methods": fitting_methods, "risk_predictions": risk_predictions, "feature_indices": idx}


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


def add_clock_age_and_acceleration(pred_df, covariate_cols, train_mask_col="split"):
    df = pred_df.copy()
    covariate_cols = [c for c in covariate_cols if c in df.columns]
    if "age_at_imaging" not in covariate_cols and "age_at_imaging" in df.columns:
        covariate_cols = ["age_at_imaging"] + covariate_cols
    if not covariate_cols:
        warnings.warn("No residualization covariates available. Acceleration set to missing.")
        df[ACCEL_Z_COL] = df[ACCEL_YEARS_COL] = df[CLOCK_AGE_COL] = np.nan
        return df, None

    train = df.loc[df[train_mask_col] == "train"].copy()
    if train.shape[0] < 10:
        warnings.warn("Too few training samples for clock-acceleration residualization.")
        df[ACCEL_Z_COL] = df[ACCEL_YEARS_COL] = df[CLOCK_AGE_COL] = np.nan
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

    X_train = prep.fit_transform(X_train_raw)
    X_all = prep.transform(X_all_raw)
    lr = LinearRegression().fit(X_train, train[RISK_COL].values)
    expected = lr.predict(X_all)
    resid_raw = df[RISK_COL].values - expected
    train_index = df[train_mask_col].values == "train"
    resid_mean = float(np.nanmean(resid_raw[train_index]))
    resid_sd = float(np.nanstd(resid_raw[train_index]))
    resid = resid_raw - resid_mean
    df[ACCEL_Z_COL] = resid / resid_sd if resid_sd > 0 else np.nan

    feat_names = []
    feat_names.extend([f"num__{c}" for c in numeric_covs])
    if categorical_covs:
        ohe = prep.named_transformers_["cat"].named_steps["onehot"]
        try:
            cat_names = ohe.get_feature_names_out(categorical_covs)
        except AttributeError:
            cat_names = ohe.get_feature_names(categorical_covs)
        feat_names.extend([f"cat__{c}" for c in cat_names])

    beta_age = np.nan
    if "num__age_at_imaging" in feat_names:
        beta_age = float(lr.coef_[feat_names.index("num__age_at_imaging")])
    if np.isfinite(beta_age) and abs(beta_age) > 1e-8:
        df[ACCEL_YEARS_COL] = resid / beta_age
        df[CLOCK_AGE_COL] = df["age_at_imaging"] + df[ACCEL_YEARS_COL]
    else:
        warnings.warn("Adjusted age coefficient is near zero or unavailable. Year-scale acceleration set to missing.")
        df[ACCEL_YEARS_COL] = np.nan
        df[CLOCK_AGE_COL] = np.nan

    info = {
        "residualization_covariates": covariate_cols,
        "numeric_residualization_covariates": numeric_covs,
        "categorical_residualization_covariates": categorical_covs,
        "risk_score_covariate_model_intercept": float(lr.intercept_),
        "risk_score_covariate_model_coef": {k: float(v) for k, v in zip(feat_names, lr.coef_)},
        "adjusted_age_coefficient_risk_score_per_year": float(beta_age) if np.isfinite(beta_age) else None,
        "risk_score_residual_mean_train": resid_mean,
        "risk_score_residual_sd_train": resid_sd,
        "note": "Clock acceleration is the residual of the Cox risk score after adjustment for retained non-brain covariates. The z-score is recommended for downstream analyses.",
    }
    return df, info


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    admin_censor_date = pd.to_datetime(args.admin_censor_date)
    l1_ratios = tuple(float(x) for x in args.l1_ratios.split(","))

    print("Loading dementia/assessment data...")
    dementia = load_dementia_data(args.dementia_xlsx, args.id_match_csv)
    dementia["admin_censor_date"] = admin_censor_date

    print("Loading brain MRI data...")
    brain = load_brain_data(args.brain_tsv, args.imaging_session_id)

    print("Loading optional covariates...")
    cov = load_covariates(args.covariate_csv)

    print("Merging data...")
    df = brain.merge(dementia, on="participant_id", how="inner")
    if cov is not None:
        df = df.merge(cov, on="participant_id", how="left", suffixes=("", "_cov"))

    print("Constructing prospective dementia survival outcome...")
    df = construct_survival_dataset(df)
    df = df.loc[df["time_days"] >= args.min_followup_days].copy()

    print("Adding brain-specific covariates...")
    df = add_basic_covariates(df)
    brain_feature_cols = infer_brain_feature_columns(df)
    print(f"Found {len(brain_feature_cols)} brain MRI features.")

    df, numeric_cols, categorical_cols, numeric_covariates, brain_feature_cols = build_design_matrix(df, brain_feature_cols)
    required_model_cols = ["participant_id", "time_years", "event", "age_at_imaging", "sex"]
    df = df.dropna(subset=required_model_cols).copy()

    print("Final prospective imaging dementia dataset:")
    print(f"  N = {df.shape[0]}")
    print(f"  Incident dementia events after imaging = {int(df['event'].sum())}")
    print(f"  Censored = {int((~df['event']).sum())}")
    print(f"  Median follow-up years = {df['time_years'].median():.2f}")
    if df["event"].sum() < 20:
        warnings.warn("Very few dementia events. The fitted brain MRI dementia clock may be unstable.")

    keep_dataset_cols = ["participant_id", "baseline_date", "imaging_date", "dementia_date", "death_date", "admin_censor_date", "censor_date", "end_date", "time_days", "time_years", "event", "age_at_imaging", "sex"]
    if "DLICV" in df.columns:
        keep_dataset_cols.append("DLICV")
    df[keep_dataset_cols].to_csv(outdir / f"{PREFIX}_survival_dataset.tsv", sep="\t", index=False)

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
    X_train_raw, other, numeric_cols_kept, categorical_cols_kept, dropped_numeric = drop_high_missing_features(
        X_train_raw,
        [X_val_raw, X_test_raw, X_trainval_raw],
        numeric_cols,
        categorical_cols,
        args.max_feature_missing,
    )
    X_val_raw, X_test_raw, X_trainval_raw = other

    residualization_covariates = [c for c in (numeric_cols_kept + categorical_cols_kept) if c not in brain_feature_cols]
    print("Residualizing clock acceleration on retained non-brain covariates:")
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
    best, penalty_factor = fit_and_select_coxnet(X_train, y_train, X_val, y_val, feature_names, l1_ratios, args.n_alphas)
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

    print("Running incremental-value analysis: does brain MRI add value beyond covariates?")
    inc = run_incremental_value_analysis(
        X_train,
        X_val,
        X_test,
        X_trainval,
        y_train,
        y_val,
        y_test,
        y_trainval,
        feature_names,
        final_model,
        args.n_bootstrap_incremental,
        args.random_state,
    )
    model_comparison_df = inc["model_comparison_df"]
    delta_cindex_df = inc["delta_cindex_df"]
    model_comparison_df.to_csv(outdir / f"{PREFIX}_model_comparison.tsv", sep="\t", index=False)
    delta_cindex_df.to_csv(outdir / f"{PREFIX}_incremental_value_delta_cindex.tsv", sep="\t", index=False)
    print(delta_cindex_df.to_string(index=False))

    def make_pred_frame(part, split, risk):
        base = ["participant_id", "imaging_date", "dementia_date", "death_date", "admin_censor_date", "censor_date", "end_date", "time_years", "event", "age_at_imaging", "sex"]
        if "DLICV" in part.columns:
            base.append("DLICV")
        extra = [c for c in residualization_covariates if c in part.columns and c not in base]
        out = part[base + extra].copy()
        out["split"] = split
        out[RISK_COL] = risk
        return out

    pred_train = make_pred_frame(df_train, "train", risk_train)
    pred_val = make_pred_frame(df_val, "validation", risk_val)
    pred_test = make_pred_frame(df_test, "test", risk_test)
    comparison_cols = {
        "M0_age_sex": "risk_score_M0_age_sex",
        "M1_covariate_baseline": "risk_score_M1_covariate_baseline",
        "M2_brain_mri_only": "risk_score_M2_brain_mri_only",
        "M3_full_covariates_plus_brain_mri": "risk_score_M3_full_covariates_plus_brain_mri",
    }
    for model_name, col in comparison_cols.items():
        pred_train[col] = inc["risk_predictions"]["train"][model_name]
        pred_val[col] = inc["risk_predictions"]["validation"][model_name]
        pred_test[col] = inc["risk_predictions"]["test"][model_name]

    risk_times = [5.0, 10.0, 15.0]
    pred_train = pd.concat([pred_train.reset_index(drop=True), predict_absolute_risk(final_model, X_train, risk_times)], axis=1)
    pred_val = pd.concat([pred_val.reset_index(drop=True), predict_absolute_risk(final_model, X_val, risk_times)], axis=1)
    pred_test = pd.concat([pred_test.reset_index(drop=True), predict_absolute_risk(final_model, X_test, risk_times)], axis=1)
    pred_all = pd.concat([pred_train, pred_val, pred_test], ignore_index=True)

    print("Adding residualized dementia-clock acceleration...")
    pred_all, clock_transform_info = add_clock_age_and_acceleration(pred_all, residualization_covariates)
    pred_all.to_csv(outdir / f"{PREFIX}_predictions.tsv", sep="\t", index=False)
    pred_test.to_csv(outdir / f"{PREFIX}_test_predictions.tsv", sep="\t", index=False)

    print("Saving coefficients...")
    coef = np.asarray(final_model.coef_).reshape(-1)
    is_brain = brain_feature_mask(feature_names)
    coef_df = pd.DataFrame({
        "feature": feature_names,
        "coefficient": coef,
        "abs_coefficient": np.abs(coef),
        "penalty_factor": penalty_factor,
        "is_nonzero": coef != 0,
        "is_brain_mri_feature": is_brain,
    }).sort_values("abs_coefficient", ascending=False)
    coef_df.to_csv(outdir / f"{PREFIX}_coefficients.tsv", sep="\t", index=False)
    nonzero_coef_df = coef_df.loc[coef_df["is_nonzero"]].copy()
    nonzero_coef_df.to_csv(outdir / f"{PREFIX}_nonzero_coefficients.tsv", sep="\t", index=False)

    model_comparison_records = json.loads(model_comparison_df.to_json(orient="records"))
    delta_cindex_records = json.loads(delta_cindex_df.to_json(orient="records"))
    model_bundle = {
        "preprocessor": preprocessor,
        "model": final_model,
        "feature_names": feature_names,
        "numeric_cols_kept": numeric_cols_kept,
        "categorical_cols_kept": categorical_cols_kept,
        "brain_feature_cols": brain_feature_cols,
        "dropped_numeric": dropped_numeric,
        "residualization_covariates": residualization_covariates,
        "best": best,
        "penalty_factor": penalty_factor,
        "clock_transform_info": clock_transform_info,
        "incremental_value_model_comparison": model_comparison_records,
        "incremental_value_delta_cindex": delta_cindex_records,
        "incremental_value_fitting_methods": inc["fitting_methods"],
        "brain_tsv_input": args.brain_tsv,
        "admin_censor_date": str(admin_censor_date.date()),
    }
    joblib.dump(model_bundle, outdir / f"{PREFIX}_model.joblib")

    delta_row = delta_cindex_df.iloc[0]
    performance = {
        "n_total": int(df.shape[0]),
        "n_events_total": int(df["event"].sum()),
        "n_censored_total": int((~df["event"]).sum()),
        "median_followup_years": float(df["time_years"].median()),
        "n_train": int(df_train.shape[0]),
        "n_events_train": int(df_train["event"].sum()),
        "n_validation": int(df_val.shape[0]),
        "n_events_validation": int(df_val["event"].sum()),
        "n_test": int(df_test.shape[0]),
        "n_events_test": int(df_test["event"].sum()),
        "cindex_train": float(cindex_train),
        "cindex_validation": float(cindex_val),
        "cindex_trainval": float(cindex_trainval),
        "cindex_test": float(cindex_test),
        "incremental_value_model_comparison": model_comparison_records,
        "incremental_value_delta_cindex": delta_cindex_records,
        "cindex_test_M1_covariate_baseline": float(model_comparison_df.loc[(model_comparison_df["model"] == "M1_covariate_baseline") & (model_comparison_df["split"] == "test"), "cindex"].iloc[0]),
        "cindex_test_M3_full_covariates_plus_brain_mri": float(model_comparison_df.loc[(model_comparison_df["model"] == "M3_full_covariates_plus_brain_mri") & (model_comparison_df["split"] == "test"), "cindex"].iloc[0]),
        "delta_cindex_test_M3_vs_M1": float(delta_row["delta_cindex"]),
        "delta_cindex_test_M3_vs_M1_ci_lower": float(delta_row["delta_cindex_ci_lower"]),
        "delta_cindex_test_M3_vs_M1_ci_upper": float(delta_row["delta_cindex_ci_upper"]),
        "delta_cindex_test_M3_vs_M1_p_two_sided": float(delta_row["empirical_p_two_sided_delta_not_equal_0"]),
        "best_l1_ratio": float(best["l1_ratio"]),
        "best_alpha": float(best["alpha"]),
        "best_validation_cindex_during_tuning": float(best["cindex"]),
        "used_penalty_factor": bool(best["used_penalty_factor"]),
        "n_original_brain_features": int(len(brain_feature_cols)),
        "n_numeric_cols_kept": int(len(numeric_cols_kept)),
        "n_categorical_cols_kept": int(len(categorical_cols_kept)),
        "n_nonzero_coefficients": int(nonzero_coef_df.shape[0]),
        "n_residualization_covariates": int(len(residualization_covariates)),
        "residualization_covariates": residualization_covariates,
        "admin_censor_date": str(admin_censor_date.date()),
        "time_zero": "UKB imaging assessment date, field 53-2.0",
        "event_date": "UKB dementia date, field 42018-0.0",
        "censoring": "earlier of death date field 40000-0.0 and endpoint-specific administrative censor date",
        "note": "Primary score is brain_mri_dementia_risk_score from elastic-net Cox. Acceleration z-score is residualized on retained non-brain covariates and is recommended for downstream analyses.",
    }
    with open(outdir / f"{PREFIX}_performance.json", "w") as f:
        json.dump(performance, f, indent=2)

    print("Done.")
    print(f"Outputs written to: {outdir}")
    print("Main output files:")
    for name in ["survival_dataset.tsv", "predictions.tsv", "test_predictions.tsv", "coefficients.tsv", "nonzero_coefficients.tsv", "model_comparison.tsv", "incremental_value_delta_cindex.tsv", "model.joblib", "performance.json"]:
        print(f"  {outdir / f'{PREFIX}_{name}'}")


if __name__ == "__main__":
    main()
