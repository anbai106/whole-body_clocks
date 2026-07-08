#!/usr/bin/env python3
"""
Build ADNI brain MRI AD L'EPOCH using elastic-net Cox survival modeling.

Time zero:
  ADNI baseline visit, preferably Visit_Code == bl.

Eligible baseline diagnosis:
  CN.

Event:
  first follow-up diagnosis of MCI or AD.

Censoring:
  last available follow-up visit if no MCI/AD conversion.

Features:
  Hard-coded MUSE gray-matter ROI volumes, matched to the UKB MUSE GM ROI set.

Model:
  Elastic-net Cox survival model using baseline MUSE GM ROI volumes plus covariates.

Important constraint:
  Model selection requires at least a minimum number of nonzero MUSE GM ROI coefficients,
  so the final L'EPOCH remains brain-imaging related and is not purely covariate driven.

Primary output:
  adni_brain_mri_ad_lepoch_risk_score

Recommended downstream phenotype:
  adni_brain_mri_ad_lepoch_acceleration_z
"""

import argparse
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
        "  conda install -c conda-forge scikit-survival\n"
        "or activate your existing survival_clock environment."
    ) from e


# ============================================================
# Constants
# ============================================================

RISK_COL = "adni_brain_mri_ad_lepoch_risk_score"
ACCEL_Z_COL = "adni_brain_mri_ad_lepoch_acceleration_z"
ACCEL_YEARS_COL = "adni_brain_mri_ad_lepoch_acceleration_years"
CLOCK_AGE_COL = "adni_brain_mri_ad_lepoch_clock_age_years"

# Hard-coded MUSE GM ROIs requested by the user.
MUSE_GM_ROIS = [
    "MUSE_Volume_23", "MUSE_Volume_30", "MUSE_Volume_31",
    "MUSE_Volume_32", "MUSE_Volume_36", "MUSE_Volume_37",
    "MUSE_Volume_38", "MUSE_Volume_39", "MUSE_Volume_47",
    "MUSE_Volume_48", "MUSE_Volume_55", "MUSE_Volume_56",
    "MUSE_Volume_57", "MUSE_Volume_58", "MUSE_Volume_59",
    "MUSE_Volume_60", "MUSE_Volume_71", "MUSE_Volume_72",
    "MUSE_Volume_73", "MUSE_Volume_75", "MUSE_Volume_76",
    "MUSE_Volume_100", "MUSE_Volume_101", "MUSE_Volume_102",
    "MUSE_Volume_103", "MUSE_Volume_104", "MUSE_Volume_105",
    "MUSE_Volume_106", "MUSE_Volume_107", "MUSE_Volume_108",
    "MUSE_Volume_109", "MUSE_Volume_112", "MUSE_Volume_113",
    "MUSE_Volume_114", "MUSE_Volume_115", "MUSE_Volume_116",
    "MUSE_Volume_117", "MUSE_Volume_118", "MUSE_Volume_119",
    "MUSE_Volume_120", "MUSE_Volume_121", "MUSE_Volume_122",
    "MUSE_Volume_123", "MUSE_Volume_124", "MUSE_Volume_125",
    "MUSE_Volume_128", "MUSE_Volume_129", "MUSE_Volume_132",
    "MUSE_Volume_133", "MUSE_Volume_134", "MUSE_Volume_135",
    "MUSE_Volume_136", "MUSE_Volume_137", "MUSE_Volume_138",
    "MUSE_Volume_139", "MUSE_Volume_140", "MUSE_Volume_141",
    "MUSE_Volume_142", "MUSE_Volume_143", "MUSE_Volume_144",
    "MUSE_Volume_145", "MUSE_Volume_146", "MUSE_Volume_147",
    "MUSE_Volume_148", "MUSE_Volume_149", "MUSE_Volume_150",
    "MUSE_Volume_151", "MUSE_Volume_152", "MUSE_Volume_153",
    "MUSE_Volume_154", "MUSE_Volume_155", "MUSE_Volume_156",
    "MUSE_Volume_157", "MUSE_Volume_160", "MUSE_Volume_161",
    "MUSE_Volume_162", "MUSE_Volume_163", "MUSE_Volume_164",
    "MUSE_Volume_165", "MUSE_Volume_166", "MUSE_Volume_167",
    "MUSE_Volume_168", "MUSE_Volume_169", "MUSE_Volume_170",
    "MUSE_Volume_171", "MUSE_Volume_172", "MUSE_Volume_173",
    "MUSE_Volume_174", "MUSE_Volume_175", "MUSE_Volume_176",
    "MUSE_Volume_177", "MUSE_Volume_178", "MUSE_Volume_179",
    "MUSE_Volume_180", "MUSE_Volume_181", "MUSE_Volume_182",
    "MUSE_Volume_183", "MUSE_Volume_184", "MUSE_Volume_185",
    "MUSE_Volume_186", "MUSE_Volume_187", "MUSE_Volume_190",
    "MUSE_Volume_191", "MUSE_Volume_192", "MUSE_Volume_193",
    "MUSE_Volume_194", "MUSE_Volume_195", "MUSE_Volume_196",
    "MUSE_Volume_197", "MUSE_Volume_198", "MUSE_Volume_199",
    "MUSE_Volume_200", "MUSE_Volume_201", "MUSE_Volume_202",
    "MUSE_Volume_203", "MUSE_Volume_204", "MUSE_Volume_205",
    "MUSE_Volume_206", "MUSE_Volume_207"
]


# ============================================================
# Argument parsing
# ============================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Build ADNI brain MRI AD L'EPOCH using baseline MUSE GM ROIs and Coxnet survival modeling."
    )

    parser.add_argument("--input-file", required=True, help="ADNI longitudinal MUSE file, TSV or CSV.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--prefix", default="adni_brain_mri_ad_lepoch")

    parser.add_argument("--id-col", default="PTID")
    parser.add_argument("--visit-col", default="Visit_Code")
    parser.add_argument("--date-col", default="Date")
    parser.add_argument("--dx-col", default="DX_Binary")

    parser.add_argument("--baseline-dx", default="CN")
    parser.add_argument("--event-dx", default="MCI,AD")

    parser.add_argument("--covariates", default="Age,Sex,DLICV,SITE")

    parser.add_argument("--test-size", type=float, default=0.20)
    parser.add_argument("--validation-size", type=float, default=0.20)
    parser.add_argument("--random-state", type=int, default=20260707)
    parser.add_argument("--stratify-age-bins", type=int, default=5)

    parser.add_argument("--max-feature-missing", type=float, default=0.30)
    parser.add_argument("--min-followup-days", type=int, default=1)

    parser.add_argument("--l1-ratios", default="0.1,0.25,0.5,0.75,1.0")
    parser.add_argument("--n-alphas", type=int, default=120)
    parser.add_argument("--alpha-min-ratio", type=float, default=0.001)

    parser.add_argument(
        "--min-nonzero-brain-features",
        type=int,
        default=5,
        help="Require at least this many nonzero MUSE GM ROI coefficients during model selection."
    )

    parser.add_argument("--n-bootstrap-incremental", type=int, default=1000)

    return parser.parse_args()


# ============================================================
# Utility functions
# ============================================================

def log(msg):
    print(msg, flush=True)


def read_table(path):
    path = Path(path)
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path, low_memory=False)
    return pd.read_csv(path, sep="\t", low_memory=False)


def parse_list_arg(x):
    if x is None or str(x).strip() == "":
        return []
    return [v.strip() for v in str(x).split(",") if v.strip()]


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


def parse_date_series(s):
    return pd.to_datetime(s, errors="coerce")


def visit_code_to_month(x):
    if pd.isna(x):
        return np.nan

    x = str(x).strip().lower()

    if x in ["bl", "base", "baseline", "m00", "m0", "screen", "screening", "sc"]:
        return 0.0

    m = re.search(r"m(?:onth)?\s*0*([0-9]+)", x)
    if m:
        return float(m.group(1))

    m = re.search(r"([0-9]+)", x)
    if m:
        return float(m.group(1))

    return np.nan


def make_onehot_encoder():
    try:
        return OneHotEncoder(handle_unknown="ignore", sparse_output=False)
    except TypeError:
        return OneHotEncoder(handle_unknown="ignore", sparse=False)


def compute_cindex(y, risk):
    return float(concordance_index_censored(
        y["event"],
        y["time"],
        np.asarray(risk).reshape(-1)
    )[0])


def brain_feature_mask(feature_names):
    feature_names = np.asarray(feature_names).astype(str)
    return np.char.startswith(feature_names, "num__MUSE_Volume_")


def make_stratify_vector(df, age_col="Age", event_col="event", age_bins=5):
    if event_col not in df.columns:
        return None

    counts = df[event_col].value_counts()
    if len(counts) < 2 or counts.min() < 2:
        return None

    if age_col in df.columns and age_bins and age_bins > 1:
        try:
            age_numeric = pd.to_numeric(df[age_col], errors="coerce")
            age_bin = pd.qcut(
                age_numeric.rank(method="first"),
                q=age_bins,
                labels=False,
                duplicates="drop"
            )
            lab = df[event_col].astype(int).astype(str) + "_age" + age_bin.astype(str)
            if lab.value_counts().min() >= 2:
                return lab
        except Exception:
            pass

    return df[event_col].astype(int)


def fit_coxph_or_ridge_fallback(X, y, model_name="cox_model"):
    try:
        try:
            model = CoxPHSurvivalAnalysis(alpha=0.0, ties="breslow")
        except TypeError:
            model = CoxPHSurvivalAnalysis(alpha=0.0)
        model.fit(X, y)
        return model, "CoxPHSurvivalAnalysis(alpha=0.0)"
    except Exception as exc:
        warnings.warn(
            f"{model_name}: unpenalized CoxPH failed ({exc}). "
            "Falling back to weakly penalized Coxnet."
        )

    last = None
    for alpha in [1e-6, 1e-5, 1e-4, 1e-3, 1e-2]:
        try:
            model = CoxnetSurvivalAnalysis(
                l1_ratio=0.01,
                alphas=[alpha],
                fit_baseline_model=True,
                max_iter=100000
            )
            model.fit(X, y)
            return model, f"CoxnetSurvivalAnalysis(l1_ratio=0.01, alpha={alpha})"
        except Exception as exc:
            last = exc

    raise RuntimeError(f"{model_name}: all Cox fitting attempts failed. Last error: {last}")


def predict_risk_score(model, X):
    return np.asarray(model.predict(X)).reshape(-1)


# ============================================================
# ADNI survival dataset construction
# ============================================================

def construct_adni_survival_dataset(
    df,
    id_col,
    visit_col,
    date_col,
    dx_col,
    baseline_dx,
    event_dx_set,
    min_followup_days
):
    d = df.copy()

    d["_dx_norm"] = d[dx_col].apply(normalize_dx)
    d["_date"] = parse_date_series(d[date_col])
    d["_visit_month"] = d[visit_col].apply(visit_code_to_month)

    d["_date_sort"] = d["_date"].map(lambda x: x.toordinal() if pd.notna(x) else np.nan)
    d["_visit_sort"] = d["_visit_month"] * 30.4375
    d["_sort_key"] = d["_date_sort"].fillna(d["_visit_sort"])

    d = d.loc[d[id_col].notna()].copy()
    d = d.loc[d["_sort_key"].notna()].copy()
    d = d.sort_values([id_col, "_sort_key"], kind="mergesort")

    baseline_rows = []

    for pid, g in d.groupby(id_col, sort=False):
        g = g.copy().sort_values("_sort_key", kind="mergesort")

        bl_mask = g[visit_col].astype(str).str.lower().isin(
            ["bl", "base", "baseline", "m00", "m0"]
        )

        if bl_mask.any():
            baseline = g.loc[bl_mask].iloc[0]
        else:
            baseline = g.iloc[0]

        baseline_dx_value = normalize_dx(baseline[dx_col])
        if baseline_dx_value != baseline_dx:
            continue

        base_sort = baseline["_sort_key"]
        base_date = baseline["_date"]
        base_month = baseline["_visit_month"]

        follow = g.loc[g["_sort_key"] > base_sort].copy()
        if follow.empty:
            continue

        follow["_is_event"] = follow["_dx_norm"].isin(event_dx_set)

        if follow["_is_event"].any():
            end_row = follow.loc[follow["_is_event"]].iloc[0]
            event = True
        else:
            end_row = follow.iloc[-1]
            event = False

        if pd.notna(base_date) and pd.notna(end_row["_date"]):
            time_days = (end_row["_date"] - base_date).days
        elif pd.notna(base_month) and pd.notna(end_row["_visit_month"]):
            time_days = int(round((end_row["_visit_month"] - base_month) * 30.4375))
        else:
            time_days = np.nan

        if not np.isfinite(time_days) or time_days < min_followup_days:
            continue

        b = baseline.copy()
        b["baseline_dx"] = baseline_dx_value
        b["event"] = bool(event)
        b["time_days"] = float(time_days)
        b["time_years"] = float(time_days) / 365.25
        b["event_or_censor_dx"] = normalize_dx(end_row[dx_col])
        b["event_or_censor_visit"] = end_row[visit_col]
        b["event_or_censor_date"] = end_row["_date"]

        baseline_rows.append(b)

    if len(baseline_rows) == 0:
        raise ValueError(
            "No eligible baseline CN participants with follow-up were found. "
            "Check PTID, Visit_Code, Date, and DX_Binary."
        )

    out = pd.DataFrame(baseline_rows).reset_index(drop=True)

    return out


# ============================================================
# Feature selection and preprocessing
# ============================================================

def select_hardcoded_muse_gm_rois(df):
    available = [c for c in MUSE_GM_ROIS if c in df.columns]
    missing = [c for c in MUSE_GM_ROIS if c not in df.columns]

    if missing:
        warnings.warn(
            f"{len(missing)} requested MUSE GM ROI columns are missing and will be skipped. "
            f"First missing columns: {missing[:10]}"
        )

    if len(available) == 0:
        raise ValueError(
            "None of the hard-coded MUSE GM ROI columns were found in the ADNI file."
        )

    return available, missing


def filter_features_by_missingness_and_variance(df, feature_cols, max_missing):
    keep = []
    rows = []

    for c in feature_cols:
        x = pd.to_numeric(df[c], errors="coerce")
        missing_rate = float(x.isna().mean())
        variance = float(x.var(skipna=True)) if x.notna().sum() > 1 else np.nan

        keep_flag = (
            missing_rate <= max_missing
            and np.isfinite(variance)
            and variance > 0
        )

        rows.append({
            "roi_feature": c,
            "missing_rate": missing_rate,
            "variance": variance,
            "kept": keep_flag
        })

        if keep_flag:
            keep.append(c)

    if len(keep) == 0:
        raise ValueError(
            "No MUSE GM ROI features passed missingness/variance QC. "
            "Consider relaxing --max-feature-missing."
        )

    return keep, pd.DataFrame(rows)


def build_design_matrix_columns(df, roi_cols, covariates):
    numeric_covariates = []
    categorical_covariates = []

    for c in covariates:
        if c not in df.columns:
            warnings.warn(f"Covariate {c} is missing and will be skipped.")
            continue

        if c in ["Sex", "SITE"]:
            categorical_covariates.append(c)
            continue

        x = pd.to_numeric(df[c], errors="coerce")
        non_missing_original = df[c].notna().sum()
        non_missing_numeric = x.notna().sum()

        if non_missing_original > 0 and non_missing_numeric / max(non_missing_original, 1) > 0.90:
            df[c] = x
            numeric_covariates.append(c)
        else:
            categorical_covariates.append(c)

    for c in roi_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    numeric_cols = numeric_covariates + roi_cols
    categorical_cols = categorical_covariates

    return df, numeric_cols, categorical_cols, numeric_covariates, categorical_covariates


def make_preprocessor(numeric_cols, categorical_cols):
    transformers = []

    if len(numeric_cols) > 0:
        num = Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler())
        ])
        transformers.append(("num", num, numeric_cols))

    if len(categorical_cols) > 0:
        cat = Pipeline([
            ("imputer", SimpleImputer(strategy="most_frequent")),
            ("onehot", make_onehot_encoder())
        ])
        transformers.append(("cat", cat, categorical_cols))

    return ColumnTransformer(transformers, remainder="drop")


def get_feature_names(preprocessor):
    names = []

    for name, transformer, cols in preprocessor.transformers_:
        if name == "num":
            names.extend([f"num__{c}" for c in cols])
        elif name == "cat":
            cat_pipe = preprocessor.named_transformers_[name]
            ohe = cat_pipe.named_steps["onehot"]
            try:
                cat_names = ohe.get_feature_names_out(cols)
            except AttributeError:
                cat_names = ohe.get_feature_names(cols)
            names.extend([f"cat__{c}" for c in cat_names])

    return np.asarray(names)


# ============================================================
# Coxnet model fitting with nonzero brain feature requirement
# ============================================================

def fit_and_select_coxnet(
    X_train,
    y_train,
    X_val,
    y_val,
    feature_names,
    l1_ratios,
    n_alphas,
    alpha_min_ratio,
    min_nonzero_brain_features
):
    is_brain = brain_feature_mask(feature_names)

    if is_brain.sum() == 0:
        raise ValueError("No brain MRI features detected in transformed feature matrix.")

    penalty_factor = np.ones(len(feature_names), dtype=float)
    penalty_factor[~is_brain] = 0.0

    best = {
        "cindex": -np.inf,
        "l1_ratio": None,
        "alpha": None,
        "coef": None,
        "n_nonzero_brain_features": 0,
        "n_nonzero_total_features": 0,
        "used_penalty_factor": True,
        "selection_rule": f"validation C-index among models with >= {min_nonzero_brain_features} nonzero brain features"
    }

    fallback = {
        "cindex": -np.inf,
        "l1_ratio": None,
        "alpha": None,
        "coef": None,
        "n_nonzero_brain_features": 0,
        "n_nonzero_total_features": 0,
        "used_penalty_factor": True,
        "selection_rule": "fallback: best validation C-index among models with the maximum available number of nonzero brain features"
    }

    for l1_ratio in l1_ratios:
        log(f"Fitting Coxnet path for l1_ratio={l1_ratio}")

        try:
            model = CoxnetSurvivalAnalysis(
                l1_ratio=l1_ratio,
                n_alphas=n_alphas,
                alpha_min_ratio=alpha_min_ratio,
                penalty_factor=penalty_factor,
                fit_baseline_model=False,
                max_iter=100000
            )
            used_penalty_factor = True
        except TypeError:
            warnings.warn(
                "Installed scikit-survival does not support penalty_factor. "
                "Covariates will also be penalized."
            )
            model = CoxnetSurvivalAnalysis(
                l1_ratio=l1_ratio,
                n_alphas=n_alphas,
                alpha_min_ratio=alpha_min_ratio,
                fit_baseline_model=False,
                max_iter=100000
            )
            used_penalty_factor = False

        model.fit(X_train, y_train)

        coefs = model.coef_
        if coefs.ndim == 1:
            coefs = coefs[:, None]

        for j, alpha in enumerate(model.alphas_):
            coef_j = coefs[:, j]
            nonzero = np.abs(coef_j) > 1e-8
            n_nonzero_brain = int(np.sum(nonzero & is_brain))
            n_nonzero_total = int(np.sum(nonzero))

            risk_val = np.dot(X_val, coef_j)
            cindex = compute_cindex(y_val, risk_val)

            if not np.isfinite(cindex):
                continue

            # Fallback tracks best model among those with the greatest number of brain features.
            if (
                n_nonzero_brain > fallback["n_nonzero_brain_features"]
                or (
                    n_nonzero_brain == fallback["n_nonzero_brain_features"]
                    and cindex > fallback["cindex"]
                )
            ):
                fallback.update({
                    "cindex": float(cindex),
                    "l1_ratio": float(l1_ratio),
                    "alpha": float(alpha),
                    "coef": coef_j.copy(),
                    "n_nonzero_brain_features": n_nonzero_brain,
                    "n_nonzero_total_features": n_nonzero_total,
                    "used_penalty_factor": bool(used_penalty_factor)
                })

            # Primary selection requires a minimum number of nonzero brain features.
            if n_nonzero_brain < min_nonzero_brain_features:
                continue

            if cindex > best["cindex"]:
                best.update({
                    "cindex": float(cindex),
                    "l1_ratio": float(l1_ratio),
                    "alpha": float(alpha),
                    "coef": coef_j.copy(),
                    "n_nonzero_brain_features": n_nonzero_brain,
                    "n_nonzero_total_features": n_nonzero_total,
                    "used_penalty_factor": bool(used_penalty_factor)
                })

        log(
            f"  best valid so far: C-index={best['cindex']:.4f}, "
            f"l1_ratio={best['l1_ratio']}, alpha={best['alpha']}, "
            f"nonzero brain={best['n_nonzero_brain_features']}"
        )
        log(
            f"  fallback so far: C-index={fallback['cindex']:.4f}, "
            f"l1_ratio={fallback['l1_ratio']}, alpha={fallback['alpha']}, "
            f"nonzero brain={fallback['n_nonzero_brain_features']}"
        )

    if best["alpha"] is None:
        warnings.warn(
            f"No candidate model retained >= {min_nonzero_brain_features} nonzero brain features. "
            "Using fallback model with the maximum available number of nonzero brain features."
        )
        best = fallback

    if best["alpha"] is None or best["n_nonzero_brain_features"] == 0:
        raise RuntimeError(
            "Failed to select a brain-imaging-related Coxnet model. "
            "All candidate models had zero nonzero MUSE GM ROI coefficients. "
            "Consider reducing --min-nonzero-brain-features, lowering --alpha-min-ratio, "
            "or checking whether MUSE features have valid variance."
        )

    return best, penalty_factor


def fit_final_model(X_trainval, y_trainval, best, penalty_factor):
    try:
        model = CoxnetSurvivalAnalysis(
            l1_ratio=best["l1_ratio"],
            alphas=[best["alpha"]],
            penalty_factor=penalty_factor,
            fit_baseline_model=True,
            max_iter=100000
        )
    except TypeError:
        model = CoxnetSurvivalAnalysis(
            l1_ratio=best["l1_ratio"],
            alphas=[best["alpha"]],
            fit_baseline_model=True,
            max_iter=100000
        )

    model.fit(X_trainval, y_trainval)
    return model


# ============================================================
# Incremental value analysis
# ============================================================

def get_incremental_model_feature_indices(feature_names):
    feature_names = np.asarray(feature_names).astype(str)

    is_brain = brain_feature_mask(feature_names)
    is_age = feature_names == "num__Age"
    is_sex = np.char.startswith(feature_names, "cat__Sex")

    return {
        "M0_age_sex": np.where(is_age | is_sex)[0],
        "M1_covariate_baseline": np.where(~is_brain)[0],
        "M2_brain_mri_only": np.where(is_brain)[0],
        "M3_full_covariates_plus_brain_mri": np.arange(len(feature_names))
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
        "n_bootstrap_successful": int(boot.size)
    }

    if boot.size == 0:
        out.update({
            "delta_cindex_ci_lower": np.nan,
            "delta_cindex_ci_upper": np.nan,
            "empirical_p_two_sided_delta_not_equal_0": np.nan,
            "empirical_p_one_sided_delta_le_0": np.nan,
            "interpretation": "Bootstrap failed."
        })
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

    out.update({
        "delta_cindex_ci_lower": float(lo),
        "delta_cindex_ci_upper": float(hi),
        "empirical_p_two_sided_delta_not_equal_0": p2,
        "empirical_p_one_sided_delta_le_0": p_le0,
        "interpretation": interp
    })

    return out


def run_incremental_value_analysis(
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
    n_bootstrap,
    random_state
):
    idx = get_incremental_model_feature_indices(feature_names)

    split_data = {
        "train": (X_train, y_train),
        "validation": (X_val, y_val),
        "test": (X_test, y_test),
        "trainval": (X_trainval, y_trainval)
    }

    model_labels = {
        "M0_age_sex": "Age + sex",
        "M1_covariate_baseline": "Covariate baseline",
        "M2_brain_mri_only": "Brain MRI only",
        "M3_full_covariates_plus_brain_mri": "Covariates + brain MRI"
    }

    fitted_models = {}
    fitting_methods = {}
    risk_predictions = {s: {} for s in split_data}

    for model_name in ["M0_age_sex", "M1_covariate_baseline", "M2_brain_mri_only"]:
        ind = idx[model_name]

        if len(ind) == 0:
            warnings.warn(f"Skipping incremental model {model_name}: zero features.")
            continue

        log(f"Fitting incremental model: {model_name} ({len(ind)} features)")

        model, method = fit_coxph_or_ridge_fallback(
            X_trainval[:, ind],
            y_trainval,
            model_name=model_name
        )

        fitted_models[model_name] = model
        fitting_methods[model_name] = method

        for split, (X, _) in split_data.items():
            risk_predictions[split][model_name] = predict_risk_score(model, X[:, ind])

    fitting_methods["M3_full_covariates_plus_brain_mri"] = "Selected elastic-net Cox model from main pipeline"

    for split, (X, _) in split_data.items():
        risk_predictions[split]["M3_full_covariates_plus_brain_mri"] = predict_risk_score(final_model, X)

    rows = []

    for model_name in ["M0_age_sex", "M1_covariate_baseline", "M2_brain_mri_only", "M3_full_covariates_plus_brain_mri"]:
        if model_name not in fitting_methods:
            continue

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
                "n_events": int(np.sum(y["event"]))
            })

    model_comparison_df = pd.DataFrame(rows)

    if "M1_covariate_baseline" in risk_predictions["test"]:
        delta_cindex_df = pd.DataFrame([
            paired_bootstrap_delta_cindex(
                y_test,
                risk_predictions["test"]["M3_full_covariates_plus_brain_mri"],
                risk_predictions["test"]["M1_covariate_baseline"],
                n_bootstrap,
                random_state
            )
        ])
    else:
        delta_cindex_df = pd.DataFrame([{
            "comparison": "M3_full_covariates_plus_brain_mri_vs_M1_covariate_baseline",
            "cindex_full": np.nan,
            "cindex_baseline": np.nan,
            "delta_cindex": np.nan,
            "interpretation": "M1 covariate baseline unavailable."
        }])

    return {
        "model_comparison_df": model_comparison_df,
        "delta_cindex_df": delta_cindex_df,
        "fitted_models": fitted_models,
        "fitting_methods": fitting_methods,
        "risk_predictions": risk_predictions,
        "feature_indices": idx
    }


# ============================================================
# Clock acceleration residualization
# ============================================================

def add_clock_age_and_acceleration(pred_df, covariate_cols, train_mask_col="split"):
    df = pred_df.copy()

    covariate_cols = [c for c in covariate_cols if c in df.columns]

    if not covariate_cols:
        warnings.warn("No residualization covariates available. Acceleration set to missing.")
        df[ACCEL_Z_COL] = np.nan
        df[ACCEL_YEARS_COL] = np.nan
        df[CLOCK_AGE_COL] = np.nan
        return df, None

    train = df.loc[df[train_mask_col] == "train"].copy()

    if train.shape[0] < 10:
        warnings.warn("Too few training samples for clock-acceleration residualization.")
        df[ACCEL_Z_COL] = np.nan
        df[ACCEL_YEARS_COL] = np.nan
        df[CLOCK_AGE_COL] = np.nan
        return df, None

    numeric_covs = []
    categorical_covs = []

    for c in covariate_cols:
        if c in ["Sex", "SITE"] or not pd.api.types.is_numeric_dtype(df[c]):
            categorical_covs.append(c)
        else:
            numeric_covs.append(c)

    transformers = []

    if numeric_covs:
        transformers.append((
            "num",
            Pipeline([("imputer", SimpleImputer(strategy="median"))]),
            numeric_covs
        ))

    if categorical_covs:
        transformers.append((
            "cat",
            Pipeline([
                ("imputer", SimpleImputer(strategy="most_frequent")),
                ("onehot", make_onehot_encoder())
            ]),
            categorical_covs
        ))

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

    lr = LinearRegression()
    lr.fit(X_train, train[RISK_COL].values)

    expected = lr.predict(X_all)
    resid_raw = df[RISK_COL].values - expected

    train_index = df[train_mask_col].values == "train"
    resid_mean = float(np.nanmean(resid_raw[train_index]))
    resid_sd = float(np.nanstd(resid_raw[train_index], ddof=1))

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
    if "num__Age" in feat_names:
        beta_age = float(lr.coef_[feat_names.index("num__Age")])

    if np.isfinite(beta_age) and abs(beta_age) > 1e-8:
        df[ACCEL_YEARS_COL] = resid / beta_age
        df[CLOCK_AGE_COL] = df["Age"] + df[ACCEL_YEARS_COL]
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
        "note": "Acceleration z-score is the residualized Cox risk score after adjustment for retained non-brain covariates."
    }

    return df, info


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


# ============================================================
# Main
# ============================================================

def main():
    args = parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    prefix = args.prefix
    l1_ratios = tuple(float(x) for x in args.l1_ratios.split(","))

    baseline_dx = normalize_dx(args.baseline_dx)
    event_dx_set = set([normalize_dx(x) for x in parse_list_arg(args.event_dx)])
    covariates = parse_list_arg(args.covariates)

    log("============================================================")
    log("Build ADNI brain MRI AD L'EPOCH")
    log("============================================================")
    log(f"Input file: {args.input_file}")
    log(f"Output dir: {outdir}")
    log(f"Baseline DX: {baseline_dx}")
    log(f"Event DX: {sorted(event_dx_set)}")
    log(f"Covariates: {covariates}")
    log(f"Hard-coded MUSE GM ROIs requested: {len(MUSE_GM_ROIS)}")
    log(f"Minimum nonzero brain features required: {args.min_nonzero_brain_features}")
    log("============================================================")

    # ------------------------------------------------------------
    # Load ADNI file
    # ------------------------------------------------------------

    df_raw = read_table(args.input_file)

    required = [
        args.id_col,
        args.visit_col,
        args.date_col,
        args.dx_col
    ]

    missing_required = [c for c in required if c not in df_raw.columns]
    if missing_required:
        raise ValueError(f"Missing required columns: {missing_required}")

    # ------------------------------------------------------------
    # Select hard-coded MUSE GM ROI features
    # ------------------------------------------------------------

    roi_available, roi_missing = select_hardcoded_muse_gm_rois(df_raw)

    pd.DataFrame({
        "roi_feature": MUSE_GM_ROIS,
        "present_in_file": [c in df_raw.columns for c in MUSE_GM_ROIS]
    }).to_csv(
        outdir / f"{prefix}_hardcoded_muse_gm_roi_presence.tsv",
        sep="\t",
        index=False
    )

    log(f"MUSE GM ROIs available in file: {len(roi_available)} / {len(MUSE_GM_ROIS)}")

    # ------------------------------------------------------------
    # Build survival dataset
    # ------------------------------------------------------------

    df = construct_adni_survival_dataset(
        df=df_raw,
        id_col=args.id_col,
        visit_col=args.visit_col,
        date_col=args.date_col,
        dx_col=args.dx_col,
        baseline_dx=baseline_dx,
        event_dx_set=event_dx_set,
        min_followup_days=args.min_followup_days
    )

    log("Constructed ADNI baseline-CN survival dataset:")
    log(f"  N baseline CN with follow-up = {df.shape[0]}")
    log(f"  Events MCI/AD = {int(df['event'].sum())}")
    log(f"  Censored = {int((~df['event']).sum())}")
    log(f"  Median follow-up years = {df['time_years'].median():.2f}")

    if int(df["event"].sum()) < 10:
        warnings.warn(
            "Very few conversion events. ADNI AD L'EPOCH model may be unstable."
        )

    # ------------------------------------------------------------
    # Basic covariate checks and types
    # ------------------------------------------------------------

    for c in covariates:
        if c not in df.columns:
            warnings.warn(f"Requested covariate {c} is missing from baseline data.")

    if "Age" in df.columns:
        df["Age"] = pd.to_numeric(df["Age"], errors="coerce")

    if "DLICV" in df.columns:
        df["DLICV"] = pd.to_numeric(df["DLICV"], errors="coerce")

    if "Sex" in df.columns:
        df["Sex"] = df["Sex"].astype(str).str.strip().replace({
            "0": "Female",
            "0.0": "Female",
            "1": "Male",
            "1.0": "Male",
            "F": "Female",
            "M": "Male",
            "female": "Female",
            "male": "Male"
        })

    if "SITE" in df.columns:
        df["SITE"] = df["SITE"].astype(str).str.strip()

    # ------------------------------------------------------------
    # Feature QC
    # ------------------------------------------------------------

    roi_cols, roi_qc = filter_features_by_missingness_and_variance(
        df,
        roi_available,
        max_missing=args.max_feature_missing
    )

    roi_qc.to_csv(
        outdir / f"{prefix}_muse_gm_roi_qc.tsv",
        sep="\t",
        index=False
    )

    pd.DataFrame({"roi_feature": roi_cols}).to_csv(
        outdir / f"{prefix}_selected_muse_gm_rois.tsv",
        sep="\t",
        index=False
    )

    log(f"MUSE GM ROIs retained after missingness/variance QC: {len(roi_cols)}")

    # ------------------------------------------------------------
    # Build model design matrix columns
    # ------------------------------------------------------------

    df, numeric_cols, categorical_cols, numeric_covariates, categorical_covariates = build_design_matrix_columns(
        df=df,
        roi_cols=roi_cols,
        covariates=covariates
    )

    required_model_cols = [args.id_col, "time_years", "event"] + [
        c for c in ["Age", "Sex"] if c in df.columns
    ]

    df = df.dropna(subset=required_model_cols).copy()

    survival_keep_cols = [
        args.id_col,
        args.visit_col,
        args.date_col,
        args.dx_col,
        "baseline_dx",
        "event_or_censor_dx",
        "event_or_censor_visit",
        "event_or_censor_date",
        "time_days",
        "time_years",
        "event"
    ]

    for c in covariates:
        if c in df.columns and c not in survival_keep_cols:
            survival_keep_cols.append(c)

    df[survival_keep_cols].to_csv(
        outdir / f"{prefix}_survival_dataset.tsv",
        sep="\t",
        index=False
    )

    log("Final model dataset:")
    log(f"  N = {df.shape[0]}")
    log(f"  Events = {int(df['event'].sum())}")
    log(f"  Censored = {int((~df['event']).sum())}")
    log(f"  Numeric covariates = {numeric_covariates}")
    log(f"  Categorical covariates = {categorical_covariates}")

    # ------------------------------------------------------------
    # Train/validation/test split
    # ------------------------------------------------------------

    log("Splitting into train/validation/test...")

    strat_all = make_stratify_vector(
        df,
        age_col="Age" if "Age" in df.columns else None,
        event_col="event",
        age_bins=args.stratify_age_bins
    )

    df_trainval, df_test = train_test_split(
        df,
        test_size=args.test_size,
        random_state=args.random_state,
        stratify=strat_all
    )

    strat_trainval = make_stratify_vector(
        df_trainval,
        age_col="Age" if "Age" in df_trainval.columns else None,
        event_col="event",
        age_bins=args.stratify_age_bins
    )

    df_train, df_val = train_test_split(
        df_trainval,
        test_size=args.validation_size,
        random_state=args.random_state,
        stratify=strat_trainval
    )

    log(f"  Train N={df_train.shape[0]}, events={int(df_train['event'].sum())}")
    log(f"  Val   N={df_val.shape[0]}, events={int(df_val['event'].sum())}")
    log(f"  Test  N={df_test.shape[0]}, events={int(df_test['event'].sum())}")

    all_cols = numeric_cols + categorical_cols

    X_train_raw = df_train[all_cols].copy()
    X_val_raw = df_val[all_cols].copy()
    X_test_raw = df_test[all_cols].copy()
    X_trainval_raw = df_trainval[all_cols].copy()

    preprocessor = make_preprocessor(numeric_cols, categorical_cols)

    X_train = preprocessor.fit_transform(X_train_raw)
    X_val = preprocessor.transform(X_val_raw)
    X_test = preprocessor.transform(X_test_raw)
    X_trainval = preprocessor.transform(X_trainval_raw)

    feature_names = get_feature_names(preprocessor)
    is_brain = brain_feature_mask(feature_names)

    log(f"Transformed total features: {len(feature_names)}")
    log(f"Transformed MUSE GM ROI features: {int(is_brain.sum())}")

    y_train = Surv.from_arrays(
        event=df_train["event"].astype(bool).values,
        time=df_train["time_years"].astype(float).values
    )

    y_val = Surv.from_arrays(
        event=df_val["event"].astype(bool).values,
        time=df_val["time_years"].astype(float).values
    )

    y_test = Surv.from_arrays(
        event=df_test["event"].astype(bool).values,
        time=df_test["time_years"].astype(float).values
    )

    y_trainval = Surv.from_arrays(
        event=df_trainval["event"].astype(bool).values,
        time=df_trainval["time_years"].astype(float).values
    )

    # ------------------------------------------------------------
    # Tune Coxnet with required nonzero brain features
    # ------------------------------------------------------------

    log("Tuning elastic-net Cox model with nonzero MUSE GM ROI requirement...")

    best, penalty_factor = fit_and_select_coxnet(
        X_train=X_train,
        y_train=y_train,
        X_val=X_val,
        y_val=y_val,
        feature_names=feature_names,
        l1_ratios=l1_ratios,
        n_alphas=args.n_alphas,
        alpha_min_ratio=args.alpha_min_ratio,
        min_nonzero_brain_features=args.min_nonzero_brain_features
    )

    log("Best validation model:")
    log(json.dumps({k: v for k, v in best.items() if k != "coef"}, indent=2))

    # ------------------------------------------------------------
    # Refit final model on train+validation
    # ------------------------------------------------------------

    log("Refitting final model on train+validation...")

    final_model = fit_final_model(
        X_trainval=X_trainval,
        y_trainval=y_trainval,
        best=best,
        penalty_factor=penalty_factor
    )

    coef_final = np.asarray(final_model.coef_).reshape(-1)
    nonzero_final = np.abs(coef_final) > 1e-8
    final_nonzero_brain = int(np.sum(nonzero_final & is_brain))

    if final_nonzero_brain == 0:
        raise RuntimeError(
            "Final model has zero nonzero MUSE GM ROI coefficients. "
            "Stopping because the L'EPOCH would not be brain-imaging related."
        )

    if final_nonzero_brain < args.min_nonzero_brain_features:
        warnings.warn(
            f"Final model retained {final_nonzero_brain} nonzero MUSE GM ROI features, "
            f"less than requested minimum {args.min_nonzero_brain_features}. "
            "Proceeding because at least one brain-imaging feature is retained."
        )

    log(f"Final nonzero MUSE GM ROI coefficients: {final_nonzero_brain}")

    # ------------------------------------------------------------
    # Predictions and C-index
    # ------------------------------------------------------------

    log("Generating risk scores...")

    risk_train = predict_risk_score(final_model, X_train)
    risk_val = predict_risk_score(final_model, X_val)
    risk_test = predict_risk_score(final_model, X_test)
    risk_trainval = predict_risk_score(final_model, X_trainval)

    cindex_train = compute_cindex(y_train, risk_train)
    cindex_val = compute_cindex(y_val, risk_val)
    cindex_test = compute_cindex(y_test, risk_test)
    cindex_trainval = compute_cindex(y_trainval, risk_trainval)

    log(f"Train C-index:      {cindex_train:.4f}")
    log(f"Validation C-index: {cindex_val:.4f}")
    log(f"Train+Val C-index:  {cindex_trainval:.4f}")
    log(f"Test C-index:       {cindex_test:.4f}")

    # ------------------------------------------------------------
    # Incremental value analysis
    # ------------------------------------------------------------

    log("Running incremental-value analysis...")

    inc = run_incremental_value_analysis(
        X_train=X_train,
        X_val=X_val,
        X_test=X_test,
        X_trainval=X_trainval,
        y_train=y_train,
        y_val=y_val,
        y_test=y_test,
        y_trainval=y_trainval,
        feature_names=feature_names,
        final_model=final_model,
        n_bootstrap=args.n_bootstrap_incremental,
        random_state=args.random_state
    )

    model_comparison_df = inc["model_comparison_df"]
    delta_cindex_df = inc["delta_cindex_df"]

    model_comparison_df.to_csv(
        outdir / f"{prefix}_model_comparison.tsv",
        sep="\t",
        index=False
    )

    delta_cindex_df.to_csv(
        outdir / f"{prefix}_incremental_value_delta_cindex.tsv",
        sep="\t",
        index=False
    )

    log(delta_cindex_df.to_string(index=False))

    # ------------------------------------------------------------
    # Participant-level predictions
    # ------------------------------------------------------------

    def make_pred_frame(part, split, risk):
        base = [
            args.id_col,
            args.visit_col,
            args.date_col,
            args.dx_col,
            "event_or_censor_dx",
            "event_or_censor_visit",
            "event_or_censor_date",
            "time_years",
            "event"
        ]

        for c in covariates:
            if c in part.columns and c not in base:
                base.append(c)

        out = part[base].copy()
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
        "M3_full_covariates_plus_brain_mri": "risk_score_M3_full_covariates_plus_brain_mri"
    }

    for model_name, col in comparison_cols.items():
        if model_name in inc["risk_predictions"]["train"]:
            pred_train[col] = inc["risk_predictions"]["train"][model_name]
            pred_val[col] = inc["risk_predictions"]["validation"][model_name]
            pred_test[col] = inc["risk_predictions"]["test"][model_name]

    risk_times = [1.0, 2.0, 3.0, 5.0]

    pred_train = pd.concat(
        [pred_train.reset_index(drop=True), predict_absolute_risk(final_model, X_train, risk_times)],
        axis=1
    )

    pred_val = pd.concat(
        [pred_val.reset_index(drop=True), predict_absolute_risk(final_model, X_val, risk_times)],
        axis=1
    )

    pred_test = pd.concat(
        [pred_test.reset_index(drop=True), predict_absolute_risk(final_model, X_test, risk_times)],
        axis=1
    )

    pred_all = pd.concat([pred_train, pred_val, pred_test], ignore_index=True)

    # ------------------------------------------------------------
    # Residualized L'EPOCH acceleration
    # ------------------------------------------------------------

    log("Adding residualized AD L'EPOCH acceleration...")

    residualization_covariates = [c for c in covariates if c in pred_all.columns]

    pred_all, clock_transform_info = add_clock_age_and_acceleration(
        pred_all,
        covariate_cols=residualization_covariates,
        train_mask_col="split"
    )

    pred_all.to_csv(
        outdir / f"{prefix}_predictions.tsv",
        sep="\t",
        index=False
    )

    pred_test.to_csv(
        outdir / f"{prefix}_test_predictions.tsv",
        sep="\t",
        index=False
    )

    # ------------------------------------------------------------
    # Coefficients
    # ------------------------------------------------------------

    log("Saving coefficients...")

    coef_df = pd.DataFrame({
        "feature": feature_names,
        "coefficient": coef_final,
        "abs_coefficient": np.abs(coef_final),
        "penalty_factor": penalty_factor,
        "is_nonzero": nonzero_final,
        "is_muse_gm_roi": is_brain
    }).sort_values("abs_coefficient", ascending=False)

    coef_df.to_csv(
        outdir / f"{prefix}_coefficients.tsv",
        sep="\t",
        index=False
    )

    nonzero_coef_df = coef_df.loc[coef_df["is_nonzero"]].copy()

    nonzero_coef_df.to_csv(
        outdir / f"{prefix}_nonzero_coefficients.tsv",
        sep="\t",
        index=False
    )

    # ------------------------------------------------------------
    # Model bundle and performance
    # ------------------------------------------------------------

    model_comparison_records = json.loads(model_comparison_df.to_json(orient="records"))
    delta_cindex_records = json.loads(delta_cindex_df.to_json(orient="records"))

    model_bundle = {
        "preprocessor": preprocessor,
        "model": final_model,
        "feature_names": feature_names,
        "hardcoded_muse_gm_rois": MUSE_GM_ROIS,
        "available_muse_gm_rois": roi_available,
        "selected_muse_gm_rois": roi_cols,
        "missing_muse_gm_rois": roi_missing,
        "numeric_cols": numeric_cols,
        "categorical_cols": categorical_cols,
        "numeric_covariates": numeric_covariates,
        "categorical_covariates": categorical_covariates,
        "covariates": covariates,
        "best": {k: v for k, v in best.items() if k != "coef"},
        "penalty_factor": penalty_factor,
        "clock_transform_info": clock_transform_info,
        "incremental_value_model_comparison": model_comparison_records,
        "incremental_value_delta_cindex": delta_cindex_records,
        "incremental_value_fitting_methods": inc["fitting_methods"],
        "input_file": args.input_file,
        "time_zero": "ADNI baseline visit, preferably Visit_Code == bl",
        "event": "First follow-up diagnosis of MCI or AD",
        "censoring": "Last available follow-up visit without MCI/AD conversion"
    }

    joblib.dump(
        model_bundle,
        outdir / f"{prefix}_model.joblib"
    )

    if delta_cindex_df.shape[0] > 0:
        delta_row = delta_cindex_df.iloc[0].to_dict()
    else:
        delta_row = {}

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
        "best_l1_ratio": float(best["l1_ratio"]),
        "best_alpha": float(best["alpha"]),
        "best_validation_cindex_during_tuning": float(best["cindex"]),
        "used_penalty_factor": bool(best["used_penalty_factor"]),
        "min_nonzero_brain_features_requested": int(args.min_nonzero_brain_features),
        "best_model_nonzero_brain_features_validation": int(best["n_nonzero_brain_features"]),
        "final_model_nonzero_brain_features": int(final_nonzero_brain),
        "n_hardcoded_muse_gm_rois_requested": int(len(MUSE_GM_ROIS)),
        "n_muse_gm_rois_available": int(len(roi_available)),
        "n_muse_gm_rois_after_qc": int(len(roi_cols)),
        "n_nonzero_coefficients_total": int(nonzero_coef_df.shape[0]),
        "n_nonzero_muse_gm_roi_coefficients": int(final_nonzero_brain),
        "covariates": covariates,
        "time_zero": "ADNI baseline visit, preferably Visit_Code == bl",
        "eligible_baseline_dx": baseline_dx,
        "event_dx": sorted(list(event_dx_set)),
        "event_definition": "First follow-up diagnosis of MCI or AD",
        "censoring": "Last available follow-up visit without MCI/AD conversion",
        "primary_score": RISK_COL,
        "recommended_downstream_score": ACCEL_Z_COL,
        "incremental_value_model_comparison": model_comparison_records,
        "incremental_value_delta_cindex": delta_cindex_records,
        "delta_cindex_test_M3_vs_M1": delta_row.get("delta_cindex", np.nan),
        "note": (
            "Primary AD L'EPOCH score is the Cox log-risk score. "
            "Acceleration z-score is residualized on retained non-brain covariates. "
            "Model selection enforces nonzero MUSE GM ROI coefficients to ensure that "
            "the clock remains brain-imaging related."
        )
    }

    with open(outdir / f"{prefix}_performance.json", "w") as f:
        json.dump(performance, f, indent=2)

    pd.DataFrame([performance]).to_csv(
        outdir / f"{prefix}_performance.tsv",
        sep="\t",
        index=False
    )

    # ------------------------------------------------------------
    # Finish
    # ------------------------------------------------------------

    log("============================================================")
    log("ADNI brain MRI AD L'EPOCH complete.")
    log("Main output files:")
    for name in [
        "survival_dataset.tsv",
        "predictions.tsv",
        "test_predictions.tsv",
        "coefficients.tsv",
        "nonzero_coefficients.tsv",
        "model_comparison.tsv",
        "incremental_value_delta_cindex.tsv",
        "performance.json",
        "model.joblib"
    ]:
        log(f"  {outdir / f'{prefix}_{name}'}")
    log("============================================================")


if __name__ == "__main__":
    main()