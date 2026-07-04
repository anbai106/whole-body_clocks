#!/usr/bin/env python3
# ============================================================
# Fix unstable immune metabolomics mortality-clock transform
#
# This does NOT retrain the Cox survival model.
# It only replaces the unstable post-hoc residualization transform:
#   immune_metabolomics_mortality_clock_acceleration_z
#   immune_metabolomics_mortality_clock_acceleration_years
#   immune_metabolomics_mortality_clock_age_years
#
# It:
#   1. Loads baseline immune predictions.
#   2. Refits a stable residualizer on the baseline TRAIN split.
#   3. Uses drop-first one-hot coding + small ridge penalty.
#   4. Updates the immune model joblib clock_transform_info.
#   5. Patches baseline and instance-1 prediction TSVs.

## This script will create backups before overwriting:
# immune_metabolomics_mortality_clock_model.joblib.before_stable_immune_clock_transform.joblib
# immune_metabolomics_mortality_clock_predictions.tsv.before_stable_immune_clock_transform.tsv
# immune_metabolomics_mortality_clock_apply_instance_1_0_predictions.tsv.before_stable_immune_clock_transform.tsv
# ============================================================

import json
import shutil
import warnings
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder


# -----------------------------
# 1. Paths
# -----------------------------

root = Path("/cbica/home/wenju/Reproducibile_paper/WholeBodyClock")

long_root = Path(
    "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/"
    "mortality_clock/longitudinal/metabolomics"
)

organ_label = "Immune"
organ = "immune"

model_joblib = (
    root
    / f"{organ_label}_metabolomics_mortality_clock"
    / f"{organ}_metabolomics_mortality_clock_model.joblib"
)

baseline_pred_file = (
    root
    / f"{organ_label}_metabolomics_mortality_clock"
    / f"{organ}_metabolomics_mortality_clock_predictions.tsv"
)

instance1_pred_file = (
    long_root
    / organ_label
    / f"{organ}_metabolomics_mortality_clock_apply_instance_1_0_predictions.tsv"
)

out_debug_json = (
    root
    / f"{organ_label}_metabolomics_mortality_clock"
    / f"{organ}_metabolomics_mortality_clock_stable_transform_debug.json"
)


# -----------------------------
# 2. Column names
# -----------------------------

risk_col = f"{organ}_metabolomics_mortality_risk_score"
z_col = f"{organ}_metabolomics_mortality_clock_acceleration_z"
yrs_col = f"{organ}_metabolomics_mortality_clock_acceleration_years"
age_col = f"{organ}_metabolomics_mortality_clock_age_years"

covariate_cols = [
    "age_at_baseline",
    "bmi_at_baseline",
    "diastolic_bp_at_baseline",
    "systolic_bp_at_baseline",
    "sex",
    "smoking_status_at_baseline",
    "uk_biobank_assessment_centre_f54_0_0",
]

numeric_covs = [
    "age_at_baseline",
    "bmi_at_baseline",
    "diastolic_bp_at_baseline",
    "systolic_bp_at_baseline",
]

categorical_covs = [
    "sex",
    "smoking_status_at_baseline",
    "uk_biobank_assessment_centre_f54_0_0",
]


# -----------------------------
# 3. Helpers
# -----------------------------

def make_onehot_drop_first():
    """
    drop='first' avoids the dummy-variable trap when an intercept is included.
    handle_unknown='ignore' allows instance-1 categories not seen in baseline train.
    """
    try:
        return OneHotEncoder(
            handle_unknown="ignore",
            drop="first",
            sparse_output=False,
        )
    except TypeError:
        return OneHotEncoder(
            handle_unknown="ignore",
            drop="first",
            sparse=False,
        )


def backup_once(path: Path, tag: str) -> Path:
    backup = path.with_name(path.name + f".before_{tag}")
    if not backup.exists():
        shutil.copy2(path, backup)
        print(f"Backup written: {backup}")
    else:
        print(f"Backup already exists: {backup}")
    return backup


def ensure_columns(df, cols):
    df = df.copy()
    for c in cols:
        if c not in df.columns:
            warnings.warn(f"Missing covariate column {c}; creating as NA.")
            df[c] = np.nan
    return df


def clean_covariates(df):
    df = ensure_columns(df, covariate_cols).copy()

    for c in numeric_covs:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    for c in categorical_covs:
        df[c] = df[c].astype("object")
        df.loc[df[c].isna(), c] = np.nan

    return df


def get_feature_names(prep):
    feat_names = []

    if "num" in prep.named_transformers_:
        feat_names.extend([f"num__{c}" for c in numeric_covs])

    if "cat" in prep.named_transformers_:
        ohe = prep.named_transformers_["cat"].named_steps["onehot"]
        try:
            cat_names = ohe.get_feature_names_out(categorical_covs)
        except AttributeError:
            cat_names = ohe.get_feature_names(categorical_covs)

        feat_names.extend([f"cat__{x}" for x in cat_names])

    return list(feat_names)


def build_preprocessor():
    num_pipe = Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
        ]
    )

    cat_pipe = Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="most_frequent")),
            ("onehot", make_onehot_drop_first()),
        ]
    )

    prep = ColumnTransformer(
        transformers=[
            ("num", num_pipe, numeric_covs),
            ("cat", cat_pipe, categorical_covs),
        ],
        remainder="drop",
    )

    return prep


def fit_stable_transform_from_baseline(baseline_df):
    if risk_col not in baseline_df.columns:
        raise ValueError(f"{risk_col} not found in baseline prediction file.")

    if "split" not in baseline_df.columns:
        warnings.warn("split column not found in baseline file. Using all baseline rows to fit residualizer.")
        train_mask = np.repeat(True, baseline_df.shape[0])
    else:
        train_mask = baseline_df["split"].astype(str).str.lower().eq("train").values

    if train_mask.sum() < 10:
        raise ValueError(f"Too few training rows for stable transform: {train_mask.sum()}")

    baseline_df = clean_covariates(baseline_df)

    train_df = baseline_df.loc[train_mask].copy()

    X_train_raw = train_df[covariate_cols].copy()
    X_all_raw = baseline_df[covariate_cols].copy()

    prep = build_preprocessor()
    X_train = prep.fit_transform(X_train_raw)
    X_all = prep.transform(X_all_raw)

    # Small ridge penalty prevents numerical instability.
    # This changes only the post-hoc acceleration calibration, not the Cox model.
    ridge = Ridge(alpha=1e-4, fit_intercept=True)
    ridge.fit(X_train, pd.to_numeric(train_df[risk_col], errors="coerce").values)

    expected_all = ridge.predict(X_all)
    risk_all = pd.to_numeric(baseline_df[risk_col], errors="coerce").values

    resid_raw_all = risk_all - expected_all
    resid_raw_train = resid_raw_all[train_mask]

    mean_train = float(np.nanmean(resid_raw_train))
    sd_train = float(np.nanstd(resid_raw_train))

    feat_names = get_feature_names(prep)
    coef = np.asarray(ridge.coef_).reshape(-1)

    if len(feat_names) != len(coef):
        raise RuntimeError(f"Feature-name length {len(feat_names)} != coefficient length {len(coef)}")

    coef_dict = {k: float(v) for k, v in zip(feat_names, coef)}

    beta_age = coef_dict.get("num__age_at_baseline", np.nan)

    max_abs_coef = float(np.nanmax(np.abs(coef))) if len(coef) else np.nan

    if not np.isfinite(sd_train) or sd_train <= 0:
        raise ValueError(f"Invalid residual SD: {sd_train}")

    if not np.isfinite(beta_age) or abs(beta_age) <= 1e-8:
        raise ValueError(f"Invalid age coefficient: {beta_age}")

    if max_abs_coef > 1000:
        raise ValueError(
            f"Stable residualizer still has very large coefficients: max_abs_coef={max_abs_coef}"
        )

    info = {
        "residualization_covariates": covariate_cols,
        "numeric_residualization_covariates": numeric_covs,
        "categorical_residualization_covariates": categorical_covs,
        "risk_score_covariate_model_intercept": float(ridge.intercept_),
        "risk_score_covariate_model_coef": coef_dict,
        "adjusted_age_coefficient_risk_score_per_year": float(beta_age),
        "risk_score_residual_mean_train": mean_train,
        "risk_score_residual_sd_train": sd_train,
        "note": (
            "Stable immune metabolomics clock transform refit from baseline predictions. "
            "Uses drop-first one-hot encoding plus Ridge(alpha=1e-4) to avoid the "
            "dummy-variable trap and numerical instability. Cox model was not retrained."
        ),
    }

    debug = {
        "n_baseline_rows": int(baseline_df.shape[0]),
        "n_train_rows": int(train_mask.sum()),
        "ridge_alpha": 1e-4,
        "intercept": float(ridge.intercept_),
        "max_abs_coef": max_abs_coef,
        "beta_age": float(beta_age),
        "residual_mean_train": mean_train,
        "residual_sd_train": sd_train,
        "top_abs_coefficients": sorted(
            [
                {"feature": k, "coefficient": float(v), "abs_coefficient": abs(float(v))}
                for k, v in coef_dict.items()
            ],
            key=lambda d: d["abs_coefficient"],
            reverse=True,
        )[:30],
    }

    return prep, ridge, info, debug


def patch_prediction_file(path, prep, ridge, info, age_for_transform_col, age_for_clock_col):
    if not path.exists():
        warnings.warn(f"Prediction file does not exist; skipping: {path}")
        return None

    backup_once(path, "stable_immune_clock_transform.tsv")

    df = pd.read_csv(path, sep="\t")

    if risk_col not in df.columns:
        raise ValueError(f"{risk_col} not found in {path}")

    df = ensure_columns(df, covariate_cols)

    # For longitudinal instance 1, the model column is still named age_at_baseline,
    # but the biologically correct sample-age covariate is age_at_imaging.
    if age_for_transform_col not in df.columns:
        raise ValueError(f"{age_for_transform_col} not found in {path}")

    df["age_at_baseline"] = pd.to_numeric(df[age_for_transform_col], errors="coerce")

    df_clean = clean_covariates(df)

    X = prep.transform(df_clean[covariate_cols])
    expected = ridge.predict(X)

    risk = pd.to_numeric(df[risk_col], errors="coerce").values
    resid_raw = risk - expected
    resid = resid_raw - float(info["risk_score_residual_mean_train"])

    sd_train = float(info["risk_score_residual_sd_train"])
    beta_age = float(info["adjusted_age_coefficient_risk_score_per_year"])

    df[z_col] = resid / sd_train
    df[yrs_col] = resid / beta_age

    if age_for_clock_col in df.columns:
        clock_base_age = pd.to_numeric(df[age_for_clock_col], errors="coerce")
    else:
        clock_base_age = pd.to_numeric(df[age_for_transform_col], errors="coerce")

    df[age_col] = clock_base_age + df[yrs_col]

    # Hard sanity check
    z_abs_max = float(np.nanmax(np.abs(pd.to_numeric(df[z_col], errors="coerce"))))
    years_abs_max = float(np.nanmax(np.abs(pd.to_numeric(df[yrs_col], errors="coerce"))))

    if z_abs_max > 50 or years_abs_max > 200:
        raise ValueError(
            f"Patched values still look invalid for {path}: "
            f"max_abs_z={z_abs_max}, max_abs_years={years_abs_max}"
        )

    df.to_csv(path, sep="\t", index=False)

    summary = {
        "file": str(path),
        "n_rows": int(df.shape[0]),
        "z_min": float(np.nanmin(df[z_col])),
        "z_median": float(np.nanmedian(df[z_col])),
        "z_mean": float(np.nanmean(df[z_col])),
        "z_max": float(np.nanmax(df[z_col])),
        "years_min": float(np.nanmin(df[yrs_col])),
        "years_median": float(np.nanmedian(df[yrs_col])),
        "years_mean": float(np.nanmean(df[yrs_col])),
        "years_max": float(np.nanmax(df[yrs_col])),
    }

    print("\nPatched:", path)
    print(json.dumps(summary, indent=2))

    return summary


# -----------------------------
# 4. Main
# -----------------------------

def main():
    print("============================================================")
    print("Fixing immune metabolomics mortality-clock transform")
    print("============================================================")
    print("Model joblib:", model_joblib)
    print("Baseline predictions:", baseline_pred_file)
    print("Instance-1 predictions:", instance1_pred_file)

    if not model_joblib.exists():
        raise FileNotFoundError(model_joblib)

    if not baseline_pred_file.exists():
        raise FileNotFoundError(baseline_pred_file)

    baseline_df = pd.read_csv(baseline_pred_file, sep="\t")

    prep, ridge, new_info, debug = fit_stable_transform_from_baseline(baseline_df)

    print("\nStable residualizer debug:")
    print(json.dumps(debug, indent=2))

    # Backup and patch model bundle
    backup_once(model_joblib, "stable_immune_clock_transform.joblib")

    bundle = joblib.load(model_joblib)
    old_info = bundle.get("clock_transform_info", None)

    bundle["clock_transform_info_original_unstable_backup"] = old_info
    bundle["clock_transform_info"] = new_info

    joblib.dump(bundle, model_joblib)

    print("\nUpdated model joblib clock_transform_info:")
    print(model_joblib)

    # Patch current baseline and instance-1 prediction files.
    # Baseline: sample age is age_at_baseline.
    baseline_summary = patch_prediction_file(
        path=baseline_pred_file,
        prep=prep,
        ridge=ridge,
        info=new_info,
        age_for_transform_col="age_at_baseline",
        age_for_clock_col="age_at_baseline",
    )

    # Instance 1: sample age is age_at_imaging.
    instance1_summary = patch_prediction_file(
        path=instance1_pred_file,
        prep=prep,
        ridge=ridge,
        info=new_info,
        age_for_transform_col="age_at_imaging",
        age_for_clock_col="age_at_imaging",
    )

    debug_out = {
        "model_joblib": str(model_joblib),
        "baseline_pred_file": str(baseline_pred_file),
        "instance1_pred_file": str(instance1_pred_file),
        "stable_transform_debug": debug,
        "baseline_summary": baseline_summary,
        "instance1_summary": instance1_summary,
    }

    with open(out_debug_json, "w") as f:
        json.dump(debug_out, f, indent=2)

    print("\nWrote debug JSON:")
    print(out_debug_json)

    print("\nFinished. The Cox model was not retrained; only the immune clock transform was stabilized.")


if __name__ == "__main__":
    main()