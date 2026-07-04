#!/usr/bin/env python3
# ============================================================
# 4_fix_pulmonary_proteomics_clock_transform.py
#
# Stabilize the Pulmonary proteomics mortality-clock transform.
#
# This DOES NOT retrain the Cox survival model.
# It only refits the post-hoc covariate residualization transform
# used to derive:
#   pulmonary_proteomics_mortality_clock_acceleration_z
#   pulmonary_proteomics_mortality_clock_acceleration_years
#   pulmonary_proteomics_mortality_clock_age_years
#
# It patches:
#   1) Baseline instance 0_0 predictions
#   2) Follow-up proteomics instance 2_0 predictions
#   3) Follow-up proteomics instance 3_0 predictions
#   4) The model joblib clock_transform_info
#
# Default reference split:
#   REFERENCE_SPLIT = "train"
#
# To use a different split, edit REFERENCE_SPLIT below or pass:
#   --reference-split test
# ============================================================

import argparse
import json
import shutil
import traceback
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
# 1. Defaults
# -----------------------------

ROOT_DEFAULT = Path("/cbica/home/wenju/Reproducibile_paper/WholeBodyClock")
LONG_ROOT_DEFAULT = Path(
    "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/"
    "mortality_clock/longitudinal/proteomics"
)

ORGAN_LABEL = "Pulmonary"
ORGAN = "pulmonary"

REFERENCE_SPLIT_DEFAULT = "train"
RIDGE_ALPHA_DEFAULT = 1e-4
BACKUP_TAG = "stable_pulmonary_proteomics_clock_transform"

MAX_ABS_Z_ALLOWED = 50
MAX_ABS_YEARS_ALLOWED = 200
MAX_ABS_COEF_ALLOWED = 1000

COVARIATE_COLS = [
    "age_at_baseline",
    "bmi_at_baseline",
    "diastolic_bp_at_baseline",
    "systolic_bp_at_baseline",
    "sex",
    "smoking_status_at_baseline",
    "uk_biobank_assessment_centre_f54_0_0",
]

NUMERIC_COVS = [
    "age_at_baseline",
    "bmi_at_baseline",
    "diastolic_bp_at_baseline",
    "systolic_bp_at_baseline",
]

CATEGORICAL_COVS = [
    "sex",
    "smoking_status_at_baseline",
    "uk_biobank_assessment_centre_f54_0_0",
]


# -----------------------------
# 2. Arguments
# -----------------------------

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--root", default=str(ROOT_DEFAULT))
    p.add_argument("--long-root", default=str(LONG_ROOT_DEFAULT))
    p.add_argument(
        "--reference-split",
        default=REFERENCE_SPLIT_DEFAULT,
        choices=["train", "test", "validation", "trainval", "all"],
    )
    p.add_argument("--ridge-alpha", type=float, default=RIDGE_ALPHA_DEFAULT)
    p.add_argument("--skip-missing-followup", action="store_true", default=True)
    return p.parse_args()


def get_paths(root: Path, long_root: Path):
    model_joblib = (
        root
        / "Pulmonary_proteomics_mortality_clock"
        / "pulmonary_proteomics_mortality_clock_model.joblib"
    )

    baseline_pred_file = (
        root
        / "Pulmonary_proteomics_mortality_clock"
        / "pulmonary_proteomics_mortality_clock_predictions.tsv"
    )

    followup_pred_files = {
        "2_0": (
            long_root
            / "Pulmonary"
            / "pulmonary_proteomics_mortality_clock_apply_instance_2_0_predictions.tsv"
        ),
        "3_0": (
            long_root
            / "Pulmonary"
            / "pulmonary_proteomics_mortality_clock_apply_instance_3_0_predictions.tsv"
        ),
    }

    outdir = long_root / "Pulmonary"
    outdir.mkdir(parents=True, exist_ok=True)

    debug_json = outdir / "pulmonary_proteomics_mortality_clock_stable_transform_debug.json"
    summary_json = outdir / "pulmonary_proteomics_stable_clock_transform_summary.json"

    return model_joblib, baseline_pred_file, followup_pred_files, debug_json, summary_json


def get_clock_cols():
    risk_col = "pulmonary_proteomics_mortality_risk_score"
    z_col = "pulmonary_proteomics_mortality_clock_acceleration_z"
    yrs_col = "pulmonary_proteomics_mortality_clock_acceleration_years"
    age_col = "pulmonary_proteomics_mortality_clock_age_years"
    return risk_col, z_col, yrs_col, age_col


# -----------------------------
# 3. Helper functions
# -----------------------------

def make_onehot_drop_first():
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
    df = ensure_columns(df, COVARIATE_COLS).copy()

    for c in NUMERIC_COVS:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    for c in CATEGORICAL_COVS:
        df[c] = df[c].astype("object")
        df.loc[df[c].isna(), c] = np.nan

    return df


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

    return ColumnTransformer(
        transformers=[
            ("num", num_pipe, NUMERIC_COVS),
            ("cat", cat_pipe, CATEGORICAL_COVS),
        ],
        remainder="drop",
    )


def get_feature_names(prep):
    feat_names = []

    if "num" in prep.named_transformers_:
        feat_names.extend([f"num__{c}" for c in NUMERIC_COVS])

    if "cat" in prep.named_transformers_:
        ohe = prep.named_transformers_["cat"].named_steps["onehot"]
        try:
            cat_names = ohe.get_feature_names_out(CATEGORICAL_COVS)
        except AttributeError:
            cat_names = ohe.get_feature_names(CATEGORICAL_COVS)
        feat_names.extend([f"cat__{x}" for x in cat_names])

    return list(feat_names)


def get_reference_mask(baseline_df, reference_split):
    reference_split = str(reference_split).lower()

    if reference_split == "all":
        return np.repeat(True, baseline_df.shape[0])

    if "split" not in baseline_df.columns:
        warnings.warn(
            "split column not found in baseline file. Using all baseline rows to fit residualizer."
        )
        return np.repeat(True, baseline_df.shape[0])

    split_values = baseline_df["split"].astype(str).str.lower()

    if reference_split == "trainval":
        mask = split_values.isin(["train", "validation"]).values
    else:
        mask = split_values.eq(reference_split).values

    if mask.sum() < 10:
        raise ValueError(
            f"Too few rows for reference_split={reference_split}: {mask.sum()}"
        )

    return mask


def first_existing_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None


# -----------------------------
# 4. Fit stable transform
# -----------------------------

def fit_stable_transform_from_baseline(
    baseline_df,
    risk_col,
    reference_split,
    ridge_alpha,
):
    if risk_col not in baseline_df.columns:
        raise ValueError(f"{risk_col} not found in baseline prediction file.")

    baseline_df = clean_covariates(baseline_df)
    ref_mask = get_reference_mask(baseline_df, reference_split)
    ref_df = baseline_df.loc[ref_mask].copy()

    prep = build_preprocessor()
    X_ref = prep.fit_transform(ref_df[COVARIATE_COLS].copy())
    X_all = prep.transform(baseline_df[COVARIATE_COLS].copy())

    y_ref = pd.to_numeric(ref_df[risk_col], errors="coerce").values
    if np.isnan(y_ref).all():
        raise ValueError(f"All reference risk scores are missing for {risk_col}.")

    ridge = Ridge(alpha=ridge_alpha, fit_intercept=True)
    ridge.fit(X_ref, y_ref)

    expected_all = ridge.predict(X_all)
    risk_all = pd.to_numeric(baseline_df[risk_col], errors="coerce").values
    resid_raw_all = risk_all - expected_all
    resid_raw_ref = resid_raw_all[ref_mask]

    mean_ref = float(np.nanmean(resid_raw_ref))
    sd_ref = float(np.nanstd(resid_raw_ref))

    feat_names = get_feature_names(prep)
    coef = np.asarray(ridge.coef_).reshape(-1)

    if len(feat_names) != len(coef):
        raise RuntimeError(
            f"Feature-name length {len(feat_names)} != coefficient length {len(coef)}"
        )

    coef_dict = {k: float(v) for k, v in zip(feat_names, coef)}
    beta_age = coef_dict.get("num__age_at_baseline", np.nan)
    max_abs_coef = float(np.nanmax(np.abs(coef))) if len(coef) else np.nan

    if not np.isfinite(sd_ref) or sd_ref <= 0:
        raise ValueError(f"Invalid residual SD: {sd_ref}")

    if not np.isfinite(beta_age) or abs(beta_age) <= 1e-8:
        raise ValueError(f"Invalid age coefficient: {beta_age}")

    if max_abs_coef > MAX_ABS_COEF_ALLOWED:
        raise ValueError(
            f"Stable residualizer still has very large coefficients: max_abs_coef={max_abs_coef}"
        )

    info = {
        "residualization_covariates": COVARIATE_COLS,
        "numeric_residualization_covariates": NUMERIC_COVS,
        "categorical_residualization_covariates": CATEGORICAL_COVS,
        "risk_score_covariate_model_intercept": float(ridge.intercept_),
        "risk_score_covariate_model_coef": coef_dict,
        "adjusted_age_coefficient_risk_score_per_year": float(beta_age),
        "risk_score_residual_mean_train": mean_ref,
        "risk_score_residual_sd_train": sd_ref,
        "reference_split_for_stable_transform": reference_split,
        "ridge_alpha": ridge_alpha,
        "note": (
            "Stable Pulmonary proteomics clock transform refit from baseline predictions. "
            "Uses drop-first one-hot encoding plus Ridge regularization to avoid the dummy-variable trap. "
            "The Cox mortality model itself was not retrained."
        ),
    }

    debug = {
        "organ_label": ORGAN_LABEL,
        "organ": ORGAN,
        "n_baseline_rows": int(baseline_df.shape[0]),
        "reference_split": reference_split,
        "n_reference_rows": int(ref_mask.sum()),
        "ridge_alpha": ridge_alpha,
        "intercept": float(ridge.intercept_),
        "max_abs_coef": max_abs_coef,
        "beta_age": float(beta_age),
        "residual_mean_reference": mean_ref,
        "residual_sd_reference": sd_ref,
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


# -----------------------------
# 5. Patch prediction files and model bundle
# -----------------------------

def patch_prediction_file(
    path: Path,
    prep,
    ridge,
    info,
    risk_col,
    z_col,
    yrs_col,
    age_col,
    instance_label,
    age_for_transform_candidates,
):
    if not path.exists():
        warnings.warn(f"Prediction file does not exist; skipping: {path}")
        return None

    backup_once(path, BACKUP_TAG + ".tsv")

    df = pd.read_csv(path, sep="\t")

    if risk_col not in df.columns:
        raise ValueError(f"{risk_col} not found in {path}")

    df = ensure_columns(df, COVARIATE_COLS)

    age_for_transform_col = first_existing_col(df, age_for_transform_candidates)
    if age_for_transform_col is None:
        raise ValueError(
            f"None of the age columns {age_for_transform_candidates} found in {path}"
        )

    # The residualizer expects age_at_baseline as the age covariate name.
    # For follow-up instances, age_at_baseline is overwritten with the sample age.
    df["age_at_baseline"] = pd.to_numeric(df[age_for_transform_col], errors="coerce")

    if "age_at_imaging" not in df.columns:
        df["age_at_imaging"] = df["age_at_baseline"]

    df_clean = clean_covariates(df)
    X = prep.transform(df_clean[COVARIATE_COLS])
    expected = ridge.predict(X)

    risk = pd.to_numeric(df[risk_col], errors="coerce").values
    resid_raw = risk - expected
    resid = resid_raw - float(info["risk_score_residual_mean_train"])

    sd_ref = float(info["risk_score_residual_sd_train"])
    beta_age = float(info["adjusted_age_coefficient_risk_score_per_year"])

    df[z_col] = resid / sd_ref
    df[yrs_col] = resid / beta_age
    df[age_col] = pd.to_numeric(df[age_for_transform_col], errors="coerce") + df[yrs_col]

    z_values = pd.to_numeric(df[z_col], errors="coerce")
    yrs_values = pd.to_numeric(df[yrs_col], errors="coerce")

    z_abs_max = float(np.nanmax(np.abs(z_values)))
    years_abs_max = float(np.nanmax(np.abs(yrs_values)))

    if z_abs_max > MAX_ABS_Z_ALLOWED or years_abs_max > MAX_ABS_YEARS_ALLOWED:
        raise ValueError(
            f"Patched values look invalid for {path}: max_abs_z={z_abs_max}, max_abs_years={years_abs_max}"
        )

    df.to_csv(path, sep="\t", index=False)

    summary = {
        "organ_label": ORGAN_LABEL,
        "organ": ORGAN,
        "instance_label": instance_label,
        "file": str(path),
        "n_rows": int(df.shape[0]),
        "age_for_transform_col": age_for_transform_col,
        "z_min": float(np.nanmin(z_values)),
        "z_median": float(np.nanmedian(z_values)),
        "z_mean": float(np.nanmean(z_values)),
        "z_max": float(np.nanmax(z_values)),
        "years_min": float(np.nanmin(yrs_values)),
        "years_median": float(np.nanmedian(yrs_values)),
        "years_mean": float(np.nanmean(yrs_values)),
        "years_max": float(np.nanmax(yrs_values)),
    }

    print("\nPatched:", path)
    print(json.dumps(summary, indent=2))
    return summary


def patch_model_joblib(model_joblib: Path, new_info, reference_split, ridge_alpha):
    backup_once(model_joblib, BACKUP_TAG + ".joblib")

    bundle = joblib.load(model_joblib)

    if "clock_transform_info_original_before_stable_backup" not in bundle:
        bundle["clock_transform_info_original_before_stable_backup"] = bundle.get(
            "clock_transform_info", None
        )

    bundle["clock_transform_info"] = new_info
    bundle["clock_transform_info_stabilized"] = True
    bundle["clock_transform_info_stabilized_organ"] = ORGAN
    bundle["clock_transform_info_stabilized_organ_label"] = ORGAN_LABEL
    bundle["clock_transform_info_stabilized_reference_split"] = reference_split
    bundle["clock_transform_info_stabilized_ridge_alpha"] = ridge_alpha

    joblib.dump(bundle, model_joblib)
    print("\nUpdated model joblib clock_transform_info:", model_joblib)


# -----------------------------
# 6. Main
# -----------------------------

def main():
    args = parse_args()
    root = Path(args.root)
    long_root = Path(args.long_root)
    reference_split = args.reference_split
    ridge_alpha = args.ridge_alpha

    (
        model_joblib,
        baseline_pred_file,
        followup_pred_files,
        debug_json,
        summary_json,
    ) = get_paths(root, long_root)

    risk_col, z_col, yrs_col, age_col = get_clock_cols()

    print("=" * 80)
    print("Stabilizing Pulmonary proteomics mortality-clock transform")
    print("=" * 80)
    print("Model joblib:", model_joblib)
    print("Baseline predictions:", baseline_pred_file)
    print("Follow-up predictions:")
    for inst, path in followup_pred_files.items():
        print(f"  {inst}: {path}")
    print("Reference split:", reference_split)
    print("Ridge alpha:", ridge_alpha)

    if not model_joblib.exists():
        raise FileNotFoundError(model_joblib)
    if not baseline_pred_file.exists():
        raise FileNotFoundError(baseline_pred_file)

    baseline_df = pd.read_csv(baseline_pred_file, sep="\t")

    prep, ridge, new_info, debug = fit_stable_transform_from_baseline(
        baseline_df=baseline_df,
        risk_col=risk_col,
        reference_split=reference_split,
        ridge_alpha=ridge_alpha,
    )

    print("\nStable residualizer debug:")
    print(json.dumps(debug, indent=2))

    patch_model_joblib(
        model_joblib=model_joblib,
        new_info=new_info,
        reference_split=reference_split,
        ridge_alpha=ridge_alpha,
    )

    results = {
        "organ_label": ORGAN_LABEL,
        "organ": ORGAN,
        "model_joblib": str(model_joblib),
        "baseline_pred_file": str(baseline_pred_file),
        "followup_pred_files": {k: str(v) for k, v in followup_pred_files.items()},
        "reference_split": reference_split,
        "ridge_alpha": ridge_alpha,
        "stable_transform_debug": debug,
        "patched_files": {},
    }

    results["patched_files"]["0_0"] = patch_prediction_file(
        path=baseline_pred_file,
        prep=prep,
        ridge=ridge,
        info=new_info,
        risk_col=risk_col,
        z_col=z_col,
        yrs_col=yrs_col,
        age_col=age_col,
        instance_label="0_0",
        age_for_transform_candidates=["age_at_baseline", "age_at_imaging"],
    )

    for instance_label, pred_file in followup_pred_files.items():
        if not pred_file.exists() and args.skip_missing_followup:
            warnings.warn(f"Skipping missing follow-up file: {pred_file}")
            results["patched_files"][instance_label] = None
            continue

        results["patched_files"][instance_label] = patch_prediction_file(
            path=pred_file,
            prep=prep,
            ridge=ridge,
            info=new_info,
            risk_col=risk_col,
            z_col=z_col,
            yrs_col=yrs_col,
            age_col=age_col,
            instance_label=instance_label,
            age_for_transform_candidates=["age_at_imaging", "age_at_baseline"],
        )

    with open(debug_json, "w") as f:
        json.dump(results, f, indent=2)

    with open(summary_json, "w") as f:
        json.dump(
            {
                "organ_label": ORGAN_LABEL,
                "organ": ORGAN,
                "reference_split": reference_split,
                "ridge_alpha": ridge_alpha,
                "patched_files": results["patched_files"],
                "debug_json": str(debug_json),
                "note": "Cox survival model was not retrained; only the post-hoc transform was stabilized.",
            },
            f,
            indent=2,
        )

    print("\nWrote debug JSON:", debug_json)
    print("Wrote summary JSON:", summary_json)
    print("\nFinished Pulmonary proteomics transform stabilization.")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        print("\nFAILED")
        traceback.print_exc()
        raise
