#!/usr/bin/env python3
# ============================================================
# Stabilize metabolomics mortality-clock transforms for
# Endocrine, Digestive, Hepatic, and Immune clocks.
#
# This DOES NOT retrain the Cox survival models.
#
# It only replaces the post-hoc residualization transform for:
#   {organ}_metabolomics_mortality_clock_acceleration_z
#   {organ}_metabolomics_mortality_clock_acceleration_years
#   {organ}_metabolomics_mortality_clock_age_years
#
# For each organ, it:
#   1. Loads baseline predictions.
#   2. Refits a stable residualizer on a baseline reference split.
#   3. Uses drop-first one-hot coding + small Ridge penalty.
#   4. Updates the model joblib clock_transform_info.
#   5. Patches baseline and instance-1 prediction TSVs.
#
# Default reference split:
#   REFERENCE_SPLIT = "train"
#
# To use the baseline test split for transform calibration, change:
#   REFERENCE_SPLIT = "test"
#
# Backups are created before overwriting.
# ============================================================

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
# 1. Global settings
# -----------------------------

ROOT = Path("/cbica/home/wenju/Reproducibile_paper/WholeBodyClock")

LONG_ROOT = Path(
    "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/"
    "mortality_clock/longitudinal/metabolomics"
)

ORGANS = {
    "Endocrine": "endocrine",
    "Digestive": "digestive",
    "Hepatic": "hepatic",
    "Immune": "immune",
}

# Choose "train", "test", "validation", "trainval", or "all".
# I keep "train" as default to match your current immune-fix script.
# Change to "test" if you want baseline-test-calibrated z-score transforms.
REFERENCE_SPLIT = "train"

RIDGE_ALPHA = 1e-4

BACKUP_TAG = "stable_metabolomics_clock_transform"

SUMMARY_JSON = (
    LONG_ROOT
    / "metabolomics_stable_clock_transform_all_organs_summary.json"
)

# Sanity thresholds.
MAX_ABS_Z_ALLOWED = 50
MAX_ABS_YEARS_ALLOWED = 200
MAX_ABS_COEF_ALLOWED = 1000


# -----------------------------
# 2. Shared covariate settings
# -----------------------------

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
# 3. Path and column helpers
# -----------------------------

def get_paths(organ_label, organ):
    model_joblib = (
        ROOT
        / f"{organ_label}_metabolomics_mortality_clock"
        / f"{organ}_metabolomics_mortality_clock_model.joblib"
    )

    baseline_pred_file = (
        ROOT
        / f"{organ_label}_metabolomics_mortality_clock"
        / f"{organ}_metabolomics_mortality_clock_predictions.tsv"
    )

    instance1_pred_file = (
        LONG_ROOT
        / organ_label
        / f"{organ}_metabolomics_mortality_clock_apply_instance_1_0_predictions.tsv"
    )

    out_debug_json = (
        ROOT
        / f"{organ_label}_metabolomics_mortality_clock"
        / f"{organ}_metabolomics_mortality_clock_stable_transform_debug.json"
    )

    return model_joblib, baseline_pred_file, instance1_pred_file, out_debug_json


def get_clock_cols(organ):
    risk_col = f"{organ}_metabolomics_mortality_risk_score"
    z_col = f"{organ}_metabolomics_mortality_clock_acceleration_z"
    yrs_col = f"{organ}_metabolomics_mortality_clock_acceleration_years"
    age_col = f"{organ}_metabolomics_mortality_clock_age_years"
    return risk_col, z_col, yrs_col, age_col


# -----------------------------
# 4. General helpers
# -----------------------------

def make_onehot_drop_first():
    """
    drop='first' avoids the dummy-variable trap when an intercept is included.
    handle_unknown='ignore' allows instance-1 categories not seen in baseline reference split.
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
    df = ensure_columns(df, COVARIATE_COLS).copy()

    for c in NUMERIC_COVS:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    for c in CATEGORICAL_COVS:
        df[c] = df[c].astype("object")
        df.loc[df[c].isna(), c] = np.nan

    return df


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
            ("num", num_pipe, NUMERIC_COVS),
            ("cat", cat_pipe, CATEGORICAL_COVS),
        ],
        remainder="drop",
    )

    return prep


def get_reference_mask(baseline_df, reference_split):
    reference_split = str(reference_split).lower()

    if reference_split == "all":
        return np.repeat(True, baseline_df.shape[0])

    if "split" not in baseline_df.columns:
        warnings.warn(
            "split column not found in baseline file. "
            "Using all baseline rows to fit residualizer."
        )
        return np.repeat(True, baseline_df.shape[0])

    split_values = baseline_df["split"].astype(str).str.lower()

    if reference_split == "trainval":
        mask = split_values.isin(["train", "validation"]).values
    else:
        mask = split_values.eq(reference_split).values

    if mask.sum() < 10:
        raise ValueError(
            f"Too few rows for REFERENCE_SPLIT={reference_split}: {mask.sum()}"
        )

    return mask


# -----------------------------
# 5. Fit stable transform
# -----------------------------

def fit_stable_transform_from_baseline(
    baseline_df,
    organ_label,
    organ,
    risk_col,
    reference_split=REFERENCE_SPLIT,
):
    if risk_col not in baseline_df.columns:
        raise ValueError(f"{risk_col} not found in baseline prediction file.")

    baseline_df = clean_covariates(baseline_df)
    ref_mask = get_reference_mask(baseline_df, reference_split)
    ref_df = baseline_df.loc[ref_mask].copy()

    X_ref_raw = ref_df[COVARIATE_COLS].copy()
    X_all_raw = baseline_df[COVARIATE_COLS].copy()

    prep = build_preprocessor()
    X_ref = prep.fit_transform(X_ref_raw)
    X_all = prep.transform(X_all_raw)

    ridge = Ridge(alpha=RIDGE_ALPHA, fit_intercept=True)
    ridge.fit(X_ref, pd.to_numeric(ref_df[risk_col], errors="coerce").values)

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
            f"Stable residualizer still has very large coefficients: "
            f"max_abs_coef={max_abs_coef}"
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
        "ridge_alpha": RIDGE_ALPHA,
        "note": (
            f"Stable {organ_label} metabolomics clock transform refit from baseline "
            f"{reference_split} split predictions. Uses drop-first one-hot encoding "
            f"plus Ridge(alpha={RIDGE_ALPHA}) to avoid the dummy-variable trap and "
            "numerical instability. Cox model was not retrained."
        ),
    }

    debug = {
        "organ_label": organ_label,
        "organ": organ,
        "n_baseline_rows": int(baseline_df.shape[0]),
        "reference_split": reference_split,
        "n_reference_rows": int(ref_mask.sum()),
        "ridge_alpha": RIDGE_ALPHA,
        "intercept": float(ridge.intercept_),
        "max_abs_coef": max_abs_coef,
        "beta_age": float(beta_age),
        "residual_mean_reference": mean_ref,
        "residual_sd_reference": sd_ref,
        "top_abs_coefficients": sorted(
            [
                {
                    "feature": k,
                    "coefficient": float(v),
                    "abs_coefficient": abs(float(v)),
                }
                for k, v in coef_dict.items()
            ],
            key=lambda d: d["abs_coefficient"],
            reverse=True,
        )[:30],
    }

    return prep, ridge, info, debug


# -----------------------------
# 6. Patch prediction file
# -----------------------------

def patch_prediction_file(
    path,
    prep,
    ridge,
    info,
    organ_label,
    organ,
    risk_col,
    z_col,
    yrs_col,
    age_col,
    age_for_transform_col,
    age_for_clock_col,
):
    if not path.exists():
        warnings.warn(f"Prediction file does not exist; skipping: {path}")
        return None

    backup_once(path, BACKUP_TAG + ".tsv")

    df = pd.read_csv(path, sep="\t")

    if risk_col not in df.columns:
        raise ValueError(f"{risk_col} not found in {path}")

    df = ensure_columns(df, COVARIATE_COLS)

    if age_for_transform_col not in df.columns:
        raise ValueError(f"{age_for_transform_col} not found in {path}")

    # The residualizer always expects a column named age_at_baseline.
    # For longitudinal instance 1, we fill it with age_at_imaging to reflect sample age.
    df["age_at_baseline"] = pd.to_numeric(df[age_for_transform_col], errors="coerce")

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

    if age_for_clock_col in df.columns:
        clock_base_age = pd.to_numeric(df[age_for_clock_col], errors="coerce")
    else:
        clock_base_age = pd.to_numeric(df[age_for_transform_col], errors="coerce")

    df[age_col] = clock_base_age + df[yrs_col]

    z_values = pd.to_numeric(df[z_col], errors="coerce")
    yrs_values = pd.to_numeric(df[yrs_col], errors="coerce")

    z_abs_max = float(np.nanmax(np.abs(z_values)))
    years_abs_max = float(np.nanmax(np.abs(yrs_values)))

    if z_abs_max > MAX_ABS_Z_ALLOWED or years_abs_max > MAX_ABS_YEARS_ALLOWED:
        raise ValueError(
            f"Patched values still look invalid for {path}: "
            f"max_abs_z={z_abs_max}, max_abs_years={years_abs_max}"
        )

    df.to_csv(path, sep="\t", index=False)

    summary = {
        "organ_label": organ_label,
        "organ": organ,
        "file": str(path),
        "n_rows": int(df.shape[0]),
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


# -----------------------------
# 7. Patch model joblib
# -----------------------------

def patch_model_joblib(model_joblib, organ_label, organ, new_info):
    backup_once(model_joblib, BACKUP_TAG + ".joblib")

    bundle = joblib.load(model_joblib)

    # Preserve the first original transform if not already preserved.
    if "clock_transform_info_original_before_stable_backup" not in bundle:
        bundle["clock_transform_info_original_before_stable_backup"] = bundle.get(
            "clock_transform_info",
            None,
        )

    # Also preserve immune-specific old backup if it exists from previous script.
    if (
        "clock_transform_info_original_unstable_backup" in bundle
        and "clock_transform_info_original_before_stable_backup" not in bundle
    ):
        bundle["clock_transform_info_original_before_stable_backup"] = bundle[
            "clock_transform_info_original_unstable_backup"
        ]

    bundle["clock_transform_info"] = new_info
    bundle["clock_transform_info_stabilized"] = True
    bundle["clock_transform_info_stabilized_organ"] = organ
    bundle["clock_transform_info_stabilized_organ_label"] = organ_label
    bundle["clock_transform_info_stabilized_reference_split"] = REFERENCE_SPLIT
    bundle["clock_transform_info_stabilized_ridge_alpha"] = RIDGE_ALPHA

    joblib.dump(bundle, model_joblib)

    print("\nUpdated model joblib clock_transform_info:")
    print(model_joblib)


# -----------------------------
# 8. Run one organ
# -----------------------------

def run_one_organ(organ_label, organ):
    print("\n" + "=" * 80)
    print(f"Stabilizing metabolomics mortality-clock transform: {organ_label}")
    print("=" * 80)

    (
        model_joblib,
        baseline_pred_file,
        instance1_pred_file,
        out_debug_json,
    ) = get_paths(organ_label, organ)

    risk_col, z_col, yrs_col, age_col = get_clock_cols(organ)

    print("Model joblib:", model_joblib)
    print("Baseline predictions:", baseline_pred_file)
    print("Instance-1 predictions:", instance1_pred_file)
    print("Reference split:", REFERENCE_SPLIT)

    if not model_joblib.exists():
        raise FileNotFoundError(model_joblib)

    if not baseline_pred_file.exists():
        raise FileNotFoundError(baseline_pred_file)

    if not instance1_pred_file.exists():
        warnings.warn(f"Instance-1 prediction file does not exist: {instance1_pred_file}")

    baseline_df = pd.read_csv(baseline_pred_file, sep="\t")

    prep, ridge, new_info, debug = fit_stable_transform_from_baseline(
        baseline_df=baseline_df,
        organ_label=organ_label,
        organ=organ,
        risk_col=risk_col,
        reference_split=REFERENCE_SPLIT,
    )

    print("\nStable residualizer debug:")
    print(json.dumps(debug, indent=2))

    patch_model_joblib(
        model_joblib=model_joblib,
        organ_label=organ_label,
        organ=organ,
        new_info=new_info,
    )

    # Baseline file: sample age is age_at_baseline.
    baseline_summary = patch_prediction_file(
        path=baseline_pred_file,
        prep=prep,
        ridge=ridge,
        info=new_info,
        organ_label=organ_label,
        organ=organ,
        risk_col=risk_col,
        z_col=z_col,
        yrs_col=yrs_col,
        age_col=age_col,
        age_for_transform_col="age_at_baseline",
        age_for_clock_col="age_at_baseline",
    )

    # Instance 1 file: sample age is age_at_imaging.
    instance1_summary = patch_prediction_file(
        path=instance1_pred_file,
        prep=prep,
        ridge=ridge,
        info=new_info,
        organ_label=organ_label,
        organ=organ,
        risk_col=risk_col,
        z_col=z_col,
        yrs_col=yrs_col,
        age_col=age_col,
        age_for_transform_col="age_at_imaging",
        age_for_clock_col="age_at_imaging",
    )

    debug_out = {
        "organ_label": organ_label,
        "organ": organ,
        "model_joblib": str(model_joblib),
        "baseline_pred_file": str(baseline_pred_file),
        "instance1_pred_file": str(instance1_pred_file),
        "reference_split": REFERENCE_SPLIT,
        "stable_transform_debug": debug,
        "baseline_summary": baseline_summary,
        "instance1_summary": instance1_summary,
    }

    with open(out_debug_json, "w") as f:
        json.dump(debug_out, f, indent=2)

    print("\nWrote organ debug JSON:")
    print(out_debug_json)

    print(
        f"\nFinished {organ_label}. Cox model was not retrained; "
        "only the post-hoc clock transform was stabilized."
    )

    return debug_out


# -----------------------------
# 9. Main
# -----------------------------

def main():
    print("=" * 80)
    print("Stabilizing metabolomics mortality-clock transforms for all organs")
    print("=" * 80)
    print("Organs:", ", ".join(ORGANS.keys()))
    print("Reference split:", REFERENCE_SPLIT)
    print("Ridge alpha:", RIDGE_ALPHA)

    all_results = {}
    failures = {}

    for organ_label, organ in ORGANS.items():
        try:
            all_results[organ_label] = run_one_organ(organ_label, organ)
        except Exception as e:
            print("\n" + "!" * 80)
            print(f"FAILED: {organ_label}")
            print("!" * 80)
            print(str(e))
            traceback.print_exc()

            failures[organ_label] = {
                "organ": organ,
                "error": str(e),
                "traceback": traceback.format_exc(),
            }

    summary = {
        "reference_split": REFERENCE_SPLIT,
        "ridge_alpha": RIDGE_ALPHA,
        "successful_organs": list(all_results.keys()),
        "failed_organs": failures,
        "results": all_results,
    }

    with open(SUMMARY_JSON, "w") as f:
        json.dump(summary, f, indent=2)

    print("\n" + "=" * 80)
    print("Finished metabolomics clock-transform stabilization.")
    print("Summary JSON:")
    print(SUMMARY_JSON)
    print("Successful organs:", ", ".join(summary["successful_organs"]))
    if failures:
        print("Failed organs:", ", ".join(failures.keys()))
    else:
        print("Failed organs: none")
    print("=" * 80)


if __name__ == "__main__":
    main()