#!/usr/bin/env python3
# ============================================================
# Apply pre-trained organ MRI mortality clocks to UKB instance 3
#
# Organs currently intended:
#   heart
#   pancreas
#
# This is application-only:
#   - no Cox model refitting
#   - no train/validation/test split
#   - loads saved *_mri_mortality_clock_model.joblib
#   - applies saved preprocessor and saved Cox model
#   - computes risk score, absolute risk, acceleration_z,
#     acceleration_years, and clock_age_years using the saved
#     clock_transform_info from the training script
#
# Final clean output:
#   participant_id
#   application_instance
#   application_source_file
#   {organ}_mri_mortality_risk_score
#   sample_date
#   death_date
#   age_at_baseline
#   age_at_imaging
#   sex
#   bmi_at_imaging
#   diastolic_bp_at_imaging
#   systolic_bp_at_imaging
#   uk_biobank_assessment_centre_f54_2_0
#   risk_5y
#   risk_10y
#   risk_15y
#   {organ}_mri_mortality_clock_acceleration_z
#   {organ}_mri_mortality_clock_acceleration_years
#   {organ}_mri_mortality_clock_age_years
#   n_model_mri_features_expected
#   n_model_mri_features_present_in_input
#   n_model_mri_features_missing_from_input
#
# IMPORTANT:
#   Original MRI features are used for model prediction but are NOT
#   written to the final prediction file.
# ============================================================

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


# -----------------------------
# 1. Arguments
# -----------------------------

def parse_args():
    p = argparse.ArgumentParser()

    p.add_argument("--organ", required=True, help="Organ name, e.g., heart or pancreas.")
    p.add_argument("--model-joblib", required=True, help="Path to pretrained *_mri_mortality_clock_model.joblib.")

    p.add_argument(
        "--input-tsv",
        action="append",
        required=True,
        help=(
            "Application TSV. Use either PATH or LABEL:PATH. "
            "Example: 3_0:/path/imaging_heart_3_0.tsv. "
            "Can be supplied multiple times."
        ),
    )

    p.add_argument(
        "--covariate-csv",
        default="/cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv",
        help="UKB covariate CSV with eid and instance-specific covariates.",
    )

    p.add_argument(
        "--death-xlsx",
        default="/cbica/home/wenju/Dataset/UKBB_UMelbourne/Death_related_var_from_Ye.xlsx",
        help="Death/assessment-date file.",
    )

    p.add_argument(
        "--id-match-csv",
        default="/cbica/home/wenju/Dataset/UKBB_UMelbourne/UKB_UMelbourne_vs_Penn_match_key.csv",
        help="ID mapping file from UMelbourne ID to Penn ID.",
    )

    p.add_argument("--admin-censor-date", default="2022-11-30")
    p.add_argument("--outdir", required=True)
    p.add_argument("--application-session-id", default="ses-M3", help="Use 'none' to disable session filtering.")
    p.add_argument("--imaging-instance", default="3")
    p.add_argument("--model-instance", default="2")
    p.add_argument("--feature-start-column", default="diagnosis")
    p.add_argument("--risk-times", default="5,10,15")

    p.add_argument(
        "--complete-case-organ-features",
        action="store_true",
        help="Drop participants missing any model-used organ MRI feature.",
    )

    p.add_argument(
        "--allow-missing-model-columns",
        action="store_true",
        help="Create missing saved model columns as NA and let the saved preprocessor impute them.",
    )

    # Kept only for backward compatibility with old Bash commands.
    # This script intentionally ignores it and always writes clean output.
    p.add_argument(
        "--include-features-in-output",
        action="store_true",
        help="Deprecated and ignored. Final output never includes original MRI features.",
    )

    return p.parse_args()


# -----------------------------
# 2. Generic helpers
# -----------------------------

def clean_name(x):
    x = re.sub(r"[^A-Za-z0-9]+", "_", str(x).strip().lower())
    x = re.sub(r"_+", "_", x).strip("_")
    if not x:
        raise ValueError("--organ is empty after sanitization.")
    return x


def output_prefix(organ):
    return f"{organ}_mri_mortality_clock"


def parse_risk_times(s):
    vals = []
    for x in str(s).split(","):
        x = x.strip()
        if x:
            vals.append(float(x))
    if not vals:
        raise ValueError("--risk-times cannot be empty.")
    return vals


def parse_input_spec(spec):
    spec = str(spec)
    if ":" in spec and not spec.startswith("/"):
        label, path = spec.split(":", 1)
        return label.strip(), path.strip()
    path = spec.strip()
    return Path(path).stem, path


def expand_paths(path):
    if any(ch in path for ch in ["*", "?", "["]):
        paths = sorted(glob.glob(path))
        if not paths:
            raise FileNotFoundError(f"No files matched: {path}")
        return paths
    return [path]


def as_numeric(s):
    return pd.to_numeric(s, errors="coerce")


def normalize_sex(s):
    return (
        s.astype(str)
        .str.strip()
        .replace({
            "0": "Female",
            "0.0": "Female",
            "1": "Male",
            "1.0": "Male",
            "F": "Female",
            "M": "Male",
            "female": "Female",
            "male": "Male",
            "Female": "Female",
            "Male": "Male",
        })
    )


def mean_existing_numeric_columns(df, cols, out_col):
    present = [c for c in cols if c in df.columns]
    if not present:
        return df

    for c in present:
        df[c] = as_numeric(df[c])

    df[out_col] = df[present].mean(axis=1, skipna=True)
    return df


def pick_first_existing(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None


# -----------------------------
# 3. Load application data
# -----------------------------

def load_application_tsv(input_tsv, organ, session_id):
    frames = []

    for spec in input_tsv:
        label, path_spec = parse_input_spec(spec)

        for path in expand_paths(path_spec):
            if not os.path.exists(path):
                raise FileNotFoundError(f"Application TSV not found: {path}")

            part = pd.read_csv(path, sep="\t")

            if "participant_id" not in part.columns:
                raise ValueError(f"participant_id is missing from {path}")

            part["application_instance"] = label
            part["application_source_file"] = Path(path).name
            frames.append(part)

            print(f"Loaded {organ} application TSV: {path}; rows={part.shape[0]}, cols={part.shape[1]}")

    df = pd.concat(frames, axis=0, ignore_index=True, sort=False)

    if session_id and str(session_id).lower() not in {"none", "null", ""}:
        if "session_id" in df.columns:
            before = df.shape[0]
            df = df.loc[df["session_id"].astype(str) == str(session_id)].copy()
            print(f"Filtered to session_id={session_id}: {before} -> {df.shape[0]} rows")
        else:
            warnings.warn("--application-session-id was provided, but session_id is missing.")

    df["participant_id"] = as_numeric(df["participant_id"])
    df = df.dropna(subset=["participant_id"]).copy()
    df["participant_id"] = df["participant_id"].astype(int)

    df = df.sort_values(["participant_id", "application_instance"], kind="mergesort")

    dup = int(df.duplicated(["participant_id", "application_instance"]).sum())
    if dup > 0:
        warnings.warn(f"Found {dup} duplicated participant_id/application_instance rows. Keeping first.")
        df = df.drop_duplicates(["participant_id", "application_instance"], keep="first")

    print("Application dataset:")
    print(f"  Rows = {df.shape[0]}")
    print(f"  Columns = {df.shape[1]}")

    return df


# -----------------------------
# 4. Load death and covariates
# -----------------------------

def load_death_data(death_xlsx, id_match_csv, imaging_instance="3", admin_censor_date="2022-11-30"):
    if not os.path.exists(death_xlsx):
        warnings.warn(f"Death file not found: {death_xlsx}. Continuing without date annotations.")
        return None

    if not os.path.exists(id_match_csv):
        warnings.warn(f"ID match file not found: {id_match_csv}. Continuing without date annotations.")
        return None

    d = pd.read_excel(death_xlsx)
    m = pd.read_csv(id_match_csv)

    if "eid" not in d.columns:
        warnings.warn("Death file does not contain eid. Continuing without date annotations.")
        return None

    d = d.rename(columns={"eid": "participant_id_umel"})

    if "id" not in m.columns or "id_upenn" not in m.columns:
        warnings.warn("ID match file must contain id and id_upenn. Continuing without date annotations.")
        return None

    m = m.rename(columns={"id": "participant_id_umel", "id_upenn": "participant_id"})
    d = m.merge(d, on="participant_id_umel", how="inner")

    baseline_date_col = "53-0.0"
    model_imaging_date_col = "53-2.0"
    target_imaging_date_col = f"53-{imaging_instance}.0"
    death_date_col = "40000-0.0"

    keep = ["participant_id"]
    for c in [baseline_date_col, model_imaging_date_col, target_imaging_date_col, death_date_col]:
        if c in d.columns:
            keep.append(c)
        else:
            warnings.warn(f"Death/assessment file is missing column: {c}")

    d = d[keep].copy()

    d["baseline_date"] = pd.to_datetime(d[baseline_date_col], errors="coerce") if baseline_date_col in d.columns else pd.NaT
    d["model_imaging_date"] = pd.to_datetime(d[model_imaging_date_col], errors="coerce") if model_imaging_date_col in d.columns else pd.NaT
    d["imaging_date"] = pd.to_datetime(d[target_imaging_date_col], errors="coerce") if target_imaging_date_col in d.columns else pd.NaT
    d["death_date"] = pd.to_datetime(d[death_date_col], errors="coerce") if death_date_col in d.columns else pd.NaT

    d["admin_censor_date"] = pd.to_datetime(admin_censor_date)

    d["event"] = (
        d["death_date"].notna()
        & d["imaging_date"].notna()
        & (d["death_date"] > d["imaging_date"])
        & (d["death_date"] <= d["admin_censor_date"])
    )

    d["end_date"] = d["admin_censor_date"]
    d.loc[d["event"], "end_date"] = d.loc[d["event"], "death_date"]

    d["time_days"] = (d["end_date"] - d["imaging_date"]).dt.days
    d["time_years"] = d["time_days"] / 365.25

    out_cols = [
        "participant_id",
        "baseline_date",
        "model_imaging_date",
        "imaging_date",
        "death_date",
        "admin_censor_date",
        "end_date",
        "event",
        "time_days",
        "time_years",
    ]

    return d[out_cols].copy()


def load_covariates(path):
    if path is None or str(path).lower() in {"none", ""}:
        return None

    if not os.path.exists(path):
        warnings.warn(f"Covariate file not found: {path}. Continuing without it.")
        return None

    cov = pd.read_csv(path)

    if "eid" in cov.columns:
        cov = cov.rename(columns={"eid": "participant_id"})
    elif "participant_id" not in cov.columns:
        warnings.warn("Covariate file does not contain eid or participant_id. Continuing without it.")
        return None

    cov["participant_id"] = as_numeric(cov["participant_id"])
    cov = cov.dropna(subset=["participant_id"]).copy()
    cov["participant_id"] = cov["participant_id"].astype(int)

    dup = int(cov["participant_id"].duplicated().sum())
    if dup > 0:
        warnings.warn(f"Found {dup} duplicated participant_id rows in covariate file. Keeping first.")
        cov = cov.drop_duplicates("participant_id", keep="first")

    return cov


def add_application_covariates(df, imaging_instance="3", model_instance="2"):
    df = df.copy()

    baseline_age_source = pick_first_existing(
        df,
        [
            "age_when_attended_assessment_centre_f21003_0_0",
            "age_at_baseline",
            "age_at_recruitment_f21022_0_0",
        ],
    )
    if baseline_age_source is not None:
        df["age_at_baseline"] = as_numeric(df[baseline_age_source])
    else:
        df["age_at_baseline"] = np.nan

    age_source = pick_first_existing(
        df,
        [
            f"age_when_attended_assessment_centre_f21003_{imaging_instance}_0",
            f"age_when_attended_assessment_centre_f21003_{model_instance}_0",
            "age_at_imaging",
            "diagnosis",
        ],
    )

    if age_source is None:
        raise ValueError(
            "Could not infer age_at_imaging. Expected f21003 instance-3/2 field, "
            "age_at_imaging, or numeric diagnosis fallback."
        )

    df["age_at_imaging"] = as_numeric(df[age_source])

    sex_source = pick_first_existing(df, ["sex_f31_0_0", "sex", "Sex"])
    if sex_source is None:
        raise ValueError("Could not infer sex. Expected sex_f31_0_0, sex, or Sex.")
    df["sex"] = normalize_sex(df[sex_source])

    bmi_source = pick_first_existing(
        df,
        [
            f"body_mass_index_bmi_f23104_{imaging_instance}_0",
            f"body_mass_index_bmi_f23104_{model_instance}_0",
            "bmi_at_imaging",
        ],
    )
    if bmi_source is not None:
        df["bmi_at_imaging"] = as_numeric(df[bmi_source])
    elif "bmi_at_imaging" not in df.columns:
        df["bmi_at_imaging"] = np.nan

    df = mean_existing_numeric_columns(
        df,
        [
            f"diastolic_blood_pressure_automated_reading_f4079_{imaging_instance}_0",
            f"diastolic_blood_pressure_automated_reading_f4079_{imaging_instance}_1",
            f"diastolic_blood_pressure_automated_reading_f4079_{model_instance}_0",
            f"diastolic_blood_pressure_automated_reading_f4079_{model_instance}_1",
            "diastolic_bp_at_imaging",
        ],
        "diastolic_bp_at_imaging",
    )

    if "diastolic_bp_at_imaging" not in df.columns:
        df["diastolic_bp_at_imaging"] = np.nan

    df = mean_existing_numeric_columns(
        df,
        [
            f"systolic_blood_pressure_automated_reading_f4080_{imaging_instance}_0",
            f"systolic_blood_pressure_automated_reading_f4080_{imaging_instance}_1",
            f"systolic_blood_pressure_automated_reading_f4080_{model_instance}_0",
            f"systolic_blood_pressure_automated_reading_f4080_{model_instance}_1",
            "systolic_bp_at_imaging",
        ],
        "systolic_bp_at_imaging",
    )

    if "systolic_bp_at_imaging" not in df.columns:
        df["systolic_bp_at_imaging"] = np.nan

    model_center_name = f"uk_biobank_assessment_centre_f54_{model_instance}_0"
    center_source = pick_first_existing(
        df,
        [
            f"uk_biobank_assessment_centre_f54_{imaging_instance}_0",
            f"uk_biobank_assessment_centre_f54_{model_instance}_0",
            model_center_name,
        ],
    )

    if center_source is not None:
        df[model_center_name] = df[center_source].astype("object")
    elif model_center_name not in df.columns:
        df[model_center_name] = pd.NA

    return df


# -----------------------------
# 5. Model application
# -----------------------------

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


def numeric_fill_for_residualization(series, name):
    x = as_numeric(series)

    if x.isna().all():
        warnings.warn(
            f"Residualization numeric covariate {name} is all missing. "
            "Using 0 for expected-risk reconstruction."
        )
        return pd.Series(np.repeat(0.0, len(x)), index=x.index)

    med = float(np.nanmedian(x))
    return x.fillna(med)


def categorical_match(series, category_value):
    s_obj = series.astype("object")
    s_str = s_obj.astype(str).str.strip()
    cat_str = str(category_value).strip()

    mask = s_str == cat_str

    try:
        cat_float = float(cat_str)
        s_num = pd.to_numeric(s_obj, errors="coerce")
        mask = mask | np.isclose(s_num.astype(float), cat_float, equal_nan=False)
    except Exception:
        pass

    return mask.astype(float)


def parse_categorical_coef_name(term, categorical_covs):
    if not term.startswith("cat__"):
        return None, None

    stem = term[len("cat__"):]

    for cov in sorted(categorical_covs, key=len, reverse=True):
        prefix = f"{cov}_"
        if stem.startswith(prefix):
            category = stem[len(prefix):]
            return cov, category

    return None, None


def compute_clock_transforms_from_saved_info(df, risk, organ, clock_info):
    risk_col = f"{organ}_mri_mortality_risk_score"
    z_col = f"{organ}_mri_mortality_clock_acceleration_z"
    yrs_col = f"{organ}_mri_mortality_clock_acceleration_years"
    age_col = f"{organ}_mri_mortality_clock_age_years"

    out = df.copy()
    out[risk_col] = risk

    if not clock_info:
        warnings.warn("clock_transform_info is missing from model bundle. Clock columns set to NA.")
        out[z_col] = np.nan
        out[yrs_col] = np.nan
        out[age_col] = np.nan
        return out

    intercept = float(clock_info.get("risk_score_covariate_model_intercept", 0.0))
    coef_dict = clock_info.get("risk_score_covariate_model_coef", {})

    covariates = list(clock_info.get("residualization_covariates", []))
    numeric_covs = list(clock_info.get("numeric_residualization_covariates", []))
    categorical_covs = list(clock_info.get("categorical_residualization_covariates", []))

    if not numeric_covs and not categorical_covs:
        for c in covariates:
            if c in out.columns and c != "sex" and pd.api.types.is_numeric_dtype(out[c]):
                numeric_covs.append(c)
            else:
                categorical_covs.append(c)

    expected = np.repeat(intercept, out.shape[0]).astype(float)

    for term, beta in coef_dict.items():
        beta = float(beta)

        if term.startswith("num__"):
            cov = term[len("num__"):]
            if cov not in out.columns:
                warnings.warn(f"Residualization numeric covariate missing: {cov}. Using 0.")
                vals = np.repeat(0.0, out.shape[0])
            else:
                vals = numeric_fill_for_residualization(out[cov], cov).values.astype(float)
            expected += beta * vals

        elif term.startswith("cat__"):
            cov, category = parse_categorical_coef_name(term, categorical_covs)

            if cov is None or category is None:
                warnings.warn(f"Could not parse categorical residualization term: {term}. Ignoring.")
                continue

            if cov not in out.columns:
                warnings.warn(f"Residualization categorical covariate missing: {cov}. Term set to 0.")
                vals = np.repeat(0.0, out.shape[0])
            else:
                vals = categorical_match(out[cov], category).values.astype(float)

            expected += beta * vals

        else:
            warnings.warn(f"Unknown residualization coefficient term: {term}. Ignoring.")

    mean_train = float(clock_info.get("risk_score_residual_mean_train", 0.0))
    sd_train = float(clock_info.get("risk_score_residual_sd_train", np.nan))
    beta_age = clock_info.get("adjusted_age_coefficient_risk_score_per_year", None)

    resid_raw = risk - expected
    resid = resid_raw - mean_train

    if np.isfinite(sd_train) and sd_train > 0:
        out[z_col] = resid / sd_train
    else:
        warnings.warn("Training residual SD is missing or zero. acceleration_z set to NA.")
        out[z_col] = np.nan

    if beta_age is not None and np.isfinite(float(beta_age)) and abs(float(beta_age)) > 1e-8:
        beta_age = float(beta_age)
        out[yrs_col] = resid / beta_age

        if "age_at_imaging" in out.columns:
            out[age_col] = as_numeric(out["age_at_imaging"]) + out[yrs_col]
        else:
            out[age_col] = np.nan
    else:
        warnings.warn("Adjusted age coefficient is missing or near zero. Year-scale acceleration set to NA.")
        out[yrs_col] = np.nan
        out[age_col] = np.nan

    return out


def prepare_model_input_columns(df, numeric_cols_kept, categorical_cols_kept):
    df = df.copy()

    expected_cols = list(numeric_cols_kept) + list(categorical_cols_kept)

    for c in expected_cols:
        if c not in df.columns:
            warnings.warn(f"Saved model expects column missing in application data: {c}. Creating as NA.")
            df[c] = np.nan

    for c in numeric_cols_kept:
        df[c] = as_numeric(df[c])

    for c in categorical_cols_kept:
        df[c] = df[c].astype("object")

    return df, expected_cols


def validate_clean_output_columns(pred_out, organ, model_used_organ_features):
    required = [
        f"{organ}_mri_mortality_risk_score",
        f"{organ}_mri_mortality_clock_acceleration_z",
        f"{organ}_mri_mortality_clock_acceleration_years",
        f"{organ}_mri_mortality_clock_age_years",
    ]

    missing = [c for c in required if c not in pred_out.columns]
    if missing:
        raise RuntimeError(f"Missing required output columns for {organ}: {missing}")

    malformed = [
        c for c in pred_out.columns
        if c.count(f"{organ}_mri_mortality_clock_acceleration_") > 1
    ]

    if malformed:
        raise RuntimeError(f"Malformed clock acceleration columns found for {organ}: {malformed}")

    leaked_features = [c for c in model_used_organ_features if c in pred_out.columns]
    if leaked_features:
        raise RuntimeError(
            f"Original MRI feature columns were accidentally written to final output for {organ}. "
            f"Examples: {leaked_features[:10]}"
        )

    return True


def apply_one_instance(
    df_app,
    bundle,
    organ,
    outdir,
    risk_times,
    complete_case_organ_features=False,
    allow_missing_model_columns=False,
):
    pref = output_prefix(organ)

    preprocessor = bundle["preprocessor"]
    model = bundle["model"]
    numeric_cols_kept = list(bundle["numeric_cols_kept"])
    categorical_cols_kept = list(bundle["categorical_cols_kept"])
    organ_feature_cols_original = list(bundle.get("organ_feature_cols", []))
    clock_info = bundle.get("clock_transform_info", None)

    model_used_organ_features = [
        c for c in organ_feature_cols_original
        if c in numeric_cols_kept
    ]

    if not model_used_organ_features:
        cov_like = {
            "age_at_imaging",
            "bmi_at_imaging",
            "diastolic_bp_at_imaging",
            "systolic_bp_at_imaging",
        }
        model_used_organ_features = [
            c for c in numeric_cols_kept
            if c not in cov_like
        ]

    if not model_used_organ_features:
        raise RuntimeError(f"No model-used organ MRI features could be identified for {organ}.")

    present_organ_features = [c for c in model_used_organ_features if c in df_app.columns]
    missing_organ_features = [c for c in model_used_organ_features if c not in df_app.columns]

    print("Saved model input summary:")
    print(f"  Numeric columns kept = {len(numeric_cols_kept)}")
    print(f"  Categorical columns kept = {len(categorical_cols_kept)}")
    print(f"  Original organ features in bundle = {len(organ_feature_cols_original)}")
    print(f"  Model-used organ features = {len(model_used_organ_features)}")
    print(f"  Present model-used organ features in application TSV = {len(present_organ_features)}")
    print(f"  Missing model-used organ features in application TSV = {len(missing_organ_features)}")

    if missing_organ_features:
        for c in missing_organ_features[:20]:
            print(f"  missing organ feature: {c}")
        if len(missing_organ_features) > 20:
            print("  ...")

        if not allow_missing_model_columns:
            raise ValueError(
                f"{len(missing_organ_features)} model-used organ features are missing from application data. "
                "Use --allow-missing-model-columns only if you want saved preprocessor imputation."
            )

    df_model, expected_cols = prepare_model_input_columns(
        df_app,
        numeric_cols_kept=numeric_cols_kept,
        categorical_cols_kept=categorical_cols_kept,
    )

    n_before = df_model.shape[0]

    df_model["n_model_mri_features_expected"] = len(model_used_organ_features)
    df_model["n_model_mri_features_missing_from_input"] = df_model[model_used_organ_features].isna().sum(axis=1)
    df_model["n_model_mri_features_present_in_input"] = (
        df_model["n_model_mri_features_expected"]
        - df_model["n_model_mri_features_missing_from_input"]
    )

    if complete_case_organ_features:
        before = df_model.shape[0]
        df_model = df_model.loc[df_model["n_model_mri_features_missing_from_input"] == 0].copy()
        print(f"Complete-case organ-feature filter: {before} -> {df_model.shape[0]} participants")

    if df_model.empty:
        raise ValueError("No participants remain after complete-case filtering.")

    X_raw = df_model[numeric_cols_kept + categorical_cols_kept].copy()
    X = preprocessor.transform(X_raw)

    risk = np.asarray(model.predict(X)).reshape(-1)

    pred = compute_clock_transforms_from_saved_info(
        df=df_model,
        risk=risk,
        organ=organ,
        clock_info=clock_info,
    )

    abs_risk = predict_absolute_risk(model, X, risk_times)
    pred = pd.concat([pred.reset_index(drop=True), abs_risk.reset_index(drop=True)], axis=1)

    risk_col = f"{organ}_mri_mortality_risk_score"
    z_col = f"{organ}_mri_mortality_clock_acceleration_z"
    yrs_col = f"{organ}_mri_mortality_clock_acceleration_years"
    age_col = f"{organ}_mri_mortality_clock_age_years"

    pred["sample_date"] = pred["imaging_date"] if "imaging_date" in pred.columns else pd.NaT

    id_cols = [
        "participant_id",
        "application_instance",
        "application_source_file",
        risk_col,
        "sample_date",
        "death_date",
        "age_at_baseline",
        "age_at_imaging",
        "sex",
        "bmi_at_imaging",
        "diastolic_bp_at_imaging",
        "systolic_bp_at_imaging",
        "uk_biobank_assessment_centre_f54_2_0",
    ]

    risk_cols = [f"risk_{t:g}y" for t in risk_times]

    clock_cols = [
        z_col,
        yrs_col,
        age_col,
    ]

    qc_cols = [
        "n_model_mri_features_expected",
        "n_model_mri_features_present_in_input",
        "n_model_mri_features_missing_from_input",
    ]

    ordered_cols = []
    for c in id_cols + risk_cols + clock_cols + qc_cols:
        if c in pred.columns and c not in ordered_cols:
            ordered_cols.append(c)

    pred_out = pred[ordered_cols].copy()

    validate_clean_output_columns(pred_out, organ, model_used_organ_features)

    instance_labels = sorted(pred_out["application_instance"].dropna().astype(str).unique().tolist())
    instance_label = instance_labels[0] if len(instance_labels) == 1 else "combined"

    out_pred = outdir / f"{pref}_apply_instance_{instance_label}_predictions.tsv"
    pred_out.to_csv(out_pred, sep="\t", index=False)

    summary = {
        "organ": organ,
        "application_instance": instance_label,
        "status": "generated",
        "n_input_rows_before_complete_case_filter": int(n_before),
        "n_output_rows": int(pred_out.shape[0]),
        "n_numeric_cols_kept": int(len(numeric_cols_kept)),
        "n_categorical_cols_kept": int(len(categorical_cols_kept)),
        "n_original_organ_features_in_bundle": int(len(organ_feature_cols_original)),
        "n_model_mri_features_expected": int(len(model_used_organ_features)),
        "n_model_mri_features_present_in_application_tsv": int(len(present_organ_features)),
        "n_model_mri_features_missing_from_application_tsv": int(len(missing_organ_features)),
        "missing_model_used_organ_features": missing_organ_features,
        "complete_case_organ_features": bool(complete_case_organ_features),
        "allow_missing_model_columns": bool(allow_missing_model_columns),
        "include_features_in_output": False,
        "risk_times_years": risk_times,
        "prediction_file": str(out_pred),
        "output_columns": ordered_cols,
    }

    out_summary = outdir / f"{pref}_apply_instance_{instance_label}_summary.json"
    with open(out_summary, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"Saved clean predictions: {out_pred}")
    print(f"Saved summary: {out_summary}")
    print("Required clock columns verified:")
    print(f"  {risk_col}")
    print(f"  {z_col}")
    print(f"  {yrs_col}")
    print(f"  {age_col}")
    print("Original MRI features were not written to final prediction file.")

    return pred_out, summary


# -----------------------------
# 6. Main
# -----------------------------

def main():
    args = parse_args()

    if args.include_features_in_output:
        warnings.warn(
            "--include-features-in-output was provided, but this revised script ignores it. "
            "Final prediction files will not include original MRI feature columns."
        )

    organ = clean_name(args.organ)
    pref = output_prefix(organ)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    risk_times = parse_risk_times(args.risk_times)

    print("============================================================")
    print("Applying pre-trained organ MRI mortality clock")
    print(f"Organ: {organ}")
    print(f"Model: {args.model_joblib}")
    print(f"Output directory: {outdir}")
    print("============================================================")

    if not os.path.exists(args.model_joblib):
        raise FileNotFoundError(f"Model joblib not found: {args.model_joblib}")

    print("Loading pretrained model bundle...")
    bundle = joblib.load(args.model_joblib)

    print("Loading application MRI data...")
    df_app = load_application_tsv(
        input_tsv=args.input_tsv,
        organ=organ,
        session_id=args.application_session_id,
    )

    print("Loading optional death/assessment-date annotations...")
    death = load_death_data(
        death_xlsx=args.death_xlsx,
        id_match_csv=args.id_match_csv,
        imaging_instance=args.imaging_instance,
        admin_censor_date=args.admin_censor_date,
    )

    if death is not None:
        df_app = df_app.merge(death, on="participant_id", how="left")

    print("Loading optional covariates...")
    cov = load_covariates(args.covariate_csv)

    if cov is not None:
        df_app = df_app.merge(cov, on="participant_id", how="left", suffixes=("", "_cov"))

    print("Constructing application covariates...")
    df_app = add_application_covariates(
        df_app,
        imaging_instance=args.imaging_instance,
        model_instance=args.model_instance,
    )

    pred, summary = apply_one_instance(
        df_app=df_app,
        bundle=bundle,
        organ=organ,
        outdir=outdir,
        risk_times=risk_times,
        complete_case_organ_features=args.complete_case_organ_features,
        allow_missing_model_columns=args.allow_missing_model_columns,
    )

    combined_pred = outdir / f"{pref}_apply_longitudinal_instances_combined_predictions.tsv"
    pred.to_csv(combined_pred, sep="\t", index=False)

    combined_summary = outdir / f"{pref}_apply_longitudinal_instances_combined_summary.json"
    with open(combined_summary, "w") as f:
        json.dump(summary, f, indent=2)

    print("============================================================")
    print("Done.")
    print(f"Combined clean prediction file: {combined_pred}")
    print(f"Combined summary file: {combined_summary}")
    print("============================================================")


if __name__ == "__main__":
    main()