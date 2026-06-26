#!/usr/bin/env python3

"""
Build a brain MRI mortality clock in UK Biobank using elastic-net Cox survival modeling.

Inputs:
  1. UKB death/assessment-date Excel file from UMelbourne
  2. UKB UMelbourne-to-Penn participant ID matching file
  3. Brain MRI MUSE gray-matter volume TSV file
  4. Optional UKB covariate CSV file
  5. Administrative censor date

Main survival outcome:
  time zero = imaging assessment date, UKB field 53-2.0
  event date = date of death, UKB field 40000-0.0
  censor date = user-provided administrative censor date

Example:
  python build_brain_mri_mortality_clock.py \
    --death-xlsx /cbica/home/wenju/Dataset/UKBB_UMelbourne/Death_related_var_from_Ye.xlsx \
    --id-match-csv /cbica/home/wenju/Dataset/UKBB_UMelbourne/UKB_UMelbourne_vs_Penn_match_key.csv \
    --brain-tsv /cbica/home/wenju/Reproducibile_paper/BrainAge/data/imaging/T1_MUSE_GM.tsv \
    --covariate-csv /cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv \
    --admin-censor-date 2024-12-31 \
    --outdir /cbica/home/wenju/Reproducibile_paper/MortalityClock/brain_mri_mortality_clock
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
    from sksurv.linear_model import CoxnetSurvivalAnalysis
    from sksurv.metrics import concordance_index_censored
    from sksurv.util import Surv
except ImportError as e:
    raise ImportError(
        "This script requires scikit-survival. Install with, for example:\n"
        "  conda install -c conda-forge scikit-survival\n"
        "or check your cluster environment/module system."
    ) from e


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--death-xlsx",
        default="/cbica/home/wenju/Dataset/UKBB_UMelbourne/Death_related_var_from_Ye.xlsx",
        help="Excel file containing UKB death and assessment-date fields.",
    )
    parser.add_argument(
        "--id-match-csv",
        default="/cbica/home/wenju/Dataset/UKBB_UMelbourne/UKB_UMelbourne_vs_Penn_match_key.csv",
        help="CSV matching UMelbourne IDs to Penn participant IDs.",
    )
    parser.add_argument(
        "--brain-tsv",
        default="/cbica/home/wenju/Reproducibile_paper/BrainAge/data/imaging/T1_MUSE_GM.tsv",
        help="Brain MRI MUSE gray-matter volume TSV.",
    )
    parser.add_argument(
        "--covariate-csv",
        default="/cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv",
        help="Optional covariate CSV.",
    )
    parser.add_argument(
        "--admin-censor-date",
        required=True,
        help=(
            "Administrative censor date for death registry follow-up, e.g. 2024-12-31. "
            "Participants alive by this date are censored at this date."
        ),
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory.",
    )
    parser.add_argument(
        "--imaging-session-id",
        type=int,
        default=1,
        help=(
            "Imaging session_id to use from the brain TSV. "
            "Default is 1, which appears to correspond to the UKB imaging visit in your file."
        ),
    )
    parser.add_argument(
        "--test-size",
        type=float,
        default=0.20,
        help="Fraction of participants reserved as held-out test set.",
    )
    parser.add_argument(
        "--validation-size",
        type=float,
        default=0.20,
        help="Fraction of training participants used for inner validation.",
    )
    parser.add_argument(
        "--random-state",
        type=int,
        default=2026,
        help="Random seed.",
    )
    parser.add_argument(
        "--max-feature-missing",
        type=float,
        default=0.20,
        help="Drop imaging/covariate features with missingness above this threshold in the training set.",
    )
    parser.add_argument(
        "--l1-ratios",
        default="0.1,0.25,0.5,0.75,1.0",
        help="Comma-separated elastic-net l1_ratio values to tune.",
    )
    parser.add_argument(
        "--n-alphas",
        type=int,
        default=100,
        help="Number of alpha values for Coxnet regularization path.",
    )
    parser.add_argument(
        "--min-followup-days",
        type=int,
        default=1,
        help="Exclude participants with follow-up time shorter than this many days.",
    )

    return parser.parse_args()


def make_onehot_encoder():
    """
    scikit-learn changed sparse -> sparse_output.
    This helper keeps the script compatible with older/newer versions.
    """
    try:
        return OneHotEncoder(handle_unknown="ignore", sparse_output=False)
    except TypeError:
        return OneHotEncoder(handle_unknown="ignore", sparse=False)


def load_death_data(death_xlsx, id_match_csv):
    df_ukb_death = pd.read_excel(death_xlsx)
    df_id_match = pd.read_csv(id_match_csv)

    df_ukb_death = df_ukb_death.rename(columns={"eid": "participant_id_umel"})
    df_id_match = df_id_match.rename(
        columns={"id": "participant_id_umel", "id_upenn": "participant_id"}
    )

    df_ukb_death = df_id_match.merge(df_ukb_death, on="participant_id_umel", how="inner")

    required_cols = ["participant_id", "53-0.0", "53-2.0", "40000-0.0"]
    missing = [c for c in required_cols if c not in df_ukb_death.columns]
    if missing:
        raise ValueError(f"Death file is missing required columns: {missing}")

    df_ukb_death = df_ukb_death[required_cols].copy()

    df_ukb_death["baseline_date"] = pd.to_datetime(df_ukb_death["53-0.0"], errors="coerce")
    df_ukb_death["imaging_date"] = pd.to_datetime(df_ukb_death["53-2.0"], errors="coerce")
    df_ukb_death["death_date"] = pd.to_datetime(df_ukb_death["40000-0.0"], errors="coerce")

    return df_ukb_death


def load_brain_data(brain_tsv, imaging_session_id):
    df_brain = pd.read_csv(brain_tsv, sep="\t")

    if "participant_id" not in df_brain.columns:
        raise ValueError("Brain TSV must contain participant_id.")

    if "session_id" in df_brain.columns and imaging_session_id is not None:
        before = df_brain.shape[0]
        df_brain = df_brain.loc[df_brain["session_id"] == imaging_session_id].copy()
        after = df_brain.shape[0]
        print(f"Filtered brain file to session_id={imaging_session_id}: {before} -> {after} rows")

    # Keep one imaging row per participant.
    if "session_id" in df_brain.columns:
        df_brain = df_brain.sort_values(["participant_id", "session_id"])
    else:
        df_brain = df_brain.sort_values(["participant_id"])

    duplicated = df_brain["participant_id"].duplicated().sum()
    if duplicated > 0:
        warnings.warn(
            f"Found {duplicated} duplicated participant_id rows in brain file. "
            "Keeping the first row per participant."
        )
        df_brain = df_brain.drop_duplicates("participant_id", keep="first")

    return df_brain


def load_covariates(covariate_csv):
    if covariate_csv is None or str(covariate_csv).lower() in ["none", ""]:
        return None

    if not os.path.exists(covariate_csv):
        warnings.warn(f"Covariate file not found: {covariate_csv}. Continuing without it.")
        return None

    df_cov = pd.read_csv(covariate_csv)

    if "eid" not in df_cov.columns:
        warnings.warn("Covariate file does not contain 'eid'. Continuing without it.")
        return None

    df_cov = df_cov.rename(columns={"eid": "participant_id"})
    return df_cov


def construct_survival_dataset(df):
    """
    Create prospective mortality survival outcome from imaging date.

    Event:
      death_date exists and death_date > imaging_date and death_date <= admin_censor_date

    Censor:
      no death observed by admin_censor_date

    Exclusions:
      missing imaging date
      death before/on imaging date
      imaging date after censor date
      non-positive follow-up
    """
    df = df.copy()

    df["death_before_or_on_imaging"] = (
        df["death_date"].notna() & df["imaging_date"].notna() & (df["death_date"] <= df["imaging_date"])
    )

    n_predeath = int(df["death_before_or_on_imaging"].sum())
    if n_predeath > 0:
        warnings.warn(
            f"Excluding {n_predeath} participants with death date before/on imaging date."
        )

    df = df.loc[df["imaging_date"].notna()].copy()
    df = df.loc[~df["death_before_or_on_imaging"]].copy()
    df = df.loc[df["imaging_date"] <= df["admin_censor_date"]].copy()

    df["event"] = (
        df["death_date"].notna()
        & (df["death_date"] > df["imaging_date"])
        & (df["death_date"] <= df["admin_censor_date"])
    )

    df["end_date"] = df["admin_censor_date"]
    df.loc[df["event"], "end_date"] = df.loc[df["event"], "death_date"]

    df["time_days"] = (df["end_date"] - df["imaging_date"]).dt.days
    df["time_years"] = df["time_days"] / 365.25

    return df


def infer_feature_columns(df):
    brain_feature_cols = [
        c for c in df.columns
        if c.startswith("MUSE_Volume_") or c.startswith("WMLS_Volume_")
    ]

    if len(brain_feature_cols) == 0:
        raise ValueError(
            "No brain MRI feature columns found. Expected names starting with "
            "'MUSE_Volume_' or 'WMLS_Volume_'."
        )

    return brain_feature_cols


def add_basic_covariates(df):
    df = df.copy()

    # Prefer Age and Sex from the brain imaging file because they appear to correspond to imaging time.
    if "Age" not in df.columns:
        raise ValueError("Brain file must contain Age column for age at imaging.")

    if "Sex" not in df.columns:
        raise ValueError("Brain file must contain Sex column.")

    df["age_at_imaging"] = pd.to_numeric(df["Age"], errors="coerce")

    # Convert Sex to a simple categorical string.
    # Your brain TSV uses F/M. This also handles 0/1 if needed.
    df["sex"] = df["Sex"].astype(str).str.strip()
    df["sex"] = df["sex"].replace({"0": "Female", "1": "Male", "F": "Female", "M": "Male"})

    # DLICV is useful to adjust brain volume measures.
    if "DLICV" in df.columns:
        df["DLICV"] = pd.to_numeric(df["DLICV"], errors="coerce")

    return df


def build_design_matrix(df, brain_feature_cols):
    """
    Define covariates and feature columns.

    The model includes age, sex, DLICV, optional scanner/assessment covariates, and brain MRI features.

    We report a clock-acceleration score by residualizing the final risk score
    against all retained non-brain covariates in the training set.
    """
    df = df.copy()

    numeric_covariates = ["age_at_imaging"]

    if "DLICV" in df.columns:
        numeric_covariates.append("DLICV")

    optional_numeric_covariates = [
        "scanner_lateral_x_brain_position_f25756_2_0",
        "scanner_transverse_y_brain_position_f25757_2_0",
        "scanner_longitudinal_z_brain_position_f25758_2_0",
        "mean_rfmri_head_motion_averaged_across_space_and_time_points_f25741_2_0",
    ]

    for c in optional_numeric_covariates:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
            numeric_covariates.append(c)

    categorical_covariates = ["sex"]

    optional_categorical_covariates = [
        "uk_biobank_assessment_centre_f54_2_0",
    ]

    for c in optional_categorical_covariates:
        if c in df.columns:
            df[c] = df[c].astype("category")
            categorical_covariates.append(c)

    # Convert brain features to numeric.
    for c in brain_feature_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    numeric_cols = numeric_covariates + brain_feature_cols
    categorical_cols = categorical_covariates

    return df, numeric_cols, categorical_cols, numeric_covariates, brain_feature_cols


def safe_stratify_vector(event_series):
    """
    Return event labels for stratified splitting only if both classes have enough samples.
    Otherwise return None.
    """
    counts = event_series.value_counts()
    if len(counts) < 2:
        return None
    if counts.min() < 2:
        return None
    return event_series


def get_feature_names(preprocessor):
    names = []

    # Numeric names.
    numeric_features = preprocessor.named_transformers_["num"].feature_names_in_
    names.extend([f"num__{c}" for c in numeric_features])

    # Categorical names.
    cat_pipeline = preprocessor.named_transformers_["cat"]
    ohe = cat_pipeline.named_steps["onehot"]
    cat_features = preprocessor.transformers_[1][2]

    try:
        cat_names = ohe.get_feature_names_out(cat_features)
    except AttributeError:
        cat_names = ohe.get_feature_names(cat_features)

    names.extend([f"cat__{c}" for c in cat_names])

    return np.array(names)


def make_preprocessor(numeric_cols, categorical_cols):
    numeric_pipeline = Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
        ]
    )

    categorical_pipeline = Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="most_frequent")),
            ("onehot", make_onehot_encoder()),
        ]
    )

    preprocessor = ColumnTransformer(
        transformers=[
            ("num", numeric_pipeline, numeric_cols),
            ("cat", categorical_pipeline, categorical_cols),
        ],
        remainder="drop",
    )

    return preprocessor


def drop_high_missing_features(X_train_raw, X_other_raw_list, numeric_cols, categorical_cols, max_missing):
    """
    Drop numeric predictors with excessive missingness in training set.
    Categorical predictors are retained.
    """
    keep_numeric = []
    dropped_numeric = []

    for c in numeric_cols:
        miss = X_train_raw[c].isna().mean()
        if miss <= max_missing:
            keep_numeric.append(c)
        else:
            dropped_numeric.append((c, float(miss)))

    if dropped_numeric:
        print(f"Dropped {len(dropped_numeric)} numeric columns with missingness > {max_missing}.")
        for c, miss in dropped_numeric[:20]:
            print(f"  dropped: {c}, missing={miss:.3f}")
        if len(dropped_numeric) > 20:
            print("  ...")

    cols = keep_numeric + categorical_cols

    X_train_raw = X_train_raw[cols].copy()
    X_other_raw_list = [x[cols].copy() for x in X_other_raw_list]

    return X_train_raw, X_other_raw_list, keep_numeric, categorical_cols, dropped_numeric


def fit_and_select_coxnet(
    X_train,
    y_train,
    X_val,
    y_val,
    feature_names,
    brain_feature_prefixes=("num__MUSE_Volume_", "num__WMLS_Volume_"),
    l1_ratios=(0.1, 0.25, 0.5, 0.75, 1.0),
    n_alphas=100,
    random_state=2026,
):
    """
    Tune l1_ratio and alpha using validation C-index.

    Covariates are assigned penalty_factor=0 when the installed scikit-survival supports it.
    Imaging features are assigned penalty_factor=1.
    """
    del random_state

    penalty_factor = np.ones(len(feature_names), dtype=float)
    is_brain_feature = np.zeros(len(feature_names), dtype=bool)

    for prefix in brain_feature_prefixes:
        is_brain_feature |= np.char.startswith(feature_names.astype(str), prefix)

    penalty_factor[~is_brain_feature] = 0.0

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
                "Installed scikit-survival does not support penalty_factor. "
                "Covariates will be penalized together with imaging features."
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
        alphas = model.alphas_

        # coefs shape is usually n_features x n_alphas.
        if coefs.ndim == 1:
            coefs = coefs[:, None]

        for j, alpha in enumerate(alphas):
            risk_val = np.dot(X_val, coefs[:, j])
            cindex = concordance_index_censored(
                y_val["event"],
                y_val["time"],
                risk_val,
            )[0]

            if np.isfinite(cindex) and cindex > best["cindex"]:
                best.update(
                    {
                        "cindex": float(cindex),
                        "l1_ratio": float(l1_ratio),
                        "alpha": float(alpha),
                        "coef": coefs[:, j].copy(),
                    }
                )

        print(
            f"  best so far: C-index={best['cindex']:.4f}, "
            f"l1_ratio={best['l1_ratio']}, alpha={best['alpha']}"
        )

    if best["alpha"] is None:
        raise RuntimeError("Failed to select a Coxnet model. Check event counts and input data.")

    return best, penalty_factor


def fit_final_model(X_trainval, y_trainval, best, penalty_factor):
    """
    Refit final Coxnet model on train+validation using selected alpha and l1_ratio.
    fit_baseline_model=True enables absolute survival/risk prediction.
    """
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
    risk = model.predict(X)
    risk = np.asarray(risk).reshape(-1)
    return risk


def predict_absolute_risk(model, X, times_years):
    """
    Predict absolute event risk = 1 - survival probability at specified times.

    Returns a DataFrame with columns risk_{t}y.
    """
    out = {}

    try:
        surv_funcs = model.predict_survival_function(X)
    except Exception as e:
        warnings.warn(f"Could not compute absolute risks from baseline survival: {e}")
        for t in times_years:
            out[f"risk_{t:g}y"] = np.repeat(np.nan, X.shape[0])
        return pd.DataFrame(out)

    for t in times_years:
        vals = []
        for sf in surv_funcs:
            try:
                s = float(sf(t))
                vals.append(1.0 - s)
            except Exception:
                vals.append(np.nan)
        out[f"risk_{t:g}y"] = vals

    return pd.DataFrame(out)


def add_clock_age_and_acceleration(pred_df, covariate_cols=None, train_mask_col="split"):
    """
    Convert risk score to a covariate-adjusted mortality-clock age and acceleration.

    This is a post-hoc interpretability transform:
      1. Fit risk_score ~ retained non-brain covariates on the training set.
         By default this includes age_at_imaging, sex, DLICV, scanner/head-motion
         variables, and assessment center if those variables were retained.
      2. Compute residual risk beyond the expected risk from those covariates.
      3. Standardize the residual using the training-set residual SD.
      4. Convert residual risk to approximate years using the adjusted age coefficient.

    Primary model output remains the survival risk score and predicted absolute risk.
    The clock acceleration variables are covariate-adjusted residual phenotypes.
    """
    df = pred_df.copy()

    if covariate_cols is None:
        covariate_cols = ["age_at_imaging", "sex"]

    # Keep only covariates that are present. Do not include brain MRI features here.
    covariate_cols = [c for c in covariate_cols if c in df.columns]

    if "age_at_imaging" not in covariate_cols and "age_at_imaging" in df.columns:
        covariate_cols = ["age_at_imaging"] + covariate_cols

    if len(covariate_cols) == 0:
        raise ValueError("No covariates available for clock-acceleration residualization.")

    train_df = df.loc[df[train_mask_col] == "train"].copy()

    if train_df.shape[0] < 10:
        warnings.warn("Too few training samples for clock-age transform.")
        df["brain_mri_mortality_clock_acceleration_years"] = np.nan
        df["brain_mri_mortality_clock_age_years"] = np.nan
        df["brain_mri_mortality_clock_acceleration_z"] = np.nan
        return df, None

    # Identify numeric versus categorical residualization covariates.
    # We intentionally do not standardize numeric covariates here because the
    # age coefficient is used to convert residual risk-score units into years.
    numeric_covs = []
    categorical_covs = []

    for c in covariate_cols:
        if c == "sex":
            categorical_covs.append(c)
        elif pd.api.types.is_numeric_dtype(df[c]):
            numeric_covs.append(c)
        else:
            categorical_covs.append(c)

    # Make a residualization preprocessor that can handle missing covariates.
    transformers = []

    if len(numeric_covs) > 0:
        num_pipe = Pipeline(
            steps=[
                ("imputer", SimpleImputer(strategy="median")),
            ]
        )
        transformers.append(("num", num_pipe, numeric_covs))

    if len(categorical_covs) > 0:
        cat_pipe = Pipeline(
            steps=[
                ("imputer", SimpleImputer(strategy="most_frequent")),
                ("onehot", make_onehot_encoder()),
            ]
        )
        transformers.append(("cat", cat_pipe, categorical_covs))

    residualization_preprocessor = ColumnTransformer(
        transformers=transformers,
        remainder="drop",
    )

    X_train_raw = train_df[covariate_cols].copy()
    X_all_raw = df[covariate_cols].copy()

    # Ensure categorical columns are consistently string/categorical-like.
    for c in categorical_covs:
        X_train_raw[c] = X_train_raw[c].astype("object")
        X_all_raw[c] = X_all_raw[c].astype("object")

    # Ensure numeric columns are numeric.
    for c in numeric_covs:
        X_train_raw[c] = pd.to_numeric(X_train_raw[c], errors="coerce")
        X_all_raw[c] = pd.to_numeric(X_all_raw[c], errors="coerce")

    X_train = residualization_preprocessor.fit_transform(X_train_raw)
    X_all = residualization_preprocessor.transform(X_all_raw)

    lr = LinearRegression()
    lr.fit(X_train, train_df["brain_mri_mortality_risk_score"].values)

    expected = lr.predict(X_all)
    residual_raw = df["brain_mri_mortality_risk_score"].values - expected

    train_index = df[train_mask_col].values == "train"
    train_resid_raw = residual_raw[train_index]
    train_resid_mean = float(np.nanmean(train_resid_raw))
    train_resid_sd = float(np.nanstd(train_resid_raw))

    # Center using the training-set residual mean. With an intercept this should
    # be approximately zero, but explicit centering improves reproducibility.
    residual = residual_raw - train_resid_mean

    df["brain_mri_mortality_clock_acceleration_z"] = (
        residual / train_resid_sd if train_resid_sd > 0 else np.nan
    )

    # Recover the adjusted age coefficient in raw risk-score units per year.
    # This only works because numeric covariates were not standardized.
    beta_age = np.nan
    residualization_feature_names = []

    if len(numeric_covs) > 0:
        residualization_feature_names.extend([f"num__{c}" for c in numeric_covs])

    if len(categorical_covs) > 0:
        cat_pipeline = residualization_preprocessor.named_transformers_["cat"]
        ohe = cat_pipeline.named_steps["onehot"]
        try:
            cat_names = ohe.get_feature_names_out(categorical_covs)
        except AttributeError:
            cat_names = ohe.get_feature_names(categorical_covs)
        residualization_feature_names.extend([f"cat__{c}" for c in cat_names])

    age_feature_name = "num__age_at_imaging"
    if age_feature_name in residualization_feature_names:
        age_idx = residualization_feature_names.index(age_feature_name)
        beta_age = float(lr.coef_[age_idx])

    if np.isfinite(beta_age) and abs(beta_age) > 1e-8:
        df["brain_mri_mortality_clock_acceleration_years"] = residual / beta_age
        df["brain_mri_mortality_clock_age_years"] = (
            df["age_at_imaging"] + df["brain_mri_mortality_clock_acceleration_years"]
        )
    else:
        warnings.warn(
            "Adjusted age coefficient in risk_score ~ covariates is near zero or unavailable. "
            "Year-scale clock acceleration will be set to missing."
        )
        df["brain_mri_mortality_clock_acceleration_years"] = np.nan
        df["brain_mri_mortality_clock_age_years"] = np.nan

    transform_info = {
        "residualization_covariates": covariate_cols,
        "numeric_residualization_covariates": numeric_covs,
        "categorical_residualization_covariates": categorical_covs,
        "risk_score_covariate_model_intercept": float(lr.intercept_),
        "risk_score_covariate_model_coef": {
            col: float(coef) for col, coef in zip(residualization_feature_names, lr.coef_)
        },
        "adjusted_age_coefficient_risk_score_per_year": (
            float(beta_age) if np.isfinite(beta_age) else None
        ),
        "risk_score_residual_mean_train": train_resid_mean,
        "risk_score_residual_sd_train": train_resid_sd,
        "note": (
            "Clock acceleration is the residual of the Cox risk score after adjustment "
            "for all retained non-brain covariates listed in residualization_covariates. "
            "The z-score is recommended for downstream analyses. The year-scale variable "
            "is an approximate transform using the adjusted age coefficient."
        ),
    }

    return df, transform_info


def main():
    args = parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    admin_censor_date = pd.to_datetime(args.admin_censor_date)
    l1_ratios = tuple(float(x) for x in args.l1_ratios.split(","))

    print("Loading death/assessment data...")
    df_death = load_death_data(args.death_xlsx, args.id_match_csv)
    df_death["admin_censor_date"] = admin_censor_date

    print("Loading brain MRI data...")
    df_brain = load_brain_data(args.brain_tsv, args.imaging_session_id)

    print("Loading optional covariates...")
    df_cov = load_covariates(args.covariate_csv)

    print("Merging data...")
    df = df_brain.merge(df_death, on="participant_id", how="inner")

    if df_cov is not None:
        df = df.merge(df_cov, on="participant_id", how="left", suffixes=("", "_cov"))

    print("Constructing survival outcome...")
    df = construct_survival_dataset(df)
    df = df.loc[df["time_days"] >= args.min_followup_days].copy()

    print("Adding basic covariates...")
    df = add_basic_covariates(df)

    brain_feature_cols = infer_feature_columns(df)
    print(f"Found {len(brain_feature_cols)} brain MRI features.")

    df, numeric_cols, categorical_cols, numeric_covariates, brain_feature_cols = build_design_matrix(
        df, brain_feature_cols
    )

    required_model_cols = ["participant_id", "time_years", "event", "age_at_imaging", "sex"]
    df = df.dropna(subset=required_model_cols).copy()

    print("Final prospective imaging mortality dataset:")
    print(f"  N = {df.shape[0]}")
    print(f"  Deaths after imaging = {int(df['event'].sum())}")
    print(f"  Censored = {int((~df['event']).sum())}")
    print(f"  Median follow-up years = {df['time_years'].median():.2f}")

    if df["event"].sum() < 20:
        warnings.warn(
            "Very few mortality events. The fitted clock may be unstable. "
            "Consider longer follow-up, broader cohort, or lower-dimensional features."
        )

    # Save analysis dataset before splitting.
    keep_dataset_cols = [
        "participant_id",
        "baseline_date",
        "imaging_date",
        "death_date",
        "admin_censor_date",
        "end_date",
        "time_days",
        "time_years",
        "event",
        "age_at_imaging",
        "sex",
    ]
    df[keep_dataset_cols].to_csv(
        outdir / "brain_mri_mortality_survival_dataset.tsv",
        sep="\t",
        index=False,
    )

    print("Splitting into train/validation/test...")
    stratify = safe_stratify_vector(df["event"])

    df_trainval, df_test = train_test_split(
        df,
        test_size=args.test_size,
        random_state=args.random_state,
        stratify=stratify,
    )

    stratify_trainval = safe_stratify_vector(df_trainval["event"])

    df_train, df_val = train_test_split(
        df_trainval,
        test_size=args.validation_size,
        random_state=args.random_state,
        stratify=stratify_trainval,
    )

    print(f"  Train N={df_train.shape[0]}, events={int(df_train['event'].sum())}")
    print(f"  Val   N={df_val.shape[0]}, events={int(df_val['event'].sum())}")
    print(f"  Test  N={df_test.shape[0]}, events={int(df_test['event'].sum())}")

    X_train_raw = df_train[numeric_cols + categorical_cols].copy()
    X_val_raw = df_val[numeric_cols + categorical_cols].copy()
    X_test_raw = df_test[numeric_cols + categorical_cols].copy()
    X_trainval_raw = df_trainval[numeric_cols + categorical_cols].copy()

    X_train_raw, other_list, numeric_cols_kept, categorical_cols_kept, dropped_numeric = (
        drop_high_missing_features(
            X_train_raw,
            [X_val_raw, X_test_raw, X_trainval_raw],
            numeric_cols,
            categorical_cols,
            args.max_feature_missing,
        )
    )
    X_val_raw, X_test_raw, X_trainval_raw = other_list

    # Covariates used for post-hoc clock-acceleration residualization.
    # These are the retained non-brain predictors. Brain MRI features are NOT
    # included here because the clock phenotype should preserve brain MRI signal.
    residualization_covariates = [
        c for c in (numeric_cols_kept + categorical_cols_kept)
        if c not in brain_feature_cols
    ]
    print("Residualizing clock acceleration on retained non-brain covariates:")
    for c in residualization_covariates:
        print(f"  {c}")

    preprocessor = make_preprocessor(numeric_cols_kept, categorical_cols_kept)

    print("Fitting preprocessing pipeline...")
    X_train = preprocessor.fit_transform(X_train_raw)
    X_val = preprocessor.transform(X_val_raw)
    X_test = preprocessor.transform(X_test_raw)
    X_trainval = preprocessor.transform(X_trainval_raw)

    feature_names = get_feature_names(preprocessor)

    y_train = Surv.from_arrays(
        event=df_train["event"].astype(bool).values,
        time=df_train["time_years"].astype(float).values,
    )
    y_val = Surv.from_arrays(
        event=df_val["event"].astype(bool).values,
        time=df_val["time_years"].astype(float).values,
    )
    y_trainval = Surv.from_arrays(
        event=df_trainval["event"].astype(bool).values,
        time=df_trainval["time_years"].astype(float).values,
    )
    y_test = Surv.from_arrays(
        event=df_test["event"].astype(bool).values,
        time=df_test["time_years"].astype(float).values,
    )

    print("Tuning elastic-net Cox model...")
    best, penalty_factor = fit_and_select_coxnet(
        X_train=X_train,
        y_train=y_train,
        X_val=X_val,
        y_val=y_val,
        feature_names=feature_names,
        l1_ratios=l1_ratios,
        n_alphas=args.n_alphas,
        random_state=args.random_state,
    )

    print("Best validation model:")
    print(json.dumps(best | {"coef": "omitted"}, indent=2))

    print("Refitting final model on train+validation...")
    final_model = fit_final_model(X_trainval, y_trainval, best, penalty_factor)

    print("Generating predictions...")
    risk_train = predict_risk_score(final_model, X_train)
    risk_val = predict_risk_score(final_model, X_val)
    risk_test = predict_risk_score(final_model, X_test)
    risk_trainval = predict_risk_score(final_model, X_trainval)

    cindex_train = concordance_index_censored(
        df_train["event"].astype(bool).values,
        df_train["time_years"].astype(float).values,
        risk_train,
    )[0]
    cindex_val = concordance_index_censored(
        df_val["event"].astype(bool).values,
        df_val["time_years"].astype(float).values,
        risk_val,
    )[0]
    cindex_test = concordance_index_censored(
        df_test["event"].astype(bool).values,
        df_test["time_years"].astype(float).values,
        risk_test,
    )[0]
    cindex_trainval = concordance_index_censored(
        df_trainval["event"].astype(bool).values,
        df_trainval["time_years"].astype(float).values,
        risk_trainval,
    )[0]

    print(f"Train C-index:     {cindex_train:.4f}")
    print(f"Validation C-index:{cindex_val:.4f}")
    print(f"Train+Val C-index: {cindex_trainval:.4f}")
    print(f"Test C-index:      {cindex_test:.4f}")

    def make_pred_frame(df_part, split_name, risk):
        base_cols = [
            "participant_id",
            "imaging_date",
            "death_date",
            "admin_censor_date",
            "end_date",
            "time_years",
            "event",
            "age_at_imaging",
            "sex",
        ]

        extra_cov_cols = [
            c for c in residualization_covariates
            if c in df_part.columns and c not in base_cols
        ]

        out = df_part[base_cols + extra_cov_cols].copy()
        out["split"] = split_name
        out["brain_mri_mortality_risk_score"] = risk
        return out

    pred_train = make_pred_frame(df_train, "train", risk_train)
    pred_val = make_pred_frame(df_val, "validation", risk_val)
    pred_test = make_pred_frame(df_test, "test", risk_test)

    # Absolute risks at fixed horizons.
    risk_times = [5.0, 10.0, 15.0]

    pred_train = pd.concat(
        [pred_train.reset_index(drop=True), predict_absolute_risk(final_model, X_train, risk_times)],
        axis=1,
    )
    pred_val = pd.concat(
        [pred_val.reset_index(drop=True), predict_absolute_risk(final_model, X_val, risk_times)],
        axis=1,
    )
    pred_test = pd.concat(
        [pred_test.reset_index(drop=True), predict_absolute_risk(final_model, X_test, risk_times)],
        axis=1,
    )

    pred_all = pd.concat([pred_train, pred_val, pred_test], axis=0, ignore_index=True)

    print("Adding approximate mortality-clock age and acceleration...")
    pred_all, clock_transform_info = add_clock_age_and_acceleration(
        pred_all,
        covariate_cols=residualization_covariates,
    )

    pred_all.to_csv(
        outdir / "brain_mri_mortality_clock_predictions.tsv",
        sep="\t",
        index=False,
    )
    pred_test.to_csv(
        outdir / "brain_mri_mortality_clock_test_predictions.tsv",
        sep="\t",
        index=False,
    )

    print("Saving coefficients...")
    coef = np.asarray(final_model.coef_).reshape(-1)

    coef_df = pd.DataFrame(
        {
            "feature": feature_names,
            "coefficient": coef,
            "abs_coefficient": np.abs(coef),
            "penalty_factor": penalty_factor,
            "is_nonzero": coef != 0,
            "is_brain_mri_feature": [
                f.startswith("num__MUSE_Volume_") or f.startswith("num__WMLS_Volume_")
                for f in feature_names
            ],
        }
    ).sort_values("abs_coefficient", ascending=False)

    coef_df.to_csv(
        outdir / "brain_mri_mortality_clock_coefficients.tsv",
        sep="\t",
        index=False,
    )

    nonzero_coef_df = coef_df.loc[coef_df["is_nonzero"]].copy()
    nonzero_coef_df.to_csv(
        outdir / "brain_mri_mortality_clock_nonzero_coefficients.tsv",
        sep="\t",
        index=False,
    )

    print("Saving model object...")
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
        "admin_censor_date": str(admin_censor_date.date()),
    }

    joblib.dump(model_bundle, outdir / "brain_mri_mortality_clock_model.joblib")

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
        "n_original_brain_features": int(len(brain_feature_cols)),
        "n_numeric_cols_kept": int(len(numeric_cols_kept)),
        "n_categorical_cols_kept": int(len(categorical_cols_kept)),
        "n_nonzero_coefficients": int(nonzero_coef_df.shape[0]),
        "n_residualization_covariates": int(len(residualization_covariates)),
        "residualization_covariates": residualization_covariates,
        "admin_censor_date": str(admin_censor_date.date()),
        "time_zero": "UKB imaging assessment date, field 53-2.0",
        "event_date": "UKB death date, field 40000-0.0",
        "note": (
            "Primary score is brain_mri_mortality_risk_score from elastic-net Cox. "
            "Clock age/acceleration are post-hoc residualized transforms adjusted for "
            "all retained non-brain covariates."
        ),
    }

    with open(outdir / "brain_mri_mortality_clock_performance.json", "w") as f:
        json.dump(performance, f, indent=2)

    print("Done.")
    print(f"Outputs written to: {outdir}")
    print("Main output files:")
    print(f"  {outdir / 'brain_mri_mortality_clock_predictions.tsv'}")
    print(f"  {outdir / 'brain_mri_mortality_clock_test_predictions.tsv'}")
    print(f"  {outdir / 'brain_mri_mortality_clock_coefficients.tsv'}")
    print(f"  {outdir / 'brain_mri_mortality_clock_nonzero_coefficients.tsv'}")
    print(f"  {outdir / 'brain_mri_mortality_clock_model.joblib'}")
    print(f"  {outdir / 'brain_mri_mortality_clock_performance.json'}")


if __name__ == "__main__":
    main()