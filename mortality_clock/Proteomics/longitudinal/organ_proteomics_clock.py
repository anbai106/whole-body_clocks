#!/usr/bin/env python3
"""
Apply a pre-trained Pulmonary proteomics mortality clock to longitudinal UKB proteomics instances.

This script does not refit the mortality clock. It loads the model bundle saved by the
baseline instance-0 training script, harmonizes follow-up pulmonary proteomics files
(instance 2_0 and/or 3_0) to the training-time feature/covariate schema, and outputs
risk scores plus the post-hoc clock acceleration statistics:

  pulmonary_proteomics_mortality_clock_acceleration_z
  pulmonary_proteomics_mortality_clock_acceleration_years
  pulmonary_proteomics_mortality_clock_age_years

Expected follow-up input format:
  participant_id  SFTPA1  SFTPA2  SCGB3A2  SFTPD  AGER  SCGB1A1  LAMP3  CCL18  MSR1

The saved model's preprocessor handles missing proteins by using the training-set
imputation values. Therefore, this script can apply the baseline-trained model even
when follow-up instances contain only a subset of the original baseline proteins,
as long as the required pulmonary model columns are present or can be imputed.
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


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--model-joblib",
        default=(
            "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/"
            "Pulmonary_proteomics_mortality_clock/"
            "pulmonary_proteomics_mortality_clock_model.joblib"
        ),
        help="Path to the pre-trained pulmonary_proteomics_mortality_clock_model.joblib file.",
    )
    p.add_argument(
        "--input-tsv",
        action="append",
        required=True,
        help=(
            "Follow-up pulmonary proteomics TSV. Use INSTANCE:PATH, e.g. "
            "2_0:/path/proteimics_Pulmonary_2_0.tsv. Can be repeated. "
            "If INSTANCE: is omitted, --instance must be provided."
        ),
    )
    p.add_argument(
        "--instance",
        default=None,
        help="Instance label used when --input-tsv is provided without INSTANCE:, e.g. 2_0 or 3_0.",
    )
    p.add_argument(
        "--covariate-csv",
        default="/cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv",
        help="UKB covariate CSV containing eid and instance-specific covariates.",
    )
    p.add_argument(
        "--assessment-xlsx",
        default="/cbica/home/wenju/Dataset/UKBB_UMelbourne/Death_related_var_from_Ye.xlsx",
        help="Optional UMelbourne file used to add assessment dates and death date to output.",
    )
    p.add_argument(
        "--id-match-csv",
        default="/cbica/home/wenju/Dataset/UKBB_UMelbourne/UKB_UMelbourne_vs_Penn_match_key.csv",
        help="UMelbourne-to-Penn participant ID mapping file.",
    )
    p.add_argument(
        "--outdir",
        required=True,
        help="Output directory.",
    )
    p.add_argument(
        "--organ",
        default="pulmonary",
        help="Sanitized organ name used in output column names. Default: pulmonary.",
    )
    p.add_argument(
        "--allow-missing-proteins",
        action="store_true",
        default=True,
        help="Allow expected model protein columns to be absent from follow-up files and impute them using the saved preprocessor.",
    )
    p.add_argument(
        "--strict-proteins",
        action="store_true",
        help="Fail if any protein feature expected by the pre-trained model is missing from a follow-up file.",
    )
    p.add_argument(
        "--save-model-input",
        action="store_true",
        help="Also save the harmonized model-input matrix before preprocessing.",
    )
    return p.parse_args()


def clean_name(x):
    x = re.sub(r"[^A-Za-z0-9]+", "_", str(x).strip().lower())
    x = re.sub(r"_+", "_", x).strip("_")
    if not x:
        raise ValueError("Organ name is empty after sanitization.")
    return x


def parse_input_specs(input_specs, default_instance=None):
    parsed = []
    for spec in input_specs:
        spec = str(spec)
        if ":" in spec and not spec.startswith("/"):
            instance, path = spec.split(":", 1)
        else:
            if default_instance is None:
                raise ValueError(
                    "When --input-tsv is given as a plain path, --instance must be provided. "
                    "Prefer INSTANCE:PATH, e.g. 2_0:/path/file.tsv."
                )
            instance, path = default_instance, spec
        instance = instance.strip()
        path = path.strip()
        if not instance:
            raise ValueError(f"Empty instance label in --input-tsv: {spec}")
        if not path:
            raise ValueError(f"Empty path in --input-tsv: {spec}")
        parsed.append((instance, path))
    return parsed


def instance_major(instance_label):
    """Return the UKB instance number from labels such as 2_0, 3_0, ses-M0, or 2."""
    s = str(instance_label)
    m = re.search(r"(\d+)", s)
    if not m:
        raise ValueError(f"Could not parse UKB instance number from instance label: {instance_label}")
    return m.group(1)


def load_covariates(path):
    if path is None or str(path).lower() in {"none", ""}:
        return None
    if not os.path.exists(path):
        raise FileNotFoundError(f"Covariate CSV not found: {path}")
    cov = pd.read_csv(path)
    if "eid" in cov.columns:
        cov = cov.rename(columns={"eid": "participant_id"})
    if "participant_id" not in cov.columns:
        raise ValueError("Covariate file must contain eid or participant_id.")
    cov["participant_id"] = pd.to_numeric(cov["participant_id"], errors="coerce").astype("Int64")
    cov = cov.dropna(subset=["participant_id"]).copy()
    cov["participant_id"] = cov["participant_id"].astype(int)
    return cov.drop_duplicates("participant_id", keep="first")


def load_assessment_dates(path, id_match_csv, instances):
    if path is None or str(path).lower() in {"none", ""}:
        return None
    if not os.path.exists(path):
        warnings.warn(f"Assessment XLSX not found: {path}. Continuing without assessment dates.")
        return None
    if not os.path.exists(id_match_csv):
        warnings.warn(f"ID match CSV not found: {id_match_csv}. Continuing without assessment dates.")
        return None

    d = pd.read_excel(path)
    m = pd.read_csv(id_match_csv)

    d = d.rename(columns={"eid": "participant_id_umel"})
    m = m.rename(columns={"id": "participant_id_umel", "id_upenn": "participant_id"})
    d = m.merge(d, on="participant_id_umel", how="inner")

    keep = ["participant_id"]
    for inst in sorted(set(instance_major(x) for x in instances)):
        c = f"53-{inst}.0"
        if c in d.columns:
            keep.append(c)
        else:
            warnings.warn(f"Assessment date field {c} not found in {path}.")
    if "40000-0.0" in d.columns:
        keep.append("40000-0.0")

    d = d[keep].copy()
    for c in keep:
        if c.startswith("53-") or c == "40000-0.0":
            d[c] = pd.to_datetime(d[c], errors="coerce")
    d["participant_id"] = pd.to_numeric(d["participant_id"], errors="coerce").astype("Int64")
    d = d.dropna(subset=["participant_id"]).copy()
    d["participant_id"] = d["participant_id"].astype(int)
    return d.drop_duplicates("participant_id", keep="first")


def mean_existing_numeric_columns(df, cols, out_col):
    present = [c for c in cols if c in df.columns]
    if not present:
        return df
    for c in present:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df[out_col] = df[present].mean(axis=1, skipna=True)
    return df


def add_application_covariates(df, instance_label, expected_categorical_cols):
    """
    Create training-time covariate names from instance-specific UKB follow-up covariates.

    The baseline model expects names such as age_at_baseline and bmi_at_baseline.
    For longitudinal application, these columns are interpreted as age/covariates at
    the follow-up sample instance, but are deliberately kept with the training-time
    names so the saved preprocessor can be reused without refitting.
    """
    df = df.copy()
    inst = instance_major(instance_label)

    age_col = f"age_when_attended_assessment_centre_f21003_{inst}_0"
    if age_col in df.columns:
        df["age_at_baseline"] = pd.to_numeric(df[age_col], errors="coerce")
    elif "diagnosis" in df.columns:
        warnings.warn(
            f"{age_col} not found. Using numeric diagnosis as age_at_baseline fallback. Please verify."
        )
        df["age_at_baseline"] = pd.to_numeric(df["diagnosis"], errors="coerce")
    else:
        warnings.warn(f"{age_col} not found. age_at_baseline will be missing and imputed by the saved preprocessor.")
        df["age_at_baseline"] = np.nan

    # Backward-compatible alias used by downstream plotting scripts.
    df["age_at_imaging"] = df["age_at_baseline"]

    if "sex_f31_0_0" in df.columns:
        df["sex"] = df["sex_f31_0_0"].astype(str).str.strip()
    elif "sex" in df.columns:
        df["sex"] = df["sex"].astype(str).str.strip()
    elif "Sex" in df.columns:
        df["sex"] = df["Sex"].astype(str).str.strip()
    else:
        warnings.warn("sex_f31_0_0 not found. sex will be missing and imputed by the saved preprocessor.")
        df["sex"] = np.nan

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
            "Female": "Female",
            "Male": "Male",
        }
    )

    bmi_col = f"body_mass_index_bmi_f23104_{inst}_0"
    if bmi_col in df.columns:
        df["bmi_at_baseline"] = pd.to_numeric(df[bmi_col], errors="coerce")

    df = mean_existing_numeric_columns(
        df,
        [
            f"diastolic_blood_pressure_automated_reading_f4079_{inst}_0",
            f"diastolic_blood_pressure_automated_reading_f4079_{inst}_1",
        ],
        "diastolic_bp_at_baseline",
    )
    df = mean_existing_numeric_columns(
        df,
        [
            f"systolic_blood_pressure_automated_reading_f4080_{inst}_0",
            f"systolic_blood_pressure_automated_reading_f4080_{inst}_1",
        ],
        "systolic_bp_at_baseline",
    )

    smoking_col = f"smoking_status_f20116_{inst}_0"
    if smoking_col in df.columns:
        df["smoking_status_at_baseline"] = df[smoking_col].astype("category")

    centre_col_followup = f"uk_biobank_assessment_centre_f54_{inst}_0"
    centre_col_training = "uk_biobank_assessment_centre_f54_0_0"
    if centre_col_followup in df.columns:
        df[centre_col_training] = df[centre_col_followup].astype("category")
    elif centre_col_training in df.columns:
        df[centre_col_training] = df[centre_col_training].astype("category")

    # Ensure any expected categorical column exists.
    for c in expected_categorical_cols:
        if c not in df.columns:
            df[c] = np.nan

    return df


def read_followup_tsv(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Follow-up TSV not found: {path}")
    df = pd.read_csv(path, sep="\t")
    if "participant_id" not in df.columns:
        raise ValueError(f"participant_id is missing from follow-up TSV: {path}")
    df["participant_id"] = pd.to_numeric(df["participant_id"], errors="coerce").astype("Int64")
    df = df.dropna(subset=["participant_id"]).copy()
    df["participant_id"] = df["participant_id"].astype(int)
    if df["participant_id"].duplicated().any():
        n_dup = int(df["participant_id"].duplicated().sum())
        warnings.warn(f"{path}: found {n_dup} duplicated participant_id rows. Keeping first.")
        df = df.drop_duplicates("participant_id", keep="first")
    return df


def add_missing_model_columns(df, numeric_cols, categorical_cols, organ_feature_cols, strict_proteins=False):
    df = df.copy()
    missing_organ = [c for c in organ_feature_cols if c not in df.columns]
    missing_numeric_cov = [c for c in numeric_cols if c not in df.columns and c not in missing_organ]
    missing_cat = [c for c in categorical_cols if c not in df.columns]

    if missing_organ:
        msg = (
            f"Missing {len(missing_organ)} protein feature(s) expected by the pre-trained model: "
            f"{missing_organ[:20]}"
        )
        if strict_proteins:
            raise ValueError(msg)
        warnings.warn(msg + ". They will be imputed by the saved preprocessor.")

    for c in missing_organ + missing_numeric_cov:
        df[c] = np.nan
    for c in missing_cat:
        df[c] = np.nan

    for c in numeric_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    for c in categorical_cols:
        df[c] = df[c].astype("object")

    return df, missing_organ, missing_numeric_cov, missing_cat


def predict_absolute_risk(model, X, times_years=(5.0, 10.0, 15.0)):
    out = {}
    try:
        surv_funcs = model.predict_survival_function(X)
    except Exception as e:
        warnings.warn(f"Could not compute absolute risks from saved model: {e}")
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


def resolve_cat_indicator_column(coef_name, categorical_cols):
    """Map a saved one-hot coefficient name like cat__sex_Male to (sex, Male)."""
    if not coef_name.startswith("cat__"):
        return None, None
    rest = coef_name[len("cat__") :]
    for col in sorted(categorical_cols, key=len, reverse=True):
        prefix = f"{col}_"
        if rest.startswith(prefix):
            return col, rest[len(prefix) :]
    return None, None


def categorical_match(series, target):
    """Robust comparison for one-hot residualization categories."""
    s = series.copy()
    target_str = str(target)
    as_str = s.astype(str)
    out = as_str == target_str

    # Handle integer-looking vs float-looking categories, e.g. 0 vs 0.0.
    try:
        target_float = float(target_str)
        num = pd.to_numeric(s, errors="coerce")
        out = out | np.isclose(num.astype(float), target_float, equal_nan=False)
    except Exception:
        pass
    return out.fillna(False).astype(float).values


def apply_clock_age_and_acceleration_from_bundle(pred_df, organ, bundle):
    """
    Reconstruct post-hoc clock acceleration using saved clock_transform_info.

    The original training script saved the linear residualization coefficients and
    training residual mean/SD, but not the fitted residualization object itself. This
    function applies the saved coefficients directly to new follow-up covariates.
    """
    df = pred_df.copy()
    risk_col = f"{organ}_proteomics_mortality_risk_score"
    z_col = f"{organ}_proteomics_mortality_clock_acceleration_z"
    yrs_col = f"{organ}_proteomics_mortality_clock_acceleration_years"
    age_col = f"{organ}_proteomics_mortality_clock_age_years"

    info = bundle.get("clock_transform_info", None)
    if not isinstance(info, dict):
        warnings.warn("Model bundle does not contain clock_transform_info. Acceleration columns set to missing.")
        df[z_col] = np.nan
        df[yrs_col] = np.nan
        df[age_col] = np.nan
        return df, {"status": "missing_clock_transform_info"}

    coef = info.get("risk_score_covariate_model_coef", {}) or {}
    intercept = float(info.get("risk_score_covariate_model_intercept", 0.0))
    mean_train = info.get("risk_score_residual_mean_train", 0.0)
    sd_train = info.get("risk_score_residual_sd_train", np.nan)
    beta_age = info.get("adjusted_age_coefficient_risk_score_per_year", np.nan)
    residualization_covariates = info.get("residualization_covariates", []) or []
    categorical_covariates = info.get("categorical_residualization_covariates", []) or []

    expected = np.repeat(intercept, df.shape[0]).astype(float)
    imputation_notes = []

    for name, beta in coef.items():
        beta = float(beta)
        if name.startswith("num__"):
            col = name[len("num__") :]
            if col not in df.columns:
                warnings.warn(f"Residualization numeric covariate {col} is missing; using 0 contribution.")
                continue
            x = pd.to_numeric(df[col], errors="coerce")
            if x.isna().any():
                med = x.median(skipna=True)
                if not np.isfinite(med):
                    med = 0.0
                imputation_notes.append({"covariate": col, "n_missing": int(x.isna().sum()), "imputed_with_application_median": float(med)})
                x = x.fillna(med)
            expected += x.values.astype(float) * beta
        elif name.startswith("cat__"):
            col, target = resolve_cat_indicator_column(name, categorical_covariates)
            if col is None or col not in df.columns:
                continue
            expected += categorical_match(df[col], target) * beta

    resid_raw = df[risk_col].values.astype(float) - expected
    mean_train = float(mean_train) if mean_train is not None and np.isfinite(mean_train) else 0.0
    sd_train = float(sd_train) if sd_train is not None and np.isfinite(sd_train) else np.nan
    resid_centered = resid_raw - mean_train

    if np.isfinite(sd_train) and sd_train > 0:
        df[z_col] = resid_centered / sd_train
    else:
        warnings.warn("Training residual SD is missing or <=0. Acceleration z-score set to missing.")
        df[z_col] = np.nan

    beta_age = float(beta_age) if beta_age is not None and np.isfinite(beta_age) else np.nan
    if np.isfinite(beta_age) and abs(beta_age) > 1e-8:
        df[yrs_col] = resid_centered / beta_age
        df[age_col] = pd.to_numeric(df.get("age_at_baseline", np.nan), errors="coerce") + df[yrs_col]
    else:
        warnings.warn("Adjusted age coefficient is missing or near zero. Year-scale acceleration set to missing.")
        df[yrs_col] = np.nan
        df[age_col] = np.nan

    transform_summary = {
        "status": "ok",
        "risk_col": risk_col,
        "z_col": z_col,
        "years_col": yrs_col,
        "clock_age_col": age_col,
        "residualization_covariates": residualization_covariates,
        "categorical_residualization_covariates": categorical_covariates,
        "risk_score_covariate_model_intercept": intercept,
        "risk_score_residual_mean_train": mean_train,
        "risk_score_residual_sd_train": sd_train,
        "adjusted_age_coefficient_risk_score_per_year": beta_age,
        "imputation_notes_for_missing_residualization_covariates": imputation_notes,
        "note": (
            "Applied saved linear residualization coefficients from the baseline training model bundle. "
            "For numeric residualization covariates with missing values, this script uses the application-sample median because the original residualization imputer was not saved separately."
        ),
    }
    return df, transform_summary


def apply_one_instance(instance_label, path, bundle, cov, dates, args):
    organ = clean_name(args.organ)
    pref = f"{organ}_proteomics_mortality_clock"
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    preprocessor = bundle["preprocessor"]
    model = bundle["model"]
    numeric_cols_kept = list(bundle.get("numeric_cols_kept", []))
    categorical_cols_kept = list(bundle.get("categorical_cols_kept", []))
    organ_feature_cols = list(bundle.get("organ_feature_cols", []))

    if not numeric_cols_kept or not categorical_cols_kept:
        raise ValueError("Model bundle is missing numeric_cols_kept or categorical_cols_kept.")
    if not organ_feature_cols:
        raise ValueError("Model bundle is missing organ_feature_cols.")

    print("============================================================")
    print(f"Applying {pref} to instance {instance_label}")
    print(f"Input TSV: {path}")

    df = read_followup_tsv(path)
    df["application_instance"] = instance_label
    df["application_source_file"] = Path(path).name
    n_input = df.shape[0]

    if cov is not None:
        df = df.merge(cov, on="participant_id", how="left", suffixes=("", "_cov"))

    if dates is not None:
        df = df.merge(dates, on="participant_id", how="left", suffixes=("", "_date"))
        inst = instance_major(instance_label)
        date_col = f"53-{inst}.0"
        if date_col in df.columns:
            df["sample_date"] = df[date_col]
        if "40000-0.0" in df.columns:
            df["death_date"] = df["40000-0.0"]

    df = add_application_covariates(df, instance_label, categorical_cols_kept)
    df, missing_organ, missing_numeric_cov, missing_cat = add_missing_model_columns(
        df,
        numeric_cols=numeric_cols_kept,
        categorical_cols=categorical_cols_kept,
        organ_feature_cols=organ_feature_cols,
        strict_proteins=args.strict_proteins,
    )

    model_cols = numeric_cols_kept + categorical_cols_kept
    X_raw = df[model_cols].copy()
    X = preprocessor.transform(X_raw)
    risk = np.asarray(model.predict(X)).reshape(-1)

    risk_col = f"{organ}_proteomics_mortality_risk_score"
    pred = pd.DataFrame(
        {
            "participant_id": df["participant_id"].values,
            "application_instance": instance_label,
            "application_source_file": df["application_source_file"].values,
            risk_col: risk,
        }
    )

    for c in [
        "sample_date",
        "death_date",
        "age_at_baseline",
        "age_at_imaging",
        "sex",
        "bmi_at_baseline",
        "diastolic_bp_at_baseline",
        "systolic_bp_at_baseline",
        "smoking_status_at_baseline",
        "uk_biobank_assessment_centre_f54_0_0",
    ]:
        if c in df.columns:
            pred[c] = df[c].values

    # Keep protein columns used by the saved model where available. This helps QC.
    present_organ_features = [c for c in organ_feature_cols if c in df.columns]
    for c in present_organ_features:
        pred[c] = pd.to_numeric(df[c], errors="coerce").values

    abs_risk_df = predict_absolute_risk(model, X, times_years=(5.0, 10.0, 15.0))
    pred = pd.concat([pred.reset_index(drop=True), abs_risk_df.reset_index(drop=True)], axis=1)

    pred, transform_summary = apply_clock_age_and_acceleration_from_bundle(pred, organ, bundle)

    pred["n_model_protein_features_expected"] = len(organ_feature_cols)
    pred["n_model_protein_features_present_in_input"] = len(present_organ_features)
    pred["n_model_protein_features_missing_from_input"] = len(missing_organ)

    safe_instance = re.sub(r"[^A-Za-z0-9]+", "_", str(instance_label)).strip("_")
    pred_file = outdir / f"{pref}_apply_instance_{safe_instance}_predictions.tsv"
    pred.to_csv(pred_file, sep="\t", index=False)

    if args.save_model_input:
        X_raw_out = pd.concat(
            [df[["participant_id", "application_instance"]].reset_index(drop=True), X_raw.reset_index(drop=True)],
            axis=1,
        )
        X_raw_out.to_csv(outdir / f"{pref}_apply_instance_{safe_instance}_model_input.tsv", sep="\t", index=False)

    summary = {
        "organ": organ,
        "application_instance": instance_label,
        "input_tsv": path,
        "output_predictions_tsv": str(pred_file),
        "n_input_rows": int(n_input),
        "n_output_rows": int(pred.shape[0]),
        "model_joblib": args.model_joblib,
        "n_model_numeric_cols_kept": int(len(numeric_cols_kept)),
        "n_model_categorical_cols_kept": int(len(categorical_cols_kept)),
        "n_model_protein_features_expected": int(len(organ_feature_cols)),
        "n_model_protein_features_present_in_input": int(len(present_organ_features)),
        "n_model_protein_features_missing_from_input": int(len(missing_organ)),
        "missing_model_protein_features": missing_organ,
        "missing_numeric_covariates_added_as_nan": missing_numeric_cov,
        "missing_categorical_covariates_added_as_nan": missing_cat,
        "risk_score_summary": {
            "mean": float(np.nanmean(pred[risk_col])),
            "sd": float(np.nanstd(pred[risk_col])),
            "min": float(np.nanmin(pred[risk_col])),
            "max": float(np.nanmax(pred[risk_col])),
        },
        "clock_transform_summary": transform_summary,
        "note": (
            "Predictions were generated by applying the pre-trained baseline instance-0 Pulmonary proteomics mortality clock. "
            "No model refitting was performed. Follow-up covariates were mapped to the training-time covariate names so the saved preprocessor could be reused."
        ),
    }

    summary_file = outdir / f"{pref}_apply_instance_{safe_instance}_summary.json"
    with open(summary_file, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"Output predictions: {pred_file}")
    print(f"Output summary:     {summary_file}")
    print(f"N rows: {pred.shape[0]}")
    print(f"Expected proteins: {len(organ_feature_cols)}; present: {len(present_organ_features)}; missing: {len(missing_organ)}")
    return pred, summary


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    input_specs = parse_input_specs(args.input_tsv, args.instance)
    instances = [x[0] for x in input_specs]

    print("Loading pre-trained model bundle...")
    bundle = joblib.load(args.model_joblib)
    print(f"Loaded model: {args.model_joblib}")
    print(f"Model organ in bundle: {bundle.get('organ', 'NA')}")
    print(f"Output prefix in bundle: {bundle.get('out_prefix', 'NA')}")

    cov = load_covariates(args.covariate_csv)
    dates = load_assessment_dates(args.assessment_xlsx, args.id_match_csv, instances)

    all_pred = []
    summaries = []
    for instance_label, path in input_specs:
        pred, summary = apply_one_instance(instance_label, path, bundle, cov, dates, args)
        all_pred.append(pred)
        summaries.append(summary)

    organ = clean_name(args.organ)
    pref = f"{organ}_proteomics_mortality_clock"
    combined = pd.concat(all_pred, ignore_index=True, sort=False)
    combined_file = outdir / f"{pref}_apply_longitudinal_instances_combined_predictions.tsv"
    combined.to_csv(combined_file, sep="\t", index=False)

    combined_summary = {
        "model_joblib": args.model_joblib,
        "outdir": str(outdir),
        "combined_predictions_tsv": str(combined_file),
        "instances": summaries,
    }
    combined_summary_file = outdir / f"{pref}_apply_longitudinal_instances_combined_summary.json"
    with open(combined_summary_file, "w") as f:
        json.dump(combined_summary, f, indent=2)

    print("============================================================")
    print("Done applying longitudinal pulmonary proteomics mortality clock.")
    print(f"Combined output:  {combined_file}")
    print(f"Combined summary: {combined_summary_file}")
    print("Main acceleration columns:")
    print(f"  {organ}_proteomics_mortality_clock_acceleration_z")
    print(f"  {organ}_proteomics_mortality_clock_acceleration_years")
    print(f"  {organ}_proteomics_mortality_clock_age_years")


if __name__ == "__main__":
    main()
