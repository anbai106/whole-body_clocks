#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
STEP III: Link NHANES mortality EPOCH acceleration to baseline disease status

Purpose
-------
Use the non-disease-input NHANES mortality EPOCH score from STEP II and test
whether mortality_epoch_acceleration_z is associated with baseline disease
status using disease modules that were excluded from model training.

Primary design
--------------
Logistic regression:
    disease status ~ mortality_epoch_acceleration_z + covariates

Primary comparison:
    disease cases vs strict healthy controls

Sensitivity comparison:
    disease cases vs disease-negative controls

Survival analysis
-----------------
The script supports optional prospective disease-onset survival analysis if
the user provides a long-format onset table:

    SEQN    disease    time    event

However, standard public NHANES linked mortality data do not provide prospective
nonfatal disease-onset dates. Self-reported age at diagnosis is not used for
prospective survival here because disease onset often occurred before baseline.

Main outputs
------------
disease_status_wide.tsv
healthy_control_definition.tsv
logistic_results_healthy_control.tsv
logistic_results_disease_negative.tsv
disease_association_summary_all.tsv
survival_results_optional.tsv
skipped_endpoints.tsv
step3_config.json
"""

import argparse
import json
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

import statsmodels.api as sm
import statsmodels.formula.api as smf


# ============================================================
# NHANES cycle configuration
# ============================================================

CYCLES = [
    {"cycle": "1999-2000", "begin_year": 1999},
    {"cycle": "2001-2002", "begin_year": 2001},
    {"cycle": "2003-2004", "begin_year": 2003},
    {"cycle": "2005-2006", "begin_year": 2005},
    {"cycle": "2007-2008", "begin_year": 2007},
    {"cycle": "2009-2010", "begin_year": 2009},
    {"cycle": "2011-2012", "begin_year": 2011},
    {"cycle": "2013-2014", "begin_year": 2013},
    {"cycle": "2015-2016", "begin_year": 2015},
    {"cycle": "2017-2018", "begin_year": 2017},
]

DISEASE_MODULE_PATTERNS = [
    r"^MCQ",   # medical conditions
    r"^DIQ",   # diabetes
    r"^BPQ",   # blood pressure / cholesterol diagnosis
    r"^CDQ",   # cardiovascular questionnaire, if present
    r"^CKQ",   # kidney conditions, if present
    r"^KIQ",   # kidney/urologic conditions
    r"^HEQ",   # hepatitis / liver-related questionnaire
]

CORE_COVARS = [
    "RIDAGEYR",
    "RIAGENDR",
    "RIDRETH1",
    "DMDEDUC2",
    "INDFMPIR",
    "cycle",
]

DEFAULT_SCORE_COL = "mortality_epoch_acceleration_z"


# ============================================================
# Disease endpoint definitions
# ============================================================
# NHANES disease variables generally use:
#   1 = Yes
#   2 = No
#   7/9 or 77/99 etc. = missing/refused/don't know
#
# Composite endpoints use OR across available variables.

DISEASE_DEFINITIONS = {
    "diabetes": {
        "description": "Self-reported doctor-diagnosed diabetes",
        "vars_any_yes": ["DIQ010"],
        "category": "metabolic",
    },
    "hypertension": {
        "description": "Self-reported doctor-diagnosed hypertension",
        "vars_any_yes": ["BPQ020"],
        "category": "cardiovascular",
    },
    "high_cholesterol": {
        "description": "Self-reported doctor-diagnosed high cholesterol",
        "vars_any_yes": ["BPQ080"],
        "category": "cardiometabolic",
    },
    "asthma_ever": {
        "description": "Ever told had asthma",
        "vars_any_yes": ["MCQ010"],
        "category": "pulmonary",
    },
    "asthma_current": {
        "description": "Still has asthma among participants ever diagnosed",
        "vars_any_yes": ["MCQ025"],
        "category": "pulmonary",
    },
    "arthritis": {
        "description": "Self-reported arthritis",
        "vars_any_yes": ["MCQ160A"],
        "category": "musculoskeletal",
    },
    "congestive_heart_failure": {
        "description": "Self-reported congestive heart failure",
        "vars_any_yes": ["MCQ160B"],
        "category": "cardiovascular",
    },
    "coronary_heart_disease": {
        "description": "Self-reported coronary heart disease",
        "vars_any_yes": ["MCQ160C"],
        "category": "cardiovascular",
    },
    "angina": {
        "description": "Self-reported angina",
        "vars_any_yes": ["MCQ160D"],
        "category": "cardiovascular",
    },
    "myocardial_infarction": {
        "description": "Self-reported heart attack / myocardial infarction",
        "vars_any_yes": ["MCQ160E"],
        "category": "cardiovascular",
    },
    "stroke": {
        "description": "Self-reported stroke",
        "vars_any_yes": ["MCQ160F"],
        "category": "brain_vascular",
    },
    "emphysema": {
        "description": "Self-reported emphysema",
        "vars_any_yes": ["MCQ160G"],
        "category": "pulmonary",
    },
    "chronic_bronchitis": {
        "description": "Self-reported chronic bronchitis",
        "vars_any_yes": ["MCQ160K"],
        "category": "pulmonary",
    },
    "copd_composite": {
        "description": "Self-reported emphysema or chronic bronchitis",
        "vars_any_yes": ["MCQ160G", "MCQ160K"],
        "category": "pulmonary",
    },
    "liver_condition": {
        "description": "Self-reported liver condition, if MCQ160L is available",
        "vars_any_yes": ["MCQ160L"],
        "category": "hepatic",
    },
    "thyroid_condition": {
        "description": "Self-reported thyroid condition, if MCQ160M is available",
        "vars_any_yes": ["MCQ160M"],
        "category": "endocrine",
    },
    "cancer_any": {
        "description": "Ever told had cancer or malignancy",
        "vars_any_yes": ["MCQ220"],
        "category": "cancer",
    },
    "kidney_disease": {
        "description": "Self-reported weak/failing kidneys, if KIQ022 is available",
        "vars_any_yes": ["KIQ022"],
        "category": "renal",
    },
    "cardiovascular_disease_composite": {
        "description": "Any self-reported CHF, CHD, angina, MI, or stroke",
        "vars_any_yes": ["MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E", "MCQ160F"],
        "category": "cardiovascular",
    },
    "heart_disease_composite": {
        "description": "Any self-reported CHF, CHD, angina, or MI",
        "vars_any_yes": ["MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E"],
        "category": "cardiovascular",
    },
    "major_disease_composite": {
        "description": "Any major self-reported disease endpoint included in STEP III",
        "vars_any_yes": [
            "DIQ010",
            "BPQ020",
            "MCQ010",
            "MCQ160A",
            "MCQ160B",
            "MCQ160C",
            "MCQ160D",
            "MCQ160E",
            "MCQ160F",
            "MCQ160G",
            "MCQ160K",
            "MCQ220",
            "KIQ022",
        ],
        "category": "global",
    },
}


# ============================================================
# Utility functions
# ============================================================

def module_name(path: Path) -> str:
    return path.stem.upper()


def matches_any_pattern(text: str, patterns) -> bool:
    text = str(text).upper()
    return any(re.search(pattern, text, flags=re.IGNORECASE) for pattern in patterns)


def decode_object_columns(df: pd.DataFrame) -> pd.DataFrame:
    object_cols = df.select_dtypes(include=["object"]).columns
    for col in object_cols:
        df[col] = df[col].map(
            lambda x: x.decode("utf-8", errors="ignore") if isinstance(x, bytes) else x
        )
    return df


def read_xpt_safe(path: Path):
    try:
        try:
            import pyreadstat
            df, meta = pyreadstat.read_xport(str(path))
            labels = dict(zip(meta.column_names, meta.column_labels))
        except Exception:
            df = pd.read_sas(path, format="xport", encoding="latin1")
            labels = {c: "" for c in df.columns}
    except Exception as exc:
        warnings.warn(f"Could not read XPT file: {path}\nError: {exc}")
        return None, {}

    df = decode_object_columns(df)
    df.columns = [str(c).upper() for c in df.columns]
    labels = {str(k).upper(): v for k, v in labels.items()}
    return df, labels


def collapse_by_seqn(df: pd.DataFrame):
    if df is None or "SEQN" not in df.columns:
        return None
    if not df["SEQN"].duplicated().any():
        return df
    return df.groupby("SEQN", as_index=False).first()


def safe_merge(left, right):
    if right is None:
        return left
    if left is None:
        return right

    new_cols = [c for c in right.columns if c == "SEQN" or c not in left.columns]
    if len(new_cols) <= 1:
        return left

    return left.merge(right[new_cols], on="SEQN", how="outer")


def normalize_dtypes(df: pd.DataFrame) -> pd.DataFrame:
    for col in df.columns:
        if col in ["cycle"]:
            continue
        if pd.api.types.is_numeric_dtype(df[col]):
            continue

        nonmissing = df[col].notna().sum()
        if nonmissing == 0:
            continue

        numeric = pd.to_numeric(df[col], errors="coerce")
        numeric_nonmissing = numeric.notna().sum()

        if numeric_nonmissing / max(nonmissing, 1) >= 0.90:
            df[col] = numeric

    return df


def clean_special_missing_codes_series(s: pd.Series) -> pd.Series:
    s = pd.to_numeric(s, errors="coerce")
    missing_codes = {7, 9, 77, 99, 777, 999, 7777, 9999, 77777, 99999}
    return s.mask(s.isin(missing_codes), np.nan)


def yes_no_to_binary(s: pd.Series) -> pd.Series:
    """
    Convert NHANES yes/no variables to 1/0/NA.
    1 = yes, 2 = no; special missing codes -> NA.
    """
    x = clean_special_missing_codes_series(s)
    out = pd.Series(np.nan, index=s.index, dtype=float)
    out.loc[x == 1] = 1.0
    out.loc[x == 2] = 0.0
    return out


def any_yes_binary(df: pd.DataFrame, columns):
    """
    Composite disease endpoint.
    Returns:
      1 if any available component == yes
      0 if at least one available component is valid and all valid components == no
      NA if no valid component is available
    """
    available = [c for c in columns if c in df.columns]
    if not available:
        return pd.Series(np.nan, index=df.index, dtype=float), []

    mats = []
    for c in available:
        mats.append(yes_no_to_binary(df[c]).rename(c))

    M = pd.concat(mats, axis=1)

    any_yes = M.eq(1).any(axis=1)
    any_valid = M.notna().any(axis=1)
    all_valid_no = M.fillna(0).eq(0).all(axis=1) & any_valid

    out = pd.Series(np.nan, index=df.index, dtype=float)
    out.loc[any_yes] = 1.0
    out.loc[all_valid_no & ~any_yes] = 0.0

    return out, available


def p_adjust_bh(pvals):
    pvals = np.asarray(pvals, dtype=float)
    qvals = np.full_like(pvals, np.nan, dtype=float)

    valid = np.isfinite(pvals)
    p = pvals[valid]
    n = len(p)

    if n == 0:
        return qvals

    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(1, n + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.minimum(q, 1.0)

    tmp = np.empty(n, dtype=float)
    tmp[order] = q
    qvals[valid] = tmp
    return qvals


# ============================================================
# Load disease modules
# ============================================================

def load_disease_modules(nhanes_root: Path):
    all_cycle = []
    source_rows = []
    label_rows = []

    for ci in CYCLES:
        cycle = ci["cycle"]
        begin_year = ci["begin_year"]
        cycle_dir = nhanes_root / cycle

        if not cycle_dir.exists():
            warnings.warn(f"Missing cycle directory: {cycle_dir}")
            continue

        print(f"Reading disease modules for cycle: {cycle}")

        xpt_files = sorted([p for p in cycle_dir.rglob("*") if p.suffix.lower() == ".xpt"])

        merged = None

        for path in xpt_files:
            if path.parent.name != "Questionnaire":
                continue

            mod = module_name(path)

            if not matches_any_pattern(mod, DISEASE_MODULE_PATTERNS):
                continue

            df, labels = read_xpt_safe(path)
            if df is None or "SEQN" not in df.columns:
                continue

            df = normalize_dtypes(df)
            df = collapse_by_seqn(df)
            if df is None:
                continue

            source_rows.append({
                "cycle": cycle,
                "component": path.parent.name,
                "file": path.name,
                "module": mod,
                "n_columns": len(df.columns),
                "columns": ";".join(df.columns),
            })

            for col in df.columns:
                if col == "SEQN":
                    continue
                label_rows.append({
                    "cycle": cycle,
                    "file": path.name,
                    "module": mod,
                    "variable": col,
                    "label": labels.get(col, ""),
                })

            merged = safe_merge(merged, df)

        if merged is None:
            warnings.warn(f"No disease modules loaded for cycle {cycle}")
            continue

        merged["cycle"] = cycle
        merged["cycle_begin_year"] = begin_year

        all_cycle.append(merged)

    if not all_cycle:
        raise RuntimeError("No disease questionnaire modules were loaded.")

    disease_raw = pd.concat(all_cycle, axis=0, ignore_index=True, sort=False)
    disease_raw = normalize_dtypes(disease_raw)

    sources = pd.DataFrame(source_rows)
    labels = pd.DataFrame(label_rows)

    return disease_raw, sources, labels


# ============================================================
# Build harmonized disease status
# ============================================================

def build_disease_status(disease_raw: pd.DataFrame):
    out = disease_raw[["SEQN", "cycle", "cycle_begin_year"]].copy()
    mapping_rows = []

    for disease, spec in DISEASE_DEFINITIONS.items():
        status, used_vars = any_yes_binary(disease_raw, spec["vars_any_yes"])
        out[disease] = status

        mapping_rows.append({
            "disease": disease,
            "description": spec["description"],
            "category": spec["category"],
            "requested_variables": ";".join(spec["vars_any_yes"]),
            "available_variables_used": ";".join(used_vars),
            "n_available_variables": len(used_vars),
            "n_nonmissing": int(status.notna().sum()),
            "n_cases": int((status == 1).sum(skipna=True)),
            "n_controls_disease_negative": int((status == 0).sum(skipna=True)),
        })

    mapping = pd.DataFrame(mapping_rows)

    # Strict healthy control:
    # no major disease among the endpoints that were actually measurable.
    healthy_components = [
        "diabetes",
        "hypertension",
        "asthma_ever",
        "arthritis",
        "cardiovascular_disease_composite",
        "copd_composite",
        "cancer_any",
        "kidney_disease",
        "liver_condition",
    ]

    healthy_components = [c for c in healthy_components if c in out.columns]

    M = out[healthy_components]
    any_disease = M.eq(1).any(axis=1)
    any_valid = M.notna().any(axis=1)

    # Strict healthy requires no positive disease and at least one measured endpoint.
    # Missing endpoints do not automatically count as healthy, but they also do not
    # automatically exclude the participant if other major endpoints are measured
    # and all measured endpoints are negative.
    out["healthy_control_strict"] = np.where(
        (~any_disease) & any_valid,
        1.0,
        0.0,
    )

    healthy_def = pd.DataFrame({
        "healthy_control_component": healthy_components,
        "definition": "Participant must have no positive value across measured major disease endpoints.",
    })

    return out, mapping, healthy_def


# ============================================================
# Logistic disease association
# ============================================================

def prepare_covariates(df: pd.DataFrame):
    df = df.copy()

    for c in ["RIAGENDR", "RIDRETH1", "DMDEDUC2", "cycle"]:
        if c in df.columns:
            df[c] = df[c].astype("object")
            df[c] = df[c].where(df[c].notna(), "Missing")
            df[c] = df[c].astype(str)

    for c in ["RIDAGEYR", "INDFMPIR"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    return df


def fit_logistic_one(
    df: pd.DataFrame,
    disease: str,
    score_col: str,
    control_mode: str,
    min_case: int,
    min_control: int,
    use_weights: bool = False,
):
    """
    control_mode:
      healthy_control:
        cases = disease == 1
        controls = healthy_control_strict == 1
      disease_negative:
        cases = disease == 1
        controls = disease == 0
    """
    if disease not in df.columns:
        return None, {
            "disease": disease,
            "control_mode": control_mode,
            "reason": "disease_column_missing",
        }

    tmp = df.copy()

    if control_mode == "healthy_control":
        tmp["analysis_group"] = np.nan
        tmp.loc[tmp[disease] == 1, "analysis_group"] = 1.0
        tmp.loc[(tmp[disease] != 1) & (tmp["healthy_control_strict"] == 1), "analysis_group"] = 0.0
    elif control_mode == "disease_negative":
        tmp["analysis_group"] = tmp[disease]
    else:
        raise ValueError(f"Unknown control_mode: {control_mode}")

    required = ["analysis_group", score_col, "RIDAGEYR", "RIAGENDR", "RIDRETH1", "DMDEDUC2", "INDFMPIR", "cycle"]
    required = [c for c in required if c in tmp.columns]

    tmp = tmp.dropna(subset=required).copy()
    tmp = prepare_covariates(tmp)

    n = tmp.shape[0]
    n_case = int((tmp["analysis_group"] == 1).sum())
    n_control = int((tmp["analysis_group"] == 0).sum())

    if n_case < min_case or n_control < min_control:
        return None, {
            "disease": disease,
            "control_mode": control_mode,
            "reason": "too_few_cases_or_controls",
            "n": n,
            "n_case": n_case,
            "n_control": n_control,
        }

    formula_terms = [score_col]

    if "RIDAGEYR" in tmp.columns:
        formula_terms.append("RIDAGEYR")
    if "RIAGENDR" in tmp.columns:
        formula_terms.append("C(RIAGENDR)")
    if "RIDRETH1" in tmp.columns:
        formula_terms.append("C(RIDRETH1)")
    if "DMDEDUC2" in tmp.columns:
        formula_terms.append("C(DMDEDUC2)")
    if "INDFMPIR" in tmp.columns:
        formula_terms.append("INDFMPIR")
    if "cycle" in tmp.columns:
        formula_terms.append("C(cycle)")

    formula = "analysis_group ~ " + " + ".join(formula_terms)

    try:
        if use_weights and "WTMEC2YR" in tmp.columns:
            weights = pd.to_numeric(tmp["WTMEC2YR"], errors="coerce")
            weights = weights / np.nanmean(weights)
            weights = weights.replace([np.inf, -np.inf], np.nan).fillna(1.0)

            model = smf.glm(
                formula=formula,
                data=tmp,
                family=sm.families.Binomial(),
                freq_weights=weights,
            )
        else:
            model = smf.glm(
                formula=formula,
                data=tmp,
                family=sm.families.Binomial(),
            )

        fit = model.fit(cov_type="HC3", maxiter=200)

        beta = fit.params[score_col]
        se = fit.bse[score_col]
        p = fit.pvalues[score_col]
        ci_low = beta - 1.96 * se
        ci_high = beta + 1.96 * se

        result = {
            "disease": disease,
            "control_mode": control_mode,
            "n": n,
            "n_case": n_case,
            "n_control": n_control,
            "case_fraction": n_case / n,
            "score_col": score_col,
            "beta_log_or_per_1sd": beta,
            "se": se,
            "or_per_1sd": float(np.exp(beta)),
            "or_lower_95": float(np.exp(ci_low)),
            "or_upper_95": float(np.exp(ci_high)),
            "p": p,
            "formula": formula,
            "weighted": bool(use_weights and "WTMEC2YR" in tmp.columns),
            "status": "ok",
        }

        return result, None

    except Exception as exc:
        return None, {
            "disease": disease,
            "control_mode": control_mode,
            "reason": "model_failed",
            "error": str(exc),
            "n": n,
            "n_case": n_case,
            "n_control": n_control,
            "formula": formula,
        }


def run_logistic_associations(
    merged: pd.DataFrame,
    disease_mapping: pd.DataFrame,
    score_col: str,
    min_case: int,
    min_control: int,
    use_weights: bool,
):
    disease_list = disease_mapping.loc[
        disease_mapping["n_available_variables"] > 0, "disease"
    ].tolist()

    result_rows = []
    skipped_rows = []

    for mode in ["healthy_control", "disease_negative"]:
        print(f"\nRunning logistic associations: {mode}")

        for disease in disease_list:
            if disease == "major_disease_composite":
                continue

            result, skipped = fit_logistic_one(
                df=merged,
                disease=disease,
                score_col=score_col,
                control_mode=mode,
                min_case=min_case,
                min_control=min_control,
                use_weights=use_weights,
            )

            if result is not None:
                result_rows.append(result)
            if skipped is not None:
                skipped_rows.append(skipped)

    results = pd.DataFrame(result_rows)
    skipped = pd.DataFrame(skipped_rows)

    if not results.empty:
        results["q_bh"] = np.nan
        for mode in results["control_mode"].unique():
            idx = results["control_mode"] == mode
            results.loc[idx, "q_bh"] = p_adjust_bh(results.loc[idx, "p"].values)

        results = results.sort_values(
            ["control_mode", "p", "disease"],
            ascending=[True, True, True],
        )

    return results, skipped


# ============================================================
# Optional survival analysis if external disease-onset table exists
# ============================================================

def run_optional_survival(onset_table_path: Path, scores: pd.DataFrame, score_col: str, min_case: int):
    """
    Optional prospective disease-onset Cox analysis.

    Required onset table columns:
      SEQN
      disease
      time
      event

    time should be prospective follow-up time from baseline to disease onset
    or censoring. event should be 1 for incident disease, 0 for censoring.
    """
    if onset_table_path is None:
        return pd.DataFrame(), pd.DataFrame([{
            "status": "not_run",
            "reason": "No onset table provided. Standard public NHANES does not include prospective nonfatal disease-onset follow-up."
        }])

    if not onset_table_path.exists():
        return pd.DataFrame(), pd.DataFrame([{
            "status": "not_run",
            "reason": f"Onset table path does not exist: {onset_table_path}"
        }])

    try:
        onset = pd.read_csv(onset_table_path, sep=None, engine="python")
    except Exception as exc:
        return pd.DataFrame(), pd.DataFrame([{
            "status": "not_run",
            "reason": "Could not read onset table",
            "error": str(exc),
        }])

    required = {"SEQN", "disease", "time", "event"}
    missing = required - set(onset.columns)

    if missing:
        return pd.DataFrame(), pd.DataFrame([{
            "status": "not_run",
            "reason": f"Onset table missing required columns: {sorted(missing)}"
        }])

    dat = onset.merge(scores, on="SEQN", how="inner")
    dat["time"] = pd.to_numeric(dat["time"], errors="coerce")
    dat["event"] = pd.to_numeric(dat["event"], errors="coerce")

    try:
        from statsmodels.duration.hazard_regression import PHReg
    except Exception as exc:
        return pd.DataFrame(), pd.DataFrame([{
            "status": "not_run",
            "reason": "statsmodels PHReg unavailable",
            "error": str(exc),
        }])

    result_rows = []
    skipped_rows = []

    for disease, d in dat.groupby("disease"):
        d = d.dropna(subset=["time", "event", score_col, "RIDAGEYR", "RIAGENDR", "RIDRETH1", "DMDEDUC2", "INDFMPIR", "cycle"]).copy()
        d = d[(d["time"] > 0) & d["event"].isin([0, 1])]

        n = d.shape[0]
        n_event = int((d["event"] == 1).sum())

        if n_event < min_case:
            skipped_rows.append({
                "disease": disease,
                "reason": "too_few_incident_events",
                "n": n,
                "n_event": n_event,
            })
            continue

        covars = d[[score_col, "RIDAGEYR", "INDFMPIR"]].copy()

        cat = pd.get_dummies(
            d[["RIAGENDR", "RIDRETH1", "DMDEDUC2", "cycle"]].astype(str),
            drop_first=True,
            dtype=float,
        )

        X = pd.concat([covars, cat], axis=1)
        X = X.replace([np.inf, -np.inf], np.nan).dropna()
        d2 = d.loc[X.index].copy()

        try:
            model = PHReg(
                endog=d2["time"].astype(float),
                exog=X.astype(float),
                status=d2["event"].astype(int),
            )
            fit = model.fit()

            param_names = list(X.columns)
            score_idx = param_names.index(score_col)

            beta = fit.params[score_idx]
            se = fit.bse[score_idx]
            p = fit.pvalues[score_idx]
            ci_low = beta - 1.96 * se
            ci_high = beta + 1.96 * se

            result_rows.append({
                "disease": disease,
                "n": n,
                "n_event": n_event,
                "beta_log_hr_per_1sd": beta,
                "se": se,
                "hr_per_1sd": float(np.exp(beta)),
                "hr_lower_95": float(np.exp(ci_low)),
                "hr_upper_95": float(np.exp(ci_high)),
                "p": p,
                "status": "ok",
            })

        except Exception as exc:
            skipped_rows.append({
                "disease": disease,
                "reason": "cox_model_failed",
                "error": str(exc),
                "n": n,
                "n_event": n_event,
            })

    results = pd.DataFrame(result_rows)
    skipped = pd.DataFrame(skipped_rows)

    if not results.empty:
        results["q_bh"] = p_adjust_bh(results["p"].values)
        results = results.sort_values(["p", "disease"])

    return results, skipped


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--nhanes_root",
        default="/Users/hao/Dropbox/NHANES",
        help="Root NHANES directory, e.g. /Users/hao/Dropbox/NHANES",
    )
    parser.add_argument(
        "--epoch_dir",
        default="/Users/hao/Dropbox/NHANES/output_dir/model2_nondisease_mortality_epoch",
        help="STEP II output directory containing nhanes_model2_epoch_scores.tsv",
    )
    parser.add_argument(
        "--outdir",
        default="/Users/hao/Dropbox/NHANES/output_dir/model2_nondisease_mortality_epoch/step3_disease_associations",
        help="Output directory for STEP III disease association analyses",
    )
    parser.add_argument(
        "--score_col",
        default=DEFAULT_SCORE_COL,
        help="Primary EPOCH score column. Default: mortality_epoch_acceleration_z",
    )
    parser.add_argument(
        "--min_case",
        type=int,
        default=20,
        help="Minimum number of cases required per disease endpoint",
    )
    parser.add_argument(
        "--min_control",
        type=int,
        default=20,
        help="Minimum number of controls required per disease endpoint",
    )
    parser.add_argument(
        "--use_weights",
        action="store_true",
        help="Use WTMEC2YR as frequency weights in logistic GLM sensitivity analysis",
    )
    parser.add_argument(
        "--onset_table",
        default=None,
        help="Optional prospective disease-onset table with columns: SEQN,disease,time,event",
    )

    args = parser.parse_args()

    nhanes_root = Path(args.nhanes_root)
    epoch_dir = Path(args.epoch_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    score_file = epoch_dir / "nhanes_model2_epoch_scores.tsv"
    if not score_file.exists():
        raise FileNotFoundError(f"Missing STEP II score file: {score_file}")

    print("=" * 80)
    print("STEP III: Link mortality EPOCH acceleration to disease status")
    print("=" * 80)
    print(f"NHANES root: {nhanes_root}")
    print(f"Epoch dir:   {epoch_dir}")
    print(f"Outdir:      {outdir}")
    print(f"Score col:   {args.score_col}")
    print("=" * 80)

    # --------------------------------------------------------
    # Load STEP II EPOCH scores
    # --------------------------------------------------------
    scores = pd.read_csv(score_file, sep="\t", low_memory=False)

    if args.score_col not in scores.columns:
        raise ValueError(f"Score column not found in STEP II score file: {args.score_col}")

    keep_score_cols = [
        "SEQN",
        "cycle",
        "split",
        "death",
        "followup_years_exm",
        "RIDAGEYR",
        "RIAGENDR",
        "RIDRETH1",
        "DMDEDUC2",
        "INDFMPIR",
        "WTMEC2YR",
        "SDMVPSU",
        "SDMVSTRA",
        args.score_col,
        "mortality_epoch_lp_total",
        "mortality_epoch_lp_feature_only",
        "mortality_epoch_year_equivalent",
        "mortality_epoch_acceleration_years",
    ]

    keep_score_cols = [c for c in keep_score_cols if c in scores.columns]
    scores = scores[keep_score_cols].copy()

    scores["SEQN"] = pd.to_numeric(scores["SEQN"], errors="coerce")
    scores = scores.dropna(subset=["SEQN", args.score_col]).copy()
    scores["SEQN"] = scores["SEQN"].astype(int)

    print(f"Loaded STEP II scores: {scores.shape[0]} participants")

    # --------------------------------------------------------
    # Load excluded disease modules from raw NHANES
    # --------------------------------------------------------
    disease_raw, source_table, label_table = load_disease_modules(nhanes_root)

    disease_raw["SEQN"] = pd.to_numeric(disease_raw["SEQN"], errors="coerce")
    disease_raw = disease_raw.dropna(subset=["SEQN"]).copy()
    disease_raw["SEQN"] = disease_raw["SEQN"].astype(int)

    source_table.to_csv(outdir / "disease_module_sources.tsv", sep="\t", index=False)
    label_table.to_csv(outdir / "disease_module_variable_labels.tsv", sep="\t", index=False)

    print(f"Loaded disease-module rows: {disease_raw.shape[0]}")
    print(f"Loaded disease-module columns: {disease_raw.shape[1]}")

    # --------------------------------------------------------
    # Build harmonized disease endpoints
    # --------------------------------------------------------
    disease_status, disease_mapping, healthy_def = build_disease_status(disease_raw)

    disease_status.to_csv(outdir / "disease_status_wide.tsv", sep="\t", index=False)
    disease_mapping.to_csv(outdir / "disease_endpoint_mapping_and_counts.tsv", sep="\t", index=False)
    healthy_def.to_csv(outdir / "healthy_control_definition.tsv", sep="\t", index=False)

    print("\nDisease endpoint counts:")
    print(disease_mapping[[
        "disease",
        "n_available_variables",
        "n_nonmissing",
        "n_cases",
        "n_controls_disease_negative",
    ]].to_string(index=False))

    # --------------------------------------------------------
    # Merge EPOCH scores with disease status
    # --------------------------------------------------------
    merged = scores.merge(
        disease_status.drop(columns=["cycle", "cycle_begin_year"], errors="ignore"),
        on="SEQN",
        how="inner",
    )

    print(f"\nMerged EPOCH + disease status: {merged.shape[0]} participants")

    merged.to_csv(
        outdir / "epoch_scores_with_disease_status.tsv",
        sep="\t",
        index=False,
    )

    # --------------------------------------------------------
    # Logistic associations
    # --------------------------------------------------------
    logistic_results, skipped_logistic = run_logistic_associations(
        merged=merged,
        disease_mapping=disease_mapping,
        score_col=args.score_col,
        min_case=args.min_case,
        min_control=args.min_control,
        use_weights=args.use_weights,
    )

    if not logistic_results.empty:
        logistic_results = logistic_results.merge(
            disease_mapping[["disease", "description", "category", "available_variables_used"]],
            on="disease",
            how="left",
        )

    logistic_healthy = logistic_results[
        logistic_results["control_mode"] == "healthy_control"
    ].copy() if not logistic_results.empty else pd.DataFrame()

    logistic_negative = logistic_results[
        logistic_results["control_mode"] == "disease_negative"
    ].copy() if not logistic_results.empty else pd.DataFrame()

    logistic_healthy.to_csv(
        outdir / "logistic_results_healthy_control.tsv",
        sep="\t",
        index=False,
    )

    logistic_negative.to_csv(
        outdir / "logistic_results_disease_negative.tsv",
        sep="\t",
        index=False,
    )

    logistic_results.to_csv(
        outdir / "disease_association_summary_all.tsv",
        sep="\t",
        index=False,
    )

    skipped_logistic.to_csv(
        outdir / "skipped_endpoints_logistic.tsv",
        sep="\t",
        index=False,
    )

    # --------------------------------------------------------
    # Optional disease-onset survival
    # --------------------------------------------------------
    onset_path = Path(args.onset_table) if args.onset_table is not None else None

    survival_results, skipped_survival = run_optional_survival(
        onset_table_path=onset_path,
        scores=scores,
        score_col=args.score_col,
        min_case=args.min_case,
    )

    survival_results.to_csv(
        outdir / "survival_results_optional.tsv",
        sep="\t",
        index=False,
    )

    skipped_survival.to_csv(
        outdir / "skipped_endpoints_survival.tsv",
        sep="\t",
        index=False,
    )

    # --------------------------------------------------------
    # Save config and compact report
    # --------------------------------------------------------
    config = {
        "analysis": "STEP III disease association analysis",
        "nhanes_root": str(nhanes_root),
        "epoch_dir": str(epoch_dir),
        "outdir": str(outdir),
        "score_col": args.score_col,
        "min_case": args.min_case,
        "min_control": args.min_control,
        "use_weights": args.use_weights,
        "disease_module_patterns_loaded": DISEASE_MODULE_PATTERNS,
        "main_model": "Logistic regression: disease status ~ mortality_epoch_acceleration_z + age + sex + race/ethnicity + education + income + cycle",
        "primary_control": "strict healthy controls",
        "sensitivity_control": "disease-negative controls",
        "survival_note": "Prospective disease-onset survival is only run if an external onset table is provided. Standard public NHANES does not provide prospective nonfatal disease-onset follow-up.",
    }

    with open(outdir / "step3_config.json", "w") as f:
        json.dump(config, f, indent=2)

    summary_rows = []

    if not logistic_healthy.empty:
        top_healthy = logistic_healthy.sort_values("p").head(10)
        summary_rows.append("Top healthy-control associations:")
        summary_rows.append(top_healthy[
            ["disease", "category", "n", "n_case", "n_control", "or_per_1sd", "or_lower_95", "or_upper_95", "p", "q_bh"]
        ].to_string(index=False))

    if not logistic_negative.empty:
        top_negative = logistic_negative.sort_values("p").head(10)
        summary_rows.append("\nTop disease-negative-control associations:")
        summary_rows.append(top_negative[
            ["disease", "category", "n", "n_case", "n_control", "or_per_1sd", "or_lower_95", "or_upper_95", "p", "q_bh"]
        ].to_string(index=False))

    with open(outdir / "step3_text_summary.txt", "w") as f:
        f.write("\n".join(summary_rows))

    print("\nDone.")
    print("Main outputs:")
    print(f"  {outdir / 'logistic_results_healthy_control.tsv'}")
    print(f"  {outdir / 'logistic_results_disease_negative.tsv'}")
    print(f"  {outdir / 'disease_association_summary_all.tsv'}")
    print(f"  {outdir / 'epoch_scores_with_disease_status.tsv'}")
    print(f"  {outdir / 'step3_text_summary.txt'}")

    if not logistic_healthy.empty:
        print("\nTop healthy-control associations:")
        print(logistic_healthy.sort_values("p").head(10)[[
            "disease",
            "category",
            "n",
            "n_case",
            "n_control",
            "or_per_1sd",
            "or_lower_95",
            "or_upper_95",
            "p",
            "q_bh",
        ]].to_string(index=False))


if __name__ == "__main__":
    main()