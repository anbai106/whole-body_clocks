#!/usr/bin/env python3
# ============================================================
# Function-based QC for L'EPOCH / mortality-clock year scaling
#
# Purpose:
#   Check whether a clock can be safely transformed into
#   acceleration years using the age coefficient from
#   *_coefficients.tsv.
#
# Main QC logic:
#   acceleration_years = centered risk_score / age_beta
#
# A clock is flagged if:
#   1. age_beta is missing
#   2. age_beta <= 0
#   3. abs(age_beta) is close to zero
#   4. recomputed acceleration-year distribution is implausibly wide
#
# General use:
#   run_lepoch_year_scale_qc(
#       base_dir="/path/to/WholeBodyClock",
#       clock_dir_glob="*_mi_clock",
#       analysis_label="mi"
#   )
# ============================================================

import os
import re
import glob
import json
import warnings
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")


# ============================================================
# 1. Basic readers and helpers
# ============================================================

def safe_read_tsv(path: str) -> Optional[pd.DataFrame]:
    if path is None or not os.path.exists(path):
        return None

    try:
        return pd.read_csv(path, sep="\t", low_memory=False)
    except Exception:
        try:
            return pd.read_csv(path, sep=None, engine="python", low_memory=False)
        except Exception as e:
            print(f"WARNING: failed to read {path}: {e}")
            return None


def safe_read_json(path: str) -> Dict:
    if path is None or not os.path.exists(path):
        return {}

    try:
        with open(path, "r") as f:
            return json.load(f)
    except Exception as e:
        print(f"WARNING: failed to read JSON {path}: {e}")
        return {}


def as_numeric(x) -> pd.Series:
    return pd.to_numeric(x, errors="coerce")


def first_existing_col(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def find_col_by_regex(df: pd.DataFrame, patterns: List[str]) -> Optional[str]:
    cols = list(df.columns)

    for pat in patterns:
        hits = [
            c for c in cols
            if re.search(pat, str(c), flags=re.IGNORECASE)
        ]
        if len(hits) > 0:
            return hits[0]

    return None


def find_unique_file(clock_dir: str, pattern: str) -> Optional[str]:
    files = sorted(glob.glob(os.path.join(clock_dir, pattern)))

    if len(files) == 0:
        return None

    nonempty = [f for f in files if os.path.getsize(f) > 0]

    if len(nonempty) > 0:
        return nonempty[0]

    return files[0]


def infer_modality_from_folder(folder_name: str) -> str:
    x = folder_name.lower()

    if "_mri_" in x:
        return "MRI"
    if "_proteomics_" in x:
        return "Proteomics"
    if "_metabolomics_" in x:
        return "Metabolomics"

    return "Unknown"


def infer_organ_from_folder(folder_name: str, analysis_label: str) -> str:
    x = folder_name

    # Remove terminal clock token.
    x = re.sub(r"_clock$", "", x, flags=re.IGNORECASE)

    # Remove analysis label, for example _mi, _copd, _dementia, _mortality.
    if analysis_label:
        x = re.sub(
            rf"_{re.escape(analysis_label)}$",
            "",
            x,
            flags=re.IGNORECASE
        )

    # Remove modality token.
    x = re.sub(r"_mri$", "", x, flags=re.IGNORECASE)
    x = re.sub(r"_proteomics$", "", x, flags=re.IGNORECASE)
    x = re.sub(r"_metabolomics$", "", x, flags=re.IGNORECASE)

    x = x.replace("_", " ").strip()

    if len(x) == 0:
        return folder_name

    return x[:1].upper() + x[1:]


# ============================================================
# 2. Coefficient table parsing
# ============================================================

def detect_coefficient_columns(
    coef_df: pd.DataFrame
) -> Tuple[Optional[str], Optional[str]]:
    """
    Detect term and coefficient columns from *_coefficients.tsv.
    """

    term_candidates = [
        "feature",
        "term",
        "variable",
        "covariate",
        "name",
        "Feature",
        "Variable",
        "parameter",
        "coef_name",
    ]

    beta_candidates = [
        "coefficient",
        "coef",
        "beta",
        "Coefficient",
        "Coef",
        "Beta",
        "estimate",
        "value",
    ]

    term_col = first_existing_col(coef_df, term_candidates)
    beta_col = first_existing_col(coef_df, beta_candidates)

    if term_col is None:
        term_col = find_col_by_regex(
            coef_df,
            patterns=[
                r"feature",
                r"term",
                r"variable",
                r"covariate",
                r"parameter",
                r"name",
            ],
        )

    if beta_col is None:
        beta_col = find_col_by_regex(
            coef_df,
            patterns=[
                r"coef",
                r"beta",
                r"estimate",
                r"value",
            ],
        )

    return term_col, beta_col


def identify_age_terms(
    coef_df: pd.DataFrame,
    term_col: str,
    beta_col: str,
    age_term_regex: str
) -> pd.DataFrame:
    d = coef_df.copy()

    d["term_string"] = d[term_col].astype(str)
    d["beta_numeric"] = as_numeric(d[beta_col])

    d["is_age_like"] = d["term_string"].str.contains(
        age_term_regex,
        case=False,
        regex=True,
        na=False,
    )

    d = d[d["is_age_like"]].copy()

    if d.empty:
        return d

    def priority(term: str) -> int:
        t = str(term).lower()

        exact_patterns = [
            r"^age$",
            r"^age_at_baseline$",
            r"^age_at_imaging$",
            r"^age_at_clock$",
            r"^chronological_age$",
            r"age_recruitment",
            r"age_at_recruitment",
        ]

        for i, pat in enumerate(exact_patterns):
            if re.search(pat, t):
                return i

        return 100

    d["age_term_priority"] = d["term_string"].apply(priority)

    d = d.sort_values(
        ["age_term_priority", "term_string"]
    ).copy()

    return d


def choose_age_beta(age_terms: pd.DataFrame) -> Tuple[float, str, str]:
    if age_terms is None or age_terms.empty:
        return np.nan, "", "no_age_term_found"

    valid = age_terms[age_terms["beta_numeric"].notna()].copy()

    if valid.empty:
        return np.nan, "", "age_terms_found_but_beta_missing"

    selected = valid.sort_values(
        ["age_term_priority", "term_string"]
    ).iloc[0]

    return (
        float(selected["beta_numeric"]),
        str(selected["term_string"]),
        "selected_highest_priority_age_term",
    )


# ============================================================
# 3. Prediction table parsing
# ============================================================

def detect_prediction_columns(
    pred_df: pd.DataFrame,
    analysis_label: str = ""
) -> Dict[str, Optional[str]]:
    cols = list(pred_df.columns)

    risk_score_col = None
    accel_year_col = None
    accel_z_col = None
    split_col = "split" if "split" in cols else None

    # Risk score column.
    risk_candidates = [
        c for c in cols
        if re.search(r"risk_score$", str(c), flags=re.IGNORECASE)
    ]

    if analysis_label:
        label_candidates = [
            c for c in risk_candidates
            if analysis_label.lower() in str(c).lower()
        ]
    else:
        label_candidates = []

    if len(label_candidates) > 0:
        risk_score_col = label_candidates[0]
    elif len(risk_candidates) > 0:
        risk_score_col = risk_candidates[0]
    else:
        risk_score_col = find_col_by_regex(
            pred_df,
            patterns=[
                r"risk.*score",
                r"linear.*predictor",
                r"lp$",
            ],
        )

    # Acceleration years.
    accel_year_candidates = [
        c for c in cols
        if re.search(
            r"acceleration.*years$",
            str(c),
            flags=re.IGNORECASE,
        )
    ]

    if analysis_label:
        label_candidates = [
            c for c in accel_year_candidates
            if analysis_label.lower() in str(c).lower()
        ]
    else:
        label_candidates = []

    if len(label_candidates) > 0:
        accel_year_col = label_candidates[0]
    elif len(accel_year_candidates) > 0:
        accel_year_col = accel_year_candidates[0]

    # Acceleration z.
    accel_z_candidates = [
        c for c in cols
        if re.search(
            r"acceleration.*z$",
            str(c),
            flags=re.IGNORECASE,
        )
    ]

    if analysis_label:
        label_candidates = [
            c for c in accel_z_candidates
            if analysis_label.lower() in str(c).lower()
        ]
    else:
        label_candidates = []

    if len(label_candidates) > 0:
        accel_z_col = label_candidates[0]
    elif len(accel_z_candidates) > 0:
        accel_z_col = accel_z_candidates[0]

    return {
        "risk_score_col": risk_score_col,
        "accel_year_col": accel_year_col,
        "accel_z_col": accel_z_col,
        "split_col": split_col,
    }


def summarize_numeric(x: Union[pd.Series, np.ndarray], prefix: str) -> Dict:
    y = as_numeric(pd.Series(x))
    y = y[np.isfinite(y)]

    if len(y) == 0:
        return {
            f"{prefix}_n": 0,
            f"{prefix}_mean": np.nan,
            f"{prefix}_sd": np.nan,
            f"{prefix}_min": np.nan,
            f"{prefix}_p001": np.nan,
            f"{prefix}_p01": np.nan,
            f"{prefix}_p05": np.nan,
            f"{prefix}_median": np.nan,
            f"{prefix}_p95": np.nan,
            f"{prefix}_p99": np.nan,
            f"{prefix}_p999": np.nan,
            f"{prefix}_max": np.nan,
            f"{prefix}_iqr": np.nan,
            f"{prefix}_p01_p99_range": np.nan,
        }

    q = y.quantile(
        [0, 0.001, 0.01, 0.05, 0.5, 0.95, 0.99, 0.999, 1.0]
    )

    return {
        f"{prefix}_n": int(len(y)),
        f"{prefix}_mean": float(y.mean()),
        f"{prefix}_sd": float(y.std(ddof=1)) if len(y) > 1 else np.nan,
        f"{prefix}_min": float(q.loc[0.0]),
        f"{prefix}_p001": float(q.loc[0.001]),
        f"{prefix}_p01": float(q.loc[0.01]),
        f"{prefix}_p05": float(q.loc[0.05]),
        f"{prefix}_median": float(q.loc[0.5]),
        f"{prefix}_p95": float(q.loc[0.95]),
        f"{prefix}_p99": float(q.loc[0.99]),
        f"{prefix}_p999": float(q.loc[0.999]),
        f"{prefix}_max": float(q.loc[1.0]),
        f"{prefix}_iqr": float(y.quantile(0.75) - y.quantile(0.25)),
        f"{prefix}_p01_p99_range": float(q.loc[0.99] - q.loc[0.01]),
    }


def prediction_scale_qc(
    pred_path: Optional[str],
    age_beta: float,
    analysis_label: str = ""
) -> Dict:
    out = {
        "prediction_file_exists": bool(pred_path and os.path.exists(pred_path)),
        "risk_score_col_used": "",
        "existing_accel_year_col_used": "",
        "risk_score_center_value_used": np.nan,
    }

    if pred_path is None or not os.path.exists(pred_path):
        return out

    pred = safe_read_tsv(pred_path)

    if pred is None or pred.empty:
        return out

    detected = detect_prediction_columns(
        pred,
        analysis_label=analysis_label,
    )

    risk_col = detected["risk_score_col"]
    accel_year_col = detected["accel_year_col"]
    split_col = detected["split_col"]

    out["risk_score_col_used"] = "" if risk_col is None else risk_col
    out["existing_accel_year_col_used"] = "" if accel_year_col is None else accel_year_col

    if risk_col is not None and risk_col in pred.columns:
        score = as_numeric(pred[risk_col]).replace([np.inf, -np.inf], np.nan)

        out.update(summarize_numeric(score, "risk_score"))

        # Center using training-set median when split exists.
        if split_col is not None:
            train_mask = pred[split_col].astype(str).str.lower() == "train"
            if train_mask.sum() > 0:
                center_value = score[train_mask].median(skipna=True)
            else:
                center_value = score.median(skipna=True)
        else:
            center_value = score.median(skipna=True)

        out["risk_score_center_value_used"] = (
            float(center_value) if np.isfinite(center_value) else np.nan
        )

        if np.isfinite(age_beta) and age_beta != 0:
            recomputed_years = (score - center_value) / age_beta
            out.update(
                summarize_numeric(
                    recomputed_years,
                    "recomputed_accel_years",
                )
            )
        else:
            out.update(
                summarize_numeric(
                    pd.Series([], dtype=float),
                    "recomputed_accel_years",
                )
            )

    if accel_year_col is not None and accel_year_col in pred.columns:
        out.update(
            summarize_numeric(
                pred[accel_year_col],
                "existing_accel_years",
            )
        )
    else:
        out.update(
            summarize_numeric(
                pd.Series([], dtype=float),
                "existing_accel_years",
            )
        )

    return out


# ============================================================
# 4. Delta-C-index, performance, and nonzero features
# ============================================================

def read_delta_cindex(delta_path: Optional[str]) -> Dict:
    out = {
        "delta_file_exists": bool(delta_path and os.path.exists(delta_path)),
        "delta_cindex": np.nan,
        "delta_cindex_ci_lower": np.nan,
        "delta_cindex_ci_upper": np.nan,
        "delta_cindex_p": np.nan,
        "delta_cindex_significant_positive": False,
    }

    if delta_path is None or not os.path.exists(delta_path):
        return out

    d = safe_read_tsv(delta_path)

    if d is None or d.empty:
        return out

    def get_first(candidates):
        for c in candidates:
            if c in d.columns:
                val = as_numeric(d[c]).dropna()
                if len(val) > 0:
                    return float(val.iloc[0])
        return np.nan

    out["delta_cindex"] = get_first(
        ["delta_cindex", "delta_cindex_test_M3_vs_M1"]
    )

    out["delta_cindex_ci_lower"] = get_first(
        ["delta_cindex_ci_lower", "delta_cindex_test_M3_vs_M1_ci_lower"]
    )

    out["delta_cindex_ci_upper"] = get_first(
        ["delta_cindex_ci_upper", "delta_cindex_test_M3_vs_M1_ci_upper"]
    )

    out["delta_cindex_p"] = get_first(
        [
            "empirical_p_two_sided_delta_not_equal_0",
            "delta_cindex_test_M3_vs_M1_p_two_sided",
            "p",
            "p_value",
        ]
    )

    if np.isfinite(out["delta_cindex_ci_lower"]):
        out["delta_cindex_significant_positive"] = (
            out["delta_cindex_ci_lower"] > 0
        )
    elif np.isfinite(out["delta_cindex"]) and np.isfinite(out["delta_cindex_p"]):
        out["delta_cindex_significant_positive"] = (
            out["delta_cindex"] > 0 and out["delta_cindex_p"] < 0.05
        )

    return out


def read_performance(perf_path: Optional[str]) -> Dict:
    perf = safe_read_json(perf_path)

    def get_num(field):
        try:
            return float(perf.get(field, np.nan))
        except Exception:
            return np.nan

    return {
        "performance_file_exists": bool(perf_path and os.path.exists(perf_path)),
        "cindex_train_json": get_num("cindex_train"),
        "cindex_validation_json": get_num("cindex_validation"),
        "cindex_test_json": get_num("cindex_test"),
        "cindex_test_M1_covariate_baseline_json": get_num(
            "cindex_test_M1_covariate_baseline"
        ),
        "admin_censor_date": perf.get("admin_censor_date", ""),
    }


def count_nonzero_features(nonzero_path: Optional[str]) -> Dict:
    out = {
        "nonzero_file_exists": bool(nonzero_path and os.path.exists(nonzero_path)),
        "n_nonzero_rows": np.nan,
        "n_nonzero_non_age_rows": np.nan,
    }

    if nonzero_path is None or not os.path.exists(nonzero_path):
        return out

    d = safe_read_tsv(nonzero_path)

    if d is None:
        return out

    out["n_nonzero_rows"] = int(d.shape[0])

    term_col, _ = detect_coefficient_columns(d)

    if term_col is not None:
        is_age = d[term_col].astype(str).str.contains(
            r"(^|[^a-zA-Z])age([^a-zA-Z]|$)|age_at_|chronological_age",
            case=False,
            regex=True,
            na=False,
        )
        out["n_nonzero_non_age_rows"] = int((~is_age).sum())
    else:
        out["n_nonzero_non_age_rows"] = int(d.shape[0])

    return out


# ============================================================
# 5. Clock-level QC
# ============================================================

def apply_year_scale_flags(
    row: Dict,
    age_beta_min_abs: float = 0.005,
    year_sd_warn: float = 20.0,
    year_iqr_warn: float = 30.0,
    year_p01_p99_warn: float = 100.0,
) -> Dict:
    age_beta = row.get("age_beta", np.nan)

    if not np.isfinite(age_beta):
        row["final_year_scale_qc_status"] = "FAIL_missing_age_beta"
        row["final_year_scale_qc_reason"] = "No numeric age beta found."
        return row

    if age_beta <= 0:
        row["final_year_scale_qc_status"] = "FAIL_nonpositive_age_beta"
        row["final_year_scale_qc_reason"] = (
            "Age beta is non-positive. Year transformation is biologically invalid "
            "if higher chronological age does not increase hazard."
        )
        return row

    if abs(age_beta) < age_beta_min_abs:
        row["final_year_scale_qc_status"] = "FAIL_age_beta_near_zero"
        row["final_year_scale_qc_reason"] = (
            f"abs(age_beta) < {age_beta_min_abs}. "
            "Dividing by this coefficient can explode the year scale."
        )
        return row

    warnings_list = []

    year_sd = row.get("recomputed_accel_years_sd", np.nan)
    year_iqr = row.get("recomputed_accel_years_iqr", np.nan)
    year_p01_p99 = row.get("recomputed_accel_years_p01_p99_range", np.nan)

    if np.isfinite(year_sd) and year_sd > year_sd_warn:
        warnings_list.append(
            f"recomputed acceleration-year SD {year_sd:.3g} > {year_sd_warn}"
        )

    if np.isfinite(year_iqr) and year_iqr > year_iqr_warn:
        warnings_list.append(
            f"recomputed acceleration-year IQR {year_iqr:.3g} > {year_iqr_warn}"
        )

    if np.isfinite(year_p01_p99) and year_p01_p99 > year_p01_p99_warn:
        warnings_list.append(
            f"recomputed acceleration-year p01-p99 range {year_p01_p99:.3g} > {year_p01_p99_warn}"
        )

    if warnings_list:
        row["final_year_scale_qc_status"] = "WARN_year_distribution_too_wide"
        row["final_year_scale_qc_reason"] = "; ".join(warnings_list)
    else:
        row["final_year_scale_qc_status"] = "PASS_year_scale_qc"
        row["final_year_scale_qc_reason"] = (
            "Age beta and recomputed year distribution pass chosen QC thresholds."
        )

    return row


def qc_one_clock_dir(
    clock_dir: str,
    analysis_label: str = "",
    age_term_regex: str = (
        r"(^|[^a-zA-Z])age([^a-zA-Z]|$)|"
        r"age_at_|"
        r"chronological_age|"
        r"Age_recruitment|"
        r"age_at_baseline|"
        r"age_at_imaging"
    ),
    age_beta_min_abs: float = 0.005,
    year_sd_warn: float = 20.0,
    year_iqr_warn: float = 30.0,
    year_p01_p99_warn: float = 100.0,
    skip_prediction_qc: bool = False,
) -> Tuple[Dict, pd.DataFrame]:
    """
    Run age-beta / year-scale QC for one clock folder.
    """

    folder_name = os.path.basename(clock_dir.rstrip("/"))

    coeff_file = find_unique_file(clock_dir, "*_coefficients.tsv")
    nonzero_file = find_unique_file(clock_dir, "*_nonzero_coefficients.tsv")
    pred_file = find_unique_file(clock_dir, "*_predictions.tsv")
    perf_file = find_unique_file(clock_dir, "*_performance.json")
    delta_file = find_unique_file(clock_dir, "*_incremental_value_delta_cindex.tsv")

    row = {
        "clock_folder": folder_name,
        "clock_dir": clock_dir,
        "analysis_label": analysis_label,
        "organ_label": infer_organ_from_folder(folder_name, analysis_label),
        "modality": infer_modality_from_folder(folder_name),
        "coeff_file": coeff_file or "",
        "nonzero_file": nonzero_file or "",
        "prediction_file": pred_file or "",
        "performance_file": perf_file or "",
        "delta_file": delta_file or "",
        "coeff_file_exists": bool(coeff_file and os.path.exists(coeff_file)),
    }

    all_age_terms = pd.DataFrame()

    if coeff_file is None or not os.path.exists(coeff_file):
        row.update({
            "age_beta": np.nan,
            "selected_age_term": "",
            "age_beta_selection_note": "coefficients_file_missing",
            "n_age_terms_found": 0,
        })
    else:
        coef_df = safe_read_tsv(coeff_file)

        if coef_df is None or coef_df.empty:
            row.update({
                "age_beta": np.nan,
                "selected_age_term": "",
                "age_beta_selection_note": "coefficients_file_empty_or_unreadable",
                "n_age_terms_found": 0,
            })
        else:
            term_col, beta_col = detect_coefficient_columns(coef_df)

            row["coeff_term_col_used"] = term_col or ""
            row["coeff_beta_col_used"] = beta_col or ""

            if term_col is None or beta_col is None:
                row.update({
                    "age_beta": np.nan,
                    "selected_age_term": "",
                    "age_beta_selection_note": "could_not_detect_term_or_beta_column",
                    "n_age_terms_found": 0,
                })
            else:
                age_terms = identify_age_terms(
                    coef_df=coef_df,
                    term_col=term_col,
                    beta_col=beta_col,
                    age_term_regex=age_term_regex,
                )

                age_beta, selected_age_term, note = choose_age_beta(age_terms)

                row.update({
                    "age_beta": age_beta,
                    "abs_age_beta": abs(age_beta) if np.isfinite(age_beta) else np.nan,
                    "selected_age_term": selected_age_term,
                    "age_beta_selection_note": note,
                    "n_age_terms_found": int(age_terms.shape[0]),
                    "all_age_terms": ";".join(age_terms["term_string"].astype(str).tolist()) if not age_terms.empty else "",
                    "all_age_betas": ";".join(age_terms["beta_numeric"].astype(str).tolist()) if not age_terms.empty else "",
                    "years_per_1_log_hazard_unit": (
                        1.0 / age_beta
                        if np.isfinite(age_beta) and age_beta != 0
                        else np.nan
                    ),
                })

                if not age_terms.empty:
                    all_age_terms = age_terms.copy()
                    all_age_terms["clock_folder"] = folder_name
                    all_age_terms["clock_dir"] = clock_dir
                    all_age_terms["analysis_label"] = analysis_label
                    all_age_terms["organ_label"] = row["organ_label"]
                    all_age_terms["modality"] = row["modality"]
                    all_age_terms["coeff_file"] = coeff_file

    row.update(read_performance(perf_file))
    row.update(read_delta_cindex(delta_file))
    row.update(count_nonzero_features(nonzero_file))

    if not skip_prediction_qc:
        row.update(
            prediction_scale_qc(
                pred_path=pred_file,
                age_beta=row.get("age_beta", np.nan),
                analysis_label=analysis_label,
            )
        )

    row = apply_year_scale_flags(
        row=row,
        age_beta_min_abs=age_beta_min_abs,
        year_sd_warn=year_sd_warn,
        year_iqr_warn=year_iqr_warn,
        year_p01_p99_warn=year_p01_p99_warn,
    )

    return row, all_age_terms


# ============================================================
# 6. Main reusable function
# ============================================================

def run_lepoch_year_scale_qc(
    base_dir: Optional[str] = None,
    clock_dir_glob: Optional[str] = None,
    clock_dirs: Optional[List[str]] = None,
    out_dir: Optional[str] = None,
    analysis_label: str = "clock",
    age_beta_min_abs: float = 0.005,
    year_sd_warn: float = 20.0,
    year_iqr_warn: float = 30.0,
    year_p01_p99_warn: float = 100.0,
    age_term_regex: str = (
        r"(^|[^a-zA-Z])age([^a-zA-Z]|$)|"
        r"age_at_|"
        r"chronological_age|"
        r"Age_recruitment|"
        r"age_at_baseline|"
        r"age_at_imaging"
    ),
    skip_prediction_qc: bool = False,
    verbose: bool = True,
) -> Dict[str, Union[str, pd.DataFrame]]:
    """
    General QC runner.

    Parameters
    ----------
    base_dir:
        Directory containing clock folders.
        Example: /cbica/home/wenju/Reproducibile_paper/WholeBodyClock

    clock_dir_glob:
        Glob pattern relative to base_dir.
        Examples:
            "*_mi_clock"
            "*_copd_clock"
            "*_dementia_clock"
            "*_mortality_clock"

    clock_dirs:
        Optional explicit list of clock directories.
        If supplied, base_dir and clock_dir_glob are not required.

    out_dir:
        Output directory. If None, defaults to:
            <base_dir>/<analysis_label>_lepoch_year_scale_qc

    analysis_label:
        Label used in output file names.
        Examples:
            "mi", "copd", "dementia", "mortality"

    Returns
    -------
    dict with:
        summary_df
        age_terms_df
        problematic_df
        summary_path
        age_terms_path
        problematic_path
        metadata_path
    """

    if clock_dirs is None:
        if base_dir is None or clock_dir_glob is None:
            raise ValueError(
                "Provide either clock_dirs, or both base_dir and clock_dir_glob."
            )

        base_dir = os.path.abspath(base_dir)
        search_pattern = os.path.join(base_dir, clock_dir_glob)

        clock_dirs = sorted([
            d for d in glob.glob(search_pattern)
            if os.path.isdir(d)
        ])
    else:
        clock_dirs = [os.path.abspath(d) for d in clock_dirs]
        if base_dir is None:
            base_dir = os.path.commonpath(clock_dirs)
        search_pattern = "explicit_clock_dirs"

    if len(clock_dirs) == 0:
        raise RuntimeError(
            f"No clock directories found. base_dir={base_dir}, glob={clock_dir_glob}"
        )

    if out_dir is None:
        out_dir = os.path.join(
            os.path.abspath(base_dir),
            f"{analysis_label}_lepoch_year_scale_qc"
        )
    else:
        out_dir = os.path.abspath(out_dir)

    os.makedirs(out_dir, exist_ok=True)

    if verbose:
        print("============================================================")
        print("L'EPOCH / clock year-scale QC")
        print("============================================================")
        print("Base directory:", base_dir)
        print("Analysis label:", analysis_label)
        print("Search pattern:", search_pattern)
        print("Number of clock folders:", len(clock_dirs))
        print("Output directory:", out_dir)
        print("age_beta_min_abs:", age_beta_min_abs)
        print("year_sd_warn:", year_sd_warn)
        print("year_iqr_warn:", year_iqr_warn)
        print("year_p01_p99_warn:", year_p01_p99_warn)
        print("============================================================")

    rows = []
    age_term_tables = []

    for i, clock_dir in enumerate(clock_dirs, start=1):
        if verbose:
            print(f"[{i}/{len(clock_dirs)}] {os.path.basename(clock_dir)}")

        try:
            row, age_terms = qc_one_clock_dir(
                clock_dir=clock_dir,
                analysis_label=analysis_label,
                age_term_regex=age_term_regex,
                age_beta_min_abs=age_beta_min_abs,
                year_sd_warn=year_sd_warn,
                year_iqr_warn=year_iqr_warn,
                year_p01_p99_warn=year_p01_p99_warn,
                skip_prediction_qc=skip_prediction_qc,
            )
        except Exception as e:
            row = {
                "clock_folder": os.path.basename(clock_dir),
                "clock_dir": clock_dir,
                "analysis_label": analysis_label,
                "final_year_scale_qc_status": "FAIL_exception",
                "final_year_scale_qc_reason": str(e),
            }
            age_terms = pd.DataFrame()

        rows.append(row)

        if age_terms is not None and not age_terms.empty:
            age_term_tables.append(age_terms)

    summary_df = pd.DataFrame(rows)

    if len(age_term_tables) > 0:
        age_terms_df = pd.concat(age_term_tables, ignore_index=True)
    else:
        age_terms_df = pd.DataFrame()

    problematic_df = summary_df[
        summary_df["final_year_scale_qc_status"].astype(str).str.startswith("FAIL")
        | summary_df["final_year_scale_qc_status"].astype(str).str.startswith("WARN")
    ].copy()

    status_order = {
        "FAIL_exception": 0,
        "FAIL_missing_age_beta": 1,
        "FAIL_nonpositive_age_beta": 2,
        "FAIL_age_beta_near_zero": 3,
        "WARN_year_distribution_too_wide": 4,
        "PASS_year_scale_qc": 5,
    }

    summary_df["_status_order"] = (
        summary_df["final_year_scale_qc_status"]
        .map(status_order)
        .fillna(99)
    )

    sort_cols = [
        "_status_order",
        "modality",
        "organ_label",
        "clock_folder",
    ]
    sort_cols = [c for c in sort_cols if c in summary_df.columns]

    summary_df = (
        summary_df
        .sort_values(sort_cols, na_position="last")
        .drop(columns=["_status_order"])
    )

    summary_path = os.path.join(
        out_dir,
        f"{analysis_label}_year_scale_qc_summary.tsv"
    )

    age_terms_path = os.path.join(
        out_dir,
        f"{analysis_label}_all_age_terms.tsv"
    )

    problematic_path = os.path.join(
        out_dir,
        f"{analysis_label}_problematic_year_scale_clocks.tsv"
    )

    metadata_path = os.path.join(
        out_dir,
        f"{analysis_label}_year_scale_qc_metadata.json"
    )

    summary_df.to_csv(summary_path, sep="\t", index=False)
    age_terms_df.to_csv(age_terms_path, sep="\t", index=False)
    problematic_df.to_csv(problematic_path, sep="\t", index=False)

    metadata = {
        "base_dir": base_dir,
        "analysis_label": analysis_label,
        "search_pattern": search_pattern,
        "n_clock_dirs": len(clock_dirs),
        "out_dir": out_dir,
        "age_beta_min_abs": age_beta_min_abs,
        "year_sd_warn": year_sd_warn,
        "year_iqr_warn": year_iqr_warn,
        "year_p01_p99_warn": year_p01_p99_warn,
        "age_term_regex": age_term_regex,
        "skip_prediction_qc": skip_prediction_qc,
        "outputs": {
            "summary": summary_path,
            "age_terms": age_terms_path,
            "problematic": problematic_path,
        },
        "qc_status_counts": (
            summary_df["final_year_scale_qc_status"]
            .value_counts(dropna=False)
            .to_dict()
        ),
    }

    with open(metadata_path, "w") as f:
        json.dump(metadata, f, indent=2, default=str)

    if verbose:
        print("\n============================================================")
        print("Finished.")
        print("Outputs:")
        print("  Summary:", summary_path)
        print("  All age terms:", age_terms_path)
        print("  Problematic clocks:", problematic_path)
        print("  Metadata:", metadata_path)
        print("============================================================")
        print("\nQC status counts:")
        print(
            summary_df["final_year_scale_qc_status"]
            .value_counts(dropna=False)
            .to_string()
        )

    return {
        "summary_df": summary_df,
        "age_terms_df": age_terms_df,
        "problematic_df": problematic_df,
        "summary_path": summary_path,
        "age_terms_path": age_terms_path,
        "problematic_path": problematic_path,
        "metadata_path": metadata_path,
    }


# ============================================================
# 7. Example calls
# ============================================================

if __name__ == "__main__":
    # Example 1: MI L'EPOCH clocks
    run_lepoch_year_scale_qc(
        base_dir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock",
        clock_dir_glob="*_mi_clock",
        analysis_label="mi",
        out_dir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mi_lepoch_year_scale_qc",
    )

    # Example 2: COPD L'EPOCH clocks
    # Uncomment to run:
    run_lepoch_year_scale_qc(
        base_dir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock",
        clock_dir_glob="*_copd_clock",
        analysis_label="copd",
        out_dir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/copd_lepoch_year_scale_qc",
    )

    # Example 3: Dementia L'EPOCH clocks
    # Uncomment to run:
    run_lepoch_year_scale_qc(
        base_dir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock",
        clock_dir_glob="*_dementia_clock",
        analysis_label="dementia",
        out_dir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/dementia_lepoch_year_scale_qc",
    )

    # Example 3: Asthma L'EPOCH clocks
    # Uncomment to run:
    run_lepoch_year_scale_qc(
        base_dir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock",
        clock_dir_glob="*_asthma_clock",
        analysis_label="asthma",
        out_dir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/asthma_lepoch_year_scale_qc",
    )

    # Example 3: Stroke L'EPOCH clocks
    # Uncomment to run:
    run_lepoch_year_scale_qc(
        base_dir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock",
        clock_dir_glob="*_stroke_clock",
        analysis_label="stroke",
        out_dir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/asthma_lepoch_year_scale_qc",
    )

    # Example 4: mortality clocks
    # Uncomment to run:
    run_lepoch_year_scale_qc(
        base_dir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock",
        clock_dir_glob="*_mortality_clock",
        analysis_label="mortality",
        out_dir="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_year_scale_qc",
    )