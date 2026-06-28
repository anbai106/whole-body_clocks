#!/usr/bin/env python3
"""
Compare mortality clocks versus corresponding chronological-age BAGs for incident disease onset.

For each ICD endpoint, this script runs one fair pairwise comparison per organ/modality:
  - mortality clock alone: covariates + mortality_clock_z
  - BAG alone:             covariates + BAG_z
  - joint model:           covariates + mortality_clock_z + BAG_z

Fairness rule:
  For each mortality-clock/BAG pair, both predictors are evaluated on exactly the same participants
  after applying the same time-origin, incident-case filtering, covariate missingness filtering,
  and non-missingness filtering for both scores.
"""

import argparse
import os
import warnings
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from scipy.stats import chi2, norm, wilcoxon

warnings.filterwarnings("ignore")


CLOCK_BAG_PAIRS: List[Dict[str, str]] = [
    # MRI clocks, measured at imaging visit
    {"pair_id": "Brain_mri", "organ": "Brain", "modality": "MRI", "time_origin": "imaging", "clock": "Brain_mri", "bag": "Brain_MRIBAG"},
    {"pair_id": "Adipose_mri", "organ": "Adipose", "modality": "MRI", "time_origin": "imaging", "clock": "Adipose_mri", "bag": "Adipose_MRIBAG"},
    {"pair_id": "Heart_mri", "organ": "Heart", "modality": "MRI", "time_origin": "imaging", "clock": "Heart_mri", "bag": "Heart_MRIBAG"},
    {"pair_id": "Kidney_mri", "organ": "Kidney", "modality": "MRI", "time_origin": "imaging", "clock": "Kidney_mri", "bag": "Kidney_MRIBAG"},
    {"pair_id": "Liver_mri", "organ": "Liver", "modality": "MRI", "time_origin": "imaging", "clock": "Liver_mri", "bag": "Liver_MRIBAG"},
    {"pair_id": "Pancreas_mri", "organ": "Pancreas", "modality": "MRI", "time_origin": "imaging", "clock": "Pancreas_mri", "bag": "Pancreas_MRIBAG"},
    {"pair_id": "Spleen_mri", "organ": "Spleen", "modality": "MRI", "time_origin": "imaging", "clock": "Spleen_mri", "bag": "Spleen_MRIBAG"},

    # Proteomic clocks, assumed measured at baseline/blood assessment
    {"pair_id": "Reproductive_female_proteomics", "organ": "Reproductive_female", "modality": "Proteomics", "time_origin": "baseline", "clock": "Reproductive_female_proteomics", "bag": "Reproductive_female_ProtBAG"},
    {"pair_id": "Pulmonary_proteomics", "organ": "Pulmonary", "modality": "Proteomics", "time_origin": "baseline", "clock": "Pulmonary_proteomics", "bag": "Pulmonary_ProtBAG"},
    {"pair_id": "Heart_proteomics", "organ": "Heart", "modality": "Proteomics", "time_origin": "baseline", "clock": "Heart_proteomics", "bag": "Heart_ProtBAG"},
    {"pair_id": "Brain_proteomics", "organ": "Brain", "modality": "Proteomics", "time_origin": "baseline", "clock": "Brain_proteomics", "bag": "Brain_ProtBAG"},
    {"pair_id": "Eye_proteomics", "organ": "Eye", "modality": "Proteomics", "time_origin": "baseline", "clock": "Eye_proteomics", "bag": "Eye_ProtBAG"},
    {"pair_id": "Hepatic_proteomics", "organ": "Hepatic", "modality": "Proteomics", "time_origin": "baseline", "clock": "Hepatic_proteomics", "bag": "Hepatic_ProtBAG"},
    {"pair_id": "Renal_proteomics", "organ": "Renal", "modality": "Proteomics", "time_origin": "baseline", "clock": "Renal_proteomics", "bag": "Renal_ProtBAG"},
    {"pair_id": "Reproductive_male_proteomics", "organ": "Reproductive_male", "modality": "Proteomics", "time_origin": "baseline", "clock": "Reproductive_male_proteomics", "bag": "Reproductive_male_ProtBAG"},
    {"pair_id": "Endocrine_proteomics", "organ": "Endocrine", "modality": "Proteomics", "time_origin": "baseline", "clock": "Endocrine_proteomics", "bag": "Endocrine_ProtBAG"},
    {"pair_id": "Immune_proteomics", "organ": "Immune", "modality": "Proteomics", "time_origin": "baseline", "clock": "Immune_proteomics", "bag": "Immune_ProtBAG"},
    {"pair_id": "Skin_proteomics", "organ": "Skin", "modality": "Proteomics", "time_origin": "baseline", "clock": "Skin_proteomics", "bag": "Skin_ProtBAG"},

    # Metabolomic clocks, assumed measured at baseline/blood assessment
    {"pair_id": "Endocrine_metabolomics", "organ": "Endocrine", "modality": "Metabolomics", "time_origin": "baseline", "clock": "Endocrine_metabolomics", "bag": "Endocrine_MetBAG"},
    {"pair_id": "Digestive_metabolomics", "organ": "Digestive", "modality": "Metabolomics", "time_origin": "baseline", "clock": "Digestive_metabolomics", "bag": "Digestive_MetBAG"},
    {"pair_id": "Hepatic_metabolomics", "organ": "Hepatic", "modality": "Metabolomics", "time_origin": "baseline", "clock": "Hepatic_metabolomics", "bag": "Hepatic_MetBAG"},
    {"pair_id": "Immune_metabolomics", "organ": "Immune", "modality": "Metabolomics", "time_origin": "baseline", "clock": "Immune_metabolomics", "bag": "Immune_MetBAG"},
]

MORTALITY_CLOCK_COLS = [p["clock"] for p in CLOCK_BAG_PAIRS]
BAG_COLS = [p["bag"] for p in CLOCK_BAG_PAIRS]

BASELINE_DATE_COL = "date_of_attending_assessment_centre_f53_0_0"
IMAGING_DATE_COL = "date_of_attending_assessment_centre_f53_2_0"
AGE_RECRUIT_COL = "age_at_recruitment_f21022_0_0"
SMOKING_COL = "smoking_status_f20116_0_0"
BMI_COL = "body_mass_index_bmi_f23104_0_0"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Fair Cox comparison of mortality clocks versus BAGs for incident disease onset."
    )
    parser.add_argument("--icd_tsv", required=True, type=str)
    parser.add_argument("--output_tsv", required=True, type=str)
    parser.add_argument(
        "--bag_tsv",
        default="/cbica/home/wenju/Reproducibile_paper/SleepAging/data/MomoBAG.tsv",
        type=str,
    )
    parser.add_argument(
        "--cov_tsv",
        default="/cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv",
        type=str,
    )
    parser.add_argument(
        "--date_tsv",
        default="/cbica/home/wenju/Reproducibile_paper/Multiorgan_Subtype/data/PWAS/UKBB_fullsample_death_variables.csv",
        type=str,
    )
    parser.add_argument("--min_case", default=20, type=int)
    parser.add_argument("--min_noncase", default=20, type=int)
    parser.add_argument("--penalizer", default=0.0, type=float)
    return parser.parse_args()


def disease_id_from_path(path: str) -> str:
    base = os.path.basename(path)
    if base.endswith("_diagnosis_clock.tsv"):
        return base.replace("_diagnosis_clock.tsv", "")
    return os.path.splitext(base)[0]


def clean_event_dates(series: pd.Series) -> pd.Series:
    """Treat non-case placeholders such as 0 as missing before date parsing."""
    x = series.copy()
    x = x.replace([0, 0.0, "0", "0.0", "", "NA", "NaN", "nan", "None", "-1", -1], np.nan)
    return pd.to_datetime(x, errors="coerce")


def read_covariates(path: str) -> pd.DataFrame:
    cov_all = pd.read_csv(path)
    if "eid" not in cov_all.columns:
        raise ValueError(f"Cannot find eid in covariate file: {path}")

    sex_candidates = ["sex_f31_0_0", "genetic_sex_f22001_0_0", "Sex", "sex"]
    sex_col = next((c for c in sex_candidates if c in cov_all.columns), None)

    keep = ["eid"]
    for col in [AGE_RECRUIT_COL, SMOKING_COL, BMI_COL, sex_col]:
        if col is not None and col in cov_all.columns:
            keep.append(col)

    cov = cov_all[keep].copy()
    rename = {
        "eid": "participant_id",
        AGE_RECRUIT_COL: "Age_baseline",
        SMOKING_COL: "Smoking",
        BMI_COL: "BMI",
    }
    if sex_col is not None:
        rename[sex_col] = "Sex"
    cov = cov.rename(columns=rename)

    for c in cov.columns:
        if c != "participant_id":
            cov[c] = pd.to_numeric(cov[c], errors="coerce")
    return cov


def read_assessment_dates(path: str) -> pd.DataFrame:
    d = pd.read_csv(path)
    if "eid" in d.columns:
        d = d.rename(columns={"eid": "participant_id"})
    keep = ["participant_id", BASELINE_DATE_COL, IMAGING_DATE_COL]
    missing = [c for c in keep if c not in d.columns]
    if missing:
        raise ValueError(f"Missing date columns in {path}: {missing}")
    d = d[keep].copy()
    d[BASELINE_DATE_COL] = pd.to_datetime(d[BASELINE_DATE_COL], errors="coerce")
    d[IMAGING_DATE_COL] = pd.to_datetime(d[IMAGING_DATE_COL], errors="coerce")
    return d


def construct_survival_data(args) -> Tuple[pd.DataFrame, str]:
    disease_id = disease_id_from_path(args.icd_tsv)

    df_clock_all = pd.read_csv(args.icd_tsv, sep="\t")
    required_clock_cols = ["participant_id", "case", "date"] + MORTALITY_CLOCK_COLS
    missing_clock = [c for c in required_clock_cols if c not in df_clock_all.columns]
    if missing_clock:
        raise ValueError(f"Missing columns in ICD clock TSV {args.icd_tsv}: {missing_clock}")
    df_clock = df_clock_all[required_clock_cols].copy()

    df_bag_all = pd.read_csv(args.bag_tsv, sep="\t")

    # MomoBAG stores the brain MRI BAG as Brain_PhenoBAG. Create the expected
    # analysis alias Brain_MRIBAG so the downstream clock/BAG pair list can stay
    # consistent with the other MRI BAG column names.
    if "Brain_MRIBAG" not in df_bag_all.columns and "Brain_PhenoBAG" in df_bag_all.columns:
        df_bag_all = df_bag_all.rename(columns={"Brain_PhenoBAG": "Brain_MRIBAG"})

    required_bag_cols = ["participant_id"] + BAG_COLS
    missing_bag = [c for c in required_bag_cols if c not in df_bag_all.columns]
    if missing_bag:
        raise ValueError(f"Missing columns in BAG TSV {args.bag_tsv}: {missing_bag}")
    df_bag = df_bag_all[required_bag_cols].copy()

    cov = read_covariates(args.cov_tsv)
    dates = read_assessment_dates(args.date_tsv)

    data = df_clock.merge(df_bag, on="participant_id", how="left")
    data = data.merge(cov, on="participant_id", how="left")
    data = data.merge(dates, on="participant_id", how="left")

    data["case"] = pd.to_numeric(data["case"], errors="coerce").fillna(0).astype(int)
    data["case"] = (data["case"] == 1).astype(int)

    data["event_date"] = clean_event_dates(data["date"])
    data.loc[data["case"] == 0, "event_date"] = pd.NaT

    if data.loc[data["case"] == 1, "event_date"].notna().sum() == 0:
        raise ValueError("No usable event dates among cases after date cleaning.")

    # Disease-specific administrative censor date. If you have a true disease-specific censor date,
    # replace this with that external censor date.
    # global_end_date = data.loc[data["event_date"].notna(), "event_date"].max() + pd.Timedelta(days=2)
    global_end_date = pd.Timestamp("2022-11-30")

    data["time_baseline"] = np.where(
        data["case"] == 1,
        (data["event_date"] - data[BASELINE_DATE_COL]).dt.days,
        (global_end_date - data[BASELINE_DATE_COL]).dt.days,
    )

    data["time_imaging"] = np.where(
        data["case"] == 1,
        (data["event_date"] - data[IMAGING_DATE_COL]).dt.days,
        (global_end_date - data[IMAGING_DATE_COL]).dt.days,
    )

    # Use age at the actual time origin. For MRI, the time origin is imaging visit.
    data["Age_imaging"] = data["Age_baseline"] + (
        (data[IMAGING_DATE_COL] - data[BASELINE_DATE_COL]).dt.days / 365.25
    )

    data["disease_id"] = disease_id
    return data, disease_id


def standardize_inplace(df: pd.DataFrame, src: str, dst: str) -> bool:
    vals = pd.to_numeric(df[src], errors="coerce")
    sd = vals.std()
    if not np.isfinite(sd) or sd == 0:
        return False
    df[dst] = (vals - vals.mean()) / sd
    return True


def base_covariates(df: pd.DataFrame, age_col: str) -> List[str]:
    candidates = [age_col, "Sex", "Smoking", "BMI"]
    covars = []
    for c in candidates:
        if c in df.columns:
            vals = pd.to_numeric(df[c], errors="coerce")
            if vals.notna().sum() > 0 and vals.nunique(dropna=True) > 1:
                df[c] = vals
                covars.append(c)
    return covars


def fit_cox(df: pd.DataFrame, time_col: str, covars: List[str], predictors: List[str], penalizer: float = 0.0):
    cols = [time_col, "case"] + covars + predictors
    fit_df = df[cols].copy()
    for c in cols:
        fit_df[c] = pd.to_numeric(fit_df[c], errors="coerce")
    fit_df = fit_df.replace([np.inf, -np.inf], np.nan).dropna().copy()
    fit_df = fit_df[fit_df[time_col] > 0].copy()

    # Drop zero-variance covariates; never drop requested predictors here.
    usable_covars = []
    for c in covars:
        if fit_df[c].nunique(dropna=True) > 1:
            usable_covars.append(c)
    model_cols = [time_col, "case"] + usable_covars + predictors
    fit_df = fit_df[model_cols].copy()

    if fit_df["case"].sum() == 0 or (fit_df["case"] == 0).sum() == 0:
        raise ValueError("No cases or no non-cases in fit_df.")

    # Retry with small penalization if the unpenalized model fails.
    last_error = None
    for pen in [penalizer, 0.001, 0.01, 0.1]:
        try:
            cph = CoxPHFitter(penalizer=pen)
            cph.fit(fit_df, duration_col=time_col, event_col="case", show_progress=False)
            return cph, fit_df, usable_covars, pen
        except Exception as e:
            last_error = e
    raise last_error


def extract_hr(cph: CoxPHFitter, var: str) -> Dict[str, float]:
    beta = float(cph.params_.loc[var])
    se = float(cph.standard_errors_.loc[var])
    hr = float(np.exp(beta))
    ci_lo = float(np.exp(beta - 1.96 * se))
    ci_hi = float(np.exp(beta + 1.96 * se))
    p = float(cph.summary.loc[var, "p"])
    return {"beta": beta, "se": se, "hr": hr, "ci_lo": ci_lo, "ci_hi": ci_hi, "p": p}


def get_cindex(cph: CoxPHFitter, fit_df: pd.DataFrame, time_col: str) -> float:
    risk = cph.predict_partial_hazard(fit_df).values.ravel()
    # concordance_index assumes larger score means longer survival; risk means shorter survival.
    return float(concordance_index(fit_df[time_col], -risk, fit_df["case"]))


def lrt_pvalue(cph_full: CoxPHFitter, cph_reduced: CoxPHFitter, df_diff: int) -> float:
    stat = 2.0 * (float(cph_full.log_likelihood_) - float(cph_reduced.log_likelihood_))
    if not np.isfinite(stat) or stat < 0:
        return np.nan
    return float(chi2.sf(stat, df_diff))


def empty_result(disease_id: str, pair: Dict[str, str], status: str, error: str = "") -> Dict[str, object]:
    out = {
        "disease_id": disease_id,
        "pair_id": pair["pair_id"],
        "organ": pair["organ"],
        "modality": pair["modality"],
        "time_origin": pair["time_origin"],
        "mortality_clock": pair["clock"],
        "bag": pair["bag"],
        "status": status,
        "error": error,
    }
    numeric_cols = [
        "N", "N_case", "N_noncase", "event_rate",
        "followup_years_min", "followup_years_max",
        "event_followup_years_min", "event_followup_years_max",
        "clock_bag_pearson",
        "clock_beta", "clock_se", "clock_hr", "clock_ci_lo", "clock_ci_hi", "clock_p",
        "bag_beta", "bag_se", "bag_hr", "bag_ci_lo", "bag_ci_hi", "bag_p",
        "clock_joint_beta", "clock_joint_se", "clock_joint_hr", "clock_joint_p",
        "bag_joint_beta", "bag_joint_se", "bag_joint_hr", "bag_joint_p",
        "joint_beta_diff_clock_minus_bag", "joint_se_diff", "joint_z_diff", "joint_p_diff",
        "base_cindex", "clock_cindex", "bag_cindex", "both_cindex",
        "delta_cindex_clock_minus_bag", "delta_cindex_clock_minus_base",
        "delta_cindex_bag_minus_base", "delta_cindex_both_minus_base",
        "lrt_p_clock_vs_base", "lrt_p_bag_vs_base", "lrt_p_both_vs_base",
        "penalizer_base", "penalizer_clock", "penalizer_bag", "penalizer_both",
        "disease_level_n_success_pairs", "disease_level_mean_delta_abs_beta_clock_minus_bag",
        "disease_level_wilcoxon_abs_beta_p", "disease_level_mean_delta_cindex_clock_minus_bag",
        "disease_level_wilcoxon_cindex_p",
    ]
    for c in numeric_cols:
        out[c] = np.nan
    return out


def analyze_pair(data: pd.DataFrame, disease_id: str, pair: Dict[str, str], args) -> Dict[str, object]:
    clock = pair["clock"]
    bag = pair["bag"]
    time_col = "time_imaging" if pair["time_origin"] == "imaging" else "time_baseline"
    age_col = "Age_imaging" if pair["time_origin"] == "imaging" else "Age_baseline"

    if clock not in data.columns or bag not in data.columns:
        return empty_result(disease_id, pair, "missing_columns", f"Missing {clock} or {bag}")

    covars_initial = base_covariates(data, age_col)
    needed = [clock, bag, time_col, "case"] + covars_initial
    df = data[needed].copy()
    for c in needed:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.replace([np.inf, -np.inf], np.nan).dropna().copy()

    # Incident disease filtering relative to the measurement time origin.
    df = df[df[time_col] > 0].copy()

    # Follow-up duration from the pair-specific time origin to event or censoring.
    # For MRI clocks this is imaging-to-event/censoring; for proteomics/metabolomics
    # clocks this is baseline-to-event/censoring. Event-specific summaries are
    # reported among incident cases only.
    followup_years = pd.to_numeric(df[time_col], errors="coerce") / 365.25
    event_followup_years = followup_years[df["case"] == 1]
    followup_summary = {
        "followup_years_min": float(followup_years.min()) if followup_years.notna().any() else np.nan,
        "followup_years_max": float(followup_years.max()) if followup_years.notna().any() else np.nan,
        "event_followup_years_min": float(event_followup_years.min()) if event_followup_years.notna().any() else np.nan,
        "event_followup_years_max": float(event_followup_years.max()) if event_followup_years.notna().any() else np.nan,
    }

    n_case = int(df["case"].sum())
    n_noncase = int((df["case"] == 0).sum())
    if n_case < args.min_case or n_noncase < args.min_noncase:
        out = empty_result(disease_id, pair, "insufficient_events")
        out.update({"N": len(df), "N_case": n_case, "N_noncase": n_noncase})
        out.update(followup_summary)
        return out

    if not standardize_inplace(df, clock, "clock_z"):
        return empty_result(disease_id, pair, "zero_variance_clock")
    if not standardize_inplace(df, bag, "bag_z"):
        return empty_result(disease_id, pair, "zero_variance_bag")

    # Rebuild covariates after pairwise complete-case filtering.
    covars = base_covariates(df, age_col)

    try:
        cph_base, df_base, covars_base, pen_base = fit_cox(df, time_col, covars, [], args.penalizer)
        cph_clock, df_clock, _, pen_clock = fit_cox(df, time_col, covars_base, ["clock_z"], args.penalizer)
        cph_bag, df_bag, _, pen_bag = fit_cox(df, time_col, covars_base, ["bag_z"], args.penalizer)
        cph_both, df_both, _, pen_both = fit_cox(df, time_col, covars_base, ["clock_z", "bag_z"], args.penalizer)
    except Exception as e:
        out = empty_result(disease_id, pair, "cox_fit_failed", str(e))
        out.update({"N": len(df), "N_case": n_case, "N_noncase": n_noncase})
        out.update(followup_summary)
        return out

    clock_stats = extract_hr(cph_clock, "clock_z")
    bag_stats = extract_hr(cph_bag, "bag_z")
    clock_joint = extract_hr(cph_both, "clock_z")
    bag_joint = extract_hr(cph_both, "bag_z")

    # Wald test in the joint model: H0 beta_clock = beta_BAG.
    try:
        var_mat = cph_both.variance_matrix_
        var_diff = (
            float(var_mat.loc["clock_z", "clock_z"])
            + float(var_mat.loc["bag_z", "bag_z"])
            - 2.0 * float(var_mat.loc["clock_z", "bag_z"])
        )
        se_diff = np.sqrt(var_diff) if var_diff > 0 else np.nan
        beta_diff = clock_joint["beta"] - bag_joint["beta"]
        z_diff = beta_diff / se_diff if np.isfinite(se_diff) and se_diff > 0 else np.nan
        p_diff = 2.0 * norm.sf(abs(z_diff)) if np.isfinite(z_diff) else np.nan
    except Exception:
        beta_diff, se_diff, z_diff, p_diff = np.nan, np.nan, np.nan, np.nan

    base_cindex = get_cindex(cph_base, df_base, time_col)
    clock_cindex = get_cindex(cph_clock, df_clock, time_col)
    bag_cindex = get_cindex(cph_bag, df_bag, time_col)
    both_cindex = get_cindex(cph_both, df_both, time_col)

    out = empty_result(disease_id, pair, "ok")
    out.update(
        {
            "N": int(len(df)),
            "N_case": n_case,
            "N_noncase": n_noncase,
            "event_rate": float(n_case / len(df)) if len(df) > 0 else np.nan,
            **followup_summary,
            "clock_bag_pearson": float(df[["clock_z", "bag_z"]].corr().iloc[0, 1]),
            "clock_beta": clock_stats["beta"],
            "clock_se": clock_stats["se"],
            "clock_hr": clock_stats["hr"],
            "clock_ci_lo": clock_stats["ci_lo"],
            "clock_ci_hi": clock_stats["ci_hi"],
            "clock_p": clock_stats["p"],
            "bag_beta": bag_stats["beta"],
            "bag_se": bag_stats["se"],
            "bag_hr": bag_stats["hr"],
            "bag_ci_lo": bag_stats["ci_lo"],
            "bag_ci_hi": bag_stats["ci_hi"],
            "bag_p": bag_stats["p"],
            "clock_joint_beta": clock_joint["beta"],
            "clock_joint_se": clock_joint["se"],
            "clock_joint_hr": clock_joint["hr"],
            "clock_joint_p": clock_joint["p"],
            "bag_joint_beta": bag_joint["beta"],
            "bag_joint_se": bag_joint["se"],
            "bag_joint_hr": bag_joint["hr"],
            "bag_joint_p": bag_joint["p"],
            "joint_beta_diff_clock_minus_bag": beta_diff,
            "joint_se_diff": se_diff,
            "joint_z_diff": z_diff,
            "joint_p_diff": p_diff,
            "base_cindex": base_cindex,
            "clock_cindex": clock_cindex,
            "bag_cindex": bag_cindex,
            "both_cindex": both_cindex,
            "delta_cindex_clock_minus_bag": clock_cindex - bag_cindex,
            "delta_cindex_clock_minus_base": clock_cindex - base_cindex,
            "delta_cindex_bag_minus_base": bag_cindex - base_cindex,
            "delta_cindex_both_minus_base": both_cindex - base_cindex,
            "lrt_p_clock_vs_base": lrt_pvalue(cph_clock, cph_base, 1),
            "lrt_p_bag_vs_base": lrt_pvalue(cph_bag, cph_base, 1),
            "lrt_p_both_vs_base": lrt_pvalue(cph_both, cph_base, 2),
            "penalizer_base": pen_base,
            "penalizer_clock": pen_clock,
            "penalizer_bag": pen_bag,
            "penalizer_both": pen_both,
        }
    )
    return out


def add_disease_level_tests(res: pd.DataFrame) -> pd.DataFrame:
    ok = res[res["status"] == "ok"].copy()
    n_ok = len(ok)
    mean_delta_abs_beta = np.nan
    p_abs_beta = np.nan
    mean_delta_cindex = np.nan
    p_cindex = np.nan

    if n_ok >= 5:
        delta_abs_beta = ok["clock_beta"].abs().values - ok["bag_beta"].abs().values
        mean_delta_abs_beta = float(np.nanmean(delta_abs_beta))
        try:
            p_abs_beta = float(wilcoxon(delta_abs_beta, zero_method="wilcox", alternative="two-sided").pvalue)
        except Exception:
            p_abs_beta = np.nan

        delta_cindex = ok["delta_cindex_clock_minus_bag"].values
        mean_delta_cindex = float(np.nanmean(delta_cindex))
        try:
            p_cindex = float(wilcoxon(delta_cindex, zero_method="wilcox", alternative="two-sided").pvalue)
        except Exception:
            p_cindex = np.nan

    res["disease_level_n_success_pairs"] = n_ok
    res["disease_level_mean_delta_abs_beta_clock_minus_bag"] = mean_delta_abs_beta
    res["disease_level_wilcoxon_abs_beta_p"] = p_abs_beta
    res["disease_level_mean_delta_cindex_clock_minus_bag"] = mean_delta_cindex
    res["disease_level_wilcoxon_cindex_p"] = p_cindex
    return res


def main():
    args = parse_args()
    data, disease_id = construct_survival_data(args)

    rows = []
    for pair in CLOCK_BAG_PAIRS:
        rows.append(analyze_pair(data, disease_id, pair, args))

    res = pd.DataFrame(rows)
    res = add_disease_level_tests(res)

    os.makedirs(os.path.dirname(args.output_tsv), exist_ok=True)
    res.to_csv(args.output_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
