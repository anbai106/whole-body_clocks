#!/usr/bin/env python3
# ============================================================
# Respiratory medication subtype analysis for immune
# metabolomics mortality-clock delta.
#
# Analyses:
#   1) Top medication codes among respiratory-cluster participants
#   2) Respiratory medication subtype flags:
#        inhaled corticosteroids
#        bronchodilators / beta-agonists
#        anticholinergics
#        leukotriene modifiers
#        systemic corticosteroids
#   3) Linear models:
#        delta_immune_clock ~ respiratory_cluster +
#          baseline_immune_clock +
#          age + sex + BMI + smoking + BP +
#          baseline asthma + baseline COPD
#
# Also fits subtype models:
#        delta_immune_clock ~ subtype + baseline_clock + covariates
#
# Input files expected from prior medication-cluster pipeline:
#   medication_instance0_long_classified.tsv
#   medication_participant_clusters.tsv
#   metabolomics_delta_clock_medication_cluster_requested5_long.tsv
#
# Output:
#   respiratory_top_medication_codes.tsv
#   respiratory_medication_code_subtype_classification.tsv
#   respiratory_participant_subtype_flags.tsv
#   immune_delta_respiratory_analysis_dataset.tsv
#   immune_delta_respiratory_lm_results.tsv
#   immune_delta_respiratory_adjusted_means.tsv
#   immune_delta_respiratory_model_summaries.tsv
# ============================================================

import argparse
import json
import os
import re
import warnings

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")


RESP_CLUSTER = "Respiratory medication cluster"
REF_CLUSTER = "No/minimal medication"

DISEASE_FIELD_ASTHMA = "42014"
DISEASE_FIELD_COPD = "42016"
BASELINE_DATE_FIELD = "53"
UMEL_BASELINE_DATE_RAW = "53-0.0"

DEFAULT_COVARIATES = [
    "baseline_accel_years",
    "chronological_age_0_0",
    "Sex",
    "Smoking",
    "BMI",
    "Diastolic",
    "Systolic",
    "baseline_asthma",
    "baseline_copd",
]


# ============================================================
# 1. Arguments
# ============================================================

def parse_args():
    p = argparse.ArgumentParser(
        description="Respiratory medication subtype analysis for immune delta clock."
    )

    p.add_argument(
        "--input_long_tsv",
        default=(
            "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/"
            "mortality_clock/SA/delta_metabolomics_algorithmic_disease_onset/"
            "medication_cluster_delta_clock_inputs/"
            "metabolomics_delta_clock_medication_cluster_requested5_long.tsv"
        ),
        help="Long-format delta-clock + medication-cluster file."
    )

    p.add_argument(
        "--med_long_tsv",
        default=(
            "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/"
            "mortality_clock/SA/delta_metabolomics_algorithmic_disease_onset/"
            "medication_cluster_delta_clock_inputs/"
            "medication_instance0_long_classified.tsv"
        ),
        help="Long classified baseline medication table from the previous medication-cluster script."
    )

    p.add_argument(
        "--participant_cluster_tsv",
        default=(
            "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/"
            "mortality_clock/SA/delta_metabolomics_algorithmic_disease_onset/"
            "medication_cluster_delta_clock_inputs/"
            "medication_participant_clusters.tsv"
        ),
        help="Participant-level medication-cluster table."
    )

    p.add_argument(
        "--umel_death_xlsx",
        default="/cbica/home/wenju/Dataset/UKBB_UMelbourne/Death_related_var_from_Ye.xlsx",
        help="UMelbourne UKBB Excel file containing baseline date and algorithmic outcome dates."
    )

    p.add_argument(
        "--umel_match_csv",
        default="/cbica/home/wenju/Dataset/UKBB_UMelbourne/UKB_UMelbourne_vs_Penn_match_key.csv",
        help="Mapping key from UMelbourne participant ID to Penn/UPenn participant ID."
    )

    p.add_argument(
        "--out_dir",
        default=(
            "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/"
            "mortality_clock/SA/delta_metabolomics_algorithmic_disease_onset/"
            "respiratory_medication_subtype_immune_delta_results"
        ),
        help="Output directory."
    )

    p.add_argument(
        "--outcome",
        default="delta_clock_age_years",
        choices=["delta_clock_age_years", "delta_accel_years"],
        help="Immune delta outcome to analyze."
    )

    p.add_argument(
        "--robust_cov",
        default="HC3",
        choices=["HC0", "HC1", "HC2", "HC3", "nonrobust"],
        help="Robust covariance estimator for OLS."
    )

    p.add_argument(
        "--min_n_exposure",
        default=20,
        type=int,
        help="Minimum exposed participants for reporting a subtype model."
    )

    return p.parse_args()


# ============================================================
# 2. Helpers
# ============================================================

def normalize_participant_id(df, col):
    if col not in df.columns:
        raise ValueError(f"Missing ID column: {col}")
    out = df.copy()
    out[col] = pd.to_numeric(out[col], errors="coerce")
    out = out[out[col].notna()].copy()
    out[col] = out[col].astype(np.int64)
    return out


def parse_ukb_date(series):
    x = series.copy()
    x = x.replace([0, 0.0, "0", "0.0", "", "NA", "NaN", "nan", "None", "-1", -1], np.nan)

    parsed = pd.to_datetime(x, errors="coerce")

    numeric = pd.to_numeric(x, errors="coerce")
    excel_mask = numeric.between(20000, 60000)

    if excel_mask.any():
        excel_dates = pd.to_datetime(
            numeric,
            unit="D",
            origin="1899-12-30",
            errors="coerce"
        )
        parsed = parsed.where(~excel_mask, excel_dates)

    return parsed


def find_ukb_field_column(df, field_id, required=True):
    field_id = str(field_id)

    exact_candidates = [
        field_id,
        f"{field_id}-0.0",
        f"{field_id}_0_0",
        f"f{field_id}",
        f"f{field_id}_0_0",
    ]

    for wanted in exact_candidates:
        for c in df.columns:
            if str(c) == wanted:
                return c

    patterns = [
        rf"(^|[^0-9]){re.escape(field_id)}([^0-9]|$)",
        rf"f{re.escape(field_id)}(_|$)",
    ]

    matches = []
    for c in df.columns:
        s = str(c)
        for pat in patterns:
            if re.search(pat, s):
                matches.append(c)
                break

    if matches:
        for c in matches:
            s = str(c)
            if s.endswith("-0.0") or s.endswith("_0_0"):
                return c
        return matches[0]

    if required:
        preview = ", ".join(list(map(str, df.columns[:50])))
        raise ValueError(f"Could not find UKB field {field_id}. First columns: {preview}")

    return ""


def safe_numeric(s):
    return pd.to_numeric(s, errors="coerce")


def contains_any(text, patterns):
    if pd.isna(text):
        return False
    text = str(text).lower()
    for pat in patterns:
        if re.search(pat, text, flags=re.IGNORECASE):
            return True
    return False


def fmt_model_formula(outcome, exposure, covars):
    return f"{outcome} ~ {exposure} + " + " + ".join(covars)


# ============================================================
# 3. Respiratory subtype keyword definitions
# ============================================================

INHALED_CORTICOSTEROID_PATTERNS = [
    r"\bbeclometasone\b",
    r"\bbeclomethasone\b",
    r"\bbudesonide\b",
    r"\bfluticasone\b",
    r"\bciclesonide\b",
    r"\bmometasone\b",
    r"\bqvar\b",
    r"\bpulmicort\b",
    r"\bflixotide\b",
    r"\basmanex\b",
    r"\binhaled corticosteroid\b",
    r"\binhaled steroid\b",
    r"\bsteroid inhaler\b",
    r"\bseretide\b",
    r"\bsymbicort\b",
    r"\bfostair\b",
    r"\brelvar\b",
    r"\badvair\b",
]

BETA_AGONIST_PATTERNS = [
    r"\bsalbutamol\b",
    r"\balbuterol\b",
    r"\bventolin\b",
    r"\bterbutaline\b",
    r"\bbricanyl\b",
    r"\bsalmeterol\b",
    r"\bformoterol\b",
    r"\bindacaterol\b",
    r"\bolodaterol\b",
    r"\bvilanterol\b",
    r"\bbeta[- ]?agonist\b",
    r"\bsaba\b",
    r"\blaba\b",
    r"\bbronchodilator\b",
    r"\bseretide\b",
    r"\bsymbicort\b",
    r"\bfostair\b",
    r"\brelvar\b",
    r"\badvair\b",
]

ANTICHOLINERGIC_PATTERNS = [
    r"\bipratropium\b",
    r"\btiotropium\b",
    r"\baclidinium\b",
    r"\bglycopyrronium\b",
    r"\bglycopyrrolate\b",
    r"\bumeclidinium\b",
    r"\batrovent\b",
    r"\bspiriva\b",
    r"\bseebri\b",
    r"\bincruse\b",
    r"\banticholinergic\b",
    r"\bmuscarinic\b",
    r"\blama\b",
    r"\bsama\b",
]

LEUKOTRIENE_PATTERNS = [
    r"\bmontelukast\b",
    r"\bzafirlukast\b",
    r"\bzileuton\b",
    r"\bleukotriene\b",
]

SYSTEMIC_CORTICOSTEROID_PATTERNS = [
    r"\bprednisolone\b",
    r"\bprednisone\b",
    r"\bhydrocortisone\b",
    r"\bdexamethasone\b",
    r"\bmethylprednisolone\b",
    r"\btriamcinolone\b",
    r"\bsystemic corticosteroid\b",
    r"\boral steroid\b",
    r"\bsteroid tablet\b",
]

GENERAL_RESPIRATORY_PATTERNS = [
    r"\binhaler\b",
    r"\basthma\b",
    r"\bcopd\b",
    r"\bchronic obstructive\b",
    r"\brespiratory\b",
    r"\btheophylline\b",
    r"\baminophylline\b",
]


SUBTYPE_PATTERNS = {
    "inhaled_corticosteroid": INHALED_CORTICOSTEROID_PATTERNS,
    "beta_agonist_bronchodilator": BETA_AGONIST_PATTERNS,
    "anticholinergic": ANTICHOLINERGIC_PATTERNS,
    "leukotriene_modifier": LEUKOTRIENE_PATTERNS,
    "systemic_corticosteroid": SYSTEMIC_CORTICOSTEROID_PATTERNS,
    "other_respiratory": GENERAL_RESPIRATORY_PATTERNS,
}


def classify_respiratory_subtypes(meaning):
    out = {}
    for subtype, patterns in SUBTYPE_PATTERNS.items():
        out[subtype] = int(contains_any(meaning, patterns))

    named_subtypes = [
        "inhaled_corticosteroid",
        "beta_agonist_bronchodilator",
        "anticholinergic",
        "leukotriene_modifier",
        "systemic_corticosteroid",
    ]

    out["any_specific_respiratory_subtype"] = int(any(out[x] == 1 for x in named_subtypes))

    # Other respiratory only if it has general respiratory wording but not a specific subtype.
    out["other_respiratory"] = int(
        out["other_respiratory"] == 1 and out["any_specific_respiratory_subtype"] == 0
    )

    out["any_respiratory_subtype"] = int(
        out["any_specific_respiratory_subtype"] == 1 or out["other_respiratory"] == 1
    )

    return out


# ============================================================
# 4. Read inputs
# ============================================================

def read_long_delta(path, outcome):
    if not os.path.exists(path):
        raise FileNotFoundError(path)

    dat = pd.read_csv(path, sep="\t", dtype={"participant_id": str}, low_memory=False)

    # Handle accidental pasted/header typo if present.
    if (
        "baseline_accel_yearsfollowup_accel_years" in dat.columns
        and "baseline_accel_years" not in dat.columns
    ):
        dat = dat.rename(
            columns={
                "baseline_accel_yearsfollowup_accel_years": "baseline_accel_years"
            }
        )

    required = [
        "participant_id",
        "medication_cluster",
        "organ_label",
        "organ_clean",
        outcome,
        "baseline_accel_years",
        "chronological_age_0_0",
        "Sex",
        "Smoking",
        "BMI",
        "Diastolic",
        "Systolic",
        "delta_chrono_age_1_minus_0",
    ]

    missing = [c for c in required if c not in dat.columns]
    if missing:
        raise ValueError(
            "Input long delta file missing required columns:\n"
            + "\n".join(missing)
        )

    dat = normalize_participant_id(dat, "participant_id")

    # Immune clock only.
    dat = dat[
        (dat["organ_clean"].astype(str).str.lower() == "immune")
        | (dat["organ_label"].astype(str).str.lower() == "immune")
    ].copy()

    for c in [outcome] + DEFAULT_COVARIATES:
        if c in dat.columns:
            dat[c] = safe_numeric(dat[c])

    for c in ["Sex", "Smoking"]:
        dat.loc[dat[c].isin([-1, -3]), c] = np.nan

    return dat


def read_medication_long(path):
    if not os.path.exists(path):
        raise FileNotFoundError(path)

    med = pd.read_csv(path, sep="\t", low_memory=False)

    required = [
        "participant_id",
        "participant_id_umel",
        "coding_str",
        "meaning",
    ]

    missing = [c for c in required if c not in med.columns]
    if missing:
        raise ValueError(
            "Medication long file missing required columns:\n"
            + "\n".join(missing)
        )

    med = normalize_participant_id(med, "participant_id")
    med = normalize_participant_id(med, "participant_id_umel")

    med["meaning"] = med["meaning"].astype(str)
    med["coding_str"] = med["coding_str"].astype(str)

    if "is_respiratory" not in med.columns:
        med["is_respiratory"] = med["meaning"].apply(
            lambda x: int(contains_any(x, GENERAL_RESPIRATORY_PATTERNS + BETA_AGONIST_PATTERNS + INHALED_CORTICOSTEROID_PATTERNS + ANTICHOLINERGIC_PATTERNS + LEUKOTRIENE_PATTERNS + SYSTEMIC_CORTICOSTEROID_PATTERNS))
        )
    else:
        med["is_respiratory"] = safe_numeric(med["is_respiratory"]).fillna(0).astype(int)

    subtype_df = med["meaning"].apply(classify_respiratory_subtypes).apply(pd.Series)
    med = pd.concat([med, subtype_df], axis=1)

    return med


def read_participant_clusters(path):
    if not os.path.exists(path):
        raise FileNotFoundError(path)

    cl = pd.read_csv(path, sep="\t", low_memory=False)

    required = [
        "participant_id",
        "participant_id_umel",
        "medication_cluster",
    ]

    missing = [c for c in required if c not in cl.columns]
    if missing:
        raise ValueError(
            "Participant cluster file missing required columns:\n"
            + "\n".join(missing)
        )

    cl = normalize_participant_id(cl, "participant_id")
    cl = normalize_participant_id(cl, "participant_id_umel")

    return cl


# ============================================================
# 5. Baseline asthma and COPD
# ============================================================

def read_id_mapping(path):
    if not os.path.exists(path):
        raise FileNotFoundError(path)

    match = pd.read_csv(path)

    required = ["id", "id_upenn"]
    missing = [c for c in required if c not in match.columns]
    if missing:
        raise ValueError(f"ID match file missing columns: {missing}")

    match = match.rename(
        columns={
            "id": "participant_id_umel",
            "id_upenn": "participant_id",
        }
    )

    match = normalize_participant_id(match, "participant_id_umel")
    match = normalize_participant_id(match, "participant_id")

    return match[["participant_id_umel", "participant_id"]].drop_duplicates()


def read_baseline_asthma_copd(umel_xlsx, match_csv):
    if not os.path.exists(umel_xlsx):
        raise FileNotFoundError(umel_xlsx)

    raw = pd.read_excel(umel_xlsx, engine="openpyxl")

    if "eid" not in raw.columns:
        raise ValueError("UMelbourne Excel must contain eid.")

    baseline_col = find_ukb_field_column(raw, BASELINE_DATE_FIELD, required=True)
    asthma_col = find_ukb_field_column(raw, DISEASE_FIELD_ASTHMA, required=True)
    copd_col = find_ukb_field_column(raw, DISEASE_FIELD_COPD, required=True)

    raw = raw.rename(columns={"eid": "participant_id_umel"})
    raw = normalize_participant_id(raw, "participant_id_umel")

    match = read_id_mapping(match_csv)

    d = match.merge(
        raw[["participant_id_umel", baseline_col, asthma_col, copd_col]],
        on="participant_id_umel",
        how="inner"
    )

    d["baseline_date_instance0"] = parse_ukb_date(d[baseline_col])
    d["asthma_date"] = parse_ukb_date(d[asthma_col])
    d["copd_date"] = parse_ukb_date(d[copd_col])

    d["baseline_asthma"] = np.where(
        d["asthma_date"].notna() & d["baseline_date_instance0"].notna(),
        (d["asthma_date"] <= d["baseline_date_instance0"]).astype(int),
        np.where(d["asthma_date"].isna(), 0, np.nan)
    )

    d["baseline_copd"] = np.where(
        d["copd_date"].notna() & d["baseline_date_instance0"].notna(),
        (d["copd_date"] <= d["baseline_date_instance0"]).astype(int),
        np.where(d["copd_date"].isna(), 0, np.nan)
    )

    out = d[
        [
            "participant_id",
            "participant_id_umel",
            "baseline_date_instance0",
            "asthma_date",
            "copd_date",
            "baseline_asthma",
            "baseline_copd",
        ]
    ].drop_duplicates("participant_id", keep="first")

    return out


# ============================================================
# 6. Respiratory medication code summaries
# ============================================================

def make_top_respiratory_codes(med, clusters):
    respiratory_ids = clusters.loc[
        clusters["medication_cluster"] == RESP_CLUSTER,
        ["participant_id"]
    ].drop_duplicates()

    med_resp_cluster = med.merge(respiratory_ids, on="participant_id", how="inner")

    # Use all codes among respiratory-cluster people, then identify respiratory-specific codes.
    top_all = (
        med_resp_cluster
        .groupby(["coding_str", "meaning"], dropna=False)
        .agg(
            n_participants=("participant_id", "nunique"),
            is_respiratory=("is_respiratory", "max"),
            inhaled_corticosteroid=("inhaled_corticosteroid", "max"),
            beta_agonist_bronchodilator=("beta_agonist_bronchodilator", "max"),
            anticholinergic=("anticholinergic", "max"),
            leukotriene_modifier=("leukotriene_modifier", "max"),
            systemic_corticosteroid=("systemic_corticosteroid", "max"),
            other_respiratory=("other_respiratory", "max"),
            any_respiratory_subtype=("any_respiratory_subtype", "max"),
        )
        .reset_index()
        .sort_values(["any_respiratory_subtype", "n_participants"], ascending=[False, False])
    )

    top_resp_only = top_all[top_all["any_respiratory_subtype"] == 1].copy()

    return top_all, top_resp_only


def make_participant_subtype_flags(med, clusters):
    respiratory_ids = clusters[
        ["participant_id", "participant_id_umel", "medication_cluster"]
    ].drop_duplicates()

    # Keep medications for everyone because subtype flags can be helpful even outside
    # the dominant respiratory cluster.
    subtypes = [
        "inhaled_corticosteroid",
        "beta_agonist_bronchodilator",
        "anticholinergic",
        "leukotriene_modifier",
        "systemic_corticosteroid",
        "other_respiratory",
        "any_specific_respiratory_subtype",
        "any_respiratory_subtype",
    ]

    flags = (
        med.groupby("participant_id", dropna=False)[subtypes]
        .max()
        .reset_index()
    )

    flags = respiratory_ids.merge(flags, on="participant_id", how="left")

    for c in subtypes:
        flags[c] = safe_numeric(flags[c]).fillna(0).astype(int)

    flags["is_respiratory_cluster"] = (
        flags["medication_cluster"] == RESP_CLUSTER
    ).astype(int)

    flags["is_no_minimal_medication"] = (
        flags["medication_cluster"] == REF_CLUSTER
    ).astype(int)

    return flags


# ============================================================
# 7. Linear modeling
# ============================================================

def available_covariates(df, covars):
    out = []
    for c in covars:
        if c not in df.columns:
            continue
        x = safe_numeric(df[c])
        if x.notna().sum() > 0 and x.nunique(dropna=True) > 1:
            out.append(c)
    return out


def fit_lm(df, outcome, exposure, covars, robust_cov):
    needed = [outcome, exposure] + covars

    dat = df[needed].copy()

    for c in needed:
        dat[c] = safe_numeric(dat[c])

    dat = dat.replace([np.inf, -np.inf], np.nan).dropna().copy()

    if dat.shape[0] < 50:
        raise ValueError(f"Insufficient complete cases: N={dat.shape[0]}")

    if dat[exposure].nunique(dropna=True) < 2:
        raise ValueError(f"Exposure has <2 levels: {exposure}")

    formula = fmt_model_formula(outcome, exposure, covars)

    model = smf.ols(formula=formula, data=dat)

    if robust_cov == "nonrobust":
        fit = model.fit()
    else:
        fit = model.fit(cov_type=robust_cov)

    return fit, dat, formula


def extract_binary_exposure_result(fit, fit_df, outcome, exposure, formula, model_name):
    term = exposure

    if term not in fit.params.index:
        # statsmodels can sometimes rename binary variables, but for numeric 0/1
        # it should remain exactly exposure.
        candidates = [x for x in fit.params.index if exposure in x]
        if len(candidates) == 1:
            term = candidates[0]
        else:
            raise ValueError(f"Could not find exposure term in model: {exposure}")

    beta = float(fit.params.loc[term])
    se = float(fit.bse.loc[term])
    p = float(fit.pvalues.loc[term])

    out = {
        "model_name": model_name,
        "outcome": outcome,
        "exposure": exposure,
        "term": term,
        "formula": formula,
        "N": int(fit_df.shape[0]),
        "N_exposed": int((fit_df[exposure] == 1).sum()),
        "N_unexposed": int((fit_df[exposure] == 0).sum()),
        "mean_outcome_exposed": float(fit_df.loc[fit_df[exposure] == 1, outcome].mean()),
        "mean_outcome_unexposed": float(fit_df.loc[fit_df[exposure] == 0, outcome].mean()),
        "beta": beta,
        "se": se,
        "ci_lo": beta - 1.96 * se,
        "ci_hi": beta + 1.96 * se,
        "p": p,
        "r_squared": float(fit.rsquared),
        "adj_r_squared": float(fit.rsquared_adj),
        "aic": float(fit.aic),
        "bic": float(fit.bic),
        "status": "ok",
        "error": "",
    }

    return out


def get_adjusted_means_binary(fit, fit_df, outcome, exposure, covars, model_name):
    mean_covars = {}
    for c in covars:
        mean_covars[c] = float(safe_numeric(fit_df[c]).mean())

    rows = []

    for val, label in [(0, "unexposed_or_reference"), (1, "exposed")]:
        row = {exposure: val}
        row.update(mean_covars)

        new_df = pd.DataFrame([row])

        try:
            pred = fit.get_prediction(new_df).summary_frame(alpha=0.05)
            rows.append({
                "model_name": model_name,
                "outcome": outcome,
                "exposure": exposure,
                "exposure_value": val,
                "exposure_label": label,
                "adjusted_mean": float(pred["mean"].iloc[0]),
                "adjusted_mean_se": float(pred["mean_se"].iloc[0]),
                "adjusted_ci_lo": float(pred["mean_ci_lower"].iloc[0]),
                "adjusted_ci_hi": float(pred["mean_ci_upper"].iloc[0]),
                "covariate_values": json.dumps(mean_covars),
            })
        except Exception as e:
            rows.append({
                "model_name": model_name,
                "outcome": outcome,
                "exposure": exposure,
                "exposure_value": val,
                "exposure_label": label,
                "adjusted_mean": np.nan,
                "adjusted_mean_se": np.nan,
                "adjusted_ci_lo": np.nan,
                "adjusted_ci_hi": np.nan,
                "prediction_error": str(e),
                "covariate_values": json.dumps(mean_covars),
            })

    return rows


def run_models(analysis_df, outcome, robust_cov, min_n_exposure):
    results = []
    adjusted_means = []
    summaries = []

    # Primary model:
    # respiratory cluster versus no/minimal medication.
    primary = analysis_df[
        analysis_df["medication_cluster"].isin([REF_CLUSTER, RESP_CLUSTER])
    ].copy()

    primary["respiratory_cluster_binary"] = (
        primary["medication_cluster"] == RESP_CLUSTER
    ).astype(int)

    covars = available_covariates(primary, DEFAULT_COVARIATES)

    try:
        fit, fit_df, formula = fit_lm(
            primary,
            outcome=outcome,
            exposure="respiratory_cluster_binary",
            covars=covars,
            robust_cov=robust_cov
        )

        res = extract_binary_exposure_result(
            fit=fit,
            fit_df=fit_df,
            outcome=outcome,
            exposure="respiratory_cluster_binary",
            formula=formula,
            model_name="primary_respiratory_cluster_vs_no_minimal"
        )
        results.append(res)

        adjusted_means.extend(
            get_adjusted_means_binary(
                fit=fit,
                fit_df=fit_df,
                outcome=outcome,
                exposure="respiratory_cluster_binary",
                covars=covars,
                model_name="primary_respiratory_cluster_vs_no_minimal"
            )
        )

        summaries.append({
            "model_name": "primary_respiratory_cluster_vs_no_minimal",
            "status": "ok",
            "N": int(fit_df.shape[0]),
            "formula": formula,
            "covariates_used": ",".join(covars),
            "r_squared": float(fit.rsquared),
            "adj_r_squared": float(fit.rsquared_adj),
        })

    except Exception as e:
        results.append({
            "model_name": "primary_respiratory_cluster_vs_no_minimal",
            "outcome": outcome,
            "exposure": "respiratory_cluster_binary",
            "status": "failed",
            "error": str(e),
        })
        summaries.append({
            "model_name": "primary_respiratory_cluster_vs_no_minimal",
            "status": "failed",
            "error": str(e),
        })

    # Subtype models:
    # Use no/minimal medication as reference plus participants in respiratory cluster.
    subtype_base = analysis_df[
        analysis_df["medication_cluster"].isin([REF_CLUSTER, RESP_CLUSTER])
    ].copy()

    subtype_exposures = [
        "inhaled_corticosteroid",
        "beta_agonist_bronchodilator",
        "anticholinergic",
        "leukotriene_modifier",
        "systemic_corticosteroid",
        "other_respiratory",
    ]

    for exposure in subtype_exposures:
        dat = subtype_base.copy()

        # Reference group should be no/minimal meds. Participants with respiratory
        # cluster but without this subtype are kept as unexposed in the model.
        dat[exposure] = safe_numeric(dat[exposure]).fillna(0).astype(int)

        n_exp = int((dat[exposure] == 1).sum())

        model_name = f"subtype_{exposure}_within_no_minimal_plus_respiratory_cluster"

        if n_exp < min_n_exposure:
            results.append({
                "model_name": model_name,
                "outcome": outcome,
                "exposure": exposure,
                "status": "small_exposed_n",
                "N_exposed": n_exp,
                "min_n_exposure": min_n_exposure,
                "error": f"N_exposed={n_exp} < {min_n_exposure}",
            })
            continue

        covars = available_covariates(dat, DEFAULT_COVARIATES)

        try:
            fit, fit_df, formula = fit_lm(
                dat,
                outcome=outcome,
                exposure=exposure,
                covars=covars,
                robust_cov=robust_cov
            )

            res = extract_binary_exposure_result(
                fit=fit,
                fit_df=fit_df,
                outcome=outcome,
                exposure=exposure,
                formula=formula,
                model_name=model_name
            )
            results.append(res)

            adjusted_means.extend(
                get_adjusted_means_binary(
                    fit=fit,
                    fit_df=fit_df,
                    outcome=outcome,
                    exposure=exposure,
                    covars=covars,
                    model_name=model_name
                )
            )

            summaries.append({
                "model_name": model_name,
                "status": "ok",
                "N": int(fit_df.shape[0]),
                "N_exposed": int((fit_df[exposure] == 1).sum()),
                "formula": formula,
                "covariates_used": ",".join(covars),
                "r_squared": float(fit.rsquared),
                "adj_r_squared": float(fit.rsquared_adj),
            })

        except Exception as e:
            results.append({
                "model_name": model_name,
                "outcome": outcome,
                "exposure": exposure,
                "status": "failed",
                "error": str(e),
            })
            summaries.append({
                "model_name": model_name,
                "status": "failed",
                "error": str(e),
            })

    # Joint subtype model:
    # Include all subtype flags together to evaluate independent subtype effects.
    joint_exposures = [
        "inhaled_corticosteroid",
        "beta_agonist_bronchodilator",
        "anticholinergic",
        "leukotriene_modifier",
        "systemic_corticosteroid",
    ]

    joint = subtype_base.copy()
    for exposure in joint_exposures:
        joint[exposure] = safe_numeric(joint[exposure]).fillna(0).astype(int)

    usable_joint_exposures = [
        x for x in joint_exposures
        if int((joint[x] == 1).sum()) >= min_n_exposure
        and joint[x].nunique(dropna=True) > 1
    ]

    if len(usable_joint_exposures) >= 1:
        covars = available_covariates(joint, DEFAULT_COVARIATES)
        formula = f"{outcome} ~ " + " + ".join(usable_joint_exposures + covars)

        needed = [outcome] + usable_joint_exposures + covars
        fit_df = joint[needed].copy()

        for c in needed:
            fit_df[c] = safe_numeric(fit_df[c])

        fit_df = fit_df.replace([np.inf, -np.inf], np.nan).dropna().copy()

        try:
            model = smf.ols(formula=formula, data=fit_df)
            if robust_cov == "nonrobust":
                fit = model.fit()
            else:
                fit = model.fit(cov_type=robust_cov)

            for exposure in usable_joint_exposures:
                beta = float(fit.params.loc[exposure])
                se = float(fit.bse.loc[exposure])
                p = float(fit.pvalues.loc[exposure])

                results.append({
                    "model_name": "joint_respiratory_subtype_model",
                    "outcome": outcome,
                    "exposure": exposure,
                    "term": exposure,
                    "formula": formula,
                    "N": int(fit_df.shape[0]),
                    "N_exposed": int((fit_df[exposure] == 1).sum()),
                    "N_unexposed": int((fit_df[exposure] == 0).sum()),
                    "mean_outcome_exposed": float(fit_df.loc[fit_df[exposure] == 1, outcome].mean()),
                    "mean_outcome_unexposed": float(fit_df.loc[fit_df[exposure] == 0, outcome].mean()),
                    "beta": beta,
                    "se": se,
                    "ci_lo": beta - 1.96 * se,
                    "ci_hi": beta + 1.96 * se,
                    "p": p,
                    "r_squared": float(fit.rsquared),
                    "adj_r_squared": float(fit.rsquared_adj),
                    "aic": float(fit.aic),
                    "bic": float(fit.bic),
                    "status": "ok",
                    "error": "",
                })

            summaries.append({
                "model_name": "joint_respiratory_subtype_model",
                "status": "ok",
                "N": int(fit_df.shape[0]),
                "formula": formula,
                "covariates_used": ",".join(covars),
                "subtypes_used": ",".join(usable_joint_exposures),
                "r_squared": float(fit.rsquared),
                "adj_r_squared": float(fit.rsquared_adj),
            })

        except Exception as e:
            summaries.append({
                "model_name": "joint_respiratory_subtype_model",
                "status": "failed",
                "error": str(e),
                "formula": formula,
            })

    return pd.DataFrame(results), pd.DataFrame(adjusted_means), pd.DataFrame(summaries)


def add_fdr(results):
    if results.empty or "p" not in results.columns:
        return results

    results = results.copy()
    results["p_fdr_bh"] = np.nan
    results["p_bonferroni"] = np.nan
    results["fdr_significant_0.05"] = False
    results["bonferroni_significant_0.05"] = False

    mask = results["p"].notna() & (results["status"] == "ok")
    if mask.sum() > 0:
        pvals = results.loc[mask, "p"].astype(float).values
        _, p_fdr, _, _ = multipletests(pvals, alpha=0.05, method="fdr_bh")
        p_bonf = np.minimum(pvals * len(pvals), 1.0)

        results.loc[mask, "p_fdr_bh"] = p_fdr
        results.loc[mask, "p_bonferroni"] = p_bonf
        results.loc[mask, "fdr_significant_0.05"] = p_fdr < 0.05
        results.loc[mask, "bonferroni_significant_0.05"] = p_bonf < 0.05

    return results


# ============================================================
# 8. Main
# ============================================================

def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print("============================================================")
    print("Respiratory medication subtype analysis for immune delta clock")
    print("============================================================")
    print("Input long delta TSV:", args.input_long_tsv)
    print("Medication long TSV:", args.med_long_tsv)
    print("Participant cluster TSV:", args.participant_cluster_tsv)
    print("UMelbourne disease/date Excel:", args.umel_death_xlsx)
    print("Output directory:", args.out_dir)
    print("Outcome:", args.outcome)
    print("============================================================")

    immune = read_long_delta(args.input_long_tsv, args.outcome)
    med = read_medication_long(args.med_long_tsv)
    clusters = read_participant_clusters(args.participant_cluster_tsv)
    baseline_disease = read_baseline_asthma_copd(
        args.umel_death_xlsx,
        args.umel_match_csv
    )

    # Top medication codes in respiratory cluster.
    top_all, top_resp_only = make_top_respiratory_codes(med, clusters)

    top_all_out = os.path.join(args.out_dir, "respiratory_cluster_top_all_medication_codes.tsv")
    top_resp_out = os.path.join(args.out_dir, "respiratory_cluster_top_respiratory_medication_codes.tsv")

    top_all.to_csv(top_all_out, sep="\t", index=False)
    top_resp_only.to_csv(top_resp_out, sep="\t", index=False)

    # Code subtype classification.
    code_subtype = (
        med.groupby(["coding_str", "meaning"], dropna=False)
        .agg(
            n_participants=("participant_id", "nunique"),
            is_respiratory=("is_respiratory", "max"),
            inhaled_corticosteroid=("inhaled_corticosteroid", "max"),
            beta_agonist_bronchodilator=("beta_agonist_bronchodilator", "max"),
            anticholinergic=("anticholinergic", "max"),
            leukotriene_modifier=("leukotriene_modifier", "max"),
            systemic_corticosteroid=("systemic_corticosteroid", "max"),
            other_respiratory=("other_respiratory", "max"),
            any_respiratory_subtype=("any_respiratory_subtype", "max"),
        )
        .reset_index()
        .sort_values(["any_respiratory_subtype", "n_participants"], ascending=[False, False])
    )

    code_subtype_out = os.path.join(args.out_dir, "respiratory_medication_code_subtype_classification.tsv")
    code_subtype.to_csv(code_subtype_out, sep="\t", index=False)

    # Participant subtype flags.
    subtype_flags = make_participant_subtype_flags(med, clusters)
    subtype_flags_out = os.path.join(args.out_dir, "respiratory_participant_subtype_flags.tsv")
    subtype_flags.to_csv(subtype_flags_out, sep="\t", index=False)

    # Build analysis dataset.
    analysis = immune.merge(
        subtype_flags.drop(columns=["participant_id_umel"], errors="ignore"),
        on="participant_id",
        how="left",
        suffixes=("", "_from_flags")
    )

    # Prefer medication_cluster from immune input if present; fill from subtype table otherwise.
    if "medication_cluster_from_flags" in analysis.columns:
        analysis["medication_cluster"] = analysis["medication_cluster"].where(
            analysis["medication_cluster"].notna(),
            analysis["medication_cluster_from_flags"]
        )
        analysis = analysis.drop(columns=["medication_cluster_from_flags"])

    analysis = analysis.merge(
        baseline_disease[
            [
                "participant_id",
                "baseline_date_instance0",
                "asthma_date",
                "copd_date",
                "baseline_asthma",
                "baseline_copd",
            ]
        ],
        on="participant_id",
        how="left"
    )

    # Fill subtype flags for no-medication people.
    subtype_cols = [
        "inhaled_corticosteroid",
        "beta_agonist_bronchodilator",
        "anticholinergic",
        "leukotriene_modifier",
        "systemic_corticosteroid",
        "other_respiratory",
        "any_specific_respiratory_subtype",
        "any_respiratory_subtype",
        "is_respiratory_cluster",
        "is_no_minimal_medication",
    ]

    for c in subtype_cols:
        if c not in analysis.columns:
            analysis[c] = 0
        analysis[c] = safe_numeric(analysis[c]).fillna(0).astype(int)

    # Baseline asthma/COPD should be 0 if no date; NA only if ambiguous due to date issues.
    analysis["baseline_asthma"] = safe_numeric(analysis["baseline_asthma"])
    analysis["baseline_copd"] = safe_numeric(analysis["baseline_copd"])

    analysis_out = os.path.join(args.out_dir, "immune_delta_respiratory_analysis_dataset.tsv")
    analysis.to_csv(analysis_out, sep="\t", index=False)

    # Subtype participant counts.
    subtype_count_rows = []
    for c in subtype_cols:
        subtype_count_rows.append({
            "variable": c,
            "N": int(analysis.shape[0]),
            "N_positive": int((analysis[c] == 1).sum()),
            "N_negative": int((analysis[c] == 0).sum()),
            "pct_positive": float((analysis[c] == 1).mean() * 100.0),
        })

    subtype_counts = pd.DataFrame(subtype_count_rows)
    subtype_counts_out = os.path.join(args.out_dir, "respiratory_subtype_counts_in_immune_delta_dataset.tsv")
    subtype_counts.to_csv(subtype_counts_out, sep="\t", index=False)

    # Run models.
    results, adjusted_means, summaries = run_models(
        analysis_df=analysis,
        outcome=args.outcome,
        robust_cov=args.robust_cov,
        min_n_exposure=args.min_n_exposure
    )

    results = add_fdr(results)

    results_out = os.path.join(args.out_dir, "immune_delta_respiratory_lm_results.tsv")
    means_out = os.path.join(args.out_dir, "immune_delta_respiratory_adjusted_means.tsv")
    summaries_out = os.path.join(args.out_dir, "immune_delta_respiratory_model_summaries.tsv")

    results.to_csv(results_out, sep="\t", index=False)
    adjusted_means.to_csv(means_out, sep="\t", index=False)
    summaries.to_csv(summaries_out, sep="\t", index=False)

    metadata = {
        "input_long_tsv": args.input_long_tsv,
        "med_long_tsv": args.med_long_tsv,
        "participant_cluster_tsv": args.participant_cluster_tsv,
        "umel_death_xlsx": args.umel_death_xlsx,
        "umel_match_csv": args.umel_match_csv,
        "out_dir": args.out_dir,
        "outcome": args.outcome,
        "robust_cov": args.robust_cov,
        "min_n_exposure": args.min_n_exposure,
        "outputs": {
            "top_all_medication_codes": top_all_out,
            "top_respiratory_medication_codes": top_resp_out,
            "code_subtype_classification": code_subtype_out,
            "participant_subtype_flags": subtype_flags_out,
            "analysis_dataset": analysis_out,
            "subtype_counts": subtype_counts_out,
            "lm_results": results_out,
            "adjusted_means": means_out,
            "model_summaries": summaries_out,
        },
        "n_immune_delta_rows": int(immune.shape[0]),
        "n_medication_long_rows": int(med.shape[0]),
        "n_analysis_rows": int(analysis.shape[0]),
    }

    metadata_out = os.path.join(args.out_dir, "run_metadata.json")
    with open(metadata_out, "w") as f:
        json.dump(metadata, f, indent=2, default=str)

    print("\n============================================================")
    print("Finished respiratory subtype immune-delta analyses.")
    print("Key outputs:")
    print("  Top respiratory codes:", top_resp_out)
    print("  Participant subtype flags:", subtype_flags_out)
    print("  Analysis dataset:", analysis_out)
    print("  LM results:", results_out)
    print("  Model summaries:", summaries_out)
    print("  Metadata:", metadata_out)
    print("============================================================\n")

    print("Top respiratory medication codes:")
    print(top_resp_only.head(25).to_string(index=False))

    print("\nLinear model results:")
    if not results.empty:
        cols = [
            "model_name",
            "exposure",
            "N",
            "N_exposed",
            "beta",
            "ci_lo",
            "ci_hi",
            "p",
            "p_fdr_bh",
            "status",
        ]
        show_cols = [c for c in cols if c in results.columns]
        print(results[show_cols].to_string(index=False))


if __name__ == "__main__":
    main()