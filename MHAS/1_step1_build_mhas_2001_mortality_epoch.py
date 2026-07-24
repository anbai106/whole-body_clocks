#!/usr/bin/env python3
"""
STEP 1: Build 2001-baseline MHAS analytic cohort for phenotype-based mortality EPOCH

Input folder expected:
  /Users/hao/Dropbox/MHAS/download/
    GatewayHarmonizedMHAS/
      H_MHAS_d.dta
      Gateway_Harmonized_MHAS_D_2001-2022.pdf
      H_MHAS_long.do
    GatewayHarmonizedMHASEndofLife/
      GH_MHAS_EOL_c.dta
      Gateway Harmonized MHAS End of Life C 2003-2022.pdf
      GH_MHAS_EOL_long.do
    GatewayHarmonizedMexCog/
      GH_MEX_COG_b2.dta
      ...

Primary design:
  - Baseline: MHAS Wave 1 / 2001
  - Inclusion: Wave-1 alive respondents (inw1 == 1 and r1iwstat == 1)
  - Outcome: all-cause mortality after baseline
  - Mortality timing: radyear/radmonth from Gateway Harmonized MHAS
  - Censoring: latest wave where respondent was alive or known alive non-response
  - Primary predictors: non-disease phenotype variables only

Important:
  Doctor-diagnosed disease labels are intentionally excluded from the primary
  mortality EPOCH feature matrix and saved separately for downstream analyses.

Run:
  cd /Users/hao/Dropbox/MHAS
  python 1_step1_build_mhas_2001_mortality_epoch.py

Outputs:
  /Users/hao/Dropbox/MHAS/step1_2001_mortality_epoch/
    mhas_2001_step1_clean_analytic_cohort.tsv
    mhas_2001_step1_model_input_primary_nondisease.tsv
    mhas_2001_step1_primary_feature_manifest.tsv
    mhas_2001_step1_downstream_disease_labels.tsv
    mhas_2001_step1_audit_summary.txt
"""

import argparse
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from pandas.io.stata import StataReader


# -----------------------------
# Utility functions
# -----------------------------
def log(msg: str) -> None:
    print(msg, flush=True)


def require_file(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing {label}: {path}")


def get_stata_columns_and_labels(path: Path):
    """Read Stata column names and variable labels without loading the whole file."""
    one = next(pd.read_stata(path, chunksize=1, convert_categoricals=False))
    cols = list(one.columns)
    reader = StataReader(path)
    labels = reader.variable_labels()
    try:
        reader.close()
    except Exception:
        pass
    return cols, labels


def dta_shape(path: Path):
    one = next(pd.read_stata(path, chunksize=1, convert_categoricals=False))
    first_col = one.columns[0]
    n = pd.read_stata(path, columns=[first_col], convert_categoricals=False).shape[0]
    return n, len(one.columns)


def make_date(year, month=None, default_month: int = 7, default_day: int = 15):
    """
    Construct dates from year/month. MHAS harmonized mortality often has year/month
    but not day. We assign day=15 to represent mid-month.
    """
    y = pd.to_numeric(year, errors="coerce")
    if month is None:
        m = pd.Series(default_month, index=y.index)
    else:
        m = pd.to_numeric(month, errors="coerce").fillna(default_month)

    m = m.where((m >= 1) & (m <= 12), default_month).astype("Int64")

    dates = []
    for yy, mm in zip(y, m):
        if pd.isna(yy):
            dates.append(pd.NaT)
        else:
            try:
                dates.append(pd.Timestamp(int(yy), int(mm), default_day))
            except Exception:
                dates.append(pd.NaT)
    return pd.Series(dates, index=y.index)


def winsorize_series(s, low: float = 0.005, high: float = 0.995):
    x = pd.to_numeric(s, errors="coerce")
    if x.notna().sum() < 10:
        return x
    lo, hi = x.quantile([low, high])
    return x.clip(lo, hi)


def stata_missing_to_nan(df: pd.DataFrame) -> pd.DataFrame:
    """
    pandas.read_stata(convert_categoricals=False) generally returns Stata numeric
    missings as NaN, but this wrapper is kept for clarity and future safety.
    """
    return df.replace({pd.NA: np.nan})


# -----------------------------
# Primary non-disease feature set
# -----------------------------
FEATURE_SPECS = [
    # IDs / design variables
    ("unhhidnp", "participant_id_raw", "id", "id", "Unique participant ID"),
    ("rahhidnp", "participant_id_char_raw", "id", "id", "Character participant ID"),
    ("r1wtresp", "wave1_person_weight", "design", "continuous", "Wave 1 respondent-level analysis weight"),

    # Demographics and SES
    ("r1agey", "age_2001", "demographic", "continuous", "Age at Wave 1 interview"),
    ("ragender", "sex", "demographic", "categorical", "Gender"),
    ("rabyear", "birth_year", "demographic", "continuous", "Birth year"),
    ("rabmonth", "birth_month", "demographic", "continuous", "Birth month"),
    ("raedyrs", "education_years", "socioeconomic", "continuous", "Years of education"),
    ("raeducel", "education_early_level", "socioeconomic", "categorical", "Harmonized early education"),
    ("raeducl", "education_level", "socioeconomic", "categorical", "Harmonized education"),
    ("raindlang", "speaks_indigenous_language", "socioeconomic", "categorical", "Speaks indigenous language"),
    ("r1mstat", "marital_status_2001", "socioeconomic", "categorical", "Marital status"),
    ("r1mrct", "number_marriages_2001", "socioeconomic", "continuous", "Number of marriages"),
    ("h1rural", "rural_urban_2001", "socioeconomic", "categorical", "Rural/urban residence"),
    ("h1rural_m", "locality_size_2001", "socioeconomic", "categorical", "Size of locality"),
    ("r1work", "currently_working_2001", "socioeconomic", "categorical", "Currently working"),
    ("r1lbrf_m", "labor_force_status_2001", "socioeconomic", "categorical", "Labor force status"),

    # General health, excluding disease diagnosis labels
    ("r1shlt", "self_rated_health_2001", "general_health", "categorical", "Self-rated health"),
    ("r1hltc", "health_change_2001", "general_health", "categorical", "Self-rated health change"),
    ("r1painfr", "frequent_pain_2001", "general_health", "categorical", "Frequent pain"),
    ("r1painlv", "pain_level_2001", "general_health", "categorical", "Pain level"),
    ("r1paina", "pain_interferes_2001", "general_health", "categorical", "Pain interferes with activities"),
    ("r1osleep", "trouble_sleeping_2001", "general_health", "categorical", "Trouble sleeping in past week"),

    # Lifestyle
    ("r1vigact", "vigorous_physical_activity_2001", "lifestyle", "categorical", "Vigorous physical activity 3+ times/week"),
    ("r1smokev", "ever_smoked_2001", "lifestyle", "categorical", "Ever smoked"),
    ("r1smoken", "current_smoking_2001", "lifestyle", "categorical", "Currently smokes"),
    ("r1smokef", "cigarettes_per_day_2001", "lifestyle", "continuous", "Cigarettes per day"),
    ("r1smokefm", "cigarettes_per_day_when_smoking_most_2001", "lifestyle", "continuous", "Cigarettes per day when smoking most"),
    ("r1strtsmok", "age_started_smoking_2001", "lifestyle", "continuous", "Age started smoking"),
    ("r1quitsmok", "age_quit_smoking_2001", "lifestyle", "continuous", "Age quit smoking"),
    ("r1drink", "drinks_alcohol_2001", "lifestyle", "categorical", "Drinks alcohol"),
    ("r1drinkd", "drinking_days_per_week_2001", "lifestyle", "continuous", "Drinking days per week"),
    ("r1drinkn", "drinks_per_day_2001", "lifestyle", "continuous", "Drinks per day"),
    ("r1drinkb", "ever_binge_drinks_2001", "lifestyle", "categorical", "Ever binge drinks"),
    ("r1binged", "binge_drinking_days_2001", "lifestyle", "continuous", "Binge drinking days"),

    # Mental health
    ("r1cesd_m", "cesd_modified_score_2001", "mental_health", "continuous", "Modified CES-D score"),
    ("r1cesdm_m", "cesd_missing_count_2001", "mental_health", "continuous", "CES-D missing count"),
    ("r1depres", "cesd_depressed_2001", "mental_health", "categorical", "CES-D felt depressed"),
    ("r1effort", "cesd_effort_2001", "mental_health", "categorical", "CES-D everything an effort"),
    ("r1sleepr", "cesd_restless_sleep_2001", "mental_health", "categorical", "CES-D restless sleep"),
    ("r1whappy", "cesd_happy_2001", "mental_health", "categorical", "CES-D felt happy"),
    ("r1flone", "cesd_lonely_2001", "mental_health", "categorical", "CES-D felt lonely"),
    ("r1enlife", "cesd_enjoyed_life_2001", "mental_health", "categorical", "CES-D enjoyed life"),
    ("r1fsad", "cesd_sad_2001", "mental_health", "categorical", "CES-D felt sad"),
    ("r1ftired", "cesd_tired_2001", "mental_health", "categorical", "CES-D felt tired"),
    ("r1energ", "cesd_energy_2001", "mental_health", "categorical", "CES-D had energy"),

    # Function
    ("r1adltot6", "adl_0_6_2001", "function", "continuous", "ADL difficulty count 0-6"),
    ("r1adltot6a", "any_adl_0_6_2001", "function", "categorical", "Any ADL difficulty 0-6"),
    ("r1iadlfour", "iadl_0_4_2001", "function", "continuous", "IADL difficulty count 0-4"),
    ("r1iadlfoura", "any_iadl_0_4_2001", "function", "categorical", "Any IADL difficulty 0-4"),
    ("r1mobila", "mobility_0_5_2001", "function", "continuous", "Mobility difficulty count 0-5"),
    ("r1mobilaa", "any_mobility_0_5_2001", "function", "categorical", "Any mobility difficulty 0-5"),
    ("r1lgmusa", "large_muscle_0_4_2001", "function", "continuous", "Large muscle difficulty count 0-4"),
    ("r1grossa", "gross_motor_0_5_2001", "function", "continuous", "Gross motor difficulty count 0-5"),
    ("r1finea", "fine_motor_0_3_2001", "function", "continuous", "Fine motor difficulty count 0-3"),
    ("r1mobilsev", "mobility_0_7_2001", "function", "continuous", "Seven-item mobility count 0-7"),
    ("r1uppermob", "upper_body_mobility_0_3_2001", "function", "continuous", "Upper-body mobility count 0-3"),
    ("r1lowermob", "lower_body_mobility_0_4_2001", "function", "continuous", "Lower-body mobility count 0-4"),
    ("r1walkra", "difficulty_walk_across_room_2001", "function", "categorical", "Difficulty walking across room"),
    ("r1walksa", "difficulty_walk_several_blocks_2001", "function", "categorical", "Difficulty walking several blocks"),
    ("r1walk1a", "difficulty_walk_one_block_2001", "function", "categorical", "Difficulty walking one block"),
    ("r1joga", "difficulty_jog_1km_2001", "function", "categorical", "Difficulty jogging 1 km"),
    ("r1dressa", "difficulty_dressing_2001", "function", "categorical", "Difficulty dressing"),
    ("r1batha", "difficulty_bathing_2001", "function", "categorical", "Difficulty bathing"),
    ("r1eata", "difficulty_eating_2001", "function", "categorical", "Difficulty eating"),
    ("r1beda", "difficulty_getting_in_out_bed_2001", "function", "categorical", "Difficulty getting in/out of bed"),
    ("r1toilta", "difficulty_using_toilet_2001", "function", "categorical", "Difficulty using toilet"),
    ("r1pusha", "difficulty_push_large_object_2001", "function", "categorical", "Difficulty pushing/pulling large objects"),

    # Cognition
    ("r1imrc8", "immediate_word_recall_2001", "cognition", "continuous", "Immediate word recall, first trial"),
    ("r1imrc8_m", "immediate_word_recall_avg_2001", "cognition", "continuous", "Immediate word recall, average three trials"),
    ("r1dlrc8", "delayed_word_recall_2001", "cognition", "continuous", "Delayed word recall"),
    ("r1tr16", "word_recall_summary_2001", "cognition", "continuous", "Word recall summary"),
    ("r1prmem", "proxy_memory_rating_2001", "cognition", "categorical", "Proxy memory rating"),
    ("r1prchmem", "proxy_memory_change_2001", "cognition", "categorical", "Proxy memory change"),

    # Anthropometrics and physical measures
    ("r1bmi", "bmi_self_report_2001", "anthropometric", "continuous", "Self-reported BMI"),
    ("r1bmicat", "bmi_category_self_report_2001", "anthropometric", "categorical", "Self-reported BMI category"),
    ("r1height", "height_self_report_m_2001", "anthropometric", "continuous", "Self-reported height in meters"),
    ("r1weight", "weight_self_report_kg_2001", "anthropometric", "continuous", "Self-reported weight in kg"),
    ("r1mheight", "height_measured_m_2001", "anthropometric", "continuous", "Measured height in meters"),
    ("r1mweight", "weight_measured_kg_2001", "anthropometric", "continuous", "Measured weight in kg"),
    ("r1mbmi", "bmi_measured_2001", "anthropometric", "continuous", "Measured BMI"),
    ("r1mbmicat", "bmi_category_measured_2001", "anthropometric", "categorical", "Measured BMI category"),
    ("r1mwaist", "waist_measured_cm_2001", "anthropometric", "continuous", "Measured waist circumference"),
    ("r1mwhratio", "waist_hip_ratio_measured_2001", "anthropometric", "continuous", "Measured waist-hip ratio"),
    ("r1kneehght", "knee_height_cm_2001", "anthropometric", "continuous", "Measured knee height"),

    # Health care utilization and insurance, not diagnoses
    ("r1hosp1y", "hospital_stay_1y_2001", "healthcare", "categorical", "Hospital stay in prior 12 months"),
    ("r1hspnit1y", "hospital_nights_1y_2001", "healthcare", "continuous", "Hospital nights in prior 12 months"),
    ("r1doctor1y", "doctor_visit_1y_2001", "healthcare", "categorical", "Doctor visit in prior 12 months"),
    ("r1doctim1y", "doctor_visits_count_1y_2001", "healthcare", "continuous", "Number of doctor visits in prior 12 months"),
    ("r1hipriv", "private_health_insurance_2001", "healthcare", "categorical", "Private health insurance"),
    ("r1htnum", "health_insurance_plan_count_2001", "healthcare", "continuous", "Number of health insurance plans"),

    # Early-life features
    ("raserchlth", "childhood_serious_health_problem", "early_life", "categorical", "Serious health problem before age 10"),
    ("rachchlt", "childhood_health_compared", "early_life", "categorical", "Health compared with other children before age 10"),
    ("rachhdinj", "childhood_head_injury", "early_life", "categorical", "Serious head injury before age 10"),
    ("rameduc_m", "mother_education", "early_life", "categorical", "Mother's education"),
    ("rafeduc_m", "father_education", "early_life", "categorical", "Father's education"),
]


# Disease labels saved for downstream prediction, not used as primary EPOCH features.
DISEASE_MAP = {
    "hypertension": ["r1hibpe", "r2hibpe", "r3hibpe", "r4hibpe", "r5hibpe", "r6hibpe"],
    "diabetes": ["r1diabe", "r2diabe", "r3diabe", "r4diabe", "r5diabe", "r6diabe"],
    "cancer": ["r1cancre", "r2cancre", "r3cancre", "r4cancre", "r5cancre", "r6cancre"],
    "respiratory_disease_incl_asthma": ["r1respe", "r2respe", "r3respe", "r4respe", "r5respe", "r6respe"],
    "heart_attack": ["r1hrtatte", "r2hrtatte", "r3hrtatte", "r4hrtatte", "r5hrtatte", "r6hrtatte"],
    "heart_problem": ["r4hearte", "r5hearte", "r6hearte"],  # available later in harmonized file
    "stroke": ["r1stroke", "r2stroke", "r3stroke", "r4stroke", "r5stroke", "r6stroke"],
    "arthritis": ["r1arthre", "r2arthre", "r3arthre", "r4arthre", "r5arthre", "r6arthre"],
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--base-dir",
        default="/Users/hao/Dropbox/MHAS",
        help="Base MHAS directory."
    )
    parser.add_argument(
        "--out-dir",
        default=None,
        help="Output directory. Default: <base-dir>/step1_2001_mortality_epoch"
    )
    parser.add_argument(
        "--min-feature-nonmissing",
        type=float,
        default=0.10,
        help="Minimum nonmissing fraction for a source feature to enter primary model input."
    )
    parser.add_argument(
        "--restrict-age50",
        action="store_true",
        help="If set, restrict analytic cohort to participants age >=50 at Wave 1."
    )
    args = parser.parse_args()

    base_dir = Path(args.base_dir)
    download_dir = base_dir / "download"
    out_dir = Path(args.out_dir) if args.out_dir else base_dir / "step1_2001_mortality_epoch"
    out_dir.mkdir(parents=True, exist_ok=True)

    h_file = download_dir / "GatewayHarmonizedMHAS" / "H_MHAS_d.dta"
    eol_file = download_dir / "GatewayHarmonizedMHASEndofLife" / "GH_MHAS_EOL_c.dta"
    mexcog_file = download_dir / "GatewayHarmonizedMexCog" / "GH_MEX_COG_b2.dta"

    require_file(h_file, "Gateway Harmonized MHAS file")
    require_file(eol_file, "Gateway Harmonized MHAS End of Life file")

    log(f"Reading metadata from: {h_file}")
    h_cols, h_labels = get_stata_columns_and_labels(h_file)
    eol_cols, _ = get_stata_columns_and_labels(eol_file)

    # Record uploaded/local file shapes.
    shape_rows = []
    for label, path in [
        ("H_MHAS_d", h_file),
        ("GH_MHAS_EOL_c", eol_file),
        ("GH_MEX_COG_b2", mexcog_file),
    ]:
        if path.exists() and path.suffix.lower() == ".dta":
            try:
                n, k = dta_shape(path)
                shape_rows.append({"name": label, "path": str(path), "n_rows": n, "n_columns": k})
            except Exception as e:
                shape_rows.append({"name": label, "path": str(path), "n_rows": np.nan, "n_columns": np.nan, "error": str(e)})

    # Select needed columns from main MHAS.
    base_cols = [
        "unhhidnp", "rahhidnp",
        "inw1", "inw2", "inw3", "inw4", "inw5", "inw6",
        "r1iwstat", "r2iwstat", "r3iwstat", "r4iwstat", "r5iwstat", "r6iwstat",
        "r1iwy", "r1iwm", "r2iwy", "r2iwm", "r3iwy", "r3iwm",
        "r4iwy", "r4iwm", "r5iwy", "r5iwm", "r6iwy", "r6iwm",
        "radyear", "radmonth",
    ]

    all_cols_needed = set([c for c in base_cols if c in h_cols])
    for src, _, _, _, _ in FEATURE_SPECS:
        if src in h_cols:
            all_cols_needed.add(src)
    for cols in DISEASE_MAP.values():
        for c in cols:
            if c in h_cols:
                all_cols_needed.add(c)

    log(f"Reading selected columns from H_MHAS_d.dta: {len(all_cols_needed)} columns")
    h = pd.read_stata(h_file, columns=sorted(all_cols_needed), convert_categoricals=False).copy()
    h = stata_missing_to_nan(h)

    # Construct baseline and mortality/censoring dates.
    h["baseline_date"] = make_date(h["r1iwy"], h["r1iwm"], default_month=7, default_day=15)
    h["death_date"] = make_date(h["radyear"], h["radmonth"], default_month=7, default_day=15)
    h["death_month_missing_flag"] = h["radyear"].notna() & h["radmonth"].isna()

    # Last known alive date: status 1 = respondent alive; status 4 = nonresponse but known/presumed alive.
    last_alive = pd.Series(pd.NaT, index=h.index, dtype="datetime64[ns]")
    last_alive_wave = pd.Series(np.nan, index=h.index)
    wave_info = [
        (1, "r1iwstat", "r1iwy", "r1iwm", 2001),
        (2, "r2iwstat", "r2iwy", "r2iwm", 2003),
        (3, "r3iwstat", "r3iwy", "r3iwm", 2012),
        (4, "r4iwstat", "r4iwy", "r4iwm", 2015),
        (5, "r5iwstat", "r5iwy", "r5iwm", 2018),
        (6, "r6iwstat", "r6iwy", "r6iwm", 2021),
    ]

    for wave, stat, year, month, nominal_year in wave_info:
        status = pd.to_numeric(h[stat], errors="coerce")
        alive_mask = status.isin([1, 4])

        y = pd.to_numeric(h[year], errors="coerce").fillna(nominal_year)
        m = pd.to_numeric(h[month], errors="coerce").fillna(7)
        wave_date = make_date(y, m, default_month=7, default_day=15)

        update = alive_mask & wave_date.notna() & (last_alive.isna() | (wave_date > last_alive))
        last_alive.loc[update] = wave_date.loc[update]
        last_alive_wave.loc[update] = wave

    h["last_alive_date"] = last_alive
    h["last_alive_wave"] = last_alive_wave

    h["event_death"] = ((h["death_date"].notna()) & (h["death_date"] > h["baseline_date"])).astype(int)
    h["followup_end_date"] = h["last_alive_date"]
    h.loc[h["event_death"] == 1, "followup_end_date"] = h.loc[h["event_death"] == 1, "death_date"]
    h["followup_years"] = (h["followup_end_date"] - h["baseline_date"]).dt.days / 365.25

    # Primary analytic cohort.
    cohort = h[(h["inw1"] == 1) & (h["r1iwstat"] == 1)].copy()
    cohort["age50plus_2001"] = (pd.to_numeric(cohort["r1agey"], errors="coerce") >= 50).astype(int)

    cohort = cohort[
        cohort["baseline_date"].notna()
        & cohort["followup_end_date"].notna()
        & (cohort["followup_years"] > 0)
    ].copy()

    if args.restrict_age50:
        cohort = cohort[cohort["age50plus_2001"] == 1].copy()

    # Cross-check mortality with EOL file.
    eol_read_cols = [
        c for c in [
            "unhhidnp", "rahhidnp", "raxyear", "raxmonth",
            "raxtiwy", "raxtiwm", "ralstcore", "ralstcorey", "radage"
        ]
        if c in eol_cols
    ]
    eol = pd.read_stata(eol_file, columns=eol_read_cols, convert_categoricals=False)
    eol = eol.rename(
        columns={
            "raxyear": "eol_raxyear",
            "raxmonth": "eol_raxmonth",
            "raxtiwy": "eol_interview_year",
            "raxtiwm": "eol_interview_month",
        }
    )

    cohort = cohort.merge(eol, on="unhhidnp", how="left", suffixes=("", "_eol"))
    cohort["eol_death_year_mismatch"] = (
        cohort["radyear"].notna()
        & cohort["eol_raxyear"].notna()
        & (cohort["radyear"] != cohort["eol_raxyear"])
    ).astype(int)

    # Clean analytic cohort file.
    clean = pd.DataFrame({
        "participant_id": cohort["unhhidnp"].astype("Int64").astype(str),
        "participant_id_char": cohort["rahhidnp"].astype(str),
        "baseline_date": cohort["baseline_date"].dt.strftime("%Y-%m-%d"),
        "death_date": cohort["death_date"].dt.strftime("%Y-%m-%d"),
        "followup_end_date": cohort["followup_end_date"].dt.strftime("%Y-%m-%d"),
        "event_death": cohort["event_death"].astype(int),
        "followup_years": cohort["followup_years"],
        "death_month_missing_flag": cohort["death_month_missing_flag"].astype(int),
        "last_alive_wave": cohort["last_alive_wave"],
        "age50plus_2001": cohort["age50plus_2001"],
        "eol_death_year_mismatch": cohort["eol_death_year_mismatch"],
    })

    for src, clean_name, group, ftype, desc in FEATURE_SPECS:
        if src in cohort.columns:
            clean[clean_name] = cohort[src].values

    # Specific QC: self-reported BMI sometimes has special/impossible high values.
    if "bmi_self_report_2001" in clean.columns:
        bmi = pd.to_numeric(clean["bmi_self_report_2001"], errors="coerce")
        clean["bmi_self_report_2001_flag_impossible_99"] = (bmi >= 90).astype(int)
        clean.loc[bmi >= 90, "bmi_self_report_2001"] = np.nan

    # Feature manifest.
    manifest_rows = []
    for src, clean_name, group, ftype, desc in FEATURE_SPECS:
        if src not in cohort.columns:
            status = "missing_from_file"
            frac = np.nan
            use = "no"
        else:
            status = "found"
            frac = float(pd.Series(clean[clean_name]).notna().mean())
            use = "yes" if group not in ["id", "design"] and frac >= args.min_feature_nonmissing else "no"

        manifest_rows.append({
            "source_variable": src,
            "clean_variable": clean_name,
            "group": group,
            "type": ftype,
            "status": status,
            "nonmissing_fraction_in_analytic_cohort": frac,
            "stata_label": h_labels.get(src, ""),
            "description": desc,
            "use_in_primary_epoch": use,
        })

    feature_manifest = pd.DataFrame(manifest_rows)

    # Model-ready file for QC/prototyping. Step 2 should redo train-fold imputation/scaling.
    id_outcome_cols = [
        "participant_id", "baseline_date", "death_date",
        "followup_end_date", "event_death", "followup_years"
    ]
    model_parts = [clean[id_outcome_cols].copy()]

    for _, row in feature_manifest.iterrows():
        if row["use_in_primary_epoch"] != "yes":
            continue

        var = row["clean_variable"]
        typ = row["type"]

        if typ == "continuous":
            x = winsorize_series(clean[var], 0.005, 0.995)
            med = x.median()
            if pd.isna(med):
                continue
            model_parts.append(pd.DataFrame({
                var: x.fillna(med),
                f"{var}__missing": x.isna().astype(int),
            }))
        elif typ == "categorical":
            s = clean[var].astype("string").fillna("MISSING")
            counts = s.value_counts(dropna=False)
            keep = set(counts.head(20).index)
            s = s.where(s.isin(keep), "OTHER")
            dummies = pd.get_dummies(s, prefix=var, dummy_na=False, dtype=int)
            if dummies.shape[1] > 1:
                dummies = dummies.iloc[:, 1:]  # drop first category as reference
            model_parts.append(dummies)

    model = pd.concat(model_parts, axis=1)
    model_feature_cols = [c for c in model.columns if c not in id_outcome_cols]
    zero_var = [c for c in model_feature_cols if model[c].nunique(dropna=False) <= 1]
    if zero_var:
        model = model.drop(columns=zero_var)

    # Downstream disease labels.
    base_ids = pd.DataFrame({"participant_id": cohort["unhhidnp"].astype("Int64").astype(str)})
    disease_tables = []

    for disease, vars_ in DISEASE_MAP.items():
        out = base_ids.copy()

        for v in vars_:
            if v in cohort.columns:
                w = int(v[1])
                out[f"{disease}_w{w}"] = cohort[v].values

        baseline_col = f"{disease}_w1"
        if baseline_col in out.columns:
            baseline = (pd.to_numeric(out[baseline_col], errors="coerce") == 1)
            out[f"{disease}_baseline_w1"] = baseline.astype(int)
        else:
            baseline = pd.Series([False] * len(out))
            out[f"{disease}_baseline_w1"] = np.nan

        future_cols = [c for c in out.columns if re.match(fr"{re.escape(disease)}_w[2-6]$", c)]
        if future_cols:
            future = out[future_cols].apply(lambda s: pd.to_numeric(s, errors="coerce") == 1).any(axis=1)
            out[f"{disease}_future_any_w2_w6"] = future.astype(int)
            out[f"{disease}_incident_after_w1"] = np.where(baseline, 0, future.astype(int))
        else:
            out[f"{disease}_future_any_w2_w6"] = np.nan
            out[f"{disease}_incident_after_w1"] = np.nan

        disease_tables.append(out.set_index("participant_id"))

    disease_labels = pd.concat(disease_tables, axis=1).reset_index()
    disease_labels = disease_labels.loc[:, ~disease_labels.columns.duplicated()]

    # Save outputs.
    shape_df = pd.DataFrame(shape_rows)
    clean_out = out_dir / "mhas_2001_step1_clean_analytic_cohort.tsv"
    model_out = out_dir / "mhas_2001_step1_model_input_primary_nondisease.tsv"
    manifest_out = out_dir / "mhas_2001_step1_primary_feature_manifest.tsv"
    disease_out = out_dir / "mhas_2001_step1_downstream_disease_labels.tsv"
    shape_out = out_dir / "mhas_step1_input_file_shapes.tsv"
    audit_out = out_dir / "mhas_2001_step1_audit_summary.txt"

    clean.to_csv(clean_out, sep="\t", index=False)
    model.to_csv(model_out, sep="\t", index=False)
    feature_manifest.to_csv(manifest_out, sep="\t", index=False)
    disease_labels.to_csv(disease_out, sep="\t", index=False)
    shape_df.to_csv(shape_out, sep="\t", index=False)

    # Audit summary.
    n_all = len(h)
    n_wave1_alive = int(((h["inw1"] == 1) & (h["r1iwstat"] == 1)).sum())
    n_clean = len(clean)
    n_deaths = int(clean["event_death"].sum())
    n_censored = n_clean - n_deaths
    event_rate = n_deaths / n_clean if n_clean else np.nan
    median_fu = clean["followup_years"].median()
    max_fu = clean["followup_years"].max()
    n_month_imp = int(clean["death_month_missing_flag"].sum())
    n_eol_mismatch = int(clean["eol_death_year_mismatch"].sum())
    n_features_used = int((feature_manifest["use_in_primary_epoch"] == "yes").sum())
    n_model_features = model.shape[1] - len(id_outcome_cols)

    feature_group_counts = (
        feature_manifest[feature_manifest["use_in_primary_epoch"] == "yes"]
        .groupby("group")
        .size()
        .to_string()
    )

    audit = f"""MHAS Step 1 cleaned analytic cohort for phenotype-based mortality EPOCH

Input file structure
--------------------
Base directory: {base_dir}
Main harmonized MHAS: {h_file}
End-of-life file: {eol_file}
Mex-Cog file present: {mexcog_file.exists()}

Input DTA shapes
----------------
{shape_df.to_string(index=False)}

Cohort definition
-----------------
Baseline wave: Wave 1 / 2001
Included participants: inw1 == 1 and r1iwstat == 1
Age restriction applied: {args.restrict_age50}
Outcome: all-cause mortality after baseline
Death date: radyear/radmonth from Gateway Harmonized MHAS; day assigned as 15
If radmonth missing but radyear present: month assigned as July and death_month_missing_flag=1
Censoring: latest wave date where rWiwstat is 1 or 4
Survival time: followup_end_date - baseline_date, in years

Cohort summary
--------------
All rows in H_MHAS_d: {n_all:,}
Wave-1 alive respondents before date/follow-up QC: {n_wave1_alive:,}
Final analytic cohort: {n_clean:,}
Deaths after baseline: {n_deaths:,}
Censored/non-deaths: {n_censored:,}
Event rate: {event_rate:.4f}
Median follow-up years: {median_fu:.2f}
Maximum follow-up years: {max_fu:.2f}
Deaths with month imputed from year only: {n_month_imp:,}
Death-year mismatches between H_MHAS radyear and EOL raxyear: {n_eol_mismatch:,}

Primary non-disease feature set
-------------------------------
Number of source features used in primary mortality EPOCH: {n_features_used}
Number of model-ready feature columns after dummy coding/missingness indicators: {n_model_features}

Feature groups used:
{feature_group_counts}

Disease-label policy
--------------------
Doctor-diagnosed disease variables are excluded from the primary mortality EPOCH features.
They are saved separately for downstream prediction in:
{disease_out}

Outputs
-------
Clean analytic cohort:
{clean_out}

Model input for Step 2:
{model_out}

Feature manifest:
{manifest_out}

Downstream disease labels:
{disease_out}

Audit:
{audit_out}

Important note for Step 2
-------------------------
The model-input file is useful for QC and initial modeling, but the final Cox elastic-net
training script should redo imputation and standardization within each training fold to
avoid leakage during cross-validation.
"""
    audit_out.write_text(audit)

    log("\n" + audit)
    log("STEP 1 finished successfully.")


if __name__ == "__main__":
    main()
