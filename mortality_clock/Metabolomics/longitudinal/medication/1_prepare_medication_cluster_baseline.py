#!/usr/bin/env python3
# ============================================================
# Prepare medication clusters and merge with longitudinal
# metabolomics mortality-clock delta biomarkers.
#
# Medication source:
#   /cbica/home/wenju/Dataset/UKBB_UMelbourne/medication.xlsx
#
# Medication coding:
#   coding4.tsv
#
# Medication field:
#   UK Biobank field 20003, instance 0:
#     20003-0.*
#
# Output:
#   1) medication_instance0_long_classified.tsv
#   2) medication_code_classification.tsv
#   3) medication_participant_clusters.tsv
#   4) medication_cluster_summary.tsv
#   5) metabolomics_delta_clock_medication_cluster_wide.tsv
#   6) metabolomics_delta_clock_medication_cluster_long.tsv
#   7) metabolomics_delta_clock_medication_cluster_requested5_wide.tsv
#   8) metabolomics_delta_clock_medication_cluster_requested5_long.tsv
#
# Downstream model:
#   delta_clock ~ medication_cluster + baseline_clock + covariates
# ============================================================

import argparse
import os
import re
import json
import warnings
from typing import Dict, List

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")


ORGANS: Dict[str, str] = {
    "Endocrine": "endocrine",
    "Digestive": "digestive",
    "Hepatic": "hepatic",
    "Immune": "immune",
}

REQUESTED_CLUSTER_LEVELS = [
    "No/minimal medication",
    "Cardiometabolic medication cluster",
    "Respiratory medication cluster",
    "Psychiatric/pain medication cluster",
    "High polypharmacy cluster",
]

ALL_CLUSTER_LEVELS = REQUESTED_CLUSTER_LEVELS + [
    "Other/uncategorized medication"
]

AGE_RECRUIT_COL = "age_at_recruitment_f21022_0_0"
SMOKING_COL = "smoking_status_f20116_0_0"
BMI_COL = "body_mass_index_bmi_f23104_0_0"
DIASTOLIC_COL = "diastolic_blood_pressure_automated_reading_f4079_0_0"
SYSTOLIC_COL = "systolic_blood_pressure_automated_reading_f4080_0_0"


# ============================================================
# 1. Argument parser
# ============================================================

def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Create baseline medication clusters and merge them with "
            "longitudinal metabolomics mortality-clock delta biomarkers."
        )
    )

    p.add_argument(
        "--medication_xlsx",
        default="/cbica/home/wenju/Dataset/UKBB_UMelbourne/medication.xlsx",
        help="UMelbourne UKBB medication Excel file."
    )

    p.add_argument(
        "--coding4_tsv",
        default="/cbica/home/wenju/Reproducibile_paper/Multiorgan_Subtype/output/medication/data/Ye_UMelbourne_100075/coding4.tsv",
        help="UKB coding 4 mapping TSV with columns coding and meaning."
    )

    p.add_argument(
        "--umel_match_csv",
        default="/cbica/home/wenju/Dataset/UKBB_UMelbourne/UKB_UMelbourne_vs_Penn_match_key.csv",
        help="Mapping key from UMelbourne ID to Penn/UPenn participant ID."
    )

    p.add_argument(
        "--delta_root",
        default=(
            "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/"
            "mortality_clock/longitudinal/metabolomics/"
            "metabolomics_delta_acceleration_years_landmark_survival_analysis"
        ),
        help="Directory containing organ-specific *_wide_delta_acceleration_years.tsv files."
    )

    p.add_argument(
        "--cov_tsv",
        default="/cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv",
        help="Optional covariate CSV with eid or participant_id. Set to NONE to skip."
    )

    p.add_argument(
        "--out_dir",
        default=(
            "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/"
            "mortality_clock/SA/delta_metabolomics_algorithmic_disease_onset/"
            "medication_cluster_delta_clock_inputs"
        ),
        help="Output directory."
    )

    p.add_argument(
        "--polypharmacy_threshold",
        type=int,
        default=5,
        help="Major medication count threshold for high polypharmacy."
    )

    p.add_argument(
        "--keep_other_cluster",
        action="store_true",
        help=(
            "Keep Other/uncategorized medication in requested downstream files. "
            "By default requested5 files include only the five requested clusters."
        )
    )

    return p.parse_args()


# ============================================================
# 2. Basic helpers
# ============================================================

def normalize_id_series(s: pd.Series) -> pd.Series:
    out = pd.to_numeric(s, errors="coerce")
    return out


def normalize_participant_id(df: pd.DataFrame, col: str) -> pd.DataFrame:
    if col not in df.columns:
        raise ValueError(f"Missing ID column: {col}")
    out = df.copy()
    out[col] = pd.to_numeric(out[col], errors="coerce")
    out = out[out[col].notna()].copy()
    out[col] = out[col].astype(np.int64)
    return out


def normalize_code_value(x):
    """
    Convert medication code values to integer strings.

    Handles:
      1140868226
      1.140868226e+09
      '1140868226'
      NaN

    Removes:
      -1, -3, 0, missing
    """
    if pd.isna(x):
        return None
    try:
        val = int(float(str(x).strip()))
    except Exception:
        return None

    if val in [-1, -3, 0]:
        return None
    return str(val)


def contains_any(text: str, patterns: List[str]) -> bool:
    if text is None or pd.isna(text):
        return False
    text = str(text).lower()
    for pat in patterns:
        if re.search(pat, text, flags=re.IGNORECASE):
            return True
    return False


# ============================================================
# 3. Medication classification keyword rules
# ============================================================

SUPPLEMENT_PATTERNS = [
    r"\bvitamin\b",
    r"\bmultivitamin\b",
    r"\bsupplement\b",
    r"\bfood supplement\b",
    r"\bherbal\b",
    r"\bplant\b",
    r"\bextract\b",
    r"\bfish oil\b",
    r"\bomega[- ]?3\b",
    r"\bevening primrose\b",
    r"\bco[- ]?enzyme q10\b",
    r"\bcoenzyme q10\b",
    r"\bubiquinone\b",
    r"\bchondroitin\b",
    r"\bglucosamine\b",
    r"\bmineral\b",
    r"\bcalcium product\b",
    r"\biron product\b",
    r"\bzinc\b",
    r"\bselenium\b",
    r"\bst john",
]

CARDIOMETABOLIC_PATTERNS = [
    # Lipid-lowering / statins
    r"\bstatin\b",
    r"\bsimvastatin\b",
    r"\batorvastatin\b",
    r"\bpravastatin\b",
    r"\brosuvastatin\b",
    r"\bfluvastatin\b",
    r"\bezeti?mibe\b",
    r"\bfibrate\b",
    r"\bfenofibrate\b",
    r"\bgemfibrozil\b",
    r"\bcholestyramine\b",
    r"\blipid\b",
    r"\bcholesterol\b",

    # Antihypertensives
    r"\bblood pressure\b",
    r"\bhypertension\b",
    r"\bantihypertensive\b",
    r"\bace inhibitor\b",
    r"\bangiotensin\b",
    r"\bramipril\b",
    r"\blisinopril\b",
    r"\benalapril\b",
    r"\bperindopril\b",
    r"\bcaptopril\b",
    r"\blosartan\b",
    r"\bcandesartan\b",
    r"\bvalsartan\b",
    r"\birbesartan\b",
    r"\btelmisartan\b",
    r"\bolmesartan\b",
    r"\bamlodipine\b",
    r"\bnifedipine\b",
    r"\bdiltiazem\b",
    r"\bverapamil\b",
    r"\bbeta[- ]?blocker\b",
    r"\batenolol\b",
    r"\bbisoprolol\b",
    r"\bmetoprolol\b",
    r"\bpropranolol\b",
    r"\bcarvedilol\b",
    r"\bnebivolol\b",
    r"\bdiuretic\b",
    r"\bfurosemide\b",
    r"\bfrusemide\b",
    r"\bbendroflumethiazide\b",
    r"\bhydrochlorothiazide\b",
    r"\bindapamide\b",
    r"\bspironolactone\b",
    r"\bamiloride\b",
    r"\bdoxazosin\b",

    # Diabetes / endocrine-metabolic
    r"\bdiabetes\b",
    r"\bdiabetic\b",
    r"\bmetformin\b",
    r"\binsulin\b",
    r"\bgliclazide\b",
    r"\bglipizide\b",
    r"\bglimepiride\b",
    r"\btolbutamide\b",
    r"\bpioglitazone\b",
    r"\bacarbose\b",
    r"\bglibenclamide\b",

    # Antiplatelet / anticoagulant / cardiovascular
    r"\baspirin\b",
    r"\bclopidogrel\b",
    r"\bdipyridamole\b",
    r"\bwarfarin\b",
    r"\bapixaban\b",
    r"\brivaroxaban\b",
    r"\bdabigatran\b",
    r"\bedoxaban\b",
    r"\banticoagulant\b",
    r"\bantiplatelet\b",
    r"\bnitrate\b",
    r"\bglyceryl trinitrate\b",
    r"\bgt[n]?\b",
    r"\bisosorbide\b",
    r"\bdigoxin\b",
    r"\bamiodarone\b",
    r"\bsotalol\b",
    r"\bangina\b",
    r"\bheart\b",
    r"\bcardiac\b",
]

RESPIRATORY_PATTERNS = [
    r"\binhaler\b",
    r"\bbronchodilator\b",
    r"\bsalbutamol\b",
    r"\balbuterol\b",
    r"\bventolin\b",
    r"\bterbutaline\b",
    r"\bsalmeterol\b",
    r"\bformoterol\b",
    r"\bipratropium\b",
    r"\btiotropium\b",
    r"\batrovent\b",
    r"\bspiriva\b",
    r"\btheophylline\b",
    r"\baminophylline\b",
    r"\bmontelukast\b",
    r"\bbeclometasone\b",
    r"\bbeclomethasone\b",
    r"\bbudesonide\b",
    r"\bfluticasone\b",
    r"\bciclesonide\b",
    r"\bseretide\b",
    r"\bsymbicort\b",
    r"\bqvar\b",
    r"\basthma\b",
    r"\bcopd\b",
    r"\bchronic obstructive\b",
    r"\brespiratory\b",
]

PSYCH_PAIN_PATTERNS = [
    # Psychiatric
    r"\bantidepressant\b",
    r"\bssri\b",
    r"\bsnri\b",
    r"\bcitalopram\b",
    r"\bescitalopram\b",
    r"\bfluoxetine\b",
    r"\bsertraline\b",
    r"\bparoxetine\b",
    r"\bvenlafaxine\b",
    r"\bduloxetine\b",
    r"\bmirtazapine\b",
    r"\bamitriptyline\b",
    r"\bnortriptyline\b",
    r"\bdosulepin\b",
    r"\bdothiepin\b",
    r"\btrazodone\b",
    r"\blithium\b",
    r"\bantipsychotic\b",
    r"\bolanzapine\b",
    r"\bquetiapine\b",
    r"\brisperidone\b",
    r"\bhaloperidol\b",
    r"\bchlorpromazine\b",
    r"\banxiolytic\b",
    r"\bdiazepam\b",
    r"\btemazepam\b",
    r"\blorazepam\b",
    r"\bbenzodiazepine\b",
    r"\bzopiclone\b",
    r"\bzolpidem\b",
    r"\bsleeping tablet\b",

    # Pain / analgesic / anti-inflammatory
    r"\banalgesic\b",
    r"\bpain\b",
    r"\bparacetamol\b",
    r"\bacetaminophen\b",
    r"\bibuprofen\b",
    r"\bnaproxen\b",
    r"\bdiclofenac\b",
    r"\bcelecoxib\b",
    r"\betoricoxib\b",
    r"\bnsaid\b",
    r"\banti[- ]?inflammatory\b",
    r"\bcodeine\b",
    r"\btramadol\b",
    r"\bmorphine\b",
    r"\boxycodone\b",
    r"\bfentanyl\b",
    r"\bdihydrocodeine\b",
    r"\bgabapentin\b",
    r"\bpregabalin\b",
]


def classify_meaning(meaning: str) -> Dict[str, int]:
    text = "" if pd.isna(meaning) else str(meaning).lower()

    is_supplement = contains_any(text, SUPPLEMENT_PATTERNS)
    is_cardiometabolic = contains_any(text, CARDIOMETABOLIC_PATTERNS)
    is_respiratory = contains_any(text, RESPIRATORY_PATTERNS)
    is_psych_pain = contains_any(text, PSYCH_PAIN_PATTERNS)

    return {
        "is_supplement": int(is_supplement),
        "is_cardiometabolic": int(is_cardiometabolic),
        "is_respiratory": int(is_respiratory),
        "is_psych_pain": int(is_psych_pain),
        "is_any_named_cluster_class": int(
            is_cardiometabolic or is_respiratory or is_psych_pain
        ),
    }


# ============================================================
# 4. Read medication coding
# ============================================================

def read_coding4(path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        raise FileNotFoundError(path)

    coding = pd.read_csv(path, sep="\t")
    if "coding" not in coding.columns or "meaning" not in coding.columns:
        raise ValueError("coding4_tsv must contain columns: coding, meaning")

    coding = coding.copy()
    coding["coding_str"] = coding["coding"].apply(normalize_code_value)
    coding = coding[coding["coding_str"].notna()].copy()
    coding["meaning"] = coding["meaning"].astype(str)

    class_df = coding["meaning"].apply(classify_meaning).apply(pd.Series)
    coding = pd.concat([coding, class_df], axis=1)

    coding = coding.drop_duplicates("coding_str", keep="first")
    return coding


# ============================================================
# 5. Read ID mapping
# ============================================================

def read_id_mapping(path: str) -> pd.DataFrame:
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

    match = match[["participant_id_umel", "participant_id"]].drop_duplicates()
    return match


# ============================================================
# 6. Read and reshape baseline medication field 20003-0.*
# ============================================================

def read_medication_instance0_long(
    medication_xlsx: str,
    coding4: pd.DataFrame,
    id_map: pd.DataFrame
) -> (pd.DataFrame, pd.DataFrame):
    if not os.path.exists(medication_xlsx):
        raise FileNotFoundError(medication_xlsx)

    print("Reading medication Excel file:")
    print("  ", medication_xlsx)

    med = pd.read_excel(medication_xlsx, engine="openpyxl")

    if "eid" not in med.columns:
        raise ValueError("Medication Excel file must contain eid.")

    med = med.rename(columns={"eid": "participant_id_umel"})
    med = normalize_participant_id(med, "participant_id_umel")

    med_cols = [
        c for c in med.columns
        if re.match(r"^20003-0\.[0-9]+$", str(c))
    ]

    if len(med_cols) == 0:
        raise ValueError(
            "No instance-0 medication columns found. Expected columns like 20003-0.0, 20003-0.1, ..."
        )

    print(f"Found {len(med_cols)} baseline medication columns matching 20003-0.*")

    # Base participant table, including people with no medication codes.
    participant_base = med[["participant_id_umel"]].drop_duplicates()
    participant_base = participant_base.merge(
        id_map,
        on="participant_id_umel",
        how="left"
    )

    # Long medication table.
    long = med[["participant_id_umel"] + med_cols].melt(
        id_vars="participant_id_umel",
        value_vars=med_cols,
        var_name="medication_field",
        value_name="raw_medication_code"
    )

    long["coding_str"] = long["raw_medication_code"].apply(normalize_code_value)
    long = long[long["coding_str"].notna()].copy()

    long = long.merge(id_map, on="participant_id_umel", how="left")
    long = long.merge(
        coding4[
            [
                "coding_str",
                "coding",
                "meaning",
                "is_supplement",
                "is_cardiometabolic",
                "is_respiratory",
                "is_psych_pain",
                "is_any_named_cluster_class",
            ]
        ],
        on="coding_str",
        how="left"
    )

    # Unmatched codes are treated as non-supplement major medications
    # but are not assigned to a named class.
    for c in [
        "is_supplement",
        "is_cardiometabolic",
        "is_respiratory",
        "is_psych_pain",
        "is_any_named_cluster_class",
    ]:
        long[c] = long[c].fillna(0).astype(int)

    long["meaning"] = long["meaning"].fillna("UNKNOWN_CODE_NOT_IN_CODING4")
    long["is_major_medication"] = (long["is_supplement"] == 0).astype(int)
    long["coding_int"] = pd.to_numeric(long["coding_str"], errors="coerce")

    long = long.drop_duplicates(
        ["participant_id_umel", "participant_id", "coding_str"],
        keep="first"
    )

    return participant_base, long


# ============================================================
# 7. Create participant-level clusters
# ============================================================

def dominant_cluster(row, polypharmacy_threshold: int) -> str:
    major_count = row["major_med_count"]

    if major_count >= polypharmacy_threshold:
        return "High polypharmacy cluster"

    if major_count == 0:
        return "No/minimal medication"

    class_counts = {
        "Cardiometabolic medication cluster": row["cardiometabolic_med_count"],
        "Respiratory medication cluster": row["respiratory_med_count"],
        "Psychiatric/pain medication cluster": row["psych_pain_med_count"],
    }

    max_count = max(class_counts.values())

    if max_count == 0:
        return "Other/uncategorized medication"

    # Tie-breaking rule:
    #   cardiometabolic > respiratory > psychiatric/pain
    # This is only for assigning a single dominant cluster.
    for cluster in [
        "Cardiometabolic medication cluster",
        "Respiratory medication cluster",
        "Psychiatric/pain medication cluster",
    ]:
        if class_counts[cluster] == max_count:
            return cluster

    return "Other/uncategorized medication"


def create_participant_clusters(
    participant_base: pd.DataFrame,
    long: pd.DataFrame,
    polypharmacy_threshold: int
) -> pd.DataFrame:
    if long.empty:
        raise ValueError("Medication long table is empty after filtering valid codes.")

    summary = (
        long
        .groupby(["participant_id_umel", "participant_id"], dropna=False)
        .agg(
            total_med_count=("coding_str", "nunique"),
            supplement_count=("is_supplement", "sum"),
            major_med_count=("is_major_medication", "sum"),
            cardiometabolic_med_count=("is_cardiometabolic", "sum"),
            respiratory_med_count=("is_respiratory", "sum"),
            psych_pain_med_count=("is_psych_pain", "sum"),
            named_cluster_med_count=("is_any_named_cluster_class", "sum"),
        )
        .reset_index()
    )

    clusters = participant_base.merge(
        summary,
        on=["participant_id_umel", "participant_id"],
        how="left"
    )

    count_cols = [
        "total_med_count",
        "supplement_count",
        "major_med_count",
        "cardiometabolic_med_count",
        "respiratory_med_count",
        "psych_pain_med_count",
        "named_cluster_med_count",
    ]

    for c in count_cols:
        clusters[c] = clusters[c].fillna(0).astype(int)

    clusters["has_cardiometabolic_med"] = (clusters["cardiometabolic_med_count"] > 0).astype(int)
    clusters["has_respiratory_med"] = (clusters["respiratory_med_count"] > 0).astype(int)
    clusters["has_psych_pain_med"] = (clusters["psych_pain_med_count"] > 0).astype(int)
    clusters["has_any_named_cluster_med"] = (clusters["named_cluster_med_count"] > 0).astype(int)
    clusters["is_high_polypharmacy"] = (
        clusters["major_med_count"] >= polypharmacy_threshold
    ).astype(int)

    class_count_mat = clusters[
        [
            "cardiometabolic_med_count",
            "respiratory_med_count",
            "psych_pain_med_count",
        ]
    ].values

    clusters["n_nonzero_named_classes"] = (class_count_mat > 0).sum(axis=1)
    clusters["has_multiple_named_classes"] = (
        clusters["n_nonzero_named_classes"] > 1
    ).astype(int)

    clusters["medication_cluster"] = clusters.apply(
        lambda r: dominant_cluster(r, polypharmacy_threshold),
        axis=1
    )

    clusters["medication_cluster"] = pd.Categorical(
        clusters["medication_cluster"],
        categories=ALL_CLUSTER_LEVELS,
        ordered=True
    )

    clusters["cluster_assignment_rule"] = np.select(
        [
            clusters["medication_cluster"] == "High polypharmacy cluster",
            clusters["medication_cluster"] == "No/minimal medication",
            clusters["medication_cluster"] == "Cardiometabolic medication cluster",
            clusters["medication_cluster"] == "Respiratory medication cluster",
            clusters["medication_cluster"] == "Psychiatric/pain medication cluster",
            clusters["medication_cluster"] == "Other/uncategorized medication",
        ],
        [
            f"major_med_count >= {polypharmacy_threshold}",
            "major_med_count == 0, allowing no medication or supplement-only use",
            "dominant named medication class among non-polypharmacy participants",
            "dominant named medication class among non-polypharmacy participants",
            "dominant named medication class among non-polypharmacy participants",
            "major medication use present but no named class matched keyword rules",
        ],
        default="not_assigned"
    )

    return clusters


# ============================================================
# 8. Read metabolomics delta-clock files
# ============================================================

def read_one_delta_file(delta_root: str, organ_label: str, organ_clean: str) -> pd.DataFrame:
    path = os.path.join(
        delta_root,
        organ_label,
        f"{organ_clean}_wide_delta_acceleration_years.tsv"
    )

    if not os.path.exists(path):
        raise FileNotFoundError(path)

    df = pd.read_csv(path, sep="\t")
    if "participant_id" not in df.columns:
        raise ValueError(f"Missing participant_id in {path}")

    df = normalize_participant_id(df, "participant_id")

    required = [
        "participant_id",
        "sample_date_instance1",
        "clock_acceleration_years_0_0",
        "clock_acceleration_years_1_0",
        "clock_age_years_0_0",
        "clock_age_years_1_0",
        "delta_accel_years_1_minus_0",
        "delta_clock_age_1_minus_0",
        "chronological_age_0_0",
        "chronological_age_1_0",
        "delta_chrono_age_1_minus_0",
    ]

    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in {path}: {missing}")

    keep = [
        "participant_id",
        "sample_date_instance1",
        "death_date_instance1",
        "admin_censor_date_instance1",
        "end_date_from_instance1",
        "time_from_instance1_years",
        "event_after_instance1",
        "case_status",
        "clock_acceleration_years_0_0",
        "clock_acceleration_years_1_0",
        "clock_acceleration_z_0_0",
        "clock_acceleration_z_1_0",
        "chronological_age_0_0",
        "chronological_age_1_0",
        "clock_age_years_0_0",
        "clock_age_years_1_0",
        "delta_accel_years_1_minus_0",
        "delta_accel_z_1_minus_0",
        "delta_chrono_age_1_minus_0",
        "delta_clock_age_1_minus_0",
    ]

    keep = [c for c in keep if c in df.columns]
    df = df[keep].copy()

    rename = {
        "sample_date_instance1": f"{organ_clean}_sample_date_instance1",
        "death_date_instance1": f"{organ_clean}_death_date_instance1",
        "admin_censor_date_instance1": f"{organ_clean}_admin_censor_date_instance1",
        "end_date_from_instance1": f"{organ_clean}_end_date_from_instance1",
        "time_from_instance1_years": f"{organ_clean}_mortality_time_from_instance1_years",
        "event_after_instance1": f"{organ_clean}_mortality_event_after_instance1",
        "case_status": f"{organ_clean}_mortality_case_status",

        "clock_acceleration_years_0_0": f"{organ_clean}_baseline_accel_years",
        "clock_acceleration_years_1_0": f"{organ_clean}_followup_accel_years",
        "clock_acceleration_z_0_0": f"{organ_clean}_baseline_accel_z",
        "clock_acceleration_z_1_0": f"{organ_clean}_followup_accel_z",

        "clock_age_years_0_0": f"{organ_clean}_clock_age_years_0_0",
        "clock_age_years_1_0": f"{organ_clean}_clock_age_years_1_0",

        "delta_accel_years_1_minus_0": f"{organ_clean}_delta_accel_years",
        "delta_accel_z_1_minus_0": f"{organ_clean}_delta_accel_z",
        "delta_clock_age_1_minus_0": f"{organ_clean}_delta_clock_age_years",

        "chronological_age_0_0": f"{organ_clean}_chronological_age_0_0",
        "chronological_age_1_0": f"{organ_clean}_chronological_age_1_0",
        "delta_chrono_age_1_minus_0": f"{organ_clean}_delta_chrono_age_1_minus_0",
    }

    df = df.rename(columns=rename)
    df = df.drop_duplicates("participant_id", keep="first")
    return df


def coalesce_first_available(df: pd.DataFrame, cols: List[str]):
    available = [c for c in cols if c in df.columns]
    if len(available) == 0:
        return np.nan
    out = df[available[0]]
    for c in available[1:]:
        out = out.where(out.notna(), df[c])
    return out


def read_all_delta_clocks(delta_root: str) -> pd.DataFrame:
    merged = None

    for organ_label, organ_clean in ORGANS.items():
        d = read_one_delta_file(delta_root, organ_label, organ_clean)
        if merged is None:
            merged = d
        else:
            merged = merged.merge(d, on="participant_id", how="outer")

    if merged is None:
        raise ValueError("No delta-clock files were read.")

    # Shared age/date columns. Use first available value across organs.
    merged["sample_date_instance1"] = coalesce_first_available(
        merged,
        [f"{o}_sample_date_instance1" for o in ORGANS.values()]
    )

    merged["admin_censor_date_instance1"] = coalesce_first_available(
        merged,
        [f"{o}_admin_censor_date_instance1" for o in ORGANS.values()]
    )

    merged["death_date_instance1"] = coalesce_first_available(
        merged,
        [f"{o}_death_date_instance1" for o in ORGANS.values()]
    )

    for new_col, suffix in [
        ("chronological_age_0_0", "chronological_age_0_0"),
        ("chronological_age_1_0", "chronological_age_1_0"),
        ("delta_chrono_age_1_minus_0", "delta_chrono_age_1_minus_0"),
    ]:
        cols = [f"{o}_{suffix}" for o in ORGANS.values()]
        available = [c for c in cols if c in merged.columns]
        if available:
            tmp = pd.concat(
                [pd.to_numeric(merged[c], errors="coerce") for c in available],
                axis=1
            )
            merged[new_col] = tmp.bfill(axis=1).iloc[:, 0]
        else:
            merged[new_col] = np.nan

    return merged


# ============================================================
# 9. Optional covariates
# ============================================================

def read_covariates(path: str) -> pd.DataFrame:
    if path is None or str(path).upper() == "NONE":
        return pd.DataFrame()

    if not os.path.exists(path):
        print(f"WARNING: covariate file not found, skipping: {path}")
        return pd.DataFrame()

    cov = pd.read_csv(path)

    if "eid" in cov.columns:
        cov = cov.rename(columns={"eid": "participant_id"})

    if "participant_id" not in cov.columns:
        print("WARNING: covariate file does not contain eid or participant_id; skipping covariates.")
        return pd.DataFrame()

    sex_candidates = ["sex_f31_0_0", "genetic_sex_f22001_0_0", "Sex", "sex"]
    sex_col = next((c for c in sex_candidates if c in cov.columns), None)

    keep = ["participant_id"]
    rename = {}

    for src, dst in [
        (AGE_RECRUIT_COL, "Age_recruitment"),
        (SMOKING_COL, "Smoking"),
        (BMI_COL, "BMI"),
        (DIASTOLIC_COL, "Diastolic"),
        (SYSTOLIC_COL, "Systolic"),
        (sex_col, "Sex"),
    ]:
        if src is not None and src in cov.columns:
            keep.append(src)
            rename[src] = dst

    cov = cov[keep].copy().rename(columns=rename)
    cov = normalize_participant_id(cov, "participant_id")

    for c in cov.columns:
        if c != "participant_id":
            cov[c] = pd.to_numeric(cov[c], errors="coerce")

    cov = cov.drop_duplicates("participant_id", keep="first")
    return cov


# ============================================================
# 10. Long-format table for downstream modeling
# ============================================================

def make_long_analysis_table(wide: pd.DataFrame) -> pd.DataFrame:
    long_list = []

    base_cols = [
        "participant_id",
        "participant_id_umel",
        "medication_cluster",
        "total_med_count",
        "supplement_count",
        "major_med_count",
        "cardiometabolic_med_count",
        "respiratory_med_count",
        "psych_pain_med_count",
        "named_cluster_med_count",
        "has_cardiometabolic_med",
        "has_respiratory_med",
        "has_psych_pain_med",
        "has_any_named_cluster_med",
        "has_multiple_named_classes",
        "is_high_polypharmacy",
        "sample_date_instance1",
        "admin_censor_date_instance1",
        "death_date_instance1",
        "chronological_age_0_0",
        "chronological_age_1_0",
        "delta_chrono_age_1_minus_0",
    ]

    cov_cols = [
        "Age_recruitment",
        "Sex",
        "Smoking",
        "BMI",
        "Diastolic",
        "Systolic",
    ]

    available_base_cols = [c for c in base_cols + cov_cols if c in wide.columns]

    for organ_label, organ_clean in ORGANS.items():
        required = [
            f"{organ_clean}_baseline_accel_years",
            f"{organ_clean}_followup_accel_years",
            f"{organ_clean}_delta_accel_years",
            f"{organ_clean}_delta_clock_age_years",
        ]

        missing = [c for c in required if c not in wide.columns]
        if missing:
            print(f"WARNING: skipping {organ_label}; missing columns: {missing}")
            continue

        tmp = wide[available_base_cols + required].copy()
        tmp["organ_label"] = organ_label
        tmp["organ_clean"] = organ_clean

        tmp = tmp.rename(
            columns={
                f"{organ_clean}_baseline_accel_years": "baseline_accel_years",
                f"{organ_clean}_followup_accel_years": "followup_accel_years",
                f"{organ_clean}_delta_accel_years": "delta_accel_years",
                f"{organ_clean}_delta_clock_age_years": "delta_clock_age_years",
            }
        )

        long_list.append(tmp)

    if len(long_list) == 0:
        raise ValueError("Could not create long table because no organ-specific clock columns were found.")

    out = pd.concat(long_list, axis=0, ignore_index=True)
    return out


# ============================================================
# 11. Main
# ============================================================

def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print("============================================================")
    print("Medication cluster + metabolomics delta-clock preparation")
    print("============================================================")
    print("Medication xlsx:", args.medication_xlsx)
    print("Coding4 TSV:", args.coding4_tsv)
    print("ID map:", args.umel_match_csv)
    print("Delta root:", args.delta_root)
    print("Covariates:", args.cov_tsv)
    print("Output dir:", args.out_dir)
    print("Polypharmacy threshold:", args.polypharmacy_threshold)
    print("============================================================")

    # Read coding and classify medication codes.
    coding4 = read_coding4(args.coding4_tsv)

    coding4_out = os.path.join(args.out_dir, "medication_code_classification.tsv")
    coding4.to_csv(coding4_out, sep="\t", index=False)

    # Read ID mapping.
    id_map = read_id_mapping(args.umel_match_csv)

    # Read medication data and create long classified table.
    participant_base, med_long = read_medication_instance0_long(
        medication_xlsx=args.medication_xlsx,
        coding4=coding4,
        id_map=id_map
    )

    med_long_out = os.path.join(args.out_dir, "medication_instance0_long_classified.tsv")
    med_long.to_csv(med_long_out, sep="\t", index=False)

    # Participant-level medication clusters.
    clusters = create_participant_clusters(
        participant_base=participant_base,
        long=med_long,
        polypharmacy_threshold=args.polypharmacy_threshold
    )

    clusters_out = os.path.join(args.out_dir, "medication_participant_clusters.tsv")
    clusters.to_csv(clusters_out, sep="\t", index=False)

    # Cluster summary.
    cluster_summary = (
        clusters
        .groupby("medication_cluster", dropna=False)
        .agg(
            N=("participant_id_umel", "size"),
            N_with_mapped_participant_id=("participant_id", lambda x: int(pd.notna(x).sum())),
            mean_total_med_count=("total_med_count", "mean"),
            mean_major_med_count=("major_med_count", "mean"),
            mean_cardiometabolic_med_count=("cardiometabolic_med_count", "mean"),
            mean_respiratory_med_count=("respiratory_med_count", "mean"),
            mean_psych_pain_med_count=("psych_pain_med_count", "mean"),
            pct_has_multiple_named_classes=("has_multiple_named_classes", "mean"),
        )
        .reset_index()
    )

    cluster_summary["pct_has_multiple_named_classes"] = (
        cluster_summary["pct_has_multiple_named_classes"] * 100.0
    )

    cluster_summary_out = os.path.join(args.out_dir, "medication_cluster_summary.tsv")
    cluster_summary.to_csv(cluster_summary_out, sep="\t", index=False)

    # Read delta-clock files.
    delta = read_all_delta_clocks(args.delta_root)

    # Merge with medication clusters.
    clusters_for_merge = clusters[clusters["participant_id"].notna()].copy()
    clusters_for_merge["participant_id"] = clusters_for_merge["participant_id"].astype(np.int64)

    wide = delta.merge(
        clusters_for_merge,
        on="participant_id",
        how="left"
    )

    # Optional covariates.
    cov = read_covariates(args.cov_tsv)
    if not cov.empty:
        wide = wide.merge(cov, on="participant_id", how="left")

    # Add a medication data availability marker.
    wide["has_medication_data"] = wide["participant_id_umel"].notna().astype(int)
    wide["medication_cluster"] = wide["medication_cluster"].astype(object)
    wide.loc[wide["medication_cluster"].isna(), "medication_cluster"] = "Medication data unavailable"

    wide_out = os.path.join(
        args.out_dir,
        "metabolomics_delta_clock_medication_cluster_wide.tsv"
    )
    wide.to_csv(wide_out, sep="\t", index=False)

    # Long table for modeling.
    long_analysis = make_long_analysis_table(wide)

    long_out = os.path.join(
        args.out_dir,
        "metabolomics_delta_clock_medication_cluster_long.tsv"
    )
    long_analysis.to_csv(long_out, sep="\t", index=False)

    # Requested five-cluster-only files.
    if args.keep_other_cluster:
        requested_keep = REQUESTED_CLUSTER_LEVELS + ["Other/uncategorized medication"]
    else:
        requested_keep = REQUESTED_CLUSTER_LEVELS

    wide_requested = wide[wide["medication_cluster"].isin(requested_keep)].copy()
    long_requested = long_analysis[long_analysis["medication_cluster"].isin(requested_keep)].copy()

    wide_requested_out = os.path.join(
        args.out_dir,
        "metabolomics_delta_clock_medication_cluster_requested5_wide.tsv"
    )
    long_requested_out = os.path.join(
        args.out_dir,
        "metabolomics_delta_clock_medication_cluster_requested5_long.tsv"
    )

    wide_requested.to_csv(wide_requested_out, sep="\t", index=False)
    long_requested.to_csv(long_requested_out, sep="\t", index=False)

    # Save run metadata.
    metadata = {
        "medication_xlsx": args.medication_xlsx,
        "coding4_tsv": args.coding4_tsv,
        "umel_match_csv": args.umel_match_csv,
        "delta_root": args.delta_root,
        "cov_tsv": args.cov_tsv,
        "out_dir": args.out_dir,
        "polypharmacy_threshold": args.polypharmacy_threshold,
        "medication_field_used": "20003-0.*",
        "requested_cluster_levels": REQUESTED_CLUSTER_LEVELS,
        "all_cluster_levels": ALL_CLUSTER_LEVELS,
        "outputs": {
            "medication_code_classification": coding4_out,
            "medication_instance0_long_classified": med_long_out,
            "medication_participant_clusters": clusters_out,
            "medication_cluster_summary": cluster_summary_out,
            "wide_all": wide_out,
            "long_all": long_out,
            "wide_requested5": wide_requested_out,
            "long_requested5": long_requested_out,
        },
        "n_medication_file_participants": int(participant_base.shape[0]),
        "n_valid_long_medication_rows": int(med_long.shape[0]),
        "n_clustered_participants": int(clusters.shape[0]),
        "n_delta_clock_participants": int(delta.shape[0]),
        "n_merged_wide_rows": int(wide.shape[0]),
        "n_merged_long_rows": int(long_analysis.shape[0]),
        "n_requested5_wide_rows": int(wide_requested.shape[0]),
        "n_requested5_long_rows": int(long_requested.shape[0]),
    }

    metadata_out = os.path.join(args.out_dir, "run_metadata.json")
    with open(metadata_out, "w") as f:
        json.dump(metadata, f, indent=2, default=str)

    print("\n============================================================")
    print("Finished.")
    print("Outputs:")
    print("  medication code classification:", coding4_out)
    print("  long classified medication table:", med_long_out)
    print("  participant clusters:", clusters_out)
    print("  cluster summary:", cluster_summary_out)
    print("  merged wide table:", wide_out)
    print("  merged long table:", long_out)
    print("  requested five-cluster wide table:", wide_requested_out)
    print("  requested five-cluster long table:", long_requested_out)
    print("  metadata:", metadata_out)
    print("============================================================\n")

    print("Medication cluster counts:")
    print(
        clusters["medication_cluster"]
        .value_counts(dropna=False)
        .reindex(ALL_CLUSTER_LEVELS)
        .fillna(0)
        .astype(int)
        .to_string()
    )

    print("\nRequested five-cluster counts among merged delta-clock participants:")
    print(
        wide_requested["medication_cluster"]
        .value_counts(dropna=False)
        .reindex(REQUESTED_CLUSTER_LEVELS)
        .fillna(0)
        .astype(int)
        .to_string()
    )


if __name__ == "__main__":
    main()