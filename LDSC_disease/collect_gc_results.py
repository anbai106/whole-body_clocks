#!/usr/bin/env python3

import os
import re
import glob
import argparse
import numpy as np
import pandas as pd


# ============================================================
# 1. Disease EPOCH settings
# ============================================================

DISEASE_ORDER = ["asthma", "copd", "dementia", "mi", "stroke"]

DISEASE_LABEL = {
    "asthma": "Asthma",
    "copd": "COPD",
    "dementia": "Dementia",
    "mi": "MI",
    "stroke": "Stroke",
}

MODALITY_LABEL = {
    "mri": "MRI",
    "proteomics": "Proteomics",
    "metabolomics": "Metabolomics",
}

MODALITY_ORDER = {
    "MRI": 1,
    "Proteomics": 2,
    "Metabolomics": 3,
    "Unknown": 99,
}


# ============================================================
# 2. Previous AI biomarker list
# ============================================================

PROTBAG_LIST = [
    "Reproductive_female",
    "Pulmonary",
    "Heart",
    "Brain",
    "Eye",
    "Hepatic",
    "Renal",
    "Reproductive_male",
    "Endocrine",
    "Immune",
    "Skin",
]

METBAG_LIST = [
    "Endocrine",
    "Digestive",
    "Hepatic",
    "Immune",
    "Metabolic",
]

MRIBAG_LIST = [
    "brain",
    "adipose",
    "heart",
    "kidney",
    "liver",
    "pancreas",
    "spleen",
]

MAE_BRAIN_LIST = ["r1", "r2", "r3", "r4", "r5", "r6"]
MAE_EYE_LIST = ["r1", "r2", "r3"]
MAE_HEART_LIST = ["r1", "r2"]

DNE_LIST = [
    "AD_SurrealGAN_1",
    "AD_SurrealGAN_2",
    "ASD_1",
    "ASD_2",
    "ASD_3",
    "LLD_1",
    "LLD_2",
    "SCZ_1",
    "SCZ_2",
]

MULTIMODAL_BRAIN_BAG_MAP = {
    "muse": "GM brain MRIBAG",
    "dmri_fullmetric_fa": "WM brain MRIBAG",
    "fmri": "FC brain MRIBAG",
}


def build_exact_ai_biomarker_table():
    rows = []

    for bag in PROTBAG_LIST:
        rows.append({
            "target_id": f"ProtBAG_{bag}",
            "target_category": "AI_biomarker",
            "target_source": "ProtBAG",
            "target_family": "Proteomics aging clock",
            "target_name": bag,
            "target_display": f"ProtBAG {bag}",
        })

    for bag in METBAG_LIST:
        rows.append({
            "target_id": f"MetBAG_{bag}",
            "target_category": "AI_biomarker",
            "target_source": "MetBAG",
            "target_family": "Metabolomics aging clock",
            "target_name": bag,
            "target_display": f"MetBAG {bag}",
        })

    for bag in MRIBAG_LIST:
        rows.append({
            "target_id": f"MRIBAG_{bag}",
            "target_category": "AI_biomarker",
            "target_source": "MRIBAG",
            "target_family": "Organ MRI aging clock",
            "target_name": bag,
            "target_display": f"MRIBAG {bag}",
        })

    for mae in MAE_BRAIN_LIST:
        rows.append({
            "target_id": f"{mae}_brain",
            "target_category": "AI_biomarker",
            "target_source": "Pan_disease_MAE",
            "target_family": "Brain pan-disease MAE",
            "target_name": mae,
            "target_display": f"{mae.upper()} brain MAE",
        })

    for mae in MAE_EYE_LIST:
        rows.append({
            "target_id": f"{mae}_eye",
            "target_category": "AI_biomarker",
            "target_source": "Pan_disease_MAE",
            "target_family": "Eye pan-disease MAE",
            "target_name": mae,
            "target_display": f"{mae.upper()} eye MAE",
        })

    for mae in MAE_HEART_LIST:
        rows.append({
            "target_id": f"{mae}_heart",
            "target_category": "AI_biomarker",
            "target_source": "Pan_disease_MAE",
            "target_family": "Heart pan-disease MAE",
            "target_name": mae,
            "target_display": f"{mae.upper()} heart MAE",
        })

    for dne in DNE_LIST:
        if dne.startswith("AD_SurrealGAN"):
            family = "AD DNE subtype"
        elif dne.startswith("ASD"):
            family = "ASD DNE subtype"
        elif dne.startswith("LLD"):
            family = "LLD DNE subtype"
        elif dne.startswith("SCZ"):
            family = "SCZ DNE subtype"
        else:
            family = "DNE subtype"

        rows.append({
            "target_id": dne,
            "target_category": "AI_biomarker",
            "target_source": "DNE",
            "target_family": family,
            "target_name": dne,
            "target_display": dne.replace("_", " "),
        })

    for target_id, display in MULTIMODAL_BRAIN_BAG_MAP.items():
        rows.append({
            "target_id": target_id,
            "target_category": "AI_biomarker",
            "target_source": "Brain_MRIBAG",
            "target_family": "Multi-modal brain MRI aging clock",
            "target_name": target_id,
            "target_display": display,
        })

    return pd.DataFrame(rows)


AI_BIOMARKER_DF = build_exact_ai_biomarker_table()
AI_BIOMARKER_LOOKUP = AI_BIOMARKER_DF.set_index("target_id").to_dict(orient="index")


# ============================================================
# 3. Optional disease endpoint manifest loading
# ============================================================

PGC_FALLBACK_TARGETS = {
    "AD",
    "ADHD",
    "ASD",
    "BIP",
    "SCZ",
    "OCD",
    "AUD",
    "AN",
    "MDD",
    "PTSD",
    "TS",
}


def load_finngen_targets(manifest_path):
    if not manifest_path or not os.path.exists(manifest_path):
        print(f"WARNING: FinnGen manifest not found: {manifest_path}")
        return set()

    try:
        df = pd.read_csv(manifest_path, sep="\t")

        if "phenocode" in df.columns:
            return set(df["phenocode"].dropna().astype(str).tolist())

        first_col = df.columns[0]
        vals = df[first_col].dropna().astype(str).tolist()
        vals = [x for x in vals if x != "phenocode"]

        return set(vals)

    except Exception as e:
        print(f"WARNING: Could not read FinnGen manifest: {e}")
        return set()


def load_pgc_targets(manifest_path):
    if not manifest_path or not os.path.exists(manifest_path):
        print(f"WARNING: PGC manifest not found: {manifest_path}")
        return set()

    try:
        df = pd.read_csv(manifest_path, sep="\t")

        if "Phenotype" in df.columns:
            return set(df["Phenotype"].dropna().astype(str).tolist())

        if df.shape[1] >= 2:
            vals = df.iloc[:, 1].dropna().astype(str).tolist()
            vals = [x for x in vals if x != "Phenotype"]
            return set(vals)

        return set()

    except Exception as e:
        print(f"WARNING: Could not read PGC manifest: {e}")
        return set()


# ============================================================
# 4. Helper functions
# ============================================================

def safe_float(x, default=np.nan):
    try:
        if x is None:
            return default

        s = str(x).strip()

        if s.lower() in ["", "nan", "na", "none", "null"]:
            return default

        return float(s)

    except Exception:
        return default


def p_to_log10(p):
    p = safe_float(p)

    if not np.isfinite(p):
        return np.nan, np.nan

    if p == 0:
        return -np.inf, np.inf

    return np.log10(p), -np.log10(p)


def clean_token(x):
    x = str(x)
    x = x.replace("-", "_").replace(" ", "_")
    x = re.sub(r"_+", "_", x)
    return x


def organ_key_from_raw(organ_raw):
    x = clean_token(organ_raw).lower()

    organ_map = {
        "brain": "brain",
        "eye": "eye",
        "heart": "heart",
        "cardiac": "heart",

        "liver": "hepatic",
        "hepatic": "hepatic",

        "kidney": "renal",
        "renal": "renal",

        "lung": "pulmonary",
        "pulmonary": "pulmonary",

        "pancreas": "pancreas",
        "pancreatic": "pancreas",

        "spleen": "spleen",
        "adipose": "adipose",

        "endocrine": "endocrine",
        "digestive": "digestive",
        "immune": "immune",
        "metabolic": "metabolic",
        "skin": "skin",

        "reproductive": "reproductive",
        "reproductive_female": "reproductive_female",
        "female_reproductive": "reproductive_female",
        "reproductive_male": "reproductive_male",
        "male_reproductive": "reproductive_male",
    }

    return organ_map.get(x, x)


def format_organ_label(organ_raw):
    if organ_raw is None or pd.isna(organ_raw):
        return "Other"

    x = clean_token(organ_raw)

    words = []

    for token in x.split("_"):
        t = token.lower()

        if t == "female":
            words.append("female")
        elif t == "male":
            words.append("male")
        else:
            words.append(t.capitalize())

    return " ".join(words)


def infer_modality_from_clock_prefix(clock_prefix):
    x = str(clock_prefix).lower()

    if x.endswith("_mri"):
        return "MRI"
    if x.endswith("_proteomics"):
        return "Proteomics"
    if x.endswith("_metabolomics"):
        return "Metabolomics"

    return "Unknown"


def parse_disease_clock_folder(clock_folder):
    """
    Parse disease EPOCH clock folder.

    Expected folder:
      <organ>_<modality>_<disease>_clock

    Examples:
      Brain_proteomics_dementia_clock
      Endocrine_proteomics_stroke_clock
      Hepatic_metabolomics_asthma_clock
      heart_mri_copd_clock
      spleen_mri_asthma_clock
      Reproductive_male_proteomics_mi_clock
    """

    if not clock_folder.endswith("_clock"):
        raise ValueError(f"Clock folder does not end with _clock: {clock_folder}")

    x = clock_folder[:-len("_clock")]

    disease_key = None
    x_without_disease = None

    for d in DISEASE_ORDER:
        if x.lower().endswith(f"_{d}"):
            disease_key = d
            x_without_disease = re.sub(f"_{d}$", "", x, flags=re.IGNORECASE)
            break

    if disease_key is None:
        raise ValueError(f"Cannot parse disease from clock folder: {clock_folder}")

    modality = None
    organ_raw = None

    for m in ["metabolomics", "proteomics", "mri"]:
        if x_without_disease.lower().endswith(f"_{m}"):
            modality = MODALITY_LABEL[m]
            organ_raw = re.sub(f"_{m}$", "", x_without_disease, flags=re.IGNORECASE)
            break

    if modality is None or organ_raw is None:
        raise ValueError(f"Cannot parse organ/modality from clock folder: {clock_folder}")

    organ_key = organ_key_from_raw(organ_raw)
    organ_label = format_organ_label(organ_raw)
    disease_label = DISEASE_LABEL[disease_key]

    return {
        "clock_folder": clock_folder,
        "disease_clock": f"{disease_key}__{organ_key}__{modality.lower()}",
        "disease_key": disease_key,
        "disease_label": disease_label,
        "modality": modality,
        "modality_key": modality.lower(),
        "organ_key": organ_key,
        "organ_label": organ_label,
        "organ_raw": organ_raw,
        "clock_prefix": x,
    }


def parse_log_filename(log_file):
    """
    Parse pairwise LDSC log filename.

    In each clock folder, logs are organ/system-specific, e.g.:
      Endocrine_vs_I9_MI_STRICT.log
      Brain_vs_ProtBAG_Brain.log
      heart_vs_mortality_epoch_heart_mri_mortality_clock.log
      Reproductive_vs_MetBAG_Immune.log

    Only target_id is required for collection/classification.
    """

    basename = os.path.basename(log_file)
    basename_no_ext = re.sub(r"\.log$", "", basename)

    if "_vs_" not in basename_no_ext:
        raise ValueError(f"Not a pairwise LDSC log: {basename}")

    exposure_name, target_id = basename_no_ext.split("_vs_", 1)

    return exposure_name, target_id


# ============================================================
# 5. Target classification
# ============================================================

def parse_mortality_epoch_target(target_id):
    """
    Parse mortality EPOCH target IDs.

    Examples:
      mortality_epoch_Brain_proteomics_mortality_clock
      mortality_epoch_adipose_mri_mortality_clock
      mortality_epoch_Endocrine_metabolomics_mortality_clock
    """

    clock_folder = str(target_id).replace("mortality_epoch_", "", 1)
    clock_prefix = clock_folder.replace("_mortality_clock", "")

    modality = infer_modality_from_clock_prefix(clock_prefix)

    organ_raw = clock_prefix

    for suffix in ["_metabolomics", "_proteomics", "_mri"]:
        if organ_raw.lower().endswith(suffix):
            organ_raw = re.sub(f"{suffix}$", "", organ_raw, flags=re.IGNORECASE)
            break

    organ_key = organ_key_from_raw(organ_raw)
    organ_label = format_organ_label(organ_raw)

    return {
        "target_category": "Mortality_EPOCH",
        "target_source": "Mortality_EPOCH",
        "target_family": "Mortality EPOCH clock",
        "target_name": clock_folder,
        "target_display": f"Mortality EPOCH {organ_label} {modality}",
        "target_mortality_clock_folder": clock_folder,
        "target_mortality_clock_prefix": clock_prefix,
        "target_mortality_modality": modality,
        "target_mortality_organ_key": organ_key,
        "target_mortality_organ_label": organ_label,
    }


def classify_target(target_id, finngen_targets, pgc_targets):
    target = str(target_id)

    # Category 3: mortality EPOCH clocks
    if target.startswith("mortality_epoch_"):
        return parse_mortality_epoch_target(target)

    # Category 2: previous AI biomarkers
    if target in AI_BIOMARKER_LOOKUP:
        info = AI_BIOMARKER_LOOKUP[target].copy()

        info.update({
            "target_mortality_clock_folder": "",
            "target_mortality_clock_prefix": "",
            "target_mortality_modality": "",
            "target_mortality_organ_key": "",
            "target_mortality_organ_label": "",
        })

        return info

    # Category 1: disease endpoints
    if target in pgc_targets or target in PGC_FALLBACK_TARGETS:
        source = "PGC"
        family = "PGC"
    elif target in finngen_targets:
        source = "FinnGen"
        family = "FinnGen"
    else:
        source = "FinnGen_or_other"
        family = "Disease endpoint"

    return {
        "target_category": "Disease_endpoint",
        "target_source": source,
        "target_family": family,
        "target_name": target,
        "target_display": target,
        "target_mortality_clock_folder": "",
        "target_mortality_clock_prefix": "",
        "target_mortality_modality": "",
        "target_mortality_organ_key": "",
        "target_mortality_organ_label": "",
    }


# ============================================================
# 6. LDSC log parser
# ============================================================

def parse_ldsc_log(log_file):
    """
    Parse LDSC genetic correlation from .log file.

    Supported formats:
      Genetic Correlation: 0.123 (0.045)
      Z-score: 2.73
      P: 0.0063

    Fallback format:
      p1 p2 rg se z p ...
    """

    with open(log_file, "r", errors="ignore") as f:
        txt = f.read()

    gc_mean = np.nan
    gc_std = np.nan
    z_score = np.nan
    p_value = np.nan

    m_gc = re.search(
        r"Genetic Correlation:\s*([+-]?(?:nan|inf|-inf|[0-9.eE+-]+))"
        r"\s*(?:\(([+-]?(?:nan|inf|-inf|[0-9.eE+-]+))\))?",
        txt,
        flags=re.IGNORECASE,
    )

    if m_gc is not None:
        gc_mean = safe_float(m_gc.group(1))
        gc_std = safe_float(m_gc.group(2))

    m_z = re.search(
        r"^Z-score:\s*([+-]?(?:nan|inf|-inf|[0-9.eE+-]+))",
        txt,
        flags=re.IGNORECASE | re.MULTILINE,
    )

    if m_z is not None:
        z_score = safe_float(m_z.group(1))

    m_p = re.search(
        r"^P:\s*([+-]?(?:nan|inf|-inf|[0-9.eE+-]+))",
        txt,
        flags=re.IGNORECASE | re.MULTILINE,
    )

    if m_p is not None:
        p_value = safe_float(m_p.group(1))

    # LDSC table fallback.
    if not np.isfinite(gc_mean) or not np.isfinite(gc_std) or not np.isfinite(p_value):
        lines = txt.splitlines()

        for i, line in enumerate(lines):
            line_clean = line.strip()

            if re.match(r"^p1\s+p2\s+rg\s+se\s+z\s+p\b", line_clean):
                header = re.split(r"\s+", line_clean)

                if i + 1 < len(lines):
                    values = re.split(r"\s+", lines[i + 1].strip())

                    if len(values) >= len(header):
                        tmp = dict(zip(header, values))

                        gc_mean = safe_float(tmp.get("rg", gc_mean))
                        gc_std = safe_float(tmp.get("se", gc_std))
                        z_score = safe_float(tmp.get("z", z_score))
                        p_value = safe_float(tmp.get("p", p_value))

                break

    p_log10, neg_log10_p = p_to_log10(p_value)

    if np.isfinite(gc_mean) and np.isfinite(gc_std) and np.isfinite(p_value):
        parse_status = "ok"
    elif "nan" in txt.lower():
        parse_status = "nan_result"
    elif "error" in txt.lower() or "traceback" in txt.lower():
        parse_status = "ldsc_error"
    else:
        parse_status = "parse_failed"

    return {
        "gc_mean": gc_mean,
        "gc_std": gc_std,
        "Z": z_score,
        "P": p_value,
        "P_log10": p_log10,
        "neg_log10_P": neg_log10_p,
        "parse_status": parse_status,
    }


# ============================================================
# 7. Discover all 47 disease EPOCH clock folders
# ============================================================

def discover_disease_epoch_clock_folders(base_dir):
    """
    Discover all 47 disease EPOCH clock folders.

    It does not expect cross-organ logs inside any folder.
    Each folder is processed independently using its own fastGWA/ldsc/*.log files.
    """

    rows = []

    for disease in DISEASE_ORDER:
        pattern = os.path.join(base_dir, f"*_{disease}_clock")

        for clock_dir in sorted(glob.glob(pattern)):
            if not os.path.isdir(clock_dir):
                continue

            clock_folder = os.path.basename(clock_dir)
            ldsc_dir = os.path.join(clock_dir, "fastGWA", "ldsc")

            try:
                meta = parse_disease_clock_folder(clock_folder)
                parse_status = "ok"
                parse_error = ""
            except Exception as e:
                meta = {
                    "clock_folder": clock_folder,
                    "disease_clock": "",
                    "disease_key": disease,
                    "disease_label": DISEASE_LABEL.get(disease, disease),
                    "modality": "Unknown",
                    "modality_key": "unknown",
                    "organ_key": "",
                    "organ_label": "",
                    "organ_raw": "",
                    "clock_prefix": "",
                }
                parse_status = "parse_failed"
                parse_error = str(e)

            n_pairwise_logs = 0
            if os.path.isdir(ldsc_dir):
                n_pairwise_logs = len(glob.glob(os.path.join(ldsc_dir, "*_vs_*.log")))

            rows.append({
                **meta,
                "clock_dir": clock_dir,
                "ldsc_dir": ldsc_dir,
                "ldsc_dir_exists": os.path.isdir(ldsc_dir),
                "n_pairwise_logs_found": n_pairwise_logs,
                "clock_folder_parse_status": parse_status,
                "clock_folder_parse_error": parse_error,
            })

    df_clock = pd.DataFrame(rows)

    if df_clock.empty:
        raise RuntimeError(f"No disease EPOCH clock folders found under: {base_dir}")

    df_clock = (
        df_clock
        .drop_duplicates(subset=["clock_folder", "clock_dir"])
        .reset_index(drop=True)
    )

    return df_clock


def get_pairwise_logs_for_clock(ldsc_dir):
    if not os.path.isdir(ldsc_dir):
        return []

    return sorted(glob.glob(os.path.join(ldsc_dir, "*_vs_*.log")))


# ============================================================
# 8. Multiple-testing correction
# ============================================================

def add_multiple_testing_columns(df_all):
    """
    Primary Bonferroni correction:
      0.05 / number of unique targets across all three target categories:
        1. disease endpoints
        2. previous AI biomarkers
        3. mortality EPOCH clocks

    This avoids multiplying by the 47 disease clocks again and instead
    corrects by the total number of unique biomarkers/endpoints tested.
    """

    df = df_all.copy()

    df["nominal_sig"] = df["P"] < 0.05

    ok_mask = df["parse_status"].eq("ok") & np.isfinite(df["P"])

    n_valid_pairwise_tests_all_3_categories = int(ok_mask.sum())
    n_unique_targets_all_3_categories = int(df.loc[ok_mask, "target_id"].nunique())

    bonf_pairwise = (
        0.05 / n_valid_pairwise_tests_all_3_categories
        if n_valid_pairwise_tests_all_3_categories > 0
        else np.nan
    )

    bonf_unique_targets = (
        0.05 / n_unique_targets_all_3_categories
        if n_unique_targets_all_3_categories > 0
        else np.nan
    )

    df["n_valid_pairwise_tests_all_3_categories"] = n_valid_pairwise_tests_all_3_categories
    df["bonferroni_p_pairwise_tests_all_3_categories"] = bonf_pairwise
    df["bonferroni_sig_pairwise_tests_all_3_categories"] = df["P"] < bonf_pairwise

    df["n_unique_targets_all_3_categories"] = n_unique_targets_all_3_categories
    df["bonferroni_p_unique_targets_all_3_categories"] = bonf_unique_targets
    df["bonferroni_sig_unique_targets_all_3_categories"] = df["P"] < bonf_unique_targets

    # Primary flag requested for downstream interpretation.
    df["primary_bonferroni_scope"] = "unique_targets_all_3_categories"
    df["primary_bonferroni_p"] = bonf_unique_targets
    df["primary_bonferroni_sig"] = df["bonferroni_sig_unique_targets_all_3_categories"]

    # Category-specific unique target thresholds, useful for QC.
    df["n_unique_targets_target_category"] = np.nan
    df["bonferroni_p_unique_targets_target_category"] = np.nan
    df["bonferroni_sig_unique_targets_target_category"] = False

    for category, idx in df.groupby("target_category").groups.items():
        idx = list(idx)
        ok_idx = [
            i for i in idx
            if df.loc[i, "parse_status"] == "ok" and np.isfinite(df.loc[i, "P"])
        ]

        n_unique = int(df.loc[ok_idx, "target_id"].nunique())
        bonf = 0.05 / n_unique if n_unique > 0 else np.nan

        df.loc[idx, "n_unique_targets_target_category"] = n_unique
        df.loc[idx, "bonferroni_p_unique_targets_target_category"] = bonf
        df.loc[idx, "bonferroni_sig_unique_targets_target_category"] = df.loc[idx, "P"] < bonf

    return df


# ============================================================
# 9. Sorting and summary helpers
# ============================================================

def sort_gc_table(df_all):
    df = df_all.copy()

    disease_order_map = {d: i + 1 for i, d in enumerate(DISEASE_ORDER)}

    category_order = {
        "Disease_endpoint": 1,
        "AI_biomarker": 2,
        "Mortality_EPOCH": 3,
    }

    target_source_order = {
        "FinnGen": 1,
        "PGC": 2,
        "FinnGen_or_other": 3,
        "MRIBAG": 4,
        "ProtBAG": 5,
        "MetBAG": 6,
        "Brain_MRIBAG": 7,
        "Pan_disease_MAE": 8,
        "DNE": 9,
        "Mortality_EPOCH": 10,
    }

    df["disease_order"] = df["disease_key"].map(disease_order_map).fillna(99)
    df["modality_order"] = df["modality"].map(MODALITY_ORDER).fillna(99)
    df["category_order"] = df["target_category"].map(category_order).fillna(99)
    df["target_source_order"] = df["target_source"].map(target_source_order).fillna(99)

    df = (
        df.sort_values(
            [
                "disease_order",
                "modality_order",
                "organ_key",
                "clock_folder",
                "category_order",
                "target_source_order",
                "target_id",
            ]
        )
        .drop(
            columns=[
                "disease_order",
                "modality_order",
                "category_order",
                "target_source_order",
            ]
        )
        .reset_index(drop=True)
    )

    return df


def nanmean_abs(x):
    arr = np.asarray(x, dtype=float)
    arr = arr[np.isfinite(arr)]

    if arr.size == 0:
        return np.nan

    return np.mean(np.abs(arr))


def nanmax_abs(x):
    arr = np.asarray(x, dtype=float)
    arr = arr[np.isfinite(arr)]

    if arr.size == 0:
        return np.nan

    return np.max(np.abs(arr))


def save_summary_tables(df_all, output_dir):
    df_ok = df_all.loc[df_all["parse_status"] == "ok"].copy()

    summary_clock = (
        df_ok.groupby(
            [
                "clock_folder",
                "disease_clock",
                "disease_label",
                "modality",
                "organ_label",
                "target_category",
                "target_source",
            ],
            as_index=False,
        )
        .agg(
            n_tests=("target_id", "count"),
            n_unique_targets=("target_id", "nunique"),
            n_nominal_sig=("nominal_sig", "sum"),
            n_primary_bonferroni_sig=("primary_bonferroni_sig", "sum"),
            mean_rg=("gc_mean", "mean"),
            mean_abs_rg=("gc_mean", nanmean_abs),
            max_abs_rg=("gc_mean", nanmax_abs),
            min_p=("P", "min"),
            max_neg_log10_p=("neg_log10_P", "max"),
        )
        .sort_values(
            [
                "disease_label",
                "modality",
                "organ_label",
                "target_category",
                "target_source",
            ]
        )
    )

    out_clock = os.path.join(
        output_dir,
        "LDSC_gc_47_disease_epoch_summary_by_clock_and_target_category.tsv",
    )

    summary_clock.to_csv(out_clock, index=False, sep="\t", encoding="utf-8")

    summary_target = (
        df_ok.groupby(
            [
                "target_id",
                "target_display",
                "target_name",
                "target_category",
                "target_source",
                "target_family",
            ],
            as_index=False,
        )
        .agg(
            n_disease_clocks=("clock_folder", "count"),
            n_nominal_sig=("nominal_sig", "sum"),
            n_primary_bonferroni_sig=("primary_bonferroni_sig", "sum"),
            mean_rg=("gc_mean", "mean"),
            mean_abs_rg=("gc_mean", nanmean_abs),
            max_abs_rg=("gc_mean", nanmax_abs),
            min_p=("P", "min"),
            max_neg_log10_p=("neg_log10_P", "max"),
        )
        .sort_values(
            [
                "target_category",
                "target_source",
                "target_family",
                "min_p",
                "target_id",
            ]
        )
    )

    out_target = os.path.join(
        output_dir,
        "LDSC_gc_47_disease_epoch_summary_by_target.tsv",
    )

    summary_target.to_csv(out_target, index=False, sep="\t", encoding="utf-8")

    top_hits = (
        df_ok.sort_values("P")
        .loc[
            :,
            [
                "clock_folder",
                "disease_clock",
                "disease_label",
                "modality",
                "organ_label",
                "target_id",
                "target_display",
                "target_category",
                "target_source",
                "target_family",
                "gc_mean",
                "gc_std",
                "Z",
                "P",
                "neg_log10_P",
                "nominal_sig",
                "primary_bonferroni_scope",
                "primary_bonferroni_p",
                "primary_bonferroni_sig",
                "n_unique_targets_all_3_categories",
                "bonferroni_sig_unique_targets_target_category",
                "rg_abs",
                "rg_out_of_bounds",
                "log_file",
            ],
        ]
    )

    out_top = os.path.join(
        output_dir,
        "LDSC_gc_47_disease_epoch_top_hits_all_targets.tsv",
    )

    top_hits.to_csv(out_top, index=False, sep="\t", encoding="utf-8")

    print(f"Saved clock/category summary: {out_clock}")
    print(f"Saved target summary:         {out_target}")
    print(f"Saved top-hit table:          {out_top}")


# ============================================================
# 10. Save output files
# ============================================================

def save_outputs(df_all, df_clock_manifest, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    out_manifest = os.path.join(
        output_dir,
        "LDSC_gc_47_disease_epoch_clock_discovery_manifest.tsv",
    )

    out_all = os.path.join(
        output_dir,
        "LDSC_gc_47_disease_epoch_all_targets.tsv",
    )

    out_de = os.path.join(
        output_dir,
        "LDSC_gc_47_disease_epoch_vs_disease_endpoints.tsv",
    )

    out_ai = os.path.join(
        output_dir,
        "LDSC_gc_47_disease_epoch_vs_AI_biomarkers.tsv",
    )

    out_mortality = os.path.join(
        output_dir,
        "LDSC_gc_47_disease_epoch_vs_mortality_EPOCH.tsv",
    )

    out_failed = os.path.join(
        output_dir,
        "LDSC_gc_47_disease_epoch_parse_failed_or_nan.tsv",
    )

    out_ai_ref = os.path.join(
        output_dir,
        "LDSC_gc_47_disease_epoch_AI_biomarker_reference_table.tsv",
    )

    df_clock_manifest.to_csv(out_manifest, index=False, sep="\t", encoding="utf-8")
    df_all.to_csv(out_all, index=False, sep="\t", encoding="utf-8")

    df_de = df_all.loc[df_all["target_category"] == "Disease_endpoint"].copy()
    df_ai = df_all.loc[df_all["target_category"] == "AI_biomarker"].copy()
    df_mort = df_all.loc[df_all["target_category"] == "Mortality_EPOCH"].copy()
    df_failed = df_all.loc[df_all["parse_status"] != "ok"].copy()

    df_de.to_csv(out_de, index=False, sep="\t", encoding="utf-8")
    df_ai.to_csv(out_ai, index=False, sep="\t", encoding="utf-8")
    df_mort.to_csv(out_mortality, index=False, sep="\t", encoding="utf-8")
    df_failed.to_csv(out_failed, index=False, sep="\t", encoding="utf-8")
    AI_BIOMARKER_DF.to_csv(out_ai_ref, index=False, sep="\t", encoding="utf-8")

    print(f"Saved clock manifest:       {out_manifest}")
    print(f"Saved all targets:          {out_all}")
    print(f"Saved disease endpoints:    {out_de}")
    print(f"Saved AI biomarkers:        {out_ai}")
    print(f"Saved mortality EPOCH:      {out_mortality}")
    print(f"Saved failed/nan parses:    {out_failed}")
    print(f"Saved AI biomarker ref:     {out_ai_ref}")


# ============================================================
# 11. Main collection
# ============================================================

def collect_all_gc_results(
    base_dir,
    output_dir,
    finngen_manifest,
    pgc_manifest,
    expected_n_clocks,
):
    os.makedirs(output_dir, exist_ok=True)

    finngen_targets = load_finngen_targets(finngen_manifest)
    pgc_targets = load_pgc_targets(pgc_manifest)

    df_clock_manifest = discover_disease_epoch_clock_folders(base_dir)

    print("============================================================")
    print("Collecting LDSC genetic-correlation results for 47 disease EPOCH clocks")
    print("============================================================")
    print(f"Base directory:             {base_dir}")
    print(f"Output directory:           {output_dir}")
    print(f"Discovered clock folders:    {df_clock_manifest['clock_folder'].nunique()}")
    print(f"Expected disease clocks:     {expected_n_clocks}")
    print(f"FinnGen targets in manifest: {len(finngen_targets)}")
    print(f"PGC targets in manifest:     {len(pgc_targets)}")

    rows = []
    skipped = []

    for _, clock_row in df_clock_manifest.iterrows():
        clock_folder = clock_row["clock_folder"]
        ldsc_dir = clock_row["ldsc_dir"]

        if clock_row["clock_folder_parse_status"] != "ok":
            skipped.append({
                "clock_folder": clock_folder,
                "ldsc_dir": ldsc_dir,
                "log_file": "",
                "reason": "clock_folder_parse_failed",
                "details": clock_row["clock_folder_parse_error"],
            })
            continue

        if not bool(clock_row["ldsc_dir_exists"]):
            skipped.append({
                "clock_folder": clock_folder,
                "ldsc_dir": ldsc_dir,
                "log_file": "",
                "reason": "missing_ldsc_dir",
                "details": "",
            })
            continue

        log_files = get_pairwise_logs_for_clock(ldsc_dir)

        if len(log_files) == 0:
            skipped.append({
                "clock_folder": clock_folder,
                "ldsc_dir": ldsc_dir,
                "log_file": "",
                "reason": "no_pairwise_logs_found",
                "details": "",
            })
            continue

        clock_meta = {
            "clock_folder": clock_row["clock_folder"],
            "disease_clock": clock_row["disease_clock"],
            "disease_key": clock_row["disease_key"],
            "disease_label": clock_row["disease_label"],
            "modality": clock_row["modality"],
            "modality_key": clock_row["modality_key"],
            "organ_key": clock_row["organ_key"],
            "organ_label": clock_row["organ_label"],
            "organ_raw": clock_row["organ_raw"],
            "clock_prefix": clock_row["clock_prefix"],
        }

        for log_file in log_files:
            try:
                exposure_name, target_id = parse_log_filename(log_file)
            except Exception as e:
                skipped.append({
                    "clock_folder": clock_folder,
                    "ldsc_dir": ldsc_dir,
                    "log_file": log_file,
                    "reason": "log_filename_parse_failed",
                    "details": str(e),
                })
                continue

            target_info = classify_target(
                target_id=target_id,
                finngen_targets=finngen_targets,
                pgc_targets=pgc_targets,
            )

            try:
                parsed = parse_ldsc_log(log_file)
                error = ""
            except Exception as e:
                parsed = {
                    "gc_mean": np.nan,
                    "gc_std": np.nan,
                    "Z": np.nan,
                    "P": np.nan,
                    "P_log10": np.nan,
                    "neg_log10_P": np.nan,
                    "parse_status": "error",
                }
                error = str(e)

            row = {
                "clock_folder": clock_folder,
                "disease_clock": clock_meta["disease_clock"],
                "disease_key": clock_meta["disease_key"],
                "disease_label": clock_meta["disease_label"],
                "modality": clock_meta["modality"],
                "modality_key": clock_meta["modality_key"],
                "organ_key": clock_meta["organ_key"],
                "organ_label": clock_meta["organ_label"],
                "organ_raw": clock_meta["organ_raw"],
                "clock_prefix": clock_meta["clock_prefix"],

                # This is just the prefix in the log filename.
                # It is not used to infer another organ or another clock.
                "file_exposure_name": exposure_name,

                "target_id": target_id,
                "target_category": target_info["target_category"],
                "target_source": target_info["target_source"],
                "target_family": target_info["target_family"],
                "target_name": target_info["target_name"],
                "target_display": target_info["target_display"],

                "target_mortality_clock_folder": target_info.get("target_mortality_clock_folder", ""),
                "target_mortality_clock_prefix": target_info.get("target_mortality_clock_prefix", ""),
                "target_mortality_modality": target_info.get("target_mortality_modality", ""),
                "target_mortality_organ_key": target_info.get("target_mortality_organ_key", ""),
                "target_mortality_organ_label": target_info.get("target_mortality_organ_label", ""),

                "log_file": log_file,
                "log_relative": os.path.relpath(log_file, base_dir),
                "log_basename": os.path.basename(log_file),
                "error": error,
            }

            row.update(parsed)
            rows.append(row)

    skipped_out = os.path.join(
        output_dir,
        "LDSC_gc_47_disease_epoch_skipped_logs.tsv",
    )

    pd.DataFrame(skipped).to_csv(
        skipped_out,
        index=False,
        sep="\t",
        encoding="utf-8",
    )

    if len(rows) == 0:
        raise RuntimeError(
            "No LDSC genetic-correlation rows were collected. "
            f"Skipped log table written to: {skipped_out}"
        )

    df_all = pd.DataFrame(rows)

    df_all["rg_abs"] = np.abs(df_all["gc_mean"])
    df_all["rg_out_of_bounds"] = (
        np.isfinite(df_all["gc_mean"]) & (df_all["rg_abs"] > 1)
    )

    df_all = add_multiple_testing_columns(df_all)
    df_all = sort_gc_table(df_all)

    save_outputs(
        df_all=df_all,
        df_clock_manifest=df_clock_manifest,
        output_dir=output_dir,
    )

    save_summary_tables(df_all, output_dir)

    print_qc(
        df_all=df_all,
        df_clock_manifest=df_clock_manifest,
        skipped_out=skipped_out,
        expected_n_clocks=expected_n_clocks,
    )

    return df_all


# ============================================================
# 12. QC printout
# ============================================================

def print_qc(df_all, df_clock_manifest, skipped_out, expected_n_clocks):
    print("\n============================================================")
    print("Quick QC")
    print("============================================================")

    print(f"Rows collected:                 {df_all.shape[0]}")
    print(f"Clock folders in manifest:      {df_clock_manifest['clock_folder'].nunique()}")
    print(f"Clock folders with rows:        {df_all['clock_folder'].nunique()}")
    print(f"Expected disease clocks:        {expected_n_clocks}")
    print(f"Unique targets:                 {df_all['target_id'].nunique()}")

    print("\nDiscovered disease clocks by disease:")
    print(
        df_clock_manifest
        .groupby("disease_label")["clock_folder"]
        .nunique()
        .sort_index()
    )

    print("\nRows by disease:")
    print(
        df_all
        .groupby("disease_label")
        .size()
        .sort_index()
    )

    print("\nTarget category counts:")
    print(df_all["target_category"].value_counts(dropna=False))

    print("\nTarget source counts:")
    print(df_all["target_source"].value_counts(dropna=False))

    print("\nParse status:")
    print(df_all["parse_status"].value_counts(dropna=False))

    ok = df_all.query("parse_status == 'ok'").copy()

    if ok.shape[0] > 0:
        print("\nPrimary Bonferroni settings:")
        print("  Scope: unique targets across disease endpoints + AI biomarkers + mortality EPOCH clocks")
        print(f"  n valid pairwise tests: {int(ok.shape[0])}")
        print(f"  n unique targets:       {int(ok['target_id'].nunique())}")
        print(f"  p threshold:            {df_all['primary_bonferroni_p'].iloc[0]:.3e}")
        print(f"  n significant rows:     {int(df_all['primary_bonferroni_sig'].sum())}")

    n_out_bounds = int(df_all["rg_out_of_bounds"].sum())

    if n_out_bounds > 0:
        print(f"\nWARNING: {n_out_bounds} rows have |rg| > 1.")
        print(
            df_all.loc[
                df_all["rg_out_of_bounds"],
                [
                    "clock_folder",
                    "target_id",
                    "target_category",
                    "gc_mean",
                    "gc_std",
                    "P",
                    "log_basename",
                ],
            ].head(30).to_string(index=False)
        )

    if df_all["clock_folder"].nunique() != expected_n_clocks:
        print(
            f"\nWARNING: Expected {expected_n_clocks} clocks, "
            f"but collected rows from {df_all['clock_folder'].nunique()} clocks."
        )

        found = set(df_all["clock_folder"].unique())
        manifest = set(df_clock_manifest["clock_folder"].unique())
        no_results = sorted(manifest - found)

        if len(no_results) > 0:
            print("\nClock folders discovered but no rows collected:")
            for x in no_results:
                print(f"  - {x}")

    print(f"\nSkipped log table: {skipped_out}")


# ============================================================
# 13. Arguments
# ============================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Collect LDSC genetic-correlation results for all 47 disease EPOCH clocks. "
            "Each disease-clock folder is processed independently and is assumed to "
            "contain only its own organ/system-specific LDSC log files."
        )
    )

    parser.add_argument(
        "--base_dir",
        default="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock",
        type=str,
        help="Local WholeBodyClock base directory containing disease EPOCH clock folders.",
    )

    parser.add_argument(
        "--outdir",
        default="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/Result",
        type=str,
        help="Output directory for collected LDSC genetic-correlation TSV files.",
    )

    parser.add_argument(
        "--finngen_manifest",
        default="/cbica/projects/MULTI/processed/FinnGen/GWAS_summary_stats/summary_stats_R9_manifest_5000_cases.tsv",
        type=str,
        help="Optional FinnGen manifest for disease endpoint annotation.",
    )

    parser.add_argument(
        "--pgc_manifest",
        default="/cbica/projects/MULTI/processed/PGC/PGC_MUTATE.tsv",
        type=str,
        help="Optional PGC manifest for disease endpoint annotation.",
    )

    parser.add_argument(
        "--expected_n_clocks",
        default=47,
        type=int,
        help="Expected number of disease EPOCH clocks.",
    )

    return parser.parse_args()


# ============================================================
# 14. Run
# ============================================================

if __name__ == "__main__":
    args = parse_args()

    collect_all_gc_results(
        base_dir=args.base_dir,
        output_dir=args.outdir,
        finngen_manifest=args.finngen_manifest,
        pgc_manifest=args.pgc_manifest,
        expected_n_clocks=args.expected_n_clocks,
    )