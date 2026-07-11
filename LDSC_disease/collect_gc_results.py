import os
import re
import glob
import numpy as np
import pandas as pd


# ============================================================
# 1. Hardcoded paths
# ============================================================

# Local Mac path after rsync
LDSC_ROOT_DIR = (
    "/Users/hao/cubic-home/Reproducibile_paper/"
    "WholeBodyClock/mortality_clock/fastGWA/output"
)

OUTPUT_DIR_RESULT = (
    "/Users/hao/cubic-home/Reproducibile_paper/"
    "WholeBodyClock/Result"
)

# Optional manifest files for better FinnGen / PGC annotation.
# If these files do not exist locally, the script still runs.
FINNGEN_MANIFEST = (
    "/Users/hao/cubic-projects/MULTI/processed/FinnGen/"
    "GWAS_summary_stats/summary_stats_R9_manifest_5000_cases.tsv"
)

PGC_MANIFEST = (
    "/Users/hao/cubic-projects/MULTI/processed/PGC/"
    "PGC_MUTATE.tsv"
)

# If running directly on CUBIC, use:
# LDSC_ROOT_DIR = (
#     "/cbica/home/wenju/Reproducibile_paper/"
#     "WholeBodyClock/mortality_clock/fastGWA/output"
# )
# OUTPUT_DIR_RESULT = (
#     "/cbica/home/wenju/Reproducibile_paper/"
#     "WholeBodyClock/Result"
# )
# FINNGEN_MANIFEST = (
#     "/cbica/projects/MULTI/processed/FinnGen/"
#     "GWAS_summary_stats/summary_stats_R9_manifest_5000_cases.tsv"
# )
# PGC_MANIFEST = (
#     "/cbica/projects/MULTI/processed/PGC/"
#     "PGC_MUTATE.tsv"
# )


# ============================================================
# 2. Exact AI biomarker list from your LDSC bash script
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

# Pan-disease MAEs from Multiorgan_Subtype.
MAE_BRAIN_LIST = ["r1", "r2", "r3", "r4", "r5", "r6"]
MAE_EYE_LIST = ["r1", "r2", "r3"]
MAE_HEART_LIST = ["r1", "r2"]

# DNE biomarkers.
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

# Multi-modal brain BAGs.
# User requested:
#   muse              -> GM brain MRIBAG
#   dmri_fullmetric_fa -> WM brain MRIBAG
#   fmri              -> FC brain MRIBAG
MULTIMODAL_BRAIN_BAG_MAP = {
    "muse": "GM brain MRIBAG",
    "dmri_fullmetric_fa": "WM brain MRIBAG",
    "fmri": "FC brain MRIBAG",
}


def build_exact_ai_biomarker_table():
    """
    Build an exact lookup table for AI biomarkers used in LDSC.
    The target_id must match the suffix after mortality_clock_vs_ in the .log file.
    """

    rows = []

    # 11 ProtBAGs
    for bag in PROTBAG_LIST:
        rows.append({
            "target_id": f"ProtBAG_{bag}",
            "analysis_group": "AI_biomarker",
            "target_source": "ProtBAG",
            "target_family": "Proteomics aging clock",
            "target_name": bag,
            "target_display": f"ProtBAG {bag}",
        })

    # 5 MetBAGs
    for bag in METBAG_LIST:
        rows.append({
            "target_id": f"MetBAG_{bag}",
            "analysis_group": "AI_biomarker",
            "target_source": "MetBAG",
            "target_family": "Metabolomics aging clock",
            "target_name": bag,
            "target_display": f"MetBAG {bag}",
        })

    # 7 MRIBAGs
    for bag in MRIBAG_LIST:
        rows.append({
            "target_id": f"MRIBAG_{bag}",
            "analysis_group": "AI_biomarker",
            "target_source": "MRIBAG",
            "target_family": "Organ MRI aging clock",
            "target_name": bag,
            "target_display": f"MRIBAG {bag}",
        })

    # Pan-disease MAEs: brain
    for mae in MAE_BRAIN_LIST:
        rows.append({
            "target_id": f"{mae}_brain",
            "analysis_group": "AI_biomarker",
            "target_source": "Pan_disease_MAE",
            "target_family": "Brain pan-disease MAE",
            "target_name": mae,
            "target_display": f"{mae.upper()} brain MAE",
        })

    # Pan-disease MAEs: eye
    for mae in MAE_EYE_LIST:
        rows.append({
            "target_id": f"{mae}_eye",
            "analysis_group": "AI_biomarker",
            "target_source": "Pan_disease_MAE",
            "target_family": "Eye pan-disease MAE",
            "target_name": mae,
            "target_display": f"{mae.upper()} eye MAE",
        })

    # Pan-disease MAEs: heart
    for mae in MAE_HEART_LIST:
        rows.append({
            "target_id": f"{mae}_heart",
            "analysis_group": "AI_biomarker",
            "target_source": "Pan_disease_MAE",
            "target_family": "Heart pan-disease MAE",
            "target_name": mae,
            "target_display": f"{mae.upper()} heart MAE",
        })

    # DNE biomarkers
    for mae in DNE_LIST:
        if mae.startswith("AD_SurrealGAN"):
            family = "AD DNE subtype"
        elif mae.startswith("ASD"):
            family = "ASD DNE subtype"
        elif mae.startswith("LLD"):
            family = "LLD DNE subtype"
        elif mae.startswith("SCZ"):
            family = "SCZ DNE subtype"
        else:
            family = "DNE subtype"

        rows.append({
            "target_id": mae,
            "analysis_group": "AI_biomarker",
            "target_source": "DNE",
            "target_family": family,
            "target_name": mae,
            "target_display": mae.replace("_", " "),
        })

    # 3 multi-modal brain BAGs
    for target_id, display in MULTIMODAL_BRAIN_BAG_MAP.items():
        rows.append({
            "target_id": target_id,
            "analysis_group": "AI_biomarker",
            "target_source": "Brain_MRIBAG",
            "target_family": "Multi-modal brain MRI aging clock",
            "target_name": target_id,
            "target_display": display,
        })

    df_ai = pd.DataFrame(rows)
    return df_ai


AI_BIOMARKER_DF = build_exact_ai_biomarker_table()
AI_BIOMARKER_LOOKUP = AI_BIOMARKER_DF.set_index("target_id").to_dict(orient="index")


# ============================================================
# 3. Optional FinnGen / PGC manifest loading
# ============================================================

def load_finngen_targets(manifest_path):
    """
    FinnGen manifest in bash:
      summary_tsv='...summary_stats_R9_manifest_5000_cases.tsv'
      variable=`awk '{print $1}' ${summary_tsv}`
      if de != phenocode
    """
    if not os.path.exists(manifest_path):
        print(f"WARNING: FinnGen manifest not found: {manifest_path}")
        return set()

    try:
        df = pd.read_csv(manifest_path, sep="\t")
        if "phenocode" in df.columns:
            return set(df["phenocode"].dropna().astype(str).tolist())

        # Fallback: first column
        first_col = df.columns[0]
        vals = df[first_col].dropna().astype(str).tolist()
        vals = [x for x in vals if x != "phenocode"]
        return set(vals)

    except Exception as e:
        print(f"WARNING: Could not read FinnGen manifest: {e}")
        return set()


def load_pgc_targets(manifest_path):
    """
    PGC manifest in bash:
      summary_tsv='...PGC_MUTATE.tsv'
      variable=`awk -F'\t' '{print $2}' ${summary_tsv}`
      if de != Phenotype
    """
    if not os.path.exists(manifest_path):
        print(f"WARNING: PGC manifest not found: {manifest_path}")
        return set()

    try:
        df = pd.read_csv(manifest_path, sep="\t")

        if "Phenotype" in df.columns:
            return set(df["Phenotype"].dropna().astype(str).tolist())

        # Fallback: second column
        if df.shape[1] >= 2:
            vals = df.iloc[:, 1].dropna().astype(str).tolist()
            vals = [x for x in vals if x != "Phenotype"]
            return set(vals)

        return set()

    except Exception as e:
        print(f"WARNING: Could not read PGC manifest: {e}")
        return set()


FINNGEN_TARGETS = load_finngen_targets(FINNGEN_MANIFEST)
PGC_TARGETS = load_pgc_targets(PGC_MANIFEST)

# Fallback PGC labels in case manifest is not available locally.
PGC_FALLBACK_TARGETS = {
    "AD",
    "ADHD",
    "ASD",
    "ASD_1",
    "ASD_2",
    "ASD_3",
    "BIP",
    "SCZ",
    "SCZ_1",
    "SCZ_2",
    "OCD",
    "AUD",
    "AN",
    "MDD",
    "PTSD",
    "TS",
}


# ============================================================
# 4. Basic helpers
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


def infer_modality_from_clock_name(clock_name):
    x = str(clock_name).lower()

    if x.endswith("_mri"):
        return "MRI"
    if x.endswith("_proteomics"):
        return "Proteomics"
    if x.endswith("_metabolomics"):
        return "Metabolomics"

    return "Unknown"


def get_clock_name_from_folder(clock_folder):
    """
    Examples:
      spleen_mri_mortality_clock -> spleen_mri
      Brain_proteomics_mortality_clock -> Brain_proteomics
      Endocrine_metabolomics_mortality_clock -> Endocrine_metabolomics
    """
    return clock_folder.replace("_mortality_clock", "")


def parse_log_filename(log_file, clock_name_from_folder):
    """
    Expected basename examples:
      spleen_mri_vs_MRIBAG_adipose.log
      spleen_mri_vs_ProtBAG_Brain.log
      spleen_mri_vs_T2D.log
      spleen_mri_h2.log should already be excluded.

    Returns:
      exposure_name = spleen_mri
      target_id     = MRIBAG_adipose
    """
    basename = os.path.basename(log_file)
    basename_no_ext = re.sub(r"\.log$", "", basename)

    if "_vs_" in basename_no_ext:
        exposure_name, target_id = basename_no_ext.split("_vs_", 1)
    else:
        exposure_name = clock_name_from_folder
        target_id = basename_no_ext

    return exposure_name, target_id


def classify_target(target_id):
    """
    Exact AI biomarker classification from your bash script.
    Everything not in exact AI list is treated as FinnGen or PGC disease endpoint.
    """

    target = str(target_id)

    # Exact AI biomarkers
    if target in AI_BIOMARKER_LOOKUP:
        return AI_BIOMARKER_LOOKUP[target].copy()

    # PGC endpoints
    if target in PGC_TARGETS or target in PGC_FALLBACK_TARGETS:
        return {
            "analysis_group": "Disease_endpoint",
            "target_source": "PGC",
            "target_family": "PGC",
            "target_name": target,
            "target_display": target,
        }

    # FinnGen endpoints
    if target in FINNGEN_TARGETS:
        return {
            "analysis_group": "Disease_endpoint",
            "target_source": "FinnGen",
            "target_family": "FinnGen",
            "target_name": target,
            "target_display": target,
        }

    # Fallback disease endpoint.
    # This covers FinnGen if the manifest is unavailable locally,
    # and any other disease endpoint logs.
    return {
        "analysis_group": "Disease_endpoint",
        "target_source": "FinnGen_or_other",
        "target_family": "Disease endpoint",
        "target_name": target,
        "target_display": target,
    }


# ============================================================
# 5. LDSC log parser
# ============================================================

def parse_ldsc_log(log_file):
    """
    Parse LDSC genetic correlation result from a .log file.

    Supports:
      Format 1:
        Genetic Correlation: 0.123 (0.045)
        Z-score: 2.73
        P: 0.0063

      Format 2, LDSC table fallback:
        p1 p2 rg se z p ...
    """
    with open(log_file, "r", errors="ignore") as f:
        txt = f.read()

    gc_mean = np.nan
    gc_std = np.nan
    z_score = np.nan
    p_value = np.nan

    # Format 1: explicit LDSC output lines.
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

    # Format 2: LDSC table fallback.
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
# 6. Collect LDSC logs
# ============================================================

def collect_ldsc_gc_logs(ldsc_root_dir):
    """
    Collect pairwise genetic-correlation logs only.

    Important:
      This uses *_vs_*.log so h2 logs such as spleen_mri_h2.log are excluded.
    """
    pattern = os.path.join(
        ldsc_root_dir,
        "*_mortality_clock",
        "ldsc",
        "*_vs_*.log",
    )

    return sorted(glob.glob(pattern))


def collect_all_gc_results(ldsc_root_dir, output_dir_result):
    """
    Collect all pairwise LDSC genetic-correlation results for mortality clocks.
    """
    os.makedirs(output_dir_result, exist_ok=True)

    log_files = collect_ldsc_gc_logs(ldsc_root_dir)

    print(f"Found {len(log_files)} pairwise LDSC genetic-correlation logs.")

    if len(log_files) == 0:
        raise RuntimeError(f"No pairwise LDSC logs found under: {ldsc_root_dir}")

    rows = []

    for log_file in log_files:
        clock_folder = os.path.basename(os.path.dirname(os.path.dirname(log_file)))
        mortality_clock = get_clock_name_from_folder(clock_folder)
        modality = infer_modality_from_clock_name(mortality_clock)

        mortality_clock_from_file, target_id = parse_log_filename(
            log_file=log_file,
            clock_name_from_folder=mortality_clock,
        )

        target_info = classify_target(target_id)

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
            "mortality_clock": mortality_clock,
            "mortality_clock_from_file": mortality_clock_from_file,
            "modality": modality,
            "target_id": target_id,
            "analysis_group": target_info["analysis_group"],
            "target_source": target_info["target_source"],
            "target_family": target_info["target_family"],
            "target_name": target_info["target_name"],
            "target_display": target_info["target_display"],
            "log_file": log_file,
            "log_basename": os.path.basename(log_file),
            "error": error,
        }

        row.update(parsed)
        rows.append(row)

    df_all = pd.DataFrame(rows)

    df_all["clock_name_mismatch"] = (
        df_all["mortality_clock"] != df_all["mortality_clock_from_file"]
    )

    df_all["rg_abs"] = np.abs(df_all["gc_mean"])
    df_all["rg_out_of_bounds"] = (
        np.isfinite(df_all["gc_mean"]) & (df_all["gc_mean"].abs() > 1)
    )

    df_all = add_bonferroni_columns(df_all)
    df_all = sort_gc_table(df_all)

    save_gc_outputs(df_all, output_dir_result)
    save_gc_summary_tables(df_all, output_dir_result)
    save_ai_biomarker_reference_table(output_dir_result)

    print_qc(df_all)

    return df_all


# ============================================================
# 7. Multiple-testing correction and sorting
# ============================================================

def add_bonferroni_columns(df_all):
    """
    Add Bonferroni thresholds at several levels:
      1. analysis_group: AI_biomarker vs Disease_endpoint
      2. target_source: MRIBAG, ProtBAG, MetBAG, Brain_MRIBAG, DNE, PGC, FinnGen
      3. mortality_clock x analysis_group
      4. mortality_clock x target_source
    """
    df = df_all.copy()

    df["n_tests_analysis_group"] = np.nan
    df["bonferroni_p_analysis_group"] = np.nan
    df["bonferroni_sig_analysis_group"] = False

    for group, idx in df.groupby("analysis_group").groups.items():
        idx = list(idx)
        n_tests = df.loc[idx].query("parse_status == 'ok'").shape[0]
        bonf = 0.05 / n_tests if n_tests > 0 else np.nan

        df.loc[idx, "n_tests_analysis_group"] = n_tests
        df.loc[idx, "bonferroni_p_analysis_group"] = bonf
        df.loc[idx, "bonferroni_sig_analysis_group"] = df.loc[idx, "P"] < bonf

    df["n_tests_target_source"] = np.nan
    df["bonferroni_p_target_source"] = np.nan
    df["bonferroni_sig_target_source"] = False

    for source, idx in df.groupby("target_source").groups.items():
        idx = list(idx)
        n_tests = df.loc[idx].query("parse_status == 'ok'").shape[0]
        bonf = 0.05 / n_tests if n_tests > 0 else np.nan

        df.loc[idx, "n_tests_target_source"] = n_tests
        df.loc[idx, "bonferroni_p_target_source"] = bonf
        df.loc[idx, "bonferroni_sig_target_source"] = df.loc[idx, "P"] < bonf

    df["n_tests_clock_analysis_group"] = np.nan
    df["bonferroni_p_clock_analysis_group"] = np.nan
    df["bonferroni_sig_clock_analysis_group"] = False

    for keys, idx in df.groupby(["mortality_clock", "analysis_group"]).groups.items():
        idx = list(idx)
        n_tests = df.loc[idx].query("parse_status == 'ok'").shape[0]
        bonf = 0.05 / n_tests if n_tests > 0 else np.nan

        df.loc[idx, "n_tests_clock_analysis_group"] = n_tests
        df.loc[idx, "bonferroni_p_clock_analysis_group"] = bonf
        df.loc[idx, "bonferroni_sig_clock_analysis_group"] = df.loc[idx, "P"] < bonf

    df["n_tests_clock_target_source"] = np.nan
    df["bonferroni_p_clock_target_source"] = np.nan
    df["bonferroni_sig_clock_target_source"] = False

    for keys, idx in df.groupby(["mortality_clock", "target_source"]).groups.items():
        idx = list(idx)
        n_tests = df.loc[idx].query("parse_status == 'ok'").shape[0]
        bonf = 0.05 / n_tests if n_tests > 0 else np.nan

        df.loc[idx, "n_tests_clock_target_source"] = n_tests
        df.loc[idx, "bonferroni_p_clock_target_source"] = bonf
        df.loc[idx, "bonferroni_sig_clock_target_source"] = df.loc[idx, "P"] < bonf

    df["nominal_sig"] = df["P"] < 0.05

    return df


def sort_gc_table(df_all):
    df = df_all.copy()

    modality_order = {
        "MRI": 1,
        "Proteomics": 2,
        "Metabolomics": 3,
        "Unknown": 99,
    }

    analysis_order = {
        "AI_biomarker": 1,
        "Disease_endpoint": 2,
    }

    target_source_order = {
        "MRIBAG": 1,
        "ProtBAG": 2,
        "MetBAG": 3,
        "Brain_MRIBAG": 4,
        "Pan_disease_MAE": 5,
        "DNE": 6,
        "PGC": 7,
        "FinnGen": 8,
        "FinnGen_or_other": 9,
    }

    df["modality_order"] = df["modality"].map(modality_order).fillna(99)
    df["analysis_order"] = df["analysis_group"].map(analysis_order).fillna(99)
    df["target_source_order"] = df["target_source"].map(target_source_order).fillna(99)

    df = (
        df.sort_values(
            [
                "analysis_order",
                "target_source_order",
                "modality_order",
                "mortality_clock",
                "target_id",
            ]
        )
        .drop(columns=["analysis_order", "target_source_order", "modality_order"])
        .reset_index(drop=True)
    )

    return df


# ============================================================
# 8. Save outputs
# ============================================================

def save_gc_outputs(df_all, output_dir_result):
    out_all = os.path.join(
        output_dir_result,
        "LDSC_gc_mortality_clocks_all_targets.tsv",
    )
    df_all.to_csv(out_all, index=False, sep="\t", encoding="utf-8")
    print(f"Saved all results: {out_all}")

    df_ai = df_all.loc[df_all["analysis_group"] == "AI_biomarker"].copy()
    df_disease = df_all.loc[df_all["analysis_group"] == "Disease_endpoint"].copy()
    df_pgc = df_disease.loc[df_disease["target_source"] == "PGC"].copy()
    df_finngen = df_disease.loc[df_disease["target_source"].isin(["FinnGen", "FinnGen_or_other"])].copy()
    df_failed = df_all.loc[df_all["parse_status"] != "ok"].copy()

    out_ai = os.path.join(
        output_dir_result,
        "LDSC_gc_mortality_clocks_vs_exact_AI_biomarkers.tsv",
    )
    out_disease = os.path.join(
        output_dir_result,
        "LDSC_gc_mortality_clocks_vs_disease_endpoints.tsv",
    )
    out_pgc = os.path.join(
        output_dir_result,
        "LDSC_gc_mortality_clocks_vs_PGC_endpoints.tsv",
    )
    out_finngen = os.path.join(
        output_dir_result,
        "LDSC_gc_mortality_clocks_vs_FinnGen_or_other_endpoints.tsv",
    )
    out_failed = os.path.join(
        output_dir_result,
        "LDSC_gc_parse_failed_or_nan.tsv",
    )

    df_ai.to_csv(out_ai, index=False, sep="\t", encoding="utf-8")
    df_disease.to_csv(out_disease, index=False, sep="\t", encoding="utf-8")
    df_pgc.to_csv(out_pgc, index=False, sep="\t", encoding="utf-8")
    df_finngen.to_csv(out_finngen, index=False, sep="\t", encoding="utf-8")
    df_failed.to_csv(out_failed, index=False, sep="\t", encoding="utf-8")

    print(f"Saved exact AI biomarker results: {out_ai}")
    print(f"Saved disease endpoint results: {out_disease}")
    print(f"Saved PGC endpoint results: {out_pgc}")
    print(f"Saved FinnGen/other endpoint results: {out_finngen}")
    print(f"Saved failed/nan parse results: {out_failed}")


def save_ai_biomarker_reference_table(output_dir_result):
    out_ref = os.path.join(
        output_dir_result,
        "LDSC_gc_exact_AI_biomarker_reference_table.tsv",
    )
    AI_BIOMARKER_DF.to_csv(out_ref, index=False, sep="\t", encoding="utf-8")
    print(f"Saved exact AI biomarker reference table: {out_ref}")


# ============================================================
# 9. Summary tables for plotting/QC
# ============================================================

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


def save_gc_summary_tables(df_all, output_dir_result):
    df_ok = df_all.loc[df_all["parse_status"] == "ok"].copy()

    summary_clock = (
        df_ok.groupby(
            ["mortality_clock", "modality", "analysis_group", "target_source"],
            as_index=False,
        )
        .agg(
            n_tests=("target_id", "count"),
            n_nominal_sig=("nominal_sig", "sum"),
            n_bonf_sig_analysis_group=("bonferroni_sig_analysis_group", "sum"),
            n_bonf_sig_target_source=("bonferroni_sig_target_source", "sum"),
            n_bonf_sig_clock_analysis_group=("bonferroni_sig_clock_analysis_group", "sum"),
            n_bonf_sig_clock_target_source=("bonferroni_sig_clock_target_source", "sum"),
            mean_rg=("gc_mean", "mean"),
            mean_abs_rg=("gc_mean", nanmean_abs),
            max_abs_rg=("gc_mean", nanmax_abs),
            min_p=("P", "min"),
            max_neg_log10_p=("neg_log10_P", "max"),
        )
        .sort_values(
            ["analysis_group", "target_source", "modality", "mortality_clock"]
        )
    )

    out_summary_clock = os.path.join(
        output_dir_result,
        "LDSC_gc_summary_by_clock_and_target_source.tsv",
    )
    summary_clock.to_csv(out_summary_clock, index=False, sep="\t", encoding="utf-8")
    print(f"Saved clock-level summary: {out_summary_clock}")

    summary_target = (
        df_ok.groupby(
            [
                "target_id",
                "target_display",
                "target_name",
                "analysis_group",
                "target_source",
                "target_family",
            ],
            as_index=False,
        )
        .agg(
            n_clocks=("mortality_clock", "count"),
            n_nominal_sig=("nominal_sig", "sum"),
            n_bonf_sig_analysis_group=("bonferroni_sig_analysis_group", "sum"),
            n_bonf_sig_target_source=("bonferroni_sig_target_source", "sum"),
            mean_rg=("gc_mean", "mean"),
            mean_abs_rg=("gc_mean", nanmean_abs),
            max_abs_rg=("gc_mean", nanmax_abs),
            min_p=("P", "min"),
            max_neg_log10_p=("neg_log10_P", "max"),
        )
        .sort_values(
            ["analysis_group", "target_source", "target_family", "min_p", "target_id"]
        )
    )

    out_summary_target = os.path.join(
        output_dir_result,
        "LDSC_gc_summary_by_target.tsv",
    )
    summary_target.to_csv(out_summary_target, index=False, sep="\t", encoding="utf-8")
    print(f"Saved target-level summary: {out_summary_target}")

    top_hits = (
        df_ok.sort_values("P")
        .loc[
            :,
            [
                "mortality_clock",
                "modality",
                "target_id",
                "target_display",
                "target_name",
                "analysis_group",
                "target_source",
                "target_family",
                "gc_mean",
                "gc_std",
                "Z",
                "P",
                "neg_log10_P",
                "nominal_sig",
                "bonferroni_sig_analysis_group",
                "bonferroni_sig_target_source",
                "bonferroni_sig_clock_analysis_group",
                "bonferroni_sig_clock_target_source",
                "rg_out_of_bounds",
                "log_file",
            ],
        ]
    )

    out_top_hits = os.path.join(
        output_dir_result,
        "LDSC_gc_top_hits_all_targets.tsv",
    )
    top_hits.to_csv(out_top_hits, index=False, sep="\t", encoding="utf-8")
    print(f"Saved top-hit table: {out_top_hits}")


# ============================================================
# 10. QC printing
# ============================================================

def print_qc(df_all):
    print("\nQuick QC:")
    print(f"Total logs collected: {df_all.shape[0]}")
    print(f"Unique mortality clocks: {df_all['mortality_clock'].nunique()}")

    print("\nAnalysis group counts:")
    print(df_all["analysis_group"].value_counts(dropna=False))

    print("\nTarget source counts:")
    print(df_all["target_source"].value_counts(dropna=False))

    print("\nTarget family counts for AI biomarkers:")
    print(
        df_all.loc[df_all["analysis_group"] == "AI_biomarker", "target_family"]
        .value_counts(dropna=False)
    )

    print("\nParse status:")
    print(df_all["parse_status"].value_counts(dropna=False))

    n_mismatch = int(df_all["clock_name_mismatch"].sum())
    if n_mismatch > 0:
        print(f"\nWARNING: {n_mismatch} rows have clock-name mismatch.")
        print(
            df_all.loc[
                df_all["clock_name_mismatch"],
                ["mortality_clock", "mortality_clock_from_file", "log_basename"],
            ].head(30)
        )

    n_out_bounds = int(df_all["rg_out_of_bounds"].sum())
    if n_out_bounds > 0:
        print(f"\nWARNING: {n_out_bounds} rows have |rg| > 1.")
        print(
            df_all.loc[
                df_all["rg_out_of_bounds"],
                ["mortality_clock", "target_id", "gc_mean", "gc_std", "P", "log_basename"],
            ].head(30)
        )


# ============================================================
# 11. Run
# ============================================================

if __name__ == "__main__":

    df_gc = collect_all_gc_results(
        ldsc_root_dir=LDSC_ROOT_DIR,
        output_dir_result=OUTPUT_DIR_RESULT,
    )