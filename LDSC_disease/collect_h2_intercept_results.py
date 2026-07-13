#!/usr/bin/env python3

# ============================================================
# Collect LDSC h2/intercept statistics for 47 disease EPOCH clocks
#
# Disease endpoints:
#   dementia, copd, asthma, mi, stroke
#
# Expected input:
#   <BASE_DIR>/*_clock/fastGWA/ldsc/*_h2.log
#
# Example:
#   /cbica/home/wenju/Reproducibile_paper/WholeBodyClock/
#     Brain_proteomics_dementia_clock/fastGWA/ldsc/Brain_h2.log
#
# Outputs:
#   Result/LDSC_h2_intercept_47_disease_clocks.tsv
#   Result/LDSC_h2_intercept_47_disease_clocks.xlsx
#   Result/LDSC_h2_intercept_47_disease_clocks_ok_only.tsv
#   Result/LDSC_h2_intercept_47_disease_clocks_failed_or_nan.tsv
#   Result/LDSC_h2_intercept_47_disease_clocks_summary_by_disease_modality.tsv
#   Result/LDSC_h2_intercept_47_disease_clocks_QC_flags.tsv
# ============================================================

import os
import re
import glob
import math
import numpy as np
import pandas as pd


# ============================================================
# 1. Paths
# ============================================================

BASE_DIR_CANDIDATES = [
    "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock",
    "/gpfs/fs001/cbica/home/wenju/Reproducibile_paper/WholeBodyClock",
    "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock",
    os.getcwd(),
]

BASE_DIR = os.environ.get("WHOLEBODYCLOCK_BASE_DIR", "")

if not BASE_DIR:
    existing = [x for x in BASE_DIR_CANDIDATES if os.path.isdir(x)]
    if not existing:
        raise RuntimeError(
            "Could not find BASE_DIR. Please set WHOLEBODYCLOCK_BASE_DIR "
            "or edit BASE_DIR manually."
        )
    BASE_DIR = existing[0]

BASE_DIR = os.path.abspath(BASE_DIR)

OUTPUT_DIR_RESULT = os.environ.get(
    "WHOLEBODYCLOCK_RESULT_DIR",
    os.path.join(BASE_DIR, "Result")
)

os.makedirs(OUTPUT_DIR_RESULT, exist_ok=True)

EXPECTED_N_LOGS = 47

DISEASE_ORDER = ["dementia", "copd", "asthma", "mi", "stroke"]

DISEASE_LABELS = {
    "dementia": "Dementia",
    "copd": "COPD",
    "asthma": "Asthma",
    "mi": "MI",
    "stroke": "Stroke",
}

MODALITY_ORDER = {
    "MRI": 1,
    "Proteomics": 2,
    "Metabolomics": 3,
    "Unknown": 99,
}

DISEASE_ORDER_MAP = {x: i + 1 for i, x in enumerate(DISEASE_ORDER)}


# ============================================================
# 2. Helper functions
# ============================================================

def safe_float(x, default=np.nan):
    try:
        if x is None:
            return default

        s = str(x).strip()

        if s.lower() in ["", "nan", "na", "none", "null"]:
            return default

        if s.lower() in ["inf", "+inf", "infinity", "+infinity"]:
            return np.inf

        if s.lower() in ["-inf", "-infinity"]:
            return -np.inf

        return float(s)

    except Exception:
        return default


def p_to_display(p):
    p = safe_float(p)

    if not np.isfinite(p):
        return ""

    if p == 0:
        return "<1e-300"

    if p < 0.001:
        return f"{p:.2e}"

    return f"{p:.3g}"


def infer_modality_from_folder(clock_folder):
    x = str(clock_folder).lower()

    if "_mri_" in x:
        return "MRI"

    if "_proteomics_" in x:
        return "Proteomics"

    if "_metabolomics_" in x:
        return "Metabolomics"

    return "Unknown"


def infer_metadata_from_folder(clock_folder):
    """
    Expected examples:
      Brain_proteomics_dementia_clock
      heart_mri_copd_clock
      Digestive_metabolomics_mi_clock
      Reproductive_female_proteomics_stroke_clock
    """

    pattern = (
        r"^(.+)_(mri|proteomics|metabolomics)_"
        r"(asthma|copd|dementia|mi|stroke)_clock$"
    )

    m = re.match(pattern, clock_folder)

    if m is None:
        modality = infer_modality_from_folder(clock_folder)

        return {
            "organ_raw": "",
            "organ_key": "",
            "organ_label": "",
            "modality_key": "",
            "modality": modality,
            "disease_key": "",
            "disease_label": "",
            "disease_clock": "",
        }

    organ_raw = m.group(1)
    modality_key = m.group(2)
    disease_key = m.group(3)

    organ_key = organ_raw.lower()
    organ_label = organ_raw.replace("_", " ").capitalize()

    modality = {
        "mri": "MRI",
        "proteomics": "Proteomics",
        "metabolomics": "Metabolomics",
    }.get(modality_key, "Unknown")

    disease_label = DISEASE_LABELS.get(disease_key, disease_key)

    disease_clock = f"{disease_key}__{organ_key}__{modality_key}"

    return {
        "organ_raw": organ_raw,
        "organ_key": organ_key,
        "organ_label": organ_label,
        "modality_key": modality_key,
        "modality": modality,
        "disease_key": disease_key,
        "disease_label": disease_label,
        "disease_clock": disease_clock,
    }


def parse_trait_from_h2_log(h2_log):
    """
    Example:
      Brain_h2.log -> brain
      Reproductive_h2.log -> reproductive
      spleen_h2.log -> spleen
    """

    basename = os.path.basename(h2_log)
    return re.sub(r"_h2\.log$", "", basename).lower()


def is_file_trait_compatible(organ_key, file_trait_key):
    if not organ_key or not file_trait_key:
        return np.nan

    # Reproductive female/male logs are often named Reproductive_h2.log.
    if organ_key.startswith("reproductive_") and file_trait_key == "reproductive":
        return True

    return organ_key == file_trait_key


def extract_value_and_se(text, label):
    """
    Extract LDSC lines such as:
      Total Observed scale h2: 0.1234 (0.0456)
      Intercept: 1.0123 (0.0078)
      Ratio: 0.0345 (0.0123)

    Returns:
      value, se, found
    """

    num_pattern = r"([+-]?(?:nan|-?inf|[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?)"

    pattern = (
        rf"{re.escape(label)}\s*:\s*"
        rf"{num_pattern}"
        rf"(?:\s*\({num_pattern}\))?"
    )

    m = re.search(pattern, text, flags=re.IGNORECASE)

    if m is None:
        return np.nan, np.nan, False

    value = safe_float(m.group(1))
    se = safe_float(m.group(2)) if m.lastindex and m.lastindex >= 2 else np.nan

    return value, se, True


def extract_single_value(text, label):
    """
    Extract LDSC lines such as:
      Mean Chi^2: 1.021
      Lambda GC: 1.012
      h2 Z: 3.21
    """

    num_pattern = r"([+-]?(?:nan|-?inf|[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?)"

    pattern = (
        rf"{re.escape(label)}\s*:\s*"
        rf"{num_pattern}"
    )

    m = re.search(pattern, text, flags=re.IGNORECASE)

    if m is None:
        return np.nan, False

    return safe_float(m.group(1)), True


# ============================================================
# 3. LDSC h2 parser
# ============================================================

def parse_ldsc_h2_log(h2_log):
    """
    Parse one LDSC h2 log file.

    Expected key lines:
      Total Observed scale h2: ...
      Intercept: ...
      Ratio: ...
      Lambda GC: ...
      Mean Chi^2: ...

    Some LDSC logs may show non-estimable ratio lines.
    """

    with open(h2_log, "r", errors="ignore") as f:
        txt = f.read()

    h2_mean, h2_std, h2_found = extract_value_and_se(
        txt,
        "Total Observed scale h2"
    )

    intercept, intercept_std, intercept_found = extract_value_and_se(
        txt,
        "Intercept"
    )

    ratio, ratio_std, ratio_found = extract_value_and_se(
        txt,
        "Ratio"
    )

    if not ratio_found:
        ratio_text_present = re.search(
            r"Ratio\s+<|Ratio\s+not|Ratio\s+undefined|Ratio\s+NA",
            txt,
            flags=re.IGNORECASE,
        )
        if ratio_text_present:
            ratio = np.nan
            ratio_std = np.nan
            ratio_found = True

    lambda_gc, lambda_found = extract_single_value(
        txt,
        "Lambda GC"
    )

    mean_chi2, mean_chi2_found = extract_single_value(
        txt,
        "Mean Chi^2"
    )

    h2_z_from_log, h2_z_found = extract_single_value(
        txt,
        "h2 Z"
    )

    if np.isfinite(h2_mean) and np.isfinite(h2_std) and h2_std > 0:
        h2_z_calc = h2_mean / h2_std
        h2_p_calc = math.erfc(abs(h2_z_calc) / math.sqrt(2.0))
    else:
        h2_z_calc = np.nan
        h2_p_calc = np.nan

    required_found = h2_found and intercept_found
    required_finite = np.isfinite(h2_mean) and np.isfinite(intercept)

    txt_lower = txt.lower()

    if required_found and required_finite:
        parse_status = "ok"
    elif "nan" in txt_lower:
        parse_status = "nan_result"
    elif "error" in txt_lower or "traceback" in txt_lower:
        parse_status = "ldsc_error"
    else:
        parse_status = "parse_failed"

    return {
        "h2_mean": h2_mean,
        "h2_std": h2_std,
        "h2_z_from_log": h2_z_from_log,
        "h2_z_calc": h2_z_calc,
        "h2_p_calc": h2_p_calc,
        "h2_p_calc_display": p_to_display(h2_p_calc),

        "intercept": intercept,
        "intercept_std": intercept_std,

        "ratio": ratio,
        "ratio_std": ratio_std,

        "lambda_gc": lambda_gc,
        "mean_chi2": mean_chi2,

        "h2_found": h2_found,
        "intercept_found": intercept_found,
        "ratio_found": ratio_found,
        "lambda_gc_found": lambda_found,
        "mean_chi2_found": mean_chi2_found,
        "h2_z_found": h2_z_found,

        "parse_status": parse_status,
    }


# ============================================================
# 4. Collect disease h2 logs
# ============================================================

def collect_disease_ldsc_h2_logs(base_dir):
    """
    Restricted glob only:
      <base_dir>/*_clock/fastGWA/ldsc/*_h2.log

    Then keep only:
      *_asthma_clock
      *_copd_clock
      *_dementia_clock
      *_mi_clock
      *_stroke_clock
    """

    pattern = os.path.join(
        base_dir,
        "*_clock",
        "fastGWA",
        "ldsc",
        "*_h2.log",
    )

    h2_logs = sorted(glob.glob(pattern))

    disease_pattern = re.compile(
        r"_(asthma|copd|dementia|mi|stroke)_clock[/\\]fastGWA[/\\]ldsc[/\\][^/\\]+_h2\.log$",
        flags=re.IGNORECASE,
    )

    h2_logs = [
        os.path.abspath(x)
        for x in h2_logs
        if disease_pattern.search(x)
    ]

    return sorted(set(h2_logs))


# ============================================================
# 5. Sorting and summary helpers
# ============================================================

def sort_disease_h2_table(df):
    df = df.copy()

    df["disease_order"] = df["disease_key"].map(DISEASE_ORDER_MAP).fillna(999)
    df["modality_order"] = df["modality"].map(MODALITY_ORDER).fillna(99)

    df = (
        df.sort_values(
            ["disease_order", "modality_order", "organ_label", "disease_clock"]
        )
        .drop(columns=["disease_order", "modality_order"])
        .reset_index(drop=True)
    )

    return df


def save_summary_tables(df_final, output_dir_result):
    df_ok = df_final.loc[df_final["parse_status"] == "ok"].copy()

    if df_ok.empty:
        print("No OK LDSC h2 results to summarize.")
        return pd.DataFrame(), pd.DataFrame()

    summary_by_disease_modality = (
        df_ok.groupby(["disease_label", "modality"], as_index=False)
        .agg(
            n_clocks=("disease_clock", "count"),

            mean_h2=("h2_mean", "mean"),
            median_h2=("h2_mean", "median"),
            min_h2=("h2_mean", "min"),
            max_h2=("h2_mean", "max"),
            mean_h2_se=("h2_std", "mean"),

            mean_intercept=("intercept", "mean"),
            median_intercept=("intercept", "median"),
            min_intercept=("intercept", "min"),
            max_intercept=("intercept", "max"),
            mean_intercept_se=("intercept_std", "mean"),

            mean_ratio=("ratio", "mean"),
            median_ratio=("ratio", "median"),
            min_ratio=("ratio", "min"),
            max_ratio=("ratio", "max"),
            mean_ratio_se=("ratio_std", "mean"),

            mean_lambda_gc=("lambda_gc", "mean"),
            median_lambda_gc=("lambda_gc", "median"),
            min_lambda_gc=("lambda_gc", "min"),
            max_lambda_gc=("lambda_gc", "max"),

            mean_chi2=("mean_chi2", "mean"),
            median_chi2=("mean_chi2", "median"),
            min_chi2=("mean_chi2", "min"),
            max_chi2=("mean_chi2", "max"),

            n_h2_p05=("h2_significant_p05", "sum"),
            n_h2_bonferroni_47=("h2_significant_bonferroni_47", "sum"),
        )
    )

    summary_by_disease_modality["disease_order"] = (
        summary_by_disease_modality["disease_label"]
        .map({DISEASE_LABELS[x]: i + 1 for i, x in enumerate(DISEASE_ORDER)})
        .fillna(999)
    )

    summary_by_disease_modality["modality_order"] = (
        summary_by_disease_modality["modality"]
        .map(MODALITY_ORDER)
        .fillna(99)
    )

    summary_by_disease_modality = (
        summary_by_disease_modality
        .sort_values(["disease_order", "modality_order"])
        .drop(columns=["disease_order", "modality_order"])
        .reset_index(drop=True)
    )

    qc_flags = df_final[
        [
            "disease_label",
            "disease_key",
            "modality",
            "organ_label",
            "disease_clock",
            "h2_mean",
            "h2_std",
            "h2_z_calc",
            "h2_p_calc",
            "intercept",
            "intercept_std",
            "ratio",
            "ratio_std",
            "lambda_gc",
            "mean_chi2",
            "parse_status",
            "file_trait_key",
            "file_trait_compatible",
            "clock_file_trait_mismatch",
            "intercept_gt_1p05",
            "intercept_gt_1p10",
            "lambda_gc_gt_1p10",
            "ratio_gt_0p20",
            "negative_h2",
            "h2_significant_p05",
            "h2_significant_bonferroni_47",
            "h2_log_relative",
        ]
    ].copy()

    out_summary = os.path.join(
        output_dir_result,
        "LDSC_h2_intercept_47_disease_clocks_summary_by_disease_modality.tsv"
    )

    out_qc = os.path.join(
        output_dir_result,
        "LDSC_h2_intercept_47_disease_clocks_QC_flags.tsv"
    )

    summary_by_disease_modality.to_csv(out_summary, sep="\t", index=False)
    qc_flags.to_csv(out_qc, sep="\t", index=False)

    print(f"Saved: {out_summary}")
    print(f"Saved: {out_qc}")

    return summary_by_disease_modality, qc_flags


def save_excel_workbook(df_final, df_ok, df_failed, summary_tbl, qc_flags, output_dir_result):
    out_xlsx = os.path.join(
        output_dir_result,
        "LDSC_h2_intercept_47_disease_clocks.xlsx"
    )

    readme = pd.DataFrame(
        {
            "field": [
                "Purpose",
                "Input pattern",
                "Disease endpoints",
                "Detected h2 logs",
                "Expected h2 logs",
                "Main h2 field",
                "QC fields",
                "Generated on",
            ],
            "description": [
                "Collect LDSC SNP-based h2, intercept, ratio, lambda GC, and mean chi-square for disease EPOCH clocks.",
                "<base_dir>/*_clock/fastGWA/ldsc/*_h2.log restricted to asthma, copd, dementia, mi, and stroke clock folders.",
                ", ".join(DISEASE_ORDER),
                str(df_final.shape[0]),
                str(EXPECTED_N_LOGS),
                "Total Observed scale h2 from LDSC.",
                "Flags include intercept > 1.05, intercept > 1.10, lambda GC > 1.10, ratio > 0.20, negative h2, and h2 significance.",
                pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"),
            ],
        }
    )

    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
        readme.to_excel(writer, index=False, sheet_name="README")
        df_final.to_excel(writer, index=False, sheet_name="Disease_LDSC_h2")
        df_ok.to_excel(writer, index=False, sheet_name="OK_only")
        df_failed.to_excel(writer, index=False, sheet_name="Failed_or_nan")
        summary_tbl.to_excel(writer, index=False, sheet_name="Summary_disease_modality")
        qc_flags.to_excel(writer, index=False, sheet_name="QC_flags")

        workbook = writer.book

        for sheet_name in writer.sheets:
            ws = writer.sheets[sheet_name]
            ws.freeze_panes = "A2"

            # Auto-width with a cap.
            for column_cells in ws.columns:
                max_length = 0
                col_letter = column_cells[0].column_letter

                for cell in column_cells:
                    try:
                        value = str(cell.value) if cell.value is not None else ""
                    except Exception:
                        value = ""
                    max_length = max(max_length, len(value))

                ws.column_dimensions[col_letter].width = min(max(max_length + 2, 10), 45)

    print(f"Saved: {out_xlsx}")


def print_qc(df_final):
    print("\nQuick QC:")
    print(f"Total h2 logs collected: {df_final.shape[0]}")
    print(f"Unique disease clocks: {df_final['disease_clock'].nunique()}")

    print("\nParse status:")
    print(df_final["parse_status"].value_counts(dropna=False))

    print("\nBy disease, modality, and parse status:")
    print(
        df_final
        .groupby(["disease_label", "modality", "parse_status"], dropna=False)
        .size()
    )

    n_mismatch = int(df_final["clock_file_trait_mismatch"].fillna(False).sum())

    if n_mismatch > 0:
        print(f"\nWARNING: {n_mismatch} rows have file-trait mismatch.")
        print(
            df_final.loc[
                df_final["clock_file_trait_mismatch"].fillna(False),
                [
                    "disease_label",
                    "modality",
                    "organ_label",
                    "file_trait_key",
                    "h2_log_relative",
                ],
            ]
        )

    print("\nFinal table preview:")
    preview_cols = [
        "disease_label",
        "modality",
        "organ_label",
        "h2_mean",
        "h2_std",
        "h2_z_calc",
        "h2_p_calc_display",
        "intercept",
        "intercept_std",
        "ratio",
        "ratio_std",
        "lambda_gc",
        "mean_chi2",
        "parse_status",
    ]

    print(df_final[preview_cols])


# ============================================================
# 6. Main collection function
# ============================================================

def collect_ldsc_h2_intercept_47_disease_clocks(base_dir, output_dir_result):
    os.makedirs(output_dir_result, exist_ok=True)

    h2_logs = collect_disease_ldsc_h2_logs(base_dir)

    print(f"Base directory: {base_dir}")
    print(f"Output directory: {output_dir_result}")
    print(f"Found {len(h2_logs)} disease LDSC h2 log files.")

    if len(h2_logs) == 0:
        raise RuntimeError(f"No disease LDSC h2 logs found under: {base_dir}")

    if len(h2_logs) != EXPECTED_N_LOGS:
        print(
            f"WARNING: Expected {EXPECTED_N_LOGS} disease LDSC h2 logs, "
            f"but found {len(h2_logs)}. Continuing with detected logs."
        )

    rows = []

    for h2_log in h2_logs:
        clock_folder = os.path.basename(
            os.path.dirname(
                os.path.dirname(
                    os.path.dirname(h2_log)
                )
            )
        )

        meta = infer_metadata_from_folder(clock_folder)

        file_trait_key = parse_trait_from_h2_log(h2_log)
        compatible = is_file_trait_compatible(
            meta["organ_key"],
            file_trait_key
        )

        print("\n===================================================")
        print(f"Clock folder: {clock_folder}")
        print(f"Disease:      {meta['disease_label']}")
        print(f"Organ:        {meta['organ_label']}")
        print(f"Modality:     {meta['modality']}")
        print(f"Log file:     {os.path.basename(h2_log)}")

        try:
            parsed = parse_ldsc_h2_log(h2_log)
            error = ""
        except Exception as e:
            parsed = {
                "h2_mean": np.nan,
                "h2_std": np.nan,
                "h2_z_from_log": np.nan,
                "h2_z_calc": np.nan,
                "h2_p_calc": np.nan,
                "h2_p_calc_display": "",

                "intercept": np.nan,
                "intercept_std": np.nan,

                "ratio": np.nan,
                "ratio_std": np.nan,

                "lambda_gc": np.nan,
                "mean_chi2": np.nan,

                "h2_found": False,
                "intercept_found": False,
                "ratio_found": False,
                "lambda_gc_found": False,
                "mean_chi2_found": False,
                "h2_z_found": False,

                "parse_status": "error",
            }
            error = str(e)

        h2_log_relative = os.path.relpath(h2_log, base_dir)

        row = {
            "clock_folder": clock_folder,
            "disease_clock": meta["disease_clock"],
            "disease_key": meta["disease_key"],
            "disease_label": meta["disease_label"],
            "modality": meta["modality"],
            "modality_key": meta["modality_key"],
            "organ_key": meta["organ_key"],
            "organ_label": meta["organ_label"],
            "organ_raw": meta["organ_raw"],

            "file_trait_key": file_trait_key,
            "file_trait_compatible": compatible,
            "clock_file_trait_mismatch": (
                not compatible if isinstance(compatible, bool) else np.nan
            ),

            "h2_log": h2_log,
            "h2_log_relative": h2_log_relative,
            "h2_log_basename": os.path.basename(h2_log),
            "error": error,
        }

        row.update(parsed)

        rows.append(row)

        print(
            f"h2={row['h2_mean']}, "
            f"se={row['h2_std']}, "
            f"intercept={row['intercept']}, "
            f"ratio={row['ratio']}, "
            f"lambda_gc={row['lambda_gc']}, "
            f"status={row['parse_status']}"
        )

    df_final = pd.DataFrame(rows)

    # QC annotations.
    df_final["intercept_gt_1p05"] = df_final["intercept"] > 1.05
    df_final["intercept_gt_1p10"] = df_final["intercept"] > 1.10
    df_final["lambda_gc_gt_1p10"] = df_final["lambda_gc"] > 1.10
    df_final["ratio_gt_0p20"] = df_final["ratio"] > 0.20
    df_final["negative_h2"] = df_final["h2_mean"] < 0
    df_final["h2_significant_p05"] = df_final["h2_p_calc"] < 0.05
    df_final["h2_significant_bonferroni_47"] = (
        df_final["h2_p_calc"] < 0.05 / EXPECTED_N_LOGS
    )

    df_final = sort_disease_h2_table(df_final)

    df_ok = df_final.loc[df_final["parse_status"] == "ok"].copy()
    df_failed = df_final.loc[df_final["parse_status"] != "ok"].copy()

    # Main outputs.
    out_tsv = os.path.join(
        output_dir_result,
        "LDSC_h2_intercept_47_disease_clocks.tsv"
    )

    out_ok = os.path.join(
        output_dir_result,
        "LDSC_h2_intercept_47_disease_clocks_ok_only.tsv"
    )

    out_failed = os.path.join(
        output_dir_result,
        "LDSC_h2_intercept_47_disease_clocks_failed_or_nan.tsv"
    )

    df_final.to_csv(out_tsv, index=False, sep="\t", encoding="utf-8")
    df_ok.to_csv(out_ok, index=False, sep="\t", encoding="utf-8")
    df_failed.to_csv(out_failed, index=False, sep="\t", encoding="utf-8")

    print(f"\nSaved: {out_tsv}")
    print(f"Saved: {out_ok}")
    print(f"Saved: {out_failed}")

    summary_tbl, qc_flags = save_summary_tables(
        df_final,
        output_dir_result
    )

    save_excel_workbook(
        df_final=df_final,
        df_ok=df_ok,
        df_failed=df_failed,
        summary_tbl=summary_tbl,
        qc_flags=qc_flags,
        output_dir_result=output_dir_result,
    )

    print_qc(df_final)

    return df_final


# ============================================================
# 7. Run
# ============================================================

if __name__ == "__main__":

    df_ldsc_h2 = collect_ldsc_h2_intercept_47_disease_clocks(
        base_dir=BASE_DIR,
        output_dir_result=OUTPUT_DIR_RESULT,
    )

    print("\nDone.")