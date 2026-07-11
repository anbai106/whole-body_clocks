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

# If running directly on CUBIC, use:
# LDSC_ROOT_DIR = (
#     "/cbica/home/wenju/Reproducibile_paper/"
#     "WholeBodyClock/mortality_clock/fastGWA/output"
# )
# OUTPUT_DIR_RESULT = (
#     "/cbica/home/wenju/Reproducibile_paper/"
#     "WholeBodyClock/Result"
# )


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
    elif x.endswith("_proteomics"):
        return "Proteomics"
    elif x.endswith("_metabolomics"):
        return "Metabolomics"
    else:
        return "Unknown"


def get_clock_name_from_folder(clock_folder):
    """
    Examples:
      adipose_mri_mortality_clock -> adipose_mri
      Brain_proteomics_mortality_clock -> Brain_proteomics
      Endocrine_metabolomics_mortality_clock -> Endocrine_metabolomics
    """
    return clock_folder.replace("_mortality_clock", "")


def parse_clock_name_from_h2_log(h2_log):
    """
    Example:
      adipose_mri_h2.log -> adipose_mri
      Brain_proteomics_h2.log -> Brain_proteomics
      Endocrine_metabolomics_h2.log -> Endocrine_metabolomics
    """
    basename = os.path.basename(h2_log)
    return re.sub(r"_h2\.log$", "", basename)


def extract_value_and_se(text, label):
    """
    Extract value and optional SE from LDSC lines such as:
      Intercept: 1.0123 (0.0078)
      Total Observed scale h2: 0.1234 (0.0456)
      Ratio: 0.0345 (0.0123)

    Returns:
      value, se, found
    """

    pattern = (
        rf"{re.escape(label)}\s*:\s*"
        r"([+-]?(?:nan|inf|-inf|[0-9.eE+-]+))"
        r"(?:\s*\(([+-]?(?:nan|inf|-inf|[0-9.eE+-]+))\))?"
    )

    m = re.search(pattern, text, flags=re.IGNORECASE)

    if m is None:
        return np.nan, np.nan, False

    value = safe_float(m.group(1))
    se = safe_float(m.group(2))

    return value, se, True

def extract_single_value(text, label):
    """
    Extract Mean Chi^2: 1.021
    """

    pattern = (
        rf"{re.escape(label)}\s*:\s*"
        r"([+-]?(?:nan|inf|-inf|[0-9.eE+-]+))"
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

    Some LDSC logs may show:
      Ratio < 0 ...
    or omit ratio. In that case ratio is set to NA.
    """

    with open(h2_log, "r", errors="ignore") as f:
        txt = f.read()

    # Main LDSC parameters.
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

    # Some LDSC outputs include a non-estimable ratio line.
    if not ratio_found:
        if re.search(r"Ratio\s+<|Ratio\s+not|Ratio\s+undefined|Ratio\s+NA", txt, flags=re.IGNORECASE):
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

    # Optional LDSC fields.
    observed_scale_h2_p, h2_p_found = extract_single_value(
        txt,
        "h2 Z"
    )

    # Parse status.
    required_found = h2_found and intercept_found

    if required_found:
        parse_status = "ok"
    elif "nan" in txt.lower():
        parse_status = "nan_result"
    elif "error" in txt.lower() or "traceback" in txt.lower():
        parse_status = "ldsc_error"
    else:
        parse_status = "parse_failed"

    return {
        "h2_mean": h2_mean,
        "h2_std": h2_std,
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
        "parse_status": parse_status,
    }


# ============================================================
# 4. Collect all LDSC h2 logs
# ============================================================

def collect_ldsc_h2_logs(ldsc_root_dir):
    """
    Collect h2 logs only.

    Expected:
      fastGWA/output/*_mortality_clock/ldsc/*_h2.log

    This excludes genetic-correlation logs such as *_vs_*.log.
    """

    pattern = os.path.join(
        ldsc_root_dir,
        "*_mortality_clock",
        "ldsc",
        "*_h2.log"
    )

    h2_logs = sorted(glob.glob(pattern))

    return h2_logs


def collect_ldsc_h2_intercept_mortality_clocks(ldsc_root_dir, output_dir_result):
    """
    Collect LDSC SNP-based h2, intercept, ratio, lambda GC, and mean chi-square
    for all 22 mortality clocks.
    """

    os.makedirs(output_dir_result, exist_ok=True)

    h2_logs = collect_ldsc_h2_logs(ldsc_root_dir)

    print(f"Found {len(h2_logs)} LDSC h2 log files.")

    if len(h2_logs) == 0:
        raise RuntimeError(f"No LDSC h2 logs found under: {ldsc_root_dir}")

    rows = []

    for h2_log in h2_logs:
        clock_folder = os.path.basename(
            os.path.dirname(os.path.dirname(h2_log))
        )

        mortality_clock = get_clock_name_from_folder(clock_folder)
        mortality_clock_from_file = parse_clock_name_from_h2_log(h2_log)

        modality = infer_modality_from_clock_name(mortality_clock)

        print("\n===================================================")
        print(f"Clock folder: {clock_folder}")
        print(f"Mortality clock: {mortality_clock}")
        print(f"Log clock name:  {mortality_clock_from_file}")
        print(f"Modality:        {modality}")

        try:
            parsed = parse_ldsc_h2_log(h2_log)
            error = ""
        except Exception as e:
            parsed = {
                "h2_mean": np.nan,
                "h2_std": np.nan,
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
                "parse_status": "error",
            }
            error = str(e)

        row = {
            "clock_folder": clock_folder,
            "mortality_clock": mortality_clock,
            "mortality_clock_from_file": mortality_clock_from_file,
            "modality": modality,
            "h2_log": h2_log,
            "h2_log_basename": os.path.basename(h2_log),
            "error": error,
        }

        row.update(parsed)

        row["clock_name_mismatch"] = (
            row["mortality_clock"] != row["mortality_clock_from_file"]
        )

        rows.append(row)

        print(
            f"h2={row['h2_mean']}, "
            f"intercept={row['intercept']}, "
            f"ratio={row['ratio']}, "
            f"lambda_gc={row['lambda_gc']}, "
            f"status={row['parse_status']}"
        )

    df_final = pd.DataFrame(rows)

    df_final = sort_ldsc_h2_table(df_final)

    # Save main outputs.
    out_tsv = os.path.join(
        output_dir_result,
        "LDSC_h2_intercept_mortality_clocks.tsv"
    )

    out_xlsx = os.path.join(
        output_dir_result,
        "LDSC_h2_intercept_mortality_clocks.xlsx"
    )

    df_final.to_csv(out_tsv, index=False, sep="\t", encoding="utf-8")
    df_final.to_excel(out_xlsx, index=False)

    print(f"\nSaved: {out_tsv}")
    print(f"Saved: {out_xlsx}")

    # Save OK-only.
    df_ok = df_final.loc[df_final["parse_status"] == "ok"].copy()

    out_ok = os.path.join(
        output_dir_result,
        "LDSC_h2_intercept_mortality_clocks_ok_only.tsv"
    )

    df_ok.to_csv(out_ok, index=False, sep="\t", encoding="utf-8")
    print(f"Saved: {out_ok}")

    # Save failed/nan logs.
    df_failed = df_final.loc[df_final["parse_status"] != "ok"].copy()

    out_failed = os.path.join(
        output_dir_result,
        "LDSC_h2_intercept_mortality_clocks_failed_or_nan.tsv"
    )

    df_failed.to_csv(out_failed, index=False, sep="\t", encoding="utf-8")
    print(f"Saved: {out_failed}")

    # Save summaries.
    save_ldsc_h2_summary_tables(df_final, output_dir_result)

    # QC print.
    print_qc(df_final)

    return df_final


# ============================================================
# 5. Sorting and summaries
# ============================================================

def sort_ldsc_h2_table(df_final):
    df = df_final.copy()

    clock_order = [
        # 7 MRI mortality clocks
        "adipose_mri",
        "brain_mri",
        "heart_mri",
        "kidney_mri",
        "liver_mri",
        "pancreas_mri",
        "spleen_mri",

        # 11 proteomics mortality clocks
        "Brain_proteomics",
        "Endocrine_proteomics",
        "Eye_proteomics",
        "Heart_proteomics",
        "Hepatic_proteomics",
        "Immune_proteomics",
        "Pulmonary_proteomics",
        "Renal_proteomics",
        "Reproductive_female_proteomics",
        "Reproductive_male_proteomics",
        "Skin_proteomics",

        # 4 metabolomics mortality clocks
        "Digestive_metabolomics",
        "Endocrine_metabolomics",
        "Hepatic_metabolomics",
        "Immune_metabolomics",
    ]

    modality_order = {
        "MRI": 1,
        "Proteomics": 2,
        "Metabolomics": 3,
        "Unknown": 99,
    }

    clock_order_map = {x: i + 1 for i, x in enumerate(clock_order)}

    df["clock_order"] = df["mortality_clock"].map(clock_order_map).fillna(999)
    df["modality_order"] = df["modality"].map(modality_order).fillna(99)

    df = (
        df.sort_values(["clock_order", "modality_order", "mortality_clock"])
        .drop(columns=["clock_order", "modality_order"])
        .reset_index(drop=True)
    )

    return df


def save_ldsc_h2_summary_tables(df_final, output_dir_result):
    df_ok = df_final.loc[df_final["parse_status"] == "ok"].copy()

    if df_ok.empty:
        print("No OK LDSC h2 results to summarize.")
        return

    summary_by_modality = (
        df_ok.groupby("modality", as_index=False)
        .agg(
            n_clocks=("mortality_clock", "count"),

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
        )
    )

    out_summary_modality = os.path.join(
        output_dir_result,
        "LDSC_h2_intercept_summary_by_modality.tsv"
    )

    summary_by_modality.to_csv(
        out_summary_modality,
        index=False,
        sep="\t",
        encoding="utf-8"
    )

    print(f"Saved: {out_summary_modality}")

    # Flag potentially problematic LDSC outputs.
    df_qc = df_ok.copy()

    df_qc["intercept_gt_1p05"] = df_qc["intercept"] > 1.05
    df_qc["intercept_gt_1p10"] = df_qc["intercept"] > 1.10
    df_qc["lambda_gc_gt_1p10"] = df_qc["lambda_gc"] > 1.10
    df_qc["ratio_gt_0p20"] = df_qc["ratio"] > 0.20
    df_qc["negative_h2"] = df_qc["h2_mean"] < 0

    out_qc = os.path.join(
        output_dir_result,
        "LDSC_h2_intercept_QC_flags.tsv"
    )

    df_qc.to_csv(
        out_qc,
        index=False,
        sep="\t",
        encoding="utf-8"
    )

    print(f"Saved: {out_qc}")


def print_qc(df_final):
    print("\nQuick QC:")
    print(f"Total h2 logs collected: {df_final.shape[0]}")
    print(f"Unique mortality clocks: {df_final['mortality_clock'].nunique()}")

    print("\nParse status:")
    print(df_final["parse_status"].value_counts(dropna=False))

    print("\nBy modality and parse status:")
    print(df_final.groupby(["modality", "parse_status"]).size())

    n_mismatch = int(df_final["clock_name_mismatch"].sum())
    if n_mismatch > 0:
        print(f"\nWARNING: {n_mismatch} rows have clock-name mismatch.")
        print(
            df_final.loc[
                df_final["clock_name_mismatch"],
                ["mortality_clock", "mortality_clock_from_file", "h2_log_basename"]
            ]
        )

    print("\nFinal table preview:")
    print(
        df_final[
            [
                "mortality_clock",
                "modality",
                "h2_mean",
                "h2_std",
                "intercept",
                "intercept_std",
                "ratio",
                "ratio_std",
                "lambda_gc",
                "mean_chi2",
                "parse_status",
            ]
        ]
    )


# ============================================================
# 6. Run
# ============================================================

if __name__ == "__main__":

    df_ldsc_h2 = collect_ldsc_h2_intercept_mortality_clocks(
        ldsc_root_dir=LDSC_ROOT_DIR,
        output_dir_result=OUTPUT_DIR_RESULT,
    )

    print('Stop...')