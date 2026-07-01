import os
import glob
import re
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt


# ============================================================
# 1. Define all input sources
# ============================================================

H2_SOURCES = [
    {
        "clock_class": "Mortality_L_EPOCH",
        "source_name": "mortality_clocks_all_modalities",
        "root_dir": "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/mortality_clock/GCTA_h2",
        "clock_name_mode": "mortality_folder",
    },
    {
        "clock_class": "Aging_BAG",
        "source_name": "MRI_BAG",
        "root_dir": "/Users/hao/cubic-home/Reproducibile_paper/AbdoImaging/h2_gcta",
        "clock_name_mode": "aging_mri_folder",
    },
    {
        "clock_class": "Aging_BAG",
        "source_name": "ProtBAG",
        "root_dir": "/Users/hao/cubic-home/Reproducibile_paper/UKBB_Proteomics/h2",
        "clock_name_mode": "aging_proteomics_folder",
    },
    {
        "clock_class": "Aging_BAG",
        "source_name": "MetBAG",
        "root_dir": "/Users/hao/cubic-home/Reproducibile_paper/UKBB_metabolomics/h2",
        "clock_name_mode": "aging_metabolomics_folder",
    },
]


# ============================================================
# 2. Helper functions
# ============================================================

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


def get_clock_name_from_folder(clock_folder, clock_name_mode):
    """
    Preserve full original naming as much as possible.

    Mortality examples:
      adipose_mri_mortality_clock -> adipose_mri
      Brain_proteomics_mortality_clock -> Brain_proteomics
      Endocrine_metabolomics_mortality_clock -> Endocrine_metabolomics

    Aging MRI examples:
      adipose -> adipose_mri
      brain -> brain_mri

    Aging proteomics examples:
      Brain -> Brain_proteomics
      Reproductive_male -> Reproductive_male_proteomics

    Aging metabolomics examples:
      Endocrine -> Endocrine_metabolomics
      Digestive -> Digestive_metabolomics
    """

    if clock_name_mode == "mortality_folder":
        return clock_folder.replace("_mortality_clock", "")

    if clock_name_mode == "aging_mri_folder":
        return f"{clock_folder}_mri"

    if clock_name_mode == "aging_proteomics_folder":
        return f"{clock_folder}_proteomics"

    if clock_name_mode == "aging_metabolomics_folder":
        return f"{clock_folder}_metabolomics"

    return clock_folder


def safe_float(x, default=np.nan):
    try:
        return float(x)
    except Exception:
        return default


def p_to_log10(pvalue):
    """
    Return log10(P) and -log10(P). If P is 0 due to numerical underflow,
    report -inf and inf.
    """
    p = safe_float(pvalue, default=np.nan)

    if not np.isfinite(p):
        return np.nan, np.nan

    if p == 0:
        return -np.inf, np.inf

    return np.log10(p), -np.log10(p)


def mixture_chisq_p_from_lrt(lrt):
    """
    GCTA REML variance-component LRT is commonly evaluated using a 50:50
    mixture of point mass at zero and chi-square_1.

    p = 0.5 * P(chi-square_1 >= LRT)
    """
    if not np.isfinite(lrt):
        return np.nan, np.nan, np.nan

    log_p = np.log(0.5) + stats.chi2.logsf(lrt, 1)
    log10_p = log_p / np.log(10)
    neg_log10_p = -log10_p

    if log10_p < -300:
        p = 0.0
    else:
        p = 10 ** log10_p

    return p, log10_p, neg_log10_p


# ============================================================
# 3. Parsers for GCTA outputs
# ============================================================

def read_gcta_hsq(hsq_file):
    """
    Read one GCTA REML .hsq file and extract h2-related statistics.

    Expected GCTA REML .hsq format:
      Source    Variance    SE
      V(G)      ...
      V(e)      ...
      Vp        ...
      V(G)/Vp   ...
      logL      ...
      logL0     ...
      LRT       ...
      df        ...
      Pval      ...
      n         ...
    """

    df = pd.read_csv(hsq_file, sep="\t")
    df.columns = [str(c).strip() for c in df.columns]

    if "Source" not in df.columns:
        df = pd.read_csv(hsq_file, sep=r"\s+", engine="python")
        df.columns = [str(c).strip() for c in df.columns]

    if "Source" not in df.columns or "Variance" not in df.columns:
        raise ValueError(f"Cannot parse .hsq file: {hsq_file}")

    df["Source"] = df["Source"].astype(str).str.strip()

    def get_value(source_name, col="Variance", default=np.nan):
        tmp = df.loc[df["Source"] == source_name]
        if tmp.empty or col not in tmp.columns:
            return default
        return safe_float(tmp.iloc[0][col], default=default)

    def get_raw(source_name, col="Variance", default="NA"):
        tmp = df.loc[df["Source"] == source_name]
        if tmp.empty or col not in tmp.columns:
            return default
        return str(tmp.iloc[0][col])

    vg = get_value("V(G)", "Variance")
    vg_se = get_value("V(G)", "SE")

    ve = get_value("V(e)", "Variance")
    ve_se = get_value("V(e)", "SE")

    vp = get_value("Vp", "Variance")
    vp_se = get_value("Vp", "SE")

    h2 = get_value("V(G)/Vp", "Variance")
    h2_se = get_value("V(G)/Vp", "SE")

    logL = get_value("logL", "Variance")
    logL0 = get_value("logL0", "Variance")
    lrt = get_value("LRT", "Variance")
    df_lrt = get_value("df", "Variance")
    n = get_value("n", "Variance")

    pval_raw = get_raw("Pval", "Variance")

    try:
        pvalue = float(pval_raw)
    except Exception:
        pvalue = np.nan

    pvalue_from_lrt = np.nan
    pvalue_log10 = np.nan
    neg_log10_pvalue = np.nan

    if np.isfinite(lrt):
        pvalue_from_lrt, pvalue_log10, neg_log10_pvalue = mixture_chisq_p_from_lrt(lrt)

    if (pvalue == 0 or pval_raw == "0.0000e+00") and np.isfinite(pvalue_from_lrt):
        pvalue = pvalue_from_lrt

    if not np.isfinite(pvalue_log10) and np.isfinite(pvalue):
        pvalue_log10, neg_log10_pvalue = p_to_log10(pvalue)

    return {
        "result_type": "GCTA_REML_hsq",
        "h2_estimator": "REML_VG_over_Vp",
        "V_G": vg,
        "V_G_SE": vg_se,
        "V_e": ve,
        "V_e_SE": ve_se,
        "Vp": vp,
        "Vp_SE": vp_se,
        "h2": h2,
        "h2_SE": h2_se,
        "logL": logL,
        "logL0": logL0,
        "LRT": lrt,
        "df": df_lrt,
        "Pvalue_raw": pval_raw,
        "Pvalue": pvalue,
        "Pvalue_from_LRT_mixture": pvalue_from_lrt,
        "Pvalue_log10": pvalue_log10,
        "neg_log10_Pvalue": neg_log10_pvalue,
        "n": n,

        # HEreg-specific columns kept for a consistent output table.
        "HE_CP_h2": np.nan,
        "HE_CP_SE_OLS": np.nan,
        "HE_CP_SE_Jackknife": np.nan,
        "HE_CP_P_OLS": np.nan,
        "HE_CP_P_Jackknife": np.nan,
        "HE_SD_h2": np.nan,
        "HE_SD_SE_OLS": np.nan,
        "HE_SD_SE_Jackknife": np.nan,
        "HE_SD_P_OLS": np.nan,
        "HE_SD_P_Jackknife": np.nan,
    }


def read_gcta_hereg(hereg_file):
    """
    Read one GCTA-HEreg .HEreg file.

    Expected structure:
      HE-CP
      Coefficient     Estimate        SE_OLS          SE_Jackknife    P_OLS           P_Jackknife
      Intercept       ...
      V(G)/Vp         ...

      HE-SD
      Coefficient     Estimate        SE_OLS          SE_Jackknife    P_OLS           P_Jackknife
      Intercept       ...
      V(G)/Vp         ...

    Primary estimate reported in the unified columns:
      h2      = HE-CP V(G)/Vp Estimate
      h2_SE   = HE-CP V(G)/Vp SE_Jackknife
      Pvalue  = HE-CP V(G)/Vp P_Jackknife
    """

    sections = {}
    current_section = None
    header = None

    with open(hereg_file, "r") as f:
        for raw_line in f:
            line = raw_line.strip()

            if line == "":
                continue

            if line in ["HE-CP", "HE-SD"]:
                current_section = line
                sections[current_section] = {}
                header = None
                continue

            if current_section is not None and line.startswith("Coefficient"):
                header = re.split(r"\s+", line)
                continue

            if current_section is not None and header is not None:
                parts = re.split(r"\s+", line)

                if len(parts) < len(header):
                    continue

                coefficient = parts[0]
                values = dict(zip(header[1:], parts[1:]))
                sections[current_section][coefficient] = values

    def get(section, coefficient, col, default=np.nan):
        try:
            return safe_float(sections[section][coefficient][col], default=default)
        except Exception:
            return default

    def get_raw(section, coefficient, col, default="NA"):
        try:
            return str(sections[section][coefficient][col])
        except Exception:
            return default

    he_cp_h2 = get("HE-CP", "V(G)/Vp", "Estimate")
    he_cp_se_ols = get("HE-CP", "V(G)/Vp", "SE_OLS")
    he_cp_se_jackknife = get("HE-CP", "V(G)/Vp", "SE_Jackknife")
    he_cp_p_ols = get("HE-CP", "V(G)/Vp", "P_OLS")
    he_cp_p_jackknife = get("HE-CP", "V(G)/Vp", "P_Jackknife")
    he_cp_p_jackknife_raw = get_raw("HE-CP", "V(G)/Vp", "P_Jackknife")

    he_sd_h2 = get("HE-SD", "V(G)/Vp", "Estimate")
    he_sd_se_ols = get("HE-SD", "V(G)/Vp", "SE_OLS")
    he_sd_se_jackknife = get("HE-SD", "V(G)/Vp", "SE_Jackknife")
    he_sd_p_ols = get("HE-SD", "V(G)/Vp", "P_OLS")
    he_sd_p_jackknife = get("HE-SD", "V(G)/Vp", "P_Jackknife")

    he_cp_intercept = get("HE-CP", "Intercept", "Estimate")
    he_cp_intercept_se_jackknife = get("HE-CP", "Intercept", "SE_Jackknife")
    he_cp_intercept_p_jackknife = get("HE-CP", "Intercept", "P_Jackknife")

    he_sd_intercept = get("HE-SD", "Intercept", "Estimate")
    he_sd_intercept_se_jackknife = get("HE-SD", "Intercept", "SE_Jackknife")
    he_sd_intercept_p_jackknife = get("HE-SD", "Intercept", "P_Jackknife")

    pvalue_log10, neg_log10_pvalue = p_to_log10(he_cp_p_jackknife)

    return {
        "result_type": "GCTA_HEreg",
        "h2_estimator": "HE_CP_VG_over_Vp_SE_Jackknife",
        "V_G": np.nan,
        "V_G_SE": np.nan,
        "V_e": np.nan,
        "V_e_SE": np.nan,
        "Vp": np.nan,
        "Vp_SE": np.nan,

        # Unified primary columns.
        "h2": he_cp_h2,
        "h2_SE": he_cp_se_jackknife,
        "logL": np.nan,
        "logL0": np.nan,
        "LRT": np.nan,
        "df": np.nan,
        "Pvalue_raw": he_cp_p_jackknife_raw,
        "Pvalue": he_cp_p_jackknife,
        "Pvalue_from_LRT_mixture": np.nan,
        "Pvalue_log10": pvalue_log10,
        "neg_log10_Pvalue": neg_log10_pvalue,
        "n": np.nan,

        # HE-CP detailed columns.
        "HE_CP_h2": he_cp_h2,
        "HE_CP_SE_OLS": he_cp_se_ols,
        "HE_CP_SE_Jackknife": he_cp_se_jackknife,
        "HE_CP_P_OLS": he_cp_p_ols,
        "HE_CP_P_Jackknife": he_cp_p_jackknife,
        "HE_CP_Intercept": he_cp_intercept,
        "HE_CP_Intercept_SE_Jackknife": he_cp_intercept_se_jackknife,
        "HE_CP_Intercept_P_Jackknife": he_cp_intercept_p_jackknife,

        # HE-SD detailed columns.
        "HE_SD_h2": he_sd_h2,
        "HE_SD_SE_OLS": he_sd_se_ols,
        "HE_SD_SE_Jackknife": he_sd_se_jackknife,
        "HE_SD_P_OLS": he_sd_p_ols,
        "HE_SD_P_Jackknife": he_sd_p_jackknife,
        "HE_SD_Intercept": he_sd_intercept,
        "HE_SD_Intercept_SE_Jackknife": he_sd_intercept_se_jackknife,
        "HE_SD_Intercept_P_Jackknife": he_sd_intercept_p_jackknife,
    }


def read_gcta_result_file(result_file):
    """
    Dispatch parser based on file extension.
    """
    if result_file.endswith(".hsq"):
        return read_gcta_hsq(result_file)

    if result_file.endswith(".HEreg"):
        return read_gcta_hereg(result_file)

    raise ValueError(f"Unsupported GCTA result file type: {result_file}")


def collect_gcta_result_files_for_source(source):
    """
    Collect both .hsq and .HEreg files from one source directory.

    This supports:
      - MRI/proteomics mortality clocks: *.hsq
      - metabolomics mortality clocks: *.HEreg
      - aging clocks: *.hsq or *.HEreg if generated later
    """

    root_dir = source["root_dir"]

    if not os.path.exists(root_dir):
        print(f"WARNING: root directory does not exist: {root_dir}")
        return []

    result_files = []

    result_files.extend(sorted(glob.glob(os.path.join(root_dir, "*", "*.hsq"))))
    result_files.extend(sorted(glob.glob(os.path.join(root_dir, "*", "*.HEreg"))))

    result_files = sorted(result_files)

    return result_files


# ============================================================
# 4. Main collection function
# ============================================================

def collect_gcta_h2_results_all_sources(output_dir_result):
    """
    Collect h2 estimates from:
      1. mortality L'EPOCH clocks
      2. MRI BAGs
      3. ProtBAGs
      4. MetBAGs

    This version supports both:
      - GCTA REML .hsq
      - GCTA-HEreg .HEreg
    """

    os.makedirs(output_dir_result, exist_ok=True)

    rows = []

    for source in H2_SOURCES:
        clock_class = source["clock_class"]
        source_name = source["source_name"]
        root_dir = source["root_dir"]
        clock_name_mode = source["clock_name_mode"]

        print("\n=================================================")
        print(f"Collecting source: {source_name}")
        print(f"Clock class: {clock_class}")
        print(f"Root: {root_dir}")

        result_files = collect_gcta_result_files_for_source(source)

        print(f"Found {len(result_files)} GCTA result files: .hsq or .HEreg.")

        for result_file in result_files:
            clock_folder = os.path.basename(os.path.dirname(result_file))
            result_basename = os.path.basename(result_file)

            clock_name = get_clock_name_from_folder(
                clock_folder=clock_folder,
                clock_name_mode=clock_name_mode,
            )

            modality = infer_modality_from_clock_name(clock_name)

            try:
                h2_info = read_gcta_result_file(result_file)
                status = "ok"
                error = ""
            except Exception as e:
                h2_info = {}
                status = "error"
                error = str(e)

            row = {
                "clock_class": clock_class,
                "source_name": source_name,
                "clock_folder": clock_folder,
                "clock_name": clock_name,
                "modality": modality,
                "result_file": result_file,
                "result_basename": result_basename,

                # Backward-compatible names, in case your downstream code expects them.
                "hsq_file": result_file,
                "hsq_basename": result_basename,

                "status": status,
                "error": error,
            }

            row.update(h2_info)
            rows.append(row)

    if len(rows) == 0:
        raise RuntimeError("No .hsq or .HEreg files were found across all configured sources.")

    df_final = pd.DataFrame(rows)

    # Multiple-testing correction across all successfully parsed h2 estimates.
    n_tests_all = df_final.loc[df_final["status"] == "ok"].shape[0]
    bonf_all = 0.05 / n_tests_all if n_tests_all > 0 else np.nan

    df_final["n_tests_all_available"] = n_tests_all
    df_final["bonferroni_p_all"] = bonf_all

    df_final["Pvalue_nominal_significant"] = df_final["Pvalue"] < 0.05
    df_final["Pvalue_bonferroni_significant_all"] = df_final["Pvalue"] < bonf_all

    # Correction within each clock class: mortality vs aging.
    df_final["n_tests_within_clock_class"] = np.nan
    df_final["bonferroni_p_within_clock_class"] = np.nan
    df_final["Pvalue_bonferroni_significant_within_clock_class"] = False

    for clock_class, sub_idx in df_final.groupby("clock_class").groups.items():
        idx = list(sub_idx)
        n_tests = df_final.loc[idx].query("status == 'ok'").shape[0]
        bonf = 0.05 / n_tests if n_tests > 0 else np.nan

        df_final.loc[idx, "n_tests_within_clock_class"] = n_tests
        df_final.loc[idx, "bonferroni_p_within_clock_class"] = bonf
        df_final.loc[idx, "Pvalue_bonferroni_significant_within_clock_class"] = (
            df_final.loc[idx, "Pvalue"] < bonf
        )

    # Sorting.
    clock_class_order = {
        "Mortality_L_EPOCH": 1,
        "Aging_BAG": 2,
    }

    source_order = {
        "mortality_clocks_all_modalities": 1,
        "MRI_BAG": 2,
        "ProtBAG": 3,
        "MetBAG": 4,
    }

    modality_order = {
        "MRI": 1,
        "Proteomics": 2,
        "Metabolomics": 3,
        "Unknown": 99,
    }

    df_final["clock_class_order"] = df_final["clock_class"].map(clock_class_order).fillna(99)
    df_final["source_order"] = df_final["source_name"].map(source_order).fillna(99)
    df_final["modality_order"] = df_final["modality"].map(modality_order).fillna(99)

    df_final = (
        df_final.sort_values(
            ["clock_class_order", "source_order", "modality_order", "clock_name", "result_type"]
        )
        .drop(columns=["clock_class_order", "source_order", "modality_order"])
        .reset_index(drop=True)
    )

    # Save full outputs.
    out_tsv = os.path.join(
        output_dir_result,
        "GCTA_h2_results_mortality_and_aging_clocks.tsv",
    )

    out_xlsx = os.path.join(
        output_dir_result,
        "GCTA_h2_results_mortality_and_aging_clocks.xlsx",
    )

    df_final.to_csv(out_tsv, index=False, sep="\t", encoding="utf-8")
    df_final.to_excel(out_xlsx, index=False)

    print(f"\nSaved: {out_tsv}")
    print(f"Saved: {out_xlsx}")

    # Save mortality and aging separately.
    mortality_out = os.path.join(
        output_dir_result,
        "GCTA_h2_results_mortality_clocks_only.tsv",
    )

    aging_out = os.path.join(
        output_dir_result,
        "GCTA_h2_results_aging_clocks_only.tsv",
    )

    df_final.loc[df_final["clock_class"] == "Mortality_L_EPOCH"].to_csv(
        mortality_out,
        index=False,
        sep="\t",
        encoding="utf-8",
    )

    df_final.loc[df_final["clock_class"] == "Aging_BAG"].to_csv(
        aging_out,
        index=False,
        sep="\t",
        encoding="utf-8",
    )

    print(f"Saved: {mortality_out}")
    print(f"Saved: {aging_out}")

    # Save summary tables.
    save_summary_tables(df_final, output_dir_result)

    # Plot.
    plot_h2_barplot_by_clock_class(df_final, output_dir_result)

    print("\nQuick check:")
    print(
        df_final[
            [
                "clock_class",
                "source_name",
                "clock_name",
                "modality",
                "result_type",
                "h2_estimator",
                "h2",
                "h2_SE",
                "Pvalue",
                "neg_log10_Pvalue",
                "n",
                "status",
            ]
        ]
    )

    return df_final


# ============================================================
# 5. Summary tables
# ============================================================

def save_summary_tables(df_final, output_dir_result):
    df_ok = df_final.loc[df_final["status"] == "ok"].copy()

    summary_by_class_source = (
        df_ok.groupby(["clock_class", "source_name"], as_index=False)
        .agg(
            n_clocks=("clock_name", "count"),
            mean_h2=("h2", "mean"),
            median_h2=("h2", "median"),
            min_h2=("h2", "min"),
            max_h2=("h2", "max"),
            mean_h2_SE=("h2_SE", "mean"),
            mean_n=("n", "mean"),
            min_p=("Pvalue", "min"),
            max_neg_log10_p=("neg_log10_Pvalue", "max"),
            n_nominal_sig=("Pvalue_nominal_significant", "sum"),
            n_bonf_sig_all=("Pvalue_bonferroni_significant_all", "sum"),
            n_bonf_sig_within_class=(
                "Pvalue_bonferroni_significant_within_clock_class",
                "sum",
            ),
        )
        .sort_values(["clock_class", "source_name"])
    )

    out1 = os.path.join(
        output_dir_result,
        "GCTA_h2_summary_by_clock_class_and_source.tsv",
    )
    summary_by_class_source.to_csv(out1, index=False, sep="\t", encoding="utf-8")
    print(f"Saved: {out1}")

    summary_by_class_modality = (
        df_ok.groupby(["clock_class", "modality"], as_index=False)
        .agg(
            n_clocks=("clock_name", "count"),
            mean_h2=("h2", "mean"),
            median_h2=("h2", "median"),
            min_h2=("h2", "min"),
            max_h2=("h2", "max"),
            mean_h2_SE=("h2_SE", "mean"),
            mean_n=("n", "mean"),
            min_p=("Pvalue", "min"),
            max_neg_log10_p=("neg_log10_Pvalue", "max"),
            n_nominal_sig=("Pvalue_nominal_significant", "sum"),
            n_bonf_sig_all=("Pvalue_bonferroni_significant_all", "sum"),
            n_bonf_sig_within_class=(
                "Pvalue_bonferroni_significant_within_clock_class",
                "sum",
            ),
        )
        .sort_values(["clock_class", "modality"])
    )

    out2 = os.path.join(
        output_dir_result,
        "GCTA_h2_summary_by_clock_class_and_modality.tsv",
    )
    summary_by_class_modality.to_csv(out2, index=False, sep="\t", encoding="utf-8")
    print(f"Saved: {out2}")

    # Wide table for direct mortality-vs-aging comparison by clock_name.
    # This works when names match, e.g. brain_mri or Endocrine_metabolomics.
    df_wide = df_ok[
        [
            "clock_class",
            "clock_name",
            "modality",
            "result_type",
            "h2_estimator",
            "h2",
            "h2_SE",
            "Pvalue",
            "neg_log10_Pvalue",
            "n",
            "source_name",
        ]
    ].copy()

    wide = df_wide.pivot_table(
        index=["clock_name", "modality"],
        columns="clock_class",
        values=["h2", "h2_SE", "Pvalue", "neg_log10_Pvalue", "n"],
        aggfunc="first",
    )

    wide.columns = ["_".join([str(x) for x in col]).strip("_") for col in wide.columns]
    wide = wide.reset_index()

    if "h2_Mortality_L_EPOCH" in wide.columns and "h2_Aging_BAG" in wide.columns:
        wide["h2_diff_mortality_minus_aging"] = (
            wide["h2_Mortality_L_EPOCH"] - wide["h2_Aging_BAG"]
        )

    out3 = os.path.join(
        output_dir_result,
        "GCTA_h2_mortality_vs_aging_wide_by_clock_name.tsv",
    )
    wide.to_csv(out3, index=False, sep="\t", encoding="utf-8")
    print(f"Saved: {out3}")


# ============================================================
# 6. Plotting
# ============================================================

def plot_h2_barplot_by_clock_class(df_final, output_dir_result):
    """
    Simple grouped barplot of h2 estimates across mortality and aging clocks.
    """

    df_ok = df_final.loc[df_final["status"] == "ok"].copy()

    if df_ok.empty:
        print("No valid h2 estimates to plot.")
        return

    df_ok["plot_label"] = df_ok["clock_name"] + "\n" + df_ok["clock_class"]

    x = np.arange(df_ok.shape[0])
    y = df_ok["h2"].astype(float).values
    yerr = df_ok["h2_SE"].astype(float).values

    fig_width = max(10, 0.45 * df_ok.shape[0])
    plt.figure(figsize=(fig_width, 5.2))

    plt.bar(x, y)
    plt.errorbar(x, y, yerr=yerr, fmt="none", capsize=3, linewidth=1)

    plt.xticks(x, df_ok["plot_label"], rotation=70, ha="right", fontsize=7)
    plt.ylabel("SNP heritability, h²")
    plt.xlabel("")
    plt.title("GCTA SNP heritability of mortality L'EPOCH and aging clocks")

    ymax = np.nanmax(y + yerr)
    plt.ylim(0, min(1.05, max(0.1, ymax * 1.15)))

    plt.tight_layout()

    out_pdf = os.path.join(
        output_dir_result,
        "GCTA_h2_mortality_and_aging_clocks_barplot.pdf",
    )

    out_png = os.path.join(
        output_dir_result,
        "GCTA_h2_mortality_and_aging_clocks_barplot.png",
    )

    plt.savefig(out_pdf)
    plt.savefig(out_png, dpi=300)
    plt.close()

    print(f"Saved: {out_pdf}")
    print(f"Saved: {out_png}")


# ============================================================
# 7. Run
# ============================================================

if __name__ == "__main__":

    output_dir_result = (
        "/Users/hao/cubic-home/Reproducibile_paper/"
        "WholeBodyClock/Result"
    )

    collect_gcta_h2_results_all_sources(
        output_dir_result=output_dir_result,
    )