import os
import glob
import numpy as np
import pandas as pd


# ============================================================
# 1. Paths
# ============================================================

GCTB_ROOT_DIR = (
    "/Users/hao/cubic-home/Reproducibile_paper/"
    "WholeBodyClock/GCTB"
)

OUTPUT_DIR_RESULT = (
    "/Users/hao/cubic-home/Reproducibile_paper/"
    "WholeBodyClock/Result"
)

# If running directly on CUBIC, use:
# GCTB_ROOT_DIR = (
#     "/cbica/home/wenju/Reproducibile_paper/"
#     "WholeBodyClock/GCTB"
# )
# OUTPUT_DIR_RESULT = (
#     "/cbica/home/wenju/Reproducibile_paper/"
#     "WholeBodyClock/Result"
# )


# ============================================================
# 2. Helpers
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


def get_clock_name_from_folder(clock_folder):
    """
    Examples:
      adipose_mri_mortality_clock -> adipose_mri
      Brain_proteomics_mortality_clock -> Brain_proteomics
      Endocrine_metabolomics_mortality_clock -> Endocrine_metabolomics
    """
    return clock_folder.replace("_mortality_clock", "")


def safe_float(x, default=np.nan):
    try:
        return float(x)
    except Exception:
        return default


def get_parres_value(df, parameter, column):
    """
    Extract values from GCTB .parRes table.

    Expected format:
                    Mean            SD
          Pi        ...
          NnzSnp    ...
          SigmaSq   ...
          S         ...
          ResVar    ...
          GenVar    ...
          SigmaSqG  ...
          hsq       ...
    """
    if parameter not in df.index:
        return np.nan

    if column not in df.columns:
        return np.nan

    return safe_float(df.loc[parameter, column])


def read_sbayess_parres(parres_file):
    """
    Read one GCTB SBayesS .parRes file and extract:
      hsq     = SNP heritability
      S       = selection signature
      Pi      = polygenicity proportion
      NnzSnp  = estimated number of nonzero SNPs

    Example .parRes rows:
      Pi
      NnzSnp
      SigmaSq
      S
      ResVar
      GenVar
      SigmaSqG
      hsq
    """

    df = pd.read_csv(
        parres_file,
        sep=r"\s+",
        engine="python",
        index_col=0,
    )

    df.columns = [str(c).strip() for c in df.columns]
    df.index = df.index.astype(str).str.strip()

    if "Mean" not in df.columns or "SD" not in df.columns:
        raise ValueError(f"Cannot find Mean and SD columns in: {parres_file}")

    out = {
        "Pi": get_parres_value(df, "Pi", "Mean"),
        "Pi_se": get_parres_value(df, "Pi", "SD"),

        "NnzSnp": get_parres_value(df, "NnzSnp", "Mean"),
        "NnzSnp_se": get_parres_value(df, "NnzSnp", "SD"),

        "SigmaSq": get_parres_value(df, "SigmaSq", "Mean"),
        "SigmaSq_se": get_parres_value(df, "SigmaSq", "SD"),

        "S": get_parres_value(df, "S", "Mean"),
        "S_se": get_parres_value(df, "S", "SD"),

        "ResVar": get_parres_value(df, "ResVar", "Mean"),
        "ResVar_se": get_parres_value(df, "ResVar", "SD"),

        "GenVar": get_parres_value(df, "GenVar", "Mean"),
        "GenVar_se": get_parres_value(df, "GenVar", "SD"),

        "SigmaSqG": get_parres_value(df, "SigmaSqG", "Mean"),
        "SigmaSqG_se": get_parres_value(df, "SigmaSqG", "SD"),

        "h2_mean": get_parres_value(df, "hsq", "Mean"),
        "h2_se": get_parres_value(df, "hsq", "SD"),

        "n_parameter_rows": df.shape[0],
        "parres_parameters": ";".join(df.index.astype(str).tolist()),
    }

    return out


def find_parres_file(clock_dir):
    """
    Find .parRes file in one mortality-clock folder.
    """
    parres_files = sorted(glob.glob(os.path.join(clock_dir, "*.parRes")))

    if len(parres_files) == 0:
        return None, []

    return parres_files[0], parres_files


# ============================================================
# 3. Main collection function
# ============================================================

def collect_sbayess_parameters_mortality_clocks(gctb_root_dir, output_dir_result):
    """
    Collect SBayesS h2, S, Pi, and NnzSnp parameters for all 22 mortality clocks.
    """

    os.makedirs(output_dir_result, exist_ok=True)

    clock_dirs = sorted(glob.glob(os.path.join(gctb_root_dir, "*_mortality_clock")))

    if len(clock_dirs) == 0:
        raise RuntimeError(f"No mortality-clock folders found under: {gctb_root_dir}")

    rows = []

    for clock_dir in clock_dirs:
        clock_folder = os.path.basename(clock_dir)
        mortality_clock = get_clock_name_from_folder(clock_folder)
        modality = infer_modality_from_clock_name(mortality_clock)

        snpres_files = sorted(glob.glob(os.path.join(clock_dir, "*.snpRes")))
        parres_file, all_parres_files = find_parres_file(clock_dir)

        row = {
            "clock_folder": clock_folder,
            "mortality_clock": mortality_clock,
            "modality": modality,
            "clock_dir": clock_dir,
            "n_snpRes_files": len(snpres_files),
            "snpRes_file": snpres_files[0] if len(snpres_files) > 0 else "",
            "n_parRes_files": len(all_parres_files),
            "parRes_file": parres_file if parres_file is not None else "",
            "status": "",
            "error": "",
        }

        print("\n===================================================")
        print(f"Clock:    {mortality_clock}")
        print(f"Modality: {modality}")

        if parres_file is None:
            row.update({
                "h2_mean": np.nan,
                "h2_se": np.nan,
                "S": np.nan,
                "S_se": np.nan,
                "Pi": np.nan,
                "Pi_se": np.nan,
                "NnzSnp": np.nan,
                "NnzSnp_se": np.nan,
                "SigmaSq": np.nan,
                "SigmaSq_se": np.nan,
                "ResVar": np.nan,
                "ResVar_se": np.nan,
                "GenVar": np.nan,
                "GenVar_se": np.nan,
                "SigmaSqG": np.nan,
                "SigmaSqG_se": np.nan,
                "n_parameter_rows": np.nan,
                "parres_parameters": "",
            })
            row["status"] = "missing_parRes"
            row["error"] = "No .parRes file found."
            print("WARNING: missing .parRes file.")

        else:
            try:
                param_info = read_sbayess_parres(parres_file)
                row.update(param_info)
                row["status"] = "ok"

                print(
                    f"h2={row['h2_mean']:.4f}, "
                    f"S={row['S']:.4f}, "
                    f"Pi={row['Pi']:.4f}, "
                    f"NnzSnp={row['NnzSnp']:.1f}"
                )

            except Exception as e:
                row.update({
                    "h2_mean": np.nan,
                    "h2_se": np.nan,
                    "S": np.nan,
                    "S_se": np.nan,
                    "Pi": np.nan,
                    "Pi_se": np.nan,
                    "NnzSnp": np.nan,
                    "NnzSnp_se": np.nan,
                    "SigmaSq": np.nan,
                    "SigmaSq_se": np.nan,
                    "ResVar": np.nan,
                    "ResVar_se": np.nan,
                    "GenVar": np.nan,
                    "GenVar_se": np.nan,
                    "SigmaSqG": np.nan,
                    "SigmaSqG_se": np.nan,
                    "n_parameter_rows": np.nan,
                    "parres_parameters": "",
                })
                row["status"] = "parse_error"
                row["error"] = str(e)
                print(f"ERROR parsing {parres_file}: {e}")

        rows.append(row)

    df_final = pd.DataFrame(rows)

    # Sorting.
    modality_order = {
        "MRI": 1,
        "Proteomics": 2,
        "Metabolomics": 3,
        "Unknown": 99,
    }

    df_final["modality_order"] = df_final["modality"].map(modality_order).fillna(99)

    df_final = (
        df_final.sort_values(["modality_order", "mortality_clock"])
        .drop(columns=["modality_order"])
        .reset_index(drop=True)
    )

    # Main output.
    out_tsv = os.path.join(
        output_dir_result,
        "GCTB_SBayesS_parameters_mortality_clocks.tsv",
    )

    out_xlsx = os.path.join(
        output_dir_result,
        "GCTB_SBayesS_parameters_mortality_clocks.xlsx",
    )

    df_final.to_csv(out_tsv, index=False, sep="\t", encoding="utf-8")
    df_final.to_excel(out_xlsx, index=False)

    print(f"\nSaved: {out_tsv}")
    print(f"Saved: {out_xlsx}")

    # OK-only output.
    df_ok = df_final.loc[df_final["status"] == "ok"].copy()

    out_ok = os.path.join(
        output_dir_result,
        "GCTB_SBayesS_parameters_mortality_clocks_ok_only.tsv",
    )

    df_ok.to_csv(out_ok, index=False, sep="\t", encoding="utf-8")
    print(f"Saved: {out_ok}")

    # Summary by modality.
    if not df_ok.empty:
        summary_by_modality = (
            df_ok.groupby("modality", as_index=False)
            .agg(
                n_clocks=("mortality_clock", "count"),

                mean_h2=("h2_mean", "mean"),
                median_h2=("h2_mean", "median"),
                min_h2=("h2_mean", "min"),
                max_h2=("h2_mean", "max"),

                mean_S=("S", "mean"),
                median_S=("S", "median"),
                min_S=("S", "min"),
                max_S=("S", "max"),

                mean_Pi=("Pi", "mean"),
                median_Pi=("Pi", "median"),
                min_Pi=("Pi", "min"),
                max_Pi=("Pi", "max"),

                mean_NnzSnp=("NnzSnp", "mean"),
                median_NnzSnp=("NnzSnp", "median"),
                min_NnzSnp=("NnzSnp", "min"),
                max_NnzSnp=("NnzSnp", "max"),
            )
        )

        out_summary = os.path.join(
            output_dir_result,
            "GCTB_SBayesS_parameters_summary_by_modality.tsv",
        )

        summary_by_modality.to_csv(out_summary, index=False, sep="\t", encoding="utf-8")
        print(f"Saved: {out_summary}")

    print("\nQuick QC:")
    print(f"Total mortality-clock folders: {df_final.shape[0]}")
    print(df_final["status"].value_counts(dropna=False))
    print("\nBy modality:")
    print(df_final.groupby(["modality", "status"]).size())

    print("\nFinal table preview:")
    print(
        df_final[
            [
                "mortality_clock",
                "modality",
                "h2_mean",
                "h2_se",
                "S",
                "S_se",
                "Pi",
                "Pi_se",
                "NnzSnp",
                "NnzSnp_se",
                "status",
            ]
        ]
    )

    return df_final


# ============================================================
# 4. Run
# ============================================================

if __name__ == "__main__":

    df_sbayess = collect_sbayess_parameters_mortality_clocks(
        gctb_root_dir=GCTB_ROOT_DIR,
        output_dir_result=OUTPUT_DIR_RESULT,
    )

    print('Stop here...')