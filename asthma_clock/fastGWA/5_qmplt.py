#!/usr/bin/env python

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from qmplot import manhattanplot, qqplot


P_MIN = 1e-300


def read_fastgwa(path):
    """Read a GCTA fastGWA output file and standardize columns for qmplot."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"fastGWA file not found: {path}")

    # fastGWA is usually whitespace- or tab-delimited.
    df = pd.read_csv(path, sep=r"\s+", engine="python")

    # Standardize chromosome column name for qmplot.
    if "#CHROM" not in df.columns:
        if "CHR" in df.columns:
            df = df.rename(columns={"CHR": "#CHROM"})
        elif "CHROM" in df.columns:
            df = df.rename(columns={"CHROM": "#CHROM"})
        else:
            raise ValueError(
                f"Cannot find chromosome column in {path}. Columns are: {list(df.columns)}"
            )

    required = ["#CHROM", "POS", "P"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns {missing} in {path}. Columns are: {list(df.columns)}")

    # Convert key columns to numeric and clean invalid p-values.
    df["#CHROM"] = (
        df["#CHROM"].astype(str)
        .str.replace("chr", "", case=False, regex=False)
        .replace({"X": "23", "Y": "24", "MT": "25", "M": "25"})
    )
    df["#CHROM"] = pd.to_numeric(df["#CHROM"], errors="coerce")
    df["POS"] = pd.to_numeric(df["POS"], errors="coerce")
    df["P"] = pd.to_numeric(df["P"], errors="coerce")

    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=required).copy()
    df = df.loc[(df["#CHROM"] > 0) & (df["POS"] > 0) & (df["P"] >= 0) & (df["P"] <= 1)].copy()

    # Avoid infinite -log10(P) in Manhattan/QQ plots.
    df.loc[df["P"] < P_MIN, "P"] = P_MIN

    # Sort for stable plotting.
    df = df.sort_values(["#CHROM", "POS"]).reset_index(drop=True)
    df["#CHROM"] = df["#CHROM"].astype(int)

    return df


def plot_qmplot(output_dir, output_result, clock_name=None):
    os.makedirs(output_dir, exist_ok=True)
    df = read_fastgwa(output_result)

    if df.empty:
        raise ValueError(f"No valid variants after cleaning: {output_result}")

    label = clock_name if clock_name else os.path.basename(os.path.dirname(output_result))
    print(f"Plotting {label}")
    print(f"Input fastGWA: {output_result}")
    print(f"Valid variants: {df.shape[0]:,}")
    print(f"Minimum P value after clipping: {df['P'].min():.3e}")

    # Manhattan plot.
    manhattan_png = os.path.join(output_dir, "manhattan_qmplot.png")
    manhattan_pdf = os.path.join(output_dir, "manhattan_qmplot.pdf")

    plt.figure(figsize=(14, 6))
    manhattanplot(
        data=df,
        xticklabel_kws={"rotation": "vertical"},
        sign_marker_p=None,
        genomewideline=5e-8,
        suggestiveline=1e-5,
    )
    plt.title(label.replace("_", " "))
    plt.tight_layout()
    plt.savefig(manhattan_png, dpi=300, bbox_inches="tight")
    plt.savefig(manhattan_pdf, bbox_inches="tight")
    plt.close()

    # QQ plot.
    qq_png = os.path.join(output_dir, "QQ_plot.png")
    qq_pdf = os.path.join(output_dir, "QQ_plot.pdf")

    plt.figure(figsize=(6, 6))
    qqplot(
        data=df["P"],
        marker="o",
        xlabel=r"Expected $-\log_{10}{(P)}$",
        ylabel=r"Observed $-\log_{10}{(P)}$",
    )
    plt.title(label.replace("_", " "))
    plt.tight_layout()
    plt.savefig(qq_png, dpi=300, bbox_inches="tight")
    plt.savefig(qq_pdf, bbox_inches="tight")
    plt.close()

    # Save a minimal cleaned summary for audit.
    summary_tsv = os.path.join(output_dir, "qmplot_input_summary.tsv")
    pd.DataFrame(
        {
            "clock_name": [label],
            "fastgwa_file": [output_result],
            "n_variants_used": [df.shape[0]],
            "min_p": [df["P"].min()],
            "max_p": [df["P"].max()],
        }
    ).to_csv(summary_tsv, sep="\t", index=False)

    print(f"Saved: {manhattan_png}")
    print(f"Saved: {qq_png}")
    print("Done.")


def main():
    parser = argparse.ArgumentParser(description="Create Manhattan and QQ plots for GCTA fastGWA output.")
    parser.add_argument("--output_dir", required=True, help="Directory to save Manhattan/QQ plots")
    parser.add_argument("--output_result", required=True, help="Path to organ_pheno_normalized_residualized.fastGWA")
    parser.add_argument("--clock_name", default=None, help="Optional clock name for plot title")
    args = parser.parse_args()

    plot_qmplot(args.output_dir, args.output_result, args.clock_name)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        raise
