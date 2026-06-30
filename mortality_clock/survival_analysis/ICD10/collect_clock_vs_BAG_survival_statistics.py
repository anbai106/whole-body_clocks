#!/usr/bin/env python3
"""
Collect and summarize comparative survival-analysis results for mortality clocks versus BAGs.

Input: one TSV per ICD endpoint, e.g.
  cox_compare_clock_vs_BAG_I10.tsv

Each TSV is expected to contain one row per clock/BAG pair with columns produced by
survival_analysis_clock_vs_bag.py, including clock_p, bag_p, joint_p_diff,
C-index columns, N_case, status, etc.

Outputs:
  1. Combined raw table across all endpoints and all statuses
  2. Powered table restricted to N_case >= threshold and status == ok
  3. Bonferroni-annotated significant results for clock, BAG, and clock-vs-BAG difference
  4. Disease-level summary
  5. Modality-level summary
  6. Pair-level summary
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd



# =========================
# Hardcoded input settings
# =========================
# Edit these paths here if you move the project between CUBIC and local Mac.
ICD_LIST = "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/mortality_clock/SA/data/included_ICD_mortality_clock.tsv"
INPUT_DIR = "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/mortality_clock/SA/output_clock_vs_BAG"
OUTPUT_DIR = "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/Result/clock_vs_BAG_survival_summary"
FILE_PREFIX = "cox_compare_clock_vs_BAG_"
MIN_CASE = 50
ALPHA = 0.05
WRITE_EXCEL = True
ALLOW_GLOB_IF_NO_ICD_LIST = False

# Optional local Mac settings. Uncomment these lines if running locally.
# ICD_LIST = "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/mortality_clock/SA/data/included_ICD_mortality_clock.tsv"
# INPUT_DIR = "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/mortality_clock/SA/output_clock_vs_BAG"
# OUTPUT_DIR = "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/Result/clock_vs_BAG_survival_summary"

def safe_read_tsv(path: Path) -> pd.DataFrame:
    """Read a TSV robustly, preserving empty strings as NaN."""
    return pd.read_csv(path, sep="\t", low_memory=False)


def read_icd_list(icd_list_path: Path, input_dir: Path, file_prefix: str, allow_glob: bool) -> pd.DataFrame:
    """Return a DataFrame with at least disease_id and optional included_count."""
    if icd_list_path.exists():
        df = pd.read_csv(icd_list_path, sep="\t")
        if "value" not in df.columns:
            raise ValueError(f"ICD list must contain a 'value' column: {icd_list_path}")
        out = df.rename(columns={"value": "disease_id", "count": "included_count"}).copy()
        keep = ["disease_id"] + (["included_count"] if "included_count" in out.columns else [])
        return out[keep]

    if not allow_glob:
        raise FileNotFoundError(
            f"ICD list not found: {icd_list_path}. Use --allow_glob_if_no_icd_list to collect by glob."
        )

    files = sorted(input_dir.glob(f"{file_prefix}*.tsv"))
    disease_ids = [f.name.replace(file_prefix, "").replace(".tsv", "") for f in files]
    return pd.DataFrame({"disease_id": disease_ids})


def bh_fdr(p: Iterable[float]) -> np.ndarray:
    """Benjamini-Hochberg FDR q-values, implemented without statsmodels."""
    p = np.asarray(list(p), dtype=float)
    q = np.full_like(p, np.nan, dtype=float)
    finite = np.isfinite(p)
    if finite.sum() == 0:
        return q
    pv = p[finite]
    order = np.argsort(pv)
    ranked = pv[order]
    m = len(ranked)
    raw_q = ranked * m / np.arange(1, m + 1)
    raw_q = np.minimum.accumulate(raw_q[::-1])[::-1]
    raw_q = np.minimum(raw_q, 1.0)
    q_finite = np.empty_like(raw_q)
    q_finite[order] = raw_q
    q[finite] = q_finite
    return q


def add_derived_statistics(df: pd.DataFrame, alpha: float) -> pd.DataFrame:
    """Add comparison metrics, Bonferroni thresholds, and FDR q-values."""
    df = df.copy()

    # Numeric conversion for key columns. Errors remain NaN.
    numeric_cols = [
        "N", "N_case", "N_noncase", "event_rate",
        "followup_years_min", "followup_years_max",
        "event_followup_years_min", "event_followup_years_max",
        "clock_bag_pearson",
        "clock_beta", "clock_se", "clock_hr", "clock_ci_lo", "clock_ci_hi", "clock_p",
        "bag_beta", "bag_se", "bag_hr", "bag_ci_lo", "bag_ci_hi", "bag_p",
        "clock_joint_beta", "clock_joint_se", "clock_joint_hr", "clock_joint_p",
        "bag_joint_beta", "bag_joint_se", "bag_joint_hr", "bag_joint_p",
        "joint_beta_diff_clock_minus_bag", "joint_se_diff", "joint_z_diff", "joint_p_diff",
        "base_cindex", "clock_cindex", "bag_cindex", "both_cindex",
        "delta_cindex_clock_minus_bag", "delta_cindex_clock_minus_base",
        "delta_cindex_bag_minus_base", "delta_cindex_both_minus_base",
        "lrt_p_clock_vs_base", "lrt_p_bag_vs_base", "lrt_p_both_vs_base",
        "included_count",
    ]
    for c in numeric_cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    # Core comparison metrics.
    if {"clock_hr", "bag_hr"}.issubset(df.columns):
        df["hr_ratio_clock_over_bag"] = df["clock_hr"] / df["bag_hr"]
        df["hr_diff_clock_minus_bag"] = df["clock_hr"] - df["bag_hr"]

    if {"clock_beta", "bag_beta"}.issubset(df.columns):
        df["abs_beta_clock"] = df["clock_beta"].abs()
        df["abs_beta_bag"] = df["bag_beta"].abs()
        df["abs_beta_diff_clock_minus_bag"] = df["abs_beta_clock"] - df["abs_beta_bag"]
        df["clock_has_larger_abs_beta"] = df["abs_beta_diff_clock_minus_bag"] > 0

    if "delta_cindex_clock_minus_bag" in df.columns:
        df["clock_has_higher_cindex"] = df["delta_cindex_clock_minus_bag"] > 0

    if {"clock_joint_p", "bag_joint_p"}.issubset(df.columns):
        df["clock_joint_sig_nominal"] = df["clock_joint_p"] <= alpha
        df["bag_joint_sig_nominal"] = df["bag_joint_p"] <= alpha

    if "joint_p_diff" in df.columns:
        df["joint_diff_sig_nominal"] = df["joint_p_diff"] <= alpha
        df["joint_diff_favors_clock"] = (df["joint_beta_diff_clock_minus_bag"] > 0) & (df["joint_p_diff"] <= alpha)
        df["joint_diff_favors_bag"] = (df["joint_beta_diff_clock_minus_bag"] < 0) & (df["joint_p_diff"] <= alpha)

    # Multiple-testing correction for powered ok rows is applied after filtering.
    ok_mask = (df.get("status", pd.Series(index=df.index, dtype=object)) == "ok")
    ok_mask &= df.get("N_case", pd.Series(index=df.index, dtype=float)).notna()

    n_diseases = df.loc[ok_mask, "disease_id"].nunique() if "disease_id" in df.columns else 0
    n_pairs = df.loc[ok_mask, "pair_id"].nunique() if "pair_id" in df.columns else 0
    n_tests = max(int(n_diseases) * int(n_pairs), 1)
    bonf = alpha / n_tests
    df["n_diseases_for_bonferroni_all_ok"] = n_diseases
    df["n_pairs_for_bonferroni_all_ok"] = n_pairs
    df["bonferroni_p_all_ok"] = bonf

    for p_col in ["clock_p", "bag_p", "clock_joint_p", "bag_joint_p", "joint_p_diff", "lrt_p_clock_vs_base", "lrt_p_bag_vs_base", "lrt_p_both_vs_base"]:
        if p_col in df.columns:
            q_col = p_col.replace("_p", "_q_bh") if p_col.endswith("_p") else f"{p_col}_q_bh"
            sig_col = p_col.replace("_p", "_sig_bonf") if p_col.endswith("_p") else f"{p_col}_sig_bonf"
            df[q_col] = np.nan
            df[sig_col] = False
            pvals = df.loc[ok_mask, p_col]
            df.loc[ok_mask, q_col] = bh_fdr(pvals)
            df.loc[ok_mask, sig_col] = pvals <= bonf

    return df


def summarize_diseases(df_powered: pd.DataFrame) -> pd.DataFrame:
    """Disease-level summary after N_case/status filtering."""
    if df_powered.empty:
        return pd.DataFrame()

    def agg_one(g: pd.DataFrame) -> pd.Series:
        return pd.Series({
            "n_pairs_ok": g["pair_id"].nunique(),
            "n_rows_ok": len(g),
            "min_N_case": g["N_case"].min(),
            "max_N_case": g["N_case"].max(),
            "median_N_case": g["N_case"].median(),
            "mean_delta_cindex_clock_minus_bag": g["delta_cindex_clock_minus_bag"].mean(),
            "median_delta_cindex_clock_minus_bag": g["delta_cindex_clock_minus_bag"].median(),
            "n_clock_higher_cindex": int((g["delta_cindex_clock_minus_bag"] > 0).sum()),
            "prop_clock_higher_cindex": float((g["delta_cindex_clock_minus_bag"] > 0).mean()),
            "mean_abs_beta_diff_clock_minus_bag": g["abs_beta_diff_clock_minus_bag"].mean() if "abs_beta_diff_clock_minus_bag" in g else np.nan,
            "n_clock_larger_abs_beta": int((g["abs_beta_diff_clock_minus_bag"] > 0).sum()) if "abs_beta_diff_clock_minus_bag" in g else np.nan,
            "prop_clock_larger_abs_beta": float((g["abs_beta_diff_clock_minus_bag"] > 0).mean()) if "abs_beta_diff_clock_minus_bag" in g else np.nan,
            "n_clock_sig_bonf": int(g.get("clock_sig_bonf", pd.Series(False, index=g.index)).sum()),
            "n_bag_sig_bonf": int(g.get("bag_sig_bonf", pd.Series(False, index=g.index)).sum()),
            "n_joint_diff_sig_bonf": int(g.get("joint_p_diff_sig_bonf_powered", g.get("joint_p_diff_sig_bonf", pd.Series(False, index=g.index))).sum()),
            "n_joint_diff_favors_clock_nominal": int(g.get("joint_diff_favors_clock", pd.Series(False, index=g.index)).sum()),
            "n_joint_diff_favors_bag_nominal": int(g.get("joint_diff_favors_bag", pd.Series(False, index=g.index)).sum()),
            "n_joint_diff_favors_clock_bonf": int(((g.get("joint_beta_diff_clock_minus_bag", pd.Series(np.nan, index=g.index)) > 0) & g.get("joint_p_diff_sig_bonf_powered", g.get("joint_p_diff_sig_bonf", pd.Series(False, index=g.index)))).sum()),
            "n_joint_diff_favors_bag_bonf": int(((g.get("joint_beta_diff_clock_minus_bag", pd.Series(np.nan, index=g.index)) < 0) & g.get("joint_p_diff_sig_bonf_powered", g.get("joint_p_diff_sig_bonf", pd.Series(False, index=g.index)))).sum()),
            "disease_level_wilcoxon_cindex_p": g["disease_level_wilcoxon_cindex_p"].dropna().iloc[0] if "disease_level_wilcoxon_cindex_p" in g and g["disease_level_wilcoxon_cindex_p"].notna().any() else np.nan,
            "disease_level_wilcoxon_abs_beta_p": g["disease_level_wilcoxon_abs_beta_p"].dropna().iloc[0] if "disease_level_wilcoxon_abs_beta_p" in g and g["disease_level_wilcoxon_abs_beta_p"].notna().any() else np.nan,
        })

    out = df_powered.groupby("disease_id", dropna=False).apply(agg_one).reset_index()
    if "included_count" in df_powered.columns:
        counts = df_powered.groupby("disease_id", dropna=False)["included_count"].first().reset_index()
        out = counts.merge(out, on="disease_id", how="right")
    return out.sort_values(["mean_delta_cindex_clock_minus_bag", "disease_id"], ascending=[False, True])


def summarize_by_group(df_powered: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    """Generic summary by modality, pair_id, organ, etc."""
    if df_powered.empty:
        return pd.DataFrame()
    for c in group_cols:
        if c not in df_powered.columns:
            return pd.DataFrame()

    g = df_powered.groupby(group_cols, dropna=False)
    out = g.agg(
        n_diseases=("disease_id", "nunique"),
        n_rows=("disease_id", "size"),
        median_N_case=("N_case", "median"),
        mean_delta_cindex_clock_minus_bag=("delta_cindex_clock_minus_bag", "mean"),
        median_delta_cindex_clock_minus_bag=("delta_cindex_clock_minus_bag", "median"),
        prop_clock_higher_cindex=("clock_has_higher_cindex", "mean"),
        mean_clock_hr=("clock_hr", "mean"),
        mean_bag_hr=("bag_hr", "mean"),
        mean_abs_beta_diff_clock_minus_bag=("abs_beta_diff_clock_minus_bag", "mean"),
        prop_clock_larger_abs_beta=("clock_has_larger_abs_beta", "mean"),
        prop_joint_diff_favors_clock_nominal=("joint_diff_favors_clock", "mean"),
        prop_joint_diff_favors_bag_nominal=("joint_diff_favors_bag", "mean"),
    ).reset_index()
    return out.sort_values("mean_delta_cindex_clock_minus_bag", ascending=False)


def write_table(df: pd.DataFrame, path_prefix: Path, write_excel: bool = False) -> None:
    path_prefix.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path_prefix.with_suffix(".tsv"), sep="\t", index=False)
    if write_excel:
        try:
            df.to_excel(path_prefix.with_suffix(".xlsx"), index=False)
        except Exception as e:
            print(f"WARNING: failed to write Excel file for {path_prefix.name}: {e}")


def collect_data(
    icd_list: str,
    input_dir: str,
    output_dir: str,
    file_prefix: str = "cox_compare_clock_vs_BAG_",
    min_case: int = 50,
    alpha: float = 0.05,
    write_excel: bool = False,
    allow_glob_if_no_icd_list: bool = False,
) -> None:
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    icd_list_path = Path(icd_list)
    output_dir.mkdir(parents=True, exist_ok=True)

    icd_df = read_icd_list(icd_list_path, input_dir, file_prefix, allow_glob_if_no_icd_list)
    icd_df["disease_id"] = icd_df["disease_id"].astype(str)

    rows: list[pd.DataFrame] = []
    missing_files: list[dict] = []

    for _, r in icd_df.iterrows():
        disease_id = str(r["disease_id"])
        tsv = input_dir / f"{file_prefix}{disease_id}.tsv"
        if not tsv.exists() or tsv.stat().st_size == 0:
            missing_files.append({
                "disease_id": disease_id,
                "included_count": r.get("included_count", np.nan),
                "expected_file": str(tsv),
                "reason": "missing_or_empty",
            })
            print(f"Missing or empty result file: {disease_id} -> {tsv}")
            continue

        df = safe_read_tsv(tsv)
        if "disease_id" not in df.columns:
            df["disease_id"] = disease_id
        if "included_count" in icd_df.columns:
            df["included_count"] = r.get("included_count", np.nan)
        df["source_file"] = str(tsv)
        rows.append(df)

    if not rows:
        raise RuntimeError("No result TSV files were collected. Check input_dir and icd_list.")

    df_all = pd.concat(rows, ignore_index=True, sort=False)
    df_all = add_derived_statistics(df_all, alpha=alpha)

    # Primary powered table: only successful Cox models with enough cases.
    df_powered = df_all.loc[(df_all["status"] == "ok") & (df_all["N_case"] >= min_case)].copy()

    # Recompute Bonferroni within the powered table, which is usually the intended analysis family.
    n_de_powered = df_powered["disease_id"].nunique()
    n_pair_powered = df_powered["pair_id"].nunique()
    n_tests_powered = max(int(n_de_powered) * int(n_pair_powered), 1)
    p_bonf_powered = alpha / n_tests_powered

    df_powered["n_diseases_powered"] = n_de_powered
    df_powered["n_pairs_powered"] = n_pair_powered
    df_powered["bonferroni_p_powered"] = p_bonf_powered

    # Bonferroni and FDR in powered set for the key p-value families.
    p_families = [
        "clock_p", "bag_p", "clock_joint_p", "bag_joint_p", "joint_p_diff",
        "lrt_p_clock_vs_base", "lrt_p_bag_vs_base", "lrt_p_both_vs_base",
    ]
    for p_col in p_families:
        if p_col in df_powered.columns:
            stem = p_col.replace("_p", "") if p_col.endswith("_p") else p_col
            df_powered[f"{stem}_q_bh_powered"] = bh_fdr(df_powered[p_col])
            df_powered[f"{stem}_sig_bonf_powered"] = df_powered[p_col] <= p_bonf_powered

    # Main significant result table: any of the key comparison criteria.
    sig_masks = []
    for col in [
        "clock_sig_bonf_powered", "bag_sig_bonf_powered",
        "clock_joint_sig_bonf_powered", "bag_joint_sig_bonf_powered",
        "joint_diff_sig_bonf_powered",
        "lrt_p_clock_vs_base_sig_bonf_powered", "lrt_p_bag_vs_base_sig_bonf_powered",
        "lrt_p_both_vs_base_sig_bonf_powered",
    ]:
        if col in df_powered.columns:
            sig_masks.append(df_powered[col].fillna(False))

    if sig_masks:
        any_sig = np.logical_or.reduce(sig_masks)
        df_sig = df_powered.loc[any_sig].copy()
    else:
        df_sig = df_powered.iloc[0:0].copy()

    # More specific table for L'EPOCH superiority: clock better than BAG by C-index and/or joint beta test.
    superiority_mask = pd.Series(False, index=df_powered.index)
    if "delta_cindex_clock_minus_bag" in df_powered.columns:
        superiority_mask |= df_powered["delta_cindex_clock_minus_bag"] > 0
    if "joint_beta_diff_clock_minus_bag" in df_powered.columns and "joint_diff_sig_bonf_powered" in df_powered.columns:
        superiority_mask |= (df_powered["joint_beta_diff_clock_minus_bag"] > 0) & df_powered["joint_diff_sig_bonf_powered"].fillna(False)
    df_clock_better = df_powered.loc[superiority_mask].copy()

    missing_df = pd.DataFrame(missing_files)
    disease_summary = summarize_diseases(df_powered)
    modality_summary = summarize_by_group(df_powered, ["modality"])
    pair_summary = summarize_by_group(df_powered, ["pair_id", "organ", "modality"])
    organ_summary = summarize_by_group(df_powered, ["organ"])

    # Write outputs.
    write_table(df_all, output_dir / "clock_vs_BAG_all_rows_all_status", write_excel=write_excel)
    write_table(df_powered, output_dir / f"clock_vs_BAG_powered_ok_Ncase_ge_{min_case}", write_excel=write_excel)
    write_table(df_sig, output_dir / f"clock_vs_BAG_significant_any_test_Ncase_ge_{min_case}", write_excel=write_excel)
    write_table(df_clock_better, output_dir / f"clock_vs_BAG_clock_better_rows_Ncase_ge_{min_case}", write_excel=write_excel)
    write_table(disease_summary, output_dir / f"clock_vs_BAG_disease_level_summary_Ncase_ge_{min_case}", write_excel=write_excel)
    write_table(modality_summary, output_dir / f"clock_vs_BAG_modality_summary_Ncase_ge_{min_case}", write_excel=write_excel)
    write_table(pair_summary, output_dir / f"clock_vs_BAG_pair_summary_Ncase_ge_{min_case}", write_excel=write_excel)
    write_table(organ_summary, output_dir / f"clock_vs_BAG_organ_summary_Ncase_ge_{min_case}", write_excel=write_excel)
    if not missing_df.empty:
        write_table(missing_df, output_dir / "clock_vs_BAG_missing_or_empty_files", write_excel=write_excel)

    # Console summary.
    print("\n=== Collection summary ===")
    print(f"Input dir: {input_dir}")
    print(f"Output dir: {output_dir}")
    print(f"Requested ICD endpoints: {icd_df['disease_id'].nunique():,}")
    print(f"Missing/empty result files: {len(missing_files):,}")
    print(f"Collected rows: {len(df_all):,}")
    print(f"Collected diseases: {df_all['disease_id'].nunique():,}")
    print(f"Powered ok rows, N_case >= {min_case}: {len(df_powered):,}")
    print(f"Powered diseases: {n_de_powered:,}")
    print(f"Powered pairs: {n_pair_powered:,}")
    print(f"Bonferroni threshold in powered analysis: {p_bonf_powered:.3e}")
    print(f"Any significant rows: {len(df_sig):,}")
    print(f"Clock-better rows: {len(df_clock_better):,}")
    print("Done.")


def main() -> None:
    print("Running hardcoded clock-vs-BAG survival-statistics collector")
    print(f"ICD_LIST: {ICD_LIST}")
    print(f"INPUT_DIR: {INPUT_DIR}")
    print(f"OUTPUT_DIR: {OUTPUT_DIR}")
    print(f"MIN_CASE: {MIN_CASE}")

    collect_data(
        icd_list=ICD_LIST,
        input_dir=INPUT_DIR,
        output_dir=OUTPUT_DIR,
        file_prefix=FILE_PREFIX,
        min_case=MIN_CASE,
        alpha=ALPHA,
        write_excel=WRITE_EXCEL,
        allow_glob_if_no_icd_list=ALLOW_GLOB_IF_NO_ICD_LIST,
    )


if __name__ == "__main__":
    main()
