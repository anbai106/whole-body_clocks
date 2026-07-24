#!/usr/bin/env python3
"""
STEP 5 CLEAN: Final MHAS mortality EPOCH report using ONLY baseline disease-burden associations

This script intentionally removes all old incident/future-disease wording and outputs.
It only reads the updated Step 4 baseline prevalent disease results.

Disease-analysis design expected from Step 4:
  baseline disease at Wave 1 ~ mortality EPOCH score + covariates

Required Step 4 input:
  /Users/hao/Dropbox/MHAS/step4_epoch_disease_associations/
    mhas_step4_disease_association_summary.tsv
    mhas_step4_disease_or_by_quartile.tsv
    mhas_step4_disease_event_rate_by_quartile.tsv

Expected Step 4 analysis_type:
  baseline_prevalent

Outputs:
  /Users/hao/Dropbox/MHAS/step5_final_report/
    mhas_step5_final_report_baseline_burden.md
    mhas_step5_table1_mortality_performance.tsv
    mhas_step5_table2_mortality_hr.tsv
    mhas_step5_table3_baseline_burden_or_per_sd.tsv
    mhas_step5_table4_baseline_burden_q4_vs_q1.tsv
    mhas_step5_table5_baseline_burden_prevalence_by_quartile.tsv
    mhas_step5_top_selected_features.tsv
    mhas_step5_fig1_mortality_cindex.pdf/png
    mhas_step5_fig2_mortality_hr.pdf/png
    mhas_step5_fig3_baseline_burden_or_per_sd.pdf/png
    mhas_step5_fig4_baseline_burden_q4_vs_q1.pdf/png
    mhas_step5_fig5_top_baseline_burden_prevalence_by_quartile.pdf/png
    mhas_step5_audit_baseline_burden.txt

Important:
  This script does NOT create legacy filenames such as:
    mhas_step5_fig_disease_forest_final.pdf
    mhas_step5_table3_top_disease_associations.tsv

Run:
  cd /Users/hao/Project/whole-body_clocks/MHAS

  /Users/hao/opt/anaconda3/envs/DNE/bin/python \
    5_step5_mhas_final_report_baseline_burden_only.py \
    --primary-split all

Recommended workflow:
  rm -rf /Users/hao/Dropbox/MHAS/step5_final_report
  /Users/hao/opt/anaconda3/envs/DNE/bin/python 5_step5_mhas_final_report_baseline_burden_only.py --primary-split all
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------------
# Utilities
# -----------------------------------------------------------------------------
def log(msg: str) -> None:
    print(msg, flush=True)


def safe_mkdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def read_tsv(path: Path, required=False):
    if path.exists():
        return pd.read_csv(path, sep="\t", low_memory=False)
    if required:
        raise FileNotFoundError(f"Required file not found: {path}")
    log(f"Optional file not found: {path}")
    return pd.DataFrame()


def fmt_num(x, digits=3):
    if pd.isna(x):
        return "NA"
    try:
        return f"{float(x):.{digits}f}"
    except Exception:
        return str(x)


def fmt_sci(x):
    if pd.isna(x):
        return "NA"
    try:
        x = float(x)
        if x == 0:
            return "0"
        if abs(x) < 0.001:
            return f"{x:.2e}"
        return f"{x:.3g}"
    except Exception:
        return str(x)


def pretty_endpoint(x):
    return (
        str(x)
        .replace("_incl_", "_including_")
        .replace("_", " ")
        .replace("respiratory disease including asthma", "respiratory disease/asthma")
    )


def prefer_adjusted(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty or "adjusted" not in df.columns:
        return df
    adjusted_bool = df["adjusted"].fillna(False).astype(bool)
    if adjusted_bool.any():
        return df[adjusted_bool].copy()
    return df[~adjusted_bool].copy()


def keep_valid_status(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    if "status" in df.columns:
        return df[df["status"].astype(str).str.startswith("ok")].copy()
    return df


def require_baseline_prevalent(df: pd.DataFrame, name: str) -> pd.DataFrame:
    """
    Keep only Step 4 baseline-prevalent rows. If analysis_type is absent, fail loudly
    because this script is meant to avoid accidentally using old incident results.
    """
    if df.empty:
        return df
    if "analysis_type" not in df.columns:
        raise ValueError(
            f"{name} does not contain an analysis_type column. "
            "This may be an old Step 4 output. Re-run updated Step 4 baseline "
            "prevalent disease script before Step 5."
        )
    out = df[df["analysis_type"].astype(str) == "baseline_prevalent"].copy()
    if out.empty:
        raise ValueError(
            f"{name} contains no analysis_type == 'baseline_prevalent' rows. "
            "This looks like old incident/future-disease output."
        )
    return out


def or_text(row, estimate_col, lo_col="ci_lower", hi_col="ci_upper", p_col="p"):
    return (
        f"{fmt_num(row.get(estimate_col), 2)} "
        f"[{fmt_num(row.get(lo_col), 2)}, {fmt_num(row.get(hi_col), 2)}], "
        f"P={fmt_sci(row.get(p_col))}"
    )


def hr_text(row, estimate_col="hr", lo_col="ci_lower", hi_col="ci_upper", p_col="p"):
    return (
        f"{fmt_num(row.get(estimate_col), 2)} "
        f"[{fmt_num(row.get(lo_col), 2)}, {fmt_num(row.get(hi_col), 2)}], "
        f"P={fmt_sci(row.get(p_col))}"
    )


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------
def save_plot(fig, out_prefix):
    fig.tight_layout()
    fig.savefig(f"{out_prefix}.pdf")
    fig.savefig(f"{out_prefix}.png", dpi=300)
    plt.close(fig)
    log(f"Saved: {out_prefix}.pdf")
    log(f"Saved: {out_prefix}.png")


def plot_mortality_cindex(cindex_df, out_prefix):
    if cindex_df.empty:
        return

    score_order = ["clinical_baseline_lp", "lp_total", "mortality_epoch_acceleration_z"]
    score_labels = {
        "clinical_baseline_lp": "Clinical baseline",
        "lp_total": "EPOCH LP",
        "mortality_epoch_acceleration_z": "EPOCH acceleration",
    }
    split_order = ["train", "validation", "test"]

    sub = cindex_df[
        cindex_df["split"].isin(split_order)
        & cindex_df["score"].isin(score_order)
    ].copy()

    if sub.empty:
        return

    rows = []
    for split in split_order:
        for score in score_order:
            r = sub[(sub["split"] == split) & (sub["score"] == score)]
            if len(r):
                rr = r.iloc[0]
                rows.append({
                    "label": f"{split}\n{score_labels.get(score, score)}",
                    "cindex": rr["cindex"],
                    "ci_lower": rr["ci_lower"],
                    "ci_upper": rr["ci_upper"],
                })

    pdat = pd.DataFrame(rows)
    if pdat.empty:
        return

    x = np.arange(len(pdat))
    y = pdat["cindex"].astype(float).values
    lo = pdat["ci_lower"].astype(float).values
    hi = pdat["ci_upper"].astype(float).values

    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    ax.errorbar(x, y, yerr=[y - lo, hi - y], fmt="o", capsize=3)
    ax.set_xticks(x)
    ax.set_xticklabels(pdat["label"], rotation=45, ha="right")
    ax.set_ylabel("Harrell C-index")
    ax.set_title("MHAS mortality EPOCH discrimination")
    ax.grid(True, axis="y", alpha=0.25)
    ymin = max(0.45, np.nanmin(lo) - 0.05)
    ymax = min(1.00, np.nanmax(hi) + 0.05)
    ax.set_ylim(ymin, ymax)
    save_plot(fig, out_prefix)


def plot_mortality_hr(hr_df, out_prefix):
    if hr_df.empty:
        return

    sub = hr_df[
        (hr_df["score"] == "mortality_epoch_acceleration_z")
        & hr_df["split"].isin(["train", "validation", "test", "all"])
    ].copy()

    sub = prefer_adjusted(sub)
    sub = keep_valid_status(sub)
    sub = sub.dropna(subset=["hr", "ci_lower", "ci_upper"])
    if sub.empty:
        return

    order = ["train", "validation", "test", "all"]
    sub["split"] = pd.Categorical(sub["split"], categories=order, ordered=True)
    sub = sub.sort_values("split")

    y = np.arange(len(sub))
    x = sub["hr"].astype(float).values
    lo = sub["ci_lower"].astype(float).values
    hi = sub["ci_upper"].astype(float).values

    fig, ax = plt.subplots(figsize=(6.8, 4.5))
    ax.errorbar(x, y, xerr=[x - lo, hi - x], fmt="o", capsize=3)
    ax.axvline(1.0, linestyle="--", linewidth=1)
    ax.set_yticks(y)
    ax.set_yticklabels(sub["split"].astype(str))
    ax.set_xlabel("Mortality HR per 1-SD higher EPOCH acceleration")
    ax.set_title("Mortality association of EPOCH acceleration")
    ax.grid(True, axis="x", alpha=0.25)
    save_plot(fig, out_prefix)


def plot_baseline_burden_or(per_sd_df, out_prefix, primary_split="all", max_endpoints=20):
    if per_sd_df.empty:
        return

    sub = per_sd_df.copy()
    sub = sub[
        (sub["split"] == primary_split)
        & (sub["score"] == "mortality_epoch_acceleration_z")
    ].copy()
    sub = keep_valid_status(sub)
    sub = prefer_adjusted(sub)
    sub = sub.dropna(subset=["or_per_sd", "ci_lower", "ci_upper", "p"])

    if sub.empty:
        return

    sub = sub.sort_values("p", na_position="last").head(max_endpoints)
    sub = sub.sort_values("or_per_sd", ascending=True)
    sub["endpoint_label"] = sub["endpoint"].map(pretty_endpoint)

    y = np.arange(len(sub))
    x = sub["or_per_sd"].astype(float).values
    lo = sub["ci_lower"].astype(float).values
    hi = sub["ci_upper"].astype(float).values

    fig_h = max(4.8, 0.45 * len(sub) + 1.5)
    fig, ax = plt.subplots(figsize=(8.2, fig_h))
    ax.errorbar(x, y, xerr=[x - lo, hi - x], fmt="o", capsize=3)
    ax.axvline(1.0, linestyle="--", linewidth=1)
    ax.set_yticks(y)
    ax.set_yticklabels(sub["endpoint_label"])
    ax.set_xlabel("Odds ratio per 1-SD higher EPOCH acceleration")
    ax.set_title("Baseline disease-burden associations")
    ax.grid(True, axis="x", alpha=0.25)
    save_plot(fig, out_prefix)


def plot_baseline_burden_q4(quartile_df, out_prefix, primary_split="all", max_endpoints=20):
    if quartile_df.empty:
        return

    sub = quartile_df.copy()
    sub = sub[
        (sub["split"] == primary_split)
        & (sub["comparison"] == "Q4_highest vs Q1_lowest")
    ].copy()
    sub = keep_valid_status(sub)
    sub = prefer_adjusted(sub)
    sub = sub.dropna(subset=["or", "ci_lower", "ci_upper", "p"])

    if sub.empty:
        return

    sub = sub.sort_values("p", na_position="last").head(max_endpoints)
    sub = sub.sort_values("or", ascending=True)
    sub["endpoint_label"] = sub["endpoint"].map(pretty_endpoint)

    y = np.arange(len(sub))
    x = sub["or"].astype(float).values
    lo = sub["ci_lower"].astype(float).values
    hi = sub["ci_upper"].astype(float).values

    fig_h = max(4.8, 0.45 * len(sub) + 1.5)
    fig, ax = plt.subplots(figsize=(8.2, fig_h))
    ax.errorbar(x, y, xerr=[x - lo, hi - x], fmt="o", capsize=3)
    ax.axvline(1.0, linestyle="--", linewidth=1)
    ax.set_yticks(y)
    ax.set_yticklabels(sub["endpoint_label"])
    ax.set_xlabel("Odds ratio for Q4 highest vs Q1 lowest")
    ax.set_title("Baseline disease burden in highest vs lowest EPOCH quartile")
    ax.grid(True, axis="x", alpha=0.25)
    save_plot(fig, out_prefix)


def plot_top_prevalence_by_quartile(rate_df, per_sd_df, out_prefix, primary_split="all"):
    if rate_df.empty or per_sd_df.empty:
        return

    assoc = per_sd_df[
        (per_sd_df["split"] == primary_split)
        & (per_sd_df["score"] == "mortality_epoch_acceleration_z")
    ].copy()
    assoc = keep_valid_status(assoc)
    assoc = prefer_adjusted(assoc)
    assoc = assoc.dropna(subset=["p"])
    if assoc.empty:
        return

    top_endpoint = assoc.sort_values("p").iloc[0]["endpoint"]
    sub = rate_df[
        (rate_df["split"] == primary_split)
        & (rate_df["endpoint"] == top_endpoint)
    ].copy()

    if sub.empty:
        return

    order = ["Q1_lowest", "Q2", "Q3", "Q4_highest"]
    sub["quartile"] = pd.Categorical(sub["quartile"], categories=order, ordered=True)
    sub = sub.sort_values("quartile")

    fig, ax = plt.subplots(figsize=(6.8, 4.8))
    ax.bar(sub["quartile"].astype(str), sub["event_rate"].astype(float))
    ax.set_xlabel("Mortality EPOCH acceleration quartile")
    ax.set_ylabel("Baseline disease prevalence")
    ax.set_title(f"Baseline {pretty_endpoint(top_endpoint)} prevalence by EPOCH quartile")
    ymax = min(1.0, max(0.05, float(sub["event_rate"].max()) * 1.25))
    ax.set_ylim(0, ymax)
    ax.grid(True, axis="y", alpha=0.25)
    save_plot(fig, out_prefix)


# -----------------------------------------------------------------------------
# Tables and report
# -----------------------------------------------------------------------------
def make_markdown_report(out_path, table1, table2, table3, table4, table5,
                         selected_features, paths, primary_split="all", top_n=10):
    lines = []
    lines.append("# MHAS phenotype-based mortality EPOCH final report\n")
    lines.append("## Workflow summary\n")
    lines.append("This report summarizes the MHAS phenotype-based mortality EPOCH analysis.\n")
    lines.append("- Step 2 trained the elastic-net Cox mortality EPOCH model.\n")
    lines.append("- Step 3 validated mortality prediction and survival associations.\n")
    lines.append("- Step 4 tested **baseline disease-burden associations** using Wave-1 prevalent disease labels.\n")
    lines.append("- No incident/future-disease association results are included in this report.\n")

    lines.append("\n## Disease-analysis design\n")
    lines.append("Model: `baseline disease at Wave 1 ~ mortality EPOCH score + covariates`.\n")
    lines.append(f"Primary split for baseline disease-burden summaries: **{primary_split}**.\n")

    if not table1.empty:
        lines.append("\n## Table 1. Mortality model performance\n")
        lines.append(table1.to_markdown(index=False))
        lines.append("\n")

    if not table2.empty:
        lines.append("\n## Table 2. Mortality association summary\n")
        lines.append(table2.to_markdown(index=False))
        lines.append("\n")

    if not table3.empty:
        lines.append("\n## Table 3. Baseline disease-burden associations per 1-SD EPOCH acceleration\n")
        lines.append(table3.head(top_n).to_markdown(index=False))
        lines.append("\n")

    if not table4.empty:
        lines.append("\n## Table 4. Baseline disease burden in highest vs lowest EPOCH quartile\n")
        lines.append(table4.head(top_n).to_markdown(index=False))
        lines.append("\n")

    if not table5.empty:
        lines.append("\n## Table 5. Baseline prevalence by EPOCH quartile for the top endpoint\n")
        lines.append(table5.to_markdown(index=False))
        lines.append("\n")

    if not selected_features.empty:
        lines.append("\n## Top selected EPOCH model features\n")
        lines.append(selected_features.head(top_n).to_markdown(index=False))
        lines.append("\n")

    lines.append("\n## Output files\n")
    for label, path in paths.items():
        lines.append(f"- **{label}:** `{path}`")

    out_path.write_text("\n".join(lines))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--base-dir",
        default="/Users/hao/Dropbox/MHAS",
        help="Base MHAS directory."
    )
    parser.add_argument(
        "--out-dir",
        default="/Users/hao/Dropbox/MHAS/step5_final_report",
        help="Output directory. Delete old output folder before running if desired."
    )
    parser.add_argument(
        "--primary-split",
        default="all",
        choices=["train", "validation", "test", "all"],
        help="Primary split for baseline disease-burden tables. Default: all."
    )
    parser.add_argument(
        "--top-disease-n",
        type=int,
        default=20,
        help="Number of disease endpoints to keep in final disease tables/figures."
    )
    parser.add_argument("--skip-plots", action="store_true")
    args = parser.parse_args()

    base_dir = Path(args.base_dir)
    out_dir = Path(args.out_dir)
    safe_mkdir(out_dir)

    step2_dir = base_dir / "step2_mortality_epoch_model"
    step3_dir = base_dir / "step3_mortality_epoch_validation"
    step4_dir = base_dir / "step4_epoch_disease_associations"

    # Step 2/3 mortality outputs.
    perf2 = read_tsv(step2_dir / "mhas_mortality_epoch_performance.tsv")
    selected = read_tsv(step2_dir / "mhas_mortality_epoch_selected_features.tsv")
    cindex3 = read_tsv(step3_dir / "mhas_step3_cindex_bootstrap_ci.tsv")
    hr3 = read_tsv(step3_dir / "mhas_step3_cox_hr_per_sd.tsv")

    # Step 4 baseline disease-burden outputs. Required and explicitly checked.
    disease_per_sd = read_tsv(step4_dir / "mhas_step4_disease_association_summary.tsv", required=True)
    disease_quartile = read_tsv(step4_dir / "mhas_step4_disease_or_by_quartile.tsv", required=True)
    disease_rate = read_tsv(step4_dir / "mhas_step4_disease_event_rate_by_quartile.tsv", required=True)

    disease_per_sd = require_baseline_prevalent(disease_per_sd, "mhas_step4_disease_association_summary.tsv")
    disease_quartile = require_baseline_prevalent(disease_quartile, "mhas_step4_disease_or_by_quartile.tsv")
    disease_rate = require_baseline_prevalent(disease_rate, "mhas_step4_disease_event_rate_by_quartile.tsv")

    # -------------------------------------------------------------------------
    # Table 1: mortality model performance
    # -------------------------------------------------------------------------
    table1_rows = []
    if not cindex3.empty:
        for split in ["train", "validation", "test", "all"]:
            row = {"split": split}
            split_rows = cindex3[cindex3["split"] == split]
            if len(split_rows):
                row["n"] = split_rows.iloc[0].get("n", np.nan)
                row["deaths"] = split_rows.iloc[0].get("deaths", np.nan)
            for score, label in [
                ("clinical_baseline_lp", "clinical_cindex"),
                ("lp_total", "epoch_lp_cindex"),
                ("mortality_epoch_acceleration_z", "epoch_accel_cindex"),
            ]:
                r = split_rows[split_rows["score"] == score]
                if len(r):
                    rr = r.iloc[0]
                    row[label] = (
                        f"{fmt_num(rr.get('cindex'), 3)} "
                        f"[{fmt_num(rr.get('ci_lower'), 3)}, {fmt_num(rr.get('ci_upper'), 3)}]"
                    )
            if len(row) > 1:
                table1_rows.append(row)
    elif not perf2.empty:
        table1_rows = perf2.to_dict("records")
    table1 = pd.DataFrame(table1_rows)

    # -------------------------------------------------------------------------
    # Table 2: mortality HR per 1-SD EPOCH acceleration
    # -------------------------------------------------------------------------
    table2 = pd.DataFrame()
    if not hr3.empty:
        sub = hr3[hr3["score"] == "mortality_epoch_acceleration_z"].copy()
        sub = keep_valid_status(prefer_adjusted(sub))
        rows = []
        for _, r in sub.iterrows():
            rows.append({
                "split": r.get("split"),
                "n": r.get("n"),
                "deaths": r.get("deaths"),
                "HR_per_SD_EPOCH_acceleration": hr_text(r),
                "status": r.get("status", ""),
            })
        table2 = pd.DataFrame(rows)

    # -------------------------------------------------------------------------
    # Table 3: baseline burden OR per 1-SD
    # -------------------------------------------------------------------------
    sub = disease_per_sd[
        (disease_per_sd["split"] == args.primary_split)
        & (disease_per_sd["score"] == "mortality_epoch_acceleration_z")
    ].copy()
    sub = keep_valid_status(prefer_adjusted(sub))
    sub = sub.sort_values("p", na_position="last").head(args.top_disease_n)

    table3_rows = []
    for _, r in sub.iterrows():
        table3_rows.append({
            "endpoint": r.get("endpoint"),
            "n": r.get("n"),
            "cases": r.get("cases"),
            "baseline_prevalence": fmt_num(r.get("event_rate"), 3),
            "OR_per_SD_EPOCH_acceleration": or_text(r, "or_per_sd"),
            "p": fmt_sci(r.get("p")),
            "FDR_BH": fmt_sci(r.get("fdr_bh")),
            "AUC_score_only": fmt_num(r.get("auc_score_only"), 3),
            "adjustment_covariates": r.get("adjustment_covariates", ""),
        })
    table3 = pd.DataFrame(table3_rows)

    # -------------------------------------------------------------------------
    # Table 4: Q4 vs Q1 baseline burden OR
    # -------------------------------------------------------------------------
    q = disease_quartile[
        (disease_quartile["split"] == args.primary_split)
        & (disease_quartile["comparison"] == "Q4_highest vs Q1_lowest")
    ].copy()
    q = keep_valid_status(prefer_adjusted(q))
    q = q.sort_values("p", na_position="last").head(args.top_disease_n)

    table4_rows = []
    for _, r in q.iterrows():
        table4_rows.append({
            "endpoint": r.get("endpoint"),
            "n": r.get("n"),
            "cases": r.get("cases"),
            "OR_Q4_vs_Q1": or_text(r, "or"),
            "p": fmt_sci(r.get("p")),
            "FDR_BH": fmt_sci(r.get("fdr_bh")),
            "adjustment_covariates": r.get("adjustment_covariates", ""),
        })
    table4 = pd.DataFrame(table4_rows)

    # -------------------------------------------------------------------------
    # Table 5: prevalence by quartile for top endpoint
    # -------------------------------------------------------------------------
    table5 = pd.DataFrame()
    if not sub.empty:
        top_endpoint = sub.iloc[0]["endpoint"]
        table5 = disease_rate[
            (disease_rate["split"] == args.primary_split)
            & (disease_rate["endpoint"] == top_endpoint)
        ].copy()
        if not table5.empty:
            order = ["Q1_lowest", "Q2", "Q3", "Q4_highest"]
            table5["quartile"] = pd.Categorical(table5["quartile"], categories=order, ordered=True)
            table5 = table5.sort_values("quartile")
            table5 = table5[["endpoint", "split", "quartile", "n", "cases", "event_rate"]].copy()
            table5 = table5.rename(columns={"event_rate": "baseline_prevalence"})

    # -------------------------------------------------------------------------
    # Top selected features
    # -------------------------------------------------------------------------
    top_selected = pd.DataFrame()
    if not selected.empty:
        cols = [c for c in ["feature", "coef", "abs_coef"] if c in selected.columns]
        top_selected = selected[cols].copy().head(50)

    # -------------------------------------------------------------------------
    # Save tables
    # -------------------------------------------------------------------------
    table1_out = out_dir / "mhas_step5_table1_mortality_performance.tsv"
    table2_out = out_dir / "mhas_step5_table2_mortality_hr.tsv"
    table3_out = out_dir / "mhas_step5_table3_baseline_burden_or_per_sd.tsv"
    table4_out = out_dir / "mhas_step5_table4_baseline_burden_q4_vs_q1.tsv"
    table5_out = out_dir / "mhas_step5_table5_baseline_burden_prevalence_by_quartile.tsv"
    selected_out = out_dir / "mhas_step5_top_selected_features.tsv"
    report_out = out_dir / "mhas_step5_final_report_baseline_burden.md"
    audit_out = out_dir / "mhas_step5_audit_baseline_burden.txt"

    table1.to_csv(table1_out, sep="\t", index=False)
    table2.to_csv(table2_out, sep="\t", index=False)
    table3.to_csv(table3_out, sep="\t", index=False)
    table4.to_csv(table4_out, sep="\t", index=False)
    table5.to_csv(table5_out, sep="\t", index=False)
    top_selected.to_csv(selected_out, sep="\t", index=False)

    # -------------------------------------------------------------------------
    # Figures
    # -------------------------------------------------------------------------
    fig1 = out_dir / "mhas_step5_fig1_mortality_cindex"
    fig2 = out_dir / "mhas_step5_fig2_mortality_hr"
    fig3 = out_dir / "mhas_step5_fig3_baseline_burden_or_per_sd"
    fig4 = out_dir / "mhas_step5_fig4_baseline_burden_q4_vs_q1"
    fig5 = out_dir / "mhas_step5_fig5_top_baseline_burden_prevalence_by_quartile"

    if not args.skip_plots:
        plot_mortality_cindex(cindex3, str(fig1))
        plot_mortality_hr(hr3, str(fig2))
        plot_baseline_burden_or(
            disease_per_sd,
            str(fig3),
            primary_split=args.primary_split,
            max_endpoints=args.top_disease_n,
        )
        plot_baseline_burden_q4(
            disease_quartile,
            str(fig4),
            primary_split=args.primary_split,
            max_endpoints=args.top_disease_n,
        )
        plot_top_prevalence_by_quartile(
            disease_rate,
            disease_per_sd,
            str(fig5),
            primary_split=args.primary_split,
        )

    # -------------------------------------------------------------------------
    # Markdown report and audit
    # -------------------------------------------------------------------------
    paths = {
        "Table 1 mortality performance": table1_out,
        "Table 2 mortality HR": table2_out,
        "Table 3 baseline burden OR per SD": table3_out,
        "Table 4 baseline burden Q4 vs Q1": table4_out,
        "Table 5 baseline prevalence by quartile": table5_out,
        "Top selected features": selected_out,
        "Figure 1 mortality C-index": str(fig1) + ".pdf",
        "Figure 2 mortality HR": str(fig2) + ".pdf",
        "Figure 3 baseline burden OR per SD": str(fig3) + ".pdf",
        "Figure 4 baseline burden Q4 vs Q1": str(fig4) + ".pdf",
        "Figure 5 top baseline prevalence by quartile": str(fig5) + ".pdf",
    }

    make_markdown_report(
        report_out,
        table1=table1,
        table2=table2,
        table3=table3,
        table4=table4,
        table5=table5,
        selected_features=top_selected,
        paths=paths,
        primary_split=args.primary_split,
        top_n=min(10, args.top_disease_n),
    )

    audit = f"""MHAS STEP 5 CLEAN: baseline disease-burden final report

Base directory
--------------
{base_dir}

Output directory
----------------
{out_dir}

Inputs detected
---------------
Step 2 performance: {not perf2.empty}
Step 2 selected features: {not selected.empty}
Step 3 C-index table: {not cindex3.empty}
Step 3 mortality HR table: {not hr3.empty}
Step 4 baseline disease OR per SD table: {not disease_per_sd.empty}
Step 4 baseline disease quartile OR table: {not disease_quartile.empty}
Step 4 baseline disease prevalence-by-quartile table: {not disease_rate.empty}

Disease analysis included
-------------------------
ONLY baseline_prevalent disease-burden associations.

Disease analysis explicitly excluded
------------------------------------
Incident/future disease associations.
Old downstream disease forest outputs.
Legacy filenames such as mhas_step5_fig_disease_forest_final.pdf.

Primary split
-------------
{args.primary_split}

Final tables
------------
{table1_out}
{table2_out}
{table3_out}
{table4_out}
{table5_out}
{selected_out}

Final report
------------
{report_out}

Figures
-------
{str(fig1) + ".pdf"}
{str(fig2) + ".pdf"}
{str(fig3) + ".pdf"}
{str(fig4) + ".pdf"}
{str(fig5) + ".pdf"}

Notes
-----
This script intentionally fails if the Step 4 disease files do not contain
analysis_type == baseline_prevalent. This prevents accidentally reusing old
incident/future-disease results.
"""
    audit_out.write_text(audit)
    log("\n" + audit)
    log("STEP 5 clean baseline disease-burden report finished successfully.")


if __name__ == "__main__":
    main()
