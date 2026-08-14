#!/usr/bin/env python3
"""
Brain proteomics BAG-EPOCH discordance / resilience analysis.

Primary purpose
---------------
Implement analytic steps 1-3 for the matched brain-proteomics clocks:

1) Build the complete overlapping population between:
     - Brain_ProtBAG
     - brain_proteomics_mortality_clock_acceleration_z

2) Standardize Brain_ProtBAG within the overlapping sample and fit the primary
   residual model:

       EPOCH_z = intercept + beta * BAG_z + error

   The participant-level error term is the EPOCH-BAG discordance residual.
   Negative residual = mortality proximity lower than expected for the person's
   BAG (resilience direction).
   Positive residual = mortality proximity higher than expected for the person's
   BAG (latent-vulnerability direction).

3) Define candidate resilient agers (CRA) and latent vulnerability agers (LVA)
   at multiple symmetric percentile thresholds. By default, 10%, 20%, and 25%:

   CRA at threshold p:
       BAG_z >= BAG (1-p) quantile AND residual_z <= residual p quantile

   LVA at threshold p:
       BAG_z <= BAG p quantile AND residual_z >= residual (1-p) quantile

Examples:
   CRA_p10 = top 10% BAG AND bottom 10% EPOCH|BAG residual
   LVA_p10 = bottom 10% BAG AND top 10% EPOCH|BAG residual

Notes
-----
- Brain_ProtBAG is described as age-bias-corrected but not z-scored, so this
  script z-scores it in the complete overlapping analysis population.
- The EPOCH column is already z-scored and is used directly.
- The primary residual model intentionally does not re-adjust for age because
  both inputs are already age-gap/acceleration-style phenotypes. Age, sex, split,
  survival fields, and other available EPOCH metadata are retained in the output
  for downstream validation and characterization.
- This script does NOT use subsequent outcomes to define the four phenotypes.
  Outcome-based validation should be performed separately.

Dependencies: pandas, numpy
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable, List

import numpy as np
import pandas as pd


DEFAULT_EPOCH_TSV = (
    "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/"
    "Brain_proteomics_mortality_clock/brain_proteomics_mortality_clock_predictions.tsv"
)
DEFAULT_BAG_TSV = "/Users/hao/cubic-home/Reproducibile_paper/SleepAging/data/MomoBAG.tsv"
DEFAULT_OUTPUT_TSV = (
    "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/"
    "Brain_proteomics_mortality_clock/brain_proteomics_BAG_EPOCH_discordance_resilience.tsv"
)
DEFAULT_SUMMARY_TSV = (
    "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/"
    "Brain_proteomics_mortality_clock/brain_proteomics_BAG_EPOCH_discordance_summary.tsv"
)

ID_COL = "participant_id"
BAG_COL = "Brain_ProtBAG"
EPOCH_COL = "brain_proteomics_mortality_clock_acceleration_z"
BAG_Z_COL = "Brain_ProtBAG_z"
EPOCH_PRED_COL = "brain_proteomics_EPOCH_z_predicted_from_BAG"
RESID_COL = "EPOCH_BAG_discordance_residual"
RESID_Z_COL = "EPOCH_BAG_discordance_residual_z"
BAG_PERCENTILE_COL = "Brain_ProtBAG_percentile"
RESID_PERCENTILE_COL = "EPOCH_BAG_residual_percentile"

# Keep these EPOCH columns if present. The script also retains all EPOCH columns
# by default, but this list determines the preferred ordering in the final TSV.
PREFERRED_EPOCH_META = [
    "sample_date",
    "death_date",
    "admin_censor_date",
    "end_date",
    "time_years",
    "event",
    "age_at_baseline",
    "age_at_imaging",
    "sex",
    "organ_source_file",
    "bmi_at_baseline",
    "diastolic_bp_at_baseline",
    "systolic_bp_at_baseline",
    "smoking_status_at_baseline",
    "uk_biobank_assessment_centre_f54_0_0",
    "split",
    "brain_proteomics_mortality_risk_score",
    "risk_score_M0_age_sex",
    "risk_score_M1_covariate_baseline",
    "risk_score_M2_brain_proteomics_only",
    "risk_score_M3_full_covariates_plus_brain_proteomics",
    "risk_5y",
    "risk_10y",
    "risk_15y",
    EPOCH_COL,
    "brain_proteomics_mortality_clock_acceleration_years",
    "brain_proteomics_mortality_clock_age_years",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Residual-based four-phenotype BAG-EPOCH analysis for Brain_ProtBAG vs brain proteomics mortality EPOCH."
    )
    p.add_argument("--epoch_tsv", default=DEFAULT_EPOCH_TSV)
    p.add_argument("--bag_tsv", default=DEFAULT_BAG_TSV)
    p.add_argument("--output_tsv", default=DEFAULT_OUTPUT_TSV)
    p.add_argument("--summary_tsv", default=DEFAULT_SUMMARY_TSV)
    p.add_argument(
        "--thresholds",
        default="0.10,0.20,0.25",
        help="Comma-separated symmetric tail fractions used for the four phenotype definitions. Default: 0.10,0.20,0.25",
    )
    p.add_argument(
        "--drop_extra_epoch_columns",
        action="store_true",
        help=(
            "Drop non-preferred extra columns from the EPOCH predictions TSV. "
            "By default all EPOCH columns are retained."
        ),
    )
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def parse_thresholds(text: str) -> List[float]:
    vals: List[float] = []
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        value = float(item)
        if not 0 < value < 0.5:
            raise ValueError(f"Each threshold must be >0 and <0.5; got {value}")
        vals.append(value)
    vals = sorted(set(vals))
    if not vals:
        raise ValueError("At least one threshold is required.")
    return vals


def ensure_unique_ids(df: pd.DataFrame, source_name: str) -> None:
    if ID_COL not in df.columns:
        raise ValueError(f"{source_name} is missing required column: {ID_COL}")
    dup = df[ID_COL].duplicated(keep=False)
    if dup.any():
        examples = df.loc[dup, ID_COL].astype(str).head(10).tolist()
        raise ValueError(
            f"{source_name} has {int(dup.sum())} rows belonging to duplicated participant_id values. "
            f"Please resolve duplicates before this one-to-one analysis. Examples: {examples}"
        )


def safe_zscore(x: pd.Series, label: str) -> pd.Series:
    vals = pd.to_numeric(x, errors="coerce")
    mean = float(vals.mean())
    sd = float(vals.std(ddof=1))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError(f"{label} has zero or undefined SD in the overlapping sample.")
    return (vals - mean) / sd


def percentile_rank(x: pd.Series) -> pd.Series:
    # Average ranks for ties; output on approximately [0, 1].
    return x.rank(method="average", pct=True)


def threshold_tag(p: float) -> str:
    # 0.10 -> p10; 0.025 -> p2p5
    pct = p * 100
    if abs(pct - round(pct)) < 1e-10:
        return f"p{int(round(pct))}"
    return "p" + str(pct).replace(".", "p")


def main() -> None:
    args = parse_args()
    thresholds = parse_thresholds(args.thresholds)

    out_path = Path(args.output_tsv)
    summary_path = Path(args.summary_tsv)
    for path in [out_path, summary_path]:
        if path.exists() and path.stat().st_size > 0 and not args.overwrite:
            raise FileExistsError(f"Output exists: {path}. Use --overwrite to replace it.")
        path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Reading EPOCH predictions: {args.epoch_tsv}", flush=True)
    epoch = pd.read_csv(args.epoch_tsv, sep="\t", low_memory=False)
    print(f"Reading BAG table: {args.bag_tsv}", flush=True)
    bag = pd.read_csv(args.bag_tsv, sep="\t", low_memory=False)

    for c in [ID_COL, EPOCH_COL]:
        if c not in epoch.columns:
            raise ValueError(f"EPOCH TSV is missing required column: {c}")
    for c in [ID_COL, BAG_COL]:
        if c not in bag.columns:
            raise ValueError(f"BAG TSV is missing required column: {c}")

    epoch = epoch.copy()
    bag = bag[[ID_COL, BAG_COL]].copy()
    epoch[ID_COL] = pd.to_numeric(epoch[ID_COL], errors="coerce").astype("Int64")
    bag[ID_COL] = pd.to_numeric(bag[ID_COL], errors="coerce").astype("Int64")
    epoch = epoch[epoch[ID_COL].notna()].copy()
    bag = bag[bag[ID_COL].notna()].copy()

    ensure_unique_ids(epoch, "EPOCH TSV")
    ensure_unique_ids(bag, "BAG TSV")

    n_epoch = len(epoch)
    n_bag = len(bag)

    # ------------------------------------------------------------------
    # STEP 1. Matched overlapping sample.
    # ------------------------------------------------------------------
    merged = epoch.merge(bag, on=ID_COL, how="inner", validate="one_to_one")
    n_id_overlap = len(merged)

    merged[BAG_COL] = pd.to_numeric(merged[BAG_COL], errors="coerce")
    merged[EPOCH_COL] = pd.to_numeric(merged[EPOCH_COL], errors="coerce")
    merged = merged.replace([np.inf, -np.inf], np.nan)
    merged = merged.dropna(subset=[BAG_COL, EPOCH_COL]).copy()
    n_complete_overlap = len(merged)

    if n_complete_overlap < 20:
        raise ValueError(
            f"Only {n_complete_overlap} participants have non-missing matched BAG and EPOCH values; "
            "too few for a stable residual analysis."
        )

    print(
        f"Matched overlap: {n_id_overlap:,} participant IDs; "
        f"{n_complete_overlap:,} complete BAG+EPOCH pairs.",
        flush=True,
    )

    # ------------------------------------------------------------------
    # STEP 2. Standardize BAG and estimate EPOCH conditional on BAG.
    # ------------------------------------------------------------------
    merged[BAG_Z_COL] = safe_zscore(merged[BAG_COL], BAG_COL)

    # EPOCH is supplied as a z score. Recenter/re-scale only if it materially
    # deviates from z scoring? No: use the supplied phenotype exactly as requested.
    y = merged[EPOCH_COL].to_numpy(dtype=float)
    x = merged[BAG_Z_COL].to_numpy(dtype=float)
    X = np.column_stack([np.ones(len(merged), dtype=float), x])
    coef, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    intercept, beta_bag = float(coef[0]), float(coef[1])
    yhat = X @ coef
    resid = y - yhat

    merged[EPOCH_PRED_COL] = yhat
    merged[RESID_COL] = resid
    merged[RESID_Z_COL] = safe_zscore(pd.Series(resid, index=merged.index), RESID_COL)

    y_mean = float(np.mean(y))
    ss_tot = float(np.sum((y - y_mean) ** 2))
    ss_res = float(np.sum(resid**2))
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    pearson_r = float(np.corrcoef(x, y)[0, 1])
    resid_bag_r = float(np.corrcoef(merged[BAG_Z_COL], merged[RESID_COL])[0, 1])

    merged[BAG_PERCENTILE_COL] = percentile_rank(merged[BAG_Z_COL])
    merged[RESID_PERCENTILE_COL] = percentile_rank(merged[RESID_Z_COL])

    # ------------------------------------------------------------------
    # STEP 3. Continuous discordance + four categorical phenotypes at multiple thresholds.
    # ------------------------------------------------------------------
    summary_rows = []
    for p in thresholds:
        tag = threshold_tag(p)
        bag_low = float(merged[BAG_Z_COL].quantile(p))
        bag_high = float(merged[BAG_Z_COL].quantile(1.0 - p))
        resid_low = float(merged[RESID_Z_COL].quantile(p))
        resid_high = float(merged[RESID_Z_COL].quantile(1.0 - p))

        cfa_col = f"CFA_{tag}"
        cra_col = f"CRA_{tag}"
        lva_col = f"LVA_{tag}"
        cua_col = f"CUA_{tag}"
        group_col = f"aging_phenotype_{tag}"

        # Four mutually exclusive extreme BAG-EPOCH phenotypes.
        # Residual direction is defined relative to the expected EPOCH given BAG.
        merged[cfa_col] = (merged[BAG_Z_COL] <= bag_low) & (merged[RESID_Z_COL] <= resid_low)
        merged[cra_col] = (merged[BAG_Z_COL] >= bag_high) & (merged[RESID_Z_COL] <= resid_low)
        merged[lva_col] = (merged[BAG_Z_COL] <= bag_low) & (merged[RESID_Z_COL] >= resid_high)
        merged[cua_col] = (merged[BAG_Z_COL] >= bag_high) & (merged[RESID_Z_COL] >= resid_high)

        merged[group_col] = np.select(
            [
                merged[cfa_col],
                merged[cra_col],
                merged[lva_col],
                merged[cua_col],
            ],
            [
                "Concordant_favorable_ager",
                "Candidate_resilient_ager",
                "Latent_vulnerability_ager",
                "Concordant_unfavorable_ager",
            ],
            default="Other",
        )

        # Defensive check: because low/high tails do not overlap for p < 0.5,
        # no participant should satisfy more than one extreme phenotype definition.
        phenotype_count = (
            merged[[cfa_col, cra_col, lva_col, cua_col]]
            .astype(int)
            .sum(axis=1)
        )
        if int((phenotype_count > 1).sum()) > 0:
            raise RuntimeError(
                f"Threshold {tag} produced overlapping phenotype assignments, which should not occur."
            )

        summary_rows.append(
            {
                "threshold_tail_fraction": p,
                "threshold_tag": tag,
                "N_complete_overlap": n_complete_overlap,
                "BAG_z_low_cutoff": bag_low,
                "BAG_z_high_cutoff": bag_high,
                "residual_z_low_cutoff": resid_low,
                "residual_z_high_cutoff": resid_high,
                "N_CFA": int(merged[cfa_col].sum()),
                "pct_CFA": float(100.0 * merged[cfa_col].mean()),
                "N_CRA": int(merged[cra_col].sum()),
                "pct_CRA": float(100.0 * merged[cra_col].mean()),
                "N_LVA": int(merged[lva_col].sum()),
                "pct_LVA": float(100.0 * merged[lva_col].mean()),
                "N_CUA": int(merged[cua_col].sum()),
                "pct_CUA": float(100.0 * merged[cua_col].mean()),
                "N_other": int((merged[group_col] == "Other").sum()),
                "pct_other": float(100.0 * (merged[group_col] == "Other").mean()),
            }
        )

    # Global model metadata repeated in summary for provenance.
    summary = pd.DataFrame(summary_rows)
    summary.insert(0, "N_epoch_rows", n_epoch)
    summary.insert(1, "N_bag_rows", n_bag)
    summary.insert(2, "N_ID_overlap", n_id_overlap)
    summary["OLS_intercept"] = intercept
    summary["OLS_beta_BAG_z"] = beta_bag
    summary["OLS_R_squared"] = r_squared
    summary["Pearson_r_BAGz_vs_EPOCHz"] = pearson_r
    summary["Pearson_r_BAGz_vs_residual"] = resid_bag_r
    summary["BAG_column"] = BAG_COL
    summary["EPOCH_column"] = EPOCH_COL
    summary["CFA_definition"] = "low BAG_z AND low EPOCH|BAG residual_z"
    summary["CRA_definition"] = "high BAG_z AND low EPOCH|BAG residual_z"
    summary["LVA_definition"] = "low BAG_z AND high EPOCH|BAG residual_z"
    summary["CUA_definition"] = "high BAG_z AND high EPOCH|BAG residual_z"
    summary["phenotype_abbreviations"] = (
        "CFA=Concordant favorable ager; "
        "CRA=Candidate resilient ager; "
        "LVA=Latent vulnerability ager; "
        "CUA=Concordant unfavorable ager"
    )

    # Add provenance/interpretation columns to participant table.
    merged["discordance_direction"] = np.select(
        [merged[RESID_COL] < 0, merged[RESID_COL] > 0],
        ["lower_EPOCH_than_expected_for_BAG", "higher_EPOCH_than_expected_for_BAG"],
        default="as_expected_for_BAG",
    )
    merged["residual_model_intercept"] = intercept
    merged["residual_model_beta_BAG_z"] = beta_bag
    merged["residual_model_R_squared"] = r_squared

    # Arrange key columns first, followed by retained EPOCH metadata/other columns.
    threshold_cols: List[str] = []
    for p in thresholds:
        tag = threshold_tag(p)
        threshold_cols.extend(
            [
                f"CFA_{tag}",
                f"CRA_{tag}",
                f"LVA_{tag}",
                f"CUA_{tag}",
                f"aging_phenotype_{tag}",
            ]
        )

    key_cols = [
        ID_COL,
        BAG_COL,
        BAG_Z_COL,
        EPOCH_COL,
        EPOCH_PRED_COL,
        RESID_COL,
        RESID_Z_COL,
        BAG_PERCENTILE_COL,
        RESID_PERCENTILE_COL,
        "discordance_direction",
        *threshold_cols,
    ]

    preferred_meta = [c for c in PREFERRED_EPOCH_META if c in merged.columns and c not in key_cols]
    model_meta = [
        "residual_model_intercept",
        "residual_model_beta_BAG_z",
        "residual_model_R_squared",
    ]

    if not args.drop_extra_epoch_columns:
        other_cols = [
            c for c in merged.columns
            if c not in key_cols and c not in preferred_meta and c not in model_meta
        ]
    else:
        other_cols = []

    final_cols = key_cols + preferred_meta + other_cols + model_meta
    # Remove any accidental duplicates while preserving order.
    final_cols = list(dict.fromkeys(final_cols))
    merged = merged.loc[:, final_cols].sort_values(ID_COL).reset_index(drop=True)

    merged.to_csv(out_path, sep="\t", index=False, na_rep="NA")
    summary.to_csv(summary_path, sep="\t", index=False, na_rep="NA")

    print("\nPrimary residual model", flush=True)
    print(f"  N complete overlap: {n_complete_overlap:,}", flush=True)
    print(f"  EPOCH_z = {intercept:.6f} + {beta_bag:.6f} * BAG_z + residual", flush=True)
    print(f"  Pearson r(BAG_z, EPOCH_z): {pearson_r:.6f}", flush=True)
    print(f"  R^2: {r_squared:.6f}", flush=True)
    print(f"  r(BAG_z, residual): {resid_bag_r:.6g}", flush=True)
    print("\nFour-phenotype counts", flush=True)
    print(
        summary[
            [
                "threshold_tag",
                "N_CFA",
                "pct_CFA",
                "N_CRA",
                "pct_CRA",
                "N_LVA",
                "pct_LVA",
                "N_CUA",
                "pct_CUA",
            ]
        ].to_string(index=False),
        flush=True,
    )
    print(f"\nWrote participant-level TSV: {out_path}", flush=True)
    print(f"Wrote summary TSV: {summary_path}", flush=True)


if __name__ == "__main__":
    main()