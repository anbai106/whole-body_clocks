#!/usr/bin/env python3
"""
Collect detailed CRA-vs-CUA brain imaging association results across:
    1) DTI / white-matter features
    2) T1 MUSE gray-matter features
    3) functional MRI connectivity features

This script is adapted for the brain-proteomics EPOCH-BAG resilience analysis.

Expected per-feature input filename
-----------------------------------
    CRA_vs_CUA_logistic_results_<IDP>.tsv

Expected default local-Mac directory structure
----------------------------------------------
/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/mortality_clock/BWAS/
    Brain_proteomics_EPOCH_BAG_resilience/
        CRA_vs_CUA/
            DTI/
            T1/
            FC/

The collector deliberately DISCOVERS result files from each directory instead
of relying on hard-coded IDP lists. This avoids several failure modes in the
older collector:
    - a mismatch between the stated number of DTI features and the actual list,
    - needing to edit the collector every time a modality changes,
    - losing newly added IDPs,
    - failure when one result file is missing.

All columns written by fit_logistic.R are preserved. Additional columns are
added for:
    - modality,
    - imaging feature family,
    - beta direction,
    - nominal significance,
    - BH-FDR and Bonferroni correction within modality,
    - BH-FDR and Bonferroni correction across all three modalities,
    - QC and collection provenance.

Outputs
-------
Per-modality full results:
    BWAS_CRA_vs_CUA_DTI_all_statistics.tsv
    BWAS_CRA_vs_CUA_T1_all_statistics.tsv
    BWAS_CRA_vs_CUA_FC_all_statistics.tsv

Combined:
    BWAS_CRA_vs_CUA_all_modalities_all_statistics.tsv

Subsets:
    BWAS_CRA_vs_CUA_nominal_P_lt_0.05.tsv
    BWAS_CRA_vs_CUA_FDR_lt_0.05_within_modality.tsv
    BWAS_CRA_vs_CUA_FDR_lt_0.05_all_modalities.tsv

QC / summary:
    BWAS_CRA_vs_CUA_collection_audit.tsv
    BWAS_CRA_vs_CUA_modality_summary.tsv

Interpretation
--------------
The regression outcome is:
    CRA = 1
    CUA = 0

Therefore:
    Beta_log_odds_per_1SD_IDP > 0
        higher IDP is associated with greater odds of CRA

    Beta_log_odds_per_1SD_IDP < 0
        higher IDP is associated with greater odds of CUA

The sign is reported separately from statistical significance. A nonsignificant
negative beta is NOT automatically labeled as evidence for CUA.

Usage
-----
Default local Mac:
    python collect_CRA_CUA_BWAS_results.py

Explicit root:
    python collect_CRA_CUA_BWAS_results.py \
        --analysis-root /path/to/Brain_proteomics_EPOCH_BAG_resilience/CRA_vs_CUA

Explicit output:
    python collect_CRA_CUA_BWAS_results.py \
        --output-dir /path/to/summary_results

Cluster example:
    python collect_CRA_CUA_BWAS_results.py \
        --analysis-root \
        /cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/BWAS/Brain_proteomics_EPOCH_BAG_resilience/CRA_vs_CUA
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd


# =============================================================================
# DEFAULT PATHS
# =============================================================================

DEFAULT_ANALYSIS_ROOT = Path(
    "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/"
    "mortality_clock/BWAS/"
    "Brain_proteomics_EPOCH_BAG_resilience/CRA_vs_CUA"
)

DEFAULT_OUTPUT_DIR = (
    DEFAULT_ANALYSIS_ROOT / "combined_results"
)

RESULT_PREFIX = "CRA_vs_CUA_logistic_results_"
RESULT_SUFFIX = ".tsv"

MODALITY_DIRS: Dict[str, str] = {
    "DTI": "DTI",
    "T1": "T1",
    "FC": "FC",
}

MODALITY_LABELS: Dict[str, str] = {
    "DTI": "Brain-WM-DTI",
    "T1": "Brain-GM-T1",
    "FC": "Brain-FC-fMRI",
}

# Optional expected numbers for a simple completeness audit.
#
# T1 and FC are fixed from the submit scripts used here.
# DTI is left as None deliberately because the current project has used both
# an FA-only 48-feature analysis and broader DTI sets containing additional
# MD/ICVF/OD features. File discovery remains the source of truth.
DEFAULT_EXPECTED_COUNTS: Dict[str, Optional[int]] = {
    "DTI": None,
    "T1": 119,
    "FC": 210,
}


# =============================================================================
# REQUIRED / PREFERRED RESULT COLUMNS
# =============================================================================

REQUIRED_COLUMNS = [
    "IDP",
    "Beta_log_odds_per_1SD_IDP",
    "SE",
    "Z",
    "P_value",
    "OR_per_1SD_IDP",
    "OR_CI95_lower",
    "OR_CI95_upper",
    "N_case",
    "N_control",
    "N_total",
]

PREFERRED_COLUMN_ORDER = [
    # Collection / modality metadata
    "Modality",
    "Organ",
    "Feature_family",
    "IDP",
    "source_file",
    "collection_status",

    # Comparison definition
    "comparison",
    "phenotype_column",
    "phenotype_creation_method",
    "case_indicator_column",
    "control_indicator_column",
    "case_label",
    "control_label",
    "case_coding",
    "control_coding",

    # Main association statistics
    "Beta_log_odds_per_1SD_IDP",
    "SE",
    "Z",
    "P_value",
    "OR_per_1SD_IDP",
    "OR_CI95_lower",
    "OR_CI95_upper",

    # Direction and multiple-testing statistics
    "Beta_direction",
    "Nominal_P_lt_0.05",
    "P_value_FDR_BH_within_modality",
    "P_value_Bonferroni_within_modality",
    "FDR_BH_lt_0.05_within_modality",
    "Bonferroni_lt_0.05_within_modality",
    "P_value_FDR_BH_all_modalities",
    "P_value_Bonferroni_all_modalities",
    "FDR_BH_lt_0.05_all_modalities",
    "Bonferroni_lt_0.05_all_modalities",
    "N_tests_within_modality",
    "N_tests_all_modalities",

    # Sample sizes and QC
    "N_case",
    "N_control",
    "N_total",
    "N_case_in_phenotype_file",
    "N_control_in_phenotype_file",
    "N_after_merge",
    "N_complete_before_outlier_filter",
    "N_IDP_outliers_removed",
    "IDP_outlier_SD_threshold",
    "IDP_mean_before_outlier_filter",
    "IDP_SD_before_outlier_filter",
    "IDP_mean_analysis_sample",
    "IDP_SD_analysis_sample",
    "model_converged",
    "possible_separation",
    "AIC",

    # Model specification
    "phenotype_covariates",
    "general_covariates",
    "idp_covariates",
    "factor_covariates",
    "model_formula",
    "glm_warnings",

    # Original R-script direction field retained for provenance
    "Direction",
]


# =============================================================================
# ARGUMENTS
# =============================================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Collect detailed CRA-vs-CUA DTI, T1, and fMRI association results."
        )
    )

    parser.add_argument(
        "--analysis-root",
        type=Path,
        default=DEFAULT_ANALYSIS_ROOT,
        help=(
            "Root directory containing DTI/, T1/, and FC/. "
            f"Default: {DEFAULT_ANALYSIS_ROOT}"
        ),
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=(
            "Output directory. Default: <analysis-root>/combined_results"
        ),
    )

    parser.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help="Significance threshold. Default: 0.05",
    )

    parser.add_argument(
        "--expected-dti",
        type=int,
        default=None,
        help=(
            "Optional expected DTI feature count for audit purposes only. "
            "Default: unspecified."
        ),
    )

    parser.add_argument(
        "--expected-t1",
        type=int,
        default=119,
        help="Expected T1 feature count for audit. Default: 119",
    )

    parser.add_argument(
        "--expected-fc",
        type=int,
        default=210,
        help="Expected FC feature count for audit. Default: 210",
    )

    parser.add_argument(
        "--strict",
        action="store_true",
        help=(
            "Stop on malformed result files. Without --strict, malformed files "
            "are recorded in the audit and collection continues."
        ),
    )

    return parser.parse_args()


# =============================================================================
# MULTIPLE-TESTING HELPERS
# =============================================================================

def bh_fdr(p_values: Iterable[float]) -> np.ndarray:
    """
    Benjamini-Hochberg FDR adjustment.

    Returns NaN for missing/non-finite p-values.
    """
    p = np.asarray(list(p_values), dtype=float)
    adjusted = np.full(p.shape, np.nan, dtype=float)

    valid = np.isfinite(p) & (p >= 0) & (p <= 1)
    p_valid = p[valid]

    if p_valid.size == 0:
        return adjusted

    order = np.argsort(p_valid)
    ranked = p_valid[order]
    m = len(ranked)

    raw_adj = ranked * m / np.arange(1, m + 1)

    # Enforce monotonicity from largest rank backward.
    monotone = np.minimum.accumulate(raw_adj[::-1])[::-1]
    monotone = np.minimum(monotone, 1.0)

    restored = np.empty_like(monotone)
    restored[order] = monotone

    adjusted[valid] = restored
    return adjusted


def bonferroni(p_values: Iterable[float]) -> np.ndarray:
    """
    Bonferroni adjustment across non-missing valid p-values.
    """
    p = np.asarray(list(p_values), dtype=float)
    adjusted = np.full(p.shape, np.nan, dtype=float)

    valid = np.isfinite(p) & (p >= 0) & (p <= 1)
    m = int(valid.sum())

    if m == 0:
        return adjusted

    adjusted[valid] = np.minimum(
        p[valid] * m,
        1.0,
    )

    return adjusted


# =============================================================================
# FEATURE ANNOTATION
# =============================================================================

def infer_feature_family(modality: str, idp: str) -> str:
    """
    Add a compact feature-family annotation without changing the original IDP.
    """
    x = str(idp).lower()

    if modality == "DTI":
        if x.startswith("mean_fa_"):
            return "FA"
        if x.startswith("mean_md_"):
            return "MD"
        if x.startswith("mean_icvf_"):
            return "ICVF"
        if x.startswith("mean_od_"):
            return "OD"
        if x.startswith("mean_isovf_"):
            return "ISOVF"
        return "DTI_other"

    if modality == "T1":
        if str(idp).startswith("MUSE_Volume_"):
            return "MUSE_GM_volume"
        return "T1_other"

    if modality == "FC":
        if re.fullmatch(r"f_\d+", str(idp)):
            return "functional_connectivity"
        return "fMRI_other"

    return "other"


def beta_direction(beta: float) -> str:
    """
    Point-estimate direction only. Significance is handled separately.
    """
    if not np.isfinite(beta):
        return "NA"

    if beta > 0:
        return "Higher_IDP_more_CRA_like"

    if beta < 0:
        return "Higher_IDP_more_CUA_like"

    return "Zero_beta"


# =============================================================================
# FILE DISCOVERY / VALIDATION
# =============================================================================

def discover_result_files(directory: Path) -> List[Path]:
    if not directory.exists():
        return []

    return sorted(
        directory.glob(
            f"{RESULT_PREFIX}*{RESULT_SUFFIX}"
        )
    )


def expected_idp_from_filename(path: Path) -> str:
    name = path.name

    if not (
        name.startswith(RESULT_PREFIX)
        and name.endswith(RESULT_SUFFIX)
    ):
        return ""

    return name[
        len(RESULT_PREFIX) :
        -len(RESULT_SUFFIX)
    ]


def read_one_result(
    path: Path,
    modality: str,
    strict: bool,
) -> Tuple[Optional[pd.DataFrame], dict]:
    """
    Read one per-IDP result TSV.

    Returns:
        dataframe or None
        audit dictionary
    """
    audit = {
        "Modality": modality,
        "source_file": str(path),
        "filename": path.name,
        "expected_IDP_from_filename": expected_idp_from_filename(path),
        "status": "unknown",
        "message": "",
        "N_rows_in_file": np.nan,
        "IDP_in_file": "",
    }

    try:
        df = pd.read_csv(
            path,
            sep="\t",
            low_memory=False,
        )
    except Exception as exc:
        audit["status"] = "read_error"
        audit["message"] = str(exc)

        if strict:
            raise

        return None, audit

    audit["N_rows_in_file"] = len(df)

    missing = [
        c for c in REQUIRED_COLUMNS
        if c not in df.columns
    ]

    if missing:
        audit["status"] = "missing_required_columns"
        audit["message"] = ", ".join(missing)

        if strict:
            raise ValueError(
                f"{path} is missing required columns: "
                + ", ".join(missing)
            )

        return None, audit

    if len(df) == 0:
        audit["status"] = "empty_file"
        audit["message"] = "No result rows."

        if strict:
            raise ValueError(f"Empty result file: {path}")

        return None, audit

    if len(df) != 1:
        audit["status"] = "unexpected_row_count"
        audit["message"] = (
            f"Expected one row per IDP but found {len(df)}."
        )

        if strict:
            raise ValueError(
                f"Expected one row in {path}, found {len(df)}."
            )

        # Keep all rows in non-strict mode but flag the issue.

    df["IDP"] = df["IDP"].astype(str)

    audit["IDP_in_file"] = ";".join(
        sorted(df["IDP"].dropna().unique())
    )

    expected_idp = audit["expected_IDP_from_filename"]

    mismatch = (
        expected_idp
        and not all(df["IDP"] == expected_idp)
    )

    if mismatch:
        audit["status"] = "idp_filename_mismatch"
        audit["message"] = (
            f"Filename implies IDP={expected_idp}; "
            f"file contains {audit['IDP_in_file']}."
        )

        if strict:
            raise ValueError(
                f"IDP mismatch in {path}: {audit['message']}"
            )
    else:
        audit["status"] = (
            "ok"
            if len(df) == 1
            else "unexpected_row_count"
        )

    # Add collection metadata while preserving every original R output column.
    df.insert(
        0,
        "Modality",
        modality,
    )

    df.insert(
        1,
        "Organ",
        MODALITY_LABELS[modality],
    )

    df.insert(
        2,
        "Feature_family",
        [
            infer_feature_family(
                modality,
                x,
            )
            for x in df["IDP"]
        ],
    )

    df.insert(
        3,
        "source_file",
        str(path),
    )

    df.insert(
        4,
        "collection_status",
        audit["status"],
    )

    return df, audit


# =============================================================================
# RESULT ANNOTATION
# =============================================================================

def coerce_result_types(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()

    numeric_cols = [
        "Beta_log_odds_per_1SD_IDP",
        "SE",
        "Z",
        "P_value",
        "OR_per_1SD_IDP",
        "OR_CI95_lower",
        "OR_CI95_upper",
        "N_case",
        "N_control",
        "N_total",
        "N_case_in_phenotype_file",
        "N_control_in_phenotype_file",
        "N_after_merge",
        "N_complete_before_outlier_filter",
        "N_IDP_outliers_removed",
        "IDP_outlier_SD_threshold",
        "IDP_mean_before_outlier_filter",
        "IDP_SD_before_outlier_filter",
        "IDP_mean_analysis_sample",
        "IDP_SD_analysis_sample",
        "AIC",
    ]

    for col in numeric_cols:
        if col in out.columns:
            out[col] = pd.to_numeric(
                out[col],
                errors="coerce",
            )

    return out


def add_within_modality_multiple_testing(
    df: pd.DataFrame,
    alpha: float,
) -> pd.DataFrame:
    out = df.copy()

    out["P_value_FDR_BH_within_modality"] = np.nan
    out["P_value_Bonferroni_within_modality"] = np.nan
    out["N_tests_within_modality"] = np.nan

    for modality, idx in out.groupby(
        "Modality",
        sort=False,
    ).groups.items():

        idx = list(idx)

        p = pd.to_numeric(
            out.loc[idx, "P_value"],
            errors="coerce",
        ).to_numpy()

        valid_n = int(
            (
                np.isfinite(p)
                & (p >= 0)
                & (p <= 1)
            ).sum()
        )

        out.loc[
            idx,
            "P_value_FDR_BH_within_modality",
        ] = bh_fdr(p)

        out.loc[
            idx,
            "P_value_Bonferroni_within_modality",
        ] = bonferroni(p)

        out.loc[
            idx,
            "N_tests_within_modality",
        ] = valid_n

    out["FDR_BH_lt_0.05_within_modality"] = (
        out["P_value_FDR_BH_within_modality"]
        < alpha
    )

    out["Bonferroni_lt_0.05_within_modality"] = (
        out["P_value_Bonferroni_within_modality"]
        < alpha
    )

    return out


def add_global_multiple_testing(
    df: pd.DataFrame,
    alpha: float,
) -> pd.DataFrame:
    out = df.copy()

    p = pd.to_numeric(
        out["P_value"],
        errors="coerce",
    ).to_numpy()

    valid_n = int(
        (
            np.isfinite(p)
            & (p >= 0)
            & (p <= 1)
        ).sum()
    )

    out["P_value_FDR_BH_all_modalities"] = bh_fdr(p)
    out["P_value_Bonferroni_all_modalities"] = bonferroni(p)

    out["FDR_BH_lt_0.05_all_modalities"] = (
        out["P_value_FDR_BH_all_modalities"]
        < alpha
    )

    out["Bonferroni_lt_0.05_all_modalities"] = (
        out["P_value_Bonferroni_all_modalities"]
        < alpha
    )

    out["N_tests_all_modalities"] = valid_n

    return out


def annotate_results(
    df: pd.DataFrame,
    alpha: float,
) -> pd.DataFrame:
    out = coerce_result_types(df)

    out["Beta_direction"] = [
        beta_direction(x)
        for x in out[
            "Beta_log_odds_per_1SD_IDP"
        ].to_numpy(dtype=float)
    ]

    out["Nominal_P_lt_0.05"] = (
        out["P_value"] < alpha
    )

    out = add_within_modality_multiple_testing(
        out,
        alpha=alpha,
    )

    out = add_global_multiple_testing(
        out,
        alpha=alpha,
    )

    return out


def reorder_columns(df: pd.DataFrame) -> pd.DataFrame:
    preferred = [
        c for c in PREFERRED_COLUMN_ORDER
        if c in df.columns
    ]

    remaining = [
        c for c in df.columns
        if c not in preferred
    ]

    return df[
        preferred + remaining
    ]


# =============================================================================
# SUMMARY TABLES
# =============================================================================

def make_modality_summary(
    df: pd.DataFrame,
    audit_df: pd.DataFrame,
    expected_counts: Dict[str, Optional[int]],
    alpha: float,
) -> pd.DataFrame:

    rows = []

    for modality in MODALITY_DIRS:

        sub = df[
            df["Modality"] == modality
        ].copy()

        audit_sub = audit_df[
            audit_df["Modality"] == modality
        ].copy()

        expected_n = expected_counts.get(
            modality
        )

        n_files_discovered = len(
            audit_sub
        )

        n_valid_files = int(
            audit_sub["status"].isin(
                [
                    "ok",
                    "unexpected_row_count",
                    "idp_filename_mismatch",
                ]
            ).sum()
        ) if len(audit_sub) else 0

        n_results = len(sub)

        if n_results:
            p_numeric = pd.to_numeric(
                sub["P_value"],
                errors="coerce",
            )

            best_idx = p_numeric.idxmin()

            best_idp = sub.loc[
                best_idx,
                "IDP",
            ]

            best_p = sub.loc[
                best_idx,
                "P_value",
            ]

            best_beta = sub.loc[
                best_idx,
                "Beta_log_odds_per_1SD_IDP",
            ]

            best_or = sub.loc[
                best_idx,
                "OR_per_1SD_IDP",
            ]

            n_nominal = int(
                sub["Nominal_P_lt_0.05"].sum()
            )

            n_fdr_modality = int(
                sub[
                    "FDR_BH_lt_0.05_within_modality"
                ].sum()
            )

            n_fdr_global = int(
                sub[
                    "FDR_BH_lt_0.05_all_modalities"
                ].sum()
            )

            n_positive_beta = int(
                (
                    sub[
                        "Beta_log_odds_per_1SD_IDP"
                    ] > 0
                ).sum()
            )

            n_negative_beta = int(
                (
                    sub[
                        "Beta_log_odds_per_1SD_IDP"
                    ] < 0
                ).sum()
            )

        else:
            best_idp = ""
            best_p = np.nan
            best_beta = np.nan
            best_or = np.nan
            n_nominal = 0
            n_fdr_modality = 0
            n_fdr_global = 0
            n_positive_beta = 0
            n_negative_beta = 0

        rows.append(
            {
                "Modality": modality,
                "Organ": MODALITY_LABELS[modality],
                "Expected_N_features": expected_n,
                "N_result_files_discovered": n_files_discovered,
                "N_valid_result_files": n_valid_files,
                "N_result_rows_collected": n_results,
                "Expected_count_difference": (
                    np.nan
                    if expected_n is None
                    else n_results - expected_n
                ),
                "N_nominal_P_lt_alpha": n_nominal,
                "N_FDR_lt_alpha_within_modality": n_fdr_modality,
                "N_FDR_lt_alpha_all_modalities": n_fdr_global,
                "N_positive_beta_CRA_direction": n_positive_beta,
                "N_negative_beta_CUA_direction": n_negative_beta,
                "Best_IDP": best_idp,
                "Best_P_value": best_p,
                "Best_Beta": best_beta,
                "Best_OR": best_or,
                "alpha": alpha,
            }
        )

    return pd.DataFrame(
        rows
    )


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    args = parse_args()

    if not (
        0 < args.alpha < 1
    ):
        raise ValueError(
            "--alpha must be between 0 and 1."
        )

    analysis_root = args.analysis_root.expanduser()

    output_dir = (
        args.output_dir.expanduser()
        if args.output_dir is not None
        else analysis_root / "combined_results"
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    expected_counts: Dict[str, Optional[int]] = {
        "DTI": args.expected_dti,
        "T1": args.expected_t1,
        "FC": args.expected_fc,
    }

    print(
        "============================================================"
    )
    print(
        "Collect CRA-vs-CUA brain imaging association results"
    )
    print(
        "============================================================"
    )
    print(
        f"Analysis root: {analysis_root}"
    )
    print(
        f"Output dir:    {output_dir}"
    )
    print()

    collected: List[pd.DataFrame] = []
    audit_rows: List[dict] = []

    # -------------------------------------------------------------------------
    # Discover and read all modality results.
    # -------------------------------------------------------------------------
    for modality, dirname in MODALITY_DIRS.items():

        modality_dir = (
            analysis_root / dirname
        )

        files = discover_result_files(
            modality_dir
        )

        print(
            f"{modality}: {len(files)} result file(s) discovered in"
        )
        print(
            f"  {modality_dir}"
        )

        if len(files) == 0:
            audit_rows.append(
                {
                    "Modality": modality,
                    "source_file": "",
                    "filename": "",
                    "expected_IDP_from_filename": "",
                    "status": "no_result_files_found",
                    "message": (
                        f"No files matching "
                        f"{RESULT_PREFIX}*{RESULT_SUFFIX}"
                    ),
                    "N_rows_in_file": np.nan,
                    "IDP_in_file": "",
                }
            )
            continue

        for path in files:

            df_one, audit = read_one_result(
                path=path,
                modality=modality,
                strict=args.strict,
            )

            audit_rows.append(
                audit
            )

            if df_one is not None:
                collected.append(
                    df_one
                )

    audit_df = pd.DataFrame(
        audit_rows
    )

    audit_path = (
        output_dir
        / "BWAS_CRA_vs_CUA_collection_audit.tsv"
    )

    audit_df.to_csv(
        audit_path,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    if len(collected) == 0:
        print()
        print(
            "No valid result rows were collected."
        )
        print(
            f"Audit written to: {audit_path}"
        )
        sys.exit(1)

    # -------------------------------------------------------------------------
    # Combine and annotate.
    # -------------------------------------------------------------------------
    combined = pd.concat(
        collected,
        ignore_index=True,
        sort=False,
    )

    combined = annotate_results(
        combined,
        alpha=args.alpha,
    )

    combined = reorder_columns(
        combined
    )

    # Sort by modality then P-value while preserving all statistics.
    modality_order = {
        "DTI": 0,
        "T1": 1,
        "FC": 2,
    }

    combined["_modality_order"] = combined[
        "Modality"
    ].map(
        modality_order
    )

    combined = combined.sort_values(
        by=[
            "_modality_order",
            "P_value",
            "IDP",
        ],
        ascending=[
            True,
            True,
            True,
        ],
        kind="mergesort",
        na_position="last",
    ).drop(
        columns="_modality_order"
    ).reset_index(
        drop=True
    )

    # -------------------------------------------------------------------------
    # Per-modality outputs.
    # -------------------------------------------------------------------------
    for modality in MODALITY_DIRS:

        sub = combined[
            combined["Modality"] == modality
        ].copy()

        out_path = (
            output_dir
            / f"BWAS_CRA_vs_CUA_{modality}_all_statistics.tsv"
        )

        sub.to_csv(
            out_path,
            sep="\t",
            index=False,
            na_rep="NA",
        )

        print(
            f"Wrote {modality}: {len(sub)} row(s)"
        )
        print(
            f"  {out_path}"
        )

    # -------------------------------------------------------------------------
    # Combined detailed output.
    # -------------------------------------------------------------------------
    combined_path = (
        output_dir
        / "BWAS_CRA_vs_CUA_all_modalities_all_statistics.tsv"
    )

    combined.to_csv(
        combined_path,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    # -------------------------------------------------------------------------
    # Significant / nominal subsets.
    # -------------------------------------------------------------------------
    nominal = combined[
        combined[
            "Nominal_P_lt_0.05"
        ].fillna(False)
    ].copy()

    fdr_modality = combined[
        combined[
            "FDR_BH_lt_0.05_within_modality"
        ].fillna(False)
    ].copy()

    fdr_global = combined[
        combined[
            "FDR_BH_lt_0.05_all_modalities"
        ].fillna(False)
    ].copy()

    nominal.to_csv(
        output_dir
        / "BWAS_CRA_vs_CUA_nominal_P_lt_0.05.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )

    fdr_modality.to_csv(
        output_dir
        / "BWAS_CRA_vs_CUA_FDR_lt_0.05_within_modality.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )

    fdr_global.to_csv(
        output_dir
        / "BWAS_CRA_vs_CUA_FDR_lt_0.05_all_modalities.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )

    # -------------------------------------------------------------------------
    # Modality summary.
    # -------------------------------------------------------------------------
    modality_summary = make_modality_summary(
        df=combined,
        audit_df=audit_df,
        expected_counts=expected_counts,
        alpha=args.alpha,
    )

    summary_path = (
        output_dir
        / "BWAS_CRA_vs_CUA_modality_summary.tsv"
    )

    modality_summary.to_csv(
        summary_path,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    # -------------------------------------------------------------------------
    # Console summary.
    # -------------------------------------------------------------------------
    print()
    print(
        "============================================================"
    )
    print(
        "Collection complete"
    )
    print(
        "============================================================"
    )

    print(
        f"Total association rows: {len(combined)}"
    )

    print(
        f"Nominal P < {args.alpha}: "
        f"{int(combined['Nominal_P_lt_0.05'].sum())}"
    )

    print(
        f"BH-FDR < {args.alpha} within modality: "
        f"{int(combined['FDR_BH_lt_0.05_within_modality'].sum())}"
    )

    print(
        f"BH-FDR < {args.alpha} across all modalities: "
        f"{int(combined['FDR_BH_lt_0.05_all_modalities'].sum())}"
    )

    print()

    with pd.option_context(
        "display.max_columns",
        None,
        "display.width",
        220,
    ):
        print(
            modality_summary.to_string(
                index=False
            )
        )

    print()
    print(
        "Main combined detailed file:"
    )
    print(
        f"  {combined_path}"
    )

    print(
        "Collection audit:"
    )
    print(
        f"  {audit_path}"
    )

    print(
        "Modality summary:"
    )
    print(
        f"  {summary_path}"
    )


if __name__ == "__main__":
    main()