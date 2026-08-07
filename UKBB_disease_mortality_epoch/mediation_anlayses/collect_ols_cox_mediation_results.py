#!/usr/bin/env python3
"""
Collect EPOCH OLS + Cox + bootstrap mediation results
=====================================================

This script is intended to be run locally from a Mac after copying the completed
SLURM result directory to your computer.

It does NOT use argparse. Edit only the USER SETTINGS section below.

Expected structure
------------------
JOB_DIR/
    single_model_results/
        model_000__....tsv
        model_001__....tsv
        ...
        model_314__....tsv

Pipeline
--------
1. Read all one-row model_*.tsv files from single_model_results/.
2. Verify that model indices 0-314 are represented once each.
3. Concatenate all model results and sort by model_index.
4. Recalculate multiple-testing statistics across the full 315-model family:
       Bonferroni threshold = 0.05 / 315
       Bonferroni adjusted p = min(indirect_delta_p * 315, 1)
       BH-FDR q-values across finite indirect_delta_p values
5. Recalculate useful QC flags:
       bootstrap CI excludes zero
       EPP < 5
       EPP < 10
       VIF > 5
       PH-test p < 0.05
       bootstrap success fraction < 0.80
       warnings present
       status != ok
6. Write combined result tables and summary files.

The script does NOT refit any mediation model.
"""

from pathlib import Path
import json
import re

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


# =============================================================================
# USER SETTINGS -- EDIT THESE
# =============================================================================

# Set this to the completed SLURM array job directory on your Mac.
#
# Example:
# JOB_DIR = Path(
#     "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/"
#     "mediation_OLS_Cox_bootstrap_full_single_models/job_16997583"
# )
#
# If this Python script is placed directly inside job_16997583, you can use:
# JOB_DIR = Path(__file__).resolve().parent

JOB_DIR = Path(
    "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/"
    "mediation_OLS_Cox_bootstrap_full_single_models/job_16997583"
)

EXPECTED_MODELS = 315
ALPHA = 0.05

# Results will be written here.
OUTPUT_DIR = JOB_DIR / "collected_results"

# If False, the script will stop when any of the 315 model indices are missing.
# Set True only if you intentionally want to collect an incomplete run.
ALLOW_INCOMPLETE = False


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def model_index_from_filename(path: Path):
    match = re.match(r"model_(\d+)__", path.name)
    if match is None:
        return None
    return int(match.group(1))


def read_one_model(path: Path) -> pd.DataFrame:
    """Read and validate one single-model result TSV."""
    df = pd.read_csv(path, sep="\t", low_memory=False)

    if len(df) != 1:
        raise ValueError(
            f"{path.name}: expected exactly one row, found {len(df)}"
        )

    filename_index = model_index_from_filename(path)
    if filename_index is None:
        raise ValueError(
            f"Could not obtain model index from filename: {path.name}"
        )

    if "model_index" not in df.columns:
        df["model_index"] = filename_index
    else:
        value = pd.to_numeric(
            df["model_index"], errors="coerce"
        ).iloc[0]

        if not np.isfinite(value):
            df.loc[df.index[0], "model_index"] = filename_index
        elif int(value) != filename_index:
            raise ValueError(
                f"{path.name}: model_index in file ({int(value)}) "
                f"does not match filename ({filename_index})"
            )

    df["source_file"] = path.name
    return df


def numeric_columns(df: pd.DataFrame, columns):
    """Convert selected columns to numeric when present."""
    for column in columns:
        if column in df.columns:
            df[column] = pd.to_numeric(
                df[column], errors="coerce"
            )


def add_global_statistics(df: pd.DataFrame) -> pd.DataFrame:
    """Recalculate global multiple-testing and QC statistics."""
    out = df.copy()

    numeric_columns(
        out,
        [
            "model_index",
            "model_number",
            "n",
            "deaths",
            "censored",
            "events_per_parameter_direct",
            "n_covariates",
            "a_beta",
            "a_se_hc3",
            "a_p_hc3",
            "mediator_model_r2",
            "b_log_hr",
            "b_hr",
            "b_se",
            "b_p",
            "direct_log_hr",
            "direct_hr",
            "direct_se",
            "direct_p",
            "total_log_hr",
            "total_hr",
            "total_se",
            "total_p",
            "indirect_log_hr",
            "indirect_hr_like",
            "indirect_delta_se",
            "indirect_delta_z",
            "indirect_delta_p",
            "boot_ci_low",
            "boot_ci_high",
            "boot_hr_like_ci_low",
            "boot_hr_like_ci_high",
            "boot_p_empirical",
            "bootstrap_requested",
            "bootstrap_successes",
            "bootstrap_success_fraction",
            "max_vif",
            "ph_exposure_p",
            "ph_mediator_p",
            "cox_direct_concordance",
            "cox_total_concordance",
        ],
    )

    # -------------------------------------------------------------------------
    # Bonferroni correction across the FULL prespecified family of 315 tests.
    # -------------------------------------------------------------------------
    bonf_threshold = ALPHA / EXPECTED_MODELS

    out["multiplicity_family_m"] = EXPECTED_MODELS
    out["bonferroni_alpha_threshold"] = bonf_threshold

    p = pd.to_numeric(
        out.get(
            "indirect_delta_p",
            pd.Series(np.nan, index=out.index),
        ),
        errors="coerce",
    )

    finite_p = np.isfinite(p)

    out["indirect_bonferroni_p"] = np.nan
    out.loc[
        finite_p, "indirect_bonferroni_p"
    ] = np.minimum(
        p.loc[finite_p] * EXPECTED_MODELS,
        1.0,
    )

    out["indirect_bonferroni_significant"] = (
        finite_p & (p < bonf_threshold)
    )

    # -------------------------------------------------------------------------
    # BH-FDR across all estimable indirect-effect p-values.
    # -------------------------------------------------------------------------
    out["indirect_bh_fdr"] = np.nan
    out["indirect_bh_fdr_significant"] = False

    if finite_p.any():
        _, qvalues, _, _ = multipletests(
            p.loc[finite_p].to_numpy(dtype=float),
            alpha=ALPHA,
            method="fdr_bh",
        )

        out.loc[
            finite_p, "indirect_bh_fdr"
        ] = qvalues

        out.loc[
            finite_p, "indirect_bh_fdr_significant"
        ] = qvalues < ALPHA

    # -------------------------------------------------------------------------
    # Bootstrap support.
    # -------------------------------------------------------------------------
    if {
        "boot_ci_low",
        "boot_ci_high",
    }.issubset(out.columns):

        out["bootstrap_ci_excludes_zero"] = (
            out["boot_ci_low"].notna()
            & out["boot_ci_high"].notna()
            & (
                (out["boot_ci_low"] > 0)
                | (out["boot_ci_high"] < 0)
            )
        )
    else:
        out["bootstrap_ci_excludes_zero"] = False

    # -------------------------------------------------------------------------
    # QC flags.
    # -------------------------------------------------------------------------
    if "events_per_parameter_direct" in out.columns:
        out["epp_lt_5"] = (
            out["events_per_parameter_direct"].notna()
            & (
                out["events_per_parameter_direct"]
                < 5
            )
        )

        out["epp_lt_10"] = (
            out["events_per_parameter_direct"].notna()
            & (
                out["events_per_parameter_direct"]
                < 10
            )
        )
    else:
        out["epp_lt_5"] = False
        out["epp_lt_10"] = False

    if "max_vif" in out.columns:
        out["vif_gt_5"] = (
            out["max_vif"].notna()
            & (out["max_vif"] > 5)
        )
    else:
        out["vif_gt_5"] = False

    if "ph_exposure_p" in out.columns:
        out[
            "ph_exposure_violation_p_lt_0_05"
        ] = (
            out["ph_exposure_p"].notna()
            & (out["ph_exposure_p"] < 0.05)
        )
    else:
        out[
            "ph_exposure_violation_p_lt_0_05"
        ] = False

    if "ph_mediator_p" in out.columns:
        out[
            "ph_mediator_violation_p_lt_0_05"
        ] = (
            out["ph_mediator_p"].notna()
            & (out["ph_mediator_p"] < 0.05)
        )
    else:
        out[
            "ph_mediator_violation_p_lt_0_05"
        ] = False

    if "bootstrap_success_fraction" in out.columns:
        out["bootstrap_success_lt_0_80"] = (
            out["bootstrap_success_fraction"].notna()
            & (
                out["bootstrap_success_fraction"]
                < 0.80
            )
        )
    else:
        out["bootstrap_success_lt_0_80"] = False

    if "warnings" in out.columns:
        out["has_warnings"] = (
            out["warnings"]
            .fillna("")
            .astype(str)
            .str.strip()
            .ne("")
        )
    else:
        out["has_warnings"] = False

    if "status" in out.columns:
        out["model_not_ok"] = (
            out["status"]
            .fillna("")
            .astype(str)
            .ne("ok")
        )
    else:
        out["model_not_ok"] = True

    # Conservative convenience flag.
    # PH p-values are reported separately rather than making a model fail this
    # QC flag, because PH violations should be interpreted scientifically.
    out["qc_clean_epp5"] = ~(
        out["model_not_ok"]
        | out["has_warnings"]
        | out["vif_gt_5"]
        | out["epp_lt_5"]
        | out["bootstrap_success_lt_0_80"]
    )

    return out


def compact_results(df: pd.DataFrame) -> pd.DataFrame:
    """Create a concise table for scientific review."""
    columns = [
        "model_index",
        "model_number",
        "model_id",
        "exposure_organ",
        "exposure_modality",
        "exposure_endpoint",
        "mediator_organ",
        "analysis_subset",
        "status",
        "n",
        "deaths",
        "censored",
        "events_per_parameter_direct",
        "n_covariates",

        "a_beta",
        "a_se_hc3",
        "a_p_hc3",

        "b_log_hr",
        "b_hr",
        "b_se",
        "b_p",

        "direct_log_hr",
        "direct_hr",
        "direct_p",

        "total_log_hr",
        "total_hr",
        "total_p",

        "indirect_log_hr",
        "indirect_hr_like",
        "indirect_delta_se",
        "indirect_delta_z",
        "indirect_delta_p",

        "indirect_bonferroni_p",
        "indirect_bonferroni_significant",

        "indirect_bh_fdr",
        "indirect_bh_fdr_significant",

        "boot_ci_low",
        "boot_ci_high",
        "boot_p_empirical",
        "bootstrap_ci_excludes_zero",

        "bootstrap_successes",
        "bootstrap_requested",
        "bootstrap_success_fraction",

        "max_vif",
        "vif_gt_5",

        "ph_exposure_p",
        "ph_mediator_p",
        "ph_exposure_violation_p_lt_0_05",
        "ph_mediator_violation_p_lt_0_05",

        "cox_direct_concordance",
        "cox_total_concordance",

        "epp_lt_5",
        "epp_lt_10",
        "has_warnings",
        "qc_clean_epp5",
    ]

    return df[
        [c for c in columns if c in df.columns]
    ].copy()


def qc_results(df: pd.DataFrame) -> pd.DataFrame:
    """Create dedicated model-QC table."""
    columns = [
        "model_index",
        "model_number",
        "model_id",
        "exposure_organ",
        "exposure_modality",
        "exposure_endpoint",
        "mediator_organ",
        "status",
        "message",

        "n",
        "deaths",
        "censored",
        "n_covariates",
        "events_per_parameter_direct",

        "epp_lt_5",
        "epp_lt_10",

        "max_vif",
        "vif_gt_5",

        "ph_exposure_p",
        "ph_mediator_p",
        "ph_exposure_violation_p_lt_0_05",
        "ph_mediator_violation_p_lt_0_05",

        "bootstrap_requested",
        "bootstrap_successes",
        "bootstrap_success_fraction",
        "bootstrap_success_lt_0_80",

        "has_warnings",
        "warnings",
        "model_not_ok",
        "qc_clean_epp5",
        "source_file",
    ]

    return df[
        [c for c in columns if c in df.columns]
    ].copy()


# =============================================================================
# MAIN COLLECTION
# =============================================================================

JOB_DIR = JOB_DIR.expanduser().resolve()
RESULT_DIR = JOB_DIR / "single_model_results"
OUTPUT_DIR = OUTPUT_DIR.expanduser().resolve()

print("=" * 80)
print("Collecting EPOCH single-model mediation results")
print("=" * 80)
print(f"Job directory:    {JOB_DIR}")
print(f"Results directory:{RESULT_DIR}")
print(f"Output directory: {OUTPUT_DIR}")
print(f"Expected models:  {EXPECTED_MODELS}")
print("=" * 80)

if not RESULT_DIR.is_dir():
    raise FileNotFoundError(
        f"single_model_results directory not found:\n{RESULT_DIR}"
    )

files = list(RESULT_DIR.glob("model_*.tsv"))

files = sorted(
    files,
    key=lambda p: (
        model_index_from_filename(p)
        if model_index_from_filename(p) is not None
        else 10**9
    ),
)

if len(files) == 0:
    raise FileNotFoundError(
        f"No model_*.tsv files found in:\n{RESULT_DIR}"
    )

print(f"Found {len(files)} model result files.")

frames = []
read_errors = []

for i, path in enumerate(files, start=1):

    try:
        frames.append(
            read_one_model(path)
        )

    except Exception as exc:
        read_errors.append(
            {
                "source_file": str(path),
                "error": str(exc),
            }
        )

        print(
            f"WARNING: failed to read "
            f"{path.name}: {exc}"
        )

    if i % 50 == 0 or i == len(files):
        print(
            f"Read {i}/{len(files)} files"
        )

if len(frames) == 0:
    raise RuntimeError(
        "No valid mediation result files were read."
    )

results = pd.concat(
    frames,
    ignore_index=True,
    sort=False,
)

results["model_index"] = pd.to_numeric(
    results["model_index"],
    errors="coerce",
).astype("Int64")


# =============================================================================
# CHECK MODEL INDICES
# =============================================================================

duplicate_mask = results[
    "model_index"
].duplicated(
    keep=False
)

duplicate_indices = sorted(
    results.loc[
        duplicate_mask,
        "model_index",
    ]
    .dropna()
    .astype(int)
    .unique()
    .tolist()
)

if duplicate_indices:
    raise ValueError(
        "Duplicate model indices detected: "
        f"{duplicate_indices}"
    )

present_indices = set(
    results[
        "model_index"
    ]
    .dropna()
    .astype(int)
    .tolist()
)

expected_indices = set(
    range(EXPECTED_MODELS)
)

missing_indices = sorted(
    expected_indices
    - present_indices
)

unexpected_indices = sorted(
    present_indices
    - expected_indices
)

if unexpected_indices:
    raise ValueError(
        "Unexpected model indices outside "
        f"0-{EXPECTED_MODELS - 1}: "
        f"{unexpected_indices}"
    )

print(
    f"Unique model indices collected: "
    f"{len(present_indices)}/{EXPECTED_MODELS}"
)

if missing_indices:
    print(
        f"Missing model indices "
        f"({len(missing_indices)}): "
        f"{missing_indices}"
    )

    if not ALLOW_INCOMPLETE:
        raise RuntimeError(
            "The 315-model run is incomplete. "
            "Set ALLOW_INCOMPLETE = True at "
            "the top of the script if you "
            "intentionally want partial results."
        )


# =============================================================================
# GLOBAL STATISTICS
# =============================================================================

results = results.sort_values(
    "model_index"
).reset_index(
    drop=True
)

results = add_global_statistics(
    results
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# =============================================================================
# WRITE TABLES
# =============================================================================

all_results_path = (
    OUTPUT_DIR
    / "OLS_Cox_bootstrap_mediation_all_results.tsv"
)

compact_path = (
    OUTPUT_DIR
    / "OLS_Cox_bootstrap_mediation_compact_results.tsv"
)

qc_path = (
    OUTPUT_DIR
    / "OLS_Cox_bootstrap_mediation_QC.tsv"
)

bonf_path = (
    OUTPUT_DIR
    / "OLS_Cox_bootstrap_mediation_Bonferroni_significant.tsv"
)

fdr_path = (
    OUTPUT_DIR
    / "OLS_Cox_bootstrap_mediation_BH_FDR_significant.tsv"
)

bonf_boot_path = (
    OUTPUT_DIR
    / "OLS_Cox_bootstrap_mediation_Bonferroni_plus_bootstrap.tsv"
)

bonf_boot_qc_path = (
    OUTPUT_DIR
    / "OLS_Cox_bootstrap_mediation_Bonferroni_bootstrap_QCclean.tsv"
)

missing_path = (
    OUTPUT_DIR
    / "missing_model_indices.tsv"
)

status_path = (
    OUTPUT_DIR
    / "model_status_summary.tsv"
)


results.to_csv(
    all_results_path,
    sep="\t",
    index=False,
)

compact_results(
    results
).to_csv(
    compact_path,
    sep="\t",
    index=False,
)

qc_results(
    results
).to_csv(
    qc_path,
    sep="\t",
    index=False,
)


bonf_mask = (
    results[
        "indirect_bonferroni_significant"
    ]
)

fdr_mask = (
    results[
        "indirect_bh_fdr_significant"
    ]
)

bootstrap_mask = (
    results[
        "bootstrap_ci_excludes_zero"
    ]
)

qc_clean_mask = (
    results[
        "qc_clean_epp5"
    ]
)


results.loc[
    bonf_mask
].to_csv(
    bonf_path,
    sep="\t",
    index=False,
)

results.loc[
    fdr_mask
].to_csv(
    fdr_path,
    sep="\t",
    index=False,
)

results.loc[
    bonf_mask
    & bootstrap_mask
].to_csv(
    bonf_boot_path,
    sep="\t",
    index=False,
)

results.loc[
    bonf_mask
    & bootstrap_mask
    & qc_clean_mask
].to_csv(
    bonf_boot_qc_path,
    sep="\t",
    index=False,
)


pd.DataFrame(
    {
        "missing_model_index":
        missing_indices
    }
).to_csv(
    missing_path,
    sep="\t",
    index=False,
)


if "status" in results.columns:

    status_summary = (
        results
        .groupby(
            "status",
            dropna=False,
        )
        .size()
        .reset_index(
            name="n_models"
        )
        .sort_values(
            "n_models",
            ascending=False,
        )
    )

else:

    status_summary = pd.DataFrame(
        {
            "status": ["unknown"],
            "n_models": [len(results)],
        }
    )


status_summary.to_csv(
    status_path,
    sep="\t",
    index=False,
)


if read_errors:

    pd.DataFrame(
        read_errors
    ).to_csv(
        OUTPUT_DIR
        / "unreadable_result_files.tsv",
        sep="\t",
        index=False,
    )


# =============================================================================
# SUMMARY
# =============================================================================

status_counts = (
    results[
        "status"
    ]
    .value_counts(
        dropna=False
    )
    .to_dict()
    if "status" in results.columns
    else {}
)

summary = {
    "job_dir": str(JOB_DIR),
    "result_dir": str(RESULT_DIR),
    "expected_models": EXPECTED_MODELS,
    "collected_models": int(len(results)),
    "unique_model_indices": int(
        results["model_index"].nunique()
    ),
    "missing_model_count": len(
        missing_indices
    ),
    "missing_model_indices": (
        missing_indices
    ),
    "alpha": ALPHA,
    "bonferroni_threshold": (
        ALPHA / EXPECTED_MODELS
    ),
    "status_counts": (
        status_counts
    ),
    "bonferroni_significant_n": int(
        bonf_mask.sum()
    ),
    "bh_fdr_significant_n": int(
        fdr_mask.sum()
    ),
    "bootstrap_ci_excludes_zero_n": int(
        bootstrap_mask.sum()
    ),
    "bonferroni_plus_bootstrap_n": int(
        (
            bonf_mask
            & bootstrap_mask
        ).sum()
    ),
    "bonferroni_bootstrap_qcclean_n": int(
        (
            bonf_mask
            & bootstrap_mask
            & qc_clean_mask
        ).sum()
    ),
    "qc_clean_epp5_n": int(
        qc_clean_mask.sum()
    ),
    "epp_lt_5_n": int(
        results[
            "epp_lt_5"
        ].sum()
    ),
    "epp_lt_10_n": int(
        results[
            "epp_lt_10"
        ].sum()
    ),
    "vif_gt_5_n": int(
        results[
            "vif_gt_5"
        ].sum()
    ),
    "models_with_warnings_n": int(
        results[
            "has_warnings"
        ].sum()
    ),
    "ph_exposure_violation_n": int(
        results[
            "ph_exposure_violation_p_lt_0_05"
        ].sum()
    ),
    "ph_mediator_violation_n": int(
        results[
            "ph_mediator_violation_p_lt_0_05"
        ].sum()
    ),
}


with open(
    OUTPUT_DIR
    / "collection_summary.json",
    "w",
) as f:

    json.dump(
        summary,
        f,
        indent=2,
        default=str,
    )


# =============================================================================
# PRINT SUMMARY
# =============================================================================

print()
print("=" * 80)
print("COLLECTION COMPLETED")
print("=" * 80)

print(
    f"Models collected: "
    f"{len(results)}/{EXPECTED_MODELS}"
)

print(
    f"Status counts: "
    f"{status_counts}"
)

print(
    f"Bonferroni threshold: "
    f"{ALPHA / EXPECTED_MODELS:.10g}"
)

print(
    "Bonferroni-significant "
    "indirect effects: "
    f"{int(bonf_mask.sum())}"
)

print(
    "BH-FDR-significant "
    "indirect effects: "
    f"{int(fdr_mask.sum())}"
)

print(
    "Bootstrap 95% CI "
    "excludes zero: "
    f"{int(bootstrap_mask.sum())}"
)

print(
    "Bonferroni significant "
    "+ bootstrap CI excludes zero: "
    f"{int((bonf_mask & bootstrap_mask).sum())}"
)

print(
    "Bonferroni significant "
    "+ bootstrap supported "
    "+ QC clean (EPP >= 5): "
    f"{int((bonf_mask & bootstrap_mask & qc_clean_mask).sum())}"
)

print(
    "Models with EPP < 5: "
    f"{int(results['epp_lt_5'].sum())}"
)

print(
    "Models with EPP < 10: "
    f"{int(results['epp_lt_10'].sum())}"
)

print(
    "Models with VIF > 5: "
    f"{int(results['vif_gt_5'].sum())}"
)

print(
    "Models with warnings: "
    f"{int(results['has_warnings'].sum())}"
)

print()
print(
    f"Results written to:\n"
    f"{OUTPUT_DIR}"
)

print("=" * 80)