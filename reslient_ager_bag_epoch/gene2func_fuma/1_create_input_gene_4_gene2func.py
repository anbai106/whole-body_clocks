#!/usr/bin/env python3
"""
GTEx v8-only directional eQTL gene sets for EPOCH-BAG residual resilience.

This script reads FUMA `eqtl.txt` directly and performs ONLY GTEx v8 analyses.

Analyses
========
PRIMARY:
    GTEx v8 brain tissues only
    -> strict directional resilience-associated expression genes
    -> matched GTEx v8 brain directional-eQTL background
    -> FUMA GENE2FUNC-ready gene lists

SECONDARY SYSTEMIC:
    All GTEx v8 tissues
    -> strict directional resilience-associated expression genes
    -> matched all-GTEx-v8 directional-eQTL background
    -> FUMA GENE2FUNC-ready gene lists

No GTEx v6/v7, BRAINEAC, CMC, PsychENCODE, BrainSeq, eQTLGen, or other
eQTL resources are used anywhere in the analysis.

Biological interpretation
=========================
EPOCH-BAG residual =
    observed mortality EPOCH - EPOCH expected from BAG

FUMA `alignedDirection` is oriented to the GWAS phenotype-increasing allele:

    alignedDirection == "-"
        GWAS phenotype-increasing allele decreases gene expression.

For the EPOCH-BAG residual, the phenotype-increasing allele is the
residual-increasing / vulnerability-like allele. Therefore:

    alignedDirection == "-"
        -> higher gene expression is aligned with LOWER EPOCH-BAG residual
        -> resilience-like expression direction

    alignedDirection == "+"
        -> higher gene expression is aligned with HIGHER EPOCH-BAG residual
        -> vulnerability-like expression direction

Direction is first collapsed at the gene x GTEx-v8-tissue level
===============================================================
Many correlated SNPs can map the same gene in the same tissue. They are not
treated as independent directional votes.

Each gene x tissue unit is classified as:

    negative_only:
        every informative GTEx v8 eQTL row in the tissue has
        alignedDirection == "-"

    positive_only:
        every informative row has alignedDirection == "+"

    mixed:
        both "-" and "+" are observed within that gene x tissue

Strict resilience-associated expression gene
=============================================
Within the analysis scope (brain-only or all GTEx v8 tissues), a gene is
included in the strict resilience foreground when:

    N_negative_only_tissues >= 1
    N_positive_only_tissues == 0
    N_mixed_tissues == 0

Thus every informative GTEx v8 tissue for that gene points in the
resilience-like expression direction.

Matched GENE2FUNC background
============================
The background is NOT all protein-coding genes.

For each analysis, the matched background is:
    all unique genes with
        db == "GTEx_v8"
        eqtlMapFilt == 1
        valid gene symbol
        alignedDirection in {"+", "-"}
    within the same tissue scope.

This is the set of genes that could have qualified for the foreground.

Default paths
=============
Input:
    /Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/mortality_clock/
    fuma/Brain_proteomics_mortality_clock/EPOCH_BAG_residual/eqtl.txt

Output:
    .../EPOCH_BAG_residual/directional_eqtl_gene_sets_GTEx_v8/

Validated on FUMA_job763388
===========================
Using the definitions above, the uploaded FUMA archive gives:

    PRIMARY GTEx v8 brain:
        matched background = 34 genes
        strict resilience foreground = 16 genes

    SECONDARY all GTEx v8 tissues:
        matched background = 131 genes
        strict resilience foreground = 40 genes

These counts are printed at runtime and are not hard-coded.

Usage
=====
    python fuma_eqtl_directional_gene_sets_GTEx_v8.py --overwrite

Optional explicit paths:
    python fuma_eqtl_directional_gene_sets_GTEx_v8.py \
        --eqtl /path/to/eqtl.txt \
        --outdir /path/to/output \
        --overwrite
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, List, Tuple

import numpy as np
import pandas as pd


# =============================================================================
# CONSTANTS
# =============================================================================

GTEX_DATABASE = "GTEx_v8"

DEFAULT_FUMA_DIR = Path(
    "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/mortality_clock/"
    "fuma/Brain_proteomics_mortality_clock/EPOCH_BAG_residual"
)

DEFAULT_EQTL = DEFAULT_FUMA_DIR / "eqtl.txt"

DEFAULT_OUTDIR = (
    DEFAULT_FUMA_DIR / "directional_eqtl_gene_sets_GTEx_v8"
)

REQUIRED_COLUMNS = {
    "uniqID",
    "db",
    "tissue",
    "gene",
    "p",
    "FDR",
    "alignedDirection",
    "symbol",
    "eqtlMapFilt",
}


# =============================================================================
# COMMAND LINE
# =============================================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create GTEx v8-only directional eQTL resilience gene sets "
            "from FUMA eqtl.txt."
        )
    )

    parser.add_argument(
        "--eqtl",
        type=Path,
        default=DEFAULT_EQTL,
        help=f"FUMA eqtl.txt. Default: {DEFAULT_EQTL}",
    )

    parser.add_argument(
        "--outdir",
        type=Path,
        default=DEFAULT_OUTDIR,
        help=f"Output directory. Default: {DEFAULT_OUTDIR}",
    )

    parser.add_argument(
        "--min-tissues-for-80pct",
        type=int,
        default=2,
        help=(
            "Minimum number of non-mixed informative tissues required for "
            "the secondary >=80%% resilience-consistency flag. Default: 2."
        ),
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing output files.",
    )

    return parser.parse_args()


# =============================================================================
# BASIC HELPERS
# =============================================================================

def clean_string_series(series: pd.Series) -> pd.Series:
    return series.astype("string").str.strip()


def valid_symbol_mask(series: pd.Series) -> pd.Series:
    x = clean_string_series(series)

    bad = x.str.lower().isin(
        {
            "",
            "na",
            "nan",
            "none",
            "null",
            ".",
        }
    )

    return x.notna() & ~bad


def safe_min_numeric(values: pd.Series) -> float:
    x = pd.to_numeric(
        values,
        errors="coerce",
    )

    x = x[np.isfinite(x)]

    if len(x) == 0:
        return np.nan

    return float(x.min())


def collapse_unique_strings(values: Iterable[object]) -> str:
    cleaned = []

    for value in values:
        if pd.isna(value):
            continue

        text = str(value).strip()

        if (
            not text
            or text.lower()
            in {
                "na",
                "nan",
                "none",
                "null",
                ".",
            }
        ):
            continue

        cleaned.append(text)

    return ";".join(
        sorted(
            set(cleaned)
        )
    )


def write_gene_list(
    genes: Iterable[str],
    path: Path,
    overwrite: bool,
) -> int:
    """
    One gene symbol per line, no header.
    This is directly suitable for FUMA GENE2FUNC input.
    """
    genes = sorted(
        {
            str(gene).strip()
            for gene in genes
            if str(gene).strip()
        }
    )

    if path.exists() and not overwrite:
        raise FileExistsError(
            f"Output exists: {path}\n"
            "Re-run with --overwrite to replace it."
        )

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    path.write_text(
        "".join(
            f"{gene}\n"
            for gene in genes
        )
    )

    return len(genes)


def write_tsv(
    df: pd.DataFrame,
    path: Path,
    overwrite: bool,
) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(
            f"Output exists: {path}\n"
            "Re-run with --overwrite to replace it."
        )

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    df.to_csv(
        path,
        sep="\t",
        index=False,
        na_rep="NA",
    )


# =============================================================================
# READ AND FILTER FUMA eqtl.txt
# =============================================================================

def read_eqtl(
    path: Path,
) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(
            f"FUMA eqtl.txt not found: {path}"
        )

    if path.stat().st_size == 0:
        raise ValueError(
            f"FUMA eqtl.txt is empty: {path}"
        )

    print(
        f"Reading FUMA eQTL file:\n  {path}",
        flush=True,
    )

    df = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        low_memory=False,
    )

    missing = sorted(
        REQUIRED_COLUMNS
        - set(df.columns)
    )

    if missing:
        raise ValueError(
            "FUMA eqtl.txt is missing required column(s): "
            + ", ".join(missing)
        )

    return df


def prepare_gtex_v8_directional_rows(
    df: pd.DataFrame,
) -> Tuple[pd.DataFrame, dict]:
    """
    Apply the ONLY allowed database filter: GTEx_v8.

    Then require:
      - eqtlMapFilt == 1
      - valid gene symbol
      - alignedDirection exactly "+" or "-"

    Rows lacking directional information are excluded from both foreground
    and background construction.
    """
    out = df.copy()

    for col in [
        "uniqID",
        "db",
        "tissue",
        "gene",
        "symbol",
        "alignedDirection",
        "eqtlMapFilt",
    ]:
        out[col] = clean_string_series(
            out[col]
        )

    provenance = {
        "database_constraint": GTEX_DATABASE,
        "N_rows_raw_eqtl": int(
            len(out)
        ),
    }

    db_mask = (
        out["db"]
        == GTEX_DATABASE
    )

    provenance["N_rows_GTEx_v8"] = int(
        db_mask.sum()
    )

    mapped_mask = out[
        "eqtlMapFilt"
    ].isin(
        {
            "1",
            "1.0",
            "TRUE",
            "True",
            "true",
        }
    )

    symbol_mask = valid_symbol_mask(
        out["symbol"]
    )

    direction_mask = out[
        "alignedDirection"
    ].isin(
        {
            "+",
            "-",
        }
    )

    keep = (
        db_mask
        & mapped_mask
        & symbol_mask
        & direction_mask
    )

    out = out.loc[
        keep
    ].copy()

    if len(out) == 0:
        available_dbs = sorted(
            df["db"]
            .dropna()
            .astype(str)
            .str.strip()
            .unique()
            .tolist()
        )

        raise ValueError(
            "No directionally informative mapped GTEx_v8 eQTL rows "
            "were found.\nAvailable FUMA db labels include:\n  "
            + "\n  ".join(
                available_dbs
            )
        )

    out["p_numeric"] = pd.to_numeric(
        out["p"],
        errors="coerce",
    )

    out["FDR_numeric"] = pd.to_numeric(
        out["FDR"],
        errors="coerce",
    )

    provenance[
        "N_directional_mapped_GTEx_v8_rows"
    ] = int(
        len(out)
    )

    provenance[
        "N_directional_GTEx_v8_genes"
    ] = int(
        out["symbol"].nunique()
    )

    provenance[
        "N_directional_GTEx_v8_tissues"
    ] = int(
        out["tissue"].nunique()
    )

    return (
        out,
        provenance,
    )


# =============================================================================
# TISSUE SCOPE
# =============================================================================

def get_primary_brain_rows(
    gtex_v8_rows: pd.DataFrame,
) -> pd.DataFrame:
    """
    Primary scope:
        GTEx_v8 rows whose FUMA tissue label begins with "Brain_".
    """
    mask = (
        gtex_v8_rows["tissue"]
        .fillna("")
        .str.startswith(
            "Brain_"
        )
    )

    brain = gtex_v8_rows.loc[
        mask
    ].copy()

    if len(brain) == 0:
        raise ValueError(
            "No GTEx_v8 Brain_* rows were found."
        )

    return brain


def get_secondary_all_tissue_rows(
    gtex_v8_rows: pd.DataFrame,
) -> pd.DataFrame:
    """
    Secondary systemic scope:
        every directionally informative mapped GTEx_v8 tissue.
    """
    return gtex_v8_rows.copy()


# =============================================================================
# GENE x TISSUE COLLAPSE
# =============================================================================

def classify_tissue_unit_direction(
    directions: pd.Series,
) -> str:
    observed = set(
        directions
        .dropna()
        .astype(str)
    )

    if observed == {"-"}:
        return "negative_only"

    if observed == {"+"}:
        return "positive_only"

    if observed == {
        "+",
        "-",
    }:
        return "mixed"

    return "uninformative"


def build_gene_tissue_units(
    rows: pd.DataFrame,
) -> pd.DataFrame:
    """
    Collapse all SNP-eQTL rows into one directional unit per gene x tissue.

    Because every row is already GTEx_v8, database is constant and does not
    contribute an additional independent evidence dimension.
    """
    columns = [
        "symbol",
        "db",
        "tissue",
        "unit_direction",
        "N_eQTL_rows",
        "N_unique_SNPs",
        "N_negative_rows",
        "N_positive_rows",
        "min_eQTL_p",
        "min_eQTL_FDR",
        "Ensembl_gene_ids",
    ]

    if len(rows) == 0:
        return pd.DataFrame(
            columns=columns
        )

    records = []

    grouped = rows.groupby(
        [
            "symbol",
            "tissue",
        ],
        dropna=False,
        sort=True,
    )

    for (
        symbol,
        tissue,
    ), group in grouped:

        directions = group[
            "alignedDirection"
        ]

        records.append(
            {
                "symbol": symbol,
                "db": GTEX_DATABASE,
                "tissue": tissue,
                "unit_direction": classify_tissue_unit_direction(
                    directions
                ),
                "N_eQTL_rows": int(
                    len(group)
                ),
                "N_unique_SNPs": int(
                    group[
                        "uniqID"
                    ].nunique(
                        dropna=True
                    )
                ),
                "N_negative_rows": int(
                    (
                        directions
                        == "-"
                    ).sum()
                ),
                "N_positive_rows": int(
                    (
                        directions
                        == "+"
                    ).sum()
                ),
                "min_eQTL_p": safe_min_numeric(
                    group[
                        "p_numeric"
                    ]
                ),
                "min_eQTL_FDR": safe_min_numeric(
                    group[
                        "FDR_numeric"
                    ]
                ),
                "Ensembl_gene_ids": collapse_unique_strings(
                    group[
                        "gene"
                    ]
                ),
            }
        )

    return pd.DataFrame.from_records(
        records,
        columns=columns,
    )


# =============================================================================
# GENE-LEVEL DIRECTION CONSISTENCY
# =============================================================================

def build_gene_summary(
    rows: pd.DataFrame,
    tissue_units: pd.DataFrame,
    scope_name: str,
    min_tissues_for_80pct: int,
) -> pd.DataFrame:
    """
    One row per gene.

    direction_score =
        (negative-only tissues - positive-only tissues)
        / (negative-only tissues + positive-only tissues)

    Positive score:
        expression direction is more resilience-like.

    Negative score:
        expression direction is more vulnerability-like.

    Mixed tissues are explicitly retained as discordant evidence and are
    excluded from the denominator of direction_score.
    """
    if len(tissue_units) == 0:
        return pd.DataFrame()

    gene_records = []

    for (
        symbol,
        group,
    ) in tissue_units.groupby(
        "symbol",
        sort=True,
    ):

        n_negative = int(
            (
                group[
                    "unit_direction"
                ]
                == "negative_only"
            ).sum()
        )

        n_positive = int(
            (
                group[
                    "unit_direction"
                ]
                == "positive_only"
            ).sum()
        )

        n_mixed = int(
            (
                group[
                    "unit_direction"
                ]
                == "mixed"
            ).sum()
        )

        n_nonmixed = (
            n_negative
            + n_positive
        )

        if n_nonmixed > 0:
            resilience_fraction = (
                n_negative
                / n_nonmixed
            )

            vulnerability_fraction = (
                n_positive
                / n_nonmixed
            )

            direction_score = (
                n_negative
                - n_positive
            ) / n_nonmixed

        else:
            resilience_fraction = np.nan
            vulnerability_fraction = np.nan
            direction_score = np.nan

        strict_resilience = (
            n_negative >= 1
            and n_positive == 0
            and n_mixed == 0
        )

        strict_vulnerability = (
            n_positive >= 1
            and n_negative == 0
            and n_mixed == 0
        )

        strict_resilience_multitissue = (
            n_negative >= 2
            and n_positive == 0
            and n_mixed == 0
        )

        predominant_resilience_80pct = (
            n_nonmixed
            >= min_tissues_for_80pct
            and resilience_fraction
            >= 0.80
        )

        predominant_vulnerability_80pct = (
            n_nonmixed
            >= min_tissues_for_80pct
            and vulnerability_fraction
            >= 0.80
        )

        gene_records.append(
            {
                "scope": scope_name,
                "database": GTEX_DATABASE,
                "symbol": symbol,
                "N_informative_tissues": int(
                    len(group)
                ),
                "N_negative_only_tissues": n_negative,
                "N_positive_only_tissues": n_positive,
                "N_mixed_tissues": n_mixed,
                "N_nonmixed_tissues": n_nonmixed,
                "resilience_fraction_nonmixed": resilience_fraction,
                "vulnerability_fraction_nonmixed": vulnerability_fraction,
                "direction_score": direction_score,
                "strict_resilience_expression_gene": strict_resilience,
                "strict_resilience_multitissue_gene": strict_resilience_multitissue,
                "strict_vulnerability_expression_gene": strict_vulnerability,
                "predominant_resilience_80pct": predominant_resilience_80pct,
                "predominant_vulnerability_80pct": predominant_vulnerability_80pct,
                "tissues": collapse_unique_strings(
                    group[
                        "tissue"
                    ]
                ),
                "Ensembl_gene_ids": collapse_unique_strings(
                    group[
                        "Ensembl_gene_ids"
                    ]
                ),
            }
        )

    gene_summary = pd.DataFrame.from_records(
        gene_records
    )

    row_records = []

    for (
        symbol,
        group,
    ) in rows.groupby(
        "symbol",
        sort=True,
    ):

        row_records.append(
            {
                "symbol": symbol,
                "N_directional_eQTL_rows": int(
                    len(group)
                ),
                "N_unique_directional_SNPs": int(
                    group[
                        "uniqID"
                    ].nunique(
                        dropna=True
                    )
                ),
                "N_negative_rows": int(
                    (
                        group[
                            "alignedDirection"
                        ]
                        == "-"
                    ).sum()
                ),
                "N_positive_rows": int(
                    (
                        group[
                            "alignedDirection"
                        ]
                        == "+"
                    ).sum()
                ),
                "min_eQTL_p": safe_min_numeric(
                    group[
                        "p_numeric"
                    ]
                ),
                "min_eQTL_FDR": safe_min_numeric(
                    group[
                        "FDR_numeric"
                    ]
                ),
            }
        )

    row_summary = pd.DataFrame.from_records(
        row_records
    )

    gene_summary = gene_summary.merge(
        row_summary,
        on="symbol",
        how="left",
        validate="one_to_one",
    )

    def classify_gene(
        row: pd.Series,
    ) -> str:

        if bool(
            row[
                "strict_resilience_expression_gene"
            ]
        ):
            return "strict_resilience_expression"

        if bool(
            row[
                "strict_vulnerability_expression_gene"
            ]
        ):
            return "strict_vulnerability_expression"

        if (
            row[
                "N_mixed_tissues"
            ]
            > 0
        ):
            return "mixed_direction"

        if (
            row[
                "direction_score"
            ]
            > 0
        ):
            return "predominantly_resilience_like"

        if (
            row[
                "direction_score"
            ]
            < 0
        ):
            return "predominantly_vulnerability_like"

        return "balanced_or_unclassified"

    gene_summary[
        "gene_direction_class"
    ] = gene_summary.apply(
        classify_gene,
        axis=1,
    )

    gene_summary = gene_summary.sort_values(
        by=[
            "strict_resilience_expression_gene",
            "direction_score",
            "N_informative_tissues",
            "N_unique_directional_SNPs",
            "symbol",
        ],
        ascending=[
            False,
            False,
            False,
            False,
            True,
        ],
        kind="mergesort",
    ).reset_index(
        drop=True
    )

    return gene_summary


# =============================================================================
# GENE LISTS
# =============================================================================

def strict_resilience_genes(
    gene_summary: pd.DataFrame,
) -> List[str]:

    if len(gene_summary) == 0:
        return []

    return (
        gene_summary.loc[
            gene_summary[
                "strict_resilience_expression_gene"
            ].astype(bool),
            "symbol",
        ]
        .dropna()
        .astype(str)
        .sort_values()
        .tolist()
    )


def matched_background_genes(
    gene_summary: pd.DataFrame,
) -> List[str]:
    """
    All directionally informative genes in the SAME GTEx v8 tissue scope.
    """
    if len(gene_summary) == 0:
        return []

    return sorted(
        gene_summary[
            "symbol"
        ]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )


# =============================================================================
# ANALYSIS WRAPPER
# =============================================================================

def analyze_scope(
    rows: pd.DataFrame,
    scope_name: str,
    min_tissues_for_80pct: int,
) -> dict:

    units = build_gene_tissue_units(
        rows
    )

    genes = build_gene_summary(
        rows=rows,
        tissue_units=units,
        scope_name=scope_name,
        min_tissues_for_80pct=min_tissues_for_80pct,
    )

    foreground = strict_resilience_genes(
        genes
    )

    background = matched_background_genes(
        genes
    )

    return {
        "rows": rows,
        "units": units,
        "genes": genes,
        "foreground": foreground,
        "background": background,
    }


# =============================================================================
# README / SUMMARY
# =============================================================================

def build_summary(
    primary: dict,
    secondary: dict,
    provenance: dict,
) -> pd.DataFrame:

    rows = [
        {
            "analysis_priority": "Primary",
            "scope": "GTEx_v8_brain_only",
            "database": GTEX_DATABASE,
            "tissue_definition": 'tissue starts with "Brain_"',
            "N_directional_eQTL_rows": len(
                primary["rows"]
            ),
            "N_tissues": primary[
                "rows"
            ][
                "tissue"
            ].nunique(),
            "N_gene_tissue_units": len(
                primary[
                    "units"
                ]
            ),
            "N_matched_background_genes": len(
                primary[
                    "background"
                ]
            ),
            "N_strict_resilience_genes": len(
                primary[
                    "foreground"
                ]
            ),
            "N_strict_resilience_multitissue_genes": int(
                primary[
                    "genes"
                ][
                    "strict_resilience_multitissue_gene"
                ].sum()
            ),
            "N_strict_vulnerability_genes": int(
                primary[
                    "genes"
                ][
                    "strict_vulnerability_expression_gene"
                ].sum()
            ),
        },
        {
            "analysis_priority": "Secondary",
            "scope": "GTEx_v8_all_tissues",
            "database": GTEX_DATABASE,
            "tissue_definition": "all informative GTEx_v8 tissues",
            "N_directional_eQTL_rows": len(
                secondary["rows"]
            ),
            "N_tissues": secondary[
                "rows"
            ][
                "tissue"
            ].nunique(),
            "N_gene_tissue_units": len(
                secondary[
                    "units"
                ]
            ),
            "N_matched_background_genes": len(
                secondary[
                    "background"
                ]
            ),
            "N_strict_resilience_genes": len(
                secondary[
                    "foreground"
                ]
            ),
            "N_strict_resilience_multitissue_genes": int(
                secondary[
                    "genes"
                ][
                    "strict_resilience_multitissue_gene"
                ].sum()
            ),
            "N_strict_vulnerability_genes": int(
                secondary[
                    "genes"
                ][
                    "strict_vulnerability_expression_gene"
                ].sum()
            ),
        },
    ]

    summary = pd.DataFrame(
        rows
    )

    for key, value in provenance.items():
        summary[key] = value

    return summary


def build_overlap_table(
    primary: dict,
    secondary: dict,
) -> pd.DataFrame:

    primary_genes = set(
        primary[
            "foreground"
        ]
    )

    secondary_genes = set(
        secondary[
            "foreground"
        ]
    )

    all_genes = sorted(
        primary_genes
        | secondary_genes
    )

    records = []

    for gene in all_genes:
        in_primary = (
            gene
            in primary_genes
        )

        in_secondary = (
            gene
            in secondary_genes
        )

        if (
            in_primary
            and in_secondary
        ):
            category = (
                "brain_and_all_tissue_strict"
            )

        elif in_primary:
            category = (
                "brain_strict_only"
            )

        else:
            category = (
                "all_tissue_strict_only"
            )

        records.append(
            {
                "symbol": gene,
                "GTEx_v8_brain_strict_resilience": in_primary,
                "GTEx_v8_all_tissues_strict_resilience": in_secondary,
                "overlap_category": category,
            }
        )

    return pd.DataFrame.from_records(
        records
    )


def write_readme(
    path: Path,
    overwrite: bool,
    summary: pd.DataFrame,
    primary: dict,
    secondary: dict,
) -> None:

    if path.exists() and not overwrite:
        raise FileExistsError(
            f"Output exists: {path}\n"
            "Re-run with --overwrite to replace it."
        )

    primary_row = summary.loc[
        summary[
            "scope"
        ]
        == "GTEx_v8_brain_only"
    ].iloc[0]

    secondary_row = summary.loc[
        summary[
            "scope"
        ]
        == "GTEx_v8_all_tissues"
    ].iloc[0]

    primary_set = set(
        primary[
            "foreground"
        ]
    )

    secondary_set = set(
        secondary[
            "foreground"
        ]
    )

    overlap = (
        primary_set
        & secondary_set
    )

    text = f"""GTEx v8-only directional eQTL resilience gene sets
==================================================

DATABASE CONSTRAINT
-------------------
Every analysis in this directory uses ONLY:
    db == "GTEx_v8"

No other eQTL database contributes to foregrounds, backgrounds,
direction scores, or tissue-consistency calculations.

EPOCH-BAG RESIDUAL
------------------
Residual = observed mortality EPOCH - EPOCH expected from BAG.

FUMA DIRECTION
--------------
alignedDirection = "-"
    The GWAS residual-increasing allele decreases gene expression.

Therefore, higher expression is aligned with LOWER EPOCH-BAG residual,
which is the resilience-like direction.

alignedDirection = "+"
    Higher expression is aligned with HIGHER EPOCH-BAG residual,
which is the vulnerability-like direction.

PRIMARY ANALYSIS: GTEx v8 BRAIN ONLY
------------------------------------
Tissue definition:
    GTEx_v8 rows with tissue names beginning "Brain_"

Number of brain tissues:
    {int(primary_row["N_tissues"])}

Matched directional-eQTL background:
    {int(primary_row["N_matched_background_genes"])} genes

Strict resilience-associated expression foreground:
    {int(primary_row["N_strict_resilience_genes"])} genes

FUMA GENE2FUNC files:
    GENE2FUNC_foreground_strict_resilience_GTEx_v8_brain.txt
    GENE2FUNC_background_directional_eqtl_GTEx_v8_brain.txt

SECONDARY SYSTEMIC ANALYSIS: ALL GTEx v8 TISSUES
------------------------------------------------
Tissue definition:
    every directionally informative mapped GTEx_v8 tissue

Number of tissues:
    {int(secondary_row["N_tissues"])}

Matched directional-eQTL background:
    {int(secondary_row["N_matched_background_genes"])} genes

Strict resilience-associated expression foreground:
    {int(secondary_row["N_strict_resilience_genes"])} genes

FUMA GENE2FUNC files:
    GENE2FUNC_foreground_strict_resilience_GTEx_v8_all_tissues.txt
    GENE2FUNC_background_directional_eqtl_GTEx_v8_all_tissues.txt

STRICT GENE DEFINITION
----------------------
The script first collapses SNP rows into one gene x GTEx-v8-tissue unit.

A gene is a strict resilience-associated expression gene if:
    N_negative_only_tissues >= 1
    N_positive_only_tissues == 0
    N_mixed_tissues == 0

Thus all informative tissues for that gene within the analysis scope
point toward the resilience-like expression direction.

MATCHED BACKGROUND
------------------
The background is not all protein-coding genes.

It contains every gene with:
    db == "GTEx_v8"
    eqtlMapFilt == 1
    valid symbol
    alignedDirection in {{"+", "-"}}

within the same tissue scope as the foreground.

This is the correct matched universe of genes that could have entered
that foreground under the directional eQTL selection rule.

PRIMARY vs SECONDARY OVERLAP
----------------------------
Brain-only strict resilience genes:
    {len(primary_set)}

All-GTEx-v8-tissue strict resilience genes:
    {len(secondary_set)}

Strict genes shared by both:
    {len(overlap)}

Interpretation:
    A brain-only strict gene can fail the all-tissue strict definition if
    an informative non-brain GTEx v8 tissue has an opposing or mixed
    direction. This does not invalidate its brain-specific signal; it
    indicates tissue-dependent directionality.

RECOMMENDED MANUSCRIPT USE
--------------------------
Primary:
    Run GENE2FUNC using the GTEx v8 brain strict foreground and its
    GTEx v8 brain matched background.

Secondary systemic:
    Run GENE2FUNC using the all-GTEx-v8 strict foreground and its
    all-GTEx-v8 matched background.

Prefer the term:
    "resilience-associated expression genes"

rather than:
    "resilience genes"

because directional eQTL alignment is associative and does not by itself
establish that increasing expression causally produces resilience.
"""

    path.write_text(
        text
    )


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    args = parse_args()

    if (
        args.min_tissues_for_80pct
        < 1
    ):
        raise ValueError(
            "--min-tissues-for-80pct must be >= 1"
        )

    args.outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    raw = read_eqtl(
        args.eqtl
    )

    gtex_v8_rows, provenance = (
        prepare_gtex_v8_directional_rows(
            raw
        )
    )

    primary_rows = (
        get_primary_brain_rows(
            gtex_v8_rows
        )
    )

    secondary_rows = (
        get_secondary_all_tissue_rows(
            gtex_v8_rows
        )
    )

    print(
        "\nGTEx v8-only input after directional filtering:",
        flush=True,
    )
    print(
        f"  Rows: {len(gtex_v8_rows):,}",
        flush=True,
    )
    print(
        f"  Genes: {gtex_v8_rows['symbol'].nunique():,}",
        flush=True,
    )
    print(
        f"  Tissues: {gtex_v8_rows['tissue'].nunique():,}",
        flush=True,
    )

    primary = analyze_scope(
        rows=primary_rows,
        scope_name="GTEx_v8_brain_only",
        min_tissues_for_80pct=args.min_tissues_for_80pct,
    )

    secondary = analyze_scope(
        rows=secondary_rows,
        scope_name="GTEx_v8_all_tissues",
        min_tissues_for_80pct=args.min_tissues_for_80pct,
    )

    summary = build_summary(
        primary=primary,
        secondary=secondary,
        provenance=provenance,
    )

    overlap = build_overlap_table(
        primary=primary,
        secondary=secondary,
    )

    # -------------------------------------------------------------------------
    # Output paths
    # -------------------------------------------------------------------------
    outputs = {
        # FUMA-ready foreground/background lists
        "primary_foreground": (
            args.outdir
            / "GENE2FUNC_foreground_strict_resilience_GTEx_v8_brain.txt"
        ),
        "primary_background": (
            args.outdir
            / "GENE2FUNC_background_directional_eqtl_GTEx_v8_brain.txt"
        ),
        "secondary_foreground": (
            args.outdir
            / "GENE2FUNC_foreground_strict_resilience_GTEx_v8_all_tissues.txt"
        ),
        "secondary_background": (
            args.outdir
            / "GENE2FUNC_background_directional_eqtl_GTEx_v8_all_tissues.txt"
        ),

        # Gene-level detail
        "primary_gene_summary": (
            args.outdir
            / "eqtl_gene_direction_consistency_GTEx_v8_brain.tsv"
        ),
        "secondary_gene_summary": (
            args.outdir
            / "eqtl_gene_direction_consistency_GTEx_v8_all_tissues.tsv"
        ),

        # Gene x tissue audit
        "primary_units": (
            args.outdir
            / "eqtl_gene_tissue_direction_units_GTEx_v8_brain.tsv"
        ),
        "secondary_units": (
            args.outdir
            / "eqtl_gene_tissue_direction_units_GTEx_v8_all_tissues.tsv"
        ),

        # Tissue audit
        "primary_tissues": (
            args.outdir
            / "GTEx_v8_brain_tissues_used.tsv"
        ),
        "secondary_tissues": (
            args.outdir
            / "GTEx_v8_all_tissues_used.tsv"
        ),

        # Cross-analysis audit
        "overlap": (
            args.outdir
            / "GTEx_v8_brain_vs_all_tissues_strict_resilience_overlap.tsv"
        ),
        "summary": (
            args.outdir
            / "GTEx_v8_directional_eqtl_gene_sets_summary.tsv"
        ),
        "readme": (
            args.outdir
            / "README_GTEx_v8_directional_eqtl_gene_sets.txt"
        ),
    }

    # -------------------------------------------------------------------------
    # Write GENE2FUNC lists
    # -------------------------------------------------------------------------
    n_primary_foreground = write_gene_list(
        primary["foreground"],
        outputs["primary_foreground"],
        args.overwrite,
    )

    n_primary_background = write_gene_list(
        primary["background"],
        outputs["primary_background"],
        args.overwrite,
    )

    n_secondary_foreground = write_gene_list(
        secondary["foreground"],
        outputs["secondary_foreground"],
        args.overwrite,
    )

    n_secondary_background = write_gene_list(
        secondary["background"],
        outputs["secondary_background"],
        args.overwrite,
    )

    # -------------------------------------------------------------------------
    # Write detailed tables
    # -------------------------------------------------------------------------
    write_tsv(
        primary["genes"],
        outputs["primary_gene_summary"],
        args.overwrite,
    )

    write_tsv(
        secondary["genes"],
        outputs["secondary_gene_summary"],
        args.overwrite,
    )

    write_tsv(
        primary["units"],
        outputs["primary_units"],
        args.overwrite,
    )

    write_tsv(
        secondary["units"],
        outputs["secondary_units"],
        args.overwrite,
    )

    primary_tissues = pd.DataFrame(
        {
            "database": GTEX_DATABASE,
            "tissue": sorted(
                primary_rows[
                    "tissue"
                ]
                .dropna()
                .astype(str)
                .unique()
                .tolist()
            ),
            "analysis": "Primary_GTEx_v8_brain_only",
        }
    )

    secondary_tissues = pd.DataFrame(
        {
            "database": GTEX_DATABASE,
            "tissue": sorted(
                secondary_rows[
                    "tissue"
                ]
                .dropna()
                .astype(str)
                .unique()
                .tolist()
            ),
            "analysis": "Secondary_GTEx_v8_all_tissues",
        }
    )

    write_tsv(
        primary_tissues,
        outputs["primary_tissues"],
        args.overwrite,
    )

    write_tsv(
        secondary_tissues,
        outputs["secondary_tissues"],
        args.overwrite,
    )

    write_tsv(
        overlap,
        outputs["overlap"],
        args.overwrite,
    )

    write_tsv(
        summary,
        outputs["summary"],
        args.overwrite,
    )

    write_readme(
        outputs["readme"],
        args.overwrite,
        summary,
        primary,
        secondary,
    )

    # -------------------------------------------------------------------------
    # Console summary
    # -------------------------------------------------------------------------
    print(
        "\n============================================================",
        flush=True,
    )
    print(
        "GTEx v8-only directional eQTL gene-set construction complete",
        flush=True,
    )
    print(
        "============================================================",
        flush=True,
    )

    print(
        "\nPRIMARY: GTEx v8 brain-only",
        flush=True,
    )
    print(
        f"  Brain tissues: {primary_rows['tissue'].nunique():,}",
        flush=True,
    )
    print(
        f"  Matched background: {n_primary_background:,} genes",
        flush=True,
    )
    print(
        f"  Strict resilience foreground: {n_primary_foreground:,} genes",
        flush=True,
    )

    print(
        "\nSECONDARY: all GTEx v8 tissues",
        flush=True,
    )
    print(
        f"  GTEx v8 tissues: {secondary_rows['tissue'].nunique():,}",
        flush=True,
    )
    print(
        f"  Matched background: {n_secondary_background:,} genes",
        flush=True,
    )
    print(
        f"  Strict resilience foreground: {n_secondary_foreground:,} genes",
        flush=True,
    )

    primary_set = set(
        primary[
            "foreground"
        ]
    )

    secondary_set = set(
        secondary[
            "foreground"
        ]
    )

    print(
        f"\nStrict-gene overlap: {len(primary_set & secondary_set):,}",
        flush=True,
    )

    print(
        "\nFUMA GENE2FUNC files:",
        flush=True,
    )

    for key in [
        "primary_foreground",
        "primary_background",
        "secondary_foreground",
        "secondary_background",
    ]:
        print(
            f"  {outputs[key]}",
            flush=True,
        )

    print(
        "\nDetailed audit files:",
        flush=True,
    )

    for key in [
        "primary_gene_summary",
        "secondary_gene_summary",
        "primary_units",
        "secondary_units",
        "primary_tissues",
        "secondary_tissues",
        "overlap",
        "summary",
        "readme",
    ]:
        print(
            f"  {outputs[key]}",
            flush=True,
        )
    print('Stop...')

if __name__ == "__main__":
    main()