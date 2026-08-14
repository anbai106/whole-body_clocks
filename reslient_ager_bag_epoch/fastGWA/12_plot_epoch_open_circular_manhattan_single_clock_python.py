#!/usr/bin/env python3

"""
Create an open circular Manhattan plot for an EPOCH GWAS.

Inputs
------
1. fastGWA summary-statistics file:
   organ_pheno_normalized_residualized.fastGWA.zip

2. FUMA genomic-risk-locus file:
   GenomicRiskLoci.txt

3. FUMA GWAS Catalog annotation file:
   gwascatalog.txt

Outputs
-------
- <out_prefix>.png
- <out_prefix>.pdf
- <out_prefix>.svg
- <out_prefix>_lead_snps.tsv
- <out_prefix>_wordcloud_frequencies.tsv

The top lead SNP for each locus is taken from the `rsID` column of
GenomicRiskLoci.txt. This is preferable to taking the first entry in
LeadSNPs because the LeadSNPs field can contain multiple variants whose
order does not necessarily identify the locus-level top SNP.

Required python package
pip install pandas numpy matplotlib wordcloud

"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
import zipfile
from collections import OrderedDict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import patheffects
from matplotlib.colors import to_hex
from matplotlib.lines import Line2D

try:
    from wordcloud import WordCloud
except ImportError as exc:
    raise ImportError(
        "The Python package 'wordcloud' is required. Install it with:\n"
        "pip install wordcloud"
    ) from exc


###############################################################################
# Global configuration
###############################################################################

CHR_ORDER = [str(i) for i in range(1, 23)] + ["X"]
CHR_INDEX = {chrom: index for index, chrom in enumerate(CHR_ORDER)}

# Soft alternating chromosome colors.
CHR_COLORS = [
    "#64B5D9",
    "#F39A38",
    "#62C58A",
    "#F05B68",
    "#A98BC7",
    "#BA8F8F",
    "#E8A6C4",
    "#A6A9AD",
    "#C7C681",
    "#72AFCF",
    "#E6B976",
    "#80C995",
    "#C1A5D5",
    "#A98BC7",
    "#EAA9C5",
    "#AAB0B5",
    "#C1C44C",
    "#B39ACB",
    "#E798B5",
    "#83B6D5",
    "#9DCA9A",
    "#B9BEC3",
    "#92979C",
]

# Ordered keyword rules for automatically grouping GWAS Catalog traits.
# The first matching category is assigned.
CATEGORY_RULES = OrderedDict(
    [
        (
            "Brain structure and function",
            [
                r"\bbrain\b",
                r"\bcortical\b",
                r"\bcortex\b",
                r"\bhippocamp",
                r"\bwhite matter\b",
                r"\bgrey matter\b",
                r"\bgray matter\b",
                r"\bintracranial\b",
                r"\bneuroimaging\b",
                r"\bbrain volume\b",
                r"\bbrain morph",
                r"\bsubcortical\b",
                r"\bconnectiv",
                r"\bneuronal\b",
            ],
        ),
        (
            "Neurological diseases",
            [
                r"\balzheimer",
                r"\bdementia\b",
                r"\bparkinson",
                r"\bmultiple sclerosis\b",
                r"\bepilep",
                r"\bstroke\b",
                r"\bmigraine\b",
                r"\bneuropath",
                r"\bamyotrophic",
                r"\bals\b",
                r"\bhuntington",
                r"\bneurological",
                r"\bneurodegener",
            ],
        ),
        (
            "Depression and psychiatric traits",
            [
                r"\bdepress",
                r"\banxiety\b",
                r"\bschizophren",
                r"\bbipolar\b",
                r"\bpsychiatr",
                r"\bmental health\b",
                r"\bneurotic",
                r"\bwell-being\b",
                r"\bwellbeing\b",
                r"\bmood\b",
                r"\bpost-traumatic\b",
                r"\bptsd\b",
            ],
        ),
        (
            "Cognitive function and education",
            [
                r"\bcognitive\b",
                r"\bintelligence\b",
                r"\beducation",
                r"\bmemory\b",
                r"\breaction time\b",
                r"\bexecutive function\b",
                r"\bmathematical ability\b",
                r"\bverbal\b",
            ],
        ),
        (
            "Sleep and circadian traits",
            [
                r"\bsleep\b",
                r"\binsomnia\b",
                r"\bchronotype\b",
                r"\bcircadian\b",
                r"\bdaytime sleepiness\b",
            ],
        ),
        (
            "Body and muscle composition",
            [
                r"\bbody mass index\b",
                r"\bbmi\b",
                r"\bbody composition\b",
                r"\bbody fat\b",
                r"\bfat mass\b",
                r"\blean mass\b",
                r"\bmuscle\b",
                r"\bwaist\b",
                r"\bhip circumference\b",
                r"\bheight\b",
                r"\bweight\b",
                r"\bobesity\b",
                r"\bsarcopen",
                r"\bgrip strength\b",
            ],
        ),
        (
            "Blood cell traits",
            [
                r"\bplatelet\b",
                r"\bred blood cell\b",
                r"\bred cell\b",
                r"\bwhite blood cell\b",
                r"\bwhite cell\b",
                r"\bleukocyte\b",
                r"\bleukocyte\b",
                r"\bneutrophil\b",
                r"\bmonocyte\b",
                r"\blymphocyte\b",
                r"\beosinophil\b",
                r"\bbasophil\b",
                r"\bhematocrit\b",
                r"\bhaematocrit\b",
                r"\bhemoglobin\b",
                r"\bhaemoglobin\b",
                r"\bcorpuscular\b",
                r"\bblood-cell\b",
                r"\bblood cell\b",
            ],
        ),
        (
            "Immune cells and inflammation",
            [
                r"\bimmune\b",
                r"\binflamm",
                r"\bcytokine\b",
                r"\binterleukin\b",
                r"\bc-reactive protein\b",
                r"\bcrp\b",
                r"\bautoimmune\b",
                r"\ballergy\b",
                r"\basthma\b",
                r"\bimmunoglobulin\b",
            ],
        ),
        (
            "Cardiovascular function and disease",
            [
                r"\bcardiovascular\b",
                r"\bcoronary\b",
                r"\bheart\b",
                r"\bcardiac\b",
                r"\bmyocard",
                r"\batrial fibrillation\b",
                r"\baortic\b",
                r"\bartery\b",
                r"\barterial\b",
                r"\bvascular\b",
                r"\bventric",
                r"\bcarotid\b",
                r"\bheart rate\b",
                r"\becg\b",
                r"\belectrocard",
            ],
        ),
        (
            "Blood pressure",
            [
                r"\bblood pressure\b",
                r"\bhypertension\b",
                r"\bsystolic\b",
                r"\bdiastolic\b",
                r"\bpulse pressure\b",
            ],
        ),
        (
            "Respiratory function",
            [
                r"\blung\b",
                r"\bpulmonary\b",
                r"\brespiratory\b",
                r"\bforced expiratory\b",
                r"\bfev1\b",
                r"\bforced vital capacity\b",
                r"\bfvc\b",
                r"\bcopd\b",
            ],
        ),
        (
            "Glucose metabolism and diabetes",
            [
                r"\bdiabetes\b",
                r"\bglucose\b",
                r"\bglycemic\b",
                r"\bglycaemic\b",
                r"\binsulin\b",
                r"\bhba1c\b",
                r"\bhemoglobin a1c\b",
                r"\bhaemoglobin a1c\b",
            ],
        ),
        (
            "Cholesterol and lipid metabolism",
            [
                r"\bcholesterol\b",
                r"\btriglyceride\b",
                r"\blipid\b",
                r"\blipoprotein\b",
                r"\bhdl\b",
                r"\bldl\b",
                r"\bapolipoprotein\b",
                r"\bfatty acid\b",
                r"\bphospholipid\b",
            ],
        ),
        (
            "Cellular energy and metabolite traits",
            [
                r"\bmetabol",
                r"\bamino acid\b",
                r"\blactate\b",
                r"\bpyruvate\b",
                r"\bcitrate\b",
                r"\bketone\b",
                r"\benergy expenditure\b",
                r"\bmitochond",
                r"\boxidative\b",
            ],
        ),
        (
            "Liver function",
            [
                r"\bliver\b",
                r"\bhepatic\b",
                r"\balt\b",
                r"\bast\b",
                r"\balanine aminotransferase\b",
                r"\baspartate aminotransferase\b",
                r"\bgamma-glutamyl\b",
                r"\bbilirubin\b",
            ],
        ),
        (
            "Renal and urinary function",
            [
                r"\bkidney\b",
                r"\brenal\b",
                r"\burinary\b",
                r"\burine\b",
                r"\bcreatinine\b",
                r"\begfr\b",
                r"\bglomerular\b",
                r"\burate\b",
                r"\buric acid\b",
            ],
        ),
        (
            "Bone and joint traits",
            [
                r"\bbone\b",
                r"\bosteoporosis\b",
                r"\bosteoarthritis\b",
                r"\bfracture\b",
                r"\bjoint\b",
                r"\bmineral density\b",
                r"\bcalcaneal\b",
            ],
        ),
        (
            "Endocrine and hormone traits",
            [
                r"\bhormone\b",
                r"\bthyroid\b",
                r"\btestosterone\b",
                r"\bestrogen\b",
                r"\boestrogen\b",
                r"\bcortisol\b",
                r"\bsex hormone\b",
                r"\bendocrine\b",
                r"\bmenopause\b",
                r"\bmenarche\b",
            ],
        ),
        (
            "Cancer",
            [
                r"\bcancer\b",
                r"\bcarcinoma\b",
                r"\bneoplasm\b",
                r"\btumou?r\b",
                r"\bleukemia\b",
                r"\bleukaemia\b",
                r"\blymphoma\b",
                r"\bmelanoma\b",
                r"\bglioma\b",
            ],
        ),
        (
            "Digestive function and disease",
            [
                r"\bdigestive\b",
                r"\bgastro",
                r"\bintestinal\b",
                r"\bbowel\b",
                r"\bcrohn",
                r"\bcolitis\b",
                r"\bceliac\b",
                r"\bcoeliac\b",
                r"\bpancrea",
                r"\bgallbladder\b",
            ],
        ),
        (
            "Reproductive traits",
            [
                r"\breproductive\b",
                r"\bfertility\b",
                r"\bpregnan",
                r"\bbirth weight\b",
                r"\bgestational\b",
                r"\bovarian\b",
                r"\bprostate\b",
                r"\bmenstrual\b",
                r"\bsperm\b",
            ],
        ),
        (
            "Skin and dermatological traits",
            [
                r"\bskin\b",
                r"\bdermat",
                r"\bpsoriasis\b",
                r"\beczema\b",
                r"\bpigmentation\b",
                r"\bhair color\b",
                r"\bhair colour\b",
                r"\bvitiligo\b",
            ],
        ),
        (
            "Eye structure and function",
            [
                r"\beye\b",
                r"\bretina\b",
                r"\bretinal\b",
                r"\bglaucoma\b",
                r"\bcataract\b",
                r"\bvision\b",
                r"\bvisual acuity\b",
                r"\bmyopia\b",
                r"\bcornea\b",
            ],
        ),
        (
            "Behavior and lifestyle factors",
            [
                r"\bsmoking\b",
                r"\balcohol\b",
                r"\bdiet\b",
                r"\bphysical activity\b",
                r"\bsedentary\b",
                r"\bcoffee\b",
                r"\bcaffeine\b",
                r"\bbehavior\b",
                r"\bbehaviour\b",
            ],
        ),
        (
            "Medication response",
            [
                r"\bdrug response\b",
                r"\bmedication\b",
                r"\bpharmacogen",
                r"\btreatment response\b",
            ],
        ),
        (
            "Longevity and aging",
            [
                r"\baging\b",
                r"\bageing\b",
                r"\blongevity\b",
                r"\blifespan\b",
                r"\bparental age\b",
                r"\bbiological age\b",
            ],
        ),
    ]
)


###############################################################################
# General helpers
###############################################################################


def normalize_chr(value: object) -> Optional[str]:
    """Normalize chromosome labels to 1-22 or X."""
    if pd.isna(value):
        return None

    text = str(value).strip()
    text = re.sub(r"^chr", "", text, flags=re.IGNORECASE)

    if text.endswith(".0"):
        text = text[:-2]

    if text == "23":
        text = "X"

    text = text.upper()

    if text in CHR_ORDER:
        return text

    return None


def detect_column(
    dataframe: pd.DataFrame,
    candidates: Sequence[str],
    required: bool = True,
) -> Optional[str]:
    """Return the first matching column, case-insensitively."""
    normalized = {str(column).strip().lower(): column for column in dataframe.columns}

    for candidate in candidates:
        match = normalized.get(candidate.lower())
        if match is not None:
            return match

    if required:
        raise ValueError(
            "Could not identify a required column. "
            f"Accepted names were: {', '.join(candidates)}. "
            f"Observed columns were: {', '.join(map(str, dataframe.columns))}"
        )

    return None


def read_table_robust(path: Path) -> pd.DataFrame:
    """
    Read a whitespace- or tab-delimited table, including ZIP/GZIP files.

    Notes
    -----
    pandas does not support ``low_memory`` when ``engine="python"``.
    Therefore, ``low_memory=False`` is used only with the default C engine.
    For ZIP files, the first non-directory member is read.
    """
    suffixes = [suffix.lower() for suffix in path.suffixes]

    if ".zip" in suffixes:
        with zipfile.ZipFile(path, "r") as archive:
            members = [
                member
                for member in archive.namelist()
                if not member.endswith("/") and not member.startswith("__MACOSX")
            ]

            if not members:
                raise ValueError(f"No readable file was found inside {path}")

            if len(members) > 1:
                print(
                    f"Warning: {path} contains multiple files. "
                    f"Reading the first member: {members[0]}",
                    file=sys.stderr,
                )

            with archive.open(members[0]) as handle:
                # The Python parser supports regex whitespace separators,
                # but it does not accept low_memory.
                return pd.read_csv(
                    handle,
                    sep=r"\s+",
                    engine="python",
                )

    compression = "gzip" if ".gz" in suffixes else "infer"

    try:
        # First try a true tab-delimited file with the default C engine.
        return pd.read_csv(
            path,
            sep="\t",
            compression=compression,
            low_memory=False,
        )
    except (pd.errors.ParserError, UnicodeDecodeError, ValueError):
        # Fall back to arbitrary whitespace. Do not pass low_memory here,
        # because it is unsupported by the Python parsing engine.
        return pd.read_csv(
            path,
            sep=r"\s+",
            engine="python",
            compression=compression,
        )


###############################################################################
# Input readers
###############################################################################


def read_fastgwa(path: Path) -> pd.DataFrame:
    """Read and standardize fastGWA summary statistics."""
    dataframe = read_table_robust(path)

    chromosome_column = detect_column(
        dataframe,
        ["CHR", "CHROM", "CHROMOSOME", "#CHROM"],
    )
    position_column = detect_column(
        dataframe,
        ["POS", "BP", "POSITION", "BASE_PAIR_LOCATION"],
    )
    pvalue_column = detect_column(
        dataframe,
        ["P", "PVAL", "PVALUE", "P_VALUE"],
    )
    snp_column = detect_column(
        dataframe,
        ["SNP", "RSID", "ID", "MARKERNAME"],
        required=False,
    )

    selected = dataframe[
        [chromosome_column, position_column, pvalue_column]
        + ([snp_column] if snp_column is not None else [])
    ].copy()

    rename_map = {
        chromosome_column: "CHR",
        position_column: "POS",
        pvalue_column: "P",
    }

    if snp_column is not None:
        rename_map[snp_column] = "SNP"

    selected = selected.rename(columns=rename_map)

    selected["CHR"] = selected["CHR"].map(normalize_chr)
    selected["POS"] = pd.to_numeric(selected["POS"], errors="coerce")
    selected["P"] = pd.to_numeric(selected["P"], errors="coerce")

    selected = selected.dropna(subset=["CHR", "POS", "P"])
    selected = selected[
        selected["CHR"].isin(CHR_ORDER)
        & np.isfinite(selected["P"])
        & (selected["P"] > 0)
        & (selected["P"] <= 1)
        & (selected["POS"] > 0)
    ].copy()

    selected["POS"] = selected["POS"].astype(np.int64)

    if "SNP" not in selected.columns:
        selected["SNP"] = (
            selected["CHR"].astype(str)
            + ":"
            + selected["POS"].astype(str)
        )
    else:
        selected["SNP"] = selected["SNP"].astype(str).str.strip()

    selected["CHR_INDEX"] = selected["CHR"].map(CHR_INDEX)
    selected = selected.sort_values(["CHR_INDEX", "POS", "P"]).reset_index(drop=True)
    selected = selected.drop(columns="CHR_INDEX")

    return selected


def read_genomic_risk_loci(path: Path) -> pd.DataFrame:
    """
    Read FUMA GenomicRiskLoci.txt.

    The `rsID` field is used as the locus-level top lead SNP label.
    """
    loci = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    loci.columns = [str(column).strip() for column in loci.columns]

    required_columns = [
        "GenomicLocus",
        "rsID",
        "chr",
        "pos",
        "p",
    ]

    missing = [column for column in required_columns if column not in loci.columns]
    if missing:
        raise ValueError(
            f"Missing required columns in {path}: {', '.join(missing)}"
        )

    loci["GenomicLocus"] = pd.to_numeric(
        loci["GenomicLocus"],
        errors="coerce",
    )
    loci["CHR"] = loci["chr"].map(normalize_chr)
    loci["POS"] = pd.to_numeric(loci["pos"], errors="coerce")
    loci["P"] = pd.to_numeric(loci["p"], errors="coerce")
    loci["LabelSNP"] = loci["rsID"].astype(str).str.strip()

    loci = loci.dropna(
        subset=["GenomicLocus", "CHR", "POS", "P", "LabelSNP"]
    ).copy()

    loci = loci[
        loci["CHR"].isin(CHR_ORDER)
        & np.isfinite(loci["P"])
        & (loci["P"] > 0)
        & (loci["P"] <= 1)
        & (loci["POS"] > 0)
        & (loci["LabelSNP"] != "")
        & (loci["LabelSNP"].str.lower() != "nan")
    ].copy()

    loci["GenomicLocus"] = loci["GenomicLocus"].astype(int)
    loci["POS"] = loci["POS"].astype(np.int64)

    # One row per genomic locus; select the smallest locus-level P value if needed.
    loci = (
        loci.sort_values(["GenomicLocus", "P"])
        .drop_duplicates("GenomicLocus", keep="first")
        .reset_index(drop=True)
    )

    return loci[
        [
            "GenomicLocus",
            "LabelSNP",
            "CHR",
            "POS",
            "P",
        ]
    ]


def read_gwas_catalog(path: Path) -> pd.DataFrame:
    """Read FUMA gwascatalog.txt and standardize key fields."""
    catalog = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        quoting=csv.QUOTE_NONE,
        on_bad_lines="warn",
        low_memory=False,
    )

    catalog.columns = [str(column).strip() for column in catalog.columns]

    required_columns = ["GenomicLocus", "Trait"]
    missing = [column for column in required_columns if column not in catalog.columns]

    if missing:
        raise ValueError(
            f"Missing required columns in {path}: {', '.join(missing)}"
        )

    catalog["GenomicLocus"] = pd.to_numeric(
        catalog["GenomicLocus"],
        errors="coerce",
    )
    catalog["Trait"] = (
        catalog["Trait"]
        .fillna("")
        .astype(str)
        .str.replace(r"\s+", " ", regex=True)
        .str.strip()
    )

    if "PMID" not in catalog.columns:
        catalog["PMID"] = ""

    if "IndSigSNP" not in catalog.columns:
        catalog["IndSigSNP"] = ""

    catalog["PMID"] = catalog["PMID"].fillna("").astype(str).str.strip()
    catalog["IndSigSNP"] = (
        catalog["IndSigSNP"].fillna("").astype(str).str.strip()
    )

    catalog = catalog.dropna(subset=["GenomicLocus"]).copy()
    catalog = catalog[
        (catalog["Trait"] != "")
        & (catalog["Trait"].str.lower() != "nan")
    ].copy()

    catalog["GenomicLocus"] = catalog["GenomicLocus"].astype(int)

    return catalog


###############################################################################
# Lead-SNP matching
###############################################################################


def match_lead_snps_to_gwas(
    gwas: pd.DataFrame,
    loci: pd.DataFrame,
) -> pd.DataFrame:
    """
    Match FUMA locus-level top SNPs to fastGWA results.

    Matching priority:
    1. Exact SNP ID;
    2. Exact chromosome and position;
    3. Use the chromosome, position, and P value reported by FUMA.
    """
    gwas_by_snp = (
        gwas.sort_values("P")
        .drop_duplicates("SNP", keep="first")
        .set_index("SNP")
    )

    gwas_by_position = (
        gwas.sort_values("P")
        .drop_duplicates(["CHR", "POS"], keep="first")
        .set_index(["CHR", "POS"])
    )

    output_rows = []

    for row in loci.itertuples(index=False):
        match_source = "FUMA locus file"
        matched_snp = row.LabelSNP
        matched_chr = row.CHR
        matched_pos = int(row.POS)
        matched_p = float(row.P)

        if row.LabelSNP in gwas_by_snp.index:
            gwas_row = gwas_by_snp.loc[row.LabelSNP]

            if isinstance(gwas_row, pd.DataFrame):
                gwas_row = gwas_row.iloc[0]

            matched_snp = row.LabelSNP
            matched_chr = str(gwas_row["CHR"])
            matched_pos = int(gwas_row["POS"])
            matched_p = float(gwas_row["P"])
            match_source = "fastGWA SNP ID"

        elif (row.CHR, int(row.POS)) in gwas_by_position.index:
            gwas_row = gwas_by_position.loc[(row.CHR, int(row.POS))]

            if isinstance(gwas_row, pd.DataFrame):
                gwas_row = gwas_row.iloc[0]

            matched_snp = row.LabelSNP
            matched_chr = str(gwas_row["CHR"])
            matched_pos = int(gwas_row["POS"])
            matched_p = float(gwas_row["P"])
            match_source = "fastGWA chromosome-position"

        output_rows.append(
            {
                "GenomicLocus": int(row.GenomicLocus),
                "LabelSNP": matched_snp,
                "CHR": matched_chr,
                "POS": matched_pos,
                "P": matched_p,
                "MatchSource": match_source,
            }
        )

    matched = pd.DataFrame(output_rows)
    matched["neglog10P"] = -np.log10(matched["P"])
    matched["CHR_INDEX"] = matched["CHR"].map(CHR_INDEX)

    matched = matched.sort_values(
        ["CHR_INDEX", "POS", "P"]
    ).drop(columns="CHR_INDEX")

    return matched.reset_index(drop=True)


###############################################################################
# GWAS Catalog word-cloud processing
###############################################################################


def normalize_trait_text(trait: str) -> str:
    """Normalize a GWAS Catalog trait string."""
    text = re.sub(r"\s+", " ", str(trait)).strip()
    text = re.sub(r"\s*\(.*?adjusted.*?\)\s*", " ", text, flags=re.IGNORECASE)
    text = re.sub(r"\s+", " ", text).strip(" ;,")
    return text


def assign_trait_category(trait: str) -> str:
    """Assign one broad phenotype category using ordered keyword rules."""
    lower_trait = normalize_trait_text(trait).lower()

    for category, patterns in CATEGORY_RULES.items():
        for pattern in patterns:
            if re.search(pattern, lower_trait, flags=re.IGNORECASE):
                return category

    return "Other phenotypic traits"


def build_wordcloud_frequencies(
    catalog: pd.DataFrame,
    valid_loci: Iterable[int],
    mode: str = "category",
    count_method: str = "locus",
    minimum_frequency: int = 1,
    maximum_terms: int = 80,
    include_other: bool = False,
) -> Tuple[Dict[str, float], pd.DataFrame]:
    """
    Build word-cloud frequencies from FUMA GWAS Catalog annotations.

    Parameters
    ----------
    mode
        "category": group traits into broad phenotype categories.
        "trait": use normalized GWAS Catalog trait names directly.

    count_method
        "locus":
            Count each term once per GenomicLocus. Recommended because it
            prevents a locus with many linked catalog SNPs from dominating.

        "association":
            Count unique GenomicLocus + PMID + Trait combinations.

        "row":
            Count every row in gwascatalog.txt.
    """
    valid_loci_set = {int(value) for value in valid_loci}

    subset = catalog[
        catalog["GenomicLocus"].isin(valid_loci_set)
    ].copy()

    if subset.empty:
        raise ValueError(
            "No GWAS Catalog annotations matched the genomic loci in "
            "GenomicRiskLoci.txt."
        )

    subset["TraitNormalized"] = subset["Trait"].map(normalize_trait_text)

    if mode == "category":
        subset["WordCloudTerm"] = subset["TraitNormalized"].map(
            assign_trait_category
        )
    elif mode == "trait":
        subset["WordCloudTerm"] = subset["TraitNormalized"]
    else:
        raise ValueError("wordcloud mode must be 'category' or 'trait'")

    subset = subset[
        (subset["WordCloudTerm"] != "")
        & subset["WordCloudTerm"].notna()
    ].copy()

    if not include_other:
        subset = subset[
            subset["WordCloudTerm"] != "Other phenotypic traits"
        ].copy()

    if count_method == "locus":
        counted = subset.drop_duplicates(
            ["GenomicLocus", "WordCloudTerm"]
        )
    elif count_method == "association":
        counted = subset.drop_duplicates(
            ["GenomicLocus", "PMID", "TraitNormalized", "WordCloudTerm"]
        )
    elif count_method == "row":
        counted = subset
    else:
        raise ValueError(
            "count_method must be 'locus', 'association', or 'row'"
        )

    frequency_table = (
        counted.groupby("WordCloudTerm", as_index=False)
        .size()
        .rename(columns={"size": "Frequency"})
        .sort_values(
            ["Frequency", "WordCloudTerm"],
            ascending=[False, True],
        )
    )

    frequency_table = frequency_table[
        frequency_table["Frequency"] >= minimum_frequency
    ].head(maximum_terms)

    if frequency_table.empty:
        raise ValueError(
            "No word-cloud terms remained after filtering. "
            "Try --wordcloud-min-frequency 1 or --include-other."
        )

    frequencies = dict(
        zip(
            frequency_table["WordCloudTerm"],
            frequency_table["Frequency"].astype(float),
        )
    )

    return frequencies, frequency_table.reset_index(drop=True)


###############################################################################
# Circular genomic layout
###############################################################################


def build_genome_layout(
    gwas: pd.DataFrame,
    gap_fraction: float,
) -> Tuple[
    List[str],
    Dict[str, float],
    Dict[str, float],
    Dict[str, float],
    float,
]:
    """
    Construct cumulative genomic coordinates using observed chromosome maxima.
    """
    chromosome_lengths = (
        gwas.groupby("CHR")["POS"]
        .max()
        .reindex(CHR_ORDER)
        .dropna()
        .astype(float)
    )

    present_chromosomes = chromosome_lengths.index.tolist()

    if not present_chromosomes:
        raise ValueError("No autosomal or chromosome-X variants were found.")

    total_bp = float(chromosome_lengths.sum())
    gap_bp = total_bp * gap_fraction

    starts: Dict[str, float] = {}
    ends: Dict[str, float] = {}
    midpoints: Dict[str, float] = {}

    cursor = 0.0

    for chromosome in present_chromosomes:
        starts[chromosome] = cursor
        ends[chromosome] = cursor + float(chromosome_lengths[chromosome])
        midpoints[chromosome] = (
            starts[chromosome] + ends[chromosome]
        ) / 2.0
        cursor = ends[chromosome] + gap_bp

    genome_span = cursor - gap_bp

    return (
        present_chromosomes,
        starts,
        ends,
        midpoints,
        genome_span,
    )


def add_angles(
    dataframe: pd.DataFrame,
    starts: Dict[str, float],
    genome_span: float,
    theta_start_degrees: float,
    sweep_degrees: float,
) -> pd.DataFrame:
    """Map chromosome-position coordinates to polar angles."""
    output = dataframe.copy()

    output["CUM_POS"] = [
        starts[chromosome] + position
        for chromosome, position in zip(output["CHR"], output["POS"])
    ]

    output["THETA_DEG"] = (
        theta_start_degrees
        + sweep_degrees * output["CUM_POS"] / genome_span
    )
    output["THETA"] = np.deg2rad(output["THETA_DEG"])

    return output


def choose_background_variants(
    gwas: pd.DataFrame,
    lead_snps: pd.DataFrame,
    keep_p_threshold: float,
    maximum_background_points: int,
    seed: int,
) -> pd.DataFrame:
    """
    Keep all stronger associations and a reproducible random background sample.
    """
    lead_positions = set(
        zip(
            lead_snps["CHR"].astype(str),
            lead_snps["POS"].astype(int),
        )
    )

    is_lead = np.fromiter(
        (
            (str(chromosome), int(position)) in lead_positions
            for chromosome, position in zip(gwas["CHR"], gwas["POS"])
        ),
        dtype=bool,
        count=len(gwas),
    )

    must_keep = (gwas["P"].to_numpy() <= keep_p_threshold) | is_lead

    retained = gwas.loc[must_keep].copy()
    background = gwas.loc[~must_keep].copy()

    if len(background) > maximum_background_points:
        background = background.sample(
            n=maximum_background_points,
            random_state=seed,
            replace=False,
        )

    return pd.concat(
        [retained, background],
        ignore_index=True,
    )


###############################################################################
# Word-cloud drawing
###############################################################################


def create_wordcloud_colors(
    frequencies: Dict[str, float],
) -> Dict[str, str]:
    """Assign reproducible colors to word-cloud terms."""
    cmap = plt.get_cmap("tab20")
    sorted_terms = sorted(
        frequencies,
        key=lambda term: (-frequencies[term], term),
    )

    return {
        term: to_hex(cmap(index % cmap.N))
        for index, term in enumerate(sorted_terms)
    }


def draw_wordcloud(
    figure: plt.Figure,
    frequencies: Dict[str, float],
    title: str,
    bounding_box: Tuple[float, float, float, float],
    random_state: int,
) -> None:
    """Draw a GWAS Catalog word cloud in a rectangular inset."""
    inset = figure.add_axes(bounding_box)

    term_colors = create_wordcloud_colors(frequencies)

    def color_function(
        word: str,
        font_size: int,
        position: Tuple[int, int],
        orientation: Optional[int],
        random_state: Optional[np.random.RandomState] = None,
        **kwargs: object,
    ) -> str:
        del font_size, position, orientation, random_state, kwargs
        return term_colors.get(word, "#555555")

    word_cloud = WordCloud(
        width=1600,
        height=1050,
        background_color="white",
        max_words=len(frequencies),
        prefer_horizontal=0.92,
        relative_scaling=0.45,
        collocations=False,
        random_state=random_state,
        margin=4,
        min_font_size=10,
        max_font_size=115,
    ).generate_from_frequencies(frequencies)

    word_cloud = word_cloud.recolor(color_func=color_function)

    inset.imshow(word_cloud, interpolation="bilinear")
    inset.set_xticks([])
    inset.set_yticks([])
    inset.set_title(title, fontsize=11, pad=7)

    for spine in inset.spines.values():
        spine.set_visible(True)
        spine.set_color("#B6B6B6")
        spine.set_linewidth(1.2)


###############################################################################
# Label placement
###############################################################################


def spread_label_radii(
    leads: pd.DataFrame,
    base_offset: float,
    step: float,
    angular_window_degrees: float,
) -> pd.DataFrame:
    """
    Assign staggered label radii to nearby lead SNPs.

    This does not fully solve all possible label collisions, but substantially
    improves dense loci without requiring a separate label-repulsion package.
    """
    output = leads.sort_values("THETA_DEG").copy()
    levels: List[int] = []
    recent_angles: List[Tuple[float, int]] = []

    for theta in output["THETA_DEG"]:
        recent_angles = [
            (angle, level)
            for angle, level in recent_angles
            if theta - angle < angular_window_degrees
        ]

        occupied_levels = {level for _, level in recent_angles}
        level = 0

        while level in occupied_levels:
            level += 1

        levels.append(level)
        recent_angles.append((theta, level))

    output["LABEL_LEVEL"] = levels
    output["LABEL_OFFSET"] = (
        base_offset + output["LABEL_LEVEL"] * step
    )

    return output


###############################################################################
# Plotting
###############################################################################


def plot_open_circular_manhattan(
    gwas: pd.DataFrame,
    lead_snps: pd.DataFrame,
    wordcloud_frequencies: Dict[str, float],
    output_prefix: Path,
    title: str,
    wordcloud_title: str,
    p_suggestive: float,
    p_genome_wide: float,
    p_cap: Optional[float],
    theta_start_degrees: float,
    sweep_degrees: float,
    maximum_background_points: int,
    seed: int,
    figure_size: float,
    point_size: float,
) -> None:
    """Create the complete open circular Manhattan figure."""
    if sweep_degrees <= 0 or sweep_degrees >= 360:
        raise ValueError("sweep_degrees must be greater than 0 and below 360.")

    selected_gwas = choose_background_variants(
        gwas=gwas,
        lead_snps=lead_snps,
        keep_p_threshold=p_suggestive,
        maximum_background_points=maximum_background_points,
        seed=seed,
    )

    (
        present_chromosomes,
        chromosome_starts,
        chromosome_ends,
        chromosome_midpoints,
        genome_span,
    ) = build_genome_layout(
        gwas=gwas,
        gap_fraction=0.004,
    )

    selected_gwas = add_angles(
        selected_gwas,
        starts=chromosome_starts,
        genome_span=genome_span,
        theta_start_degrees=theta_start_degrees,
        sweep_degrees=sweep_degrees,
    )

    lead_snps = lead_snps[
        lead_snps["CHR"].isin(present_chromosomes)
    ].copy()

    lead_snps = add_angles(
        lead_snps,
        starts=chromosome_starts,
        genome_span=genome_span,
        theta_start_degrees=theta_start_degrees,
        sweep_degrees=sweep_degrees,
    )

    selected_gwas["NEG_LOG10_P"] = -np.log10(selected_gwas["P"])
    lead_snps["NEG_LOG10_P"] = -np.log10(lead_snps["P"])

    observed_maximum = max(
        float(selected_gwas["NEG_LOG10_P"].max()),
        float(lead_snps["NEG_LOG10_P"].max()),
        -math.log10(p_genome_wide),
    )

    if p_cap is None:
        radial_cap = max(
            10.0,
            math.ceil(observed_maximum + 1.0),
        )
    else:
        radial_cap = float(p_cap)

    radial_inner = 0.47
    radial_outer = 1.04
    chromosome_radius = 1.19
    chromosome_label_radius = 1.285

    def pvalue_height_to_radius(values):
        clipped = np.minimum(values, radial_cap)
        return radial_inner + (
            np.asarray(clipped) / radial_cap
        ) * (radial_outer - radial_inner)

    selected_gwas["RADIUS"] = pvalue_height_to_radius(
        selected_gwas["NEG_LOG10_P"]
    )
    lead_snps["RADIUS"] = pvalue_height_to_radius(
        lead_snps["NEG_LOG10_P"]
    )

    suggestive_radius = float(
        pvalue_height_to_radius(-math.log10(p_suggestive))
    )
    genome_wide_radius = float(
        pvalue_height_to_radius(-math.log10(p_genome_wide))
    )

    chromosome_colors = {
        chromosome: CHR_COLORS[index % len(CHR_COLORS)]
        for index, chromosome in enumerate(present_chromosomes)
    }

    figure = plt.figure(
        figsize=(figure_size, figure_size),
        facecolor="white",
    )

    axis = figure.add_axes(
        [0.035, 0.035, 0.93, 0.93],
        projection="polar",
    )

    axis.set_theta_zero_location("N")
    axis.set_theta_direction(-1)
    axis.set_ylim(0, chromosome_label_radius + 0.08)
    axis.set_xticks([])
    axis.set_yticks([])
    axis.grid(False)
    axis.spines["polar"].set_visible(False)
    axis.set_facecolor("white")

    # Threshold arcs.
    theta_arc = np.linspace(
        np.deg2rad(theta_start_degrees),
        np.deg2rad(theta_start_degrees + sweep_degrees),
        1800,
    )

    axis.plot(
        theta_arc,
        np.full_like(theta_arc, suggestive_radius),
        color="#075EB6",
        linewidth=2.5,
        linestyle=(0, (7, 5)),
        zorder=3,
    )

    axis.plot(
        theta_arc,
        np.full_like(theta_arc, genome_wide_radius),
        color="#E31A1C",
        linewidth=2.5,
        linestyle=(0, (7, 5)),
        zorder=3,
    )

    # SNP points by chromosome.
    for chromosome in present_chromosomes:
        chromosome_data = selected_gwas[
            selected_gwas["CHR"] == chromosome
        ]

        if chromosome_data.empty:
            continue

        axis.scatter(
            chromosome_data["THETA"],
            chromosome_data["RADIUS"],
            s=point_size,
            facecolors="none",
            edgecolors=chromosome_colors[chromosome],
            linewidths=0.55,
            alpha=0.64,
            rasterized=True,
            zorder=2,
        )

    # Radial stems for suggestive variants.
    suggestive_variants = selected_gwas[
        selected_gwas["P"] <= p_suggestive
    ]

    for chromosome in present_chromosomes:
        chromosome_data = suggestive_variants[
            suggestive_variants["CHR"] == chromosome
        ]

        if chromosome_data.empty:
            continue

        color = chromosome_colors[chromosome]

        for row in chromosome_data.itertuples(index=False):
            axis.plot(
                [row.THETA, row.THETA],
                [radial_inner, row.RADIUS],
                color=color,
                linewidth=0.65,
                alpha=0.26,
                zorder=1,
            )

    # Emphasized locus-level top lead SNPs.
    for row in lead_snps.itertuples(index=False):
        color = chromosome_colors[row.CHR]

        axis.plot(
            [row.THETA, row.THETA],
            [radial_inner, row.RADIUS],
            color=color,
            linewidth=1.8,
            alpha=0.9,
            zorder=5,
        )

        axis.scatter(
            [row.THETA],
            [row.RADIUS],
            s=46,
            facecolors="white",
            edgecolors=color,
            linewidths=1.8,
            zorder=6,
        )

    # Outer chromosome arcs and labels.
    for chromosome in present_chromosomes:
        chromosome_theta_start = np.deg2rad(
            theta_start_degrees
            + sweep_degrees
            * chromosome_starts[chromosome]
            / genome_span
        )
        chromosome_theta_end = np.deg2rad(
            theta_start_degrees
            + sweep_degrees
            * chromosome_ends[chromosome]
            / genome_span
        )
        chromosome_theta_midpoint = np.deg2rad(
            theta_start_degrees
            + sweep_degrees
            * chromosome_midpoints[chromosome]
            / genome_span
        )

        chromosome_arc = np.linspace(
            chromosome_theta_start,
            chromosome_theta_end,
            200,
        )

        axis.plot(
            chromosome_arc,
            np.full_like(chromosome_arc, chromosome_radius),
            color="#A6A9AD",
            linewidth=8.5,
            solid_capstyle="butt",
            zorder=4,
        )

        axis.text(
            chromosome_theta_midpoint,
            chromosome_label_radius,
            f"Chr{chromosome}",
            ha="center",
            va="center",
            fontsize=12.5,
            fontweight="bold",
            color="black",
            zorder=10,
        )

    # Stagger nearby lead-SNP labels.
    labeled_leads = spread_label_radii(
        lead_snps,
        base_offset=0.035,
        step=0.033,
        angular_window_degrees=3.2,
    )

    for row in labeled_leads.itertuples(index=False):
        degrees = row.THETA_DEG % 360.0

        # Place text on the side away from the center.
        if 0 <= degrees < 180:
            horizontal_alignment = "left"
            text_offset = 7
        else:
            horizontal_alignment = "right"
            text_offset = -7

        label_radius = min(
            row.RADIUS + row.LABEL_OFFSET,
            chromosome_radius - 0.035,
        )

        annotation = axis.annotate(
            row.LabelSNP,
            xy=(row.THETA, row.RADIUS),
            xytext=(row.THETA, label_radius),
            textcoords="data",
            ha=horizontal_alignment,
            va="center",
            fontsize=7.6,
            fontweight="bold",
            color="black",
            arrowprops={
                "arrowstyle": "-",
                "color": "#4C4C4C",
                "linewidth": 0.45,
                "alpha": 0.65,
            },
            annotation_clip=False,
            zorder=9,
        )

        annotation.set_path_effects(
            [
                patheffects.withStroke(
                    linewidth=2.3,
                    foreground="white",
                )
            ]
        )

        # Small horizontal adjustment in display coordinates.
        annotation.set_position(
            (
                annotation.get_position()[0],
                annotation.get_position()[1],
            )
        )

    # Panel label and title.
    figure.text(
        0.047,
        0.955,
        "a",
        fontsize=31,
        fontweight="bold",
        ha="left",
        va="top",
    )

    figure.text(
        0.105,
        0.952,
        title,
        fontsize=16,
        fontweight="bold",
        ha="left",
        va="top",
    )

    # P-value scale label in the opening.
    scale_text = (
        r"$-\log_{10}(P)$"
        f"\nmaximum shown = {radial_cap:g}"
    )

    figure.text(
        0.485,
        0.835,
        scale_text,
        fontsize=10.5,
        ha="center",
        va="top",
    )

    # Legend.
    legend_handles = [
        Line2D(
            [0],
            [0],
            color="#E31A1C",
            linewidth=2.5,
            linestyle=(0, (7, 5)),
            label=r"Genome-wide: $P<5\times10^{-8}$",
        ),
        Line2D(
            [0],
            [0],
            color="#075EB6",
            linewidth=2.5,
            linestyle=(0, (7, 5)),
            label=rf"Suggestive: $P<{p_suggestive:g}$",
        ),
    ]

    figure.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.500, 0.765),
        frameon=False,
        fontsize=9,
        ncol=1,
        handlelength=3.2,
    )

    # GWAS Catalog word cloud.
    draw_wordcloud(
        figure=figure,
        frequencies=wordcloud_frequencies,
        title=wordcloud_title,
        bounding_box=(0.548, 0.548, 0.405, 0.355),
        random_state=seed,
    )

    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    output_files = {
        "PNG": output_prefix.with_suffix(".png"),
        "PDF": output_prefix.with_suffix(".pdf"),
        "SVG": output_prefix.with_suffix(".svg"),
    }

    figure.savefig(
        output_files["PNG"],
        dpi=500,
        bbox_inches="tight",
        facecolor="white",
    )
    figure.savefig(
        output_files["PDF"],
        dpi=500,
        bbox_inches="tight",
        facecolor="white",
    )
    figure.savefig(
        output_files["SVG"],
        dpi=500,
        bbox_inches="tight",
        facecolor="white",
    )

    plt.close(figure)

    for label, path in output_files.items():
        print(f"Wrote {label}: {path}")


###############################################################################
# Command-line interface
###############################################################################


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create an open circular Manhattan plot from fastGWA results, "
            "FUMA genomic loci, and FUMA GWAS Catalog annotations."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    default_clock = "Brain_proteomics_mortality_clock"
    default_title = "Brain proteomics mortality EPOCH"

    default_root = Path(
        "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/"
        "mortality_clock"
    )

    parser.add_argument(
        "--gwas",
        type=Path,
        default=(
                default_root
                / "fastGWA"
                / "output"
                / default_clock
                / "organ_pheno_normalized_residualized.fastGWA.zip"
        ),
        help="fastGWA summary-statistics file.",
    )

    parser.add_argument(
        "--loci",
        type=Path,
        default=(
                default_root
                / "fuma"
                / default_clock
                / "GenomicRiskLoci.txt"
        ),
        help="FUMA GenomicRiskLoci.txt file.",
    )

    parser.add_argument(
        "--gwascatalog",
        type=Path,
        default=(
                default_root
                / "fuma"
                / default_clock
                / "gwascatalog.txt"
        ),
        help="FUMA gwascatalog.txt file.",
    )

    parser.add_argument(
        "--out-prefix",
        type=Path,
        default=(
                default_root
                / "fuma"
                / default_clock
                / "EPOCH_open_circular"
        ),
        help="Output prefix without an extension.",
    )

    parser.add_argument(
        "--title",
        default=default_title,
        help="Figure title.",
    )

    parser.add_argument(
        "--wordcloud-mode",
        choices=["category", "trait"],
        default="category",
        help=(
            "'category' groups GWAS Catalog traits into broad phenotype "
            "categories; 'trait' plots the original normalized Trait names."
        ),
    )
    parser.add_argument(
        "--wordcloud-count",
        choices=["locus", "association", "row"],
        default="locus",
        help=(
            "How word-cloud frequencies are counted. 'locus' counts a term "
            "once per genomic locus and is recommended."
        ),
    )
    parser.add_argument(
        "--wordcloud-min-frequency",
        type=int,
        default=1,
        help="Minimum term frequency retained in the word cloud.",
    )
    parser.add_argument(
        "--wordcloud-max-terms",
        type=int,
        default=60,
        help="Maximum number of word-cloud terms.",
    )
    parser.add_argument(
        "--include-other",
        action="store_true",
        help="Include the residual 'Other phenotypic traits' category.",
    )

    parser.add_argument(
        "--p-suggestive",
        type=float,
        default=1e-5,
        help="Suggestive-association threshold.",
    )
    parser.add_argument(
        "--p-genome-wide",
        type=float,
        default=5e-8,
        help="Genome-wide-significance threshold.",
    )
    parser.add_argument(
        "--p-cap",
        type=float,
        default=None,
        help=(
            "Optional maximum displayed -log10(P). When omitted, the script "
            "chooses a value from the observed data."
        ),
    )
    parser.add_argument(
        "--theta-start",
        type=float,
        default=82.0,
        help="Starting angle of the open genomic arc, in degrees.",
    )
    parser.add_argument(
        "--sweep",
        type=float,
        default=296.0,
        help="Angular span of the genome, leaving an opening for the inset.",
    )
    parser.add_argument(
        "--max-background-points",
        type=int,
        default=350000,
        help="Maximum number of randomly sampled background SNPs.",
    )
    parser.add_argument(
        "--point-size",
        type=float,
        default=7.0,
        help="Background SNP point size.",
    )
    parser.add_argument(
        "--figure-size",
        type=float,
        default=13.0,
        help="Square figure size in inches.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=2026,
        help="Random seed for background sampling and word-cloud layout.",
    )

    return parser.parse_args()


def validate_arguments(arguments: argparse.Namespace) -> None:
    for path_argument in [
        arguments.gwas,
        arguments.loci,
        arguments.gwascatalog,
    ]:
        if not path_argument.exists():
            raise FileNotFoundError(f"Input file not found: {path_argument}")

    if not (0 < arguments.p_genome_wide <= 1):
        raise ValueError("--p-genome-wide must be between 0 and 1.")

    if not (0 < arguments.p_suggestive <= 1):
        raise ValueError("--p-suggestive must be between 0 and 1.")

    if arguments.p_genome_wide >= arguments.p_suggestive:
        raise ValueError(
            "--p-genome-wide should be smaller than --p-suggestive."
        )

    if arguments.wordcloud_min_frequency < 1:
        raise ValueError("--wordcloud-min-frequency must be at least 1.")

    if arguments.wordcloud_max_terms < 1:
        raise ValueError("--wordcloud-max-terms must be at least 1.")

    if arguments.max_background_points < 1:
        raise ValueError("--max-background-points must be at least 1.")


def main() -> None:
    arguments = parse_arguments()
    validate_arguments(arguments)

    print(f"Reading fastGWA results: {arguments.gwas}")
    gwas = read_fastgwa(arguments.gwas)
    print(f"Retained {len(gwas):,} valid GWAS variants.")

    print(f"Reading FUMA genomic loci: {arguments.loci}")
    loci = read_genomic_risk_loci(arguments.loci)
    print(f"Retained {len(loci):,} genomic loci.")

    print("Matching locus-level top SNPs to fastGWA results.")
    lead_snps = match_lead_snps_to_gwas(gwas, loci)

    lead_output = Path(
        f"{arguments.out_prefix}_lead_snps.tsv"
    )
    lead_output.parent.mkdir(parents=True, exist_ok=True)
    lead_snps.to_csv(lead_output, sep="\t", index=False)
    print(f"Wrote lead-SNP table: {lead_output}")

    match_summary = (
        lead_snps["MatchSource"]
        .value_counts(dropna=False)
        .rename_axis("MatchSource")
        .reset_index(name="N")
    )
    print("\nLead-SNP matching summary:")
    print(match_summary.to_string(index=False))

    print(f"\nReading FUMA GWAS Catalog annotations: {arguments.gwascatalog}")
    catalog = read_gwas_catalog(arguments.gwascatalog)
    print(f"Retained {len(catalog):,} valid GWAS Catalog rows.")

    wordcloud_frequencies, wordcloud_table = build_wordcloud_frequencies(
        catalog=catalog,
        valid_loci=loci["GenomicLocus"],
        mode=arguments.wordcloud_mode,
        count_method=arguments.wordcloud_count,
        minimum_frequency=arguments.wordcloud_min_frequency,
        maximum_terms=arguments.wordcloud_max_terms,
        include_other=arguments.include_other,
    )

    wordcloud_output = Path(
        f"{arguments.out_prefix}_wordcloud_frequencies.tsv"
    )
    wordcloud_table.to_csv(
        wordcloud_output,
        sep="\t",
        index=False,
    )
    print(f"Wrote word-cloud frequency table: {wordcloud_output}")

    print("\nTop word-cloud terms:")
    print(wordcloud_table.head(20).to_string(index=False))

    if arguments.wordcloud_mode == "category":
        wordcloud_title = (
            "Phenotypic categories previously associated with lead loci"
        )
    else:
        wordcloud_title = (
            "GWAS Catalog traits previously associated with lead loci"
        )

    plot_open_circular_manhattan(
        gwas=gwas,
        lead_snps=lead_snps,
        wordcloud_frequencies=wordcloud_frequencies,
        output_prefix=arguments.out_prefix,
        title=arguments.title,
        wordcloud_title=wordcloud_title,
        p_suggestive=arguments.p_suggestive,
        p_genome_wide=arguments.p_genome_wide,
        p_cap=arguments.p_cap,
        theta_start_degrees=arguments.theta_start,
        sweep_degrees=arguments.sweep,
        maximum_background_points=arguments.max_background_points,
        seed=arguments.seed,
        figure_size=arguments.figure_size,
        point_size=arguments.point_size,
    )


if __name__ == "__main__":
    main()