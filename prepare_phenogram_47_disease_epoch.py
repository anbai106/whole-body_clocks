#!/usr/bin/env python3

import argparse
import os
import re
from glob import glob
from typing import Dict, List, Optional

import numpy as np
import pandas as pd


# ============================================================
# Disease-clock settings
# ============================================================

DISEASE_ORDER = ["asthma", "copd", "dementia", "mi", "stroke"]

DISEASE_LABEL = {
    "asthma": "Asthma",
    "copd": "COPD",
    "dementia": "Dementia",
    "mi": "MI",
    "stroke": "Stroke",
}

DISEASE_ABBREV = {
    "asthma": "AST",
    "copd": "COPD",
    "dementia": "DEM",
    "mi": "MI",
    "stroke": "STR",
}

MODALITY_LABEL = {
    "mri": "MRI",
    "proteomics": "Proteomics",
    "metabolomics": "Metabolomics",
}

MODALITY_ORDER = {
    "MRI": 1,
    "Proteomics": 2,
    "Metabolomics": 3,
}

# PhenoGram shape file:
# GROUP = modality.
SHAPE_TABLE = pd.DataFrame(
    {
        "group": ["MRI", "Proteomics", "Metabolomics"],
        "shape": ["circle", "triangle", "square"],
    }
)

# Consistent organ/system color order across all 5 disease panels.
COLOR_GROUP_ORDER = [
    "Brain",
    "Eye",
    "Heart",
    "Liver",
    "Kidney",
    "Pulmonary",
    "Pancreas",
    "Spleen",
    "Adipose",
    "Endocrine",
    "Digestive",
    "Immune",
    "Metabolic",
    "Skin",
    "Reproductive Female",
    "Reproductive Male",
    "Other",
]

# R color names accepted by PhenoGram/base R.
COLOR_PALETTE = [
    "royalblue",
    "darkorange",
    "firebrick",
    "forestgreen",
    "purple",
    "deepskyblue3",
    "goldenrod",
    "sienna4",
    "gray55",
    "magenta",
    "turquoise3",
    "olivedrab3",
    "navy",
    "hotpink",
    "darkorchid",
    "black",
    "brown",
    "cyan3",
    "deeppink",
    "steelblue",
    "darkgoldenrod",
    "dodgerblue3",
]


# ============================================================
# General helpers
# ============================================================

def normalize_locus_id(x):
    if pd.isna(x):
        return np.nan

    s = str(x)

    if ":" in s:
        s = s.split(":")[0]

    s = re.sub(r"[^0-9]", "", s)

    if s == "":
        return np.nan

    return int(s)


def normalize_chr(x):
    if pd.isna(x):
        return np.nan

    s = str(x).strip()
    s = s.replace("chr", "").replace("CHR", "")
    s = s.replace("23", "X").replace("24", "Y")

    return s


def chr_sort_value(x):
    x = str(x).replace("chr", "").replace("CHR", "")

    if x == "X":
        return 23
    if x == "Y":
        return 24
    if x in ["M", "MT"]:
        return 25

    try:
        return int(x)
    except Exception:
        return 99


def find_first_existing_col(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    exact = {c: c for c in df.columns}
    lower = {c.lower(): c for c in df.columns}

    for c in candidates:
        if c in exact:
            return exact[c]

    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]

    return None


def format_original_organ_name(organ_raw):
    if pd.isna(organ_raw):
        return "Other"

    x = str(organ_raw).strip()
    x = x.replace("-", "_").replace(" ", "_")
    x = re.sub(r"_+", "_", x)

    words = []

    for token in x.split("_"):
        token_lower = token.lower()

        if token_lower == "mri":
            words.append("MRI")
        elif token_lower == "female":
            words.append("Female")
        elif token_lower == "male":
            words.append("Male")
        else:
            words.append(token_lower.capitalize())

    return " ".join(words)


def canonical_color_group_from_organ_raw(organ_raw):
    """
    Canonical organ/system group used only for color grouping.

    Disease is NOT encoded by color here because disease is encoded by panel.
    """

    if pd.isna(organ_raw):
        return "Other"

    x = str(organ_raw).strip().replace("-", "_").replace(" ", "_")
    x_lower = x.lower()

    color_group_map = {
        "brain": "Brain",
        "eye": "Eye",
        "heart": "Heart",
        "cardiac": "Heart",

        "liver": "Liver",
        "hepatic": "Liver",

        "kidney": "Kidney",
        "renal": "Kidney",

        "lung": "Pulmonary",
        "pulmonary": "Pulmonary",

        "pancreas": "Pancreas",
        "pancreatic": "Pancreas",

        "spleen": "Spleen",
        "adipose": "Adipose",

        "endocrine": "Endocrine",
        "digestive": "Digestive",
        "immune": "Immune",
        "metabolic": "Metabolic",
        "skin": "Skin",

        "reproductive_female": "Reproductive Female",
        "female_reproductive": "Reproductive Female",
        "reproductive_male": "Reproductive Male",
        "male_reproductive": "Reproductive Male",
    }

    return color_group_map.get(x_lower, format_original_organ_name(x))


def parse_disease_clock_folder(clock_name: str) -> Dict[str, str]:
    """
    Parse disease EPOCH folder names.

    Expected pattern:
      <organ>_<modality>_<disease>_clock

    Examples:
      Brain_proteomics_dementia_clock
      heart_mri_copd_clock
      Endocrine_metabolomics_asthma_clock
      Reproductive_male_proteomics_mi_clock
    """

    if not clock_name.endswith("_clock"):
        raise ValueError(f"Clock folder does not end with _clock: {clock_name}")

    x = clock_name[:-len("_clock")]

    disease = None
    x_without_disease = None

    for d in DISEASE_ORDER:
        if x.lower().endswith(f"_{d}"):
            disease = d
            x_without_disease = re.sub(f"_{d}$", "", x, flags=re.IGNORECASE)
            break

    if disease is None:
        raise ValueError(f"Cannot parse disease from folder: {clock_name}")

    modality = None
    organ_raw = None

    for m in ["metabolomics", "proteomics", "mri"]:
        if x_without_disease.lower().endswith(f"_{m}"):
            modality = MODALITY_LABEL[m]
            organ_raw = re.sub(f"_{m}$", "", x_without_disease, flags=re.IGNORECASE)
            break

    if modality is None or organ_raw is None:
        raise ValueError(f"Cannot parse modality/organ from folder: {clock_name}")

    organ_display = format_original_organ_name(organ_raw)
    color_group = canonical_color_group_from_organ_raw(organ_raw)
    disease_display = DISEASE_LABEL[disease]
    disease_abbrev = DISEASE_ABBREV[disease]

    # Combined file phenotype includes disease.
    phenotype_combined = f"{disease_display} {organ_display} {modality}"

    # Disease-panel phenotype excludes disease because disease is encoded by panel.
    phenotype_panel = f"{organ_display} {modality}"

    return {
        "clock": clock_name,
        "clock_prefix": x,
        "organ_raw": organ_raw,
        "organ_display": organ_display,
        "color_group": color_group,
        "modality": modality,
        "disease": disease,
        "disease_display": disease_display,
        "disease_abbrev": disease_abbrev,
        "phenotype_combined": phenotype_combined,
        "phenotype_panel": phenotype_panel,
    }


# ============================================================
# Cytoband annotation
# ============================================================

def read_cytoband(cyto_path: str) -> pd.DataFrame:
    df_cyto = pd.read_csv(
        cyto_path,
        sep="\t",
        header=None,
        compression="infer",
    )

    df_cyto.columns = ["chrom", "start", "end", "band", "stain"]

    df_cyto["chrom_clean"] = (
        df_cyto["chrom"]
        .astype(str)
        .str.replace("chr", "", regex=False)
        .str.replace("CHR", "", regex=False)
    )

    return df_cyto


def get_ucsc_cytoband(df_cyto: pd.DataFrame, chrom, pos):
    if pd.isna(chrom) or pd.isna(pos):
        return "NA"

    chrom = normalize_chr(chrom)
    pos = int(pos)

    pos_row = df_cyto.loc[
        (df_cyto["chrom_clean"].astype(str) == str(chrom))
        & (df_cyto["start"] <= pos)
        & (df_cyto["end"] >= pos)
    ]

    if pos_row.shape[0] >= 1:
        return f"{chrom}{pos_row['band'].values[0]}"

    return "NA"


# ============================================================
# FUMA gene annotation
# ============================================================

def is_nonempty_mapping_value(x):
    if pd.isna(x):
        return False

    s = str(x).strip()

    return s not in ["", ".", "NA", "NaN", "nan", "None", "0"]


def read_physical_position_genes(genes_file: str) -> pd.DataFrame:
    empty = pd.DataFrame(
        columns=["GenomicLocus_norm", "MappedGene_PhysicalPosition"]
    )

    if not os.path.exists(genes_file) or os.path.getsize(genes_file) == 0:
        return empty

    try:
        df_gene = pd.read_csv(genes_file, sep="\t")
    except Exception:
        return empty

    if df_gene.empty:
        return empty

    locus_col = find_first_existing_col(df_gene, ["GenomicLocus", "GenomicLocusID"])
    symbol_col = find_first_existing_col(df_gene, ["symbol", "Symbol", "gene", "Gene"])

    if locus_col is None or symbol_col is None:
        return empty

    df_gene = df_gene.copy()
    df_gene["GenomicLocus_norm"] = df_gene[locus_col].apply(normalize_locus_id)
    df_gene = df_gene.dropna(subset=["GenomicLocus_norm", symbol_col])

    if df_gene.empty:
        return empty

    # Prefer FUMA physical-position mapped genes if available.
    pos_col = find_first_existing_col(df_gene, ["posMapSNPs", "posMapSNP", "posMap"])

    if pos_col is not None:
        df_gene_pos = df_gene.loc[
            df_gene[pos_col].apply(is_nonempty_mapping_value)
        ].copy()

        if not df_gene_pos.empty:
            df_gene = df_gene_pos

    df_gene[symbol_col] = df_gene[symbol_col].astype(str)

    df_gene = (
        df_gene.groupby("GenomicLocus_norm", as_index=False)[symbol_col]
        .agg(lambda x: ";".join(sorted(set(x))))
        .rename(columns={symbol_col: "MappedGene_PhysicalPosition"})
    )

    return df_gene


# ============================================================
# Color mapping and PhenoGram writing
# ============================================================

def build_global_organ_color_map(df_final: pd.DataFrame) -> Dict[str, str]:
    observed = sorted(df_final["ColorGroup"].dropna().unique().tolist())

    ordered = [g for g in COLOR_GROUP_ORDER if g in observed]
    ordered += [g for g in observed if g not in ordered]

    color_map = {}

    for i, group in enumerate(ordered):
        color_map[group] = COLOR_PALETTE[i % len(COLOR_PALETTE)]

    return color_map


def write_color_spec(
    df_phenogram: pd.DataFrame,
    color_map: Dict[str, str],
    out_path: str,
):
    """
    PhenoGram color-spec file:
      phenotype    color

    We map each phenotype to the organ/system color through COLORGROUP.
    """

    color_spec = (
        df_phenogram[["PHENOTYPE", "COLORGROUP"]]
        .drop_duplicates()
        .sort_values(["COLORGROUP", "PHENOTYPE"])
        .copy()
    )

    color_spec["color"] = color_spec["COLORGROUP"].map(color_map)

    color_spec = (
        color_spec[["PHENOTYPE", "color"]]
        .rename(columns={"PHENOTYPE": "phenotype"})
    )

    color_spec.to_csv(out_path, index=False, sep="\t", encoding="utf-8")


def make_phenogram_table(
    df: pd.DataFrame,
    phenotype_col: str,
    annotation_col: str = "Cytoband_UCSC",
) -> pd.DataFrame:
    out = df[
        [
            "TopLeadSNP",
            "Chromosome",
            "Position",
            phenotype_col,
            annotation_col,
            "Modality",
            "ColorGroup",
        ]
    ].copy()

    out = out.rename(
        columns={
            "TopLeadSNP": "SNP",
            "Chromosome": "CHR",
            "Position": "POS",
            phenotype_col: "PHENOTYPE",
            annotation_col: "ANNOTATION",
            "Modality": "GROUP",
            "ColorGroup": "COLORGROUP",
        }
    )

    out["CHR"] = out["CHR"].apply(normalize_chr)
    out["POS"] = out["POS"].astype(int)

    return out


def write_phenogram_outputs(
    df_phenogram: pd.DataFrame,
    out_path: str,
):
    df_phenogram.to_csv(out_path, index=False, sep="\t", encoding="utf-8")

    legacy_out = out_path.replace(".tsv", "_legacy_lowercase.tsv")

    df_legacy = df_phenogram.rename(
        columns={
            "SNP": "snp",
            "CHR": "chr",
            "POS": "pos",
            "PHENOTYPE": "phenotype",
            "ANNOTATION": "annotation",
            "GROUP": "group",
            "COLORGROUP": "colorgroup",
        }
    )

    df_legacy.to_csv(legacy_out, index=False, sep="\t", encoding="utf-8")

    return legacy_out


# ============================================================
# Main annotation
# ============================================================

def discover_fuma_locus_files(base_dir: str) -> List[str]:
    files = []

    for disease in DISEASE_ORDER:
        pat = os.path.join(
            base_dir,
            f"*_{disease}_clock",
            "fuma",
            "GenomicRiskLoci.txt",
        )
        files.extend(glob(pat))

    return sorted(set(files))


def read_fuma_genomic_risk_loci(genomic_loci_file: str) -> pd.DataFrame:
    df_loci = pd.read_csv(genomic_loci_file, sep="\t")

    if df_loci.empty:
        return df_loci

    col_map = {
        "GenomicLocus": find_first_existing_col(
            df_loci,
            ["GenomicLocus", "GenomicLocusID", "locus"],
        ),
        "rsID": find_first_existing_col(
            df_loci,
            ["rsID", "rsid", "rsIDuniq", "leadSNP", "LeadSNP"],
        ),
        "chr": find_first_existing_col(
            df_loci,
            ["chr", "CHR", "chrom", "Chromosome"],
        ),
        "pos": find_first_existing_col(
            df_loci,
            ["pos", "POS", "bp", "BP", "position", "Position"],
        ),
        "p": find_first_existing_col(
            df_loci,
            ["p", "P", "P-value", "P_value", "pvalue", "PVAL", "pval"],
        ),
    }

    missing = [k for k, v in col_map.items() if v is None]

    if missing:
        raise ValueError(f"Missing required FUMA columns: {missing}")

    out = df_loci[
        [
            col_map["GenomicLocus"],
            col_map["rsID"],
            col_map["chr"],
            col_map["pos"],
            col_map["p"],
        ]
    ].copy()

    out.columns = ["GenomicLocus", "rsID", "chr", "pos", "p"]

    return out


def annotate_genomic_loci_disease_epoch_clocks_for_phenogram_panels(
    base_dir: str,
    output_dir_result: str,
    cyto_path: str,
    p_threshold: float,
    expected_n_clocks: int = 47,
):
    """
    Prepare PhenoGram input for 47 disease-specific EPOCH clocks.

    Final visual encoding:
      Disease       = separate PhenoGram panel/file
      Modality      = GROUP / shape
      Organ-system  = COLORGROUP / color
      Annotation    = cytoband
    """

    os.makedirs(output_dir_result, exist_ok=True)

    print("============================================================")
    print("Preparing PhenoGram input for 47 disease EPOCH clock panels")
    print("============================================================")
    print(f"Base directory:          {base_dir}")
    print(f"Output directory:        {output_dir_result}")
    print(f"Cytoband file:           {cyto_path}")
    print(f"Significance threshold:  {p_threshold:.3e}")
    print(f"Expected clocks:         {expected_n_clocks}")

    if not os.path.isdir(base_dir):
        raise FileNotFoundError(f"Cannot find base directory: {base_dir}")

    if not os.path.exists(cyto_path):
        raise FileNotFoundError(f"Cannot find cytoband file: {cyto_path}")

    df_cyto = read_cytoband(cyto_path)
    genomic_locus_files = discover_fuma_locus_files(base_dir)

    print(
        f"\nNumber of disease EPOCH FUMA GenomicRiskLoci files found: "
        f"{len(genomic_locus_files)}"
    )

    for f in genomic_locus_files:
        print("  -", os.path.relpath(f, base_dir))

    if len(genomic_locus_files) != expected_n_clocks:
        print(
            f"\nWARNING: Expected {expected_n_clocks} disease EPOCH clocks, "
            f"but found {len(genomic_locus_files)} GenomicRiskLoci files."
        )

    all_loci = []
    skipped = []

    for genomic_loci_file in genomic_locus_files:
        fuma_dir = os.path.dirname(genomic_loci_file)
        clock_dir = os.path.dirname(fuma_dir)
        clock_name = os.path.basename(clock_dir)
        genes_file = os.path.join(fuma_dir, "genes.txt")

        try:
            meta = parse_disease_clock_folder(clock_name)
        except Exception as e:
            skipped.append((clock_name, "parse_failed", str(e)))
            continue

        try:
            df_loci = read_fuma_genomic_risk_loci(genomic_loci_file)
        except Exception as e:
            skipped.append((clock_name, "cannot_read_GenomicRiskLoci", str(e)))
            continue

        if df_loci.empty:
            skipped.append((clock_name, "empty_GenomicRiskLoci", ""))
            continue

        df_loci["p"] = pd.to_numeric(df_loci["p"], errors="coerce")
        df_loci["pos"] = pd.to_numeric(df_loci["pos"], errors="coerce")
        df_loci["chr"] = df_loci["chr"].apply(normalize_chr)

        df_loci = df_loci.dropna(subset=["p", "pos", "chr", "rsID"]).copy()
        df_loci = df_loci.loc[df_loci["p"] < p_threshold].copy()

        if df_loci.empty:
            skipped.append(
                (
                    clock_name,
                    "no_loci_below_threshold",
                    f"p<{p_threshold:.3e}",
                )
            )
            continue

        df_loci["GenomicLocus_norm"] = df_loci["GenomicLocus"].apply(
            normalize_locus_id
        )

        df_gene = read_physical_position_genes(genes_file)

        df_loci = df_loci.merge(df_gene, how="left", on="GenomicLocus_norm")

        df_loci["MappedGene_PhysicalPosition"] = (
            df_loci["MappedGene_PhysicalPosition"].fillna("NA")
        )

        df_loci = df_loci.rename(
            columns={
                "rsID": "TopLeadSNP",
                "chr": "Chromosome",
                "pos": "Position",
                "p": "P-value",
            }
        )

        df_loci["Clock"] = meta["clock"]
        df_loci["ClockPrefix"] = meta["clock_prefix"]
        df_loci["Disease"] = meta["disease"]
        df_loci["DiseaseDisplay"] = meta["disease_display"]
        df_loci["DiseaseAbbrev"] = meta["disease_abbrev"]
        df_loci["OrganRaw"] = meta["organ_raw"]
        df_loci["OrganDisplay"] = meta["organ_display"]
        df_loci["ColorGroup"] = meta["color_group"]
        df_loci["Modality"] = meta["modality"]

        df_loci["PhenotypeCombined"] = meta["phenotype_combined"]
        df_loci["PhenotypePanel"] = meta["phenotype_panel"]

        df_loci["Cytoband_UCSC"] = df_loci.apply(
            lambda row: get_ucsc_cytoband(
                df_cyto=df_cyto,
                chrom=row["Chromosome"],
                pos=row["Position"],
            ),
            axis=1,
        )

        # Combined plot can optionally display disease prefix in annotation.
        df_loci["AnnotationCombined"] = (
            df_loci["DiseaseAbbrev"].astype(str)
            + ":"
            + df_loci["Cytoband_UCSC"].astype(str)
        )

        keep_cols = [
            "Clock",
            "ClockPrefix",
            "Disease",
            "DiseaseDisplay",
            "DiseaseAbbrev",
            "GenomicLocus",
            "TopLeadSNP",
            "Chromosome",
            "Position",
            "P-value",
            "MappedGene_PhysicalPosition",
            "Cytoband_UCSC",
            "AnnotationCombined",
            "PhenotypeCombined",
            "PhenotypePanel",
            "OrganRaw",
            "OrganDisplay",
            "ColorGroup",
            "Modality",
        ]

        all_loci.append(df_loci[keep_cols].copy())

    skipped_out = os.path.join(
        output_dir_result,
        "Fuma_loci_annotation_47_disease_epoch_panels_skipped.tsv",
    )

    pd.DataFrame(
        skipped,
        columns=["Clock", "Reason", "Details"],
    ).to_csv(skipped_out, index=False, sep="\t", encoding="utf-8")

    if len(all_loci) == 0:
        raise RuntimeError(
            "No significant loci found across disease EPOCH clocks. "
            f"Skipped log written to: {skipped_out}"
        )

    df_final = pd.concat(all_loci, ignore_index=True)

    disease_order_map = {d: i + 1 for i, d in enumerate(DISEASE_ORDER)}

    df_final["disease_order"] = df_final["Disease"].map(disease_order_map).fillna(99)
    df_final["modality_order"] = df_final["Modality"].map(MODALITY_ORDER).fillna(99)
    df_final["chr_order"] = df_final["Chromosome"].apply(chr_sort_value)

    df_final = (
        df_final.sort_values(
            [
                "disease_order",
                "modality_order",
                "ColorGroup",
                "PhenotypePanel",
                "chr_order",
                "Position",
                "P-value",
            ]
        )
        .drop(columns=["disease_order", "modality_order", "chr_order"])
        .reset_index(drop=True)
    )

    # ========================================================
    # Detailed all-loci output
    # ========================================================

    detailed_out = os.path.join(
        output_dir_result,
        "Fuma_loci_annotation_47_disease_epoch_panels.tsv",
    )

    df_final.to_csv(detailed_out, index=False, sep="\t", encoding="utf-8")

    # ========================================================
    # Global shape and color files
    # ========================================================

    shape_out = os.path.join(
        output_dir_result,
        "PhenoGram_47_disease_epoch_group_shape_modality.tsv",
    )

    SHAPE_TABLE.to_csv(shape_out, index=False, sep="\t", encoding="utf-8")

    color_map = build_global_organ_color_map(df_final)

    color_map_out = os.path.join(
        output_dir_result,
        "PhenoGram_47_disease_epoch_global_organ_color_map.tsv",
    )

    pd.DataFrame(
        [{"COLORGROUP": k, "color": v} for k, v in color_map.items()]
    ).to_csv(color_map_out, index=False, sep="\t", encoding="utf-8")

    # ========================================================
    # Optional combined file
    # Disease is included in PHENOTYPE and annotation.
    # This is less recommended than panel files, but useful for QC.
    # ========================================================

    df_combined = make_phenogram_table(
        df_final,
        phenotype_col="PhenotypeCombined",
        annotation_col="AnnotationCombined",
    )

    combined_out = os.path.join(
        output_dir_result,
        "PhenoGram_47_disease_epoch_combined_qc.tsv",
    )

    combined_legacy_out = write_phenogram_outputs(
        df_phenogram=df_combined,
        out_path=combined_out,
    )

    combined_color_spec_out = os.path.join(
        output_dir_result,
        "PhenoGram_47_disease_epoch_combined_qc_color_spec_same_organ_color.tsv",
    )

    write_color_spec(
        df_phenogram=df_combined,
        color_map=color_map,
        out_path=combined_color_spec_out,
    )

    # ========================================================
    # Recommended disease-specific panel files
    # Disease = panel
    # Modality = shape
    # Organ/system = color
    # ========================================================

    panel_manifest_rows = []

    for disease in DISEASE_ORDER:
        disease_label = DISEASE_LABEL[disease]

        df_disease = df_final.loc[df_final["Disease"] == disease].copy()

        if df_disease.empty:
            panel_manifest_rows.append(
                {
                    "disease": disease,
                    "disease_label": disease_label,
                    "n_loci": 0,
                    "n_clocks": 0,
                    "phenogram_file": "NA",
                    "legacy_file": "NA",
                    "color_spec_file": "NA",
                    "status": "no_loci_retained",
                }
            )
            continue

        df_panel = make_phenogram_table(
            df_disease,
            phenotype_col="PhenotypePanel",
            annotation_col="Cytoband_UCSC",
        )

        panel_out = os.path.join(
            output_dir_result,
            f"PhenoGram_panel_{disease}_disease_epoch.tsv",
        )

        panel_legacy_out = write_phenogram_outputs(
            df_phenogram=df_panel,
            out_path=panel_out,
        )

        panel_color_spec_out = os.path.join(
            output_dir_result,
            f"PhenoGram_panel_{disease}_disease_epoch_color_spec_same_organ_color.tsv",
        )

        write_color_spec(
            df_phenogram=df_panel,
            color_map=color_map,
            out_path=panel_color_spec_out,
        )

        panel_manifest_rows.append(
            {
                "disease": disease,
                "disease_label": disease_label,
                "n_loci": df_panel.shape[0],
                "n_clocks": df_disease["Clock"].nunique(),
                "phenogram_file": panel_out,
                "legacy_file": panel_legacy_out,
                "color_spec_file": panel_color_spec_out,
                "status": "ok",
            }
        )

    panel_manifest = pd.DataFrame(panel_manifest_rows)

    panel_manifest_out = os.path.join(
        output_dir_result,
        "PhenoGram_47_disease_epoch_panel_manifest.tsv",
    )

    panel_manifest.to_csv(
        panel_manifest_out,
        index=False,
        sep="\t",
        encoding="utf-8",
    )

    # ========================================================
    # Summary files
    # ========================================================

    summary = (
        df_final.groupby(
            [
                "DiseaseDisplay",
                "Modality",
                "PhenotypePanel",
                "ColorGroup",
                "Clock",
            ],
            as_index=False,
        )
        .agg(
            n_loci=("TopLeadSNP", "count"),
            min_p=("P-value", "min"),
            n_loci_with_physical_gene=(
                "MappedGene_PhysicalPosition",
                lambda x: (x.astype(str) != "NA").sum(),
            ),
        )
        .sort_values(
            [
                "DiseaseDisplay",
                "Modality",
                "ColorGroup",
                "PhenotypePanel",
                "Clock",
            ]
        )
    )

    summary_out = os.path.join(
        output_dir_result,
        "Fuma_loci_annotation_47_disease_epoch_panels_summary.tsv",
    )

    summary.to_csv(summary_out, index=False, sep="\t", encoding="utf-8")

    # ========================================================
    # Final report
    # ========================================================

    print("\nFinished.")
    print(f"Detailed loci file:              {detailed_out}")
    print(f"Combined QC PhenoGram input:     {combined_out}")
    print(f"Combined QC legacy input:        {combined_legacy_out}")
    print(f"Combined QC color spec:          {combined_color_spec_out}")
    print(f"Panel manifest:                  {panel_manifest_out}")
    print(f"Group shape file:                {shape_out}")
    print(f"Global organ color map:          {color_map_out}")
    print(f"Summary file:                    {summary_out}")
    print(f"Skipped log:                     {skipped_out}")
    print(f"Total significant loci retained: {df_final.shape[0]}")
    print(f"Total clocks with loci retained: {df_final['Clock'].nunique()}")
    print(f"Total disease endpoints:         {df_final['DiseaseDisplay'].nunique()}")

    print("\nRecommended PhenoGram setup for the paper figure:")
    print("  Use one PhenoGram input file per disease panel:")
    print("    PhenoGram_panel_asthma_disease_epoch.tsv")
    print("    PhenoGram_panel_copd_disease_epoch.tsv")
    print("    PhenoGram_panel_dementia_disease_epoch.tsv")
    print("    PhenoGram_panel_mi_disease_epoch.tsv")
    print("    PhenoGram_panel_stroke_disease_epoch.tsv")
    print("")
    print("  Encoding:")
    print("    Disease       = separate panel")
    print("    Modality      = GROUP / shape")
    print("    Organ-system  = COLORGROUP / color")
    print("    Annotation    = cytoband")
    print("")
    print("  Upload the same group-shape file for every disease panel:")
    print(f"    {shape_out}")
    print("")
    print("  Upload the corresponding disease-specific color-spec file")
    print("  listed in the panel manifest if you want identical organ colors")
    print("  across all 5 panels:")
    print(f"    {panel_manifest_out}")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Prepare disease-specific PhenoGram panel input files from FUMA "
            "GenomicRiskLoci.txt files for 47 disease-specific EPOCH clocks."
        )
    )

    parser.add_argument(
        "--base_dir",
        default="/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock",
        type=str,
        help="WholeBodyClock base directory containing disease EPOCH clock folders.",
    )

    parser.add_argument(
        "--outdir",
        default="/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/Result",
        type=str,
        help="Output directory for PhenoGram input files.",
    )

    parser.add_argument(
        "--cyto_path",
        default="/Users/hao/cubic-home/Dataset/GRch37_cytoband/cytoBand.txt.gz",
        type=str,
        help="UCSC GRCh37/hg19 cytoBand.txt or cytoBand.txt.gz file.",
    )

    parser.add_argument(
        "--p_threshold",
        default=5e-8 / 47,
        type=float,
        help="P-value threshold for retaining FUMA loci. Default: 5e-8/47.",
    )

    parser.add_argument(
        "--expected_n_clocks",
        default=47,
        type=int,
        help="Expected number of disease EPOCH clocks with FUMA outputs.",
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    annotate_genomic_loci_disease_epoch_clocks_for_phenogram_panels(
        base_dir=args.base_dir,
        output_dir_result=args.outdir,
        cyto_path=args.cyto_path,
        p_threshold=args.p_threshold,
        expected_n_clocks=args.expected_n_clocks,
    )