import os
import re
import numpy as np
import pandas as pd


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


def format_original_organ_name(organ_raw):
    """
    Preserve original clock wording, but make it figure-ready.

    Examples:
      brain -> Brain
      adipose -> Adipose
      pulmonary -> Pulmonary
      hepatic -> Hepatic
      renal -> Renal
      reproductive_female -> Reproductive Female
      reproductive_male -> Reproductive Male

    This function does NOT map:
      hepatic -> liver
      renal -> kidney
      pulmonary -> lung
      reproductive_female / reproductive_male -> reproductive
    """

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
    Canonical organ/system group used ONLY for color grouping.

    This ensures:
      Liver and Hepatic share the same color group: Liver
      Kidney and Renal share the same color group: Kidney
      Lung and Pulmonary share the same color group: Pulmonary

    Reproductive Female and Reproductive Male are kept separate.
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


def parse_clock_name(clock_name):
    """
    Convert FUMA folder name into:
      original organ label for legend,
      modality for shape,
      original phenotype label for PhenoGram legend,
      canonical organ color group for shared colors.

    Examples:
      liver_mri_mortality_clock
        -> organ_display = Liver
        -> modality = MRI
        -> phenotype = Liver MRI
        -> color_group = Liver

      Hepatic_metabolomics_mortality_clock
        -> organ_display = Hepatic
        -> modality = Metabolomics
        -> phenotype = Hepatic Metabolomics
        -> color_group = Liver

      Reproductive_female_proteomics_mortality_clock
        -> organ_display = Reproductive Female
        -> modality = Proteomics
        -> phenotype = Reproductive Female Proteomics
        -> color_group = Reproductive Female
    """

    x = clock_name.replace("_mortality_clock", "")

    if x.lower().endswith("_mri"):
        modality = "MRI"
        organ_raw = re.sub("_mri$", "", x, flags=re.IGNORECASE)

    elif x.lower().endswith("_proteomics"):
        modality = "Proteomics"
        organ_raw = re.sub("_proteomics$", "", x, flags=re.IGNORECASE)

    elif x.lower().endswith("_metabolomics"):
        modality = "Metabolomics"
        organ_raw = re.sub("_metabolomics$", "", x, flags=re.IGNORECASE)

    else:
        modality = "Unknown"
        organ_raw = x

    organ_display = format_original_organ_name(organ_raw)
    color_group = canonical_color_group_from_organ_raw(organ_raw)
    phenotype = f"{organ_display} {modality}"

    return organ_raw, organ_display, modality, phenotype, color_group


def read_cytoband(cyto_path):
    """
    UCSC cytoband file columns:
      chrom, chromStart, chromEnd, name, gieStain
    """

    df_cyto = pd.read_csv(cyto_path, sep="\t", header=None)
    df_cyto.columns = ["chrom", "start", "end", "band", "stain"]

    df_cyto["chrom_clean"] = (
        df_cyto["chrom"].astype(str).str.replace("chr", "", regex=False)
    )

    return df_cyto


def get_ucsc_cytoband(df_cyto, chrom, pos):
    """
    Return cytoband annotation in PhenoGram-friendly format:
      chr1 + p36.33 -> 1p36.33
      chrX + q28    -> Xq28
    """

    chrom = str(chrom).replace("chr", "")
    pos = int(pos)

    pos_row = df_cyto.loc[
        (df_cyto["chrom_clean"] == chrom)
        & (df_cyto["start"] <= pos)
        & (df_cyto["end"] >= pos)
    ]

    if pos_row.shape[0] >= 1:
        return f"{chrom}{pos_row['band'].values[0]}"

    return "NA"


def is_nonempty_mapping_value(x):
    if pd.isna(x):
        return False

    s = str(x).strip()

    return s not in ["", ".", "NA", "NaN", "nan", "None", "0"]


def read_physical_position_genes(genes_file):
    """
    Read FUMA genes.txt and return physical-position mapped genes by GenomicLocus.

    These genes are saved in the detailed annotation file, but cytobands are used
    as PhenoGram ANNOTATION to satisfy the 10-character label limit.
    """

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

    if "symbol" not in df_gene.columns or "GenomicLocus" not in df_gene.columns:
        return empty

    df_gene = df_gene.copy()
    df_gene["GenomicLocus_norm"] = df_gene["GenomicLocus"].apply(normalize_locus_id)
    df_gene = df_gene.dropna(subset=["GenomicLocus_norm", "symbol"])

    if df_gene.empty:
        return empty

    # Prefer FUMA physical-position mapped genes if available.
    if "posMapSNPs" in df_gene.columns:
        df_gene_pos = df_gene.loc[
            df_gene["posMapSNPs"].apply(is_nonempty_mapping_value)
        ].copy()

        if not df_gene_pos.empty:
            df_gene = df_gene_pos

    df_gene["symbol"] = df_gene["symbol"].astype(str)

    df_gene = (
        df_gene.groupby("GenomicLocus_norm", as_index=False)["symbol"]
        .agg(lambda x: ";".join(sorted(set(x))))
        .rename(columns={"symbol": "MappedGene_PhysicalPosition"})
    )

    return df_gene


def build_phenogram_color_spec(df_phenogram, output_dir_result):
    """
    Build optional PhenoGram color-spec file.

    This guarantees that phenotypes sharing a canonical COLORGROUP have the
    exact same color while keeping original phenotype labels in the legend.

    The file format follows PhenoGram's color-spec example:
      phenotype    color
      Brain MRI    blue
      Brain Proteomics    blue
    """

    # R color names accepted by base plotting. Keep enough distinct colors.
    palette = [
        "blue",
        "red",
        "green3",
        "navy",
        "magenta",
        "darkgreen",
        "orange",
        "gold",
        "turquoise3",
        "sienna4",
        "hotpink",
        "olivedrab3",
        "gray55",
        "purple",
        "cyan3",
        "brown",
        "black",
        "deeppink",
        "darkorange",
        "steelblue",
    ]

    color_groups = sorted(df_phenogram["COLORGROUP"].dropna().unique().tolist())
    color_map = {
        group: palette[i % len(palette)]
        for i, group in enumerate(color_groups)
    }

    color_spec = (
        df_phenogram[["PHENOTYPE", "COLORGROUP"]]
        .drop_duplicates()
        .sort_values(["COLORGROUP", "PHENOTYPE"])
        .copy()
    )

    color_spec["color"] = color_spec["COLORGROUP"].map(color_map)
    color_spec = color_spec[["PHENOTYPE", "color"]]
    color_spec = color_spec.rename(columns={"PHENOTYPE": "phenotype"})

    color_spec_out = os.path.join(
        output_dir_result,
        "PhenoGram_mortality_clocks_color_spec_same_organ_color.tsv",
    )

    color_spec.to_csv(color_spec_out, index=False, sep="\t", encoding="utf-8")

    color_map_out = os.path.join(
        output_dir_result,
        "PhenoGram_mortality_clocks_color_group_map.tsv",
    )

    pd.DataFrame(
        [{"COLORGROUP": k, "color": v} for k, v in color_map.items()]
    ).to_csv(color_map_out, index=False, sep="\t", encoding="utf-8")

    return color_spec_out, color_map_out


def annotate_genomic_loci_mortality_clocks_for_phenogram(
    output_dir_fuma,
    output_dir_result,
    cyto_path,
    p_threshold=5e-8 / 22,
    expected_n_clocks=22,
):
    """
    Prepare PhenoGram input for significant loci across mortality clocks.

    Main PhenoGram columns:
      CHR         = chromosome
      POS         = lead SNP position
      PHENOTYPE   = original display label, e.g., Hepatic Metabolomics
      ANNOTATION  = cytoband, e.g., 1p36.33
      GROUP       = modality, e.g., MRI / Proteomics / Metabolomics
      COLORGROUP  = canonical organ color group, e.g., Liver

    This means:
      - Legend uses original clock names.
      - Hepatic and Liver can share the same color via COLORGROUP or color-spec file.
      - Reproductive Female and Reproductive Male are kept separate.
      - Shape is controlled by modality through GROUP.
    """

    os.makedirs(output_dir_result, exist_ok=True)

    print(f"Using significance threshold: {p_threshold:.3e}")

    df_cyto = read_cytoband(cyto_path)

    clock_dirs = [
        d for d in os.listdir(output_dir_fuma)
        if os.path.isdir(os.path.join(output_dir_fuma, d))
        and d.endswith("_mortality_clock")
    ]
    clock_dirs = sorted(clock_dirs)

    print(f"Number of mortality-clock FUMA folders found: {len(clock_dirs)}")
    for d in clock_dirs:
        print("  -", d)

    if len(clock_dirs) != expected_n_clocks:
        print(
            f"WARNING: Expected {expected_n_clocks} mortality clocks, "
            f"but found {len(clock_dirs)} folders."
        )

    all_loci = []
    skipped = []

    for clock_name in clock_dirs:
        clock_dir = os.path.join(output_dir_fuma, clock_name)

        genomic_loci_file = os.path.join(clock_dir, "GenomicRiskLoci.txt")
        genes_file = os.path.join(clock_dir, "genes.txt")

        organ_raw, organ_display, modality, phenotype, color_group = parse_clock_name(clock_name)

        if not os.path.exists(genomic_loci_file):
            skipped.append((clock_name, "missing GenomicRiskLoci.txt"))
            continue

        try:
            df_loci = pd.read_csv(genomic_loci_file, sep="\t")
        except Exception as e:
            skipped.append((clock_name, f"cannot read GenomicRiskLoci.txt: {e}"))
            continue

        if df_loci.empty:
            skipped.append((clock_name, "empty GenomicRiskLoci.txt"))
            continue

        required_cols = ["GenomicLocus", "rsID", "chr", "pos", "p"]
        missing_cols = [c for c in required_cols if c not in df_loci.columns]

        if missing_cols:
            skipped.append((clock_name, f"missing columns: {missing_cols}"))
            continue

        df_loci = df_loci[required_cols].copy()
        df_loci["p"] = pd.to_numeric(df_loci["p"], errors="coerce")
        df_loci["pos"] = pd.to_numeric(df_loci["pos"], errors="coerce")

        df_loci = df_loci.dropna(subset=["p", "pos", "chr", "rsID"])
        df_loci = df_loci.loc[df_loci["p"] < p_threshold].copy()

        if df_loci.empty:
            skipped.append((clock_name, "no significant loci"))
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

        df_loci["Clock"] = clock_name
        df_loci["ClockPrefix"] = clock_name.replace("_mortality_clock", "")
        df_loci["OrganRaw"] = organ_raw
        df_loci["OrganDisplay"] = organ_display
        df_loci["ColorGroup"] = color_group
        df_loci["Modality"] = modality
        df_loci["Phenotype"] = phenotype

        df_loci["Cytoband_UCSC"] = df_loci.apply(
            lambda row: get_ucsc_cytoband(
                df_cyto=df_cyto,
                chrom=row["Chromosome"],
                pos=row["Position"],
            ),
            axis=1,
        )

        df_loci["PHENOGRAM_CHR"] = (
            df_loci["Chromosome"].astype(str).str.replace("chr", "", regex=False)
        )
        df_loci["PHENOGRAM_POS"] = df_loci["Position"].astype(int)
        df_loci["PHENOGRAM_PHENOTYPE"] = df_loci["Phenotype"]
        df_loci["PHENOGRAM_ANNOTATION"] = df_loci["Cytoband_UCSC"]
        df_loci["PHENOGRAM_GROUP"] = df_loci["Modality"]
        df_loci["PHENOGRAM_COLORGROUP"] = df_loci["ColorGroup"]

        keep_cols = [
            "Clock",
            "ClockPrefix",
            "GenomicLocus",
            "TopLeadSNP",
            "Chromosome",
            "Position",
            "P-value",
            "MappedGene_PhysicalPosition",
            "Cytoband_UCSC",
            "Phenotype",
            "OrganRaw",
            "OrganDisplay",
            "ColorGroup",
            "Modality",
            "PHENOGRAM_CHR",
            "PHENOGRAM_POS",
            "PHENOGRAM_PHENOTYPE",
            "PHENOGRAM_ANNOTATION",
            "PHENOGRAM_GROUP",
            "PHENOGRAM_COLORGROUP",
        ]

        df_loci = df_loci[keep_cols]
        all_loci.append(df_loci)

    if len(all_loci) == 0:
        raise RuntimeError("No significant loci found across mortality clocks.")

    df_final = pd.concat(all_loci, ignore_index=True)

    modality_order = {"MRI": 1, "Proteomics": 2, "Metabolomics": 3}
    df_final["modality_order"] = df_final["Modality"].map(modality_order).fillna(99)
    df_final["Chromosome_numeric"] = pd.to_numeric(
        df_final["PHENOGRAM_CHR"], errors="coerce"
    )

    df_final = (
        df_final.sort_values(
            [
                "modality_order",
                "ColorGroup",
                "Phenotype",
                "Chromosome_numeric",
                "PHENOGRAM_POS",
                "P-value",
            ]
        )
        .drop(columns=["modality_order", "Chromosome_numeric"])
        .reset_index(drop=True)
    )

    detailed_out = os.path.join(
        output_dir_result,
        "Fuma_loci_annotation_mortality_clocks.tsv",
    )
    df_final.to_csv(detailed_out, index=False, sep="\t", encoding="utf-8")

    df_phenogram = df_final[
        [
            "TopLeadSNP",
            "PHENOGRAM_CHR",
            "PHENOGRAM_POS",
            "PHENOGRAM_PHENOTYPE",
            "PHENOGRAM_ANNOTATION",
            "PHENOGRAM_GROUP",
            "PHENOGRAM_COLORGROUP",
        ]
    ].copy()

    df_phenogram = df_phenogram.rename(
        columns={
            "TopLeadSNP": "SNP",
            "PHENOGRAM_CHR": "CHR",
            "PHENOGRAM_POS": "POS",
            "PHENOGRAM_PHENOTYPE": "PHENOTYPE",
            "PHENOGRAM_ANNOTATION": "ANNOTATION",
            "PHENOGRAM_GROUP": "GROUP",
            "PHENOGRAM_COLORGROUP": "COLORGROUP",
        }
    )

    phenogram_out = os.path.join(
        output_dir_result,
        "PhenoGram_mortality_clocks_original_labels_same_organ_colors.tsv",
    )
    df_phenogram.to_csv(phenogram_out, index=False, sep="\t", encoding="utf-8")

    df_phenogram_legacy = df_phenogram.rename(
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

    phenogram_legacy_out = os.path.join(
        output_dir_result,
        "PhenoGram_mortality_clocks_original_labels_same_organ_colors_legacy_lowercase.tsv",
    )
    df_phenogram_legacy.to_csv(
        phenogram_legacy_out,
        index=False,
        sep="\t",
        encoding="utf-8",
    )

    shape_df = pd.DataFrame(
        {
            "group": ["MRI", "Proteomics", "Metabolomics"],
            "shape": ["circle", "triangle", "square"],
        }
    )

    shape_out = os.path.join(
        output_dir_result,
        "PhenoGram_mortality_clocks_group_shape.tsv",
    )
    shape_df.to_csv(shape_out, index=False, sep="\t", encoding="utf-8")

    color_spec_out, color_map_out = build_phenogram_color_spec(
        df_phenogram=df_phenogram,
        output_dir_result=output_dir_result,
    )

    summary = (
        df_final.groupby(["Modality", "Phenotype", "ColorGroup", "Clock"], as_index=False)
        .agg(
            n_loci=("TopLeadSNP", "count"),
            min_p=("P-value", "min"),
            n_loci_with_physical_gene=(
                "MappedGene_PhysicalPosition",
                lambda x: (x.astype(str) != "NA").sum(),
            ),
        )
        .sort_values(["Modality", "ColorGroup", "Phenotype", "Clock"])
    )

    summary_out = os.path.join(
        output_dir_result,
        "Fuma_loci_annotation_mortality_clocks_summary.tsv",
    )
    summary.to_csv(summary_out, index=False, sep="\t", encoding="utf-8")

    df_skipped = pd.DataFrame(skipped, columns=["Clock", "Reason"])
    skipped_out = os.path.join(
        output_dir_result,
        "Fuma_loci_annotation_mortality_clocks_skipped.tsv",
    )
    df_skipped.to_csv(skipped_out, index=False, sep="\t", encoding="utf-8")

    print("\nFinished.")
    print(f"Detailed loci file: {detailed_out}")
    print(f"Main PhenoGram input: {phenogram_out}")
    print(f"Legacy lowercase PhenoGram input: {phenogram_legacy_out}")
    print(f"Group shape file: {shape_out}")
    print(f"Optional exact color spec file: {color_spec_out}")
    print(f"Color group map: {color_map_out}")
    print(f"Summary file: {summary_out}")
    print(f"Skipped log: {skipped_out}")
    print(f"Total significant loci retained: {df_final.shape[0]}")
    print(f"Total clocks with significant loci: {df_final['Clock'].nunique()}")

    print("\nRecommended PhenoGram setup:")
    print("  Input file:       PhenoGram_mortality_clocks_original_labels_same_organ_colors.tsv")
    print("  Genome:           Human GRCh37/hg19")
    print("  Phenogram color:  Grouped colors")
    print("  Group shape file: PhenoGram_mortality_clocks_group_shape.tsv")
    print("  Annotation?:      Select this if you want cytoband labels displayed")
    print("  GROUP:            modality, used for shape")
    print("  COLORGROUP:       canonical organ color group, used for color")
    print("")
    print("If you want Liver MRI and Hepatic clocks to have the EXACT same color,")
    print("choose 'File for specifying colors' and upload:")
    print(f"  {color_spec_out}")


if __name__ == "__main__":

    output_dir_fuma = (
        "/Users/hao/cubic-home/Reproducibile_paper/"
        "WholeBodyClock/mortality_clock/fuma"
    )
    output_dir_result = (
        "/Users/hao/cubic-home/Reproducibile_paper/"
        "WholeBodyClock/Result"
    )
    cyto_path = (
        "/Users/hao/cubic-home/Dataset/GRch37_cytoband/"
        "cytoBand.txt.gz"
    )

    # If running directly on CUBIC:
    # output_dir_fuma = (
    #     "/cbica/home/wenju/Reproducibile_paper/"
    #     "WholeBodyClock/mortality_clock/fuma"
    # )
    # output_dir_result = (
    #     "/cbica/home/wenju/Reproducibile_paper/"
    #     "WholeBodyClock/Result"
    # )
    # cyto_path = (
    #     "/cbica/home/wenju/Dataset/GRch37_cytoband/"
    #     "cytoBand.txt.gz"
    # )

    annotate_genomic_loci_mortality_clocks_for_phenogram(
        output_dir_fuma=output_dir_fuma,
        output_dir_result=output_dir_result,
        cyto_path=cyto_path,
        p_threshold=5e-8 / 22,
        expected_n_clocks=22,
    )