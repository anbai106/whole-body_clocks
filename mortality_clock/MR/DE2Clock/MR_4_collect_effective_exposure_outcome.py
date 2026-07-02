#!/usr/bin/env python3

import os
import re
from pathlib import Path
import pandas as pd


# ============================================================
# Collect MR results for DE2Clock
#
# Exposure: disease endpoints, FinnGen or PGC
# Outcome: 22 mortality clocks
#
# Expected structure:
#   DE2Clock/FinnGen/<clock_folder>/MR/MR_<FinnGen_phenocode>_2_<clock_folder>_OR.tsv
#   DE2Clock/PGC/<clock_folder>/MR/MR_<PGC_code>_2_<clock_folder>_OR.tsv
# ============================================================

# ----------------------------
# 1. Paths
# ----------------------------
output_dir_results = Path(
    "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/Result"
)

output_dir_finngen = Path(
    "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/mortality_clock/MR/DE2Clock/FinnGen"
)

output_dir_pgc = Path(
    "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/mortality_clock/MR/DE2Clock/PGC"
)

finngen_manifest = Path(
    "/Users/hao/cubic-projects/MULTI/processed/FinnGen/GWAS_summary_stats/summary_stats_R9_manifest_5000_cases.tsv"
)

output_dir_results.mkdir(parents=True, exist_ok=True)


# ----------------------------
# 2. Expected mortality-clock folders
# ----------------------------
expected_clocks = [
    "adipose_mri_mortality_clock",
    "brain_mri_mortality_clock",
    "Brain_proteomics_mortality_clock",
    "Digestive_metabolomics_mortality_clock",
    "Endocrine_metabolomics_mortality_clock",
    "Endocrine_proteomics_mortality_clock",
    "Eye_proteomics_mortality_clock",
    "heart_mri_mortality_clock",
    "Heart_proteomics_mortality_clock",
    "Hepatic_metabolomics_mortality_clock",
    "Hepatic_proteomics_mortality_clock",
    "Immune_metabolomics_mortality_clock",
    "Immune_proteomics_mortality_clock",
    "kidney_mri_mortality_clock",
    "liver_mri_mortality_clock",
    "pancreas_mri_mortality_clock",
    "Pulmonary_proteomics_mortality_clock",
    "Renal_proteomics_mortality_clock",
    "Reproductive_female_proteomics_mortality_clock",
    "Reproductive_male_proteomics_mortality_clock",
    "Skin_proteomics_mortality_clock",
    "spleen_mri_mortality_clock",
]

pgc_outcomes = ["AD", "ADHD", "BIP", "SCZ"]
primary_method = "Inverse variance weighted"


# ----------------------------
# 3. Helper functions
# ----------------------------
def clean_clock_name(clock_folder_name: str) -> str:
    """Remove mortality-clock suffix for a cleaner clock ID."""
    return re.sub(r"_mortality_clock$", "", clock_folder_name)


def discover_clock_folders(base_dirs, expected=None):
    """
    Discover clock folders that contain MR subdirectories.
    Keeps expected clock order first, then adds extra discovered folders.
    """
    discovered = set()

    for base_dir in base_dirs:
        if not base_dir.exists():
            print(f"[WARNING] Base directory does not exist: {base_dir}")
            continue

        for p in base_dir.iterdir():
            if p.is_dir() and (p / "MR").is_dir():
                discovered.add(p.name)

    if expected is None:
        return sorted(discovered)

    ordered = [x for x in expected if x in discovered]
    extra = sorted(discovered.difference(expected))
    missing = [x for x in expected if x not in discovered]

    if missing:
        print("[WARNING] Expected clock folders not found:")
        for x in missing:
            print(f"  - {x}")

    if extra:
        print("[INFO] Extra discovered clock folders not in expected list:")
        for x in extra:
            print(f"  - {x}")

    return ordered + extra


def load_finngen_phenocodes():
    """Load FinnGen phenocodes from manifest."""
    if not finngen_manifest.exists():
        print(f"[WARNING] FinnGen manifest not found: {finngen_manifest}")
        return []

    manifest = pd.read_csv(finngen_manifest, sep="\t")

    if "phenocode" not in manifest.columns:
        raise ValueError(
            f"Manifest does not contain a 'phenocode' column: {finngen_manifest}"
        )

    phenos = (
        manifest["phenocode"]
        .dropna()
        .astype(str)
        .drop_duplicates()
        .tolist()
    )

    print(f"[INFO] Loaded {len(phenos)} FinnGen phenocodes from manifest.")

    return phenos


def parse_de2clock_filename(path: Path, clock_folder: str):
    """
    Parse disease code from DE2Clock filename:
      MR_<disease_code>_2_<clock_folder>_OR.tsv

    Because FinnGen phenocodes can contain underscores, parse using
    the exact suffix _2_<clock_folder>_OR.tsv.
    """
    name = path.name
    prefix = "MR_"
    suffix = f"_2_{clock_folder}_OR.tsv"

    if name.startswith(prefix) and name.endswith(suffix):
        return name[len(prefix):-len(suffix)]

    # Flexible fallback
    m = re.match(r"^MR_(.+)_2_.+_OR\.tsv$", name)
    if m:
        return m.group(1)

    return None


def read_mr_file(tsv_path: Path, source: str, clock_folder: str, disease_code: str):
    """Read one MR result file and add metadata."""
    try:
        df = pd.read_csv(tsv_path, sep="\t")
    except Exception as e:
        print(f"[WARNING] Failed to read {tsv_path}: {e}")
        return None

    if df.empty:
        print(f"[WARNING] Empty file: {tsv_path}")
        return None

    # Preserve original 2SampleMR labels before standardizing
    if "exposure" in df.columns:
        df["exposure_original"] = df["exposure"]
    if "outcome" in df.columns:
        df["outcome_original"] = df["outcome"]

    # Drop 2SampleMR internal IDs if present
    for col in ["id.exposure", "id.outcome"]:
        if col in df.columns:
            df = df.drop(columns=[col])

    clock_id = clean_clock_name(clock_folder)

    # Add standardized metadata
    df["MR_direction"] = "DE2Clock"
    df["target_source"] = source
    df["disease_code"] = disease_code
    df["clock_folder"] = clock_folder
    df["clock_id"] = clock_id
    df["file_path"] = str(tsv_path)

    # Standardize exposure/outcome for downstream counting
    # Exposure = disease endpoint; outcome = mortality clock
    df["exposure"] = disease_code
    df["outcome"] = clock_id

    # Numeric conversion
    if "pval" in df.columns:
        df["pval"] = pd.to_numeric(df["pval"], errors="coerce")
    else:
        df["pval"] = pd.NA
        print(f"[WARNING] No pval column in {tsv_path}")

    for col in ["b", "se", "or", "or_lci95", "or_uci95", "nsnp"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


def collect_finngen_results(clock_folders, finngen_phenos):
    """Collect FinnGen DE2Clock MR results."""
    dfs = []

    for clock in clock_folders:
        mr_dir = output_dir_finngen / clock / "MR"

        if not mr_dir.exists():
            print(f"[INFO] FinnGen MR directory missing for {clock}: {mr_dir}")
            continue

        found_files = set()

        # First use manifest phenocodes for exact expected files
        for pheno in finngen_phenos:
            expected_file = mr_dir / f"MR_{pheno}_2_{clock}_OR.tsv"
            if expected_file.exists():
                found_files.add(expected_file)

        # Then glob all MR files to catch any files not in manifest
        for tsv in mr_dir.glob("MR_*_2_*_OR.tsv"):
            found_files.add(tsv)

        for tsv in sorted(found_files):
            disease_code = parse_de2clock_filename(tsv, clock)

            if disease_code is None:
                print(f"[WARNING] Could not parse disease code from: {tsv}")
                disease_code = "UNKNOWN"

            df = read_mr_file(
                tsv_path=tsv,
                source="FinnGen",
                clock_folder=clock,
                disease_code=disease_code,
            )

            if df is not None:
                dfs.append(df)

    return dfs


def collect_pgc_results(clock_folders):
    """Collect PGC DE2Clock MR results."""
    dfs = []

    for clock in clock_folders:
        mr_dir = output_dir_pgc / clock / "MR"

        if not mr_dir.exists():
            print(f"[INFO] PGC MR directory missing for {clock}: {mr_dir}")
            continue

        found_files = set()

        # Exact expected PGC files
        for pgc in pgc_outcomes:
            expected_file = mr_dir / f"MR_{pgc}_2_{clock}_OR.tsv"
            if expected_file.exists():
                found_files.add(expected_file)

        # Flexible fallback
        for tsv in mr_dir.glob("MR_*_2_*_OR.tsv"):
            found_files.add(tsv)

        for tsv in sorted(found_files):
            disease_code = parse_de2clock_filename(tsv, clock)

            if disease_code is None:
                print(f"[WARNING] Could not parse PGC code from: {tsv}")
                disease_code = "UNKNOWN"

            df = read_mr_file(
                tsv_path=tsv,
                source="PGC",
                clock_folder=clock,
                disease_code=disease_code,
            )

            if df is not None:
                dfs.append(df)

    return dfs


# ----------------------------
# 4. Collect all DE2Clock results
# ----------------------------
clock_folders = discover_clock_folders(
    base_dirs=[output_dir_finngen, output_dir_pgc],
    expected=expected_clocks,
)

print(f"[INFO] Clock folders to scan: {len(clock_folders)}")
for c in clock_folders:
    print(f"  - {c}")

finngen_phenos = load_finngen_phenocodes()

all_dfs = []
all_dfs.extend(collect_finngen_results(clock_folders, finngen_phenos))
all_dfs.extend(collect_pgc_results(clock_folders))

if len(all_dfs) == 0:
    raise RuntimeError(
        "No DE2Clock MR result files were found. "
        "Please check the input directories and file naming."
    )

df_final = pd.concat(all_dfs, ignore_index=True)

# Drop exact duplicated rows from exact + glob collection
dedup_cols = [
    "target_source",
    "disease_code",
    "clock_id",
    "method",
    "exposure",
    "outcome",
]

dedup_cols = [x for x in dedup_cols if x in df_final.columns]

df_final = (
    df_final
    .drop_duplicates(subset=dedup_cols, keep="first")
    .reset_index(drop=True)
)


# ----------------------------
# 5. Primary IVW results and multiple testing
# ----------------------------
df_final_ivw = df_final.loc[
    df_final["method"].eq(primary_method)
].copy()

if df_final_ivw.empty:
    available_methods = sorted(df_final["method"].dropna().unique().tolist())
    raise RuntimeError(
        f"No rows with method == '{primary_method}'. "
        f"Available methods are: {available_methods}"
    )

n_disease = df_final_ivw["disease_code"].nunique()
n_clock = df_final_ivw["clock_id"].nunique()
n_tests = n_disease * n_clock

print(f"Unique mortality clocks tested, IVW: {n_clock}")
print(f"Unique disease endpoints tested, IVW: {n_disease}")
print(f"Total disease x clock tests, IVW: {n_tests}")

# Add thresholds to both all-method and IVW tables
for df in [df_final, df_final_ivw]:
    df["P_bon_n_disease"] = 0.05 / n_disease
    df["P_bon_n_clock"] = 0.05 / n_clock
    df["P_bon_n_clock_de"] = 0.05 / n_tests

    df["sig_nominal"] = df["pval"] <= 0.05
    df["sig_by_disease"] = df["pval"] <= df["P_bon_n_disease"]
    df["sig_by_clock"] = df["pval"] <= df["P_bon_n_clock"]
    df["sig_by_clock_de"] = df["pval"] <= df["P_bon_n_clock_de"]


# ----------------------------
# 6. Significant result tables
# ----------------------------
# Recommended primary significance calls: IVW only
df_ivw_sig_by_de = df_final_ivw.loc[df_final_ivw["sig_by_disease"]].copy()
df_ivw_sig_by_clock = df_final_ivw.loc[df_final_ivw["sig_by_clock"]].copy()
df_ivw_sig_by_clock_de = df_final_ivw.loc[df_final_ivw["sig_by_clock_de"]].copy()

# Optional all-method significant tables
df_all_sig_by_de = df_final.loc[df_final["sig_by_disease"]].copy()
df_all_sig_by_clock = df_final.loc[df_final["sig_by_clock"]].copy()
df_all_sig_by_clock_de = df_final.loc[df_final["sig_by_clock_de"]].copy()


# ----------------------------
# 7. Diagnostics
# ----------------------------
summary_by_source = (
    df_final_ivw
    .groupby("target_source", dropna=False)
    .agg(
        n_rows=("pval", "size"),
        n_clocks=("clock_id", "nunique"),
        n_diseases=("disease_code", "nunique"),
        n_nominal=("sig_nominal", "sum"),
        n_sig_by_disease=("sig_by_disease", "sum"),
        n_sig_by_clock=("sig_by_clock", "sum"),
        n_sig_by_clock_de=("sig_by_clock_de", "sum"),
        min_pval=("pval", "min"),
    )
    .reset_index()
)

summary_by_clock = (
    df_final_ivw
    .groupby(["clock_id", "target_source"], dropna=False)
    .agg(
        n_diseases=("disease_code", "nunique"),
        n_nominal=("sig_nominal", "sum"),
        n_sig_by_disease=("sig_by_disease", "sum"),
        n_sig_by_clock=("sig_by_clock", "sum"),
        n_sig_by_clock_de=("sig_by_clock_de", "sum"),
        min_pval=("pval", "min"),
    )
    .reset_index()
)

summary_by_disease = (
    df_final_ivw
    .groupby(["target_source", "disease_code"], dropna=False)
    .agg(
        n_clocks=("clock_id", "nunique"),
        n_nominal=("sig_nominal", "sum"),
        n_sig_by_disease=("sig_by_disease", "sum"),
        n_sig_by_clock=("sig_by_clock", "sum"),
        n_sig_by_clock_de=("sig_by_clock_de", "sum"),
        min_pval=("pval", "min"),
    )
    .reset_index()
)

thresholds = pd.DataFrame(
    {
        "MR_direction": ["DE2Clock"],
        "primary_method": [primary_method],
        "n_clock": [n_clock],
        "n_disease": [n_disease],
        "n_tests": [n_tests],
        "P_nominal": [0.05],
        "P_bon_n_disease": [0.05 / n_disease],
        "P_bon_n_clock": [0.05 / n_clock],
        "P_bon_n_clock_de": [0.05 / n_tests],
    }
)

print("\n[INFO] IVW summary by source")
print(summary_by_source)

print("\n[INFO] Multiple-testing thresholds")
print(thresholds)


# ----------------------------
# 8. Save outputs
# ----------------------------
prefix = "2SampleMR_DE2clock"

df_final.to_csv(
    output_dir_results / f"{prefix}_all_methods.tsv",
    index=False,
    sep="\t",
    encoding="utf-8",
)

df_final_ivw.to_csv(
    output_dir_results / f"{prefix}_ivw.tsv",
    index=False,
    sep="\t",
    encoding="utf-8",
)

df_ivw_sig_by_de.to_csv(
    output_dir_results / f"{prefix}_ivw_sig_by_DE.tsv",
    index=False,
    sep="\t",
    encoding="utf-8",
)

df_ivw_sig_by_clock.to_csv(
    output_dir_results / f"{prefix}_ivw_sig_by_clock.tsv",
    index=False,
    sep="\t",
    encoding="utf-8",
)

df_ivw_sig_by_clock_de.to_csv(
    output_dir_results / f"{prefix}_ivw_sig_by_clock_DE.tsv",
    index=False,
    sep="\t",
    encoding="utf-8",
)

df_all_sig_by_de.to_csv(
    output_dir_results / f"{prefix}_all_methods_sig_by_DE.tsv",
    index=False,
    sep="\t",
    encoding="utf-8",
)

df_all_sig_by_clock.to_csv(
    output_dir_results / f"{prefix}_all_methods_sig_by_clock.tsv",
    index=False,
    sep="\t",
    encoding="utf-8",
)

df_all_sig_by_clock_de.to_csv(
    output_dir_results / f"{prefix}_all_methods_sig_by_clock_DE.tsv",
    index=False,
    sep="\t",
    encoding="utf-8",
)

summary_by_source.to_csv(
    output_dir_results / f"{prefix}_ivw_summary_by_source.tsv",
    index=False,
    sep="\t",
    encoding="utf-8",
)

summary_by_clock.to_csv(
    output_dir_results / f"{prefix}_ivw_summary_by_clock.tsv",
    index=False,
    sep="\t",
    encoding="utf-8",
)

summary_by_disease.to_csv(
    output_dir_results / f"{prefix}_ivw_summary_by_disease.tsv",
    index=False,
    sep="\t",
    encoding="utf-8",
)

thresholds.to_csv(
    output_dir_results / f"{prefix}_multiple_testing_thresholds.tsv",
    index=False,
    sep="\t",
    encoding="utf-8",
)

print("\n[INFO] Saved outputs to:")
print(output_dir_results)

print("\nSTOP here...")