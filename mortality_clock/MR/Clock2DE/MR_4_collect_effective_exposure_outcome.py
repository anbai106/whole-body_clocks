#!/usr/bin/env python3

import os
import re
from pathlib import Path
import pandas as pd


# ============================================================
# Collect MR results for mortality Clock2DE
# Exposure: 22 mortality clocks
# Outcomes: FinnGen + PGC disease endpoints
# ============================================================

# ----------------------------
# 1. Paths
# ----------------------------
output_dir_finngen = Path(
    "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/mortality_clock/MR/Clock2DE/FinnGen"
)

output_dir_pgc = Path(
    "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/mortality_clock/MR/Clock2DE/PGC"
)

output_dir_results = Path(
    "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/Result"
)

output_dir_results.mkdir(parents=True, exist_ok=True)

finngen_manifest = Path(
    "/Users/hao/cubic-projects/MULTI/processed/FinnGen/GWAS_summary_stats/summary_stats_R9_manifest_5000_cases.tsv"
)


# ----------------------------
# 2. Expected mortality clocks
# ----------------------------
# These should match folder names under:
# Clock2DE/FinnGen/<clock_name>/MR/
# Clock2DE/PGC/<clock_name>/MR/
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


# ----------------------------
# 3. Helper functions
# ----------------------------
def clean_clock_name(clock_folder_name: str) -> str:
    """Remove mortality-clock suffix for cleaner downstream labels."""
    return re.sub(r"_mortality_clock$", "", clock_folder_name)


def discover_clock_folders(base_dirs, expected=None):
    """
    Discover clock folders that contain MR subdirectories.
    Keeps expected order first, then adds any extra discovered folders.
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
        print("[WARNING] Expected clock folders not found in FinnGen/PGC directories:")
        for x in missing:
            print(f"  - {x}")

    if extra:
        print("[INFO] Extra discovered clock folders not in expected list:")
        for x in extra:
            print(f"  - {x}")

    return ordered + extra


def read_mr_file(tsv_path: Path, source: str, clock_folder: str, outcome_code: str):
    """Read one MR result file and add metadata."""
    try:
        df = pd.read_csv(tsv_path, sep="\t")
    except Exception as e:
        print(f"[WARNING] Failed to read {tsv_path}: {e}")
        return None

    if df.empty:
        print(f"[WARNING] Empty file: {tsv_path}")
        return None

    # Drop 2SampleMR internal IDs if present
    for col in ["id.exposure", "id.outcome"]:
        if col in df.columns:
            df = df.drop(columns=[col])

    df["target_source"] = source
    df["outcome_code"] = outcome_code
    df["clock_folder"] = clock_folder
    df["clock_id"] = clean_clock_name(clock_folder)
    df["file_path"] = str(tsv_path)

    # Robust exposure/outcome fields
    if "exposure" not in df.columns:
        df["exposure"] = df["clock_id"]

    if "outcome" not in df.columns:
        df["outcome"] = outcome_code

    # Make sure pval is numeric
    if "pval" in df.columns:
        df["pval"] = pd.to_numeric(df["pval"], errors="coerce")
    else:
        df["pval"] = pd.NA
        print(f"[WARNING] No pval column in {tsv_path}")

    # Standardize beta/se/or columns if present
    for col in ["b", "se", "or", "or_lci95", "or_uci95", "nsnp"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


def collect_pgc_results(clock_folders):
    """Collect PGC MR results."""
    dfs = []

    for clock in clock_folders:
        mr_dir = output_dir_pgc / clock / "MR"

        if not mr_dir.exists():
            print(f"[INFO] PGC MR directory missing for {clock}: {mr_dir}")
            continue

        for pgc in pgc_outcomes:
            expected_file = mr_dir / f"MR_{clock}_2_{pgc}_OR.tsv"

            if expected_file.exists():
                df = read_mr_file(expected_file, "PGC", clock, pgc)
                if df is not None:
                    dfs.append(df)
            else:
                # Fallback: search flexibly in case naming differs slightly
                matches = list(mr_dir.glob(f"MR_*_2_{pgc}_OR.tsv"))
                if len(matches) == 0:
                    continue

                for tsv in matches:
                    df = read_mr_file(tsv, "PGC", clock, pgc)
                    if df is not None:
                        dfs.append(df)

    return dfs


def load_finngen_phenocodes():
    """Load FinnGen phenocodes from manifest."""
    if not finngen_manifest.exists():
        print(f"[WARNING] FinnGen manifest not found: {finngen_manifest}")
        return []

    manifest = pd.read_csv(finngen_manifest, sep="\t")

    if "phenocode" not in manifest.columns:
        raise ValueError(f"Manifest does not contain a 'phenocode' column: {finngen_manifest}")

    phenos = (
        manifest["phenocode"]
        .dropna()
        .astype(str)
        .drop_duplicates()
        .tolist()
    )

    print(f"[INFO] Loaded {len(phenos)} FinnGen phenocodes from manifest.")

    return phenos


def parse_finngen_outcome_from_filename(path: Path, clock: str):
    """
    Parse FinnGen outcome from filename:
    MR_<clock>_2_<phenocode>_OR.tsv
    """
    prefix = f"MR_{clock}_2_"
    suffix = "_OR.tsv"

    name = path.name

    if name.startswith(prefix) and name.endswith(suffix):
        return name[len(prefix):-len(suffix)]

    # More flexible fallback
    m = re.match(r"^MR_.+?_2_(.+)_OR\.tsv$", name)
    if m:
        return m.group(1)

    return None


def collect_finngen_results(clock_folders, finngen_phenos):
    """Collect FinnGen MR results."""
    dfs = []

    for clock in clock_folders:
        mr_dir = output_dir_finngen / clock / "MR"

        if not mr_dir.exists():
            print(f"[INFO] FinnGen MR directory missing for {clock}: {mr_dir}")
            continue

        found_files = set()

        # First, use manifest phenocodes to look for exact expected files
        for pheno in finngen_phenos:
            expected_file = mr_dir / f"MR_{clock}_2_{pheno}_OR.tsv"

            if expected_file.exists():
                found_files.add(expected_file)

        # Second, glob all MR files to catch outcomes not in the manifest
        for tsv in mr_dir.glob("MR_*_2_*_OR.tsv"):
            found_files.add(tsv)

        for tsv in sorted(found_files):
            outcome_code = parse_finngen_outcome_from_filename(tsv, clock)

            if outcome_code is None:
                print(f"[WARNING] Could not parse FinnGen outcome from filename: {tsv}")
                outcome_code = "UNKNOWN"

            df = read_mr_file(tsv, "FinnGen", clock, outcome_code)
            if df is not None:
                dfs.append(df)

    return dfs


# ----------------------------
# 4. Collect all results
# ----------------------------
clock_folders = discover_clock_folders(
    base_dirs=[output_dir_finngen, output_dir_pgc],
    expected=expected_clocks
)

print(f"[INFO] Clock folders to scan: {len(clock_folders)}")
for c in clock_folders:
    print(f"  - {c}")

finngen_phenos = load_finngen_phenocodes()

all_dfs = []
all_dfs.extend(collect_pgc_results(clock_folders))
all_dfs.extend(collect_finngen_results(clock_folders, finngen_phenos))

if len(all_dfs) == 0:
    raise RuntimeError("No MR result files were found. Please check the input directories and file naming.")

df_final = pd.concat(all_dfs, ignore_index=True)

# Drop exact duplicated rows if any arose from exact + glob fallback
dedup_cols = [
    "target_source",
    "clock_id",
    "outcome_code",
    "method",
    "exposure",
    "outcome",
]

dedup_cols = [x for x in dedup_cols if x in df_final.columns]
df_final = df_final.drop_duplicates(subset=dedup_cols, keep="first").reset_index(drop=True)


# ----------------------------
# 5. Multiple testing thresholds
# ----------------------------
# Primary unit is clock x disease endpoint.
# Use IVW as the primary MR method for significance calls.
primary_method = "Inverse variance weighted"

df_final_ivw = df_final.loc[df_final["method"].eq(primary_method)].copy()

if df_final_ivw.empty:
    raise RuntimeError(
        f"No rows with method == '{primary_method}'. "
        "Please check method names in the MR result files."
    )

n_disease = df_final_ivw["outcome_code"].nunique()
n_clock = df_final_ivw["clock_id"].nunique()
n_tests = n_disease * n_clock

print(f"Unique clocks tested, IVW: {n_clock}")
print(f"Unique disease endpoints tested, IVW: {n_disease}")
print(f"Total clock x disease tests, IVW: {n_tests}")

df_final["P_bon_n_disease"] = 0.05 / n_disease
df_final["P_bon_n_clock"] = 0.05 / n_clock
df_final["P_bon_n_clock_de"] = 0.05 / n_tests

df_final_ivw["P_bon_n_disease"] = 0.05 / n_disease
df_final_ivw["P_bon_n_clock"] = 0.05 / n_clock
df_final_ivw["P_bon_n_clock_de"] = 0.05 / n_tests

# Significance columns
for df in [df_final, df_final_ivw]:
    df["sig_nominal"] = df["pval"] <= 0.05
    df["sig_by_disease"] = df["pval"] <= df["P_bon_n_disease"]
    df["sig_by_clock"] = df["pval"] <= df["P_bon_n_clock"]
    df["sig_by_clock_de"] = df["pval"] <= df["P_bon_n_clock_de"]


# ----------------------------
# 6. Significant result tables
# ----------------------------
# Recommended primary significant tables should use IVW only.
df_ivw_sig_by_de = df_final_ivw.loc[df_final_ivw["sig_by_disease"]].copy()
df_ivw_sig_by_clock = df_final_ivw.loc[df_final_ivw["sig_by_clock"]].copy()
df_ivw_sig_by_clock_de = df_final_ivw.loc[df_final_ivw["sig_by_clock_de"]].copy()

# Optional all-method significant tables, useful for sensitivity checking
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
        n_outcomes=("outcome_code", "nunique"),
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
        n_outcomes=("outcome_code", "nunique"),
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
        "n_clock": [n_clock],
        "n_disease": [n_disease],
        "n_tests": [n_tests],
        "P_nominal": [0.05],
        "P_bon_n_disease": [0.05 / n_disease],
        "P_bon_n_clock": [0.05 / n_clock],
        "P_bon_n_clock_de": [0.05 / n_tests],
        "primary_method": [primary_method],
    }
)

print("\n[INFO] IVW summary by source")
print(summary_by_source)

print("\n[INFO] Multiple-testing thresholds")
print(thresholds)


# ----------------------------
# 8. Save outputs
# ----------------------------
prefix = "2SampleMR_MortalityClock2DE"

df_final.to_csv(
    output_dir_results / f"{prefix}_all_methods.tsv",
    index=False,
    sep="\t",
    encoding="utf-8"
)

df_final_ivw.to_csv(
    output_dir_results / f"{prefix}_ivw.tsv",
    index=False,
    sep="\t",
    encoding="utf-8"
)

df_ivw_sig_by_de.to_csv(
    output_dir_results / f"{prefix}_ivw_sig_by_DE.tsv",
    index=False,
    sep="\t",
    encoding="utf-8"
)

df_ivw_sig_by_clock.to_csv(
    output_dir_results / f"{prefix}_ivw_sig_by_clock.tsv",
    index=False,
    sep="\t",
    encoding="utf-8"
)

df_ivw_sig_by_clock_de.to_csv(
    output_dir_results / f"{prefix}_ivw_sig_by_clock_DE.tsv",
    index=False,
    sep="\t",
    encoding="utf-8"
)

df_all_sig_by_de.to_csv(
    output_dir_results / f"{prefix}_all_methods_sig_by_DE.tsv",
    index=False,
    sep="\t",
    encoding="utf-8"
)

df_all_sig_by_clock.to_csv(
    output_dir_results / f"{prefix}_all_methods_sig_by_clock.tsv",
    index=False,
    sep="\t",
    encoding="utf-8"
)

df_all_sig_by_clock_de.to_csv(
    output_dir_results / f"{prefix}_all_methods_sig_by_clock_DE.tsv",
    index=False,
    sep="\t",
    encoding="utf-8"
)

summary_by_source.to_csv(
    output_dir_results / f"{prefix}_ivw_summary_by_source.tsv",
    index=False,
    sep="\t",
    encoding="utf-8"
)

summary_by_clock.to_csv(
    output_dir_results / f"{prefix}_ivw_summary_by_clock.tsv",
    index=False,
    sep="\t",
    encoding="utf-8"
)

thresholds.to_csv(
    output_dir_results / f"{prefix}_multiple_testing_thresholds.tsv",
    index=False,
    sep="\t",
    encoding="utf-8"
)

print("\n[INFO] Saved outputs to:")
print(output_dir_results)

print("\nSTOP here...")