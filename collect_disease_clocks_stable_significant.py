#!/usr/bin/env python3
# ============================================================
# Collect stable significant disease-clock acceleration-z scores
#
# Python 3.7/3.8/3.9 compatible version
#
# Main fix:
#   - Output score columns now use single underscores:
#       asthma_spleen_mri_clock_acceleration_z
#     instead of:
#       asthma__spleen_mri__clock_acceleration_z
#
# Input:
#   1. Good-clock table from the R plotting script:
#      <base_dir>/all_disease_lepoch_incremental_value_scale_qc/
#        all_disease_lepoch_main_text_good_clocks.tsv
#
#      or fallback:
#        all_disease_lepoch_incremental_value_scale_qc_plot_table.tsv
#
#   2. Per-clock prediction files:
#      <base_dir>/<clock_folder>/<prefix>_predictions.tsv
#
# Output:
#   <base_dir>/all_disease_lepoch_incremental_value_scale_qc/
#      stable_significant_disease_clock_acceleration_z_wide.tsv
#      stable_significant_disease_clock_acceleration_z_long_format.tsv
#      stable_significant_disease_clock_acceleration_z_metadata.tsv
#
# Purpose:
#   Collect *_clock_acceleration_z for clocks that are:
#      significant positive M3-M1 AND stable year-scale QC
#
# Recommended run:
#   python collect_stable_significant_disease_clock_scores.py \
#     --base_dir /cbica/home/wenju/Reproducibile_paper/WholeBodyClock \
#     --no_long
# ============================================================

from __future__ import print_function

import argparse
import re
import sys
from pathlib import Path
from typing import Any, List, Optional

import pandas as pd


# ============================================================
# 1. Helpers
# ============================================================

def info(msg):
    print(msg, flush=True)


def warn(msg):
    print("WARNING: {}".format(msg), file=sys.stderr, flush=True)


def find_base_dir(user_base_dir=None):
    # type: (Optional[str]) -> Path
    candidates = []

    if user_base_dir:
        candidates.append(Path(user_base_dir))

    candidates.extend([
        Path("/cbica/home/wenju/Reproducibile_paper/WholeBodyClock"),
        Path("/gpfs/fs001/cbica/home/wenju/Reproducibile_paper/WholeBodyClock"),
        Path("/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock"),
        Path.cwd(),
    ])

    for p in candidates:
        if p.exists():
            return p.resolve()

    raise FileNotFoundError("Could not detect base_dir. Please provide --base_dir.")


def read_tsv(path):
    # type: (Path) -> pd.DataFrame
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", low_memory=False)


def clean_for_column_name(x):
    # type: (Any) -> str
    x = str(x)
    x = x.strip()
    x = x.replace(" ", "_")
    x = re.sub(r"[^A-Za-z0-9_]+", "_", x)
    x = re.sub(r"_+", "_", x)
    x = x.strip("_")
    return x.lower()


def infer_prefix_from_folder(folder):
    # type: (str) -> str
    """
    Folder examples:
      Brain_proteomics_copd_clock
      Reproductive_female_proteomics_mi_clock
      adipose_mri_asthma_clock

    Prediction prefix is usually lowercase:
      brain_proteomics_copd_clock
      reproductive_female_proteomics_mi_clock
      adipose_mri_asthma_clock
    """
    return Path(str(folder)).name.lower()


def find_prediction_file(base_dir, folder, prefix=None):
    # type: (Path, str, Optional[str]) -> Optional[Path]
    folder_path = base_dir / str(folder)

    if prefix:
        candidate = folder_path / "{}_predictions.tsv".format(prefix)
        if candidate.exists():
            return candidate

    inferred_prefix = infer_prefix_from_folder(folder)
    candidate = folder_path / "{}_predictions.tsv".format(inferred_prefix)
    if candidate.exists():
        return candidate

    if not folder_path.exists():
        return None

    hits = sorted([
        p for p in folder_path.glob("*_predictions.tsv")
        if not p.name.endswith("_test_predictions.tsv")
    ])

    if len(hits) > 0:
        return hits[0]

    return None


def detect_id_col(columns):
    # type: (List[str]) -> str
    candidates = [
        "participant_id",
        "eid",
        "subject_id",
        "IID",
        "id",
        "ID",
    ]

    for c in candidates:
        if c in columns:
            return c

    for c in columns:
        cl = c.lower()
        if cl in {"participant", "participantid", "subject", "subjectid"}:
            return c

    raise ValueError(
        "Could not detect participant ID column. "
        "Available columns include: {}".format(columns[:30])
    )


def detect_accel_z_col(columns, disease, folder, clock_id=None):
    # type: (List[str], str, str, Optional[str]) -> str
    """
    Expected examples:
      brain_mri_asthma_clock_acceleration_z
      pulmonary_proteomics_copd_clock_acceleration_z
      endocrine_metabolomics_mi_clock_acceleration_z
    """

    cols = list(columns)
    disease = str(disease).lower()
    folder_lower = infer_prefix_from_folder(folder)

    exact_candidates = [
        "{}_acceleration_z".format(folder_lower),
        "{}_clock_acceleration_z".format(folder_lower),
    ]

    for c in exact_candidates:
        if c in cols:
            return c

    accel_cols = [
        c for c in cols
        if re.search(r"clock_acceleration_z$", c, flags=re.IGNORECASE)
    ]

    disease_hits = [
        c for c in accel_cols
        if disease in c.lower()
    ]

    if len(disease_hits) == 1:
        return disease_hits[0]

    if len(disease_hits) > 1:
        folder_hits = [
            c for c in disease_hits
            if folder_lower.replace("_clock", "") in c.lower()
        ]
        if len(folder_hits) > 0:
            return folder_hits[0]
        return disease_hits[0]

    if len(accel_cols) == 1:
        return accel_cols[0]

    if len(accel_cols) > 1:
        folder_hits = [
            c for c in accel_cols
            if folder_lower.replace("_clock", "") in c.lower()
        ]
        if len(folder_hits) > 0:
            return folder_hits[0]

    raise ValueError(
        "Could not detect *_clock_acceleration_z column for "
        "disease={}, folder={}. Candidate acceleration columns: {}".format(
            disease,
            folder,
            accel_cols[:20],
        )
    )


def get_good_clock_table(base_dir, analysis_dir):
    # type: (Path, Path) -> pd.DataFrame
    """
    Prefer the main-text good-clock table. If unavailable, use the full plot table
    and filter to filled-circle class:
      Good: significant + stable scale
    """

    good_file = analysis_dir / "all_disease_lepoch_main_text_good_clocks.tsv"
    plot_file = analysis_dir / "all_disease_lepoch_incremental_value_scale_qc_plot_table.tsv"

    if good_file.exists():
        info("Reading good-clock table:\n  {}".format(good_file))
        df = read_tsv(good_file)
        df["source_good_clock_file"] = str(good_file)
        return df

    if plot_file.exists():
        info("Good-clock table not found. Filtering plot table:\n  {}".format(plot_file))
        df = read_tsv(plot_file)

        if "plot_status" not in df.columns:
            raise ValueError("{} does not contain column 'plot_status'.".format(plot_file))

        df = df[df["plot_status"] == "Good: significant + stable scale"].copy()
        df["source_good_clock_file"] = str(plot_file)
        return df

    raise FileNotFoundError(
        "Could not find either:\n"
        "  {}\n"
        "or\n"
        "  {}\n"
        "Please run the R figure/QC script first.".format(good_file, plot_file)
    )


def get_required_col(df, candidates, name):
    # type: (pd.DataFrame, List[str], str) -> str
    for c in candidates:
        if c in df.columns:
            return c

    raise ValueError(
        "Could not find required {} column. Tried: {}. Available columns: {}".format(
            name,
            candidates,
            list(df.columns),
        )
    )


def make_safe_disease_value(x):
    # type: (Any) -> str
    return clean_for_column_name(str(x).strip().lower())


def make_score_output_column(disease, clock_label):
    # type: (str, str) -> str
    """
    Fixed naming convention:
      asthma_spleen_mri_clock_acceleration_z
      dementia_brain_proteomics_clock_acceleration_z
      mi_endocrine_metabolomics_clock_acceleration_z

    Previous version used:
      asthma__spleen_mri__clock_acceleration_z
    """
    disease_clean = clean_for_column_name(disease)
    clock_clean = clean_for_column_name(clock_label)

    return "{}_{}_clock_acceleration_z".format(
        disease_clean,
        clock_clean,
    )


def drop_duplicate_ids(df, id_col, score_col, prediction_file):
    # type: (pd.DataFrame, str, str, Path) -> pd.DataFrame
    """
    If duplicate participant_id rows exist, keep the first non-missing score.
    This should rarely happen, but it prevents wide merge problems.
    """
    if not df[id_col].duplicated().any():
        return df

    warn(
        "Duplicate participant IDs detected in {}. "
        "Keeping the first non-missing score per participant.".format(prediction_file)
    )

    df["_score_is_missing_tmp"] = df[score_col].isna().astype(int)
    df = df.sort_values(by=[id_col, "_score_is_missing_tmp"])
    df = df.drop_duplicates(subset=[id_col], keep="first")
    df = df.drop(columns=["_score_is_missing_tmp"])
    return df


# ============================================================
# 2. Main collection
# ============================================================

def collect_scores(base_dir, analysis_dir, output_prefix, save_long=True, save_wide=True):
    # type: (Path, Path, str, bool, bool) -> None

    good = get_good_clock_table(base_dir, analysis_dir)

    if good.empty:
        raise ValueError("No good clocks found. Nothing to collect.")

    disease_col = get_required_col(good, ["disease"], "disease")
    folder_col = get_required_col(good, ["folder"], "folder")

    clock_label_col = "clock_label" if "clock_label" in good.columns else None
    clock_id_col = "clock_id" if "clock_id" in good.columns else None
    modality_col = "modality" if "modality" in good.columns else None
    organ_label_col = "organ_label" if "organ_label" in good.columns else None

    good = good.drop_duplicates(subset=[disease_col, folder_col]).copy()

    info("Number of stable significant disease clocks to collect: {}".format(len(good)))

    long_chunks = []
    wide_df = None
    metadata_rows = []

    for idx, row in good.iterrows():
        disease = make_safe_disease_value(row[disease_col])
        folder = str(row[folder_col])
        prefix = infer_prefix_from_folder(folder)

        prediction_file = find_prediction_file(base_dir, folder, prefix=prefix)

        if prediction_file is None:
            warn("Prediction file not found for disease={}, folder={}. Skipping.".format(disease, folder))
            metadata_rows.append({
                "disease": disease,
                "folder": folder,
                "prefix": prefix,
                "prediction_file": None,
                "status": "missing_prediction_file",
            })
            continue

        try:
            header = pd.read_csv(prediction_file, sep="\t", nrows=0)
            columns = list(header.columns)
        except Exception as e:
            warn("Failed to read header from {}: {}. Skipping.".format(prediction_file, e))
            metadata_rows.append({
                "disease": disease,
                "folder": folder,
                "prefix": prefix,
                "prediction_file": str(prediction_file),
                "status": "header_read_failed: {}".format(e),
            })
            continue

        try:
            id_col = detect_id_col(columns)
            score_col = detect_accel_z_col(
                columns=columns,
                disease=disease,
                folder=folder,
                clock_id=str(row[clock_id_col]) if clock_id_col else None,
            )
        except Exception as e:
            warn("{}: {}. Skipping.".format(prediction_file, e))
            metadata_rows.append({
                "disease": disease,
                "folder": folder,
                "prefix": prefix,
                "prediction_file": str(prediction_file),
                "status": "column_detection_failed: {}".format(e),
            })
            continue

        info("[{}] {} | {}".format(len(metadata_rows) + 1, disease, folder))
        info("    file: {}".format(prediction_file.name))
        info("    id_col: {}".format(id_col))
        info("    score_col: {}".format(score_col))

        usecols = [id_col, score_col]

        try:
            tmp = pd.read_csv(
                prediction_file,
                sep="\t",
                usecols=usecols,
                dtype={id_col: "str"},
                low_memory=False,
            )
        except Exception as e:
            warn("Failed to read {}: {}. Skipping.".format(prediction_file, e))
            metadata_rows.append({
                "disease": disease,
                "folder": folder,
                "prefix": prefix,
                "prediction_file": str(prediction_file),
                "id_col": id_col,
                "score_col": score_col,
                "status": "read_failed: {}".format(e),
            })
            continue

        tmp = tmp.rename(columns={id_col: "participant_id"})
        tmp["participant_id"] = tmp["participant_id"].astype(str)
        tmp[score_col] = pd.to_numeric(tmp[score_col], errors="coerce")
        tmp = drop_duplicate_ids(tmp, "participant_id", score_col, prediction_file)

        clock_label = str(row[clock_label_col]) if clock_label_col else folder
        clock_id = str(row[clock_id_col]) if clock_id_col else clean_for_column_name(clock_label)
        modality = str(row[modality_col]) if modality_col else ""
        organ_label = str(row[organ_label_col]) if organ_label_col else ""

        out_score_col = make_score_output_column(disease, clock_label)

        tmp_wide = tmp[["participant_id", score_col]].rename(
            columns={score_col: out_score_col}
        )

        if save_wide:
            if wide_df is None:
                wide_df = tmp_wide
            else:
                wide_df = wide_df.merge(tmp_wide, on="participant_id", how="outer")

        if save_long:
            tmp_long = tmp[["participant_id", score_col]].rename(
                columns={score_col: "clock_acceleration_z"}
            )
            tmp_long.insert(1, "disease", disease)
            tmp_long.insert(2, "clock_label", clock_label)
            tmp_long.insert(3, "clock_id", clock_id)
            tmp_long.insert(4, "modality", modality)
            tmp_long.insert(5, "organ_label", organ_label)
            tmp_long.insert(6, "folder", folder)
            tmp_long.insert(7, "score_col_original", score_col)
            tmp_long.insert(8, "score_col_wide", out_score_col)
            long_chunks.append(tmp_long)

        n_nonmissing = int(tmp_wide[out_score_col].notna().sum())
        n_total = int(tmp_wide.shape[0])

        metadata_rows.append({
            "disease": disease,
            "clock_label": clock_label,
            "clock_id": clock_id,
            "modality": modality,
            "organ_label": organ_label,
            "folder": folder,
            "prefix": prefix,
            "prediction_file": str(prediction_file),
            "id_col_original": id_col,
            "score_col_original": score_col,
            "score_col_wide": out_score_col,
            "n_rows": n_total,
            "n_nonmissing_scores": n_nonmissing,
            "status": "collected",
        })

    metadata = pd.DataFrame(metadata_rows)

    metadata_file = analysis_dir / "{}_metadata.tsv".format(output_prefix)
    metadata.to_csv(metadata_file, sep="\t", index=False)

    info("\nMetadata saved:\n  {}".format(metadata_file))

    if save_wide:
        if wide_df is None:
            warn("No wide score table was created.")
        else:
            wide_file = analysis_dir / "{}_wide.tsv".format(output_prefix)
            wide_df.to_csv(wide_file, sep="\t", index=False)
            info("Wide score matrix saved:\n  {}".format(wide_file))
            info("  rows: {:,}".format(wide_df.shape[0]))
            info("  columns: {:,}".format(wide_df.shape[1]))

    if save_long:
        if len(long_chunks) == 0:
            warn("No long-format score table was created.")
        else:
            long_df = pd.concat(long_chunks, axis=0, ignore_index=True)
            long_file = analysis_dir / "{}_long_format.tsv".format(output_prefix)
            long_df.to_csv(long_file, sep="\t", index=False)
            info("Long-format score table saved:\n  {}".format(long_file))
            info("  rows: {:,}".format(long_df.shape[0]))
            info("  columns: {:,}".format(long_df.shape[1]))

    if not metadata.empty and "status" in metadata.columns:
        collected = metadata[metadata["status"] == "collected"].shape[0]
    else:
        collected = 0

    info("\nFinished.")
    info("Collected clocks: {}".format(collected))
    info("Requested good clocks: {}".format(len(good)))


# ============================================================
# 3. CLI
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Collect *_clock_acceleration_z scores for stable significant disease L'EPOCH clocks."
    )

    parser.add_argument(
        "--base_dir",
        type=str,
        default=None,
        help="WholeBodyClock base directory.",
    )

    parser.add_argument(
        "--analysis_dir",
        type=str,
        default=None,
        help=(
            "Directory containing all_disease_lepoch_main_text_good_clocks.tsv. "
            "Default: <base_dir>/all_disease_lepoch_incremental_value_scale_qc"
        ),
    )

    parser.add_argument(
        "--output_prefix",
        type=str,
        default="stable_significant_disease_clock_acceleration_z",
        help="Output file prefix.",
    )

    parser.add_argument(
        "--no_long",
        action="store_true",
        help="Do not save long-format output.",
    )

    parser.add_argument(
        "--no_wide",
        action="store_true",
        help="Do not save wide-format output.",
    )

    args = parser.parse_args()

    base_dir = find_base_dir(args.base_dir)

    if args.analysis_dir:
        analysis_dir = Path(args.analysis_dir).resolve()
    else:
        analysis_dir = base_dir / "all_disease_lepoch_incremental_value_scale_qc"

    if not analysis_dir.exists():
        raise FileNotFoundError("analysis_dir does not exist: {}".format(analysis_dir))

    info("============================================================")
    info("Collecting stable significant disease-clock acceleration-z scores")
    info("============================================================")
    info("Base directory: {}".format(base_dir))
    info("Analysis directory: {}".format(analysis_dir))
    info("Output prefix: {}".format(args.output_prefix))
    info("============================================================")

    collect_scores(
        base_dir=base_dir,
        analysis_dir=analysis_dir,
        output_prefix=args.output_prefix,
        save_long=not args.no_long,
        save_wide=not args.no_wide,
    )


if __name__ == "__main__":
    main()