#!/usr/bin/env python3
# ============================================================
# Prepare complete-case longitudinal abdominal imaging data
# for applying pre-trained abdominal MRI mortality clocks to
# UKB instance 3.
#
# Goal:
#   For each abdominal organ clock, create a TSV with the same
#   feature-column names as the baseline instance-2 training file,
#   but filled with instance-3 imaging values from abdo_imaging.csv.
#
# Example:
#   Baseline training column:
#     Visceral_fat_volume_21085-2.0
#
#   Raw longitudinal column:
#     21085-3.0
#
#   Output longitudinal column:
#     Visceral_fat_volume_21085-2.0
#
# Important:
#   This script keeps only participants with ALL required organ
#   features non-missing at instance 3. If no complete-case
#   participants remain for an organ, the organ is skipped with
#   a warning.
# ============================================================

import os
import re
import json
import warnings
from pathlib import Path

import pandas as pd


# -----------------------------
# 1. Input/output paths
# -----------------------------

RAW_ABDO_FILE = (
    "/Users/hao/cubic-home/Reproducibile_paper/AbdoImaging/"
    "Data_raw_before_Huizi/Multiorgan_imaging_BAG/abdo_imaging.csv"
)

BASELINE_SPLIT_DIR = (
    "/Users/hao/cubic-home/Reproducibile_paper/AbdoImaging/"
    "Data_split_by_Huizi/MLNI_dataset_splits"
)

OUTPUT_DIR = (
    "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/"
    "mortality_clock/longitudinal/imaging/data"
)

QC_DIR = os.path.join(OUTPUT_DIR, "qc")

ORGAN_LIST = ["adipose", "kidney", "liver", "pancreas", "spleen"]

TARGET_INSTANCE = "3"
MODEL_INSTANCE = "2"
OUTPUT_SESSION_ID = "ses-M3"


# -----------------------------
# 2. Helper functions
# -----------------------------

def ensure_dir(path: str) -> None:
    Path(path).mkdir(parents=True, exist_ok=True)


def extract_ukb_field_suffix(colname: str):
    """
    Extract UKB field-instance-array suffix from a baseline MLNI feature name.

    Examples:
      Visceral_fat_volume_21085-2.0
      Posterior_thigh_muscle_fat_infiltration_(MFI)_(left)_23355-2.0
      Some_feature_23359-2.1

    Returns:
      field_id, instance, array
      e.g., ("21085", "2", "0")
    """
    m = re.search(r"(\d+)-(\d+)\.(\d+)$", str(colname))
    if m is None:
        return None

    return {
        "field_id": m.group(1),
        "instance": m.group(2),
        "array": m.group(3),
        "full_suffix": f"{m.group(1)}-{m.group(2)}.{m.group(3)}",
    }


def make_target_raw_col(field_id: str, array: str, target_instance: str = TARGET_INSTANCE) -> str:
    """
    Convert baseline field/array to raw longitudinal field name.

    Example:
      field_id=21085, array=0, target_instance=3
      -> 21085-3.0
    """
    return f"{field_id}-{target_instance}.{array}"


def read_abdominal_raw_longitudinal(raw_file: str, target_instance: str = TARGET_INSTANCE) -> pd.DataFrame:
    """
    Read raw abdo_imaging.csv and keep eid plus all target-instance columns.

    This keeps all arrays, not only -3.0. For example, it keeps -3.1 and -3.2.
    """
    print(f"Reading raw abdominal imaging file:\n  {raw_file}")

    df = pd.read_csv(raw_file)

    if "eid" not in df.columns:
        raise ValueError("Raw abdominal imaging file must contain column 'eid'.")

    target_pattern = re.compile(rf"^\d+-{re.escape(target_instance)}\.\d+$")
    keep_cols = ["eid"] + [c for c in df.columns if target_pattern.match(str(c))]

    if len(keep_cols) == 1:
        raise ValueError(
            f"No instance-{target_instance} columns found in raw abdominal imaging file."
        )

    df = df[keep_cols].copy()
    df = df.rename(columns={"eid": "participant_id"})

    df["participant_id"] = pd.to_numeric(df["participant_id"], errors="coerce")
    df = df.dropna(subset=["participant_id"]).copy()
    df["participant_id"] = df["participant_id"].astype(int)

    duplicated = int(df["participant_id"].duplicated().sum())
    if duplicated > 0:
        warnings.warn(
            f"Found {duplicated} duplicated participant_id rows in raw abdominal file. "
            "Keeping the first occurrence."
        )
        df = df.drop_duplicates("participant_id", keep="first")

    print("Raw longitudinal abdominal data:")
    print(f"  N participants = {df.shape[0]}")
    print(f"  N instance-{target_instance} imaging columns = {df.shape[1] - 1}")

    return df


def read_baseline_organ_template(organ: str, baseline_split_dir: str = BASELINE_SPLIT_DIR) -> pd.DataFrame:
    """
    Read the baseline MLNI organ split file used to train the clock.

    Expected columns:
      participant_id, session_id, diagnosis, feature1, feature2, ...
    """
    path = os.path.join(baseline_split_dir, f"{organ}_train_val_test.tsv")

    if not os.path.exists(path):
        raise FileNotFoundError(f"Baseline organ file not found: {path}")

    print(f"\nReading baseline template for {organ}:\n  {path}")
    df = pd.read_csv(path, sep="\t")

    required = ["participant_id"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{organ} baseline file is missing required columns: {missing}")

    return df


def get_baseline_feature_columns(baseline_df: pd.DataFrame, organ: str):
    """
    Robustly infer imaging feature columns from the baseline training file.

    Works for either:
      participant_id, session_id, diagnosis, feature1, feature2, ...
    or:
      participant_id, feature1, feature2, ...
    """
    non_feature_cols = {
        "participant_id",
        "eid",
        "session_id",
        "diagnosis",
        "organ_source_file",
        "organ_source_order",
        "organ_source_row",
    }

    candidate_cols = [c for c in baseline_df.columns if c not in non_feature_cols]

    feature_cols = []
    skipped_cols = []

    for c in candidate_cols:
        if extract_ukb_field_suffix(c) is not None:
            feature_cols.append(c)
        else:
            skipped_cols.append(c)

    if skipped_cols:
        print(f"  {organ}: skipped {len(skipped_cols)} non-feature columns.")
        for c in skipped_cols[:20]:
            print(f"    skipped: {c}")
        if len(skipped_cols) > 20:
            print("    ...")

    if len(feature_cols) == 0:
        raise ValueError(
            f"No baseline imaging features found for organ={organ}. "
            "Expected feature names ending with a suffix like 21085-2.0."
        )

    return feature_cols


def build_feature_mapping(baseline_feature_cols, raw_long_cols, organ: str):
    """
    Build mapping from baseline model-expected feature names to raw longitudinal columns.

    Output keeps the baseline feature name but uses values from the target-instance column.
    """
    raw_long_set = set(raw_long_cols)

    rows = []
    missing_rows = []
    malformed_rows = []

    for baseline_col in baseline_feature_cols:
        parsed = extract_ukb_field_suffix(baseline_col)

        if parsed is None:
            malformed_rows.append({
                "organ": organ,
                "baseline_feature_col": baseline_col,
                "reason": "Could not extract trailing UKB suffix like 21085-2.0",
            })
            continue

        raw_target_col = make_target_raw_col(
            field_id=parsed["field_id"],
            array=parsed["array"],
            target_instance=TARGET_INSTANCE,
        )

        row = {
            "organ": organ,
            "baseline_feature_col": baseline_col,
            "baseline_suffix": parsed["full_suffix"],
            "field_id": parsed["field_id"],
            "baseline_instance": parsed["instance"],
            "array": parsed["array"],
            "raw_target_col": raw_target_col,
            "output_feature_col": baseline_col,
            "found_in_raw_longitudinal": raw_target_col in raw_long_set,
        }

        rows.append(row)

        if raw_target_col not in raw_long_set:
            missing_rows.append(row)

    mapping = pd.DataFrame(rows)
    missing = pd.DataFrame(missing_rows)
    malformed = pd.DataFrame(malformed_rows)

    return mapping, missing, malformed


def write_skip_summary(
    organ: str,
    reason: str,
    qc_dir: str,
    n_expected: int = None,
    n_mapped: int = None,
    n_missing: int = None,
    n_complete_case: int = None,
):
    """
    Write a JSON summary when an organ is skipped.
    """
    summary = {
        "organ": organ,
        "target_instance": TARGET_INSTANCE,
        "model_instance": MODEL_INSTANCE,
        "status": "skipped",
        "reason": reason,
        "n_expected_model_features": n_expected,
        "n_mapped_instance3_features": n_mapped,
        "n_missing_instance3_features": n_missing,
        "n_complete_case_participants": n_complete_case,
    }

    out_json = os.path.join(qc_dir, f"{organ}_instance3_prepare_summary.json")
    with open(out_json, "w") as f:
        json.dump(summary, f, indent=2)

    warnings.warn(f"Skipping organ={organ}: {reason}")
    return summary


def create_organ_longitudinal_file(
    organ: str,
    raw_long: pd.DataFrame,
    output_dir: str = OUTPUT_DIR,
    qc_dir: str = QC_DIR,
):
    """
    Create one organ-specific longitudinal complete-case input file.

    Output format:
      participant_id, session_id, diagnosis, baseline-style feature columns...

    Complete-case rule:
      Keep only participants with all expected organ features non-missing.
    """
    baseline_df = read_baseline_organ_template(organ)
    baseline_feature_cols = get_baseline_feature_columns(baseline_df, organ)

    mapping, missing, malformed = build_feature_mapping(
        baseline_feature_cols=baseline_feature_cols,
        raw_long_cols=list(raw_long.columns),
        organ=organ,
    )

    mapping_file = os.path.join(qc_dir, f"{organ}_instance3_feature_mapping.tsv")
    missing_file = os.path.join(qc_dir, f"{organ}_instance3_missing_features.tsv")
    malformed_file = os.path.join(qc_dir, f"{organ}_instance3_malformed_baseline_features.tsv")

    mapping.to_csv(mapping_file, sep="\t", index=False)

    if not missing.empty:
        missing.to_csv(missing_file, sep="\t", index=False)

    if not malformed.empty:
        malformed.to_csv(malformed_file, sep="\t", index=False)

    n_expected = len(baseline_feature_cols)
    n_mapped = int(mapping["found_in_raw_longitudinal"].sum()) if not mapping.empty else 0
    n_missing = n_expected - n_mapped

    print(f"Organ: {organ}")
    print(f"  Expected model features = {n_expected}")
    print(f"  Mapped instance-{TARGET_INSTANCE} features = {n_mapped}")
    print(f"  Missing instance-{TARGET_INSTANCE} features = {n_missing}")

    if not malformed.empty:
        print(f"  WARNING: {malformed.shape[0]} baseline feature names were malformed.")
        print(f"  See: {malformed_file}")

    # If any expected feature is not available as a raw instance-3 column,
    # complete-case output is impossible for this organ.
    if n_mapped == 0:
        return write_skip_summary(
            organ=organ,
            reason="No expected model features could be mapped to instance-3 raw columns.",
            qc_dir=qc_dir,
            n_expected=n_expected,
            n_mapped=n_mapped,
            n_missing=n_missing,
            n_complete_case=0,
        )

    if n_missing > 0:
        return write_skip_summary(
            organ=organ,
            reason=(
                f"{n_missing} expected model features are not available in instance-3 raw data. "
                "Complete-case organ file was not generated."
            ),
            qc_dir=qc_dir,
            n_expected=n_expected,
            n_mapped=n_mapped,
            n_missing=n_missing,
            n_complete_case=0,
        )

    # Build output table
    out = raw_long[["participant_id"]].copy()

    # Keep these columns so downstream scripts using --feature-start-column diagnosis still work.
    out["session_id"] = OUTPUT_SESSION_ID
    out["diagnosis"] = pd.NA

    # Add feature columns in the exact same order and names as baseline training file.
    for _, row in mapping.iterrows():
        baseline_col = row["baseline_feature_col"]
        raw_target_col = row["raw_target_col"]
        out[baseline_col] = pd.to_numeric(raw_long[raw_target_col], errors="coerce")

    # Complete-case filter: keep participants with ALL organ features present.
    n_before = out.shape[0]
    missing_count_per_person = out[baseline_feature_cols].isna().sum(axis=1)
    out["n_missing_model_features"] = missing_count_per_person
    out["n_observed_model_features"] = len(baseline_feature_cols) - missing_count_per_person

    complete_mask = out["n_missing_model_features"] == 0
    n_complete = int(complete_mask.sum())

    print(f"  Participants before complete-case filter = {n_before}")
    print(f"  Participants with all {n_expected} features observed = {n_complete}")
    print(f"  Participants removed for partial/missing features = {n_before - n_complete}")

    missingness_by_feature = pd.DataFrame({
        "organ": organ,
        "feature": baseline_feature_cols,
        "n_missing": out[baseline_feature_cols].isna().sum().astype(int).values,
        "missing_rate": out[baseline_feature_cols].isna().mean().values,
    }).sort_values(["n_missing", "feature"], ascending=[False, True])

    missingness_file = os.path.join(qc_dir, f"{organ}_instance3_feature_missingness.tsv")
    missingness_by_feature.to_csv(missingness_file, sep="\t", index=False)

    person_missingness_file = os.path.join(qc_dir, f"{organ}_instance3_person_feature_missingness.tsv")
    out[["participant_id", "n_observed_model_features", "n_missing_model_features"]].to_csv(
        person_missingness_file,
        sep="\t",
        index=False,
    )

    if n_complete == 0:
        return write_skip_summary(
            organ=organ,
            reason=(
                "No participants have complete instance-3 data for all expected model features. "
                "Complete-case organ file was not generated."
            ),
            qc_dir=qc_dir,
            n_expected=n_expected,
            n_mapped=n_mapped,
            n_missing=n_missing,
            n_complete_case=0,
        )

    out_complete = out.loc[complete_mask].copy()

    # Do not write helper missingness columns into the model-input file.
    out_complete = out_complete[
        ["participant_id", "session_id", "diagnosis"] + baseline_feature_cols
    ]

    out_file = os.path.join(output_dir, f"imaging_{organ}_3_0.tsv")
    out_complete.to_csv(out_file, sep="\t", index=False)

    qc_summary = {
        "organ": organ,
        "target_instance": TARGET_INSTANCE,
        "model_instance": MODEL_INSTANCE,
        "status": "generated",
        "output_session_id": OUTPUT_SESSION_ID,
        "n_participants_raw_instance3": int(raw_long.shape[0]),
        "n_participants_before_complete_case_filter": int(n_before),
        "n_participants_output_complete_case": int(out_complete.shape[0]),
        "n_expected_model_features": int(n_expected),
        "n_mapped_instance3_features": int(n_mapped),
        "n_missing_instance3_features": int(n_missing),
        "n_malformed_baseline_feature_names": int(malformed.shape[0]),
        "output_file": out_file,
        "mapping_file": mapping_file,
        "feature_missingness_file": missingness_file,
        "person_missingness_file": person_missingness_file,
        "missing_feature_mapping_file": missing_file if not missing.empty else None,
        "malformed_file": malformed_file if not malformed.empty else None,
    }

    qc_json = os.path.join(qc_dir, f"{organ}_instance3_prepare_summary.json")
    with open(qc_json, "w") as f:
        json.dump(qc_summary, f, indent=2)

    print(f"  Saved complete-case file: {out_file}")
    print(f"  QC mapping: {mapping_file}")
    print(f"  Feature missingness: {missingness_file}")
    print(f"  Participant missingness: {person_missingness_file}")
    print(f"  QC summary: {qc_json}")

    return qc_summary


# -----------------------------
# 3. Main
# -----------------------------

def main():
    ensure_dir(OUTPUT_DIR)
    ensure_dir(QC_DIR)

    raw_long = read_abdominal_raw_longitudinal(
        RAW_ABDO_FILE,
        target_instance=TARGET_INSTANCE,
    )

    all_summaries = []

    for organ in ORGAN_LIST:
        try:
            summary = create_organ_longitudinal_file(
                organ=organ,
                raw_long=raw_long,
                output_dir=OUTPUT_DIR,
                qc_dir=QC_DIR,
            )
            all_summaries.append(summary)

        except Exception as e:
            warnings.warn(f"Failed for organ={organ}: {e}")
            summary = write_skip_summary(
                organ=organ,
                reason=f"Unexpected error: {e}",
                qc_dir=QC_DIR,
            )
            all_summaries.append(summary)

    summary_df = pd.DataFrame(all_summaries)
    summary_tsv = os.path.join(QC_DIR, "all_abdominal_organs_instance3_prepare_summary.tsv")
    summary_df.to_csv(summary_tsv, sep="\t", index=False)

    generated = summary_df.loc[summary_df["status"] == "generated", "organ"].tolist()
    skipped = summary_df.loc[summary_df["status"] == "skipped", "organ"].tolist()

    print("\nFinished preparing longitudinal abdominal imaging instance-3 complete-case files.")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"QC directory: {QC_DIR}")
    print(f"Combined QC summary: {summary_tsv}")

    print("\nGenerated organs:")
    if generated:
        for organ in generated:
            print(f"  {organ}: {os.path.join(OUTPUT_DIR, f'imaging_{organ}_3_0.tsv')}")
    else:
        print("  None")

    print("\nSkipped organs:")
    if skipped:
        for organ in skipped:
            print(f"  {organ}")
    else:
        print("  None")


if __name__ == "__main__":
    main()