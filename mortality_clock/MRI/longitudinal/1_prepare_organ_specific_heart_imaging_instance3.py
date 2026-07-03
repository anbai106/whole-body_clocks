#!/usr/bin/env python3
# ============================================================
# Prepare complete-case longitudinal heart MRI data for applying
# the pre-trained heart MRI mortality clock to UKB instance 3.
#
# Goal:
#   Create a TSV with the same feature-column names as the
#   baseline instance-2 heart training file, but filled with
#   instance-3 heart MRI values from UKBB_fullsample_heart.csv.
#
# Example:
#   Baseline training column:
#     lv_end_diastolic_volume_f24100_2_0
#
#   Raw longitudinal column:
#     lv_end_diastolic_volume_f24100_3_0
#
#   Output longitudinal column:
#     lv_end_diastolic_volume_f24100_2_0
#
# This preserves the exact model-expected feature names while using
# longitudinal instance-3 values.
#
# Complete-case rule:
#   Keep only participants with ALL required heart MRI features
#   non-missing at instance 3.
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

RAW_HEART_FILE = (
    "/Users/hao/cubic-home/Reproducibile_paper/BrainHeart/data/UKBB_fullsample_heart.csv"
)

BASELINE_HEART_FILE = (
    "/Users/hao/cubic-home/Reproducibile_paper/AbdoImaging/"
    "Data_split_by_Huizi/MLNI_dataset_splits/heart_train_val_test.tsv"
)

OUTPUT_DIR = (
    "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/"
    "mortality_clock/longitudinal/imaging/data"
)

QC_DIR = os.path.join(OUTPUT_DIR, "qc_heart_instance3")

TARGET_INSTANCE = "3"
MODEL_INSTANCE = "2"
OUTPUT_SESSION_ID = "ses-M3"

OUTPUT_FILE = os.path.join(OUTPUT_DIR, "imaging_heart_3_0.tsv")


# -----------------------------
# 2. Helper functions
# -----------------------------

def ensure_dir(path: str) -> None:
    Path(path).mkdir(parents=True, exist_ok=True)


def extract_heart_feature_suffix(colname: str):
    """
    Extract UKB field-instance-array suffix from a heart MRI feature name.

    Examples:
      lv_end_diastolic_volume_f24100_2_0
      lv_circumferential_strain_aha_1_f24141_2_0
      lv_longitudinal_strain_global_f24181_3_0

    Returns:
      field_id, instance, array
      e.g., ("24100", "2", "0")
    """
    m = re.search(r"_f(\d+)_(\d+)_(\d+)$", str(colname))
    if m is None:
        return None

    return {
        "field_id": m.group(1),
        "instance": m.group(2),
        "array": m.group(3),
        "full_suffix": f"f{m.group(1)}_{m.group(2)}_{m.group(3)}",
    }


def make_target_raw_col_from_baseline_col(
    baseline_col: str,
    target_instance: str = TARGET_INSTANCE,
):
    """
    Convert a baseline model feature column name to the matching
    raw instance-3 column name.

    Example:
      lv_end_diastolic_volume_f24100_2_0
      -> lv_end_diastolic_volume_f24100_3_0
    """
    parsed = extract_heart_feature_suffix(baseline_col)
    if parsed is None:
        return None

    old_suffix = f"_f{parsed['field_id']}_{parsed['instance']}_{parsed['array']}"
    new_suffix = f"_f{parsed['field_id']}_{target_instance}_{parsed['array']}"

    if not baseline_col.endswith(old_suffix):
        return None

    return baseline_col[: -len(old_suffix)] + new_suffix


def read_raw_heart_longitudinal(
    raw_file: str,
    target_instance: str = TARGET_INSTANCE,
) -> pd.DataFrame:
    """
    Read raw heart MRI file and keep eid plus all target-instance heart columns.

    Raw heart columns look like:
      lv_end_diastolic_volume_f24100_2_0
      lv_end_diastolic_volume_f24100_3_0
    """
    print(f"Reading raw heart MRI file:\n  {raw_file}")

    df = pd.read_csv(raw_file)

    if "eid" not in df.columns:
        raise ValueError("Raw heart MRI file must contain column 'eid'.")

    target_pattern = re.compile(rf"_f\d+_{re.escape(target_instance)}_\d+$")
    keep_cols = ["eid"] + [c for c in df.columns if target_pattern.search(str(c))]

    if len(keep_cols) == 1:
        raise ValueError(
            f"No instance-{target_instance} heart MRI columns found in raw file."
        )

    df = df[keep_cols].copy()
    df = df.rename(columns={"eid": "participant_id"})

    df["participant_id"] = pd.to_numeric(df["participant_id"], errors="coerce")
    df = df.dropna(subset=["participant_id"]).copy()
    df["participant_id"] = df["participant_id"].astype(int)

    n_dup = int(df["participant_id"].duplicated().sum())
    if n_dup > 0:
        warnings.warn(
            f"Found {n_dup} duplicated participant_id rows in raw heart file. "
            "Keeping the first occurrence."
        )
        df = df.drop_duplicates("participant_id", keep="first")

    print("Raw longitudinal heart MRI data:")
    print(f"  N participants = {df.shape[0]}")
    print(f"  N instance-{target_instance} heart MRI columns = {df.shape[1] - 1}")

    return df


def read_baseline_heart_template(baseline_file: str = BASELINE_HEART_FILE) -> pd.DataFrame:
    """
    Read the baseline heart split file used to train the heart MRI mortality clock.

    Expected columns:
      participant_id, session_id, diagnosis, feature1, feature2, ...
    """
    print(f"\nReading baseline heart template:\n  {baseline_file}")

    if not os.path.exists(baseline_file):
        raise FileNotFoundError(f"Baseline heart file not found: {baseline_file}")

    df = pd.read_csv(baseline_file, sep="\t")

    if "participant_id" not in df.columns:
        raise ValueError("Baseline heart file is missing participant_id.")

    return df


def get_baseline_heart_feature_columns(baseline_df: pd.DataFrame):
    """
    Robustly infer heart MRI feature columns from the baseline training file.

    This works for either:
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
        if extract_heart_feature_suffix(c) is not None:
            feature_cols.append(c)
        else:
            skipped_cols.append(c)

    if skipped_cols:
        print(f"Skipped {len(skipped_cols)} non-feature columns from baseline heart file.")
        for c in skipped_cols[:20]:
            print(f"  skipped: {c}")
        if len(skipped_cols) > 20:
            print("  ...")

    if len(feature_cols) == 0:
        raise ValueError(
            "No baseline heart MRI features found. Expected feature names ending "
            "with suffix like _f24100_2_0."
        )

    return feature_cols


def build_heart_feature_mapping(
    baseline_feature_cols,
    raw_long_cols,
    target_instance: str = TARGET_INSTANCE,
):
    """
    Build mapping from baseline model-expected heart feature names
    to raw instance-3 heart feature names.

    Output keeps baseline feature names but uses raw instance-3 values.
    """
    raw_long_set = set(raw_long_cols)

    rows = []
    missing_rows = []
    malformed_rows = []

    for baseline_col in baseline_feature_cols:
        parsed = extract_heart_feature_suffix(baseline_col)

        if parsed is None:
            malformed_rows.append({
                "baseline_feature_col": baseline_col,
                "reason": "Could not extract suffix like _f24100_2_0",
            })
            continue

        raw_target_col = make_target_raw_col_from_baseline_col(
            baseline_col,
            target_instance=target_instance,
        )

        row = {
            "baseline_feature_col": baseline_col,
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
    reason: str,
    qc_dir: str = QC_DIR,
    n_expected: int = None,
    n_mapped: int = None,
    n_missing: int = None,
    n_complete_case: int = None,
):
    summary = {
        "organ": "heart",
        "target_instance": TARGET_INSTANCE,
        "model_instance": MODEL_INSTANCE,
        "status": "skipped",
        "reason": reason,
        "n_expected_model_features": n_expected,
        "n_mapped_instance3_features": n_mapped,
        "n_missing_instance3_features": n_missing,
        "n_complete_case_participants": n_complete_case,
    }

    out_json = os.path.join(qc_dir, "heart_instance3_prepare_summary.json")
    with open(out_json, "w") as f:
        json.dump(summary, f, indent=2)

    warnings.warn(f"Skipping heart: {reason}")
    return summary


def prepare_heart_instance3_complete_case():
    ensure_dir(OUTPUT_DIR)
    ensure_dir(QC_DIR)

    raw_long = read_raw_heart_longitudinal(
        RAW_HEART_FILE,
        target_instance=TARGET_INSTANCE,
    )

    baseline_df = read_baseline_heart_template(BASELINE_HEART_FILE)
    baseline_feature_cols = get_baseline_heart_feature_columns(baseline_df)

    mapping, missing, malformed = build_heart_feature_mapping(
        baseline_feature_cols=baseline_feature_cols,
        raw_long_cols=list(raw_long.columns),
        target_instance=TARGET_INSTANCE,
    )

    mapping_file = os.path.join(QC_DIR, "heart_instance3_feature_mapping.tsv")
    missing_file = os.path.join(QC_DIR, "heart_instance3_missing_features.tsv")
    malformed_file = os.path.join(QC_DIR, "heart_instance3_malformed_baseline_features.tsv")

    mapping.to_csv(mapping_file, sep="\t", index=False)

    if not missing.empty:
        missing.to_csv(missing_file, sep="\t", index=False)

    if not malformed.empty:
        malformed.to_csv(malformed_file, sep="\t", index=False)

    n_expected = len(baseline_feature_cols)
    n_mapped = int(mapping["found_in_raw_longitudinal"].sum()) if not mapping.empty else 0
    n_missing = n_expected - n_mapped

    print("\nHeart MRI feature mapping:")
    print(f"  Expected model features = {n_expected}")
    print(f"  Mapped instance-{TARGET_INSTANCE} features = {n_mapped}")
    print(f"  Missing instance-{TARGET_INSTANCE} features = {n_missing}")

    if not malformed.empty:
        print(f"  WARNING: {malformed.shape[0]} baseline feature names were malformed.")
        print(f"  See: {malformed_file}")

    if n_mapped == 0:
        return write_skip_summary(
            reason="No expected heart MRI features could be mapped to instance-3 raw columns.",
            n_expected=n_expected,
            n_mapped=n_mapped,
            n_missing=n_missing,
            n_complete_case=0,
        )

    if n_missing > 0:
        return write_skip_summary(
            reason=(
                f"{n_missing} expected heart MRI model features are not available "
                "in the raw instance-3 file. Complete-case file was not generated."
            ),
            n_expected=n_expected,
            n_mapped=n_mapped,
            n_missing=n_missing,
            n_complete_case=0,
        )

    # Build output table.
    # Keep session_id and diagnosis because many application scripts use
    # --feature-start-column diagnosis.
    out = raw_long[["participant_id"]].copy()
    out["session_id"] = OUTPUT_SESSION_ID
    out["diagnosis"] = pd.NA

    # Add features in the exact same order and names as baseline training file.
    for _, row in mapping.iterrows():
        baseline_col = row["baseline_feature_col"]
        raw_target_col = row["raw_target_col"]
        out[baseline_col] = pd.to_numeric(raw_long[raw_target_col], errors="coerce")

    # Complete-case filter: participant must have all expected heart features.
    n_before = out.shape[0]

    missing_count_per_person = out[baseline_feature_cols].isna().sum(axis=1)
    out["n_missing_model_features"] = missing_count_per_person
    out["n_observed_model_features"] = len(baseline_feature_cols) - missing_count_per_person

    complete_mask = out["n_missing_model_features"] == 0
    n_complete = int(complete_mask.sum())

    print("\nComplete-case filtering:")
    print(f"  Participants before complete-case filter = {n_before}")
    print(f"  Participants with all {n_expected} heart features observed = {n_complete}")
    print(f"  Participants removed for partial/missing heart features = {n_before - n_complete}")

    # Feature-level missingness QC
    missingness_by_feature = pd.DataFrame({
        "organ": "heart",
        "feature": baseline_feature_cols,
        "n_missing": out[baseline_feature_cols].isna().sum().astype(int).values,
        "missing_rate": out[baseline_feature_cols].isna().mean().values,
    }).sort_values(["n_missing", "feature"], ascending=[False, True])

    missingness_file = os.path.join(QC_DIR, "heart_instance3_feature_missingness.tsv")
    missingness_by_feature.to_csv(missingness_file, sep="\t", index=False)

    # Participant-level missingness QC
    person_missingness_file = os.path.join(QC_DIR, "heart_instance3_person_feature_missingness.tsv")
    out[["participant_id", "n_observed_model_features", "n_missing_model_features"]].to_csv(
        person_missingness_file,
        sep="\t",
        index=False,
    )

    if n_complete == 0:
        return write_skip_summary(
            reason=(
                "No participants have complete instance-3 data for all expected "
                "heart MRI model features. Complete-case file was not generated."
            ),
            n_expected=n_expected,
            n_mapped=n_mapped,
            n_missing=n_missing,
            n_complete_case=0,
        )

    out_complete = out.loc[complete_mask].copy()

    # Do not write helper missingness columns into model-input file.
    out_complete = out_complete[
        ["participant_id", "session_id", "diagnosis"] + baseline_feature_cols
    ]

    out_complete.to_csv(OUTPUT_FILE, sep="\t", index=False)

    qc_summary = {
        "organ": "heart",
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
        "output_file": OUTPUT_FILE,
        "mapping_file": mapping_file,
        "feature_missingness_file": missingness_file,
        "person_missingness_file": person_missingness_file,
        "missing_feature_mapping_file": missing_file if not missing.empty else None,
        "malformed_file": malformed_file if not malformed.empty else None,
    }

    qc_json = os.path.join(QC_DIR, "heart_instance3_prepare_summary.json")
    with open(qc_json, "w") as f:
        json.dump(qc_summary, f, indent=2)

    print("\nSaved complete-case heart instance-3 file:")
    print(f"  {OUTPUT_FILE}")
    print("\nQC files:")
    print(f"  Mapping: {mapping_file}")
    print(f"  Feature missingness: {missingness_file}")
    print(f"  Participant missingness: {person_missingness_file}")
    print(f"  Summary: {qc_json}")

    return qc_summary


# -----------------------------
# 3. Main
# -----------------------------

def main():
    summary = prepare_heart_instance3_complete_case()

    summary_tsv = os.path.join(QC_DIR, "heart_instance3_prepare_summary.tsv")
    pd.DataFrame([summary]).to_csv(summary_tsv, sep="\t", index=False)

    print("\nFinished preparing longitudinal heart MRI instance-3 complete-case file.")
    print(f"QC summary TSV: {summary_tsv}")

    if summary.get("status") == "generated":
        print(f"Generated file: {OUTPUT_FILE}")
    else:
        print("No heart instance-3 file was generated.")
        print(f"Reason: {summary.get('reason')}")


if __name__ == "__main__":
    main()