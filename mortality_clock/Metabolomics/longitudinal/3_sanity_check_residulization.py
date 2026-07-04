#!/usr/bin/env python3
# ============================================================
# Sanity check for 4 metabolomics mortality clocks
#
# Checks both:
#   1) baseline instance 0 predictions
#   2) longitudinal instance 1 predictions
#
# Organs:
#   Endocrine, Digestive, Hepatic, Immune
#
# Main columns checked:
#   {organ}_metabolomics_mortality_risk_score
#   {organ}_metabolomics_mortality_clock_acceleration_z
#   {organ}_metabolomics_mortality_clock_acceleration_years
#   {organ}_metabolomics_mortality_clock_age_years
#
# Outputs:
#   printed summaries
#   metabolomics_clock_sanity_check_summary.tsv
#   metabolomics_clock_sanity_check_extreme_rows.tsv
# ============================================================

import numpy as np
import pandas as pd
from pathlib import Path


# -----------------------------
# 1. Paths
# -----------------------------

root = Path("/cbica/home/wenju/Reproducibile_paper/WholeBodyClock")

long_root = Path(
    "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/"
    "mortality_clock/longitudinal/metabolomics"
)

out_summary = long_root / "metabolomics_clock_sanity_check_summary.tsv"
out_extreme = long_root / "metabolomics_clock_sanity_check_extreme_rows.tsv"


# -----------------------------
# 2. Organs
# -----------------------------

organs = {
    "Endocrine": "endocrine",
    "Digestive": "digestive",
    "Hepatic": "hepatic",
    "Immune": "immune",
}


# -----------------------------
# 3. Sanity thresholds
# -----------------------------

z_abs_warn = 50
years_abs_warn = 200
clock_age_min = 0
clock_age_max = 150


# -----------------------------
# 4. Helper functions
# -----------------------------

def describe_numeric(x):
    x = pd.to_numeric(x, errors="coerce")
    desc = x.describe(percentiles=[0.01, 0.05, 0.5, 0.95, 0.99])
    return {
        "n": int(x.notna().sum()),
        "n_missing": int(x.isna().sum()),
        "mean": desc.get("mean", np.nan),
        "std": desc.get("std", np.nan),
        "min": desc.get("min", np.nan),
        "p01": desc.get("1%", np.nan),
        "p05": desc.get("5%", np.nan),
        "median": desc.get("50%", np.nan),
        "p95": desc.get("95%", np.nan),
        "p99": desc.get("99%", np.nan),
        "max": desc.get("max", np.nan),
        "max_abs": float(np.nanmax(np.abs(x))) if x.notna().any() else np.nan,
    }


def get_files(organ_label, organ):
    baseline_file = (
        root
        / f"{organ_label}_metabolomics_mortality_clock"
        / f"{organ}_metabolomics_mortality_clock_predictions.tsv"
    )

    instance1_file = (
        long_root
        / organ_label
        / f"{organ}_metabolomics_mortality_clock_apply_instance_1_0_predictions.tsv"
    )

    return {
        "baseline_0_0": baseline_file,
        "instance_1_0": instance1_file,
    }


def get_columns(organ):
    return [
        f"{organ}_metabolomics_mortality_risk_score",
        f"{organ}_metabolomics_mortality_clock_acceleration_z",
        f"{organ}_metabolomics_mortality_clock_acceleration_years",
        f"{organ}_metabolomics_mortality_clock_age_years",
    ]


def classify_column(col):
    if col.endswith("_mortality_risk_score"):
        return "risk_score"
    if col.endswith("_clock_acceleration_z"):
        return "acceleration_z"
    if col.endswith("_clock_acceleration_years"):
        return "acceleration_years"
    if col.endswith("_clock_age_years"):
        return "clock_age_years"
    return "other"


def flag_column_values(df, col, organ_label, organ, file_label, path):
    x = pd.to_numeric(df[col], errors="coerce")
    col_type = classify_column(col)

    if col_type == "acceleration_z":
        bad = x.abs() > z_abs_warn
        reason = f"abs(z) > {z_abs_warn}"
    elif col_type == "acceleration_years":
        bad = x.abs() > years_abs_warn
        reason = f"abs(years) > {years_abs_warn}"
    elif col_type == "clock_age_years":
        bad = (x < clock_age_min) | (x > clock_age_max)
        reason = f"clock age outside [{clock_age_min}, {clock_age_max}]"
    else:
        bad = pd.Series(False, index=df.index)
        reason = ""

    if not bad.any():
        return pd.DataFrame()

    keep_cols = ["participant_id"]
    for c in [
        "application_instance",
        "application_source_file",
        "sample_date",
        "death_date",
        "age_at_baseline",
        "age_at_imaging",
        "sex",
    ]:
        if c in df.columns:
            keep_cols.append(c)

    out = df.loc[bad, keep_cols].copy()
    out["organ_label"] = organ_label
    out["organ"] = organ
    out["file_label"] = file_label
    out["file"] = str(path)
    out["column"] = col
    out["value"] = x.loc[bad].values
    out["reason"] = reason

    return out


# -----------------------------
# 5. Main sanity check
# -----------------------------

summary_rows = []
extreme_rows = []

for organ_label, organ in organs.items():
    print("\n" + "=" * 90)
    print(f"Organ: {organ_label} ({organ})")
    print("=" * 90)

    files = get_files(organ_label, organ)
    cols = get_columns(organ)

    for file_label, path in files.items():
        print("\n" + "-" * 90)
        print(f"{file_label}: {path}")
        print("-" * 90)

        if not path.exists():
            print(f"ERROR: file not found: {path}")
            summary_rows.append({
                "organ_label": organ_label,
                "organ": organ,
                "file_label": file_label,
                "file": str(path),
                "column": None,
                "column_type": None,
                "status": "file_missing",
            })
            continue

        df = pd.read_csv(path, sep="\t")
        print(f"N rows: {df.shape[0]}")
        print(f"N columns: {df.shape[1]}")

        for col in cols:
            col_type = classify_column(col)

            if col not in df.columns:
                print(f"\nMISSING COLUMN: {col}")
                summary_rows.append({
                    "organ_label": organ_label,
                    "organ": organ,
                    "file_label": file_label,
                    "file": str(path),
                    "column": col,
                    "column_type": col_type,
                    "status": "column_missing",
                })
                continue

            x = pd.to_numeric(df[col], errors="coerce")
            stats = describe_numeric(x)

            status = "ok"
            warning_notes = []

            if col_type == "acceleration_z" and stats["max_abs"] > z_abs_warn:
                status = "warning"
                warning_notes.append(f"max_abs_z>{z_abs_warn}")

            if col_type == "acceleration_years" and stats["max_abs"] > years_abs_warn:
                status = "warning"
                warning_notes.append(f"max_abs_years>{years_abs_warn}")

            if col_type == "clock_age_years":
                if stats["min"] < clock_age_min or stats["max"] > clock_age_max:
                    status = "warning"
                    warning_notes.append(f"clock_age_outside_{clock_age_min}_{clock_age_max}")

            print(f"\n{col}")
            print(x.describe(percentiles=[0.01, 0.05, 0.5, 0.95, 0.99]))
            print(f"max_abs = {stats['max_abs']}")
            print(f"status = {status}" + (f" ({'; '.join(warning_notes)})" if warning_notes else ""))

            summary_rows.append({
                "organ_label": organ_label,
                "organ": organ,
                "file_label": file_label,
                "file": str(path),
                "column": col,
                "column_type": col_type,
                "status": status,
                "warning_notes": ";".join(warning_notes),
                **stats,
            })

            bad_df = flag_column_values(
                df=df,
                col=col,
                organ_label=organ_label,
                organ=organ,
                file_label=file_label,
                path=path,
            )

            if bad_df.shape[0] > 0:
                extreme_rows.append(bad_df)


# -----------------------------
# 6. Save outputs
# -----------------------------

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(out_summary, sep="\t", index=False)

if extreme_rows:
    extreme_df = pd.concat(extreme_rows, axis=0, ignore_index=True)
else:
    extreme_df = pd.DataFrame()

extreme_df.to_csv(out_extreme, sep="\t", index=False)

print("\n" + "=" * 90)
print("Finished metabolomics clock sanity check.")
print("=" * 90)
print(f"Summary saved to: {out_summary}")
print(f"Extreme rows saved to: {out_extreme}")

print("\nWarning summary:")
if "status" in summary_df.columns:
    warn = summary_df.loc[summary_df["status"] == "warning"].copy()
    if warn.empty:
        print("No warning-level values detected.")
    else:
        print(
            warn[
                [
                    "organ_label",
                    "file_label",
                    "column",
                    "warning_notes",
                    "min",
                    "median",
                    "max",
                    "max_abs",
                ]
            ].to_string(index=False)
        )