#!/usr/bin/env python3

"""
Summarize baseline demographics for AIBL, OASIS, and BLSA.

Primary outputs
---------------
1. Study-level baseline demographic summary:
       Study
       N
       Age_mean
       Age_sd
       Age_mean_sd
       Female_n
       Female_percent
       Female_n_percent

2. Site-level baseline demographic summary.

3. Participant-level selected baseline records for QC.

Baseline selection
------------------
One harmonized AD EPOCH scan is selected per Study + participant_id using:

    1. explicit baseline scan flag;
    2. longitudinal scan number == 1;
    3. smallest absolute years-since-baseline;
    4. explicit baseline visit code;
    5. earliest scan date;
    6. youngest scan age;
    7. first available row.

Sex source
----------
The script first uses a Sex column from the harmonized prediction file when
available. If sex is absent or missing, it matches the selected baseline scan
to the closest record in the full iSTAGING sample-information file by:

    1. closest Date;
    2. closest Age;
    3. earliest available sample-information record.

Running
-------
The script contains project-specific default paths and can therefore be run
directly in PyCharm without command-line arguments.

Optional command-line overrides are still supported:

python summarize_external_baseline_demographics.py \
  --prediction-file /path/to/predictions.tsv \
  --sample-file /path/to/external_5_studies_istaging.tsv \
  --outdir /path/to/output
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd


STUDY_ORDER = ["AIBL", "OASIS", "BLSA"]

# Country in which each cohort collected its data.
STUDY_COUNTRY = {
    "AIBL": "Australia",
    "OASIS": "USA",
    "BLSA": "USA",
}


def log(message: str) -> None:
    print(message, flush=True)


def clean_string(series: pd.Series) -> pd.Series:
    output = series.astype("string").str.strip()
    output = output.replace(
        {
            "": pd.NA,
            "NA": pd.NA,
            "NaN": pd.NA,
            "nan": pd.NA,
            "None": pd.NA,
            "null": pd.NA,
            "<NA>": pd.NA,
        }
    )
    return output


def normalize_study(series: pd.Series) -> pd.Series:
    raw = clean_string(series)
    upper = raw.str.upper()

    output = raw.copy()
    output.loc[upper.str.startswith("AIBL", na=False)] = "AIBL"
    output.loc[upper.str.startswith("OASIS", na=False)] = "OASIS"
    output.loc[upper.str.startswith("BLSA", na=False)] = "BLSA"
    return output


def normalize_sex(series: pd.Series) -> pd.Series:
    """
    Return Female, Male, or NA for common sex encodings.
    """
    raw = clean_string(series).str.upper()

    output = pd.Series(pd.NA, index=series.index, dtype="string")

    female_values = {
        "F",
        "FEMALE",
        "WOMAN",
        "WOMEN",
        "0",
        "2",
    }

    male_values = {
        "M",
        "MALE",
        "MAN",
        "MEN",
        "1",
    }

    output.loc[raw.isin(female_values)] = "Female"
    output.loc[raw.isin(male_values)] = "Male"

    output.loc[
        output.isna()
        & raw.str.contains(r"^F", regex=True, na=False)
    ] = "Female"

    output.loc[
        output.isna()
        & raw.str.contains(r"^M", regex=True, na=False)
    ] = "Male"

    return output


def parse_date(series: pd.Series) -> pd.Series:
    return pd.to_datetime(series, errors="coerce")


def find_column(
    columns: Sequence[str],
    preferred: Sequence[str],
    regex: Optional[str],
    label: str,
    required: bool = True,
) -> Optional[str]:
    for candidate in preferred:
        if candidate in columns:
            return candidate

    if regex is not None:
        matches = [
            column
            for column in columns
            if re.search(regex, column, flags=re.IGNORECASE)
        ]

        if len(matches) == 1:
            return matches[0]

        if len(matches) > 1:
            raise ValueError(
                f"Multiple candidate columns found for {label}: "
                + ", ".join(matches)
            )

    if required:
        raise ValueError(
            f"Could not identify {label}. Available columns include:\n"
            + ", ".join(list(columns)[:180])
        )

    return None


def visit_to_month(series: pd.Series) -> pd.Series:
    raw = clean_string(series).str.lower()
    output = pd.Series(np.nan, index=series.index, dtype=float)

    baseline_mask = raw.isin(
        [
            "bl",
            "base",
            "baseline",
            "m00",
            "m0",
            "screen",
            "screening",
            "sc",
        ]
    )
    output.loc[baseline_mask] = 0.0

    unresolved = output.isna() & raw.notna()

    if unresolved.any():
        extracted = raw.loc[unresolved].str.extract(
            r"m(?:onth)?\s*0*([0-9]+)",
            expand=False,
        )
        output.loc[unresolved] = pd.to_numeric(
            extracted,
            errors="coerce",
        )

    unresolved = output.isna() & raw.notna()

    if unresolved.any():
        extracted = raw.loc[unresolved].str.extract(
            r"([0-9]+)",
            expand=False,
        )
        output.loc[unresolved] = pd.to_numeric(
            extracted,
            errors="coerce",
        )

    return output


def coerce_boolean(series: pd.Series) -> pd.Series:
    raw = clean_string(series).str.upper()

    output = pd.Series(pd.NA, index=series.index, dtype="boolean")

    output.loc[
        raw.isin(
            [
                "TRUE",
                "T",
                "YES",
                "Y",
                "1",
            ]
        )
    ] = True

    output.loc[
        raw.isin(
            [
                "FALSE",
                "F",
                "NO",
                "N",
                "0",
            ]
        )
    ] = False

    return output


def read_prediction_file(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(
        path,
        sep="\t",
        low_memory=False,
    )

    columns = raw.columns.tolist()

    id_col = find_column(
        columns,
        ["PTID", "participant_id", "IID", "eid"],
        r"(^PTID$|participant.*id|^IID$|^eid$)",
        "participant ID",
    )

    study_col = find_column(
        columns,
        ["external_Study", "Study", "STUDY"],
        r"(^|_)study$",
        "Study",
    )

    site_col = find_column(
        columns,
        ["external_SITE", "SITE", "Site"],
        r"(^|_)site$",
        "SITE",
        required=False,
    )

    age_col = find_column(
        columns,
        ["Age", "age_at_scan_used_for_model", "AGE"],
        r"(^|_)age($|_at_scan)",
        "Age",
    )

    sex_col = find_column(
        columns,
        ["Sex", "SEX", "sex"],
        r"^sex$",
        "Sex",
        required=False,
    )

    date_col = find_column(
        columns,
        ["Date", "scan_date", "MRI_Date"],
        r"(^|_)date$",
        "scan Date",
        required=False,
    )

    visit_col = find_column(
        columns,
        ["Visit_Code", "VISCODE", "visit_code", "Visit"],
        r"(visit|viscode)",
        "visit code",
        required=False,
    )

    baseline_flag_col = find_column(
        columns,
        ["is_external_baseline_scan", "is_baseline_scan"],
        r"baseline.*scan",
        "baseline scan indicator",
        required=False,
    )

    scan_number_col = find_column(
        columns,
        ["longitudinal_scan_number", "scan_number"],
        r"scan.*number",
        "longitudinal scan number",
        required=False,
    )

    years_col = find_column(
        columns,
        [
            "years_since_external_baseline",
            "Delta_Baseline",
            "years_since_baseline",
        ],
        r"(years.*baseline|delta_baseline)",
        "years since baseline",
        required=False,
    )

    prediction = pd.DataFrame(
        {
            "prediction_source_row": np.arange(
                raw.shape[0],
                dtype=int,
            ),
            "participant_id": clean_string(
                raw[id_col]
            ),
            "study": normalize_study(
                raw[study_col]
            ),
            "site": (
                clean_string(raw[site_col])
                if site_col is not None
                else pd.Series(
                    "Unknown",
                    index=raw.index,
                    dtype="string",
                )
            ),
            "age": pd.to_numeric(
                raw[age_col],
                errors="coerce",
            ),
            "sex_prediction": (
                normalize_sex(raw[sex_col])
                if sex_col is not None
                else pd.Series(
                    pd.NA,
                    index=raw.index,
                    dtype="string",
                )
            ),
            "scan_date": (
                parse_date(raw[date_col])
                if date_col is not None
                else pd.Series(
                    pd.NaT,
                    index=raw.index,
                )
            ),
            "visit_code": (
                clean_string(raw[visit_col])
                if visit_col is not None
                else pd.Series(
                    pd.NA,
                    index=raw.index,
                    dtype="string",
                )
            ),
            "baseline_flag": (
                coerce_boolean(raw[baseline_flag_col])
                if baseline_flag_col is not None
                else pd.Series(
                    pd.NA,
                    index=raw.index,
                    dtype="boolean",
                )
            ),
            "scan_number": (
                pd.to_numeric(
                    raw[scan_number_col],
                    errors="coerce",
                )
                if scan_number_col is not None
                else pd.Series(
                    np.nan,
                    index=raw.index,
                    dtype=float,
                )
            ),
            "years_since_baseline": (
                pd.to_numeric(
                    raw[years_col],
                    errors="coerce",
                )
                if years_col is not None
                else pd.Series(
                    np.nan,
                    index=raw.index,
                    dtype=float,
                )
            ),
        }
    )

    prediction = prediction.loc[
        prediction["study"].isin(STUDY_ORDER)
        & prediction["participant_id"].notna()
        & prediction["age"].notna()
    ].copy()

    return prediction


def select_one_baseline(prediction: pd.DataFrame) -> pd.DataFrame:
    data = prediction.copy()

    data["visit_month"] = visit_to_month(
        data["visit_code"]
    )

    data["explicit_baseline_visit"] = (
        data["visit_month"] == 0
    )

    # np.select requires ordinary boolean ndarrays. Pandas nullable
    # BooleanArray objects may contain pd.NA and therefore must be converted.
    condition_baseline_flag = (
        data["baseline_flag"]
        .fillna(False)
        .astype(bool)
        .to_numpy()
    )

    condition_scan_number = (
        data["scan_number"]
        .eq(1)
        .fillna(False)
        .to_numpy(dtype=bool)
    )

    condition_years_available = (
        data["years_since_baseline"]
        .notna()
        .to_numpy(dtype=bool)
    )

    condition_baseline_visit = (
        data["explicit_baseline_visit"]
        .fillna(False)
        .astype(bool)
        .to_numpy()
    )

    condition_date_available = (
        data["scan_date"]
        .notna()
        .to_numpy(dtype=bool)
    )

    condition_age_available = (
        data["age"]
        .notna()
        .to_numpy(dtype=bool)
    )

    data["baseline_priority"] = np.select(
        [
            condition_baseline_flag,
            condition_scan_number,
            condition_years_available,
            condition_baseline_visit,
            condition_date_available,
            condition_age_available,
        ],
        [
            1,
            2,
            3,
            4,
            5,
            6,
        ],
        default=7,
    )

    data["baseline_distance"] = np.select(
        [
            data["baseline_priority"] == 1,
            data["baseline_priority"] == 2,
            data["baseline_priority"] == 3,
            data["baseline_priority"] == 4,
            data["baseline_priority"] == 5,
            data["baseline_priority"] == 6,
        ],
        [
            0.0,
            0.0,
            data["years_since_baseline"].abs(),
            data["visit_month"].abs(),
            data["scan_date"].map(
                lambda value: (
                    value.toordinal()
                    if pd.notna(value)
                    else np.nan
                )
            ),
            data["age"],
        ],
        default=data["prediction_source_row"].astype(float),
    )

    data = data.sort_values(
        [
            "study",
            "participant_id",
            "baseline_priority",
            "baseline_distance",
            "scan_date",
            "age",
            "prediction_source_row",
        ],
        kind="mergesort",
    )

    baseline = (
        data.groupby(
            [
                "study",
                "participant_id",
            ],
            as_index=False,
            sort=False,
        )
        .head(1)
        .copy()
    )

    baseline["baseline_selection_source"] = np.select(
        [
            (baseline["baseline_priority"] == 1).to_numpy(dtype=bool),
            (baseline["baseline_priority"] == 2).to_numpy(dtype=bool),
            (baseline["baseline_priority"] == 3).to_numpy(dtype=bool),
            (baseline["baseline_priority"] == 4).to_numpy(dtype=bool),
            (baseline["baseline_priority"] == 5).to_numpy(dtype=bool),
            (baseline["baseline_priority"] == 6).to_numpy(dtype=bool),
        ],
        [
            "explicit baseline flag",
            "scan number 1",
            "closest years-since-baseline to zero",
            "explicit baseline visit code",
            "earliest scan Date",
            "youngest scan Age",
        ],
        default="first available row",
    )

    return baseline



def calculate_longitudinal_followup(
    prediction: pd.DataFrame,
    baseline: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Calculate participant-level longitudinal follow-up from the selected
    baseline scan to the latest available harmonized prediction scan.

    Time is calculated in the following order:
        1. Scan Date difference when baseline and follow-up dates are available.
        2. Scan Age difference otherwise.

    Negative values, which can arise from imperfect metadata alignment, are
    excluded. Participants with only one qualified scan have zero follow-up.
    """
    baseline_reference = baseline[
        [
            "study",
            "participant_id",
            "scan_date",
            "age",
        ]
    ].rename(
        columns={
            "scan_date": "selected_baseline_date",
            "age": "selected_baseline_age",
        }
    )

    longitudinal = prediction.merge(
        baseline_reference,
        on=["study", "participant_id"],
        how="inner",
        validate="many_to_one",
    )

    longitudinal["followup_years_by_date"] = (
        longitudinal["scan_date"]
        - longitudinal["selected_baseline_date"]
    ).dt.total_seconds() / (365.25 * 24.0 * 60.0 * 60.0)

    longitudinal["followup_years_by_age"] = (
        longitudinal["age"]
        - longitudinal["selected_baseline_age"]
    )

    longitudinal["followup_years"] = (
        longitudinal["followup_years_by_date"]
        .where(
            longitudinal["followup_years_by_date"].notna(),
            longitudinal["followup_years_by_age"],
        )
    )

    longitudinal["followup_time_source"] = np.select(
        [
            longitudinal["followup_years_by_date"].notna().to_numpy(dtype=bool),
            longitudinal["followup_years_by_age"].notna().to_numpy(dtype=bool),
        ],
        [
            "Date",
            "Age",
        ],
        default="Unavailable",
    )

    # Small negative values can reflect rounding between dates and ages.
    longitudinal.loc[
        longitudinal["followup_years"].between(-0.05, 0, inclusive="left"),
        "followup_years",
    ] = 0.0

    longitudinal = longitudinal.loc[
        longitudinal["followup_years"].notna()
        & (longitudinal["followup_years"] >= 0)
    ].copy()

    participant_followup = (
        longitudinal.groupby(
            ["study", "participant_id"],
            as_index=False,
            observed=False,
        )
        .agg(
            n_longitudinal_scans=("prediction_source_row", "nunique"),
            followup_years=("followup_years", "max"),
            followup_time_source=(
                "followup_time_source",
                lambda values: (
                    "Date"
                    if (values == "Date").any()
                    else (
                        "Age"
                        if (values == "Age").any()
                        else "Unavailable"
                    )
                ),
            ),
        )
    )

    study_followup = (
        participant_followup.groupby(
            "study",
            as_index=False,
            observed=False,
        )
        .agg(
            N_with_followup=("participant_id", "nunique"),
            N_with_at_least_2_scans=(
                "n_longitudinal_scans",
                lambda values: int((values >= 2).sum()),
            ),
            Followup_mean_years=("followup_years", "mean"),
            Followup_sd_years=("followup_years", "std"),
            Followup_median_years=("followup_years", "median"),
            Followup_q25_years=(
                "followup_years",
                lambda values: values.quantile(0.25),
            ),
            Followup_q75_years=(
                "followup_years",
                lambda values: values.quantile(0.75),
            ),
            Followup_min_years=("followup_years", "min"),
            Followup_max_years=("followup_years", "max"),
        )
    )

    study_followup["Followup_median_IQR"] = study_followup.apply(
        lambda row: (
            f"{row['Followup_median_years']:.2f} "
            f"({row['Followup_q25_years']:.2f}-"
            f"{row['Followup_q75_years']:.2f})"
        ),
        axis=1,
    )

    study_followup["Event_followup_years"] = study_followup.apply(
        lambda row: (
            f"0-{row['Followup_max_years']:.1f} years"
            if pd.notna(row["Followup_max_years"])
            else "NA"
        ),
        axis=1,
    )

    return participant_followup, study_followup

def read_sample_sex_file(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(
        path,
        sep="\t",
        low_memory=False,
    )

    columns = raw.columns.tolist()

    id_col = find_column(
        columns,
        ["PTID", "participant_id", "IID", "eid"],
        r"(^PTID$|participant.*id|^IID$|^eid$)",
        "sample participant ID",
    )

    study_col = find_column(
        columns,
        ["Study", "STUDY"],
        r"(^|_)study$",
        "sample Study",
    )

    sex_col = find_column(
        columns,
        ["Sex", "SEX", "sex"],
        r"^sex$",
        "sample Sex",
    )

    age_col = find_column(
        columns,
        ["Age", "AGE"],
        r"^age$",
        "sample Age",
        required=False,
    )

    date_col = find_column(
        columns,
        ["Date", "scan_date", "MRI_Date"],
        r"(^|_)date$",
        "sample Date",
        required=False,
    )

    sample = pd.DataFrame(
        {
            "sample_source_row": np.arange(
                raw.shape[0],
                dtype=int,
            ),
            "participant_id": clean_string(
                raw[id_col]
            ),
            "study": normalize_study(
                raw[study_col]
            ),
            "sex_sample": normalize_sex(
                raw[sex_col]
            ),
            "sample_age": (
                pd.to_numeric(
                    raw[age_col],
                    errors="coerce",
                )
                if age_col is not None
                else pd.Series(
                    np.nan,
                    index=raw.index,
                    dtype=float,
                )
            ),
            "sample_date": (
                parse_date(raw[date_col])
                if date_col is not None
                else pd.Series(
                    pd.NaT,
                    index=raw.index,
                )
            ),
        }
    )

    sample = sample.loc[
        sample["study"].isin(STUDY_ORDER)
        & sample["participant_id"].notna()
    ].copy()

    return sample


def match_sex_to_baseline(
    baseline: pd.DataFrame,
    sample: pd.DataFrame,
) -> pd.DataFrame:
    candidates = baseline.merge(
        sample,
        on=[
            "study",
            "participant_id",
        ],
        how="left",
        validate="one_to_many",
    )

    candidates["date_difference_days"] = (
        candidates["scan_date"]
        - candidates["sample_date"]
    ).abs().dt.total_seconds() / 86400.0

    candidates["age_difference_years"] = (
        candidates["age"]
        - candidates["sample_age"]
    ).abs()

    candidates["sex_match_priority"] = np.select(
        [
            candidates["date_difference_days"].notna(),
            candidates["age_difference_years"].notna(),
            candidates["sample_date"].notna(),
            candidates["sample_age"].notna(),
        ],
        [
            1,
            2,
            3,
            4,
        ],
        default=5,
    )

    candidates["sex_match_distance"] = np.select(
        [
            candidates["sex_match_priority"] == 1,
            candidates["sex_match_priority"] == 2,
            candidates["sex_match_priority"] == 3,
            candidates["sex_match_priority"] == 4,
        ],
        [
            candidates["date_difference_days"],
            candidates["age_difference_years"],
            candidates["sample_date"].map(
                lambda value: (
                    value.toordinal()
                    if pd.notna(value)
                    else np.nan
                )
            ),
            candidates["sample_age"],
        ],
        default=candidates["sample_source_row"].astype(float),
    )

    candidates = candidates.sort_values(
        [
            "study",
            "participant_id",
            "sex_match_priority",
            "sex_match_distance",
            "sample_date",
            "sample_age",
            "sample_source_row",
        ],
        kind="mergesort",
    )

    matched = (
        candidates.groupby(
            [
                "study",
                "participant_id",
            ],
            as_index=False,
            sort=False,
        )
        .head(1)
        .copy()
    )

    matched["sex"] = matched["sex_prediction"].fillna(
        matched["sex_sample"]
    )

    matched["sex_source"] = np.where(
        matched["sex_prediction"].notna(),
        "harmonized prediction file",
        np.where(
            matched["sex_sample"].notna(),
            "matched iSTAGING sample file",
            "missing",
        ),
    )

    return matched


def summarize_group(
    data: pd.DataFrame,
    group_columns: list[str],
) -> pd.DataFrame:
    summary = (
        data.groupby(
            group_columns,
            dropna=False,
            observed=False,
        )
        .agg(
            N=(
                "participant_id",
                "nunique",
            ),
            Age_mean=(
                "age",
                "mean",
            ),
            Age_sd=(
                "age",
                "std",
            ),
            Age_median=(
                "age",
                "median",
            ),
            Age_min=(
                "age",
                "min",
            ),
            Age_max=(
                "age",
                "max",
            ),
            Sex_nonmissing_n=(
                "sex",
                lambda values: values.notna().sum(),
            ),
            Female_n=(
                "sex",
                lambda values: (
                    values == "Female"
                ).sum(),
            ),
            Male_n=(
                "sex",
                lambda values: (
                    values == "Male"
                ).sum(),
            ),
        )
        .reset_index()
    )

    summary["Female_percent"] = np.where(
        summary["Sex_nonmissing_n"] > 0,
        100.0
        * summary["Female_n"]
        / summary["Sex_nonmissing_n"],
        np.nan,
    )

    summary["Age_mean_sd"] = summary.apply(
        lambda row: (
            f"{row['Age_mean']:.2f} \u00b1 "
            f"{row['Age_sd']:.2f}"
            if pd.notna(row["Age_mean"])
            and pd.notna(row["Age_sd"])
            else "NA"
        ),
        axis=1,
    )

    summary["Female_n_percent"] = summary.apply(
        lambda row: (
            f"{int(row['Female_n']):,}/"
            f"{row['Female_percent']:.1f}%"
            if pd.notna(row["Female_percent"])
            else "NA"
        ),
        axis=1,
    )

    return summary


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Summarize baseline age, sex, and sample size for "
            "AIBL, OASIS, and BLSA."
        )
    )

    parser.add_argument(
        "--prediction-file",
        default=(
            "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/"
            "adni_lepoch/results_external_longitudinal_ad_epoch_harmonized/"
            "external_5_studies_adni_brain_mri_ad_epoch_"
            "harmonized_scan_level_predictions.tsv"
        ),
        help=(
            "Harmonized external AD EPOCH scan-level prediction TSV. "
            "A project-specific default is provided so the script can be "
            "run directly in PyCharm."
        ),
    )

    parser.add_argument(
        "--sample-file",
        default=(
            "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/"
            "adni_lepoch/external_5_studies_istaging.tsv"
        ),
        help=(
            "Full external iSTAGING sample-information TSV, used as a "
            "fallback source for sex. A project-specific default is provided."
        ),
    )

    parser.add_argument(
        "--outdir",
        default=(
            "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/"
            "adni_lepoch/results_external_longitudinal_ad_epoch_comparison/"
            "baseline_demographics"
        ),
        help=(
            "Output directory. A project-specific default is provided so the "
            "script can be run directly in PyCharm."
        ),
    )

    parser.add_argument(
        "--prefix",
        default="AIBL_OASIS_BLSA_baseline_demographics",
    )

    args = parser.parse_args()

    prediction_file = Path(
        args.prediction_file
    )
    sample_file = Path(
        args.sample_file
    )
    outdir = Path(
        args.outdir
    )
    outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    if not prediction_file.exists():
        raise FileNotFoundError(
            f"Prediction file does not exist: {prediction_file}"
        )

    if not sample_file.exists():
        raise FileNotFoundError(
            f"Sample file does not exist: {sample_file}"
        )

    log("=" * 72)
    log("Baseline demographic summary for AIBL, OASIS, and BLSA")
    log(f"Prediction file: {prediction_file}")
    log(f"Sample file:     {sample_file}")
    log(f"Output folder:   {outdir}")
    log("=" * 72)

    prediction = read_prediction_file(
        prediction_file
    )

    baseline = select_one_baseline(
        prediction
    )

    participant_followup, study_followup = calculate_longitudinal_followup(
        prediction=prediction,
        baseline=baseline,
    )

    sample = read_sample_sex_file(
        sample_file
    )

    baseline_with_sex = match_sex_to_baseline(
        baseline,
        sample,
    )

    baseline_with_sex["study"] = pd.Categorical(
        baseline_with_sex["study"],
        categories=STUDY_ORDER,
        ordered=True,
    )

    baseline_with_sex = baseline_with_sex.sort_values(
        [
            "study",
            "site",
            "participant_id",
        ],
        kind="mergesort",
    )

    study_summary = summarize_group(
        baseline_with_sex,
        ["study"],
    )

    site_summary = summarize_group(
        baseline_with_sex,
        [
            "study",
            "site",
        ],
    )

    study_summary["study"] = study_summary[
        "study"
    ].astype("string")

    site_summary["study"] = site_summary[
        "study"
    ].astype("string")

    study_summary = study_summary.loc[
        study_summary["study"].isin(
            STUDY_ORDER
        )
    ].copy()

    study_followup["study"] = study_followup["study"].astype("string")

    study_summary["Country"] = study_summary["study"].map(
        STUDY_COUNTRY
    )

    study_summary = study_summary.merge(
        study_followup,
        on="study",
        how="left",
        validate="one_to_one",
    )

    # Put table-ready columns first.
    preferred_study_columns = [
        "study",
        "Country",
        "N",
        "Age_mean_sd",
        "Female_n_percent",
        "Event_followup_years",
        "Age_mean",
        "Age_sd",
        "Age_median",
        "Age_min",
        "Age_max",
        "Sex_nonmissing_n",
        "Female_n",
        "Female_percent",
        "Male_n",
        "N_with_followup",
        "N_with_at_least_2_scans",
        "Followup_mean_years",
        "Followup_sd_years",
        "Followup_median_years",
        "Followup_q25_years",
        "Followup_q75_years",
        "Followup_min_years",
        "Followup_max_years",
        "Followup_median_IQR",
    ]

    study_summary = study_summary[
        [
            column
            for column in preferred_study_columns
            if column in study_summary.columns
        ]
        + [
            column
            for column in study_summary.columns
            if column not in preferred_study_columns
        ]
    ]

    site_summary = site_summary.loc[
        site_summary["study"].isin(
            STUDY_ORDER
        )
    ].copy()

    study_file = outdir / (
        f"{args.prefix}_study_summary.tsv"
    )

    site_file = outdir / (
        f"{args.prefix}_site_summary.tsv"
    )

    selected_file = outdir / (
        f"{args.prefix}_selected_baseline_participants.tsv"
    )

    selection_qc_file = outdir / (
        f"{args.prefix}_baseline_selection_QC.tsv"
    )

    sex_source_file = outdir / (
        f"{args.prefix}_sex_source_QC.tsv"
    )

    participant_followup_file = outdir / (
        f"{args.prefix}_participant_followup.tsv"
    )

    study_followup_file = outdir / (
        f"{args.prefix}_study_followup_summary.tsv"
    )

    study_summary.to_csv(
        study_file,
        sep="\t",
        index=False,
    )

    site_summary.to_csv(
        site_file,
        sep="\t",
        index=False,
    )

    baseline_with_sex.to_csv(
        selected_file,
        sep="\t",
        index=False,
    )

    (
        baseline_with_sex.groupby(
            [
                "study",
                "baseline_selection_source",
            ],
            dropna=False,
            observed=False,
        )
        .size()
        .reset_index(
            name="N"
        )
        .to_csv(
            selection_qc_file,
            sep="\t",
            index=False,
        )
    )

    (
        baseline_with_sex.groupby(
            [
                "study",
                "sex_source",
            ],
            dropna=False,
            observed=False,
        )
        .size()
        .reset_index(
            name="N"
        )
        .to_csv(
            sex_source_file,
            sep="\t",
            index=False,
        )
    )

    participant_followup.to_csv(
        participant_followup_file,
        sep="\t",
        index=False,
    )

    study_followup.to_csv(
        study_followup_file,
        sep="\t",
        index=False,
    )

    log("")
    log("Study-level baseline demographics:")
    display_columns = [
        "study",
        "Country",
        "N",
        "Age_mean_sd",
        "Female_n_percent",
        "Event_followup_years",
        "Followup_median_IQR",
        "N_with_at_least_2_scans",
    ]
    log(
        study_summary[
            display_columns
        ].to_string(
            index=False
        )
    )

    log("")
    log(f"Wrote study summary:       {study_file}")
    log(f"Wrote site summary:        {site_file}")
    log(f"Wrote selected data:       {selected_file}")
    log(f"Wrote participant followup:{participant_followup_file}")
    log(f"Wrote study followup:      {study_followup_file}")
    log("=" * 72)


if __name__ == "__main__":
    main()