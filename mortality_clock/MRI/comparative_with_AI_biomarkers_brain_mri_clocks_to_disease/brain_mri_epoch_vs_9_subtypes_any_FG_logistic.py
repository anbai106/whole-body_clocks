#!/usr/bin/env python3
"""
Case-control logistic-regression comparison of the brain MRI mortality EPOCH
versus nine AI-derived disease-subtype biomarkers for future ICD-10 F/G disease.

Why this analysis?
------------------
The mortality EPOCH was trained with a Cox survival objective. A reviewer could
therefore argue that a survival-model comparison favors EPOCH because the
downstream endpoint is also analyzed with survival methods. This script removes
time-to-event modeling from the downstream comparison and asks a simpler binary
question: can each biomarker distinguish participants who subsequently develop
an ICD-10 F/G disorder from participants who remain F/G-disease-free?

The script supports two complementary binary endpoints:

1) observed_followup
   - case: first recorded inpatient ICD-10 F* or G* diagnosis occurs after MRI
           and on/before censoring;
   - control: no recorded F/G diagnosis on/before censoring.
   This most directly mirrors the case/non-case groups used in the survival
   analysis, but follow-up duration can differ among controls.

2) fixed_horizon
   - case: first recorded F/G diagnosis occurs within a fixed number of years
           after MRI;
   - control: remains F/G-disease-free through that same horizon and is observed
              for at least the full horizon.
   This removes differential follow-up opportunity and is therefore an important
   sensitivity analysis for a binary logistic design.

Primary fairness rules
----------------------
- Use the held-out EPOCH test split by default.
- Use brain MRI date as the common index date.
- Exclude any participant with an F/G diagnosis on/before MRI.
- Strictly exclude unresolved F/G codes with missing paired 41280 dates by default.
- Use one common complete-case population for EPOCH + all nine subtype scores.
- Standardize all biomarkers in the exact same analysis sample.
- Use the same covariates for every adjusted model:
      age at imaging, sex, smoking, BMI.
- Compare both:
      (a) biomarker-only logistic models, and
      (b) covariate-adjusted logistic models.
- For predictive discrimination, use repeated stratified K-fold out-of-fold
  predictions, not apparent/in-sample AUC.
- Use paired participant-level bootstrap CIs for AUC and AUPRC differences.

Biomarkers
----------
Mortality EPOCH:
  brain_mri_mortality_clock_acceleration_z

AI disease subtypes:
  AD1, AD2,
  ASD1, ASD2, ASD3,
  LLD1, LLD2,
  SCZ1, SCZ2

Outputs
-------
cohort_qc.tsv
sample_flow.tsv
analysis_common_sample.tsv
marker_standardization.tsv
single_marker_unadjusted.tsv
single_marker_adjusted.tsv
epoch_vs_subtype_adjusted_joint.tsv
combined_adjusted_models.tsv
cv_model_metrics.tsv
cv_out_of_fold_predictions.tsv
bootstrap_discrimination_comparisons.tsv
run_metadata.txt
"""

import argparse
import os
import sys
import warnings
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy.stats import chi2

try:
    import statsmodels.api as sm
except ImportError as exc:
    raise ImportError(
        "statsmodels is required. Install it in the survival environment "
        "(e.g., conda install statsmodels)."
    ) from exc

try:
    from sklearn.metrics import (
        roc_auc_score,
        average_precision_score,
        brier_score_loss,
        log_loss,
    )
    from sklearn.model_selection import RepeatedStratifiedKFold
except ImportError as exc:
    raise ImportError(
        "scikit-learn is required. Install it in the survival environment."
    ) from exc

warnings.filterwarnings("ignore")

SUBTYPES: List[str] = [
    "AD1", "AD2",
    "ASD1", "ASD2", "ASD3",
    "LLD1", "LLD2",
    "SCZ1", "SCZ2",
]

EPOCH_RAW = "brain_mri_mortality_clock_acceleration_z"
EPOCH_NAME = "Brain_MRI_mortality_EPOCH"

BASE_COVARS: List[str] = [
    "Age_imaging",
    "Sex",
    "Smoking",
    "BMI",
]


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Fair logistic-regression comparison of brain MRI mortality EPOCH "
            "versus nine AI-derived disease subtypes for future ICD-10 F/G disease."
        )
    )

    p.add_argument(
        "--epoch_tsv",
        default=(
            "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/"
            "brain_mri_mortality_clock/"
            "brain_mri_mortality_clock_predictions.tsv"
        ),
    )
    p.add_argument(
        "--subtype_tsv",
        default=(
            "/cbica/projects/MULTI/processed/UKBB/"
            "derived_AI_biomakers_across_projects/"
            "UKBB_487894_participant_58_biomarker_matched_ID.tsv"
        ),
    )
    p.add_argument(
        "--subtype_id_col",
        default="id_upenn",
    )
    p.add_argument(
        "--icd10_csv",
        default=(
            "/cbica/home/wenju/Reproducibile_paper/BrainEye/data/"
            "UKBB_fullsample_ICD10.csv"
        ),
    )
    p.add_argument(
        "--umel_death_xlsx",
        default=(
            "/cbica/home/wenju/Dataset/UKBB_UMelbourne/"
            "Death_related_var_from_Ye.xlsx"
        ),
    )
    p.add_argument(
        "--umel_match_csv",
        default=(
            "/cbica/home/wenju/Dataset/UKBB_UMelbourne/"
            "UKB_UMelbourne_vs_Penn_match_key.csv"
        ),
    )
    p.add_argument(
        "--cov_tsv",
        default=(
            "/cbica/home/wenju/Reproducibile_paper/PRS_UKBB/"
            "prediction/data/UKBB_fullsample_covariate.csv"
        ),
    )
    p.add_argument(
        "--output_dir",
        required=True,
    )
    p.add_argument(
        "--split",
        choices=["test", "train", "all"],
        default="test",
        help=(
            "Use test for the primary analysis so the mortality EPOCH is evaluated "
            "only in participants held out from its own training."
        ),
    )
    p.add_argument(
        "--endpoint_mode",
        choices=["observed_followup", "fixed_horizon"],
        default="observed_followup",
    )
    p.add_argument(
        "--horizon_years",
        type=float,
        default=3.0,
        help=(
            "Used only for --endpoint_mode fixed_horizon. Controls must be "
            "observed and F/G-disease-free through this horizon."
        ),
    )
    p.add_argument(
        "--admin_censor_date",
        default="2022-11-30",
    )
    p.add_argument(
        "--min_cases",
        type=int,
        default=20,
    )
    p.add_argument(
        "--min_controls",
        type=int,
        default=20,
    )
    p.add_argument(
        "--cv_folds",
        type=int,
        default=5,
    )
    p.add_argument(
        "--cv_repeats",
        type=int,
        default=5,
    )
    p.add_argument(
        "--bootstrap",
        type=int,
        default=500,
        help="Paired participant bootstrap replicates for AUC/AUPRC differences.",
    )
    p.add_argument(
        "--seed",
        type=int,
        default=20260818,
    )
    p.add_argument(
        "--allow_unresolved_fg_dates",
        action="store_true",
        help=(
            "By default, an F/G code whose paired first-diagnosis date is missing "
            "causes exclusion because prevalent versus incident status is unresolved."
        ),
    )

    return p.parse_args()


def norm_id(x: pd.Series) -> pd.Series:
    return pd.to_numeric(
        x,
        errors="coerce",
    ).astype("Int64")


def parse_date(x: pd.Series) -> pd.Series:
    s = x.copy().replace(
        [
            0,
            0.0,
            "0",
            "0.0",
            "",
            "NA",
            "NaN",
            "nan",
            "None",
            "-1",
            -1,
        ],
        np.nan,
    )

    out = pd.to_datetime(
        s,
        errors="coerce",
    )

    num = pd.to_numeric(
        s,
        errors="coerce",
    )

    excel_mask = num.between(
        20000,
        60000,
    )

    if excel_mask.any():
        excel = pd.to_datetime(
            num,
            unit="D",
            origin="1899-12-30",
            errors="coerce",
        )
        out = out.where(
            ~excel_mask,
            excel,
        )

    return out


def bool_numpy(mask) -> np.ndarray:
    """Safely convert pandas nullable/object masks to strict NumPy bool arrays."""
    if isinstance(mask, pd.Series):
        return (
            mask
            .fillna(False)
            .astype(bool)
            .to_numpy(dtype=bool)
        )

    arr = np.asarray(mask)

    if arr.dtype == bool:
        return arr

    return (
        pd.Series(arr)
        .fillna(False)
        .astype(bool)
        .to_numpy(dtype=bool)
    )


def read_epoch(args) -> Tuple[pd.DataFrame, Dict[str, object]]:
    d = pd.read_csv(
        args.epoch_tsv,
        sep="\t",
    )

    required = [
        "participant_id",
        "imaging_date",
        "age_at_imaging",
        EPOCH_RAW,
    ]

    missing = [
        c for c in required
        if c not in d.columns
    ]

    if missing:
        raise ValueError(
            f"Missing EPOCH columns: {missing}"
        )

    qc = {
        "epoch_rows_raw": int(len(d)),
    }

    if args.split != "all":
        if "split" not in d.columns:
            raise ValueError(
                f"--split {args.split} requested, but EPOCH file has no split column."
            )

        d = d[
            d["split"]
            .astype(str)
            .str.lower()
            .eq(args.split.lower())
        ].copy()

    qc["epoch_rows_after_split"] = int(len(d))

    keep = required + [
        c for c in [
            "split",
            "sex",
            "death_date",
            "admin_censor_date",
            "end_date",
        ]
        if c in d.columns
    ]

    d = d[keep].copy()

    d["participant_id"] = norm_id(
        d["participant_id"]
    )

    d = d[
        d["participant_id"].notna()
    ].copy()

    d["imaging_date"] = parse_date(
        d["imaging_date"]
    )

    d["Age_imaging"] = pd.to_numeric(
        d["age_at_imaging"],
        errors="coerce",
    )

    d[EPOCH_RAW] = pd.to_numeric(
        d[EPOCH_RAW],
        errors="coerce",
    )

    admin_default = pd.Timestamp(
        args.admin_censor_date
    )

    if "admin_censor_date" in d.columns:
        d["admin_censor_date"] = parse_date(
            d["admin_censor_date"]
        ).fillna(admin_default)
    else:
        d["admin_censor_date"] = admin_default

    if "death_date" in d.columns:
        d["death_date"] = parse_date(
            d["death_date"]
        )
    else:
        d["death_date"] = pd.NaT

    fallback_end = d[
        "admin_censor_date"
    ].copy()

    death_before_admin = (
        d["death_date"].notna()
        & (
            d["death_date"]
            < fallback_end
        )
    )

    fallback_end.loc[
        death_before_admin
    ] = d.loc[
        death_before_admin,
        "death_date",
    ]

    if "end_date" in d.columns:
        d["censor_date"] = parse_date(
            d["end_date"]
        ).fillna(
            fallback_end
        )
    else:
        d["censor_date"] = fallback_end

    if "sex" in d.columns:
        s = (
            d["sex"]
            .astype(str)
            .str.strip()
            .str.lower()
        )

        d["Sex_epoch_fallback"] = np.where(
            s.isin(
                [
                    "male",
                    "m",
                    "1",
                    "1.0",
                ]
            ),
            1.0,
            np.where(
                s.isin(
                    [
                        "female",
                        "f",
                        "0",
                        "0.0",
                    ]
                ),
                0.0,
                np.nan,
            ),
        )
    else:
        d["Sex_epoch_fallback"] = np.nan

    d = d[
        d["imaging_date"].notna()
        & d["censor_date"].notna()
        & (
            d["censor_date"]
            > d["imaging_date"]
        )
    ].copy()

    qc["epoch_duplicate_rows_before_dedup"] = int(
        d["participant_id"]
        .duplicated(keep=False)
        .sum()
    )

    d = (
        d
        .sort_values(
            [
                "participant_id",
                "imaging_date",
            ]
        )
        .drop_duplicates(
            "participant_id",
            keep="first",
        )
        .copy()
    )

    qc["epoch_unique_ids_usable"] = int(
        d["participant_id"].nunique()
    )

    return d, qc


def read_subtypes(
    args,
    epoch_ids: set,
) -> Tuple[pd.DataFrame, Dict[str, object]]:
    d = pd.read_csv(
        args.subtype_tsv,
        sep="\t",
    )

    if args.subtype_id_col not in d.columns:
        raise ValueError(
            f"Missing subtype ID column: {args.subtype_id_col}"
        )

    missing = [
        c for c in SUBTYPES
        if c not in d.columns
    ]

    if missing:
        raise ValueError(
            f"Missing subtype columns: {missing}"
        )

    d = d[
        [
            args.subtype_id_col
        ] + SUBTYPES
    ].copy()

    d = d.rename(
        columns={
            args.subtype_id_col: "participant_id"
        }
    )

    d["participant_id"] = norm_id(
        d["participant_id"]
    )

    d = d[
        d["participant_id"].notna()
    ].copy()

    for c in SUBTYPES:
        d[c] = pd.to_numeric(
            d[c],
            errors="coerce",
        )

    qc = {
        "subtype_rows_raw": int(len(d)),
        "subtype_duplicate_rows_before_dedup": int(
            d["participant_id"]
            .duplicated(keep=False)
            .sum()
        ),
    }

    d["_n_observed"] = (
        d[SUBTYPES]
        .notna()
        .sum(axis=1)
    )

    d = (
        d
        .sort_values(
            [
                "participant_id",
                "_n_observed",
            ],
            ascending=[
                True,
                False,
            ],
        )
        .drop_duplicates(
            "participant_id",
            keep="first",
        )
        .drop(
            columns="_n_observed"
        )
    )

    qc["subtype_unique_ids"] = int(
        d["participant_id"].nunique()
    )

    qc["subtype_overlap_with_epoch"] = int(
        d["participant_id"]
        .isin(epoch_ids)
        .sum()
    )

    return d, qc


def read_covariates(args) -> pd.DataFrame:
    d = pd.read_csv(
        args.cov_tsv
    )

    if "eid" not in d.columns:
        raise ValueError(
            "Covariate file must contain eid."
        )

    sex_col = next(
        (
            c for c in [
                "sex_f31_0_0",
                "genetic_sex_f22001_0_0",
                "Sex",
                "sex",
            ]
            if c in d.columns
        ),
        None,
    )

    smoking_col = next(
        (
            c for c in [
                "smoking_status_f20116_0_0",
                "Smoking",
                "smoking",
            ]
            if c in d.columns
        ),
        None,
    )

    bmi_col = next(
        (
            c for c in [
                "body_mass_index_bmi_f23104_0_0",
                "BMI",
                "bmi",
            ]
            if c in d.columns
        ),
        None,
    )

    keep = ["eid"] + [
        c for c in [
            sex_col,
            smoking_col,
            bmi_col,
        ]
        if c is not None
    ]

    d = d[
        keep
    ].copy()

    rename = {
        "eid": "participant_id",
    }

    if sex_col is not None:
        rename[sex_col] = "Sex"

    if smoking_col is not None:
        rename[smoking_col] = "Smoking"

    if bmi_col is not None:
        rename[bmi_col] = "BMI"

    d = d.rename(
        columns=rename
    )

    d["participant_id"] = norm_id(
        d["participant_id"]
    )

    d = d[
        d["participant_id"].notna()
    ].copy()

    for c in [
        "Sex",
        "Smoking",
        "BMI",
    ]:
        if c not in d.columns:
            d[c] = np.nan

        d[c] = pd.to_numeric(
            d[c],
            errors="coerce",
        )

    return (
        d[
            [
                "participant_id",
                "Sex",
                "Smoking",
                "BMI",
            ]
        ]
        .drop_duplicates(
            "participant_id",
            keep="first",
        )
        .copy()
    )


def read_icd_codes(
    args,
    ids: set,
) -> Tuple[pd.DataFrame, List[str]]:
    d = pd.read_csv(
        args.icd10_csv
    )

    if "eid" not in d.columns:
        raise ValueError(
            "ICD file must contain eid."
        )

    code_cols = [
        c for c in d.columns
        if c.startswith(
            "diagnoses_icd10_f41270_"
        )
    ]

    if not code_cols:
        raise ValueError(
            "No diagnoses_icd10_f41270_* columns found."
        )

    d = d[
        [
            "eid"
        ] + code_cols
    ].copy()

    d = d.rename(
        columns={
            "eid": "participant_id"
        }
    )

    d["participant_id"] = norm_id(
        d["participant_id"]
    )

    d = d[
        d["participant_id"].isin(ids)
    ].copy()

    d = d.drop_duplicates(
        "participant_id",
        keep="first",
    )

    return d, code_cols


def read_icd_dates(
    args,
    ids: set,
) -> pd.DataFrame:
    death = pd.read_excel(
        args.umel_death_xlsx,
        engine="openpyxl",
    )

    match = pd.read_csv(
        args.umel_match_csv
    )

    if "eid" not in death.columns:
        raise ValueError(
            "UMelbourne diagnosis-date file must contain eid."
        )

    if not {
        "id",
        "id_upenn",
    }.issubset(
        match.columns
    ):
        raise ValueError(
            "UMelbourne match key must contain id and id_upenn."
        )

    date_cols = [
        c for c in death.columns
        if c.startswith("41280-")
    ]

    if not date_cols:
        raise ValueError(
            "No 41280-* diagnosis date columns found."
        )

    death = death[
        [
            "eid"
        ] + date_cols
    ].copy()

    death = death.rename(
        columns={
            "eid": "participant_id_umel"
        }
    )

    match = match[
        [
            "id",
            "id_upenn",
        ]
    ].copy()

    match = match.rename(
        columns={
            "id": "participant_id_umel",
            "id_upenn": "participant_id",
        }
    )

    death["participant_id_umel"] = norm_id(
        death["participant_id_umel"]
    )

    match["participant_id_umel"] = norm_id(
        match["participant_id_umel"]
    )

    match["participant_id"] = norm_id(
        match["participant_id"]
    )

    match = match[
        match["participant_id"].isin(ids)
    ].copy()

    d = match.merge(
        death,
        on="participant_id_umel",
        how="inner",
    )

    d = d[
        [
            "participant_id"
        ] + date_cols
    ].copy()

    d = d.rename(
        columns={
            c: (
                c
                .replace(
                    "41280-",
                    "",
                )
                .replace(
                    ".",
                    "_",
                )
            )
            for c in date_cols
        }
    )

    return d.drop_duplicates(
        "participant_id",
        keep="first",
    )


def derive_fg_dates(
    base: pd.DataFrame,
    diag: pd.DataFrame,
    code_cols: List[str],
    dates: pd.DataFrame,
) -> Tuple[pd.DataFrame, Dict[str, object]]:
    """Derive earliest valid inpatient ICD-10 F/G diagnosis date."""

    x = (
        base
        .merge(
            diag,
            on="participant_id",
            how="left",
        )
        .merge(
            dates,
            on="participant_id",
            how="left",
        )
    )

    earliest = pd.Series(
        pd.NaT,
        index=x.index,
        dtype="datetime64[ns]",
    )

    earliest_code = pd.Series(
        pd.NA,
        index=x.index,
        dtype="string",
    )

    fg_count = np.zeros(
        len(x),
        dtype=np.int32,
    )

    unresolved = np.zeros(
        len(x),
        dtype=bool,
    )

    matched_slots = 0
    unmatched_slots = 0

    for col in code_cols:
        suffix = col.replace(
            "diagnoses_icd10_f41270_",
            "",
        )

        codes = (
            x[col]
            .astype("string")
            .str.upper()
            .str.strip()
            .str.replace(
                r"[^A-Z0-9\.]",
                "",
                regex=True,
            )
        )

        is_fg_np = bool_numpy(
            codes.str.match(
                r"^[FG][0-9]",
                na=False,
            )
        )

        if not is_fg_np.any():
            continue

        fg_count += is_fg_np.astype(
            np.int32
        )

        if suffix not in x.columns:
            unresolved |= is_fg_np
            unmatched_slots += 1
            continue

        matched_slots += 1

        dt = parse_date(
            x[suffix]
        )

        missing_date_np = bool_numpy(
            dt.isna()
        )

        present_date_np = bool_numpy(
            dt.notna()
        )

        unresolved |= (
            is_fg_np
            & missing_date_np
        )

        valid_np = (
            is_fg_np
            & present_date_np
        )

        if not valid_np.any():
            continue

        valid = pd.Series(
            valid_np,
            index=x.index,
            dtype=bool,
        )

        earlier = valid & (
            earliest.isna()
            | (
                dt
                < earliest
            )
        )

        earlier = pd.Series(
            bool_numpy(
                earlier
            ),
            index=x.index,
            dtype=bool,
        )

        earliest.loc[
            earlier
        ] = dt.loc[
            earlier
        ]

        earliest_code.loc[
            earlier
        ] = codes.loc[
            earlier
        ]

    out = base.copy()

    out["earliest_FG_date"] = earliest.values
    out["earliest_FG_code"] = earliest_code.values
    out["n_FG_codes_recorded"] = fg_count
    out[
        "has_FG_code_missing_paired_date"
    ] = unresolved

    qc = {
        "icd_41270_slots": int(
            len(code_cols)
        ),
        "icd_slots_with_matching_41280": int(
            matched_slots
        ),
        "icd_slots_without_matching_41280": int(
            unmatched_slots
        ),
        "participants_with_any_FG_code": int(
            (
                out[
                    "n_FG_codes_recorded"
                ]
                > 0
            ).sum()
        ),
        "participants_with_unresolved_FG_date": int(
            out[
                "has_FG_code_missing_paired_date"
            ].sum()
        ),
    }

    return out, qc


def define_binary_endpoint(
    d: pd.DataFrame,
    args,
) -> Tuple[pd.DataFrame, Dict[str, object]]:
    out = d.copy()

    out["prevalent_FG_at_imaging"] = (
        out[
            "earliest_FG_date"
        ].notna()
        & (
            out[
                "earliest_FG_date"
            ]
            <= out[
                "imaging_date"
            ]
        )
    )

    if args.allow_unresolved_fg_dates:
        out[
            "exclude_unresolved_FG_date"
        ] = False
    else:
        out[
            "exclude_unresolved_FG_date"
        ] = out[
            "has_FG_code_missing_paired_date"
        ].astype(bool)

    basic_eligible = (
        ~out[
            "prevalent_FG_at_imaging"
        ].astype(bool)
        & ~out[
            "exclude_unresolved_FG_date"
        ].astype(bool)
    )

    incident_after_mri_before_censor = (
        out[
            "earliest_FG_date"
        ].notna()
        & (
            out[
                "earliest_FG_date"
            ]
            > out[
                "imaging_date"
            ]
        )
        & (
            out[
                "earliest_FG_date"
            ]
            <= out[
                "censor_date"
            ]
        )
    )

    if args.endpoint_mode == "observed_followup":
        out["case"] = (
            basic_eligible
            & incident_after_mri_before_censor
        ).astype(int)

        out["binary_endpoint_eligible"] = (
            basic_eligible
            & (
                out[
                    "censor_date"
                ]
                > out[
                    "imaging_date"
                ]
            )
        )

        out["followup_years_for_binary_endpoint"] = (
            (
                out[
                    "censor_date"
                ]
                - out[
                    "imaging_date"
                ]
            ).dt.days
            / 365.25
        )

    elif args.endpoint_mode == "fixed_horizon":
        if (
            not np.isfinite(
                args.horizon_years
            )
            or args.horizon_years <= 0
        ):
            raise ValueError(
                "--horizon_years must be > 0 for fixed_horizon."
            )

        horizon_days = int(
            round(
                args.horizon_years
                * 365.25
            )
        )

        out["horizon_date"] = (
            out[
                "imaging_date"
            ]
            + pd.to_timedelta(
                horizon_days,
                unit="D",
            )
        )

        event_by_horizon = (
            incident_after_mri_before_censor
            & (
                out[
                    "earliest_FG_date"
                ]
                <= out[
                    "horizon_date"
                ]
            )
        )

        fully_observed_through_horizon = (
            out[
                "censor_date"
            ]
            >= out[
                "horizon_date"
            ]
        )

        # Cases do not need to remain observable after the event.
        # Controls must be observable through the full fixed horizon.
        eligible_case = (
            basic_eligible
            & event_by_horizon
        )

        eligible_control = (
            basic_eligible
            & ~event_by_horizon
            & fully_observed_through_horizon
        )

        out[
            "binary_endpoint_eligible"
        ] = (
            eligible_case
            | eligible_control
        )

        out["case"] = (
            eligible_case
        ).astype(int)

        out[
            "followup_years_for_binary_endpoint"
        ] = args.horizon_years

    else:
        raise ValueError(
            f"Unknown endpoint_mode: {args.endpoint_mode}"
        )

    qc = {
        "participants_prevalent_FG_at_imaging": int(
            out[
                "prevalent_FG_at_imaging"
            ].sum()
        ),
        "participants_excluded_unresolved_FG_date": int(
            out[
                "exclude_unresolved_FG_date"
            ].sum()
        ),
        "binary_endpoint_eligible_N_before_complete_case": int(
            out[
                "binary_endpoint_eligible"
            ].sum()
        ),
        "binary_endpoint_cases_before_complete_case": int(
            out.loc[
                out[
                    "binary_endpoint_eligible"
                ],
                "case",
            ].sum()
        ),
    }

    out = out[
        out[
            "binary_endpoint_eligible"
        ]
    ].copy()

    return out, qc


def build_common_sample(
    endpoint: pd.DataFrame,
    subtypes: pd.DataFrame,
    cov: pd.DataFrame,
) -> Tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
]:
    ep = endpoint[
        [
            c for c in endpoint.columns
            if c not in SUBTYPES
        ]
    ].copy()

    d = (
        ep
        .merge(
            subtypes,
            on="participant_id",
            how="inner",
        )
        .merge(
            cov,
            on="participant_id",
            how="left",
        )
    )

    n_before_cc = int(
        len(d)
    )

    if (
        "Sex_epoch_fallback"
        in d.columns
    ):
        d["Sex"] = d[
            "Sex"
        ].where(
            d[
                "Sex"
            ].notna(),
            d[
                "Sex_epoch_fallback"
            ],
        )

    needed = [
        "participant_id",
        "case",
    ] + BASE_COVARS + [
        EPOCH_RAW
    ] + SUBTYPES

    for c in needed:
        if c != "participant_id":
            d[c] = pd.to_numeric(
                d[c],
                errors="coerce",
            )

    d = (
        d
        .replace(
            [
                np.inf,
                -np.inf,
            ],
            np.nan,
        )
        .dropna(
            subset=needed
        )
        .copy()
    )

    if d["case"].nunique() < 2:
        raise RuntimeError(
            "Common complete-case sample does not contain both cases and controls."
        )

    std_rows = []

    for raw in [
        EPOCH_RAW
    ] + SUBTYPES:
        zcol = (
            "EPOCH_z"
            if raw == EPOCH_RAW
            else f"{raw}_z"
        )

        vals = pd.to_numeric(
            d[raw],
            errors="coerce",
        )

        mean = float(
            vals.mean()
        )

        sd = float(
            vals.std(
                ddof=1
            )
        )

        if (
            not np.isfinite(sd)
            or sd <= 0
        ):
            raise ValueError(
                f"Cannot standardize {raw}; SD={sd}"
            )

        d[zcol] = (
            vals
            - mean
        ) / sd

        std_rows.append(
            {
                "marker": (
                    EPOCH_NAME
                    if raw == EPOCH_RAW
                    else raw
                ),
                "raw_column": raw,
                "z_column": zcol,
                "mean_in_common_sample": mean,
                "sd_in_common_sample": sd,
            }
        )

    flow = pd.DataFrame(
        [
            {
                "stage": "Binary-endpoint eligible before biomarker/covariate complete-case filter",
                "N": int(
                    len(endpoint)
                ),
                "N_case": int(
                    endpoint[
                        "case"
                    ].sum()
                ),
                "N_control": int(
                    (
                        endpoint[
                            "case"
                        ]
                        == 0
                    ).sum()
                ),
            },
            {
                "stage": "After EPOCH-subtype-covariate merge, before complete-case filter",
                "N": n_before_cc,
                "N_case": np.nan,
                "N_control": np.nan,
            },
            {
                "stage": "Strict common complete-case logistic sample",
                "N": int(
                    len(d)
                ),
                "N_case": int(
                    d[
                        "case"
                    ].sum()
                ),
                "N_control": int(
                    (
                        d[
                            "case"
                        ]
                        == 0
                    ).sum()
                ),
            },
        ]
    )

    return (
        d,
        pd.DataFrame(
            std_rows
        ),
        flow,
    )


def design_matrix(
    d: pd.DataFrame,
    predictors: List[str],
) -> pd.DataFrame:
    x = d[
        predictors
    ].copy()

    for c in predictors:
        x[c] = pd.to_numeric(
            x[c],
            errors="coerce",
        )

    x = sm.add_constant(
        x,
        has_constant="add",
    )

    return x


def fit_glm(
    d: pd.DataFrame,
    predictors: List[str],
    label: str,
):
    y = pd.to_numeric(
        d["case"],
        errors="coerce",
    ).astype(int)

    x = design_matrix(
        d,
        predictors,
    )

    if (
        x.isna().any().any()
        or y.isna().any()
    ):
        raise ValueError(
            f"{label}: missing values entered logistic model."
        )

    try:
        model = sm.GLM(
            y,
            x,
            family=sm.families.Binomial(),
        )

        fit = model.fit(
            maxiter=200,
            disp=0,
        )

        return fit

    except Exception as exc:
        raise RuntimeError(
            f"{label}: logistic GLM failed: {exc}"
        ) from exc


def extract_or(
    fit,
    var: str,
) -> Dict[str, float]:
    beta = float(
        fit.params.loc[
            var
        ]
    )

    se = float(
        fit.bse.loc[
            var
        ]
    )

    p = float(
        fit.pvalues.loc[
            var
        ]
    )

    return {
        "beta": beta,
        "se": se,
        "or": float(
            np.exp(beta)
        ),
        "ci_low": float(
            np.exp(
                beta
                - 1.96
                * se
            )
        ),
        "ci_high": float(
            np.exp(
                beta
                + 1.96
                * se
            )
        ),
        "p": p,
    }


def lrt(
    full,
    reduced,
    df_diff: int,
) -> Tuple[float, float]:
    stat = 2.0 * (
        float(
            full.llf
        )
        - float(
            reduced.llf
        )
    )

    if (
        not np.isfinite(stat)
        or stat < 0
    ):
        return (
            np.nan,
            np.nan,
        )

    return (
        float(stat),
        float(
            chi2.sf(
                stat,
                df_diff,
            )
        ),
    )


def single_marker_inference(
    d: pd.DataFrame,
    adjusted: bool,
) -> Tuple[pd.DataFrame, Dict[str, object]]:
    rows = []

    base_predictors = (
        BASE_COVARS
        if adjusted
        else []
    )

    base_fit = fit_glm(
        d,
        base_predictors,
        (
            "Adjusted base"
            if adjusted
            else "Intercept-only base"
        ),
    )

    models = {
        "BASE": base_fit,
    }

    marker_specs = [
        (
            EPOCH_NAME,
            "EPOCH_z",
        )
    ] + [
        (
            s,
            f"{s}_z",
        )
        for s in SUBTYPES
    ]

    for marker_name, zcol in marker_specs:
        predictors = (
            base_predictors
            + [
                zcol
            ]
        )

        fit = fit_glm(
            d,
            predictors,
            marker_name,
        )

        stats = extract_or(
            fit,
            zcol,
        )

        lr_stat, lr_p = lrt(
            fit,
            base_fit,
            1,
        )

        rows.append(
            {
                "analysis": (
                    "covariate_adjusted"
                    if adjusted
                    else "biomarker_only"
                ),
                "marker": marker_name,
                "N": int(
                    len(d)
                ),
                "N_case": int(
                    d[
                        "case"
                    ].sum()
                ),
                "N_control": int(
                    (
                        d[
                            "case"
                        ]
                        == 0
                    ).sum()
                ),
                "OR_per_1SD": stats[
                    "or"
                ],
                "OR_CI_low": stats[
                    "ci_low"
                ],
                "OR_CI_high": stats[
                    "ci_high"
                ],
                "beta": stats[
                    "beta"
                ],
                "se": stats[
                    "se"
                ],
                "p_marker": stats[
                    "p"
                ],
                "lrt_chisq_vs_base": lr_stat,
                "lrt_p_vs_base": lr_p,
                "AIC": float(
                    fit.aic
                ),
                "log_likelihood": float(
                    fit.llf
                ),
            }
        )

        models[
            marker_name
        ] = fit

    return (
        pd.DataFrame(
            rows
        ),
        models,
    )


def adjusted_pairwise_inference(
    d: pd.DataFrame,
    adjusted_single_models: Dict[str, object],
) -> pd.DataFrame:
    rows = []

    epoch_single = adjusted_single_models[
        EPOCH_NAME
    ]

    for subtype in SUBTYPES:
        subtype_single = adjusted_single_models[
            subtype
        ]

        fit = fit_glm(
            d,
            BASE_COVARS
            + [
                "EPOCH_z",
                f"{subtype}_z",
            ],
            f"Adjusted EPOCH + {subtype}",
        )

        epoch_stats = extract_or(
            fit,
            "EPOCH_z",
        )

        subtype_stats = extract_or(
            fit,
            f"{subtype}_z",
        )

        stat_epoch, p_epoch = lrt(
            fit,
            subtype_single,
            1,
        )

        stat_subtype, p_subtype = lrt(
            fit,
            epoch_single,
            1,
        )

        pearson = float(
            d[
                [
                    "EPOCH_z",
                    f"{subtype}_z",
                ]
            ]
            .corr()
            .iloc[
                0,
                1,
            ]
        )

        rows.append(
            {
                "subtype": subtype,
                "N": int(
                    len(d)
                ),
                "N_case": int(
                    d[
                        "case"
                    ].sum()
                ),
                "pearson_EPOCH_subtype": pearson,
                "EPOCH_conditional_OR": epoch_stats[
                    "or"
                ],
                "EPOCH_conditional_CI_low": epoch_stats[
                    "ci_low"
                ],
                "EPOCH_conditional_CI_high": epoch_stats[
                    "ci_high"
                ],
                "EPOCH_conditional_p": epoch_stats[
                    "p"
                ],
                "subtype_conditional_OR": subtype_stats[
                    "or"
                ],
                "subtype_conditional_CI_low": subtype_stats[
                    "ci_low"
                ],
                "subtype_conditional_CI_high": subtype_stats[
                    "ci_high"
                ],
                "subtype_conditional_p": subtype_stats[
                    "p"
                ],
                "lrt_chisq_EPOCH_beyond_subtype": stat_epoch,
                "lrt_p_EPOCH_beyond_subtype": p_epoch,
                "lrt_chisq_subtype_beyond_EPOCH": stat_subtype,
                "lrt_p_subtype_beyond_EPOCH": p_subtype,
                "joint_AIC": float(
                    fit.aic
                ),
            }
        )

    return pd.DataFrame(
        rows
    )


def combined_adjusted_inference(
    d: pd.DataFrame,
    adjusted_single_models: Dict[str, object],
) -> pd.DataFrame:
    subtype_z = [
        f"{s}_z"
        for s in SUBTYPES
    ]

    base_fit = adjusted_single_models[
        "BASE"
    ]

    epoch_fit = adjusted_single_models[
        EPOCH_NAME
    ]

    all9_fit = fit_glm(
        d,
        BASE_COVARS
        + subtype_z,
        "Adjusted all 9 subtypes",
    )

    all10_fit = fit_glm(
        d,
        BASE_COVARS
        + subtype_z
        + [
            "EPOCH_z"
        ],
        "Adjusted all 9 subtypes + EPOCH",
    )

    stat_all9_vs_base, p_all9_vs_base = lrt(
        all9_fit,
        base_fit,
        9,
    )

    stat_epoch_beyond_all9, p_epoch_beyond_all9 = lrt(
        all10_fit,
        all9_fit,
        1,
    )

    stat_all9_beyond_epoch, p_all9_beyond_epoch = lrt(
        all10_fit,
        epoch_fit,
        9,
    )

    epoch_conditional = extract_or(
        all10_fit,
        "EPOCH_z",
    )

    return pd.DataFrame(
        [
            {
                "N": int(
                    len(d)
                ),
                "N_case": int(
                    d[
                        "case"
                    ].sum()
                ),
                "N_control": int(
                    (
                        d[
                            "case"
                        ]
                        == 0
                    ).sum()
                ),
                "AIC_covariates_only": float(
                    base_fit.aic
                ),
                "AIC_EPOCH_adjusted": float(
                    epoch_fit.aic
                ),
                "AIC_all9_adjusted": float(
                    all9_fit.aic
                ),
                "AIC_all9_plus_EPOCH_adjusted": float(
                    all10_fit.aic
                ),
                "lrt_chisq_all9_vs_covariates": stat_all9_vs_base,
                "lrt_p_all9_vs_covariates": p_all9_vs_base,
                "lrt_chisq_EPOCH_beyond_all9": stat_epoch_beyond_all9,
                "lrt_p_EPOCH_beyond_all9": p_epoch_beyond_all9,
                "lrt_chisq_all9_beyond_EPOCH": stat_all9_beyond_epoch,
                "lrt_p_all9_beyond_EPOCH": p_all9_beyond_epoch,
                "EPOCH_conditional_OR_given_all9": epoch_conditional[
                    "or"
                ],
                "EPOCH_conditional_CI_low": epoch_conditional[
                    "ci_low"
                ],
                "EPOCH_conditional_CI_high": epoch_conditional[
                    "ci_high"
                ],
                "EPOCH_conditional_p_given_all9": epoch_conditional[
                    "p"
                ],
            }
        ]
    )


def get_cv_model_specs() -> Dict[str, List[str]]:
    specs: Dict[str, List[str]] = {}

    specs["Covariates_only"] = BASE_COVARS.copy()

    specs[
        "EPOCH_adjusted"
    ] = BASE_COVARS + [
        "EPOCH_z"
    ]

    for s in SUBTYPES:
        specs[
            f"{s}_adjusted"
        ] = BASE_COVARS + [
            f"{s}_z"
        ]

    specs[
        "All9_adjusted"
    ] = BASE_COVARS + [
        f"{s}_z"
        for s in SUBTYPES
    ]

    specs[
        "All9_plus_EPOCH_adjusted"
    ] = BASE_COVARS + [
        f"{s}_z"
        for s in SUBTYPES
    ] + [
        "EPOCH_z"
    ]

    specs[
        "EPOCH_only"
    ] = [
        "EPOCH_z"
    ]

    for s in SUBTYPES:
        specs[
            f"{s}_only"
        ] = [
            f"{s}_z"
        ]

    specs[
        "All9_only"
    ] = [
        f"{s}_z"
        for s in SUBTYPES
    ]

    specs[
        "All9_plus_EPOCH_only"
    ] = [
        f"{s}_z"
        for s in SUBTYPES
    ] + [
        "EPOCH_z"
    ]

    return specs


def fit_predict_fold(
    train_df: pd.DataFrame,
    test_df: pd.DataFrame,
    predictors: List[str],
    label: str,
) -> np.ndarray:
    y_train = (
        pd.to_numeric(
            train_df[
                "case"
            ],
            errors="coerce",
        )
        .astype(int)
    )

    x_train = design_matrix(
        train_df,
        predictors,
    )

    x_test = design_matrix(
        test_df,
        predictors,
    )

    x_test = x_test.reindex(
        columns=x_train.columns,
        fill_value=0.0,
    )

    try:
        model = sm.GLM(
            y_train,
            x_train,
            family=sm.families.Binomial(),
        )

        fit = model.fit(
            maxiter=200,
            disp=0,
        )

        pred = np.asarray(
            fit.predict(
                x_test
            ),
            dtype=float,
        )

    except Exception:
        # Prediction-only numerical fallback for rare separation/convergence cases.
        model = sm.GLM(
            y_train,
            x_train,
            family=sm.families.Binomial(),
        )

        fit = model.fit_regularized(
            alpha=1e-6,
            L1_wt=0.0,
            maxiter=1000,
        )

        pred = np.asarray(
            fit.predict(
                x_test
            ),
            dtype=float,
        )

    pred = np.clip(
        pred,
        1e-8,
        1.0 - 1e-8,
    )

    return pred


def repeated_cv_predictions(
    d: pd.DataFrame,
    args,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    y = (
        pd.to_numeric(
            d[
                "case"
            ],
            errors="coerce",
        )
        .astype(int)
        .to_numpy()
    )

    specs = get_cv_model_specs()

    n = len(d)

    pred_sum = {
        name: np.zeros(
            n,
            dtype=float,
        )
        for name in specs
    }

    pred_count = {
        name: np.zeros(
            n,
            dtype=np.int32,
        )
        for name in specs
    }

    cv = RepeatedStratifiedKFold(
        n_splits=args.cv_folds,
        n_repeats=args.cv_repeats,
        random_state=args.seed,
    )

    total_splits = (
        args.cv_folds
        * args.cv_repeats
    )

    print(
        f"Starting repeated stratified CV: "
        f"{args.cv_folds} folds x {args.cv_repeats} repeats "
        f"= {total_splits} held-out splits.",
        flush=True,
    )

    for split_i, (
        train_idx,
        test_idx,
    ) in enumerate(
        cv.split(
            np.zeros(
                n
            ),
            y,
        ),
        start=1,
    ):
        train_df = d.iloc[
            train_idx
        ].copy()

        test_df = d.iloc[
            test_idx
        ].copy()

        for model_name, predictors in specs.items():
            pred = fit_predict_fold(
                train_df,
                test_df,
                predictors,
                model_name,
            )

            pred_sum[
                model_name
            ][
                test_idx
            ] += pred

            pred_count[
                model_name
            ][
                test_idx
            ] += 1

        if (
            split_i % args.cv_folds == 0
            or split_i == total_splits
        ):
            print(
                f"  completed CV split "
                f"{split_i}/{total_splits}",
                flush=True,
            )

    oof = pd.DataFrame(
        {
            "participant_id": d[
                "participant_id"
            ].astype(
                "Int64"
            ).to_numpy(),
            "case": y,
        }
    )

    metric_rows = []

    prevalence = float(
        y.mean()
    )

    for model_name in specs:
        if np.any(
            pred_count[
                model_name
            ] == 0
        ):
            raise RuntimeError(
                f"CV internal error: some participants have no prediction for {model_name}."
            )

        pred = (
            pred_sum[
                model_name
            ]
            / pred_count[
                model_name
            ]
        )

        oof[
            f"pred_{model_name}"
        ] = pred

        auc = float(
            roc_auc_score(
                y,
                pred,
            )
        )

        auprc = float(
            average_precision_score(
                y,
                pred,
            )
        )

        brier = float(
            brier_score_loss(
                y,
                pred,
            )
        )

        ll = float(
            log_loss(
                y,
                pred,
                labels=[
                    0,
                    1,
                ],
            )
        )

        metric_rows.append(
            {
                "model": model_name,
                "N": int(
                    n
                ),
                "N_case": int(
                    y.sum()
                ),
                "N_control": int(
                    (
                        y == 0
                    ).sum()
                ),
                "case_prevalence": prevalence,
                "ROC_AUC_repeated_OOF": auc,
                "AUPRC_repeated_OOF": auprc,
                "Brier_repeated_OOF": brier,
                "log_loss_repeated_OOF": ll,
                "cv_folds": int(
                    args.cv_folds
                ),
                "cv_repeats": int(
                    args.cv_repeats
                ),
            }
        )

    return (
        oof,
        pd.DataFrame(
            metric_rows
        ),
    )


def bootstrap_difference(
    y: np.ndarray,
    pred_a: np.ndarray,
    pred_b: np.ndarray,
    B: int,
    rng,
) -> Dict[str, float]:
    observed_auc = float(
        roc_auc_score(
            y,
            pred_a,
        )
        - roc_auc_score(
            y,
            pred_b,
        )
    )

    observed_auprc = float(
        average_precision_score(
            y,
            pred_a,
        )
        - average_precision_score(
            y,
            pred_b,
        )
    )

    auc_diffs = []
    auprc_diffs = []

    n = len(y)

    for _ in range(B):
        idx = rng.integers(
            0,
            n,
            size=n,
        )

        yb = y[
            idx
        ]

        if np.unique(
            yb
        ).size < 2:
            continue

        pa = pred_a[
            idx
        ]

        pb = pred_b[
            idx
        ]

        auc_diffs.append(
            roc_auc_score(
                yb,
                pa,
            )
            - roc_auc_score(
                yb,
                pb,
            )
        )

        auprc_diffs.append(
            average_precision_score(
                yb,
                pa,
            )
            - average_precision_score(
                yb,
                pb,
            )
        )

    def summarize(
        vals,
        observed,
    ):
        v = np.asarray(
            vals,
            dtype=float,
        )

        v = v[
            np.isfinite(
                v
            )
        ]

        if len(v) < 20:
            return (
                observed,
                np.nan,
                np.nan,
                np.nan,
                int(
                    len(v)
                ),
            )

        lo, hi = np.quantile(
            v,
            [
                0.025,
                0.975,
            ],
        )

        p_left = (
            np.sum(
                v <= 0
            )
            + 1.0
        ) / (
            len(v)
            + 1.0
        )

        p_right = (
            np.sum(
                v >= 0
            )
            + 1.0
        ) / (
            len(v)
            + 1.0
        )

        p = min(
            1.0,
            2.0
            * min(
                p_left,
                p_right,
            ),
        )

        return (
            float(
                observed
            ),
            float(
                lo
            ),
            float(
                hi
            ),
            float(
                p
            ),
            int(
                len(v)
            ),
        )

    auc_summary = summarize(
        auc_diffs,
        observed_auc,
    )

    auprc_summary = summarize(
        auprc_diffs,
        observed_auprc,
    )

    return {
        "delta_ROC_AUC": auc_summary[
            0
        ],
        "ROC_AUC_CI_low": auc_summary[
            1
        ],
        "ROC_AUC_CI_high": auc_summary[
            2
        ],
        "ROC_AUC_bootstrap_p": auc_summary[
            3
        ],
        "delta_AUPRC": auprc_summary[
            0
        ],
        "AUPRC_CI_low": auprc_summary[
            1
        ],
        "AUPRC_CI_high": auprc_summary[
            2
        ],
        "AUPRC_bootstrap_p": auprc_summary[
            3
        ],
        "successful_bootstrap_replicates": min(
            auc_summary[
                4
            ],
            auprc_summary[
                4
            ],
        ),
    }


def holm_adjust(
    pvalues: pd.Series,
) -> np.ndarray:
    p = pd.to_numeric(
        pvalues,
        errors="coerce",
    ).to_numpy(
        dtype=float
    )

    out = np.full(
        len(p),
        np.nan,
        dtype=float,
    )

    valid_idx = np.where(
        np.isfinite(
            p
        )
    )[0]

    if len(
        valid_idx
    ) == 0:
        return out

    p_valid = p[
        valid_idx
    ]

    order = np.argsort(
        p_valid
    )

    m = len(
        p_valid
    )

    adjusted_ordered = np.empty(
        m,
        dtype=float,
    )

    running_max = 0.0

    for rank, ord_idx in enumerate(
        order
    ):
        raw = p_valid[
            ord_idx
        ]

        multiplier = (
            m
            - rank
        )

        adj = min(
            1.0,
            raw
            * multiplier,
        )

        running_max = max(
            running_max,
            adj,
        )

        adjusted_ordered[
            ord_idx
        ] = min(
            1.0,
            running_max,
        )

    out[
        valid_idx
    ] = adjusted_ordered

    return out


def bootstrap_comparisons(
    oof: pd.DataFrame,
    args,
) -> pd.DataFrame:
    y = (
        pd.to_numeric(
            oof[
                "case"
            ],
            errors="coerce",
        )
        .astype(int)
        .to_numpy()
    )

    comparisons = []

    # Direct EPOCH versus each disease subtype, adjusted and biomarker-only.
    for context in [
        "adjusted",
        "only",
    ]:
        epoch_model = (
            "EPOCH_adjusted"
            if context == "adjusted"
            else "EPOCH_only"
        )

        for s in SUBTYPES:
            subtype_model = (
                f"{s}_adjusted"
                if context == "adjusted"
                else f"{s}_only"
            )

            comparisons.append(
                (
                    f"EPOCH_vs_{s}",
                    context,
                    epoch_model,
                    subtype_model,
                )
            )

    # EPOCH versus the combined nine-subtype panel.
    comparisons.extend(
        [
            (
                "EPOCH_vs_all9",
                "adjusted",
                "EPOCH_adjusted",
                "All9_adjusted",
            ),
            (
                "EPOCH_vs_all9",
                "biomarker_only",
                "EPOCH_only",
                "All9_only",
            ),
            (
                "EPOCH_increment_beyond_all9",
                "adjusted",
                "All9_plus_EPOCH_adjusted",
                "All9_adjusted",
            ),
            (
                "EPOCH_increment_beyond_all9",
                "biomarker_only",
                "All9_plus_EPOCH_only",
                "All9_only",
            ),
        ]
    )

    rng = np.random.default_rng(
        args.seed
        + 1009
    )

    rows = []

    print(
        f"Starting paired bootstrap discrimination comparisons: "
        f"B={args.bootstrap}",
        flush=True,
    )

    for i, (
        comparison,
        context,
        model_a,
        model_b,
    ) in enumerate(
        comparisons,
        start=1,
    ):
        pa = pd.to_numeric(
            oof[
                f"pred_{model_a}"
            ],
            errors="coerce",
        ).to_numpy(
            dtype=float
        )

        pb = pd.to_numeric(
            oof[
                f"pred_{model_b}"
            ],
            errors="coerce",
        ).to_numpy(
            dtype=float
        )

        summary = bootstrap_difference(
            y,
            pa,
            pb,
            args.bootstrap,
            rng,
        )

        rows.append(
            {
                "comparison": comparison,
                "context": context,
                "model_A": model_a,
                "model_B": model_b,
                "positive_difference_favors": model_a,
                **summary,
            }
        )

        if (
            i % 5 == 0
            or i == len(
                comparisons
            )
        ):
            print(
                f"  completed bootstrap comparison "
                f"{i}/{len(comparisons)}",
                flush=True,
            )

    res = pd.DataFrame(
        rows
    )

    # Holm adjustment across the nine direct EPOCH-vs-subtype comparisons
    # separately for each modeling context.
    res[
        "ROC_AUC_bootstrap_p_Holm_9subtypes"
    ] = np.nan

    res[
        "AUPRC_bootstrap_p_Holm_9subtypes"
    ] = np.nan

    for context in [
        "adjusted",
        "only",
    ]:
        mask = (
            res[
                "context"
            ].eq(
                context
            )
            & res[
                "comparison"
            ].str.match(
                r"^EPOCH_vs_(AD1|AD2|ASD1|ASD2|ASD3|LLD1|LLD2|SCZ1|SCZ2)$"
            )
        )

        res.loc[
            mask,
            "ROC_AUC_bootstrap_p_Holm_9subtypes",
        ] = holm_adjust(
            res.loc[
                mask,
                "ROC_AUC_bootstrap_p",
            ]
        )

        res.loc[
            mask,
            "AUPRC_bootstrap_p_Holm_9subtypes",
        ] = holm_adjust(
            res.loc[
                mask,
                "AUPRC_bootstrap_p",
            ]
        )

    return res


def print_table(
    title: str,
    df: pd.DataFrame,
):
    print(
        "\n"
        + "=" * 120
    )
    print(
        title
    )
    print(
        "=" * 120
    )

    with pd.option_context(
        "display.max_columns",
        None,
        "display.width",
        260,
        "display.max_rows",
        300,
    ):
        print(
            df.to_string(
                index=False
            )
            if len(
                df
            )
            else "<empty>"
        )


def main():
    args = parse_args()

    os.makedirs(
        args.output_dir,
        exist_ok=True,
    )

    qc_rows = []

    epoch, q1 = read_epoch(
        args
    )

    qc_rows.extend(
        [
            {
                "metric": k,
                "value": v,
            }
            for k, v in q1.items()
        ]
    )

    epoch_ids = set(
        epoch[
            "participant_id"
        ]
        .dropna()
        .astype(int)
    )

    subtypes, q2 = read_subtypes(
        args,
        epoch_ids,
    )

    qc_rows.extend(
        [
            {
                "metric": k,
                "value": v,
            }
            for k, v in q2.items()
        ]
    )

    base = epoch.merge(
        subtypes,
        on="participant_id",
        how="inner",
    )

    qc_rows.append(
        {
            "metric": "epoch_subtype_candidate_overlap_N",
            "value": int(
                len(base)
            ),
        }
    )

    candidate_ids = set(
        base[
            "participant_id"
        ]
        .dropna()
        .astype(int)
    )

    cov = read_covariates(
        args
    )

    diag, code_cols = read_icd_codes(
        args,
        candidate_ids,
    )

    dates = read_icd_dates(
        args,
        candidate_ids,
    )

    with_fg_dates, q3 = derive_fg_dates(
        base,
        diag,
        code_cols,
        dates,
    )

    qc_rows.extend(
        [
            {
                "metric": k,
                "value": v,
            }
            for k, v in q3.items()
        ]
    )

    endpoint, q4 = define_binary_endpoint(
        with_fg_dates,
        args,
    )

    qc_rows.extend(
        [
            {
                "metric": k,
                "value": v,
            }
            for k, v in q4.items()
        ]
    )

    analysis, std_tbl, flow = build_common_sample(
        endpoint,
        subtypes,
        cov,
    )

    n_case = int(
        analysis[
            "case"
        ].sum()
    )

    n_control = int(
        (
            analysis[
                "case"
            ]
            == 0
        ).sum()
    )

    if n_case < args.min_cases:
        raise RuntimeError(
            f"Only {n_case} cases remain in the strict common sample; "
            f"min_cases={args.min_cases}."
        )

    if n_control < args.min_controls:
        raise RuntimeError(
            f"Only {n_control} controls remain in the strict common sample; "
            f"min_controls={args.min_controls}."
        )

    if n_case < args.cv_folds:
        raise RuntimeError(
            f"N_case={n_case} is smaller than cv_folds={args.cv_folds}."
        )

    unadjusted_tbl, unadjusted_models = single_marker_inference(
        analysis,
        adjusted=False,
    )

    adjusted_tbl, adjusted_models = single_marker_inference(
        analysis,
        adjusted=True,
    )

    pairwise_tbl = adjusted_pairwise_inference(
        analysis,
        adjusted_models,
    )

    combined_tbl = combined_adjusted_inference(
        analysis,
        adjusted_models,
    )

    oof_tbl, cv_metrics_tbl = repeated_cv_predictions(
        analysis,
        args,
    )

    bootstrap_tbl = bootstrap_comparisons(
        oof_tbl,
        args,
    )

    # --------------------------------------------------------------------------
    # Save outputs
    # --------------------------------------------------------------------------

    pd.DataFrame(
        qc_rows
    ).to_csv(
        os.path.join(
            args.output_dir,
            "cohort_qc.tsv",
        ),
        sep="\t",
        index=False,
    )

    flow.to_csv(
        os.path.join(
            args.output_dir,
            "sample_flow.tsv",
        ),
        sep="\t",
        index=False,
    )

    std_tbl.to_csv(
        os.path.join(
            args.output_dir,
            "marker_standardization.tsv",
        ),
        sep="\t",
        index=False,
    )

    save_cols = [
        "participant_id",
        "imaging_date",
        "censor_date",
        "earliest_FG_date",
        "earliest_FG_code",
        "n_FG_codes_recorded",
        "has_FG_code_missing_paired_date",
        "prevalent_FG_at_imaging",
        "case",
        "followup_years_for_binary_endpoint",
    ] + BASE_COVARS + [
        EPOCH_RAW,
        "EPOCH_z",
    ] + SUBTYPES + [
        f"{s}_z"
        for s in SUBTYPES
    ]

    if (
        "horizon_date"
        in analysis.columns
    ):
        save_cols.append(
            "horizon_date"
        )

    if (
        "split"
        in analysis.columns
    ):
        save_cols.append(
            "split"
        )

    save_cols = [
        c for c in save_cols
        if c in analysis.columns
    ]

    analysis[
        save_cols
    ].to_csv(
        os.path.join(
            args.output_dir,
            "analysis_common_sample.tsv",
        ),
        sep="\t",
        index=False,
    )

    unadjusted_tbl.to_csv(
        os.path.join(
            args.output_dir,
            "single_marker_unadjusted.tsv",
        ),
        sep="\t",
        index=False,
    )

    adjusted_tbl.to_csv(
        os.path.join(
            args.output_dir,
            "single_marker_adjusted.tsv",
        ),
        sep="\t",
        index=False,
    )

    pairwise_tbl.to_csv(
        os.path.join(
            args.output_dir,
            "epoch_vs_subtype_adjusted_joint.tsv",
        ),
        sep="\t",
        index=False,
    )

    combined_tbl.to_csv(
        os.path.join(
            args.output_dir,
            "combined_adjusted_models.tsv",
        ),
        sep="\t",
        index=False,
    )

    cv_metrics_tbl.to_csv(
        os.path.join(
            args.output_dir,
            "cv_model_metrics.tsv",
        ),
        sep="\t",
        index=False,
    )

    oof_tbl.to_csv(
        os.path.join(
            args.output_dir,
            "cv_out_of_fold_predictions.tsv",
        ),
        sep="\t",
        index=False,
    )

    bootstrap_tbl.to_csv(
        os.path.join(
            args.output_dir,
            "bootstrap_discrimination_comparisons.tsv",
        ),
        sep="\t",
        index=False,
    )

    with open(
        os.path.join(
            args.output_dir,
            "run_metadata.txt",
        ),
        "w",
    ) as f:
        f.write(
            f"split={args.split}\n"
        )
        f.write(
            f"endpoint_mode={args.endpoint_mode}\n"
        )
        f.write(
            f"horizon_years={args.horizon_years}\n"
        )
        f.write(
            "endpoint_codes=ICD-10 F* or G* inpatient diagnoses\n"
        )
        f.write(
            "prevalent_exclusion=any valid F/G diagnosis on_or_before MRI\n"
        )
        f.write(
            f"strict_unresolved_FG_exclusion={not args.allow_unresolved_fg_dates}\n"
        )
        f.write(
            "common_complete_case=True\n"
        )
        f.write(
            f"base_covariates={','.join(BASE_COVARS)}\n"
        )
        f.write(
            "biomarker_standardization=within exact common logistic sample\n"
        )
        f.write(
            f"cv_folds={args.cv_folds}\n"
        )
        f.write(
            f"cv_repeats={args.cv_repeats}\n"
        )
        f.write(
            f"bootstrap={args.bootstrap}\n"
        )
        f.write(
            f"seed={args.seed}\n"
        )
        f.write(
            f"analysis_N={len(analysis)}\n"
        )
        f.write(
            f"analysis_cases={n_case}\n"
        )
        f.write(
            f"analysis_controls={n_control}\n"
        )

    # --------------------------------------------------------------------------
    # Print main results
    # --------------------------------------------------------------------------

    print_table(
        "SAMPLE FLOW",
        flow,
    )

    print_table(
        "BIOMARKER-ONLY LOGISTIC MODELS",
        unadjusted_tbl,
    )

    print_table(
        "COVARIATE-ADJUSTED SINGLE-MARKER LOGISTIC MODELS",
        adjusted_tbl,
    )

    print_table(
        "ADJUSTED JOINT LOGISTIC MODELS: EPOCH + EACH SUBTYPE",
        pairwise_tbl,
    )

    print_table(
        "ADJUSTED COMBINED MODELS: ALL 9 SUBTYPES VS EPOCH",
        combined_tbl,
    )

    print_table(
        "REPEATED OUT-OF-FOLD PREDICTIVE PERFORMANCE",
        cv_metrics_tbl,
    )

    print_table(
        "PAIRED BOOTSTRAP DISCRIMINATION COMPARISONS",
        bootstrap_tbl,
    )

    print(
        "\nAnalysis complete.",
        flush=True,
    )

    print(
        f"Output directory: {args.output_dir}",
        flush=True,
    )


if __name__ == "__main__":
    main()
