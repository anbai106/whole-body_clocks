#!/usr/bin/env python3
"""
Post-hoc proportional-hazards (PH) diagnostics for the 23 existing mortality
EPOCH clocks.

IMPORTANT
---------
This script DOES NOT retrain or retune any mortality EPOCH model.

It reads the existing participant-level:
    *_mortality_clock_predictions.tsv

files, restricts the primary analysis to the held-out TEST split, and evaluates
whether the EPOCH predictor satisfies the proportional-hazards assumption.

Primary predictor
-----------------
By default, the script uses the saved mortality-clock acceleration z-score:
    *_mortality_clock_acceleration_z

This is appropriate when the reported downstream Cox hazard ratios use EPOCH
acceleration. If the manuscript instead reports the original Cox linear
predictor/risk score, run with:
    --predictor-mode risk_score

If acceleration is requested but unavailable, the script falls back to the
saved mortality risk score, standardized within the held-out test set, and
records that fallback explicitly.

Primary Cox model
-----------------
The default post-hoc model is:
    Surv(time_years, event) ~ EPOCH + age + sex

Age and sex are included as covariates in the primary analysis. Importantly,
the clock names "Reproductive female proteomics" and "Reproductive male
proteomics" refer to organ/proteomic feature definitions; they do NOT imply
female-only or male-only analysis samples. Participants are never filtered by
clock name. Therefore, both males and females are retained for these clocks,
and sex is included in the Cox model exactly as for the other clocks.

By default, --covariate-mode age_sex requires both age and a nonconstant
male/female sex variable to be available in the held-out test set. This strict
behavior prevents accidental omission of sex from any clock. Use
--covariate-mode none only for an intentionally unadjusted sensitivity analysis.

PH diagnostics
--------------
For each clock, the script reports:

1. Standard Cox association for the EPOCH predictor
   - coefficient
   - HR and 95% CI
   - Wald P value

2. Schoenfeld-residual PH diagnostic
   - lifelines proportional_hazard_test
   - primary time transform = Kaplan-Meier ("km"), analogous in spirit to
     cox.zph(..., transform="km")
   - EPOCH-specific PH chi-square statistic and P value
   - BH-FDR across the 23 EPOCH-specific PH tests
   - covariate-specific PH tests saved in a separate long table

3. Descriptive scaled-Schoenfeld residual trend
   - Spearman rho between the EPOCH scaled Schoenfeld residual and log event
     time (descriptive; NOT used as the formal PH test)
   - residual diagnostic plot for each clock

4. Piecewise time-varying EPOCH sensitivity analysis
   - by default splits follow-up at 5 years
   - compares:
       reduced: EPOCH + age + sex
       full:    EPOCH + EPOCH x I(t > 5y) + age + sex
     using CoxTimeVaryingFitter on start-stop data
   - reports HR during 0-5 years
   - reports HR after 5 years
   - reports interaction Wald P
   - reports a 1-df likelihood-ratio P when both models fit without penalization

The 5-year split is a prespecified sensitivity analysis, not a replacement for
the Schoenfeld test. Change with --time-split-years if another clinically
meaningful split is preferred.

Outputs
-------
Master directory:
    <base-dir>/mortality_EPOCH_PH_diagnostics

Master files:
    mortality_EPOCH_23_PH_diagnostics.tsv
    mortality_EPOCH_23_PH_covariate_tests.tsv
    mortality_EPOCH_23_PH_manuscript_table.tsv
    mortality_EPOCH_23_PH_run_manifest.tsv

Per-clock directory:
    <clock folder>/ph_diagnostics/

Per-clock files include:
    *_PH_summary.tsv
    *_PH_covariate_tests.tsv
    *_scaled_schoenfeld_epoch.tsv
    *_scaled_schoenfeld_epoch.png
    *_scaled_schoenfeld_epoch.pdf
    *_time_varying_5y.tsv

Dependencies
------------
    pandas
    numpy
    scipy
    lifelines
    matplotlib

Example
-------
python compute_23_mortality_EPOCH_PH_diagnostics.py \
    --base-dir /cbica/home/wenju/Reproducibile_paper/WholeBodyClock \
    --test-split test \
    --predictor-mode acceleration \
    --covariate-mode age_sex \
    --time-transform km \
    --time-split-years 5
"""

from __future__ import print_function

import argparse
import math
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from scipy.stats import chi2, spearmanr
except ImportError as exc:
    raise ImportError("This script requires scipy.") from exc

try:
    from lifelines import CoxPHFitter, CoxTimeVaryingFitter
    from lifelines.statistics import proportional_hazard_test
except ImportError as exc:
    raise ImportError(
        "This script requires lifelines. In the survival_clock environment, "
        "install with: conda install -c conda-forge lifelines"
    ) from exc

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError as exc:
    raise ImportError("This script requires matplotlib.") from exc


# =============================================================================
# 1. Exact 23-clock manifest copied from the existing calibration pipeline
# =============================================================================

CLOCKS = [
    {"folder": "adipose_mri_mortality_clock", "clock": "Adipose MRI", "modality": "MRI"},
    {"folder": "brain_mri_mortality_clock", "clock": "Brain MRI", "modality": "MRI"},
    {"folder": "Brain_proteomics_mortality_clock", "clock": "Brain proteomics", "modality": "Proteomics"},
    {"folder": "Digestive_metabolomics_mortality_clock", "clock": "Digestive metabolomics", "modality": "Metabolomics"},
    {"folder": "Endocrine_metabolomics_mortality_clock", "clock": "Endocrine metabolomics", "modality": "Metabolomics"},
    {"folder": "Endocrine_proteomics_mortality_clock", "clock": "Endocrine proteomics", "modality": "Proteomics"},
    {"folder": "Eye_proteomics_mortality_clock", "clock": "Eye proteomics", "modality": "Proteomics"},
    {"folder": "heart_mri_mortality_clock", "clock": "Heart MRI", "modality": "MRI"},
    {"folder": "Heart_proteomics_mortality_clock", "clock": "Heart proteomics", "modality": "Proteomics"},
    {"folder": "Hepatic_metabolomics_mortality_clock", "clock": "Hepatic metabolomics", "modality": "Metabolomics"},
    {"folder": "Hepatic_proteomics_mortality_clock", "clock": "Hepatic proteomics", "modality": "Proteomics"},
    {"folder": "Immune_metabolomics_mortality_clock", "clock": "Immune metabolomics", "modality": "Metabolomics"},
    {"folder": "Immune_proteomics_mortality_clock", "clock": "Immune proteomics", "modality": "Proteomics"},
    {"folder": "kidney_mri_mortality_clock", "clock": "Kidney MRI", "modality": "MRI"},
    {"folder": "liver_mri_mortality_clock", "clock": "Liver MRI", "modality": "MRI"},
    {"folder": "Metabolic_metabolomics_mortality_clock", "clock": "Metabolic metabolomics", "modality": "Metabolomics"},
    {"folder": "pancreas_mri_mortality_clock", "clock": "Pancreas MRI", "modality": "MRI"},
    {"folder": "Pulmonary_proteomics_mortality_clock", "clock": "Pulmonary proteomics", "modality": "Proteomics"},
    {"folder": "Renal_proteomics_mortality_clock", "clock": "Renal proteomics", "modality": "Proteomics"},
    {"folder": "Reproductive_female_proteomics_mortality_clock", "clock": "Reproductive female proteomics", "modality": "Proteomics"},
    {"folder": "Reproductive_male_proteomics_mortality_clock", "clock": "Reproductive male proteomics", "modality": "Proteomics"},
    {"folder": "Skin_proteomics_mortality_clock", "clock": "Skin proteomics", "modality": "Proteomics"},
    {"folder": "spleen_mri_mortality_clock", "clock": "Spleen MRI", "modality": "MRI"},
]


# =============================================================================
# 2. CLI
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Compute post-hoc proportional-hazards diagnostics for the 23 "
            "existing mortality EPOCH clocks from saved prediction TSV files."
        )
    )

    p.add_argument(
        "--base-dir",
        default="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock",
    )

    p.add_argument(
        "--output-dir",
        default=None,
        help=(
            "Master output directory. Default: "
            "<base-dir>/mortality_EPOCH_PH_diagnostics"
        ),
    )

    p.add_argument(
        "--test-split",
        default="test",
        help="Held-out split used for the primary PH diagnostics.",
    )

    p.add_argument(
        "--predictor-mode",
        choices=["acceleration", "risk_score"],
        default="acceleration",
        help=(
            "Primary EPOCH predictor to test. 'acceleration' uses the saved "
            "mortality-clock acceleration z-score when available; 'risk_score' "
            "uses the original mortality risk score standardized within test."
        ),
    )

    p.add_argument(
        "--covariate-mode",
        choices=["age_sex", "none"],
        default="age_sex",
        help=(
            "Primary PH Cox model adjustment. 'age_sex' REQUIRES both age "
            "and a nonconstant male/female sex variable and includes both in "
            "every clock. Clock names never trigger sex-based filtering. "
            "'none' tests EPOCH alone."
        ),
    )

    p.add_argument(
        "--time-transform",
        choices=["km", "rank", "log", "identity"],
        default="km",
        help=(
            "Time transform passed to lifelines proportional_hazard_test. "
            "Default 'km' is closest to the usual cox.zph KM transform."
        ),
    )

    p.add_argument(
        "--time-split-years",
        type=float,
        default=5.0,
        help=(
            "Prespecified split for the piecewise time-varying sensitivity "
            "analysis. Default: 5 years."
        ),
    )

    p.add_argument(
        "--min-events-per-time-band",
        type=int,
        default=20,
        help=(
            "Minimum deaths required both on/before and after the time split "
            "to fit the piecewise time-varying sensitivity model."
        ),
    )

    p.add_argument(
        "--skip-time-varying",
        action="store_true",
        help="Run Schoenfeld diagnostics only and skip the piecewise model.",
    )

    p.add_argument(
        "--fail-on-missing-clock",
        action="store_true",
        help="Stop instead of recording a failed manifest row.",
    )

    return p.parse_args()


# =============================================================================
# 3. General helpers
# =============================================================================

def info(msg):
    print(msg, flush=True)


def warn(msg):
    print("WARNING: {}".format(msg), file=sys.stderr, flush=True)


def sanitize_prefix(x):
    x = re.sub(r"_predictions\.tsv$", "", str(x), flags=re.IGNORECASE)
    x = re.sub(r"[^A-Za-z0-9._-]+", "_", x)
    return x.strip("_")


def parse_event(series):
    if pd.api.types.is_bool_dtype(series):
        return series.astype(bool)

    numeric = pd.to_numeric(series, errors="coerce")

    if numeric.notna().mean() >= 0.90:
        return numeric.fillna(0).astype(int).astype(bool)

    s = series.astype(str).str.strip().str.lower()

    true_vals = {
        "1", "true", "t", "yes", "y", "event", "dead", "death"
    }

    false_vals = {
        "0", "false", "f", "no", "n", "censored", "alive"
    }

    recognized = s.isin(true_vals | false_vals)

    if recognized.mean() < 0.90:
        bad = sorted(
            s[~recognized].dropna().unique()
        )[:10]

        raise ValueError(
            "Could not reliably parse event values. Examples: {}".format(bad)
        )

    return s.isin(true_vals)


def normalize_sex(series):
    s = series.astype(str).str.strip().str.lower()

    out = pd.Series(
        np.nan,
        index=series.index,
        dtype="float64",
    )

    female = s.isin(
        {
            "0", "0.0", "f", "female", "woman", "women"
        }
    )

    male = s.isin(
        {
            "1", "1.0", "m", "male", "man", "men"
        }
    )

    out.loc[female] = 0.0
    out.loc[male] = 1.0

    return out


def resolve_prediction_file(clock_dir):
    hits = sorted(
        clock_dir.glob("*_mortality_clock_predictions.tsv")
    )

    hits = [
        x
        for x in hits
        if ".before_" not in x.name
        and x.name.endswith("_mortality_clock_predictions.tsv")
    ]

    if len(hits) == 1:
        return hits[0]

    if len(hits) == 0:
        raise FileNotFoundError(
            "No primary *_mortality_clock_predictions.tsv found in {}".format(
                clock_dir
            )
        )

    raise RuntimeError(
        "Multiple primary prediction files found in {}: {}".format(
            clock_dir,
            [x.name for x in hits],
        )
    )


def detect_primary_risk_score(columns):
    candidates = [
        str(c)
        for c in columns
        if str(c).endswith("_mortality_risk_score")
        and not str(c).startswith("risk_score_")
    ]

    candidates = [
        c
        for c in candidates
        if "M0_" not in c
        and "M1_" not in c
        and "M2_" not in c
        and "M3_" not in c
    ]

    if len(candidates) == 1:
        return candidates[0]

    if len(candidates) == 0:
        return None

    candidates = sorted(
        candidates,
        key=lambda x: (len(x), x),
    )

    warn(
        "Multiple primary risk-score candidates found; using {} from {}".format(
            candidates[0],
            candidates,
        )
    )

    return candidates[0]


def detect_acceleration_z(columns):
    cols = [str(c) for c in columns]

    preferred = [
        c
        for c in cols
        if re.search(
            r"_mortality_(?:clock|epoch)_acceleration_z$",
            c,
            flags=re.IGNORECASE,
        )
    ]

    if len(preferred) == 1:
        return preferred[0]

    if len(preferred) > 1:
        preferred = sorted(
            preferred,
            key=lambda x: (
                "clock" not in x.lower(),
                len(x),
                x,
            ),
        )

        warn(
            "Multiple mortality acceleration-z candidates found; using {} "
            "from {}".format(
                preferred[0],
                preferred,
            )
        )

        return preferred[0]

    broader = [
        c
        for c in cols
        if "mortality" in c.lower()
        and "acceleration" in c.lower()
        and c.lower().endswith("_z")
    ]

    if len(broader) == 1:
        return broader[0]

    return None


def detect_age_column(columns):
    cols = list(columns)

    preferred = [
        "age_at_baseline",
        "age_at_imaging",
        "age_when_attended_assessment_centre_f21003_0_0",
        "Age",
        "age",
    ]

    for c in preferred:
        if c in cols:
            return c

    hits = [
        c
        for c in cols
        if re.search(
            r"(^|_)age($|_at_|_when_)",
            str(c),
            flags=re.IGNORECASE,
        )
    ]

    if len(hits) == 1:
        return hits[0]

    return None


def detect_sex_column(columns):
    cols = list(columns)

    preferred = [
        "sex",
        "Sex",
        "SEX",
        "sex_f31_0_0",
    ]

    for c in preferred:
        if c in cols:
            return c

    hits = [
        c
        for c in cols
        if re.search(
            r"(^|_)(sex|gender)($|_)",
            str(c),
            flags=re.IGNORECASE,
        )
    ]

    if len(hits) == 1:
        return hits[0]

    return None


def standardize_numeric(series):
    x = pd.to_numeric(
        series,
        errors="coerce",
    )

    mean = float(
        x.mean()
    )

    sd = float(
        x.std(ddof=1)
    )

    if not np.isfinite(sd) or sd <= 0:
        return pd.Series(
            np.nan,
            index=series.index,
            dtype="float64",
        ), mean, sd

    return (
        (x - mean) / sd
    ), mean, sd


def bh_fdr(pvalues):
    p = np.asarray(
        pvalues,
        dtype=float,
    )

    out = np.full(
        p.shape,
        np.nan,
        dtype=float,
    )

    ok = np.isfinite(p)

    if not np.any(ok):
        return out

    pv = p[ok]
    m = len(pv)

    order = np.argsort(
        pv
    )

    ranked = pv[order]

    adjusted = ranked * m / np.arange(
        1,
        m + 1,
        dtype=float,
    )

    adjusted = np.minimum.accumulate(
        adjusted[::-1]
    )[::-1]

    adjusted = np.clip(
        adjusted,
        0.0,
        1.0,
    )

    unsorted = np.empty_like(
        adjusted
    )

    unsorted[order] = adjusted
    out[ok] = unsorted

    return out


def safe_float(x):
    try:
        value = float(x)
    except Exception:
        return np.nan

    if not np.isfinite(value):
        return np.nan

    return value


# =============================================================================
# 4. Prepare one held-out test-set Cox dataset
# =============================================================================

def prepare_analysis_data(df, args):
    required = {
        "participant_id",
        "split",
        "time_years",
        "event",
    }

    missing = sorted(
        required -
        set(df.columns)
    )

    if missing:
        raise ValueError(
            "Missing required prediction columns: {}".format(missing)
        )

    d = df.copy()

    d["participant_id"] = d[
        "participant_id"
    ].astype(str)

    d["split"] = d[
        "split"
    ].astype(str).str.lower()

    d["time_years"] = pd.to_numeric(
        d["time_years"],
        errors="coerce",
    )

    d["event_bool"] = parse_event(
        d["event"]
    )

    d = d.loc[
        d["split"] ==
        str(args.test_split).lower()
    ].copy()

    d = d.loc[
        np.isfinite(
            d["time_years"]
        )
        & (
            d["time_years"] >
            0
        )
    ].copy()

    if d.empty:
        raise ValueError(
            "No usable held-out test rows."
        )

    acceleration_col = detect_acceleration_z(
        d.columns
    )

    risk_score_col = detect_primary_risk_score(
        d.columns
    )

    predictor_source = None
    predictor_scale_note = None

    if args.predictor_mode == "acceleration":
        if acceleration_col is not None:
            predictor_source = acceleration_col

            d["epoch_predictor"] = pd.to_numeric(
                d[acceleration_col],
                errors="coerce",
            )

            predictor_scale_note = (
                "Existing saved EPOCH acceleration z-score; not re-standardized."
            )

        elif risk_score_col is not None:
            warn(
                "Acceleration z-score was not found; falling back to the "
                "saved mortality risk score standardized within the test set."
            )

            predictor_source = risk_score_col

            (
                d["epoch_predictor"],
                risk_mean,
                risk_sd,
            ) = standardize_numeric(
                d[risk_score_col]
            )

            predictor_scale_note = (
                "Fallback: saved mortality risk score standardized within "
                "held-out test set (mean={:.6g}, SD={:.6g})."
            ).format(
                risk_mean,
                risk_sd,
            )

        else:
            raise ValueError(
                "Neither a mortality acceleration z-score nor a primary "
                "mortality risk score could be identified."
            )

    else:
        if risk_score_col is None:
            raise ValueError(
                "--predictor-mode risk_score was requested, but no primary "
                "mortality risk score could be identified."
            )

        predictor_source = risk_score_col

        (
            d["epoch_predictor"],
            risk_mean,
            risk_sd,
        ) = standardize_numeric(
            d[risk_score_col]
        )

        predictor_scale_note = (
            "Saved mortality risk score standardized within held-out test set "
            "(mean={:.6g}, SD={:.6g})."
        ).format(
            risk_mean,
            risk_sd,
        )

    age_col = detect_age_column(
        d.columns
    )

    sex_col = detect_sex_column(
        d.columns
    )

    covariates = []

    if args.covariate_mode == "age_sex":
        # IMPORTANT:
        # Clock names such as "Reproductive female proteomics" and
        # "Reproductive male proteomics" NEVER trigger sex-based filtering.
        # These are organ/proteomic clock labels, not sex-restricted cohorts.
        # All available test participants are retained before complete-case QC.

        if age_col is None:
            raise ValueError(
                "covariate-mode age_sex requires an age column, but none "
                "could be identified in the saved prediction file."
            )

        (
            d["age_covariate_z"],
            age_mean,
            age_sd,
        ) = standardize_numeric(
            d[age_col]
        )

        if not np.isfinite(age_sd) or age_sd <= 0:
            raise ValueError(
                "covariate-mode age_sex requires a nonconstant age variable."
            )

        covariates.append(
            "age_covariate_z"
        )

        if sex_col is None:
            raise ValueError(
                "covariate-mode age_sex requires a sex column, but none "
                "could be identified in the saved prediction file."
            )

        d["sex_male"] = normalize_sex(
            d[sex_col]
        )

        n_unique_sex = int(
            d["sex_male"].dropna().nunique()
        )

        if n_unique_sex < 2:
            raise ValueError(
                "covariate-mode age_sex requires both mapped female and male "
                "participants in the held-out test data. Found {} mapped sex "
                "level(s). No clock is automatically treated as sex-restricted "
                "based on its name.".format(n_unique_sex)
            )

        covariates.append(
            "sex_male"
        )

    else:
        age_mean = np.nan
        age_sd = np.nan
        n_unique_sex = 0

    model_cols = [
        "participant_id",
        "time_years",
        "event_bool",
        "epoch_predictor",
    ] + covariates

    model_df = d[
        model_cols
    ].copy()

    for c in [
        "time_years",
        "epoch_predictor",
    ] + covariates:
        model_df[c] = pd.to_numeric(
            model_df[c],
            errors="coerce",
        )

    model_df = model_df.replace(
        [np.inf, -np.inf],
        np.nan,
    )

    model_df = model_df.dropna().copy()

    model_df = model_df.loc[
        model_df["time_years"] >
        0
    ].copy()

    # lifelines residual routines are easier to use with a clean integer index.
    model_df = model_df.reset_index(
        drop=True
    )

    if model_df.shape[0] < 20:
        raise ValueError(
            "Fewer than 20 complete test participants for PH analysis."
        )

    if int(
        model_df[
            "event_bool"
        ].sum()
    ) < 5:
        raise ValueError(
            "Fewer than 5 deaths in the held-out PH analysis sample."
        )

    if args.covariate_mode == "age_sex":
        n_female_complete = int(
            (model_df["sex_male"] == 0.0).sum()
        )
        n_male_complete = int(
            (model_df["sex_male"] == 1.0).sum()
        )
        age_covariate_used = True
        sex_covariate_used = True

        if n_female_complete == 0 or n_male_complete == 0:
            raise ValueError(
                "After complete-case filtering, both female and male "
                "participants must remain for the age+sex PH model. "
                "Female N={}, male N={}.".format(
                    n_female_complete,
                    n_male_complete,
                )
            )
    else:
        n_female_complete = np.nan
        n_male_complete = np.nan
        age_covariate_used = False
        sex_covariate_used = False

    meta = {
        "predictor_source_col": predictor_source,
        "predictor_mode_requested": args.predictor_mode,
        "predictor_scale_note": predictor_scale_note,
        "acceleration_col_detected": acceleration_col,
        "risk_score_col_detected": risk_score_col,
        "age_col_detected": age_col,
        "sex_col_detected": sex_col,
        "covariates_used": ";".join(covariates),
        "age_mean_before_z": age_mean,
        "age_sd_before_z": age_sd,
        "n_unique_sex_mapped": n_unique_sex,
        "n_female_test_complete_case": n_female_complete,
        "n_male_test_complete_case": n_male_complete,
        "age_covariate_used": age_covariate_used,
        "sex_covariate_used": sex_covariate_used,
        "sex_handling_note": (
            "No participant is filtered by clock name. Reproductive female "
            "and reproductive male proteomics clocks retain both sexes; sex "
            "is modeled as a covariate exactly as for the other clocks."
            if args.covariate_mode == "age_sex"
            else "Unadjusted sensitivity analysis; age/sex covariates disabled."
        ),
    }

    return model_df, meta


# =============================================================================
# 5. Standard Cox + Schoenfeld PH diagnostics
# =============================================================================

def fit_cox_with_fallback(model_df):
    covariate_cols = [
        c
        for c in model_df.columns
        if c not in {
            "participant_id",
            "time_years",
            "event_bool",
        }
    ]

    fit_df = model_df[
        [
            "time_years",
            "event_bool",
        ] + covariate_cols
    ].copy()

    last_error = None

    for penalizer in [
        0.0,
        1e-8,
        1e-6,
        1e-5,
        1e-4,
    ]:
        try:
            cph = CoxPHFitter(
                penalizer=penalizer
            )

            cph.fit(
                fit_df,
                duration_col="time_years",
                event_col="event_bool",
                show_progress=False,
            )

            return (
                cph,
                fit_df,
                penalizer,
            )

        except Exception as exc:
            last_error = exc

    raise RuntimeError(
        "CoxPHFitter failed for all penalizers. Last error: {}".format(
            last_error
        )
    )


def extract_epoch_cox_stats(cph):
    term = "epoch_predictor"

    if term not in cph.summary.index:
        raise KeyError(
            "epoch_predictor is missing from CoxPHFitter summary."
        )

    row = cph.summary.loc[
        term
    ]

    return {
        "epoch_coef": safe_float(
            row.get("coef")
        ),
        "epoch_se": safe_float(
            row.get("se(coef)")
        ),
        "epoch_hr": safe_float(
            row.get("exp(coef)")
        ),
        "epoch_hr_ci_lower": safe_float(
            row.get(
                "exp(coef) lower 95%"
            )
        ),
        "epoch_hr_ci_upper": safe_float(
            row.get(
                "exp(coef) upper 95%"
            )
        ),
        "epoch_wald_p": safe_float(
            row.get("p")
        ),
    }


def schoenfeld_diagnostics(
    cph,
    fit_df,
    time_transform,
):
    ph = proportional_hazard_test(
        cph,
        fit_df,
        time_transform=time_transform,
    )

    ph_summary = ph.summary.copy()

    if "epoch_predictor" not in ph_summary.index:
        raise KeyError(
            "epoch_predictor is missing from proportional_hazard_test output."
        )

    epoch_row = ph_summary.loc[
        "epoch_predictor"
    ]

    covariate_rows = ph_summary.reset_index().rename(
        columns={
            "index": "term",
            "test_statistic": "ph_chisq",
            "p": "ph_p",
        }
    )

    covariate_rows["time_transform"] = time_transform

    # Descriptive residual-time trend. This is NOT the formal PH test.
    residuals = cph.compute_residuals(
        fit_df,
        kind="scaled_schoenfeld",
    )

    if "epoch_predictor" not in residuals.columns:
        raise KeyError(
            "epoch_predictor is missing from scaled Schoenfeld residuals."
        )

    event_index = residuals.index

    event_times = pd.to_numeric(
        fit_df.loc[
            event_index,
            "time_years",
        ],
        errors="coerce",
    ).to_numpy(
        dtype=float
    )

    epoch_residual = pd.to_numeric(
        residuals[
            "epoch_predictor"
        ],
        errors="coerce",
    ).to_numpy(
        dtype=float
    )

    ok = (
        np.isfinite(event_times)
        & (
            event_times >
            0
        )
        & np.isfinite(
            epoch_residual
        )
    )

    if int(
        np.sum(ok)
    ) >= 3:
        rho, rho_p = spearmanr(
            np.log(
                event_times[ok]
            ),
            epoch_residual[ok],
        )

        rho = safe_float(
            rho
        )

        rho_p = safe_float(
            rho_p
        )

    else:
        rho = np.nan
        rho_p = np.nan

    residual_df = pd.DataFrame(
        {
            "event_time_years": event_times,
            "log_event_time": np.where(
                event_times > 0,
                np.log(event_times),
                np.nan,
            ),
            "scaled_schoenfeld_residual_epoch": epoch_residual,
        }
    )

    result = {
        "ph_time_transform": time_transform,
        "epoch_ph_chisq": safe_float(
            epoch_row.get(
                "test_statistic"
            )
        ),
        "epoch_ph_p": safe_float(
            epoch_row.get("p")
        ),
        "epoch_scaled_schoenfeld_spearman_rho_log_time": rho,
        "epoch_scaled_schoenfeld_spearman_rho_p_descriptive": rho_p,
    }

    return (
        result,
        covariate_rows,
        residual_df,
    )


def save_schoenfeld_plot(
    residual_df,
    clock_name,
    p_value,
    png_path,
    pdf_path,
):
    d = residual_df.loc[
        np.isfinite(
            residual_df[
                "event_time_years"
            ]
        )
        & np.isfinite(
            residual_df[
                "scaled_schoenfeld_residual_epoch"
            ]
        )
    ].copy()

    if d.shape[0] < 3:
        return

    fig, ax = plt.subplots(
        figsize=(6.4, 4.8)
    )

    ax.axhline(
        0.0,
        linestyle="--",
        linewidth=1.0,
        color="0.45",
    )

    ax.scatter(
        d[
            "event_time_years"
        ].values,
        d[
            "scaled_schoenfeld_residual_epoch"
        ].values,
        s=14,
        alpha=0.28,
        edgecolors="none",
        label="Scaled Schoenfeld residuals",
    )

    # Quantile-bin mean trend: dependency-free visual aid only.
    q = min(
        10,
        max(
            3,
            int(
                d.shape[0] // 20
            ),
        ),
    )

    try:
        d["time_bin"] = pd.qcut(
            d[
                "event_time_years"
            ].rank(
                method="first"
            ),
            q=q,
            labels=False,
            duplicates="drop",
        )

        trend = (
            d.groupby(
                "time_bin",
                observed=True,
            )
            .agg(
                event_time_years=(
                    "event_time_years",
                    "median",
                ),
                residual=(
                    "scaled_schoenfeld_residual_epoch",
                    "mean",
                ),
            )
            .reset_index(
                drop=True
            )
        )

        if trend.shape[0] >= 2:
            ax.plot(
                trend[
                    "event_time_years"
                ].values,
                trend[
                    "residual"
                ].values,
                marker="o",
                linewidth=1.7,
                label="Binned mean trend",
            )

    except Exception:
        pass

    p_label = (
        "NA"
        if not np.isfinite(
            p_value
        )
        else "{:.3g}".format(
            p_value
        )
    )

    ax.set_xlabel(
        "Event time (years)"
    )

    ax.set_ylabel(
        "Scaled Schoenfeld residual for EPOCH"
    )

    ax.set_title(
        "{} mortality EPOCH\nSchoenfeld PH test P={}".format(
            clock_name,
            p_label,
        )
    )

    ax.legend(
        frameon=False
    )

    fig.tight_layout()

    fig.savefig(
        png_path,
        dpi=300,
        bbox_inches="tight",
    )

    fig.savefig(
        pdf_path,
        bbox_inches="tight",
    )

    plt.close(
        fig
    )


# =============================================================================
# 6. Piecewise time-varying EPOCH sensitivity analysis
# =============================================================================

def make_two_band_start_stop(
    model_df,
    split_years,
):
    rows = []

    static_covariates = [
        c
        for c in model_df.columns
        if c not in {
            "participant_id",
            "time_years",
            "event_bool",
        }
    ]

    for row in model_df.itertuples(
        index=False
    ):
        record = row._asdict()

        pid = str(
            record[
                "participant_id"
            ]
        )

        duration = float(
            record[
                "time_years"
            ]
        )

        event = bool(
            record[
                "event_bool"
            ]
        )

        base_cov = {
            c: float(
                record[c]
            )
            for c in static_covariates
        }

        if duration <= split_years:
            first = {
                "participant_id": pid,
                "start": 0.0,
                "stop": duration,
                "event_bool": event,
                **base_cov,
            }

            first[
                "epoch_after_split"
            ] = 0.0

            rows.append(
                first
            )

        else:
            first = {
                "participant_id": pid,
                "start": 0.0,
                "stop": float(
                    split_years
                ),
                "event_bool": False,
                **base_cov,
            }

            first[
                "epoch_after_split"
            ] = 0.0

            rows.append(
                first
            )

            second = {
                "participant_id": pid,
                "start": float(
                    split_years
                ),
                "stop": duration,
                "event_bool": event,
                **base_cov,
            }

            second[
                "epoch_after_split"
            ] = float(
                record[
                    "epoch_predictor"
                ]
            )

            rows.append(
                second
            )

    long_df = pd.DataFrame(
        rows
    )

    long_df = long_df.loc[
        long_df["stop"] >
        long_df["start"]
    ].copy()

    return long_df


def fit_ctv_with_fallback(
    long_df,
    include_interaction,
):
    base_cols = [
        "participant_id",
        "start",
        "stop",
        "event_bool",
        "epoch_predictor",
    ]

    covariates = [
        c
        for c in [
            "age_covariate_z",
            "sex_male",
        ]
        if c in long_df.columns
    ]

    cols = (
        base_cols
        + covariates
        + (
            [
                "epoch_after_split"
            ]
            if include_interaction
            else []
        )
    )

    d = long_df[
        cols
    ].copy()

    last_error = None

    for penalizer in [
        0.0,
        1e-8,
        1e-6,
        1e-5,
        1e-4,
    ]:
        try:
            model = CoxTimeVaryingFitter(
                penalizer=penalizer
            )

            model.fit(
                d,
                id_col="participant_id",
                event_col="event_bool",
                start_col="start",
                stop_col="stop",
                show_progress=False,
            )

            return (
                model,
                penalizer,
            )

        except Exception as exc:
            last_error = exc

    raise RuntimeError(
        "CoxTimeVaryingFitter failed for all penalizers. Last error: {}".format(
            last_error
        )
    )


def linear_combination_hr(
    params,
    variance,
    weights,
):
    names = list(
        weights.keys()
    )

    beta = sum(
        float(
            weights[name]
        )
        * float(
            params.loc[
                name
            ]
        )
        for name in names
    )

    w = np.asarray(
        [
            float(
                weights[name]
            )
            for name in names
        ],
        dtype=float,
    )

    V = variance.loc[
        names,
        names,
    ].to_numpy(
        dtype=float
    )

    var_beta = float(
        w @ V @ w
    )

    se = (
        math.sqrt(
            max(
                var_beta,
                0.0,
            )
        )
        if np.isfinite(
            var_beta
        )
        else np.nan
    )

    z = 1.959963984540054

    return {
        "coef": beta,
        "se": se,
        "hr": math.exp(
            beta
        ),
        "hr_ci_lower": (
            math.exp(
                beta -
                z * se
            )
            if np.isfinite(
                se
            )
            else np.nan
        ),
        "hr_ci_upper": (
            math.exp(
                beta +
                z * se
            )
            if np.isfinite(
                se
            )
            else np.nan
        ),
    }


def time_varying_sensitivity(
    model_df,
    split_years,
    min_events_per_band,
):
    event = model_df[
        "event_bool"
    ].astype(bool)

    time = model_df[
        "time_years"
    ].astype(float)

    n_early = int(
        (
            event
            & (
                time <=
                split_years
            )
        ).sum()
    )

    n_late = int(
        (
            event
            & (
                time >
                split_years
            )
        ).sum()
    )

    base = {
        "time_split_years": float(
            split_years
        ),
        "n_events_on_or_before_time_split": n_early,
        "n_events_after_time_split": n_late,
    }

    if (
        n_early <
        min_events_per_band
        or n_late <
        min_events_per_band
    ):
        return {
            **base,
            "time_varying_status": (
                "skipped: insufficient events in one or both time bands"
            ),
            "time_varying_reduced_penalizer": np.nan,
            "time_varying_full_penalizer": np.nan,
            "epoch_hr_early": np.nan,
            "epoch_hr_early_ci_lower": np.nan,
            "epoch_hr_early_ci_upper": np.nan,
            "epoch_hr_late": np.nan,
            "epoch_hr_late_ci_lower": np.nan,
            "epoch_hr_late_ci_upper": np.nan,
            "epoch_after_split_interaction_coef": np.nan,
            "epoch_after_split_interaction_hr_ratio": np.nan,
            "epoch_after_split_interaction_p": np.nan,
            "time_varying_lrt_chisq": np.nan,
            "time_varying_lrt_p": np.nan,
            "time_varying_lrt_note": (
                "Not fit because one or both time bands had fewer than "
                "{} deaths.".format(
                    min_events_per_band
                )
            ),
        }

    long_df = make_two_band_start_stop(
        model_df,
        split_years,
    )

    reduced, reduced_pen = fit_ctv_with_fallback(
        long_df,
        include_interaction=False,
    )

    full, full_pen = fit_ctv_with_fallback(
        long_df,
        include_interaction=True,
    )

    early = linear_combination_hr(
        full.params_,
        full.variance_matrix_,
        {
            "epoch_predictor": 1.0,
        },
    )

    late = linear_combination_hr(
        full.params_,
        full.variance_matrix_,
        {
            "epoch_predictor": 1.0,
            "epoch_after_split": 1.0,
        },
    )

    interaction_coef = safe_float(
        full.params_.loc[
            "epoch_after_split"
        ]
    )

    interaction_p = safe_float(
        full.summary.loc[
            "epoch_after_split",
            "p",
        ]
    )

    if (
        reduced_pen == 0.0
        and full_pen == 0.0
    ):
        lrt = max(
            0.0,
            2.0 * (
                float(
                    full.log_likelihood_
                )
                -
                float(
                    reduced.log_likelihood_
                )
            ),
        )

        lrt_p = float(
            chi2.sf(
                lrt,
                1,
            )
        )

        lrt_note = (
            "Valid 1-df likelihood-ratio comparison of unpenalized "
            "piecewise time-varying versus proportional model."
        )

    else:
        lrt = np.nan
        lrt_p = np.nan

        lrt_note = (
            "LRT not reported because at least one CoxTimeVaryingFitter "
            "required penalization; use the interaction Wald P instead."
        )

    return {
        **base,
        "time_varying_status": "success",
        "time_varying_reduced_penalizer": reduced_pen,
        "time_varying_full_penalizer": full_pen,
        "epoch_hr_early": early["hr"],
        "epoch_hr_early_ci_lower": early["hr_ci_lower"],
        "epoch_hr_early_ci_upper": early["hr_ci_upper"],
        "epoch_hr_late": late["hr"],
        "epoch_hr_late_ci_lower": late["hr_ci_lower"],
        "epoch_hr_late_ci_upper": late["hr_ci_upper"],
        "epoch_after_split_interaction_coef": interaction_coef,
        "epoch_after_split_interaction_hr_ratio": (
            math.exp(
                interaction_coef
            )
            if np.isfinite(
                interaction_coef
            )
            else np.nan
        ),
        "epoch_after_split_interaction_p": interaction_p,
        "time_varying_lrt_chisq": lrt,
        "time_varying_lrt_p": lrt_p,
        "time_varying_lrt_note": lrt_note,
    }


# =============================================================================
# 7. Analyze one clock
# =============================================================================

def analyze_clock(
    clock_info,
    base_dir,
    master_dir,
    args,
):
    clock_dir = (
        base_dir /
        clock_info[
            "folder"
        ]
    )

    if not clock_dir.exists():
        raise FileNotFoundError(
            "Clock directory does not exist: {}".format(
                clock_dir
            )
        )

    pred_file = resolve_prediction_file(
        clock_dir
    )

    prefix = sanitize_prefix(
        pred_file.name
    )

    per_clock_dir = (
        clock_dir /
        "ph_diagnostics"
    )

    per_clock_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    info("")
    info("=" * 96)
    info(
        "Clock: {}".format(
            clock_info[
                "clock"
            ]
        )
    )
    info(
        "Modality: {}".format(
            clock_info[
                "modality"
            ]
        )
    )
    info(
        "Predictions: {}".format(
            pred_file
        )
    )
    info("=" * 96)

    df = pd.read_csv(
        pred_file,
        sep="\t",
        low_memory=False,
    )

    model_df, meta = prepare_analysis_data(
        df,
        args,
    )

    n_test = int(
        model_df.shape[0]
    )

    n_events = int(
        model_df[
            "event_bool"
        ].sum()
    )

    median_followup = float(
        model_df[
            "time_years"
        ].median()
    )

    max_followup = float(
        model_df[
            "time_years"
        ].max()
    )

    cph, fit_df, cox_penalizer = fit_cox_with_fallback(
        model_df
    )

    cox_stats = extract_epoch_cox_stats(
        cph
    )

    (
        ph_stats,
        ph_covariate_df,
        residual_df,
    ) = schoenfeld_diagnostics(
        cph,
        fit_df,
        args.time_transform,
    )

    save_schoenfeld_plot(
        residual_df,
        clock_name=clock_info[
            "clock"
        ],
        p_value=ph_stats[
            "epoch_ph_p"
        ],
        png_path=(
            per_clock_dir /
            "{}_scaled_schoenfeld_epoch.png".format(
                prefix
            )
        ),
        pdf_path=(
            per_clock_dir /
            "{}_scaled_schoenfeld_epoch.pdf".format(
                prefix
            )
        ),
    )

    residual_df.to_csv(
        per_clock_dir /
        "{}_scaled_schoenfeld_epoch.tsv".format(
            prefix
        ),
        sep="\t",
        index=False,
        na_rep="NA",
    )

    ph_covariate_df.insert(
        0,
        "clock",
        clock_info[
            "clock"
        ],
    )

    ph_covariate_df.insert(
        1,
        "modality",
        clock_info[
            "modality"
        ],
    )

    ph_covariate_df.insert(
        2,
        "folder",
        clock_info[
            "folder"
        ],
    )

    ph_covariate_df.to_csv(
        per_clock_dir /
        "{}_PH_covariate_tests.tsv".format(
            prefix
        ),
        sep="\t",
        index=False,
        na_rep="NA",
    )

    if args.skip_time_varying:
        tv_stats = {
            "time_split_years": float(
                args.time_split_years
            ),
            "n_events_on_or_before_time_split": np.nan,
            "n_events_after_time_split": np.nan,
            "time_varying_status": "skipped by --skip-time-varying",
            "time_varying_reduced_penalizer": np.nan,
            "time_varying_full_penalizer": np.nan,
            "epoch_hr_early": np.nan,
            "epoch_hr_early_ci_lower": np.nan,
            "epoch_hr_early_ci_upper": np.nan,
            "epoch_hr_late": np.nan,
            "epoch_hr_late_ci_lower": np.nan,
            "epoch_hr_late_ci_upper": np.nan,
            "epoch_after_split_interaction_coef": np.nan,
            "epoch_after_split_interaction_hr_ratio": np.nan,
            "epoch_after_split_interaction_p": np.nan,
            "time_varying_lrt_chisq": np.nan,
            "time_varying_lrt_p": np.nan,
            "time_varying_lrt_note": "Skipped by user request.",
        }

    else:
        tv_stats = time_varying_sensitivity(
            model_df,
            split_years=args.time_split_years,
            min_events_per_band=args.min_events_per_time_band,
        )

    tv_df = pd.DataFrame(
        [
            {
                "clock": clock_info[
                    "clock"
                ],
                "modality": clock_info[
                    "modality"
                ],
                "folder": clock_info[
                    "folder"
                ],
                **tv_stats,
            }
        ]
    )

    tv_df.to_csv(
        per_clock_dir /
        "{}_time_varying_{:g}y.tsv".format(
            prefix,
            float(
                args.time_split_years
            ),
        ),
        sep="\t",
        index=False,
        na_rep="NA",
    )

    row = {
        "clock": clock_info[
            "clock"
        ],
        "modality": clock_info[
            "modality"
        ],
        "folder": clock_info[
            "folder"
        ],
        "prediction_file": str(
            pred_file
        ),
        "test_split": args.test_split,
        "n_test_complete_case": n_test,
        "n_events_test": n_events,
        "median_followup_test_years": median_followup,
        "max_followup_test_years": max_followup,
        "predictor_mode_requested": meta[
            "predictor_mode_requested"
        ],
        "predictor_source_col": meta[
            "predictor_source_col"
        ],
        "predictor_scale_note": meta[
            "predictor_scale_note"
        ],
        "acceleration_col_detected": meta[
            "acceleration_col_detected"
        ],
        "risk_score_col_detected": meta[
            "risk_score_col_detected"
        ],
        "age_col_detected": meta[
            "age_col_detected"
        ],
        "sex_col_detected": meta[
            "sex_col_detected"
        ],
        "covariate_mode": args.covariate_mode,
        "covariates_used": meta[
            "covariates_used"
        ],
        "age_covariate_used": meta[
            "age_covariate_used"
        ],
        "sex_covariate_used": meta[
            "sex_covariate_used"
        ],
        "n_female_test_complete_case": meta[
            "n_female_test_complete_case"
        ],
        "n_male_test_complete_case": meta[
            "n_male_test_complete_case"
        ],
        "sex_handling_note": meta[
            "sex_handling_note"
        ],
        "cox_penalizer_used": cox_penalizer,
        **cox_stats,
        **ph_stats,
        **tv_stats,
    }

    summary_df = pd.DataFrame(
        [
            row
        ]
    )

    summary_df.to_csv(
        per_clock_dir /
        "{}_PH_summary.tsv".format(
            prefix
        ),
        sep="\t",
        index=False,
        na_rep="NA",
    )

    info(
        "N={:,}; deaths={:,}; HR={:.3f} [{:.3f}, {:.3f}]; "
        "Schoenfeld P={}; time-varying interaction P={}".format(
            n_test,
            n_events,
            row[
                "epoch_hr"
            ]
            if np.isfinite(
                row[
                    "epoch_hr"
                ]
            )
            else np.nan,
            row[
                "epoch_hr_ci_lower"
            ]
            if np.isfinite(
                row[
                    "epoch_hr_ci_lower"
                ]
            )
            else np.nan,
            row[
                "epoch_hr_ci_upper"
            ]
            if np.isfinite(
                row[
                    "epoch_hr_ci_upper"
                ]
            )
            else np.nan,
            "{:.4g}".format(
                row[
                    "epoch_ph_p"
                ]
            )
            if np.isfinite(
                row[
                    "epoch_ph_p"
                ]
            )
            else "NA",
            "{:.4g}".format(
                row[
                    "epoch_after_split_interaction_p"
                ]
            )
            if np.isfinite(
                row[
                    "epoch_after_split_interaction_p"
                ]
            )
            else "NA",
        )
    )

    manifest = {
        "clock": clock_info[
            "clock"
        ],
        "modality": clock_info[
            "modality"
        ],
        "folder": clock_info[
            "folder"
        ],
        "status": "success",
        "prediction_file": str(
            pred_file
        ),
        "n_test": n_test,
        "n_events": n_events,
        "output_dir": str(
            per_clock_dir
        ),
        "error": "",
    }

    return (
        summary_df,
        ph_covariate_df,
        manifest,
    )


# =============================================================================
# 8. Main batch
# =============================================================================

def main():
    args = parse_args()

    base_dir = Path(
        args.base_dir
    ).resolve()

    if args.output_dir:
        master_dir = Path(
            args.output_dir
        ).resolve()
    else:
        master_dir = (
            base_dir /
            "mortality_EPOCH_PH_diagnostics"
        )

    master_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    if len(CLOCKS) != 23:
        raise RuntimeError(
            "Internal manifest does not contain exactly 23 clocks."
        )

    if (
        not np.isfinite(
            args.time_split_years
        )
        or args.time_split_years <= 0
    ):
        raise ValueError(
            "--time-split-years must be > 0."
        )

    all_summary = []
    all_covariate_tests = []
    manifest = []

    info("=" * 96)
    info("POST-HOC PH DIAGNOSTICS FOR 23 MORTALITY EPOCH CLOCKS")
    info("=" * 96)
    info(
        "Base directory: {}".format(
            base_dir
        )
    )
    info(
        "Master output: {}".format(
            master_dir
        )
    )
    info(
        "Test split: {}".format(
            args.test_split
        )
    )
    info(
        "Primary predictor mode: {}".format(
            args.predictor_mode
        )
    )
    info(
        "Covariate mode: {}".format(
            args.covariate_mode
        )
    )
    if args.covariate_mode == "age_sex":
        info(
            "Sex handling: both female and male participants are required; "
            "clock names never trigger sex-based filtering."
        )
    info(
        "Schoenfeld time transform: {}".format(
            args.time_transform
        )
    )
    info(
        "Time-varying split: {:g} years".format(
            args.time_split_years
        )
    )
    info(
        "Time-varying sensitivity: {}".format(
            "SKIPPED"
            if args.skip_time_varying
            else "ENABLED"
        )
    )
    info("=" * 96)

    for i, clock_info in enumerate(
        CLOCKS,
        start=1,
    ):
        info("")
        info(
            "[{:02d}/23] {}".format(
                i,
                clock_info[
                    "clock"
                ],
            )
        )

        try:
            summary_df, cov_df, m = analyze_clock(
                clock_info,
                base_dir,
                master_dir,
                args,
            )

            all_summary.append(
                summary_df
            )

            all_covariate_tests.append(
                cov_df
            )

            manifest.append(
                m
            )

        except Exception as exc:
            msg = "{}".format(
                exc
            )

            warn(
                "{} FAILED: {}".format(
                    clock_info[
                        "clock"
                    ],
                    msg,
                )
            )

            manifest.append(
                {
                    "clock": clock_info[
                        "clock"
                    ],
                    "modality": clock_info[
                        "modality"
                    ],
                    "folder": clock_info[
                        "folder"
                    ],
                    "status": "failed",
                    "prediction_file": "",
                    "n_test": np.nan,
                    "n_events": np.nan,
                    "output_dir": "",
                    "error": msg,
                }
            )

            if args.fail_on_missing_clock:
                raise

    manifest_df = pd.DataFrame(
        manifest
    )

    manifest_path = (
        master_dir /
        "mortality_EPOCH_23_PH_run_manifest.tsv"
    )

    manifest_df.to_csv(
        manifest_path,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    if not all_summary:
        raise RuntimeError(
            "PH diagnostics failed for all clocks. See run manifest."
        )

    master = pd.concat(
        all_summary,
        ignore_index=True,
    )

    # BH-FDR across the 23 EPOCH-specific Schoenfeld tests.
    master[
        "epoch_ph_p_fdr_bh_23"
    ] = bh_fdr(
        master[
            "epoch_ph_p"
        ].values
    )

    master[
        "ph_flag_raw_p_lt_0_05"
    ] = (
        master[
            "epoch_ph_p"
        ] < 0.05
    )

    master[
        "ph_flag_fdr_lt_0_05"
    ] = (
        master[
            "epoch_ph_p_fdr_bh_23"
        ] < 0.05
    )

    master_path = (
        master_dir /
        "mortality_EPOCH_23_PH_diagnostics.tsv"
    )

    master.to_csv(
        master_path,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    if all_covariate_tests:
        covariate_master = pd.concat(
            all_covariate_tests,
            ignore_index=True,
        )

        covariate_master_path = (
            master_dir /
            "mortality_EPOCH_23_PH_covariate_tests.tsv"
        )

        covariate_master.to_csv(
            covariate_master_path,
            sep="\t",
            index=False,
            na_rep="NA",
        )
    else:
        covariate_master_path = None

    report_cols = [
        "clock",
        "modality",
        "n_test_complete_case",
        "n_events_test",
        "median_followup_test_years",
        "predictor_source_col",
        "covariates_used",
        "age_covariate_used",
        "sex_covariate_used",
        "n_female_test_complete_case",
        "n_male_test_complete_case",
        "epoch_hr",
        "epoch_hr_ci_lower",
        "epoch_hr_ci_upper",
        "epoch_wald_p",
        "ph_time_transform",
        "epoch_ph_chisq",
        "epoch_ph_p",
        "epoch_ph_p_fdr_bh_23",
        "ph_flag_raw_p_lt_0_05",
        "ph_flag_fdr_lt_0_05",
        "epoch_scaled_schoenfeld_spearman_rho_log_time",
        "time_split_years",
        "n_events_on_or_before_time_split",
        "n_events_after_time_split",
        "epoch_hr_early",
        "epoch_hr_early_ci_lower",
        "epoch_hr_early_ci_upper",
        "epoch_hr_late",
        "epoch_hr_late_ci_lower",
        "epoch_hr_late_ci_upper",
        "epoch_after_split_interaction_hr_ratio",
        "epoch_after_split_interaction_p",
        "time_varying_lrt_chisq",
        "time_varying_lrt_p",
        "time_varying_status",
    ]

    manuscript_table = master[
        [
            c
            for c in report_cols
            if c in master.columns
        ]
    ].copy()

    manuscript_path = (
        master_dir /
        "mortality_EPOCH_23_PH_manuscript_table.tsv"
    )

    manuscript_table.to_csv(
        manuscript_path,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    n_success = int(
        (
            manifest_df[
                "status"
            ] == "success"
        ).sum()
    )

    n_failed = int(
        (
            manifest_df[
                "status"
            ] == "failed"
        ).sum()
    )

    n_raw = int(
        master[
            "ph_flag_raw_p_lt_0_05"
        ].fillna(
            False
        ).sum()
    )

    n_fdr = int(
        master[
            "ph_flag_fdr_lt_0_05"
        ].fillna(
            False
        ).sum()
    )

    info("")
    info("=" * 96)
    info("PH DIAGNOSTIC ANALYSIS FINISHED")
    info("=" * 96)
    info(
        "Successful clocks: {}/23".format(
            n_success
        )
    )
    info(
        "Failed clocks: {}".format(
            n_failed
        )
    )
    info(
        "EPOCH PH raw P<0.05: {}/{}".format(
            n_raw,
            master.shape[0],
        )
    )
    info(
        "EPOCH PH BH-FDR<0.05: {}/{}".format(
            n_fdr,
            master.shape[0],
        )
    )
    info("Master files:")
    info(
        "  {}".format(
            master_path
        )
    )
    info(
        "  {}".format(
            manuscript_path
        )
    )

    if covariate_master_path is not None:
        info(
            "  {}".format(
                covariate_master_path
            )
        )

    info(
        "  {}".format(
            manifest_path
        )
    )
    info("=" * 96)


if __name__ == "__main__":
    main()
