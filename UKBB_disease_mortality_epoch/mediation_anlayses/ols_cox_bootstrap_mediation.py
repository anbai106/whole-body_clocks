#!/usr/bin/env python3
"""
EPOCH statistical mediation: OLS + Cox proportional hazards + bootstrap
=====================================================================

METHOD PIPELINE
---------------
Scientific question
    For each baseline molecular disease EPOCH acceleration-z score, test whether
    its association with later all-cause mortality is statistically mediated by
    an MRI mortality EPOCH acceleration-z score measured at the imaging visit.

Temporal ordering
    X (exposure): baseline proteomics/metabolomics disease EPOCH acceleration-z
        -> M (mediator): MRI mortality EPOCH acceleration-z at imaging visit 2
        -> Y (outcome): observed all-cause mortality after the MRI landmark.

Model grid
    45 molecular disease EPOCH exposures x 7 MRI mortality EPOCH mediators =
    315 exposure-mediator models when the current wide table contains the same
    significant clocks used in the SEM analysis.

Population
    FULL available population is used by default. The optional --test-only flag
    restricts each MRI mortality prediction file to rows with split == "test".
    The split restriction is applied before construction of the mortality
    landmark dataset. Models then use complete cases after merging exposure,
    mediator, mortality follow-up, and covariates.

Landmark survival outcome
    MRI imaging_date is preferred as time zero; sample_date is accepted as a
    fallback for older files. A death is counted only when death_date is after
    the MRI landmark and on/before admin_censor_date. Survival time is measured
    from MRI landmark to death or administrative censoring. Unlike the previous
    fixed-5-year WLSMV SEM, participants with shorter valid follow-up remain in
    the Cox analysis and contribute their observed censored follow-up.

Covariate adjustment
    Default covariates match the prior SEM as closely as possible:
      - imaging-visit age
      - imaging-visit BMI
      - sex
      - imaging-visit smoking status
      - imaging-visit assessment centre
      - genetic PCs 1-10
    Categorical covariates are one-hot encoded. Zero-variance dummy variables are
    removed within each exposure-mediator model. Complete-case filtering occurs
    before model fitting.

Standardization
    Exposure X and mediator M are re-standardized to mean 0 and SD 1 within the
    final complete-case model sample. Thus a and b correspond to one-SD changes
    in X and M within that model's analysis sample.

Path a: OLS mediator model
    M = alpha + a*X + covariates + error
    OLS is fitted with HC3 heteroskedasticity-robust standard errors.

Paths b and c': Cox survival model
    h(t) = h0(t) * exp(b*M + c'*X + covariates)
    b is the conditional MRI-mediator log-hazard coefficient and c' is the
    direct exposure log-hazard coefficient conditional on M.

Total effect c: Cox survival model
    h(t) = h0(t) * exp(c*X + covariates)

Indirect association
    The statistical indirect association is approximated on the Cox log-hazard
    scale by a*b. exp(a*b) is also reported as an indirect-association HR-like
    transform. This product-of-coefficients quantity is a statistical mediation
    approximation, NOT a formally identified causal natural indirect effect.
    Cox coefficients are non-collapsible, so c need not equal c' + a*b.
    Proportion mediated is intentionally NOT calculated.

Uncertainty for a*b
    1) Delta/Sobel approximation:
         SE(a*b) = sqrt(b^2*SE(a)^2 + a^2*SE(b)^2)
       This provides a continuous two-sided p-value suitable for stringent
       multiple-testing correction.
    2) Participant-level nonparametric bootstrap:
       participants are resampled with replacement; the OLS a-path and the Cox
       b/c' model are refitted; a*b is saved. Percentile 95% CIs and an empirical
       two-sided bootstrap p-value are reported when enough replicates succeed.
       Bootstrap p-values have finite resolution, so Bonferroni inference uses
       the delta-method indirect p-value by default, while bootstrap CIs are used
       as the main robustness/uncertainty check.

Multiple testing
    The primary multiplicity family is all requested exposure x mediator models
    (normally 315). Bonferroni-adjusted p-values are reported for the indirect
    delta-method p-values using m = number of requested models, including models
    that later fail or are skipped. The script also reports the exact familywise
    threshold alpha/m. BH-FDR values are saved as secondary information only.

Diagnostics and safeguards
    - minimum complete-case N and minimum number of deaths
    - Cox convergence/error capture
    - maximum VIF for X, M, and covariates
    - Schoenfeld-residual proportional-hazards tests for X and M when available
    - bootstrap success count and success fraction
    - per-model QC/status table and captured warnings
    - no proportion-mediated statistic

Primary interpretation
    Prioritize models that:
      1) fit without Cox convergence problems,
      2) have adequate numbers of deaths,
      3) have Bonferroni-significant indirect delta p-values,
      4) have bootstrap CIs excluding zero,
      5) have acceptable VIF/PH diagnostics,
      6) show compatible direction in the optional held-out test analysis.

Required packages
    numpy pandas scipy statsmodels lifelines

Typical full-population run
    python run_epoch_ols_cox_bootstrap_mediation.py \
      --bootstrap 1000

Held-out test-only sensitivity run
    python run_epoch_ols_cox_bootstrap_mediation.py \
      --test-only \
      --bootstrap 1000 \
      --output-dir /path/to/test_only_output
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
import traceback
import warnings
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import statsmodels.api as sm
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test
from scipy import stats
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.outliers_influence import variance_inflation_factor


# -----------------------------------------------------------------------------
# Defaults matching the current CUBIC project layout
# -----------------------------------------------------------------------------
DEFAULT_ROOT = Path("/cbica/home/wenju/Reproducibile_paper/WholeBodyClock")
DEFAULT_EPOCH_WIDE = DEFAULT_ROOT / "collected_significant_epoch_clocks" / "significant_epoch_clocks_wide.tsv"
DEFAULT_COVARIATES = Path(
    "/cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv"
)
DEFAULT_OUTPUT_FULL = DEFAULT_ROOT / "mediation_OLS_Cox_bootstrap_full"
DEFAULT_OUTPUT_TEST = DEFAULT_ROOT / "mediation_OLS_Cox_bootstrap_test"

MRI_MORTALITY_ORGANS = ["adipose", "brain", "heart", "kidney", "liver", "pancreas", "spleen"]

CONTINUOUS_COVARIATE_MAP = {
    "age_mri": "age_when_attended_assessment_centre_f21003_2_0",
    "bmi_mri": "body_mass_index_bmi_f23104_2_0",
}
CATEGORICAL_COVARIATE_MAP = {
    "sex": "sex_f31_0_0",
    "smoking_mri": "smoking_status_f20116_2_0",
    "centre_mri": "uk_biobank_assessment_centre_f54_2_0",
}


@dataclass
class ModelResult:
    model_id: str
    exposure_column: str
    exposure_organ: str
    exposure_modality: str
    exposure_endpoint: str
    mediator_column: str
    mediator_organ: str
    analysis_subset: str
    status: str
    message: str
    n: int
    deaths: int
    censored: int
    events_per_parameter_direct: float
    n_covariates: int
    a_beta: float = np.nan
    a_se_hc3: float = np.nan
    a_p_hc3: float = np.nan
    mediator_model_r2: float = np.nan
    b_log_hr: float = np.nan
    b_hr: float = np.nan
    b_se: float = np.nan
    b_p: float = np.nan
    direct_log_hr: float = np.nan
    direct_hr: float = np.nan
    direct_se: float = np.nan
    direct_p: float = np.nan
    total_log_hr: float = np.nan
    total_hr: float = np.nan
    total_se: float = np.nan
    total_p: float = np.nan
    indirect_log_hr: float = np.nan
    indirect_hr_like: float = np.nan
    indirect_delta_se: float = np.nan
    indirect_delta_z: float = np.nan
    indirect_delta_p: float = np.nan
    boot_ci_low: float = np.nan
    boot_ci_high: float = np.nan
    boot_hr_like_ci_low: float = np.nan
    boot_hr_like_ci_high: float = np.nan
    boot_p_empirical: float = np.nan
    bootstrap_requested: int = 0
    bootstrap_successes: int = 0
    bootstrap_success_fraction: float = np.nan
    max_vif: float = np.nan
    ph_exposure_p: float = np.nan
    ph_mediator_p: float = np.nan
    cox_direct_concordance: float = np.nan
    cox_total_concordance: float = np.nan
    warnings: str = ""


def clean_message(value: object) -> str:
    if value is None:
        return ""
    text = str(value)
    text = re.sub(r"[\r\n\t]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def normalize_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(value).lower()).strip("_")


def safe_exp(x: float) -> float:
    if not np.isfinite(x):
        return np.nan
    return float(np.exp(np.clip(x, -700, 700)))


def read_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        gz = Path(str(path) + ".gz")
        if gz.exists():
            path = gz
        else:
            raise FileNotFoundError(f"File not found: {path} (or {gz})")
    low = "".join(path.suffixes).lower()
    if ".csv" in low:
        return pd.read_csv(path, low_memory=False)
    return pd.read_csv(path, sep="\t", low_memory=False)


def standardize_id(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "participant_id" in out.columns:
        source = "participant_id"
    elif "eid" in out.columns:
        source = "eid"
    else:
        raise ValueError("No participant_id or eid column found.")
    ids = pd.to_numeric(out[source], errors="coerce")
    out = out.loc[ids.notna()].copy()
    out["participant_id"] = ids.loc[ids.notna()].astype("int64").astype(str)
    if source != "participant_id":
        out = out.drop(columns=[source])
    return out


def parse_epoch_metadata(column_name: str) -> dict[str, str] | None:
    # Expected collector form: epoch__brain__proteomics__dementia__acceleration_z
    pieces = column_name.split("__")
    if len(pieces) != 5 or pieces[0] != "epoch" or pieces[4] != "acceleration_z":
        return None
    return {
        "column": column_name,
        "organ": pieces[1],
        "modality": pieces[2],
        "endpoint": pieces[3],
        "measure": pieces[4],
    }


def classify_epoch_columns(epoch: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for col in epoch.columns:
        meta = parse_epoch_metadata(col)
        if meta is not None:
            rows.append(meta)
    metadata = pd.DataFrame(rows)
    if metadata.empty:
        raise ValueError("No standardized epoch__...__acceleration_z columns found.")

    exposures = metadata[
        metadata["modality"].isin(["proteomics", "metabolomics"])
        & metadata["endpoint"].ne("mortality")
    ].copy()
    mediators = metadata[
        metadata["modality"].eq("mri")
        & metadata["endpoint"].eq("mortality")
    ].copy()
    exposures = exposures.sort_values(["organ", "modality", "endpoint"]).reset_index(drop=True)
    mediators = mediators.sort_values(["organ"]).reset_index(drop=True)
    if exposures.empty:
        raise ValueError("No molecular disease EPOCH exposures found.")
    if mediators.empty:
        raise ValueError("No MRI mortality EPOCH mediators found.")
    return exposures, mediators


def find_prediction_file(root: Path, organ: str) -> Path:
    direct_dir = root / f"{organ}_mri_mortality_clock"
    if not direct_dir.is_dir():
        target = normalize_name(direct_dir.name)
        candidates = [p for p in root.iterdir() if p.is_dir() and normalize_name(p.name) == target]
        if len(candidates) != 1:
            raise FileNotFoundError(
                f"Could not uniquely resolve MRI mortality clock directory for {organ}: {candidates}"
            )
        direct_dir = candidates[0]
    files = sorted(direct_dir.glob("*_clock_predictions.tsv"))
    if len(files) != 1:
        raise FileNotFoundError(f"Expected one *_clock_predictions.tsv in {direct_dir}; found {files}")
    return files[0]


def build_landmark_survival(prediction_file: Path, test_only: bool) -> pd.DataFrame:
    header = pd.read_csv(prediction_file, sep="\t", nrows=0)
    landmark_candidates = ["imaging_date", "sample_date"]
    landmark = next((c for c in landmark_candidates if c in header.columns), None)
    if landmark is None:
        raise ValueError(
            f"{prediction_file} has neither imaging_date nor sample_date."
        )

    required = ["participant_id", landmark, "death_date", "admin_censor_date"]
    missing = [c for c in required if c not in header.columns]
    if missing:
        raise ValueError(f"{prediction_file} missing required fields: {missing}")
    usecols = required + (["split"] if "split" in header.columns else [])
    dat = standardize_id(pd.read_csv(prediction_file, sep="\t", usecols=usecols, low_memory=False))

    subset_label = "test" if test_only else "full"
    if test_only:
        if "split" not in dat.columns:
            raise ValueError(
                f"--test-only requested but {prediction_file} has no split column."
            )
        dat = dat[dat["split"].astype(str).str.strip().str.lower().eq("test")].copy()

    dat["landmark_date"] = pd.to_datetime(dat[landmark], errors="coerce")
    dat["death_date"] = pd.to_datetime(dat["death_date"], errors="coerce")
    dat["admin_censor_date"] = pd.to_datetime(dat["admin_censor_date"], errors="coerce")
    dat = dat.dropna(subset=["participant_id", "landmark_date", "admin_censor_date"]).copy()

    death_valid = (
        dat["death_date"].notna()
        & (dat["death_date"] > dat["landmark_date"])
        & (dat["death_date"] <= dat["admin_censor_date"])
    )
    dat["mortality_event"] = death_valid.astype(int)
    dat["observed_end_date"] = dat["admin_censor_date"]
    dat.loc[death_valid, "observed_end_date"] = dat.loc[death_valid, "death_date"]
    dat["mortality_time_years"] = (
        dat["observed_end_date"] - dat["landmark_date"]
    ).dt.total_seconds() / (365.25 * 24 * 3600)
    dat = dat[
        np.isfinite(dat["mortality_time_years"])
        & dat["mortality_time_years"].gt(0)
    ].copy()

    if dat["participant_id"].duplicated().any():
        dat = dat.sort_values(["participant_id", "landmark_date"]).drop_duplicates(
            "participant_id", keep="first"
        )

    print(
        f"Mortality source: {prediction_file} | landmark={landmark} | subset={subset_label} "
        f"| N={len(dat)} | deaths={int(dat['mortality_event'].sum())}",
        flush=True,
    )
    return dat[["participant_id", "landmark_date", "mortality_time_years", "mortality_event"]]


def prepare_covariates(path: Path, n_pcs: int) -> tuple[pd.DataFrame, list[str], list[str]]:
    cov = standardize_id(pd.read_csv(path, low_memory=False))
    rename_map: dict[str, str] = {}
    continuous: list[str] = []
    categorical: list[str] = []

    for new, old in CONTINUOUS_COVARIATE_MAP.items():
        if old in cov.columns:
            rename_map[old] = new
            continuous.append(new)
        else:
            warnings.warn(f"Requested continuous covariate not found and will be omitted: {old}")
    for new, old in CATEGORICAL_COVARIATE_MAP.items():
        if old in cov.columns:
            rename_map[old] = new
            categorical.append(new)
        else:
            warnings.warn(f"Requested categorical covariate not found and will be omitted: {old}")

    cov = cov.rename(columns=rename_map)
    for new in continuous:
        cov[new] = pd.to_numeric(cov[new], errors="coerce")
    for new in categorical:
        text = cov[new].astype("string").str.strip()
        cov[new] = text.replace({"": pd.NA, "NA": pd.NA, "NaN": pd.NA, ".": pd.NA, "-9999": pd.NA})

    for pc in range(1, n_pcs + 1):
        source = f"genetic_principal_components_f22009_0_{pc}"
        if source in cov.columns:
            new = f"PC{pc}"
            cov[new] = pd.to_numeric(cov[source], errors="coerce")
            continuous.append(new)
        else:
            warnings.warn(f"Requested PC not found and will be omitted: {source}")

    keep = ["participant_id"] + continuous + categorical
    cov = cov[keep].copy()
    if cov["participant_id"].duplicated().any():
        warnings.warn("Duplicate participant IDs in covariate file; keeping first record.")
        cov = cov.drop_duplicates("participant_id", keep="first")
    return cov, continuous, categorical


def zscore(values: pd.Series) -> pd.Series:
    x = pd.to_numeric(values, errors="coerce")
    sd = x.std(ddof=1)
    if not np.isfinite(sd) or sd <= 0:
        return pd.Series(np.nan, index=x.index)
    return (x - x.mean()) / sd


def make_model_dataset(
    merged: pd.DataFrame,
    exposure_column: str,
    mediator_column: str,
    continuous_covariates: list[str],
    categorical_covariates: list[str],
) -> tuple[pd.DataFrame, list[str]]:
    columns = [
        "participant_id",
        exposure_column,
        mediator_column,
        "mortality_time_years",
        "mortality_event",
    ] + [c for c in continuous_covariates + categorical_covariates if c in merged.columns]
    data = merged[columns].copy()
    data = data.rename(columns={exposure_column: "exposure", mediator_column: "mediator"})

    numeric = ["exposure", "mediator", "mortality_time_years", "mortality_event"] + [
        c for c in continuous_covariates if c in data.columns
    ]
    for c in numeric:
        data[c] = pd.to_numeric(data[c], errors="coerce")
    for c in categorical_covariates:
        if c in data.columns:
            text = data[c].astype("string").str.strip()
            data[c] = text.replace({"": pd.NA, "NA": pd.NA, "NaN": pd.NA, ".": pd.NA, "-9999": pd.NA})

    required = ["exposure", "mediator", "mortality_time_years", "mortality_event"] + [
        c for c in continuous_covariates + categorical_covariates if c in data.columns
    ]
    data = data.dropna(subset=required).copy()
    data = data[np.isfinite(data["mortality_time_years"]) & data["mortality_time_years"].gt(0)].copy()
    data["exposure"] = zscore(data["exposure"])
    data["mediator"] = zscore(data["mediator"])
    data = data.dropna(subset=["exposure", "mediator"]).copy()

    present_cat = [c for c in categorical_covariates if c in data.columns]
    informative_cat = [c for c in present_cat if data[c].nunique(dropna=True) >= 2]
    if informative_cat:
        data = pd.get_dummies(data, columns=informative_cat, drop_first=True, dtype=float)
    for c in set(present_cat) - set(informative_cat):
        data = data.drop(columns=[c])

    excluded = {"participant_id", "exposure", "mediator", "mortality_time_years", "mortality_event", "landmark_date"}
    covariates = [c for c in data.columns if c not in excluded]
    valid_covariates = []
    for c in covariates:
        data[c] = pd.to_numeric(data[c], errors="coerce")
        vals = data[c].replace([np.inf, -np.inf], np.nan)
        if vals.notna().sum() > 1 and vals.std(ddof=1) > 0:
            valid_covariates.append(c)
    data = data[["participant_id", "exposure", "mediator", "mortality_time_years", "mortality_event"] + valid_covariates]
    data = data.replace([np.inf, -np.inf], np.nan).dropna().copy()
    return data, valid_covariates


def calculate_vif(data: pd.DataFrame, columns: list[str]) -> float:
    if not columns:
        return np.nan
    x = data[columns].astype(float)
    x = sm.add_constant(x, has_constant="add")
    values = []
    for i, col in enumerate(x.columns):
        if col == "const":
            continue
        try:
            v = float(variance_inflation_factor(x.values, i))
            if np.isfinite(v):
                values.append(v)
        except Exception:
            pass
    return float(max(values)) if values else np.nan


def fit_ols_and_cox(
    data: pd.DataFrame,
    covariates: list[str],
    cox_penalizer: float,
) -> tuple[object, CoxPHFitter, CoxPHFitter, pd.DataFrame, pd.DataFrame]:
    x_a = sm.add_constant(data[["exposure"] + covariates].astype(float), has_constant="add")
    med_model = sm.OLS(data["mediator"].astype(float), x_a).fit(cov_type="HC3")

    direct_cols = ["mortality_time_years", "mortality_event", "exposure", "mediator"] + covariates
    direct_data = data[direct_cols].astype(float)
    cph_direct = CoxPHFitter(penalizer=cox_penalizer)
    cph_direct.fit(
        direct_data,
        duration_col="mortality_time_years",
        event_col="mortality_event",
        show_progress=False,
    )

    total_cols = ["mortality_time_years", "mortality_event", "exposure"] + covariates
    total_data = data[total_cols].astype(float)
    cph_total = CoxPHFitter(penalizer=cox_penalizer)
    cph_total.fit(
        total_data,
        duration_col="mortality_time_years",
        event_col="mortality_event",
        show_progress=False,
    )
    return med_model, cph_direct, cph_total, direct_data, total_data


def indirect_delta(a: float, se_a: float, b: float, se_b: float) -> tuple[float, float, float, float]:
    indirect = a * b
    variance = (b * b * se_a * se_a) + (a * a * se_b * se_b)
    if not np.isfinite(variance) or variance <= 0:
        return indirect, np.nan, np.nan, np.nan
    se = math.sqrt(variance)
    z = indirect / se
    p = 2.0 * stats.norm.sf(abs(z))
    return indirect, se, z, float(p)


def empirical_two_sided_p(samples: np.ndarray) -> float:
    if samples.size == 0:
        return np.nan
    lower = (np.sum(samples <= 0) + 1.0) / (samples.size + 1.0)
    upper = (np.sum(samples >= 0) + 1.0) / (samples.size + 1.0)
    return float(min(1.0, 2.0 * min(lower, upper)))


def bootstrap_indirect(
    data: pd.DataFrame,
    covariates: list[str],
    n_bootstrap: int,
    seed: int,
    cox_penalizer: float,
) -> tuple[np.ndarray, list[str]]:
    if n_bootstrap <= 0:
        return np.array([], dtype=float), []
    rng = np.random.default_rng(seed)
    n = len(data)
    estimates: list[float] = []
    warning_messages: list[str] = []

    for _ in range(n_bootstrap):
        idx = rng.integers(0, n, size=n)
        boot = data.iloc[idx].reset_index(drop=True)
        # A bootstrap resample with zero deaths cannot fit a Cox model.
        if boot["mortality_event"].sum() == 0:
            continue
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                x_a = sm.add_constant(boot[["exposure"] + covariates].astype(float), has_constant="add")
                med = sm.OLS(boot["mediator"].astype(float), x_a).fit()
                direct_cols = ["mortality_time_years", "mortality_event", "exposure", "mediator"] + covariates
                cph = CoxPHFitter(penalizer=cox_penalizer)
                cph.fit(
                    boot[direct_cols].astype(float),
                    duration_col="mortality_time_years",
                    event_col="mortality_event",
                    show_progress=False,
                )
                a = float(med.params["exposure"])
                b = float(cph.params_["mediator"])
                value = a * b
                if np.isfinite(value):
                    estimates.append(value)
                if caught and len(warning_messages) < 20:
                    warning_messages.extend(clean_message(w.message) for w in caught[:2])
        except Exception as exc:
            if len(warning_messages) < 20:
                warning_messages.append(f"bootstrap fit failure: {clean_message(exc)}")
            continue
    return np.asarray(estimates, dtype=float), warning_messages


def ph_pvalues(cph: CoxPHFitter, direct_data: pd.DataFrame) -> tuple[float, float]:
    try:
        test = proportional_hazard_test(cph, direct_data, time_transform="rank")
        summary = test.summary
        px = float(summary.loc["exposure", "p"]) if "exposure" in summary.index else np.nan
        pm = float(summary.loc["mediator", "p"]) if "mediator" in summary.index else np.nan
        return px, pm
    except Exception:
        return np.nan, np.nan


def analyze_model(
    merged: pd.DataFrame,
    exposure_meta: pd.Series,
    mediator_meta: pd.Series,
    continuous_covariates: list[str],
    categorical_covariates: list[str],
    minimum_n: int,
    minimum_deaths: int,
    n_bootstrap: int,
    min_bootstrap_success_fraction: float,
    seed: int,
    cox_penalizer: float,
    subset_label: str,
    bootstrap_dir: Path | None,
) -> ModelResult:
    exposure_col = str(exposure_meta["column"])
    mediator_col = str(mediator_meta["column"])
    model_id = "__".join([
        str(exposure_meta["organ"]),
        str(exposure_meta["modality"]),
        str(exposure_meta["endpoint"]),
        "to",
        str(mediator_meta["organ"]),
        "mri_mortality",
    ])

    try:
        data, covariates = make_model_dataset(
            merged,
            exposure_column=exposure_col,
            mediator_column=mediator_col,
            continuous_covariates=continuous_covariates,
            categorical_covariates=categorical_covariates,
        )
    except Exception as exc:
        return ModelResult(
            model_id=model_id, exposure_column=exposure_col,
            exposure_organ=str(exposure_meta["organ"]), exposure_modality=str(exposure_meta["modality"]),
            exposure_endpoint=str(exposure_meta["endpoint"]), mediator_column=mediator_col,
            mediator_organ=str(mediator_meta["organ"]), analysis_subset=subset_label,
            status="data_error", message=clean_message(exc), n=0, deaths=0, censored=0,
            events_per_parameter_direct=np.nan, n_covariates=0, bootstrap_requested=n_bootstrap,
        )

    n = len(data)
    deaths = int(data["mortality_event"].sum())
    censored = n - deaths
    n_parameters_direct = 2 + len(covariates)
    epp = deaths / n_parameters_direct if n_parameters_direct > 0 else np.nan

    base = dict(
        model_id=model_id, exposure_column=exposure_col,
        exposure_organ=str(exposure_meta["organ"]), exposure_modality=str(exposure_meta["modality"]),
        exposure_endpoint=str(exposure_meta["endpoint"]), mediator_column=mediator_col,
        mediator_organ=str(mediator_meta["organ"]), analysis_subset=subset_label,
        n=n, deaths=deaths, censored=censored, events_per_parameter_direct=epp,
        n_covariates=len(covariates), bootstrap_requested=n_bootstrap,
    )

    if n < minimum_n or deaths < minimum_deaths:
        return ModelResult(
            **base,
            status="skipped_low_sample",
            message=f"N={n}, deaths={deaths}; thresholds are N>={minimum_n}, deaths>={minimum_deaths}",
        )

    warning_messages: list[str] = []
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            med_model, cph_direct, cph_total, direct_data, _ = fit_ols_and_cox(
                data, covariates, cox_penalizer=cox_penalizer
            )
            warning_messages.extend(clean_message(w.message) for w in caught)
    except Exception as exc:
        return ModelResult(
            **base,
            status="model_error",
            message=clean_message(exc),
            warnings=" | ".join(dict.fromkeys(warning_messages)),
        )

    try:
        a = float(med_model.params["exposure"])
        a_se = float(med_model.bse["exposure"])
        a_p = float(med_model.pvalues["exposure"])
        b = float(cph_direct.params_["mediator"])
        b_se = float(cph_direct.standard_errors_["mediator"])
        b_p = float(cph_direct.summary.loc["mediator", "p"])
        direct = float(cph_direct.params_["exposure"])
        direct_se = float(cph_direct.standard_errors_["exposure"])
        direct_p = float(cph_direct.summary.loc["exposure", "p"])
        total = float(cph_total.params_["exposure"])
        total_se = float(cph_total.standard_errors_["exposure"])
        total_p = float(cph_total.summary.loc["exposure", "p"])
        indirect, ind_se, ind_z, ind_p = indirect_delta(a, a_se, b, b_se)
    except Exception as exc:
        return ModelResult(
            **base,
            status="extraction_error",
            message=clean_message(exc),
            warnings=" | ".join(dict.fromkeys(warning_messages)),
        )

    max_vif = calculate_vif(data, ["exposure", "mediator"] + covariates)
    ph_x, ph_m = ph_pvalues(cph_direct, direct_data)

    boot, boot_warnings = bootstrap_indirect(
        data=data,
        covariates=covariates,
        n_bootstrap=n_bootstrap,
        seed=seed,
        cox_penalizer=cox_penalizer,
    )
    warning_messages.extend(boot_warnings)
    success_fraction = boot.size / n_bootstrap if n_bootstrap > 0 else np.nan
    enough_boot = n_bootstrap > 0 and success_fraction >= min_bootstrap_success_fraction and boot.size >= 100
    if enough_boot:
        ci_low, ci_high = np.quantile(boot, [0.025, 0.975])
        boot_p = empirical_two_sided_p(boot)
    else:
        ci_low = ci_high = boot_p = np.nan
        if n_bootstrap > 0:
            warning_messages.append(
                f"insufficient bootstrap success: {boot.size}/{n_bootstrap} "
                f"(< required fraction {min_bootstrap_success_fraction:.2f})"
            )

    if bootstrap_dir is not None:
        bootstrap_dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"indirect_log_hr": boot}).to_csv(
            bootstrap_dir / f"bootstrap__{normalize_name(model_id)}.tsv.gz",
            sep="\t", index=False, compression="gzip"
        )

    status = "ok"
    message = ""
    if warning_messages:
        message = "warnings captured; inspect warnings column"

    return ModelResult(
        **base,
        status=status,
        message=message,
        a_beta=a,
        a_se_hc3=a_se,
        a_p_hc3=a_p,
        mediator_model_r2=float(med_model.rsquared),
        b_log_hr=b,
        b_hr=safe_exp(b),
        b_se=b_se,
        b_p=b_p,
        direct_log_hr=direct,
        direct_hr=safe_exp(direct),
        direct_se=direct_se,
        direct_p=direct_p,
        total_log_hr=total,
        total_hr=safe_exp(total),
        total_se=total_se,
        total_p=total_p,
        indirect_log_hr=indirect,
        indirect_hr_like=safe_exp(indirect),
        indirect_delta_se=ind_se,
        indirect_delta_z=ind_z,
        indirect_delta_p=ind_p,
        boot_ci_low=float(ci_low),
        boot_ci_high=float(ci_high),
        boot_hr_like_ci_low=safe_exp(float(ci_low)),
        boot_hr_like_ci_high=safe_exp(float(ci_high)),
        boot_p_empirical=boot_p,
        bootstrap_successes=int(boot.size),
        bootstrap_success_fraction=float(success_fraction) if np.isfinite(success_fraction) else np.nan,
        max_vif=max_vif,
        ph_exposure_p=ph_x,
        ph_mediator_p=ph_m,
        cox_direct_concordance=float(cph_direct.concordance_index_),
        cox_total_concordance=float(cph_total.concordance_index_),
        warnings=" | ".join(dict.fromkeys(x for x in warning_messages if x)),
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run 45x7 EPOCH OLS+Cox statistical mediation with bootstrap uncertainty."
    )
    parser.add_argument("--epoch-wide", type=Path, default=DEFAULT_EPOCH_WIDE)
    parser.add_argument("--wholebody-root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--covariates", type=Path, default=DEFAULT_COVARIATES)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument(
        "--test-only", action="store_true",
        help="Restrict MRI mortality prediction files to split == test. Default is full population."
    )
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--min-bootstrap-success-fraction", type=float, default=0.80)
    parser.add_argument("--seed", type=int, default=20260806)
    parser.add_argument("--minimum-n", type=int, default=500)
    parser.add_argument("--minimum-deaths", type=int, default=20)
    parser.add_argument("--n-pcs", type=int, default=10)
    parser.add_argument(
        "--cox-penalizer", type=float, default=0.0,
        help="L2 penalizer for Cox models. Default 0.0 (unpenalized primary analysis)."
    )
    parser.add_argument(
        "--save-bootstrap-samples", action="store_true",
        help="Save per-model bootstrap a*b samples as compressed TSV files."
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.bootstrap < 0:
        raise ValueError("--bootstrap must be >= 0")
    if not 0 < args.min_bootstrap_success_fraction <= 1:
        raise ValueError("--min-bootstrap-success-fraction must be in (0, 1].")

    subset_label = "test" if args.test_only else "full"
    output_dir = args.output_dir
    if output_dir is None:
        output_dir = DEFAULT_OUTPUT_TEST if args.test_only else DEFAULT_OUTPUT_FULL
    output_dir.mkdir(parents=True, exist_ok=True)
    bootstrap_dir = output_dir / "bootstrap_samples" if args.save_bootstrap_samples else None

    print("=" * 80, flush=True)
    print("EPOCH OLS + Cox + bootstrap statistical mediation", flush=True)
    print(f"Analysis subset: {subset_label}", flush=True)
    print(f"Bootstrap replicates per fitted model: {args.bootstrap}", flush=True)
    print(f"Output directory: {output_dir}", flush=True)
    print("=" * 80, flush=True)

    epoch = standardize_id(read_table(args.epoch_wide))
    if epoch["participant_id"].duplicated().any():
        raise ValueError("EPOCH wide table contains duplicated participant_id values.")
    exposures, mediators = classify_epoch_columns(epoch)
    print(f"EPOCH table: {len(epoch)} participants x {epoch.shape[1]} columns", flush=True)
    print(f"Molecular disease exposures: {len(exposures)}", flush=True)
    print(f"MRI mortality mediators: {len(mediators)}", flush=True)

    model_grid = exposures.assign(_key=1).merge(mediators.assign(_key=1), on="_key", suffixes=("_x", "_m")).drop(columns="_key")
    n_requested = len(model_grid)
    print(f"Models requested: {n_requested}", flush=True)

    cov, continuous_covariates, categorical_covariates = prepare_covariates(args.covariates, args.n_pcs)
    print(f"Continuous covariates: {continuous_covariates}", flush=True)
    print(f"Categorical covariates: {categorical_covariates}", flush=True)

    mortality_by_mediator: dict[str, pd.DataFrame] = {}
    for _, med in mediators.iterrows():
        organ = str(med["organ"])
        prediction = find_prediction_file(args.wholebody_root, organ)
        mortality_by_mediator[str(med["column"])] = build_landmark_survival(prediction, args.test_only)

    results: list[dict] = []
    model_counter = 0
    for _, exp in exposures.iterrows():
        for _, med in mediators.iterrows():
            model_counter += 1
            exp_col = str(exp["column"])
            med_col = str(med["column"])
            print(
                f"[{model_counter}/{n_requested}] {exp['organ']} {exp['modality']} {exp['endpoint']} "
                f"-> {med['organ']} MRI mortality EPOCH -> post-MRI mortality",
                flush=True,
            )
            outcome = mortality_by_mediator[med_col]
            merged = (
                epoch[["participant_id", exp_col, med_col]]
                .merge(outcome, on="participant_id", how="inner")
                .merge(cov, on="participant_id", how="left")
            )
            result = analyze_model(
                merged=merged,
                exposure_meta=exp,
                mediator_meta=med,
                continuous_covariates=continuous_covariates,
                categorical_covariates=categorical_covariates,
                minimum_n=args.minimum_n,
                minimum_deaths=args.minimum_deaths,
                n_bootstrap=args.bootstrap,
                min_bootstrap_success_fraction=args.min_bootstrap_success_fraction,
                seed=args.seed + model_counter,
                cox_penalizer=args.cox_penalizer,
                subset_label=subset_label,
                bootstrap_dir=bootstrap_dir,
            )
            results.append(asdict(result))
            print(
                f"    status={result.status} N={result.n} deaths={result.deaths} "
                f"indirect_p={result.indirect_delta_p if np.isfinite(result.indirect_delta_p) else 'NA'} "
                f"bootstrap={result.bootstrap_successes}/{result.bootstrap_requested}",
                flush=True,
            )

    table = pd.DataFrame(results)

    # Primary family-wise correction. Use the number of ALL requested models, not
    # merely the models that happened to fit, so the declared family is explicit
    # and conservative.
    family_m = n_requested
    bonf_threshold = 0.05 / family_m if family_m > 0 else np.nan
    table["multiplicity_family_m"] = family_m
    table["bonferroni_alpha_threshold"] = bonf_threshold
    table["indirect_bonferroni_p"] = np.nan
    table["indirect_bonferroni_significant"] = False
    table["indirect_bh_fdr"] = np.nan

    valid = table["indirect_delta_p"].notna() & np.isfinite(table["indirect_delta_p"])
    if valid.any():
        raw = table.loc[valid, "indirect_delta_p"].astype(float).to_numpy()
        # Bonferroni uses the predeclared full family m=315 (or whatever the grid is).
        bonf = np.minimum(raw * family_m, 1.0)
        table.loc[valid, "indirect_bonferroni_p"] = bonf
        table.loc[valid, "indirect_bonferroni_significant"] = raw < bonf_threshold
        # BH-FDR is secondary and is calculated among estimable p-values.
        _, qvals, _, _ = multipletests(raw, alpha=0.05, method="fdr_bh")
        table.loc[valid, "indirect_bh_fdr"] = qvals

    table["bootstrap_ci_excludes_zero"] = (
        table["boot_ci_low"].notna()
        & table["boot_ci_high"].notna()
        & ((table["boot_ci_low"] > 0) | (table["boot_ci_high"] < 0))
    )
    table["ph_exposure_violation_p_lt_0_05"] = table["ph_exposure_p"].notna() & table["ph_exposure_p"].lt(0.05)
    table["ph_mediator_violation_p_lt_0_05"] = table["ph_mediator_p"].notna() & table["ph_mediator_p"].lt(0.05)
    table["vif_gt_5"] = table["max_vif"].notna() & table["max_vif"].gt(5)

    result_path = output_dir / "OLS_Cox_bootstrap_mediation_results.tsv"
    qc_path = output_dir / "OLS_Cox_bootstrap_model_QC.tsv"
    table.to_csv(result_path, sep="\t", index=False)

    qc_columns = [
        "model_id", "exposure_column", "exposure_organ", "exposure_modality", "exposure_endpoint",
        "mediator_column", "mediator_organ", "analysis_subset", "status", "message", "n", "deaths",
        "censored", "events_per_parameter_direct", "n_covariates", "bootstrap_requested",
        "bootstrap_successes", "bootstrap_success_fraction", "max_vif", "ph_exposure_p", "ph_mediator_p",
        "warnings",
    ]
    table[qc_columns].to_csv(qc_path, sep="\t", index=False)

    exposures.to_csv(output_dir / "exposure_metadata.tsv", sep="\t", index=False)
    mediators.to_csv(output_dir / "mediator_metadata.tsv", sep="\t", index=False)

    settings = {
        "analysis_subset": subset_label,
        "epoch_wide": str(args.epoch_wide),
        "wholebody_root": str(args.wholebody_root),
        "covariates": str(args.covariates),
        "output_dir": str(output_dir),
        "bootstrap": args.bootstrap,
        "minimum_n": args.minimum_n,
        "minimum_deaths": args.minimum_deaths,
        "n_pcs": args.n_pcs,
        "cox_penalizer": args.cox_penalizer,
        "n_requested_models": n_requested,
        "bonferroni_family_m": family_m,
        "bonferroni_alpha_threshold": bonf_threshold,
        "continuous_covariates": continuous_covariates,
        "categorical_covariates": categorical_covariates,
        "indirect_inference": (
            "a*b on Cox log-hazard scale; delta-method p-value is primary for Bonferroni; "
            "participant bootstrap percentile CI is robustness inference."
        ),
        "causal_warning": (
            "Statistical mediation only. This is not a formally identified causal natural indirect effect; "
            "Cox non-collapsibility means total need not equal direct plus indirect."
        ),
    }
    (output_dir / "analysis_settings.json").write_text(json.dumps(settings, indent=2) + "\n")

    status_counts = table["status"].value_counts(dropna=False).to_dict()
    n_ok = int((table["status"] == "ok").sum())
    n_bonf = int(table["indirect_bonferroni_significant"].fillna(False).sum())
    n_boot_ci = int(table["bootstrap_ci_excludes_zero"].fillna(False).sum())

    print("\nCompleted.", flush=True)
    print(f"Output directory: {output_dir}", flush=True)
    print(f"Models requested: {n_requested}", flush=True)
    print(f"Models fitted successfully: {n_ok}", flush=True)
    print(f"Status counts: {status_counts}", flush=True)
    print(f"Bonferroni threshold: {bonf_threshold:.8g} (0.05/{family_m})", flush=True)
    print(f"Bonferroni-significant indirect effects: {n_bonf}", flush=True)
    print(f"Bootstrap 95% CI excludes zero: {n_boot_ci}", flush=True)
    print(f"Main results: {result_path}", flush=True)
    print(f"QC results: {qc_path}", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("Interrupted by user.", file=sys.stderr)
        raise
    except Exception as exc:
        print(f"FATAL ERROR: {clean_message(exc)}", file=sys.stderr)
        traceback.print_exc()
        raise
