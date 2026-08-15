#!/usr/bin/env python3
"""
Individual Cox survival analyses for incident Alzheimer's disease onset (G309):
brain-proteomics mortality EPOCH versus each underlying brain-enriched protein.

The three organ feature tables are combined to recover the full brain-proteomics
sample:
    training/training_4589.tsv
    PT/patient_pop.tsv
    test/ind_test_500.tsv

One separate Cox model is fitted for each predictor:
    covariates + predictor_z

There is NO joint EPOCH + protein model.

Primary effect:
    hazard ratio (HR) per 1-SD higher EPOCH or protein.

Default covariates:
    age_at_baseline, sex, bmi_at_baseline, smoking_status_at_baseline

Time origin:
    proteomics/EPOCH sample_date.

Incident cases:
    sample_date -> G309 diagnosis date.

Non-cases:
    sample_date -> earlier of death or administrative censoring.

By default, all predictors are evaluated in a common complete-case sample so
the forest-plot HRs are based on identical participants.
"""

from __future__ import annotations

import argparse
import os
import warnings
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test
from lifelines.utils import concordance_index
from scipy.stats import chi2

warnings.filterwarnings("ignore")

DEFAULT_EPOCH_COL = "brain_proteomics_mortality_clock_acceleration_z"
DEFAULT_COVARIATES = [
    "age_at_baseline",
    "sex",
    "bmi_at_baseline",
    "smoking_status_at_baseline",
]
DEFAULT_ADMIN_CENSOR_DATE = "2022-11-30"

# Non-protein columns expected in organ feature tables.
PROTEIN_METADATA_COLUMNS = {
    "participant_id",
    "eid",
    "session_id",
    "diagnosis",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Compare brain-proteomics mortality EPOCH with each underlying "
            "brain-enriched protein for incident G309 onset."
        )
    )
    p.add_argument("--icd-tsv", required=True)
    p.add_argument("--epoch-tsv", required=True)
    p.add_argument(
        "--organ-tsv",
        required=True,
        help=(
            "Comma-separated full-sample organ feature TSVs, e.g. "
            "training_4589.tsv,patient_pop.tsv,ind_test_500.tsv"
        ),
    )
    p.add_argument("--output-tsv", required=True)
    p.add_argument("--epoch-col", default=DEFAULT_EPOCH_COL)
    p.add_argument(
        "--sample-mode",
        choices=["common", "predictor_specific"],
        default="common",
        help=(
            "common: all EPOCH/protein models use identical complete cases; "
            "predictor_specific: each predictor uses its own complete cases."
        ),
    )
    p.add_argument(
        "--split",
        choices=["all", "train", "validation", "test"],
        default="all",
        help="Optional filter on the EPOCH prediction split column.",
    )
    p.add_argument(
        "--session-id",
        default="ses-M0",
        help="Session retained from organ TSVs when session_id exists.",
    )
    p.add_argument(
        "--covariates",
        default=",".join(DEFAULT_COVARIATES),
        help="Comma-separated covariates taken from the EPOCH prediction file.",
    )
    p.add_argument(
        "--protein-exclude",
        default="",
        help="Additional comma-separated organ-TSV columns to exclude.",
    )
    p.add_argument(
        "--admin-censor-date",
        default=DEFAULT_ADMIN_CENSOR_DATE,
    )
    p.add_argument("--min-case", type=int, default=20)
    p.add_argument("--min-noncase", type=int, default=20)
    p.add_argument("--penalizer", type=float, default=0.0)
    return p.parse_args()


def split_csv(text: str) -> List[str]:
    if text is None or not str(text).strip():
        return []
    return [x.strip() for x in str(text).split(",") if x.strip()]


def normalize_participant_id(
    df: pd.DataFrame,
    col: str = "participant_id",
) -> pd.DataFrame:
    out = df.copy()
    if col not in out.columns and "eid" in out.columns:
        out = out.rename(columns={"eid": col})
    if col not in out.columns:
        raise ValueError(f"Missing participant identifier: {col}")
    out[col] = pd.to_numeric(out[col], errors="coerce").astype("Int64")
    out = out.loc[out[col].notna()].copy()
    return out


def ensure_unique(
    df: pd.DataFrame,
    cols: Sequence[str],
    label: str,
) -> None:
    dup = df.duplicated(list(cols), keep=False)
    if dup.any():
        example = df.loc[dup, list(cols)].head(1).to_dict("records")[0]
        raise ValueError(
            f"{label} has duplicated key(s) {list(cols)}. Example: {example}"
        )


def clean_date(series: pd.Series) -> pd.Series:
    x = series.copy()
    x = x.replace(
        [0, 0.0, "0", "0.0", "", "NA", "NaN", "nan", "None", "-1", -1],
        np.nan,
    )
    return pd.to_datetime(x, errors="coerce")


def encode_sex(series: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    n_original = int(series.notna().sum())
    if n_original > 0 and numeric.notna().sum() / n_original >= 0.9:
        return numeric.astype(float)

    text = series.astype(str).str.strip().str.lower()
    return text.map(
        {
            "female": 0.0,
            "f": 0.0,
            "male": 1.0,
            "m": 1.0,
        }
    ).astype(float)


def bh_fdr(values: Sequence[float]) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    out = np.full(p.shape, np.nan, dtype=float)
    valid = np.isfinite(p) & (p >= 0) & (p <= 1)
    pv = p[valid]
    if pv.size == 0:
        return out

    order = np.argsort(pv)
    ranked = pv[order]
    m = ranked.size
    adj = ranked * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.minimum(adj, 1.0)

    restored = np.empty_like(adj)
    restored[order] = adj
    out[valid] = restored
    return out


def bonferroni(values: Sequence[float]) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    out = np.full(p.shape, np.nan, dtype=float)
    valid = np.isfinite(p) & (p >= 0) & (p <= 1)
    n = int(valid.sum())
    if n > 0:
        out[valid] = np.minimum(p[valid] * n, 1.0)
    return out


def read_disease_file(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    required = ["participant_id", "case", "date"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Disease TSV missing columns: {missing}")

    df = normalize_participant_id(df[required].copy())
    ensure_unique(df, ["participant_id"], "Disease TSV")

    df["case"] = (
        pd.to_numeric(df["case"], errors="coerce").fillna(0) == 1
    ).astype(int)
    df["event_date"] = clean_date(df["date"])
    df.loc[df["case"] == 0, "event_date"] = pd.NaT
    return df.drop(columns=["date"])


def read_epoch_file(
    path: str,
    epoch_col: str,
    covariates: List[str],
    split: str,
    default_admin_censor_date: str,
) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)

    required = ["participant_id", "sample_date", epoch_col] + covariates
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"EPOCH TSV missing columns: {missing}")

    if split != "all":
        if "split" not in df.columns:
            raise ValueError("--split requested but EPOCH TSV has no split column.")
        df = df.loc[
            df["split"].astype(str).str.strip().str.lower() == split.lower()
        ].copy()

    keep = ["participant_id", "sample_date", epoch_col] + covariates
    for optional in [
        "death_date",
        "admin_censor_date",
        "split",
        "organ_source_file",
    ]:
        if optional in df.columns:
            keep.append(optional)
    keep = list(dict.fromkeys(keep))

    df = normalize_participant_id(df[keep].copy())
    ensure_unique(df, ["participant_id"], "EPOCH TSV")

    df["sample_date"] = clean_date(df["sample_date"])

    if "death_date" in df.columns:
        df["death_date"] = clean_date(df["death_date"])
    else:
        df["death_date"] = pd.NaT

    if "admin_censor_date" in df.columns:
        df["admin_censor_date"] = clean_date(df["admin_censor_date"])
    else:
        df["admin_censor_date"] = pd.NaT

    df["admin_censor_date"] = df["admin_censor_date"].fillna(
        pd.Timestamp(default_admin_censor_date)
    )

    if "organ_source_file" in df.columns:
        df["organ_source_file"] = (
            df["organ_source_file"]
            .astype(str)
            .map(lambda x: os.path.basename(x.strip()) if x.strip() else x)
        )

    df[epoch_col] = pd.to_numeric(df[epoch_col], errors="coerce")

    for c in covariates:
        if c == "sex":
            df[c] = encode_sex(df[c])
        else:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    return df


def read_full_organ_sample(
    comma_paths: str,
    session_id: str,
    extra_exclude: List[str],
) -> Tuple[pd.DataFrame, List[str], List[str]]:
    paths = split_csv(comma_paths)
    if not paths:
        raise ValueError("--organ-tsv did not contain any paths.")

    for path in paths:
        if not os.path.isfile(path):
            raise FileNotFoundError(f"Organ feature TSV not found: {path}")

    frames: List[pd.DataFrame] = []
    protein_cols: Optional[List[str]] = None
    source_names: List[str] = []

    exclude = set(PROTEIN_METADATA_COLUMNS)
    exclude.update(extra_exclude)

    for index, path in enumerate(paths):
        source_name = os.path.basename(path)
        source_names.append(source_name)

        df = pd.read_csv(path, sep="\t", low_memory=False)
        df = normalize_participant_id(df)

        if session_id and "session_id" in df.columns:
            df = df.loc[
                df["session_id"].astype(str).str.strip() == session_id
            ].copy()

        if index == 0:
            protein_cols = [c for c in df.columns if c not in exclude]
            if not protein_cols:
                raise ValueError(
                    f"No protein columns detected in first organ TSV: {path}"
                )
        else:
            assert protein_cols is not None
            missing = [c for c in protein_cols if c not in df.columns]
            if missing:
                raise ValueError(
                    f"Organ TSV {path} is missing protein columns defined by "
                    f"the first file: {missing}"
                )

        assert protein_cols is not None
        keep = ["participant_id"] + protein_cols
        df = df[keep].copy()

        for c in protein_cols:
            df[c] = pd.to_numeric(df[c], errors="coerce")

        df["_protein_source_file"] = source_name
        frames.append(df)

    assert protein_cols is not None

    full = pd.concat(frames, axis=0, ignore_index=True, sort=False)

    # Duplicates across source files are allowed because source-file matching
    # can disambiguate them using EPOCH organ_source_file. Duplicates within
    # the same source are not allowed.
    ensure_unique(
        full,
        ["participant_id", "_protein_source_file"],
        "Combined organ TSVs",
    )

    usable = []
    for c in protein_cols:
        x = pd.to_numeric(full[c], errors="coerce")
        if x.notna().sum() > 1 and x.nunique(dropna=True) > 1:
            usable.append(c)

    if not usable:
        raise ValueError("No usable protein predictors after combining organ TSVs.")

    return (
        full[["participant_id", "_protein_source_file"] + usable].copy(),
        usable,
        source_names,
    )


def merge_epoch_and_proteins(
    epoch: pd.DataFrame,
    proteins: pd.DataFrame,
) -> Tuple[pd.DataFrame, str]:
    # Preferred: match by participant and exact source file because the EPOCH
    # prediction file records organ_source_file.
    if "organ_source_file" in epoch.columns:
        merged = epoch.merge(
            proteins,
            left_on=["participant_id", "organ_source_file"],
            right_on=["participant_id", "_protein_source_file"],
            how="inner",
            validate="one_to_one",
        )
        method = "participant_id_plus_organ_source_file"
        return merged, method

    # Fallback: participant-only merge is valid only if participant IDs are
    # globally unique across the concatenated organ tables.
    ensure_unique(
        proteins,
        ["participant_id"],
        "Combined organ TSVs without EPOCH source-file matching",
    )
    merged = epoch.merge(
        proteins,
        on="participant_id",
        how="inner",
        validate="one_to_one",
    )
    return merged, "participant_id_only"


def build_survival_table(
    disease: pd.DataFrame,
    epoch: pd.DataFrame,
    proteins: pd.DataFrame,
) -> Tuple[pd.DataFrame, Dict[str, object]]:
    epoch_protein, merge_method = merge_epoch_and_proteins(epoch, proteins)

    data = disease.merge(
        epoch_protein,
        on="participant_id",
        how="inner",
        validate="one_to_one",
    )

    qc: Dict[str, object] = {
        "protein_merge_method": merge_method,
        "N_epoch_protein_overlap": int(len(epoch_protein)),
        "N_after_disease_epoch_protein_merge": int(len(data)),
    }

    if "_protein_source_file" in data.columns:
        counts = data["_protein_source_file"].value_counts(dropna=False)
        qc["source_file_counts_before_survival_qc"] = ";".join(
            f"{k}:{int(v)}" for k, v in counts.items()
        )

    death = data["death_date"]
    admin = data["admin_censor_date"]

    data["censor_date"] = admin
    earlier_death = death.notna() & admin.notna() & (death < admin)
    data.loc[earlier_death, "censor_date"] = death.loc[earlier_death]
    death_only = death.notna() & admin.isna()
    data.loc[death_only, "censor_date"] = death.loc[death_only]

    case_missing_date = (data["case"] == 1) & data["event_date"].isna()
    prevalent_case = (
        (data["case"] == 1)
        & data["event_date"].notna()
        & data["sample_date"].notna()
        & (data["event_date"] <= data["sample_date"])
    )
    case_after_censor = (
        (data["case"] == 1)
        & data["event_date"].notna()
        & data["censor_date"].notna()
        & (data["event_date"] > data["censor_date"])
    )

    qc["N_excluded_case_missing_event_date"] = int(case_missing_date.sum())
    qc["N_excluded_prevalent_case"] = int(prevalent_case.sum())
    qc["N_excluded_case_after_censor"] = int(case_after_censor.sum())

    invalid_case = case_missing_date | prevalent_case | case_after_censor
    data = data.loc[~invalid_case].copy()

    data["end_date_survival"] = data["censor_date"]
    data.loc[
        data["case"] == 1,
        "end_date_survival",
    ] = data.loc[data["case"] == 1, "event_date"]

    data["time_days"] = (
        data["end_date_survival"] - data["sample_date"]
    ).dt.days

    invalid_followup = (
        data["sample_date"].isna()
        | data["end_date_survival"].isna()
        | data["time_days"].isna()
        | (data["time_days"] <= 0)
    )

    qc["N_excluded_invalid_followup"] = int(invalid_followup.sum())
    data = data.loc[~invalid_followup].copy()
    data["time_years"] = data["time_days"] / 365.25

    qc["N_valid_incident_survival_before_model_missingness"] = int(len(data))

    if "_protein_source_file" in data.columns:
        counts = data["_protein_source_file"].value_counts(dropna=False)
        qc["source_file_counts_after_survival_qc"] = ";".join(
            f"{k}:{int(v)}" for k, v in counts.items()
        )

    return data, qc


def usable_covariates(
    df: pd.DataFrame,
    requested: List[str],
) -> List[str]:
    out = []
    for c in requested:
        if c not in df.columns:
            continue
        x = pd.to_numeric(df[c], errors="coerce")
        if x.notna().sum() > 1 and x.nunique(dropna=True) > 1:
            out.append(c)
    return out


def standardize(
    series: pd.Series,
) -> Tuple[pd.Series, float, float]:
    x = pd.to_numeric(series, errors="coerce")
    mean = float(x.mean())
    sd = float(x.std(ddof=1))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("Predictor has zero or invalid SD.")
    return (x - mean) / sd, mean, sd


def fit_cox_with_retry(
    fit_df: pd.DataFrame,
    duration_col: str,
    event_col: str,
    initial_penalizer: float,
) -> Tuple[CoxPHFitter, float]:
    candidates = []
    for pen in [initial_penalizer, 0.001, 0.01, 0.1]:
        if pen not in candidates:
            candidates.append(pen)

    last_error = None
    for pen in candidates:
        try:
            cph = CoxPHFitter(penalizer=pen)
            cph.fit(
                fit_df,
                duration_col=duration_col,
                event_col=event_col,
                show_progress=False,
            )
            return cph, pen
        except Exception as exc:
            last_error = exc

    raise last_error


def get_cindex(
    cph: CoxPHFitter,
    fit_df: pd.DataFrame,
    duration_col: str,
    event_col: str,
) -> float:
    risk = cph.predict_partial_hazard(fit_df).values.ravel()
    return float(
        concordance_index(
            fit_df[duration_col],
            -risk,
            fit_df[event_col],
        )
    )


def lrt_pvalue(
    full: CoxPHFitter,
    reduced: CoxPHFitter,
) -> float:
    stat = 2.0 * (
        float(full.log_likelihood_) - float(reduced.log_likelihood_)
    )
    if not np.isfinite(stat) or stat < 0:
        return np.nan
    return float(chi2.sf(stat, 1))


def ph_test_pvalue(
    cph: CoxPHFitter,
    fit_df: pd.DataFrame,
    predictor: str,
) -> float:
    try:
        result = proportional_hazard_test(
            cph,
            fit_df,
            time_transform="rank",
        )
        return float(result.summary.loc[predictor, "p"])
    except Exception:
        return np.nan


def empty_result(
    disease_id: str,
    predictor: str,
    predictor_type: str,
    status: str,
    error: str = "",
) -> Dict[str, object]:
    return {
        "disease_id": disease_id,
        "predictor_type": predictor_type,
        "predictor": predictor,
        "status": status,
        "error": error,
    }


def analyze_predictor(
    data: pd.DataFrame,
    disease_id: str,
    predictor: str,
    predictor_type: str,
    covariates: List[str],
    args: argparse.Namespace,
    common_sample: bool,
) -> Dict[str, object]:
    needed = ["time_days", "time_years", "case", predictor] + covariates
    df = data[needed].copy()

    for c in needed:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    if not common_sample:
        df = df.replace([np.inf, -np.inf], np.nan).dropna().copy()

    df = df.loc[df["time_days"] > 0].copy()

    n_case = int(df["case"].sum())
    n_noncase = int((df["case"] == 0).sum())

    if n_case < args.min_case or n_noncase < args.min_noncase:
        out = empty_result(
            disease_id,
            predictor,
            predictor_type,
            "insufficient_events",
        )
        out.update(
            {
                "N": len(df),
                "N_case": n_case,
                "N_noncase": n_noncase,
            }
        )
        return out

    try:
        predictor_z, raw_mean, raw_sd = standardize(df[predictor])
    except Exception as exc:
        return empty_result(
            disease_id,
            predictor,
            predictor_type,
            "standardization_failed",
            str(exc),
        )

    df["predictor_z"] = predictor_z
    covars = usable_covariates(df, covariates)

    fit_cols = ["time_days", "case"] + covars + ["predictor_z"]
    fit_df = (
        df[fit_cols]
        .replace([np.inf, -np.inf], np.nan)
        .dropna()
        .copy()
    )

    n_case = int(fit_df["case"].sum())
    n_noncase = int((fit_df["case"] == 0).sum())

    if n_case < args.min_case or n_noncase < args.min_noncase:
        out = empty_result(
            disease_id,
            predictor,
            predictor_type,
            "insufficient_events_after_model_missingness",
        )
        out.update(
            {
                "N": len(fit_df),
                "N_case": n_case,
                "N_noncase": n_noncase,
            }
        )
        return out

    base_fit_df = fit_df[["time_days", "case"] + covars].copy()

    try:
        cph_base, pen_base = fit_cox_with_retry(
            base_fit_df,
            "time_days",
            "case",
            args.penalizer,
        )
        cph_pred, pen_pred = fit_cox_with_retry(
            fit_df,
            "time_days",
            "case",
            args.penalizer,
        )
    except Exception as exc:
        out = empty_result(
            disease_id,
            predictor,
            predictor_type,
            "cox_fit_failed",
            str(exc),
        )
        out.update(
            {
                "N": len(fit_df),
                "N_case": n_case,
                "N_noncase": n_noncase,
            }
        )
        return out

    beta = float(cph_pred.params_.loc["predictor_z"])
    se = float(cph_pred.standard_errors_.loc["predictor_z"])
    z = beta / se if np.isfinite(se) and se > 0 else np.nan
    hr = float(np.exp(beta))
    ci_lo = float(np.exp(beta - 1.96 * se))
    ci_hi = float(np.exp(beta + 1.96 * se))
    p_value = float(cph_pred.summary.loc["predictor_z", "p"])

    base_c = get_cindex(cph_base, base_fit_df, "time_days", "case")
    pred_c = get_cindex(cph_pred, fit_df, "time_days", "case")

    followup = df["time_years"]
    event_followup = df.loc[df["case"] == 1, "time_years"]

    out = empty_result(disease_id, predictor, predictor_type, "ok")
    out.update(
        {
            "N": int(len(fit_df)),
            "N_case": n_case,
            "N_noncase": n_noncase,
            "event_rate": float(n_case / len(fit_df)),
            "predictor_raw_mean": raw_mean,
            "predictor_raw_sd": raw_sd,
            "beta_per_1SD": beta,
            "se": se,
            "z": z,
            "hr_per_1SD": hr,
            "ci95_lower": ci_lo,
            "ci95_upper": ci_hi,
            "p_value": p_value,
            "base_cindex": base_c,
            "predictor_cindex": pred_c,
            "delta_cindex_vs_base": pred_c - base_c,
            "lrt_p_vs_base": lrt_pvalue(cph_pred, cph_base),
            "ph_test_p_value": ph_test_pvalue(
                cph_pred,
                fit_df,
                "predictor_z",
            ),
            "followup_years_min": float(followup.min()),
            "followup_years_median": float(followup.median()),
            "followup_years_max": float(followup.max()),
            "event_followup_years_min": (
                float(event_followup.min()) if len(event_followup) else np.nan
            ),
            "event_followup_years_median": (
                float(event_followup.median()) if len(event_followup) else np.nan
            ),
            "event_followup_years_max": (
                float(event_followup.max()) if len(event_followup) else np.nan
            ),
            "covariates": ",".join(covars),
            "penalizer_base": pen_base,
            "penalizer_predictor": pen_pred,
        }
    )
    return out


def add_multiple_testing(results: pd.DataFrame) -> pd.DataFrame:
    out = results.copy()

    if "p_value" not in out.columns:
        out["p_value"] = np.nan

    p = pd.to_numeric(out["p_value"], errors="coerce").to_numpy(float)
    out["p_fdr_bh_all_predictors"] = bh_fdr(p)
    out["p_bonferroni_all_predictors"] = bonferroni(p)

    out["p_fdr_bh_proteins_only"] = np.nan
    out["p_bonferroni_proteins_only"] = np.nan

    mask = out["predictor_type"].eq("protein")
    protein_p = pd.to_numeric(
        out.loc[mask, "p_value"],
        errors="coerce",
    ).to_numpy(float)

    out.loc[mask, "p_fdr_bh_proteins_only"] = bh_fdr(protein_p)
    out.loc[mask, "p_bonferroni_proteins_only"] = bonferroni(protein_p)
    return out


def main() -> None:
    args = parse_args()

    covariates = split_csv(args.covariates)
    protein_exclude = split_csv(args.protein_exclude)

    disease_id = Path(args.icd_tsv).name
    for suffix in [
        "_diagnosis_clock_disease_free.tsv",
        "_diagnosis_clock.tsv",
        ".tsv",
    ]:
        if disease_id.endswith(suffix):
            disease_id = disease_id[: -len(suffix)]
            break

    print("=" * 84)
    print("Brain-proteomics EPOCH vs underlying proteins: G309 incident survival")
    print("=" * 84)
    print(f"Disease endpoint: {disease_id}")
    print(f"EPOCH column:     {args.epoch_col}")
    print(f"Sample mode:      {args.sample_mode}")
    print(f"EPOCH split:      {args.split}")
    print(f"Covariates:       {', '.join(covariates)}")
    print()

    disease = read_disease_file(args.icd_tsv)

    epoch = read_epoch_file(
        path=args.epoch_tsv,
        epoch_col=args.epoch_col,
        covariates=covariates,
        split=args.split,
        default_admin_censor_date=args.admin_censor_date,
    )

    proteins, protein_cols, source_names = read_full_organ_sample(
        comma_paths=args.organ_tsv,
        session_id=args.session_id,
        extra_exclude=protein_exclude,
    )

    print("Full organ feature sample:")
    for source in source_names:
        n = int((proteins["_protein_source_file"] == source).sum())
        print(f"  {source}: {n:,}")
    print(f"Detected protein predictors: {len(protein_cols)}")
    print("  " + ", ".join(protein_cols))
    print()

    data, survival_qc = build_survival_table(
        disease=disease,
        epoch=epoch,
        proteins=proteins,
    )

    print("Survival-data QC:")
    for key, value in survival_qc.items():
        print(f"  {key}: {value}")
    print()

    predictors = [(args.epoch_col, "EPOCH")] + [
        (protein, "protein") for protein in protein_cols
    ]

    if args.sample_mode == "common":
        common_needed = (
            ["participant_id", "time_days", "time_years", "case"]
            + covariates
            + [p for p, _ in predictors]
        )

        missing = [c for c in common_needed if c not in data.columns]
        if missing:
            raise ValueError(
                "Missing columns for common-sample analysis: "
                + ", ".join(missing)
            )

        common = data[common_needed].copy()

        for c in (
            ["time_days", "time_years", "case"]
            + covariates
            + [p for p, _ in predictors]
        ):
            common[c] = pd.to_numeric(common[c], errors="coerce")

        before = len(common)
        common = (
            common.replace([np.inf, -np.inf], np.nan)
            .dropna()
            .copy()
        )
        common = common.loc[common["time_days"] > 0].copy()

        print("Common complete-case sample:")
        print(f"  Before model missingness: {before:,}")
        print(f"  Final N:                  {len(common):,}")
        print(f"  G309 cases:               {int(common['case'].sum()):,}")
        print(f"  Non-cases:                {int((common['case'] == 0).sum()):,}")
        print()

        analysis_data = common
        common_sample = True
    else:
        analysis_data = data
        common_sample = False

    rows = []

    for i, (predictor, predictor_type) in enumerate(predictors, start=1):
        print(
            f"[{i:02d}/{len(predictors):02d}] "
            f"{predictor_type}: {predictor}"
        )

        result = analyze_predictor(
            data=analysis_data,
            disease_id=disease_id,
            predictor=predictor,
            predictor_type=predictor_type,
            covariates=covariates,
            args=args,
            common_sample=common_sample,
        )

        result["sample_mode"] = args.sample_mode
        result["epoch_split"] = args.split
        result["time_origin"] = "proteomics_sample_date"
        result["epoch_column"] = args.epoch_col
        result["N_underlying_proteins"] = len(protein_cols)
        result["organ_tsv_files"] = ",".join(source_names)

        for key, value in survival_qc.items():
            result[key] = value

        rows.append(result)

        if result["status"] == "ok":
            print(
                f"    HR={result['hr_per_1SD']:.4f} "
                f"[{result['ci95_lower']:.4f}, {result['ci95_upper']:.4f}] "
                f"P={result['p_value']:.3g} "
                f"delta-C={result['delta_cindex_vs_base']:+.4f}"
            )
        else:
            print(
                f"    status={result['status']} "
                f"error={result.get('error', '')}"
            )

    results = pd.DataFrame(rows)
    results = add_multiple_testing(results)

    type_order = results["predictor_type"].map(
        {"EPOCH": 0, "protein": 1}
    ).fillna(2)
    results["_type_order"] = type_order
    results = results.sort_values(
        ["_type_order", "p_value", "predictor"],
        ascending=[True, True, True],
        na_position="last",
        kind="mergesort",
    ).drop(columns="_type_order")

    preferred = [
        "disease_id",
        "predictor_type",
        "predictor",
        "status",
        "error",
        "sample_mode",
        "epoch_split",
        "time_origin",
        "N",
        "N_case",
        "N_noncase",
        "event_rate",
        "beta_per_1SD",
        "se",
        "z",
        "hr_per_1SD",
        "ci95_lower",
        "ci95_upper",
        "p_value",
        "p_fdr_bh_all_predictors",
        "p_bonferroni_all_predictors",
        "p_fdr_bh_proteins_only",
        "p_bonferroni_proteins_only",
        "base_cindex",
        "predictor_cindex",
        "delta_cindex_vs_base",
        "lrt_p_vs_base",
        "ph_test_p_value",
        "predictor_raw_mean",
        "predictor_raw_sd",
        "followup_years_min",
        "followup_years_median",
        "followup_years_max",
        "event_followup_years_min",
        "event_followup_years_median",
        "event_followup_years_max",
        "covariates",
        "penalizer_base",
        "penalizer_predictor",
        "N_underlying_proteins",
        "epoch_column",
        "organ_tsv_files",
        "protein_merge_method",
        "N_epoch_protein_overlap",
        "N_after_disease_epoch_protein_merge",
        "source_file_counts_before_survival_qc",
        "N_excluded_case_missing_event_date",
        "N_excluded_prevalent_case",
        "N_excluded_case_after_censor",
        "N_excluded_invalid_followup",
        "N_valid_incident_survival_before_model_missingness",
        "source_file_counts_after_survival_qc",
    ]

    ordered = [c for c in preferred if c in results.columns]
    remaining = [c for c in results.columns if c not in ordered]
    results = results[ordered + remaining]

    output = Path(args.output_tsv)
    output.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(
        output,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    print()
    print("=" * 84)
    print("Finished")
    print("=" * 84)
    print(f"Output: {output}")
    print(f"Rows:   {len(results)} (1 EPOCH + {len(protein_cols)} proteins)")
    print(f"Successful Cox models: {int((results['status'] == 'ok').sum())}")

    epoch_row = results.loc[
        (results["predictor_type"] == "EPOCH")
        & (results["status"] == "ok")
    ]
    if not epoch_row.empty:
        row = epoch_row.iloc[0]
        print()
        print("EPOCH:")
        print(
            f"  HR per SD = {row['hr_per_1SD']:.4f} "
            f"[{row['ci95_lower']:.4f}, {row['ci95_upper']:.4f}]"
        )
        print(f"  P = {row['p_value']:.3g}")
        print(f"  C-index = {row['predictor_cindex']:.4f}")


if __name__ == "__main__":
    main()