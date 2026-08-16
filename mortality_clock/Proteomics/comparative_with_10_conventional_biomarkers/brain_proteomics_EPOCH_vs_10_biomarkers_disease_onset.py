#!/usr/bin/env python3
"""
Apple-to-apple incident-disease survival comparison:

    Brain proteomics-based mortality EPOCH
    versus
    10 conventional UK Biobank biomarkers

Primary endpoints are supplied as disease-free TSV files containing:
    participant_id, case, date
where case indicates incident disease and date is the disease-onset date for cases.

The script is intended for the two analyses requested here:
    G309 = Alzheimer's disease (G30.9)
    I500 = heart failure (I50.0)

Key design features
-------------------
1. Use the full brain-proteomics mortality-EPOCH prediction file by default (train + validation + test) to maximize the disease-onset sample size.
2. Use the EPOCH acceleration z-score, not the raw mortality risk score.
3. Define time zero as the baseline sample/assessment date used by the EPOCH.
4. Exclude disease onset on/before baseline.
5. Censor noncases at the earlier of death or administrative censoring.
6. Build ONE strict common sample per disease: EPOCH + all 10 biomarkers + valid
   disease follow-up. Every individual predictor is evaluated in that same sample.
7. Default baseline covariates are age, sex, genetic ethnic grouping, assessment
   centre, smoking status, and BMI. Blood pressure is not included because systolic
   BP is itself one of the 10 benchmark biomarkers.
8. Standardize EPOCH and each biomarker within the common sample. HRs are per 1 SD.
9. Fit, in addition to individual predictor models, a 10-biomarker panel and a
   10-biomarker-panel + EPOCH model. This directly tests whether EPOCH adds
   prognostic information beyond the complete conventional-biomarker panel.
10. Direct C-index differences are evaluated with participant-level paired bootstrap
    using fixed fitted-model risk scores.

Conventional biomarkers (same definitions as the prior comparison script)
---------------------------------------------------------------------------
1. Overall health rating          UKB field 2178, baseline array 0
2. Left grip strength             field 46, baseline array 0
3. Right grip strength            field 47, baseline array 0
4. Systolic blood pressure        mean of field 4080 arrays 0 and 1
5. Peak expiratory flow           maximum of field 3064 arrays 0,1,2
6. C-reactive protein             log1p(field 30710, array 0)
7. Albumin                        field 30600, array 0
8. Cystatin C                     field 30720, array 0
9. HbA1c                          field 30750, array 0
10. Red-cell distribution width   field 30070, array 0

Main outputs
------------
*_individual_predictor_summary.tsv
*_model_panel_summary.tsv
*_epoch_beyond_10_biomarker_panel.tsv
*_paired_delta_cindex.tsv
*_common_sample_ids.tsv
*_qc.tsv
*_baseline_model_coefficients.tsv
*_biomarker_panel_coefficients.tsv
*_biomarker_panel_plus_epoch_coefficients.tsv
"""

from __future__ import print_function

import argparse
import os
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

try:
    from lifelines import CoxPHFitter
    from lifelines.statistics import proportional_hazard_test
    from lifelines.utils import concordance_index
except ImportError as exc:
    raise ImportError(
        "Missing lifelines. Activate the survival environment used for the existing analyses."
    ) from exc

try:
    from scipy.stats import chi2
except ImportError as exc:
    raise ImportError("Missing scipy.") from exc


# =============================================================================
# General helpers
# =============================================================================

def info(msg):
    print(msg, flush=True)


def warn(msg):
    print("WARNING: {}".format(msg), file=sys.stderr, flush=True)


def clean_id(x):
    if pd.isna(x):
        return None
    s = str(x).strip()
    if s == "":
        return None
    if re.match(r"^\d+\.0$", s):
        s = s[:-2]
    return s


def clean_col_name(x):
    s = str(x).strip().lower()
    s = re.sub(r"[^a-z0-9]+", "_", s)
    s = re.sub(r"_+", "_", s)
    return s.strip("_")


def split_csv(x):
    if x is None or str(x).strip() == "":
        return []
    return [z.strip() for z in str(x).split(",") if z.strip()]


def parse_date_series(s):
    out = pd.to_datetime(s, errors="coerce")
    numeric = pd.to_numeric(s, errors="coerce")
    if numeric.notna().sum() > out.notna().sum():
        out2 = pd.to_datetime(
            numeric, unit="D", origin="1899-12-30", errors="coerce"
        )
        if out2.notna().sum() > out.notna().sum():
            out = out2
    return out


def bh_fdr(p_values):
    p = np.asarray(p_values, dtype=float)
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


def bonferroni(p_values):
    p = np.asarray(p_values, dtype=float)
    out = np.full(p.shape, np.nan, dtype=float)
    valid = np.isfinite(p) & (p >= 0) & (p <= 1)
    n = int(valid.sum())
    if n:
        out[valid] = np.minimum(p[valid] * n, 1.0)
    return out


def choose_column_by_id_overlap(
    df,
    target_ids,
    min_overlap=20,
    prefer_patterns=None,
    exclude_cols=None,
):
    if exclude_cols is None:
        exclude_cols = set()
    best_col = None
    best_overlap = -1
    best_bonus = -1
    for col in df.columns:
        if col in exclude_cols:
            continue
        try:
            vals = set(v for v in df[col].map(clean_id).dropna().unique() if v is not None)
        except Exception:
            continue
        if not vals:
            continue
        overlap = len(vals.intersection(target_ids))
        bonus = 0
        if prefer_patterns:
            low = str(col).lower()
            for pat in prefer_patterns:
                if re.search(pat, low):
                    bonus += 1
        if overlap > best_overlap or (overlap == best_overlap and bonus > best_bonus):
            best_col = col
            best_overlap = overlap
            best_bonus = bonus
    if best_overlap < min_overlap:
        return None, best_overlap
    return best_col, best_overlap


def make_earliest_date(df, cols):
    parsed = [parse_date_series(df[c]) for c in cols]
    return pd.concat(parsed, axis=1).min(axis=1)


def save_cox_coefficients(cph, path):
    out = cph.summary.copy()
    out.insert(0, "term", out.index.astype(str))
    out = out.reset_index(drop=True)
    out.to_csv(path, sep="\t", index=False)


# =============================================================================
# Brain proteomics mortality EPOCH
# =============================================================================

def resolve_epoch_column(columns, override=None):
    cols = list(columns)
    if override:
        if override not in cols:
            raise ValueError("Requested --epoch-col not found: {}".format(override))
        return override

    preferred = [
        "brain_proteomics_mortality_clock_acceleration_z",
        "Brain_proteomics_mortality_clock_acceleration_z",
        "brain_proteomics_mortality_epoch_acceleration_z",
    ]
    for c in preferred:
        if c in cols:
            return c

    hits = []
    for c in cols:
        low = clean_col_name(c)
        if (
            "brain" in low
            and "proteom" in low
            and "mortality" in low
            and "acceleration" in low
            and (low.endswith("_z") or "acceleration_z" in low)
            and "year" not in low
        ):
            hits.append(c)

    if len(hits) == 1:
        return hits[0]
    if len(hits) > 1:
        raise ValueError(
            "Multiple candidate Brain proteomics mortality EPOCH acceleration columns found: {}. "
            "Pass --epoch-col explicitly.".format(hits)
        )

    candidates = [
        c for c in cols
        if "brain" in clean_col_name(c)
        and "proteom" in clean_col_name(c)
        and "mortality" in clean_col_name(c)
    ]
    raise ValueError(
        "Could not resolve the Brain proteomics mortality EPOCH acceleration z-score. "
        "Pass --epoch-col explicitly. Brain/proteomics/mortality candidate columns: {}".format(
            candidates[:30]
        )
    )


def read_epoch_predictions(path, epoch_col_override=None, analysis_split="all"):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(str(path))

    header = list(pd.read_csv(path, sep="\t", nrows=0).columns)
    if "participant_id" not in header:
        raise ValueError("participant_id is missing from EPOCH predictions file.")
    epoch_col = resolve_epoch_column(header, epoch_col_override)

    optional = [
        "split",
        "baseline_date",
        "sample_date",
        "death_date",
        "admin_censor_date",
        "age_at_baseline",
    ]
    usecols = ["participant_id", epoch_col] + [c for c in optional if c in header]
    df = pd.read_csv(path, sep="\t", usecols=usecols, dtype={"participant_id": "str"})
    df["participant_id"] = df["participant_id"].map(clean_id)
    df = df.dropna(subset=["participant_id"]).drop_duplicates("participant_id").copy()

    if analysis_split == "test" and "split" in df.columns:
        before = len(df)
        df = df.loc[df["split"].astype(str).str.lower() == "test"].copy()
        info("Filtered EPOCH predictions to held-out test split: {:,} -> {:,}".format(before, len(df)))
    elif analysis_split == "test" and "split" not in df.columns:
        raise ValueError(
            "--analysis-split test was requested, but the EPOCH prediction file has no 'split' column. "
            "Use brain_proteomics_mortality_clock_predictions.tsv, which contains the split column, "
            "or run with --analysis-split all."
        )
    else:
        info("Using FULL available EPOCH prediction sample (train + validation + test): N={:,}".format(len(df)))

    df = df.rename(columns={epoch_col: "epoch_raw"})
    df["epoch_raw"] = pd.to_numeric(df["epoch_raw"], errors="coerce")

    for c in ["baseline_date", "sample_date", "death_date", "admin_censor_date"]:
        if c in df.columns:
            df[c] = parse_date_series(df[c])

    # EPOCH measurement date: prefer sample_date, then baseline_date.
    df["epoch_baseline_date"] = pd.NaT
    if "sample_date" in df.columns:
        df["epoch_baseline_date"] = df["sample_date"]
    if "baseline_date" in df.columns:
        df["epoch_baseline_date"] = df["epoch_baseline_date"].where(
            df["epoch_baseline_date"].notna(), df["baseline_date"]
        )

    keep = ["participant_id", "epoch_raw", "epoch_baseline_date"]
    if "death_date" in df.columns:
        df = df.rename(columns={"death_date": "epoch_death_date"})
        keep.append("epoch_death_date")
    if "admin_censor_date" in df.columns:
        df = df.rename(columns={"admin_censor_date": "epoch_admin_censor_date"})
        keep.append("epoch_admin_censor_date")
    if "age_at_baseline" in df.columns:
        keep.append("age_at_baseline")
    if "split" in df.columns:
        keep.append("split")

    return df[keep].copy(), epoch_col


# =============================================================================
# ID linkage and baseline/death data
# =============================================================================

def detect_field53_col(df, instance=0, user_col=None):
    if user_col:
        if user_col not in df.columns:
            raise ValueError("Requested Field 53 column not found: {}".format(user_col))
        return user_col

    exact = [
        "53-{}.0".format(instance),
        "53_{}_0".format(instance),
        "date_of_attending_assessment_centre_f53_{}_0".format(instance),
    ]
    for c in exact:
        if c in df.columns:
            return c

    candidates = []
    patterns = [
        r"(^|[^0-9])53([^0-9]+){}([^0-9]+)0([^0-9]|$)".format(instance),
        r"f[._-]?53[._-]?{}[._-]?0".format(instance),
        r"date.*attend.*{}.*0".format(instance),
        r"assessment.*date.*{}.*0".format(instance),
    ]
    for col in df.columns:
        low = str(col).lower()
        score = sum(25 for pat in patterns if re.search(pat, low))
        if "53" in low:
            score += 5
        if "date" in low:
            score += 4
        if score == 0:
            continue
        parsed = parse_date_series(df[col])
        n_dates = int(parsed.notna().sum())
        if n_dates > 100:
            score += 20
        elif n_dates > 10:
            score += 10
        candidates.append((score, n_dates, col))
    if not candidates:
        raise ValueError("Could not auto-detect baseline Field 53 instance 0_0.")
    candidates.sort(key=lambda x: (x[0], x[1]), reverse=True)
    return candidates[0][2]


def detect_death_date_cols(df, user_col=None):
    if user_col:
        if user_col not in df.columns:
            raise ValueError("Requested death-date column not found: {}".format(user_col))
        return [user_col]
    candidates = []
    for col in df.columns:
        low = str(col).lower()
        score = 0
        if "40000" in low:
            score += 30
        if "death" in low:
            score += 15
        if "date" in low:
            score += 5
        if score == 0:
            continue
        parsed = parse_date_series(df[col])
        n_dates = int(parsed.notna().sum())
        if n_dates > 0:
            score += min(20, n_dates // 100)
        candidates.append((score, n_dates, col))
    candidates.sort(key=lambda x: (x[0], x[1]), reverse=True)
    out = [x[2] for x in candidates if x[1] > 0]
    if not out:
        raise ValueError("Could not auto-detect a death-date column.")
    return out


def detect_score_id_col_in_match(id_match_df, score_ids, override=None):
    if override:
        if override not in id_match_df.columns:
            raise ValueError("Requested ID-match score-ID column not found: {}".format(override))
        return override
    col, overlap = choose_column_by_id_overlap(
        id_match_df,
        score_ids,
        min_overlap=20,
        prefer_patterns=[r"penn", r"eid", r"ukb", r"participant", r"id"],
    )
    if col is None:
        raise ValueError("Could not detect the score-ID column in ID-match CSV.")
    info("ID-match score-ID column: {} overlap={}".format(col, overlap))
    return col


def map_death_to_score_ids(
    death_df,
    id_match_df,
    score_ids,
    death_id_col_arg=None,
    idmatch_score_col_arg=None,
    idmatch_death_col_arg=None,
):
    death_df = death_df.copy()
    id_match_df = id_match_df.copy()

    if death_id_col_arg:
        if death_id_col_arg not in death_df.columns:
            raise ValueError("death_id_col not found: {}".format(death_id_col_arg))
        death_id_col = death_id_col_arg
        direct_overlap = len(set(death_df[death_id_col].map(clean_id).dropna()).intersection(score_ids))
    else:
        death_id_col, direct_overlap = choose_column_by_id_overlap(
            death_df, score_ids, min_overlap=20, prefer_patterns=[r"eid", r"id", r"participant"]
        )

    if death_id_col is not None and direct_overlap >= 20:
        info("Death file uses score IDs directly: {} overlap={}".format(death_id_col, direct_overlap))
        death_df["participant_id"] = death_df[death_id_col].map(clean_id)
        return death_df

    info("Death file does not directly match score IDs. Using ID-match CSV.")
    idmatch_score_col = detect_score_id_col_in_match(id_match_df, score_ids, idmatch_score_col_arg)

    if idmatch_death_col_arg:
        if idmatch_death_col_arg not in id_match_df.columns:
            raise ValueError("idmatch_death_col not found: {}".format(idmatch_death_col_arg))
        idmatch_death_col = idmatch_death_col_arg
    else:
        idmatch_death_col = None

    if death_id_col_arg:
        death_id_col = death_id_col_arg
    else:
        death_id_col = None

    if death_id_col is None or idmatch_death_col is None:
        best = None
        for dcol in death_df.columns:
            dvals = set(death_df[dcol].map(clean_id).dropna().unique())
            if not dvals:
                continue
            for mcol in id_match_df.columns:
                if mcol == idmatch_score_col:
                    continue
                mvals = set(id_match_df[mcol].map(clean_id).dropna().unique())
                overlap = len(dvals.intersection(mvals))
                if best is None or overlap > best[0]:
                    best = (overlap, dcol, mcol)
        if best is None or best[0] < 20:
            raise ValueError("Could not detect death-ID mapping through ID-match CSV.")
        if death_id_col is None:
            death_id_col = best[1]
        if idmatch_death_col is None:
            idmatch_death_col = best[2]

    map_df = id_match_df[[idmatch_score_col, idmatch_death_col]].copy()
    map_df["participant_id"] = map_df[idmatch_score_col].map(clean_id)
    map_df["_death_merge_id"] = map_df[idmatch_death_col].map(clean_id)
    map_df = map_df[["participant_id", "_death_merge_id"]].dropna().drop_duplicates()

    death_df["_death_merge_id"] = death_df[death_id_col].map(clean_id)
    out = death_df.merge(map_df, on="_death_merge_id", how="left")
    info("Death rows mapped to score IDs: {:,}".format(out["participant_id"].notna().sum()))
    return out


def map_generic_table_to_score_ids(
    df,
    id_match_df,
    score_ids,
    table_label,
    table_id_col_arg=None,
    idmatch_score_col_arg=None,
):
    out = df.copy()
    if table_id_col_arg:
        if table_id_col_arg not in out.columns:
            raise ValueError("{} ID column not found: {}".format(table_label, table_id_col_arg))
        source_id_col = table_id_col_arg
        direct_overlap = len(set(out[source_id_col].map(clean_id).dropna()).intersection(score_ids))
    else:
        source_id_col, direct_overlap = choose_column_by_id_overlap(
            out,
            score_ids,
            min_overlap=20,
            prefer_patterns=[r"^eid$", r"participant", r"ukb", r"id"],
        )

    if source_id_col is not None and direct_overlap >= 20:
        info("{} uses score IDs directly: {} overlap={}".format(table_label, source_id_col, direct_overlap))
        out["participant_id"] = out[source_id_col].map(clean_id)
        return out, source_id_col, "direct"

    info("{} does not directly match score IDs. Using ID-match CSV.".format(table_label))
    idmatch_score_col = detect_score_id_col_in_match(id_match_df, score_ids, idmatch_score_col_arg)

    preferred = [c for c in out.columns if re.search(r"(^eid$|participant|ukb|id)", str(c).lower())]
    if source_id_col is not None:
        source_candidates = [(source_id_col, set(out[source_id_col].map(clean_id).dropna().unique()))]
    else:
        if not preferred:
            preferred = list(out.columns)
        source_candidates = []
        for c in preferred:
            vals = set(out[c].map(clean_id).dropna().unique())
            if len(vals) >= 20:
                source_candidates.append((c, vals))

    best = None
    for source_col, source_vals in source_candidates:
        for match_col in id_match_df.columns:
            if match_col == idmatch_score_col:
                continue
            match_vals = set(id_match_df[match_col].map(clean_id).dropna().unique())
            overlap = len(source_vals.intersection(match_vals))
            if best is None or overlap > best[0]:
                best = (overlap, source_col, match_col)
    if best is None or best[0] < 20:
        raise ValueError("Could not link {} to score IDs through ID-match CSV.".format(table_label))

    overlap, source_id_col, match_source_col = best
    info("{} source ID column: {}".format(table_label, source_id_col))
    info("ID-match column corresponding to {}: {} overlap={}".format(table_label, match_source_col, overlap))

    map_df = id_match_df[[idmatch_score_col, match_source_col]].copy()
    map_df["participant_id"] = map_df[idmatch_score_col].map(clean_id)
    map_df["_source_merge_id"] = map_df[match_source_col].map(clean_id)
    map_df = map_df[["participant_id", "_source_merge_id"]].dropna().drop_duplicates()

    # If the source table already has a participant_id column but it is not in
    # the EPOCH/score ID system, preserve it under a temporary provenance name
    # so the mapped score-system participant_id can be created unambiguously.
    if "participant_id" in out.columns:
        out = out.rename(columns={"participant_id": "_source_participant_id_original"})
        source_id_col_for_merge = (
            "_source_participant_id_original" if source_id_col == "participant_id" else source_id_col
        )
    else:
        source_id_col_for_merge = source_id_col

    out["_source_merge_id"] = out[source_id_col_for_merge].map(clean_id)
    out = out.merge(map_df, on="_source_merge_id", how="left")
    return out, source_id_col, "id_match_csv"


# =============================================================================
# Covariates
# =============================================================================

def detect_first_matching_col(df, patterns):
    for pat in patterns:
        for col in df.columns:
            if re.search(pat, str(col).lower()):
                return col
    return None


def detect_age_col(cov_df, instance=0):
    preferred = "age_when_attended_assessment_centre_f21003_{}_0".format(instance)
    if preferred in cov_df.columns:
        return preferred
    for patterns in [
        [r"21003.*{}.*0".format(instance), r"age.*{}.*0".format(instance)],
        [r"age_at_recruitment", r"age_when_attended_assessment", r"^age$", r"age"],
    ]:
        hits = []
        for col in cov_df.columns:
            low = str(col).lower()
            if not any(re.search(pat, low) for pat in patterns):
                continue
            vals = pd.to_numeric(cov_df[col], errors="coerce")
            if vals.notna().sum() > 100 and vals.between(20, 100).mean() > 0.50:
                hits.append((vals.notna().sum(), col))
        if hits:
            hits.sort(reverse=True)
            return hits[0][1]
    return None


def select_comparison_covariates(cov_df, instance=0, covariate_cols_arg=None):
    if covariate_cols_arg:
        cols = split_csv(covariate_cols_arg)
        missing = [x for x in cols if x not in cov_df.columns]
        if missing:
            raise ValueError("Requested covariates not found: {}".format(missing))
        return cols, "User-specified exact comparison covariates: {}".format(";".join(cols))

    selected = []
    source_records = []

    age_col = detect_age_col(cov_df, instance)
    if age_col:
        selected.append(age_col)
        source_records.append("Age={}".format(age_col))

    candidates = [
        ("Sex", ["sex_f31_0_0", "genetic_sex_f22001_0_0"], [r"^sex$", r"genetic_sex", r"reported_sex"]),
        ("Ethnicity", ["genetic_ethnic_grouping_f22006_0_0"], [r"genetic_ethnic_grouping", r"ethnic", r"race", r"21000"]),
        ("Assessment_center", ["uk_biobank_assessment_centre_f54_{}_0".format(instance)], [r"54.*{}.*0".format(instance), r"assessment.*center", r"assessment.*centre"]),
        ("Smoking", ["smoking_status_f20116_{}_0".format(instance)], [r"smoking_status.*{}.*0".format(instance), r"20116.*{}.*0".format(instance), r"smoking_status"]),
        ("BMI", ["body_mass_index_bmi_f23104_{}_0".format(instance), "body_mass_index_bmi_f21001_{}_0".format(instance)], [r"body_mass_index.*{}.*0".format(instance), r"bmi.*{}.*0".format(instance), r"body_mass_index", r"bmi"]),
    ]

    for label, exacts, patterns in candidates:
        col = next((x for x in exacts if x in cov_df.columns), None)
        if col is None:
            col = detect_first_matching_col(cov_df, patterns)
        if col:
            selected.append(col)
            source_records.append("{}={}".format(label, col))
        else:
            warn("{} covariate was not auto-detected.".format(label))

    out = []
    for c in selected:
        if c in cov_df.columns and c not in out:
            out.append(c)
    return out, "; ".join(source_records)


def force_categorical_covariate(col):
    low = str(col).lower()
    return any(
        p in low
        for p in [
            "sex",
            "genetic_ethnic_grouping",
            "ethnic_background",
            "assessment_centre",
            "assessment_center",
            "smoking_status",
        ]
    )


def build_covariate_design(df, duration_col, event_col, covariate_cols):
    d = df[[duration_col, event_col] + list(covariate_cols)].copy()
    d = d.replace([np.inf, -np.inf], np.nan)
    y = d[[duration_col, event_col]].reset_index(drop=True)
    x = d[list(covariate_cols)].copy()
    processed = []

    for col in x.columns:
        s = x[col]
        if force_categorical_covariate(col):
            cat = s.astype("object").where(s.notna(), "Missing").astype(str)
            if cat.nunique(dropna=False) <= 1:
                continue
            if cat.nunique(dropna=False) > 80:
                warn("Skipping high-cardinality categorical covariate {}.".format(col))
                continue
            processed.append(pd.get_dummies(cat, prefix=clean_col_name(col), drop_first=True).astype(float))
            continue

        numeric = pd.to_numeric(s, errors="coerce")
        if numeric.notna().mean() >= 0.90 and numeric.nunique(dropna=True) > 5:
            med = numeric.median()
            if pd.isna(med):
                continue
            processed.append(pd.DataFrame({clean_col_name(col): numeric.fillna(med).astype(float)}, index=x.index))
        else:
            cat = s.astype("object").where(s.notna(), "Missing").astype(str)
            if cat.nunique(dropna=False) <= 1:
                continue
            if cat.nunique(dropna=False) > 80:
                warn("Skipping high-cardinality categorical covariate {}.".format(col))
                continue
            processed.append(pd.get_dummies(cat, prefix=clean_col_name(col), drop_first=True).astype(float))

    if processed:
        xp = pd.concat(processed, axis=1)
        keep = [c for c in xp.columns if xp[c].nunique(dropna=True) > 1]
        xp = xp[keep].reset_index(drop=True)
    else:
        xp = pd.DataFrame(index=np.arange(len(y)))
    return pd.concat([y, xp], axis=1)


# =============================================================================
# Ten conventional biomarkers
# =============================================================================

def exact_or_prefix(df, exact, prefix=None):
    if exact in df.columns:
        return exact
    prefix = exact if prefix is None else prefix
    hits = [c for c in df.columns if str(c).startswith(prefix)]
    if not hits:
        return None
    hits = sorted(hits, key=lambda c: (0 if re.search(r"_0_0$", str(c)) else 1, str(c)))
    return hits[0]


def positive_numeric(s):
    x = pd.to_numeric(s, errors="coerce")
    x = x.where(np.isfinite(x), np.nan)
    return x.where(x > 0, np.nan)


def build_conventional_biomarkers(bio_df):
    df = bio_df.copy()
    source = {}

    c = exact_or_prefix(df, "overall_health_rating_f2178_0_0", "overall_health_rating_f2178_0_")
    if c is None:
        raise ValueError("Overall health rating field 2178 was not found.")
    health = pd.to_numeric(df[c], errors="coerce").where(lambda x: x.isin([1, 2, 3, 4]), np.nan)
    df["biomarker_overall_health_rating"] = health
    source["biomarker_overall_health_rating"] = c

    c = exact_or_prefix(df, "hand_grip_strength_left_f46_0_0", "hand_grip_strength_left_f46_0_")
    if c is None:
        raise ValueError("Left grip-strength field 46 was not found.")
    df["biomarker_grip_strength_left"] = positive_numeric(df[c])
    source["biomarker_grip_strength_left"] = c

    c = exact_or_prefix(df, "hand_grip_strength_right_f47_0_0", "hand_grip_strength_right_f47_0_")
    if c is None:
        raise ValueError("Right grip-strength field 47 was not found.")
    df["biomarker_grip_strength_right"] = positive_numeric(df[c])
    source["biomarker_grip_strength_right"] = c

    sbp_cols = [c for c in df.columns if re.match(r"^systolic_blood_pressure_automated_reading_f4080_0_[01]$", str(c))]
    if not sbp_cols:
        raise ValueError("Baseline automated systolic-BP field 4080 was not found.")
    sbp_mat = pd.concat([positive_numeric(df[c]) for c in sorted(sbp_cols)], axis=1)
    df["biomarker_systolic_bp"] = sbp_mat.mean(axis=1, skipna=True)
    source["biomarker_systolic_bp"] = ";".join(sorted(sbp_cols))

    pef_cols = [c for c in df.columns if re.match(r"^peak_expiratory_flow_pef_f3064_0_[012]$", str(c))]
    if not pef_cols:
        raise ValueError("Baseline peak-expiratory-flow field 3064 was not found.")
    pef_mat = pd.concat([positive_numeric(df[c]) for c in sorted(pef_cols)], axis=1)
    df["biomarker_peak_expiratory_flow"] = pef_mat.max(axis=1, skipna=True)
    source["biomarker_peak_expiratory_flow"] = ";".join(sorted(pef_cols))

    c = exact_or_prefix(df, "creactive_protein_f30710_0_0", "creactive_protein_f30710_0_")
    if c is None:
        raise ValueError("C-reactive-protein field 30710 was not found.")
    crp = pd.to_numeric(df[c], errors="coerce")
    crp = crp.where(np.isfinite(crp), np.nan).where(crp >= 0, np.nan)
    df["biomarker_log1p_crp"] = np.log1p(crp)
    source["biomarker_log1p_crp"] = "{} [log1p]".format(c)

    for out_col, exact, prefix, label in [
        ("biomarker_albumin", "albumin_f30600_0_0", "albumin_f30600_0_", "Albumin"),
        ("biomarker_cystatin_c", "cystatin_c_f30720_0_0", "cystatin_c_f30720_0_", "Cystatin C"),
        ("biomarker_hba1c", "glycated_haemoglobin_hba1c_f30750_0_0", "glycated_haemoglobin_hba1c_f30750_0_", "HbA1c"),
        ("biomarker_rdw", "red_blood_cell_erythrocyte_distribution_width_f30070_0_0", "red_blood_cell_erythrocyte_distribution_width_f30070_0_", "RDW"),
    ]:
        c = exact_or_prefix(df, exact, prefix)
        if c is None:
            raise ValueError("{} field was not found.".format(label))
        df[out_col] = positive_numeric(df[c])
        source[out_col] = c

    cols = [
        "biomarker_overall_health_rating",
        "biomarker_grip_strength_left",
        "biomarker_grip_strength_right",
        "biomarker_systolic_bp",
        "biomarker_peak_expiratory_flow",
        "biomarker_log1p_crp",
        "biomarker_albumin",
        "biomarker_cystatin_c",
        "biomarker_hba1c",
        "biomarker_rdw",
    ]
    names = {
        "biomarker_overall_health_rating": "Overall health rating",
        "biomarker_grip_strength_left": "Left grip strength",
        "biomarker_grip_strength_right": "Right grip strength",
        "biomarker_systolic_bp": "Systolic blood pressure",
        "biomarker_peak_expiratory_flow": "Peak expiratory flow",
        "biomarker_log1p_crp": "C-reactive protein (log1p)",
        "biomarker_albumin": "Albumin",
        "biomarker_cystatin_c": "Cystatin C",
        "biomarker_hba1c": "HbA1c",
        "biomarker_rdw": "Red-cell distribution width",
    }
    return df, cols, names, source


# =============================================================================
# Disease outcome
# =============================================================================

def parse_case_series(s):
    numeric = pd.to_numeric(s, errors="coerce")
    out = pd.Series(np.nan, index=s.index, dtype=float)
    out.loc[numeric == 0] = 0.0
    out.loc[numeric == 1] = 1.0
    text = s.astype(str).str.strip().str.lower()
    out.loc[text.isin(["true", "yes", "case", "event"])] = 1.0
    out.loc[text.isin(["false", "no", "control", "noncase", "non-case"])] = 0.0
    return out


def read_disease_file(path, disease_code, case_col="case", date_col="date"):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(str(path))
    header = list(pd.read_csv(path, sep="\t", nrows=0).columns)
    required = ["participant_id", case_col, date_col]
    missing = [c for c in required if c not in header]
    if missing:
        raise ValueError("Disease TSV missing required columns {}: {}".format(missing, path))
    keep = required + ([disease_code] if disease_code in header else [])
    d = pd.read_csv(path, sep="\t", usecols=keep, dtype={"participant_id": "str"}, low_memory=False)
    d["participant_id"] = d["participant_id"].map(clean_id)
    d["reported_case"] = parse_case_series(d[case_col])
    d["disease_date"] = parse_date_series(d[date_col])
    d = d.dropna(subset=["participant_id"]).drop_duplicates("participant_id", keep="first")
    return d[["participant_id", "reported_case", "disease_date"]]


# =============================================================================
# Cox helpers
# =============================================================================

def standardize_predictor(s):
    x = pd.to_numeric(s, errors="coerce")
    mean = float(x.mean())
    sd = float(x.std(ddof=1))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("Predictor has zero/invalid SD.")
    return (x - mean) / sd, mean, sd


def fit_cox_with_retry(design_df, duration_col, event_col, initial_penalizer=0.0):
    candidates = []
    for p in [initial_penalizer, 0.0, 0.001, 0.01, 0.05, 0.1]:
        if p not in candidates:
            candidates.append(p)
    last_error = None
    for penalizer in candidates:
        try:
            cph = CoxPHFitter(penalizer=penalizer)
            cph.fit(design_df, duration_col=duration_col, event_col=event_col, show_progress=False)
            risk = cph.predict_log_partial_hazard(design_df).values.reshape(-1)
            cindex = concordance_index(
                event_times=design_df[duration_col].values,
                predicted_scores=-risk,
                event_observed=design_df[event_col].values,
            )
            return cph, float(cindex), float(penalizer), risk
        except Exception as exc:
            last_error = exc
    raise last_error


def valid_lrt(full_cph, reduced_cph, full_penalizer, reduced_penalizer, df_diff=1):
    if float(full_penalizer) != 0.0 or float(reduced_penalizer) != 0.0:
        return np.nan, np.nan
    stat = 2.0 * (float(full_cph.log_likelihood_) - float(reduced_cph.log_likelihood_))
    if not np.isfinite(stat) or stat < 0:
        return np.nan, np.nan
    return float(stat), float(chi2.sf(stat, df_diff))


def predictor_ph_pvalue(cph, design_df, predictor_col="predictor_z"):
    try:
        test = proportional_hazard_test(cph, design_df, time_transform="rank")
        return float(test.summary.loc[predictor_col, "p"])
    except Exception:
        return np.nan


def paired_bootstrap_delta_cindex(time, event, risk_a, risk_b, label_a, label_b, n_boot, seed):
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=int)
    risk_a = np.asarray(risk_a, dtype=float)
    risk_b = np.asarray(risk_b, dtype=float)
    n = len(time)
    c_a = float(concordance_index(time, -risk_a, event))
    c_b = float(concordance_index(time, -risk_b, event))
    delta = c_a - c_b
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(int(n_boot)):
        idx = rng.integers(0, n, size=n)
        if np.sum(event[idx] == 1) < 2 or np.sum(event[idx] == 0) < 2:
            continue
        try:
            ca = float(concordance_index(time[idx], -risk_a[idx], event[idx]))
            cb = float(concordance_index(time[idx], -risk_b[idx], event[idx]))
            d = ca - cb
            if np.isfinite(d):
                boots.append(d)
        except Exception:
            continue
    boots = np.asarray(boots, dtype=float)
    if boots.size:
        lo, hi = np.quantile(boots, [0.025, 0.975])
        p_le0 = float(np.mean(boots <= 0.0))
        p_ge0 = float(np.mean(boots >= 0.0))
        p2 = float(min(1.0, 2.0 * min(p_le0, p_ge0)))
        # With finite resamples, report a non-zero lower bound later in manuscripts.
    else:
        lo = hi = p2 = np.nan
    return {
        "model_a": label_a,
        "model_b": label_b,
        "cindex_a": c_a,
        "cindex_b": c_b,
        "delta_cindex_a_minus_b": delta,
        "delta_cindex_ci_lower": float(lo) if np.isfinite(lo) else np.nan,
        "delta_cindex_ci_upper": float(hi) if np.isfinite(hi) else np.nan,
        "empirical_p_two_sided": p2,
        "n_bootstrap_requested": int(n_boot),
        "n_bootstrap_successful": int(boots.size),
        "bootstrap_note": "Participant bootstrap of C-index difference using fixed fitted-model risk scores; Cox models are not refit within each bootstrap.",
    }


# =============================================================================
# Main analysis
# =============================================================================

def run(args):
    outdir = Path(args.output_dir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = "brain_proteomics_mortality_EPOCH_vs_10_biomarkers_{}_{}_disease_onset".format(
        clean_col_name(args.disease_code), clean_col_name(args.analysis_split)
    )

    summary_out = outdir / "{}_individual_predictor_summary.tsv".format(prefix)
    model_out = outdir / "{}_model_panel_summary.tsv".format(prefix)
    incremental_out = outdir / "{}_epoch_beyond_10_biomarker_panel.tsv".format(prefix)
    pairwise_out = outdir / "{}_paired_delta_cindex.tsv".format(prefix)
    ids_out = outdir / "{}_common_sample_ids.tsv".format(prefix)
    qc_out = outdir / "{}_qc.tsv".format(prefix)
    base_coef_out = outdir / "{}_baseline_model_coefficients.tsv".format(prefix)
    panel_coef_out = outdir / "{}_biomarker_panel_coefficients.tsv".format(prefix)
    panel_epoch_coef_out = outdir / "{}_biomarker_panel_plus_epoch_coefficients.tsv".format(prefix)

    info("=" * 90)
    info("Brain proteomics mortality EPOCH vs 10 conventional biomarkers")
    info("Incident disease: {} ({})".format(args.disease_label, args.disease_code))
    info("Analysis split: {}".format(args.analysis_split))
    info("=" * 90)

    epoch_df, epoch_source_col = read_epoch_predictions(
        args.epoch_predictions,
        epoch_col_override=args.epoch_col,
        analysis_split=args.analysis_split,
    )
    score_ids = set(epoch_df["participant_id"].dropna().unique())
    info("EPOCH rows: {:,}; non-missing acceleration: {:,}".format(
        len(epoch_df), int(epoch_df["epoch_raw"].notna().sum())
    ))
    info("EPOCH source column: {}".format(epoch_source_col))

    id_match_df = pd.read_csv(args.id_match_csv, low_memory=False)

    # Disease file may already use Penn/score IDs; map only if needed.
    disease_raw = read_disease_file(
        args.disease_tsv, args.disease_code, case_col=args.case_col, date_col=args.date_col
    )
    disease_df, disease_id_col, disease_link_method = map_generic_table_to_score_ids(
        disease_raw,
        id_match_df,
        score_ids,
        table_label="Disease outcome file",
        table_id_col_arg="participant_id",
        idmatch_score_col_arg=args.idmatch_score_col,
    )
    disease_keep = disease_df[["participant_id", "reported_case", "disease_date"]].copy()
    disease_keep = disease_keep.dropna(subset=["participant_id"]).drop_duplicates("participant_id")

    # Death / baseline-date file.
    death_raw = pd.read_excel(args.death_xlsx, sheet_name=0, engine="openpyxl")
    death_df = map_death_to_score_ids(
        death_raw,
        id_match_df,
        score_ids,
        death_id_col_arg=args.death_id_col,
        idmatch_score_col_arg=args.idmatch_score_col,
        idmatch_death_col_arg=args.idmatch_death_col,
    )
    field53_col = detect_field53_col(death_df, 0, args.field53_0_col)
    death_cols = detect_death_date_cols(death_df, args.death_date_col)
    death_df["ukb_baseline_date"] = parse_date_series(death_df[field53_col])
    death_df["ukb_death_date"] = make_earliest_date(death_df, death_cols)
    death_keep = (
        death_df[["participant_id", "ukb_baseline_date", "ukb_death_date"]]
        .dropna(subset=["participant_id"])
        .groupby("participant_id", as_index=False)
        .agg({"ukb_baseline_date": "min", "ukb_death_date": "min"})
    )

    # Covariates.
    cov_raw = pd.read_csv(args.covariate_csv, low_memory=False)
    cov_df, cov_id_col, cov_link_method = map_generic_table_to_score_ids(
        cov_raw,
        id_match_df,
        score_ids,
        table_label="Covariate file",
        table_id_col_arg=args.covariate_id_col,
        idmatch_score_col_arg=args.idmatch_score_col,
    )
    covariate_cols, covariate_source_desc = select_comparison_covariates(
        cov_df, instance=0, covariate_cols_arg=args.covariate_cols
    )
    if not covariate_cols:
        raise ValueError("No baseline comparison covariates were selected.")
    cov_keep = cov_df[["participant_id"] + covariate_cols].dropna(subset=["participant_id"]).drop_duplicates("participant_id")
    info("Baseline comparison covariates: {}".format("; ".join(covariate_cols)))

    # 10 conventional biomarkers.
    bio_raw = pd.read_csv(args.biomarker_csv, low_memory=False)
    bio_df, bio_id_col, bio_link_method = map_generic_table_to_score_ids(
        bio_raw,
        id_match_df,
        score_ids,
        table_label="10-biomarker file",
        table_id_col_arg=args.biomarker_id_col,
        idmatch_score_col_arg=args.idmatch_score_col,
    )
    bio_df, biomarker_cols, display_names, biomarker_sources = build_conventional_biomarkers(bio_df)
    bio_keep = bio_df[["participant_id"] + biomarker_cols].dropna(subset=["participant_id"]).drop_duplicates("participant_id")

    # Merge all sources before defining the analysis sample.
    dat = epoch_df.merge(disease_keep, on="participant_id", how="left")
    dat = dat.merge(death_keep, on="participant_id", how="left")
    dat = dat.merge(cov_keep, on="participant_id", how="left")
    dat = dat.merge(bio_keep, on="participant_id", how="left")
    n_after_merge = int(len(dat))

    # Baseline date: EPOCH sample date is primary, Field 53 is fallback.
    dat["baseline_date"] = dat["epoch_baseline_date"].where(
        dat["epoch_baseline_date"].notna(), dat["ukb_baseline_date"]
    )

    admin_const = pd.Timestamp(args.admin_censor_date)
    dat["analysis_admin_censor_date"] = admin_const
    if "epoch_admin_censor_date" in dat.columns:
        dat["analysis_admin_censor_date"] = dat["epoch_admin_censor_date"].where(
            dat["epoch_admin_censor_date"].notna(), admin_const
        )

    # Death date: take earliest available death record.
    death_parts = [dat["ukb_death_date"]]
    if "epoch_death_date" in dat.columns:
        death_parts.append(dat["epoch_death_date"])
    dat["analysis_death_date"] = pd.concat(death_parts, axis=1).min(axis=1)

    dat = dat.loc[dat["baseline_date"].notna()].copy()
    dat = dat.loc[dat["baseline_date"] <= dat["analysis_admin_censor_date"]].copy()
    dat = dat.loc[
        dat["analysis_death_date"].isna()
        | (dat["analysis_death_date"] > dat["baseline_date"])
    ].copy()

    dat["censor_date"] = dat["analysis_admin_censor_date"]
    death_before_admin = (
        dat["analysis_death_date"].notna()
        & (dat["analysis_death_date"] < dat["analysis_admin_censor_date"])
    )
    dat.loc[death_before_admin, "censor_date"] = dat.loc[death_before_admin, "analysis_death_date"]

    # Require interpretable case status. A reported case must have an event date.
    n_missing_case_status = int(dat["reported_case"].isna().sum())
    dat = dat.loc[dat["reported_case"].notna()].copy()
    invalid_case_missing_date = (dat["reported_case"] == 1) & dat["disease_date"].isna()
    n_case_missing_date = int(invalid_case_missing_date.sum())
    dat = dat.loc[~invalid_case_missing_date].copy()

    prevalent = (
        (dat["reported_case"] == 1)
        & dat["disease_date"].notna()
        & (dat["disease_date"] <= dat["baseline_date"])
    )
    n_prevalent = int(prevalent.sum())
    dat = dat.loc[~prevalent].copy()

    dat["event"] = (
        (dat["reported_case"] == 1)
        & dat["disease_date"].notna()
        & (dat["disease_date"] > dat["baseline_date"])
        & (dat["disease_date"] <= dat["censor_date"])
    ).astype(int)

    dat["end_date"] = dat["censor_date"]
    dat.loc[dat["event"] == 1, "end_date"] = dat.loc[dat["event"] == 1, "disease_date"]
    dat["followup_time_years"] = (dat["end_date"] - dat["baseline_date"]).dt.days / 365.25
    dat = dat.loc[dat["followup_time_years"].notna() & (dat["followup_time_years"] > 0)].copy()
    n_valid_followup = int(len(dat))

    # Strict apple-to-apple sample: EPOCH + every biomarker observed.
    predictor_cols = ["epoch_raw"] + biomarker_cols
    common_mask = np.ones(len(dat), dtype=bool)
    for c in predictor_cols:
        x = pd.to_numeric(dat[c], errors="coerce")
        common_mask &= x.notna().to_numpy() & np.isfinite(x).to_numpy()
    common = dat.loc[common_mask].copy()

    if common.empty:
        raise ValueError("Common EPOCH + 10-biomarker analysis sample is empty.")
    n_common = int(len(common))
    n_cases = int(common["event"].sum())
    n_noncases = int(n_common - n_cases)
    if n_cases < args.min_events:
        warn("Only {} incident events in the common sample.".format(n_cases))

    info("=" * 90)
    info("COMMON APPLE-TO-APPLE DISEASE-ONSET SAMPLE")
    info("Rows after source merge: {:,}".format(n_after_merge))
    info("Rows with valid disease follow-up: {:,}".format(n_valid_followup))
    info("Rows complete for EPOCH + all 10 biomarkers: {:,}".format(n_common))
    info("Incident cases: {:,}; noncases: {:,}".format(n_cases, n_noncases))
    info("Median follow-up: {:.3f} years".format(float(common["followup_time_years"].median())))
    info("=" * 90)

    common[[
        "participant_id", "baseline_date", "disease_date", "analysis_death_date",
        "censor_date", "end_date", "followup_time_years", "event"
    ]].to_csv(ids_out, sep="\t", index=False)

    duration_col = "followup_time_years"
    event_col = "event"

    # Baseline covariate model.
    base_design = build_covariate_design(common, duration_col, event_col, covariate_cols)
    if len(base_design) != n_common:
        raise ValueError("Baseline design unexpectedly changed common-sample N.")
    base_cph, base_cindex, base_pen, base_risk = fit_cox_with_retry(
        base_design, duration_col, event_col, initial_penalizer=args.penalizer
    )
    save_cox_coefficients(base_cph, base_coef_out)

    predictor_specs = [{
        "column": "epoch_raw",
        "display": "Brain proteomics mortality EPOCH",
        "type": "EPOCH",
        "source": epoch_source_col,
        "transformation": "EPOCH acceleration z-score; re-standardized within common disease sample",
    }]
    transformations = {
        "biomarker_overall_health_rating": "ordinal 1-4; standardized within common sample",
        "biomarker_grip_strength_left": "raw baseline value; standardized within common sample",
        "biomarker_grip_strength_right": "raw baseline value; standardized within common sample",
        "biomarker_systolic_bp": "mean baseline automated readings; standardized within common sample",
        "biomarker_peak_expiratory_flow": "maximum baseline trial; standardized within common sample",
        "biomarker_log1p_crp": "log1p baseline CRP; standardized within common sample",
        "biomarker_albumin": "raw baseline value; standardized within common sample",
        "biomarker_cystatin_c": "raw baseline value; standardized within common sample",
        "biomarker_hba1c": "raw baseline value; standardized within common sample",
        "biomarker_rdw": "raw baseline value; standardized within common sample",
    }
    for c in biomarker_cols:
        predictor_specs.append({
            "column": c,
            "display": display_names[c],
            "type": "Conventional biomarker",
            "source": biomarker_sources[c],
            "transformation": transformations[c],
        })

    rows = []
    risks = {"Baseline covariates": base_risk}
    fitted_models = {}
    standardized = {}

    for i, spec in enumerate(predictor_specs, start=1):
        col = spec["column"]
        z, raw_mean, raw_sd = standardize_predictor(common[col])
        standardized[col] = np.asarray(z, dtype=float)
        design = base_design.copy()
        design["predictor_z"] = standardized[col]
        cph, cindex, used_pen, risk = fit_cox_with_retry(
            design, duration_col, event_col, initial_penalizer=args.penalizer
        )
        fitted_models[col] = cph
        risks[spec["display"]] = risk
        row = cph.summary.loc["predictor_z"]
        coef = float(row["coef"])
        se = float(row["se(coef)"])
        hr = float(np.exp(coef))
        ci_lo = float(np.exp(float(row["coef lower 95%"])))
        ci_hi = float(np.exp(float(row["coef upper 95%"])))
        p = float(row["p"])
        lrt_chi, lrt_p = valid_lrt(cph, base_cph, used_pen, base_pen, 1)
        result = {
            "disease_code": args.disease_code,
            "disease_label": args.disease_label,
            "analysis_split": args.analysis_split,
            "predictor_order": i,
            "predictor_type": spec["type"],
            "predictor": spec["display"],
            "predictor_internal_col": col,
            "source_variable": spec["source"],
            "transformation": spec["transformation"],
            "n_analysis_rows": n_common,
            "n_cases": n_cases,
            "n_noncases": n_noncases,
            "event_rate": float(n_cases / float(n_common)),
            "median_followup_years": float(common[duration_col].median()),
            "predictor_raw_mean_common_sample": raw_mean,
            "predictor_raw_sd_common_sample": raw_sd,
            "coef_per_1sd": coef,
            "coef_se": se,
            "hr_per_1sd": hr,
            "hr_ci_lower": ci_lo,
            "hr_ci_upper": ci_hi,
            "p_value": p,
            "cindex_baseline_covariates": base_cindex,
            "cindex_baseline_plus_predictor": cindex,
            "delta_cindex_vs_baseline": cindex - base_cindex,
            "lrt_chisq_vs_baseline": lrt_chi,
            "lrt_p_vs_baseline": lrt_p,
            "ph_test_p_predictor": predictor_ph_pvalue(cph, design, "predictor_z"),
            "baseline_penalizer": base_pen,
            "used_penalizer": used_pen,
            "baseline_covariates": ";".join(covariate_cols),
        }
        rows.append(result)
        info("[{:02d}/11] {:38s} HR={:.3f} [{:.3f}, {:.3f}] P={:.3g} C={:.4f} delta-C={:+.4f}".format(
            i, spec["display"][:38], hr, ci_lo, ci_hi, p, cindex, cindex - base_cindex
        ))

    results = pd.DataFrame(rows)
    results["p_fdr_bh_all_11_predictors"] = bh_fdr(results["p_value"].values)
    results["p_bonferroni_all_11_predictors"] = bonferroni(results["p_value"].values)
    conventional = results["predictor_type"] == "Conventional biomarker"
    results["p_fdr_bh_10_conventional_only"] = np.nan
    results["p_bonferroni_10_conventional_only"] = np.nan
    results.loc[conventional, "p_fdr_bh_10_conventional_only"] = bh_fdr(results.loc[conventional, "p_value"].values)
    results.loc[conventional, "p_bonferroni_10_conventional_only"] = bonferroni(results.loc[conventional, "p_value"].values)
    results.to_csv(summary_out, sep="\t", index=False, na_rep="NA")

    # -------------------------------------------------------------------------
    # Joint 10-biomarker panel and incremental EPOCH beyond the entire panel.
    # -------------------------------------------------------------------------
    panel_design = base_design.copy()
    biomarker_z_names = []
    for j, c in enumerate(biomarker_cols, start=1):
        zname = "bio{:02d}_z".format(j)
        biomarker_z_names.append(zname)
        panel_design[zname] = standardized[c]

    panel_cph, panel_c, panel_pen, panel_risk = fit_cox_with_retry(
        panel_design, duration_col, event_col, initial_penalizer=args.penalizer
    )
    save_cox_coefficients(panel_cph, panel_coef_out)
    risks["10-biomarker panel"] = panel_risk

    panel_epoch_design = panel_design.copy()
    panel_epoch_design["epoch_z"] = standardized["epoch_raw"]
    panel_epoch_cph, panel_epoch_c, panel_epoch_pen, panel_epoch_risk = fit_cox_with_retry(
        panel_epoch_design, duration_col, event_col, initial_penalizer=args.penalizer
    )
    save_cox_coefficients(panel_epoch_cph, panel_epoch_coef_out)
    risks["10-biomarker panel + Brain EPOCH"] = panel_epoch_risk

    epoch_row = results.loc[results["predictor_type"] == "EPOCH"].iloc[0]
    panel_lrt_chi, panel_lrt_p = valid_lrt(panel_cph, base_cph, panel_pen, base_pen, len(biomarker_z_names))
    add_epoch_chi, add_epoch_p = valid_lrt(panel_epoch_cph, panel_cph, panel_epoch_pen, panel_pen, 1)
    panel_epoch_stats = panel_epoch_cph.summary.loc["epoch_z"]

    model_rows = [
        {
            "disease_code": args.disease_code,
            "disease_label": args.disease_label,
            "model": "Baseline covariates",
            "cindex": base_cindex,
            "delta_cindex_vs_baseline": 0.0,
            "penalizer": base_pen,
            "lrt_chisq_vs_reduced": np.nan,
            "lrt_p_vs_reduced": np.nan,
            "reduced_model": "",
        },
        {
            "disease_code": args.disease_code,
            "disease_label": args.disease_label,
            "model": "Baseline covariates + Brain EPOCH",
            "cindex": float(epoch_row["cindex_baseline_plus_predictor"]),
            "delta_cindex_vs_baseline": float(epoch_row["delta_cindex_vs_baseline"]),
            "penalizer": float(epoch_row["used_penalizer"]),
            "lrt_chisq_vs_reduced": float(epoch_row["lrt_chisq_vs_baseline"]) if pd.notna(epoch_row["lrt_chisq_vs_baseline"]) else np.nan,
            "lrt_p_vs_reduced": float(epoch_row["lrt_p_vs_baseline"]) if pd.notna(epoch_row["lrt_p_vs_baseline"]) else np.nan,
            "reduced_model": "Baseline covariates",
        },
        {
            "disease_code": args.disease_code,
            "disease_label": args.disease_label,
            "model": "Baseline covariates + 10-biomarker panel",
            "cindex": panel_c,
            "delta_cindex_vs_baseline": panel_c - base_cindex,
            "penalizer": panel_pen,
            "lrt_chisq_vs_reduced": panel_lrt_chi,
            "lrt_p_vs_reduced": panel_lrt_p,
            "reduced_model": "Baseline covariates",
        },
        {
            "disease_code": args.disease_code,
            "disease_label": args.disease_label,
            "model": "Baseline covariates + 10-biomarker panel + Brain EPOCH",
            "cindex": panel_epoch_c,
            "delta_cindex_vs_baseline": panel_epoch_c - base_cindex,
            "penalizer": panel_epoch_pen,
            "lrt_chisq_vs_reduced": add_epoch_chi,
            "lrt_p_vs_reduced": add_epoch_p,
            "reduced_model": "Baseline covariates + 10-biomarker panel",
        },
    ]
    pd.DataFrame(model_rows).to_csv(model_out, sep="\t", index=False, na_rep="NA")

    epoch_joint_beta = float(panel_epoch_stats["coef"])
    epoch_joint_se = float(panel_epoch_stats["se(coef)"])
    epoch_joint_hr = float(np.exp(epoch_joint_beta))
    epoch_joint_ci_lo = float(np.exp(float(panel_epoch_stats["coef lower 95%"])))
    epoch_joint_ci_hi = float(np.exp(float(panel_epoch_stats["coef upper 95%"])))
    epoch_joint_p = float(panel_epoch_stats["p"])

    # Direct paired C-index comparisons.
    time = common[duration_col].to_numpy(dtype=float)
    event = common[event_col].to_numpy(dtype=int)
    pairwise_rows = []
    epoch_model_name = "Brain proteomics mortality EPOCH"
    for j, c in enumerate(biomarker_cols, start=1):
        biomarker_name = display_names[c]
        row = paired_bootstrap_delta_cindex(
            time, event,
            risks[epoch_model_name], risks[biomarker_name],
            "Baseline + Brain EPOCH", "Baseline + {}".format(biomarker_name),
            args.n_bootstrap, args.random_state + j,
        )
        row.update({"comparison_type": "EPOCH_vs_individual_biomarker", "disease_code": args.disease_code, "disease_label": args.disease_label})
        pairwise_rows.append(row)

    row = paired_bootstrap_delta_cindex(
        time, event,
        risks[epoch_model_name], risks["10-biomarker panel"],
        "Baseline + Brain EPOCH", "Baseline + 10-biomarker panel",
        args.n_bootstrap, args.random_state + 100,
    )
    row.update({"comparison_type": "EPOCH_vs_10_biomarker_panel", "disease_code": args.disease_code, "disease_label": args.disease_label})
    pairwise_rows.append(row)

    row_panel = paired_bootstrap_delta_cindex(
        time, event,
        risks["10-biomarker panel + Brain EPOCH"], risks["10-biomarker panel"],
        "Baseline + 10-biomarker panel + Brain EPOCH", "Baseline + 10-biomarker panel",
        args.n_bootstrap, args.random_state + 200,
    )
    row_panel.update({"comparison_type": "incremental_EPOCH_beyond_10_biomarker_panel", "disease_code": args.disease_code, "disease_label": args.disease_label})
    pairwise_rows.append(row_panel)
    pairwise_df = pd.DataFrame(pairwise_rows)
    pairwise_df.to_csv(pairwise_out, sep="\t", index=False, na_rep="NA")

    incremental = pd.DataFrame([{
        "disease_code": args.disease_code,
        "disease_label": args.disease_label,
        "analysis_split": args.analysis_split,
        "N": n_common,
        "N_cases": n_cases,
        "cindex_10_biomarker_panel": panel_c,
        "cindex_10_biomarker_panel_plus_epoch": panel_epoch_c,
        "delta_cindex_epoch_beyond_panel": panel_epoch_c - panel_c,
        "delta_cindex_ci_lower": row_panel["delta_cindex_ci_lower"],
        "delta_cindex_ci_upper": row_panel["delta_cindex_ci_upper"],
        "delta_cindex_bootstrap_p": row_panel["empirical_p_two_sided"],
        "epoch_beta_in_joint_panel_model": epoch_joint_beta,
        "epoch_se_in_joint_panel_model": epoch_joint_se,
        "epoch_hr_per_1sd_in_joint_panel_model": epoch_joint_hr,
        "epoch_hr_ci_lower": epoch_joint_ci_lo,
        "epoch_hr_ci_upper": epoch_joint_ci_hi,
        "epoch_wald_p_in_joint_panel_model": epoch_joint_p,
        "lrt_chisq_adding_epoch_to_10_biomarker_panel": add_epoch_chi,
        "lrt_p_adding_epoch_to_10_biomarker_panel": add_epoch_p,
        "panel_penalizer": panel_pen,
        "panel_plus_epoch_penalizer": panel_epoch_pen,
        "lrt_note": "Formal LRT reported only when both nested Cox fits are unpenalized.",
    }])
    incremental.to_csv(incremental_out, sep="\t", index=False, na_rep="NA")

    # QC / provenance.
    baseline_date_discrepancy_days = np.nan
    both_dates = common["epoch_baseline_date"].notna() & common["ukb_baseline_date"].notna()
    if both_dates.any():
        baseline_date_discrepancy_days = float(
            np.nanmedian(np.abs((common.loc[both_dates, "epoch_baseline_date"] - common.loc[both_dates, "ukb_baseline_date"]).dt.days))
        )

    qc_rows = [
        ("disease_code", args.disease_code),
        ("disease_label", args.disease_label),
        ("disease_tsv", str(Path(args.disease_tsv).resolve())),
        ("epoch_predictions", str(Path(args.epoch_predictions).resolve())),
        ("epoch_source_column", epoch_source_col),
        ("analysis_split", args.analysis_split),
        ("field53_col", field53_col),
        ("death_date_cols", ";".join(death_cols)),
        ("rows_after_all_source_merges", n_after_merge),
        ("rows_missing_reported_case_status_excluded", n_missing_case_status),
        ("reported_cases_missing_event_date_excluded", n_case_missing_date),
        ("prevalent_or_on_baseline_cases_excluded", n_prevalent),
        ("rows_valid_followup", n_valid_followup),
        ("common_sample_rows", n_common),
        ("common_sample_cases", n_cases),
        ("common_sample_noncases", n_noncases),
        ("common_sample_event_rate", float(n_cases / float(n_common))),
        ("common_sample_median_followup_years", float(common[duration_col].median())),
        ("median_absolute_epoch_vs_field53_baseline_date_difference_days", baseline_date_discrepancy_days),
        ("disease_link_method", disease_link_method),
        ("disease_id_col", disease_id_col),
        ("covariate_link_method", cov_link_method),
        ("covariate_id_col", cov_id_col),
        ("biomarker_link_method", bio_link_method),
        ("biomarker_id_col", bio_id_col),
        ("baseline_covariates", ";".join(covariate_cols)),
        ("baseline_covariate_source_summary", covariate_source_desc),
        ("baseline_cindex", base_cindex),
        ("baseline_penalizer_used", base_pen),
        ("common_sample_rule", "complete Brain EPOCH + all 10 conventional biomarkers + valid disease follow-up"),
        ("blood_pressure_covariates_in_default_baseline", "No; systolic BP is a benchmark predictor"),
        ("primary_comparison_population", "held-out mortality-EPOCH test split" if args.analysis_split == "test" else "full EPOCH sample: mortality-clock train + validation + test participants"),
        ("out_of_sample_interpretation", "Yes: held-out from mortality-EPOCH development" if args.analysis_split == "test" else "No: full-sample analysis includes mortality-EPOCH development participants; use for maximum-power association/incremental comparison"),
    ]
    for c in biomarker_cols:
        qc_rows.append(("source_{}".format(c), biomarker_sources[c]))
    pd.DataFrame(qc_rows, columns=["metric", "value"]).to_csv(qc_out, sep="\t", index=False)

    info("=" * 90)
    info("Finished disease-onset comparison")
    info("Individual predictors: {}".format(summary_out))
    info("Model panel summary: {}".format(model_out))
    info("EPOCH beyond 10-biomarker panel: {}".format(incremental_out))
    info("Paired C-index comparisons: {}".format(pairwise_out))
    info("QC: {}".format(qc_out))
    info("=" * 90)


# =============================================================================
# CLI
# =============================================================================

def main():
    p = argparse.ArgumentParser(
        description=(
            "Compare brain proteomics mortality EPOCH with 10 conventional biomarkers "
            "for incident disease onset in one identical comparison population. Full sample is the default to maximize N."
        )
    )
    p.add_argument("--epoch-predictions", required=True)
    p.add_argument("--epoch-col", default=None)
    p.add_argument(
        "--analysis-split", choices=["test", "all"], default="all",
        help=(
            "Default 'all' uses train + validation + test participants from the full EPOCH predictions file "
            "to maximize sample size. Use 'test' for a held-out sensitivity analysis."
        ),
    )

    p.add_argument("--disease-tsv", required=True)
    p.add_argument("--disease-code", required=True)
    p.add_argument("--disease-label", required=True)
    p.add_argument("--case-col", default="case")
    p.add_argument("--date-col", default="date")

    p.add_argument("--death-xlsx", required=True)
    p.add_argument("--id-match-csv", required=True)
    p.add_argument("--covariate-csv", required=True)
    p.add_argument("--biomarker-csv", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--admin-censor-date", default="2022-11-30")

    p.add_argument("--field53-0-col", default=None)
    p.add_argument("--death-date-col", default=None)
    p.add_argument("--death-id-col", default=None)
    p.add_argument("--idmatch-score-col", default=None)
    p.add_argument("--idmatch-death-col", default=None)
    p.add_argument("--covariate-id-col", default=None)
    p.add_argument("--biomarker-id-col", default=None)
    p.add_argument("--covariate-cols", default=None)

    p.add_argument("--penalizer", type=float, default=0.0)
    p.add_argument("--min-events", type=int, default=20)
    p.add_argument("--n-bootstrap", type=int, default=1000)
    p.add_argument("--random-state", type=int, default=2026)

    args = p.parse_args()
    run(args)


if __name__ == "__main__":
    main()