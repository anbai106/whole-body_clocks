#!/usr/bin/env python3
"""
Apple-to-apple landmark survival comparison for all-cause mortality:

    Stroke hepatic-proteomics disease-specific EPOCH
    versus
    10 conventional UK Biobank mortality biomarkers

The key design feature is a SINGLE COMMON ANALYSIS SAMPLE. A participant is
included only if all of the following are available:

    - stroke hepatic-proteomics EPOCH
    - all 10 conventional biomarkers
    - valid landmark date / mortality follow-up
    - participant linkage to the covariate/death data

Every predictor is then evaluated in a separate Cox model on exactly this
same set of participants and exactly the same deaths.

Landmark definition follows the existing disease-clock mortality analysis:
    proteomics EPOCH -> UKB Field 53 instance 0_0

All-cause mortality:
    event = death after landmark and on/before administrative censoring
    non-event = censored at the administrative censor date

Primary comparison model:
    baseline covariates + one standardized predictor

Default baseline covariates:
    age at assessment
    sex
    genetic ethnic grouping
    assessment centre
    smoking status
    BMI

Systolic and diastolic BP are deliberately NOT included in the default
comparison baseline because systolic BP is itself one of the 10 benchmark
predictors. This avoids adjusting the systolic-BP benchmark for itself and
avoids asymmetric over-adjustment across predictors.

The original mortality script can still be emulated with --covariate-cols
if a different exact covariate set is desired.

Biomarker definitions at UKB baseline (instance 0):
    1. Overall health rating          field 2178, array 0
    2. Left grip strength             field 46, array 0
    3. Right grip strength            field 47, array 0
    4. Systolic blood pressure        mean of field 4080 arrays 0 and 1
    5. Peak expiratory flow           maximum of field 3064 arrays 0,1,2
    6. C-reactive protein             log1p(field 30710, array 0)
    7. Albumin                        field 30600, array 0
    8. Cystatin C                     field 30720, array 0
    9. HbA1c                          field 30750, array 0
   10. Red-cell distribution width    field 30070, array 0

All predictors, including the EPOCH, are standardized within the common
analysis sample. HRs therefore represent the hazard ratio per 1-SD higher
predictor.

Outputs:
    *_summary.tsv
    *_common_sample_ids.tsv
    *_qc.tsv
    *_baseline_model_coefficients.tsv
"""

from __future__ import print_function

import argparse
import math
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
except ImportError:
    raise ImportError(
        "Missing lifelines. Activate the same environment used for the "
        "existing mortality survival analysis."
    )

try:
    from scipy.stats import chi2
except ImportError:
    raise ImportError("Missing scipy.")


# =============================================================================
# 1. General helpers
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
    return [z.strip() for z in str(x).split(",") if z.strip() != ""]


def read_tsv_header(path):
    return list(pd.read_csv(path, sep="\t", nrows=0).columns)


def parse_date_series(s):
    out = pd.to_datetime(s, errors="coerce")

    numeric = pd.to_numeric(s, errors="coerce")
    if numeric.notna().sum() > out.notna().sum():
        out2 = pd.to_datetime(
            numeric,
            unit="D",
            origin="1899-12-30",
            errors="coerce",
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
    if n > 0:
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
            vals = set(
                v
                for v in df[col].map(clean_id).dropna().unique()
                if v is not None
            )
        except Exception:
            continue

        if len(vals) == 0:
            continue

        overlap = len(vals.intersection(target_ids))

        bonus = 0
        if prefer_patterns is not None:
            low = str(col).lower()
            for pat in prefer_patterns:
                if re.search(pat, low):
                    bonus += 1

        if overlap > best_overlap or (
            overlap == best_overlap and bonus > best_bonus
        ):
            best_col = col
            best_overlap = overlap
            best_bonus = bonus

    if best_overlap < min_overlap:
        return None, best_overlap

    return best_col, best_overlap


def make_earliest_date(df, cols):
    parsed = [parse_date_series(df[c]) for c in cols]
    return pd.concat(parsed, axis=1).min(axis=1)


# =============================================================================
# 2. Resolve the target stroke hepatic-proteomics EPOCH
# =============================================================================

def infer_modality(score_col, folder, modality):
    text = "{} {} {}".format(score_col, folder, modality).lower()
    if "mri" in text:
        return "MRI"
    if "proteomics" in text:
        return "Proteomics"
    if "metabolomics" in text:
        return "Metabolomics"
    return str(modality)


def resolve_target_clock(
    metadata_tsv,
    wide_tsv,
    good_clock_tsv,
    disease_key,
    organ_key,
    modality_key,
    score_col_override=None,
):
    wide_cols = read_tsv_header(wide_tsv)

    if "participant_id" not in wide_cols:
        raise ValueError("participant_id is missing from score-wide TSV.")

    if score_col_override:
        if score_col_override not in wide_cols:
            alt = score_col_override.replace("__", "_")
            if alt in wide_cols:
                score_col_override = alt
            else:
                raise ValueError(
                    "Requested --score-col not found in score-wide TSV: {}".format(
                        score_col_override
                    )
                )

        return {
            "disease": disease_key,
            "folder": "",
            "clock_label": "stroke hepatic proteomics EPOCH",
            "modality": "Proteomics",
            "organ_label": organ_key,
            "score_col": score_col_override,
        }

    meta = pd.read_csv(metadata_tsv, sep="\t", low_memory=False)

    if "status" in meta.columns:
        meta = meta[
            meta["status"].astype(str).str.lower() == "collected"
        ].copy()

    if good_clock_tsv and Path(good_clock_tsv).exists():
        good = pd.read_csv(good_clock_tsv, sep="\t", low_memory=False)
        if {"disease", "folder"}.issubset(set(good.columns)):
            good_key = good[["disease", "folder"]].drop_duplicates().copy()
            good_key["disease"] = good_key["disease"].astype(str).str.lower()
            good_key["folder"] = good_key["folder"].astype(str)

            if "disease" in meta.columns and "folder" in meta.columns:
                meta["disease"] = meta["disease"].astype(str).str.lower()
                meta["folder"] = meta["folder"].astype(str)
                meta = meta.merge(
                    good_key,
                    on=["disease", "folder"],
                    how="inner",
                )

    required = ["disease", "folder"]
    missing = [c for c in required if c not in meta.columns]
    if missing:
        raise ValueError(
            "Metadata TSV missing required column(s): {}".format(missing)
        )

    disease_norm = str(disease_key).strip().lower()
    candidates = meta[
        meta["disease"].astype(str).str.strip().str.lower() == disease_norm
    ].copy()

    if candidates.empty:
        raise ValueError(
            "No metadata rows found for disease='{}'.".format(disease_key)
        )

    # Search the combined metadata text. For "hepatic", also accept "liver".
    organ_terms = [str(organ_key).lower()]
    if str(organ_key).lower() == "hepatic":
        organ_terms.append("liver")

    def combined_text(row):
        vals = []
        for col in [
            "folder",
            "clock_label",
            "clock_id",
            "modality",
            "organ_label",
            "score_col_wide",
        ]:
            if col in candidates.columns and pd.notna(row.get(col, np.nan)):
                vals.append(str(row[col]))
        return " ".join(vals).lower()

    keep_rows = []
    for idx, row in candidates.iterrows():
        text = combined_text(row)

        organ_hit = any(term in text for term in organ_terms)
        modality_hit = str(modality_key).lower() in text

        if organ_hit and modality_hit:
            keep_rows.append(idx)

    candidates = candidates.loc[keep_rows].copy()

    if candidates.empty:
        raise ValueError(
            "Could not resolve the stroke hepatic-proteomics clock from metadata. "
            "Use --score-col to provide the exact score column."
        )

    resolved_rows = []

    for _, r in candidates.iterrows():
        score_col = None

        if "score_col_wide" in candidates.columns and pd.notna(
            r.get("score_col_wide", np.nan)
        ):
            candidate = str(r["score_col_wide"])
            if candidate in wide_cols:
                score_col = candidate
            elif candidate.replace("__", "_") in wide_cols:
                score_col = candidate.replace("__", "_")

        if score_col is None:
            clock_label = (
                str(r["clock_label"])
                if "clock_label" in candidates.columns
                else str(r["folder"])
            )
            candidate = "{}_{}_clock_acceleration_z".format(
                clean_col_name(disease_key),
                clean_col_name(clock_label),
            )
            if candidate in wide_cols:
                score_col = candidate

        if score_col is not None:
            resolved_rows.append(
                {
                    "disease": str(r["disease"]),
                    "folder": str(r["folder"]),
                    "clock_label": (
                        str(r["clock_label"])
                        if "clock_label" in candidates.columns
                        else str(r["folder"])
                    ),
                    "modality": (
                        str(r["modality"])
                        if "modality" in candidates.columns
                        else infer_modality(score_col, str(r["folder"]), "")
                    ),
                    "organ_label": (
                        str(r["organ_label"])
                        if "organ_label" in candidates.columns
                        else organ_key
                    ),
                    "score_col": score_col,
                }
            )

    # Remove exact duplicate resolved score columns.
    unique = {}
    for r in resolved_rows:
        unique[r["score_col"]] = r
    resolved_rows = list(unique.values())

    if len(resolved_rows) == 0:
        raise ValueError(
            "Candidate stroke hepatic-proteomics metadata row(s) were found, "
            "but no score column could be resolved in the wide score TSV."
        )

    if len(resolved_rows) > 1:
        msg = [
            "Multiple stroke hepatic-proteomics score columns were resolved:",
        ]
        for r in resolved_rows:
            msg.append(
                "  score_col={} folder={} clock_label={}".format(
                    r["score_col"],
                    r["folder"],
                    r["clock_label"],
                )
            )
        msg.append("Pass --score-col to select one explicitly.")
        raise ValueError("\n".join(msg))

    return resolved_rows[0]


# =============================================================================
# 3. Date/death-field detection
# =============================================================================

def detect_field53_col(df, instance, user_col=None):
    if user_col:
        if user_col not in df.columns:
            raise ValueError(
                "Requested Field 53 column not found: {}".format(user_col)
            )
        return user_col

    candidates = []

    patterns = [
        r"(^|[^0-9])53([^0-9]+){}([^0-9]+)0([^0-9]|$)".format(instance),
        r"f[._-]?53[._-]?{}[._-]?0".format(instance),
        r"53_{}_0".format(instance),
        r"53-{}\.0".format(instance),
        r"53\.{}\.0".format(instance),
        r"date.*attend.*{}.*0".format(instance),
        r"assessment.*date.*{}.*0".format(instance),
    ]

    for col in df.columns:
        low = str(col).lower()
        score = 0

        for pat in patterns:
            if re.search(pat, low):
                score += 25

        if "53" in low:
            score += 5
        if "date" in low:
            score += 4
        if (
            "attend" in low
            or "assessment" in low
            or "centre" in low
            or "center" in low
        ):
            score += 4
        if str(instance) in low:
            score += 2

        if score == 0:
            continue

        parsed = parse_date_series(df[col])
        n_dates = int(parsed.notna().sum())

        if n_dates > 100:
            score += 20
        elif n_dates > 10:
            score += 10

        candidates.append((score, n_dates, col))

    candidates = sorted(
        candidates,
        key=lambda x: (x[0], x[1]),
        reverse=True,
    )

    if len(candidates) == 0:
        raise ValueError(
            "Could not auto-detect Field 53 instance {}_0. "
            "Pass --field53-0-col.".format(instance)
        )

    best = candidates[0][2]
    info(
        "Detected Field 53 instance {}_0 column: {}".format(
            instance,
            best,
        )
    )
    return best


def detect_death_date_cols(df, user_col=None):
    if user_col:
        if user_col not in df.columns:
            raise ValueError(
                "Requested death-date column not found: {}".format(user_col)
            )
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

    candidates = sorted(
        candidates,
        key=lambda x: (x[0], x[1]),
        reverse=True,
    )

    if len(candidates) == 0:
        raise ValueError(
            "Could not auto-detect death-date column. "
            "Pass --death-date-col."
        )

    out = [x[2] for x in candidates if x[1] > 0]
    info("Detected death-date column(s): {}".format("; ".join(out)))
    return out


# =============================================================================
# 4. ID linkage
# =============================================================================

def detect_score_id_col_in_match(id_match_df, score_ids, override=None):
    if override:
        if override not in id_match_df.columns:
            raise ValueError(
                "Requested ID-match score-ID column not found: {}".format(
                    override
                )
            )
        return override

    col, overlap = choose_column_by_id_overlap(
        id_match_df,
        score_ids,
        min_overlap=20,
        prefer_patterns=[
            r"penn",
            r"eid",
            r"ukb",
            r"participant",
            r"id",
        ],
    )

    if col is None:
        raise ValueError(
            "Could not detect the score-ID column in ID-match CSV."
        )

    info(
        "ID-match score-ID column: {} overlap={}".format(
            col,
            overlap,
        )
    )
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
            raise ValueError(
                "death_id_col not found: {}".format(death_id_col_arg)
            )
        death_id_col = death_id_col_arg
        direct_overlap = len(
            set(
                death_df[death_id_col].map(clean_id).dropna()
            ).intersection(score_ids)
        )
    else:
        death_id_col, direct_overlap = choose_column_by_id_overlap(
            death_df,
            score_ids,
            min_overlap=20,
            prefer_patterns=[r"eid", r"id", r"participant"],
        )

    if death_id_col is not None and direct_overlap >= 20:
        info("Death file uses score IDs directly.")
        info(
            "Death ID column: {} overlap={}".format(
                death_id_col,
                direct_overlap,
            )
        )
        death_df["participant_id"] = death_df[death_id_col].map(clean_id)
        return death_df

    info("Death file does not directly match score IDs. Using ID-match CSV.")

    idmatch_score_col = detect_score_id_col_in_match(
        id_match_df,
        score_ids,
        override=idmatch_score_col_arg,
    )

    if death_id_col_arg:
        death_id_col = death_id_col_arg
    else:
        death_id_col = None

    if idmatch_death_col_arg:
        if idmatch_death_col_arg not in id_match_df.columns:
            raise ValueError(
                "idmatch_death_col not found: {}".format(
                    idmatch_death_col_arg
                )
            )
        idmatch_death_col = idmatch_death_col_arg
    else:
        idmatch_death_col = None

    if death_id_col is None or idmatch_death_col is None:
        best_overlap = -1
        best_death_col = None
        best_match_col = None

        for dcol in death_df.columns:
            death_vals = set(
                death_df[dcol].map(clean_id).dropna().unique()
            )
            if len(death_vals) == 0:
                continue

            for mcol in id_match_df.columns:
                if mcol == idmatch_score_col:
                    continue

                match_vals = set(
                    id_match_df[mcol].map(clean_id).dropna().unique()
                )
                overlap = len(death_vals.intersection(match_vals))

                if overlap > best_overlap:
                    best_overlap = overlap
                    best_death_col = dcol
                    best_match_col = mcol

        if best_overlap < 20:
            raise ValueError(
                "Could not detect death-ID mapping through ID-match CSV."
            )

        if death_id_col is None:
            death_id_col = best_death_col
        if idmatch_death_col is None:
            idmatch_death_col = best_match_col

    info("Death-file ID column: {}".format(death_id_col))
    info(
        "ID-match death/Melbourne-ID column: {}".format(
            idmatch_death_col
        )
    )

    map_df = id_match_df[
        [idmatch_score_col, idmatch_death_col]
    ].copy()
    map_df["participant_id"] = map_df[idmatch_score_col].map(clean_id)
    map_df["_death_merge_id"] = map_df[idmatch_death_col].map(clean_id)
    map_df = (
        map_df[
            ["participant_id", "_death_merge_id"]
        ]
        .dropna()
        .drop_duplicates()
    )

    death_df["_death_merge_id"] = death_df[death_id_col].map(clean_id)
    out = death_df.merge(
        map_df,
        on="_death_merge_id",
        how="left",
    )

    info(
        "Death rows mapped to score IDs: {:,}".format(
            out["participant_id"].notna().sum()
        )
    )
    return out


def map_generic_table_to_score_ids(
    df,
    id_match_df,
    score_ids,
    table_label,
    table_id_col_arg=None,
    idmatch_score_col_arg=None,
):
    """
    Map a covariate/biomarker table to the clock participant IDs.

    First preference is a direct ID overlap. If that fails, the function uses
    the ID-match CSV and automatically finds the match-table column that
    overlaps the source table.
    """
    out = df.copy()

    if table_id_col_arg:
        if table_id_col_arg not in out.columns:
            raise ValueError(
                "{} ID column not found: {}".format(
                    table_label,
                    table_id_col_arg,
                )
            )
        source_id_col = table_id_col_arg
        direct_overlap = len(
            set(
                out[source_id_col].map(clean_id).dropna()
            ).intersection(score_ids)
        )
    else:
        source_id_col, direct_overlap = choose_column_by_id_overlap(
            out,
            score_ids,
            min_overlap=20,
            prefer_patterns=[
                r"^eid$",
                r"participant",
                r"ukb",
                r"id",
            ],
        )

    if source_id_col is not None and direct_overlap >= 20:
        info(
            "{} uses score IDs directly: {} overlap={}".format(
                table_label,
                source_id_col,
                direct_overlap,
            )
        )
        out["participant_id"] = out[source_id_col].map(clean_id)
        return out, source_id_col, "direct"

    # No direct match: resolve through ID-match CSV.
    info(
        "{} does not directly match score IDs. Using ID-match CSV.".format(
            table_label
        )
    )

    # Choose a plausible source ID column if one was not explicitly supplied.
    if source_id_col is None:
        preferred = [
            c
            for c in out.columns
            if re.search(r"(^eid$|participant|ukb|id)", str(c).lower())
        ]
        if len(preferred) == 0:
            preferred = list(out.columns)

        # Keep source columns that have enough non-missing ID-like values.
        source_candidates = []
        for c in preferred:
            vals = set(out[c].map(clean_id).dropna().unique())
            if len(vals) >= 20:
                source_candidates.append((c, vals))
    else:
        source_candidates = [
            (
                source_id_col,
                set(out[source_id_col].map(clean_id).dropna().unique()),
            )
        ]

    if len(source_candidates) == 0:
        raise ValueError(
            "Could not identify an ID column in {}.".format(table_label)
        )

    idmatch_score_col = detect_score_id_col_in_match(
        id_match_df,
        score_ids,
        override=idmatch_score_col_arg,
    )

    best = None

    for source_col, source_vals in source_candidates:
        for match_col in id_match_df.columns:
            if match_col == idmatch_score_col:
                continue

            match_vals = set(
                id_match_df[match_col].map(clean_id).dropna().unique()
            )
            overlap = len(source_vals.intersection(match_vals))

            if best is None or overlap > best[0]:
                best = (
                    overlap,
                    source_col,
                    match_col,
                )

    if best is None or best[0] < 20:
        raise ValueError(
            "Could not link {} to score IDs through ID-match CSV.".format(
                table_label
            )
        )

    overlap, source_id_col, match_source_col = best

    info(
        "{} source ID column: {}".format(
            table_label,
            source_id_col,
        )
    )
    info(
        "ID-match column corresponding to {}: {} overlap={}".format(
            table_label,
            match_source_col,
            overlap,
        )
    )

    map_df = id_match_df[
        [idmatch_score_col, match_source_col]
    ].copy()
    map_df["participant_id"] = map_df[idmatch_score_col].map(clean_id)
    map_df["_source_merge_id"] = map_df[match_source_col].map(clean_id)
    map_df = (
        map_df[
            ["participant_id", "_source_merge_id"]
        ]
        .dropna()
        .drop_duplicates()
    )

    out["_source_merge_id"] = out[source_id_col].map(clean_id)
    out = out.merge(
        map_df,
        on="_source_merge_id",
        how="left",
    )

    return out, source_id_col, "id_match_csv"


# =============================================================================
# 5. Covariate selection
# =============================================================================

def detect_first_matching_col(df, patterns):
    for pat in patterns:
        for col in df.columns:
            if re.search(pat, str(col).lower()):
                return col
    return None


def exact_or_none(df, col):
    return col if col in df.columns else None


def detect_age_col(cov_df, instance):
    preferred = "age_when_attended_assessment_centre_f21003_{}_0".format(
        instance
    )
    if preferred in cov_df.columns:
        return preferred

    primary_patterns = [
        r"21003.*{}.*0".format(instance),
        r"f[._-]?21003[._-]?{}[._-]?0".format(instance),
        r"age.*{}.*0".format(instance),
    ]

    fallback_patterns = [
        r"age_at_recruitment",
        r"age_when_attended_assessment",
        r"^age$",
        r"age",
    ]

    for patterns in [primary_patterns, fallback_patterns]:
        hits = []

        for col in cov_df.columns:
            low = str(col).lower()
            if not any(re.search(pat, low) for pat in patterns):
                continue

            vals = pd.to_numeric(cov_df[col], errors="coerce")
            if (
                vals.notna().sum() > 100
                and vals.between(20, 100).mean() > 0.50
            ):
                hits.append((vals.notna().sum(), col))

        if len(hits) > 0:
            hits = sorted(hits, reverse=True)
            return hits[0][1]

    return None


def select_comparison_covariates(
    cov_df,
    instance,
    covariate_cols_arg=None,
):
    """
    Default apple-to-apple comparison baseline.

    We intentionally exclude systolic and diastolic BP from the default
    baseline because systolic BP is one of the benchmark predictors.
    """
    cov_df = cov_df.copy()

    if covariate_cols_arg:
        cols = split_csv(covariate_cols_arg)
        missing = [x for x in cols if x not in cov_df.columns]
        if missing:
            raise ValueError(
                "Requested covariates not found: {}".format(missing)
            )
        return (
            cov_df,
            cols,
            "User-specified exact comparison covariates: {}".format(
                ";".join(cols)
            ),
        )

    selected = []
    source_records = []

    # 1. Age at assessment.
    age_col = detect_age_col(cov_df, instance)
    if age_col is not None:
        selected.append(age_col)
        source_records.append("Age={}".format(age_col))
    else:
        warn("Age covariate was not auto-detected.")

    # 2. Sex.
    sex_col = exact_or_none(cov_df, "sex_f31_0_0")
    if sex_col is None:
        sex_col = detect_first_matching_col(
            cov_df,
            [
                r"^sex$",
                r"31.*0.*0",
                r"genetic_sex",
                r"reported_sex",
            ],
        )
    if sex_col is not None:
        selected.append(sex_col)
        source_records.append("Sex={}".format(sex_col))
    else:
        warn("Sex covariate was not auto-detected.")

    # 3. Genetic ethnic grouping.
    ethnicity_col = exact_or_none(
        cov_df,
        "genetic_ethnic_grouping_f22006_0_0",
    )
    if ethnicity_col is None:
        ethnicity_col = detect_first_matching_col(
            cov_df,
            [
                r"genetic_ethnic_grouping",
                r"ethnic",
                r"race",
                r"21000",
            ],
        )
    if ethnicity_col is not None:
        selected.append(ethnicity_col)
        source_records.append("Ethnicity={}".format(ethnicity_col))
    else:
        warn("Ethnicity covariate was not auto-detected.")

    # 4. Assessment centre.
    assessment_center_col = exact_or_none(
        cov_df,
        "uk_biobank_assessment_centre_f54_{}_0".format(instance),
    )
    if assessment_center_col is None:
        assessment_center_col = detect_first_matching_col(
            cov_df,
            [
                r"54.*{}.*0".format(instance),
                r"assessment.*center",
                r"assessment.*centre",
            ],
        )
    if assessment_center_col is not None:
        selected.append(assessment_center_col)
        source_records.append(
            "Assessment_center={}".format(assessment_center_col)
        )
    else:
        warn("Assessment-center covariate was not auto-detected.")

    # 5. Smoking.
    smoking_col = exact_or_none(
        cov_df,
        "smoking_status_f20116_{}_0".format(instance),
    )
    if smoking_col is None:
        smoking_col = detect_first_matching_col(
            cov_df,
            [
                r"smoking_status.*{}.*0".format(instance),
                r"20116.*{}.*0".format(instance),
                r"smoking_status",
            ],
        )
    if smoking_col is not None:
        selected.append(smoking_col)
        source_records.append("Smoking={}".format(smoking_col))
    else:
        warn("Smoking-status covariate was not auto-detected.")

    # 6. BMI.
    bmi_candidates = [
        "body_mass_index_bmi_f23104_{}_0".format(instance),
        "body_mass_index_bmi_f21001_{}_0".format(instance),
    ]
    bmi_col = None
    for c in bmi_candidates:
        if c in cov_df.columns:
            bmi_col = c
            break

    if bmi_col is None:
        bmi_col = detect_first_matching_col(
            cov_df,
            [
                r"body_mass_index.*{}.*0".format(instance),
                r"bmi.*{}.*0".format(instance),
                r"body_mass_index",
                r"bmi",
            ],
        )

    if bmi_col is not None:
        selected.append(bmi_col)
        source_records.append("BMI={}".format(bmi_col))
    else:
        warn("BMI covariate was not auto-detected.")

    # De-duplicate.
    out = []
    for c in selected:
        if c in cov_df.columns and c not in out:
            out.append(c)

    return (
        cov_df,
        out,
        "; ".join(source_records),
    )


def force_categorical_covariate(col):
    low = str(col).lower()
    patterns = [
        "sex",
        "genetic_ethnic_grouping",
        "ethnic_background",
        "assessment_centre",
        "assessment_center",
        "smoking_status",
    ]
    return any(p in low for p in patterns)


def build_covariate_design(
    df,
    duration_col,
    event_col,
    covariate_cols,
):
    """
    Match the existing mortality pipeline:
      - numeric covariates -> median imputation
      - categorical covariates -> Missing category + dummy coding

    The common participant set is fixed BEFORE this function is called.
    """
    d = df[
        [duration_col, event_col] + list(covariate_cols)
    ].copy()

    d = d.replace([np.inf, -np.inf], np.nan)

    y = d[[duration_col, event_col]].reset_index(drop=True)
    x = d[list(covariate_cols)].copy()

    processed = []

    for col in x.columns:
        s = x[col]

        if force_categorical_covariate(col):
            cat = (
                s.astype("object")
                .where(s.notna(), "Missing")
                .astype(str)
            )
            n_unique = cat.nunique(dropna=False)

            if n_unique <= 1:
                continue
            if n_unique > 80:
                warn(
                    "Skipping high-cardinality categorical covariate "
                    "{} with {} levels.".format(col, n_unique)
                )
                continue

            dummy = pd.get_dummies(
                cat,
                prefix=clean_col_name(col),
                drop_first=True,
            )
            processed.append(dummy.astype(float))
            continue

        numeric = pd.to_numeric(s, errors="coerce")
        numeric_fraction = numeric.notna().mean()

        if (
            numeric_fraction >= 0.90
            and numeric.nunique(dropna=True) > 5
        ):
            med = numeric.median()
            if pd.isna(med):
                continue
            numeric = numeric.fillna(med).astype(float)
            processed.append(
                pd.DataFrame(
                    {clean_col_name(col): numeric},
                    index=x.index,
                )
            )
        else:
            cat = (
                s.astype("object")
                .where(s.notna(), "Missing")
                .astype(str)
            )
            n_unique = cat.nunique(dropna=False)

            if n_unique <= 1:
                continue
            if n_unique > 80:
                warn(
                    "Skipping high-cardinality categorical covariate "
                    "{} with {} levels.".format(col, n_unique)
                )
                continue

            dummy = pd.get_dummies(
                cat,
                prefix=clean_col_name(col),
                drop_first=True,
            )
            processed.append(dummy.astype(float))

    if len(processed) == 0:
        x_processed = pd.DataFrame(index=x.index)
    else:
        x_processed = pd.concat(processed, axis=1)

    keep = [
        c
        for c in x_processed.columns
        if x_processed[c].nunique(dropna=True) > 1
    ]
    x_processed = x_processed[keep].reset_index(drop=True)

    return pd.concat([y, x_processed], axis=1)


# =============================================================================
# 6. Construct the 10 conventional biomarkers
# =============================================================================

def exact_or_prefix(df, exact, prefix=None):
    if exact in df.columns:
        return exact

    if prefix is None:
        prefix = exact

    hits = [c for c in df.columns if str(c).startswith(prefix)]
    if len(hits) == 0:
        return None

    # Prefer instance 0 array 0 when exact naming varies slightly.
    hits = sorted(
        hits,
        key=lambda c: (
            0 if re.search(r"_0_0$", str(c)) else 1,
            str(c),
        ),
    )
    return hits[0]


def positive_numeric(s):
    x = pd.to_numeric(s, errors="coerce")
    x = x.where(np.isfinite(x), np.nan)
    return x.where(x > 0, np.nan)


def build_conventional_biomarkers(bio_df):
    df = bio_df.copy()
    source = {}

    # 1. Overall health rating, UKB 2178:
    # typically 1=Excellent, 2=Good, 3=Fair, 4=Poor.
    c = exact_or_prefix(
        df,
        "overall_health_rating_f2178_0_0",
        "overall_health_rating_f2178_0_",
    )
    if c is None:
        raise ValueError("Overall health rating field 2178 was not found.")
    health = pd.to_numeric(df[c], errors="coerce")
    health = health.where(health.isin([1, 2, 3, 4]), np.nan)
    df["biomarker_overall_health_rating"] = health
    source["biomarker_overall_health_rating"] = c

    # 2. Left grip.
    c = exact_or_prefix(
        df,
        "hand_grip_strength_left_f46_0_0",
        "hand_grip_strength_left_f46_0_",
    )
    if c is None:
        raise ValueError("Left grip-strength field 46 was not found.")
    df["biomarker_grip_strength_left"] = positive_numeric(df[c])
    source["biomarker_grip_strength_left"] = c

    # 3. Right grip.
    c = exact_or_prefix(
        df,
        "hand_grip_strength_right_f47_0_0",
        "hand_grip_strength_right_f47_0_",
    )
    if c is None:
        raise ValueError("Right grip-strength field 47 was not found.")
    df["biomarker_grip_strength_right"] = positive_numeric(df[c])
    source["biomarker_grip_strength_right"] = c

    # 4. Mean systolic BP from automated readings, field 4080 instance 0.
    sbp_cols = [
        c
        for c in df.columns
        if re.match(
            r"^systolic_blood_pressure_automated_reading_f4080_0_[01]$",
            str(c),
        )
    ]
    if len(sbp_cols) == 0:
        raise ValueError(
            "Baseline automated systolic-BP field 4080 was not found."
        )
    sbp_mat = pd.concat(
        [positive_numeric(df[c]) for c in sorted(sbp_cols)],
        axis=1,
    )
    df["biomarker_systolic_bp"] = sbp_mat.mean(axis=1, skipna=True)
    source["biomarker_systolic_bp"] = ";".join(sorted(sbp_cols))

    # 5. Peak expiratory flow: best/max baseline trial, field 3064.
    pef_cols = [
        c
        for c in df.columns
        if re.match(
            r"^peak_expiratory_flow_pef_f3064_0_[012]$",
            str(c),
        )
    ]
    if len(pef_cols) == 0:
        raise ValueError(
            "Baseline peak-expiratory-flow field 3064 was not found."
        )
    pef_mat = pd.concat(
        [positive_numeric(df[c]) for c in sorted(pef_cols)],
        axis=1,
    )
    df["biomarker_peak_expiratory_flow"] = pef_mat.max(
        axis=1,
        skipna=True,
    )
    source["biomarker_peak_expiratory_flow"] = ";".join(
        sorted(pef_cols)
    )

    # 6. CRP, log1p transformed because of strong right skew.
    c = exact_or_prefix(
        df,
        "creactive_protein_f30710_0_0",
        "creactive_protein_f30710_0_",
    )
    if c is None:
        raise ValueError("C-reactive-protein field 30710 was not found.")
    crp = pd.to_numeric(df[c], errors="coerce")
    crp = crp.where(np.isfinite(crp), np.nan)
    crp = crp.where(crp >= 0, np.nan)
    df["biomarker_log1p_crp"] = np.log1p(crp)
    source["biomarker_log1p_crp"] = "{} [log1p]".format(c)

    # 7. Albumin.
    c = exact_or_prefix(
        df,
        "albumin_f30600_0_0",
        "albumin_f30600_0_",
    )
    if c is None:
        raise ValueError("Albumin field 30600 was not found.")
    df["biomarker_albumin"] = positive_numeric(df[c])
    source["biomarker_albumin"] = c

    # 8. Cystatin C.
    c = exact_or_prefix(
        df,
        "cystatin_c_f30720_0_0",
        "cystatin_c_f30720_0_",
    )
    if c is None:
        raise ValueError("Cystatin-C field 30720 was not found.")
    df["biomarker_cystatin_c"] = positive_numeric(df[c])
    source["biomarker_cystatin_c"] = c

    # 9. HbA1c.
    c = exact_or_prefix(
        df,
        "glycated_haemoglobin_hba1c_f30750_0_0",
        "glycated_haemoglobin_hba1c_f30750_0_",
    )
    if c is None:
        raise ValueError("HbA1c field 30750 was not found.")
    df["biomarker_hba1c"] = positive_numeric(df[c])
    source["biomarker_hba1c"] = c

    # 10. RDW.
    c = exact_or_prefix(
        df,
        "red_blood_cell_erythrocyte_distribution_width_f30070_0_0",
        "red_blood_cell_erythrocyte_distribution_width_f30070_0_",
    )
    if c is None:
        raise ValueError(
            "Red-cell distribution width field 30070 was not found."
        )
    df["biomarker_rdw"] = positive_numeric(df[c])
    source["biomarker_rdw"] = c

    biomarker_cols = [
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

    display_names = {
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

    return df, biomarker_cols, display_names, source


# =============================================================================
# 7. Cox model helpers
# =============================================================================

def standardize_predictor(s):
    x = pd.to_numeric(s, errors="coerce")
    mean = float(x.mean())
    sd = float(x.std(ddof=1))

    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("Predictor has zero/invalid SD.")

    return (x - mean) / sd, mean, sd


def fit_cox_with_retry(
    design_df,
    duration_col,
    event_col,
    initial_penalizer,
):
    candidates = []
    for p in [
        initial_penalizer,
        0.001,
        0.01,
        0.05,
        0.1,
    ]:
        if p not in candidates:
            candidates.append(p)

    last_error = None

    for penalizer in candidates:
        try:
            cph = CoxPHFitter(penalizer=penalizer)
            cph.fit(
                design_df,
                duration_col=duration_col,
                event_col=event_col,
                show_progress=False,
            )

            risk = cph.predict_partial_hazard(
                design_df
            ).values.reshape(-1)

            cindex = concordance_index(
                event_times=design_df[duration_col].values,
                predicted_scores=-risk,
                event_observed=design_df[event_col].values,
            )

            return cph, float(cindex), float(penalizer)

        except Exception as exc:
            last_error = exc

    raise last_error


def lrt_pvalue(full_cph, base_cph):
    stat = 2.0 * (
        float(full_cph.log_likelihood_)
        - float(base_cph.log_likelihood_)
    )
    if not np.isfinite(stat) or stat < 0:
        return np.nan
    return float(chi2.sf(stat, 1))


def predictor_ph_pvalue(cph, design_df, predictor_col="predictor_z"):
    try:
        test = proportional_hazard_test(
            cph,
            design_df,
            time_transform="rank",
        )
        return float(test.summary.loc[predictor_col, "p"])
    except Exception:
        return np.nan


# =============================================================================
# 8. Main analysis
# =============================================================================

def run(args):
    score_wide_tsv = Path(args.score_wide_tsv).resolve()
    metadata_tsv = Path(args.score_metadata_tsv).resolve()
    good_clock_tsv = (
        Path(args.good_clock_tsv).resolve()
        if args.good_clock_tsv
        else None
    )
    death_xlsx = Path(args.death_xlsx).resolve()
    id_match_csv = Path(args.id_match_csv).resolve()
    covariate_csv = Path(args.covariate_csv).resolve()
    biomarker_csv = Path(args.biomarker_csv).resolve()
    outdir = Path(args.output_dir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    target = resolve_target_clock(
        metadata_tsv=metadata_tsv,
        wide_tsv=score_wide_tsv,
        good_clock_tsv=good_clock_tsv,
        disease_key=args.disease_key,
        organ_key=args.organ_key,
        modality_key=args.modality_key,
        score_col_override=args.score_col,
    )

    score_col = target["score_col"]
    modality = infer_modality(
        score_col,
        target["folder"],
        target["modality"],
    )

    if modality.lower() != "proteomics":
        raise ValueError(
            "Resolved target is not proteomics: modality={}".format(
                modality
            )
        )

    field53_instance = 0

    prefix = "stroke_hepatic_proteomics_EPOCH_vs_10_biomarkers_mortality"
    summary_out = outdir / "{}_summary.tsv".format(prefix)
    ids_out = outdir / "{}_common_sample_ids.tsv".format(prefix)
    qc_out = outdir / "{}_qc.tsv".format(prefix)
    base_coef_out = outdir / "{}_baseline_model_coefficients.tsv".format(
        prefix
    )

    info("=" * 88)
    info("Apple-to-apple mortality comparison")
    info("=" * 88)
    info("Disease: {}".format(target["disease"]))
    info("Folder: {}".format(target["folder"]))
    info("Clock label: {}".format(target["clock_label"]))
    info("Organ label: {}".format(target["organ_label"]))
    info("Modality: {}".format(modality))
    info("Score column: {}".format(score_col))
    info("Landmark Field 53 instance: 0_0")
    info("Output directory: {}".format(outdir))
    info("=" * 88)

    # -------------------------------------------------------------------------
    # Read target EPOCH only.
    # -------------------------------------------------------------------------
    score_df = pd.read_csv(
        score_wide_tsv,
        sep="\t",
        usecols=["participant_id", score_col],
        dtype={"participant_id": "str"},
        low_memory=False,
    )

    score_df["participant_id"] = score_df["participant_id"].map(clean_id)
    score_df = score_df.dropna(subset=["participant_id"]).copy()
    score_df = score_df.drop_duplicates(subset=["participant_id"])
    score_df = score_df.rename(columns={score_col: "epoch_raw"})
    score_df["epoch_raw"] = pd.to_numeric(
        score_df["epoch_raw"],
        errors="coerce",
    )

    score_ids = set(score_df["participant_id"].dropna().unique())

    info("Score rows: {:,}".format(score_df.shape[0]))
    info(
        "Non-missing stroke hepatic-proteomics EPOCH: {:,}".format(
            score_df["epoch_raw"].notna().sum()
        )
    )

    # -------------------------------------------------------------------------
    # Read mortality/linkage/covariate/biomarker data.
    # -------------------------------------------------------------------------
    death_df_raw = pd.read_excel(
        death_xlsx,
        sheet_name=0,
        engine="openpyxl",
    )
    id_match_df = pd.read_csv(
        id_match_csv,
        low_memory=False,
    )
    cov_df_raw = pd.read_csv(
        covariate_csv,
        low_memory=False,
    )
    bio_df_raw = pd.read_csv(
        biomarker_csv,
        low_memory=False,
    )

    # -------------------------------------------------------------------------
    # Death mapping, exactly following the existing mortality pipeline logic.
    # -------------------------------------------------------------------------
    death_df = map_death_to_score_ids(
        death_df=death_df_raw,
        id_match_df=id_match_df,
        score_ids=score_ids,
        death_id_col_arg=args.death_id_col,
        idmatch_score_col_arg=args.idmatch_score_col,
        idmatch_death_col_arg=args.idmatch_death_col,
    )

    field53_col = detect_field53_col(
        death_df,
        instance=field53_instance,
        user_col=args.field53_0_col,
    )

    death_date_cols = detect_death_date_cols(
        death_df,
        user_col=args.death_date_col,
    )

    death_df["baseline_date"] = parse_date_series(
        death_df[field53_col]
    )
    death_df["death_date"] = make_earliest_date(
        death_df,
        death_date_cols,
    )

    death_keep = (
        death_df[
            ["participant_id", "baseline_date", "death_date"]
        ]
        .dropna(subset=["participant_id"])
        .copy()
    )

    death_keep = death_keep.groupby(
        "participant_id",
        as_index=False,
    ).agg(
        {
            "baseline_date": "min",
            "death_date": "min",
        }
    )

    # -------------------------------------------------------------------------
    # Covariate linkage + selection.
    # -------------------------------------------------------------------------
    cov_df, cov_id_col, cov_link_method = map_generic_table_to_score_ids(
        df=cov_df_raw,
        id_match_df=id_match_df,
        score_ids=score_ids,
        table_label="Covariate file",
        table_id_col_arg=args.covariate_id_col,
        idmatch_score_col_arg=args.idmatch_score_col,
    )

    cov_df, covariate_cols, covariate_source_desc = (
        select_comparison_covariates(
            cov_df,
            instance=field53_instance,
            covariate_cols_arg=args.covariate_cols,
        )
    )

    if len(covariate_cols) == 0:
        raise ValueError(
            "No baseline comparison covariates were selected."
        )

    info("Comparison baseline covariates:")
    for c in covariate_cols:
        info("  - {}".format(c))
    info("Covariate source summary:")
    info("  {}".format(covariate_source_desc))

    cov_keep = (
        cov_df[
            ["participant_id"] + covariate_cols
        ]
        .dropna(subset=["participant_id"])
        .drop_duplicates(subset=["participant_id"])
        .copy()
    )

    # -------------------------------------------------------------------------
    # Biomarker linkage + construction.
    # -------------------------------------------------------------------------
    bio_df, bio_id_col, bio_link_method = map_generic_table_to_score_ids(
        df=bio_df_raw,
        id_match_df=id_match_df,
        score_ids=score_ids,
        table_label="10-biomarker file",
        table_id_col_arg=args.biomarker_id_col,
        idmatch_score_col_arg=args.idmatch_score_col,
    )

    bio_df, biomarker_cols, display_names, biomarker_sources = (
        build_conventional_biomarkers(bio_df)
    )

    bio_keep = (
        bio_df[
            ["participant_id"] + biomarker_cols
        ]
        .dropna(subset=["participant_id"])
        .drop_duplicates(subset=["participant_id"])
        .copy()
    )

    info("Constructed conventional biomarkers:")
    for c in biomarker_cols:
        info(
            "  - {} <= {}".format(
                display_names[c],
                biomarker_sources[c],
            )
        )

    # -------------------------------------------------------------------------
    # Merge all sources BEFORE defining the common comparison population.
    # -------------------------------------------------------------------------
    dat = score_df.merge(
        death_keep,
        on="participant_id",
        how="left",
    )
    dat = dat.merge(
        cov_keep,
        on="participant_id",
        how="left",
    )
    dat = dat.merge(
        bio_keep,
        on="participant_id",
        how="left",
    )

    n_after_all_merges = int(dat.shape[0])

    # -------------------------------------------------------------------------
    # Mortality follow-up: same landmark rule as existing proteomics analysis.
    # -------------------------------------------------------------------------
    admin_censor_date = pd.to_datetime(args.admin_censor_date)

    dat = dat[dat["baseline_date"].notna()].copy()
    dat = dat[dat["baseline_date"] <= admin_censor_date].copy()

    # Exclude deaths before/on landmark.
    dat = dat[
        dat["death_date"].isna()
        | (dat["death_date"] > dat["baseline_date"])
    ].copy()

    dat["event"] = (
        dat["death_date"].notna()
        & (dat["death_date"] > dat["baseline_date"])
        & (dat["death_date"] <= admin_censor_date)
    ).astype(int)

    dat["end_date"] = dat["death_date"].where(
        dat["event"] == 1,
        admin_censor_date,
    )

    dat["followup_time_years"] = (
        dat["end_date"] - dat["baseline_date"]
    ).dt.days / 365.25

    dat = dat[
        dat["followup_time_years"].notna()
        & (dat["followup_time_years"] > 0)
    ].copy()

    n_valid_followup = int(dat.shape[0])

    # -------------------------------------------------------------------------
    # Strict common predictor sample.
    #
    # This is the apple-to-apple step:
    # EPOCH + ALL 10 conventional biomarkers must be observed.
    # -------------------------------------------------------------------------
    predictor_cols = ["epoch_raw"] + biomarker_cols

    common_mask = np.ones(dat.shape[0], dtype=bool)

    for c in predictor_cols:
        x = pd.to_numeric(dat[c], errors="coerce")
        common_mask &= x.notna().values & np.isfinite(x).values

    common = dat.loc[common_mask].copy()

    if common.empty:
        raise ValueError(
            "The common EPOCH + 10-biomarker analysis sample is empty."
        )

    n_common = int(common.shape[0])
    n_deaths = int(common["event"].sum())

    if n_deaths < args.min_events:
        warn(
            "Only {} mortality events in common sample.".format(n_deaths)
        )

    info("=" * 88)
    info("COMMON APPLE-TO-APPLE ANALYSIS SAMPLE")
    info("=" * 88)
    info("Rows after all source merges: {:,}".format(n_after_all_merges))
    info("Rows with valid mortality follow-up: {:,}".format(n_valid_followup))
    info("Rows complete for EPOCH + all 10 biomarkers: {:,}".format(n_common))
    info("Deaths in common sample: {:,}".format(n_deaths))
    info(
        "Event rate: {:.3f}%".format(
            100.0 * n_deaths / float(n_common)
        )
    )
    info(
        "Median follow-up: {:.3f} years".format(
            float(np.nanmedian(common["followup_time_years"]))
        )
    )
    info("=" * 88)

    # Save exact analysis IDs so the population is auditable.
    common[
        [
            "participant_id",
            "baseline_date",
            "death_date",
            "end_date",
            "followup_time_years",
            "event",
        ]
    ].to_csv(
        ids_out,
        sep="\t",
        index=False,
    )

    # -------------------------------------------------------------------------
    # Build and fit the baseline covariate model ONCE.
    # -------------------------------------------------------------------------
    duration_col = "followup_time_years"
    event_col = "event"

    base_design = build_covariate_design(
        common,
        duration_col=duration_col,
        event_col=event_col,
        covariate_cols=covariate_cols,
    )

    if int(base_design.shape[0]) != n_common:
        raise ValueError(
            "Baseline design unexpectedly changed common-sample N."
        )

    base_cph, base_cindex, base_penalizer = fit_cox_with_retry(
        base_design,
        duration_col=duration_col,
        event_col=event_col,
        initial_penalizer=args.penalizer,
    )

    base_coef = base_cph.summary.copy()
    base_coef.insert(0, "term", base_coef.index)
    base_coef = base_coef.reset_index(drop=True)
    base_coef.to_csv(
        base_coef_out,
        sep="\t",
        index=False,
    )

    # -------------------------------------------------------------------------
    # Fit one predictor at a time on EXACTLY the same rows.
    # -------------------------------------------------------------------------
    predictor_specs = [
        {
            "column": "epoch_raw",
            "display": "Stroke hepatic-proteomics EPOCH",
            "type": "EPOCH",
            "source": score_col,
            "transformation": "standardized within common sample",
        }
    ]

    transformations = {
        "biomarker_overall_health_rating": (
            "ordinal 1-4; standardized within common sample"
        ),
        "biomarker_grip_strength_left": (
            "raw baseline value; standardized within common sample"
        ),
        "biomarker_grip_strength_right": (
            "raw baseline value; standardized within common sample"
        ),
        "biomarker_systolic_bp": (
            "mean of baseline automated readings; standardized within common sample"
        ),
        "biomarker_peak_expiratory_flow": (
            "maximum baseline trial; standardized within common sample"
        ),
        "biomarker_log1p_crp": (
            "log1p baseline CRP; standardized within common sample"
        ),
        "biomarker_albumin": (
            "raw baseline value; standardized within common sample"
        ),
        "biomarker_cystatin_c": (
            "raw baseline value; standardized within common sample"
        ),
        "biomarker_hba1c": (
            "raw baseline value; standardized within common sample"
        ),
        "biomarker_rdw": (
            "raw baseline value; standardized within common sample"
        ),
    }

    for c in biomarker_cols:
        predictor_specs.append(
            {
                "column": c,
                "display": display_names[c],
                "type": "Conventional biomarker",
                "source": biomarker_sources[c],
                "transformation": transformations[c],
            }
        )

    rows = []

    for i, spec in enumerate(predictor_specs, start=1):
        col = spec["column"]

        z, raw_mean, raw_sd = standardize_predictor(common[col])

        full_design = base_design.copy()
        full_design["predictor_z"] = np.asarray(z, dtype=float)

        if int(full_design.shape[0]) != n_common:
            raise ValueError(
                "Predictor {} changed analysis N.".format(spec["display"])
            )

        cph, cindex, used_penalizer = fit_cox_with_retry(
            full_design,
            duration_col=duration_col,
            event_col=event_col,
            initial_penalizer=args.penalizer,
        )

        row = cph.summary.loc["predictor_z"]

        coef = float(row["coef"])
        se = float(row["se(coef)"])
        hr = float(np.exp(coef))
        ci_lower = float(np.exp(float(row["coef lower 95%"])))
        ci_upper = float(np.exp(float(row["coef upper 95%"])))
        p_value = float(row["p"])

        result = {
            "predictor_order": i,
            "predictor_type": spec["type"],
            "predictor": spec["display"],
            "predictor_internal_col": col,
            "source_variable": spec["source"],
            "transformation": spec["transformation"],
            "n_analysis_rows": n_common,
            "n_deaths": n_deaths,
            "event_rate": float(n_deaths / float(n_common)),
            "median_followup_years": float(
                np.nanmedian(common[duration_col])
            ),
            "predictor_raw_mean_common_sample": raw_mean,
            "predictor_raw_sd_common_sample": raw_sd,
            "coef_per_1sd": coef,
            "coef_se": se,
            "hr_per_1sd": hr,
            "hr_ci_lower": ci_lower,
            "hr_ci_upper": ci_upper,
            "p_value": p_value,
            "cindex_baseline_covariates": base_cindex,
            "cindex_baseline_plus_predictor": cindex,
            "delta_cindex_vs_baseline": cindex - base_cindex,
            "lrt_p_vs_baseline": lrt_pvalue(cph, base_cph),
            "ph_test_p_predictor": predictor_ph_pvalue(
                cph,
                full_design,
                predictor_col="predictor_z",
            ),
            "landmark_field53_instance": "0_0",
            "landmark_field53_column": field53_col,
            "admin_censor_date": args.admin_censor_date,
            "baseline_covariates": ";".join(covariate_cols),
            "baseline_covariate_source_summary": covariate_source_desc,
            "initial_penalizer": args.penalizer,
            "used_penalizer": used_penalizer,
            "score_col": score_col,
            "clock_folder": target["folder"],
            "clock_label": target["clock_label"],
        }

        rows.append(result)

        info(
            "[{:02d}/{:02d}] {:38s} "
            "HR={:.3f} [{:.3f}, {:.3f}] "
            "P={:.3g} C={:.4f} delta-C={:+.4f}".format(
                i,
                len(predictor_specs),
                spec["display"][:38],
                hr,
                ci_lower,
                ci_upper,
                p_value,
                cindex,
                cindex - base_cindex,
            )
        )

    results = pd.DataFrame(rows)

    results["p_fdr_bh_all_11_predictors"] = bh_fdr(
        results["p_value"].values
    )
    results["p_bonferroni_all_11_predictors"] = bonferroni(
        results["p_value"].values
    )

    results["p_fdr_bh_10_conventional_only"] = np.nan
    results["p_bonferroni_10_conventional_only"] = np.nan

    conventional_mask = (
        results["predictor_type"] == "Conventional biomarker"
    )

    conventional_p = results.loc[
        conventional_mask,
        "p_value",
    ].values

    results.loc[
        conventional_mask,
        "p_fdr_bh_10_conventional_only",
    ] = bh_fdr(conventional_p)

    results.loc[
        conventional_mask,
        "p_bonferroni_10_conventional_only",
    ] = bonferroni(conventional_p)

    results.to_csv(
        summary_out,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    # -------------------------------------------------------------------------
    # QC.
    # -------------------------------------------------------------------------
    qc_rows = [
        ("target_disease", target["disease"]),
        ("target_folder", target["folder"]),
        ("target_clock_label", target["clock_label"]),
        ("target_modality", modality),
        ("target_score_col", score_col),
        ("field53_instance", "0_0"),
        ("field53_col", field53_col),
        ("death_date_cols", ";".join(death_date_cols)),
        ("score_rows", int(score_df.shape[0])),
        (
            "score_nonmissing_epoch",
            int(score_df["epoch_raw"].notna().sum()),
        ),
        ("rows_after_all_source_merges", n_after_all_merges),
        ("rows_valid_followup", n_valid_followup),
        ("common_sample_rows", n_common),
        ("common_sample_deaths", n_deaths),
        (
            "common_sample_event_rate",
            float(n_deaths / float(n_common)),
        ),
        (
            "common_sample_median_followup_years",
            float(np.nanmedian(common[duration_col])),
        ),
        ("covariate_link_method", cov_link_method),
        ("covariate_id_col", cov_id_col),
        ("biomarker_link_method", bio_link_method),
        ("biomarker_id_col", bio_id_col),
        ("baseline_covariates", ";".join(covariate_cols)),
        (
            "baseline_covariate_source_summary",
            covariate_source_desc,
        ),
        ("baseline_cindex", base_cindex),
        ("baseline_penalizer_used", base_penalizer),
        (
            "common_sample_rule",
            "complete EPOCH + all 10 conventional biomarkers",
        ),
        (
            "blood_pressure_covariates_in_default_baseline",
            "No; SBP is a benchmark predictor",
        ),
    ]

    for c in biomarker_cols:
        qc_rows.append(
            (
                "source_{}".format(c),
                biomarker_sources[c],
            )
        )

    pd.DataFrame(
        qc_rows,
        columns=["metric", "value"],
    ).to_csv(
        qc_out,
        sep="\t",
        index=False,
    )

    info("=" * 88)
    info("Finished apple-to-apple mortality comparison")
    info("=" * 88)
    info("Summary:")
    info("  {}".format(summary_out))
    info("Common sample IDs:")
    info("  {}".format(ids_out))
    info("QC:")
    info("  {}".format(qc_out))
    info("Baseline model coefficients:")
    info("  {}".format(base_coef_out))
    info("=" * 88)


# =============================================================================
# 9. CLI
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compare stroke hepatic-proteomics EPOCH with 10 conventional "
            "UKB mortality biomarkers in one identical landmark population."
        )
    )

    parser.add_argument("--score-wide-tsv", required=True)
    parser.add_argument("--score-metadata-tsv", required=True)
    parser.add_argument("--good-clock-tsv", default=None)

    parser.add_argument("--death-xlsx", required=True)
    parser.add_argument("--id-match-csv", required=True)
    parser.add_argument("--covariate-csv", required=True)
    parser.add_argument("--biomarker-csv", required=True)
    parser.add_argument("--output-dir", required=True)

    parser.add_argument("--disease-key", default="stroke")
    parser.add_argument("--organ-key", default="hepatic")
    parser.add_argument("--modality-key", default="proteomics")
    parser.add_argument(
        "--score-col",
        default=None,
        help=(
            "Optional exact stroke hepatic-proteomics EPOCH score column. "
            "If omitted, it is resolved from the metadata."
        ),
    )

    parser.add_argument(
        "--admin-censor-date",
        default="2022-11-30",
    )

    parser.add_argument("--field53-0-col", default=None)
    parser.add_argument("--death-date-col", default=None)
    parser.add_argument("--death-id-col", default=None)
    parser.add_argument("--idmatch-score-col", default=None)
    parser.add_argument("--idmatch-death-col", default=None)
    parser.add_argument("--covariate-id-col", default=None)
    parser.add_argument("--biomarker-id-col", default=None)

    parser.add_argument(
        "--covariate-cols",
        default=None,
        help=(
            "Optional comma-separated exact baseline covariates. "
            "Default comparison baseline is age, sex, genetic ethnic "
            "grouping, assessment centre, smoking and BMI. BP is excluded "
            "because systolic BP is a benchmark predictor."
        ),
    )

    parser.add_argument("--penalizer", type=float, default=0.01)
    parser.add_argument("--min-events", type=int, default=20)

    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
