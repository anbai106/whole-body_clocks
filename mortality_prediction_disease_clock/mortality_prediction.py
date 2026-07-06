#!/usr/bin/env python3
# ============================================================
# Run mortality Cox survival analysis for one stable/significant
# disease-clock acceleration-z score.
#
# Time zero:
#   MRI clocks                  -> UKB Field 53 instance 2_0
#   Proteomics/metabolomics     -> UKB Field 53 instance 0_0
#
# Field 53 is read from the Melbourne death-related file.
#
# Default covariates:
#   age at assessment, sex, genetic ethnic grouping,
#   assessment center, smoking status, BMI,
#   mean diastolic BP, mean systolic BP
#
# Important revision:
#   - Genetic PCs are NOT included by default.
#   - Smoking, BMI, diastolic BP, and systolic BP are added.
#   - Instance-specific covariates are used:
#       MRI: 2_0
#       Proteomics/metabolomics: 0_0
#
# Output:
#   <BASE_DIR>/<clock_folder>/survival_analysis_mortality/
# ============================================================

from __future__ import print_function

import argparse
import re
import sys
from pathlib import Path
from typing import Any, List, Optional

import numpy as np
import pandas as pd

try:
    from lifelines import CoxPHFitter
    from lifelines.utils import concordance_index
except ImportError:
    raise ImportError("Missing lifelines. Install with: pip install lifelines")


# ============================================================
# 1. Helpers
# ============================================================

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
    s = s.strip("_")
    return s


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


def choose_column_by_id_overlap(df, target_ids, min_overlap=20, prefer_patterns=None, exclude_cols=None):
    if exclude_cols is None:
        exclude_cols = set()

    best_col = None
    best_overlap = -1
    best_bonus = -1

    for col in df.columns:
        if col in exclude_cols:
            continue

        try:
            vals = df[col].map(clean_id)
        except Exception:
            continue

        vals = set(v for v in vals.dropna().unique() if v is not None)

        if len(vals) == 0:
            continue

        overlap = len(vals.intersection(target_ids))

        bonus = 0
        if prefer_patterns is not None:
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


def get_task_row(tasks_tsv, task_index):
    tasks = pd.read_csv(tasks_tsv, sep="\t")

    if "array_id" not in tasks.columns:
        raise ValueError("tasks TSV must contain array_id.")

    hit = tasks[tasks["array_id"].astype(int) == int(task_index)]

    if hit.empty:
        raise ValueError("No task found for array_id={}".format(task_index))

    return hit.iloc[0]


def infer_modality(score_col, folder, modality):
    text = "{} {} {}".format(score_col, folder, modality).lower()

    if "mri" in text:
        return "MRI"

    if "proteomics" in text:
        return "Proteomics"

    if "metabolomics" in text:
        return "Metabolomics"

    return str(modality)


def infer_field53_instance(score_col, folder, modality):
    modality2 = infer_modality(score_col, folder, modality)

    if modality2 == "MRI":
        return 2

    return 0


# ============================================================
# 2. Date and field detection
# ============================================================

def detect_field53_col(df, instance, user_col=None):
    if user_col is not None and user_col != "":
        if user_col not in df.columns:
            raise ValueError("Requested Field 53 column not found: {}".format(user_col))
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

        if "attend" in low or "assessment" in low or "centre" in low or "center" in low:
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

    candidates = sorted(candidates, key=lambda x: (x[0], x[1]), reverse=True)

    if len(candidates) == 0:
        raise ValueError(
            "Could not auto-detect Field 53 instance {}_0. "
            "Please pass --field53-{}-col.".format(instance, instance)
        )

    best = candidates[0][2]

    info("Detected Field 53 instance {}_0 column: {}".format(instance, best))

    return best


def detect_death_date_cols(df, user_col=None):
    if user_col is not None and user_col != "":
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

    candidates = sorted(candidates, key=lambda x: (x[0], x[1]), reverse=True)

    if len(candidates) == 0:
        raise ValueError("Could not auto-detect death-date column. Pass --death-date-col.")

    out = [x[2] for x in candidates if x[1] > 0]

    info("Detected death-date column(s): {}".format("; ".join(out)))

    return out


def make_earliest_date(df, cols):
    parsed = []

    for col in cols:
        parsed.append(parse_date_series(df[col]))

    mat = pd.concat(parsed, axis=1)

    return mat.min(axis=1)


# ============================================================
# 3. ID mapping
# ============================================================

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

    if death_id_col_arg is not None and death_id_col_arg != "":
        if death_id_col_arg not in death_df.columns:
            raise ValueError("death_id_col not found: {}".format(death_id_col_arg))
        death_id_col = death_id_col_arg
        direct_overlap = len(set(death_df[death_id_col].map(clean_id).dropna()).intersection(score_ids))
    else:
        death_id_col, direct_overlap = choose_column_by_id_overlap(
            death_df,
            score_ids,
            min_overlap=20,
            prefer_patterns=[r"eid", r"id", r"participant"],
        )

    if death_id_col is not None and direct_overlap >= 20:
        info("Death file uses score IDs directly.")
        info("Death ID column: {} overlap={}".format(death_id_col, direct_overlap))
        death_df["participant_id"] = death_df[death_id_col].map(clean_id)
        return death_df

    info("Death file does not directly match score IDs. Using ID-match CSV.")

    if idmatch_score_col_arg is not None and idmatch_score_col_arg != "":
        if idmatch_score_col_arg not in id_match_df.columns:
            raise ValueError("idmatch_score_col not found: {}".format(idmatch_score_col_arg))
        idmatch_score_col = idmatch_score_col_arg
        score_overlap = len(set(id_match_df[idmatch_score_col].map(clean_id).dropna()).intersection(score_ids))
    else:
        idmatch_score_col, score_overlap = choose_column_by_id_overlap(
            id_match_df,
            score_ids,
            min_overlap=20,
            prefer_patterns=[r"penn", r"eid", r"ukb", r"participant", r"id"],
        )

    if idmatch_score_col is None:
        raise ValueError("Could not detect score-ID column in ID-match CSV.")

    info("ID-match score-ID column: {} overlap={}".format(idmatch_score_col, score_overlap))

    if death_id_col_arg is not None and death_id_col_arg != "":
        death_id_col = death_id_col_arg
    else:
        death_id_col = None

    if idmatch_death_col_arg is not None and idmatch_death_col_arg != "":
        if idmatch_death_col_arg not in id_match_df.columns:
            raise ValueError("idmatch_death_col not found: {}".format(idmatch_death_col_arg))
        idmatch_death_col = idmatch_death_col_arg
    else:
        idmatch_death_col = None

    if death_id_col is None or idmatch_death_col is None:
        best_overlap = -1
        best_death_col = None
        best_match_col = None

        for dcol in death_df.columns:
            death_vals = set(death_df[dcol].map(clean_id).dropna().unique())
            if len(death_vals) == 0:
                continue

            for mcol in id_match_df.columns:
                if mcol == idmatch_score_col:
                    continue

                match_vals = set(id_match_df[mcol].map(clean_id).dropna().unique())
                overlap = len(death_vals.intersection(match_vals))

                if overlap > best_overlap:
                    best_overlap = overlap
                    best_death_col = dcol
                    best_match_col = mcol

        if best_overlap < 20:
            raise ValueError(
                "Could not detect death-ID/Melbourne-ID mapping. "
                "Pass --death-id-col and --idmatch-death-col."
            )

        if death_id_col is None:
            death_id_col = best_death_col

        if idmatch_death_col is None:
            idmatch_death_col = best_match_col

    info("Death-file ID column: {}".format(death_id_col))
    info("ID-match death/Melbourne-ID column: {}".format(idmatch_death_col))

    map_df = id_match_df[[idmatch_score_col, idmatch_death_col]].copy()
    map_df["participant_id"] = map_df[idmatch_score_col].map(clean_id)
    map_df["_death_merge_id"] = map_df[idmatch_death_col].map(clean_id)
    map_df = map_df[["participant_id", "_death_merge_id"]].dropna().drop_duplicates()

    death_df["_death_merge_id"] = death_df[death_id_col].map(clean_id)

    out = death_df.merge(map_df, on="_death_merge_id", how="left")

    info("Death rows mapped to score IDs: {:,}".format(out["participant_id"].notna().sum()))

    return out


def map_covariates_to_score_ids(cov_df, score_ids, covariate_id_col_arg=None):
    cov_df = cov_df.copy()

    if covariate_id_col_arg is not None and covariate_id_col_arg != "":
        if covariate_id_col_arg not in cov_df.columns:
            raise ValueError("covariate_id_col not found: {}".format(covariate_id_col_arg))
        cov_id_col = covariate_id_col_arg
        overlap = len(set(cov_df[cov_id_col].map(clean_id).dropna()).intersection(score_ids))
    else:
        cov_id_col, overlap = choose_column_by_id_overlap(
            cov_df,
            score_ids,
            min_overlap=20,
            prefer_patterns=[r"eid", r"ukb", r"participant", r"id"],
        )

    if cov_id_col is None:
        raise ValueError("Could not detect covariate ID column. Pass --covariate-id-col.")

    info("Covariate ID column: {} overlap={}".format(cov_id_col, overlap))

    cov_df["participant_id"] = cov_df[cov_id_col].map(clean_id)

    return cov_df, cov_id_col


# ============================================================
# 4. Covariate selection and preprocessing
# ============================================================

def detect_first_matching_col(df, patterns):
    for pat in patterns:
        for col in df.columns:
            if re.search(pat, str(col).lower()):
                return col
    return None


def exact_or_none(df, col):
    if col in df.columns:
        return col
    return None


def detect_age_col(cov_df, instance):
    preferred = "age_when_attended_assessment_centre_f21003_{}_0".format(instance)
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
            if vals.notna().sum() > 100 and vals.between(20, 100).mean() > 0.50:
                hits.append((vals.notna().sum(), col))

        if len(hits) > 0:
            hits = sorted(hits, reverse=True)
            return hits[0][1]

    return None


def add_mean_column(cov_df, new_col, source_cols):
    existing = [c for c in source_cols if c in cov_df.columns]

    if len(existing) == 0:
        return cov_df, None, ""

    vals = []
    for c in existing:
        vals.append(pd.to_numeric(cov_df[c], errors="coerce"))

    mat = pd.concat(vals, axis=1)
    cov_df[new_col] = mat.mean(axis=1, skipna=True)

    return cov_df, new_col, ";".join(existing)


def select_common_covariates(cov_df, instance, covariate_cols_arg=None):
    """
    Default clinical covariates for mortality prediction.

    Revised default:
      - age at assessment, instance-specific
      - sex
      - genetic ethnic grouping
      - assessment center, instance-specific
      - smoking status, instance-specific
      - BMI, instance-specific
      - mean diastolic BP, instance-specific
      - mean systolic BP, instance-specific

    Genetic PCs are intentionally excluded by default.

    If --covariate-cols is supplied, the user-specified list is used exactly.
    """

    cov_df = cov_df.copy()

    if covariate_cols_arg is not None and covariate_cols_arg.strip() != "":
        cols = [x.strip() for x in covariate_cols_arg.split(",") if x.strip() != ""]
        missing = [x for x in cols if x not in cov_df.columns]
        if len(missing) > 0:
            raise ValueError("Requested covariates not found: {}".format(missing))

        source_desc = "User-specified covariates: {}".format(";".join(cols))
        return cov_df, cols, source_desc

    selected = []
    source_records = []

    # 1. Age at assessment, instance-specific.
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
    ethnicity_col = exact_or_none(cov_df, "genetic_ethnic_grouping_f22006_0_0")
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

    # 4. Assessment center, instance-specific.
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
                r"assessment_center",
                r"assessment_centre",
            ],
        )

    if assessment_center_col is not None:
        selected.append(assessment_center_col)
        source_records.append("Assessment_center={}".format(assessment_center_col))
    else:
        warn("Assessment center covariate was not auto-detected.")

    # 5. Smoking status, instance-specific.
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
                r"smoking",
            ],
        )

    if smoking_col is not None:
        selected.append(smoking_col)
        source_records.append("Smoking={}".format(smoking_col))
    else:
        warn("Smoking-status covariate was not auto-detected.")

    # 6. BMI, instance-specific.
    bmi_col = exact_or_none(
        cov_df,
        "body_mass_index_bmi_f23104_{}_0".format(instance),
    )
    if bmi_col is None:
        bmi_col = detect_first_matching_col(
            cov_df,
            [
                r"body_mass_index.*{}.*0".format(instance),
                r"bmi.*{}.*0".format(instance),
                r"23104.*{}.*0".format(instance),
                r"body_mass_index",
                r"bmi",
            ],
        )

    if bmi_col is not None:
        selected.append(bmi_col)
        source_records.append("BMI={}".format(bmi_col))
    else:
        warn("BMI covariate was not auto-detected.")

    # 7. Mean diastolic BP, instance-specific.
    diastolic_source_cols = [
        "diastolic_blood_pressure_automated_reading_f4079_{}_0".format(instance),
        "diastolic_blood_pressure_automated_reading_f4079_{}_1".format(instance),
    ]
    diastolic_mean_col = "diastolic_blood_pressure_mean_f4079_{}".format(instance)

    cov_df, diastolic_col, diastolic_sources = add_mean_column(
        cov_df,
        diastolic_mean_col,
        diastolic_source_cols,
    )

    if diastolic_col is not None:
        selected.append(diastolic_col)
        source_records.append("Diastolic_mean={}".format(diastolic_sources))
    else:
        warn(
            "Diastolic BP covariate was not detected for instance {}. "
            "Expected one of: {}".format(instance, diastolic_source_cols)
        )

    # 8. Mean systolic BP, instance-specific.
    systolic_source_cols = [
        "systolic_blood_pressure_automated_reading_f4080_{}_0".format(instance),
        "systolic_blood_pressure_automated_reading_f4080_{}_1".format(instance),
    ]
    systolic_mean_col = "systolic_blood_pressure_mean_f4080_{}".format(instance)

    cov_df, systolic_col, systolic_sources = add_mean_column(
        cov_df,
        systolic_mean_col,
        systolic_source_cols,
    )

    if systolic_col is not None:
        selected.append(systolic_col)
        source_records.append("Systolic_mean={}".format(systolic_sources))
    else:
        warn(
            "Systolic BP covariate was not detected for instance {}. "
            "Expected one of: {}".format(instance, systolic_source_cols)
        )

    # De-duplicate while preserving order.
    out = []
    for c in selected:
        if c in cov_df.columns and c not in out:
            out.append(c)

    source_desc = "; ".join(source_records)

    return cov_df, out, source_desc


def force_categorical_covariate(col):
    low = str(col).lower()

    categorical_patterns = [
        "sex",
        "genetic_ethnic_grouping",
        "ethnic_background",
        "assessment_centre",
        "assessment_center",
        "smoking_status",
    ]

    return any(pat in low for pat in categorical_patterns)


def preprocess_design_matrix(df, duration_col, event_col, clock_col, covariate_cols):
    cols = [duration_col, event_col] + list(covariate_cols)

    if clock_col is not None:
        cols.append(clock_col)

    d = df[cols].copy()
    d = d.replace([np.inf, -np.inf], np.nan)
    d = d[d[duration_col].notna() & d[event_col].notna()].copy()
    d = d[d[duration_col] > 0].copy()

    if clock_col is not None:
        d = d[d[clock_col].notna()].copy()

    y = d[[duration_col, event_col]].copy()

    x_cols = list(covariate_cols)
    if clock_col is not None:
        x_cols = [clock_col] + x_cols

    x = d[x_cols].copy()

    processed = []

    for col in x.columns:
        s = x[col]

        # Force selected categorical covariates to dummy variables even if coded numerically.
        if force_categorical_covariate(col):
            cat = s.astype("object").where(s.notna(), "Missing").astype(str)
            n_unique = cat.nunique(dropna=False)

            if n_unique <= 1:
                continue

            if n_unique > 80:
                warn("Skipping high-cardinality categorical covariate {} with {} levels.".format(col, n_unique))
                continue

            dummy = pd.get_dummies(cat, prefix=clean_col_name(col), drop_first=True)
            processed.append(dummy.astype(float))
            continue

        numeric = pd.to_numeric(s, errors="coerce")
        numeric_fraction = numeric.notna().mean()

        if numeric_fraction >= 0.90 and numeric.nunique(dropna=True) > 5:
            med = numeric.median()
            if pd.isna(med):
                continue
            numeric = numeric.fillna(med).astype(float)
            processed.append(pd.DataFrame({clean_col_name(col): numeric}, index=x.index))
        else:
            cat = s.astype("object").where(s.notna(), "Missing").astype(str)
            n_unique = cat.nunique(dropna=False)

            if n_unique <= 1:
                continue

            if n_unique > 80:
                warn("Skipping high-cardinality categorical covariate {} with {} levels.".format(col, n_unique))
                continue

            dummy = pd.get_dummies(cat, prefix=clean_col_name(col), drop_first=True)
            processed.append(dummy.astype(float))

    if len(processed) == 0:
        x_processed = pd.DataFrame(index=d.index)
    else:
        x_processed = pd.concat(processed, axis=1)

    keep_cols = []
    for col in x_processed.columns:
        if x_processed[col].nunique(dropna=True) > 1:
            keep_cols.append(col)

    x_processed = x_processed[keep_cols].copy()

    out = pd.concat(
        [
            y.reset_index(drop=True),
            x_processed.reset_index(drop=True),
        ],
        axis=1,
    )

    return out


def fit_cox(design_df, duration_col, event_col, penalizer):
    cph = CoxPHFitter(penalizer=penalizer)

    cph.fit(
        design_df,
        duration_col=duration_col,
        event_col=event_col,
        show_progress=False,
    )

    risk = cph.predict_partial_hazard(design_df).values.reshape(-1)

    cindex = concordance_index(
        event_times=design_df[duration_col].values,
        predicted_scores=-risk,
        event_observed=design_df[event_col].values,
    )

    return cph, float(cindex)


# ============================================================
# 5. Main analysis
# ============================================================

def run_one_clock(args):
    base_dir = Path(args.base_dir).resolve()
    score_wide_tsv = Path(args.score_wide_tsv).resolve()
    tasks_tsv = Path(args.tasks_tsv).resolve()
    death_xlsx = Path(args.death_xlsx).resolve()
    id_match_csv = Path(args.id_match_csv).resolve()
    covariate_csv = Path(args.covariate_csv).resolve()

    task = get_task_row(tasks_tsv, args.task_index)

    disease = str(task["disease"])
    folder = str(task["folder"])
    clock_label = str(task["clock_label"])
    modality = str(task["modality"])
    score_col = str(task["score_col_wide"])

    field53_instance = infer_field53_instance(score_col, folder, modality)
    modality_inferred = infer_modality(score_col, folder, modality)

    outdir = base_dir / folder / "survival_analysis_mortality"
    outdir.mkdir(parents=True, exist_ok=True)

    safe_score = clean_col_name(score_col)

    summary_out = outdir / "{}_mortality_survival_summary.tsv".format(safe_score)
    full_coef_out = outdir / "{}_mortality_survival_full_model_coefficients.tsv".format(safe_score)
    covar_coef_out = outdir / "{}_mortality_survival_covariate_model_coefficients.tsv".format(safe_score)
    qc_out = outdir / "{}_mortality_survival_qc.tsv".format(safe_score)

    info("============================================================")
    info("Mortality survival analysis for one disease clock")
    info("============================================================")
    info("Task index: {}".format(args.task_index))
    info("Disease: {}".format(disease))
    info("Folder: {}".format(folder))
    info("Clock label: {}".format(clock_label))
    info("Modality: {}".format(modality_inferred))
    info("Score column: {}".format(score_col))
    info("Field 53 instance: {}_0".format(field53_instance))
    info("Output directory: {}".format(outdir))
    info("============================================================")

    wide_header = read_tsv_header(score_wide_tsv)

    if score_col not in wide_header:
        alt = score_col.replace("__", "_")
        if alt in wide_header:
            score_col = alt
        else:
            raise ValueError("Score column not found in wide score TSV: {}".format(score_col))

    score_df = pd.read_csv(
        score_wide_tsv,
        sep="\t",
        usecols=["participant_id", score_col],
        dtype={"participant_id": "str"},
        low_memory=False,
    )

    score_df["participant_id"] = score_df["participant_id"].map(clean_id)
    score_df = score_df.dropna(subset=["participant_id"]).copy()
    score_df = score_df.rename(columns={score_col: "clock_z"})
    score_df["clock_z"] = pd.to_numeric(score_df["clock_z"], errors="coerce")

    score_ids = set(score_df["participant_id"].dropna().unique())

    info("Score rows: {:,}".format(score_df.shape[0]))
    info("Non-missing clock_z: {:,}".format(score_df["clock_z"].notna().sum()))

    death_df_raw = pd.read_excel(death_xlsx, sheet_name=0, engine="openpyxl")
    id_match_df = pd.read_csv(id_match_csv, low_memory=False)
    cov_df_raw = pd.read_csv(covariate_csv, low_memory=False)

    death_df = map_death_to_score_ids(
        death_df=death_df_raw,
        id_match_df=id_match_df,
        score_ids=score_ids,
        death_id_col_arg=args.death_id_col,
        idmatch_score_col_arg=args.idmatch_score_col,
        idmatch_death_col_arg=args.idmatch_death_col,
    )

    if field53_instance == 2:
        field53_col = detect_field53_col(
            death_df,
            instance=2,
            user_col=args.field53_2_col,
        )
    else:
        field53_col = detect_field53_col(
            death_df,
            instance=0,
            user_col=args.field53_0_col,
        )

    death_date_cols = detect_death_date_cols(
        death_df,
        user_col=args.death_date_col,
    )

    death_df["baseline_date"] = parse_date_series(death_df[field53_col])
    death_df["death_date"] = make_earliest_date(death_df, death_date_cols)

    death_keep = death_df[
        ["participant_id", "baseline_date", "death_date"]
    ].dropna(subset=["participant_id"]).copy()

    death_keep = death_keep.groupby("participant_id", as_index=False).agg({
        "baseline_date": "min",
        "death_date": "min",
    })

    cov_df, cov_id_col = map_covariates_to_score_ids(
        cov_df_raw,
        score_ids=score_ids,
        covariate_id_col_arg=args.covariate_id_col,
    )

    cov_df, covariate_cols, covariate_source_desc = select_common_covariates(
        cov_df,
        instance=field53_instance,
        covariate_cols_arg=args.covariate_cols,
    )

    if len(covariate_cols) == 0:
        raise ValueError("No covariates selected. Pass --covariate-cols.")

    info("Selected covariates:")
    for c in covariate_cols:
        info("  - {}".format(c))

    info("Covariate source summary:")
    info("  {}".format(covariate_source_desc))

    cov_keep = cov_df[
        ["participant_id"] + covariate_cols
    ].dropna(subset=["participant_id"]).drop_duplicates(subset=["participant_id"]).copy()

    dat = score_df.merge(death_keep, on="participant_id", how="left")
    dat = dat.merge(cov_keep, on="participant_id", how="left")

    admin_censor_date = pd.to_datetime(args.admin_censor_date)

    dat = dat[dat["baseline_date"].notna()].copy()
    dat = dat[dat["baseline_date"] <= admin_censor_date].copy()

    dat["death_after_baseline"] = (
        dat["death_date"].notna()
        & (dat["death_date"] > dat["baseline_date"])
        & (dat["death_date"] <= admin_censor_date)
    )

    dat["event"] = dat["death_after_baseline"].astype(int)

    dat["end_date"] = dat["death_date"].where(dat["event"] == 1, admin_censor_date)

    dat["followup_time_years"] = (
        dat["end_date"] - dat["baseline_date"]
    ).dt.days / 365.25

    dat = dat[dat["followup_time_years"].notna()].copy()
    dat = dat[dat["followup_time_years"] > 0].copy()

    # Exclude deaths before or on baseline date.
    dat = dat[
        dat["death_date"].isna()
        | (dat["death_date"] > dat["baseline_date"])
    ].copy()

    n_before_model = int(dat.shape[0])
    n_with_clock = int(dat["clock_z"].notna().sum())
    n_deaths_with_clock = int(dat.loc[dat["clock_z"].notna(), "event"].sum())

    info("Rows with valid follow-up: {:,}".format(n_before_model))
    info("Rows with non-missing clock: {:,}".format(n_with_clock))
    info("Deaths with non-missing clock: {:,}".format(n_deaths_with_clock))

    duration_col = "followup_time_years"
    event_col = "event"

    cov_design = preprocess_design_matrix(
        df=dat,
        duration_col=duration_col,
        event_col=event_col,
        clock_col=None,
        covariate_cols=covariate_cols,
    )

    full_design = preprocess_design_matrix(
        df=dat,
        duration_col=duration_col,
        event_col=event_col,
        clock_col="clock_z",
        covariate_cols=covariate_cols,
    )

    if "clock_z" not in full_design.columns:
        raise ValueError("clock_z was lost during preprocessing.")

    if int(full_design[event_col].sum()) < args.min_events:
        warn(
            "Few mortality events: {} deaths. Model may be unstable.".format(
                int(full_design[event_col].sum())
            )
        )

    cph_cov, cindex_cov = fit_cox(
        cov_design,
        duration_col=duration_col,
        event_col=event_col,
        penalizer=args.penalizer,
    )

    cph_full, cindex_full = fit_cox(
        full_design,
        duration_col=duration_col,
        event_col=event_col,
        penalizer=args.penalizer,
    )

    full_coef = cph_full.summary.copy()
    full_coef.insert(0, "term", full_coef.index)
    full_coef = full_coef.reset_index(drop=True)

    cov_coef = cph_cov.summary.copy()
    cov_coef.insert(0, "term", cov_coef.index)
    cov_coef = cov_coef.reset_index(drop=True)

    full_coef.to_csv(full_coef_out, sep="\t", index=False)
    cov_coef.to_csv(covar_coef_out, sep="\t", index=False)

    clock_row = full_coef[full_coef["term"] == "clock_z"]

    if clock_row.empty:
        raise ValueError("clock_z coefficient not found in full Cox model.")

    clock_row = clock_row.iloc[0]

    coef = float(clock_row["coef"])
    hr = float(np.exp(coef))
    hr_ci_lower = float(np.exp(float(clock_row["coef lower 95%"])))
    hr_ci_upper = float(np.exp(float(clock_row["coef upper 95%"])))
    p_value = float(clock_row["p"])

    summary = pd.DataFrame([{
        "disease": disease,
        "clock_label": clock_label,
        "folder": folder,
        "modality": modality_inferred,
        "score_col": score_col,
        "field53_instance": "{}_0".format(field53_instance),
        "field53_col_from_melbourne": field53_col,
        "death_date_cols": ";".join(death_date_cols),
        "n_analysis_rows": int(full_design.shape[0]),
        "n_deaths": int(full_design[event_col].sum()),
        "event_rate": float(full_design[event_col].mean()),
        "median_followup_years": float(np.nanmedian(full_design[duration_col])),
        "clock_hr_per_1sd": hr,
        "clock_hr_ci_lower": hr_ci_lower,
        "clock_hr_ci_upper": hr_ci_upper,
        "clock_coef": coef,
        "clock_coef_se": float(clock_row["se(coef)"]),
        "clock_p": p_value,
        "cindex_covariates": cindex_cov,
        "cindex_covariates_plus_clock": cindex_full,
        "delta_cindex_clock_vs_covariates": cindex_full - cindex_cov,
        "admin_censor_date": args.admin_censor_date,
        "penalizer": args.penalizer,
        "death_xlsx": str(death_xlsx),
        "id_match_csv": str(id_match_csv),
        "covariate_csv": str(covariate_csv),
        "covariate_id_col": cov_id_col,
        "covariates_selected": ";".join(covariate_cols),
        "covariate_source_summary": covariate_source_desc,
    }])

    summary.to_csv(summary_out, sep="\t", index=False)

    qc = pd.DataFrame([
        {"metric": "score_rows", "value": score_df.shape[0]},
        {"metric": "score_nonmissing_clock_z", "value": int(score_df["clock_z"].notna().sum())},
        {"metric": "rows_with_valid_followup_before_model", "value": n_before_model},
        {"metric": "rows_with_nonmissing_clock_before_model", "value": n_with_clock},
        {"metric": "deaths_with_nonmissing_clock_before_model", "value": n_deaths_with_clock},
        {"metric": "covariate_model_rows", "value": int(cov_design.shape[0])},
        {"metric": "covariate_model_cols", "value": int(cov_design.shape[1])},
        {"metric": "full_model_rows", "value": int(full_design.shape[0])},
        {"metric": "full_model_cols", "value": int(full_design.shape[1])},
        {"metric": "field53_instance", "value": "{}_0".format(field53_instance)},
        {"metric": "field53_col_from_melbourne", "value": field53_col},
        {"metric": "death_date_cols", "value": ";".join(death_date_cols)},
        {"metric": "covariate_id_col", "value": cov_id_col},
        {"metric": "covariates_selected", "value": ";".join(covariate_cols)},
        {"metric": "covariate_source_summary", "value": covariate_source_desc},
        {"metric": "genetic_pcs_included_by_default", "value": "No"},
    ])

    qc.to_csv(qc_out, sep="\t", index=False)

    info("============================================================")
    info("Finished mortality survival analysis")
    info("Summary:")
    info("  {}".format(summary_out))
    info("Full model coefficients:")
    info("  {}".format(full_coef_out))
    info("Covariate model coefficients:")
    info("  {}".format(covar_coef_out))
    info("QC:")
    info("  {}".format(qc_out))
    info("============================================================")


# ============================================================
# 6. CLI
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Run mortality Cox model for one stable/significant disease clock."
    )

    parser.add_argument("--base-dir", required=True)
    parser.add_argument("--score-wide-tsv", required=True)
    parser.add_argument("--tasks-tsv", required=True)
    parser.add_argument("--task-index", required=True, type=int)

    parser.add_argument("--death-xlsx", required=True)
    parser.add_argument("--id-match-csv", required=True)
    parser.add_argument("--covariate-csv", required=True)
    parser.add_argument("--admin-censor-date", default="2022-11-30")

    parser.add_argument("--death-id-col", default=None)
    parser.add_argument("--death-date-col", default=None)
    parser.add_argument("--idmatch-score-col", default=None)
    parser.add_argument("--idmatch-death-col", default=None)
    parser.add_argument("--covariate-id-col", default=None)

    parser.add_argument("--field53-0-col", default=None)
    parser.add_argument("--field53-2-col", default=None)

    parser.add_argument(
        "--covariate-cols",
        default=None,
        help=(
            "Comma-separated exact covariate columns. "
            "If omitted, default covariates are: age at assessment, sex, "
            "genetic ethnic grouping, assessment center, smoking, BMI, "
            "mean diastolic BP, mean systolic BP. Genetic PCs are excluded by default."
        ),
    )

    parser.add_argument("--penalizer", type=float, default=0.01)
    parser.add_argument("--min-events", type=int, default=20)

    args = parser.parse_args()

    run_one_clock(args)


if __name__ == "__main__":
    main()