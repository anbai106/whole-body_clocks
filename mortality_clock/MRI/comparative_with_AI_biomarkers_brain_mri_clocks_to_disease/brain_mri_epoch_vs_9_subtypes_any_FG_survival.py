#!/usr/bin/env python3
"""Fair survival comparison: brain MRI mortality EPOCH vs 9 AI disease subtypes.

BUGFIX VERSION: robust nullable-boolean handling in F/G endpoint construction.

Endpoint: first recorded inpatient ICD-10 F* or G* diagnosis after the brain MRI
visit, among participants with no valid F/G diagnosis on or before that visit.
Primary analysis uses the held-out EPOCH test split and one strict common
complete-case sample for all 10 predictors.

Pipeline logic
Start from the brain MRI mortality EPOCH population. The script reads brain_mri_mortality_clock_predictions.tsv, using imaging_date as time zero and brain_mri_mortality_clock_acceleration_z as the EPOCH predictor. Importantly, the primary analysis defaults to split == "test", rather than mixing the EPOCH training participants into its evaluation. This is a stricter comparison than your previous full-sample clock-vs-BAG analysis and avoids favoring EPOCH through training-set reuse. You can later run --split all as a sensitivity analysis.
Merge the nine AI-derived disease subtype scores using id_upenn → participant_id: AD1, AD2, ASD1–3, LLD1–2, and SCZ1–2. The script then restricts all analyses to a single common complete-case population containing EPOCH plus all nine subtype scores. Thus, every reported EPOCH-versus-subtype C-index comes from exactly the same participants.
Construct a new composite incident F/G endpoint. Your existing population-preparation code uses UKBB inpatient ICD-10 field 41270 and the corresponding first-diagnosis-date field 41280. The new script scans every 41270 slot and identifies any code beginning with F or G, then pairs that code with the corresponding 41280 date. For each participant it determines the earliest recorded F/G diagnosis date. Note that G* technically represents the entire ICD-10 nervous-system chapter, not only CNS diseases, which is exactly the prefix-based definition you requested.
Define an incident, disease-free-at-MRI population. A participant is excluded if their earliest valid F/G diagnosis occurs on or before their brain MRI date. Among the remaining participants, case = 1 if the first F/G diagnosis occurs after MRI and before censoring; otherwise case = 0. This is conceptually similar to your existing disease-free logic, which distinguishes disease-free participants and disease-specific patients before constructing survival dates. I deliberately do not require participants to be free of every other ICD disease—only F/G disease—because the outcome here is first future brain/mental/nervous-system diagnosis.
Apply conservative date QC. If an F/G code exists but its paired 41280 date is missing, the participant is excluded by default because it is impossible to determine whether that diagnosis was prevalent or incident. The output cohort_qc.tsv tells you how many participants this affects. You can override this with --allow_unresolved_fg_dates, although I would keep the strict default for the manuscript.
Use the same baseline covariates for every model: age at imaging, sex, smoking status, and BMI. These are the same core adjustment variables used by your current comparative Cox script. Age is anchored to the MRI visit rather than recruitment, consistent with your previous imaging survival logic.
Standardize all ten predictors within the exact same common sample. Although EPOCH is already expressed as a z-score, it is re-standardized along with all nine subtype scores. Therefore, every HR is interpretable as the hazard ratio per 1-SD higher score in the common analysis population.
Fit the primary single-marker models: covariates, covariates + EPOCH, and covariates + each subtype separately. For every marker, the output gives N, cases, HR/95% CI/P, base C-index, marker C-index, ΔC-index over the base model, and LRT P value versus the base model.
Fit direct EPOCH-versus-subtype joint models. For each subtype, the model is covariates + EPOCH + subtype. This gives conditional HRs for both predictors, the joint C-index, correlation between EPOCH and the subtype, and two especially useful nested tests: whether EPOCH adds information beyond that subtype, and whether that subtype adds information beyond EPOCH.
Test EPOCH against all nine disease subtypes jointly. The script fits covariates + all 9 subtypes and then covariates + all 9 subtypes + EPOCH. This directly answers whether mortality EPOCH retains incremental information after accounting for the full set of disease-specific AI signatures.
Use paired bootstrap inference for discrimination. The SLURM script requests 200 bootstrap replicates by default. Each replicate re-fits the base, EPOCH, and nine subtype Cox models on the same resampled participants. It provides 95% bootstrap CIs and two-sided bootstrap P values for C-index(EPOCH) − C-index(subtype) and each marker's ΔC-index versus the base model. This is considerably stronger than simply comparing point C-indices.

"""

import argparse
import os
import warnings
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from scipy.stats import chi2

warnings.filterwarnings("ignore")

SUBTYPES = ["AD1", "AD2", "ASD1", "ASD2", "ASD3", "LLD1", "LLD2", "SCZ1", "SCZ2"]
EPOCH_RAW = "brain_mri_mortality_clock_acceleration_z"
EPOCH_NAME = "Brain_MRI_mortality_EPOCH"
BASE_COVARS = ["Age_imaging", "Sex", "Smoking", "BMI"]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--epoch_tsv", default="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/brain_mri_mortality_clock/brain_mri_mortality_clock_predictions.tsv")
    p.add_argument("--subtype_tsv", default="/cbica/projects/MULTI/processed/UKBB/derived_AI_biomakers_across_projects/UKBB_487894_participant_58_biomarker_matched_ID.tsv")
    p.add_argument("--subtype_id_col", default="id_upenn")
    p.add_argument("--icd10_csv", default="/cbica/home/wenju/Reproducibile_paper/BrainEye/data/UKBB_fullsample_ICD10.csv")
    p.add_argument("--umel_death_xlsx", default="/cbica/home/wenju/Dataset/UKBB_UMelbourne/Death_related_var_from_Ye.xlsx")
    p.add_argument("--umel_match_csv", default="/cbica/home/wenju/Dataset/UKBB_UMelbourne/UKB_UMelbourne_vs_Penn_match_key.csv")
    p.add_argument("--cov_tsv", default="/cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv")
    p.add_argument("--output_dir", default="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/brain_mri_mortality_clock/comparative_with_9_disease_subtypes_any_FG")
    p.add_argument("--split", choices=["test", "train", "all"], default="test")
    p.add_argument("--admin_censor_date", default="2022-11-30")
    p.add_argument("--min_events", type=int, default=20)
    p.add_argument("--bootstrap", type=int, default=200)
    p.add_argument("--seed", type=int, default=20260818)
    p.add_argument("--allow_unresolved_fg_dates", action="store_true")
    return p.parse_args()


def norm_id(x):
    return pd.to_numeric(x, errors="coerce").astype("Int64")


def parse_date(x):
    s = x.copy().replace([0, 0.0, "0", "0.0", "", "NA", "NaN", "nan", "None", "-1", -1], np.nan)
    out = pd.to_datetime(s, errors="coerce")
    num = pd.to_numeric(s, errors="coerce")
    mask = num.between(20000, 60000)
    if mask.any():
        excel = pd.to_datetime(num, unit="D", origin="1899-12-30", errors="coerce")
        out = out.where(~mask, excel)
    return out


def bool_numpy(mask):
    """Return a strict NumPy bool array from pandas/NumPy nullable masks.

    Pandas nullable boolean/string operations can produce object-dtype arrays
    after to_numpy(). Missing values are treated as False so NumPy in-place
    logical operations such as |= are always type-safe.
    """
    if isinstance(mask, pd.Series):
        return mask.fillna(False).astype(bool).to_numpy(dtype=bool)

    arr = np.asarray(mask)
    if arr.dtype == bool:
        return arr

    return (
        pd.Series(arr)
        .fillna(False)
        .astype(bool)
        .to_numpy(dtype=bool)
    )


def read_epoch(args):
    d = pd.read_csv(args.epoch_tsv, sep="\t")
    req = ["participant_id", "imaging_date", "age_at_imaging", EPOCH_RAW]
    miss = [c for c in req if c not in d.columns]
    if miss:
        raise ValueError(f"Missing EPOCH columns: {miss}")
    qc = {"epoch_rows_raw": len(d)}
    if args.split != "all":
        if "split" not in d.columns:
            raise ValueError("EPOCH predictions file has no split column")
        d = d[d["split"].astype(str).str.lower().eq(args.split)].copy()
    qc["epoch_rows_after_split"] = len(d)
    keep = req + [c for c in ["split", "sex", "death_date", "admin_censor_date", "end_date"] if c in d.columns]
    d = d[keep].copy()
    d["participant_id"] = norm_id(d["participant_id"])
    d = d[d["participant_id"].notna()].copy()
    d["imaging_date"] = parse_date(d["imaging_date"])
    d["Age_imaging"] = pd.to_numeric(d["age_at_imaging"], errors="coerce")
    d[EPOCH_RAW] = pd.to_numeric(d[EPOCH_RAW], errors="coerce")
    admin_default = pd.Timestamp(args.admin_censor_date)
    d["admin_censor_date"] = parse_date(d["admin_censor_date"]) if "admin_censor_date" in d else admin_default
    if not isinstance(d["admin_censor_date"], pd.Series):
        d["admin_censor_date"] = admin_default
    d["admin_censor_date"] = d["admin_censor_date"].fillna(admin_default)
    d["death_date"] = parse_date(d["death_date"]) if "death_date" in d else pd.NaT
    fallback_end = d["admin_censor_date"].copy()
    death_mask = d["death_date"].notna() & (d["death_date"] < fallback_end)
    fallback_end.loc[death_mask] = d.loc[death_mask, "death_date"]
    if "end_date" in d.columns:
        d["censor_date"] = parse_date(d["end_date"]).fillna(fallback_end)
    else:
        d["censor_date"] = fallback_end
    if "sex" in d.columns:
        s = d["sex"].astype(str).str.strip().str.lower()
        d["Sex_epoch_fallback"] = np.where(s.isin(["male", "m", "1", "1.0"]), 1.0, np.where(s.isin(["female", "f", "0", "0.0"]), 0.0, np.nan))
    else:
        d["Sex_epoch_fallback"] = np.nan
    d = d[d["imaging_date"].notna() & d["censor_date"].notna() & (d["censor_date"] > d["imaging_date"])].copy()
    qc["epoch_duplicate_rows_before_dedup"] = int(d["participant_id"].duplicated(keep=False).sum())
    d = d.sort_values(["participant_id", "imaging_date"]).drop_duplicates("participant_id", keep="first")
    qc["epoch_unique_ids_usable"] = d["participant_id"].nunique()
    return d, qc


def read_subtypes(args, epoch_ids):
    d = pd.read_csv(args.subtype_tsv, sep="\t")
    if args.subtype_id_col not in d.columns:
        raise ValueError(f"Missing subtype ID column: {args.subtype_id_col}")
    miss = [c for c in SUBTYPES if c not in d.columns]
    if miss:
        raise ValueError(f"Missing subtype columns: {miss}")
    d = d[[args.subtype_id_col] + SUBTYPES].rename(columns={args.subtype_id_col: "participant_id"}).copy()
    d["participant_id"] = norm_id(d["participant_id"])
    d = d[d["participant_id"].notna()].copy()
    for c in SUBTYPES:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    qc = {"subtype_rows_raw": len(d), "subtype_duplicate_rows_before_dedup": int(d["participant_id"].duplicated(keep=False).sum())}
    d["_n"] = d[SUBTYPES].notna().sum(axis=1)
    d = d.sort_values(["participant_id", "_n"], ascending=[True, False]).drop_duplicates("participant_id").drop(columns="_n")
    qc["subtype_unique_ids"] = d["participant_id"].nunique()
    qc["subtype_overlap_with_epoch"] = int(d["participant_id"].isin(epoch_ids).sum())
    return d, qc


def read_covariates(args):
    d = pd.read_csv(args.cov_tsv)
    if "eid" not in d.columns:
        raise ValueError("Covariate file must contain eid")
    sex = next((c for c in ["sex_f31_0_0", "genetic_sex_f22001_0_0", "Sex", "sex"] if c in d.columns), None)
    smoke = next((c for c in ["smoking_status_f20116_0_0", "Smoking", "smoking"] if c in d.columns), None)
    bmi = next((c for c in ["body_mass_index_bmi_f23104_0_0", "BMI", "bmi"] if c in d.columns), None)
    keep = ["eid"] + [c for c in [sex, smoke, bmi] if c is not None]
    d = d[keep].copy().rename(columns={"eid": "participant_id", **({sex: "Sex"} if sex else {}), **({smoke: "Smoking"} if smoke else {}), **({bmi: "BMI"} if bmi else {})})
    d["participant_id"] = norm_id(d["participant_id"])
    d = d[d["participant_id"].notna()].copy()
    for c in ["Sex", "Smoking", "BMI"]:
        if c not in d:
            d[c] = np.nan
        d[c] = pd.to_numeric(d[c], errors="coerce")
    return d[["participant_id", "Sex", "Smoking", "BMI"]].drop_duplicates("participant_id")


def read_icd_codes(args, ids):
    d = pd.read_csv(args.icd10_csv)
    if "eid" not in d.columns:
        raise ValueError("ICD file must contain eid")
    code_cols = [c for c in d.columns if c.startswith("diagnoses_icd10_f41270_")]
    if not code_cols:
        raise ValueError("No diagnoses_icd10_f41270_* columns found")
    d = d[["eid"] + code_cols].rename(columns={"eid": "participant_id"}).copy()
    d["participant_id"] = norm_id(d["participant_id"])
    d = d[d["participant_id"].isin(ids)].drop_duplicates("participant_id")
    return d, code_cols


def read_icd_dates(args, ids):
    death = pd.read_excel(args.umel_death_xlsx, engine="openpyxl")
    match = pd.read_csv(args.umel_match_csv)
    if "eid" not in death.columns or not {"id", "id_upenn"}.issubset(match.columns):
        raise ValueError("UMelbourne diagnosis-date inputs have unexpected ID columns")
    date_cols = [c for c in death.columns if c.startswith("41280-")]
    if not date_cols:
        raise ValueError("No 41280-* diagnosis date columns found")
    death = death[["eid"] + date_cols].rename(columns={"eid": "participant_id_umel"})
    match = match[["id", "id_upenn"]].rename(columns={"id": "participant_id_umel", "id_upenn": "participant_id"})
    death["participant_id_umel"] = norm_id(death["participant_id_umel"])
    match["participant_id_umel"] = norm_id(match["participant_id_umel"])
    match["participant_id"] = norm_id(match["participant_id"])
    match = match[match["participant_id"].isin(ids)]
    d = match.merge(death, on="participant_id_umel", how="inner")[["participant_id"] + date_cols]
    d = d.rename(columns={c: c.replace("41280-", "").replace(".", "_") for c in date_cols})
    return d.drop_duplicates("participant_id")


def derive_fg_endpoint(base, diag, code_cols, dates, allow_unresolved):
    """Construct the first incident ICD-10 F/G endpoint after brain MRI.

    Bug fix:
    pandas string operations may return nullable BooleanDtype masks. Plain
    .to_numpy() can then yield dtype=object, which causes NumPy UFuncTypeError
    during in-place operations against dtype=bool arrays. Every mask used in
    NumPy logical operations is explicitly converted with bool_numpy().
    """
    x = (
        base
        .merge(diag, on="participant_id", how="left")
        .merge(dates, on="participant_id", how="left")
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

    fg_count = np.zeros(len(x), dtype=np.int32)
    unresolved = np.zeros(len(x), dtype=bool)

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

        # Critical bug fix: always use a true NumPy bool array here.
        is_fg_np = bool_numpy(
            codes.str.match(
                r"^[FG][0-9]",
                na=False,
            )
        )

        if not is_fg_np.any():
            continue

        fg_count += is_fg_np.astype(np.int32)

        # A code without its paired date cannot be classified as prevalent
        # versus incident and is therefore unresolved under the strict default.
        if suffix not in x.columns:
            unresolved |= is_fg_np
            unmatched_slots += 1
            continue

        matched_slots += 1
        dt = parse_date(x[suffix])

        date_missing_np = bool_numpy(dt.isna())
        date_present_np = bool_numpy(dt.notna())

        unresolved |= (
            is_fg_np
            & date_missing_np
        )

        valid_np = (
            is_fg_np
            & date_present_np
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
            | (dt < earliest)
        )

        earlier = pd.Series(
            bool_numpy(earlier),
            index=x.index,
            dtype=bool,
        )

        earliest.loc[earlier] = dt.loc[earlier]
        earliest_code.loc[earlier] = codes.loc[earlier]

    out = base.copy()

    out["earliest_FG_date"] = earliest.values
    out["earliest_FG_code"] = earliest_code.values
    out["n_FG_codes_recorded"] = fg_count
    out["has_FG_code_missing_paired_date"] = unresolved

    out["prevalent_FG_at_imaging"] = (
        out["earliest_FG_date"].notna()
        & (
            out["earliest_FG_date"]
            <= out["imaging_date"]
        )
    )

    if allow_unresolved:
        out["exclude_unresolved_FG_date"] = False
    else:
        out["exclude_unresolved_FG_date"] = (
            out["has_FG_code_missing_paired_date"]
            .astype(bool)
        )

    eligible = (
        ~out["prevalent_FG_at_imaging"].astype(bool)
        & ~out["exclude_unresolved_FG_date"].astype(bool)
    )

    out["case"] = (
        eligible
        & out["earliest_FG_date"].notna()
        & (
            out["earliest_FG_date"]
            > out["imaging_date"]
        )
        & (
            out["earliest_FG_date"]
            <= out["censor_date"]
        )
    ).astype(int)

    out["analysis_end_date"] = out["censor_date"]

    case_mask = out["case"].eq(1)
    out.loc[
        case_mask,
        "analysis_end_date",
    ] = out.loc[
        case_mask,
        "earliest_FG_date",
    ]

    out["time_years"] = (
        (
            out["analysis_end_date"]
            - out["imaging_date"]
        ).dt.days
        / 365.25
    )

    qc = {
        "icd_41270_slots": len(code_cols),
        "icd_slots_with_matching_41280": matched_slots,
        "icd_slots_without_matching_41280": unmatched_slots,
        "participants_with_any_FG_code": int(
            (
                out["n_FG_codes_recorded"]
                > 0
            ).sum()
        ),
        "participants_with_unresolved_FG_date": int(
            out["has_FG_code_missing_paired_date"].sum()
        ),
        "participants_prevalent_FG_at_imaging": int(
            out["prevalent_FG_at_imaging"].sum()
        ),
        "incident_FG_events_before_endpoint_filter": int(
            out["case"].sum()
        ),
    }

    time_numeric = pd.to_numeric(
        out["time_years"],
        errors="coerce",
    )

    time_ok = (
        np.isfinite(time_numeric)
        & (time_numeric > 0)
    )

    out = out[
        eligible
        & time_ok
    ].copy()

    qc["endpoint_eligible_N"] = len(out)
    qc["endpoint_eligible_events"] = int(
        out["case"].sum()
    )

    return out, qc


def standardize(df, src, dst):
    v = pd.to_numeric(df[src], errors="coerce")
    mu, sd = float(v.mean()), float(v.std(ddof=1))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError(f"Cannot standardize {src}; SD={sd}")
    df[dst] = (v - mu) / sd
    return mu, sd


def fit_cox(df, predictors, label):
    cols = ["time_years", "case"] + BASE_COVARS + predictors
    x = df[cols].copy()
    for c in cols:
        x[c] = pd.to_numeric(x[c], errors="coerce")
    x = x.replace([np.inf, -np.inf], np.nan).dropna()
    if x["case"].nunique() < 2:
        raise ValueError(f"{label}: need cases and noncases")
    last = None
    for pen in [0.0, 0.001, 0.01, 0.1]:
        try:
            cph = CoxPHFitter(penalizer=pen)
            cph.fit(x, duration_col="time_years", event_col="case", show_progress=False)
            return cph, x, pen
        except Exception as e:
            last = e
    raise RuntimeError(f"{label} failed: {last}")


def get_cindex(cph, x):
    risk = cph.predict_partial_hazard(x).values.ravel()
    return float(concordance_index(x["time_years"], -risk, x["case"]))


def hr_stats(cph, var):
    b = float(cph.params_.loc[var]); se = float(cph.standard_errors_.loc[var]); p = float(cph.summary.loc[var, "p"])
    return {"beta": b, "se": se, "hr": float(np.exp(b)), "lo": float(np.exp(b - 1.96*se)), "hi": float(np.exp(b + 1.96*se)), "p": p}


def lrt(full, reduced, df_diff, pen_full, pen_reduced):
    if pen_full != 0.0 or pen_reduced != 0.0:
        return np.nan, np.nan
    stat = 2.0 * (float(full.log_likelihood_) - float(reduced.log_likelihood_))
    if not np.isfinite(stat) or stat < 0:
        return np.nan, np.nan
    return stat, float(chi2.sf(stat, df_diff))


def build_common_sample(endpoint, subtypes, cov):
    ep = endpoint[[c for c in endpoint.columns if c not in SUBTYPES]].copy()
    d = ep.merge(subtypes, on="participant_id", how="inner").merge(cov, on="participant_id", how="left")
    if "Sex_epoch_fallback" in d.columns:
        d["Sex"] = d["Sex"].where(d["Sex"].notna(), d["Sex_epoch_fallback"])
    needed = ["participant_id", "time_years", "case"] + BASE_COVARS + [EPOCH_RAW] + SUBTYPES
    before = len(d)
    for c in needed:
        if c != "participant_id":
            d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d.replace([np.inf, -np.inf], np.nan).dropna(subset=needed).copy()
    std_rows = []
    for raw in [EPOCH_RAW] + SUBTYPES:
        z = "EPOCH_z" if raw == EPOCH_RAW else raw + "_z"
        mu, sd = standardize(d, raw, z)
        std_rows.append({"marker": EPOCH_NAME if raw == EPOCH_RAW else raw, "raw_column": raw, "z_column": z, "mean": mu, "sd": sd})
    flow = pd.DataFrame([
        {"stage": "Endpoint-eligible brain MRI population", "N": len(endpoint), "events": int(endpoint["case"].sum())},
        {"stage": "Before strict common complete-case filter", "N": before, "events": np.nan},
        {"stage": "Strict common complete-case sample", "N": len(d), "events": int(d["case"].sum())},
    ])
    return d, pd.DataFrame(std_rows), flow


def run_models(d):
    base, xb, pbase = fit_cox(d, [], "base")
    base_c = get_cindex(base, xb)
    models = {"BASE": {"fit": base, "x": xb, "pen": pbase, "c": base_c}}
    rows = []
    specs = [(EPOCH_NAME, "EPOCH_z")] + [(s, s + "_z") for s in SUBTYPES]
    for name, z in specs:
        fit, x, pen = fit_cox(d, [z], name)
        c = get_cindex(fit, x); h = hr_stats(fit, z); stat, p = lrt(fit, base, 1, pen, pbase)
        rows.append({"marker": name, "N": len(x), "N_case": int(x["case"].sum()), "HR_per_1SD": h["hr"], "HR_CI_low": h["lo"], "HR_CI_high": h["hi"], "beta": h["beta"], "se": h["se"], "p_marker": h["p"], "base_cindex": base_c, "marker_cindex": c, "delta_cindex_vs_base": c-base_c, "lrt_chisq_vs_base": stat, "lrt_p_vs_base": p, "penalizer": pen})
        models[name] = {"fit": fit, "x": x, "pen": pen, "c": c}
    return pd.DataFrame(rows), models


def pairwise_models(d, models):
    rows = []; e = models[EPOCH_NAME]
    for s in SUBTYPES:
        m = models[s]
        fit, x, pen = fit_cox(d, ["EPOCH_z", s + "_z"], f"EPOCH+{s}")
        c = get_cindex(fit, x); he = hr_stats(fit, "EPOCH_z"); hs = hr_stats(fit, s + "_z")
        st_e, p_e = lrt(fit, m["fit"], 1, pen, m["pen"])
        st_s, p_s = lrt(fit, e["fit"], 1, pen, e["pen"])
        r = float(d[["EPOCH_z", s + "_z"]].corr().iloc[0,1])
        rows.append({"subtype": s, "N": len(x), "N_case": int(x["case"].sum()), "pearson_EPOCH_subtype": r, "EPOCH_standalone_cindex": e["c"], "subtype_standalone_cindex": m["c"], "delta_cindex_EPOCH_minus_subtype": e["c"]-m["c"], "joint_cindex": c, "delta_cindex_joint_minus_EPOCH": c-e["c"], "delta_cindex_joint_minus_subtype": c-m["c"], "EPOCH_conditional_HR": he["hr"], "EPOCH_conditional_CI_low": he["lo"], "EPOCH_conditional_CI_high": he["hi"], "EPOCH_conditional_p": he["p"], "subtype_conditional_HR": hs["hr"], "subtype_conditional_CI_low": hs["lo"], "subtype_conditional_CI_high": hs["hi"], "subtype_conditional_p": hs["p"], "lrt_chisq_EPOCH_beyond_subtype": st_e, "lrt_p_EPOCH_beyond_subtype": p_e, "lrt_chisq_subtype_beyond_EPOCH": st_s, "lrt_p_subtype_beyond_EPOCH": p_s, "joint_penalizer": pen})
    return pd.DataFrame(rows)


def combined_models(d, models):
    zs = [s + "_z" for s in SUBTYPES]
    f9, x9, p9 = fit_cox(d, zs, "all9")
    f10, x10, p10 = fit_cox(d, zs + ["EPOCH_z"], "all9+EPOCH")
    c9, c10 = get_cindex(f9, x9), get_cindex(f10, x10)
    stat, p = lrt(f10, f9, 1, p10, p9); h = hr_stats(f10, "EPOCH_z")
    ce = models[EPOCH_NAME]["c"]
    return pd.DataFrame([{ "N": len(x10), "N_case": int(x10["case"].sum()), "base_cindex": models["BASE"]["c"], "EPOCH_only_cindex": ce, "all_9_subtypes_cindex": c9, "all_9_subtypes_plus_EPOCH_cindex": c10, "delta_cindex_EPOCH_only_minus_all9": ce-c9, "delta_cindex_EPOCH_beyond_all9": c10-c9, "delta_cindex_all9_beyond_EPOCH": c10-ce, "lrt_chisq_EPOCH_beyond_all9": stat, "lrt_p_EPOCH_beyond_all9": p, "EPOCH_conditional_HR_given_all9": h["hr"], "EPOCH_conditional_CI_low": h["lo"], "EPOCH_conditional_CI_high": h["hi"], "EPOCH_conditional_p_given_all9": h["p"], "penalizer_all9": p9, "penalizer_all10": p10 }])


def bootstrap_summary(d, marker_tbl, B, seed):
    if B <= 0:
        return pd.DataFrame()
    obs = dict(zip(marker_tbl["marker"], marker_tbl["marker_cindex"])); base_obs = float(marker_tbl["base_cindex"].iloc[0])
    rng = np.random.default_rng(seed); n = len(d); rec = []
    print(f"Starting paired bootstrap: B={B}", flush=True)
    for b in range(B):
        boot = d.iloc[rng.integers(0, n, size=n)].reset_index(drop=True)
        if boot["case"].nunique() < 2:
            continue
        try:
            fb, xb, _ = fit_cox(boot, [], f"boot{b} base"); cb = get_cindex(fb, xb)
            fe, xe, _ = fit_cox(boot, ["EPOCH_z"], f"boot{b} EPOCH"); ce = get_cindex(fe, xe)
            row = {"EPOCH_vs_base": ce-cb}
            for s in SUBTYPES:
                fs, xs, _ = fit_cox(boot, [s + "_z"], f"boot{b} {s}"); cs = get_cindex(fs, xs)
                row["EPOCH_minus_" + s] = ce-cs
                row[s + "_vs_base"] = cs-cb
            rec.append(row)
        except Exception:
            continue
        if (b+1) % 10 == 0 or b == B-1:
            print(f"  {b+1}/{B}, successful={len(rec)}", flush=True)
    br = pd.DataFrame(rec)
    def summarize(name, vals, observed):
        v = pd.to_numeric(vals, errors="coerce").dropna().to_numpy()
        if len(v) < 20:
            return {"comparison": name, "observed_difference": observed, "bootstrap_CI_low": np.nan, "bootstrap_CI_high": np.nan, "bootstrap_p_two_sided": np.nan, "successful_bootstrap_replicates": len(v)}
        lo, hi = np.quantile(v, [0.025, 0.975]); p0 = (np.sum(v <= 0)+1)/(len(v)+1); p1 = (np.sum(v >= 0)+1)/(len(v)+1)
        return {"comparison": name, "observed_difference": observed, "bootstrap_CI_low": float(lo), "bootstrap_CI_high": float(hi), "bootstrap_p_two_sided": float(min(1, 2*min(p0,p1))), "successful_bootstrap_replicates": len(v)}
    if br.empty:
        return pd.DataFrame([summarize("bootstrap_failed", pd.Series(dtype=float), np.nan)])
    out = [summarize("EPOCH vs base", br["EPOCH_vs_base"], obs[EPOCH_NAME]-base_obs)]
    for s in SUBTYPES:
        out.append(summarize("EPOCH - " + s, br["EPOCH_minus_"+s], obs[EPOCH_NAME]-obs[s]))
        out.append(summarize(s + " vs base", br[s+"_vs_base"], obs[s]-base_obs))
    return pd.DataFrame(out)


def print_table(title, df):
    print("\n" + "="*110); print(title); print("="*110)
    with pd.option_context("display.max_columns", None, "display.width", 240, "display.max_rows", 200):
        print(df.to_string(index=False) if len(df) else "<empty>")


def main():
    args = parse_args(); os.makedirs(args.output_dir, exist_ok=True)
    qc_rows = []
    epoch, q1 = read_epoch(args); qc_rows += [{"metric": k, "value": v} for k,v in q1.items()]
    sub, q2 = read_subtypes(args, set(epoch["participant_id"].dropna().astype(int))); qc_rows += [{"metric": k, "value": v} for k,v in q2.items()]
    base = epoch.merge(sub, on="participant_id", how="inner")
    qc_rows.append({"metric": "epoch_subtype_candidate_overlap_N", "value": len(base)})
    ids = set(base["participant_id"].dropna().astype(int))
    cov = read_covariates(args); diag, code_cols = read_icd_codes(args, ids); dates = read_icd_dates(args, ids)
    endpoint, q3 = derive_fg_endpoint(base, diag, code_cols, dates, args.allow_unresolved_fg_dates); qc_rows += [{"metric": k, "value": v} for k,v in q3.items()]
    analysis, std_tbl, flow = build_common_sample(endpoint, sub, cov)
    if int(analysis["case"].sum()) < args.min_events:
        raise RuntimeError(f"Only {int(analysis['case'].sum())} events remain; min_events={args.min_events}")
    marker_tbl, models = run_models(analysis)
    pair_tbl = pairwise_models(analysis, models)
    combined_tbl = combined_models(analysis, models)
    boot_tbl = bootstrap_summary(analysis, marker_tbl, args.bootstrap, args.seed)
    pd.DataFrame(qc_rows).to_csv(os.path.join(args.output_dir, "cohort_qc.tsv"), sep="\t", index=False)
    flow.to_csv(os.path.join(args.output_dir, "sample_flow.tsv"), sep="\t", index=False)
    std_tbl.to_csv(os.path.join(args.output_dir, "marker_standardization.tsv"), sep="\t", index=False)
    keep = ["participant_id", "imaging_date", "censor_date", "earliest_FG_date", "earliest_FG_code", "n_FG_codes_recorded", "case", "time_years"] + BASE_COVARS + [EPOCH_RAW, "EPOCH_z"] + SUBTYPES + [s+"_z" for s in SUBTYPES] + (["split"] if "split" in analysis.columns else [])
    analysis[[c for c in keep if c in analysis.columns]].to_csv(os.path.join(args.output_dir, "analysis_common_sample.tsv"), sep="\t", index=False)
    marker_tbl.to_csv(os.path.join(args.output_dir, "marker_models_common_sample.tsv"), sep="\t", index=False)
    pair_tbl.to_csv(os.path.join(args.output_dir, "epoch_vs_subtype_pairwise.tsv"), sep="\t", index=False)
    combined_tbl.to_csv(os.path.join(args.output_dir, "combined_subtypes_vs_epoch.tsv"), sep="\t", index=False)
    boot_tbl.to_csv(os.path.join(args.output_dir, "bootstrap_cindex_common_sample.tsv"), sep="\t", index=False)
    with open(os.path.join(args.output_dir, "run_metadata.txt"), "w") as f:
        f.write(f"split={args.split}\nendpoint=first recorded inpatient ICD-10 F* or G* diagnosis after MRI\n")
        f.write(f"strict_common_complete_case=True\nbase_covariates={','.join(BASE_COVARS)}\n")
        f.write(f"analysis_N={len(analysis)}\nanalysis_events={int(analysis['case'].sum())}\nbootstrap={args.bootstrap}\nseed={args.seed}\n")
    print_table("SAMPLE FLOW", flow)
    print_table("MARKER STANDARDIZATION", std_tbl)
    print_table("PRIMARY SINGLE-MARKER COX MODELS", marker_tbl)
    print_table("PAIRWISE EPOCH VS EACH SUBTYPE", pair_tbl)
    print_table("ALL NINE SUBTYPES VS EPOCH", combined_tbl)
    print_table("PAIRED BOOTSTRAP C-INDEX COMPARISONS", boot_tbl)
    print("\nAnalysis complete. Output directory:", args.output_dir)


if __name__ == "__main__":
    main()