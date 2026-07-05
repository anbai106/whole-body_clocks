#!/usr/bin/env python3
# ============================================================
# Linear models:
#   delta_clock_age_years ~ medication_cluster +
#     baseline_accel_years +
#     chronological_age_0_0 +
#     Sex + Smoking + BMI + Diastolic + Systolic +
#     delta_chrono_age_1_minus_0
#
# Input:
#   metabolomics_delta_clock_medication_cluster_requested5_long.tsv
#
# Output:
#   medication_cluster_delta_clock_lm_cluster_effects.tsv
#   medication_cluster_delta_clock_lm_adjusted_means.tsv
#   medication_cluster_delta_clock_lm_cluster_counts.tsv
#   medication_cluster_delta_clock_lm_model_summaries.tsv
#   medication_cluster_delta_clock_lm_complete_cases.tsv
# ============================================================

import argparse
import os
import re
import json
import warnings

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")


CLUSTER_LEVELS = [
    "No/minimal medication",
    "Cardiometabolic medication cluster",
    "Respiratory medication cluster",
    "Psychiatric/pain medication cluster",
    "High polypharmacy cluster",
]

REFERENCE_CLUSTER = "No/minimal medication"

ORGAN_ORDER = ["Endocrine", "Digestive", "Hepatic", "Immune"]

DEFAULT_COVARIATES = [
    "baseline_accel_years",
    "chronological_age_0_0",
    "Sex",
    "Smoking",
    "BMI",
    "Diastolic",
    "Systolic",
    "delta_chrono_age_1_minus_0",
]


def parse_args():
    p = argparse.ArgumentParser(
        description="Run medication-cluster association models for metabolomics delta mortality clocks."
    )

    p.add_argument(
        "--input_tsv",
        default=(
            "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/"
            "mortality_clock/SA/delta_metabolomics_algorithmic_disease_onset/"
            "medication_cluster_delta_clock_inputs/"
            "metabolomics_delta_clock_medication_cluster_requested5_long.tsv"
        ),
        help="Long-format medication-cluster + delta-clock TSV."
    )

    p.add_argument(
        "--out_dir",
        default=(
            "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/"
            "mortality_clock/SA/delta_metabolomics_algorithmic_disease_onset/"
            "medication_cluster_delta_clock_lm_results"
        ),
        help="Output directory."
    )

    p.add_argument(
        "--outcome",
        default="delta_clock_age_years",
        choices=["delta_clock_age_years", "delta_accel_years"],
        help="Outcome to model. Default is delta_clock_age_years."
    )

    p.add_argument(
        "--reference_cluster",
        default=REFERENCE_CLUSTER,
        help="Reference medication cluster."
    )

    p.add_argument(
        "--min_n_cluster",
        default=20,
        type=int,
        help="Minimum N for reporting a cluster contrast. Model still fits if possible."
    )

    p.add_argument(
        "--robust_cov",
        default="HC3",
        choices=["HC0", "HC1", "HC2", "HC3", "nonrobust"],
        help="Robust covariance estimator for OLS. Default HC3."
    )

    p.add_argument(
        "--save_complete_cases",
        action="store_true",
        help="Save complete-case rows used in each organ-specific model."
    )

    return p.parse_args()


def safe_numeric(s):
    return pd.to_numeric(s, errors="coerce")


def repair_common_header_issues(df):
    """
    This handles a possible pasted/header typo:
      baseline_accel_yearsfollowup_accel_years
    If that appears, rename it to baseline_accel_years because the model only
    needs baseline_accel_years.
    """
    if (
        "baseline_accel_yearsfollowup_accel_years" in df.columns
        and "baseline_accel_years" not in df.columns
    ):
        df = df.rename(
            columns={
                "baseline_accel_yearsfollowup_accel_years": "baseline_accel_years"
            }
        )

    # Trim accidental spaces in column names.
    df = df.rename(columns={c: str(c).strip() for c in df.columns})
    return df


def read_input(path):
    if not os.path.exists(path):
        raise FileNotFoundError(path)

    print("Reading input TSV:")
    print("  ", path)

    dat = pd.read_csv(path, sep="\t", dtype={"participant_id": str}, low_memory=False)
    dat = repair_common_header_issues(dat)

    required = [
        "participant_id",
        "medication_cluster",
        "organ_label",
        "organ_clean",
        "baseline_accel_years",
        "chronological_age_0_0",
        "Sex",
        "Smoking",
        "BMI",
        "Diastolic",
        "Systolic",
        "delta_chrono_age_1_minus_0",
    ]

    missing = [c for c in required if c not in dat.columns]
    if missing:
        raise ValueError(
            "Input file is missing required columns:\n"
            + "\n".join(missing)
            + "\n\nAvailable columns:\n"
            + "\n".join(dat.columns)
        )

    return dat


def clean_analysis_data(dat, outcome):
    dat = dat.copy()

    if outcome not in dat.columns:
        raise ValueError(f"Outcome column not found: {outcome}")

    # Keep the five requested medication clusters only.
    dat = dat[dat["medication_cluster"].isin(CLUSTER_LEVELS)].copy()

    # Categorical ordering.
    dat["medication_cluster"] = pd.Categorical(
        dat["medication_cluster"],
        categories=CLUSTER_LEVELS,
        ordered=True
    )

    dat["organ_label"] = pd.Categorical(
        dat["organ_label"],
        categories=ORGAN_ORDER,
        ordered=True
    )

    numeric_cols = list(set(DEFAULT_COVARIATES + [outcome]))

    for c in numeric_cols:
        if c in dat.columns:
            dat[c] = safe_numeric(dat[c])

    # UKB missing codes occasionally appear as negative values.
    for c in ["Sex", "Smoking"]:
        if c in dat.columns:
            dat.loc[dat[c].isin([-1, -3]), c] = np.nan

    # Exclude impossible or uninformative values only lightly.
    dat = dat.replace([np.inf, -np.inf], np.nan)

    return dat


def available_covariates(df, requested_covariates):
    covars = []
    for c in requested_covariates:
        if c not in df.columns:
            continue
        x = safe_numeric(df[c])
        if x.notna().sum() > 0 and x.nunique(dropna=True) > 1:
            covars.append(c)
    return covars


def build_formula(outcome, reference_cluster, covars):
    covar_text = " + ".join(covars)
    formula = (
        f"{outcome} ~ "
        f"C(medication_cluster, Treatment(reference='{reference_cluster}'))"
    )
    if covar_text:
        formula += " + " + covar_text
    return formula


def term_name_for_cluster(reference_cluster, cluster):
    return (
        f"C(medication_cluster, Treatment(reference='{reference_cluster}'))"
        f"[T.{cluster}]"
    )


def fit_one_organ(dat, organ_label, args):
    df = dat[dat["organ_label"] == organ_label].copy()

    if df.empty:
        return None, pd.DataFrame(), pd.DataFrame(), {
            "organ_label": organ_label,
            "status": "no_rows",
            "error": "No rows for organ."
        }

    covars = available_covariates(df, DEFAULT_COVARIATES)

    needed = [
        args.outcome,
        "medication_cluster",
        "organ_label",
        "organ_clean",
        "participant_id",
    ] + covars

    df = df[needed].copy()

    for c in covars + [args.outcome]:
        df[c] = safe_numeric(df[c])

    df = df.dropna(subset=[args.outcome, "medication_cluster"] + covars).copy()

    # Drop unused categories after complete-case filtering.
    df["medication_cluster"] = pd.Categorical(
        df["medication_cluster"],
        categories=CLUSTER_LEVELS,
        ordered=True
    )

    cluster_counts = (
        df.groupby("medication_cluster", observed=False)
        .size()
        .reset_index(name="N")
    )
    cluster_counts["organ_label"] = organ_label

    n_clusters_present = int((cluster_counts["N"] > 0).sum())
    n_total = int(df.shape[0])

    if n_total < 50:
        summary = {
            "organ_label": organ_label,
            "status": "insufficient_total_n",
            "error": f"N={n_total} after complete-case filtering.",
            "N": n_total,
            "n_clusters_present": n_clusters_present,
        }
        return None, cluster_counts, df, summary

    if n_clusters_present < 2:
        summary = {
            "organ_label": organ_label,
            "status": "insufficient_clusters",
            "error": "Fewer than 2 medication clusters present after complete-case filtering.",
            "N": n_total,
            "n_clusters_present": n_clusters_present,
        }
        return None, cluster_counts, df, summary

    formula = build_formula(args.outcome, args.reference_cluster, covars)

    try:
        model = smf.ols(formula=formula, data=df)

        if args.robust_cov == "nonrobust":
            fit = model.fit()
        else:
            fit = model.fit(cov_type=args.robust_cov)

        summary = {
            "organ_label": organ_label,
            "status": "ok",
            "error": "",
            "N": n_total,
            "n_clusters_present": n_clusters_present,
            "formula": formula,
            "covariates_used": ",".join(covars),
            "r_squared": float(fit.rsquared),
            "adj_r_squared": float(fit.rsquared_adj),
            "aic": float(fit.aic),
            "bic": float(fit.bic),
            "robust_cov": args.robust_cov,
        }

        return fit, cluster_counts, df, summary

    except Exception as e:
        summary = {
            "organ_label": organ_label,
            "status": "fit_failed",
            "error": str(e),
            "N": n_total,
            "n_clusters_present": n_clusters_present,
            "formula": formula,
            "covariates_used": ",".join(covars),
            "robust_cov": args.robust_cov,
        }
        return None, cluster_counts, df, summary


def extract_cluster_effects(fit, cluster_counts, complete_df, organ_label, args):
    rows = []

    if fit is None:
        return pd.DataFrame()

    n_total = int(complete_df.shape[0])

    n_ref = int(
        cluster_counts.loc[
            cluster_counts["medication_cluster"] == args.reference_cluster,
            "N"
        ].sum()
    )

    for cluster in CLUSTER_LEVELS:
        if cluster == args.reference_cluster:
            continue

        n_cluster = int(
            cluster_counts.loc[
                cluster_counts["medication_cluster"] == cluster,
                "N"
            ].sum()
        )

        term = term_name_for_cluster(args.reference_cluster, cluster)

        row = {
            "organ_label": organ_label,
            "outcome": args.outcome,
            "exposure_cluster": cluster,
            "reference_cluster": args.reference_cluster,
            "N": n_total,
            "N_reference": n_ref,
            "N_exposure": n_cluster,
            "min_n_cluster": args.min_n_cluster,
            "status": "ok",
            "term": term,
        }

        if n_cluster < args.min_n_cluster:
            row.update({
                "status": "small_cluster_n",
                "beta": np.nan,
                "se": np.nan,
                "ci_lo": np.nan,
                "ci_hi": np.nan,
                "p": np.nan,
            })
        elif term not in fit.params.index:
            row.update({
                "status": "term_not_in_model",
                "beta": np.nan,
                "se": np.nan,
                "ci_lo": np.nan,
                "ci_hi": np.nan,
                "p": np.nan,
            })
        else:
            beta = float(fit.params.loc[term])
            se = float(fit.bse.loc[term])
            p = float(fit.pvalues.loc[term])
            row.update({
                "beta": beta,
                "se": se,
                "ci_lo": beta - 1.96 * se,
                "ci_hi": beta + 1.96 * se,
                "p": p,
            })

        rows.append(row)

    return pd.DataFrame(rows)


def adjusted_means(fit, cluster_counts, complete_df, organ_label, args):
    if fit is None:
        return pd.DataFrame()

    covars = available_covariates(complete_df, DEFAULT_COVARIATES)

    mean_covars = {}
    for c in covars:
        mean_covars[c] = float(pd.to_numeric(complete_df[c], errors="coerce").mean())

    pred_rows = []

    for cluster in CLUSTER_LEVELS:
        n_cluster = int(
            cluster_counts.loc[
                cluster_counts["medication_cluster"] == cluster,
                "N"
            ].sum()
        )

        if n_cluster == 0:
            continue

        row = {
            "medication_cluster": cluster,
            **mean_covars,
        }

        new_df = pd.DataFrame([row])
        new_df["medication_cluster"] = pd.Categorical(
            new_df["medication_cluster"],
            categories=CLUSTER_LEVELS,
            ordered=True
        )

        try:
            pred = fit.get_prediction(new_df).summary_frame(alpha=0.05)
            pred_rows.append({
                "organ_label": organ_label,
                "outcome": args.outcome,
                "medication_cluster": cluster,
                "N_cluster": n_cluster,
                "adjusted_mean": float(pred["mean"].iloc[0]),
                "adjusted_mean_se": float(pred["mean_se"].iloc[0]),
                "adjusted_ci_lo": float(pred["mean_ci_lower"].iloc[0]),
                "adjusted_ci_hi": float(pred["mean_ci_upper"].iloc[0]),
                "covariate_values": json.dumps(mean_covars),
            })
        except Exception as e:
            pred_rows.append({
                "organ_label": organ_label,
                "outcome": args.outcome,
                "medication_cluster": cluster,
                "N_cluster": n_cluster,
                "adjusted_mean": np.nan,
                "adjusted_mean_se": np.nan,
                "adjusted_ci_lo": np.nan,
                "adjusted_ci_hi": np.nan,
                "prediction_error": str(e),
                "covariate_values": json.dumps(mean_covars),
            })

    return pd.DataFrame(pred_rows)


def add_multiple_testing(effects):
    if effects.empty:
        return effects

    effects = effects.copy()
    effects["p_fdr_bh"] = np.nan
    effects["p_bonferroni"] = np.nan
    effects["fdr_significant_0.05"] = False
    effects["bonferroni_significant_0.05"] = False

    mask = effects["p"].notna()
    if mask.sum() > 0:
        pvals = effects.loc[mask, "p"].values
        _, p_fdr, _, _ = multipletests(pvals, alpha=0.05, method="fdr_bh")
        p_bonf = np.minimum(pvals * len(pvals), 1.0)

        effects.loc[mask, "p_fdr_bh"] = p_fdr
        effects.loc[mask, "p_bonferroni"] = p_bonf
        effects.loc[mask, "fdr_significant_0.05"] = p_fdr < 0.05
        effects.loc[mask, "bonferroni_significant_0.05"] = p_bonf < 0.05

    return effects


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    dat_raw = read_input(args.input_tsv)
    dat = clean_analysis_data(dat_raw, args.outcome)

    all_effects = []
    all_means = []
    all_counts = []
    all_summaries = []
    all_complete = []

    for organ_label in ORGAN_ORDER:
        print(f"\nRunning model for {organ_label}...")

        fit, counts, complete_df, summary = fit_one_organ(dat, organ_label, args)

        all_counts.append(counts)
        all_summaries.append(summary)

        if complete_df is not None and not complete_df.empty:
            tmp_complete = complete_df.copy()
            tmp_complete["organ_label_model"] = organ_label
            all_complete.append(tmp_complete)

        if fit is not None:
            effects = extract_cluster_effects(
                fit=fit,
                cluster_counts=counts,
                complete_df=complete_df,
                organ_label=organ_label,
                args=args
            )
            means = adjusted_means(
                fit=fit,
                cluster_counts=counts,
                complete_df=complete_df,
                organ_label=organ_label,
                args=args
            )
            all_effects.append(effects)
            all_means.append(means)

            print(f"  N complete cases: {complete_df.shape[0]}")
            print(f"  R2: {fit.rsquared:.4f}; adjusted R2: {fit.rsquared_adj:.4f}")
        else:
            print("  Model did not fit:", summary.get("error", ""))

    effects = pd.concat(all_effects, ignore_index=True) if all_effects else pd.DataFrame()
    effects = add_multiple_testing(effects)

    means = pd.concat(all_means, ignore_index=True) if all_means else pd.DataFrame()
    counts = pd.concat(all_counts, ignore_index=True) if all_counts else pd.DataFrame()
    summaries = pd.DataFrame(all_summaries)
    complete = pd.concat(all_complete, ignore_index=True) if all_complete else pd.DataFrame()

    prefix = os.path.join(args.out_dir, "medication_cluster_delta_clock_lm")

    effects_out = prefix + "_cluster_effects.tsv"
    means_out = prefix + "_adjusted_means.tsv"
    counts_out = prefix + "_cluster_counts.tsv"
    summaries_out = prefix + "_model_summaries.tsv"
    complete_out = prefix + "_complete_cases.tsv"
    metadata_out = prefix + "_metadata.json"

    effects.to_csv(effects_out, sep="\t", index=False)
    means.to_csv(means_out, sep="\t", index=False)
    counts.to_csv(counts_out, sep="\t", index=False)
    summaries.to_csv(summaries_out, sep="\t", index=False)

    if args.save_complete_cases:
        complete.to_csv(complete_out, sep="\t", index=False)

    metadata = {
        "input_tsv": args.input_tsv,
        "out_dir": args.out_dir,
        "outcome": args.outcome,
        "reference_cluster": args.reference_cluster,
        "cluster_levels": CLUSTER_LEVELS,
        "covariates_requested": DEFAULT_COVARIATES,
        "robust_cov": args.robust_cov,
        "min_n_cluster": args.min_n_cluster,
        "outputs": {
            "cluster_effects": effects_out,
            "adjusted_means": means_out,
            "cluster_counts": counts_out,
            "model_summaries": summaries_out,
            "complete_cases": complete_out if args.save_complete_cases else None,
        },
        "n_input_rows": int(dat_raw.shape[0]),
        "n_rows_after_requested_cluster_filter": int(dat.shape[0]),
    }

    with open(metadata_out, "w") as f:
        json.dump(metadata, f, indent=2)

    print("\n============================================================")
    print("Finished medication-cluster delta-clock linear models.")
    print("Outputs:")
    print("  Cluster effects:", effects_out)
    print("  Adjusted means:", means_out)
    print("  Cluster counts:", counts_out)
    print("  Model summaries:", summaries_out)
    if args.save_complete_cases:
        print("  Complete cases:", complete_out)
    print("  Metadata:", metadata_out)
    print("============================================================")


if __name__ == "__main__":
    main()