#!/usr/bin/env python3

import argparse
import math
import re
import warnings
from pathlib import Path
import numpy as np
import pandas as pd

try:
    from scipy import stats
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False


# ============================================================
# Helper functions
# ============================================================

def log(msg):
    print(msg, flush=True)


def safe_name(x):
    x = re.sub(r"[^A-Za-z0-9_.-]+", "_", x)
    x = x.replace("__", "_")
    return x.strip("_")


def pvalue_from_t(tval, df):
    if not np.isfinite(tval) or not np.isfinite(df) or df <= 0:
        return np.nan

    if HAVE_SCIPY:
        return 2.0 * stats.t.sf(abs(tval), df)

    # Fallback: normal approximation
    return math.erfc(abs(tval) / math.sqrt(2.0))


def zscore(x):
    x = pd.to_numeric(x, errors="coerce")
    mu = x.mean(skipna=True)
    sd = x.std(skipna=True, ddof=1)

    if not np.isfinite(sd) or sd == 0:
        return x * np.nan

    return (x - mu) / sd


def fit_ols(df, outcome, exposure, covariates, standardize_outcome=False):
    """
    Complete-case OLS using numpy.

    Model:
        outcome ~ exposure + covariates

    exposure:
        APOE_e4e4_vs_e2e2, where e2/e2 = 0 and e4/e4 = 1

    Returns beta/se/t/p for the exposure.
    """

    cols = [outcome, exposure] + covariates
    d = df[cols].copy()

    for c in cols:
        d[c] = pd.to_numeric(d[c], errors="coerce")

    d = d.replace([np.inf, -np.inf], np.nan).dropna()

    if d.shape[0] < len(covariates) + 5:
        return {
            "outcome": outcome,
            "n": d.shape[0],
            "n_e2e2": np.nan,
            "n_e4e4": np.nan,
            "mean_e2e2": np.nan,
            "mean_e4e4": np.nan,
            "mean_diff_e4e4_minus_e2e2": np.nan,
            "beta": np.nan,
            "se": np.nan,
            "t": np.nan,
            "p": np.nan,
            "r2": np.nan,
            "df_resid": np.nan,
        }

    y_raw = d[outcome].astype(float)

    n_e2e2 = int((d[exposure] == 0).sum())
    n_e4e4 = int((d[exposure] == 1).sum())

    mean_e2e2 = y_raw[d[exposure] == 0].mean()
    mean_e4e4 = y_raw[d[exposure] == 1].mean()
    mean_diff = mean_e4e4 - mean_e2e2

    y = y_raw.copy()

    if standardize_outcome:
        sd = y.std(ddof=1)
        if not np.isfinite(sd) or sd == 0:
            return {
                "outcome": outcome,
                "n": d.shape[0],
                "n_e2e2": n_e2e2,
                "n_e4e4": n_e4e4,
                "mean_e2e2": mean_e2e2,
                "mean_e4e4": mean_e4e4,
                "mean_diff_e4e4_minus_e2e2": mean_diff,
                "beta": np.nan,
                "se": np.nan,
                "t": np.nan,
                "p": np.nan,
                "r2": np.nan,
                "df_resid": np.nan,
            }

        y = (y - y.mean()) / sd

    X = d[[exposure] + covariates].astype(float).to_numpy()
    X = np.column_stack([np.ones(X.shape[0]), X])
    y_np = y.to_numpy(dtype=float)

    n, p = X.shape
    df_resid = n - p

    try:
        beta_hat = np.linalg.lstsq(X, y_np, rcond=None)[0]
        resid = y_np - X @ beta_hat

        rss = float(resid.T @ resid)
        tss = float(((y_np - y_np.mean()) ** 2).sum())

        sigma2 = rss / df_resid
        xtx_inv = np.linalg.pinv(X.T @ X)
        cov_beta = sigma2 * xtx_inv
        se = np.sqrt(np.diag(cov_beta))

        # exposure is column 1 after intercept
        beta_exposure = beta_hat[1]
        se_exposure = se[1]
        tval = beta_exposure / se_exposure if se_exposure > 0 else np.nan
        pval = pvalue_from_t(tval, df_resid)
        r2 = 1.0 - rss / tss if tss > 0 else np.nan

    except Exception as e:
        warnings.warn(f"OLS failed for {outcome}: {e}")
        beta_exposure = np.nan
        se_exposure = np.nan
        tval = np.nan
        pval = np.nan
        r2 = np.nan

    return {
        "outcome": outcome,
        "n": n,
        "n_e2e2": n_e2e2,
        "n_e4e4": n_e4e4,
        "mean_e2e2": mean_e2e2,
        "mean_e4e4": mean_e4e4,
        "mean_diff_e4e4_minus_e2e2": mean_diff,
        "beta": beta_exposure,
        "se": se_exposure,
        "t": tval,
        "p": pval,
        "r2": r2,
        "df_resid": df_resid,
    }


def rename_result(prefix, res):
    return {
        f"{prefix}_{k}": v
        for k, v in res.items()
    }


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="APOE e4/e4 vs e2/e2 association for one dementia L'EPOCH clock and matched BAG."
    )

    parser.add_argument("--apoe-file", required=True)
    parser.add_argument("--lepoch-file", required=True)
    parser.add_argument("--bag-file", required=True)
    parser.add_argument("--cov-file", required=True)
    parser.add_argument("--outdir", required=True)

    parser.add_argument("--lepoch-col", required=True)
    parser.add_argument("--bag-col", required=True)

    parser.add_argument("--apoe-id-col", default="FID")
    parser.add_argument("--cov-id-col", default="eid")
    parser.add_argument("--lepoch-id-col", default="participant_id")
    parser.add_argument("--bag-id-col", default="participant_id")

    args = parser.parse_args()

    apoe_file = Path(args.apoe_file)
    lepoch_file = Path(args.lepoch_file)
    bag_file = Path(args.bag_file)
    cov_file = Path(args.cov_file)
    outdir = Path(args.outdir)

    outdir.mkdir(parents=True, exist_ok=True)

    lepoch_col = args.lepoch_col
    bag_col = args.bag_col

    log("============================================================")
    log("APOE e4/e4 vs e2/e2 association")
    log(f"L'EPOCH clock: {lepoch_col}")
    log(f"Matched BAG:   {bag_col}")
    log("============================================================")

    # ------------------------------------------------------------
    # Covariates
    # ------------------------------------------------------------

    age_col = "age_at_recruitment_f21022_0_0"
    sex_col = "sex_f31_0_0"
    bmi_col = "body_mass_index_bmi_f23104_0_0"

    dbp_cols = [
        "diastolic_blood_pressure_automated_reading_f4079_0_0",
        "diastolic_blood_pressure_automated_reading_f4079_0_1",
    ]

    sbp_cols = [
        "systolic_blood_pressure_automated_reading_f4080_0_0",
        "systolic_blood_pressure_automated_reading_f4080_0_1",
    ]

    pc_cols = [
        f"genetic_principal_components_f22009_0_{i}"
        for i in range(1, 6)
    ]

    cov_usecols = [
        args.cov_id_col,
        age_col,
        sex_col,
        bmi_col,
    ] + dbp_cols + sbp_cols + pc_cols

    # ------------------------------------------------------------
    # Load APOE
    # ------------------------------------------------------------

    log("Loading APOE file...")

    apoe_usecols = [
        args.apoe_id_col,
        "APOE_genotype",
        "APOE_group",
        "APOE_e2_count",
        "APOE_e4_count",
    ]

    apoe = pd.read_csv(
        apoe_file,
        sep="\t",
        usecols=lambda c: c in apoe_usecols,
        low_memory=False,
    )

    if args.apoe_id_col not in apoe.columns:
        raise ValueError(f"APOE ID column not found: {args.apoe_id_col}")

    if "APOE_genotype" not in apoe.columns:
        raise ValueError("APOE file must contain APOE_genotype.")

    apoe = apoe.rename(columns={args.apoe_id_col: "participant_id"})
    apoe["participant_id"] = pd.to_numeric(apoe["participant_id"], errors="coerce").astype("Int64")

    # Keep only homozygous APOE2 and APOE4 carriers
    apoe = apoe.loc[apoe["APOE_genotype"].isin(["e2/e2", "e4/e4"])].copy()

    apoe["APOE_e4e4_vs_e2e2"] = np.where(
        apoe["APOE_genotype"] == "e4/e4",
        1,
        0,
    )

    apoe_counts = (
        apoe["APOE_genotype"]
        .value_counts(dropna=False)
        .rename_axis("APOE_genotype")
        .reset_index(name="N_APOE_homozygous_before_merge")
    )

    log("APOE homozygous counts before merge:")
    log(apoe_counts.to_string(index=False))

    # Helpful warning because negative FID often means sample IDs may not be real UKB eids.
    n_nonpositive = int((apoe["participant_id"] <= 0).sum())
    if n_nonpositive > 0:
        warnings.warn(
            f"{n_nonpositive} APOE participant_id values are <= 0 after using "
            f"{args.apoe_id_col} as participant_id. If merge size is zero or tiny, "
            "the APOE sample IDs may not be UKB eids. Check whether IID or another "
            "sample column should be used instead."
        )

    # ------------------------------------------------------------
    # Load L'EPOCH clock
    # ------------------------------------------------------------

    log("Loading L'EPOCH file...")

    lepoch = pd.read_csv(
        lepoch_file,
        sep="\t",
        usecols=[args.lepoch_id_col, lepoch_col],
        low_memory=False,
    )

    lepoch = lepoch.rename(columns={args.lepoch_id_col: "participant_id"})
    lepoch["participant_id"] = pd.to_numeric(lepoch["participant_id"], errors="coerce").astype("Int64")
    lepoch[lepoch_col] = pd.to_numeric(lepoch[lepoch_col], errors="coerce")

    # ------------------------------------------------------------
    # Load BAG file
    # ------------------------------------------------------------

    log("Loading BAG file...")

    bag = pd.read_csv(
        bag_file,
        sep="\t",
        usecols=[args.bag_id_col, bag_col],
        low_memory=False,
    )

    bag = bag.rename(columns={args.bag_id_col: "participant_id"})
    bag["participant_id"] = pd.to_numeric(bag["participant_id"], errors="coerce").astype("Int64")
    bag[bag_col] = pd.to_numeric(bag[bag_col], errors="coerce")

    # ------------------------------------------------------------
    # Load covariates
    # ------------------------------------------------------------

    log("Loading covariate file...")

    cov = pd.read_csv(
        cov_file,
        sep=",",
        usecols=cov_usecols,
        low_memory=False,
    )

    cov = cov.rename(columns={args.cov_id_col: "participant_id"})
    cov["participant_id"] = pd.to_numeric(cov["participant_id"], errors="coerce").astype("Int64")

    cov2 = pd.DataFrame()
    cov2["participant_id"] = cov["participant_id"]

    cov2["age"] = pd.to_numeric(cov[age_col], errors="coerce")
    cov2["sex"] = pd.to_numeric(cov[sex_col], errors="coerce")
    cov2["bmi"] = pd.to_numeric(cov[bmi_col], errors="coerce")

    dbp_numeric = cov[dbp_cols].apply(pd.to_numeric, errors="coerce")
    sbp_numeric = cov[sbp_cols].apply(pd.to_numeric, errors="coerce")

    cov2["diastolic_bp"] = dbp_numeric.mean(axis=1, skipna=True)
    cov2["systolic_bp"] = sbp_numeric.mean(axis=1, skipna=True)

    for i, pc in enumerate(pc_cols, start=1):
        cov2[f"PC{i}"] = pd.to_numeric(cov[pc], errors="coerce")

    covariates = [
        "age",
        "sex",
        "bmi",
        "systolic_bp",
        "diastolic_bp",
        "PC1",
        "PC2",
        "PC3",
        "PC4",
        "PC5",
    ]

    # ------------------------------------------------------------
    # Merge
    # ------------------------------------------------------------

    log("Merging files...")

    df = (
        apoe[["participant_id", "APOE_genotype", "APOE_e4e4_vs_e2e2"]]
        .merge(cov2, on="participant_id", how="inner")
        .merge(lepoch[["participant_id", lepoch_col]], on="participant_id", how="inner")
        .merge(bag[["participant_id", bag_col]], on="participant_id", how="inner")
    )

    log(f"Merged N before complete-case filtering: {df.shape[0]}")

    if df.shape[0] == 0:
        raise ValueError(
            "No merged rows. Check APOE participant IDs. "
            "You specified FID as participant_id; however, if FID values look like -1, -2, etc., "
            "they may not match UKB eid. In that case, regenerate APOE status with real UKB sample IDs."
        )

    merge_counts = (
        df["APOE_genotype"]
        .value_counts(dropna=False)
        .rename_axis("APOE_genotype")
        .reset_index(name="N_after_merge_before_complete_case")
    )

    log("APOE counts after merge:")
    log(merge_counts.to_string(index=False))

    # ------------------------------------------------------------
    # Association model for L'EPOCH
    # ------------------------------------------------------------

    exposure = "APOE_e4e4_vs_e2e2"

    log("Running L'EPOCH association model...")

    lepoch_raw = fit_ols(
        df=df,
        outcome=lepoch_col,
        exposure=exposure,
        covariates=covariates,
        standardize_outcome=False,
    )

    lepoch_std = fit_ols(
        df=df,
        outcome=lepoch_col,
        exposure=exposure,
        covariates=covariates,
        standardize_outcome=True,
    )

    # ------------------------------------------------------------
    # Association model for matched BAG
    # ------------------------------------------------------------

    log("Running matched BAG association model...")

    bag_raw = fit_ols(
        df=df,
        outcome=bag_col,
        exposure=exposure,
        covariates=covariates,
        standardize_outcome=False,
    )

    bag_std = fit_ols(
        df=df,
        outcome=bag_col,
        exposure=exposure,
        covariates=covariates,
        standardize_outcome=True,
    )

    # ------------------------------------------------------------
    # Fold comparison
    # ------------------------------------------------------------

    lepoch_beta_std = lepoch_std["beta"]
    bag_beta_std = bag_std["beta"]

    abs_lepoch_beta_std = abs(lepoch_beta_std) if np.isfinite(lepoch_beta_std) else np.nan
    abs_bag_beta_std = abs(bag_beta_std) if np.isfinite(bag_beta_std) else np.nan

    if np.isfinite(abs_bag_beta_std) and abs_bag_beta_std != 0:
        lepoch_vs_bag_abs_effect_fold = abs_lepoch_beta_std / abs_bag_beta_std
    else:
        lepoch_vs_bag_abs_effect_fold = np.nan

    if np.isfinite(bag_beta_std) and bag_beta_std != 0:
        lepoch_vs_bag_signed_effect_ratio = lepoch_beta_std / bag_beta_std
    else:
        lepoch_vs_bag_signed_effect_ratio = np.nan

    if np.isfinite(lepoch_beta_std) and np.isfinite(bag_beta_std):
        same_direction = np.sign(lepoch_beta_std) == np.sign(bag_beta_std)
    else:
        same_direction = np.nan

    if np.isfinite(abs_lepoch_beta_std) and np.isfinite(abs_bag_beta_std):
        stronger_marker = "L_EPOCH" if abs_lepoch_beta_std > abs_bag_beta_std else "BAG"
    else:
        stronger_marker = "NA"

    # ------------------------------------------------------------
    # Output one-row result
    # ------------------------------------------------------------

    result = {
        "lepoch_clock": lepoch_col,
        "matched_bag": bag_col,
        "exposure": "APOE_e4e4_vs_e2e2",
        "reference_group": "APOE_e2e2",
        "comparison_group": "APOE_e4e4",
        "model_covariates": "age, sex, bmi, systolic_bp, diastolic_bp, PC1, PC2, PC3, PC4, PC5",
    }

    result.update(rename_result("lepoch_raw", lepoch_raw))
    result.update(rename_result("lepoch_std", lepoch_std))
    result.update(rename_result("bag_raw", bag_raw))
    result.update(rename_result("bag_std", bag_std))

    result.update({
        "abs_lepoch_beta_std": abs_lepoch_beta_std,
        "abs_bag_beta_std": abs_bag_beta_std,
        "lepoch_vs_bag_abs_effect_fold": lepoch_vs_bag_abs_effect_fold,
        "lepoch_vs_bag_signed_effect_ratio": lepoch_vs_bag_signed_effect_ratio,
        "same_direction": same_direction,
        "stronger_marker": stronger_marker,
    })

    out_prefix = safe_name(lepoch_col)
    out_file = outdir / f"{out_prefix}_association_and_fold.tsv"
    qc_file = outdir / f"{out_prefix}_qc_counts.tsv"

    pd.DataFrame([result]).to_csv(out_file, sep="\t", index=False)

    qc_rows = []
    qc_rows.append({
        "stage": "APOE_homozygous_before_merge",
        "N_total": apoe.shape[0],
        "N_e2e2": int((apoe["APOE_genotype"] == "e2/e2").sum()),
        "N_e4e4": int((apoe["APOE_genotype"] == "e4/e4").sum()),
    })
    qc_rows.append({
        "stage": "after_merge_before_complete_case",
        "N_total": df.shape[0],
        "N_e2e2": int((df["APOE_genotype"] == "e2/e2").sum()),
        "N_e4e4": int((df["APOE_genotype"] == "e4/e4").sum()),
    })
    qc_rows.append({
        "stage": "lepoch_complete_case",
        "N_total": lepoch_raw["n"],
        "N_e2e2": lepoch_raw["n_e2e2"],
        "N_e4e4": lepoch_raw["n_e4e4"],
    })
    qc_rows.append({
        "stage": "bag_complete_case",
        "N_total": bag_raw["n"],
        "N_e2e2": bag_raw["n_e2e2"],
        "N_e4e4": bag_raw["n_e4e4"],
    })

    pd.DataFrame(qc_rows).to_csv(qc_file, sep="\t", index=False)

    log("============================================================")
    log(f"Wrote result: {out_file}")
    log(f"Wrote QC:     {qc_file}")
    log("Key interpretation:")
    log("  lepoch_std_beta = adjusted standardized difference in L'EPOCH, e4/e4 vs e2/e2")
    log("  bag_std_beta    = adjusted standardized difference in matched BAG, e4/e4 vs e2/e2")
    log("  lepoch_vs_bag_abs_effect_fold = abs(lepoch_std_beta) / abs(bag_std_beta)")
    log("============================================================")


if __name__ == "__main__":
    main()