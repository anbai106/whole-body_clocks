import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
import warnings
warnings.filterwarnings("ignore")
import argparse

parser = argparse.ArgumentParser(description="Cox survival analysis with BAGs")

# Mandatory-ish arguments
parser.add_argument(
    "--icd_tsv",
    type=str,
    help="Input TSV with participant_id, BAGs, case, and event date",
)
parser.add_argument(
    "--output_tsv",
    type=str,
    help="Output TSV to store hazard ratio results",
)


def construct_survival_data(icd_tsv: str):
    """
    Construct survival data for BAG–disease Cox models.

    - BAGs ending with _ProtBAG or _MetBAG (and Female_BAG) are measured at baseline.
    - BAGs ending with _MRIBAG are measured at the second (imaging) visit.
    - We compute:
        - time_baseline: follow-up time from baseline assessment
        - time_imaging : follow-up time from imaging assessment

    Returns
    -------
    data : pd.DataFrame
        Merged dataset with BAGs, covariates, event indicator, dates, and time variables.
    baseline_bags : list of str
        BAG variables measured at baseline (use time_baseline).
    mri_bags : list of str
        BAG variables measured at imaging (use time_imaging, and filter time_imaging >= 0).
    n_case : int
        Number of cases (after removing negative baseline times).
    n_noncase : int
        Number of non-cases (after removing negative baseline times).
    """

    # ---------- 1) Read disease + BAG data ----------
    df_disease = pd.read_csv(icd_tsv, sep="\t")[
        [
            "participant_id",
            "Brain_ProtBAG", "Brain_MRIBAG", "Eye_ProtBAG",
            "Heart_ProtBAG", "Heart_MRIBAG", "Pulmonary_ProtBAG",
            "Hepatic_ProtBAG", "Liver_MRIBAG",
            "Digestive_MetBAG", "Hepatic_MetBAG", "Immune_MetBAG",
            "Pancreas_MRIBAG", "Spleen_MRIBAG", "Metabolic_MetBAG",
            "Renal_ProtBAG",
            "Endocrine_ProtBAG", "Immune_ProtBAG", "Adipose_MRIBAG", "Skin_ProtBAG",
            "Reproductive_female_ProtBAG",
            "case",
            "date",  # event or censoring date (ICD endpoint)
        ]
    ]

    # Full BAG list
    bag_var = [
        "Brain_ProtBAG", "Brain_MRIBAG", "Eye_ProtBAG",
        "Heart_ProtBAG", "Heart_MRIBAG", "Pulmonary_ProtBAG",
        "Hepatic_ProtBAG", "Liver_MRIBAG",
        "Digestive_MetBAG", "Hepatic_MetBAG", "Immune_MetBAG",
        "Pancreas_MRIBAG", "Spleen_MRIBAG", "Metabolic_MetBAG",
        "Renal_ProtBAG",
        "Endocrine_ProtBAG", "Immune_ProtBAG", "Adipose_MRIBAG", "Skin_ProtBAG",
        "Reproductive_female_ProtBAG"
    ]

    # Split into baseline vs MRI BAGs
    mri_bags = [b for b in bag_var if b.endswith("MRIBAG")]
    baseline_bags = [b for b in bag_var if b not in mri_bags]

    # ---------- 2) Covariates (age at recruitment) ----------
    cov = pd.read_csv(
        "/cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv"
    )[["eid", "age_at_recruitment_f21022_0_0", "smoking_status_f20116_0_0",
                         'body_mass_index_bmi_f23104_0_0']]
    cov = cov.rename(
        columns={
            "eid": "participant_id",
            "age_at_recruitment_f21022_0_0": "Age",
            "smoking_status_f20116_0_0": "Smoking",
        }
    )

    # ---------- 3) Assessment centre dates (baseline + imaging) ----------
    df_date_center = pd.read_csv(
        "/cbica/home/wenju/Reproducibile_paper/Multiorgan_Subtype/data/PWAS/UKBB_fullsample_death_variables.csv"
    )
    df_date_center = df_date_center.rename(columns={"eid": "participant_id"})
    df_date_center = df_date_center[
        [
            "participant_id",
            "date_of_attending_assessment_centre_f53_0_0",  # baseline
            "date_of_attending_assessment_centre_f53_2_0",  # imaging
        ]
    ]

    # ---------- 4) Merge all ----------
    data = cov.merge(df_disease, on="participant_id")
    data = data.merge(df_date_center, on="participant_id", how="left")

    # Convert dates
    data["date"] = pd.to_datetime(data["date"])
    data["date_of_attending_assessment_centre_f53_0_0"] = pd.to_datetime(
        data["date_of_attending_assessment_centre_f53_0_0"]
    )
    data["date_of_attending_assessment_centre_f53_2_0"] = pd.to_datetime(
        data["date_of_attending_assessment_centre_f53_2_0"]
    )

    # Global censor date (same logic as original: max event date + 2 days)
    global_end_date = data["date"].max() + pd.Timedelta(days=2)

    # ---------- 5) time_baseline: from baseline assessment ----------
    # Initialize with censoring time
    data["time_baseline"] = (
        global_end_date - data["date_of_attending_assessment_centre_f53_0_0"]
    ).dt.days

    # For those with an event, use event date instead of censoring
    mask_event = data["date"].notna()
    data.loc[mask_event, "time_baseline"] = (
        data.loc[mask_event, "date"]
        - data.loc[mask_event, "date_of_attending_assessment_centre_f53_0_0"]
    ).dt.days

    # ---------- 6) time_imaging: from imaging assessment ----------
    # If imaging date is missing, time_imaging will be NaN
    data["time_imaging"] = np.nan
    has_imaging = data["date_of_attending_assessment_centre_f53_2_0"].notna()

    # Start with censoring
    data.loc[has_imaging, "time_imaging"] = (
        global_end_date - data.loc[has_imaging, "date_of_attending_assessment_centre_f53_2_0"]
    ).dt.days

    # For those with an event, use event date
    mask_event_imaging = has_imaging & data["date"].notna()
    data.loc[mask_event_imaging, "time_imaging"] = (
        data.loc[mask_event_imaging, "date"]
        - data.loc[mask_event_imaging, "date_of_attending_assessment_centre_f53_2_0"]
    ).dt.days

    # ---------- 7) Remove negative baseline times ----------
    # (Events before baseline assessment are not valid for incident analysis)
    data = data.loc[data["time_baseline"] >= 0].copy()

    # Recompute case / noncase counts after filtering
    n_case = int((data["case"] == 1).sum())
    n_noncase = int((data["case"] == 0).sum())

    if n_case == 0 or n_noncase == 0:
        raise Exception("There is no cases or non-cases after filtering by baseline time >= 0.")

    return data, baseline_bags, mri_bags, n_case, n_noncase


def sa_hazard_ratio(icd_tsv: str, output_tsv: str):
    """
    Run Cox models for each BAG.

    - For baseline BAGs (ProtBAG, MetBAG, Female_BAG): use time_baseline.
    - For MRIBAGs: use time_imaging, restricted to rows with time_imaging >= 0.
    """
    data, baseline_bags, mri_bags, _, _ = construct_survival_data(icd_tsv)
    all_var = baseline_bags + mri_bags

    results = []
    for v in all_var:
        # Choose appropriate time column
        if v in baseline_bags:
            time_col = "time_baseline"
        else:
            time_col = "time_imaging"

        cols = [v, "Age", "Smoking",
                         'body_mass_index_bmi_f23104_0_0', time_col, "case"]

        if not set(cols).issubset(data.columns):
            results.append([np.nan, np.nan, np.nan, np.nan, v, np.nan, np.nan])
            continue

        df = data[cols].dropna().copy()

        # For MRIBAGs, exclude negative imaging times (events before imaging)
        if v in mri_bags:
            df = df[df[time_col] >= 0]

        # per-variable event counts (after dropna / filtering)
        n_case = int((df["case"] == 1).sum())
        n_noncase = int((df["case"] == 0).sum())

        # require >10 cases and at least 1 noncase
        if n_case <= 10 or n_noncase == 0:
            results.append([np.nan, np.nan, np.nan, np.nan, v, n_case, n_noncase])
            continue

        # standardize predictor if variance > 0
        sd = float(df[v].std())
        if not np.isfinite(sd) or sd == 0.0:
            results.append([np.nan, np.nan, np.nan, np.nan, v, n_case, n_noncase])
            continue

        df[v] = (df[v] - df[v].mean()) / sd

        # fit Cox model
        try:
            cph = CoxPHFitter()
            cph.fit(df, duration_col=time_col, event_col="case")
            hr = float(cph.hazard_ratios_.loc[v])
            ci_lo = float(np.exp(cph.confidence_intervals_.loc[v, "95% lower-bound"]))
            ci_hi = float(np.exp(cph.confidence_intervals_.loc[v, "95% upper-bound"]))
            pval = float(cph.summary.loc[v, "p"])
            results.append([hr, ci_lo, ci_hi, pval, v, n_case, n_noncase])
        except Exception:
            # on any fitting error, record NA
            results.append([np.nan, np.nan, np.nan, np.nan, v, n_case, n_noncase])

    test_result = pd.DataFrame(
        results,
        columns=[
            "hazard_ratio",
            "CI_lower_bound",
            "CI_upper_bound",
            "p_value",
            "var",
            "N_case",
            "N_noncase",
        ],
    )
    test_result.to_csv(output_tsv, index=False, sep="\t", encoding="utf-8")
    return test_result


def main(options):
    sa_hazard_ratio(options.icd_tsv, options.output_tsv)


if __name__ == "__main__":
    commandline = parser.parse_known_args()
    options = commandline[0]
    if commandline[1]:
        raise Exception("unknown arguments: %s" % parser.parse_known_args()[1])
    main(options)