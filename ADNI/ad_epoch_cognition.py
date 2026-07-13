#!/usr/bin/env python3

import argparse
import glob
import os
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats


# ============================================================
# Exact ADNI iSTAGING column names
# ============================================================

ID_COL = "PTID"
DATE_COL = "Date"
ID_CANDIDATES = ["PTID", "participant_id", "RID"]

COGNITIVE_COLS = [
    "Animal_Fluency",
    "RAVLT",
    "TMT_A",
    "TMT_B",
]

COMPARATOR_BIOMARKER_COLS = [
    "SPARE_BA",
    "SPARE_AD",
    "Abeta_CSF",
    "Tau_CSF",
    "PTau_CSF",
]

COVARIATE_COLS = [
    "Age",
    "Sex",
    "Education_Years",
    "APOE4_Alleles",
    "DLICV",
    "SITE",
]

DEFAULT_EPOCH_COL_CANDIDATES = [
    "adni_brain_mri_ad_lepoch_acceleration_z",
    "adni_brain_mri_ad_lepoch_clock_acceleration_z",
    "adni_brain_mri_ad_epoch_acceleration_z",
    "adni_brain_mri_ad_epoch_clock_acceleration_z",
    "adni_brain_mri_ad_lepoch_acceleration_years",
    "adni_brain_mri_ad_lepoch_clock_acceleration_years",
    "adni_brain_mri_ad_epoch_acceleration_years",
    "adni_brain_mri_ad_epoch_clock_acceleration_years",
    "adni_brain_mri_ad_lepoch_risk_score",
    "adni_brain_mri_ad_epoch_risk_score",
]


# ============================================================
# Utility functions
# ============================================================

def read_table(path: str) -> pd.DataFrame:
    if path.endswith(".tsv") or path.endswith(".txt"):
        return pd.read_csv(path, sep="\t", low_memory=False)
    if path.endswith(".csv"):
        return pd.read_csv(path, low_memory=False)
    if path.endswith(".xlsx"):
        return pd.read_excel(path)
    return pd.read_csv(path, sep="\t", low_memory=False)


def normalize_id_series(s: pd.Series) -> pd.Series:
    out = s.astype(str).str.strip()
    out = out.str.replace(r"\.0$", "", regex=True)
    return out


def require_columns(df: pd.DataFrame, cols: List[str], label: str):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in {label}: {missing}")


def find_id_col(df: pd.DataFrame, label: str) -> str:
    for c in ID_CANDIDATES:
        if c in df.columns:
            return c
    raise ValueError(
        f"Could not find an ADNI ID column in {label}. "
        f"Tried: {ID_CANDIDATES}. Available columns include: {list(df.columns)[:40]}"
    )


def parse_dates(s: pd.Series) -> pd.Series:
    x = s.copy()
    x = x.replace([0, 0.0, "0", "0.0", "", "NA", "NaN", "nan", "None", "-1", -1], np.nan)
    parsed = pd.to_datetime(x, errors="coerce")

    numeric = pd.to_numeric(x, errors="coerce")
    excel_mask = numeric.between(20000, 60000)
    if excel_mask.any():
        excel_dates = pd.to_datetime(numeric, unit="D", origin="1899-12-30", errors="coerce")
        parsed = parsed.where(~excel_mask, excel_dates)

    return parsed


def numeric_clean(s: pd.Series) -> pd.Series:
    if pd.api.types.is_numeric_dtype(s):
        return pd.to_numeric(s, errors="coerce")

    x = s.astype(str).str.strip()
    x = x.str.replace(",", "", regex=False)
    x = x.str.replace("<", "", regex=False)
    x = x.str.replace(">", "", regex=False)
    x = x.replace(["", "NA", "NaN", "nan", "None", "-1"], np.nan)
    return pd.to_numeric(x, errors="coerce")


def zscore_array(x: np.ndarray) -> Optional[np.ndarray]:
    x = np.asarray(x, dtype=float)
    sd = np.nanstd(x, ddof=1)
    if not np.isfinite(sd) or sd <= 0:
        return None
    return (x - np.nanmean(x)) / sd


def resolve_prediction_file(path_or_dir: str, prefer_longitudinal: bool = False) -> str:
    if os.path.isfile(path_or_dir):
        return path_or_dir

    if not os.path.isdir(path_or_dir):
        raise FileNotFoundError(f"Cannot find AD EPOCH path: {path_or_dir}")

    patterns = [
        "*predictions*.tsv",
        "*prediction*.tsv",
        "*longitudinal*.tsv",
        "*.tsv",
    ]

    files = []
    for pat in patterns:
        files.extend(glob.glob(os.path.join(path_or_dir, pat)))

    files = sorted(set(files))

    if len(files) == 0:
        raise FileNotFoundError(f"No TSV prediction files found under: {path_or_dir}")

    def score_file(p: str) -> int:
        b = os.path.basename(p).lower()
        score = 0
        if "predictions" in b:
            score += 20
        if "prediction" in b:
            score += 10
        if prefer_longitudinal and "longitudinal" in b:
            score += 20
        if "test" in b:
            score -= 20
        if "summary" in b or "manifest" in b or "coefficient" in b or "performance" in b:
            score -= 50
        return score

    files = sorted(files, key=score_file, reverse=True)
    return files[0]


def detect_epoch_col(df: pd.DataFrame, requested_col: str = "auto") -> str:
    if requested_col and requested_col != "auto":
        if requested_col not in df.columns:
            raise ValueError(f"Requested AD EPOCH column not found: {requested_col}")
        return requested_col

    for c in DEFAULT_EPOCH_COL_CANDIDATES:
        if c in df.columns:
            return c

    lower_cols = [(c, c.lower()) for c in df.columns]

    for key in [
        "acceleration_z",
        "clock_acceleration_z",
        "acceleration_year",
        "clock_acceleration_year",
        "risk_score",
    ]:
        hits = [
            c for c, cl in lower_cols
            if key in cl and ("ad" in cl or "epoch" in cl or "lepoch" in cl)
        ]
        if len(hits) > 0:
            return hits[0]

    hits = [c for c, cl in lower_cols if "risk_score" in cl]
    if len(hits) > 0:
        return hits[0]

    raise ValueError(
        "Could not detect AD EPOCH score column. "
        "Please pass --baseline_epoch_col or --longitudinal_epoch_col explicitly."
    )


def detect_date_col(df: pd.DataFrame) -> Optional[str]:
    for c in [
        "Date",
        "EXAMDATE",
        "ExamDate",
        "scan_date",
        "Scan_Date",
        "MRI_Date",
        "MRI_date",
        "epoch_date",
        "baseline_date",
        "selected_baseline_date",
    ]:
        if c in df.columns:
            return c
    return None


# ============================================================
# Load ADNI iSTAGING and AD EPOCH scores
# ============================================================

def load_adni_istaging(path: str) -> pd.DataFrame:
    df = read_table(path)

    exact_needed = [ID_COL, DATE_COL] + COGNITIVE_COLS + COMPARATOR_BIOMARKER_COLS + COVARIATE_COLS
    require_columns(df, exact_needed, "ADNI iSTAGING table")

    df = df.copy()
    df[ID_COL] = normalize_id_series(df[ID_COL])
    df[DATE_COL] = parse_dates(df[DATE_COL])
    df = df[df[DATE_COL].notna()].copy()

    numeric_cols = (
        COGNITIVE_COLS
        + COMPARATOR_BIOMARKER_COLS
        + ["Age", "Education_Years", "APOE4_Alleles", "DLICV"]
    )

    for c in numeric_cols:
        df[c] = numeric_clean(df[c])

    print("[INFO] ADNI iSTAGING loaded")
    print(f"       rows with valid Date: {len(df):,}")
    print(f"       unique PTID: {df[ID_COL].nunique():,}")

    return df


def load_epoch_scores(
    path_or_dir: str,
    value_name: str,
    requested_col: str = "auto",
    prefer_longitudinal: bool = False
) -> pd.DataFrame:
    path = resolve_prediction_file(path_or_dir, prefer_longitudinal=prefer_longitudinal)
    df = read_table(path)

    id_col = find_id_col(df, f"AD EPOCH prediction file: {path}")
    epoch_col = detect_epoch_col(df, requested_col=requested_col)
    date_col = detect_date_col(df)

    keep = [id_col, epoch_col]
    if date_col is not None and date_col not in keep:
        keep.append(date_col)

    out = df[keep].copy()
    out = out.rename(columns={id_col: ID_COL})
    out[ID_COL] = normalize_id_series(out[ID_COL])

    out[value_name] = numeric_clean(out[epoch_col])

    if epoch_col != value_name:
        out = out.drop(columns=[epoch_col])

    if date_col is not None:
        out = out.rename(columns={date_col: "epoch_date"})
        out["epoch_date"] = parse_dates(out["epoch_date"])
    else:
        out["epoch_date"] = pd.NaT

    out = out.dropna(subset=[ID_COL, value_name]).copy()

    print(f"[INFO] Loaded {value_name}")
    print(f"       file: {path}")
    print(f"       ID column: {id_col} -> {ID_COL}")
    print(f"       score column: {epoch_col}")
    print(f"       date column: {date_col if date_col is not None else 'none'}")
    print(f"       rows: {len(out):,}")
    print(f"       unique PTID: {out[ID_COL].nunique():,}")

    return out


def add_anchor_date(epoch_df: pd.DataFrame, adni: pd.DataFrame) -> pd.DataFrame:
    first_dates = (
        adni[[ID_COL, DATE_COL]]
        .dropna()
        .sort_values(DATE_COL)
        .groupby(ID_COL, as_index=False)
        .first()
        .rename(columns={DATE_COL: "first_adni_date"})
    )

    out = epoch_df.merge(first_dates, on=ID_COL, how="left")
    out["anchor_date"] = out["epoch_date"].where(out["epoch_date"].notna(), out["first_adni_date"])
    out = out.drop(columns=["first_adni_date"])

    return out


def compute_slope_per_year(
    long_epoch: pd.DataFrame,
    value_col: str,
    out_col: str,
    min_scans: int,
    min_followup_years: float
) -> pd.DataFrame:
    rows = []

    if "epoch_date" not in long_epoch.columns:
        raise ValueError("Longitudinal EPOCH file must contain a usable date column to compute slope.")

    for pid, g in long_epoch.groupby(ID_COL):
        g = g.dropna(subset=[value_col, "epoch_date"]).copy()
        g = g.sort_values("epoch_date")

        if len(g) < min_scans:
            continue

        t = (g["epoch_date"] - g["epoch_date"].min()).dt.days.values.astype(float) / 365.25
        y = g[value_col].values.astype(float)

        followup = np.nanmax(t) - np.nanmin(t)
        if not np.isfinite(followup) or followup < min_followup_years:
            continue

        if np.nanstd(y, ddof=1) <= 0:
            continue

        X = np.column_stack([np.ones(len(t)), t])
        coef = np.linalg.lstsq(X, y, rcond=None)[0]

        rows.append({
            ID_COL: pid,
            out_col: float(coef[1]),
            "AD_EPOCH_slope_intercept": float(coef[0]),
            "n_epoch_scans_for_slope": int(len(g)),
            "epoch_slope_followup_years": float(followup),
            "anchor_date": g["epoch_date"].max(),
            "baseline_epoch_date": g["epoch_date"].min(),
        })

    out = pd.DataFrame(rows)
    print(f"[INFO] Computed AD EPOCH slope per year for {len(out):,} participants")
    return out


# ============================================================
# Extract nearest ADNI variables within date window
# ============================================================

def nearest_visit_values(
    adni: pd.DataFrame,
    anchors: pd.DataFrame,
    variables: List[str],
    window_days: int
) -> pd.DataFrame:
    grouped = {
        pid: g.sort_values(DATE_COL).copy()
        for pid, g in adni.groupby(ID_COL)
    }

    rows = []

    for _, a in anchors.iterrows():
        pid = a[ID_COL]
        anchor_date = a["anchor_date"]

        row = {
            ID_COL: pid,
            "anchor_date": anchor_date,
        }

        g = grouped.get(pid, None)

        for var in variables:
            row[var] = np.nan
            row[f"{var}_date"] = pd.NaT
            row[f"{var}_days_from_anchor"] = np.nan

        if g is None or pd.isna(anchor_date):
            rows.append(row)
            continue

        for var in variables:
            gg = g[[DATE_COL, var]].copy()
            gg = gg[gg[var].notna()].copy()

            if len(gg) == 0:
                continue

            gg["abs_days"] = (gg[DATE_COL] - anchor_date).abs().dt.days
            gg = gg[gg["abs_days"] <= window_days].copy()

            if len(gg) == 0:
                continue

            gg = gg.sort_values(["abs_days", DATE_COL])
            best = gg.iloc[0]

            row[var] = best[var]
            row[f"{var}_date"] = best[DATE_COL]
            row[f"{var}_days_from_anchor"] = float(best["abs_days"])

        rows.append(row)

    out = pd.DataFrame(rows)
    return out


# ============================================================
# Regression and paired permutation test
# ============================================================

def build_covariate_matrix(df: pd.DataFrame, covariates: List[str]) -> Tuple[np.ndarray, List[str]]:
    mats = []
    names = []

    for c in covariates:
        if c not in df.columns:
            continue

        if c in ["Sex", "SITE"]:
            d = pd.get_dummies(df[c].astype(str), prefix=c, drop_first=True)
            if d.shape[1] > 0:
                mats.append(d.values.astype(float))
                names.extend(list(d.columns))
        else:
            x = numeric_clean(df[c])
            if x.notna().all() and x.nunique(dropna=True) > 1:
                z = zscore_array(x.values.astype(float))
                if z is not None:
                    mats.append(z.reshape(-1, 1))
                    names.append(c)

    if len(mats) == 0:
        return np.zeros((len(df), 0)), []

    return np.column_stack(mats), names


def ols_standardized_beta(
    data: pd.DataFrame,
    y_col: str,
    x_col: str,
    covariates: List[str],
    min_n: int
) -> Dict[str, object]:
    needed = [y_col, x_col] + covariates
    needed = [c for c in needed if c in data.columns]

    df = data[needed].copy()
    df[y_col] = numeric_clean(df[y_col])
    df[x_col] = numeric_clean(df[x_col])
    df = df.replace([np.inf, -np.inf], np.nan).dropna().copy()

    if len(df) < min_n:
        return {"status": "insufficient_n", "n": len(df)}

    y = zscore_array(df[y_col].values.astype(float))
    x = zscore_array(df[x_col].values.astype(float))

    if y is None or x is None:
        return {"status": "zero_variance", "n": len(df)}

    Xcov, cov_names = build_covariate_matrix(df, covariates)
    X = np.column_stack([np.ones(len(df)), x, Xcov])

    try:
        coef = np.linalg.lstsq(X, y, rcond=None)[0]
        yhat = X @ coef
        resid = y - yhat

        rank = np.linalg.matrix_rank(X)
        dfree = len(y) - rank

        if dfree <= 0:
            return {"status": "no_residual_df", "n": len(df)}

        sigma2 = float(np.sum(resid ** 2) / dfree)
        xtx_inv = np.linalg.pinv(X.T @ X)
        se = float(np.sqrt(sigma2 * xtx_inv[1, 1]))

        beta = float(coef[1])
        tval = beta / se if se > 0 else np.nan
        pval = float(2.0 * stats.t.sf(abs(tval), dfree)) if np.isfinite(tval) else np.nan

        ss_tot = float(np.sum((y - np.mean(y)) ** 2))
        ss_res = float(np.sum(resid ** 2))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

        return {
            "status": "ok",
            "n": int(len(df)),
            "standardized_beta": beta,
            "se": se,
            "t": float(tval),
            "p": pval,
            "r2": float(r2),
            "df": int(dfree),
            "covariate_terms": ",".join(cov_names),
        }

    except Exception as e:
        return {
            "status": f"ols_failed: {e}",
            "n": len(df),
        }


def residualize(df: pd.DataFrame, var_col: str, covariates: List[str]) -> Optional[np.ndarray]:
    y = numeric_clean(df[var_col]).values.astype(float)

    Xcov, _ = build_covariate_matrix(df, covariates)
    X = np.column_stack([np.ones(len(df)), Xcov])

    try:
        coef = np.linalg.lstsq(X, y, rcond=None)[0]
        resid = y - X @ coef
        return zscore_array(resid)
    except Exception:
        return None


def paired_permutation_abs_beta(
    data: pd.DataFrame,
    y_col: str,
    epoch_col: str,
    comparator_col: str,
    covariates: List[str],
    min_n: int,
    n_perm: int,
    seed: int
) -> Dict[str, object]:
    needed = [y_col, epoch_col, comparator_col] + covariates
    needed = [c for c in needed if c in data.columns]

    df = data[needed].copy()
    df[y_col] = numeric_clean(df[y_col])
    df[epoch_col] = numeric_clean(df[epoch_col])
    df[comparator_col] = numeric_clean(df[comparator_col])
    df = df.replace([np.inf, -np.inf], np.nan).dropna().copy()

    if len(df) < min_n:
        return {"status": "insufficient_n", "n": len(df)}

    y = residualize(df, y_col, covariates)
    xe = residualize(df, epoch_col, covariates)
    xm = residualize(df, comparator_col, covariates)

    if y is None or xe is None or xm is None:
        return {"status": "residualization_failed", "n": len(df)}

    def beta_1d(x, yy):
        denom = float(np.dot(x, x))
        if denom <= 0:
            return np.nan
        return float(np.dot(x, yy) / denom)

    beta_epoch = beta_1d(xe, y)
    beta_marker = beta_1d(xm, y)
    obs = abs(beta_epoch) - abs(beta_marker)

    rng = np.random.default_rng(seed)
    perm_stats = np.zeros(n_perm, dtype=float)

    for b in range(n_perm):
        swap = rng.random(len(df)) < 0.5

        xe_p = xe.copy()
        xm_p = xm.copy()

        xe_p[swap] = xm[swap]
        xm_p[swap] = xe[swap]

        be = beta_1d(xe_p, y)
        bm = beta_1d(xm_p, y)

        perm_stats[b] = abs(be) - abs(bm)

    p_perm = (np.sum(np.abs(perm_stats) >= abs(obs)) + 1.0) / (n_perm + 1.0)

    return {
        "status": "ok",
        "n": int(len(df)),
        "beta_epoch": float(beta_epoch),
        "beta_comparator": float(beta_marker),
        "delta_abs_beta_epoch_minus_comparator": float(obs),
        "p_perm_two_sided": float(p_perm),
    }


# ============================================================
# Main
# ============================================================

def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Compare ADNI AD EPOCH baseline/slope associations with cognition "
            "against SPARE and CSF biomarkers."
        )
    )

    p.add_argument("--adni_tsv", required=True, type=str)
    p.add_argument("--baseline_epoch", required=True, type=str)
    p.add_argument("--longitudinal_epoch", required=True, type=str)
    p.add_argument("--outdir", required=True, type=str)

    p.add_argument("--baseline_epoch_col", default="auto", type=str)
    p.add_argument("--longitudinal_epoch_col", default="auto", type=str)

    p.add_argument("--window_days", default=365, type=int)
    p.add_argument("--n_perm", default=10000, type=int)
    p.add_argument("--min_n", default=30, type=int)
    p.add_argument("--seed", default=2026, type=int)

    p.add_argument("--min_slope_scans", default=2, type=int)
    p.add_argument("--min_slope_followup_years", default=0.5, type=float)

    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print("[INFO] Reading ADNI iSTAGING table")
    adni = load_adni_istaging(args.adni_tsv)

    variables_to_extract = COGNITIVE_COLS + COMPARATOR_BIOMARKER_COLS + COVARIATE_COLS

    # ------------------------------------------------------------
    # Baseline AD EPOCH
    # ------------------------------------------------------------
    print("[INFO] Loading baseline AD EPOCH")
    baseline_epoch = load_epoch_scores(
        args.baseline_epoch,
        value_name="AD_EPOCH_baseline",
        requested_col=args.baseline_epoch_col,
        prefer_longitudinal=False,
    )

    baseline_epoch = add_anchor_date(baseline_epoch, adni)
    baseline_anchors = baseline_epoch[[ID_COL, "anchor_date", "AD_EPOCH_baseline"]].copy()

    print("[INFO] Extracting nearest ADNI measures for baseline analysis")
    baseline_nearest = nearest_visit_values(
        adni=adni,
        anchors=baseline_anchors[[ID_COL, "anchor_date"]],
        variables=variables_to_extract,
        window_days=args.window_days,
    )

    baseline_data = baseline_anchors.merge(
        baseline_nearest.drop(columns=["anchor_date"]),
        on=ID_COL,
        how="left",
    )

    baseline_data["analysis_type"] = "baseline_AD_EPOCH"

    # ------------------------------------------------------------
    # Longitudinal AD EPOCH slope
    # ------------------------------------------------------------
    print("[INFO] Loading longitudinal AD EPOCH")
    long_epoch = load_epoch_scores(
        args.longitudinal_epoch,
        value_name="AD_EPOCH_longitudinal",
        requested_col=args.longitudinal_epoch_col,
        prefer_longitudinal=True,
    )

    print("[INFO] Computing AD EPOCH slope per year")
    slope_epoch = compute_slope_per_year(
        long_epoch=long_epoch,
        value_col="AD_EPOCH_longitudinal",
        out_col="AD_EPOCH_slope_per_year",
        min_scans=args.min_slope_scans,
        min_followup_years=args.min_slope_followup_years,
    )

    print(f"[INFO] Baseline AD EPOCH anchor rows: {len(baseline_data):,}")
    print(f"[INFO] AD EPOCH slope rows: {len(slope_epoch):,}")

    if len(slope_epoch) > 0:
        print("[INFO] Extracting nearest ADNI measures for slope analysis")
        slope_nearest = nearest_visit_values(
            adni=adni,
            anchors=slope_epoch[[ID_COL, "anchor_date"]],
            variables=variables_to_extract,
            window_days=args.window_days,
        )

        slope_data = slope_epoch.merge(
            slope_nearest.drop(columns=["anchor_date"]),
            on=ID_COL,
            how="left",
        )

        slope_data["analysis_type"] = "AD_EPOCH_slope_per_year"
    else:
        slope_data = pd.DataFrame()

    # ------------------------------------------------------------
    # Save analysis-ready datasets
    # ------------------------------------------------------------
    baseline_dataset_path = os.path.join(args.outdir, "analysis_dataset_baseline_ad_epoch.tsv")
    baseline_data.to_csv(baseline_dataset_path, sep="\t", index=False)

    if len(slope_data) > 0:
        slope_dataset_path = os.path.join(args.outdir, "analysis_dataset_ad_epoch_slope.tsv")
        slope_data.to_csv(slope_dataset_path, sep="\t", index=False)

    # ------------------------------------------------------------
    # Association and permutation analyses
    # ------------------------------------------------------------
    association_rows = []
    permutation_rows = []

    analysis_sets = [
        ("baseline_AD_EPOCH", baseline_data, "AD_EPOCH_baseline"),
    ]

    if len(slope_data) > 0:
        analysis_sets.append(
            ("AD_EPOCH_slope_per_year", slope_data, "AD_EPOCH_slope_per_year")
        )

    for analysis_type, data, epoch_predictor in analysis_sets:
        predictors = [epoch_predictor] + COMPARATOR_BIOMARKER_COLS

        for cognitive_score in COGNITIVE_COLS:
            for predictor in predictors:
                res = ols_standardized_beta(
                    data=data,
                    y_col=cognitive_score,
                    x_col=predictor,
                    covariates=COVARIATE_COLS,
                    min_n=args.min_n,
                )

                association_rows.append({
                    "analysis_type": analysis_type,
                    "cognitive_score": cognitive_score,
                    "predictor": predictor,
                    "predictor_class": (
                        "AD_EPOCH" if predictor == epoch_predictor else "Comparator_biomarker"
                    ),
                    "n": res.get("n", np.nan),
                    "standardized_beta": res.get("standardized_beta", np.nan),
                    "se": res.get("se", np.nan),
                    "t": res.get("t", np.nan),
                    "p": res.get("p", np.nan),
                    "r2": res.get("r2", np.nan),
                    "df": res.get("df", np.nan),
                    "status": res.get("status", "unknown"),
                    "covariates": ",".join(COVARIATE_COLS),
                    "window_days": args.window_days,
                })

            for comparator in COMPARATOR_BIOMARKER_COLS:
                res = paired_permutation_abs_beta(
                    data=data,
                    y_col=cognitive_score,
                    epoch_col=epoch_predictor,
                    comparator_col=comparator,
                    covariates=COVARIATE_COLS,
                    min_n=args.min_n,
                    n_perm=args.n_perm,
                    seed=args.seed,
                )

                permutation_rows.append({
                    "analysis_type": analysis_type,
                    "cognitive_score": cognitive_score,
                    "epoch_predictor": epoch_predictor,
                    "comparator_biomarker": comparator,
                    "n": res.get("n", np.nan),
                    "beta_epoch": res.get("beta_epoch", np.nan),
                    "beta_comparator": res.get("beta_comparator", np.nan),
                    "delta_abs_beta_epoch_minus_comparator": res.get(
                        "delta_abs_beta_epoch_minus_comparator", np.nan
                    ),
                    "p_perm_two_sided": res.get("p_perm_two_sided", np.nan),
                    "n_perm": args.n_perm,
                    "status": res.get("status", "unknown"),
                    "covariates": ",".join(COVARIATE_COLS),
                    "window_days": args.window_days,
                })

    assoc = pd.DataFrame(association_rows)
    perm = pd.DataFrame(permutation_rows)

    assoc_path = os.path.join(args.outdir, "ad_epoch_biomarker_cognition_standardized_betas.tsv")
    perm_path = os.path.join(args.outdir, "ad_epoch_vs_biomarker_permutation_tests.tsv")

    assoc.to_csv(assoc_path, sep="\t", index=False)
    perm.to_csv(perm_path, sep="\t", index=False)

    # ------------------------------------------------------------
    # Column manifest
    # ------------------------------------------------------------
    manifest = pd.DataFrame({
        "column_type": (
            ["id"]
            + ["date"]
            + ["cognitive_score"] * len(COGNITIVE_COLS)
            + ["comparator_biomarker"] * len(COMPARATOR_BIOMARKER_COLS)
            + ["covariate"] * len(COVARIATE_COLS)
        ),
        "column_name": (
            [ID_COL]
            + [DATE_COL]
            + COGNITIVE_COLS
            + COMPARATOR_BIOMARKER_COLS
            + COVARIATE_COLS
        ),
    })

    manifest_path = os.path.join(args.outdir, "exact_column_manifest.tsv")
    manifest.to_csv(manifest_path, sep="\t", index=False)

    print("[DONE] Analysis-ready baseline dataset:")
    print(baseline_dataset_path)

    if len(slope_data) > 0:
        print("[DONE] Analysis-ready slope dataset:")
        print(slope_dataset_path)

    print("[DONE] Standardized beta table:")
    print(assoc_path)

    print("[DONE] Paired permutation comparison table:")
    print(perm_path)

    print("[DONE] Column manifest:")
    print(manifest_path)

    print("[DONE] Output directory:")
    print(args.outdir)


if __name__ == "__main__":
    main()