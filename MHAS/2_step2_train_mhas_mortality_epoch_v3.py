#!/usr/bin/env python3
"""
STEP 2 v3: Train phenotype-based MHAS mortality EPOCH clock

This version fixes the singular-matrix error observed in v1/v2:
  "Convergence halted due to matrix inversion problems. Suspicion is high collinearity."

Main change:
  The train-fold preprocessing now removes:
    1) columns that become constant after median imputation,
    2) exact duplicate columns,
    3) highly correlated columns above --correlation-threshold.

Input:
  /Users/hao/Dropbox/MHAS/step1_2001_mortality_epoch/
    mhas_2001_step1_model_input_primary_nondisease.tsv

Run:
  cd /Users/hao/Project/whole-body_clocks/MHAS
  /Users/hao/opt/anaconda3/envs/DNE/bin/python 2_step2_train_mhas_mortality_epoch_v3.py

Outputs:
  /Users/hao/Dropbox/MHAS/step2_mortality_epoch_model/
"""

import argparse
import pickle
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from sklearn.impute import SimpleImputer
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

try:
    from lifelines import CoxPHFitter
    from lifelines.utils import concordance_index
except Exception as e:
    raise ImportError(
        "This script requires lifelines. Install it with:\n"
        "  pip install lifelines\n"
        "or:\n"
        "  conda install -c conda-forge lifelines"
    ) from e


def log(msg: str) -> None:
    print(msg, flush=True)


def safe_mkdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def to_numeric_df(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for c in out.columns:
        out[c] = pd.to_numeric(out[c], errors="coerce")
    return out.astype("float64")


def finite_check(df: pd.DataFrame, label: str) -> None:
    arr = df.to_numpy(dtype=float)
    if not np.isfinite(arr).all():
        n_bad = int((~np.isfinite(arr)).sum())
        raise ValueError(f"{label} contains {n_bad} non-finite values after preprocessing.")


def harrell_cindex(time, event, score):
    """
    Larger score = higher risk. lifelines concordance_index expects larger predicted
    value = longer survival, so pass -score.
    """
    time = pd.Series(time).astype(float)
    event = pd.Series(event).astype(int)
    score = pd.Series(score).astype(float)
    mask = time.notna() & event.notna() & score.notna()
    if int(mask.sum()) < 10 or event[mask].nunique() < 2:
        return np.nan
    return float(concordance_index(time[mask], -score[mask], event[mask]))


def split_data(df, seed=20260721, train_frac=0.60, val_frac=0.20, test_frac=0.20):
    if not np.isclose(train_frac + val_frac + test_frac, 1.0):
        raise ValueError("train_frac + val_frac + test_frac must equal 1.")

    y = df["event_death"].astype(int)
    idx_all = np.arange(len(df))

    idx_train, idx_temp = train_test_split(
        idx_all,
        train_size=train_frac,
        random_state=seed,
        stratify=y,
    )

    temp = df.iloc[idx_temp]
    y_temp = temp["event_death"].astype(int)
    rel_val_frac = val_frac / (val_frac + test_frac)

    idx_val_rel, idx_test_rel = train_test_split(
        np.arange(len(temp)),
        train_size=rel_val_frac,
        random_state=seed + 1,
        stratify=y_temp,
    )

    idx_val = idx_temp[idx_val_rel]
    idx_test = idx_temp[idx_test_rel]

    split = pd.Series("unassigned", index=df.index)
    split.iloc[idx_train] = "train"
    split.iloc[idx_val] = "validation"
    split.iloc[idx_test] = "test"
    return split


def drop_highly_correlated_columns(X: pd.DataFrame, threshold: float = 0.90):
    """
    Greedy correlation filter using train data only.
    Keeps earlier columns and drops later columns with absolute correlation above threshold.
    """
    if X.shape[1] <= 1:
        return list(X.columns), []

    corr = X.corr().abs()
    upper = corr.where(np.triu(np.ones(corr.shape), k=1).astype(bool))

    to_drop = []
    reasons = []

    for col in upper.columns:
        high_corr = upper[col][upper[col] > threshold]
        if len(high_corr) > 0:
            ref = high_corr.idxmax()
            max_corr = float(high_corr.max())
            to_drop.append(col)
            reasons.append({
                "feature": col,
                "reason": f"high_correlation_abs_gt_{threshold}",
                "reference_feature": ref,
                "abs_correlation": max_corr,
            })

    keep = [c for c in X.columns if c not in set(to_drop)]
    return keep, reasons


def fit_preprocess(X_train, X_other, clip_quantiles=(0.005, 0.995), correlation_threshold=0.99):
    """
    Train-fold preprocessing:
      1) numeric conversion
      2) train-only all-missing / observed zero-variance removal
      3) train-only winsorization bounds
      4) train-only median imputation
      5) remove post-imputation constants
      6) remove exact duplicate columns
      7) remove highly correlated columns
      8) train-only z-scaling

    This is the key fix for lifelines singular-matrix convergence errors.
    """
    X_train = to_numeric_df(X_train)
    X_other = to_numeric_df(X_other)

    drop_records = []

    # Drop all-missing and observed zero-variance columns before imputation.
    nonmissing = X_train.notna().sum(axis=0)
    keep = nonmissing[nonmissing > 0].index.tolist()
    for c in X_train.columns:
        if c not in keep:
            drop_records.append({"feature": c, "reason": "all_missing_in_train", "reference_feature": "", "abs_correlation": np.nan})

    X_train = X_train[keep]
    X_other = X_other[keep]

    nunique = X_train.nunique(dropna=True)
    keep2 = nunique[nunique > 1].index.tolist()
    for c in X_train.columns:
        if c not in keep2:
            drop_records.append({"feature": c, "reason": "zero_observed_variance_in_train", "reference_feature": "", "abs_correlation": np.nan})

    X_train = X_train[keep2]
    X_other = X_other[keep2]

    initial_cols = list(X_train.columns)

    q_low, q_high = clip_quantiles
    lo = X_train.quantile(q_low, numeric_only=True)
    hi = X_train.quantile(q_high, numeric_only=True)

    X_train_clip = X_train.clip(lower=lo, upper=hi, axis=1).astype("float64")
    X_other_clip = X_other.clip(lower=lo, upper=hi, axis=1).astype("float64")

    imputer = SimpleImputer(strategy="median")
    X_train_imp = imputer.fit_transform(X_train_clip)
    X_other_imp = imputer.transform(X_other_clip)

    X_train_imp = pd.DataFrame(X_train_imp, columns=initial_cols, index=X_train.index).astype("float64")
    X_other_imp = pd.DataFrame(X_other_imp, columns=initial_cols, index=X_other.index).astype("float64")

    # Drop columns that become constant after imputation/winsorization.
    post_std = X_train_imp.std(axis=0, ddof=0)
    keep3 = post_std[post_std > 1e-10].index.tolist()
    for c in X_train_imp.columns:
        if c not in keep3:
            drop_records.append({"feature": c, "reason": "constant_after_imputation_or_winsorization", "reference_feature": "", "abs_correlation": np.nan})

    X_train_imp = X_train_imp[keep3]
    X_other_imp = X_other_imp[keep3]

    # Drop exact duplicates after imputation.
    duplicate_mask = X_train_imp.T.duplicated()
    duplicate_cols = list(X_train_imp.columns[duplicate_mask])
    if duplicate_cols:
        # Find a reference duplicate for reporting.
        for c in duplicate_cols:
            ref = ""
            for k in X_train_imp.columns:
                if k == c:
                    break
                if np.array_equal(X_train_imp[c].to_numpy(), X_train_imp[k].to_numpy()):
                    ref = k
                    break
            drop_records.append({"feature": c, "reason": "exact_duplicate_after_imputation", "reference_feature": ref, "abs_correlation": 1.0})

        X_train_imp = X_train_imp.loc[:, ~duplicate_mask]
        X_other_imp = X_other_imp[X_train_imp.columns]

    # Drop highly correlated columns.
    keep_corr, corr_reasons = drop_highly_correlated_columns(X_train_imp, threshold=correlation_threshold)
    drop_records.extend(corr_reasons)
    X_train_imp = X_train_imp[keep_corr]
    X_other_imp = X_other_imp[keep_corr]

    selected_cols = list(X_train_imp.columns)
    if len(selected_cols) == 0:
        raise ValueError("No feature columns remain after preprocessing.")

    scaler = StandardScaler()
    X_train_scaled = pd.DataFrame(
        scaler.fit_transform(X_train_imp),
        columns=selected_cols,
        index=X_train_imp.index,
    ).astype("float64")
    X_other_scaled = pd.DataFrame(
        scaler.transform(X_other_imp),
        columns=selected_cols,
        index=X_other_imp.index,
    ).astype("float64")

    finite_check(X_train_scaled, "X_train_scaled")
    finite_check(X_other_scaled, "X_other_scaled")

    preprocess = {
        "initial_feature_columns": initial_cols,
        "selected_feature_columns": selected_cols,
        "winsor_low": lo.to_dict(),
        "winsor_high": hi.to_dict(),
        "imputer": imputer,
        "scaler": scaler,
        "clip_quantiles": clip_quantiles,
        "correlation_threshold": correlation_threshold,
        "dropped_features": drop_records,
    }

    return X_train_scaled, X_other_scaled, preprocess


def apply_preprocess(X, preprocess):
    initial_cols = preprocess["initial_feature_columns"]
    selected_cols = preprocess["selected_feature_columns"]

    X = to_numeric_df(X.copy())

    for c in initial_cols:
        if c not in X.columns:
            X[c] = np.nan

    X = X[initial_cols]
    lo = pd.Series(preprocess["winsor_low"])
    hi = pd.Series(preprocess["winsor_high"])
    X_clip = X.clip(lower=lo, upper=hi, axis=1).astype("float64")

    X_imp = preprocess["imputer"].transform(X_clip)
    X_imp = pd.DataFrame(X_imp, columns=initial_cols, index=X.index).astype("float64")
    X_imp = X_imp[selected_cols]

    X_scaled = pd.DataFrame(
        preprocess["scaler"].transform(X_imp),
        columns=selected_cols,
        index=X.index,
    ).astype("float64")

    finite_check(X_scaled, "X_scaled")
    return X_scaled


def fit_cox_elastic_net(X, time, event, penalizer, l1_ratio, max_steps=200):
    dat = X.copy().astype("float64")
    dat["followup_years"] = np.asarray(time, dtype=float)
    dat["event_death"] = np.asarray(event, dtype=int)

    if not np.isfinite(dat.drop(columns=["event_death"]).to_numpy(dtype=float)).all():
        raise ValueError("Cox input contains non-finite values.")

    cph = CoxPHFitter(penalizer=penalizer, l1_ratio=l1_ratio)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            cph.fit(
                dat,
                duration_col="followup_years",
                event_col="event_death",
                show_progress=False,
                fit_options={"max_steps": max_steps},
            )
        except TypeError as e:
            # Compatibility with older lifelines that do not accept fit_options.
            if "fit_options" not in str(e):
                raise
            cph.fit(
                dat,
                duration_col="followup_years",
                event_col="event_death",
                show_progress=False,
            )

    return cph


def predict_lp(cph, X):
    ph = np.asarray(cph.predict_partial_hazard(X), dtype=float).reshape(-1)
    ph = np.clip(ph, 1e-300, 1e300)
    return np.log(ph)


def make_clinical_baseline_features(df):
    candidate_prefixes = [
        "age_2001",
        "sex_",
        "current_smoking_2001_",
        "ever_smoked_2001_",
        "bmi_self_report_2001",
        "bmi_measured_2001",
        "self_rated_health_2001_",
        "cesd_modified_score_2001",
        "adl_0_6_2001",
        "iadl_0_4_2001",
    ]

    cols = []
    exclude = {
        "participant_id", "baseline_date", "death_date",
        "followup_end_date", "event_death", "followup_years"
    }

    for c in df.columns:
        if c in exclude:
            continue
        for p in candidate_prefixes:
            if c == p or c.startswith(p):
                cols.append(c)
                break

    seen = set()
    return [c for c in cols if not (c in seen or seen.add(c))]


def fit_age_sex_residualizer(df_pred, train_mask):
    covars = []
    if "age_2001" in df_pred.columns:
        covars.append("age_2001")
    covars += [c for c in df_pred.columns if c.startswith("sex_")]

    if len(covars) == 0:
        mean_lp = float(df_pred.loc[train_mask, "lp_total"].astype(float).mean())
        return {"type": "intercept_only", "mean_lp": mean_lp, "covars": []}

    X_train = df_pred.loc[train_mask, covars].apply(pd.to_numeric, errors="coerce")
    y_train = df_pred.loc[train_mask, "lp_total"].astype(float)

    imp = SimpleImputer(strategy="median")
    X_train_imp = imp.fit_transform(X_train)

    lr = LinearRegression()
    lr.fit(X_train_imp, y_train)

    return {"type": "linear", "covars": covars, "imputer": imp, "model": lr}


def apply_residualizer(df_pred, residualizer):
    if residualizer["type"] == "intercept_only":
        return df_pred["lp_total"].astype(float) - residualizer["mean_lp"]

    covars = residualizer["covars"]
    X = df_pred[covars].apply(pd.to_numeric, errors="coerce")
    X_imp = residualizer["imputer"].transform(X)
    expected = residualizer["model"].predict(X_imp)
    return df_pred["lp_total"].astype(float) - expected


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        default="/Users/hao/Dropbox/MHAS/step1_2001_mortality_epoch/mhas_2001_step1_model_input_primary_nondisease.tsv",
        help="Step 1 model-input TSV."
    )
    parser.add_argument(
        "--out-dir",
        default="/Users/hao/Dropbox/MHAS/step2_mortality_epoch_model",
        help="Output directory."
    )
    parser.add_argument("--seed", type=int, default=20260721)
    parser.add_argument("--train-frac", type=float, default=0.60)
    parser.add_argument("--val-frac", type=float, default=0.20)
    parser.add_argument("--test-frac", type=float, default=0.20)
    parser.add_argument(
        "--penalizer-grid",
        default="0.001,0.003,0.01,0.03,0.1,0.3,1.0,3.0,10.0",
        help="Comma-separated lifelines penalizer grid."
    )
    parser.add_argument(
        "--l1-ratio-grid",
        default="0.0,0.25,0.5,0.75",
        help="Comma-separated l1_ratio grid. 0=ridge, 1=lasso."
    )
    parser.add_argument(
        "--min-feature-nonmissing",
        type=float,
        default=0.01,
        help="Drop feature columns with < this fraction nonmissing before train split."
    )
    parser.add_argument(
        "--correlation-threshold",
        type=float,
        default=0.99,
        help="Drop later columns with train-set absolute correlation above this threshold."
    )
    parser.add_argument("--min-events", type=int, default=50)
    parser.add_argument("--max-steps", type=int, default=200)
    parser.add_argument("--no-clinical-baseline", action="store_true")
    args = parser.parse_args()

    input_path = Path(args.input)
    out_dir = Path(args.out_dir)
    safe_mkdir(out_dir)

    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    log(f"Reading Step 1 model input: {input_path}")
    df = pd.read_csv(input_path, sep="\t", low_memory=False)

    required = ["participant_id", "event_death", "followup_years"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df["event_death"] = pd.to_numeric(df["event_death"], errors="coerce").astype("Int64")
    df["followup_years"] = pd.to_numeric(df["followup_years"], errors="coerce")
    df = df[df["event_death"].notna() & df["followup_years"].notna()].copy()
    df["event_death"] = df["event_death"].astype(int)
    df = df[df["followup_years"] > 0].copy().reset_index(drop=True)

    n_events = int(df["event_death"].sum())
    if n_events < args.min_events:
        raise ValueError(f"Too few mortality events for training: {n_events} < {args.min_events}")

    id_cols = [
        "participant_id", "baseline_date", "death_date",
        "followup_end_date", "event_death", "followup_years"
    ]
    feature_cols = [c for c in df.columns if c not in id_cols]

    X_all_raw = to_numeric_df(df[feature_cols])
    nonmiss_frac = X_all_raw.notna().mean(axis=0)
    feature_cols = nonmiss_frac[nonmiss_frac >= args.min_feature_nonmissing].index.tolist()
    X_all_raw = X_all_raw[feature_cols]

    log(f"Analytic N: {len(df):,}")
    log(f"Deaths: {n_events:,}")
    log(f"Candidate feature columns after missingness filter: {len(feature_cols):,}")

    df["split"] = split_data(
        df,
        seed=args.seed,
        train_frac=args.train_frac,
        val_frac=args.val_frac,
        test_frac=args.test_frac,
    )

    train_mask = df["split"] == "train"
    val_mask = df["split"] == "validation"
    trainval_mask = train_mask | val_mask

    split_summary = (
        df.groupby("split")
        .agg(
            n=("participant_id", "size"),
            deaths=("event_death", "sum"),
            median_followup_years=("followup_years", "median"),
        )
        .reset_index()
    )
    split_summary["event_rate"] = split_summary["deaths"] / split_summary["n"]

    X_train, _, preprocess = fit_preprocess(
        X_all_raw.loc[train_mask],
        X_all_raw.loc[val_mask],
        clip_quantiles=(0.005, 0.995),
        correlation_threshold=args.correlation_threshold,
    )
    X_val = apply_preprocess(X_all_raw.loc[val_mask], preprocess)
    X_all = apply_preprocess(X_all_raw, preprocess)

    dropped_features = pd.DataFrame(preprocess["dropped_features"])
    dropped_out = out_dir / "mhas_mortality_epoch_dropped_features_preprocessing.tsv"
    dropped_features.to_csv(dropped_out, sep="\t", index=False)

    log(f"Feature columns after collinearity-safe preprocessing: {len(preprocess['selected_feature_columns']):,}")
    log(f"Dropped features during preprocessing: {len(dropped_features):,}")

    time_train = df.loc[train_mask, "followup_years"].astype(float)
    event_train = df.loc[train_mask, "event_death"].astype(int)
    time_val = df.loc[val_mask, "followup_years"].astype(float)
    event_val = df.loc[val_mask, "event_death"].astype(int)

    penalizer_grid = [float(x) for x in args.penalizer_grid.split(",") if x.strip()]
    l1_ratio_grid = [float(x) for x in args.l1_ratio_grid.split(",") if x.strip()]

    tuning_rows = []
    best = None
    tuning_out = out_dir / "mhas_mortality_epoch_hyperparameter_tuning.tsv"

    log("Tuning elastic-net Cox hyperparameters...")
    for penalizer in penalizer_grid:
        for l1_ratio in l1_ratio_grid:
            row = {
                "penalizer": penalizer,
                "l1_ratio": l1_ratio,
                "train_cindex": np.nan,
                "validation_cindex": np.nan,
                "n_nonzero": np.nan,
                "status": "error",
                "error": "",
            }

            try:
                cph = fit_cox_elastic_net(
                    X_train,
                    time_train,
                    event_train,
                    penalizer=penalizer,
                    l1_ratio=l1_ratio,
                    max_steps=args.max_steps,
                )

                lp_train = predict_lp(cph, X_train)
                lp_val = predict_lp(cph, X_val)

                row["train_cindex"] = harrell_cindex(time_train, event_train, lp_train)
                row["validation_cindex"] = harrell_cindex(time_val, event_val, lp_val)
                row["n_nonzero"] = int((np.abs(cph.params_.values) > 1e-8).sum())
                row["status"] = "ok"

                if best is None or (
                    pd.notna(row["validation_cindex"])
                    and row["validation_cindex"] > best["validation_cindex"]
                ):
                    best = {
                        "penalizer": penalizer,
                        "l1_ratio": l1_ratio,
                        "validation_cindex": row["validation_cindex"],
                        "model": cph,
                    }

            except Exception as e:
                row["error"] = str(e).replace("\n", " ")[:1000]

            tuning_rows.append(row)
            err_short = f" error={row['error'][:180]}" if row["status"] != "ok" else ""
            log(
                f"  penalizer={penalizer:<7g} l1_ratio={l1_ratio:<4g} "
                f"val_C={row['validation_cindex']} status={row['status']}{err_short}"
            )
            pd.DataFrame(tuning_rows).to_csv(tuning_out, sep="\t", index=False)

    tuning = pd.DataFrame(tuning_rows)
    if best is None:
        first_err = tuning.loc[tuning["error"].astype(str).str.len() > 0, "error"]
        first_err = first_err.iloc[0] if len(first_err) else "No error message captured."
        raise RuntimeError(
            "All Cox hyperparameter fits failed even after collinearity filtering. "
            f"First captured error: {first_err}\n"
            f"Tuning table saved to: {tuning_out}\n"
            f"Dropped-feature table saved to: {dropped_out}"
        )

    best_penalizer = best["penalizer"]
    best_l1_ratio = best["l1_ratio"]

    log(
        f"Best hyperparameters: penalizer={best_penalizer}, "
        f"l1_ratio={best_l1_ratio}, validation C={best['validation_cindex']:.4f}"
    )

    # Final model on train + validation, using train-only preprocessing for conservative test evaluation.
    X_trainval = X_all.loc[trainval_mask]

    final_cph = fit_cox_elastic_net(
        X_trainval,
        df.loc[trainval_mask, "followup_years"].astype(float),
        df.loc[trainval_mask, "event_death"].astype(int),
        penalizer=best_penalizer,
        l1_ratio=best_l1_ratio,
        max_steps=args.max_steps,
    )

    pred = df[id_cols + ["split"]].copy()
    pred["lp_total"] = predict_lp(final_cph, X_all)

    # Acceleration residualizes lp_total on age and sex terms using training split only.
    pred_for_resid = pd.concat([pred.reset_index(drop=True), X_all_raw.reset_index(drop=True)], axis=1)
    residualizer = fit_age_sex_residualizer(pred_for_resid, train_mask.reset_index(drop=True))
    pred["mortality_epoch_acceleration"] = apply_residualizer(pred_for_resid, residualizer)

    train_acc = pred.loc[train_mask, "mortality_epoch_acceleration"].astype(float)
    acc_mean = float(train_acc.mean())
    acc_sd = float(train_acc.std(ddof=0))
    if acc_sd == 0 or pd.isna(acc_sd):
        acc_sd = 1.0
    pred["mortality_epoch_acceleration_z"] = (
        pred["mortality_epoch_acceleration"] - acc_mean
    ) / acc_sd

    q = pred.loc[train_mask, "mortality_epoch_acceleration_z"].quantile([0.25, 0.50, 0.75]).values
    pred["mortality_epoch_quartile"] = pd.cut(
        pred["mortality_epoch_acceleration_z"],
        bins=[-np.inf, q[0], q[1], q[2], np.inf],
        labels=["Q1_lowest", "Q2", "Q3", "Q4_highest"],
        include_lowest=True,
    ).astype(str)

    performance_rows = []
    for split_name in ["train", "validation", "test", "all"]:
        mask = np.ones(len(df), dtype=bool) if split_name == "all" else (df["split"] == split_name).values
        performance_rows.append({
            "split": split_name,
            "n": int(mask.sum()),
            "deaths": int(df.loc[mask, "event_death"].sum()),
            "event_rate": float(df.loc[mask, "event_death"].mean()),
            "median_followup_years": float(df.loc[mask, "followup_years"].median()),
            "cindex_lp_total": harrell_cindex(
                df.loc[mask, "followup_years"],
                df.loc[mask, "event_death"],
                pred.loc[mask, "lp_total"],
            ),
            "cindex_acceleration_z": harrell_cindex(
                df.loc[mask, "followup_years"],
                df.loc[mask, "event_death"],
                pred.loc[mask, "mortality_epoch_acceleration_z"],
            ),
        })

    performance = pd.DataFrame(performance_rows)

    clinical_cols = []
    if not args.no_clinical_baseline:
        clinical_cols = make_clinical_baseline_features(X_all_raw)
        if len(clinical_cols) > 0:
            log(f"Fitting clinical baseline comparator with {len(clinical_cols)} columns...")

            Xc_train, _, c_preprocess = fit_preprocess(
                X_all_raw.loc[train_mask, clinical_cols],
                X_all_raw.loc[val_mask, clinical_cols],
                clip_quantiles=(0.005, 0.995),
                correlation_threshold=args.correlation_threshold,
            )
            Xc_all = apply_preprocess(X_all_raw[clinical_cols], c_preprocess)
            Xc_trainval = Xc_all.loc[trainval_mask]

            clinical_cph = fit_cox_elastic_net(
                Xc_trainval,
                df.loc[trainval_mask, "followup_years"].astype(float),
                df.loc[trainval_mask, "event_death"].astype(int),
                penalizer=0.01,
                l1_ratio=0.0,
                max_steps=args.max_steps,
            )
            pred["clinical_baseline_lp"] = predict_lp(clinical_cph, Xc_all)

            clinical_rows = []
            for split_name in ["train", "validation", "test", "all"]:
                mask = np.ones(len(df), dtype=bool) if split_name == "all" else (df["split"] == split_name).values
                clinical_rows.append({
                    "split": split_name,
                    "cindex_clinical_baseline": harrell_cindex(
                        df.loc[mask, "followup_years"],
                        df.loc[mask, "event_death"],
                        pred.loc[mask, "clinical_baseline_lp"],
                    )
                })
            performance = performance.merge(pd.DataFrame(clinical_rows), on="split", how="left")
        else:
            log("No clinical baseline columns found; skipping comparator.")

    coef = final_cph.params_.reset_index()
    coef.columns = ["feature", "coef"]
    coef["abs_coef"] = coef["coef"].abs()
    coef["selected_nonzero"] = coef["abs_coef"] > 1e-8
    coef = coef.sort_values("abs_coef", ascending=False)
    selected = coef[coef["selected_nonzero"]].copy()

    pred_out = out_dir / "mhas_mortality_epoch_predictions.tsv"
    coef_out = out_dir / "mhas_mortality_epoch_coefficients.tsv"
    selected_out = out_dir / "mhas_mortality_epoch_selected_features.tsv"
    perf_out = out_dir / "mhas_mortality_epoch_performance.tsv"
    split_out = out_dir / "mhas_mortality_epoch_split_assignments.tsv"
    model_out = out_dir / "mhas_mortality_epoch_model.pkl"
    prep_out = out_dir / "mhas_mortality_epoch_preprocessing.pkl"
    audit_out = out_dir / "mhas_mortality_epoch_audit.txt"

    pred.to_csv(pred_out, sep="\t", index=False)
    coef.to_csv(coef_out, sep="\t", index=False)
    selected.to_csv(selected_out, sep="\t", index=False)
    performance.to_csv(perf_out, sep="\t", index=False)
    df[["participant_id", "split", "event_death", "followup_years"]].to_csv(split_out, sep="\t", index=False)

    with open(model_out, "wb") as f:
        pickle.dump(final_cph, f)

    with open(prep_out, "wb") as f:
        pickle.dump({
            "preprocess": preprocess,
            "residualizer": residualizer,
            "acceleration_train_mean": acc_mean,
            "acceleration_train_sd": acc_sd,
            "quartile_cutpoints_train": q.tolist(),
            "feature_columns_raw": feature_cols,
            "clinical_baseline_columns": clinical_cols,
            "best_penalizer": best_penalizer,
            "best_l1_ratio": best_l1_ratio,
            "seed": args.seed,
        }, f)

    audit = f"""MHAS STEP 2 v3: phenotype-based mortality EPOCH training

Input
-----
{input_path}

Output directory
----------------
{out_dir}

Analytic cohort
---------------
N: {len(df):,}
Deaths: {n_events:,}
Event rate: {n_events / len(df):.4f}
Median follow-up years: {df["followup_years"].median():.2f}

Split summary
-------------
{split_summary.to_string(index=False)}

Feature preprocessing
---------------------
Raw candidate feature columns after missingness filter: {len(feature_cols):,}
Feature columns after collinearity-safe preprocessing: {len(preprocess["selected_feature_columns"]):,}
Dropped features during preprocessing: {len(dropped_features):,}
Correlation threshold: {args.correlation_threshold}

Model
-----
Elastic-net Cox proportional hazards model
Selected penalizer: {best_penalizer}
Selected l1_ratio: {best_l1_ratio}
Nonzero coefficients: {int(selected.shape[0])}

Score definitions
-----------------
lp_total:
  Cox linear predictor using the selected primary non-disease phenotype features.

mortality_epoch_acceleration:
  residual from lp_total ~ age_2001 + sex terms.
  The residualization model is fit in the training split only.

mortality_epoch_acceleration_z:
  mortality_epoch_acceleration standardized using training-split mean/SD.

Performance
-----------
{performance.to_string(index=False)}

Clinical baseline comparator
----------------------------
Clinical baseline columns used: {len(clinical_cols)}
{", ".join(clinical_cols[:80]) if clinical_cols else "None"}

Output files
------------
Predictions:
{pred_out}

Coefficients:
{coef_out}

Selected features:
{selected_out}

Performance:
{perf_out}

Hyperparameter tuning:
{tuning_out}

Dropped features:
{dropped_out}

Model pickle:
{model_out}

Preprocessing pickle:
{prep_out}
"""
    audit_out.write_text(audit)
    log("\n" + audit)
    log("STEP 2 finished successfully.")


if __name__ == "__main__":
    main()
