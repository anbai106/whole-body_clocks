#!/usr/bin/env python3
"""
6_demographic_table.py

Create Supplementary Table 3 demographic rows for ADNI, NHANES, and MHAS only.

Rows generated:
  1. ADNI   - Brain MRI
  2. NHANES - Survey and lab biomarkers & mortality
  3. MHAS   - Survey and lab biomarkers & mortality

For each dataset, the script reports:
  N
  Age mean+/-SD
  Female count / female percentage

Default sources:
  ADNI remote:
    wenju@cubic-login5:~/Reproducibile_paper/WholeBodyClock/adni_lepoch/results_brain_mri_ad_lepoch/adni_brain_mri_ad_lepoch_predictions.tsv

  NHANES local:
    /Users/hao/Dropbox/NHANES/output_dir/model2_nondisease_mortality_epoch/nhanes_model2_analysis_table.tsv.gz

  MHAS local:
    /Users/hao/Dropbox/MHAS/step1_2001_mortality_epoch/mhas_2001_step1_clean_analytic_cohort.tsv

Recommended run:
  cd /Users/hao/Project/whole-body_clocks
  python 6_demographic_table.py --out-dir /Users/hao/Dropbox/2026_EPOCH/Supplementary_Tables
"""

import argparse
import os
import re
import subprocess
from pathlib import Path
from typing import Optional, List, Tuple, Dict

import numpy as np
import pandas as pd

PM = "\u00b1"


def log(msg: str) -> None:
    print(msg, flush=True)


def expand_path(x: Optional[str]) -> Optional[Path]:
    if x is None or str(x).strip() == "":
        return None
    return Path(os.path.expanduser(str(x)))


def read_table(path: Path, nrows=None) -> pd.DataFrame:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    name = path.name.lower()
    suffix = path.suffix.lower()
    if name.endswith(".tsv.gz") or name.endswith(".txt.gz"):
        return pd.read_csv(path, sep="\t", compression="gzip", low_memory=False, nrows=nrows)
    if name.endswith(".csv.gz"):
        return pd.read_csv(path, compression="gzip", low_memory=False, nrows=nrows)
    if suffix == ".csv":
        return pd.read_csv(path, low_memory=False, nrows=nrows)
    if suffix in [".tsv", ".txt"]:
        return pd.read_csv(path, sep="\t", low_memory=False, nrows=nrows)
    if suffix == ".dta":
        return pd.read_stata(path, convert_categoricals=False)
    if suffix == ".xpt":
        return pd.read_sas(path, format="xport")
    raise ValueError(f"Unsupported file type: {path}")


def safe_num(s) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def clean_col(x: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(x).lower())


def format_n(n) -> str:
    return "TODO" if pd.isna(n) else f"{int(n):,}"


def format_age(mean, sd) -> str:
    if pd.isna(mean) or pd.isna(sd):
        return f"TODO{PM}TODO"
    return f"{mean:.2f}{PM}{sd:.2f}"


def format_female(n_female, pct_female) -> str:
    if pd.isna(n_female) or pd.isna(pct_female):
        return "TODO/TODO %"
    return f"{int(n_female):,}/{pct_female:.0f}%"


def first_existing_col(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    cols = list(df.columns)
    col_map = {clean_col(c): c for c in cols}
    for cand in candidates:
        if cand in cols:
            return cand
        cc = clean_col(cand)
        if cc in col_map:
            return col_map[cc]
    return None


def detect_id_col(df: pd.DataFrame, dataset: str, user_col: Optional[str] = None) -> Optional[str]:
    if user_col and user_col in df.columns:
        return user_col
    dataset = dataset.lower()
    if dataset == "adni":
        candidates = ["PTID", "RID", "participant_id", "subject_id", "id", "ID"]
    elif dataset == "nhanes":
        candidates = ["SEQN", "participant_id", "id", "ID"]
    elif dataset == "mhas":
        candidates = ["participant_id", "unhhidnp", "rahhidnp", "participant_id_raw", "id", "ID"]
    else:
        candidates = ["participant_id", "PTID", "RID", "SEQN", "id", "ID"]
    return first_existing_col(df, candidates)


def detect_age_col(df: pd.DataFrame, dataset: str, user_col: Optional[str] = None) -> Optional[str]:
    if user_col and user_col in df.columns:
        return user_col
    dataset = dataset.lower()
    if dataset == "adni":
        candidates = ["Age", "AGE", "age", "age_years", "baseline_age"]
    elif dataset == "nhanes":
        candidates = ["RIDAGEYR", "age_at_exam", "Age", "AGE", "age", "age_years", "baseline_age"]
    elif dataset == "mhas":
        candidates = ["age_2001", "r1agey", "Age", "AGE", "age", "age_years"]
    else:
        candidates = ["Age", "AGE", "age", "age_years"]
    return first_existing_col(df, candidates)


def detect_sex_col(df: pd.DataFrame, dataset: str, user_col: Optional[str] = None) -> Optional[str]:
    if user_col and user_col in df.columns:
        return user_col
    dataset = dataset.lower()
    if dataset == "adni":
        candidates = ["Sex", "SEX", "sex", "gender", "Gender"]
    elif dataset == "nhanes":
        candidates = ["RIAGENDR", "sex", "Sex", "gender", "Gender"]
    elif dataset == "mhas":
        candidates = ["sex", "ragender", "gender", "Gender", "Sex"]
    else:
        candidates = ["sex", "Sex", "gender", "Gender"]

    col = first_existing_col(df, candidates)
    if col is not None:
        return col

    for c in df.columns:
        lc = str(c).lower()
        if lc in [
            "female", "sex_female", "gender_female",
            "riagendr_2", "riagendr_2.0", "sex_2", "sex_2.0",
            "ragender_2", "ragender_2.0"
        ]:
            return c
        if "female" in lc and ("sex" in lc or "gender" in lc):
            return c
    return None


def female_indicator_from_column(s: pd.Series, dataset: str, col_name: str) -> pd.Series:
    dataset = dataset.lower()
    col_l = str(col_name).lower()

    if (
        col_l in [
            "female", "sex_female", "gender_female",
            "riagendr_2", "riagendr_2.0", "sex_2", "sex_2.0",
            "ragender_2", "ragender_2.0"
        ]
        or ("female" in col_l and ("sex" in col_l or "gender" in col_l))
    ):
        x = safe_num(s)
        out = pd.Series(np.nan, index=s.index, dtype="float64")
        out[x == 1] = 1
        out[x == 0] = 0
        return out

    txt = s.astype("object").astype(str).str.strip().str.lower()
    txt = txt.replace({"nan": "", "none": "", "na": "", "n/a": "", "<na>": "", "missing": ""})
    out = pd.Series(np.nan, index=s.index, dtype="float64")

    if dataset == "adni":
        female_values = {"female", "f", "woman", "women", "0", "0.0"}
        male_values = {"male", "m", "man", "men", "1", "1.0"}
    else:
        female_values = {"female", "f", "woman", "women", "2", "2.0"}
        male_values = {"male", "m", "man", "men", "1", "1.0"}

    out[txt.isin(female_values)] = 1
    out[txt.isin(male_values)] = 0

    x = safe_num(s)
    if dataset == "adni":
        out[(out.isna()) & (x == 0)] = 1
        out[(out.isna()) & (x == 1)] = 0
    else:
        out[(out.isna()) & (x == 2)] = 1
        out[(out.isna()) & (x == 1)] = 0
    return out


def deduplicate_participants(df: pd.DataFrame, id_col: Optional[str]) -> pd.DataFrame:
    if id_col is None or id_col not in df.columns:
        return df.copy()
    out = df.copy()
    out[id_col] = out[id_col].astype(str)
    out = out[out[id_col].notna() & (out[id_col] != "nan")].copy()
    return out.drop_duplicates(subset=[id_col], keep="first")


def summarize_demographics(
    file_path: Optional[Path],
    dataset: str,
    age_col: Optional[str] = None,
    sex_col: Optional[str] = None,
    id_col: Optional[str] = None,
) -> Tuple[Dict, str]:
    if file_path is None:
        summary = {
            "N_value": np.nan, "age_mean": np.nan, "age_sd": np.nan,
            "female_n": np.nan, "female_pct": np.nan, "source_file": "",
            "status": "TODO: input file not found or not provided",
            "id_col": "", "age_col": "", "sex_col": "",
        }
        return summary, f"{dataset}: no input file found/provided."

    try:
        df = read_table(file_path)
    except Exception as e:
        summary = {
            "N_value": np.nan, "age_mean": np.nan, "age_sd": np.nan,
            "female_n": np.nan, "female_pct": np.nan, "source_file": str(file_path),
            "status": f"ERROR: failed to read file: {e}",
            "id_col": "", "age_col": "", "sex_col": "",
        }
        return summary, f"{dataset}: failed to read {file_path}: {e}"

    id_c = detect_id_col(df, dataset, id_col)
    df0 = deduplicate_participants(df, id_c)
    age_c = detect_age_col(df0, dataset, age_col)
    sex_c = detect_sex_col(df0, dataset, sex_col)

    status_parts = []
    if age_c is None:
        status_parts.append("age column not detected")
    if sex_c is None:
        status_parts.append("sex/female column not detected")

    age = safe_num(df0[age_c]) if age_c else pd.Series(np.nan, index=df0.index)
    female = female_indicator_from_column(df0[sex_c], dataset, sex_c) if sex_c else pd.Series(np.nan, index=df0.index)

    n = len(df0)
    age_mean = float(age.mean()) if age.notna().sum() > 0 else np.nan
    age_sd = float(age.std(ddof=1)) if age.notna().sum() > 1 else np.nan
    female_n = int((female == 1).sum()) if female.notna().sum() > 0 else np.nan
    female_pct = float((female == 1).mean() * 100.0) if female.notna().sum() > 0 else np.nan

    if not status_parts:
        status_parts.append("ok")

    summary = {
        "N_value": n,
        "age_mean": age_mean,
        "age_sd": age_sd,
        "female_n": female_n,
        "female_pct": female_pct,
        "source_file": str(file_path),
        "status": "; ".join(status_parts),
        "id_col": id_c or "",
        "age_col": age_c or "",
        "sex_col": sex_c or "",
    }

    audit = (
        f"{dataset}: file={file_path}\n"
        f"  rows after participant deduplication: {n}\n"
        f"  id_col={id_c}\n"
        f"  age_col={age_c}, nonmissing_age={int(age.notna().sum())}\n"
        f"  sex_col={sex_c}, nonmissing_sex={int(female.notna().sum())}\n"
        f"  female_n={female_n}, female_pct={female_pct}\n"
        f"  status={summary['status']}"
    )
    return summary, audit


def choose_first_existing(paths: List[str]) -> Optional[Path]:
    for x in paths:
        p = expand_path(x)
        if p and p.exists():
            return p
    return None


def quick_header(path: Path) -> List[str]:
    try:
        return list(read_table(path, nrows=0).columns)
    except Exception:
        return []


def file_has_demo_columns(path: Path, dataset: str) -> int:
    cols = quick_header(path)
    if not cols:
        return -999
    dummy = pd.DataFrame(columns=cols)
    score = 0
    if detect_id_col(dummy, dataset) is not None:
        score += 5
    if detect_age_col(dummy, dataset) is not None:
        score += 10
    if detect_sex_col(dummy, dataset) is not None:
        score += 10

    name = path.name.lower()
    if "analysis_table" in name:
        score += 20
    if "predictions" in name:
        score += 12
    if "survival_dataset" in name:
        score += 10
    if "epoch_scores" in name:
        score += 6
    if "model_input" in name or "clean_analytic" in name:
        score += 6
    if "performance" in name or "summary" in name or "coefficient" in name:
        score -= 30
    if "old" in name or "backup" in name:
        score -= 20
    return score


def find_best_file(root: str, dataset: str, patterns: List[str]) -> Optional[Path]:
    root_path = expand_path(root)
    if root_path is None or not root_path.exists():
        return None

    candidates = []
    for pat in patterns:
        candidates.extend(root_path.glob(pat))

    candidates = [
        p for p in candidates
        if p.is_file()
        and (
            p.name.lower().endswith(".tsv.gz")
            or p.name.lower().endswith(".csv.gz")
            or p.suffix.lower() in [".tsv", ".csv", ".txt", ".dta", ".xpt"]
        )
    ]

    scored = []
    for p in candidates:
        score = file_has_demo_columns(p, dataset)
        if score > 0:
            scored.append((score, str(p), p))

    if not scored:
        return None
    scored.sort(key=lambda z: (-z[0], z[1]))
    return scored[0][2]


def scp_from_remote(remote: str, local: Path) -> bool:
    local.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["scp", remote, str(local)]
    log("Running: " + " ".join(cmd))
    try:
        subprocess.run(cmd, check=True)
        return local.exists()
    except Exception as e:
        log(f"WARNING: scp failed for {remote}: {e}")
        return False


def get_adni_file(args) -> Optional[Path]:
    if args.adni_file:
        p = expand_path(args.adni_file)
        if p and p.exists():
            return p
        raise FileNotFoundError(f"--adni-file does not exist: {args.adni_file}")

    local_cache = expand_path(args.adni_local_cache)
    if local_cache and local_cache.exists():
        return local_cache

    direct = choose_first_existing([
        "~/Reproducibile_paper/WholeBodyClock/adni_lepoch/results_brain_mri_ad_lepoch/adni_brain_mri_ad_lepoch_predictions.tsv",
        "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch/results_brain_mri_ad_lepoch/adni_brain_mri_ad_lepoch_predictions.tsv",
    ])
    if direct:
        return direct

    if not args.no_download_adni_from_cubic:
        remotes = [args.adni_remote]
        if args.adni_remote_fallback:
            remotes.append(args.adni_remote_fallback)
        for remote in remotes:
            if scp_from_remote(remote, local_cache):
                return local_cache

    return None


def get_nhanes_file(args) -> Optional[Path]:
    if args.nhanes_file:
        p = expand_path(args.nhanes_file)
        if p and p.exists():
            return p
        raise FileNotFoundError(f"--nhanes-file does not exist: {args.nhanes_file}")

    direct = choose_first_existing([
        "/Users/hao/Dropbox/NHANES/output_dir/model2_nondisease_mortality_epoch/nhanes_model2_analysis_table.tsv.gz",
        "/Users/hao/Dropbox/NHANES/output_dir/model2_nondisease_mortality_epoch/nhanes_model2_epoch_scores.tsv",
    ])
    if direct and file_has_demo_columns(direct, "nhanes") > 0:
        return direct

    return find_best_file(
        "/Users/hao/Dropbox/NHANES/output_dir/model2_nondisease_mortality_epoch",
        "nhanes",
        [
            "nhanes_model2_analysis_table.tsv.gz",
            "*analysis_table*.tsv.gz",
            "*analysis_table*.tsv",
            "*epoch_scores*.tsv",
            "*predictions*.tsv",
            "*.csv",
            "*.csv.gz",
        ],
    )


def get_mhas_file(args) -> Optional[Path]:
    if args.mhas_file:
        p = expand_path(args.mhas_file)
        if p and p.exists():
            return p
        raise FileNotFoundError(f"--mhas-file does not exist: {args.mhas_file}")

    direct = choose_first_existing([
        "/Users/hao/Dropbox/MHAS/step1_2001_mortality_epoch/mhas_2001_step1_clean_analytic_cohort.tsv",
        "/Users/hao/Dropbox/MHAS/step1_2001_mortality_epoch/mhas_2001_step1_model_input_primary_nondisease.tsv",
        "/Users/hao/Dropbox/MHAS/step2_mortality_epoch_model/mhas_mortality_epoch_predictions.tsv",
    ])
    if direct:
        return direct

    return find_best_file(
        "/Users/hao/Dropbox/MHAS",
        "mhas",
        ["**/*clean_analytic*.tsv", "**/*model_input*.tsv", "**/*predictions*.tsv", "**/*.csv"],
    )


def build_table(adni, nhanes, mhas) -> pd.DataFrame:
    return pd.DataFrame([
        {
            "Data type": "Individual",
            "BAG/Omics": "Brain MRI",
            "Study": "ADNI",
            "Country": "USA/Canada",
            "N": format_n(adni["N_value"]),
            "Age": format_age(adni["age_mean"], adni["age_sd"]),
            "Sex (female)": format_female(adni["female_n"], adni["female_pct"]),
        },
        {
            "Data type": "Individual",
            "BAG/Omics": "Survey and lab biomarkers & mortality",
            "Study": "NHANES",
            "Country": "USA",
            "N": format_n(nhanes["N_value"]),
            "Age": format_age(nhanes["age_mean"], nhanes["age_sd"]),
            "Sex (female)": format_female(nhanes["female_n"], nhanes["female_pct"]),
        },
        {
            "Data type": "Individual",
            "BAG/Omics": "Survey and lab biomarkers & mortality",
            "Study": "MHAS",
            "Country": "Mexico",
            "N": format_n(mhas["N_value"]),
            "Age": format_age(mhas["age_mean"], mhas["age_sd"]),
            "Sex (female)": format_female(mhas["female_n"], mhas["female_pct"]),
        },
    ])


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--adni-file", default=None, help="Local ADNI predictions/survival_dataset TSV/CSV.")
    parser.add_argument(
        "--adni-remote",
        default="wenju@cubic-login5:~/Reproducibile_paper/WholeBodyClock/adni_lepoch/results_brain_mri_ad_lepoch/adni_brain_mri_ad_lepoch_predictions.tsv",
        help="Primary remote CUBIC ADNI predictions file for scp."
    )
    parser.add_argument(
        "--adni-remote-fallback",
        default="wenju@cubic-login.uphs.upenn.edu:/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/adni_lepoch/results_brain_mri_ad_lepoch/adni_brain_mri_ad_lepoch_predictions.tsv",
        help="Fallback remote CUBIC ADNI predictions file for scp."
    )
    parser.add_argument(
        "--adni-local-cache",
        default="/Users/hao/Dropbox/2026_EPOCH/Supplementary_Tables/adni_brain_mri_ad_lepoch_predictions.tsv",
        help="Local cache location for ADNI file copied from CUBIC."
    )
    parser.add_argument(
        "--no-download-adni-from-cubic",
        action="store_true",
        help="Do not attempt to copy ADNI from CUBIC if local ADNI file is missing."
    )

    parser.add_argument("--nhanes-file", default=None, help="Local NHANES row-level TSV/CSV/GZ.")
    parser.add_argument("--mhas-file", default=None, help="Local MHAS row-level TSV/CSV.")

    parser.add_argument("--adni-age-col", default=None)
    parser.add_argument("--adni-sex-col", default=None)
    parser.add_argument("--adni-id-col", default=None)

    parser.add_argument("--nhanes-age-col", default=None)
    parser.add_argument("--nhanes-sex-col", default=None)
    parser.add_argument("--nhanes-id-col", default=None)

    parser.add_argument("--mhas-age-col", default=None)
    parser.add_argument("--mhas-sex-col", default=None)
    parser.add_argument("--mhas-id-col", default=None)

    parser.add_argument("--out-dir", default="/Users/hao/Dropbox/2026_EPOCH/Supplementary_Tables")
    parser.add_argument("--strict", action="store_true")

    args = parser.parse_args()

    out_dir = expand_path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    adni_path = get_adni_file(args)
    nhanes_path = get_nhanes_file(args)
    mhas_path = get_mhas_file(args)

    log(f"Resolved ADNI file: {adni_path}")
    log(f"Resolved NHANES file: {nhanes_path}")
    log(f"Resolved MHAS file: {mhas_path}")

    adni, adni_audit = summarize_demographics(
        adni_path, "adni",
        age_col=args.adni_age_col,
        sex_col=args.adni_sex_col,
        id_col=args.adni_id_col,
    )
    nhanes, nhanes_audit = summarize_demographics(
        nhanes_path, "nhanes",
        age_col=args.nhanes_age_col,
        sex_col=args.nhanes_sex_col,
        id_col=args.nhanes_id_col,
    )
    mhas, mhas_audit = summarize_demographics(
        mhas_path, "mhas",
        age_col=args.mhas_age_col,
        sex_col=args.mhas_sex_col,
        id_col=args.mhas_id_col,
    )

    if args.strict:
        bad = []
        for name, summary in [("ADNI", adni), ("NHANES", nhanes), ("MHAS", mhas)]:
            if summary["status"] != "ok":
                bad.append(f"{name}: {summary['status']}")
        if bad:
            raise RuntimeError("Strict mode failed:\n" + "\n".join(bad))

    table = build_table(adni, nhanes, mhas)

    tsv_out = out_dir / "Supplementary_Table_3_ADNI_NHANES_MHAS_demographics.tsv"
    csv_out = out_dir / "Supplementary_Table_3_ADNI_NHANES_MHAS_demographics.csv"
    md_out = out_dir / "Supplementary_Table_3_ADNI_NHANES_MHAS_demographics.md"
    audit_out = out_dir / "Supplementary_Table_3_ADNI_NHANES_MHAS_demographics_audit.txt"

    table.to_csv(tsv_out, sep="\t", index=False)
    table.to_csv(csv_out, index=False)

    try:
        md = table.to_markdown(index=False)
    except Exception:
        md = table.to_string(index=False)
    md_out.write_text(md + "\n", encoding="utf-8")

    audit = f"""Supplementary Table 3 demographic extraction audit

ADNI
----
{adni_audit}

NHANES
------
{nhanes_audit}

MHAS
----
{mhas_audit}

Resolved input paths
--------------------
ADNI: {adni_path}
NHANES: {nhanes_path}
MHAS: {mhas_path}

Output files
------------
TSV: {tsv_out}
CSV: {csv_out}
Markdown: {md_out}

Notes
-----
- This script outputs only ADNI, NHANES, and MHAS rows.
- N is computed after deduplicating participants by the detected ID column.
- Age is reported as mean plus/minus SD.
- Sex is reported as female count/female percentage.
- ADNI default source is copied from the CUBIC AD EPOCH predictions file.
- NHANES default source is:
  /Users/hao/Dropbox/NHANES/output_dir/model2_nondisease_mortality_epoch/nhanes_model2_analysis_table.tsv.gz
- MHAS default source is:
  /Users/hao/Dropbox/MHAS/step1_2001_mortality_epoch/mhas_2001_step1_clean_analytic_cohort.tsv
"""
    audit_out.write_text(audit, encoding="utf-8")

    log("\n" + md + "\n")
    log(audit)
    log("Done.")


if __name__ == "__main__":
    main()
