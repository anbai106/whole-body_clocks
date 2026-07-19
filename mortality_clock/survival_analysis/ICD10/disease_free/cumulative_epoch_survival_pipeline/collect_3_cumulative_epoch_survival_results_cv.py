#!/usr/bin/env python3
"""
Collect cross-validated cumulative EPOCH survival outputs across ICD endpoints.

Outputs
-------
1. *_all_steps.tsv
   All cumulative model steps for all collected diseases.

2. *_final_step.tsv
   One row per disease corresponding to the final cumulative model.
   By default, the expected final model contains all 15 clocks.

3. *_final_step_ranked_by_delta_gain.tsv
   Final models ranked from largest to smallest cross-validated C-index gain.

4. *_top2_delta_gain.tsv
   The two diseases with the largest cross-validated final-model gain.

5. *_best_cv_cindex_step.tsv
   The cumulative step with the largest cross-validated C-index for each disease.

6. *_rank_order.tsv
   Combined clock-ranking files, when available.

7. *_collection_summary.tsv
   Collection and quality-control summary.

Definitions
-----------
delta_gain:
    final model cross-validated C-index minus baseline cross-validated C-index.

delta_gain_apparent:
    final model apparent/in-sample C-index minus baseline apparent C-index.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


REQUIRED_RESULT_COLUMNS = {
    "disease_id",
    "cumulative_step",
    "n_clocks",
    "status",
    "c_index",
    "base_c_index",
    "cv_c_index",
    "cv_base_c_index",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Collect cross-validated cumulative mortality EPOCH survival "
            "outputs across ICD endpoints."
        )
    )

    parser.add_argument(
        "--disease_file",
        default=(
            "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/"
            "mortality_clock/SA/data/included_ICD_mortality_clock.tsv"
        ),
        help=(
            "TSV containing disease IDs. The first column is used unless "
            "--disease_column is supplied."
        ),
    )

    parser.add_argument(
        "--disease_column",
        default=None,
        help=(
            "Optional disease-ID column name in --disease_file. "
            "Default: use the first column."
        ),
    )

    parser.add_argument(
        "--input_dir",
        default=(
            "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/"
            "mortality_clock/SA/output_cumulative_EPOCH_PM/disease_free_cv"
        ),
        help="Directory containing per-disease cumulative result TSV files.",
    )

    parser.add_argument(
        "--output_prefix",
        default=(
            "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/"
            "mortality_clock/SA/output_cumulative_EPOCH_PM/"
            "combined_cumulative_EPOCH_PM_cv"
        ),
        help="Prefix for combined output files.",
    )

    parser.add_argument(
        "--min_case",
        type=int,
        default=20,
        help="Minimum number of incident cases required. Default: 20.",
    )

    parser.add_argument(
        "--expected_final_n_clocks",
        type=int,
        default=15,
        help=(
            "Expected number of clocks in the complete final model. "
            "Default: 15."
        ),
    )

    parser.add_argument(
        "--top_n",
        type=int,
        default=2,
        help=(
            "Number of diseases with the largest final cross-validated gain "
            "to save separately. Default: 2."
        ),
    )

    parser.add_argument(
        "--allow_incomplete_final",
        action="store_true",
        help=(
            "Allow the largest available cumulative step to be treated as "
            "the final model when the expected number of clocks is absent. "
            "By default, diseases without all expected clocks are excluded "
            "from final-model outputs."
        ),
    )

    return parser.parse_args()


def load_diseases(path: str, disease_column: str | None = None) -> list[str]:
    disease_path = Path(path)

    if not disease_path.is_file():
        raise FileNotFoundError(f"Disease file does not exist: {disease_path}")

    df = pd.read_csv(disease_path, sep="\t", dtype=str)

    if df.empty or df.shape[1] == 0:
        raise ValueError(f"Disease file is empty: {disease_path}")

    if disease_column is None:
        series = df.iloc[:, 0]
    else:
        if disease_column not in df.columns:
            raise ValueError(
                f"Column '{disease_column}' was not found in {disease_path}. "
                f"Available columns: {list(df.columns)}"
            )
        series = df[disease_column]

    diseases = (
        series.dropna()
        .astype(str)
        .str.strip()
        .loc[lambda x: x.ne("")]
        .drop_duplicates()
        .tolist()
    )

    if not diseases:
        raise ValueError(f"No disease IDs were found in {disease_path}")

    return diseases


def safe_numeric(series: pd.Series) -> pd.Series:
    """Convert a pandas Series to numeric, coercing invalid values to NaN."""
    return pd.to_numeric(series, errors="coerce")


def validate_result_columns(df: pd.DataFrame, path: Path) -> None:
    missing = sorted(REQUIRED_RESULT_COLUMNS.difference(df.columns))

    if missing:
        raise ValueError(
            f"{path} is missing required columns: {', '.join(missing)}"
        )


def read_result_file(path: Path, expected_disease: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    validate_result_columns(df, path)

    if df.empty:
        raise ValueError(f"Result file is empty: {path}")

    # Normalize numeric columns used by the collector.
    numeric_columns = [
        "N",
        "N_case",
        "N_noncase",
        "event_rate",
        "cumulative_step",
        "n_clocks",
        "c_index",
        "base_c_index",
        "delta_c_index_vs_base",
        "delta_c_index_vs_previous",
        "cv_folds",
        "cv_c_index",
        "cv_base_c_index",
        "delta_cv_c_index_vs_base",
        "delta_cv_c_index_vs_previous",
    ]

    for column in numeric_columns:
        if column in df.columns:
            df[column] = safe_numeric(df[column])

    # Protect against mismatches between filename and embedded disease ID.
    observed_ids = (
        df["disease_id"]
        .dropna()
        .astype(str)
        .str.strip()
        .unique()
        .tolist()
    )

    if observed_ids and set(observed_ids) != {expected_disease}:
        raise ValueError(
            f"Disease-ID mismatch in {path}. "
            f"Expected '{expected_disease}', observed {observed_ids}."
        )

    # Ensure disease_id is present even if malformed/missing in some rows.
    df["disease_id"] = expected_disease
    df["source_result_file"] = str(path)

    return df


def select_baseline_row(disease_df: pd.DataFrame) -> pd.Series:
    """
    Select the baseline model row.

    Primary definition:
      cumulative_step == 0

    Secondary checks:
      added_pair_id == BASE or added_clock == BASE
    """
    baseline = disease_df.loc[disease_df["cumulative_step"].eq(0)].copy()

    if baseline.empty and "added_pair_id" in disease_df.columns:
        baseline = disease_df.loc[
            disease_df["added_pair_id"].astype(str).str.upper().eq("BASE")
        ].copy()

    if baseline.empty and "added_clock" in disease_df.columns:
        baseline = disease_df.loc[
            disease_df["added_clock"].astype(str).str.upper().eq("BASE")
        ].copy()

    if baseline.empty:
        raise ValueError(
            f"No baseline row found for disease "
            f"{disease_df['disease_id'].iloc[0]}."
        )

    # If duplicates exist, use the first after sorting by cumulative step.
    baseline = baseline.sort_values("cumulative_step", kind="stable")
    return baseline.iloc[0]


def select_final_row(
    disease_df: pd.DataFrame,
    expected_final_n_clocks: int,
    allow_incomplete_final: bool,
) -> tuple[pd.Series | None, bool]:
    """
    Select the final cumulative model.

    Returns
    -------
    final_row
        Selected final row, or None if no complete final model exists and
        incomplete models are not allowed.

    complete_final_model
        True when n_clocks equals expected_final_n_clocks.
    """
    usable = disease_df.copy()

    # Prefer rows whose model fitting succeeded.
    if "status" in usable.columns:
        ok_rows = usable.loc[
            usable["status"].astype(str).str.lower().eq("ok")
        ].copy()
        if not ok_rows.empty:
            usable = ok_rows

    complete = usable.loc[
        usable["n_clocks"].eq(expected_final_n_clocks)
    ].copy()

    if not complete.empty:
        complete = complete.sort_values(
            ["cumulative_step", "n_clocks"],
            ascending=[False, False],
            kind="stable",
        )
        return complete.iloc[0], True

    if not allow_incomplete_final:
        return None, False

    usable = usable.sort_values(
        ["n_clocks", "cumulative_step"],
        ascending=[False, False],
        kind="stable",
    )

    if usable.empty:
        return None, False

    return usable.iloc[0], False


def build_final_summary(
    all_res: pd.DataFrame,
    expected_final_n_clocks: int,
    allow_incomplete_final: bool,
) -> tuple[pd.DataFrame, list[str], list[str]]:
    final_rows: list[dict] = []
    missing_baseline: list[str] = []
    incomplete_final: list[str] = []

    for disease_id, disease_df in all_res.groupby(
        "disease_id", sort=False, dropna=False
    ):
        disease_id = str(disease_id)

        try:
            baseline = select_baseline_row(disease_df)
        except ValueError:
            missing_baseline.append(disease_id)
            continue

        final, complete_final_model = select_final_row(
            disease_df=disease_df,
            expected_final_n_clocks=expected_final_n_clocks,
            allow_incomplete_final=allow_incomplete_final,
        )

        if final is None:
            incomplete_final.append(disease_id)
            continue

        if not complete_final_model:
            incomplete_final.append(disease_id)

        row = final.to_dict()

        baseline_c_index = pd.to_numeric(
            baseline.get("c_index", np.nan), errors="coerce"
        )
        baseline_cv_c_index = pd.to_numeric(
            baseline.get("cv_c_index", np.nan), errors="coerce"
        )

        final_c_index = pd.to_numeric(
            final.get("c_index", np.nan), errors="coerce"
        )
        final_cv_c_index = pd.to_numeric(
            final.get("cv_c_index", np.nan), errors="coerce"
        )

        row.update(
            {
                "baseline_c_index_from_step0": baseline_c_index,
                "baseline_cv_c_index_from_step0": baseline_cv_c_index,
                "final_c_index": final_c_index,
                "final_cv_c_index": final_cv_c_index,
                "delta_gain_apparent": (
                    final_c_index - baseline_c_index
                    if pd.notna(final_c_index)
                    and pd.notna(baseline_c_index)
                    else np.nan
                ),
                "delta_gain": (
                    final_cv_c_index - baseline_cv_c_index
                    if pd.notna(final_cv_c_index)
                    and pd.notna(baseline_cv_c_index)
                    else np.nan
                ),
                "complete_final_model": complete_final_model,
                "expected_final_n_clocks": expected_final_n_clocks,
            }
        )

        # Internal consistency checks against values already emitted by
        # the endpoint-level analysis script.
        reported_apparent_gain = pd.to_numeric(
            final.get("delta_c_index_vs_base", np.nan),
            errors="coerce",
        )
        reported_cv_gain = pd.to_numeric(
            final.get("delta_cv_c_index_vs_base", np.nan),
            errors="coerce",
        )

        row["reported_delta_c_index_vs_base"] = reported_apparent_gain
        row["reported_delta_cv_c_index_vs_base"] = reported_cv_gain

        row["delta_gain_apparent_difference_from_reported"] = (
            row["delta_gain_apparent"] - reported_apparent_gain
            if pd.notna(row["delta_gain_apparent"])
            and pd.notna(reported_apparent_gain)
            else np.nan
        )

        row["delta_gain_difference_from_reported"] = (
            row["delta_gain"] - reported_cv_gain
            if pd.notna(row["delta_gain"])
            and pd.notna(reported_cv_gain)
            else np.nan
        )

        final_rows.append(row)

    final_df = pd.DataFrame(final_rows)

    if not final_df.empty:
        final_df = final_df.sort_values(
            ["delta_gain", "disease_id"],
            ascending=[False, True],
            na_position="last",
            kind="stable",
        ).reset_index(drop=True)

        final_df.insert(
            0,
            "delta_gain_rank",
            np.arange(1, len(final_df) + 1),
        )

    return final_df, missing_baseline, incomplete_final


def select_best_cv_step(all_res: pd.DataFrame) -> pd.DataFrame:
    """
    Select the step with the highest valid cross-validated C-index for
    every disease.

    Ties are resolved in favor of the smaller cumulative step, yielding
    the more parsimonious model.
    """
    best = all_res.copy()

    if "status" in best.columns:
        best = best.loc[
            best["status"].astype(str).str.lower().eq("ok")
        ].copy()

    if "cv_status" in best.columns:
        best = best.loc[
            best["cv_status"].astype(str).str.lower().eq("ok")
        ].copy()

    best = best.loc[best["cv_c_index"].notna()].copy()

    if best.empty:
        return best

    best = (
        best.sort_values(
            ["disease_id", "cv_c_index", "cumulative_step"],
            ascending=[True, False, True],
            kind="stable",
        )
        .groupby("disease_id", as_index=False, sort=False)
        .head(1)
        .reset_index(drop=True)
    )

    return best


def main() -> int:
    args = parse_args()

    input_dir = Path(args.input_dir)
    out_prefix = Path(args.output_prefix)

    if not input_dir.is_dir():
        raise NotADirectoryError(
            f"Input directory does not exist: {input_dir}"
        )

    if args.expected_final_n_clocks < 1:
        raise ValueError("--expected_final_n_clocks must be at least 1.")

    if args.top_n < 1:
        raise ValueError("--top_n must be at least 1.")

    diseases = load_diseases(
        path=args.disease_file,
        disease_column=args.disease_column,
    )

    frames: list[pd.DataFrame] = []
    rank_frames: list[pd.DataFrame] = []

    missing_outputs: list[str] = []
    unreadable_outputs: list[str] = []
    below_min_case: list[str] = []

    for disease in diseases:
        result_path = (
            input_dir / f"cox_cumulative_EPOCH_PM_{disease}.tsv"
        )
        rank_path = input_dir / f"rank_order_EPOCH_PM_{disease}.tsv"

        if not result_path.is_file():
            missing_outputs.append(disease)
            continue

        try:
            result_df = read_result_file(
                path=result_path,
                expected_disease=disease,
            )
        except Exception as exc:
            print(
                f"WARNING: Could not read {result_path}: {exc}",
                file=sys.stderr,
            )
            unreadable_outputs.append(disease)
            continue

        n_case_values = result_df["N_case"].dropna()

        if not n_case_values.empty:
            n_case = int(n_case_values.iloc[0])
            if n_case < args.min_case:
                below_min_case.append(disease)
                continue

        frames.append(result_df)

        if rank_path.is_file():
            try:
                rank_df = pd.read_csv(
                    rank_path,
                    sep="\t",
                    low_memory=False,
                )

                if "disease_id" not in rank_df.columns:
                    rank_df.insert(0, "disease_id", disease)

                rank_df["source_rank_file"] = str(rank_path)
                rank_frames.append(rank_df)

            except Exception as exc:
                print(
                    f"WARNING: Could not read {rank_path}: {exc}",
                    file=sys.stderr,
                )

    if not frames:
        raise RuntimeError(
            f"No cumulative result files were collected from {input_dir}"
        )

    all_res = pd.concat(
        frames,
        ignore_index=True,
        sort=False,
    )

    all_res = all_res.sort_values(
        ["disease_id", "cumulative_step"],
        ascending=[True, True],
        kind="stable",
    ).reset_index(drop=True)

    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    all_steps_path = Path(str(out_prefix) + "_all_steps.tsv")
    all_res.to_csv(
        all_steps_path,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    final_df, missing_baseline, incomplete_final = build_final_summary(
        all_res=all_res,
        expected_final_n_clocks=args.expected_final_n_clocks,
        allow_incomplete_final=args.allow_incomplete_final,
    )

    if final_df.empty:
        raise RuntimeError(
            "No valid final models were identified. Check "
            "--expected_final_n_clocks or use --allow_incomplete_final."
        )

    final_path = Path(str(out_prefix) + "_final_step.tsv")
    final_df.to_csv(
        final_path,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    ranked_path = Path(
        str(out_prefix) + "_final_step_ranked_by_delta_gain.tsv"
    )
    final_df.to_csv(
        ranked_path,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    top_n_df = final_df.loc[
        final_df["delta_gain"].notna()
    ].head(args.top_n).copy()

    top_n_path = Path(
        str(out_prefix) + f"_top{args.top_n}_delta_gain.tsv"
    )
    top_n_df.to_csv(
        top_n_path,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    best_cv = select_best_cv_step(all_res)
    best_cv_path = Path(
        str(out_prefix) + "_best_cv_cindex_step.tsv"
    )
    best_cv.to_csv(
        best_cv_path,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    if rank_frames:
        rank_all = pd.concat(
            rank_frames,
            ignore_index=True,
            sort=False,
        )

        rank_path = Path(str(out_prefix) + "_rank_order.tsv")
        rank_all.to_csv(
            rank_path,
            sep="\t",
            index=False,
            na_rep="NA",
        )
    else:
        rank_path = None

    summary = pd.DataFrame(
        {
            "n_diseases_in_list": [len(diseases)],
            "n_diseases_with_result_file": [
                len(diseases) - len(missing_outputs)
            ],
            "n_diseases_collected_after_qc": [
                all_res["disease_id"].nunique()
            ],
            "n_diseases_in_final_output": [
                final_df["disease_id"].nunique()
            ],
            "n_missing_outputs": [len(missing_outputs)],
            "n_unreadable_outputs": [len(unreadable_outputs)],
            "n_below_min_case": [len(below_min_case)],
            "n_missing_baseline": [len(missing_baseline)],
            "n_incomplete_final_model": [len(set(incomplete_final))],
            "min_case": [args.min_case],
            "expected_final_n_clocks": [
                args.expected_final_n_clocks
            ],
            "allow_incomplete_final": [
                args.allow_incomplete_final
            ],
            "missing_outputs_first_20": [
                ",".join(missing_outputs[:20])
            ],
            "unreadable_outputs_first_20": [
                ",".join(unreadable_outputs[:20])
            ],
            "below_min_case_first_20": [
                ",".join(below_min_case[:20])
            ],
            "missing_baseline_first_20": [
                ",".join(missing_baseline[:20])
            ],
            "incomplete_final_first_20": [
                ",".join(sorted(set(incomplete_final))[:20])
            ],
        }
    )

    summary_path = Path(
        str(out_prefix) + "_collection_summary.tsv"
    )
    summary.to_csv(
        summary_path,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    print(
        f"Collected {all_res['disease_id'].nunique()} disease endpoints."
    )
    print(
        f"Final models available for "
        f"{final_df['disease_id'].nunique()} endpoints."
    )
    print(f"Wrote: {all_steps_path}")
    print(f"Wrote: {final_path}")
    print(f"Wrote: {ranked_path}")
    print(f"Wrote: {top_n_path}")
    print(f"Wrote: {best_cv_path}")

    if rank_path is not None:
        print(f"Wrote: {rank_path}")

    print(f"Wrote: {summary_path}")

    print("\nTop diseases ranked by final cross-validated C-index gain:")
    display_columns = [
        "delta_gain_rank",
        "disease_id",
        "N",
        "N_case",
        "baseline_cv_c_index_from_step0",
        "final_cv_c_index",
        "delta_gain",
        "n_clocks",
        "complete_final_model",
    ]

    existing_columns = [
        column
        for column in display_columns
        if column in top_n_df.columns
    ]

    print(
        top_n_df[existing_columns].to_string(index=False)
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())