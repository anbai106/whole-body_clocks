#!/usr/bin/env python3
"""Collect STEP 3 cumulative EPOCH survival outputs across ICD endpoints."""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Collect cumulative mortality EPOCH survival outputs.")
    p.add_argument("--disease_file", default="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/SA/data/included_ICD_mortality_clock.tsv")
    p.add_argument("--input_dir", default="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/SA/output_cumulative_EPOCH_PM/disease_free")
    p.add_argument("--output_prefix", default="/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/SA/output_cumulative_EPOCH_PM/combined_cumulative_EPOCH_PM")
    p.add_argument("--min_case", type=int, default=20)
    return p.parse_args()


def load_diseases(path: str) -> list[str]:
    df = pd.read_csv(path, sep="\t", dtype=str)
    return df.iloc[:, 0].dropna().astype(str).tolist()


def main() -> int:
    args = parse_args()
    input_dir = Path(args.input_dir)
    diseases = load_diseases(args.disease_file)
    frames = []
    rank_frames = []
    missing = []

    for disease in diseases:
        res_path = input_dir / f"cox_cumulative_EPOCH_PM_{disease}.tsv"
        rank_path = input_dir / f"rank_order_EPOCH_PM_{disease}.tsv"
        if not res_path.is_file():
            missing.append(disease)
            continue
        res = pd.read_csv(res_path, sep="\t")
        if "N_case" in res.columns and len(res) and int(res["N_case"].iloc[0]) < args.min_case:
            continue
        frames.append(res)
        if rank_path.is_file():
            rank_frames.append(pd.read_csv(rank_path, sep="\t"))

    if not frames:
        raise RuntimeError(f"No cumulative result files found in {input_dir}")

    all_res = pd.concat(frames, ignore_index=True)
    out_prefix = Path(args.output_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    all_res.to_csv(str(out_prefix) + "_all_steps.tsv", sep="\t", index=False, na_rep="NA")

    final = all_res.sort_values(["disease_id", "cumulative_step"]).groupby("disease_id", as_index=False).tail(1)
    final.to_csv(str(out_prefix) + "_final_step.tsv", sep="\t", index=False, na_rep="NA")

    best = all_res[all_res["status"] == "ok"].copy()
    best = best.sort_values(["disease_id", "c_index", "cumulative_step"], ascending=[True, False, True]).groupby("disease_id", as_index=False).head(1)
    best.to_csv(str(out_prefix) + "_best_cindex_step.tsv", sep="\t", index=False, na_rep="NA")

    if rank_frames:
        rank_all = pd.concat(rank_frames, ignore_index=True)
        rank_all.to_csv(str(out_prefix) + "_rank_order.tsv", sep="\t", index=False, na_rep="NA")

    summary = pd.DataFrame({
        "n_diseases_in_list": [len(diseases)],
        "n_diseases_collected": [all_res["disease_id"].nunique()],
        "n_missing_outputs": [len(missing)],
        "missing_outputs_first_20": [",".join(missing[:20])],
    })
    summary.to_csv(str(out_prefix) + "_collection_summary.tsv", sep="\t", index=False)

    print(f"Collected {all_res['disease_id'].nunique()} disease endpoints.")
    print(f"Wrote {out_prefix}_all_steps.tsv")
    print(f"Wrote {out_prefix}_final_step.tsv")
    print(f"Wrote {out_prefix}_best_cindex_step.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
