from __future__ import annotations
import pandas as pd
#!/usr/bin/env python3
"""
Collect 47 significant disease-specific EPOCH clocks and 22 mortality EPOCH
clocks into one participant-level table.

Input structure
---------------
WholeBodyClock/
    Brain_proteomics_dementia_clock/
        brain_proteomics_dementia_clock_predictions.tsv
    brain_mri_mortality_clock/
        brain_mri_mortality_clock_predictions.tsv
    ...

Each prediction file must contain:
    participant_id
    *_clock_acceleration_z

Outputs
-------
1. significant_epoch_clocks_wide.tsv
   One row per participant and 69 EPOCH acceleration-z columns.

2. significant_epoch_clock_manifest.tsv
   Clock metadata, source paths, sample sizes and missingness.

3. significant_epoch_collection_qc.txt
   Collection summary and any failures.

Default behavior
----------------
An outer join is used so that participants are retained even when they do not
have measurements for every organ or omics modality.
"""

import argparse
import re
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd


DEFAULT_ROOT = Path(
    "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock"
)


# ============================================================================
# 47 significant disease-specific EPOCH clocks
# ============================================================================

DISEASE_CLOCKS = [
    "Brain_proteomics_dementia_clock",
    "Heart_proteomics_dementia_clock",
    "Immune_proteomics_stroke_clock",
    "Brain_proteomics_mi_clock",
    "Heart_proteomics_mi_clock",
    "Metabolic_metabolomics_asthma_clock",
    "Brain_proteomics_stroke_clock",
    "Hepatic_metabolomics_asthma_clock",
    "Metabolic_metabolomics_dementia_clock",
    "Digestive_metabolomics_asthma_clock",
    "Hepatic_metabolomics_dementia_clock",
    "Metabolic_metabolomics_stroke_clock",
    "Digestive_metabolomics_dementia_clock",
    "Hepatic_metabolomics_mi_clock",
    "Pulmonary_proteomics_mi_clock",
    "Digestive_metabolomics_mi_clock",
    "Hepatic_metabolomics_stroke_clock",
    "Pulmonary_proteomics_stroke_clock",
    "Digestive_metabolomics_stroke_clock",
    "Hepatic_proteomics_asthma_clock",
    "Reproductive_female_proteomics_asthma_clock",
    "Endocrine_metabolomics_asthma_clock",
    "Hepatic_proteomics_dementia_clock",
    "Reproductive_female_proteomics_copd_clock",
    "Endocrine_metabolomics_dementia_clock",
    "Hepatic_proteomics_mi_clock",
    "Reproductive_female_proteomics_dementia_clock",
    "Endocrine_metabolomics_mi_clock",
    "Hepatic_proteomics_stroke_clock",
    "Reproductive_female_proteomics_mi_clock",
    "Endocrine_metabolomics_stroke_clock",
    "Immune_metabolomics_asthma_clock",
    "Reproductive_female_proteomics_stroke_clock",
    "Endocrine_proteomics_asthma_clock",
    "Immune_metabolomics_mi_clock",
    "Reproductive_male_proteomics_dementia_clock",
    "Endocrine_proteomics_dementia_clock",
    "Immune_metabolomics_stroke_clock",
    "Reproductive_male_proteomics_mi_clock",
    "Endocrine_proteomics_mi_clock",
    "Immune_proteomics_asthma_clock",
    "Reproductive_male_proteomics_stroke_clock",
    "Endocrine_proteomics_stroke_clock",
    "Immune_proteomics_dementia_clock",
    "spleen_mri_asthma_clock",
    "heart_mri_copd_clock",
    "Immune_proteomics_mi_clock",
]


# ============================================================================
# 22 mortality EPOCH clocks
# ============================================================================

MORTALITY_CLOCKS = [
    # Seven MRI mortality clocks
    "brain_mri_mortality_clock",
    "adipose_mri_mortality_clock",
    "heart_mri_mortality_clock",
    "kidney_mri_mortality_clock",
    "liver_mri_mortality_clock",
    "pancreas_mri_mortality_clock",
    "spleen_mri_mortality_clock",

    # Eleven proteomics mortality clocks
    "Reproductive_female_proteomics_mortality_clock",
    "Pulmonary_proteomics_mortality_clock",
    "Heart_proteomics_mortality_clock",
    "Brain_proteomics_mortality_clock",
    "Eye_proteomics_mortality_clock",
    "Hepatic_proteomics_mortality_clock",
    "Renal_proteomics_mortality_clock",
    "Reproductive_male_proteomics_mortality_clock",
    "Endocrine_proteomics_mortality_clock",
    "Immune_proteomics_mortality_clock",
    "Skin_proteomics_mortality_clock",

    # Four metabolomics mortality clocks
    "Endocrine_metabolomics_mortality_clock",
    "Digestive_metabolomics_mortality_clock",
    "Hepatic_metabolomics_mortality_clock",
    "Immune_metabolomics_mortality_clock",
]


@dataclass
class ClockRecord:
    clock_name: str
    clock_class: str
    organ: str
    modality: str
    endpoint: str
    source_file: str
    source_acceleration_column: str
    output_column: str
    n_rows: int
    n_unique_participants: int
    n_nonmissing: int
    n_missing: int


def normalize_name(value: str) -> str:
    """Normalize names for case-insensitive file and directory matching."""
    return re.sub(
        r"[^a-z0-9]+",
        "_",
        value.lower(),
    ).strip("_")


def parse_clock_name(clock_name: str) -> tuple[str, str, str]:
    """
    Parse:
        Reproductive_female_proteomics_asthma_clock

    into:
        organ = reproductive_female
        modality = proteomics
        endpoint = asthma
    """

    stem = re.sub(
        r"_clock$",
        "",
        clock_name,
        flags=re.IGNORECASE,
    )

    parts = stem.split("_")

    modality_index = None

    for index, token in enumerate(parts):
        if token.lower() in {
            "mri",
            "proteomics",
            "metabolomics",
        }:
            modality_index = index
            break

    if modality_index is None:
        raise ValueError(
            f"Could not identify modality in clock name: {clock_name}"
        )

    if modality_index == 0 or modality_index == len(parts) - 1:
        raise ValueError(
            f"Could not parse organ/modality/endpoint: {clock_name}"
        )

    organ = "_".join(parts[:modality_index]).lower()
    modality = parts[modality_index].lower()
    endpoint = "_".join(parts[modality_index + 1:]).lower()

    return organ, modality, endpoint


def resolve_clock_directory(
    root: Path,
    requested_name: str,
) -> Path:
    """
    Resolve the clock directory.

    First tries the supplied capitalization exactly. If that fails, searches
    for a case-insensitive normalized match.
    """

    exact_path = root / requested_name

    if exact_path.is_dir():
        return exact_path

    requested_normalized = normalize_name(requested_name)

    matches = [
        path
        for path in root.iterdir()
        if path.is_dir()
        and normalize_name(path.name) == requested_normalized
    ]

    if len(matches) == 1:
        return matches[0]

    if len(matches) == 0:
        raise FileNotFoundError(
            f"Clock directory not found: {requested_name}"
        )

    raise RuntimeError(
        f"Multiple directories matched {requested_name}: "
        + ", ".join(str(path) for path in matches)
    )


def resolve_prediction_file(
    clock_directory: Path,
    clock_name: str,
) -> Path:
    """
    Find the *_clock_predictions.tsv file inside a clock directory.
    """

    candidates = sorted(
        clock_directory.glob("*_clock_predictions.tsv")
    )

    if len(candidates) == 1:
        return candidates[0]

    if len(candidates) == 0:
        raise FileNotFoundError(
            f"No *_clock_predictions.tsv file found in "
            f"{clock_directory}"
        )

    expected_stem = re.sub(
        r"_clock$",
        "",
        clock_name,
        flags=re.IGNORECASE,
    )

    expected_normalized = normalize_name(expected_stem)

    matched = []

    for candidate in candidates:
        candidate_stem = re.sub(
            r"_clock_predictions\.tsv$",
            "",
            candidate.name,
            flags=re.IGNORECASE,
        )

        if normalize_name(candidate_stem) == expected_normalized:
            matched.append(candidate)

    if len(matched) == 1:
        return matched[0]

    raise RuntimeError(
        f"Could not uniquely resolve prediction file in "
        f"{clock_directory}. Candidates: "
        + ", ".join(path.name for path in candidates)
    )


def identify_acceleration_column(
    columns: Iterable[str],
    clock_name: str,
) -> str:
    """
    Identify the clock-specific acceleration-z column.
    """

    columns = list(columns)

    candidates = [
        column
        for column in columns
        if column.lower().endswith(
            "_clock_acceleration_z"
        )
    ]

    if len(candidates) == 1:
        return candidates[0]

    clock_stem = re.sub(
        r"_clock$",
        "",
        clock_name,
        flags=re.IGNORECASE,
    )

    expected_column = (
        f"{clock_stem.lower()}_clock_acceleration_z"
    )

    exact_case_insensitive = [
        column
        for column in columns
        if column.lower() == expected_column
    ]

    if len(exact_case_insensitive) == 1:
        return exact_case_insensitive[0]

    normalized_expected = normalize_name(expected_column)

    normalized_matches = [
        column
        for column in candidates
        if normalize_name(column) == normalized_expected
    ]

    if len(normalized_matches) == 1:
        return normalized_matches[0]

    raise ValueError(
        f"Expected one acceleration-z column for "
        f"{clock_name}, but found: {candidates}"
    )


def collect_one_clock(
    root: Path,
    clock_name: str,
    clock_class: str,
    duplicate_policy: str,
) -> tuple[pd.Series, ClockRecord]:
    """
    Read one clock file and return its participant-indexed acceleration-z
    values and QC metadata.
    """

    clock_directory = resolve_clock_directory(
        root=root,
        requested_name=clock_name,
    )

    prediction_file = resolve_prediction_file(
        clock_directory=clock_directory,
        clock_name=clock_name,
    )

    header = pd.read_csv(
        prediction_file,
        sep="\t",
        nrows=0,
    )

    if "participant_id" not in header.columns:
        raise ValueError(
            f"participant_id is missing from {prediction_file}"
        )

    acceleration_column = identify_acceleration_column(
        columns=header.columns,
        clock_name=clock_name,
    )

    data = pd.read_csv(
        prediction_file,
        sep="\t",
        usecols=[
            "participant_id",
            acceleration_column,
        ],
        low_memory=False,
    )

    data["participant_id"] = pd.to_numeric(
        data["participant_id"],
        errors="raise",
    ).astype("Int64")

    data[acceleration_column] = pd.to_numeric(
        data[acceleration_column],
        errors="coerce",
    )

    duplicate_mask = data["participant_id"].duplicated(
        keep=False
    )

    if duplicate_mask.any():
        number_duplicate_rows = int(
            duplicate_mask.sum()
        )

        duplicate_examples = (
            data.loc[
                duplicate_mask,
                "participant_id",
            ]
            .head(10)
            .tolist()
        )

        if duplicate_policy == "error":
            raise ValueError(
                f"{prediction_file} contains "
                f"{number_duplicate_rows} rows associated with "
                f"duplicated participant IDs. "
                f"Examples: {duplicate_examples}"
            )

        if duplicate_policy == "first":
            data = data.drop_duplicates(
                subset="participant_id",
                keep="first",
            )

        elif duplicate_policy == "mean":
            data = (
                data.groupby(
                    "participant_id",
                    as_index=False,
                    dropna=False,
                )[acceleration_column]
                .mean()
            )

    organ, modality, endpoint = parse_clock_name(
        clock_name
    )

    output_column = (
        f"epoch__{organ}__{modality}__"
        f"{endpoint}__acceleration_z"
    )

    values = (
        data.set_index("participant_id")[
            acceleration_column
        ]
        .rename(output_column)
    )

    record = ClockRecord(
        clock_name=clock_name,
        clock_class=clock_class,
        organ=organ,
        modality=modality,
        endpoint=endpoint,
        source_file=str(prediction_file),
        source_acceleration_column=acceleration_column,
        output_column=output_column,
        n_rows=int(len(data)),
        n_unique_participants=int(
            data["participant_id"].nunique(
                dropna=True
            )
        ),
        n_nonmissing=int(
            data[acceleration_column].notna().sum()
        ),
        n_missing=int(
            data[acceleration_column].isna().sum()
        ),
    )

    return values, record


def validate_clock_lists() -> None:
    """Ensure that the hard-coded lists have the intended sizes."""

    if len(DISEASE_CLOCKS) != 47:
        raise RuntimeError(
            f"Expected 47 disease clocks, "
            f"found {len(DISEASE_CLOCKS)}"
        )

    if len(MORTALITY_CLOCKS) != 22:
        raise RuntimeError(
            f"Expected 22 mortality clocks, "
            f"found {len(MORTALITY_CLOCKS)}"
        )

    all_clock_names = (
        DISEASE_CLOCKS + MORTALITY_CLOCKS
    )

    normalized_names = [
        normalize_name(name)
        for name in all_clock_names
    ]

    duplicated_names = sorted(
        {
            name
            for name in normalized_names
            if normalized_names.count(name) > 1
        }
    )

    if duplicated_names:
        raise RuntimeError(
            "Duplicated clock names detected: "
            + ", ".join(duplicated_names)
        )


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Collect 47 significant disease EPOCH clocks "
            "and 22 mortality EPOCH clocks."
        )
    )

    parser.add_argument(
        "--root",
        type=Path,
        default=DEFAULT_ROOT,
        help=(
            "WholeBodyClock root directory. "
            f"Default: {DEFAULT_ROOT}"
        ),
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=(
            "Output directory. Default: "
            "<root>/collected_significant_epoch_clocks"
        ),
    )

    parser.add_argument(
        "--join",
        choices=[
            "outer",
            "inner",
        ],
        default="outer",
        help=(
            "Join strategy across clock files. "
            "Default: outer."
        ),
    )

    parser.add_argument(
        "--duplicate-policy",
        choices=[
            "error",
            "first",
            "mean",
        ],
        default="error",
        help=(
            "Handling of duplicate participant IDs "
            "within a prediction file."
        ),
    )

    parser.add_argument(
        "--allow-missing-files",
        action="store_true",
        help=(
            "Skip clocks that cannot be collected "
            "instead of stopping immediately."
        ),
    )

    return parser


def main() -> int:
    arguments = build_argument_parser().parse_args()

    validate_clock_lists()

    root = arguments.root.expanduser().resolve()

    if not root.is_dir():
        raise NotADirectoryError(
            f"WholeBodyClock directory does not exist: "
            f"{root}"
        )

    if arguments.output_dir is None:
        output_directory = (
            root
            / "collected_significant_epoch_clocks"
        )
    else:
        output_directory = (
            arguments.output_dir
            .expanduser()
            .resolve()
        )

    output_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    requested_clocks = [
        *[
            (clock_name, "disease")
            for clock_name in DISEASE_CLOCKS
        ],
        *[
            (clock_name, "mortality")
            for clock_name in MORTALITY_CLOCKS
        ],
    ]

    collected_series: list[pd.Series] = []
    records: list[ClockRecord] = []
    failures: list[str] = []

    print(f"WholeBodyClock root: {root}")

    print(
        f"Collecting {len(DISEASE_CLOCKS)} disease "
        f"and {len(MORTALITY_CLOCKS)} mortality clocks"
    )

    for index, (
        clock_name,
        clock_class,
    ) in enumerate(
        requested_clocks,
        start=1,
    ):
        print(
            f"[{index:02d}/{len(requested_clocks)}] "
            f"{clock_name}"
        )

        try:
            values, record = collect_one_clock(
                root=root,
                clock_name=clock_name,
                clock_class=clock_class,
                duplicate_policy=(
                    arguments.duplicate_policy
                ),
            )

            collected_series.append(values)
            records.append(record)

        except Exception as error:
            failure_message = (
                f"{clock_name}\t"
                f"{type(error).__name__}: {error}"
            )

            failures.append(failure_message)

            print(
                f"  ERROR: {failure_message}",
                file=sys.stderr,
            )

            if not arguments.allow_missing_files:
                raise

    if not collected_series:
        raise RuntimeError(
            "No clock files were successfully collected."
        )

    combined = pd.concat(
        collected_series,
        axis=1,
        join=arguments.join,
        sort=True,
    )

    combined.index.name = "participant_id"

    combined = (
        combined
        .reset_index()
        .sort_values(
            by="participant_id",
            kind="stable",
        )
    )

    epoch_columns = [
        column
        for column in combined.columns
        if column.startswith("epoch__")
    ]

    if not failures and len(epoch_columns) != 69:
        raise RuntimeError(
            f"Expected 69 EPOCH columns, "
            f"but collected {len(epoch_columns)}"
        )

    wide_output = (
        output_directory
        / "significant_epoch_clocks_wide.tsv"
    )

    manifest_output = (
        output_directory
        / "significant_epoch_clock_manifest.tsv"
    )

    qc_output = (
        output_directory
        / "significant_epoch_collection_qc.txt"
    )

    combined.to_csv(
        wide_output,
        sep="\t",
        index=False,
        compression="gzip",
    )

    manifest = pd.DataFrame(
        [
            asdict(record)
            for record in records
        ]
    )

    manifest.to_csv(
        manifest_output,
        sep="\t",
        index=False,
    )

    number_disease_collected = int(
        (
            manifest["clock_class"] == "disease"
        ).sum()
    )

    number_mortality_collected = int(
        (
            manifest["clock_class"] == "mortality"
        ).sum()
    )

    qc_lines = [
        "Significant EPOCH clock collection QC",
        "====================================",
        f"Root directory: {root}",
        f"Join strategy: {arguments.join}",
        (
            "Duplicate participant policy: "
            f"{arguments.duplicate_policy}"
        ),
        "",
        (
            "Disease clocks requested: "
            f"{len(DISEASE_CLOCKS)}"
        ),
        (
            "Disease clocks collected: "
            f"{number_disease_collected}"
        ),
        (
            "Mortality clocks requested: "
            f"{len(MORTALITY_CLOCKS)}"
        ),
        (
            "Mortality clocks collected: "
            f"{number_mortality_collected}"
        ),
        (
            "Total clocks collected: "
            f"{len(records)}"
        ),
        (
            "Participants in wide output: "
            f"{len(combined)}"
        ),
        (
            "EPOCH columns in wide output: "
            f"{len(epoch_columns)}"
        ),
        "",
        f"Wide output: {wide_output}",
        f"Manifest: {manifest_output}",
        "",
        "Failures",
        "--------",
    ]

    if failures:
        qc_lines.extend(failures)
    else:
        qc_lines.append("None")

    qc_output.write_text(
        "\n".join(qc_lines) + "\n",
        encoding="utf-8",
    )

    print("\nCollection completed")

    print(
        f"  Disease clocks: "
        f"{number_disease_collected}"
    )

    print(
        f"  Mortality clocks: "
        f"{number_mortality_collected}"
    )

    print(
        f"  Total EPOCH columns: "
        f"{len(epoch_columns)}"
    )

    print(
        f"  Participants: "
        f"{len(combined):,}"
    )

    print(f"  Wide table: {wide_output}")
    print(f"  Manifest: {manifest_output}")
    print(f"  QC report: {qc_output}")

    if failures:
        print(
            f"Warning: {len(failures)} clocks failed. "
            f"See {qc_output}",
            file=sys.stderr,
        )

        return 2

    return 0


if __name__ == "__main__":
    raise SystemExit(main())