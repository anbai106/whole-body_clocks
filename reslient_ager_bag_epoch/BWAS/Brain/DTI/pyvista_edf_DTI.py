#!/usr/bin/env python3
"""
Visualize within-modality BH-FDR-significant DTI FA associations for
brain-proteomics EPOCH-BAG candidate resilient agers (CRA) versus
concordant unfavorable agers (CUA).

SIGNIFICANCE
------------
Only tracts satisfying

    P_value_FDR_BH_within_modality < 0.05

are rendered.

The exact same significant-tract mesh is used in every view.

COLOR
-----
The visualization uses the signed logistic-regression coefficient

    Beta_log_odds_per_1SD_IDP

with outcome coding

    CRA = 1
    CUA = 0

Therefore

    beta > 0  -> higher FA is associated with greater odds of CRA
    beta < 0  -> higher FA is associated with greater odds of CUA

A symmetric zero-centered diverging color scale is used.

OBLIQUE FULL-BRAIN VIEWS
------------------------
The pure lateral views are replaced by elevated 45-degree superior-lateral
oblique views, which generally expose deep white-matter structures more clearly
while retaining whole-brain context.

Default views:
    1. Superior
    2. Left superior-lateral oblique (45 degrees)
    3. Right superior-lateral oblique (45 degrees)

The oblique angle can be changed without editing the script:

    --oblique-angle 35
    --oblique-angle 45
    --oblique-angle 55

Angle definition:
    0 degrees  = pure lateral
    45 degrees = equal lateral and superior components
    90 degrees = pure superior

All views use:
    - parallel / orthographic projection
    - automatic projected-bound framing
    - the union of both cortical hemisphere bounds
    - the same framing margin

This prevents cropping and keeps the whole brain visible in every panel.

JHU / UK BIOBANK MAPPING
------------------------
The 48 UK Biobank FA fields f25056-f25103 map in order to JHU tract numbers
1-48:

    f25056 -> tract 1
    ...
    f25103 -> tract 48

DEFAULT INPUT
-------------
/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/mortality_clock/BWAS/
Brain_proteomics_EPOCH_BAG_resilience/CRA_vs_CUA/combined_results/
BWAS_CRA_vs_CUA_DTI_all_statistics.tsv
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyvista as pv
from matplotlib.image import imread


# =============================================================================
# DEFAULT PATHS
# =============================================================================

DEFAULT_INPUT_TSV = Path(
    "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/"
    "mortality_clock/BWAS/"
    "Brain_proteomics_EPOCH_BAG_resilience/"
    "CRA_vs_CUA/combined_results/"
    "BWAS_CRA_vs_CUA_DTI_all_statistics.tsv"
)

DEFAULT_DICTIONARY_TSV = Path(
    "/Users/hao/cubic-home/Dataset/Atlas/WM_JHU/"
    "dictionary_JHU_labels.tsv"
)

DEFAULT_TRACT_DIR = Path(
    "/Users/hao/cubic-home/Dataset/Atlas/WM_JHU/"
    "atlas-JHULabel"
)

DEFAULT_BG_LEFT = Path(
    "/Users/hao/cubic-home/Dataset/Atlas/WM_JHU/"
    "FsAverage_hemi-left_pial.vtk"
)

DEFAULT_BG_RIGHT = Path(
    "/Users/hao/cubic-home/Dataset/Atlas/WM_JHU/"
    "FsAverage_hemi-right_pial.vtk"
)

DEFAULT_OUTPUT_DIR = Path(
    "/Users/hao/Dropbox/2026_EPOCH/Fig/orig/BWAS/DTI"
)


# =============================================================================
# ANALYSIS SETTINGS
# =============================================================================

MODALITY = "DTI"
DEFAULT_FEATURE_FAMILY = "FA"

FDR_COL = "P_value_FDR_BH_within_modality"
BETA_COL = "Beta_log_odds_per_1SD_IDP"
RAW_P_COL = "P_value"
OR_COL = "OR_per_1SD_IDP"

SCALAR_NAME = "CRA_vs_CUA_beta"

DEFAULT_FDR_THRESHOLD = 0.05
DEFAULT_OBLIQUE_ANGLE = 45.0

FA_FIRST_FIELD = 25056
FA_LAST_FIELD = 25103
FA_TRACT_OFFSET = 25055


# =============================================================================
# ARGUMENTS
# =============================================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Visualize within-modality BH-FDR-significant DTI FA tracts "
            "using full-brain superior and 45-degree oblique views."
        )
    )

    parser.add_argument(
        "--input-tsv",
        type=Path,
        default=DEFAULT_INPUT_TSV,
    )

    parser.add_argument(
        "--dictionary-tsv",
        type=Path,
        default=DEFAULT_DICTIONARY_TSV,
    )

    parser.add_argument(
        "--tract-dir",
        type=Path,
        default=DEFAULT_TRACT_DIR,
    )

    parser.add_argument(
        "--bg-left",
        type=Path,
        default=DEFAULT_BG_LEFT,
    )

    parser.add_argument(
        "--bg-right",
        type=Path,
        default=DEFAULT_BG_RIGHT,
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
    )

    parser.add_argument(
        "--fdr-threshold",
        type=float,
        default=DEFAULT_FDR_THRESHOLD,
        help="Threshold for P_value_FDR_BH_within_modality.",
    )

    parser.add_argument(
        "--feature-family",
        default=DEFAULT_FEATURE_FAMILY,
    )

    parser.add_argument(
        "--oblique-angle",
        type=float,
        default=DEFAULT_OBLIQUE_ANGLE,
        help=(
            "Elevation above the lateral plane in degrees. "
            "0=pure lateral, 45=equal lateral/superior, 90=pure superior. "
            "Default: 45."
        ),
    )

    parser.add_argument(
        "--background-opacity",
        type=float,
        default=0.055,
        help=(
            "Cortical surface opacity. Lower than previous version so deep "
            "tracts remain visible. Default: 0.055."
        ),
    )

    parser.add_argument(
        "--tract-opacity",
        type=float,
        default=1.0,
    )

    parser.add_argument(
        "--frame-margin",
        type=float,
        default=1.10,
        help="Full-brain framing margin. Default: 1.10.",
    )

    parser.add_argument(
        "--window-width",
        type=int,
        default=1800,
    )

    parser.add_argument(
        "--window-height",
        type=int,
        default=1400,
    )

    parser.add_argument(
        "--title",
        default="CRA vs CUA: FDR-significant white-matter FA associations",
    )

    parser.add_argument(
        "--scalar-bar-title",
        default="Beta (CRA vs CUA)",
    )

    parser.add_argument(
        "--no-composite",
        action="store_true",
    )

    return parser.parse_args()


# =============================================================================
# VALIDATION
# =============================================================================

def require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"{label} not found:\n  {path}")


def require_dir(path: Path, label: str) -> None:
    if not path.is_dir():
        raise FileNotFoundError(f"{label} not found:\n  {path}")


def require_columns(
    df: pd.DataFrame,
    columns: Sequence[str],
    label: str,
) -> None:
    missing = [col for col in columns if col not in df.columns]
    if missing:
        raise ValueError(
            f"{label} is missing required column(s): "
            + ", ".join(missing)
        )


# =============================================================================
# IDP -> JHU TRACT NUMBER
# =============================================================================

def parse_ukb_field_number(idp: str) -> Optional[int]:
    match = re.search(r"_f(\d+)_\d+_\d+$", str(idp))
    return None if match is None else int(match.group(1))


def fa_idp_to_tract_number(idp: str) -> Optional[int]:
    field_number = parse_ukb_field_number(idp)

    if field_number is None:
        return None

    if not (FA_FIRST_FIELD <= field_number <= FA_LAST_FIELD):
        return None

    tract_number = field_number - FA_TRACT_OFFSET

    if not (1 <= tract_number <= 48):
        return None

    return tract_number


# =============================================================================
# READ / FILTER RESULTS
# =============================================================================

def prepare_results(
    input_tsv: Path,
    dictionary_tsv: Path,
    feature_family: str,
    fdr_threshold: float,
) -> Tuple[pd.DataFrame, pd.DataFrame]:

    stats = pd.read_csv(
        input_tsv,
        sep="\t",
        low_memory=False,
    )

    require_columns(
        stats,
        [
            "Modality",
            "Feature_family",
            "IDP",
            BETA_COL,
            FDR_COL,
            RAW_P_COL,
            OR_COL,
        ],
        "DTI statistics TSV",
    )

    stats = stats.loc[
        stats["Modality"].astype(str).str.upper().eq(MODALITY)
        & stats["Feature_family"].astype(str).str.upper().eq(feature_family.upper())
    ].copy()

    if stats.empty:
        raise ValueError(
            f"No DTI rows with Feature_family={feature_family} were found."
        )

    numeric_columns = [
        BETA_COL,
        FDR_COL,
        RAW_P_COL,
        OR_COL,
        "SE",
        "Z",
        "OR_CI95_lower",
        "OR_CI95_upper",
        "N_case",
        "N_control",
        "N_total",
    ]

    for col in numeric_columns:
        if col in stats.columns:
            stats[col] = pd.to_numeric(stats[col], errors="coerce")

    stats["UKB_field_number"] = stats["IDP"].map(parse_ukb_field_number)
    stats["Tract_number"] = stats["IDP"].map(fa_idp_to_tract_number)

    unmapped = stats.loc[
        stats["Tract_number"].isna(),
        "IDP",
    ]

    if len(unmapped) > 0:
        raise ValueError(
            "Some FA IDPs could not be mapped to JHU tract numbers 1-48:\n  "
            + "\n  ".join(unmapped.astype(str))
        )

    stats["Tract_number"] = stats["Tract_number"].astype(int)

    dictionary = pd.read_csv(
        dictionary_tsv,
        sep="\t",
        low_memory=False,
    )

    require_columns(
        dictionary,
        ["Tract_number", "Tract_abb"],
        "JHU dictionary",
    )

    dictionary["Tract_number"] = pd.to_numeric(
        dictionary["Tract_number"],
        errors="coerce",
    )

    dictionary = dictionary.dropna(
        subset=["Tract_number", "Tract_abb"]
    ).copy()

    dictionary["Tract_number"] = dictionary["Tract_number"].astype(int)

    if dictionary["Tract_number"].duplicated().any():
        raise ValueError(
            "JHU dictionary contains duplicated Tract_number values."
        )

    stats = stats.merge(
        dictionary[["Tract_number", "Tract_abb"]],
        on="Tract_number",
        how="left",
        validate="many_to_one",
    )

    if stats["Tract_abb"].isna().any():
        bad = stats.loc[
            stats["Tract_abb"].isna(),
            ["IDP", "Tract_number"],
        ]
        raise ValueError(
            "Some tract numbers are absent from the JHU dictionary:\n"
            + bad.to_string(index=False)
        )

    # ONLY inferential significance filter.
    significant = stats.loc[
        stats[FDR_COL].notna()
        & (stats[FDR_COL] < fdr_threshold)
    ].copy()

    significant = significant.sort_values(
        by=[FDR_COL, RAW_P_COL, "Tract_number"],
        ascending=[True, True, True],
        kind="mergesort",
    ).reset_index(drop=True)

    return stats, significant


# =============================================================================
# TRACT MESH
# =============================================================================

def tract_vtk_path(
    tract_dir: Path,
    tract_abb: str,
    tract_number: int,
) -> Path:
    return tract_dir / f"JHU_{tract_abb}_{tract_number}.vtk"


def build_significant_mesh(
    significant: pd.DataFrame,
    tract_dir: Path,
) -> Tuple[Optional[pv.DataSet], List[Path]]:

    if significant.empty:
        return None, []

    meshes: List[pv.DataSet] = []
    used_paths: List[Path] = []

    for _, row in significant.iterrows():
        tract_number = int(row["Tract_number"])
        tract_abb = str(row["Tract_abb"])
        beta = float(row[BETA_COL])
        raw_p = float(row[RAW_P_COL])
        fdr = float(row[FDR_COL])
        odds_ratio = float(row[OR_COL])

        vtk_file = tract_vtk_path(
            tract_dir,
            tract_abb,
            tract_number,
        )

        require_file(
            vtk_file,
            f"JHU VTK for significant tract {tract_abb}",
        )

        mesh = pv.read(str(vtk_file)).copy()

        if mesh.n_cells <= 0:
            raise ValueError(
                f"JHU tract mesh has zero cells: {vtk_file}"
            )

        mesh.cell_data[SCALAR_NAME] = np.full(
            mesh.n_cells,
            beta,
            dtype=float,
        )

        meshes.append(mesh)
        used_paths.append(vtk_file)

        print(
            f"{tract_abb:<12} "
            f"tract={tract_number:02d} "
            f"beta={beta:+.4f} "
            f"OR={odds_ratio:.3f} "
            f"P={raw_p:.3g} "
            f"FDR={fdr:.3g}"
        )

    mesh_final = meshes[0]

    for mesh in meshes[1:]:
        mesh_final = mesh_final.merge(
            mesh,
            merge_points=False,
        )

    return mesh_final, used_paths


# =============================================================================
# COLOR LIMITS
# =============================================================================

def symmetric_beta_limits(
    significant: pd.DataFrame,
) -> Tuple[float, float]:

    beta = pd.to_numeric(
        significant[BETA_COL],
        errors="coerce",
    ).to_numpy(dtype=float)

    beta = beta[np.isfinite(beta)]

    if beta.size == 0:
        return -1.0, 1.0

    max_abs = float(np.max(np.abs(beta)))

    if max_abs <= 0:
        max_abs = 1.0

    max_abs *= 1.05

    return -max_abs, max_abs


# =============================================================================
# OBLIQUE VIEW DEFINITIONS
# =============================================================================

def build_view_specs(
    oblique_angle_degrees: float,
) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    """
    Construct superior plus bilateral elevated lateral-oblique views.

    Coordinate convention follows the user's existing JHU rendering:
        lateral axis  = +/-Y
        superior axis = +Z

    angle:
        0 degrees  -> pure lateral
        45 degrees -> equal lateral and superior
        90 degrees -> superior
    """

    if not (0.0 < oblique_angle_degrees < 90.0):
        raise ValueError(
            "--oblique-angle must be > 0 and < 90 degrees."
        )

    theta = math.radians(
        oblique_angle_degrees
    )

    lateral_component = math.cos(theta)
    superior_component = math.sin(theta)

    # Full-brain elevated lateral views.
    left_oblique = np.array(
        [
            0.0,
            lateral_component,
            superior_component,
        ],
        dtype=float,
    )

    right_oblique = np.array(
        [
            0.0,
            -lateral_component,
            superior_component,
        ],
        dtype=float,
    )

    return {
        "superior": (
            np.array([0.0, 0.0, 1.0], dtype=float),
            np.array([0.0, 1.0, 0.0], dtype=float),
        ),
        f"left_oblique_{int(round(oblique_angle_degrees))}deg": (
            left_oblique,
            np.array([0.0, 0.0, 1.0], dtype=float),
        ),
        f"right_oblique_{int(round(oblique_angle_degrees))}deg": (
            right_oblique,
            np.array([0.0, 0.0, 1.0], dtype=float),
        ),
    }


# =============================================================================
# AUTOMATIC PROJECTED-BOUND FULL-BRAIN FRAMING
# =============================================================================

def normalize_vector(
    vector: np.ndarray,
) -> np.ndarray:

    vector = np.asarray(
        vector,
        dtype=float,
    )

    norm = np.linalg.norm(vector)

    if not np.isfinite(norm) or norm == 0:
        raise ValueError(
            f"Invalid camera vector: {vector}"
        )

    return vector / norm


def union_bounds(
    *bounds_list: Sequence[float],
) -> Tuple[float, float, float, float, float, float]:

    arrays = [
        np.asarray(bounds, dtype=float)
        for bounds in bounds_list
    ]

    if len(arrays) == 0:
        raise ValueError(
            "At least one bounds tuple is required."
        )

    return (
        min(float(b[0]) for b in arrays),
        max(float(b[1]) for b in arrays),
        min(float(b[2]) for b in arrays),
        max(float(b[3]) for b in arrays),
        min(float(b[4]) for b in arrays),
        max(float(b[5]) for b in arrays),
    )


def bounds_corners(
    bounds: Sequence[float],
) -> np.ndarray:

    xmin, xmax, ymin, ymax, zmin, zmax = [
        float(x)
        for x in bounds
    ]

    return np.asarray(
        [
            [x, y, z]
            for x in [xmin, xmax]
            for y in [ymin, ymax]
            for z in [zmin, zmax]
        ],
        dtype=float,
    )


def compute_full_brain_camera(
    bounds: Sequence[float],
    camera_direction: np.ndarray,
    preferred_view_up: np.ndarray,
    aspect_ratio: float,
    frame_margin: float,
) -> Tuple[
    Tuple[float, float, float],
    Tuple[float, float, float],
    Tuple[float, float, float],
    float,
]:
    """
    Construct an orthographic camera that contains the complete brain bounds.
    """

    if frame_margin <= 1.0:
        raise ValueError(
            "frame_margin must be > 1.0."
        )

    if aspect_ratio <= 0:
        raise ValueError(
            "aspect_ratio must be positive."
        )

    corners = bounds_corners(bounds)
    center = np.mean(corners, axis=0)

    direction = normalize_vector(
        camera_direction
    )

    up0 = normalize_vector(
        preferred_view_up
    )

    # Orthogonalize screen-up against viewing direction.
    up = up0 - np.dot(up0, direction) * direction

    # If almost parallel, choose a stable fallback.
    if np.linalg.norm(up) < 1e-8:
        fallback = np.array([1.0, 0.0, 0.0], dtype=float)
        up = fallback - np.dot(fallback, direction) * direction

    up = normalize_vector(up)

    right = normalize_vector(
        np.cross(direction, up)
    )

    centered = corners - center

    horizontal = centered @ right
    vertical = centered @ up

    half_width = float(
        np.max(np.abs(horizontal))
    )

    half_height = float(
        np.max(np.abs(vertical))
    )

    # Parallel-projection viewport:
    # vertical half extent   = parallel_scale
    # horizontal half extent = parallel_scale * aspect ratio
    parallel_scale = max(
        half_height,
        half_width / aspect_ratio,
    )

    parallel_scale *= frame_margin

    brain_radius = float(
        np.max(
            np.linalg.norm(
                centered,
                axis=1,
            )
        )
    )

    if not np.isfinite(brain_radius) or brain_radius <= 0:
        brain_radius = 100.0

    camera_distance = brain_radius * 4.0

    camera_position = (
        center
        + direction * camera_distance
    )

    return (
        tuple(camera_position.tolist()),
        tuple(center.tolist()),
        tuple(up.tolist()),
        float(parallel_scale),
    )


# =============================================================================
# PLOTTING
# =============================================================================

def scalar_bar_args(
    title: str,
) -> dict:

    return {
        "title": title,
        "title_font_size": 22,
        "label_font_size": 18,
        "n_labels": 5,
        "fmt": "%.2f",
        "font_family": "arial",
        "vertical": True,
        "position_x": 0.86,
        "position_y": 0.12,
        "height": 0.72,
        "width": 0.08,
    }


def add_background_surfaces(
    plotter: pv.Plotter,
    mesh_bg_left: pv.DataSet,
    mesh_bg_right: pv.DataSet,
    opacity: float,
) -> None:

    plotter.add_mesh(
        mesh_bg_left,
        color="lightgray",
        opacity=opacity,
        smooth_shading=True,
        show_scalar_bar=False,
    )

    plotter.add_mesh(
        mesh_bg_right,
        color="lightgray",
        opacity=opacity,
        smooth_shading=True,
        show_scalar_bar=False,
    )


def render_view(
    mesh_stat: Optional[pv.DataSet],
    mesh_bg_left: pv.DataSet,
    mesh_bg_right: pv.DataSet,
    full_brain_bounds: Sequence[float],
    camera_direction: np.ndarray,
    preferred_view_up: np.ndarray,
    output_png: Path,
    clim: Tuple[float, float],
    title: str,
    scalar_bar_title: str,
    background_opacity: float,
    tract_opacity: float,
    frame_margin: float,
    window_size: Tuple[int, int],
) -> None:

    pv.set_plot_theme(
        "document"
    )

    plotter = pv.Plotter(
        off_screen=True,
        window_size=window_size,
    )

    plotter.set_background(
        "white"
    )

    if mesh_stat is not None:
        plotter.add_mesh(
            mesh_stat,
            scalars=SCALAR_NAME,
            preference="cell",
            clim=list(clim),
            cmap="coolwarm",
            opacity=tract_opacity,
            smooth_shading=True,
            show_scalar_bar=True,
            scalar_bar_args=scalar_bar_args(
                scalar_bar_title
            ),
        )

    add_background_surfaces(
        plotter,
        mesh_bg_left,
        mesh_bg_right,
        background_opacity,
    )

    if title:
        plotter.add_text(
            title,
            position="upper_left",
            font_size=13,
            color="black",
            font="arial",
        )

    if mesh_stat is not None:
        plotter.add_text(
            "beta > 0: higher FA associated with CRA",
            position="lower_left",
            font_size=10,
            color="black",
            font="arial",
        )

    aspect_ratio = (
        float(window_size[0])
        / float(window_size[1])
    )

    (
        camera_position,
        focal_point,
        view_up,
        parallel_scale,
    ) = compute_full_brain_camera(
        bounds=full_brain_bounds,
        camera_direction=camera_direction,
        preferred_view_up=preferred_view_up,
        aspect_ratio=aspect_ratio,
        frame_margin=frame_margin,
    )

    plotter.camera.parallel_projection = True
    plotter.camera.position = camera_position
    plotter.camera.focal_point = focal_point
    plotter.camera.up = view_up
    plotter.camera.parallel_scale = parallel_scale

    plotter.reset_camera_clipping_range()

    print(
        f"  camera direction={tuple(np.round(camera_direction, 4))}"
    )
    print(
        f"  camera position={tuple(round(x, 3) for x in camera_position)}"
    )
    print(
        f"  parallel scale={parallel_scale:.3f}"
    )

    plotter.screenshot(
        str(output_png),
        transparent_background=False,
        return_img=False,
    )

    plotter.close()


# =============================================================================
# COMPOSITE
# =============================================================================

def make_composite(
    png_paths: List[Path],
    panel_titles: List[str],
    output_png: Path,
) -> None:

    fig, axes = plt.subplots(
        1,
        len(png_paths),
        figsize=(18, 6),
    )

    if len(png_paths) == 1:
        axes = [axes]

    for ax, path, label in zip(
        axes,
        png_paths,
        panel_titles,
    ):
        ax.imshow(
            imread(path)
        )
        ax.axis(
            "off"
        )
        ax.set_title(
            label,
            fontsize=13,
        )

    fig.tight_layout(
        pad=0.3
    )

    fig.savefig(
        output_png,
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )

    plt.close(fig)


# =============================================================================
# OUTPUT TABLE / SUMMARY
# =============================================================================

def write_significant_table(
    significant: pd.DataFrame,
    output_tsv: Path,
) -> None:

    preferred = [
        "Tract_number",
        "Tract_abb",
        "IDP",
        "Feature_family",
        BETA_COL,
        "SE",
        "Z",
        RAW_P_COL,
        FDR_COL,
        OR_COL,
        "OR_CI95_lower",
        "OR_CI95_upper",
        "Beta_direction",
        "N_case",
        "N_control",
        "N_total",
        "phenotype_column",
        "case_indicator_column",
        "control_indicator_column",
        "model_formula",
    ]

    ordered = [
        col
        for col in preferred
        if col in significant.columns
    ]

    remaining = [
        col
        for col in significant.columns
        if col not in ordered
    ]

    significant[
        ordered + remaining
    ].to_csv(
        output_tsv,
        sep="\t",
        index=False,
        na_rep="NA",
    )


def write_summary(
    output_txt: Path,
    input_tsv: Path,
    all_fa: pd.DataFrame,
    significant: pd.DataFrame,
    fdr_threshold: float,
    clim: Tuple[float, float],
    frame_margin: float,
    oblique_angle: float,
) -> None:

    lines = [
        "CRA vs CUA DTI FA FDR visualization",
        "=" * 76,
        f"Input TSV: {input_tsv}",
        f"Significance: {FDR_COL} < {fdr_threshold}",
        f"N FA tests: {len(all_fa)}",
        f"N FDR-significant tracts: {len(significant)}",
        f"Visualization scalar: {BETA_COL}",
        "Outcome coding: CRA=1, CUA=0",
        "Positive beta: higher FA associated with CRA",
        "Negative beta: higher FA associated with CUA",
        f"Oblique elevation angle: {oblique_angle} degrees",
        f"Full-brain frame margin: {frame_margin}",
        f"Symmetric beta scale: [{clim[0]:.6f}, {clim[1]:.6f}]",
        "Views: superior + bilateral superior-lateral oblique",
        "All views use the same FDR-significant tract mesh.",
        "All views use automatic orthographic projected-bound framing.",
        "",
    ]

    if significant.empty:
        lines.append(
            "No FA tract survived within-modality BH-FDR."
        )
    else:
        lines.append(
            "FDR-significant tracts:"
        )

        for _, row in significant.iterrows():
            lines.append(
                f"  {int(row['Tract_number']):02d} "
                f"{row['Tract_abb']}: "
                f"beta={row[BETA_COL]:+.6f}, "
                f"OR={row[OR_COL]:.4f}, "
                f"P={row[RAW_P_COL]:.6g}, "
                f"FDR={row[FDR_COL]:.6g}"
            )

    output_txt.write_text(
        "\n".join(lines)
        + "\n"
    )


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:

    args = parse_args()

    if not (
        0 < args.fdr_threshold < 1
    ):
        raise ValueError(
            "--fdr-threshold must be between 0 and 1."
        )

    if not (
        0 < args.oblique_angle < 90
    ):
        raise ValueError(
            "--oblique-angle must be > 0 and < 90."
        )

    if not (
        0 <= args.background_opacity <= 1
    ):
        raise ValueError(
            "--background-opacity must be between 0 and 1."
        )

    if not (
        0 <= args.tract_opacity <= 1
    ):
        raise ValueError(
            "--tract-opacity must be between 0 and 1."
        )

    if args.frame_margin <= 1.0:
        raise ValueError(
            "--frame-margin must be > 1.0."
        )

    input_tsv = args.input_tsv.expanduser()
    dictionary_tsv = args.dictionary_tsv.expanduser()
    tract_dir = args.tract_dir.expanduser()
    bg_left = args.bg_left.expanduser()
    bg_right = args.bg_right.expanduser()
    output_dir = args.output_dir.expanduser()

    require_file(
        input_tsv,
        "Input DTI statistics TSV",
    )

    require_file(
        dictionary_tsv,
        "JHU dictionary",
    )

    require_dir(
        tract_dir,
        "JHU tract directory",
    )

    require_file(
        bg_left,
        "Left cortical background VTK",
    )

    require_file(
        bg_right,
        "Right cortical background VTK",
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    print("=" * 80)
    print("CRA vs CUA DTI FA FDR oblique brain visualization")
    print("=" * 80)
    print(f"Input TSV:       {input_tsv}")
    print(f"Significance:    {FDR_COL} < {args.fdr_threshold}")
    print(f"Effect plotted:  {BETA_COL}")
    print(f"Oblique angle:   {args.oblique_angle} degrees")
    print(f"Output dir:      {output_dir}")
    print()

    all_fa, significant = prepare_results(
        input_tsv=input_tsv,
        dictionary_tsv=dictionary_tsv,
        feature_family=args.feature_family,
        fdr_threshold=args.fdr_threshold,
    )

    print(
        f"Eligible DTI features: {len(all_fa)}"
    )
    print(
        f"FDR-significant tracts: {len(significant)}"
    )
    print()

    significant_tsv = (
        output_dir
        / "CRA_vs_CUA_DTI_FA_FDR_significant_tracts.tsv"
    )

    write_significant_table(
        significant,
        significant_tsv,
    )

    mesh_bg_left = pv.read(
        str(bg_left)
    )

    mesh_bg_right = pv.read(
        str(bg_right)
    )

    full_brain_bounds = union_bounds(
        mesh_bg_left.bounds,
        mesh_bg_right.bounds,
    )

    mesh_stat, _ = build_significant_mesh(
        significant,
        tract_dir,
    )

    clim = symmetric_beta_limits(
        significant
    )

    if mesh_stat is not None:
        combined_vtk = (
            output_dir
            / "CRA_vs_CUA_DTI_FA_FDR_significant_beta.vtk"
        )

        mesh_stat.save(
            str(combined_vtk)
        )

        print(
            f"Saved combined mesh: {combined_vtk}"
        )

    view_specs = build_view_specs(
        args.oblique_angle
    )

    screenshot_paths: List[Path] = []
    panel_titles: List[str] = []

    for view_name, (
        camera_direction,
        preferred_view_up,
    ) in view_specs.items():

        output_png = (
            output_dir
            / f"CRA_vs_CUA_DTI_FA_FDR_beta_{view_name}.png"
        )

        print()
        print(
            f"Rendering {view_name}:"
        )

        render_view(
            mesh_stat=mesh_stat,
            mesh_bg_left=mesh_bg_left,
            mesh_bg_right=mesh_bg_right,
            full_brain_bounds=full_brain_bounds,
            camera_direction=camera_direction,
            preferred_view_up=preferred_view_up,
            output_png=output_png,
            clim=clim,
            title=args.title,
            scalar_bar_title=args.scalar_bar_title,
            background_opacity=args.background_opacity,
            tract_opacity=args.tract_opacity,
            frame_margin=args.frame_margin,
            window_size=(
                args.window_width,
                args.window_height,
            ),
        )

        screenshot_paths.append(
            output_png
        )

        if view_name == "superior":
            panel_titles.append(
                "Superior"
            )
        elif view_name.startswith("left_"):
            panel_titles.append(
                f"Left oblique ({args.oblique_angle:g} degrees)"
            )
        else:
            panel_titles.append(
                f"Right oblique ({args.oblique_angle:g} degrees)"
            )

        print(
            f"Saved: {output_png}"
        )

    if not args.no_composite:
        composite_png = (
            output_dir
            / "CRA_vs_CUA_DTI_FA_FDR_beta_oblique_composite.png"
        )

        make_composite(
            screenshot_paths,
            panel_titles,
            composite_png,
        )

        print()
        print(
            f"Saved: {composite_png}"
        )

    summary_txt = (
        output_dir
        / "CRA_vs_CUA_DTI_FA_FDR_summary.txt"
    )

    write_summary(
        output_txt=summary_txt,
        input_tsv=input_tsv,
        all_fa=all_fa,
        significant=significant,
        fdr_threshold=args.fdr_threshold,
        clim=clim,
        frame_margin=args.frame_margin,
        oblique_angle=args.oblique_angle,
    )

    print(
        f"Saved: {summary_txt}"
    )

    print()
    print("Interpretation:")
    print("  beta > 0: higher FA associated with CRA")
    print("  beta < 0: higher FA associated with CUA")
    print("  only within-modality BH-FDR significant tracts are rendered")
    print("  the same significant tracts are used in all three views")
    print("  oblique views retain full-brain automatic framing")


if __name__ == "__main__":
    main()