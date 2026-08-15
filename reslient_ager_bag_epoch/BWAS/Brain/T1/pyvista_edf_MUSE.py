#!/usr/bin/env python3
"""
Visualize within-modality BH-FDR-significant T1 gray-matter associations
for brain-proteomics EPOCH-BAG candidate resilient agers (CRA) versus
concordant unfavorable agers (CUA).

SIGNIFICANCE
------------
Only regions satisfying

    P_value_FDR_BH_within_modality < 0.05

are rendered.

The exact same significant-region mesh is used for every brain view.

COLOR
-----
The visualization uses the signed logistic-regression coefficient

    Beta_log_odds_per_1SD_IDP

with outcome coding

    CRA = 1
    CUA = 0

Therefore

    beta > 0  -> larger regional GM volume is associated with greater odds of CRA
    beta < 0  -> larger regional GM volume is associated with greater odds of CUA

A symmetric zero-centered diverging color scale is used.

FULL-BRAIN CAMERA FRAMING
-------------------------
This revision does NOT use fixed camera coordinates. Instead, every view is
framed automatically from the bounds of MUSE_GM_all.vtk using parallel
projection. The projected brain bounds are calculated for each camera
orientation and the parallel scale is set with a margin.

This prevents one lateral view from becoming excessively zoomed/cropped and
ensures that superior, left-lateral, and right-lateral views all contain the
full brain at a comparable scale.

MUSE ROI MAPPING
----------------
An IDP such as

    MUSE_Volume_140

is mapped to

    MUSE_GM_140.vtk

in the MUSE VTK directory.

DEFAULT INPUT
-------------
/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/mortality_clock/BWAS/
Brain_proteomics_EPOCH_BAG_resilience/CRA_vs_CUA/combined_results/
BWAS_CRA_vs_CUA_T1_all_statistics.tsv

OUTPUTS
-------
CRA_vs_CUA_T1_GM_FDR_significant_regions.tsv
CRA_vs_CUA_T1_GM_FDR_significant_beta.vtk

CRA_vs_CUA_T1_GM_FDR_beta_superior.png
CRA_vs_CUA_T1_GM_FDR_beta_left_lateral.png
CRA_vs_CUA_T1_GM_FDR_beta_right_lateral.png

CRA_vs_CUA_T1_GM_FDR_beta_composite.png
CRA_vs_CUA_T1_GM_FDR_summary.txt

USAGE
-----
python plot_CRA_CUA_T1_GM_FDR_brain_visualization.py

Optional:
python plot_CRA_CUA_T1_GM_FDR_brain_visualization.py \
    --input-tsv /path/to/BWAS_CRA_vs_CUA_T1_all_statistics.tsv \
    --output-dir /path/to/output
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
    "BWAS_CRA_vs_CUA_T1_all_statistics.tsv"
)

DEFAULT_MUSE_VTK_DIR = Path(
    "/Users/hao/cubic-home/Dataset/Atlas/MUSE/VTK_GM"
)

DEFAULT_BACKGROUND_VTK = (
    DEFAULT_MUSE_VTK_DIR / "MUSE_GM_all.vtk"
)

DEFAULT_OUTPUT_DIR = Path(
    "/Users/hao/Dropbox/2026_EPOCH/Fig/orig/BWAS/T1"
)


# =============================================================================
# ANALYSIS COLUMNS
# =============================================================================

MODALITY = "T1"
FEATURE_FAMILY = "MUSE_GM_volume"

FDR_COL = "P_value_FDR_BH_within_modality"
BETA_COL = "Beta_log_odds_per_1SD_IDP"
RAW_P_COL = "P_value"
OR_COL = "OR_per_1SD_IDP"

SCALAR_NAME = "CRA_vs_CUA_beta"

DEFAULT_FDR_THRESHOLD = 0.05


# =============================================================================
# STANDARD BRAIN VIEWS
# =============================================================================
#
# MUSE is assumed to use a standard anatomical coordinate system in which:
#     +X / -X = lateral views
#     +Z       = superior view
#
# The camera distance and zoom are NOT hard-coded. They are derived from the
# full-brain background bounds below.
#
# Each entry:
#     camera_direction = vector from brain center TO camera
#     preferred_view_up = desired vertical direction on screen
# =============================================================================

VIEW_SPECS: Dict[str, Tuple[np.ndarray, np.ndarray]] = {
    "superior": (
        np.array([0.0, 0.0, 1.0], dtype=float),
        np.array([0.0, 1.0, 0.0], dtype=float),
    ),
    "left_lateral": (
        np.array([-1.0, 0.0, 0.0], dtype=float),
        np.array([0.0, 0.0, 1.0], dtype=float),
    ),
    "right_lateral": (
        np.array([1.0, 0.0, 0.0], dtype=float),
        np.array([0.0, 0.0, 1.0], dtype=float),
    ),
}


# =============================================================================
# ARGUMENTS
# =============================================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Visualize within-modality BH-FDR-significant T1 MUSE GM "
            "associations for CRA vs CUA."
        )
    )

    parser.add_argument(
        "--input-tsv",
        type=Path,
        default=DEFAULT_INPUT_TSV,
        help="Collected T1 statistics TSV.",
    )

    parser.add_argument(
        "--muse-vtk-dir",
        type=Path,
        default=DEFAULT_MUSE_VTK_DIR,
        help="Directory containing MUSE_GM_<ROI>.vtk meshes.",
    )

    parser.add_argument(
        "--background-vtk",
        type=Path,
        default=DEFAULT_BACKGROUND_VTK,
        help="Full-brain MUSE background VTK.",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Output directory.",
    )

    parser.add_argument(
        "--fdr-threshold",
        type=float,
        default=DEFAULT_FDR_THRESHOLD,
        help=(
            "Threshold for P_value_FDR_BH_within_modality. "
            "Default: 0.05."
        ),
    )

    parser.add_argument(
        "--background-opacity",
        type=float,
        default=0.08,
        help="Opacity of MUSE_GM_all.vtk. Default: 0.08.",
    )

    parser.add_argument(
        "--roi-opacity",
        type=float,
        default=1.0,
        help="Opacity of significant MUSE regions. Default: 1.0.",
    )

    parser.add_argument(
        "--frame-margin",
        type=float,
        default=1.12,
        help=(
            "Multiplicative margin around projected full-brain bounds. "
            "Values >1 zoom out. Default: 1.12."
        ),
    )

    parser.add_argument(
        "--window-width",
        type=int,
        default=1800,
        help="Screenshot width in pixels.",
    )

    parser.add_argument(
        "--window-height",
        type=int,
        default=1400,
        help="Screenshot height in pixels.",
    )

    parser.add_argument(
        "--title",
        default="CRA vs CUA: FDR-significant T1 gray-matter associations",
        help="Figure title. Pass an empty string to suppress.",
    )

    parser.add_argument(
        "--scalar-bar-title",
        default="Beta (CRA vs CUA)",
        help="Color-bar title.",
    )

    parser.add_argument(
        "--no-composite",
        action="store_true",
        help="Do not create the three-panel composite PNG.",
    )

    return parser.parse_args()


# =============================================================================
# BASIC VALIDATION
# =============================================================================

def require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(
            f"{label} not found:\n  {path}"
        )


def require_dir(path: Path, label: str) -> None:
    if not path.is_dir():
        raise FileNotFoundError(
            f"{label} not found:\n  {path}"
        )


def require_columns(
    df: pd.DataFrame,
    columns: Sequence[str],
    label: str,
) -> None:
    missing = [
        col
        for col in columns
        if col not in df.columns
    ]

    if missing:
        raise ValueError(
            f"{label} is missing required column(s): "
            + ", ".join(missing)
        )


# =============================================================================
# MUSE IDP -> VTK
# =============================================================================

def muse_roi_number(idp: str) -> Optional[int]:
    """
    MUSE_Volume_140 -> 140
    """

    match = re.fullmatch(
        r"MUSE_Volume_(\d+)",
        str(idp).strip(),
    )

    if match is None:
        return None

    return int(
        match.group(1)
    )


def muse_vtk_path(
    muse_vtk_dir: Path,
    idp: str,
) -> Path:
    roi = muse_roi_number(
        idp
    )

    if roi is None:
        raise ValueError(
            f"Cannot parse MUSE ROI number from IDP: {idp}"
        )

    return (
        muse_vtk_dir
        / f"MUSE_GM_{roi}.vtk"
    )


# =============================================================================
# READ AND FILTER STATISTICS
# =============================================================================

def prepare_results(
    input_tsv: Path,
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
        "T1 statistics TSV",
    )

    # -------------------------------------------------------------------------
    # Restrict exactly to the T1 MUSE GM analysis.
    # -------------------------------------------------------------------------

    stats = stats.loc[
        stats["Modality"]
        .astype(str)
        .str.upper()
        .eq(MODALITY)
    ].copy()

    stats = stats.loc[
        stats["Feature_family"]
        .astype(str)
        .eq(FEATURE_FAMILY)
    ].copy()

    if stats.empty:
        raise ValueError(
            "No rows remain after filtering to "
            f"Modality={MODALITY} and Feature_family={FEATURE_FAMILY}."
        )

    # -------------------------------------------------------------------------
    # Numeric conversion.
    # -------------------------------------------------------------------------

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
            stats[col] = pd.to_numeric(
                stats[col],
                errors="coerce",
            )

    # -------------------------------------------------------------------------
    # Parse MUSE ROI numbers and ensure every requested T1 row is mappable.
    # -------------------------------------------------------------------------

    stats["MUSE_ROI_number"] = stats[
        "IDP"
    ].map(
        muse_roi_number
    )

    unmapped = stats.loc[
        stats["MUSE_ROI_number"].isna(),
        "IDP",
    ]

    if len(unmapped) > 0:
        raise ValueError(
            "Some T1 IDPs are not valid MUSE_Volume_<number> names:\n  "
            + "\n  ".join(
                unmapped.astype(str)
            )
        )

    stats["MUSE_ROI_number"] = stats[
        "MUSE_ROI_number"
    ].astype(int)

    # -------------------------------------------------------------------------
    # Requested significance definition.
    #
    # This is the ONLY inferential significance filter.
    # There is NO:
    #   raw P threshold
    #   OR > 1 filter
    #   beta > 0 filter
    #   BAG filter
    # -------------------------------------------------------------------------

    significant = stats.loc[
        stats[FDR_COL].notna()
        & (stats[FDR_COL] < fdr_threshold)
    ].copy()

    significant = significant.sort_values(
        by=[
            FDR_COL,
            RAW_P_COL,
            "MUSE_ROI_number",
        ],
        ascending=[
            True,
            True,
            True,
        ],
        kind="mergesort",
    ).reset_index(
        drop=True
    )

    return stats, significant


# =============================================================================
# BUILD SIGNIFICANT REGION MESH
# =============================================================================

def build_significant_mesh(
    significant: pd.DataFrame,
    muse_vtk_dir: Path,
) -> Tuple[Optional[pv.DataSet], List[Path]]:

    if significant.empty:
        return None, []

    meshes: List[pv.DataSet] = []
    used_paths: List[Path] = []

    for _, row in significant.iterrows():

        idp = str(
            row["IDP"]
        )

        roi_number = int(
            row["MUSE_ROI_number"]
        )

        beta = float(
            row[BETA_COL]
        )

        raw_p = float(
            row[RAW_P_COL]
        )

        fdr = float(
            row[FDR_COL]
        )

        odds_ratio = float(
            row[OR_COL]
        )

        vtk_file = muse_vtk_path(
            muse_vtk_dir,
            idp,
        )

        require_file(
            vtk_file,
            f"MUSE VTK for {idp}",
        )

        mesh = pv.read(
            str(vtk_file)
        ).copy()

        if mesh.n_cells <= 0:
            raise ValueError(
                f"MUSE mesh has zero cells: {vtk_file}"
            )

        mesh.cell_data[
            SCALAR_NAME
        ] = np.full(
            mesh.n_cells,
            beta,
            dtype=float,
        )

        meshes.append(
            mesh
        )

        used_paths.append(
            vtk_file
        )

        print(
            f"{idp:<20} "
            f"ROI={roi_number:03d} "
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

    return (
        mesh_final,
        used_paths,
    )


# =============================================================================
# ZERO-CENTERED COLOR LIMITS
# =============================================================================

def symmetric_beta_limits(
    significant: pd.DataFrame,
) -> Tuple[float, float]:

    beta = pd.to_numeric(
        significant[BETA_COL],
        errors="coerce",
    ).to_numpy(
        dtype=float
    )

    beta = beta[
        np.isfinite(beta)
    ]

    if beta.size == 0:
        return (
            -1.0,
            1.0,
        )

    max_abs = float(
        np.max(
            np.abs(beta)
        )
    )

    if max_abs <= 0:
        max_abs = 1.0

    max_abs *= 1.05

    return (
        -max_abs,
        max_abs,
    )


# =============================================================================
# AUTOMATIC FULL-BRAIN CAMERA
# =============================================================================

def normalize_vector(
    vector: np.ndarray,
) -> np.ndarray:

    vector = np.asarray(
        vector,
        dtype=float,
    )

    norm = np.linalg.norm(
        vector
    )

    if not np.isfinite(norm) or norm == 0:
        raise ValueError(
            f"Invalid camera vector: {vector}"
        )

    return (
        vector / norm
    )


def bounds_corners(
    bounds: Sequence[float],
) -> np.ndarray:
    """
    Convert PyVista bounds

        (xmin, xmax, ymin, ymax, zmin, zmax)

    to the eight bounding-box corners.
    """

    xmin, xmax, ymin, ymax, zmin, zmax = [
        float(x)
        for x in bounds
    ]

    corners = []

    for x in [
        xmin,
        xmax,
    ]:
        for y in [
            ymin,
            ymax,
        ]:
            for z in [
                zmin,
                zmax,
            ]:
                corners.append(
                    [
                        x,
                        y,
                        z,
                    ]
                )

    return np.asarray(
        corners,
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
    Construct an orthographic camera that contains the entire supplied bounds.

    Returns
    -------
    camera_position
    focal_point
    view_up
    parallel_scale

    The parallel_scale explicitly accounts for both projected vertical and
    projected horizontal brain extent, so the whole brain remains visible
    regardless of screenshot aspect ratio.
    """

    if frame_margin <= 1.0:
        raise ValueError(
            "frame_margin should be > 1.0 to leave visible space "
            "around the full brain."
        )

    if aspect_ratio <= 0:
        raise ValueError(
            "aspect_ratio must be positive."
        )

    corners = bounds_corners(
        bounds
    )

    center = np.mean(
        corners,
        axis=0,
    )

    # Direction from brain center toward camera.
    d = normalize_vector(
        camera_direction
    )

    # Make view-up orthogonal to the viewing direction.
    up0 = normalize_vector(
        preferred_view_up
    )

    up = (
        up0
        - np.dot(
            up0,
            d,
        )
        * d
    )

    up = normalize_vector(
        up
    )

    # Horizontal screen axis.
    right = normalize_vector(
        np.cross(
            d,
            up,
        )
    )

    centered = (
        corners
        - center
    )

    horizontal_projection = (
        centered
        @ right
    )

    vertical_projection = (
        centered
        @ up
    )

    half_width = float(
        np.max(
            np.abs(
                horizontal_projection
            )
        )
    )

    half_height = float(
        np.max(
            np.abs(
                vertical_projection
            )
        )
    )

    # In parallel projection:
    #
    #     visible vertical half-height = parallel_scale
    #     visible horizontal half-width = parallel_scale * aspect_ratio
    #
    # Choose scale large enough for BOTH dimensions.
    parallel_scale = max(
        half_height,
        half_width / aspect_ratio,
    )

    parallel_scale *= frame_margin

    # Camera distance does not determine scale in parallel projection, but it
    # must be sufficiently far away to avoid clipping.
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

    camera_distance = (
        brain_radius
        * 4.0
    )

    camera_position = (
        center
        + d
        * camera_distance
    )

    return (
        tuple(
            camera_position.tolist()
        ),
        tuple(
            center.tolist()
        ),
        tuple(
            up.tolist()
        ),
        float(
            parallel_scale
        ),
    )


# =============================================================================
# PYVISTA PLOTTING
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


def render_view(
    mesh_stat: Optional[pv.DataSet],
    mesh_bg: pv.DataSet,
    camera_direction: np.ndarray,
    view_up: np.ndarray,
    output_png: Path,
    clim: Tuple[float, float],
    title: str,
    scalar_bar_title: str,
    background_opacity: float,
    roi_opacity: float,
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

    # Add significant regions first.
    if mesh_stat is not None:

        plotter.add_mesh(
            mesh_stat,
            scalars=SCALAR_NAME,
            preference="cell",
            clim=list(
                clim
            ),
            cmap="coolwarm",
            opacity=roi_opacity,
            smooth_shading=True,
            show_scalar_bar=True,
            scalar_bar_args=scalar_bar_args(
                scalar_bar_title
            ),
        )

    # Full MUSE brain background.
    plotter.add_mesh(
        mesh_bg,
        color="lightgray",
        opacity=background_opacity,
        smooth_shading=True,
        show_scalar_bar=False,
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
            "beta > 0: larger GM volume associated with CRA",
            position="lower_left",
            font_size=10,
            color="black",
            font="arial",
        )

    aspect_ratio = (
        float(
            window_size[0]
        )
        / float(
            window_size[1]
        )
    )

    (
        camera_position,
        focal_point,
        camera_view_up,
        parallel_scale,
    ) = compute_full_brain_camera(
        bounds=mesh_bg.bounds,
        camera_direction=camera_direction,
        preferred_view_up=view_up,
        aspect_ratio=aspect_ratio,
        frame_margin=frame_margin,
    )

    # Orthographic framing gives identical geometry-based zoom logic to every
    # view and avoids perspective-dependent lateral cropping.
    plotter.camera.parallel_projection = True

    plotter.camera.position = (
        camera_position
    )

    plotter.camera.focal_point = (
        focal_point
    )

    plotter.camera.up = (
        camera_view_up
    )

    plotter.camera.parallel_scale = (
        parallel_scale
    )

    plotter.reset_camera_clipping_range()

    plotter.screenshot(
        str(
            output_png
        ),
        transparent_background=False,
        return_img=False,
    )

    plotter.close()


# =============================================================================
# COMPOSITE FIGURE
# =============================================================================

def make_composite(
    png_paths: List[Path],
    panel_titles: List[str],
    output_png: Path,
) -> None:

    if len(png_paths) == 0:
        return

    fig, axes = plt.subplots(
        1,
        len(
            png_paths
        ),
        figsize=(
            18,
            6,
        ),
    )

    if len(
        png_paths
    ) == 1:
        axes = [
            axes
        ]

    for ax, png, label in zip(
        axes,
        png_paths,
        panel_titles,
    ):
        image = imread(
            png
        )

        ax.imshow(
            image
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

    plt.close(
        fig
    )


# =============================================================================
# TABLE / SUMMARY OUTPUT
# =============================================================================

def write_significant_table(
    significant: pd.DataFrame,
    output_tsv: Path,
) -> None:

    preferred = [
        "MUSE_ROI_number",
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
        ordered
        + remaining
    ].to_csv(
        output_tsv,
        sep="\t",
        index=False,
        na_rep="NA",
    )


def write_summary(
    output_txt: Path,
    input_tsv: Path,
    all_t1: pd.DataFrame,
    significant: pd.DataFrame,
    fdr_threshold: float,
    clim: Tuple[float, float],
    frame_margin: float,
) -> None:

    lines: List[str] = [
        "CRA vs CUA T1 gray-matter FDR visualization",
        "=" * 76,
        f"Input TSV: {input_tsv}",
        f"Modality filter: {MODALITY}",
        f"Feature family: {FEATURE_FAMILY}",
        f"Significance definition: {FDR_COL} < {fdr_threshold}",
        f"N T1 MUSE tests: {len(all_t1)}",
        f"N FDR-significant regions: {len(significant)}",
        f"Visualization scalar: {BETA_COL}",
        "Outcome coding: CRA=1, CUA=0",
        "Positive beta: larger GM volume associated with greater odds of CRA",
        "Negative beta: larger GM volume associated with greater odds of CUA",
        f"Symmetric beta color scale: [{clim[0]:.6f}, {clim[1]:.6f}]",
        f"Full-brain frame margin: {frame_margin}",
        "All rendered views use the exact same FDR-significant region mesh.",
        "All views use automatic orthographic full-brain framing.",
        "",
    ]

    if significant.empty:

        lines.append(
            "No T1 MUSE region survived the requested within-modality BH-FDR threshold."
        )

    else:

        lines.append(
            "FDR-significant MUSE regions:"
        )

        for _, row in significant.iterrows():

            lines.append(
                (
                    f"  ROI {int(row['MUSE_ROI_number']):03d} "
                    f"{row['IDP']}: "
                    f"beta={row[BETA_COL]:+.6f}, "
                    f"OR={row[OR_COL]:.6f}, "
                    f"P={row[RAW_P_COL]:.6g}, "
                    f"FDR={row[FDR_COL]:.6g}"
                )
            )

    output_txt.write_text(
        "\n".join(
            lines
        )
        + "\n"
    )


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:

    args = parse_args()

    if not (
        0
        < args.fdr_threshold
        < 1
    ):
        raise ValueError(
            "--fdr-threshold must be between 0 and 1."
        )

    if not (
        0
        <= args.background_opacity
        <= 1
    ):
        raise ValueError(
            "--background-opacity must be between 0 and 1."
        )

    if not (
        0
        <= args.roi_opacity
        <= 1
    ):
        raise ValueError(
            "--roi-opacity must be between 0 and 1."
        )

    if args.frame_margin <= 1.0:
        raise ValueError(
            "--frame-margin must be > 1.0."
        )

    input_tsv = (
        args.input_tsv
        .expanduser()
    )

    muse_vtk_dir = (
        args.muse_vtk_dir
        .expanduser()
    )

    background_vtk = (
        args.background_vtk
        .expanduser()
    )

    output_dir = (
        args.output_dir
        .expanduser()
    )

    require_file(
        input_tsv,
        "T1 statistics TSV",
    )

    require_dir(
        muse_vtk_dir,
        "MUSE VTK directory",
    )

    require_file(
        background_vtk,
        "MUSE full-brain background VTK",
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    print(
        "=" * 80
    )

    print(
        "CRA vs CUA T1 MUSE GM FDR brain visualization"
    )

    print(
        "=" * 80
    )

    print(
        f"Input TSV:       {input_tsv}"
    )

    print(
        f"Significance:    {FDR_COL} < {args.fdr_threshold}"
    )

    print(
        f"Effect plotted:  {BETA_COL}"
    )

    print(
        f"MUSE VTK dir:    {muse_vtk_dir}"
    )

    print(
        f"Background:      {background_vtk}"
    )

    print(
        f"Output dir:      {output_dir}"
    )

    print()

    all_t1, significant = prepare_results(
        input_tsv=input_tsv,
        fdr_threshold=args.fdr_threshold,
    )

    print(
        f"Eligible T1 MUSE features: {len(all_t1)}"
    )

    print(
        f"FDR-significant regions: {len(significant)}"
    )

    print()

    # -------------------------------------------------------------------------
    # Save the exact FDR-significant region set used in every view.
    # -------------------------------------------------------------------------

    significant_tsv = (
        output_dir
        / "CRA_vs_CUA_T1_GM_FDR_significant_regions.tsv"
    )

    write_significant_table(
        significant=significant,
        output_tsv=significant_tsv,
    )

    # -------------------------------------------------------------------------
    # Full-brain background.
    # -------------------------------------------------------------------------

    mesh_bg = pv.read(
        str(
            background_vtk
        )
    )

    # -------------------------------------------------------------------------
    # Build ONE significant-region mesh. Reuse it identically in all views.
    # -------------------------------------------------------------------------

    mesh_stat, used_vtks = build_significant_mesh(
        significant=significant,
        muse_vtk_dir=muse_vtk_dir,
    )

    clim = symmetric_beta_limits(
        significant
    )

    if mesh_stat is not None:

        combined_vtk = (
            output_dir
            / "CRA_vs_CUA_T1_GM_FDR_significant_beta.vtk"
        )

        mesh_stat.save(
            str(
                combined_vtk
            )
        )

        print()

        print(
            "Saved combined FDR-significant region mesh:"
        )

        print(
            f"  {combined_vtk}"
        )

    # -------------------------------------------------------------------------
    # Render standard full-brain views.
    # -------------------------------------------------------------------------

    screenshot_paths: List[Path] = []

    panel_titles: List[str] = []

    for view_name, (
        camera_direction,
        preferred_view_up,
    ) in VIEW_SPECS.items():

        output_png = (
            output_dir
            / f"CRA_vs_CUA_T1_GM_FDR_beta_{view_name}.png"
        )

        render_view(
            mesh_stat=mesh_stat,
            mesh_bg=mesh_bg,
            camera_direction=camera_direction,
            view_up=preferred_view_up,
            output_png=output_png,
            clim=clim,
            title=args.title,
            scalar_bar_title=args.scalar_bar_title,
            background_opacity=args.background_opacity,
            roi_opacity=args.roi_opacity,
            frame_margin=args.frame_margin,
            window_size=(
                args.window_width,
                args.window_height,
            ),
        )

        screenshot_paths.append(
            output_png
        )

        panel_titles.append(
            view_name
            .replace(
                "_",
                " ",
            )
            .title()
        )

        print(
            f"Saved: {output_png}"
        )

    # -------------------------------------------------------------------------
    # Composite.
    # -------------------------------------------------------------------------

    if not args.no_composite:

        composite_png = (
            output_dir
            / "CRA_vs_CUA_T1_GM_FDR_beta_composite.png"
        )

        make_composite(
            png_paths=screenshot_paths,
            panel_titles=panel_titles,
            output_png=composite_png,
        )

        print(
            f"Saved: {composite_png}"
        )

    # -------------------------------------------------------------------------
    # Summary.
    # -------------------------------------------------------------------------

    summary_txt = (
        output_dir
        / "CRA_vs_CUA_T1_GM_FDR_summary.txt"
    )

    write_summary(
        output_txt=summary_txt,
        input_tsv=input_tsv,
        all_t1=all_t1,
        significant=significant,
        fdr_threshold=args.fdr_threshold,
        clim=clim,
        frame_margin=args.frame_margin,
    )

    print(
        f"Saved: {summary_txt}"
    )

    print(
        f"Saved exact significant-region table: {significant_tsv}"
    )

    print()

    print(
        "Interpretation:"
    )

    print(
        "  beta > 0: larger GM volume associated with greater odds of CRA"
    )

    print(
        "  beta < 0: larger GM volume associated with greater odds of CUA"
    )

    print(
        "  significance is defined only by within-modality BH-FDR < threshold"
    )

    print(
        "  superior, left-lateral, and right-lateral views use the SAME significant mesh"
    )

    print(
        "  all views are automatically framed to contain the full brain"
    )


if __name__ == "__main__":
    main()