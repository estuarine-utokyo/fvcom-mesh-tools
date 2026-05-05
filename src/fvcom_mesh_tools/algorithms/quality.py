"""Per-element mesh-quality metrics.

The metrics here are CRS-agnostic: they treat ``mesh.nodes`` as a 2-D
point cloud in whatever units it happens to be in. For lon/lat meshes
the absolute edge-length numbers are degrees-of-arc rather than metres,
but the *ratios* used by ``alpha_quality`` and ``min_interior_angle``
are scale-free, so before/after comparisons remain meaningful.

Use the dedicated haversine / projected helpers in notebooks 06 / 07
when an absolute metric distance is required.
"""

from __future__ import annotations

import numpy as np

from fvcom_mesh_tools.io import Fort14Mesh


def edge_lengths_planar(mesh: Fort14Mesh) -> np.ndarray:
    """Per-element ``(NE, 3)`` array of edge lengths in node-coordinate units.

    Column order: ``[l_01, l_12, l_20]`` matching ``elements`` columns.
    """
    p0 = mesh.nodes[mesh.elements[:, 0]]
    p1 = mesh.nodes[mesh.elements[:, 1]]
    p2 = mesh.nodes[mesh.elements[:, 2]]
    return np.column_stack([
        np.linalg.norm(p1 - p0, axis=1),
        np.linalg.norm(p2 - p1, axis=1),
        np.linalg.norm(p0 - p2, axis=1),
    ])


def alpha_quality(mesh: Fort14Mesh) -> np.ndarray:
    """Per-triangle alpha-quality ``4*sqrt(3)*A / (l01^2 + l12^2 + l20^2)``.

    1 = equilateral, 0 = degenerate. Scale-free so it works on lon/lat
    meshes without unit conversion.
    """
    p0 = mesh.nodes[mesh.elements[:, 0]]
    p1 = mesh.nodes[mesh.elements[:, 1]]
    p2 = mesh.nodes[mesh.elements[:, 2]]
    twice_signed = (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )
    area = 0.5 * np.abs(twice_signed)
    ll = edge_lengths_planar(mesh)
    denom = (ll ** 2).sum(axis=1)
    return 4.0 * np.sqrt(3.0) * area / np.where(denom == 0, 1.0, denom)


def min_interior_angle(mesh: Fort14Mesh, in_degrees: bool = True) -> np.ndarray:
    """Per-triangle smallest interior angle.

    With ``in_degrees=True`` (default) returns degrees; otherwise radians.
    Scale-free.
    """
    ll = edge_lengths_planar(mesh)
    a = ll[:, 1]  # opposite of vertex 0
    b = ll[:, 2]  # opposite of vertex 1
    c = ll[:, 0]  # opposite of vertex 2

    def _angle(opp: np.ndarray, e1: np.ndarray, e2: np.ndarray) -> np.ndarray:
        cos = (e1 ** 2 + e2 ** 2 - opp ** 2) / np.where(
            e1 * e2 == 0, 1.0, 2.0 * e1 * e2
        )
        return np.arccos(np.clip(cos, -1.0, 1.0))

    A = _angle(a, b, c)
    B = _angle(b, c, a)
    C = _angle(c, a, b)
    rad = np.minimum(np.minimum(A, B), C)
    return np.degrees(rad) if in_degrees else rad
