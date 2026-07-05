"""Constraint-respecting finishing of CDT-boundary meshes.

Extraction of the PoC #93 stages into a reusable function: the mesh
arrives with its boundary ON the engineered shoreline (CDT egfix)
plus the sliver population the disabled boundary-deleters would have
removed by retreating. Finishing removes the slivers WITHOUT giving
up the boundary:

1. cap-triangle deletion — slivers whose three vertices all lie on
   the mesh boundary; deleting them exposes the constrained edges,
   so conformity is untouched;
2. weld of near-coincident nodes (chain corners);
3. ONE budgeted `phase_h_optimize` round with the shoreline
   projector (boundary nodes may slide ALONG the line only).

No convergence loops: each stage runs once (user discipline,
docs/DESIGN_HISTORY.md).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable

import numpy as np

from fvcom_mesh_tools.algorithms.perp_local import _tri_quality
from fvcom_mesh_tools.io import Fort14Mesh
from fvcom_mesh_tools.mesh_clean import (
    compact_nodes,
    keep_components,
    remove_elements,
    weld_close_nodes,
)
from fvcom_mesh_tools.mesh_clean_phase_h import (
    build_coastline_projector,
    phase_h_optimize,
)

__all__ = ["finish_constrained_mesh", "quality_counts"]

_DEG_PER_M = 1.0 / 111_194.92664455873


def quality_counts(mesh: Fort14Mesh) -> dict[str, Any]:
    """C1/C2/flip counts and the minimum interior angle."""
    mn, mx, tw = _tri_quality(mesh.nodes[mesh.elements])
    return {
        "c1": int((mn < 30).sum()),
        "c2": int((mx > 130).sum()),
        "flipped": int((tw <= 0).sum()),
        "min_angle_deg": float(mn.min()) if mn.size else 90.0,
    }


def _boundary_node_set(mesh: Fort14Mesh) -> set[int]:
    els = mesh.elements
    raw = np.vstack([els[:, [0, 1]], els[:, [1, 2]], els[:, [2, 0]]])
    raw.sort(axis=1)
    codes = raw[:, 0].astype(np.int64) * mesh.n_nodes + raw[:, 1]
    uniq, counts = np.unique(codes, return_counts=True)
    b = uniq[counts == 1]
    return set(np.column_stack([b // mesh.n_nodes,
                                b % mesh.n_nodes]).ravel().tolist())


def finish_constrained_mesh(
    mesh: Fort14Mesh,
    shoreline_shp: Path,
    work_dir: Path,
    *,
    min_angle: float = 30.0,
    max_angle: float = 130.0,
    cap_min_angle: float = 20.0,
    cap_max_angle: float = 150.0,
    cap_rounds: int = 3,
    weld_tol_m: float = 2.0,
    optimize_budget_s: float = 600.0,
    projector_snap_m: float = 200.0,
    utm_epsg: int = 32654,
    log: Callable[[str], None] = print,
) -> tuple[Fort14Mesh, dict[str, Any]]:
    """Finish a CDT-boundary mesh (UTM coordinates) as in PoC #93.

    ``shoreline_shp`` is the engineered shoreline (EPSG:4326
    polygons) the boundary is constrained to; a no-CRS UTM copy for
    the projector is cached under ``work_dir``.
    """
    import geopandas as gpd

    info: dict[str, Any] = {"input": quality_counts(mesh)}
    log(f"[finish] input NP={mesh.n_nodes:,} quality={info['input']}")

    # 1. cap-triangle deletion (bounded rounds, monotone).
    removed = 0
    for _ in range(cap_rounds):
        bset = _boundary_node_set(mesh)
        mn, mx, _tw = _tri_quality(mesh.nodes[mesh.elements])
        bad = (mn < cap_min_angle) | (mx > cap_max_angle)
        allb = np.isin(mesh.elements, list(bset)).all(axis=1)
        kill = bad & allb
        if not kill.any():
            break
        mesh = remove_elements(mesh, ~kill)
        mesh, _ = keep_components(mesh)
        removed += int(kill.sum())
    mesh, _ = compact_nodes(mesh)
    info["caps_removed"] = removed
    log(f"[finish] caps removed: {removed} quality={quality_counts(mesh)}")

    # 2. weld.
    mesh, winfo = weld_close_nodes(mesh, tol=weld_tol_m)
    info["weld"] = winfo
    log(f"[finish] weld: {winfo}")

    # 3. one budgeted optimize round, boundary sliding on-line only.
    work_dir.mkdir(parents=True, exist_ok=True)
    shp_utm = work_dir / "finish_shoreline_utm_nocrs.shp"
    if not shp_utm.exists():
        gdf = gpd.read_file(shoreline_shp).to_crs(utm_epsg)
        gdf = gdf.set_crs(None, allow_override=True)
        gdf.to_file(shp_utm)
    projector = build_coastline_projector(
        [shp_utm],
        max_snap_distance_m=projector_snap_m / _DEG_PER_M,
        mean_latitude_deg=0.0,
    )
    mesh, oinfo = phase_h_optimize(
        mesh,
        min_angle_target=min_angle,
        max_angle_target=max_angle,
        pass_f_enabled=True,
        pass_g_enabled=True,
        pass_g_min_angle_target=min_angle,
        max_outer_rounds=1,
        coastline_projector=projector,
        time_budget_s=optimize_budget_s,
    )
    mesh, cinfo = compact_nodes(mesh)
    info["optimize"] = {
        "budget_exhausted": bool(oinfo.get("budget_exhausted")),
        "compact": cinfo,
    }
    info["output"] = quality_counts(mesh)
    log(f"[finish] done NP={mesh.n_nodes:,} quality={info['output']}")
    return mesh, info
