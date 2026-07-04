"""PoC #93 — constraint-respecting cleanup of the CDT mesh.

#92 v3: boundary p50 7 mm ON the engineered shoreline (CDT egfix),
but the slivers the old boundary-deleters used to remove (by
retreating) are still there: C1=1,725 / C2=526. This pass removes
them WITHOUT giving up the boundary: cap-triangle deletion (all
three nodes on the boundary chain — deleting exposes the two
constrained edges, conformity intact), interior-apex vertex removal
for needle slivers, weld, then ONE budgeted optimize round with the
engineered-shoreline projector (boundary nodes slide along the
line only).
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
from pyproj import Transformer

from fvcom_mesh_tools.algorithms.perp_local import _tri_quality
from fvcom_mesh_tools.io import read_fort14, write_fort14
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
from fvcom_mesh_tools.plotting import plot_mesh_overview

REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "outputs" / "92_cdt_hmin300.14"
OUT = REPO / "outputs" / "93_cdt_cleaned.14"
SHORE = REPO / "outputs" / "osm_shoreline" / "osm_land_opened_land150.shp"
UTM_DIR = REPO / "outputs" / "osm_shoreline_utm_nocrs"
_TO_M = Transformer.from_crs("EPSG:4326", "EPSG:32654", always_xy=True)
_TO_LL = Transformer.from_crs("EPSG:32654", "EPSG:4326", always_xy=True)
_DEG_PER_M = 1.0 / 111_194.92664455873


def _fig(mesh_utm, tag):
    import copy

    m2 = copy.deepcopy(mesh_utm)
    lon, lat = _TO_LL.transform(m2.nodes[:, 0], m2.nodes[:, 1])
    m2.nodes = np.column_stack([lon, lat])
    for view, zb in {
        "tokyoport": (139.72, 35.53, 139.90, 35.68),
        "full": None,
    }.items():
        plot_mesh_overview(
            m2, REPO / "outputs" / "figures" / f"93_{tag}_{view}.png",
            crs="EPSG:4326", cell_m=None,
            coast=(139.4, 34.9, 140.4, 35.85), zoom=zb, dpi=220,
            title=f"93 {tag} - {view}",
        )


def _quality(mesh):
    mn, mx, tw = _tri_quality(mesh.nodes[mesh.elements])
    return (int((mn < 30).sum()), int((mx > 130).sum()),
            int((tw <= 0).sum()), float(mn.min()))


def main() -> int:
    import logging

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        datefmt="%H:%M:%S")
    t0 = time.perf_counter()
    mesh = read_fort14(SRC)
    x, y = _TO_M.transform(mesh.nodes[:, 0], mesh.nodes[:, 1])
    mesh.nodes = np.column_stack([x, y])
    print(f"[93] input NP={mesh.n_nodes:,} NE={mesh.n_elements:,} "
          f"quality C1/C2/flip/min={_quality(mesh)}", flush=True)
    _fig(mesh, "s0_input")

    # boundary node set
    def _bnodes(m):
        els = m.elements
        raw = np.vstack([els[:, [0, 1]], els[:, [1, 2]], els[:, [2, 0]]])
        raw.sort(axis=1)
        codes = raw[:, 0].astype(np.int64) * m.n_nodes + raw[:, 1]
        uniq, cnt = np.unique(codes, return_counts=True)
        bb = uniq[cnt == 1]
        return set(np.column_stack([bb // m.n_nodes,
                                    bb % m.n_nodes]).ravel().tolist())

    # S1: cap-triangle deletion — sliver whose 3 nodes are all
    # boundary nodes: delete (exposes the constrained edges).
    for rnd in range(3):
        bset = _bnodes(mesh)
        mn, mx, tw = _tri_quality(mesh.nodes[mesh.elements])
        bad = (mn < 20.0) | (mx > 150.0)
        allb = np.isin(mesh.elements, list(bset)).all(axis=1)
        kill = bad & allb
        if not kill.any():
            break
        mesh = remove_elements(mesh, ~kill)
        mesh, _ = keep_components(mesh)
        print(f"[93] S1 r{rnd}: removed {int(kill.sum())} cap slivers",
              flush=True)
    mesh, _ = compact_nodes(mesh)
    print(f"[93] S1 done: NP={mesh.n_nodes:,} "
          f"quality={_quality(mesh)}", flush=True)
    _fig(mesh, "s1_caps")

    # S2: weld near-coincident nodes (chain corners).
    mesh, winfo = weld_close_nodes(mesh, tol=2.0)
    print(f"[93] S2 weld: {winfo}", flush=True)

    # S3: ONE budgeted optimize with the engineered-shoreline
    # projector (boundary nodes slide ALONG the line only).
    import geopandas as gpd

    UTM_DIR.mkdir(parents=True, exist_ok=True)
    shp = UTM_DIR / "land_opened_utm54.shp"
    if not shp.exists():
        gdf = gpd.read_file(SHORE).to_crs(32654)
        gdf = gdf.set_crs(None, allow_override=True)
        gdf.to_file(shp)
    projector = build_coastline_projector(
        [shp], max_snap_distance_m=200.0 / _DEG_PER_M,
        mean_latitude_deg=0.0,
    )
    mesh, oinfo = phase_h_optimize(
        mesh,
        min_angle_target=30.0,
        max_angle_target=130.0,
        pass_f_enabled=True,
        pass_g_enabled=True,
        pass_g_min_angle_target=30.0,
        max_outer_rounds=1,
        coastline_projector=projector,
        time_budget_s=600.0,
    )
    mesh, cinfo = compact_nodes(mesh)
    print(f"[93] S3 optimize: NP={mesh.n_nodes:,} "
          f"quality={_quality(mesh)} compact={cinfo}", flush=True)
    _fig(mesh, "s3_opt")

    write_fort14(mesh, OUT)
    print(f"[93] wrote {OUT}  wall={time.perf_counter() - t0:.0f} s",
          flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
