"""PoC #59k-2 stage A — structural finalize + projection, checkpointed.

The single-job #59k hit the elapse limit inside ``phase_h_finish``:
the finish stage is designed for a ~100-violation residual, and the
v3 raw mesh was fed to it without the deterministic optimize passes
(a chain-order mistake — the v1 lineage had Phase A-G applied long
before 58l). This stage reruns the fast part, SAVES the checkpoint,
and measures the metric-space C1/C2/C4/C5 counts so the next stage
can be sized with data instead of guesses.
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
from pyproj import Transformer

from fvcom_mesh_tools.diagnostics import node_valence
from fvcom_mesh_tools.io import read_fort14, write_fort14
from fvcom_mesh_tools.mesh_clean import keep_components, rebuild_boundaries, remove_elements
from fvcom_mesh_tools.mesh_clean_phase_h import (
    _per_edge_area_change,
    _per_element_quality,
)
from fvcom_mesh_tools.qa import fvcom_boundary_element_flags

REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "outputs" / "59j_v3_raw_100m.14"
OUT = REPO / "outputs" / "59k_stageA_utm.14"

BBOX = (139.565, 35.10, 140.172, 35.86)
TOL_DEG = 150.0 / 111_195.0
MIN_DEPTH_M = 2.0


def main() -> int:
    t0 = time.perf_counter()
    mesh = read_fort14(SRC)
    mesh, kinfo = keep_components(mesh)
    mesh = rebuild_boundaries(mesh, bbox=BBOX, tol_deg=TOL_DEG,
                              land_ibtype=20, open_merge_coast_gap=50)
    for rnd in range(30):
        flags = fvcom_boundary_element_flags(mesh)
        bad = flags["r4_mask"] | flags["fake_open_mask"]
        if not bad.any():
            break
        mesh = remove_elements(mesh, ~bad)
        mesh, _ = keep_components(mesh)
        mesh = rebuild_boundaries(mesh, bbox=BBOX, tol_deg=TOL_DEG,
                                  land_ibtype=20, open_merge_coast_gap=50)
    mesh.depths[:] = np.maximum(mesh.depths, MIN_DEPTH_M)
    tr = Transformer.from_crs("EPSG:4326", "EPSG:32654", always_xy=True)
    x, y = tr.transform(mesh.nodes[:, 0], mesh.nodes[:, 1])
    mesh.nodes = np.column_stack([x, y])
    mesh.title = "59k stage A: v3 structural + UTM"
    write_fort14(mesh, OUT)
    print(f"[59k2A] wrote {OUT}  NP={mesh.n_nodes:,} NE={mesh.n_elements:,}",
          flush=True)

    _alpha, min_ang, max_ang = _per_element_quality(mesh.nodes, mesh.elements)
    _uv, _pair, ac = _per_edge_area_change(mesh.nodes, mesh.elements)
    val = node_valence(mesh.elements, mesh.n_nodes)
    print(
        f"[59k2A] metric-space violations: "
        f"C1={int((min_ang < 30.0).sum())} "
        f"C2={int((max_ang > 130.0).sum())} "
        f"C4={int((ac > 0.5).sum())} "
        f"C5={int((val > 8).sum())} (max valence {int(val.max())})",
        flush=True,
    )
    print(f"[59k2A] wall: {time.perf_counter() - t0:.1f} s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
