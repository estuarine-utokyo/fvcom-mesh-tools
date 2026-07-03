"""PoC #59a — structural finalize of the Phase-H quality endpoint.

The first ``fmesh-mesh-qa`` run (2026-07-04) failed
``outputs/58l_chained.14`` on 10 gates. This PoC clears the
*structural* subset in lon/lat space:

* disjoint components (44 comps / 4,878 elems)  -> keep largest
* orphan nodes (2)                              -> node compaction
  (inside ``remove_elements``)
* R4 / fake-ISBCE=2 open-boundary elements      -> iterative delete +
  bbox boundary rebuild until both masks are empty
* depth < 2 m (9,688 nodes, min 0.00 m)         -> clip to 2.0 m
  (kickoff §5)

The metric-space quality residual (~55 C1) is left for PoC #59b
(UTM 54N projection + perpfix + metric-space ``phase_h_finish``).

Boundary rebuild mirrors the 54a lineage parameters: bbox
(139.46, 34.99, 140.10, 35.74), tol 150 m, land ibtype 20,
open_merge_coast_gap 50.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from fvcom_mesh_tools.io import read_fort14, write_fort14
from fvcom_mesh_tools.mesh_clean import (
    keep_components,
    rebuild_boundaries,
    remove_elements,
)
from fvcom_mesh_tools.qa import (
    format_report,
    fvcom_boundary_element_flags,
    run_qa,
)

REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "outputs" / "58l_chained.14"
OUT = REPO / "outputs" / "59a_structural_finalized.14"
QA_JSON = REPO / "outputs" / "59a_structural_finalized_qa.json"

BBOX = (139.46, 34.99, 140.10, 35.74)
TOL_DEG = 150.0 / 111_195.0
LAND_IBTYPE = 20
OPEN_MERGE_COAST_GAP = 50
MIN_DEPTH_M = 2.0
MAX_ROUNDS = 30


def _rebuild(mesh):
    return rebuild_boundaries(
        mesh,
        bbox=BBOX,
        tol_deg=TOL_DEG,
        land_ibtype=LAND_IBTYPE,
        open_merge_coast_gap=OPEN_MERGE_COAST_GAP,
    )


def main() -> int:
    t0 = time.perf_counter()
    mesh = read_fort14(SRC)
    print(
        f"[59a] input:  NP={mesh.n_nodes:,} NE={mesh.n_elements:,} "
        f"open={[len(s) for s in mesh.open_boundaries]} "
        f"land={len(mesh.land_boundaries)}",
        flush=True,
    )

    mesh, kinfo = keep_components(mesh)
    print(
        f"[59a] keep_components: {kinfo['n_components_before']} comps -> "
        f"largest ({kinfo['n_elements_removed']:,} elems, "
        f"{kinfo['n_nodes_removed']:,} nodes removed)",
        flush=True,
    )
    mesh = _rebuild(mesh)

    removed_total = 0
    for rnd in range(1, MAX_ROUNDS + 1):
        flags = fvcom_boundary_element_flags(mesh)
        n_r4 = int(flags["r4_mask"].sum())
        n_fake = int(flags["fake_open_mask"].sum())
        print(f"[59a] round {rnd}: R4={n_r4} fake_open={n_fake}", flush=True)
        bad = flags["r4_mask"] | flags["fake_open_mask"]
        if not bad.any():
            break
        removed_total += int(bad.sum())
        mesh = remove_elements(mesh, ~bad)
        # Deletions can split off small pools; keep the largest again
        # before re-deriving the boundary.
        mesh, kinfo = keep_components(mesh)
        removed_total += int(kinfo["n_elements_removed"])
        mesh = _rebuild(mesh)
    else:
        raise SystemExit("[59a] R4/fake-open deletion did not converge")

    n_shallow = int((mesh.depths < MIN_DEPTH_M).sum())
    mesh.depths[:] = np.maximum(mesh.depths, MIN_DEPTH_M)
    print(
        f"[59a] depth clip: {n_shallow:,} nodes raised to {MIN_DEPTH_M} m",
        flush=True,
    )

    mesh.title = "PoC 59a structural finalize of 58l_chained"
    write_fort14(mesh, OUT)
    print(
        f"[59a] wrote {OUT}  NP={mesh.n_nodes:,} NE={mesh.n_elements:,} "
        f"open={[len(s) for s in mesh.open_boundaries]} "
        f"land={len(mesh.land_boundaries)} "
        f"(elements removed total: {removed_total:,})",
        flush=True,
    )

    report = run_qa(mesh, name=OUT.name, path=OUT)
    print(format_report(report, lang="ja"), flush=True)
    QA_JSON.write_text(
        json.dumps(report.to_dict(), indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print(f"[59a] wall: {time.perf_counter() - t0:.1f} s", flush=True)
    # Structural gates must pass here; metric-space quality (C1 etc.)
    # is PoC #59b's job, so do not fail the job on them.
    structural = {
        "ccw_all_elements", "no_isolated_elements", "r4_mixed_boundary",
        "manifold_boundary", "no_duplicate_nodes", "no_orphan_nodes",
        "no_tiny_area", "isbce2_authentic", "obc_on_boundary",
        "obc_chain_adjacency", "obc_interior_neighbor", "obc_ordering",
        "single_component", "obc_reachable", "min_depth_clip",
    }
    bad_structural = [
        c.check_id for c in report.checks
        if c.check_id in structural and c.gate and not c.skipped and not c.passed
    ]
    if bad_structural:
        print(f"[59a] STRUCTURAL GATES STILL FAILING: {bad_structural}", flush=True)
        return 1
    print("[59a] all structural gates pass", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
