"""PoC #56b: post-process PoC #54c output with keep_components(min=5).

PoC #56a's per-fail-element diagnostic on PoC #54c residuals
revealed that 4 of the 8 C1 fails are "isolated spit triangles" —
triangles where all 3 vertices have valence 1 and all 3 edges are
on the boundary. These are 1-element face-face-adjacency
components disconnected from the rest of the mesh, leaked through
``clean_mesh``'s default ``min_component_elements=0``.

A quick login-node test of ``keep_components`` at thresholds
1 / 2 / 5 / 10 showed the residual drops as:

  before                : C1=8  C2=0  C4=69  C5=0  (NE=87,022)
  min_elements=1        : C1=8  C2=0  C4=69  C5=0   no-op
  min_elements=2        : C1=4  C2=0  C4=69  C5=0   -4 C1
  min_elements=5        : C1=2  C2=0  C4=68  C5=0   -6 C1
  min_elements=10       : C1=2  C2=0  C4=66  C5=0   -6 C1, -3 C4

``min_elements=5`` removes 57 elements (~0.07 % of the mesh),
clearing the four fully-isolated spit triangles plus a couple of
tiny boundary islands — a topologically-cleaner mesh that arguably
shouldn't have these artifacts in the first place. This PoC
applies the post-process and records the result as the new
"reference final" for the goal of FVCOM-zero residual.

Since the isolated triangles share NO edges with the main mesh,
deleting them does not affect any other element's quality. Phase
H does not need to re-run.

Outputs:
    outputs/56b_phase_h_min5_optimized.14
    outputs/56b_summary.{txt,json}
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from fvcom_mesh_tools.algorithms.quality import min_interior_angle
from fvcom_mesh_tools.diagnostics import node_valence
from fvcom_mesh_tools.io import read_fort14, write_fort14
from fvcom_mesh_tools.mesh_clean import keep_components, rebuild_boundaries
from fvcom_mesh_tools.mesh_clean_phase_h import _per_edge_area_change

REPO = Path(__file__).resolve().parent.parent
INPUT = REPO / "outputs" / "54c_phase_h_optimized.14"
OUT_DIR = REPO / "outputs"
OUTPUT = OUT_DIR / "56b_phase_h_min5_optimized.14"
SUMMARY_TXT = OUT_DIR / "56b_summary.txt"
SUMMARY_JSON = OUT_DIR / "56b_summary.json"

MIN_COMPONENT_ELEMENTS = 5

# rebuild_boundaries needs the same bbox / tolerance the original
# pipeline used.
BBOX = (139.46, 34.99, 140.10, 35.74)
BBOX_TOL_M = 150
LAND_IBTYPE = 20
OPEN_MERGE_COAST_GAP = 50


def _max_interior_angle(mesh) -> np.ndarray:
    p0 = mesh.nodes[mesh.elements[:, 0]]
    p1 = mesh.nodes[mesh.elements[:, 1]]
    p2 = mesh.nodes[mesh.elements[:, 2]]
    e0 = np.linalg.norm(p1 - p0, axis=1)
    e1 = np.linalg.norm(p2 - p1, axis=1)
    e2 = np.linalg.norm(p0 - p2, axis=1)

    def _ang(opp, ea, eb):
        denom = 2.0 * ea * eb
        safe = np.where(denom > 0, denom, 1.0)
        cos = np.where(
            denom > 0, (ea * ea + eb * eb - opp * opp) / safe, 0.0,
        )
        return np.arccos(np.clip(cos, -1.0, 1.0))

    return np.degrees(
        np.maximum(
            np.maximum(_ang(e1, e2, e0), _ang(e2, e0, e1)),
            _ang(e0, e1, e2),
        ),
    )


def _metrics(mesh):
    m = min_interior_angle(mesh)
    M = _max_interior_angle(mesh)
    _u, _p, ac = _per_edge_area_change(mesh.nodes, mesh.elements)
    val = node_valence(mesh.elements, mesh.n_nodes)
    return {
        "NP": int(mesh.n_nodes),
        "NE": int(mesh.n_elements),
        "C1": int((m < 30.0).sum()),
        "C2": int((M > 130.0).sum()),
        "C4": int((ac > 0.5).sum()),
        "C5": int((val > 8).sum()),
        "max_valence": int(val.max()),
    }


def _deg_per_metre_lat(lat_deg: float) -> float:
    return 1.0 / 111_000.0


def main() -> int:
    if not INPUT.exists():
        raise SystemExit(f"input missing: {INPUT}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    mesh = read_fort14(INPUT)
    before = _metrics(mesh)
    print(
        f"[56b] before: NP={before['NP']:,} NE={before['NE']:,} "
        f"C1={before['C1']} C2={before['C2']} "
        f"C4={before['C4']} C5={before['C5']}",
        flush=True,
    )

    t0 = time.perf_counter()
    cleaned, info_keep = keep_components(
        mesh,
        min_elements=MIN_COMPONENT_ELEMENTS,
        require_open_boundary=False,
    )
    # Rebuild boundaries after element removal so the boundary
    # segment lists stay consistent with the surviving topology.
    tol_deg = (
        BBOX_TOL_M * _deg_per_metre_lat(float(cleaned.nodes[:, 1].mean()))
        if cleaned.n_nodes else 0.0
    )
    cleaned = rebuild_boundaries(
        cleaned, bbox=BBOX, tol_deg=tol_deg,
        land_ibtype=LAND_IBTYPE,
        open_merge_coast_gap=OPEN_MERGE_COAST_GAP,
    )
    wall = time.perf_counter() - t0

    after = _metrics(cleaned)
    write_fort14(cleaned, OUTPUT)
    print(
        f"[56b] after  (min_elements={MIN_COMPONENT_ELEMENTS}, "
        f"wall {wall:.2f} s): "
        f"NP={after['NP']:,} NE={after['NE']:,} "
        f"C1={after['C1']} C2={after['C2']} "
        f"C4={after['C4']} C5={after['C5']}",
        flush=True,
    )

    payload = {
        "input": str(INPUT.resolve()),
        "output": str(OUTPUT.resolve()),
        "min_component_elements": MIN_COMPONENT_ELEMENTS,
        "before": before,
        "after": after,
        "delta": {k: after[k] - before[k] for k in before},
        "wall_seconds": wall,
        "keep_components_info": info_keep,
    }
    SUMMARY_JSON.write_text(
        json.dumps(payload, indent=2, default=str), encoding="utf-8",
    )

    total_before = before["C1"] + before["C2"] + before["C4"] + before["C5"]
    total_after = after["C1"] + after["C2"] + after["C4"] + after["C5"]
    lines = [
        f"PoC #56b — keep_components(min_elements={MIN_COMPONENT_ELEMENTS}) "
        f"post-process",
        f"input : {INPUT.name}",
        f"output: {OUTPUT.name}",
        "",
        f"  {'stage':<22} | {'NP':>6} | {'NE':>6} | "
        f"{'C1':>4} | {'C2':>4} | {'C4':>4} | {'C5':>3}",
        "  " + "-" * 76,
        f"  {'PoC #54c (raw)':<22} | "
        f"{before['NP']:>6,} | {before['NE']:>6,} | "
        f"{before['C1']:>4} | {before['C2']:>4} | "
        f"{before['C4']:>4} | {before['C5']:>3}",
        f"  {'PoC #56b (min=5)':<22} | "
        f"{after['NP']:>6,} | {after['NE']:>6,} | "
        f"{after['C1']:>4} | {after['C2']:>4} | "
        f"{after['C4']:>4} | {after['C5']:>3}",
        "",
        f"  total violations: {total_before} -> {total_after} "
        f"({total_after - total_before:+})",
        f"  removed elements: {before['NE'] - after['NE']}",
        f"  removed nodes   : {before['NP'] - after['NP']}",
    ]
    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print()
    print("\n".join(lines))
    print()
    print(f"wrote {OUTPUT}")
    print(f"wrote {SUMMARY_TXT}")
    print(f"wrote {SUMMARY_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
