"""PoC #59d — relocate the open boundary to the Uraga Channel transect.

Domain forensics on the #59c mesh showed the lineage's open boundary
sits on the EASTERN 140.10E data-clip line near the Chiba coast
(shallow, zigzag — the source of every remaining perpendicularity
failure and of the #59c perpfix/finish ping-pong), while the true
southern cut — a 2.8 km, 50-node transect across the Uraga Channel at
35.10N with depths to 338 m — was classified LAND (the 54a-era bbox
south edge 34.99N misses the DEM edge 35.10N by ~11 km, far beyond
the 150 m tolerance).

This PoC:

1. re-derives the outer ring and takes the longest contiguous run of
   ring nodes with ``lat <= lat_min + 0.002 deg`` as THE open
   boundary (the Uraga transect — deep, smooth water per kickoff
   §7.2);
2. reclassifies the eastern 140.10E cut as a land wall (**domain
   design decision to confirm with the user**: alternative is a
   secondary clamped OBC);
3. runs the #59c convergence loop (quality fails -> phase_h_finish +
   compact; perp fails -> damped perpfix), plus end-trimming of the
   open segment if junction-element checks (R4 / fake-ISBCE=2) fail.

Output: outputs/59d_gate_passed.14 (+ _qa.json, cycle log).
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
from pyproj import Transformer

from fvcom_mesh_tools.algorithms import align_open_boundary_first_ring
from fvcom_mesh_tools.algorithms.boundary import (
    boundary_edges_from_tris,
    chain_edges_to_loops,
    outer_loop,
)
from fvcom_mesh_tools.io import Fort14Mesh, read_fort14, write_fort14
from fvcom_mesh_tools.mesh_clean import compact_nodes
from fvcom_mesh_tools.mesh_clean_phase_h import phase_h_finish
from fvcom_mesh_tools.qa import format_report, run_qa

REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "outputs" / "59c_gate_passed.14"
OUT = REPO / "outputs" / "59d_gate_passed.14"
QA_JSON = REPO / "outputs" / "59d_gate_passed_qa.json"
CYCLES_JSON = REPO / "outputs" / "59d_cycle_log.json"

LAT_TOL_DEG = 0.002       # ~220 m band above the domain's southern edge
LAND_IBTYPE = 20
MAX_CYCLES = 10
MAX_END_TRIMS = 4
QUALITY_IDS = {"c1_min_angle", "c2_max_angle", "c4_area_change", "c5_valence"}
JUNCTION_IDS = {"r4_mixed_boundary", "isbce2_authentic", "obc_interior_neighbor"}


def _south_run(mesh: Fort14Mesh) -> tuple[np.ndarray, list[np.ndarray], int]:
    """Longest contiguous outer-ring run within the southern lat band,
    plus the island loops (closing duplicate dropped)."""
    tr = Transformer.from_crs("EPSG:32654", "EPSG:4326", always_xy=True)
    _lon, lat = tr.transform(mesh.nodes[:, 0], mesh.nodes[:, 1])
    loops = chain_edges_to_loops(boundary_edges_from_tris(mesh.elements))
    outer = outer_loop(loops, mesh.nodes)
    ring = outer[:-1]
    rlat = lat[ring]
    mask = rlat <= float(rlat.min()) + LAT_TOL_DEG
    idx = np.where(mask)[0]
    if idx.size == 0:
        raise SystemExit("[59d] no southern-band ring nodes found")
    runs: list[tuple[int, int]] = []
    start = prev = int(idx[0])
    for p in idx[1:]:
        p = int(p)
        if p == prev + 1:
            prev = p
        else:
            runs.append((start, prev))
            start = prev = p
    runs.append((start, prev))
    a, b = max(runs, key=lambda r: r[1] - r[0])
    # Rotate the ring so the open run does not wrap the array end.
    ring = np.roll(ring, -a)
    b -= a
    islands = [lp[:-1].copy() for lp in loops if lp is not outer]
    return ring, islands, b  # open run = ring[0 : b + 1]


def _apply_obc(mesh: Fort14Mesh, ring: np.ndarray, islands, open_end: int,
               trim: int) -> Fort14Mesh:
    """Set open boundary = ring[trim : open_end + 1 - trim]; the rest of
    the outer ring plus the islands become land (ibtype 20). Open and
    land share their junction endpoints per fort.14 convention."""
    lo, hi = trim, open_end - trim
    if hi - lo < 2:
        raise SystemExit("[59d] open segment trimmed away entirely")
    open_seg = ring[lo : hi + 1].copy()
    land_seg = np.concatenate([ring[hi:], ring[: lo + 1]])
    land = [(LAND_IBTYPE, land_seg)]
    land += [(LAND_IBTYPE, isl.copy()) for isl in islands]
    return Fort14Mesh(
        title=mesh.title,
        nodes=mesh.nodes,
        depths=mesh.depths,
        elements=mesh.elements,
        open_boundaries=[open_seg],
        land_boundaries=land,
    )


def _perpfix_all(mesh, *, alpha: float, n_iters: int):
    for k in range(len(mesh.open_boundaries)):
        mesh, _ = align_open_boundary_first_ring(
            mesh, alpha=alpha, n_iters=n_iters,
            smooth_iters=2, smooth_alpha=0.3, segment_index=k,
        )
    return mesh


def _failed_gate_ids(report) -> set[str]:
    return {
        c.check_id for c in report.checks
        if c.gate and not c.skipped and not c.passed
    }


def main() -> int:
    t0 = time.perf_counter()
    mesh = read_fort14(SRC)
    ring, islands, open_end = _south_run(mesh)
    trim = 0
    mesh = _apply_obc(mesh, ring, islands, open_end, trim)
    seg = mesh.open_boundaries[0]
    d = mesh.depths[seg]
    print(
        f"[59d] Uraga OBC: {len(seg)} nodes, depth "
        f"{d.min():.1f}-{d.max():.1f} m (mean {d.mean():.1f}); "
        f"eastern 140.10E cut reclassified as land wall",
        flush=True,
    )

    log: list[dict] = []
    report = run_qa(mesh, name=OUT.name, path=OUT)
    for cycle in range(MAX_CYCLES):
        failed = _failed_gate_ids(report)
        counts = {
            c.check_id: int(c.n_violations)
            for c in report.checks if c.gate and not c.skipped and c.n_violations
        }
        print(f"[59d] cycle {cycle}: failed = {sorted(failed)} {counts}", flush=True)
        log.append({"cycle": cycle, "failed_gates": sorted(failed), "counts": counts})
        if not failed:
            break
        if failed & JUNCTION_IDS and trim < MAX_END_TRIMS:
            trim += 1
            mesh = _apply_obc(mesh, ring, islands, open_end, trim)
            print(f"[59d] cycle {cycle}: junction issue — trimmed OBC ends "
                  f"(trim={trim}, {len(mesh.open_boundaries[0])} nodes)", flush=True)
        elif failed & QUALITY_IDS:
            mesh, finfo = phase_h_finish(mesh, seed=60 + cycle)
            mesh, cinfo = compact_nodes(mesh)
            print(
                f"[59d] cycle {cycle}: finish {finfo.get('before')} -> "
                f"{finfo.get('after')}; compacted {cinfo['n_orphans_removed']}",
                flush=True,
            )
            if cinfo["n_orphans_removed"]:
                # Node ids changed; re-derive the ring for later trims.
                ring, islands, open_end = _south_run(mesh)
                trim = 0
                mesh = _apply_obc(mesh, ring, islands, open_end, trim)
        elif "obc_perpendicularity" in failed:
            mesh = _perpfix_all(mesh, alpha=0.6, n_iters=2)
            print(f"[59d] cycle {cycle}: damped perpfix applied", flush=True)
        else:
            print(f"[59d] cycle {cycle}: unhandled residual — stopping", flush=True)
            break
        report = run_qa(mesh, name=OUT.name, path=OUT)

    mesh.title = "PoC 59d UTM54N Tokyo Bay, OBC at Uraga Channel transect"
    write_fort14(mesh, OUT)
    QA_JSON.write_text(
        json.dumps(report.to_dict(), indent=2, ensure_ascii=False), encoding="utf-8",
    )
    CYCLES_JSON.write_text(
        json.dumps(log, indent=2, ensure_ascii=False), encoding="utf-8",
    )
    print(format_report(report, lang="ja"), flush=True)
    print(f"[59d] wrote {OUT}", flush=True)
    print(f"[59d] wall: {time.perf_counter() - t0:.1f} s", flush=True)
    return 0 if report.passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
