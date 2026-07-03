"""PoC #59h — open the FULL west+south artificial arc (tokyo_bay_v2).

The v1 OBC covered only the deep 2.8 km Uraga run of the southern
cut; the rest of the southern 35.10N line and the entire western
139.565E line (Sagami side, depths to 626 m) were left as land
walls — physically wrong for tidal forcing (user 2026-07-04). This
PoC reclassifies the whole contiguous artificial arc (~308 nodes) as
ONE open segment, re-converges every QA gate (junction trims,
R4/fake-open deletions, metric finish, the toolkit-promoted local
perp fixer), and exports the ``tokyo_bay_v2`` FVCOM input set for
the M2 tidal test (#59i).

The eastern 140.10E clip stays a land wall in v2; it disappears
entirely in the v3 east-extended rebuild (#59j).
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
from pyproj import Transformer

from fvcom_mesh_tools.algorithms import align_open_boundary_local
from fvcom_mesh_tools.algorithms.boundary import (
    boundary_edges_from_tris,
    chain_edges_to_loops,
    outer_loop,
)
from fvcom_mesh_tools.io import Fort14Mesh, read_fort14, write_fort14
from fvcom_mesh_tools.io.fvcom_native import export_fvcom_case
from fvcom_mesh_tools.mesh_clean import compact_nodes, keep_components, remove_elements
from fvcom_mesh_tools.mesh_clean_phase_h import phase_h_finish
from fvcom_mesh_tools.qa import format_report, fvcom_boundary_element_flags, run_qa

REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "outputs" / "59e_gate_passed.14"
OUT = REPO / "outputs" / "59h_gate_passed.14"
QA_JSON = REPO / "outputs" / "59h_gate_passed_qa.json"
LOG_JSON = REPO / "outputs" / "59h_cycle_log.json"
EXPORT_DIR = REPO / "outputs" / "fvcom_inputs_v2"
CASENAME = "tokyo_bay_v2"

BAND_DEG = 0.008          # west/south artificial-cut band (~890 m)
LAND_IBTYPE = 20
MAX_CYCLES = 14
MAX_END_TRIMS = 4
QUALITY_IDS = {"c1_min_angle", "c2_max_angle", "c4_area_change", "c5_valence"}
JUNCTION_IDS = {"r4_mixed_boundary", "isbce2_authentic", "obc_interior_neighbor"}
FATAL_OR_QUALITY = {
    "node_index_valid", "ccw_all_elements", "no_isolated_elements",
    "r4_mixed_boundary", "manifold_boundary", "no_duplicate_nodes",
    "no_orphan_nodes", "no_tiny_area", "isbce2_authentic",
    "obc_on_boundary", "obc_chain_adjacency", "obc_interior_neighbor",
    "obc_ordering", "c1_min_angle", "c2_max_angle", "c4_area_change",
    "c5_valence", "single_component", "obc_reachable", "min_depth_clip",
}

_TR = Transformer.from_crs("EPSG:32654", "EPSG:4326", always_xy=True)


def _arc(mesh: Fort14Mesh):
    """Longest contiguous outer-ring run inside the west/south band."""
    lon, lat = _TR.transform(mesh.nodes[:, 0], mesh.nodes[:, 1])
    loops = chain_edges_to_loops(boundary_edges_from_tris(mesh.elements))
    outer = outer_loop(loops, mesh.nodes)
    ring = outer[:-1]
    rlon, rlat = lon[ring], lat[ring]
    mask = (rlat <= rlat.min() + BAND_DEG) | (rlon <= rlon.min() + BAND_DEG)
    idx = np.where(mask)[0]
    if idx.size == 0:
        raise SystemExit("[59h] no artificial-arc nodes found")
    runs = []
    s = p = int(idx[0])
    for q in idx[1:]:
        q = int(q)
        if q == p + 1:
            p = q
        else:
            runs.append((s, p))
            s = p = q
    runs.append((s, p))
    a, b = max(runs, key=lambda r: r[1] - r[0])
    ring = np.roll(ring, -a)
    islands = [lp[:-1].copy() for lp in loops if lp is not outer]
    return ring, islands, b - a


def _apply_obc(mesh, ring, islands, open_end, trim):
    lo, hi = trim, open_end - trim
    if hi - lo < 2:
        raise SystemExit("[59h] open arc trimmed away")
    open_seg = ring[lo : hi + 1].copy()
    land_seg = np.concatenate([ring[hi:], ring[: lo + 1]])
    land = [(LAND_IBTYPE, land_seg)] + [(LAND_IBTYPE, i.copy()) for i in islands]
    return Fort14Mesh(
        title=mesh.title, nodes=mesh.nodes, depths=mesh.depths,
        elements=mesh.elements, open_boundaries=[open_seg],
        land_boundaries=land,
    )


def _failed(report):
    return {c.check_id for c in report.checks
            if c.gate and not c.skipped and not c.passed}


def main() -> int:
    t0 = time.perf_counter()
    mesh = read_fort14(SRC)
    ring, islands, open_end = _arc(mesh)
    trim = 0
    mesh = _apply_obc(mesh, ring, islands, open_end, trim)
    seg = mesh.open_boundaries[0]
    d = mesh.depths[seg]
    print(
        f"[59h] arc OBC: {len(seg)} nodes, depth {d.min():.1f}-{d.max():.1f} m "
        f"(mean {d.mean():.1f})",
        flush=True,
    )

    log = []
    report = run_qa(mesh, name=OUT.name, path=OUT)
    for cycle in range(MAX_CYCLES):
        failed = _failed(report)
        counts = {c.check_id: int(c.n_violations) for c in report.checks
                  if c.gate and not c.skipped and c.n_violations}
        print(f"[59h] cycle {cycle}: failed = {sorted(failed)} {counts}", flush=True)
        log.append({"cycle": cycle, "failed": sorted(failed), "counts": counts})
        if not failed:
            break
        if failed & JUNCTION_IDS:
            if trim < MAX_END_TRIMS:
                trim += 1
                mesh = _apply_obc(mesh, ring, islands, open_end, trim)
                print(f"[59h] cycle {cycle}: trimmed arc ends (trim={trim})",
                      flush=True)
            else:
                flags = fvcom_boundary_element_flags(mesh)
                bad = flags["r4_mask"] | flags["fake_open_mask"]
                if not bad.any():
                    print(f"[59h] cycle {cycle}: junction residual without "
                          "R4/fake elements — stopping", flush=True)
                    break
                mesh = remove_elements(mesh, ~bad)
                mesh, _ = keep_components(mesh)
                ring, islands, open_end = _arc(mesh)
                trim = 0
                mesh = _apply_obc(mesh, ring, islands, open_end, trim)
                print(f"[59h] cycle {cycle}: deleted {int(bad.sum())} "
                      "R4/fake elements, arc re-derived", flush=True)
        elif failed & QUALITY_IDS:
            mesh, finfo = phase_h_finish(mesh, seed=70 + cycle)
            mesh, cinfo = compact_nodes(mesh)
            print(f"[59h] cycle {cycle}: finish {finfo.get('before')} -> "
                  f"{finfo.get('after')}; compacted "
                  f"{cinfo['n_orphans_removed']}", flush=True)
            if cinfo["n_orphans_removed"]:
                ring, islands, open_end = _arc(mesh)
                trim = 0
                mesh = _apply_obc(mesh, ring, islands, open_end, trim)
        elif "obc_perpendicularity" in failed:
            mesh, pinfo = align_open_boundary_local(mesh, seed=4242 + cycle)
            print(f"[59h] cycle {cycle}: local perp fixer "
                  f"accepted={pinfo['accepted_total']} "
                  f"remaining={len(pinfo['remaining'])}", flush=True)
            if pinfo["accepted_total"] == 0:
                report = run_qa(mesh, name=OUT.name, path=OUT)
                log.append({"cycle": f"{cycle}-final",
                            "failed": sorted(_failed(report))})
                print("[59h] perp fixer exhausted — stopping", flush=True)
                break
        else:
            print(f"[59h] cycle {cycle}: unhandled residual — stopping", flush=True)
            break
        report = run_qa(mesh, name=OUT.name, path=OUT)

    mesh.title = "PoC 59h UTM54N Tokyo Bay, full west+south arc OBC"
    write_fort14(mesh, OUT)
    QA_JSON.write_text(json.dumps(report.to_dict(), indent=2, ensure_ascii=False),
                       encoding="utf-8")
    LOG_JSON.write_text(json.dumps(log, indent=2, ensure_ascii=False),
                        encoding="utf-8")
    print(format_report(report, lang="ja"), flush=True)

    failed = _failed(report)
    # Export for the tidal test as long as nothing startup-fatal or
    # quality-level remains (a small perp residue on the jagged cut is
    # advisory-grade and must not block the forcing-machinery test).
    if not (failed & FATAL_OR_QUALITY):
        lon, lat = _TR.transform(mesh.nodes[:, 0], mesh.nodes[:, 1])
        written = export_fvcom_case(
            mesh, EXPORT_DIR, CASENAME,
            obc_type=1, cor=lat, write_empty_spg=True,
        )
        for k, p in written.items():
            print(f"[59h] export {k}: {p}", flush=True)
    else:
        print(f"[59h] EXPORT SKIPPED — fatal/quality gates failing: "
              f"{sorted(failed & FATAL_OR_QUALITY)}", flush=True)

    print(f"[59h] wall: {time.perf_counter() - t0:.1f} s", flush=True)
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
