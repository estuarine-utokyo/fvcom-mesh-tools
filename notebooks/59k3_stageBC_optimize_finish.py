"""PoC #59k-3 — stages B+C of the v3 finishing chain (checkpointed).

Stage A (59k2) measured the raw metric-space residual of the v3 mesh
at C1=615 / C2=38 / C4=628 / C5=3 — an order of magnitude above what
``phase_h_finish`` is designed for (the #59k single-job attempt died
in finish's super-linear vertex-remove stage). This job restores the
validated chain order:

  B. deterministic ``phase_h_optimize`` (Pass A/B + F + G, FVCOM
     targets 30/130/0.5/8) to grind the bulk down — checkpointed;
  C. ``phase_h_finish`` on the residual, then the west+south arc
     OBC, the QA convergence loop, and the tokyo_bay_v3 export
     (with OBC sponge).
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
from fvcom_mesh_tools.mesh_clean_phase_h import (
    phase_h_finish,
    phase_h_optimize,
)
from fvcom_mesh_tools.qa import format_report, fvcom_boundary_element_flags, run_qa

REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "outputs" / "59k_stageA_utm.14"
STAGE_B = REPO / "outputs" / "59k_stageB_optimized.14"
OUT = REPO / "outputs" / "59k_v3_gate_passed.14"
QA_JSON = REPO / "outputs" / "59k_v3_gate_passed_qa.json"
LOG_JSON = REPO / "outputs" / "59k_cycle_log.json"
EXPORT_DIR = REPO / "outputs" / "fvcom_inputs_v3"
CASENAME = "tokyo_bay_v3"

BAND_DEG = 0.008
LAND_IBTYPE = 20
MAX_CYCLES = 16
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

_TO_LL = Transformer.from_crs("EPSG:32654", "EPSG:4326", always_xy=True)


def _arc(mesh):
    lon, lat = _TO_LL.transform(mesh.nodes[:, 0], mesh.nodes[:, 1])
    loops = chain_edges_to_loops(boundary_edges_from_tris(mesh.elements))
    outer = outer_loop(loops, mesh.nodes)
    ring = outer[:-1]
    rlon, rlat = lon[ring], lat[ring]
    mask = (rlat <= rlat.min() + BAND_DEG) | (rlon <= rlon.min() + BAND_DEG)
    idx = np.where(mask)[0]
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
        raise SystemExit("[59k3] open arc trimmed away")
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
    print(f"[59k3] stage B input: NP={mesh.n_nodes:,} NE={mesh.n_elements:,}",
          flush=True)

    mesh, oinfo = phase_h_optimize(
        mesh,
        min_angle_target=30.0,
        max_angle_target=130.0,
        pass_f_enabled=True,
        pass_g_enabled=True,
        pass_g_min_angle_target=30.0,
    )
    print(f"[59k3] optimize done in {time.perf_counter() - t0:.0f} s: "
          f"rounds={oinfo.get('rounds', '?')}", flush=True)
    mesh, cinfo = compact_nodes(mesh)
    mesh.title = "59k stage B: v3 optimized (Pass A/B+F+G, 30/130)"
    write_fort14(mesh, STAGE_B)
    print(f"[59k3] wrote {STAGE_B} (compacted {cinfo['n_orphans_removed']})",
          flush=True)

    t1 = time.perf_counter()
    mesh, finfo = phase_h_finish(mesh, seed=90)
    mesh, cinfo = compact_nodes(mesh)
    print(f"[59k3] finish: {finfo.get('before')} -> {finfo.get('after')} "
          f"({time.perf_counter() - t1:.0f} s; compacted "
          f"{cinfo['n_orphans_removed']})", flush=True)

    ring, islands, open_end = _arc(mesh)
    trim = 0
    mesh = _apply_obc(mesh, ring, islands, open_end, trim)
    seg = mesh.open_boundaries[0]
    d = mesh.depths[seg]
    print(f"[59k3] arc OBC: {len(seg)} nodes, depth {d.min():.1f}-{d.max():.1f} m",
          flush=True)

    log = []
    report = run_qa(mesh, name=OUT.name, path=OUT)
    for cycle in range(MAX_CYCLES):
        failed = _failed(report)
        counts = {c.check_id: int(c.n_violations) for c in report.checks
                  if c.gate and not c.skipped and c.n_violations}
        print(f"[59k3] cycle {cycle}: failed = {sorted(failed)} {counts}",
              flush=True)
        log.append({"cycle": cycle, "failed": sorted(failed), "counts": counts})
        if not failed:
            break
        if failed & JUNCTION_IDS:
            if trim < MAX_END_TRIMS:
                trim += 1
                mesh = _apply_obc(mesh, ring, islands, open_end, trim)
                print(f"[59k3] cycle {cycle}: trimmed arc (trim={trim})",
                      flush=True)
            else:
                flags = fvcom_boundary_element_flags(mesh)
                bad = flags["r4_mask"] | flags["fake_open_mask"]
                if not bad.any():
                    break
                mesh = remove_elements(mesh, ~bad)
                mesh, _ = keep_components(mesh)
                ring, islands, open_end = _arc(mesh)
                trim = 0
                mesh = _apply_obc(mesh, ring, islands, open_end, trim)
                print(f"[59k3] cycle {cycle}: deleted {int(bad.sum())} R4/fake",
                      flush=True)
        elif failed & QUALITY_IDS:
            mesh, finfo = phase_h_finish(mesh, seed=91 + cycle)
            mesh, cinfo = compact_nodes(mesh)
            print(f"[59k3] cycle {cycle}: finish {finfo.get('before')} -> "
                  f"{finfo.get('after')}; compacted "
                  f"{cinfo['n_orphans_removed']}", flush=True)
            if cinfo["n_orphans_removed"]:
                ring, islands, open_end = _arc(mesh)
                trim = 0
                mesh = _apply_obc(mesh, ring, islands, open_end, trim)
        elif "obc_perpendicularity" in failed:
            mesh, pinfo = align_open_boundary_local(mesh, seed=9000 + cycle)
            print(f"[59k3] cycle {cycle}: perp fixer accepted="
                  f"{pinfo['accepted_total']} remaining="
                  f"{len(pinfo['remaining'])}", flush=True)
            if pinfo["accepted_total"] == 0:
                report = run_qa(mesh, name=OUT.name, path=OUT)
                break
        else:
            break
        report = run_qa(mesh, name=OUT.name, path=OUT)

    mesh.title = "PoC 59k UTM54N Tokyo Bay v3 (100 m coastline-fidelity)"
    write_fort14(mesh, OUT)
    QA_JSON.write_text(json.dumps(report.to_dict(), indent=2, ensure_ascii=False),
                       encoding="utf-8")
    LOG_JSON.write_text(json.dumps(log, indent=2, ensure_ascii=False),
                        encoding="utf-8")
    print(format_report(report, lang="ja"), flush=True)

    failed = _failed(report)
    if not (failed & FATAL_OR_QUALITY):
        lon, lat = _TO_LL.transform(mesh.nodes[:, 0], mesh.nodes[:, 1])
        obc = mesh.open_boundaries[0]
        sponge = [(int(v), 3000.0, 0.001) for v in obc]
        written = export_fvcom_case(
            mesh, EXPORT_DIR, CASENAME, obc_type=1, cor=lat, sponge=sponge,
        )
        for k, p in written.items():
            print(f"[59k3] export {k}: {p}", flush=True)
    else:
        print(f"[59k3] EXPORT SKIPPED: {sorted(failed & FATAL_OR_QUALITY)}",
              flush=True)
    print(f"[59k3] wall: {time.perf_counter() - t0:.1f} s", flush=True)
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
