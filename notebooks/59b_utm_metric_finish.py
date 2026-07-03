"""PoC #59b — UTM 54N projection + perpfix + metric-space Phase H finish.

Consumes the PoC #59a structural-finalized mesh (lon/lat) and produces
the first FVCOM-CARTESIAN-space Tokyo Bay mesh intended to pass the
full ``fmesh-mesh-qa`` gate:

  1. project EPSG:4326 -> EPSG:32654 (UTM 54N, WGS84; the kickoff
     fixes "UTM 54N" without a datum — WGS84 assumed and recorded)
  2. perpfix every open segment (first-ring perpendicular projection
     + second-ring Laplacian absorption)
  3. ``phase_h_finish`` in metric space — the first QA run showed
     ~55 metric-space C1 violations that are invisible in lon/lat
     (Phase H optimized in lon/lat; FVCOM runs CARTESIAN)
  4. re-gate; if only perpendicularity regressed, alternate a light
     perpfix / finish cycle (bounded)

Outputs: outputs/59b_utm54n_finished.14 (+ _qa.json, cycle log).
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
from pyproj import Transformer

from fvcom_mesh_tools.algorithms import align_open_boundary_first_ring
from fvcom_mesh_tools.io import read_fort14, write_fort14
from fvcom_mesh_tools.mesh_clean_phase_h import phase_h_finish
from fvcom_mesh_tools.qa import format_report, run_qa

REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "outputs" / "59a_structural_finalized.14"
OUT = REPO / "outputs" / "59b_utm54n_finished.14"
QA_JSON = REPO / "outputs" / "59b_utm54n_finished_qa.json"
CYCLES_JSON = REPO / "outputs" / "59b_cycle_log.json"

MAX_CYCLES = 3
QUALITY_IDS = {"c1_min_angle", "c2_max_angle", "c4_area_change", "c5_valence"}


def _failed_gate_ids(report) -> set[str]:
    return {
        c.check_id for c in report.checks
        if c.gate and not c.skipped and not c.passed
    }


def _qa_counts(report) -> dict[str, int]:
    return {
        c.check_id: int(c.n_violations)
        for c in report.checks if not c.skipped
    }


def _perpfix_all(mesh, *, alpha: float, n_iters: int):
    infos = []
    for k in range(len(mesh.open_boundaries)):
        mesh, info = align_open_boundary_first_ring(
            mesh, alpha=alpha, n_iters=n_iters,
            smooth_iters=2, smooth_alpha=0.3, segment_index=k,
        )
        infos.append(info)
    return mesh, infos


def main() -> int:
    t0 = time.perf_counter()
    mesh = read_fort14(SRC)
    print(
        f"[59b] input: NP={mesh.n_nodes:,} NE={mesh.n_elements:,} "
        f"open={[len(s) for s in mesh.open_boundaries]}",
        flush=True,
    )

    tr = Transformer.from_crs("EPSG:4326", "EPSG:32654", always_xy=True)
    x, y = tr.transform(mesh.nodes[:, 0], mesh.nodes[:, 1])
    mesh.nodes = np.column_stack([x, y])
    print(
        f"[59b] projected to EPSG:32654; bbox "
        f"({x.min():.0f}, {y.min():.0f}) - ({x.max():.0f}, {y.max():.0f})",
        flush=True,
    )

    mesh, pinfos = _perpfix_all(mesh, alpha=1.0, n_iters=3)
    print(
        f"[59b] perpfix: moved per segment = {[i.get('moved') for i in pinfos]}",
        flush=True,
    )

    cycle_log: list[dict] = []
    report = None
    for cycle in range(MAX_CYCLES):
        t1 = time.perf_counter()
        mesh, finfo = phase_h_finish(mesh, seed=42 + cycle)
        finish_wall = time.perf_counter() - t1
        print(
            f"[59b] cycle {cycle}: phase_h_finish "
            f"{finfo.get('before')} -> {finfo.get('after')} "
            f"({finish_wall:.0f} s)",
            flush=True,
        )

        report = run_qa(mesh, name=OUT.name, path=OUT)
        failed = _failed_gate_ids(report)
        cycle_log.append({
            "cycle": cycle,
            "finish_before": finfo.get("before"),
            "finish_after": finfo.get("after"),
            "finish_wall_s": finish_wall,
            "failed_gates": sorted(failed),
            "counts": _qa_counts(report),
        })
        print(f"[59b] cycle {cycle}: failed gates = {sorted(failed)}", flush=True)

        if not failed:
            break
        if failed == {"obc_perpendicularity"}:
            # Quality holds; only the ring drifted. Re-align lightly and
            # re-gate without another full finish.
            mesh, _ = _perpfix_all(mesh, alpha=0.5, n_iters=2)
            report = run_qa(mesh, name=OUT.name, path=OUT)
            failed = _failed_gate_ids(report)
            cycle_log.append({
                "cycle": f"{cycle}-perpfix",
                "failed_gates": sorted(failed),
                "counts": _qa_counts(report),
            })
            print(
                f"[59b] cycle {cycle}: after light perpfix, "
                f"failed gates = {sorted(failed)}",
                flush=True,
            )
            if not failed:
                break
        if not (failed & QUALITY_IDS) and cycle_log:
            # Non-quality residual that another finish pass cannot help;
            # stop iterating and surface it.
            break

    mesh.title = "PoC 59b UTM54N (EPSG:32654) perpfix + metric phase_h_finish"
    write_fort14(mesh, OUT)
    QA_JSON.write_text(
        json.dumps(report.to_dict(), indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    CYCLES_JSON.write_text(
        json.dumps(cycle_log, indent=2, ensure_ascii=False), encoding="utf-8",
    )
    print(format_report(report, lang="ja"), flush=True)
    print(f"[59b] wrote {OUT}", flush=True)
    print(f"[59b] wall: {time.perf_counter() - t0:.1f} s", flush=True)
    return 0 if report.passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
