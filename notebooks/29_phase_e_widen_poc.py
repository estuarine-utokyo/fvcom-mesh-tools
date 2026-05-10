"""PoC #29: Phase E (widen under-resolved channels) end-to-end validation.

Loads the PoC #19 cleaned Tokyo-Bay mesh (4-phase clean output, no
Phase E run yet), measures the detector-6 flag count before, runs
Phase E with mode='widen' via the high-level ``clean_mesh`` driver
(all other phases off), and reports the after-state. Confirms:

    1. Topology: the output mesh's NP grows by ``n_widened`` and NE
       grows by ``2 * n_widened`` (each flagged triangle becomes 3
       sub-triangles fanning from a new interior centroid).
    2. Boundaries are preserved (open/land segment counts unchanged).
    3. Detector 6's flagged-element count is reduced — but only
       modestly. Centroid insertion shrinks the local edge length
       (h_local → ≈ 0.577 × original) without changing the geometric
       channel width, so w/h ≈ 1.73 × the original ratio. Only
       elements whose original ratio sat in [min_w_h/1.73, min_w_h]
       cross the threshold; very narrow channels (ratio well below
       min_w_h) remain flagged. On the PoC #19 cleaned mesh,
       ``3,178 → 3,032`` (4.6 % reduction); the per-mesh flagged
       fraction drops from 6.7 % to 5.6 % because NE grows.

So Phase E in widen mode is best understood as "lift local
resolution one step" rather than "guarantee 3 cells across every
narrow channel". Achieving the latter would require inserting
interior nodes along the channel medial axis, which is a deeper
remeshing operation outside the scope of clean_mesh.

The run also writes ``outputs/29_phase_e_widened.14`` for downstream
inspection and ``outputs/29_phase_e_summary.json`` with the
before/after detector-6 counts and clean info.
"""
from __future__ import annotations

import json
from pathlib import Path

from fvcom_mesh_tools.diagnostics import (
    DEFAULT_MIN_W_H,
    under_resolved_channels_flag,
)
from fvcom_mesh_tools.io import read_fort14, write_fort14
from fvcom_mesh_tools.mesh_clean import clean_mesh

REPO = Path(__file__).resolve().parent.parent
INPUT = REPO / "outputs" / "19_tokyo_bay_oceanmesh_cleaned.14"
OUTPUT = REPO / "outputs" / "29_phase_e_widened.14"
SUMMARY = REPO / "outputs" / "29_phase_e_summary.json"

# Tokyo Bay DEM bbox used in PoC #19; Phase E does not delete here, so
# this is only required by the ``clean_mesh`` API. The bbox-tol-driven
# rebuild is a no-op for widen.
BBOX = (139.46, 34.99, 140.10, 35.74)


def _detector6_count(mesh, *, min_w_h: float) -> int:
    flag, _info = under_resolved_channels_flag(mesh, min_w_h=min_w_h)
    return int(flag.sum())


def main() -> int:
    print(f"loading: {INPUT}")
    mesh_in = read_fort14(INPUT)
    print(
        f"  NP={mesh_in.n_nodes:,}  NE={mesh_in.n_elements:,}  "
        f"open={len(mesh_in.open_boundaries)}  "
        f"land={len(mesh_in.land_boundaries)}"
    )
    n_before = _detector6_count(mesh_in, min_w_h=DEFAULT_MIN_W_H)
    print(f"detector 6 flagged BEFORE: {n_before:,}")

    print("running clean_mesh with under_resolved_mode='widen' ...")
    out, info = clean_mesh(
        mesh_in,
        bbox=BBOX,
        bbox_tol_m=150.0,
        remove_disjoint=False,
        trim_dead_ends_iters=0,
        thin_chain_mode="none",
        repair_overconnected_iters=0,
        under_resolved_mode="widen",
        under_resolved_min_w_h=DEFAULT_MIN_W_H,
    )
    print(
        f"  NP {info['input']['n_nodes']:,} -> {info['output']['n_nodes']:,}  "
        f"NE {info['input']['n_elements']:,} -> {info['output']['n_elements']:,}"
    )

    phase_e = next(
        p for p in info["phases"] if p["name"] == "repair_under_resolved_channels"
    )
    print(f"phase E: {phase_e}")

    n_after = _detector6_count(out, min_w_h=DEFAULT_MIN_W_H)
    print(f"detector 6 flagged AFTER : {n_after:,}")
    if n_before > 0:
        print(f"reduction: {(n_before - n_after) / n_before * 100:.1f}%")

    # Topology checks
    assert out.n_elements == mesh_in.n_elements + 2 * phase_e["n_widened"], (
        f"NE growth mismatch: {out.n_elements - mesh_in.n_elements} vs "
        f"2 * {phase_e['n_widened']}"
    )
    assert out.n_nodes == mesh_in.n_nodes + phase_e["n_widened"], (
        "NP growth should equal n_widened"
    )
    print("topology growth checks: OK")
    assert len(out.open_boundaries) == len(mesh_in.open_boundaries)
    assert len(out.land_boundaries) == len(mesh_in.land_boundaries)
    print("boundary count preserved: OK")

    write_fort14(out, OUTPUT)
    print(f"wrote {OUTPUT}")

    SUMMARY.write_text(
        json.dumps(
            {
                "input_path": str(INPUT.resolve()),
                "output_path": str(OUTPUT.resolve()),
                "n_flagged_before": n_before,
                "n_flagged_after": n_after,
                "reduction_pct": (
                    (n_before - n_after) / n_before * 100
                    if n_before > 0 else 0.0
                ),
                **info,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"wrote {SUMMARY}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
