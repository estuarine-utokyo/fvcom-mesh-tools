"""PoC #31: Phase F (skewed-element removal) end-to-end validation.

Loads the PoC #19 4-phase-cleaned Tokyo-Bay mesh and runs Phase F
(``ocsmesh.utils.cleanup_skewed_el`` wrapped by
``fvcom_mesh_tools.mesh_clean.repair_skewed_elements``) over a few
threshold configurations. Reports the count of triangles deleted at
each setting, plus the change in the global angle distribution
(min interior angle p01/p05/p50, max p95/p99).

The point is to characterise how aggressive the angle thresholds
need to be to make a measurable difference on a real cleaned mesh,
so users have a starting point for ``--repair-skewed-min-angle-deg``
and ``--repair-skewed-max-angle-deg``.

Settings tested:

    * ocsmesh defaults: ``[1.0, 175.0]`` (very permissive — only
      catches near-degenerate slivers).
    * Conservative FVCOM-friendly: ``[5.0, 170.0]``.
    * Aggressive: ``[10.0, 160.0]``.

Outputs:
    outputs/31_phase_f_summary.txt
    outputs/31_phase_f_summary.json
    outputs/31_tokyo_bay_phase_f_default.14
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from fvcom_mesh_tools.algorithms import min_interior_angle
from fvcom_mesh_tools.io import read_fort14, write_fort14
from fvcom_mesh_tools.mesh_clean import (
    DEFAULT_BBOX_TOL_M,
    repair_skewed_elements,
)

REPO = Path(__file__).resolve().parent.parent
INPUT = REPO / "outputs" / "19_tokyo_bay_oceanmesh_cleaned.14"
OUTPUT_F14 = REPO / "outputs" / "31_tokyo_bay_phase_f_default.14"
SUMMARY_TXT = REPO / "outputs" / "31_phase_f_summary.txt"
SUMMARY_JSON = REPO / "outputs" / "31_phase_f_summary.json"

# PoC #19 DEM bbox.
BBOX = (139.46, 34.99, 140.10, 35.74)


def _angle_stats(mesh) -> dict[str, float]:
    """Return percentile stats over all triangles' (min interior, max
    interior) angle pairs."""
    NP = mesh.n_nodes
    elems = mesh.elements
    nodes = mesh.nodes
    if elems.size == 0:
        return {}
    # Per-element three angles via dot-product trick.
    p0 = nodes[elems[:, 0]]
    p1 = nodes[elems[:, 1]]
    p2 = nodes[elems[:, 2]]

    def _ang(a, b, c):
        u = b - a
        v = c - a
        cos_t = (u * v).sum(axis=1) / (
            np.linalg.norm(u, axis=1) * np.linalg.norm(v, axis=1) + 1e-30
        )
        return np.degrees(np.arccos(np.clip(cos_t, -1.0, 1.0)))

    a0 = _ang(p0, p1, p2)
    a1 = _ang(p1, p0, p2)
    a2 = _ang(p2, p0, p1)
    angs = np.column_stack([a0, a1, a2])
    min_a = angs.min(axis=1)
    max_a = angs.max(axis=1)
    # Also recompute with the package helper for cross-check.
    # min_interior_angle returns degrees by default.
    pkg_min = min_interior_angle(mesh)
    return {
        "n_elements": int(elems.shape[0]),
        "n_nodes": int(NP),
        "min_angle_p01_deg": float(np.percentile(min_a, 1)),
        "min_angle_p05_deg": float(np.percentile(min_a, 5)),
        "min_angle_p50_deg": float(np.median(min_a)),
        "max_angle_p50_deg": float(np.median(max_a)),
        "max_angle_p95_deg": float(np.percentile(max_a, 95)),
        "max_angle_p99_deg": float(np.percentile(max_a, 99)),
        "frac_lt_5deg": float((min_a < 5.0).sum() / len(min_a)),
        "frac_lt_10deg": float((min_a < 10.0).sum() / len(min_a)),
        "frac_lt_20deg": float((pkg_min < 20.0).sum() / len(min_a)),
        "frac_gt_170deg": float((max_a > 170.0).sum() / len(max_a)),
        "frac_gt_175deg": float((max_a > 175.0).sum() / len(max_a)),
    }


def _bbox_tol_deg(mesh) -> float:
    """Convert DEFAULT_BBOX_TOL_M to degrees at the mesh's mid-latitude."""
    EARTH_R_M = 6_371_000.0
    lat0 = float(mesh.nodes[:, 1].mean())
    return DEFAULT_BBOX_TOL_M / (
        EARTH_R_M * np.cos(np.deg2rad(lat0)) * np.pi / 180.0
    )


def main() -> int:
    print(f"loading: {INPUT}")
    mesh_in = read_fort14(INPUT)
    base = _angle_stats(mesh_in)
    print(
        f"  NP={mesh_in.n_nodes:,}  NE={mesh_in.n_elements:,}  "
        f"open={len(mesh_in.open_boundaries)}  "
        f"land={len(mesh_in.land_boundaries)}"
    )
    print(f"  before: {base}")

    runs = [
        {"label": "ocsmesh-defaults", "min_deg": 1.0, "max_deg": 175.0},
        {"label": "conservative",     "min_deg": 5.0, "max_deg": 170.0},
        {"label": "aggressive",       "min_deg": 10.0, "max_deg": 160.0},
    ]
    results = {"input": {"path": str(INPUT.resolve()), **base}, "runs": []}

    tol = _bbox_tol_deg(mesh_in)
    default_out = None
    for r in runs:
        cleaned, info = repair_skewed_elements(
            mesh_in,
            min_angle_deg=r["min_deg"],
            max_angle_deg=r["max_deg"],
            bbox=BBOX, tol_deg=tol, land_ibtype=0,
        )
        after = _angle_stats(cleaned)
        entry = {
            "label": r["label"],
            "min_angle_deg": r["min_deg"],
            "max_angle_deg": r["max_deg"],
            "phase_f_info": info,
            "after": after,
        }
        results["runs"].append(entry)
        print(
            f"  [{r['label']}] removed={info['n_elements_removed']:,}  "
            f"NE {mesh_in.n_elements:,} -> {cleaned.n_elements:,}"
        )
        if r["label"] == "ocsmesh-defaults":
            default_out = cleaned

    if default_out is not None:
        write_fort14(default_out, OUTPUT_F14)
        print(f"wrote {OUTPUT_F14}")

    SUMMARY_JSON.write_text(json.dumps(results, indent=2), encoding="utf-8")

    lines = [f"PoC #31 Phase F sweep on {INPUT}"]
    lines.append(f"input: NP={base['n_elements']:,} (NE), "
                 f"min_p05={base['min_angle_p05_deg']:.2f}°, "
                 f"max_p99={base['max_angle_p99_deg']:.2f}°, "
                 f"frac<5°={base['frac_lt_5deg']:.4%}")
    for entry in results["runs"]:
        a = entry["after"]
        lines.append(
            f"  {entry['label']:<18} thresholds=[{entry['min_angle_deg']:>5g}°, "
            f"{entry['max_angle_deg']:>5g}°]  removed="
            f"{entry['phase_f_info']['n_elements_removed']:>6,}  "
            f"NE_after={a['n_elements']:>7,}  "
            f"min_p05={a['min_angle_p05_deg']:>5.2f}°  "
            f"max_p99={a['max_angle_p99_deg']:>5.2f}°  "
            f"frac<20°={a['frac_lt_20deg']:.4%}"
        )
    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print()
    print("\n".join(lines))
    print(f"wrote {SUMMARY_TXT}")
    print(f"wrote {SUMMARY_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
