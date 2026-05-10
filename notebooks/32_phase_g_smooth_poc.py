"""PoC #32: Phase G (Laplacian smoothing) end-to-end validation.

Loads the PoC #19 4-phase-cleaned Tokyo-Bay mesh and runs Phase G
(``oceanmesh.laplacian2`` wrapped by
``fvcom_mesh_tools.mesh_clean.smooth_mesh_laplacian``) at three
iteration / tolerance presets. Reports:

    * Maximum and mean node displacement (degrees, ≈ metres at the
      mesh latitude).
    * Element-quality changes: alpha mean, frac<20°, min-angle p05.
    * Topology preservation (NP / NE / boundary counts unchanged).

The point is to characterise how aggressive the Laplacian sweep
needs to be on a real cleaned mesh, so users have starting points
for ``--smooth-laplacian-iters`` and ``--smooth-laplacian-tol``.

Settings tested:

    * default: ``max_iter=20, tol=0.01`` — oceanmesh defaults.
    * gentle:  ``max_iter=5,  tol=0.01`` — early-exit.
    * deep:    ``max_iter=50, tol=1e-4`` — let it converge.

Outputs:
    outputs/32_phase_g_summary.txt
    outputs/32_phase_g_summary.json
    outputs/32_tokyo_bay_phase_g_default.14
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from fvcom_mesh_tools.algorithms import alpha_quality, min_interior_angle
from fvcom_mesh_tools.io import read_fort14, write_fort14
from fvcom_mesh_tools.mesh_clean import smooth_mesh_laplacian

REPO = Path(__file__).resolve().parent.parent
INPUT = REPO / "outputs" / "19_tokyo_bay_oceanmesh_cleaned.14"
OUTPUT_F14 = REPO / "outputs" / "32_tokyo_bay_phase_g_default.14"
SUMMARY_TXT = REPO / "outputs" / "32_phase_g_summary.txt"
SUMMARY_JSON = REPO / "outputs" / "32_phase_g_summary.json"


def _quality(mesh) -> dict[str, float]:
    if mesh.n_elements == 0:
        return {"n_elements": 0}
    a = alpha_quality(mesh)
    # min_interior_angle returns degrees by default.
    ang = min_interior_angle(mesh)
    return {
        "n_elements": int(mesh.n_elements),
        "n_nodes": int(mesh.n_nodes),
        "alpha_mean": float(a.mean()),
        "alpha_p05": float(np.percentile(a, 5)),
        "min_angle_p05_deg": float(np.percentile(ang, 5)),
        "min_angle_p50_deg": float(np.median(ang)),
        "frac_lt_20deg": float((ang < 20.0).sum() / len(ang)),
    }


def main() -> int:
    print(f"loading: {INPUT}")
    mesh_in = read_fort14(INPUT)
    base = _quality(mesh_in)
    print(
        f"  NP={mesh_in.n_nodes:,}  NE={mesh_in.n_elements:,}  "
        f"open={len(mesh_in.open_boundaries)}  "
        f"land={len(mesh_in.land_boundaries)}"
    )
    print(f"  before: {base}")

    presets = [
        {"label": "default", "max_iter": 20, "tol": 0.01},
        {"label": "gentle",  "max_iter":  5, "tol": 0.01},
        {"label": "deep",    "max_iter": 50, "tol": 1e-4},
    ]
    results = {"input": {"path": str(INPUT.resolve()), **base}, "runs": []}

    default_out = None
    for r in presets:
        out, info = smooth_mesh_laplacian(
            mesh_in, max_iter=r["max_iter"], tol=r["tol"],
        )
        after = _quality(out)
        results["runs"].append({
            "label": r["label"],
            "max_iter": r["max_iter"],
            "tol": r["tol"],
            "phase_g_info": info,
            "after": after,
            # Topology preservation checks.
            "topology_preserved": (
                out.n_nodes == mesh_in.n_nodes
                and out.n_elements == mesh_in.n_elements
                and len(out.open_boundaries) == len(mesh_in.open_boundaries)
                and len(out.land_boundaries) == len(mesh_in.land_boundaries)
            ),
        })
        print(
            f"  [{r['label']:<8}] iters={r['max_iter']:>2}  tol={r['tol']:.0e}  "
            f"moved={info['n_nodes_moved']:,}  "
            f"max_disp={info['displacement_max']:.3e}°  "
            f"alpha {base['alpha_mean']:.4f} -> {after['alpha_mean']:.4f}  "
            f"frac<20° {base['frac_lt_20deg']:.4%} -> {after['frac_lt_20deg']:.4%}"
        )
        if r["label"] == "default":
            default_out = out

    if default_out is not None:
        write_fort14(default_out, OUTPUT_F14)
        print(f"wrote {OUTPUT_F14}")

    SUMMARY_JSON.write_text(json.dumps(results, indent=2), encoding="utf-8")

    lines = [f"PoC #32 Phase G sweep on {INPUT}"]
    lines.append(
        f"input: NP={base['n_nodes']:,}  NE={base['n_elements']:,}  "
        f"alpha={base['alpha_mean']:.4f}  "
        f"min_p05={base['min_angle_p05_deg']:.2f}°  "
        f"frac<20°={base['frac_lt_20deg']:.4%}"
    )
    for entry in results["runs"]:
        a = entry["after"]
        gi = entry["phase_g_info"]
        lines.append(
            f"  {entry['label']:<8} iter={entry['max_iter']:>2}  "
            f"tol={entry['tol']:>7.0e}  "
            f"moved={gi['n_nodes_moved']:>6,}  "
            f"max_disp={gi['displacement_max']:.3e}°  "
            f"alpha={a['alpha_mean']:.4f}  "
            f"min_p05={a['min_angle_p05_deg']:>5.2f}°  "
            f"frac<20°={a['frac_lt_20deg']:.4%}  "
            f"topo_ok={entry['topology_preserved']}"
        )
    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print()
    print("\n".join(lines))
    print(f"wrote {SUMMARY_TXT}")
    print(f"wrote {SUMMARY_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
