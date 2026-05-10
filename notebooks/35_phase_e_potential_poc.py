"""PoC #35: Stage 1 of "true medial-axis Phase E" — potential analysis.

Loads two real meshes:

    * ``outputs/19_tokyo_bay_oceanmesh_cleaned.14`` — A+B+C+D-cleaned
      Tokyo Bay (3,178 detector-6-flagged elements, the Phase E
      input).
    * ``outputs/29_phase_e_widened.14``               — same mesh
      after the existing centroid-widen Phase E (still 3,032
      flagged elements; 4.6 % reduction).

Runs ``analyze_under_resolved_channels`` on both and reports:

    * how many face-face-connected channels detector 6 sees;
    * the per-channel ``n_elements`` distribution;
    * the existing centroid-widen new-node cost (= n_elements);
    * the medial-axis estimate of new nodes needed to reach
      ``target_cells_across`` cells across (default 3) — the upper
      bound on Stage 2's cost.

This is *not* a new repair. Stage 1 only **measures** the gap.
Stage 2 (real medial-axis insertion + local CDT re-meshing) is
deferred. The report is the input to that decision.

Outputs:
    outputs/35_phase_e_potential_summary.txt
    outputs/35_phase_e_potential_summary.json
"""
from __future__ import annotations

import json
from pathlib import Path

from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.mesh_clean import analyze_under_resolved_channels

REPO = Path(__file__).resolve().parent.parent
OUT_DIR = REPO / "outputs"
SUMMARY_TXT = OUT_DIR / "35_phase_e_potential_summary.txt"
SUMMARY_JSON = OUT_DIR / "35_phase_e_potential_summary.json"

INPUTS = [
    ("cleaned_pre_E", OUT_DIR / "19_tokyo_bay_oceanmesh_cleaned.14"),
    ("after_centroid_widen", OUT_DIR / "29_phase_e_widened.14"),
]

TARGET_CELLS_ACROSS = 3


def _summarize(report: dict) -> str:
    comps = report["components"]
    if not comps:
        return "  (no flagged elements)"
    n_e = [c["n_elements"] for c in comps]
    n_e_sorted = sorted(n_e, reverse=True)
    head = ", ".join(str(x) for x in n_e_sorted[:8])
    return (
        f"  flagged    : {report['total_flagged_elements']:,}\n"
        f"  components : {report['n_components']:,}  "
        f"(top-8 sizes: [{head}{', ...' if len(n_e_sorted) > 8 else ''}])\n"
        f"  current Phase E new nodes        : "
        f"{report['current_phase_e_new_nodes']:,}\n"
        f"  medial-axis estimate (target={TARGET_CELLS_ACROSS} cells across) : "
        f"{report['medial_axis_new_nodes_estimate']:,}\n"
        f"  delta vs current                 : "
        f"{report['delta_nodes_vs_current']:+,}"
    )


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    payload: dict = {"target_cells_across": TARGET_CELLS_ACROSS, "runs": []}
    out_lines: list[str] = [
        "PoC #35 Phase E Stage 1 — medial-axis potential analysis",
        f"target_cells_across = {TARGET_CELLS_ACROSS}",
        "",
    ]
    for label, path in INPUTS:
        if not path.exists():
            print(f"[35] skipping missing {path}")
            continue
        mesh = read_fort14(path)
        print(
            f"[35] {label}: NP={mesh.n_nodes:,}  NE={mesh.n_elements:,}"
        )
        report = analyze_under_resolved_channels(
            mesh, target_cells_across=TARGET_CELLS_ACROSS,
        )
        # Drop per-component records from the printed summary; keep
        # them in the JSON for downstream plotting.
        payload["runs"].append({
            "label": label,
            "path": str(path.resolve()),
            "n_nodes": int(mesh.n_nodes),
            "n_elements": int(mesh.n_elements),
            "report": report,
        })
        out_lines.append(f"=== {label} ({path.name}) ===")
        out_lines.append(f"  NP={mesh.n_nodes:,}  NE={mesh.n_elements:,}")
        out_lines.append(_summarize(report))
        out_lines.append("")

    SUMMARY_JSON.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    SUMMARY_TXT.write_text("\n".join(out_lines), encoding="utf-8")
    print()
    print("\n".join(out_lines))
    print(f"wrote {SUMMARY_TXT}")
    print(f"wrote {SUMMARY_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
