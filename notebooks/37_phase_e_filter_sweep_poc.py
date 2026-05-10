"""PoC #37: ``--min-channel-elements`` sweep for Phase E Stage 2 Go/No-Go.

PoC #35 measured the gap between the existing centroid-widen Phase E
and a hypothetical medial-axis Stage 2: 3,178 detector-6-flagged
elements split into 1,010 components with mean ~3 elements / component
on the cleaned PoC #19 Tokyo Bay mesh. The conclusion was that the
flagged elements are mostly small isolated clusters at river-mouth
corners / jetty tips, *not* the long ribbon-like inlets a real
medial-axis insertion would refine — so Stage 2 was deferred and the
``--min-channel-elements`` filter was landed (commit ``f48d137``) as
the immediate follow-up.

This PoC reruns the analysis with that filter applied at several
thresholds (1, 2, 3, 5, 10, 20) and reports, per threshold:

    * post-filter total_flagged_elements / n_components;
    * the per-component long_axis_m / n_elements distribution
      (top-10 components and aggregate stats);
    * existing centroid-widen new-node cost vs medial-axis estimate;
    * the medial-axis-vs-centroid ratio — the leading indicator of
      whether Stage 2 (real CDT re-meshing) actually wins on the
      *interesting* channels that survive the filter.

The decision criterion documented in the project's task tracker:
**Stage 2 is justified if filtered (e.g. min_channel_elements >= 3)
channels still total > ~500 elements with mean long_axis_m > ~3 ×
h_local_median.** The post-filter "channel landscape" is the input
to that decision.

Outputs:
    outputs/37_phase_e_filter_sweep_summary.txt
    outputs/37_phase_e_filter_sweep_summary.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.mesh_clean import analyze_under_resolved_channels

REPO = Path(__file__).resolve().parent.parent
OUT_DIR = REPO / "outputs"
SUMMARY_TXT = OUT_DIR / "37_phase_e_filter_sweep_summary.txt"
SUMMARY_JSON = OUT_DIR / "37_phase_e_filter_sweep_summary.json"

INPUTS = [
    ("cleaned_pre_E", OUT_DIR / "19_tokyo_bay_oceanmesh_cleaned.14"),
    ("after_centroid_widen", OUT_DIR / "29_phase_e_widened.14"),
]

FILTER_VALUES = [1, 2, 3, 5, 10, 20]
TARGET_CELLS_ACROSS = 3


def _component_stats(report: dict) -> dict:
    """Aggregate component-level numbers used for the Stage-2 decision."""
    comps = report["components"]
    if not comps:
        return {
            "n_components": 0,
            "mean_n_elements": 0.0,
            "max_n_elements": 0,
            "mean_long_axis_m": 0.0,
            "max_long_axis_m": 0.0,
            "mean_long_axis_over_h": 0.0,
            "max_long_axis_over_h": 0.0,
            "top_n_elements": [],
        }
    n_e = np.asarray([c["n_elements"] for c in comps], dtype=float)
    long_axis = np.asarray([c["long_axis_m"] for c in comps], dtype=float)
    h_med = np.asarray([c["h_local_median_m"] for c in comps], dtype=float)
    valid = h_med > 0
    ratio = np.where(valid, long_axis / np.where(valid, h_med, 1.0), 0.0)
    return {
        "n_components": int(n_e.size),
        "mean_n_elements": float(n_e.mean()),
        "max_n_elements": int(n_e.max()),
        "mean_long_axis_m": float(long_axis.mean()),
        "max_long_axis_m": float(long_axis.max()),
        "mean_long_axis_over_h": float(ratio[valid].mean()) if valid.any() else 0.0,
        "max_long_axis_over_h": float(ratio.max()),
        "top_n_elements": sorted([int(x) for x in n_e], reverse=True)[:10],
    }


def _format_table(rows: list[dict]) -> list[str]:
    headers = [
        "min_n", "flagged", "comps", "mean_n_e", "max_n_e",
        "mean_la_m", "mean_la/h", "centroid_n", "medial_n",
        "ratio_m/c",
    ]
    out = ["  " + "  ".join(h.rjust(11) for h in headers),
           "  " + "  ".join("-" * 11 for _ in headers)]
    for r in rows:
        s = r["stats"]
        ratio = (
            r["medial_axis_new_nodes_estimate"]
            / r["current_phase_e_new_nodes"]
            if r["current_phase_e_new_nodes"] > 0 else 0.0
        )
        out.append("  " + "  ".join([
            f"{r['min_channel_elements']:>11d}",
            f"{r['total_flagged_elements']:>11,}",
            f"{s['n_components']:>11,}",
            f"{s['mean_n_elements']:>11.2f}",
            f"{s['max_n_elements']:>11,}",
            f"{s['mean_long_axis_m']:>11.0f}",
            f"{s['mean_long_axis_over_h']:>11.2f}",
            f"{r['current_phase_e_new_nodes']:>11,}",
            f"{r['medial_axis_new_nodes_estimate']:>11,}",
            f"{ratio:>11.2f}",
        ]))
    return out


def _analyse_one(label: str, path: Path) -> dict:
    print(f"\n[37] === {label} ({path.name}) ===")
    mesh = read_fort14(path)
    print(f"[37]   NP={mesh.n_nodes:,}  NE={mesh.n_elements:,}")
    rows: list[dict] = []
    for n in FILTER_VALUES:
        report = analyze_under_resolved_channels(
            mesh,
            target_cells_across=TARGET_CELLS_ACROSS,
            min_channel_elements=n,
        )
        stats = _component_stats(report)
        cur = report["current_phase_e_new_nodes"]
        med = report["medial_axis_new_nodes_estimate"]
        ratio_txt = f"{med / cur:.2f}" if cur > 0 else "-"
        rows.append({
            "min_channel_elements": n,
            "total_flagged_elements": report["total_flagged_elements"],
            "current_phase_e_new_nodes": cur,
            "medial_axis_new_nodes_estimate": med,
            "delta_nodes_vs_current": report["delta_nodes_vs_current"],
            "stats": stats,
        })
        print(
            f"[37]   min_n={n:>3d}  "
            f"flagged={report['total_flagged_elements']:>5,}  "
            f"comps={stats['n_components']:>5,}  "
            f"mean_la/h={stats['mean_long_axis_over_h']:>5.2f}  "
            f"medial/centroid={ratio_txt}"
        )
    return {
        "label": label,
        "path": str(path.resolve()),
        "n_nodes": int(mesh.n_nodes),
        "n_elements": int(mesh.n_elements),
        "rows": rows,
    }


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    payload: dict = {
        "target_cells_across": TARGET_CELLS_ACROSS,
        "filter_values": FILTER_VALUES,
        "runs": [],
    }
    out_lines: list[str] = [
        "PoC #37 --min-channel-elements sweep on detector 6 / Phase E",
        f"target_cells_across = {TARGET_CELLS_ACROSS}",
        "",
    ]
    for label, path in INPUTS:
        if not path.exists():
            print(f"[37] skipping missing {path}")
            continue
        run = _analyse_one(label, path)
        payload["runs"].append(run)
        out_lines.append(f"=== {run['label']} ({Path(run['path']).name}) ===")
        out_lines.append(
            f"  NP={run['n_nodes']:,}  NE={run['n_elements']:,}"
        )
        out_lines.append("")
        out_lines.extend(_format_table(run["rows"]))
        out_lines.append("")

    SUMMARY_JSON.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    SUMMARY_TXT.write_text("\n".join(out_lines), encoding="utf-8")
    print()
    print("\n".join(out_lines))
    print(f"\nwrote {SUMMARY_TXT}")
    print(f"wrote {SUMMARY_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
