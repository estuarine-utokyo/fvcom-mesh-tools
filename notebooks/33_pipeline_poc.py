"""PoC #33: ``fmesh-mesh-pipeline`` end-to-end on the PoC #19 raw mesh.

Drives the pipeline through its three rungs against a realistic
"messy" starting point (the raw oceanmesh+rivers Tokyo Bay mesh —
144 dual-graph components, 5,496 disjoint elements, 3 over-connected
nodes) with FVCOM-friendly thresholds. Records which rung satisfies
the gate, the per-rung metrics, and confirms the n_flipped safety
net behaves correctly.

Inputs:
    outputs/19_tokyo_bay_oceanmesh.14

Outputs:
    outputs/33_pipeline_passing.14         — final fort.14 from the
                                             rung that passed
    outputs/33_pipeline_summary.json       — full per-rung history

Threshold preset (FVCOM 4.x conservative):
    --min-alpha 0.95
    --max-frac-lt-20deg 0.005
    --max-valence 8
    --max-flipped 0
    --max-disjoint-elems 0
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

from fvcom_mesh_tools.cli.meshpipeline import main as pipeline_main

REPO = Path(__file__).resolve().parent.parent
INPUT = REPO / "outputs" / "19_tokyo_bay_oceanmesh.14"
OUTPUT = REPO / "outputs" / "33_pipeline_passing.14"
SUMMARY = REPO / "outputs" / "33_pipeline_summary.json"


def main() -> int:
    if not INPUT.exists():
        raise SystemExit(f"missing input: {INPUT}")
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    argv = [
        str(INPUT), str(OUTPUT),
        "--bbox", "139.46", "34.99", "140.10", "35.74",
        "--bbox-tol-m", "150",
        "--open-merge-coast-gap", "50",
        "--min-alpha", "0.95",
        "--max-frac-lt-20deg", "0.005",
        "--max-valence", "8",
        "--max-flipped", "0",
        "--max-disjoint-elems", "0",
        "--summary", str(SUMMARY),
    ]
    print(f"[33] running fmesh-mesh-pipeline on {INPUT}")
    rc = pipeline_main(argv)
    print(f"[33] pipeline exit code: {rc}")

    payload = json.loads(SUMMARY.read_text())
    print()
    print(f"[33] final rung:    {payload['final']['rung_label']}")
    print(f"[33] final passed:  {payload['final']['passed']}")
    print(f"[33] rungs run:     {len(payload['history'])} of "
          f"{payload['max_iters']}")
    for h in payload["history"]:
        m = h["metrics"]
        print(
            f"[33]   {h['rung_label']:<14} "
            f"alpha={m['alpha_mean']:.4f}  "
            f"frac<20°={m['frac_lt_20deg']:.4%}  "
            f"max_v={m['max_valence']}  "
            f"flipped={m['n_flipped']}  "
            f"disjoint={m['n_disjoint_elems']}  "
            f"passed={h['passed']}"
        )
    print(f"\n[33] wrote {OUTPUT}")
    print(f"[33] wrote {SUMMARY}")
    return rc


if __name__ == "__main__":
    sys.exit(main())
