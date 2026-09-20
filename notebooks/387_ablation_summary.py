"""Summarise the per-edit ablation (jobs/octopus/387_edit_ablation_each.sh).

One row per variant: what the chain produces when exactly one of the four
human-judgment edits is excluded, against the certified mesh (all edits on)
and, when present, the all-edits-off run of job 115210.

Usage: python notebooks/387_ablation_summary.py outputs/ablation_387
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
COLLECT = Path(sys.argv[1] if len(sys.argv) > 1 else ROOT / "outputs/ablation_387")

# Certified reference and the all-off run, measured earlier on OCTOPUS.
REFERENCE = {
    "all edits ON (certified, 115144)": dict(
        log=ROOT / "logs/380_recert.115144.log",
        mesh=ROOT / "outputs/sample_repro/sample_repro_final.14",
    ),
    "all edits OFF (115210)": dict(
        log=ROOT / "logs/385_no_edits.115210.log",
        mesh=ROOT / "outputs/sample_repro_noedits/sample_repro_final.14",
    ),
}

PATTERNS = {
    "QA": (r"総合判定: (\S+?)\s", str),
    "dt_s": (r"min dt = ([\d.]+) s", float),
    "nodes": (r"節点数 ([\d,]+)", lambda s: int(s.replace(",", ""))),
    "uncovered": (r"sample elements NOT covered by our mesh: (\d+)", int),
    "on_land_unintended": (r"UNINTENDED: (\d+)", int),
    "one_wide": (r"confirmed one-wide: (\d+)", int),
    "chokes": (r"real chokes: (\d+) sites", int),
    "stray": (r"STRAY over-resolution: (\d+)", int),
    "wall": (r"WALL crossings: (\d+)", int),
}


def measure(log: Path) -> dict:
    out = {}
    if not log.exists():
        return {k: None for k in PATTERNS}
    text = log.read_text(errors="replace")
    for key, (pattern, cast) in PATTERNS.items():
        hits = re.findall(pattern, text)
        # First occurrence: job 115210 appends a re-measurement of the
        # certified mesh after the edit-free one.
        out[key] = cast(hits[0]) if hits else None
    axis = re.findall(r"\[372\].*?offset[^\d-]*(-?[\d.]+)", text)
    out["axis_hits"] = len(axis)
    return out


def haneda_gap(mesh: Path) -> int | None:
    """Arc points of edit_001 with no element centroid within 200 m.

    The Haneda D-runway passage is the one place the all-off run lost a
    through path; this counts how badly a variant loses it.
    """
    edit = ROOT / "recipes/edits/sample_repro/edit_001_haneda_d_runway.json"
    if not mesh.exists() or not edit.exists():
        return None
    from pyproj import Transformer

    fw = Transformer.from_crs(4326, 32654, always_xy=True)
    lines = mesh.read_text().splitlines()
    ne, nn = map(int, lines[1].split()[:2])
    nod = np.array([lines[2 + i].split()[1:3] for i in range(nn)], float)
    tri = np.array([lines[2 + nn + i].split()[2:5] for i in range(ne)], int) - 1
    cen = nod[tri].mean(1)
    arc = np.array([fw.transform(x, y) for x, y in json.loads(edit.read_text())["arc"]])
    gaps = [p for p in arc if not (np.hypot(cen[:, 0] - p[0], cen[:, 1] - p[1]) < 200).any()]
    return len(gaps)


rows = {}
for label, ref in REFERENCE.items():
    rows[label] = measure(ref["log"]) | {"haneda_gaps": haneda_gap(ref["mesh"])}
for d in sorted(COLLECT.glob("excl_*")):
    rows[f"without {d.name[5:]}"] = measure(d / "run.log") | {
        "haneda_gaps": haneda_gap(d / "sample_repro_final.14")
    }

cols = ["QA", "nodes", "dt_s", "uncovered", "on_land_unintended", "one_wide",
        "chokes", "stray", "wall", "haneda_gaps"]
width = max(len(k) for k in rows) + 2
print(f"{'variant':{width}}" + "".join(f"{c:>20}" for c in cols))
print("-" * (width + 20 * len(cols)))
for label, r in rows.items():
    print(f"{label:{width}}" + "".join(f"{str(r.get(c)):>20}" for c in cols))
print(
    "\nhaneda_gaps = points of the edit_001 arc (16 total) with no element "
    "centroid within 200 m: 0 means the pier passage is meshed."
)
