"""Compare the coastal-resolution sweep variants against the certified mesh.

One row per variant: cost (nodes/elements), the achieved coastal resolution,
the time step, the QA verdict and the defect counts. Reads the collected
outputs of jobs/octopus/392_sizing_sweep.sh.

Usage: python notebooks/392_sizing_summary.py outputs/sizing_392
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
COLLECT = Path(sys.argv[1] if len(sys.argv) > 1 else ROOT / "outputs/sizing_392")

CERT = ROOT / "outputs/sample_repro/sample_repro_final.14"
# The certified numbers come from two logs: 388 regenerated the certified
# mesh (QA, connectivity, one-wide) and 390 re-measured it with the corrected
# detectors (narrow water, walls, coverage gaps).
CERT_LOGS = (sorted((ROOT / "logs").glob("388_cert_regen.*.log"))
             + sorted((ROOT / "logs").glob("390_recheck.*.log"))[-1:])


def mesh_stats(path: Path) -> dict:
    lines = Path(path).read_text().splitlines()
    ne, nn = map(int, lines[1].split()[:2])
    nod = np.array([lines[2 + i].split()[1:4] for i in range(nn)], float)
    tri = np.array([lines[2 + nn + i].split()[2:5] for i in range(ne)], int) - 1
    xy, dep = nod[:, :2], nod[:, 2]
    e = np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]])
    e.sort(axis=1)
    uniq, cnt = np.unique(e, axis=0, return_counts=True)
    length = np.hypot(*(xy[uniq[:, 0]] - xy[uniq[:, 1]]).T)
    bnd = uniq[cnt == 1]
    lb = np.hypot(*(xy[bnd[:, 0]] - xy[bnd[:, 1]]).T)
    dmid = dep[uniq].mean(1)
    shallow = dmid <= 10
    # implied external time step, as fmesh-mesh-qa computes it
    p = xy[tri]
    ar = 0.5 * np.abs((p[:, 1, 0] - p[:, 0, 0]) * (p[:, 2, 1] - p[:, 0, 1])
                      - (p[:, 1, 1] - p[:, 0, 1]) * (p[:, 2, 0] - p[:, 0, 0]))
    side = np.stack([np.hypot(*(p[:, 1] - p[:, 0]).T), np.hypot(*(p[:, 2] - p[:, 1]).T),
                     np.hypot(*(p[:, 0] - p[:, 2]).T)], axis=1)
    alt = 2 * ar[:, None] / np.maximum(side, 1e-9)
    hmax = np.maximum(dep[tri].max(1), 0.1)
    dt = float(np.min(alt.min(1) / np.sqrt(9.81 * hmax)))
    return dict(
        nodes=nn, elements=ne,
        edge_median=float(np.median(length)),
        coast_median=float(np.median(lb)), coast_p05=float(np.percentile(lb, 5)),
        shallow_median=float(np.median(length[shallow])),
        dt_s=dt,
    )


def log_stats(log) -> dict:
    paths = [log] if isinstance(log, (str, Path)) else list(log)
    text = "\n".join(Path(p).read_text(errors="replace") for p in paths if Path(p).exists())

    def first(pattern, cast=int, default=None):
        hits = re.findall(pattern, text)
        return cast(hits[0]) if hits else default

    return dict(
        qa=first(r"総合判定: (\S+?)\s", str, "?"),
        narrow=first(r"NARROW water meshed \(w/h<1, outside kept corridors\): (\d+)"),
        severe=first(r"of which severe \(w/h<0\.5\): (\d+)"),
        wall=first(r"WALL crossings: (\d+)"),
        one_wide=first(r"confirmed one-wide: (\d+)"),
        chokes=first(r"real chokes: (\d+) sites"),
        gaps=first(r"clear defects=(\d+)"),
        critical=sum(int(x) for x in re.findall(r"CRITICAL [a-z -]+:\s+(\d+)", text)) or 0,
    )


rows = {"certified (H0=290, grade=0.165)": mesh_stats(CERT) | log_stats(CERT_LOGS)}
for d in sorted(x for x in COLLECT.iterdir() if x.is_dir()):
    mesh = d / "sample_repro_final.14"
    if mesh.exists():
        rows[f"variant {d.name}"] = mesh_stats(mesh) | log_stats(d / "run.log")

cols = ["nodes", "elements", "coast_median", "coast_p05", "shallow_median", "edge_median",
        "dt_s", "qa", "critical", "one_wide", "chokes", "narrow", "severe", "wall", "gaps"]
w = max(len(k) for k in rows) + 2
print(f"{'variant':{w}}" + "".join(f"{c:>15}" for c in cols))
print("-" * (w + 15 * len(cols)))
for label, r in rows.items():
    out = []
    for c in cols:
        v = r.get(c)
        out.append(f"{v:15.1f}" if isinstance(v, float) else f"{str(v):>15}")
    print(f"{label:{w}}" + "".join(out))
print("\ncoast_* are mesh boundary edge lengths (m); shallow_median is the median edge"
      " length where the mean node depth is <= 10 m; dt_s is the implied external"
      " time step (min altitude / sqrt(g*max depth)).")

summary = {k: v for k, v in rows.items()}
COLLECT.mkdir(parents=True, exist_ok=True)
(COLLECT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
