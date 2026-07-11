# Tokyo Bay varres-3r translation: pre-faithful vs faithful-stack
# comparison (NP/quality/band stats + atlas side-by-side + bay zoom).
import os, sys
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
from fvcom_mesh_tools.io import read_fort14
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DEG = 1.0 / 111e3
OUT = Path("outputs/tb_varres_3r")
old = read_fort14(str(OUT / "tb_varres_3regions_prefaithful.14"))
new = read_fort14(str(OUT / "tb_varres_3regions.14"))
meshes = {
    "pre-faithful": (old.nodes[:, :2], old.elements),
    "faithful": (new.nodes[:, :2], new.elements),
}
for tag, (p, t) in meshes.items():
    print(f"[tbcmp] {tag}: NP={len(p):,} NE={len(t):,}", flush=True)


def edge_stats(p, t):
    e = np.vstack([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]])
    e = np.unique(np.sort(e, axis=1), axis=0)
    dlon = (p[e[:, 0], 0] - p[e[:, 1], 0])
    lat = 0.5 * (p[e[:, 0], 1] + p[e[:, 1], 1])
    L = np.hypot(dlon * np.cos(np.deg2rad(lat)),
                 p[e[:, 0], 1] - p[e[:, 1], 1]) / DEG
    mid = 0.5 * (p[e[:, 0]] + p[e[:, 1]])
    return e, L, mid


print("[tbcmp] bay-window band stats (lon 139.55-140.2):", flush=True)
for tag, (p, t) in meshes.items():
    e, L, mid = edge_stats(p, t)
    inbay = (mid[:, 0] > 139.55) & (mid[:, 0] < 140.2)
    row = []
    for lo, hi in [(34.95, 35.10), (35.10, 35.20), (35.20, 35.30),
                   (35.30, 35.40), (35.40, 35.55), (35.55, 35.70)]:
        s = inbay & (mid[:, 1] >= lo) & (mid[:, 1] < hi)
        row.append(f"{np.percentile(L[s], 50):.0f}" if s.sum() else "-")
    print(f"[tbcmp]  {tag:13s} p50 by lat band: {row}", flush=True)

# figures: full domain + bay zoom, side by side, atlas grid
windows = [
    ("full", (115.0, 168.0, 10.0, 52.0)),
    ("kanto", (138.5, 141.5, 34.0, 36.2)),
    ("bay", (139.55, 140.2, 34.95, 35.75)),
    ("baynorth", (139.7, 140.15, 35.45, 35.72)),
]
for wname, wb in windows:
    fig, axs = plt.subplots(1, 2, figsize=(17, 9))
    for ax, (tag, (p, t)) in zip(axs, meshes.items()):
        cen = p[t].mean(axis=1)
        m = ((cen[:, 0] > wb[0]) & (cen[:, 0] < wb[1])
             & (cen[:, 1] > wb[2]) & (cen[:, 1] < wb[3]))
        lw = 0.08 if wname == "full" else 0.2
        ax.triplot(p[:, 0], p[:, 1], t[m], lw=lw, color="k")
        ax.set_xlim(wb[0], wb[1]); ax.set_ylim(wb[2], wb[3])
        ax.set_aspect(1 / np.cos(np.deg2rad(0.5 * (wb[2] + wb[3]))))
        # atlas grid A-T / 1-20
        xe = np.linspace(wb[0], wb[1], 21)
        ye = np.linspace(wb[2], wb[3], 21)
        for k in range(21):
            ax.axvline(xe[k], color="tab:blue", lw=0.25, alpha=0.5)
            ax.axhline(ye[k], color="tab:blue", lw=0.25, alpha=0.5)
        for k in range(20):
            ax.text(0.5 * (xe[k] + xe[k + 1]), wb[2] -
                    0.012 * (wb[3] - wb[2]), chr(65 + k), ha="center",
                    va="top", fontsize=5, color="tab:blue")
            ax.text(wb[0] - 0.012 * (wb[1] - wb[0]),
                    0.5 * (ye[k] + ye[k + 1]), str(k + 1), ha="right",
                    va="center", fontsize=5, color="tab:blue")
        _, _, _p = tag, wb, p
        ax.set_title(f"{tag}", fontsize=10)
    fig.suptitle(f"Tokyo Bay varres-3r translation — {wname}",
                 fontsize=11)
    fig.savefig(OUT / f"tbcmp_{wname}.png", dpi=180,
                bbox_inches="tight")
    plt.close(fig)
    print(f"[tbcmp] saved tbcmp_{wname}.png", flush=True)
print("[tbcmp] done", flush=True)
