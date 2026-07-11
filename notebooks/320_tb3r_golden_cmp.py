# Tokyo Bay varres-3r: MATLAB golden vs faithful vs pre-faithful.
# Per-nest node counts, bay band stats, res-ratio at golden nodes,
# deficit tiles, and atlas side-by-side figures.
import os, sys
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
import h5py
from scipy.spatial import cKDTree
from fvcom_mesh_tools.io import read_fort14
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path as MPath

DEG = 1.0 / 111e3
OUT = Path("outputs/tb_varres_3r")

with h5py.File(OUT / "matlab_tb3r.mat", "r") as f:
    gp = np.array(f["p"]).T
    gt = np.array(f["t"], dtype=int).T - 1
new = read_fort14(str(OUT / "tb_varres_3regions.14"))
old = read_fort14(str(OUT / "tb_varres_3regions_prefaithful.14"))
meshes = {
    "matlab": (gp, gt),
    "faithful": (new.nodes[:, :2], new.elements),
    "pre-faithful": (old.nodes[:, :2], old.elements),
}

# nest polygons (from the .m); nest1 = everything else
x2 = [138.196907945361, 139.784782123012, 141.294666567004,
      139.718998032166, 138.285893837899, 138.196907945361]
y2 = [34.6003740729489, 33.5207179317026, 35.2698037998702,
      35.9357551145740, 35.2841759386230, 34.6003740729489]
p3_ring = None  # nest3 box approx: use its bbox rectangle
b3 = (139.142156, 140.392541, 34.745556, 35.859592)
nest2 = MPath(np.column_stack([x2, y2]))

for tag, (p, t) in meshes.items():
    inb3 = ((p[:, 0] >= b3[0]) & (p[:, 0] <= b3[1])
            & (p[:, 1] >= b3[2]) & (p[:, 1] <= b3[3]))
    in2 = nest2.contains_points(p) & ~inb3
    n3, n2 = int(inb3.sum()), int(in2.sum())
    n1 = len(p) - n2 - n3
    print(f"[tb3] {tag:13s}: NP={len(p):7,}  nest1={n1:7,} "
          f"nest2={n2:6,} nest3={n3:7,}", flush=True)


def edge_stats(p, t):
    e = np.vstack([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]])
    e = np.unique(np.sort(e, axis=1), axis=0)
    dlon = (p[e[:, 0], 0] - p[e[:, 1], 0])
    lat = 0.5 * (p[e[:, 0], 1] + p[e[:, 1], 1])
    L = np.hypot(dlon * np.cos(np.deg2rad(lat)),
                 p[e[:, 0], 1] - p[e[:, 1], 1]) / DEG
    mid = 0.5 * (p[e[:, 0]] + p[e[:, 1]])
    r = np.full(len(p), np.inf)
    np.minimum.at(r, e[:, 0], L)
    np.minimum.at(r, e[:, 1], L)
    return e, L, mid, r


print("[tb3] bay-window p50 by lat band (lon 139.55-140.2):",
      flush=True)
stats = {}
for tag, (p, t) in meshes.items():
    e, L, mid, r = edge_stats(p, t)
    stats[tag] = (p, t, r)
    inbay = (mid[:, 0] > 139.55) & (mid[:, 0] < 140.2)
    row = []
    for lo, hi in [(34.95, 35.10), (35.10, 35.20), (35.20, 35.30),
                   (35.30, 35.40), (35.40, 35.55), (35.55, 35.70)]:
        s = inbay & (mid[:, 1] >= lo) & (mid[:, 1] < hi)
        row.append(f"{np.percentile(L[s], 50):.0f}" if s.sum() else "-")
    print(f"[tb3]  {tag:13s}: {row}", flush=True)

# res ratio at golden nodes (faithful vs matlab)
gp_, gt_, r_g = stats["matlab"]
p_, t_, r_o = stats["faithful"]
_, nn = cKDTree(p_).query(gp_, workers=-1)
ratio = r_o[nn] / r_g
print(f"[tb3] res ratio (faithful/matlab) p10/50/90 = "
      f"{np.percentile(ratio, [10, 50, 90]).round(3)}", flush=True)
for lo, hi, bt in [(0, 200, "<200m"), (200, 500, "200-500"),
                   (500, 2000, "0.5-2k"), (2000, 1e9, ">2k")]:
    s = (r_g >= lo) & (r_g < hi)
    if s.sum() > 100:
        print(f"[tb3]  golden {bt:8s}: n={s.sum():7,} "
              f"ratio p50={np.percentile(ratio[s], 50):.3f}",
              flush=True)

# figures: kanto + bay side-by-side (matlab vs faithful)
for wname, wb in [("kanto", (138.2, 141.3, 33.6, 36.0)),
                  ("bay", (139.55, 140.2, 34.95, 35.75)),
                  ("shelf", (138.4, 140.0, 34.0, 35.1))]:
    fig, axs = plt.subplots(1, 2, figsize=(17, 9))
    for ax, tag in zip(axs, ("matlab", "faithful")):
        p, t = meshes[tag]
        cen = p[t].mean(axis=1)
        m = ((cen[:, 0] > wb[0]) & (cen[:, 0] < wb[1])
             & (cen[:, 1] > wb[2]) & (cen[:, 1] < wb[3]))
        ax.triplot(p[:, 0], p[:, 1], t[m], lw=0.12, color="k")
        ax.set_xlim(wb[0], wb[1]); ax.set_ylim(wb[2], wb[3])
        ax.set_aspect(1 / np.cos(np.deg2rad(0.5 * (wb[2] + wb[3]))))
        xe = np.linspace(wb[0], wb[1], 21)
        ye = np.linspace(wb[2], wb[3], 21)
        for k in range(21):
            ax.axvline(xe[k], color="tab:blue", lw=0.25, alpha=0.5)
            ax.axhline(ye[k], color="tab:blue", lw=0.25, alpha=0.5)
        for k in range(20):
            ax.text(0.5 * (xe[k] + xe[k + 1]),
                    wb[2] - 0.012 * (wb[3] - wb[2]), chr(65 + k),
                    ha="center", va="top", fontsize=5,
                    color="tab:blue")
            ax.text(wb[0] - 0.012 * (wb[1] - wb[0]),
                    0.5 * (ye[k] + ye[k + 1]), str(k + 1),
                    ha="right", va="center", fontsize=5,
                    color="tab:blue")
        ax.set_title(f"{tag} NP={len(p):,}", fontsize=10)
    fig.suptitle(f"TB varres-3r — {wname} (MATLAB vs faithful)",
                 fontsize=11)
    fig.savefig(OUT / f"tb3_{wname}.png", dpi=180,
                bbox_inches="tight")
    plt.close(fig)
    print(f"[tb3] saved tb3_{wname}.png", flush=True)
print("[tb3] done", flush=True)
