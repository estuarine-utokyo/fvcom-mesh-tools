# Ex6b: compare our two-stage meshes against the MATLAB goldens
# (matlab_ex6b.mat: puw/tuw stage-1, pc/tc stage-2, pf/ef constraints).
# Produces: NP/res-ratio stats, node-density deficit tiles (stage 1),
# and side-by-side atlas figures for both stages.
import os, sys
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
import h5py
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DEG = 1.0 / 111e3
OUT = Path("outputs/om2d_examples/ex6bfp")
bbox = (-95.40, -94.4, 29.14, 30.09)

with h5py.File(OUT/"matlab_ex6b.mat", "r") as f:
    g_puw = np.array(f["puw"]).T; g_tuw = np.array(f["tuw"], dtype=int).T - 1
    g_pc = np.array(f["pc"]).T;  g_tc = np.array(f["tc"], dtype=int).T - 1
    g_pf = np.array(f["pf"]).T
puw = np.load(OUT/"puw.npy"); tuw = np.load(OUT/"tuw.npy")
pc = np.load(OUT/"p.npy");   tc = np.load(OUT/"t.npy")
pf = np.load(OUT/"pfix.npy")


def nodal_res(p, t):
    e = np.vstack([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]])
    e = np.unique(np.sort(e, axis=1), axis=0)
    L = np.linalg.norm(p[e[:, 0]] - p[e[:, 1]], axis=1) / DEG
    r = np.full(len(p), np.inf)
    np.minimum.at(r, e[:, 0], L)
    np.minimum.at(r, e[:, 1], L)
    return r


for tag, (op, ot, gp, gt) in {
        "stage1": (puw, tuw, g_puw, g_tuw),
        "stage2": (pc, tc, g_pc, g_tc)}.items():
    print(f"[cmp6b] {tag}: NP ours={len(op):,} golden={len(gp):,} "
          f"({100*(len(op)-len(gp))/len(gp):+.1f}%)", flush=True)
    r_o = nodal_res(op, ot); r_g = nodal_res(gp, gt)
    _, nn = cKDTree(op).query(gp)
    ratio = r_o[nn] / r_g
    print(f"[cmp6b] {tag}: res ratio p10/50/90 = "
          f"{np.percentile(ratio, [10, 50, 90]).round(3)}", flush=True)
    for lo, hi, bt in [(0, 100, "<100m"), (100, 300, "100-300"),
                       (300, 700, "300-700"), (700, 1e9, ">700")]:
        s = (r_g >= lo) & (r_g < hi)
        if s.sum() > 50:
            print(f"[cmp6b] {tag} golden {bt:8s}: n={s.sum():7,} "
                  f"ratio p50={np.percentile(ratio[s], 50):.3f}", flush=True)

# node-density deficit tiles for stage 1 (where do we lack nodes?)
NB = 20
xe = np.linspace(bbox[0], bbox[1], NB + 1)
ye = np.linspace(bbox[2], bbox[3], NB + 1)
Ho, *_ = np.histogram2d(puw[:, 0], puw[:, 1], bins=[xe, ye])
Hg, *_ = np.histogram2d(g_puw[:, 0], g_puw[:, 1], bins=[xe, ye])
diff = Ho - Hg
worst = np.dstack(np.unravel_index(np.argsort(diff.ravel())[:10],
                                   diff.shape))[0]
print("[cmp6b] stage1 largest node deficits (tile col A-T, row 1-20):",
      flush=True)
for i, j in worst:
    print(f"[cmp6b]   {chr(65+i)}{j+1}: ours {int(Ho[i,j]):5d} "
          f"golden {int(Hg[i,j]):5d} diff {int(diff[i,j]):+6d}", flush=True)

# atlas side-by-side figures
def atlas(ax, p, t, title):
    ax.triplot(p[:, 0], p[:, 1], t, lw=0.08, color="k")
    ax.set_xlim(bbox[0], bbox[1]); ax.set_ylim(bbox[2], bbox[3])
    ax.set_aspect(1/np.cos(np.deg2rad(29.6)))
    for k in range(NB + 1):
        ax.axvline(xe[k], color="tab:blue", lw=0.25, alpha=0.5)
        ax.axhline(ye[k], color="tab:blue", lw=0.25, alpha=0.5)
    for k in range(NB):
        ax.text(0.5*(xe[k]+xe[k+1]), ye[0]-0.012, chr(65+k),
                ha="center", va="top", fontsize=5, color="tab:blue")
        ax.text(xe[0]-0.012, 0.5*(ye[k]+ye[k+1]), str(k+1),
                ha="right", va="center", fontsize=5, color="tab:blue")
    ax.set_title(title, fontsize=9)


for tag, op, ot, gp, gt in [("stage1", puw, tuw, g_puw, g_tuw),
                            ("stage2", pc, tc, g_pc, g_tc)]:
    fig, axs = plt.subplots(1, 2, figsize=(17, 9))
    atlas(axs[0], op, ot, f"ours {tag} NP={len(op):,}")
    atlas(axs[1], gp, gt, f"OM2D golden {tag} NP={len(gp):,}")
    fig.savefig(OUT/f"cmp6b_{tag}_full.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[cmp6b] saved cmp6b_{tag}_full.png", flush=True)

# constraint overlay: our pfix vs golden pf
fig, ax = plt.subplots(figsize=(10, 9))
ax.plot(g_pf[:, 0], g_pf[:, 1], ".", ms=0.8, color="k", label=f"golden pfix {len(g_pf):,}")
ax.plot(pf[:, 0], pf[:, 1], ".", ms=0.8, color="crimson", alpha=0.6,
        label=f"ours pfix {len(pf):,}")
ax.set_xlim(bbox[0], bbox[1]); ax.set_ylim(bbox[2], bbox[3])
ax.set_aspect(1/np.cos(np.deg2rad(29.6)))
ax.legend(markerscale=12, fontsize=8)
ax.set_title("Ex6b fixed constraints: ours vs golden")
fig.savefig(OUT/"cmp6b_pfix.png", dpi=180, bbox_inches="tight")
print("[cmp6b] saved cmp6b_pfix.png", flush=True)

# zoomed panels: bay mouth + floodplain transition
for tag, zb in [("zoom_baymouth", (-94.85, -94.6, 29.25, 29.45)),
                ("zoom_floodplain", (-95.2, -94.95, 29.55, 29.80))]:
    fig, axs = plt.subplots(1, 2, figsize=(16, 8))
    for ax, (p_, t_, ti) in zip(
            axs, [(pc, tc, f"ours stage2"), (g_pc, g_tc, "OM2D golden stage2")]):
        ax.triplot(p_[:, 0], p_[:, 1], t_, lw=0.15, color="k")
        ax.set_xlim(zb[0], zb[1]); ax.set_ylim(zb[2], zb[3])
        ax.set_aspect(1/np.cos(np.deg2rad(29.6)))
        ax.set_title(ti, fontsize=9)
    fig.savefig(OUT/f"cmp6b_{tag}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[cmp6b] saved cmp6b_{tag}.png", flush=True)
print("[cmp6b] done", flush=True)
