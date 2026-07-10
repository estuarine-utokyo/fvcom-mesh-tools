# Ex7 Global: compare our mesh (lonlat p.npy/t.npy) against the
# MATLAB golden (matlab_ex7.mat pc/tc). Stats + atlas-gridded
# figures: global resolution maps, res histogram, coastal zooms.
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
OUT = Path("outputs/om2d_examples/ex7glob")
bbox = (-180.0, 180.0, -89.0, 90.0)

with h5py.File(OUT/"matlab_ex7.mat", "r") as f:
    gp = np.array(f["pc"]).T
    gt = np.array(f["tc"], dtype=int).T - 1
p = np.load(OUT/"p.npy")
t = np.load(OUT/"t.npy")
print(f"[cmp7] NP ours={len(p):,} golden={len(gp):,} "
      f"({100*(len(p)-len(gp))/len(gp):+.1f}%)", flush=True)
print(f"[cmp7] NT ours={len(t):,} golden={len(gt):,}", flush=True)


def nodal_res_metric(pp, tt):
    """Min incident PHYSICAL edge length (km) per node."""
    e = np.vstack([tt[:, [0, 1]], tt[:, [1, 2]], tt[:, [2, 0]]])
    e = np.unique(np.sort(e, axis=1), axis=0)
    dlon = np.abs(pp[e[:, 0], 0] - pp[e[:, 1], 0])
    dlon = np.minimum(dlon, 360 - dlon)  # dateline
    lat = 0.5 * (pp[e[:, 0], 1] + pp[e[:, 1], 1])
    dx = dlon * np.cos(np.deg2rad(lat))
    dy = pp[e[:, 0], 1] - pp[e[:, 1], 1]
    L = np.hypot(dx, dy) * 111.0  # km
    r = np.full(len(pp), np.inf)
    np.minimum.at(r, e[:, 0], L)
    np.minimum.at(r, e[:, 1], L)
    return r


r_o = nodal_res_metric(p, t)
r_g = nodal_res_metric(gp, gt)
_, nn = cKDTree(p).query(gp, workers=-1)
ratio = r_o[nn] / r_g
print(f"[cmp7] res ratio p10/50/90 = "
      f"{np.percentile(ratio, [10, 50, 90]).round(3)}", flush=True)
for lo, hi, tag in [(0, 5, "<5km"), (5, 10, "5-10km"),
                    (10, 19, "10-19km"), (19, 1e9, ">=19km")]:
    s = (r_g >= lo) & (r_g < hi)
    if s.sum() > 100:
        print(f"[cmp7] golden {tag:8s}: n={s.sum():9,} "
              f"ratio p50={np.percentile(ratio[s], 50):.3f}", flush=True)

# global resolution maps (node scatter coloured by log resolution)
NBx, NBy = 20, 20
xe = np.linspace(bbox[0], bbox[1], NBx + 1)
ye = np.linspace(bbox[2], bbox[3], NBy + 1)


def resmap(ax, pp, rr, title):
    o = np.argsort(-rr)  # draw fine (small) last
    sc = ax.scatter(pp[o, 0], pp[o, 1], c=np.log10(np.clip(rr[o], 1, 50)),
                    s=0.05, cmap="turbo", vmin=0, vmax=np.log10(30))
    ax.set_xlim(bbox[0], bbox[1]); ax.set_ylim(bbox[2], bbox[3])
    for k in range(NBx + 1):
        ax.axvline(xe[k], color="w", lw=0.2, alpha=0.4)
    for k in range(NBy + 1):
        ax.axhline(ye[k], color="w", lw=0.2, alpha=0.4)
    for k in range(NBx):
        ax.text(0.5*(xe[k]+xe[k+1]), ye[0]-4, chr(65+k), ha="center",
                va="top", fontsize=6, color="tab:blue")
    for k in range(NBy):
        ax.text(xe[0]-4, 0.5*(ye[k]+ye[k+1]), str(k+1), ha="right",
                va="center", fontsize=6, color="tab:blue")
    ax.set_title(title, fontsize=10)
    return sc


fig, axs = plt.subplots(2, 1, figsize=(16, 15))
sc = resmap(axs[0], p, r_o, f"ours NP={len(p):,} (log10 min-edge km)")
resmap(axs[1], gp, r_g, f"OM2D golden NP={len(gp):,}")
fig.colorbar(sc, ax=axs, shrink=0.5, label="log10 edge length (km)")
fig.savefig(OUT/"cmp7_global_res.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("[cmp7] saved cmp7_global_res.png", flush=True)

# resolution histograms
fig, ax = plt.subplots(figsize=(9, 5))
bins = np.logspace(np.log10(2), np.log10(60), 60)
ax.hist(r_o[np.isfinite(r_o)], bins=bins, histtype="step",
        color="crimson", label=f"ours ({len(p):,})")
ax.hist(r_g[np.isfinite(r_g)], bins=bins, histtype="step",
        color="k", label=f"golden ({len(gp):,})")
ax.set_xscale("log"); ax.set_xlabel("nodal min edge (km)")
ax.set_ylabel("nodes"); ax.legend()
ax.set_title("Ex7 global: nodal resolution distribution")
fig.savefig(OUT/"cmp7_res_hist.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("[cmp7] saved cmp7_res_hist.png", flush=True)

# coastal zooms with actual triangles
zooms = [("nz", (166, 180, -48, -34)),
         ("carib", (-70, -58, 10, 20)),
         ("japan", (128, 146, 30, 46))]
for tag, zb in zooms:
    fig, axs = plt.subplots(1, 2, figsize=(16, 8))
    for ax, (pp, tt, ti) in zip(
            axs, [(p, t, f"ours"), (gp, gt, "OM2D golden")]):
        cen_x = pp[tt].mean(axis=1)
        m = ((cen_x[:, 0] > zb[0]) & (cen_x[:, 0] < zb[1]) &
             (cen_x[:, 1] > zb[2]) & (cen_x[:, 1] < zb[3]))
        ax.triplot(pp[:, 0], pp[:, 1], tt[m], lw=0.12, color="k")
        ax.set_xlim(zb[0], zb[1]); ax.set_ylim(zb[2], zb[3])
        ax.set_aspect(1/np.cos(np.deg2rad(0.5*(zb[2]+zb[3]))))
        ax.set_title(f"{ti} — {tag}", fontsize=10)
    fig.savefig(OUT/f"cmp7_zoom_{tag}.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"[cmp7] saved cmp7_zoom_{tag}.png", flush=True)
print("[cmp7] done", flush=True)
