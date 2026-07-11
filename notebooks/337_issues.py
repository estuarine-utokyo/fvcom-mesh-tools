# Issue atlas for the final rebuilt mesh: EVERY remaining concern,
# zoomed, with the offending elements/edges highlighted.
import os
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
from pyproj import Transformer
from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.diagnostics import channel_width_metric
from fvcom_mesh_tools.plotting import _add_coast, add_atlas_grid
from fvcom_mesh_tools.algorithms.perp_local import _tri_quality

OUT = "outputs/figures"
m = read_fort14("outputs/sample_repro/sample_repro_final.14")
tr = Transformer.from_crs("EPSG:32654", "EPSG:4326", always_xy=True)
lon, lat = tr.transform(m.nodes[:, 0], m.nodes[:, 1])
P = np.column_stack([lon, lat])
T = m.elements
cent = P[T].mean(axis=1)

issues = []   # (category, elem_ids_or_None, center_ll, label)

# 1. w/h < 1 channel throats
ch = channel_width_metric(m, coords="metric")
r = ch["w_h_ratio"]
bad = np.where(np.isfinite(r) & (r < 1.0))[0]
print(f"[iss] w/h<1 elements: {len(bad)}", flush=True)
used = np.zeros(len(bad), bool)
for k in np.argsort(r[bad]):
    if used[k]:
        continue
    c = cent[bad[k]]
    dd = np.hypot((cent[bad][:, 0] - c[0]) * 91e3,
                  (cent[bad][:, 1] - c[1]) * 111e3)
    grp = dd < 2500
    used |= grp
    ids = bad[grp]
    issues.append(("w/h<1 throat", ids, c,
                   f"w/h min={r[ids].min():.2f} x{len(ids)}"))

# 2. worst Courant nodes (dt=18)
b = np.maximum(m.depths, 1.0)
e = np.vstack([T[:, [0, 1]], T[:, [1, 2]], T[:, [2, 0]]])
e = np.unique(np.sort(e, axis=1), axis=0)
L = np.linalg.norm(m.nodes[e[:, 0]] - m.nodes[e[:, 1]], axis=1)
dx = np.full(len(P), np.inf)
np.minimum.at(dx, e[:, 0], L)
np.minimum.at(dx, e[:, 1], L)
u = np.sqrt(9.81 * b) + np.sqrt(9.81 / b)
cr = 18.0 * u / dx
top = np.argsort(-cr)[:12]
usedn = set()
for v in top:
    if v in usedn or cr[v] < 1.0:
        continue
    dd = np.hypot((P[top, 0] - P[v, 0]) * 91e3,
                  (P[top, 1] - P[v, 1]) * 111e3)
    for w, d2 in zip(top, dd):
        if d2 < 2500:
            usedn.add(w)
    issues.append(("Cr>1 (dt=18)", None, P[v],
                   f"Cr={cr[v]:.2f} H={m.depths[v]:.0f}m dx={dx[v]:.0f}m"))
print(f"[iss] Cr>1 nodes: {(cr > 1.0).sum()}", flush=True)

# 3. worst-quality elements
q = _tri_quality(m.nodes[T])[0]
wq = np.argsort(q)[:3]
for ei in wq:
    issues.append(("min-angle low", np.array([ei]), cent[ei],
                   f"min angle {q[ei]:.1f} deg"))

print(f"[iss] total zoom sites: {len(issues)}", flush=True)
COAST = (139.0, 34.5, 141.3, 36.2)
n = len(issues)
ncol = 4
nrow = (n + ncol - 1) // ncol
fig, axes = plt.subplots(nrow, ncol, figsize=(22, 5.6 * nrow))
for ax in np.ravel(axes):
    ax.set_visible(False)
for ax, (catg, ids, c, lab) in zip(np.ravel(axes), issues):
    ax.set_visible(True)
    R = 0.018
    _add_coast(ax, COAST, "EPSG:4326")
    ax.triplot(lon, lat, T, lw=0.7, color="steelblue")
    if ids is not None:
        for ei in ids:
            tri = T[ei]
            ax.fill(P[tri, 0], P[tri, 1], color="crimson",
                    alpha=0.55, zorder=5)
    else:
        ax.plot([c[0]], [c[1]], marker="o", ms=13, mfc="none",
                mec="crimson", mew=2.2, zorder=6)
    ax.set_xlim(c[0] - R, c[0] + R)
    ax.set_ylim(c[1] - R * 0.82, c[1] + R * 0.82)
    add_atlas_grid(ax, crs="EPSG:4326")
    ax.set_aspect(1 / np.cos(np.deg2rad(c[1])))
    ax.set_title(f"{catg}: {lab}\n({c[0]:.3f}, {c[1]:.3f})",
                 fontsize=9)
fig.suptitle("remaining issues atlas — final mesh", y=0.995)
fig.savefig(f"{OUT}/issues_atlas.png", dpi=170, bbox_inches="tight")
print("saved issues atlas", flush=True)
