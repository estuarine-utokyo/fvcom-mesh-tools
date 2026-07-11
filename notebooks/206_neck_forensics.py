# Forensics: what kills the C5/D5 lobe during cleanup?
# - compare raw-mesh element quality (ours vs OM2D Precleaned)
# - instrument the clean stages; tag every deleted element
# - symmetric test: our clean applied to THEIR raw
import os, sys, logging
import numpy as np
logging.basicConfig(level=logging.WARNING, stream=sys.stdout)
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.io import loadmat
from pyproj import Transformer
from oceanmesh.fix_mesh import fix_mesh, simp_qual
from oceanmesh.clean import (_external_topology,
                             make_mesh_boundaries_traversable,
                             om2d_default_clean)
from oceanmesh.mesh_improve import collapse_thin_triangles

OUT = Path("outputs/om2d_examples/ex2_ny")
LOBE = (-74.118, -74.082, 40.7625, 40.790)
NECK = (-74.104, -74.080, 40.7520, 40.7725)

_tr = Transformer.from_crs(
    "EPSG:4326", "+proj=tmerc +lon_0=-74.0 +lat_0=40.75 "
    "+ellps=WGS84 +units=m", always_xy=True)

def proj(p):
    x, y = _tr.transform(p[:, 0], p[:, 1])
    return np.column_stack([x, y])

def unproj(pp):
    x, y = _tr.transform(pp[:, 0], pp[:, 1], direction="INVERSE")
    return np.column_stack([x, y])

def win_mask(c, w):
    return ((c[:, 0] > w[0]) & (c[:, 0] < w[1])
            & (c[:, 1] > w[2]) & (c[:, 1] < w[3]))

def qstats(p, t, label):
    pp = proj(p)
    q = simp_qual(pp, t)
    c = p[t].mean(axis=1)
    be, _ = _external_topology(pp, t)
    bidx = np.unique(np.asarray(be, dtype=int).reshape(-1))
    btouch = np.isin(t, bidx).any(axis=1)
    for wname, w in (("lobe", LOBE), ("neck", NECK)):
        m = win_mask(c, w)
        if not m.any():
            print(f"[fx] {label} {wname}: 0 elements", flush=True)
            continue
        qm = q[m]
        nb = int((m & btouch).sum())
        nbad = int((m & btouch & (q < 0.25)).sum())
        print(f"[fx] {label} {wname}: NT={int(m.sum())} "
              f"q10/50={np.percentile(qm,10):.3f}/"
              f"{np.percentile(qm,50):.3f} min={qm.min():.3f} "
              f"bnd-elems={nb} bnd&q<0.25={nbad}", flush=True)

theirs = loadmat(os.path.expanduser(
    "~/Github/OceanMesh2D/Examples/Precleaned_grid.mat"))
tp = np.asarray(theirs["p"], dtype=float)
tt = np.asarray(theirs["t"], dtype=int) - 1
print(f"[fx] their raw NP={len(tp):,} NT={len(tt):,}", flush=True)
op = np.load(OUT / "p_raw.npy"); ot = np.load(OUT / "t_raw.npy")
print(f"[fx] our raw NP={len(op):,} NT={len(ot):,}", flush=True)
qstats(tp, tt, "THEIR-raw")
qstats(op, ot, "OUR-raw")

def traced_clean(p, t, label, max_passes=20):
    """om2d_default_clean pass structure with kill-stage tags
    (projected frame), no pfix, smooth=True."""
    from oceanmesh.mesh_improve import bound_connectivity
    from oceanmesh.smooth2d import smooth2d
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    pp = proj(p)
    kill = {}   # centroid-key -> stage string

    def keys(pp_, t_):
        c = unproj(pp_[t_].mean(axis=1) if t_.ndim == 2 else pp_)
        return [tuple(k) for k in np.round(c, 7).tolist()]

    def mark(pp_, t_, mask, stage):
        c = unproj(pp_[t_[mask]].mean(axis=1))
        for k in np.round(c, 7).tolist():
            kill.setdefault(tuple(k), stage)

    for pas in range(max_passes):
        nt_in = len(t)
        for it in range(25):
            be, _ = _external_topology(pp, t)
            bidx = np.unique(np.asarray(be, dtype=int).reshape(-1))
            btouch = np.isin(t, bidx).any(axis=1)
            q = simp_qual(pp, t)
            bad = btouch & (q < 0.25)
            if not bad.any():
                break
            mark(pp, t, bad, f"p{pas+1}-db{it+1}")
            t = t[~bad]
            pp, t, _ = fix_mesh(pp, t, delete_unused=True)
        before = {k: True for k in keys(pp, t)}
        pp2, t2 = collapse_thin_triangles(pp, t, min_qual=0.25)
        after = set(keys(pp2, t2))
        for k in before:
            if k not in after:
                kill.setdefault(k, f"p{pas+1}-collapse")
        pp, t = pp2, t2
        # component census before MMBT
        e = np.sort(np.concatenate(
            [t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]]), axis=1)
        eu, inv, cnt = np.unique(e, axis=0, return_inverse=True,
                                 return_counts=True)
        nt = len(t)
        pair = {}
        rows = []; cols = []
        eidx = inv.reshape(3, nt).T
        from collections import defaultdict
        owner = defaultdict(list)
        for ei_row, ti in zip(eidx, range(nt)):
            for ei in ei_row:
                owner[ei].append(ti)
        for ei, tis in owner.items():
            if len(tis) == 2:
                rows.append(tis[0]); cols.append(tis[1])
        A = coo_matrix((np.ones(len(rows)), (rows, cols)),
                       shape=(nt, nt))
        ncomp, lab = connected_components(A, directed=False)
        pre = {k: l for k, l in zip(keys(pp, t), lab)}
        sizes = np.bincount(lab)
        pp3, t3 = make_mesh_boundaries_traversable(
            pp, t, min_disconnected_area=0.25)
        after3 = set(keys(pp3, t3))
        gone_by_comp = defaultdict(int)
        for k, l in pre.items():
            if k not in after3:
                gone_by_comp[l] += 1
        for l, n in sorted(gone_by_comp.items(),
                           key=lambda x: -x[1])[:6]:
            if n > 20:
                print(f"[fx] {label} pass{pas+1} MMBT killed "
                      f"{n}/{sizes[l]} of component {l} "
                      f"(ncomp={ncomp})", flush=True)
        for k in pre:
            if k not in after3:
                kill.setdefault(k, f"p{pas+1}-mmbt")
        pp, t = pp3, t3
        try:
            pp, t = bound_connectivity(pp, t, max_valence=9)
        except Exception:
            pass
        pp, t = smooth2d(pp, t)
        q = simp_qual(pp, t)
        if q.min() >= 0.25 or len(t) == nt_in:
            break
    print(f"[fx] {label} cleaned NT={len(t):,} "
          f"passes={pas+1}", flush=True)
    return unproj(pp), t, kill

p2, t2, kill_ours = traced_clean(op, ot, "OUR-raw")
p3, t3, kill_theirs = traced_clean(tp, tt, "THEIR-raw")
c3 = p3[t3].mean(axis=1)
print(f"[fx] our-clean(THEIR-raw): lobe NT="
      f"{int(win_mask(c3, LOBE).sum())}", flush=True)
c2 = p2[t2].mean(axis=1)
print(f"[fx] our-clean(OUR-raw): lobe NT="
      f"{int(win_mask(c2, LOBE).sum())}", flush=True)

# stage-kill map figure over a wide window
STAGES = ["db1", "db2", "db3", "db", "collapse", "mmbt", "later"]
CMAP = {"db1": "tab:orange", "db2": "tab:red", "db3": "firebrick",
        "db": "darkred", "collapse": "tab:purple",
        "mmbt": "tab:blue", "later": "tab:gray"}
def stage_of(s):
    if s.startswith("p1-db"):
        n = int(s.split("db")[1])
        return "db1" if n == 1 else ("db2" if n == 2 else "db3"
                                     if n == 3 else "db")
    if s == "p1-collapse":
        return "collapse"
    if s == "p1-mmbt":
        return "mmbt"
    return "later"

fig, axes = plt.subplots(1, 2, figsize=(19, 9))
for ax, (rp, rt, kills, ttl) in zip(
    axes,
    [(op, ot, kill_ours, "OUR raw: deletion stages"),
     (tp, tt, kill_theirs, "THEIR raw: deletion stages")],
):
    c = rp[rt].mean(axis=1)
    ax.triplot(rp[:, 0], rp[:, 1], rt, lw=0.15, color="0.75")
    pts = {s: [] for s in STAGES}
    for k, s in kills.items():
        pts[stage_of(s)].append(k)
    for s in STAGES:
        if pts[s]:
            a = np.asarray(pts[s])
            ax.plot(a[:, 0], a[:, 1], ".", ms=2.5, color=CMAP[s],
                    label=f"{s} ({len(a):,})")
    ax.set_xlim(-74.16, -74.06); ax.set_ylim(40.73, 40.81)
    ax.set_aspect(1 / np.cos(np.deg2rad(40.77)))
    ax.legend(loc="upper left", fontsize=8, markerscale=4)
    ax.set_title(ttl)
fig.suptitle("Cleanup kill-stage map (C4/C5-D5 window)")
fig.tight_layout()
fig.savefig(OUT / "kill_stage_map.png", dpi=170,
            bbox_inches="tight")
print("[fx] saved", flush=True)
