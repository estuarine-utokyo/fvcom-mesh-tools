# One-wide channel cell reporter (owner spec 2026-07-12): where a
# waterway is carried by cells spanning bank to bank (all three
# nodes on the land boundary), mark each such cell IN COLOUR WITH
# ITS CELL NUMBER so the next manual-editing round can address
# them by ID. This is a REPORT, not a fixer -- sub-cell-width
# channels meshed 1-wide are accepted by design (minimum mesh size
# is inviolable; banks are only widened by explicit arc edits).
import json
import os
from pathlib import Path

import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import geopandas as gpd
import matplotlib.pyplot as plt
import shapely
from pyproj import Transformer
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from shapely.geometry import LineString
from shapely.ops import unary_union
from shapely.strtree import STRtree

from fvcom_mesh_tools.gridref import TOKYO_BAY_GRID
from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.plotting import (
    add_atlas_grid,
    use_readable_style,
)

use_readable_style()
OUT = Path("outputs/sample_repro")
COSW = float(np.cos(np.deg2rad(35.35)))

m = read_fort14(str(OUT / "sample_repro_final.14"))
tr = Transformer.from_crs("EPSG:32654", "EPSG:4326", always_xy=True)
lon, lat = tr.transform(m.nodes[:, 0], m.nodes[:, 1])
P = np.column_stack([lon, lat])
T = m.elements

# boundary nodes = nodes of edges used by exactly one element
keys, eids = [], np.tile(np.arange(len(T)), 3)
for a, b in ((0, 1), (1, 2), (2, 0)):
    lo = np.minimum(T[:, a], T[:, b]).astype(np.int64)
    hi = np.maximum(T[:, a], T[:, b]).astype(np.int64)
    keys.append(lo * m.n_nodes + hi)
keys = np.concatenate(keys)
o = np.argsort(keys, kind="stable")
ks, es = keys[o], eids[o]
same = np.zeros(len(ks), bool)
same[1:] |= ks[1:] == ks[:-1]
same[:-1] |= ks[1:] == ks[:-1]
bkeys = ks[~same]
bnode = np.zeros(m.n_nodes, bool)
bnode[(bkeys // m.n_nodes)] = True
bnode[(bkeys % m.n_nodes)] = True

obc = np.zeros(m.n_nodes, bool)
for ob in m.open_boundaries:
    obc[np.asarray(ob)] = True      # already 0-based (io/fort14)

allb = bnode[T].all(axis=1)
touches_obc = obc[T].any(axis=1)

# boundary-EDGE count per cell: a true bank-to-bank cell has all
# 3 nodes on the shoreline but AT MOST ONE boundary edge (its two
# cross-channel edges are shared with neighbours). Cove-corner
# cells (3 nodes on the SAME bank) have 2-3 boundary edges and are
# NOT one-wide (owner 2026-07-12: only cell 3290 is truly 1-wide
# at the Haneda NW channel; the earlier looser flags overcounted).
bset = set(bkeys.tolist())
nbedge = np.zeros(len(T), int)
for a, b in ((0, 1), (1, 2), (2, 0)):
    lo = np.minimum(T[:, a], T[:, b]).astype(np.int64)
    hi = np.maximum(T[:, a], T[:, b]).astype(np.int64)
    k = lo * m.n_nodes + hi
    nbedge += np.fromiter((kk in bset for kk in k), bool,
                          len(T)).astype(int)
flag = np.where(allb & (nbedge <= 1) & ~touches_obc)[0]
print(f"[1wide] CONFIRMED bank-to-bank cells (all 3 nodes on "
      f"shoreline, <=1 boundary edge, no OBC): {len(flag)}",
      flush=True)

# CHOKE-EDGE detector (owner 2026-07-12: cells at Haneda
# I8-a2/I8-b2 have an interior EDGE spanning the whole channel --
# the third node is interior, so the all-3-boundary criterion
# misses them). An interior edge whose two endpoints are both
# shoreline nodes that are NOT near-neighbours ALONG the boundary
# is a bank-to-bank choke; both sharing cells are flagged.
from collections import defaultdict
badj = defaultdict(set)
for k in bkeys:
    a, b = int(k // m.n_nodes), int(k % m.n_nodes)
    badj[a].add(b)
    badj[b].add(a)


def _boundary_hops(a, b, max_hops=3):
    seen = {a}
    front = {a}
    for hop in range(1, max_hops + 1):
        front = set().union(*(badj[v] for v in front)) - seen
        if b in front:
            return hop
        seen |= front
    return None


ikeys = np.unique(ks[same])
choke = set()
edge_cells = defaultdict(list)
for kk, ee in zip(ks, eids[o]):
    edge_cells[int(kk)].append(int(ee))
for k in ikeys:
    a, b = int(k // m.n_nodes), int(k % m.n_nodes)
    if not (bnode[a] and bnode[b]) or obc[a] or obc[b]:
        continue
    if _boundary_hops(a, b) is None:      # opposite banks
        for e in edge_cells[int(k)]:
            choke.add(e)
choke -= set(np.where(touches_obc)[0].tolist())
print(f"[1wide] choke-edge cells (interior bank-to-bank edge): "
      f"{len(choke)}", flush=True)
if choke:
    flag = np.unique(np.concatenate(
        [flag, np.fromiter(choke, dtype=int)]))

# arc cross-section detector (owner 2026-07-12: the boundary-node
# criterion MISSED the 1-wide cells of the constrained Haneda
# band). Along each edit's medial arc, a station where a single
# cell carries >35 % of the water crossing = channel is one cell
# wide there. This measures width in the user's sense directly.
polys = [shapely.Polygon(P[t]) for t in T]
tree = STRtree(polys)
scale = 0.5 * (111e3 * COSW + 111e3)
xflag: set[int] = set()
arc_sources = []
for ef in sorted(Path("recipes/edits/sample_repro")
                 .glob("*.json")):
    ed = json.loads(ef.read_text())
    if "check_arcs" in ed:
        arc_sources += list(zip(ed["check_arcs"],
                                ed["check_widths_m"]))
    elif "arc" in ed and "widths_m" in ed:
        arc_sources.append((ed["arc"], ed["widths_m"]))
# sweep EVERY kept detected waterway too (owner 2026-07-12: the
# Haneda NW channel had no checker arc after edit rev 5)
wpath = OUT / "waterways.json"
if wpath.exists():
    for wrec in json.loads(wpath.read_text()):
        if wrec["action"] != "keep":
            continue
        if wrec.get("arcs_done"):
            # every carved branch arc, side branches included
            for a2, w2 in wrec["arcs_done"]:
                arc_sources.append((a2, w2))
        elif wrec["arc"]:
            arc_sources.append((wrec["arc"], wrec["widths_m"]))
print(f"[1wide] cross-section sweep over {len(arc_sources)} "
      f"arcs", flush=True)
for pairs in [arc_sources]:
    for a_raw, w_raw in pairs:
        a = np.asarray(a_raw, float)
        wprof = np.asarray(w_raw, float)
        line = LineString(a)
        n_st = max(int(line.length * scale / 100.0) + 1, 5)
        for s in np.linspace(0.0, 1.0, n_st):
            p = line.interpolate(s, normalized=True)
            pa = line.interpolate(max(s - 0.02, 0.0),
                                  normalized=True)
            pb = line.interpolate(min(s + 0.02, 1.0),
                                  normalized=True)
            tv = np.array([pb.x - pa.x, pb.y - pa.y])
            tv /= np.hypot(*tv) + 1e-15
            nv = np.array([-tv[1], tv[0]])
            wi = float(np.interp(s, np.linspace(0, 1, len(wprof)),
                                 wprof)) * 0.6 / scale
            sec = LineString([(p.x - nv[0] * wi,
                               p.y - nv[1] * wi),
                              (p.x + nv[0] * wi,
                               p.y + nv[1] * wi)])
            fr, tot = [], 0.0
            for j in tree.query(sec):
                seg = polys[j].intersection(sec)
                if not seg.is_empty and seg.length > 0:
                    fr.append((int(j), seg.length))
                    tot += seg.length
            if tot <= 0:
                continue
            sub = [j for j, ln in fr if ln / tot > 0.35]
            if len(sub) == 1 and not touches_obc[sub[0]]:
                xflag.add(sub[0])
# cross-section dominance alone is ADVISORY (tier 2): a dominant
# cell with an interior node is one of TWO rows, not a bank-to-
# bank cell (owner correction: 725/3291 etc. are not one-wide)
xonly = np.array(sorted(xflag - set(flag.tolist())), int)
print(f"[1wide] confirmed one-wide: {len(flag)}; advisory "
      f"narrow (section-dominant only): {len(xonly)}", flush=True)

# cluster flagged cells by node-sharing so nearby cells report as
# one site
lab = -np.ones(len(T), int)
if len(flag):
    nid = {}
    rows, cols = [], []
    for fi, e in enumerate(flag):
        for v in T[e]:
            nid.setdefault(v, []).append(fi)
    for lst in nid.values():
        for a2 in lst[1:]:
            rows.append(lst[0])
            cols.append(a2)
    g = coo_matrix((np.ones(len(rows)), (rows, cols)),
                   shape=(len(flag), len(flag)))
    _, cl = connected_components(g + g.T, directed=False)
    lab[flag] = cl

cent = P[T].mean(axis=1)

# TERMINAL vs CHOKE (owner 2026-07-13): the last cell of a CLOSED
# channel mouth inevitably touches both banks -- that is the
# design, not a defect (and wedge splits double-counted them).
# A site is a real flow CHOKE only if removing its cells cuts the
# local dual graph; otherwise it is a terminal mouth (grey tier).
_eq = ks[1:] == ks[:-1]
pairs_all = np.column_stack([es[:-1][_eq], es[1:][_eq]])
from collections import defaultdict as _dd, deque as _dq
adj_c = _dd(list)
for a2, b2 in pairs_all:
    adj_c[int(a2)].append(int(b2))
    adj_c[int(b2)].append(int(a2))


def _is_choke(els):
    drop = set(int(e) for e in els)
    nbrs = sorted({j for e in drop for j in adj_c[e]
                   if j not in drop})
    if len(nbrs) < 2:
        return False
    src = nbrs[0]
    seen = {src}
    q = _dq([src])
    hops = 0
    limit = max(30, 10 * len(drop))
    while q and hops < 200000:
        v = q.popleft()
        for w2 in adj_c[v]:
            if w2 in drop or w2 in seen:
                continue
            seen.add(w2)
            q.append(w2)
        hops += 1
        if len(seen) > limit + len(nbrs):
            break
    return any(n2 not in seen for n2 in nbrs[1:])


sites = []
terminals = []
for c in sorted(set(lab[flag])) if len(flag) else []:
    els = flag[lab[flag] == c]
    cc = cent[els].mean(axis=0)
    rec = {
        "cell_ids_fort14": (els + 1).tolist(),
        "center_lonlat": [round(float(cc[0]), 4),
                          round(float(cc[1]), 4)],
        "gridref": TOKYO_BAY_GRID.point_to_subcell(
            float(cc[0]), float(cc[1])),
    }
    if _is_choke(els):
        sites.append(rec)
    else:
        terminals.append(rec)

# STABLE site IDs across runs: a persistent registry maps each
# site number to a canonical location; a detected site within
# 500 m of a registered one inherits that number, new sites get
# the next free number, and numbers are never reused (a fixed
# site simply stops appearing in the ledger).
REG = OUT / "ow_registry.json"
reg = (json.loads(REG.read_text()) if REG.exists()
       else {"OW": [], "TM": []})


def _assign_ids(recs, kind):
    entries = reg.setdefault(kind, [])
    used = set()
    for rec in recs:
        cx, cy = rec["center_lonlat"]
        best, bd = None, 500.0
        for ent in entries:
            if ent["id"] in used:
                continue
            d = float(np.hypot(
                (ent["lonlat"][0] - cx) * 111e3 * 0.816,
                (ent["lonlat"][1] - cy) * 111e3))
            if d < bd:
                best, bd = ent, d
        if best is None:
            best = {"id": len(entries) + 1, "lonlat": [cx, cy]}
            entries.append(best)
        used.add(best["id"])
        rec["site"] = f"{kind}{best['id']:02d}"
    recs.sort(key=lambda r: r["site"])


_assign_ids(sites, "OW")
_assign_ids(terminals, "TM")
REG.write_text(json.dumps(reg, indent=1))
advisory = [{
    "cell_id_fort14": int(e + 1),
    "center_lonlat": [round(float(cent[e][0]), 4),
                      round(float(cent[e][1]), 4)],
    "gridref": TOKYO_BAY_GRID.point_to_subcell(
        float(cent[e][0]), float(cent[e][1])),
} for e in xonly]
(OUT / "one_wide_cells.json").write_text(json.dumps(
    {"confirmed_sites": sites, "terminal_mouths": terminals,
     "advisory_narrow": advisory}, indent=1))
print(f"[1wide] real chokes: {len(sites)} sites; terminal "
      f"mouths (closed-channel ends, by design): "
      f"{len(terminals)} sites", flush=True)
for s in sites:
    print(f"[1wide] {s['site']} [{s['gridref']}]: cells "
          f"{s['cell_ids_fort14']} at "
          f"({s['center_lonlat'][0]}, {s['center_lonlat'][1]})",
          flush=True)
for a2 in advisory:
    print(f"[1wide]   advisory narrow: cell "
          f"{a2['cell_id_fort14']} [{a2['gridref']}]", flush=True)

land = unary_union(list(gpd.read_file(
    "outputs/tb_varres_3r/land_osm_wide.shp").geometry))
gser = gpd.GeoSeries([land], crs="EPSG:4326")
show = sites          # every confirmed site gets a panel
if show:
    ncol = min(3, len(show))
    nrow = int(np.ceil(len(show) / ncol))
    fig, axes = plt.subplots(nrow, ncol,
                             figsize=(5.6 * ncol, 5.6 * nrow),
                             squeeze=False)
    axes = axes.ravel()
    for ax, s in zip(axes, show):
        cx, cy = s["center_lonlat"]
        half = 0.012
        gser.plot(ax=ax, color="0.9", edgecolor="0.6", lw=0.5,
                  zorder=1)
        els = np.asarray(s["cell_ids_fort14"]) - 1
        # light fill + strong OUTLINE so the mesh lines stay
        # visible through the highlight (owner 2026-07-12)
        for e in els:
            ring = P[np.append(T[e], T[e][0])]
            ax.fill(ring[:, 0], ring[:, 1], color="crimson",
                    alpha=0.18, zorder=2)
            ax.plot(ring[:, 0], ring[:, 1], color="crimson",
                    lw=1.6, zorder=4)
        for e in xonly:
            ring = P[np.append(T[e], T[e][0])]
            ax.plot(ring[:, 0], ring[:, 1], color="darkorange",
                    lw=1.2, ls="--", zorder=3.5)
            ax.annotate(str(e + 1), cent[e], ha="center",
                        va="center", fontsize=10,
                        color="darkorange")
        ax.triplot(lon, lat, T, lw=0.6, color="steelblue",
                   zorder=3)
        for e in els:
            ax.annotate(str(e + 1), cent[e],
                        ha="center", va="center",
                        fontsize=12, fontweight="bold",
                        color="darkred", zorder=5)
        ax.set_xlim(cx - half / COSW, cx + half / COSW)
        ax.set_ylim(cy - half, cy + half)
        ax.set_aspect(1 / COSW)
        ax.set_xticks([])
        ax.set_yticks([])
        add_atlas_grid(ax, crs="EPSG:4326")
        ax.set_title(f"{s['site']} [{s['gridref']}]  "
                     f"{len(els)} cell(s)  ({cx:.3f}, {cy:.3f})")
    for ax in axes[len(show):]:
        ax.axis("off")
    fig.suptitle("one-wide channel cells: red outline = CONFIRMED "
                 "bank-to-bank (all 3 nodes on shoreline);\n"
                 "orange dashed = advisory narrow (one cell "
                 "dominates the cross-section but 2 rows exist)")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig("outputs/figures/one_wide_cells.png", dpi=170,
                bbox_inches="tight")
    print("[1wide] saved outputs/figures/one_wide_cells.png",
          flush=True)
else:
    print("[1wide] no one-wide cells -- no figure", flush=True)
