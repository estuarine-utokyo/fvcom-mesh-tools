# Connectivity comparator (owner mandate 2026-07-12): the sample's
# water connectivity is the reference truth for Tokyo Bay.
#   SEVERANCE: sample-meshed water we failed to mesh, clustered;
#     clusters whose removal disconnects the SAMPLE dual graph are
#     severed passages -> CRITICAL.
#   BREACH: our elements sitting on ORIGINAL land (widen artifacts),
#     clustered; clusters whose removal disconnects OUR dual graph
#     are fabricated passages through land barriers -> CRITICAL.
# Exit code 1 on any critical finding.
import os
import sys
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
import geopandas as gpd
import shapely
from pyproj import Transformer
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from shapely.ops import unary_union
from shapely.strtree import STRtree
from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.plotting import (
    _add_coast,
    add_atlas_grid,
    use_readable_style,
)

use_readable_style()

tr = Transformer.from_crs("EPSG:32654", "EPSG:4326", always_xy=True)
G = os.path.expanduser('~/Github/TB-FVCOM/input/goto2023/grid/')
gd = open(G + 'TokyoBay_grd.dat').read().split('\n')
nn = int(gd[0].split('=')[1]); ne = int(gd[1].split('=')[1])
Ts = np.array([[int(w) for w in gd[2 + i].split()[1:4]]
               for i in range(ne)]) - 1
Ps = np.array([[float(w) for w in gd[2 + ne + i].split()[1:3]]
               for i in range(nn)])
lon_s, lat_s = tr.transform(Ps[:, 0], Ps[:, 1])
Pll_s = np.column_stack([lon_s, lat_s])

m = read_fort14("outputs/sample_repro/sample_repro_final.14")
lon_o, lat_o = tr.transform(m.nodes[:, 0], m.nodes[:, 1])
Pll_o = np.column_stack([lon_o, lat_o])
T_o = m.elements


def dual_pairs(T, n_nodes):
    keys = []
    for a, b in ((0, 1), (1, 2), (2, 0)):
        lo = np.minimum(T[:, a], T[:, b]).astype(np.int64)
        hi = np.maximum(T[:, a], T[:, b]).astype(np.int64)
        keys.append(lo * n_nodes + hi)
    keys = np.concatenate(keys)
    eids = np.tile(np.arange(len(T)), 3)
    o = np.argsort(keys, kind="stable")
    keys, eids = keys[o], eids[o]
    same = keys[1:] == keys[:-1]
    return np.column_stack([eids[:-1][same], eids[1:][same]])


def components_without(pairs, n, drop):
    keep = np.ones(n, bool)
    keep[list(drop)] = False
    pr = pairs[keep[pairs[:, 0]] & keep[pairs[:, 1]]]
    g = coo_matrix((np.ones(len(pr)), (pr[:, 0], pr[:, 1])),
                   shape=(n, n))
    ncomp, lab = connected_components(g + g.T, directed=False)
    return ncomp - int((~keep).sum()), lab


# ---------- SEVERANCE: sample water not covered by our mesh ------
our_polys = [shapely.Polygon(Pll_o[t]) for t in T_o]
tree = STRtree(our_polys)
cent_s = Pll_s[Ts].mean(axis=1)
pts = shapely.points(cent_s[:, 0], cent_s[:, 1])
covered = np.zeros(len(Ts), bool)
for i, pt in enumerate(pts):
    for j in tree.query(pt):
        if our_polys[j].covers(pt):
            covered[i] = True
            break
missing = np.where(~covered)[0]
print(f"[conn] sample elements NOT covered by our mesh: "
      f"{len(missing)} / {len(Ts)}", flush=True)

pairs_s = dual_pairs(Ts, nn)
n0, _ = components_without(pairs_s, len(Ts), [])
severed = []
if len(missing):
    # cluster missing elements (face adjacency)
    sub = set(missing.tolist())
    prs = pairs_s[np.isin(pairs_s, missing).all(axis=1)]
    g = coo_matrix((np.ones(len(prs)),
                    (prs[:, 0], prs[:, 1])),
                   shape=(len(Ts), len(Ts)))
    _, lab = connected_components(g + g.T, directed=False)
    # adjacency lists of the sample graph for the DETOUR test
    from collections import defaultdict, deque
    adj = defaultdict(list)
    for a, b in pairs_s:
        adj[int(a)].append(int(b))
        adj[int(b)].append(int(a))
    mset = set(int(i) for i in missing)

    def detour_cut(cl):
        """A cluster on a THROUGH waterway does not disconnect
        the sample graph when removed (the loop closes around),
        so the component test is blind to it. Instead: if the
        cluster's covered neighbours end up farther apart than
        max(20, 8 x cluster size) hops in the graph WITHOUT the
        cluster, the waterway is functionally cut."""
        clset = set(cl)
        nbrs = sorted({j for i in cl for j in adj[i]
                       if j not in mset})
        if len(nbrs) < 2:
            return False
        limit = max(20, 8 * len(cl))
        src = nbrs[0]
        dist = {src: 0}
        q = deque([src])
        while q:
            v = q.popleft()
            if dist[v] >= limit:
                continue
            for w2 in adj[v]:
                if w2 in clset or w2 in dist:
                    continue
                dist[w2] = dist[v] + 1
                q.append(w2)
        return any(n2 not in dist for n2 in nbrs[1:])

    for c in set(int(lab[i]) for i in missing):
        cl = [i for i in missing if lab[i] == c]
        nafter, _ = components_without(pairs_s, len(Ts), cl)
        cc = cent_s[cl].mean(axis=0)
        cut = nafter > n0 or detour_cut(cl)
        rec = (len(cl), float(cc[0]), float(cc[1]), cut)
        if cut:
            severed.append(rec)
        print(f"[conn]   missing x{rec[0]} at ({rec[1]:.4f}, "
              f"{rec[2]:.4f}) severs_sample={rec[3]}", flush=True)

# ---------- KEPT-NETWORK THROUGH PROBE (owner 2026-07-16) --------
# The F9-c4 lesson: the coverage+detour tests above MISSED a real
# junction severance -- our wide cells covered part of the sample
# pocket (shrinking the missing cluster to 4) and the canal loop
# around the island closed within the max(20, 8n) hop bound. The
# direct test: for every KEPT waterway network, chain its carved
# branch arcs (bridging endpoint gaps <= ~3h) and require the
# corridor-restricted dual graph of OUR mesh to connect the two
# farthest arc ends whenever the SAMPLE's corridor does.
from scipy.spatial import cKDTree as _KDT

_SXm = 111e3 * float(np.cos(np.deg2rad(35.35)))
_SYm = 111e3


def _dense_chain(_arcs, bridge_max_m):
    segs = []
    for _a in _arcs:
        _a = np.asarray(_a, float)
        for _i in range(len(_a) - 1):
            segs.append((_a[_i], _a[_i + 1]))
    endpts = [np.asarray(_a, float)[k]
              for _a in _arcs for k in (0, -1)]
    for _i in range(len(endpts)):
        for _j in range(_i + 1, len(endpts)):
            _g = np.hypot((endpts[_i][0] - endpts[_j][0]) * _SXm,
                          (endpts[_i][1] - endpts[_j][1]) * _SYm)
            if 0 < _g <= bridge_max_m:
                segs.append((endpts[_i], endpts[_j]))
    pts = []
    for _p, _q in segs:
        _n = max(2, int(np.hypot((_q[0] - _p[0]) * _SXm,
                                 (_q[1] - _p[1]) * _SYm) / 50))
        pts.append(np.linspace(_p, _q, _n))
    return np.vstack(pts), endpts


def _corridor_probe(_P, _T, dense, _A, _B, rad=450.0):
    """True/False = corridor connects the probe ends; None = not
    probeable (no cells at an end)."""
    from collections import defaultdict as _dd, deque as _dq
    _cent = _P[_T].mean(axis=1)
    _tree = _KDT(np.column_stack([dense[:, 0] * _SXm,
                                  dense[:, 1] * _SYm]))
    _d, _ = _tree.query(
        np.column_stack([_cent[:, 0] * _SXm, _cent[:, 1] * _SYm]),
        distance_upper_bound=rad)
    _sel = np.where(np.isfinite(_d))[0]
    if len(_sel) == 0:
        return None
    _edge = _dd(list)
    for _j in _sel:
        _t = _T[_j]
        for _k in range(3):
            _e = tuple(sorted((int(_t[_k]), int(_t[(_k + 1) % 3]))))
            _edge[_e].append(int(_j))
    _adj = _dd(list)
    for _cells in _edge.values():
        if len(_cells) == 2:
            _adj[_cells[0]].append(_cells[1])
            _adj[_cells[1]].append(_cells[0])

    def _near(pt):
        _dd2 = np.hypot((_cent[_sel, 0] - pt[0]) * _SXm,
                        (_cent[_sel, 1] - pt[1]) * _SYm)
        _k2 = int(np.argmin(_dd2))
        return int(_sel[_k2]), float(_dd2[_k2])

    _a, _da = _near(_A)
    _b, _db = _near(_B)
    if max(_da, _db) > 500.0:
        return None
    _seen = {_a}
    _q = _dq([_a])
    while _q:
        _u = _q.popleft()
        if _u == _b:
            return True
        for _v in _adj[_u]:
            if _v not in _seen:
                _seen.add(_v)
                _q.append(_v)
    return False


severed_mesh = []
_wpath = "outputs/sample_repro/waterways.json"
if os.path.exists(_wpath):
    import json as _json2
    for _r in _json2.loads(open(_wpath).read()):
        if _r.get("action") != "keep" or not _r.get("arcs_done"):
            continue
        _arcs = [_a for _a, _w in _r["arcs_done"]]
        _dense, _endpts = _dense_chain(_arcs, bridge_max_m=1050.0)
        _best = None
        for _i in range(len(_endpts)):
            for _j in range(_i + 1, len(_endpts)):
                _g = np.hypot(
                    (_endpts[_i][0] - _endpts[_j][0]) * _SXm,
                    (_endpts[_i][1] - _endpts[_j][1]) * _SYm)
                if _best is None or _g > _best[0]:
                    _best = (_g, _endpts[_i], _endpts[_j])
        if _best is None or _best[0] < 700.0:
            continue
        _rs = _corridor_probe(Pll_s, Ts, _dense, _best[1], _best[2])
        _ro = _corridor_probe(Pll_o, T_o, _dense, _best[1], _best[2])
        if _rs is True and _ro is False:
            severed_mesh.append(_r["center"])
            print(f"[conn] CRITICAL mesh-severed kept network at "
                  f"({_r['center'][0]:.4f}, {_r['center'][1]:.4f})"
                  f": sample corridor CONNECTED, ours "
                  f"DISCONNECTED", flush=True)
        elif _rs is True and _ro is None:
            print(f"[conn]   note: kept network at "
                  f"({_r['center'][0]:.4f}, {_r['center'][1]:.4f})"
                  f" not probeable on our mesh (no cells at an "
                  f"arc end)", flush=True)
print(f"[conn] kept-network mesh probes: {len(severed_mesh)} "
      f"CRITICAL mesh-severed", flush=True)

# ---------- BREACH: our elements on ORIGINAL land ----------------
land = unary_union(list(gpd.read_file(
    "outputs/tb_varres_3r/land_osm_wide.shp").geometry))
# Applied channel-arc edits are owner-approved design geometry
# (e.g. the pier-supported D-runway drawn as land by OSM): re-carve
# them out of the ORIGINAL land before the breach test, with the
# exact same guarded operation the runner used.
import json as _json
from pathlib import Path as _Path
from shapely.geometry import Polygon as _Poly
from fvcom_mesh_tools.channel_arcs import carve_channel_corridor
_OBC = [[139.6713, 35.1396], [139.6737, 35.1288],
        [139.6772, 35.1168], [139.6816, 35.1031],
        [139.6871, 35.0877], [139.6946, 35.0705],
        [139.7000, 35.0576], [139.7069, 35.0445],
        [139.7134, 35.0327], [139.7216, 35.0184],
        [139.7289, 35.0047], [139.7373, 34.9916],
        [139.7497, 34.9750]]
_dom = _Poly([[139.83, 34.973], [140.12, 34.973], [140.12, 35.75],
              [139.60, 35.75], [139.60, 35.20],
              [139.6642, 35.1546]] + _OBC + [[139.83, 34.973]])
_cosw = float(np.cos(np.deg2rad(35.35)))
# EDIT COMPONENT-CONNECTION GATE (owner 2026-07-15: edit_005
# rev 2 punched through the Urayasu harbour wall and merged the
# ENCLOSED basin with the sea, and no gate saw it -- the breach
# reference trusts applied edits by construction). Rule: an edit
# corridor/water_patch that touches TWO OR MORE distinct
# connected components of the ORIGINAL water is CRITICAL unless
# it declares "allow_connect": true with its evidence.
_ow_parts = [g for g in getattr(
    _dom.difference(land), "geoms",
    [_dom.difference(land)]) if not g.is_empty]
from shapely.strtree import STRtree as _ST3
_ow_tree = _ST3(_ow_parts)
_conn_gate_fail = 0


def _edit_component_check(_geom, _eid, _allow):
    global _conn_gate_fail
    _hit = set()
    for _k in _ow_tree.query(_geom, predicate="intersects"):
        _pp = _ow_parts[int(_k)]
        _ix = _geom.intersection(_pp)
        if _ix.area * (111e3 * _cosw) * 111e3 > 110.0 * 110.0:
            _hit.add(int(_k))
    if len(_hit) >= 2 and not _allow:
        _conn_gate_fail += 1
        print(f"[conn] CRITICAL edit-connect: {_eid} touches "
              f"{len(_hit)} separate ORIGINAL water components "
              f"without allow_connect", flush=True)


_SR_EXCL = {s.strip() for s in os.environ.get(
    "SR_EDITS_EXCLUDE", "").split(",") if s.strip()}
for _ef in sorted(
        _Path("recipes/edits/sample_repro").glob("*.json")):
    if _ef.stem in _SR_EXCL:
        print(f"[conn] edit {_ef.stem}: EXCLUDED "
              f"(SR_EDITS_EXCLUDE)", flush=True)
        continue
    _ed = _json.loads(_ef.read_text())
    if _ed.get("type") == "land_patch":
        import shapely.geometry as _sg
        land = unary_union([land, _sg.shape(_ed["geometry"])])
        print(f"[conn] breach reference: applied land_patch "
              f"{_ed.get('id', _ef.stem)}", flush=True)
        continue
    if _ed.get("type") == "water_patch":
        import shapely.geometry as _sg
        _wp = _sg.shape(_ed["geometry"])
        _edit_component_check(_wp, _ed.get("id", _ef.stem),
                              _ed.get("allow_connect", False))
        land = land.difference(_wp)
        print(f"[conn] breach reference: applied water_patch "
              f"{_ed.get('id', _ef.stem)}", flush=True)
        continue
    _tol = _ed.get("arc_on_land_tol_m")
    _w = (np.asarray(_ed["widths_m"], float)
          if "widths_m" in _ed else float(_ed["width_m"]))
    import shapely as _shp
    _wmax = float(np.max(_w))
    _corr = _shp.LineString(
        np.asarray(_ed["arc"], float)).buffer(
        0.5 * _wmax / 111e3)
    _edit_component_check(_corr, _ed.get("id", _ef.stem),
                          _ed.get("allow_connect", False))
    land, _ei = carve_channel_corridor(
        land, np.asarray(_ed["arc"], float), _w,
        min_gap_m=float(_ed.get("min_gap_m", 150.0)),
        metric_scale=(111e3 * _cosw, 111e3), domain_poly=_dom,
        arc_on_land_tol_m=None if _tol is None else float(_tol))
    print(f"[conn] breach reference: carved applied edit "
          f"{_ed.get('id', _ef.stem)}", flush=True)
cent_o = Pll_o[T_o].mean(axis=1)
opts = shapely.points(cent_o[:, 0], cent_o[:, 1])
# tolerance: centroid deeper than ~100 m into original land
land_shrunk = land.buffer(-100.0 / 111e3)
onland = np.array([land_shrunk.covers(p) for p in opts])
breach_idx = np.where(onland)[0]
# split INTENDED (the pipeline itself carved water there --
# waterway widening / manual edits; runner saves its adjusted
# land as land_channel_adj.shp) from UNINTENDED erosion
_pipe = unary_union(list(gpd.read_file(
    "outputs/sample_repro/land_channel_adj.shp").geometry))
unint = np.array([_pipe.covers(opts[i]) for i in breach_idx])
# choke widen-ops (finish stage, AFTER land_channel_adj): pushed
# banks are INTENDED widening, ledgered in widen_ops.json
_wops_f = _Path("outputs/sample_repro/widen_ops.json")
if _wops_f.exists() and len(breach_idx):
    _wops = _json.loads(_wops_f.read_text())
    if _wops:
        from pyproj import Transformer as _Tr
        _tr2 = _Tr.from_crs(32654, 4326, always_xy=True)
        _wc = []
        for _op in _wops:
            for _o, _n in zip(_op["old"], _op["new"]):
                _ox, _oy = _tr2.transform(*_o)
                _nx, _ny = _tr2.transform(*_n)
                _wc.append(shapely.LineString(
                    [(_ox, _oy), (_nx, _ny)]).buffer(
                    0.7 * _op["h_loc"] / 111e3))
        _wz = unary_union(_wc)
        for _k, _i in enumerate(breach_idx):
            if unint[_k] and _wz.covers(opts[_i]):
                unint[_k] = False
# Element identities for the issue map (386), 0-indexed into our mesh.
_Path("outputs/sample_repro/land_breaches.json").write_text(_json.dumps({
    "elements": [int(i) for i in breach_idx],
    "unintended": [bool(u) for u in unint],
}, indent=1) + "\n")
print(f"[conn] our elements on ORIGINAL land: {len(breach_idx)} "
      f"(intended widening: {int((~unint).sum())}, UNINTENDED: "
      f"{int(unint.sum())})", flush=True)
pairs_o = dual_pairs(T_o, m.n_nodes)
m0, _ = components_without(pairs_o, len(T_o), [])
breaches = []
if len(breach_idx):
    prs = pairs_o[np.isin(pairs_o, breach_idx).all(axis=1)]
    g = coo_matrix((np.ones(len(prs)),
                    (prs[:, 0], prs[:, 1])),
                   shape=(len(T_o), len(T_o)))
    _, lab = connected_components(g + g.T, directed=False)
    for c in set(int(lab[i]) for i in breach_idx):
        cl = [i for i in breach_idx if lab[i] == c]
        nafter, _ = components_without(pairs_o, len(T_o), cl)
        cc = cent_o[cl].mean(axis=0)
        rec = (len(cl), float(cc[0]), float(cc[1]), nafter > m0)
        if nafter > m0:
            breaches.append(rec)
        print(f"[conn]   on-land x{rec[0]} at ({rec[1]:.4f}, "
              f"{rec[2]:.4f}) bridges={rec[3]}", flush=True)

print(f"[conn] CRITICAL edit-connects:    {_conn_gate_fail}",
      flush=True)
print(f"[conn] CRITICAL severed passages: {len(severed)}",
      flush=True)
print(f"[conn] CRITICAL land breaches:    {len(breaches)}",
      flush=True)
print(f"[conn] CRITICAL mesh-severed nets: {len(severed_mesh)}",
      flush=True)

# map figure -- markers must NOT hide the mesh (owner
# 2026-07-12): hollow, translucent markers; legend outside the
# axes so no text overlaps mesh or grid labels
fig, ax = plt.subplots(figsize=(10.5, 12))
_add_coast(ax, (139.0, 34.5, 141.3, 36.2), "EPSG:4326")
ax.triplot(lon_o, lat_o, T_o, lw=0.3, color="steelblue",
           zorder=3)
if len(missing):
    ax.scatter(cent_s[missing, 0], cent_s[missing, 1], s=42,
               facecolors="none", edgecolors="darkorange",
               linewidths=1.0, alpha=0.75, zorder=4,
               label=f"sample water we closed ({len(missing)})")
for nrec, x, y, crit in severed:
    ax.plot([x], [y], marker="x", ms=16, mew=3, color="red",
            alpha=0.9, zorder=6)
for _x, _y in severed_mesh:
    ax.plot([_x], [_y], marker="X", ms=16, mew=2.4, mec="red",
            mfc="none", alpha=0.9, zorder=6)
if len(breach_idx):
    ax.scatter(cent_o[breach_idx, 0], cent_o[breach_idx, 1],
               s=55, marker="s", facecolors="none",
               edgecolors="purple", linewidths=1.2, alpha=0.75,
               zorder=4,
               label=f"our mesh on land ({len(breach_idx)})")
ax.set_xlim(139.57, 140.15); ax.set_ylim(34.93, 35.78)
add_atlas_grid(ax, crs="EPSG:4326")
ax.set_aspect(1 / np.cos(np.deg2rad(35.35)))
ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0),
          borderaxespad=0.0, frameon=True)
ax.set_title("connectivity comparison vs sample + original land\n"
             "red X = severed sample passage (critical)")
fig.savefig("outputs/figures/connectivity_check.png", dpi=190,
            bbox_inches="tight")
print("[conn] saved map", flush=True)
sys.exit(1 if (severed or breaches or severed_mesh) else 0)
