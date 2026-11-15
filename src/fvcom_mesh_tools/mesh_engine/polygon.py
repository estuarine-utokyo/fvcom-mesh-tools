"""Projected-CRS polygon engine: a region-supplied water polygon, a metre size
field and an explicit open-boundary line -> a raw FVCOM mesh.

The ``oceanmesh`` engine of ``fmesh-buildmesh`` meshes a DEM bounding box in
lon/lat and classifies its open boundary by bbox proximity, which only fits a
domain whose open sea is a bbox edge (Tokyo Bay). A new region usually brings
its own domain polygon (coast + an offshore arc + an inland closure), so this
engine meshes exactly that polygon, in metres, in the recipe CRS:

1. size field on a regular lattice in the CRS (:func:`build_size_field`): a
   size raster in metres and/or a background size, ``sizing``-schema regions
   (:mod:`fvcom_mesh_tools.sizing`; lon/lat polygons, priority, transition),
   a floor, a cap and a gradient limit;
2. the open-boundary line becomes CGAL-constrained edges with nodes at the
   local size, plus (``ladder_k > 0``) an inner guide line ``k`` x local size
   inside, so every open-boundary node gets a near-perpendicular interior edge;
3. oceanmesh (lab fork) DistMesh with a polygon signed-distance function,
   then the OM2D default clean with the constraints pinned;
4. components: those touching the open boundary are kept (the largest when
   there is no open boundary); a dropped component larger than
   ``max_dropped_component_km2`` stops the build instead of vanishing, unless
   most of it lies in ``expected_drop`` (water the region already knows is
   not connected); every component is reported (``in_keep``: it holds water
   of ``keep_polygons``). Small pools a mesh of this size cannot connect are
   dropped; a region checks the coverage of its water bodies afterwards;
5. boundary repair until stable: single-element nodes off the open boundary
   (ears), FVCOM ``{1,2,1}`` elements, zero-area boundary slivers, and R4
   corners at the open-string ends (the end moves one node inward).

Coordinates are shifted to a local origin while meshing (DistMesh's numerical
gradients lose precision at UTM magnitudes) and shifted back on output.

Ported from VR-FVCOM ``tools/mesh/build_mesh.py`` (2026-11-15/16, same owner).
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np

__all__ = ["build_size_field", "generate_polygon_mesh", "limgrad"]


def limgrad(h: np.ndarray, grade: float, dx: float, itmax: int = 5000):
    """Limit ``|grad h| <= grade`` on a regular grid of spacing ``dx``
    (eight-neighbour relaxation until no change). Returns (h, sweeps)."""
    h = np.array(h, dtype=float, copy=True)
    d1, d2 = grade * dx, grade * dx * np.sqrt(2.0)
    for it in range(itmax):
        h0 = h.copy()
        h[1:, :] = np.minimum(h[1:, :], h[:-1, :] + d1)
        h[:-1, :] = np.minimum(h[:-1, :], h[1:, :] + d1)
        h[:, 1:] = np.minimum(h[:, 1:], h[:, :-1] + d1)
        h[:, :-1] = np.minimum(h[:, :-1], h[:, 1:] + d1)
        h[1:, 1:] = np.minimum(h[1:, 1:], h[:-1, :-1] + d2)
        h[:-1, :-1] = np.minimum(h[:-1, :-1], h[1:, 1:] + d2)
        h[1:, :-1] = np.minimum(h[1:, :-1], h[:-1, 1:] + d2)
        h[:-1, 1:] = np.minimum(h[:-1, 1:], h[1:, :-1] + d2)
        if np.array_equal(h, h0):
            return h, it + 1
    return h, itmax


def _read_size_raster(path, crs_epsg: int):
    """(xs ascending, ys ascending, h[ny, nx]) of a single-band raster in
    metres whose CRS must be ``crs_epsg``."""
    import rasterio
    from pyproj import CRS

    with rasterio.open(path) as src:
        if src.crs is None or CRS.from_user_input(src.crs).to_epsg() != int(crs_epsg):
            raise ValueError(f"size raster {path}: CRS {src.crs} is not EPSG:{crs_epsg}")
        t = src.transform
        if abs(t.b) > 0 or abs(t.d) > 0 or abs(abs(t.a) - abs(t.e)) > 1e-9 * abs(t.a):
            raise ValueError(f"size raster {path}: needs a north-up grid with square cells")
        a = src.read(1).astype(float)
        if src.nodata is not None:
            a[a == src.nodata] = np.nan
        dx = abs(t.a)
        xs = t.c + dx * (np.arange(src.width) + 0.5)
        ys = t.f + t.e * (np.arange(src.height) + 0.5)
    if t.e < 0:  # north-up: flip to ascending y
        ys, a = ys[::-1], a[::-1]
    return xs, ys, a, dx


def build_size_field(
    domain,
    *,
    crs_epsg: int,
    raster=None,
    background_m: float | None = None,
    raster_m: float = 50.0,
    regions: list | None = None,
    hmin_m: float | None = None,
    max_m: float | None = None,
    gradation: float | None = None,
    log: Callable[[str], None] = print,
) -> dict[str, Any]:
    """Size field (target edge length, metres) on a lattice in EPSG ``crs_epsg``.

    ``domain`` is the shapely (Multi)Polygon to mesh, in that CRS. Either
    ``raster`` (a GeoTIFF of sizes in metres, same CRS; NaN/nodata cells take
    ``background_m`` or else ``max_m``) or ``background_m`` must be given.
    ``regions`` follow the :mod:`fvcom_mesh_tools.sizing` schema (lon/lat
    geometry, ``target_h_m``, ``transition_m``, ``priority``) and replace the
    field inside their footprint, weakest first. Then the floor ``hmin_m``,
    the cap ``max_m`` and the gradient limit ``gradation`` apply.

    Returns ``{"xs", "ys", "h", "dx", "report"}`` (ys ascending, h[ny, nx]).
    """
    if raster is None and background_m is None:
        raise ValueError("size field: give a size raster or background_m")
    if raster is not None:
        xs, ys, h, dx = _read_size_raster(raster, crs_epsg)
        fill = background_m if background_m is not None else max_m
        if np.isnan(h).any():
            if fill is None:
                raise ValueError("size raster has nodata cells: give background_m or max_m")
            h = np.where(np.isnan(h), float(fill), h)
        b = domain.bounds
        if xs[0] - dx / 2 > b[0] or xs[-1] + dx / 2 < b[2] \
                or ys[0] - dx / 2 > b[1] or ys[-1] + dx / 2 < b[3]:
            raise ValueError(f"size raster {raster} does not cover the domain {b}")
    else:
        dx = float(raster_m)
        b = domain.bounds
        pad = 10 * dx
        xs = np.arange(b[0] - pad, b[2] + pad, dx)
        ys = np.arange(b[1] - pad, b[3] + pad, dx)
        h = np.full((len(ys), len(xs)), float(background_m))
    if (h <= 0).any():
        raise ValueError("size field has non-positive sizes")
    report: dict[str, Any] = {"lattice": [int(len(xs)), int(len(ys))], "dx_m": float(dx),
                              "regions": []}
    if regions:
        import shapely
        from pyproj import Transformer
        from shapely.ops import transform as shp_transform

        from fvcom_mesh_tools.sizing import _geometry, _regions

        _regions(regions)
        fwd = Transformer.from_crs("EPSG:4326", f"EPSG:{crs_epsg}", always_xy=True)
        ordered = sorted(regions, key=lambda r: (r.get("priority", 0), -r["target_h_m"],
                                                 r["name"]))
        for r in ordered:
            poly = shp_transform(fwd.transform, _geometry(r["geometry"]))
            tr = float(r.get("transition_m", 0) or 0)
            x0, y0, x1, y1 = poly.buffer(tr).bounds
            ix = np.nonzero((xs >= x0) & (xs <= x1))[0]
            iy = np.nonzero((ys >= y0) & (ys <= y1))[0]
            if not len(ix) or not len(iy):
                report["regions"].append({"name": r["name"], "cells": 0})
                continue
            gx, gy = np.meshgrid(xs[ix], ys[iy])
            pts = shapely.points(gx, gy)
            if tr > 0:
                w = np.maximum(0.0, 1.0 - shapely.distance(pts, poly) / tr)
            else:
                w = shapely.covers(poly, pts).astype(float)
            sub = h[np.ix_(iy, ix)]
            h[np.ix_(iy, ix)] = np.where(w > 0, sub * (1 - w) + r["target_h_m"] * w, sub)
            report["regions"].append({"name": r["name"], "cells": int((w >= 1).sum()),
                                      "target_h_m": float(r["target_h_m"])})
    if hmin_m is not None:
        h = np.maximum(h, float(hmin_m))
    if max_m is not None:
        h = np.minimum(h, float(max_m))
    if gradation is not None:
        h, nit = limgrad(h, float(gradation), dx)
        report["limgrad_sweeps"] = int(nit)
    log(f"[polygon] size lattice {len(xs)}x{len(ys)} at {dx:g} m: "
        f"h {h.min():.0f}-{h.max():.0f} m")
    return {"xs": xs, "ys": ys, "h": h, "dx": float(dx), "report": report}


def _face_components(t):
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    e = np.vstack([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]])
    e.sort(axis=1)
    fid = np.tile(np.arange(len(t)), 3)
    key = e[:, 0].astype(np.int64) * (t.max() + 1) + e[:, 1]
    order = np.argsort(key)
    k, f = key[order], fid[order]
    same = np.nonzero(k[1:] == k[:-1])[0]
    a, b = f[same], f[same + 1]
    adj = coo_matrix((np.ones(len(a)), (a, b)), shape=(len(t), len(t)))
    return connected_components(adj, directed=False)


def _tri_area(p, t):
    a, b, c = p[t[:, 0]], p[t[:, 1]], p[t[:, 2]]
    return 0.5 * ((b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1])
                  - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0]))


def generate_polygon_mesh(
    water,
    size: dict[str, Any],
    *,
    obc_line=None,
    ladder_k: float = 0.0,
    dm_scale: float = 1.0,
    max_iter: int = 60,
    seed: int = 0,
    clean_dj_cutoff: float = 0.25,
    max_dropped_component_km2: float = 0.5,
    keep_polygons=None,
    expected_drop=None,
    end_trim_nodes: int = 6,
    log: Callable[[str], None] = print,
) -> dict[str, Any]:
    """Mesh the shapely (Multi)Polygon ``water`` (CRS metres).

    ``size`` is :func:`build_size_field`'s result; DistMesh is driven with
    ``h / dm_scale`` (DistMesh returns edges ~``dm_scale`` x its field).
    ``obc_line`` (shapely LineString on the domain boundary, same CRS) gives
    the open boundary; without it the mesh has land boundaries only.
    ``keep_polygons`` (shapely geometry) is reported per component;
    ``expected_drop`` (shapely geometry): a large component lying mostly
    inside it is an expected drop (logged), not a failure.

    Returns ``{"points", "cells", "open_string", "land_string", "islands",
    "components", "n_ladder"}`` (points in the CRS, 0-based ids).
    """
    import oceanmesh as om
    import shapely
    from oceanmesh.clean import om2d_default_clean
    from scipy.interpolate import RegularGridInterpolator
    from scipy.spatial import cKDTree
    from shapely.affinity import translate

    xs, ys, h = size["xs"], size["ys"], size["h"]
    x0, y0 = np.floor(np.array(water.bounds[:2]) / 1000.0) * 1000.0 - 5000.0

    def sh(g):
        return translate(g, -x0, -y0)

    water_l = sh(water)
    obc_l = sh(obc_line) if obc_line is not None else None
    keep_l = sh(keep_polygons) if keep_polygons is not None and not keep_polygons.is_empty \
        else None
    exp_l = sh(expected_drop) if expected_drop is not None and not expected_drop.is_empty \
        else None
    hmax = float(h.max())
    hreal = RegularGridInterpolator((ys - y0, xs - x0), h, bounds_error=False,
                                    fill_value=hmax)

    def hat(q):
        q = np.atleast_2d(q)
        return hreal(np.column_stack([q[:, 1], q[:, 0]]))

    def fh(q):
        return hat(q) / dm_scale

    hmin_field = float(h.min()) / dm_scale
    ds = hmin_field / 8.0
    rings = [r for q in getattr(water_l, "geoms", [water_l])
             for r in [q.exterior] + list(q.interiors)]
    bpts = np.vstack([shapely.get_coordinates(shapely.segmentize(r, ds)) for r in rings])
    tree = cKDTree(bpts)
    shapely.prepare(water_l)

    def fd(q):
        q = np.asarray(q, float)
        d, _ = tree.query(q, workers=-1)
        inside = shapely.contains_xy(water_l, q[:, 0], q[:, 1])
        return np.where(inside, -d, d)

    log(f"[polygon] signed distance: {len(bpts)} boundary samples at {ds:.1f} m, "
        f"{len(rings)} rings; local origin ({x0:.0f}, {y0:.0f})")

    pfix_l, egfix_l, arc, n_inner = [], [], None, 0
    if obc_l is not None:
        L = obc_l.length
        ss = np.linspace(0, L, 4001)
        pts = np.array([obc_l.interpolate(s).coords[0] for s in ss])
        hs = hat(pts)
        cum = np.concatenate([[0], np.cumsum(0.5 * (1 / hs[1:] + 1 / hs[:-1]) * np.diff(ss))])
        nseg = max(int(round(cum[-1])), 3)
        sk = np.interp(np.linspace(0, cum[-1], nseg + 1), cum, ss)
        arc = np.array([obc_l.interpolate(s).coords[0] for s in sk])
        h_arc = hat(arc)
        pfix_l.append(arc)
        egfix_l.append(np.column_stack([np.arange(nseg), np.arange(1, nseg + 1)]))
        if ladder_k > 0:
            tg = np.gradient(arc, axis=0)
            tg /= np.linalg.norm(tg, axis=1)[:, None]
            nrm = np.column_stack([-tg[:, 1], tg[:, 0]])
            mid = len(arc) // 2
            if not water_l.contains(shapely.Point(arc[mid] + nrm[mid] * h_arc[mid])):
                nrm = -nrm
            inner = arc + nrm * (ladder_k * h_arc)[:, None]
            coast = water_l.boundary.difference(obc_l.buffer(1.0))
            chord = shapely.LineString(arc)
            # an inner point near the chord between two open-boundary nodes made
            # zero-area elements (VR v1 build 7589989): keep >= 0.5 h from it
            ok = np.array([water_l.contains(shapely.Point(q))
                           and coast.distance(shapely.Point(q)) >= 1.0 * hh
                           and chord.distance(shapely.Point(q)) >= 0.5 * hh
                           for q, hh in zip(inner, h_arc)])
            ok[[0, -1]] = False
            idx = np.nonzero(ok)[0]
            keep = [idx[0]] if len(idx) else []
            for i in idx[1:]:
                if np.hypot(*(inner[i] - inner[keep[-1]])) >= 0.5 * h_arc[i]:
                    keep.append(i)
            keep = np.array(keep, int)
            if len(keep) >= 2:
                base = len(arc)
                pfix_l.append(inner[keep])
                adj = np.nonzero(np.diff(keep) == 1)[0]
                egfix_l.append(np.column_stack([adj, adj + 1]) + base)
                n_inner = len(keep)
        log(f"[polygon] open boundary: {len(arc)} nodes on {L / 1e3:.1f} km; "
            f"ladder k={ladder_k:g}: {n_inner} guide nodes")
    PFIX = np.vstack(pfix_l) if pfix_l else None
    SEGS = np.vstack(egfix_l) if egfix_l else None

    b = water_l.bounds
    bbox = (float(b[0]), float(b[2]), float(b[1]), float(b[3]))
    kw = dict(bbox=bbox, min_edge_length=hmin_field, max_iter=int(max_iter), seed=int(seed),
              cleanup="none")
    if PFIX is not None:
        kw.update(pfix=PFIX, egfix=SEGS)
    p, t = om.generate_mesh(fd, fh, **kw)
    log(f"[polygon] generated NP={len(p)} NE={len(t)}")

    def comp_rows(stage, p, t, keep_c=None):
        ncomp, lab = _face_components(t)
        ar = np.abs(_tri_area(p, t))
        cen = p[t].mean(axis=1)
        rows = []
        for c in range(ncomp):
            sel = lab == c
            cxy = cen[sel].mean(axis=0) + [x0, y0]
            in_keep = bool(keep_l is not None and shapely.contains_xy(
                keep_l, cen[sel, 0], cen[sel, 1]).any())
            expected = bool(exp_l is not None and shapely.contains_xy(
                exp_l, cen[sel, 0], cen[sel, 1]).mean() > 0.5)
            rows.append(dict(stage=stage, comp=c, n_elem=int(sel.sum()),
                             area_km2=float(ar[sel].sum() / 1e6),
                             x=float(cxy[0]), y=float(cxy[1]), in_keep=in_keep,
                             expected_drop=expected,
                             kept=None if keep_c is None else bool(c in keep_c)))
        return ncomp, lab, rows

    def anchor_components(p, t, lab):
        if arc is None:
            sizes = np.bincount(lab)
            return {int(sizes.argmax())}
        on_arc = cKDTree(p).query(arc)[1]
        return set(lab[np.isin(t, on_arc).any(axis=1)].tolist())

    _n, lab, rows = comp_rows("pre-clean", p, t)
    keep_c = anchor_components(p, t, lab)
    for r in rows:
        if r["comp"] not in keep_c and r["expected_drop"] \
                and r["area_km2"] > max_dropped_component_km2:
            log(f"[polygon] expected drop: {r['n_elem']} el / {r['area_km2']:.2f} km2 at "
                f"({r['x']:.0f}, {r['y']:.0f})")
    lost = [r for r in rows if r["comp"] not in keep_c and not r["expected_drop"]
            and r["area_km2"] > max_dropped_component_km2]
    if lost:
        raise RuntimeError(
            "generator components that would be dropped (no open boundary) larger than "
            f"{max_dropped_component_km2} km2: "
            + "; ".join(f"{r['n_elem']} el / {r['area_km2']:.2f} km2 at "
                        f"({r['x']:.0f}, {r['y']:.0f})" for r in lost))
    if PFIX is not None:
        pairs = np.stack([PFIX[SEGS[:, 0]], PFIX[SEGS[:, 1]]], axis=1)
        p, t = om2d_default_clean(p, t, pfix=PFIX, egfix_pairs=pairs,
                                  dj_cutoff=float(clean_dj_cutoff))
    else:
        p, t = om2d_default_clean(p, t, dj_cutoff=float(clean_dj_cutoff))
    log(f"[polygon] after OM2D default clean NP={len(p)} NE={len(t)}")
    _n, lab, rows2 = comp_rows("post-clean", p, t)
    keep_c = anchor_components(p, t, lab)
    for r in rows2:
        r["kept"] = r["comp"] in keep_c
    lost = [r for r in rows2 if not r["kept"] and not r["expected_drop"]
            and r["area_km2"] > max_dropped_component_km2]
    if lost:
        raise RuntimeError(
            "components lost in the OM2D clean (no open boundary) larger than "
            f"{max_dropped_component_km2} km2: "
            + "; ".join(f"{r['n_elem']} el / {r['area_km2']:.2f} km2 at "
                        f"({r['x']:.0f}, {r['y']:.0f})" for r in lost))
    rows += rows2
    dropped = [r for r in rows2 if not r["kept"]]
    if dropped:
        log(f"[polygon] dropped {len(dropped)} component(s) without open boundary, "
            f"{sum(r['area_km2'] for r in dropped):.2f} km2 "
            f"({sum(r['in_keep'] for r in dropped)} holding keep_polygons water)")
    t = t[np.isin(lab, list(keep_c))]
    used = np.unique(t)
    remap = -np.ones(len(p), int)
    remap[used] = np.arange(len(used))
    p, t = p[used], remap[t]
    neg = _tri_area(p, t) < 0
    t[neg] = t[neg][:, [0, 2, 1]]

    def loop_area(lp):
        x, y = p[lp, 0], p[lp, 1]
        return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))

    open_str = np.empty(0, int)
    if arc is not None:
        d_arc, arc_idx = cKDTree(p).query(arc)
        if (d_arc > 1e-6).any():
            raise RuntimeError(f"open-boundary nodes lost in the clean: "
                               f"{np.nonzero(d_arc > 1e-6)[0].tolist()[:10]}")
        ends = [int(arc_idx[0]), int(arc_idx[-1])]
    n_ears = n_fake = n_trim = 0
    for it in range(50):
        cnt = np.bincount(t.ravel(), minlength=len(p))
        loops = om.boundary_loops(t)
        if arc is not None:
            outer = np.asarray(max(loops, key=lambda lp: abs(loop_area(lp))))
            if ends[0] not in set(outer.tolist()) or ends[1] not in set(outer.tolist()):
                raise RuntimeError("an open-boundary end node is not on the outer boundary loop")
            ring = np.roll(outer, -int(np.nonzero(outer == ends[0])[0][0]))
            stop = int(np.nonzero(ring == ends[1])[0][0])
            if np.any(shapely.distance(shapely.points(p[ring[1:stop]]), obc_l) > 1.0):
                ring = np.roll(ring[::-1], 1)
                stop = int(np.nonzero(ring == ends[1])[0][0])
            open_str = ring[:stop + 1]
            off = shapely.distance(shapely.points(p[open_str]), obc_l)
            if (off > 1.0).any():
                raise RuntimeError("the outer loop between the open-boundary ends leaves the "
                                   f"line (max offset {off.max():.1f} m)")
        lone = cnt == 1
        lone[open_str] = False
        ear = lone[t].any(axis=1)
        e3 = np.sort(np.stack([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]], axis=1), axis=2)
        code = e3[..., 0].astype(np.int64) * len(p) + e3[..., 1]
        uc, cc = np.unique(code.ravel(), return_counts=True)
        bcode = uc[cc == 1]
        bnd_nodes = np.unique(np.concatenate([bcode // len(p), bcode % len(p)]))
        isonb = np.zeros(len(p), int)
        isonb[bnd_nodes] = 1
        isonb[open_str] = 2
        sm = isonb[t].sum(axis=1)
        if len(open_str) > 1:
            oc = np.sort(np.column_stack([open_str[:-1], open_str[1:]]), axis=1)
            ocode = oc[:, 0].astype(np.int64) * len(p) + oc[:, 1]
            has_open = np.isin(code, ocode).any(axis=1)
        else:
            has_open = np.zeros(len(t), bool)
        fake = (sm == 4) & ~has_open
        # sum(ISONB) > 4 WITHOUT an open edge (two non-consecutive open nodes and a
        # boundary node): deleting it cannot break the open string (VR v3 job 7590918)
        fake |= (sm > 4) & ~has_open
        r4 = (sm > 4) & has_open
        flat = (np.abs(_tri_area(p, t)) < 1.0) & np.isin(code, bcode).any(axis=1)
        ear |= flat
        if not (ear.any() or fake.any() or r4.any()):
            break
        if r4.any():
            # an R4 element near an end (land closing in on the line within
            # ``end_trim_nodes`` nodes, VR v3 job 7590924): the open string ends
            # beyond it, the trimmed nodes become land
            pos = {int(v): i for i, v in enumerate(open_str)}
            nlast = len(open_str) - 1
            i0, i1 = 0, nlast
            for el in np.nonzero(r4)[0]:
                ip = [pos[int(v)] for v in t[el] if int(v) in pos]
                if max(ip) < end_trim_nodes:
                    i0 = max(i0, max(ip) + 1)
                elif min(ip) > nlast - end_trim_nodes:
                    i1 = min(i1, min(ip) - 1)
            moved = (i0, i1) != (0, nlast)
            if moved:
                if i1 - i0 < 2:
                    raise RuntimeError("open string trimmed to nothing at the R4 ends")
                ends = [int(open_str[i0]), int(open_str[i1])]
                n_trim += 1
            bn = set(t[r4].ravel().tolist()) if not moved else {int(open_str[0])}
            if not ({open_str[0], open_str[-1]} & bn):
                where = p[t[r4]].mean(axis=1)[:5] + [x0, y0]
                raise RuntimeError("R4 violation away from the open-boundary ends at "
                                   + "; ".join(f"({a:.0f}, {b:.0f})" for a, b in where))
        n_ears += int(ear.sum())
        n_fake += int((fake & ~ear).sum())
        t = t[~(ear | fake)]
    else:
        raise RuntimeError("boundary repair did not converge in 50 passes")
    log(f"[polygon] boundary repair ({it} passes): {n_ears} ear / flat elements and {n_fake} "
        f"{{1,2,1}} elements deleted; open-string ends moved inward {n_trim} time(s)")
    used = np.unique(t)
    remap = -np.ones(len(p), int)
    remap[used] = np.arange(len(used))
    p, t = p[used], remap[t]
    loops = [np.asarray(lp) for lp in om.boundary_loops(t)]
    outer = max(loops, key=lambda lp: abs(loop_area(lp)))
    if arc is not None:
        open_str = remap[open_str]
        ring = np.roll(outer, -int(np.nonzero(outer == open_str[0])[0][0]))
        if ring[1] != open_str[1]:
            ring = np.roll(ring[::-1], 1)
        stop = int(np.nonzero(ring == open_str[-1])[0][0])
        if not np.array_equal(ring[:stop + 1], open_str):
            raise RuntimeError("open string changed during compaction")
        land_str = np.append(ring[stop:], ring[0])
    else:
        land_str = np.append(outer, outer[0])
    islands = [lp for lp in loops if not np.array_equal(lp, outer)]
    ta = _tri_area(p, t)
    if (np.abs(ta) < 1.0).any():
        raise RuntimeError(f"{int((np.abs(ta) < 1.0).sum())} degenerate elements (|area| < 1 m2)")
    log(f"[polygon] FINAL NP={len(p)} NE={len(t)}; open {len(open_str)} nodes, "
        f"{len(islands)} islands")
    return {"points": p + [x0, y0], "cells": t, "open_string": open_str,
            "land_string": land_str, "islands": islands, "components": rows,
            "n_ladder": int(n_inner)}
