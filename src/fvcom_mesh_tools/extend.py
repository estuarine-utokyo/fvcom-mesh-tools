"""Extend a finished base mesh outward: sizing, merge and the frozen contract.

The owner's decision (2026-09-29): the Tokyo Bay base mesh is adopted **as
it is** and a wider mesh is made by meshing only the sea outside it.  The
base's open boundary becomes an interior line: its nodes are fixed in the
outer generation (pfix/egfix), the outer mesh is merged onto them, and the
base's nodes and elements come through byte for byte.

Sizing (all in metres on a lon/lat lattice with metric ``x``, ``y``):

* the ambient field (coast distance, max edge) is limited to the gradation;
* the time-step floor ``dt * sqrt(g * depth) / Cr`` is raised into a field
  that is itself gradation-feasible (a graded dilation), and the ambient
  field may not go below it -- the aim is that the extension does not limit
  dt; where bands or the final depths undercut it, the build reports a
  warning (owner, 2026-10-01);
* along each constrained line (the base interface, the new open boundary)
  the size is **set** to the line's own spacing on a band, and graded away
  from it both upward and downward.

A maximum or minimum of gradation-feasible fields is gradation-feasible, so
the result needs no further limiting.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np

from fvcom_mesh_tools.io.fort14 import Fort14Mesh
from fvcom_mesh_tools.sizing import _limit

__all__ = [
    "band_field",
    "check_island_holes",
    "check_land_cover",
    "check_no_overlap",
    "land_segments",
    "compose_sizing",
    "graded_up",
    "merge_outer",
    "rfactor_smooth_free",
    "round_depths_inside",
    "trim_lone_corners",
    "verify_frozen_base",
]


#: Largest relative departure of a band from its target sizes that is accepted.
BAND_TOLERANCE = 0.05


def _grade(grade) -> float:
    """A finite, non-negative real gradation (NumPy scalars too, not a bool)
    as a float (review rounds 13 F9, 14 F4)."""
    # real numbers only: a complex one would lose its imaginary part to
    # float() (review round 15 F4)
    if isinstance(grade, (bool, np.bool_)) or not isinstance(
            grade, (int, float, np.integer, np.floating)):
        raise ValueError(f"the gradation must be a real number, not {grade!r}")
    g = float(grade)
    if not (np.isfinite(g) and g >= 0):
        raise ValueError(f"the gradation must be finite and non-negative, not {grade!r}")
    return g


def _lattice(values, x, y):
    """``values``, ``x`` and ``y`` as float arrays of one 2-D shape, with a
    finite lattice (review round 14 F5)."""
    v, x, y = (np.asarray(a, float) for a in (values, x, y))
    if v.ndim != 2 or x.shape != v.shape or y.shape != v.shape:
        raise ValueError(f"values, x and y must share one 2-D shape, not {v.shape}, "
                         f"{x.shape}, {y.shape}")
    if not (np.isfinite(x).all() and np.isfinite(y).all()):
        raise ValueError("the lattice coordinates must be finite")
    return v, x, y


def graded_up(values, x, y, grade):
    """The smallest gradation-feasible field that is >= ``values``.

    ``values`` are finite, or -inf where nothing is imposed."""
    v, x, y = _lattice(values, x, y)
    if np.isnan(v).any() or np.isposinf(v).any():
        raise ValueError("values must be finite, or -inf where nothing is imposed")
    return -_limit(-v, x, y, _grade(grade))


def band_field(x, y, line_xy, targets, half_width_m):
    """``targets`` along a polyline, spread over a band; NaN elsewhere.

    ``line_xy`` is the line in the lattice's metric coordinates and
    ``targets`` one size per vertex; a lattice point within ``half_width_m``
    of the line takes the target interpolated at its projection.
    """
    import shapely

    line_xy = np.asarray(line_xy, float)
    targets = np.asarray(targets, float)
    if len(line_xy) < 2 or len(targets) != len(line_xy):
        raise ValueError("a band needs a line of two points or more and one target per point")
    line = shapely.LineString(line_xy)
    pts = shapely.points(np.ravel(x), np.ravel(y))
    d = shapely.distance(pts, line)
    near = d <= half_width_m
    s_vertex = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(line_xy, axis=0), axis=1))]
    out = np.full(np.size(x), np.nan)
    s = shapely.line_locate_point(line, pts[near])
    out[near] = np.interp(s, s_vertex, targets)
    return out.reshape(np.shape(x))


def _check_sizing_inputs(ambient, x, y, floor, bands) -> None:
    """Sizes the limiter can work with (review round 13 F9): finite positive
    ambient sizes, a finite non-negative floor, and bands that are positive
    where set (NaN marks off-band). The lattice and the gradation are
    checked by ``_lattice`` and ``_grade``."""
    shape = np.shape(ambient)
    if not (np.isfinite(ambient).all() and (np.asarray(ambient) > 0).all()):
        raise ValueError("ambient sizes must be finite and positive")
    if floor is not None and (np.shape(floor) != shape or not np.isfinite(floor).all()
                              or (np.asarray(floor) < 0).any()):
        raise ValueError("the floor must be finite, non-negative, with ambient's shape")
    for k, b in enumerate(bands):
        on = np.isfinite(b)
        if np.shape(b) != shape or np.isinf(b).any() or (b[on] <= 0).any():
            raise ValueError(f"band {k}: positive sizes where set, NaN elsewhere, "
                             f"ambient's shape")


def compose_sizing(ambient, x, y, *, grade, floor=None, bands=()):
    """The final size field and a report; see the module docstring."""
    bands = [np.asarray(b, float) for b in bands]     # traversed twice (review r3 F3)
    # normalised once, and these arrays are the ones used (round 14 F5)
    ambient, x, y = _lattice(ambient, x, y)
    grade = _grade(grade)
    _check_sizing_inputs(ambient, x, y, floor, bands)
    h = _limit(ambient, x, y, grade)
    report = {}
    if floor is not None:
        up = graded_up(floor, x, y, grade)
        report["raised_by_floor"] = int((up > h).sum())
        h = np.maximum(h, up)
    for k, band in enumerate(bands):
        on = np.isfinite(band)
        if not on.any():
            raise ValueError(f"band {k} covers no lattice point")
        lo = _limit(np.where(on, band, np.inf), x, y, grade)
        hi = graded_up(np.where(on, band, -np.inf), x, y, grade)
        h = np.minimum(np.maximum(h, hi), lo)
        report[f"band_{k}_cells"] = int(on.sum())
    # every band must come out at its own target: two bands closer than their
    # sizes allow under the gradation cannot both hold, and the later one
    # would silently win (review round 2 F17). A band whose own sizes change
    # faster along the line than the gradation allows is smoothed a little
    # (0.9 % on the Tokyo Bay interface, 2026-10-01); up to BAND_TOLERANCE
    # that is accepted and reported, beyond it refused.
    for k, band in enumerate(bands):
        on = np.isfinite(band)
        if floor is not None:
            # where a band leaves the final field below the time-step floor
            # (review F4; counted on the field that comes out, round 3 F1).
            # The caller decides what to do: meshes are made from the real
            # depths, and the time step is the depth stage's business (owner,
            # 2026-10-01), so 446 reports it
            report[f"band_{k}_below_floor_cells"] = int(
                (on & (h < np.asarray(floor) - 1e-6)).sum())
        dev = float(np.max(np.abs(h[on] - band[on]) / band[on]))
        report[f"band_{k}_max_rel_deviation"] = dev
        if dev > BAND_TOLERANCE:
            raise ValueError(f"band {k} cannot hold its sizes (off by up to {dev:.1%}): "
                             "its own sizes vary faster than the gradation allows, or "
                             "another band is too close")
    if floor is not None:
        report["below_floor_fraction"] = float(np.mean(h < np.asarray(floor) - 1e-6))
    return h, report


def _signed_areas(nodes, elements):
    a = nodes[elements]
    return 0.5 * ((a[:, 1, 0] - a[:, 0, 0]) * (a[:, 2, 1] - a[:, 0, 1])
                  - (a[:, 2, 0] - a[:, 0, 0]) * (a[:, 1, 1] - a[:, 0, 1]))


def merge_outer(base: Fort14Mesh, outer_nodes, outer_elements, interface_outer,
                interface_base, outer_open, *, tol_m: float = 0.5) -> Fort14Mesh:
    """Append an outer mesh to ``base`` through their shared interface nodes.

    ``interface_outer[i]`` (an outer node) and ``interface_base[i]`` (a base
    node) are the same point, at most ``tol_m`` apart; the base coordinate
    is kept.  The base's nodes and elements stay first and unchanged; outer
    elements are turned to the base's orientation.  ``outer_open`` is the new
    open boundary in outer indices; it becomes the only open boundary.
    Depths of the outer nodes are NaN until a depth stage fills them.

    The base must hold float64 coordinates and depths and int64 elements
    (what the readers give): the merged arrays have those types, and the
    base is carried bit for bit (review rounds 9 F6, 10 F10). Every index
    is checked to be a whole number in range before it is used (round 9 F4).
    """
    from fvcom_mesh_tools.io.fvcom_native import _indices

    # a NaN or infinite tolerance would accept any interface (round 17 F4)
    if not (isinstance(tol_m, (int, float, np.integer, np.floating)) and np.isfinite(tol_m)
            and tol_m >= 0):
        raise ValueError(f"tol_m must be finite and non-negative, not {tol_m!r}")
    if (np.asarray(base.nodes).dtype != np.float64 or np.asarray(base.depths).dtype != np.float64
            or np.asarray(base.elements).dtype != np.int64):
        raise ValueError(f"the base must be float64 nodes and depths and int64 elements "
                         f"(nodes {np.asarray(base.nodes).dtype}, depths "
                         f"{np.asarray(base.depths).dtype}, elements "
                         f"{np.asarray(base.elements).dtype}; review round 10 F10)")
    outer_nodes = np.asarray(outer_nodes, float)
    if outer_nodes.ndim != 2 or outer_nodes.shape[1] < 2 or not np.isfinite(outer_nodes).all():
        raise ValueError("outer nodes must be finite (N, 2) coordinates")
    n_out = len(outer_nodes)
    outer_elements = _indices(outer_elements, n_out, "outer elements", ndim=2)
    if outer_elements.shape[1:] != (3,):
        raise ValueError(f"outer elements must be (NE, 3), not {outer_elements.shape}")
    io_ = _indices(interface_outer, n_out, "interface_outer")
    ib = _indices(interface_base, base.n_nodes, "interface_base")
    outer_open = _indices(outer_open, n_out, "outer_open")
    if io_.shape != ib.shape or len(io_) < 2:
        raise ValueError("the interface needs two nodes or more, paired one to one")
    gap = np.hypot(*(outer_nodes[io_] - base.nodes[ib, :2]).T)
    if gap.max() > tol_m:
        raise ValueError(f"interface nodes differ from the base by up to {gap.max():.3f} m")
    nb = base.n_nodes
    new = np.setdiff1d(np.arange(len(outer_nodes)), io_)
    index = np.full(len(outer_nodes), -1, np.int64)
    index[io_] = ib
    index[new] = nb + np.arange(len(new))
    elems = index[outer_elements]
    sign_base = np.sign(np.median(_signed_areas(base.nodes[:, :2], base.elements)))
    nodes = np.vstack([base.nodes[:, :2], outer_nodes[new]])
    flip = np.sign(_signed_areas(nodes, elems)) != sign_base
    elems[flip] = elems[flip][:, [0, 2, 1]]
    depths = np.r_[np.asarray(base.depths, float), np.full(len(new), np.nan)]
    merged = replace(base, nodes=nodes, depths=depths,
                     elements=np.vstack([base.elements, elems]),
                     open_boundaries=[index[np.asarray(outer_open, np.int64)]],
                     land_boundaries=[])
    return merged


def verify_frozen_base(merged: Fort14Mesh, base: Fort14Mesh, interface_base) -> dict:
    """Is the base inside ``merged`` exactly as it was?  Raises if not.

    Checks that the first ``NP`` nodes and ``NE`` elements are the base's,
    bit for bit, and that every base open-boundary edge is now interior
    (one base element and one outer element on it).
    """
    nb, eb = base.n_nodes, base.n_elements

    def same_bits(a, b):
        # dtype and bit pattern: numerical equality lets -0.0 pass for +0.0
        # (review round 8 F12)
        a, b = np.ascontiguousarray(a), np.ascontiguousarray(b)
        return a.dtype == b.dtype and a.shape == b.shape and a.tobytes() == b.tobytes()

    if not same_bits(merged.nodes[:nb, :2], base.nodes[:, :2]):
        raise ValueError("base node coordinates changed")
    if not same_bits(merged.elements[:eb], base.elements):   # round 10 F10
        raise ValueError("base elements changed")
    if not same_bits(merged.depths[:nb], base.depths):
        raise ValueError("base depths changed")
    # The interface given must be the base's open boundary, every edge of
    # it, as whole indices: an empty or partial one let an unmerged base pass
    # (review round 16 F5).
    from fvcom_mesh_tools.io.fvcom_native import _indices

    ib = _indices(interface_base, nb, "interface_base")

    def _edges(chain):
        c = np.asarray(chain, np.int64)
        return {frozenset((int(a), int(b))) for a, b in zip(c[:-1], c[1:])}

    want = set().union(*(_edges(c) for c in base.open_boundaries)) if base.open_boundaries \
        else set()
    if _edges(ib) != want:
        raise ValueError(f"interface_base has {len(_edges(ib))} edge(s); the base's open "
                         f"boundary has {len(want)}, and they must be the same")
    e = np.sort(np.vstack([merged.elements[:, [0, 1]], merged.elements[:, [1, 2]],
                           merged.elements[:, [2, 0]]]), axis=1)
    owner = np.tile(np.arange(merged.n_elements), 3)
    keys = e[:, 0] * merged.n_nodes + e[:, 1]
    xy = merged.nodes[:, :2]
    for a, b in zip(ib[:-1], ib[1:]):
        k = min(a, b) * merged.n_nodes + max(a, b)
        who = owner[keys == k]
        if len(who) != 2 or (who < eb).sum() != 1:
            raise ValueError(f"interface edge {a}-{b} is not shared by one base and one "
                             f"outer element (elements {who.tolist()})")
        # the two elements must lie on opposite sides of the edge, or the outer
        # one overlaps the base (review F10)
        side = []
        for e_ in who:
            c = [v for v in merged.elements[e_] if v not in (a, b)][0]
            side.append(np.sign((xy[b, 0] - xy[a, 0]) * (xy[c, 1] - xy[a, 1])
                                - (xy[b, 1] - xy[a, 1]) * (xy[c, 0] - xy[a, 0])))
        if side[0] * side[1] >= 0:
            raise ValueError(f"interface edge {a}-{b}: the base and the outer element are "
                             "on the same side (they overlap)")
    return {"n_base_nodes": nb, "n_base_elements": eb,
            "n_interface_edges": len(want),        # edges verified (round 17 F8)
            "n_nodes": merged.n_nodes, "n_elements": merged.n_elements}


def check_no_overlap(merged: Fort14Mesh, n_base_elements: int, rel_tol: float = 1e-9) -> dict:
    """No new element may cover any part of the base (review round 2 F9).

    The shared interface is a line, so an outer element may touch the base
    footprint only along it: an intersection with positive area (beyond
    ``rel_tol`` of the element's area, for round-off) is an overlap, wherever
    it is -- not only at the interface edges ``verify_frozen_base`` checks.
    Raises on the first overlap; returns counts.
    """
    import shapely

    if not (isinstance(rel_tol, (int, float, np.integer, np.floating)) and np.isfinite(rel_tol)
            and rel_tol >= 0):
        raise ValueError(f"rel_tol must be finite and non-negative, not {rel_tol!r}")
    if not (isinstance(n_base_elements, (int, np.integer))
            and 0 <= n_base_elements <= merged.n_elements):
        raise ValueError(f"n_base_elements must be in [0, {merged.n_elements}]")
    xy = merged.nodes[:, :2]
    base = shapely.union_all(shapely.polygons(xy[merged.elements[:n_base_elements]]))
    outer = shapely.polygons(xy[merged.elements[n_base_elements:]])
    tree = shapely.STRtree(outer)
    cand = np.unique(tree.query(base, predicate="intersects"))
    if len(cand):
        area = shapely.area(shapely.intersection(outer[cand], base))
        own = shapely.area(outer[cand])
        bad = cand[area > rel_tol * own]
        if len(bad):
            k = int(bad[0])
            raise ValueError(f"{len(bad)} new element(s) overlap the base, e.g. element "
                             f"{n_base_elements + k} ({float(area[cand == k][0]):.3g} m2)")
    return {"n_outer_touching_base": int(len(cand))}


def check_island_holes(mesh: Fort14Mesh, land, n_base_nodes: int) -> dict:
    """Every hole in the new part of the mesh must hold some real land.

    A closed boundary loop without an open-boundary node, touching a node
    beyond the base's ``n_base_nodes``, is an island; the polygon it bounds
    must intersect ``land`` (the supplied land, in the mesh's coordinates)
    over a positive area. A hole in open water -- elements lost in finishing
    or repair -- passed every other check (review round 11 F7). Raises on
    the first such hole; returns counts.
    """
    import shapely

    from fvcom_mesh_tools.io.fvcom_native import boundary_loops

    xy = np.asarray(mesh.nodes)[:, :2]
    obc = set(int(v) for c in mesh.open_boundaries for v in np.asarray(c).tolist())
    n_islands, bad = 0, []
    for loop in boundary_loops(np.asarray(mesh.elements)):
        loop = [int(v) for v in loop]
        if obc & set(loop) or max(loop) < n_base_nodes:
            continue
        n_islands += 1
        hole = shapely.Polygon(xy[loop])
        if not hole.is_valid:
            hole = hole.buffer(0)
        if land.intersection(hole).area <= 0:
            bad.append((loop[0], float(hole.area)))
    if bad:
        v, a = bad[0]
        raise ValueError(f"{len(bad)} hole(s) in the new mesh hold no land, e.g. at node {v} "
                         f"({a:.0f} m2): sea is missing there")
    return {"n_new_islands": n_islands}


#: Covered land is a defect where it is wider than one element: each patch of
#: land under the new elements is shrunk by this fraction of the local edge
#: length; whatever survives is land the mesh could have resolved. Narrower
#: patches -- the coast's approximation, islets below the element size --
#: are what the resolution principle drops. Tokyo Bay - Enshu, 2026-10-02:
#: 3,681 covered patches (34 km2), none survives 0.5.
LAND_COVER_ERODE = 0.5


def check_land_cover(mesh: Fort14Mesh, land, n_base_elements: int,
                     erode: float = LAND_COVER_ERODE) -> dict:
    """No land the mesh could have resolved may lie under the new elements.

    ``land`` is the true land in the mesh's coordinates, without the base's
    footprint. Each connected patch of it under the new elements is shrunk
    inward by ``erode`` times the median edge length of the elements holding
    it; a patch with anything left is wider than an element there, and
    fails. Patches are judged on their own: land outside the mesh, such as a
    mainland a covered peninsula belongs to, does not change the verdict
    (review rounds 15 F2, 16 F3, 19 F1). Raises on the first such patch;
    returns counts.
    """
    import shapely

    ne = len(np.asarray(mesh.elements))
    if not (isinstance(n_base_elements, (int, np.integer)) and 0 <= n_base_elements <= ne):
        raise ValueError(f"n_base_elements must be an integer in [0, {ne}], not "
                         f"{n_base_elements!r}")
    if not (isinstance(erode, (int, float, np.integer, np.floating)) and np.isfinite(erode)
            and erode > 0):
        raise ValueError(f"erode must be finite and positive, not {erode!r}")
    xy = np.asarray(mesh.nodes)[:, :2]
    tri = shapely.polygons(xy[np.asarray(mesh.elements)[n_base_elements:]])
    area = shapely.area(tri)
    tree = shapely.STRtree(tri)
    covered = land.intersection(shapely.union_all(tri))
    patches = [g for g in getattr(covered, "geoms", [covered])
               if g.geom_type == "Polygon" and g.area > 0]
    worst, bad = 0.0, []
    for patch in patches:
        under = tree.query(patch, predicate="intersects")
        held = shapely.area(shapely.intersection(tri[under], patch))
        under = under[held > 1e-9 * area[under]]   # not those touching an edge
        if not len(under):
            continue
        edge = float(np.median(np.sqrt(4.0 * area[under] / np.sqrt(3.0))))
        left = patch.buffer(-erode * edge).area
        worst = max(worst, left)
        if left > 0:
            c = patch.representative_point()
            bad.append((c.x, c.y, patch.area, edge))
    if bad:
        x, y, a, e = bad[0]
        raise ValueError(f"{len(bad)} patch(es) of land the mesh could resolve lie under new "
                         f"elements, e.g. {a:.0f} m2 at ({x:.0f}, {y:.0f}) under elements of "
                         f"about {e:.0f} m")
    return {"n_covered_patches": len(patches),
            "covered_area_m2": float(sum(p.area for p in patches)),
            "max_area_left_after_erosion_m2": worst, "erode": erode}


def land_segments(elements, open_chains) -> list[tuple[int, np.ndarray]]:
    """The land boundary runs: every boundary loop minus its open-boundary edges.

    Returns ``(ibtype, nodes)`` runs as ``Fort14Mesh.land_boundaries`` wants
    them: 0 for a run on a loop that carries an open boundary (the mainland),
    1 for a closed loop without one (an island).  An open chain's end nodes
    appear in the land runs too, as in ADCIRC.
    """
    from fvcom_mesh_tools.io.fvcom_native import boundary_loops

    open_edges = set()
    for c in open_chains:
        c = np.asarray(c, np.int64)
        open_edges |= {frozenset((int(a), int(b))) for a, b in zip(c[:-1], c[1:])}
    out = []
    for loop in boundary_loops(np.asarray(elements)):
        loop = [int(v) for v in loop]
        n = len(loop)
        is_open = [frozenset((loop[i], loop[(i + 1) % n])) in open_edges for i in range(n)]
        if not any(is_open):
            out.append((1, np.array(loop + [loop[0]], np.int64)))
            continue
        start = is_open.index(True)
        run: list[int] = []
        for k in range(1, n + 1):
            i = (start + k) % n
            if is_open[i]:
                if run:      # one land edge is a run too (review F11)
                    out.append((0, np.array(run + [loop[i]], np.int64)))
                run = []
            else:
                run.append(loop[i])
        if run:
            out.append((0, np.array(run + [loop[(start + n) % n]], np.int64)))
    return out


def round_depths_inside(h, hmin, hmax=None, decimals: int = 6) -> np.ndarray:
    """``h`` rounded to ``decimals``, kept inside ``[hmin, hmax]``.

    Rounding after the limiter could cross a bound that is not itself on the
    grid of ``decimals`` (3.0000004 -> 3.0 below a 3.0000004 floor); the
    bounds are moved inward to that grid first (review round 10 F12).
    """
    q = 10.0 ** decimals
    lo = np.ceil(hmin * q) / q
    hi = np.inf if hmax is None else np.floor(hmax * q) / q
    if lo > hi:
        raise ValueError(f"no {decimals}-decimal depth lies in [{hmin}, {hmax}]")
    return np.clip(np.round(np.asarray(h, float), decimals), lo, hi)


def rfactor_smooth_free(h0, ei, ej, free, *, rmax, hmin, hmax=None, max_iter=5000):
    """r-factor limiter that moves only the ``free`` nodes.

    As ``dem.m7001.rfactor_smooth`` (Beckmann-Haidvogel), but a node that is
    not free keeps its depth: on an edge with one fixed end the free end
    takes the whole correction, and an edge with two fixed ends is left as
    it is (its r is the base's own).  ``hmin`` and ``hmax`` bound the free
    depths *during* the smoothing, so the result satisfies both (a cap
    applied afterwards could break the r-factor; review F3). Returns
    ``(depth, iterations, max r over edges with a free end)``; raises when the
    limit is not reached -- infeasible (e.g. a fixed 1 m node beside a free
    node held at 3 m) or not converged within ``max_iter``.
    """
    h = np.asarray(h0, float).copy()
    free = np.asarray(free, bool)
    ei, ej = np.asarray(ei), np.asarray(ej)
    # finite positive depths and sane controls, or NaN slips through the
    # r > limit test (review round 2 F15)
    if not (0 < rmax < 1 and np.isfinite(hmin) and hmin > 0
            and (hmax is None or (np.isfinite(hmax) and hmax >= hmin))
            and int(max_iter) >= 1):
        raise ValueError(f"bad controls: rmax {rmax}, hmin {hmin}, hmax {hmax}, "
                         f"max_iter {max_iter}")
    live = free[ei] | free[ej]
    ei, ej = ei[live], ej[live]
    fi, fj = free[ei], free[ej]
    used = np.unique(np.r_[ei, ej, np.flatnonzero(free)])
    if not (np.isfinite(h[used]).all() and (h[used] > 0).all()):
        raise ValueError("depths on the limited edges and free nodes must be finite and "
                         "positive")
    # the bounds hold for every free node, on an edge or not (review round 9
    # F10)
    if hmax is not None:
        h = np.where(free, np.minimum(h, hmax), h)
    h = np.where(free, np.maximum(h, hmin), h)
    if not len(ei):
        return h, 0, 0.0
    for it in range(int(max_iter)):
        hi, hj = h[ei], h[ej]
        r = np.abs(hi - hj) / (hi + hj)
        bad = r > rmax + 1e-9
        if not bad.any():
            return h, it, float(r.max())
        excess = np.where(bad, np.abs(hi - hj) - rmax * (hi + hj), 0.0)
        sgn = np.sign(hi - hj)                       # +1: i is the deeper end
        both = fi & fj
        # a free end moves by half the excess when the other end moves too,
        # by all of it (scaled for the r-denominator) when the other is fixed
        di = np.where(both, excess / 2, excess / (1 + rmax))
        dj = np.where(both, excess / 2, excess / (1 + rmax))
        add = np.zeros_like(h)
        cnt = np.zeros_like(h)
        np.add.at(add, ei[fi & bad], -(sgn * di)[fi & bad])
        np.add.at(cnt, ei[fi & bad], 1.0)
        np.add.at(add, ej[fj & bad], (sgn * dj)[fj & bad])
        np.add.at(cnt, ej[fj & bad], 1.0)
        step = add / np.where(cnt > 0, cnt, 1.0)
        h = np.where(free, np.clip(h + step, hmin, np.inf if hmax is None else hmax), h)
    hi, hj = h[ei], h[ej]
    r = np.abs(hi - hj) / (hi + hj)
    if np.isfinite(r).all() and r.max() <= rmax + 1e-9:
        return h, int(max_iter), float(r.max())     # met on the last pass (round 2 F14)
    k = int(np.argmax(r))
    raise ValueError(f"r-factor limit {rmax} not reached in {max_iter} iterations: "
                     f"r = {r[k]:.4f} on edge {int(ei[k])}-{int(ej[k])} "
                     f"(depths {hi[k]:.3f}, {hj[k]:.3f} m); the depth bounds may make "
                     "it infeasible")


def trim_lone_corners(elements, mutable, keep_nodes=(), max_rounds=20):
    """Drop the element under a node that no other element shares.

    A cape one element wide ends in a node that sits in a single element;
    FVCOM never updates such a node. Bisecting the element (as
    ``walls.open_lone_corners`` does for wall bends) leaves two slivers at a
    cape tip, so here the tip is not resolved: the element goes, as the
    resolution principle says for what the element size cannot carry. Only
    ``mutable`` elements are dropped, never one holding a node in
    ``keep_nodes`` (the open boundary). Repeats while new lone nodes appear.

    Returns ``(elements, mutable, report)``.
    """
    t = np.asarray(elements, np.int64)
    mut = np.asarray(mutable, bool)
    keep = set(int(v) for v in keep_nodes)
    dropped, left, limited = 0, [], False
    for _ in range(max_rounds):
        if len(t) == 0:
            break
        count = np.bincount(t.ravel(), minlength=int(t.max()) + 1)
        edge_count: dict[tuple[int, int], int] = {}
        for a, b in np.sort(np.vstack([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]]), axis=1).tolist():
            edge_count[(a, b)] = edge_count.get((a, b), 0) + 1
        drop = np.zeros(len(t), bool)
        left = []
        for v in np.flatnonzero(count == 1).tolist():
            k = int(np.flatnonzero((t == v).any(axis=1))[0])
            a, b = sorted(int(x) for x in t[k] if x != v)
            # only a spike: the side facing the lone node is shared, so the
            # element's removal leaves no new lone node behind
            # Removals are judged against what survives: each one takes its
            # edges out of the count, so two elements that each lean on the
            # other's shared side are not both dropped (review round 9 F9).
            if (not drop[k] and mut[k] and not (set(t[k].tolist()) & keep)
                    and edge_count[(a, b)] == 2):
                drop[k] = True
                for e in ((t[k, 0], t[k, 1]), (t[k, 1], t[k, 2]), (t[k, 2], t[k, 0])):
                    edge_count[tuple(sorted(int(x) for x in e))] -= 1
            elif not drop[k]:
                left.append(v)
        if not drop.any():
            break
        t, mut = t[~drop], mut[~drop]
        dropped += int(drop.sum())
    else:
        limited = True
    # the report is of the mesh returned, not of the last round's start
    # (review round 10 F11)
    if len(t):
        count = np.bincount(t.ravel(), minlength=int(t.max()) + 1)
        left = np.flatnonzero(count == 1).tolist()
    else:
        left = []
    return t, mut, {"n_elements_dropped": dropped, "lone_nodes_left": left,
                    "round_limit_reached": limited}
