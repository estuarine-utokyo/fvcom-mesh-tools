"""Design an open boundary: straight sides, rounded corners, coast-normal ends.

The owner's rules for an open boundary (2026-09-29):

* it meets the coast **at a right angle** at both ends;
* between the ends it is a polygon of straight sides whose corners are
  rounded with circular arcs (sharp corners at sea leave the boundary normal,
  and the sponge and nesting zones built on it, undefined at the corner);
* its node spacing never falls below what the time step allows
  (``dt * sqrt(g * depth) / Cr``), so the boundary is not what limits dt.

Everything here works in a metric plane (the caller projects); nothing reads
data.  ``notebooks/444_design_obc.py`` applies it to a recipe and writes the
boundary's nodes, which then enter mesh generation as an input.
"""

from __future__ import annotations

import numpy as np

from fvcom_mesh_tools._checks import no_complex

__all__ = [
    "bearing_vector",
    "coast_normal",
    "fillet",
    "resample",
    "ray_intersection",
]


def bearing_vector(bearing_deg: float) -> np.ndarray:
    """Unit vector for a compass bearing (0 = north, 90 = east)."""
    b = np.radians(bearing_deg)
    return np.array([np.sin(b), np.cos(b)])


def coast_normal(land, x: float, y: float, chord_m: float = 3000.0):
    """Seaward normal of the coastline nearest ``(x, y)``.

    ``land`` is a (Multi)Polygon in the same metric plane.  The coast
    direction is the chord from ``chord_m`` before to ``chord_m`` after the
    nearest coast point, measured along the coastline -- the same measure
    the boundary's first edge is later checked against.  Returns
    ``(bearing_deg, coast_point)``.
    """
    import shapely

    from fvcom_mesh_tools._checks import real_scalar

    if real_scalar(chord_m, "chord_m") <= 0:         # not a bool, finite, positive (round 56 F2)
        raise ValueError(f"chord_m must be positive, not {chord_m!r}")
    p = shapely.Point(x, y)
    polys = list(getattr(land, "geoms", [land]))
    ring = min((g.exterior for g in polys), key=lambda r: r.distance(p))
    s = ring.project(p)
    q = np.array(ring.interpolate(s).coords[0])
    c = (np.array(ring.interpolate((s + chord_m) % ring.length).coords[0])
         - np.array(ring.interpolate((s - chord_m) % ring.length).coords[0]))
    if not np.any(c):
        raise ValueError(f"no coastline direction near ({x:.0f}, {y:.0f})")
    b = np.degrees(np.arctan2(c[0], c[1]))
    # the whole departure must stay off land, as a continuous segment from just
    # off the coast point: a thin strip of land can be crossed between probe
    # points (review F7, round 2 F5)
    for nb in ((b + 90) % 360, (b - 90) % 360):
        v = bearing_vector(nb)
        seg = shapely.LineString([q + 1e-3 * v, q + 2 * chord_m / 3 * v])
        if not land.intersects(seg):
            return float(nb), q
    raise ValueError(f"both normals at ({q[0]:.0f}, {q[1]:.0f}) point onto land")


def ray_intersection(p, u, q, w) -> np.ndarray:
    """Where the ray ``p + s u`` meets the ray ``q + t w`` (both s, t > 0)."""
    # known finite 2-vectors, and a finite positive meeting point (review
    # round 26 F6): NaN slipped past the s, t > 0 test
    if any(np.ma.is_masked(v) for v in (p, u, q, w)):
        raise ValueError("a ray has masked values")
    no_complex(p=p, u=u, q=q, w=w)
    p, u, q, w = (np.asarray(v, float) for v in (p, u, q, w))
    if any(v.shape != (2,) or not np.isfinite(v).all() for v in (p, u, q, w)):
        raise ValueError("rays need finite (2,) points and directions")
    if not (u.any() and w.any()):
        raise ValueError("a ray has no direction")
    # scaled by the largest component first, so the norm neither overflows
    # nor underflows (review round 28 F4)
    u, w = u / np.abs(u).max(), w / np.abs(w).max()
    nu, nw = np.linalg.norm(u), np.linalg.norm(w)
    # unit directions: whether two sides are parallel does not depend on how
    # long their direction vectors are (review round 27 F5)
    u, w = u / nu, w / nw
    a = np.column_stack([u, -w])
    if abs(np.linalg.det(a)) < 1e-12:
        raise ValueError("the two sides are parallel")
    with np.errstate(over="ignore", invalid="ignore"):
        s, t = np.linalg.solve(a, q - p)
        point = p + s * u
    if not (np.isfinite(s) and np.isfinite(t)) or s <= 0 or t <= 0:
        raise ValueError("the two sides meet behind one of their start points")
    if not np.isfinite(point).all():
        raise ValueError("the meeting point is not finite")
    return point


def fillet(vertices, radii, n_arc: int = 60) -> np.ndarray:
    """A polyline through ``vertices`` with every interior corner rounded.

    ``radii[i]`` is the arc radius at interior vertex ``i + 1``.  The arc is
    tangent to both sides, so the straight parts keep their direction; a
    radius too large for its sides is refused rather than overlapped.
    """
    if np.ma.is_masked(vertices) or np.ma.is_masked(radii):     # review round 25 F4
        raise ValueError("vertices or radii have masked values")
    radii = list(radii)                    # read once, even from an iterator (round 38 F1)
    no_complex(vertices=vertices, radii=radii)
    x = np.asarray(vertices, float)
    # planar (N >= 2, 2): a third column would enter some lengths and not
    # others (review round 34 F5)
    if x.ndim != 2 or x.shape[1] != 2 or len(x) < 2:
        raise ValueError(f"vertices must be (N >= 2, 2), not {x.shape}")
    if len(radii) != len(x) - 2:
        raise ValueError(f"{len(x) - 2} interior corner(s), {len(radii)} radii")
    # finite positive radii, real sides, a real arc (review r2 F20)
    if not np.isfinite(x).all() or np.any(np.linalg.norm(np.diff(x, axis=0), axis=1) <= 0):
        raise ValueError("the vertices must be finite and distinct")
    if not all(isinstance(r, (int, float, np.integer, np.floating))
               and not isinstance(r, (bool, np.bool_)) and np.isfinite(float(r)) and float(r) > 0
               for r in radii):
        raise ValueError(f"radii must be finite and positive numbers (not booleans): {radii}")
    if int(n_arc) < 2:
        raise ValueError("n_arc must be 2 or more")
    out = [x[0]]
    last_tangent = 0.0                       # length already used on the next side
    for i in range(1, len(x) - 1):
        a, b, c = x[i - 1], x[i], x[i + 1]
        u = (a - b) / np.linalg.norm(a - b)
        w = (c - b) / np.linalg.norm(c - b)
        th = np.arccos(np.clip(u @ w, -1, 1))
        if th < 1e-6:
            # the path doubles back on itself: no tangent arc exists (review
            # of the extend tools, round 3 F13)
            raise ValueError(f"the path reverses at corner {i}")
        if np.pi - th < 1e-6:                # straight on: nothing to round
            out.append(b)
            last_tangent = 0.0
            continue
        r = float(radii[i - 1])
        k = r / np.tan(th / 2)               # tangent length on each side
        la, lc = np.linalg.norm(a - b), np.linalg.norm(c - b)
        # a radius that fits exactly may miss by round-off (10 / tan(pi/4) is
        # 10.000000000000002): a relative tolerance, the excess clamped (round 56 F3)
        if k + last_tangent > la * (1 + 1e-9) or k > lc * (1 + 1e-9):
            raise ValueError(f"radius {r:.0f} m does not fit at corner {i}")
        k = max(min(k, lc, la - last_tangent), 0.0)
        p1, p2 = b + u * k, b + w * k
        bis = (u + w) / np.linalg.norm(u + w)
        cen = b + bis * r / np.sin(th / 2)
        a1 = np.arctan2(*(p1 - cen)[::-1])
        a2 = np.arctan2(*(p2 - cen)[::-1])
        da = (a2 - a1 + np.pi) % (2 * np.pi) - np.pi
        out += [cen + r * np.array([np.cos(a1 + da * f), np.sin(a1 + da * f)])
                for f in np.linspace(0, 1, n_arc)]
        last_tangent = k
    out.append(x[-1])
    line = np.array(out)
    if not np.isfinite(line).all():           # finite inputs can overflow in the arc maths
        raise ValueError("the rounded line is not finite: the coordinates are too large")
    return line


def resample(line, spacing) -> np.ndarray:
    """Nodes along ``line`` from its first to its last point.

    ``spacing`` is a number (metres) or a callable ``f(xy) -> metres`` giving
    the wanted spacing at a point. It is a floor: every edge -- the straight
    chord between consecutive nodes, which is shorter than the arc on a bend
    -- is at least the spacing at both of its ends (review F5, round 2 F4).
    Each step is lengthened until that holds; the last node is the line's
    end, and the remainder is added to the last step (a node is dropped
    when the last step would be short). A line shorter than its spacing,
    non-finite spacing or coordinates, and a step that does not advance are
    refused (review F19, round 2 F12).
    """
    if np.ma.is_masked(line):                   # unknown, not its fill (round 25 F4)
        raise ValueError("the line has masked coordinates")
    no_complex(line=line)
    xy = np.asarray(line, float)
    if xy.ndim != 2 or xy.shape[1] != 2 or len(xy) < 2:
        raise ValueError(f"the line must be (N >= 2, 2), not {xy.shape}")
    if not np.isfinite(xy).all():
        raise ValueError("the line has non-finite coordinates")
    seg = np.linalg.norm(np.diff(xy, axis=0), axis=1)
    s = np.r_[0.0, np.cumsum(seg)]
    total = s[-1]
    if total <= 0:
        raise ValueError("the line has no length")

    def at(t):
        return np.column_stack([np.interp(t, s, xy[:, 0]), np.interp(t, s, xy[:, 1])])

    if not callable(spacing):
        no_complex(spacing=spacing)
        if np.asarray(spacing).dtype.kind == "b":
            raise ValueError("spacing must be a number, not a boolean")
    f = spacing if callable(spacing) else (lambda p, h=float(spacing): np.full(len(p), h))

    def h_at(t):
        got = f(at([min(t, total)]))
        no_complex(spacing=got)            # a callable's result too (review round 38 F3)
        if np.asarray(got).dtype.kind == "b":
            raise ValueError("spacing must be numbers, not booleans")
        h = float(got[0])
        if not np.isfinite(h) or h <= 0:
            raise ValueError(f"spacing must be positive and finite (got {h} at {t:.1f} m)")
        return h

    def chord(t0, t1):
        return float(np.linalg.norm(at([t1])[0] - at([t0])[0]))

    def ok(t0, t1):
        return chord(t0, t1) >= max(h_at(t0), h_at(t1)) * (1 - 1e-12)

    t = [0.0]
    while True:
        t0 = t[-1]
        t1 = t0 + h_at(t0)
        for _ in range(64):                 # lengthen until the chord meets both ends
            if t1 >= total or ok(t0, t1):
                break
            t1 = t0 + max(t1 - t0, max(h_at(t0), h_at(t1))) * 1.01
        if not (np.isfinite(t1) and t1 > t0):
            raise ValueError(f"resampling does not advance at {t0:.1f} m")
        if t1 >= total:
            break
        if not ok(t0, t1):
            raise ValueError(f"no step from {t0:.1f} m meets the spacing")
        t.append(t1)
    # the end: the last node is the line's end; drop nodes until the last
    # step meets the floor too
    while len(t) > 1 and not ok(t[-1], total):
        t.pop()
    if not ok(t[-1], total):
        raise ValueError(f"the line ({total:.0f} m) is shorter than its spacing")
    return at(np.r_[t, total])
