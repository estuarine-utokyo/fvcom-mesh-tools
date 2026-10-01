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
    # the whole departure must stay off land, not only a far test point: a
    # thin strip of land can be crossed and left behind (review F7)
    reach = np.linspace(chord_m / 50, 2 * chord_m / 3, 25)
    for nb in ((b + 90) % 360, (b - 90) % 360):
        ray = shapely.points(q + reach[:, None] * bearing_vector(nb))
        if not shapely.contains(land, ray).any():
            return float(nb), q
    raise ValueError(f"both normals at ({q[0]:.0f}, {q[1]:.0f}) point onto land")


def ray_intersection(p, u, q, w) -> np.ndarray:
    """Where the ray ``p + s u`` meets the ray ``q + t w`` (both s, t > 0)."""
    p, u, q, w = (np.asarray(v, float) for v in (p, u, q, w))
    a = np.column_stack([u, -w])
    if abs(np.linalg.det(a)) < 1e-12:
        raise ValueError("the two sides are parallel")
    s, t = np.linalg.solve(a, q - p)
    if s <= 0 or t <= 0:
        raise ValueError("the two sides meet behind one of their start points")
    return p + s * u


def fillet(vertices, radii, n_arc: int = 60) -> np.ndarray:
    """A polyline through ``vertices`` with every interior corner rounded.

    ``radii[i]`` is the arc radius at interior vertex ``i + 1``.  The arc is
    tangent to both sides, so the straight parts keep their direction; a
    radius too large for its sides is refused rather than overlapped.
    """
    x = np.asarray(vertices, float)
    radii = list(radii)
    if len(radii) != len(x) - 2:
        raise ValueError(f"{len(x) - 2} interior corner(s), {len(radii)} radii")
    out = [x[0]]
    last_tangent = 0.0                       # length already used on the next side
    for i in range(1, len(x) - 1):
        a, b, c = x[i - 1], x[i], x[i + 1]
        u = (a - b) / np.linalg.norm(a - b)
        w = (c - b) / np.linalg.norm(c - b)
        th = np.arccos(np.clip(u @ w, -1, 1))
        if th < 1e-6 or np.pi - th < 1e-6:
            out.append(b)
            last_tangent = 0.0
            continue
        r = float(radii[i - 1])
        k = r / np.tan(th / 2)               # tangent length on each side
        if k + last_tangent > np.linalg.norm(a - b) or k > np.linalg.norm(c - b):
            raise ValueError(f"radius {r:.0f} m does not fit at corner {i}")
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
    return np.array(out)


def resample(line, spacing) -> np.ndarray:
    """Nodes along ``line`` from its first to its last point.

    ``spacing`` is a number (metres) or a callable ``f(xy) -> metres`` giving
    the wanted spacing at a point.  The first and last nodes are the line's
    ends; a last step shorter than half the local spacing is merged into the
    one before it.

    ``spacing`` is a floor: every step is at least the spacing at both of
    its ends (a step is lengthened until it is), and the remainder at the end
    is spread over all steps, which only lengthens them (review F5); a line
    shorter than its spacing is refused. Non-finite spacing is refused
    (review F19).
    """
    xy = np.asarray(line, float)
    seg = np.linalg.norm(np.diff(xy, axis=0), axis=1)
    s = np.r_[0.0, np.cumsum(seg)]
    total = s[-1]
    if total <= 0:
        raise ValueError("the line has no length")

    def at(t):
        return np.column_stack([np.interp(t, s, xy[:, 0]), np.interp(t, s, xy[:, 1])])

    if not np.isfinite(xy).all():
        raise ValueError("the line has non-finite coordinates")
    f = spacing if callable(spacing) else (lambda p, h=float(spacing): np.full(len(p), h))

    def h_at(t):
        h = float(f(at([min(t, total)]))[0])
        if not np.isfinite(h) or h <= 0:
            raise ValueError(f"spacing must be positive and finite (got {h} at {t:.1f} m)")
        return h

    t = [0.0]
    while True:
        h = h_at(t[-1])
        for _ in range(8):                  # long enough for both ends
            h2 = max(h, h_at(t[-1] + h))
            if h2 <= h:
                break
            h = h2
        if t[-1] + h > total:
            break
        t.append(t[-1] + h)
    t = np.asarray(t)
    if len(t) < 2:
        raise ValueError(f"the line ({total:.0f} m) is shorter than its spacing ({h:.0f} m)")
    # the remainder is shared out: every step grows by the same factor
    t = t * (total / t[-1])
    return at(t)
