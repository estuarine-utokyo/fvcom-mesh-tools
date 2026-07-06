"""Unified read-only QA gate for FVCOM fort.14 meshes (``fmesh-mesh-qa``).

One command evaluates the full acceptance battery from the project
kickoff (§9) as a single pass/fail report:

* FVCOM startup-fatal constraints derived from the model source
  (``docs/fvcom_source_constraints.md``): CCW orientation, no isolated
  elements, the R4 mixed-boundary-element rule (``tge.F:558-581``),
  OBC adjacency chains, manifold boundary topology, plus the hazards
  FVCOM itself never checks (duplicate nodes, orphan nodes, tiny areas,
  the ``{2,1,1}`` ISBCE=2 mis-classification, OBC interior neighbours).
* Manual quality criteria C1 (min angle), C2 (max angle),
  C4 (adjacent-element area change), C5 (valence).
* Open-boundary quality: per-node best-edge perpendicularity and
  node-list ordering.
* Connectivity: single component, everything reachable from the OBC.
* Depth: minimum-depth clip compliance.
* Informational (non-gating by default): narrow-channel w/h, local
  Delaunay fraction, implied external-mode CFL Δt.

Design notes
------------
Angles, areas, and lengths are evaluated in a **local metric space**:
lon/lat meshes are projected with the same flat-earth transform the
diagnostics module uses. FVCOM production builds are CARTESIAN, so the
run-time geometry is metric; raw lon/lat angles can differ by several
degrees at mid-latitudes (the ``fmesh-mesh-quality`` numbers remain
available for raw-space comparisons).

``ISONB`` is recomputed here exactly as ``TRIANGLE_GRID_EDGE`` does:
solid marking from topological boundary edges first, then the
open-node overwrite. The R4 gate therefore matches the model's own
PSTOP condition, independent of the fort.14 land-boundary lists.

The perpendicularity gate uses the **best** (smallest-deviation)
interior edge per OBC node, matching the manual's "one interior edge
normal to the open boundary": in any triangulation most boundary nodes
carry additional diagonal edges that can never be perpendicular.

Repair is out of scope; this module never modifies the mesh.
"""

from __future__ import annotations

import unicodedata
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from fvcom_mesh_tools.algorithms.perpendicularity import boundary_tangents
from fvcom_mesh_tools.diagnostics import (
    _to_metric,
    channel_width_metric,
    face_face_adjacency,
    node_valence,
    unreachable_elements_flag,
)
from fvcom_mesh_tools.io import Fort14Mesh

GRAVITY_M_S2: float = 9.81

#: Depth floor used only inside the implied-Δt estimate so that dry /
#: negative depths do not produce NaN; unrelated to the min-depth gate.
_DT_DEPTH_FLOOR_M: float = 0.1

#: Tolerance (degrees) for the informational local-Delaunay test:
#: an internal edge fails when the two opposite angles sum to more
#: than 180° + this. Right-triangle grids sit exactly at 180°.
_DELAUNAY_TOL_DEG: float = 0.1

# Display labels. Keys are check ids; the CLI defaults to Japanese per
# the kickoff §9 ("Japanese pass/fail table"); ``lang="en"`` is always
# available.
_CATEGORY_LABELS: dict[str, dict[str, str]] = {
    "fvcom": {"ja": "FVCOM必須", "en": "FVCOM-fatal"},
    "obc": {"ja": "開境界", "en": "open bdy"},
    "quality": {"ja": "品質C1-C5", "en": "quality"},
    "connect": {"ja": "連結性", "en": "connectivity"},
    "depth": {"ja": "水深", "en": "depth"},
    "info": {"ja": "情報", "en": "info"},
}

_CHECK_LABELS: dict[str, dict[str, str]] = {
    "node_index_valid": {"ja": "節点参照の整合", "en": "node index integrity"},
    "ccw_all_elements": {"ja": "全要素CCW(反転なし)", "en": "all elements CCW"},
    "no_isolated_elements": {"ja": "孤立要素なし", "en": "no isolated elements"},
    "r4_mixed_boundary": {
        "ja": "R4: 開境界エッジ要素に第3境界節点なし",
        "en": "R4: no extra bdy node on open-edge element",
    },
    "manifold_boundary": {"ja": "境界トポロジ多様体性", "en": "manifold boundary"},
    "no_duplicate_nodes": {"ja": "重複節点なし", "en": "no duplicate nodes"},
    "no_orphan_nodes": {"ja": "孤立節点なし", "en": "no orphan nodes"},
    "no_tiny_area": {"ja": "微小面積要素なし", "en": "no tiny-area elements"},
    "isbce2_authentic": {
        "ja": "ISBCE=2要素の開境界エッジ実在",
        "en": "ISBCE=2 elements have true open edge",
    },
    "obc_on_boundary": {"ja": "開境界節点が境界上", "en": "OBC nodes on boundary"},
    "obc_chain_adjacency": {"ja": "開境界節点の隣接チェーン", "en": "OBC adjacency chains"},
    "obc_interior_neighbor": {"ja": "開境界節点の内部隣接", "en": "OBC interior neighbour"},
    "obc_ordering": {"ja": "開境界節点列の順序", "en": "OBC node ordering"},
    "obc_perpendicularity": {
        "ja": "開境界直交性(節点毎最良エッジ)",
        "en": "OBC perpendicularity (best edge/node)",
    },
    "c1_min_angle": {"ja": "C1: 最小内角", "en": "C1: min interior angle"},
    "c2_max_angle": {"ja": "C2: 最大内角", "en": "C2: max interior angle"},
    "c4_area_change": {"ja": "C4: 隣接要素面積変化", "en": "C4: adjacent area change"},
    "c5_valence": {"ja": "C5: 節点次数", "en": "C5: node valence"},
    "single_component": {"ja": "単一連結成分", "en": "single component"},
    "obc_reachable": {"ja": "全要素が開境界到達可能", "en": "all reachable from OBC"},
    "min_depth_clip": {"ja": "最小水深クリップ", "en": "min-depth clip"},
    "channel_wh": {"ja": "狭水路 w/h", "en": "channel w/h"},
    "delaunay_local": {"ja": "局所Delaunay充足", "en": "local Delaunay"},
    "implied_dt": {"ja": "外部モード安定Δt", "en": "implied external Δt"},
}

_STATUS_LABELS: dict[str, dict[str, str]] = {
    "pass": {"ja": "PASS", "en": "PASS"},
    "fail": {"ja": "FAIL", "en": "FAIL"},
    "info": {"ja": "INFO", "en": "INFO"},
    "skip": {"ja": "SKIP", "en": "SKIP"},
}


@dataclass
class QACheck:
    """One row of the QA report."""

    check_id: str
    category: str
    gate: bool
    passed: bool
    requirement: str
    observed: str
    n_violations: int = 0
    offenders: list[dict[str, Any]] = field(default_factory=list)
    skipped: bool = False
    note: str = ""
    data: dict[str, Any] = field(default_factory=dict)

    @property
    def status(self) -> str:
        if self.skipped:
            return "skip"
        if not self.gate:
            return "info"
        return "pass" if self.passed else "fail"

    def to_dict(self) -> dict[str, Any]:
        return {
            "check_id": self.check_id,
            "category": self.category,
            "gate": bool(self.gate),
            "passed": bool(self.passed),
            "skipped": bool(self.skipped),
            "status": self.status,
            "requirement": self.requirement,
            "observed": self.observed,
            "n_violations": int(self.n_violations),
            "offenders": self.offenders,
            "note": self.note,
            "data": self.data,
        }


@dataclass
class QAReport:
    """Aggregated QA result for one mesh."""

    mesh_name: str
    mesh_path: Path | None
    n_nodes: int
    n_elements: int
    n_open_segments: int
    n_obc_nodes: int
    n_land_segments: int
    coords: str
    params: dict[str, Any]
    checks: list[QACheck]

    @property
    def passed(self) -> bool:
        return all(c.passed for c in self.checks if c.gate and not c.skipped)

    @property
    def n_gate_failed(self) -> int:
        return sum(1 for c in self.checks if c.gate and not c.skipped and not c.passed)

    @property
    def n_gate_total(self) -> int:
        return sum(1 for c in self.checks if c.gate and not c.skipped)

    def to_dict(self) -> dict[str, Any]:
        return {
            "mesh": {
                "name": self.mesh_name,
                "path": str(self.mesh_path) if self.mesh_path else None,
                "n_nodes": int(self.n_nodes),
                "n_elements": int(self.n_elements),
                "n_open_segments": int(self.n_open_segments),
                "n_obc_nodes": int(self.n_obc_nodes),
                "n_land_segments": int(self.n_land_segments),
                "coords": self.coords,
            },
            "params": self.params,
            "checks": [c.to_dict() for c in self.checks],
            "passed": bool(self.passed),
            "n_gate_failed": int(self.n_gate_failed),
            "n_gate_total": int(self.n_gate_total),
        }


# ---------------------------------------------------------------------------
# Geometry / topology primitives
# ---------------------------------------------------------------------------


def detect_coords(nodes: np.ndarray) -> str:
    """``"lonlat"`` when every coordinate fits geographic bounds, else
    ``"metric"``. Projected coordinates (UTM metres) are far outside
    [-360, 360] x [-90, 90].
    """
    if nodes.size == 0:
        return "metric"
    x_ok = np.abs(nodes[:, 0]).max() <= 360.0
    y_ok = np.abs(nodes[:, 1]).max() <= 90.0
    return "lonlat" if (x_ok and y_ok) else "metric"


def _metric_nodes(nodes: np.ndarray, coords: str) -> np.ndarray:
    if coords == "metric":
        return np.asarray(nodes, dtype=np.float64)
    lat0 = float(nodes[:, 1].mean())
    lon0 = float(nodes[:, 0].mean())
    return _to_metric(nodes, lat0=lat0, lon0=lon0)


@dataclass
class _EdgeTopo:
    """Unique undirected edges with use counts and internal-edge pairs."""

    codes_sorted: np.ndarray   # (E,) sorted edge codes u * n_nodes + v
    uv: np.ndarray             # (E, 2) decoded (u, v), aligned to codes_sorted
    counts: np.ndarray         # (E,) number of incident elements
    internal_uv: np.ndarray    # (I, 2)
    internal_pair: np.ndarray  # (I, 2) element ids sharing each internal edge


def _edge_topology(elements: np.ndarray, n_nodes: int) -> _EdgeTopo:
    ne = elements.shape[0]
    raw = np.vstack([
        elements[:, [0, 1]],
        elements[:, [1, 2]],
        elements[:, [2, 0]],
    ])
    raw.sort(axis=1)
    codes = raw[:, 0].astype(np.int64) * n_nodes + raw[:, 1]
    elem_of = np.tile(np.arange(ne, dtype=np.int64), 3)
    order = np.argsort(codes, kind="stable")
    codes_s = codes[order]
    elem_s = elem_of[order]
    uniq, first_idx, counts = np.unique(codes_s, return_index=True, return_counts=True)
    uv = np.column_stack([uniq // n_nodes, uniq % n_nodes]).astype(np.int64)
    internal_mask = counts == 2
    fi = first_idx[internal_mask]
    internal_pair = np.column_stack([elem_s[fi], elem_s[fi + 1]])
    return _EdgeTopo(
        codes_sorted=uniq,
        uv=uv,
        counts=counts,
        internal_uv=uv[internal_mask],
        internal_pair=internal_pair,
    )


def _isin_sorted(values: np.ndarray, sorted_codes: np.ndarray) -> np.ndarray:
    """Membership of ``values`` in the sorted 1-D array ``sorted_codes``."""
    if sorted_codes.size == 0:
        return np.zeros(values.shape, dtype=bool)
    idx = np.searchsorted(sorted_codes, values)
    idx = np.minimum(idx, sorted_codes.size - 1)
    return sorted_codes[idx] == values


@dataclass
class _TriGeometry:
    signed_area: np.ndarray  # (NE,) metric m² (positive = CCW)
    min_angle: np.ndarray    # (NE,) degrees
    max_angle: np.ndarray    # (NE,) degrees
    edge_len: np.ndarray     # (NE, 3) metres


def _tri_geometry(nodes_m: np.ndarray, elements: np.ndarray) -> _TriGeometry:
    p0 = nodes_m[elements[:, 0]]
    p1 = nodes_m[elements[:, 1]]
    p2 = nodes_m[elements[:, 2]]
    twice = (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )
    e0 = np.linalg.norm(p1 - p0, axis=1)
    e1 = np.linalg.norm(p2 - p1, axis=1)
    e2 = np.linalg.norm(p0 - p2, axis=1)

    def _ang(opp: np.ndarray, ea: np.ndarray, eb: np.ndarray) -> np.ndarray:
        denom = 2.0 * ea * eb
        cos = np.where(denom > 0, (ea**2 + eb**2 - opp**2) / np.where(denom > 0, denom, 1.0), 1.0)
        return np.degrees(np.arccos(np.clip(cos, -1.0, 1.0)))

    a0 = _ang(e1, e2, e0)
    a1 = _ang(e2, e0, e1)
    a2 = _ang(e0, e1, e2)
    return _TriGeometry(
        signed_area=0.5 * twice,
        min_angle=np.minimum(np.minimum(a0, a1), a2),
        max_angle=np.maximum(np.maximum(a0, a1), a2),
        edge_len=np.column_stack([e0, e1, e2]),
    )


def compute_isonb(
    n_nodes: int, boundary_uv: np.ndarray, open_boundaries: list[np.ndarray],
) -> np.ndarray:
    """Recompute FVCOM's ISONB node marker (0 interior, 1 solid, 2 open)
    exactly as ``tge.F`` does: solid from topological boundary edges,
    then the OBC-list overwrite.
    """
    isonb = np.zeros(n_nodes, dtype=np.int8)
    if boundary_uv.size:
        isonb[boundary_uv.ravel()] = 1
    for seg in open_boundaries:
        isonb[np.asarray(seg, dtype=np.int64)] = 2
    return isonb


def fvcom_boundary_element_flags(
    mesh: Fort14Mesh, *, topo: _EdgeTopo | None = None,
) -> dict[str, np.ndarray]:
    """ISONB plus the per-element FVCOM boundary-classification hazards.

    Returns a dict with:

    * ``"isonb"`` — ``(NP,)`` int8 node marker (0 interior, 1 solid,
      2 open), recomputed exactly as ``tge.F`` does.
    * ``"isonb_sum"`` — ``(NE,)`` per-element sum of the three node
      markers.
    * ``"r4_mask"`` — ``(NE,)`` bool, True where ``sum > 4``: the
      elements FVCOM PSTOPs on (``tge.F:558-581``).
    * ``"fake_open_mask"`` — ``(NE,)`` bool, True where ``sum == 4``
      but the element has no true open external edge — FVCOM silently
      mis-classifies these as ``ISBCE=2`` / ``EPOR=0``.

    Shared by :func:`run_qa` and by repair drivers that delete or fix
    the offending elements.
    """
    n_nodes = mesh.n_nodes
    ne = mesh.n_elements
    if topo is None:
        topo = _edge_topology(mesh.elements, n_nodes)
    boundary_uv = topo.uv[topo.counts == 1]
    isonb = compute_isonb(n_nodes, boundary_uv, mesh.open_boundaries)
    if ne == 0:
        empty = np.empty(0, dtype=bool)
        return {
            "isonb": isonb,
            "isonb_sum": np.empty(0, dtype=np.int64),
            "r4_mask": empty,
            "fake_open_mask": empty.copy(),
        }
    isonb_sum = isonb[mesh.elements].astype(np.int64).sum(axis=1)
    r4_mask = isonb_sum > 4

    open_edge_codes = np.empty(0, dtype=np.int64)
    if boundary_uv.size:
        both_open = (isonb[boundary_uv[:, 0]] == 2) & (isonb[boundary_uv[:, 1]] == 2)
        oe = boundary_uv[both_open]
        open_edge_codes = np.sort(oe[:, 0] * n_nodes + oe[:, 1])
    cand = isonb_sum == 4
    fake_open_mask = np.zeros(ne, dtype=bool)
    if cand.any():
        ee = np.stack([
            np.sort(mesh.elements[:, [0, 1]], axis=1),
            np.sort(mesh.elements[:, [1, 2]], axis=1),
            np.sort(mesh.elements[:, [2, 0]], axis=1),
        ], axis=1)
        elem_codes = ee[..., 0].astype(np.int64) * n_nodes + ee[..., 1]
        has_open_edge = _isin_sorted(
            elem_codes.ravel(), open_edge_codes,
        ).reshape(ne, 3).any(axis=1)
        fake_open_mask = cand & ~has_open_edge
    return {
        "isonb": isonb,
        "isonb_sum": isonb_sum,
        "r4_mask": r4_mask,
        "fake_open_mask": fake_open_mask,
    }


# ---------------------------------------------------------------------------
# Offender record helpers (locations are reported in file coordinates)
# ---------------------------------------------------------------------------


def _elem_offenders(
    idx: np.ndarray, mesh: Fort14Mesh, values: np.ndarray | None,
    *, limit: int, value_key: str = "value",
) -> list[dict[str, Any]]:
    take = idx[:limit]
    cent = mesh.nodes[mesh.elements[take]].mean(axis=1)
    out = []
    for k, i in enumerate(take):
        rec: dict[str, Any] = {
            "kind": "element", "id": int(i),
            "x": float(cent[k, 0]), "y": float(cent[k, 1]),
        }
        if values is not None:
            rec[value_key] = float(values[i])
        out.append(rec)
    return out


def _node_offenders(
    idx: np.ndarray, mesh: Fort14Mesh, values: np.ndarray | None,
    *, limit: int, value_key: str = "value",
) -> list[dict[str, Any]]:
    take = idx[:limit]
    out = []
    for i in take:
        rec: dict[str, Any] = {
            "kind": "node", "id": int(i),
            "x": float(mesh.nodes[i, 0]), "y": float(mesh.nodes[i, 1]),
        }
        if values is not None:
            rec[value_key] = float(values[i])
        out.append(rec)
    return out


def _sort_desc(idx: np.ndarray, values: np.ndarray) -> np.ndarray:
    return idx[np.argsort(values[idx])[::-1]]


def _sort_asc(idx: np.ndarray, values: np.ndarray) -> np.ndarray:
    return idx[np.argsort(values[idx])]


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def run_qa(
    mesh: Fort14Mesh,
    *,
    name: str | None = None,
    path: Path | None = None,
    min_angle_deg: float = 30.0,
    max_angle_deg: float = 130.0,
    max_area_change: float = 0.5,
    max_valence: int = 8,
    min_depth_m: float = 2.0,
    max_obc_perp_dev_deg: float = 20.0,
    coords: str = "auto",
    channel_check: bool = True,
    min_channel_wh_gate: float | None = None,
    min_dt_s: float | None = None,
    tiny_area_m2: float = 1e-6,
    duplicate_tol_m: float = 1e-3,
    max_offenders: int = 5,
    land_solid_shp: Path | None = None,
    utm_epsg: int = 32654,
    land_interior_m: float = 200.0,
) -> QAReport:
    """Run the full QA battery and return a :class:`QAReport`.

    ``min_channel_wh_gate`` and ``min_dt_s`` turn the corresponding
    informational checks into gates when supplied. ``coords`` may be
    ``"auto"`` (default), ``"lonlat"``, or ``"metric"``.
    """
    n_nodes = mesh.n_nodes
    ne = mesh.n_elements
    obc_nodes_all = (
        np.unique(np.concatenate([np.asarray(s, dtype=np.int64)
                                  for s in mesh.open_boundaries]))
        if mesh.open_boundaries else np.empty(0, dtype=np.int64)
    )
    params: dict[str, Any] = {
        "min_angle_deg": min_angle_deg,
        "max_angle_deg": max_angle_deg,
        "max_area_change": max_area_change,
        "max_valence": max_valence,
        "min_depth_m": min_depth_m,
        "max_obc_perp_dev_deg": max_obc_perp_dev_deg,
        "min_channel_wh_gate": min_channel_wh_gate,
        "min_dt_s": min_dt_s,
        "tiny_area_m2": tiny_area_m2,
        "duplicate_tol_m": duplicate_tol_m,
    }
    checks: list[QACheck] = []

    def _report(coords_resolved: str) -> QAReport:
        return QAReport(
            mesh_name=name or (path.name if path else "mesh"),
            mesh_path=path,
            n_nodes=n_nodes,
            n_elements=ne,
            n_open_segments=len(mesh.open_boundaries),
            n_obc_nodes=int(obc_nodes_all.size),
            n_land_segments=len(mesh.land_boundaries),
            coords=coords_resolved,
            params=params,
            checks=checks,
        )

    # -- 0. node index integrity (must hold before any geometry work) ------
    bad_idx = int(((mesh.elements < 0) | (mesh.elements >= n_nodes)).sum()) if ne else 0
    seg_ids = [np.asarray(s, dtype=np.int64) for s in mesh.open_boundaries] + [
        np.asarray(s, dtype=np.int64) for _ib, s in mesh.land_boundaries
    ]
    for s in seg_ids:
        if s.size:
            bad_idx += int(((s < 0) | (s >= n_nodes)).sum())
    checks.append(QACheck(
        "node_index_valid", "fvcom", True, bad_idx == 0,
        f"0 <= id < {n_nodes}", f"out-of-range refs = {bad_idx}", bad_idx,
    ))
    if bad_idx > 0 or ne == 0 or n_nodes == 0:
        if ne == 0 or n_nodes == 0:
            checks.append(QACheck(
                "no_isolated_elements", "fvcom", True, False,
                "NE > 0", f"NE = {ne}, NP = {n_nodes}", 1,
                note="empty mesh — remaining checks skipped",
            ))
        return _report(coords if coords != "auto" else "metric")

    coords_resolved = coords if coords != "auto" else detect_coords(mesh.nodes)
    nodes_m = _metric_nodes(mesh.nodes, coords_resolved)
    geom = _tri_geometry(nodes_m, mesh.elements)
    topo = _edge_topology(mesh.elements, n_nodes)
    boundary_uv = topo.uv[topo.counts == 1]
    valence = node_valence(mesh.elements, n_nodes=n_nodes)
    bflags = fvcom_boundary_element_flags(mesh, topo=topo)
    isonb_sum = bflags["isonb_sum"]
    adj = face_face_adjacency(mesh.elements)

    # -- FVCOM-fatal geometry / topology -----------------------------------

    # CCW orientation (FVCOM checks element #1 only; we check all).
    flipped = np.where(geom.signed_area < 0)[0]
    checks.append(QACheck(
        "ccw_all_elements", "fvcom", True, flipped.size == 0,
        "flipped = 0", f"flipped = {flipped.size}", int(flipped.size),
        offenders=_elem_offenders(
            _sort_asc(flipped, geom.signed_area), mesh, geom.signed_area,
            limit=max_offenders, value_key="signed_area_m2"),
    ))

    # Isolated elements (tge.F PSTOP: element with no neighbours).
    deg = np.asarray(adj.sum(axis=1)).ravel()
    isolated = np.where(deg == 0)[0] if ne > 1 else np.empty(0, dtype=np.int64)
    checks.append(QACheck(
        "no_isolated_elements", "fvcom", True, isolated.size == 0,
        "isolated = 0", f"isolated = {isolated.size}", int(isolated.size),
        offenders=_elem_offenders(isolated, mesh, None, limit=max_offenders),
    ))

    # R4: SUM(ISONB over element nodes) > 4 is a PSTOP (tge.F:558-581).
    r4_viol = np.where(bflags["r4_mask"])[0]
    checks.append(QACheck(
        "r4_mixed_boundary", "fvcom", True, r4_viol.size == 0,
        "sum(ISONB) <= 4", f"violations = {r4_viol.size}", int(r4_viol.size),
        offenders=_elem_offenders(
            _sort_desc(r4_viol, isonb_sum.astype(np.float64)), mesh,
            isonb_sum.astype(np.float64), limit=max_offenders,
            value_key="isonb_sum"),
    ))

    # Manifold boundary: no node with >2 boundary edges (tge.F ISONB
    # ordering FATAL_ERROR), no edge shared by >2 elements.
    bnd_edge_count = np.zeros(n_nodes, dtype=np.int64)
    if boundary_uv.size:
        np.add.at(bnd_edge_count, boundary_uv.ravel(), 1)
    pinch_nodes = np.where(bnd_edge_count > 2)[0]
    over_edges = topo.uv[topo.counts > 2]
    n_manifold_viol = int(pinch_nodes.size + over_edges.shape[0])
    manifold_off = _node_offenders(
        pinch_nodes, mesh, bnd_edge_count.astype(np.float64),
        limit=max_offenders, value_key="boundary_edge_count")
    for u, v in over_edges[: max(0, max_offenders - len(manifold_off))]:
        pos = int(np.searchsorted(topo.codes_sorted, int(u) * n_nodes + int(v)))
        manifold_off.append({
            "kind": "edge", "id": [int(u), int(v)],
            "x": float(mesh.nodes[[u, v], 0].mean()),
            "y": float(mesh.nodes[[u, v], 1].mean()),
            "n_incident_elements": int(topo.counts[pos]),
        })
    checks.append(QACheck(
        "manifold_boundary", "fvcom", True, n_manifold_viol == 0,
        "pinch nodes = 0, >2-elem edges = 0",
        f"pinch nodes = {pinch_nodes.size}, over-shared edges = {over_edges.shape[0]}",
        n_manifold_viol, offenders=manifold_off,
    ))

    # Duplicate nodes (silent in FVCOM).
    dup_pairs = sorted(cKDTree(nodes_m).query_pairs(r=duplicate_tol_m))
    dup_off = [
        {
            "kind": "node_pair", "id": [int(a), int(b)],
            "x": float(mesh.nodes[a, 0]), "y": float(mesh.nodes[a, 1]),
        }
        for a, b in dup_pairs[:max_offenders]
    ]
    checks.append(QACheck(
        "no_duplicate_nodes", "fvcom", True, len(dup_pairs) == 0,
        f"pairs closer than {duplicate_tol_m} m = 0",
        f"coincident pairs = {len(dup_pairs)}", len(dup_pairs), offenders=dup_off,
    ))

    # Orphan nodes (silent in FVCOM; degenerate control volumes).
    orphan = np.where(valence == 0)[0]
    checks.append(QACheck(
        "no_orphan_nodes", "fvcom", True, orphan.size == 0,
        "unreferenced nodes = 0", f"orphans = {orphan.size}", int(orphan.size),
        offenders=_node_offenders(orphan, mesh, None, limit=max_offenders),
    ))

    # Tiny / zero areas (only a WARNING in cell_area.F; gate here).
    tiny = np.where(np.abs(geom.signed_area) < tiny_area_m2)[0]
    checks.append(QACheck(
        "no_tiny_area", "fvcom", True, tiny.size == 0,
        f"|area| >= {tiny_area_m2:g} m2", f"tiny = {tiny.size}", int(tiny.size),
        offenders=_elem_offenders(
            _sort_asc(tiny, np.abs(geom.signed_area)), mesh,
            np.abs(geom.signed_area), limit=max_offenders, value_key="area_m2"),
    ))

    # ISBCE=2 authenticity: every element FVCOM would classify open
    # (sum == 4) must own a true boundary edge with both endpoints open.
    fake_open = np.where(bflags["fake_open_mask"])[0]
    checks.append(QACheck(
        "isbce2_authentic", "fvcom", True, fake_open.size == 0,
        "each ISBCE=2 element has an open external edge",
        f"mis-classified = {fake_open.size}", int(fake_open.size),
        offenders=_elem_offenders(fake_open, mesh, None, limit=max_offenders),
    ))

    # -- Open-boundary checks ----------------------------------------------

    has_obc = obc_nodes_all.size > 0
    obc_mask = np.zeros(n_nodes, dtype=bool)
    obc_mask[obc_nodes_all] = True

    boundary_node_set = np.zeros(n_nodes, dtype=bool)
    if boundary_uv.size:
        boundary_node_set[boundary_uv.ravel()] = True

    if not has_obc:
        for cid in ("obc_on_boundary", "obc_chain_adjacency",
                    "obc_interior_neighbor", "obc_ordering",
                    "obc_perpendicularity"):
            checks.append(QACheck(
                cid, "obc", True, True, "-", "-", 0, skipped=True,
                note="no open boundary in mesh",
            ))
    else:
        off_bdy = obc_nodes_all[~boundary_node_set[obc_nodes_all]]
        checks.append(QACheck(
            "obc_on_boundary", "obc", True, off_bdy.size == 0,
            "every OBC node on mesh boundary",
            f"interior OBC nodes = {off_bdy.size}", int(off_bdy.size),
            offenders=_node_offenders(off_bdy, mesh, None, limit=max_offenders),
        ))

        # Chain adjacency (mod_obcs.F PSTOP): each OBC node needs >= 1
        # OBC neighbour through a mesh edge.
        u_obc = obc_mask[topo.uv[:, 0]]
        v_obc = obc_mask[topo.uv[:, 1]]
        both = u_obc & v_obc
        has_obc_nbr = np.zeros(n_nodes, dtype=bool)
        if both.any():
            has_obc_nbr[topo.uv[both].ravel()] = True
        lonely = obc_nodes_all[~has_obc_nbr[obc_nodes_all]]
        checks.append(QACheck(
            "obc_chain_adjacency", "obc", True, lonely.size == 0,
            "every OBC node adjacent to another OBC node",
            f"lonely OBC nodes = {lonely.size}", int(lonely.size),
            offenders=_node_offenders(lonely, mesh, None, limit=max_offenders),
        ))

        # Interior neighbour (NEXT_OBC = 0 hazard, silent in FVCOM).
        one_side = u_obc ^ v_obc
        has_int_nbr = np.zeros(n_nodes, dtype=bool)
        if one_side.any():
            side_uv = topo.uv[one_side]
            obc_end = np.where(obc_mask[side_uv[:, 0]], side_uv[:, 0], side_uv[:, 1])
            has_int_nbr[obc_end] = True
        necked = obc_nodes_all[~has_int_nbr[obc_nodes_all]]
        checks.append(QACheck(
            "obc_interior_neighbor", "obc", True, necked.size == 0,
            "every OBC node has a non-OBC neighbour",
            f"necked OBC nodes = {necked.size}", int(necked.size),
            offenders=_node_offenders(necked, mesh, None, limit=max_offenders),
        ))

        # Ordering: each open segment must be a duplicate-free walk
        # along existing mesh edges (required for a clean _obc.dat).
        n_break = 0
        n_dup = 0
        order_off: list[dict[str, Any]] = []
        for k, seg in enumerate(mesh.open_boundaries):
            seg = np.asarray(seg, dtype=np.int64)
            n_dup += int(seg.size - np.unique(seg).size)
            if seg.size < 2:
                continue
            a = np.minimum(seg[:-1], seg[1:])
            b = np.maximum(seg[:-1], seg[1:])
            present = _isin_sorted(a * n_nodes + b, topo.codes_sorted)
            for j in np.where(~present)[0]:
                n_break += 1
                if len(order_off) < max_offenders:
                    order_off.append({
                        "kind": "obc_pair", "segment": int(k),
                        "id": [int(seg[j]), int(seg[j + 1])],
                        "x": float(mesh.nodes[seg[j], 0]),
                        "y": float(mesh.nodes[seg[j], 1]),
                    })
        checks.append(QACheck(
            "obc_ordering", "obc", True, (n_break + n_dup) == 0,
            "consecutive nodes share a mesh edge, no repeats",
            f"non-adjacent pairs = {n_break}, duplicates = {n_dup}",
            n_break + n_dup, offenders=order_off,
        ))

        # Perpendicularity: best incident edge per OBC node, in metric
        # space, against the local along-boundary tangent.
        min_dev = np.full(n_nodes, np.inf)
        for seg in mesh.open_boundaries:
            seg = np.asarray(seg, dtype=np.int64)
            if seg.size < 2:
                continue
            tang = boundary_tangents(nodes_m[seg])
            in_seg = np.zeros(n_nodes, dtype=bool)
            in_seg[seg] = True
            inv = np.full(n_nodes, -1, dtype=np.int64)
            inv[seg] = np.arange(seg.size)
            a_in = in_seg[topo.uv[:, 0]]
            b_in = in_seg[topo.uv[:, 1]]
            inc = topo.uv[a_in ^ b_in]
            if inc.size == 0:
                continue
            inc_a = in_seg[inc[:, 0]]
            bnode = np.where(inc_a, inc[:, 0], inc[:, 1])
            onode = np.where(inc_a, inc[:, 1], inc[:, 0])
            vec = nodes_m[onode] - nodes_m[bnode]
            norm = np.linalg.norm(vec, axis=1, keepdims=True)
            vec = vec / np.where(norm == 0, 1.0, norm)
            cos = np.clip(np.abs((vec * tang[inv[bnode]]).sum(axis=1)), 0.0, 1.0)
            dev = 90.0 - np.degrees(np.arccos(cos))
            np.minimum.at(min_dev, bnode, dev)
        obc_dev = min_dev[obc_nodes_all]
        finite = np.isfinite(obc_dev)
        worst = float(obc_dev[finite].max()) if finite.any() else float("inf")
        perp_viol = obc_nodes_all[(~finite) | (obc_dev > max_obc_perp_dev_deg)]
        checks.append(QACheck(
            "obc_perpendicularity", "obc", True, perp_viol.size == 0,
            f"best-edge deviation <= {max_obc_perp_dev_deg:g} deg / node",
            f"worst = {worst:.1f} deg, mean = "
            f"{float(obc_dev[finite].mean()) if finite.any() else float('nan'):.1f} deg",
            int(perp_viol.size),
            offenders=_node_offenders(
                _sort_desc(perp_viol, np.where(np.isfinite(min_dev), min_dev, 1e9)),
                mesh, min_dev, limit=max_offenders, value_key="deviation_deg"),
            data={"worst_deviation_deg": worst},
        ))

    # -- Manual quality criteria C1 / C2 / C4 / C5 -------------------------

    c1 = np.where(geom.min_angle < min_angle_deg)[0]
    checks.append(QACheck(
        "c1_min_angle", "quality", True, c1.size == 0,
        f">= {min_angle_deg:g} deg",
        f"min = {float(geom.min_angle.min()):.2f} deg, violations = {c1.size}",
        int(c1.size),
        offenders=_elem_offenders(
            _sort_asc(c1, geom.min_angle), mesh, geom.min_angle,
            limit=max_offenders, value_key="min_angle_deg"),
        data={"min_angle_deg": float(geom.min_angle.min())},
    ))

    c2 = np.where(geom.max_angle > max_angle_deg)[0]
    checks.append(QACheck(
        "c2_max_angle", "quality", True, c2.size == 0,
        f"<= {max_angle_deg:g} deg",
        f"max = {float(geom.max_angle.max()):.2f} deg, violations = {c2.size}",
        int(c2.size),
        offenders=_elem_offenders(
            _sort_desc(c2, geom.max_angle), mesh, geom.max_angle,
            limit=max_offenders, value_key="max_angle_deg"),
        data={"max_angle_deg": float(geom.max_angle.max())},
    ))

    areas = np.abs(geom.signed_area)
    if topo.internal_pair.shape[0]:
        a_i = areas[topo.internal_pair[:, 0]]
        a_j = areas[topo.internal_pair[:, 1]]
        larger = np.maximum(a_i, a_j)
        area_change = (larger - np.minimum(a_i, a_j)) / np.maximum(larger, 1e-30)
    else:
        area_change = np.empty(0)
    c4 = np.where(area_change > max_area_change)[0]
    c4_worst = float(area_change.max()) if area_change.size else 0.0
    c4_off = []
    for j in _sort_desc(c4, area_change)[:max_offenders]:
        u, v = topo.internal_uv[j]
        c4_off.append({
            "kind": "edge", "id": [int(u), int(v)],
            "elements": [int(topo.internal_pair[j, 0]), int(topo.internal_pair[j, 1])],
            "x": float(mesh.nodes[[u, v], 0].mean()),
            "y": float(mesh.nodes[[u, v], 1].mean()),
            "area_change": float(area_change[j]),
        })
    checks.append(QACheck(
        "c4_area_change", "quality", True, c4.size == 0,
        f"<= {max_area_change:g}",
        f"max = {c4_worst:.3f}, violations = {c4.size}", int(c4.size),
        offenders=c4_off, data={"max_area_change": c4_worst},
    ))

    c5 = np.where(valence > max_valence)[0]
    checks.append(QACheck(
        "c5_valence", "quality", True, c5.size == 0,
        f"<= {max_valence}",
        f"max = {int(valence.max())}, violations = {c5.size}", int(c5.size),
        offenders=_node_offenders(
            _sort_desc(c5, valence.astype(np.float64)), mesh,
            valence.astype(np.float64), limit=max_offenders, value_key="valence"),
        data={"max_valence": int(valence.max())},
    ))

    # -- Connectivity -------------------------------------------------------

    n_comp, labels = connected_components(adj, directed=False, return_labels=True)
    sizes = np.bincount(labels, minlength=int(n_comp))
    n_disjoint = ne - int(sizes.max()) if sizes.size else 0
    disjoint_idx = np.where(labels != int(sizes.argmax()))[0] if n_disjoint else np.empty(0, int)
    checks.append(QACheck(
        "single_component", "connect", True, n_disjoint == 0,
        "disjoint elements = 0",
        f"components = {int(n_comp)}, disjoint = {n_disjoint}", int(n_disjoint),
        offenders=_elem_offenders(disjoint_idx, mesh, None, limit=max_offenders),
    ))

    if has_obc:
        unreach = np.where(unreachable_elements_flag(adj, mesh, labels))[0]
        checks.append(QACheck(
            "obc_reachable", "connect", True, unreach.size == 0,
            "unreachable elements = 0",
            f"unreachable = {unreach.size}", int(unreach.size),
            offenders=_elem_offenders(unreach, mesh, None, limit=max_offenders),
        ))
    else:
        checks.append(QACheck(
            "obc_reachable", "connect", True, True, "-", "-", 0,
            skipped=True, note="no open boundary in mesh",
        ))

    # -- Depth ---------------------------------------------------------------

    shallow = np.where(mesh.depths < min_depth_m)[0]
    checks.append(QACheck(
        "min_depth_clip", "depth", True, shallow.size == 0,
        f">= {min_depth_m:g} m",
        f"min = {float(mesh.depths.min()):.2f} m, violations = {shallow.size}",
        int(shallow.size),
        offenders=_node_offenders(
            _sort_asc(shallow, mesh.depths), mesh, mesh.depths,
            limit=max_offenders, value_key="depth_m"),
        data={"min_depth_m": float(mesh.depths.min())},
    ))

    # -- Informational -------------------------------------------------------

    # Narrow channels (advisory; connectivity is the hard gate). The
    # cross-polyline width estimate is unreliable near open/land
    # junctions, so this only gates when the caller asks.
    if channel_check and (mesh.open_boundaries or mesh.land_boundaries):
        ch = channel_width_metric(mesh, coords=coords_resolved)
        ratio = ch["w_h_ratio"]
        finite_r = ratio[np.isfinite(ratio)]
        n_lt1 = int((finite_r < 1.0).sum())
        n_lt3 = int((finite_r < 3.0).sum())
        gate_ch = min_channel_wh_gate is not None
        thr = float(min_channel_wh_gate) if gate_ch else 1.0
        ch_viol = np.where(np.isfinite(ratio) & (ratio < thr))[0]
        checks.append(QACheck(
            "channel_wh", "info", gate_ch,
            (ch_viol.size == 0) if gate_ch else True,
            f">= {thr:g} (gate)" if gate_ch else "advisory",
            f"w/h < 1: {n_lt1}, w/h < 3: {n_lt3}",
            int(ch_viol.size) if gate_ch else n_lt1,
            offenders=_elem_offenders(
                _sort_asc(ch_viol, ratio), mesh, ratio,
                limit=max_offenders, value_key="w_h_ratio"),
            data={"n_wh_lt_1": n_lt1, "n_wh_lt_3": n_lt3},
        ))
    else:
        checks.append(QACheck(
            "channel_wh", "info", False, True, "advisory", "-", 0,
            skipped=True,
            note="skipped (disabled or no boundary polylines)",
        ))

    # Local Delaunay (advisory): opposite angles of an internal edge
    # must not sum past 180 deg.
    if topo.internal_pair.shape[0]:
        tri_sum = mesh.elements.sum(axis=1).astype(np.int64)
        uv_sum = topo.internal_uv.sum(axis=1)
        opp1 = tri_sum[topo.internal_pair[:, 0]] - uv_sum
        opp2 = tri_sum[topo.internal_pair[:, 1]] - uv_sum

        def _opp_angle(opp: np.ndarray) -> np.ndarray:
            vu = nodes_m[topo.internal_uv[:, 0]] - nodes_m[opp]
            vv = nodes_m[topo.internal_uv[:, 1]] - nodes_m[opp]
            nu = np.linalg.norm(vu, axis=1)
            nv = np.linalg.norm(vv, axis=1)
            cos = (vu * vv).sum(axis=1) / np.where(nu * nv > 0, nu * nv, 1.0)
            return np.degrees(np.arccos(np.clip(cos, -1.0, 1.0)))

        ang_sum = _opp_angle(opp1) + _opp_angle(opp2)
        n_nondel = int((ang_sum > 180.0 + _DELAUNAY_TOL_DEG).sum())
        frac = n_nondel / topo.internal_pair.shape[0]
        checks.append(QACheck(
            "delaunay_local", "info", False, True, "advisory",
            f"non-Delaunay internal edges = {n_nondel} ({frac * 100:.2f}%)",
            n_nondel, data={"n_non_delaunay": n_nondel, "fraction": frac},
        ))
    else:
        checks.append(QACheck(
            "delaunay_local", "info", False, True, "advisory", "-", 0,
            skipped=True, note="no internal edges",
        ))

    # Implied external-mode CFL dt (L_min / sqrt(g H_max), no safety
    # factor). Gates only when min_dt_s is supplied.
    l_min = geom.edge_len.min(axis=1)
    h_elem = np.maximum(mesh.depths[mesh.elements].max(axis=1), _DT_DEPTH_FLOOR_M)
    dt_elem = l_min / np.sqrt(GRAVITY_M_S2 * h_elem)
    dt_min = float(dt_elem.min())
    gate_dt = min_dt_s is not None
    dt_viol = np.where(dt_elem < float(min_dt_s))[0] if gate_dt else np.empty(0, int)
    checks.append(QACheck(
        "implied_dt", "info", gate_dt,
        (dt_viol.size == 0) if gate_dt else True,
        f">= {float(min_dt_s):g} s (gate)" if gate_dt else "advisory",
        f"min dt = {dt_min:.2f} s (elem {int(dt_elem.argmin())})",
        int(dt_viol.size),
        offenders=_elem_offenders(
            _sort_asc(dt_viol, dt_elem), mesh, dt_elem,
            limit=max_offenders, value_key="dt_s"),
        data={"dt_min_s": dt_min, "worst_element": int(dt_elem.argmin())},
    ))

    if land_solid_shp is not None and Path(land_solid_shp).exists():
        # mesh-over-land tripwire (I11/J11 wetland incident): no
        # element centroid may lie in the DEEP interior (more than
        # land_interior_m inside) of the solid pre-prep land;
        # near-boundary overlap from coastline conformity and
        # r_open-scale river openings is tolerated by construction.
        import geopandas as _gpd
        import shapely as _sh
        from pyproj import Transformer as _Tr
        from shapely.ops import unary_union as _uu

        _g = _gpd.read_file(land_solid_shp)
        if _g.crs is None:
            _g = _g.set_crs(4326)
        _solid = _uu(list(_g.to_crs(utm_epsg).geometry)).buffer(
            -float(land_interior_m))
        _cen = mesh.nodes[mesh.elements].mean(axis=1)
        if coords_resolved == "lonlat":
            _tr = _Tr.from_crs("EPSG:4326", f"EPSG:{utm_epsg}",
                               always_xy=True)
            _cx, _cy = _tr.transform(_cen[:, 0], _cen[:, 1])
        else:
            _cx, _cy = _cen[:, 0], _cen[:, 1]
        _bad = np.where(_sh.contains_xy(_solid, _cx, _cy))[0]
        checks.append(QACheck(
            "land_overlap", "structure", True, len(_bad) == 0,
            f"no element centroid > {land_interior_m:g} m inside "
            "solid pre-prep land",
            f"{len(_bad)} elements over solid land",
            int(len(_bad)),
            offenders=[{"element": int(v)} for v in
                       _bad[:max_offenders]],
        ))

    return _report(coords_resolved)


# ---------------------------------------------------------------------------
# Report rendering (Japanese default per kickoff §9; English available)
# ---------------------------------------------------------------------------


def _disp_width(s: str) -> int:
    return sum(2 if unicodedata.east_asian_width(ch) in ("W", "F") else 1 for ch in s)


def _pad(s: str, width: int, align: str = "left") -> str:
    gap = max(0, width - _disp_width(s))
    return (" " * gap + s) if align == "right" else (s + " " * gap)


def format_report(report: QAReport, *, lang: str = "ja") -> str:
    """Render the pass/fail table plus worst-offender details."""
    if lang not in ("ja", "en"):
        raise ValueError(f"unsupported lang: {lang!r}")

    if lang == "ja":
        head = (
            f"=== fmesh-mesh-qa: {report.mesh_name} ===",
            f"節点数 {report.n_nodes:,}  要素数 {report.n_elements:,}  "
            f"開境界 {report.n_open_segments} 区間 ({report.n_obc_nodes} 節点)  "
            f"陸境界 {report.n_land_segments} 区間  座標系 {report.coords}",
        )
        cols = ("判定", "区分", "チェック", "基準", "実測", "違反")
    else:
        head = (
            f"=== fmesh-mesh-qa: {report.mesh_name} ===",
            f"nodes {report.n_nodes:,}  elements {report.n_elements:,}  "
            f"open segs {report.n_open_segments} ({report.n_obc_nodes} nodes)  "
            f"land segs {report.n_land_segments}  coords {report.coords}",
        )
        cols = ("status", "category", "check", "requirement", "observed", "viol")

    rows: list[tuple[str, ...]] = []
    for c in report.checks:
        rows.append((
            _STATUS_LABELS[c.status][lang],
            _CATEGORY_LABELS[c.category][lang],
            _CHECK_LABELS[c.check_id][lang],
            c.requirement if not c.skipped else "-",
            c.observed if not c.skipped else (c.note or "-"),
            str(c.n_violations) if not c.skipped else "-",
        ))

    widths = [
        max(_disp_width(cols[i]), max(_disp_width(r[i]) for r in rows))
        for i in range(len(cols))
    ]
    lines = list(head)
    lines.append("")
    lines.append("  ".join(_pad(cols[i], widths[i]) for i in range(len(cols))))
    lines.append("  ".join("-" * widths[i] for i in range(len(cols))))
    for r in rows:
        aligned = [
            _pad(r[i], widths[i], "right" if i == len(cols) - 1 else "left")
            for i in range(len(cols))
        ]
        lines.append("  ".join(aligned))

    n_info = sum(1 for c in report.checks if not c.gate)
    if lang == "ja":
        verdict = "PASS(全合格)" if report.passed else "FAIL(不合格あり)"
        lines.append("")
        lines.append(
            f"総合判定: {verdict} — 不合格 {report.n_gate_failed} / "
            f"ゲート対象 {report.n_gate_total}(参考情報 {n_info} 件は判定対象外)"
        )
    else:
        verdict = "PASS" if report.passed else "FAIL"
        lines.append("")
        lines.append(
            f"overall: {verdict} — failed {report.n_gate_failed} / "
            f"{report.n_gate_total} gated checks ({n_info} informational)"
        )

    failed = [c for c in report.checks if c.gate and not c.skipped
              and not c.passed and c.offenders]
    if failed:
        lines.append("")
        lines.append("worst offenders:" if lang == "en" else "違反箇所(最悪順):")
        for c in failed:
            lines.append(f"  [{_CHECK_LABELS[c.check_id][lang]}]")
            for o in c.offenders:
                extra = ", ".join(
                    f"{k}={v:.4g}" if isinstance(v, float) else f"{k}={v}"
                    for k, v in o.items() if k not in ("kind", "id", "x", "y")
                )
                lines.append(
                    f"    - {o['kind']} {o['id']} @ ({o['x']:.6g}, {o['y']:.6g})"
                    + (f": {extra}" if extra else "")
                )
    return "\n".join(lines)


__all__ = [
    "GRAVITY_M_S2",
    "QACheck",
    "QAReport",
    "compute_isonb",
    "detect_coords",
    "format_report",
    "fvcom_boundary_element_flags",
    "run_qa",
]
