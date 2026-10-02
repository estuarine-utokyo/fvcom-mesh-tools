"""Writers for FVCOM's native ASCII grid-input files and SMS ``.2dm``.

Formats follow the FVCOM 5.1 cold-start readers exactly as documented
in ``docs/fvcom_source_constraints.md`` (derived from
``mod_input.F``); all are free-format, 1-indexed, keyword-header
files:

* ``casename_grd.dat``  — ``Node Number`` / ``Cell Number`` headers,
  connectivity rows ``CELL# N1 N2 N3`` (**CCW in the file** — FVCOM
  swaps columns to its internal CW convention on read), then
  coordinate rows ``NODE# X Y``.
* ``casename_dep.dat``  — ``Node Number`` header, rows ``X Y H``
  (H positive down).
* ``casename_obc.dat``  — ``OBC Node Number`` header, rows
  ``OBCNODE# GLOBALNODE# TYPE`` (type 1-10; file order is free but
  this writer emits the along-boundary order, which downstream tools
  expect).
* ``casename_cor.dat``  — ``Node Number`` header, rows ``X Y COR``
  (for CARTESIAN builds COR is the latitude in degrees used to
  compute the Coriolis parameter).
* ``casename_spg.dat``  — ``Sponge Node Number`` header, rows
  ``GLOBALNODE# RADIUS DAMPING`` (0 nodes = no sponge).
* ``casename.2dm``      — SMS interoperability export (``E3T`` /
  ``ND`` cards plus one ``NS`` nodestring per open-boundary segment).

Writers validate the format-level invariants the FVCOM reader
enforces or silently mis-handles: all-CCW connectivity (the reader
checks element #1 only; a mixed-orientation mesh would run with
corrupt geometry) and the absence of unreferenced nodes (a trailing
orphan breaks the reader's max-index check; any orphan yields a
degenerate control volume). Quality/topology acceptance is
``fmesh-mesh-qa``'s job, not the writers'.
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import Sequence

import numpy as np

from fvcom_mesh_tools.io.fort14 import Fort14Mesh


def _num(v) -> str:
    """The shortest text that reads back as the same double.

    Fixed 8/6 decimals cut 1000.123456789 to 1000.12345679, so a base mesh
    did not survive export bit for bit (review of the extend tools, round 2
    F8); a value read from such a file prints back unchanged.
    """
    return repr(float(v))


def _signed_areas(mesh: Fort14Mesh) -> np.ndarray:
    p0 = mesh.nodes[mesh.elements[:, 0]]
    p1 = mesh.nodes[mesh.elements[:, 1]]
    p2 = mesh.nodes[mesh.elements[:, 2]]
    return 0.5 * (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )


def _validate_for_export(mesh: Fort14Mesh) -> None:
    if mesh.n_elements == 0 or mesh.n_nodes == 0:
        raise ValueError("cannot export an empty mesh")
    _check_exportable(mesh)
    sa = _signed_areas(mesh)
    n_bad = int((sa <= 0).sum())
    if n_bad:
        raise ValueError(
            f"{n_bad} element(s) are not CCW (flipped or zero-area); "
            "FVCOM requires counter-clockwise connectivity in _grd.dat "
            "— repair before export"
        )
    used = np.unique(mesh.elements)
    if used.size != mesh.n_nodes:
        raise ValueError(
            f"{mesh.n_nodes - used.size} unreferenced node(s) present; "
            "run fvcom_mesh_tools.mesh_clean.compact_nodes before export"
        )


def write_grd(mesh: Fort14Mesh, path: str | Path) -> Path:
    """Write ``casename_grd.dat`` (connectivity block, then coords)."""
    _validate_for_export(mesh)
    path = Path(path).resolve()
    with path.open("w") as f:
        f.write(f"Node Number = {mesh.n_nodes}\n")
        f.write(f"Cell Number = {mesh.n_elements}\n")
        for i in range(mesh.n_elements):
            n0, n1, n2 = mesh.elements[i]
            f.write(f"{i + 1} {int(n0) + 1} {int(n1) + 1} {int(n2) + 1}\n")
        for i in range(mesh.n_nodes):
            x = _num(mesh.nodes[i, 0])
            y = _num(mesh.nodes[i, 1])
            f.write(f"{i + 1} {x} {y}\n")
    return path


def write_dep(mesh: Fort14Mesh, path: str | Path) -> Path:
    """Write ``casename_dep.dat`` (rows ``X Y H``, node order)."""
    _check_exportable(mesh)
    path = Path(path).resolve()
    with path.open("w") as f:
        f.write(f"Node Number = {mesh.n_nodes}\n")
        for i in range(mesh.n_nodes):
            x = _num(mesh.nodes[i, 0])
            y = _num(mesh.nodes[i, 1])
            f.write(f"{x} {y} {_num(mesh.depths[i])}\n")
    return path


def write_obc(
    mesh: Fort14Mesh,
    path: str | Path,
    *,
    obc_type: int | Sequence[int] = 1,
) -> Path:
    """Write ``casename_obc.dat`` from the mesh's open-boundary segments.

    ``obc_type`` is either one type for every node or one type per
    open segment (FVCOM types 1-10; odd = elevation only, even = plus
    nonlinear flux — see ``mod_obcs.F``). Zero open boundaries writes
    a valid 0-node file.
    """
    n_segs = len(mesh.open_boundaries)

    def _type(t):
        # a whole number, not a bool: int(1.9) would quietly be type 1
        # (review round 8 F6)
        if isinstance(t, (bool, np.bool_)) or not isinstance(t, (int, float, np.integer,
                                                                  np.floating)):
            raise ValueError(f"OBC type must be a whole number, not {t!r}")
        if not (np.isfinite(t) and float(t) == int(t)):
            raise ValueError(f"OBC type must be a whole number, not {t!r}")
        return int(t)

    if isinstance(obc_type, (int, float, np.integer, np.floating, bool, np.bool_)):
        seg_types = [_type(obc_type)] * n_segs
    else:
        seg_types = [_type(t) for t in obc_type]
        if len(seg_types) != n_segs:
            raise ValueError(
                f"obc_type has {len(seg_types)} entries for {n_segs} open segments"
            )
    for t in seg_types:
        if not 1 <= t <= 10:
            raise ValueError(f"OBC type {t} outside FVCOM's valid range 1-10")

    _check_exportable(mesh)
    path = Path(path).resolve()
    rows: list[tuple[int, int]] = []
    for seg, t in zip(mesh.open_boundaries, seg_types):
        for node in np.asarray(seg, dtype=np.int64):
            rows.append((int(node) + 1, t))
    with path.open("w") as f:
        f.write(f"OBC Node Number = {len(rows)}\n")
        for i, (gid, t) in enumerate(rows, start=1):
            f.write(f"{i} {gid} {t}\n")
    return path


def write_cor(
    mesh: Fort14Mesh, path: str | Path, cor: np.ndarray | Sequence[float],
) -> Path:
    """Write ``casename_cor.dat``; ``cor`` is the per-node Coriolis
    column (latitude in degrees for CARTESIAN builds).
    """
    _check_exportable(mesh)
    cor = _check_cor(mesh, cor)
    path = Path(path).resolve()
    with path.open("w") as f:
        f.write(f"Node Number = {mesh.n_nodes}\n")
        for i in range(mesh.n_nodes):
            x = _num(mesh.nodes[i, 0])
            y = _num(mesh.nodes[i, 1])
            f.write(f"{x} {y} {cor[i]:.6f}\n")
    return path


def write_spg(
    mesh: Fort14Mesh,
    path: str | Path,
    sponge: Sequence[tuple[int, float, float]] | None = None,
) -> Path:
    """Write ``casename_spg.dat``. ``sponge`` rows are 0-indexed
    ``(node_id, radius_m, damping)``; ``None`` writes the valid
    "no sponge" file (count 0).
    """
    rows = _check_sponge(mesh, sponge)
    path = Path(path).resolve()
    with path.open("w") as f:
        f.write(f"Sponge Node Number = {len(rows)}\n")
        for node, radius, damping in rows:
            # round-trip exact: a small radius must not be written as 0, which
            # FVCOM divides by (mod_setup.F; review round 8 F11)
            f.write(f"{int(node) + 1} {_num(radius)} {_num(damping)}\n")
    return path


def write_2dm(
    mesh: Fort14Mesh,
    path: str | Path,
    *,
    z_convention: str = "depth",
) -> Path:
    """Write an SMS ``.2dm`` interoperability export.

    ``z_convention``: ``"depth"`` stores the fort.14 positive-down
    depth as-is in the ND z column (the convention of this project's
    SMS-for-FVCOM workflows); ``"elevation"`` stores ``-depth``.
    One ``NS`` nodestring is emitted per open-boundary segment
    (1-indexed ids, last id negated, 10 ids per line).
    """
    if z_convention not in ("depth", "elevation"):
        raise ValueError(f"z_convention must be 'depth' or 'elevation', got {z_convention!r}")
    _validate_for_export(mesh)
    z = mesh.depths if z_convention == "depth" else -mesh.depths
    path = Path(path).resolve()
    with path.open("w") as f:
        f.write("MESH2D\n")
        for i in range(mesh.n_elements):
            n0, n1, n2 = (int(v) + 1 for v in mesh.elements[i])
            f.write(f"E3T {i + 1} {n0} {n1} {n2} 1\n")
        for i in range(mesh.n_nodes):
            x = _num(mesh.nodes[i, 0])
            y = _num(mesh.nodes[i, 1])
            f.write(f"ND {i + 1} {x} {y} {_num(z[i])}\n")
        for seg in mesh.open_boundaries:
            ids = [int(v) + 1 for v in np.asarray(seg, dtype=np.int64)]
            if not ids:
                continue
            ids[-1] = -ids[-1]
            for k in range(0, len(ids), 10):
                chunk = " ".join(str(v) for v in ids[k : k + 10])
                f.write(f"NS {chunk}\n")
    return path


def fvcom_next_obc(
    nodes: np.ndarray, elements: np.ndarray, obc_nodes: Sequence[int],
) -> tuple[np.ndarray, np.ndarray]:
    """Interior neighbour FVCOM pairs with each open-boundary node.

    Mirrors ``mod_obcs.F`` (Cartesian build): the inward normal of an OBC
    node is the (normalised sum of the) inward normals of the edges to its
    OBC neighbours, each oriented toward the centroid of the element(s)
    sharing that edge; ``NEXT_OBC`` is the edge-connected non-OBC node whose
    unit direction has the largest dot product with that normal.

    Returns ``(next_obc, margin)``: 0-indexed neighbour per OBC node and the
    gap between the best and second-best dot product. FVCOM breaks exact
    ties by its neighbour ordering, so a margin near zero means the choice
    is not reproducible from geometry alone.
    """
    nodes = np.asarray(nodes, float)
    tri = np.asarray(elements, int)
    obc = np.asarray(obc_nodes, int)
    is_obc = np.zeros(len(nodes), bool)
    is_obc[obc] = True
    nbrs: list[set[int]] = [set() for _ in range(len(nodes))]
    node_elems: list[list[int]] = [[] for _ in range(len(nodes))]
    for e, (a, b, c) in enumerate(tri):
        nbrs[a].update((b, c))
        nbrs[b].update((a, c))
        nbrs[c].update((a, b))
        for v in (a, b, c):
            node_elems[v].append(e)
    centroid = nodes[tri].mean(axis=1)

    next_obc = np.empty(len(obc), int)
    margin = np.empty(len(obc))
    for i, n in enumerate(obc):
        normal = np.zeros(2)
        for m in sorted(j for j in nbrs[n] if is_obc[j]):
            d = nodes[m] - nodes[n]
            unit = np.array([d[1], -d[0]]) / np.hypot(*d)
            for e in node_elems[n]:
                if m in tri[e]:
                    c = centroid[e] - nodes[n]
                    # FVCOM: CROSS = SIGN(1, DXC*DYN - DYC*DXN)
                    normal += np.sign(c[0] * d[1] - c[1] * d[0]) * unit
        if not normal.any():
            raise ValueError(f"OBC node {n} has no OBC neighbour; cannot define its normal")
        normal /= np.hypot(*normal)
        cand = sorted(j for j in nbrs[n] if not is_obc[j])
        if not cand:
            raise ValueError(f"OBC node {n} has no interior neighbour")
        vec = nodes[cand] - nodes[n]
        dots = (vec @ normal) / np.hypot(vec[:, 0], vec[:, 1])
        order = np.argsort(-dots)
        next_obc[i] = cand[order[0]]
        margin[i] = dots[order[0]] - dots[order[1]] if len(cand) > 1 else np.inf
    return next_obc, margin


def apply_obc_depth_control(mesh: Fort14Mesh) -> tuple[Fort14Mesh, np.ndarray]:
    """Copy of ``mesh`` with each OBC node's depth set to its NEXT_OBC depth.

    FVCOM does this at start-up (``OBC_DEPTH_CONTROL_ON``, default true,
    ``mod_startup.F``), so a mesh that does not already satisfy it runs
    with different bathymetry than the one written. Returns the new mesh
    and the per-OBC-node depth change (new - old), in open-boundary order.
    """
    obc = np.concatenate([np.asarray(b, int) for b in mesh.open_boundaries])
    nxt, _ = fvcom_next_obc(mesh.nodes, mesh.elements, obc)
    depths = mesh.depths.copy()
    change = depths[nxt] - depths[obc]
    depths[obc] = depths[nxt]
    return replace(mesh, depths=depths), change


def _indices(a, n: int, what: str, ndim: int = 1) -> np.ndarray:
    """``a`` as int64 node indices: finite whole numbers in ``[0, n)``.

    Checked before any cast, so neither -1 (NumPy's last node) nor 0.9
    (truncated to 0) reaches a file (review of the extend tools, rounds 6 F9
    and 7 F4).
    """
    arr = np.asarray(a)
    if arr.ndim != ndim:
        raise ValueError(f"{what} must be {ndim}-dimensional, not shape {arr.shape}")
    if arr.size == 0:
        return arr.astype(np.int64)
    if arr.dtype.kind not in "iuf" or arr.dtype.kind == "b":
        raise ValueError(f"{what} must be numeric node indices, not {arr.dtype}")
    if arr.dtype.kind == "f" and not (np.isfinite(arr).all() and (arr == np.round(arr)).all()):
        raise ValueError(f"{what} holds an index that is not a whole number")
    out = arr.astype(np.int64)
    if out.min() < 0 or out.max() >= n:
        raise ValueError(f"{what}: node index outside [0, {n}): {int(out.min())}..{int(out.max())}")
    return out


def _check_exportable(mesh: Fort14Mesh) -> None:
    """Refuse, before any file is written, a mesh FVCOM cannot read: element
    or boundary indices that are not whole numbers in ``[0, NP)``, or
    coordinates and depths that are not finite (review of the extend tools,
    rounds 6 F9 and 7 F4)."""
    nodes = np.asarray(mesh.nodes)
    depths = np.asarray(mesh.depths)
    n = len(nodes)
    if nodes.ndim != 2 or nodes.shape[1] < 2 or depths.shape != (n,):
        raise ValueError(f"nodes {nodes.shape} and depths {depths.shape} do not match")
    if np.asarray(mesh.elements).dtype.kind not in "iu":
        # the elements index arrays later on; whole floats would be accepted
        # here and fail there (review round 8 F4)
        raise ValueError(f"elements must be integers, not {np.asarray(mesh.elements).dtype}")
    els = _indices(mesh.elements, n, "elements", ndim=2)
    if els.size and els.shape[1] != 3:
        raise ValueError(f"elements must be (NE, 3), not {els.shape}")
    for seg in mesh.open_boundaries:
        _indices(seg, n, "an open boundary")
    for _t, seg in mesh.land_boundaries:
        _indices(seg, n, "a land boundary")
    if not np.isfinite(nodes[:, :2]).all():
        raise ValueError(f"{int((~np.isfinite(nodes[:, :2])).any(axis=1).sum())} node(s) "
                         "with a coordinate that is not finite")
    if not np.isfinite(depths).all():
        raise ValueError(f"{int((~np.isfinite(depths)).sum())} depth(s) not finite")


def _check_cor(mesh: Fort14Mesh, cor) -> np.ndarray:
    cor = np.asarray(cor, dtype=np.float64)
    if cor.shape != (mesh.n_nodes,):
        raise ValueError(f"cor shape {cor.shape} does not match n_nodes = {mesh.n_nodes}")
    if not np.isfinite(cor).all():
        raise ValueError(f"{int((~np.isfinite(cor)).sum())} Coriolis value(s) not finite")
    return cor


def _check_sponge(mesh: Fort14Mesh, sponge) -> list[tuple[int, float, float]]:
    # ``is None``, not truthiness: a NumPy row array has no truth value
    # (review round 9 F8)
    rows = [] if sponge is None else list(sponge)
    if not rows:
        return []
    arr = np.asarray(rows, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 3:
        raise ValueError(f"sponge rows are (node, radius, damping), not shape {arr.shape}")
    nodes = _indices(arr[:, 0], mesh.n_nodes, "sponge nodes")
    radius, damping = arr[:, 1], arr[:, 2]
    if not (np.isfinite(radius).all() and (radius > 0).all()):
        raise ValueError("sponge radii must be finite and positive")
    if not (np.isfinite(damping).all() and (damping >= 0).all()):
        raise ValueError("sponge damping must be finite and not negative")
    return [(int(k), float(r), float(c)) for k, r, c in zip(nodes, radius, damping)]


def export_fvcom_case(
    mesh: Fort14Mesh,
    outdir: str | Path,
    casename: str,
    *,
    obc_type: int | Sequence[int] | None = None,
    cor: np.ndarray | Sequence[float] | None = None,
    sponge: Sequence[tuple[int, float, float]] | None = None,
    write_empty_spg: bool = False,
    twodm: bool = True,
    z_convention: str = "depth",
    obc_depth_control: bool = True,
) -> dict[str, Path]:
    """Write the full FVCOM input set for ``casename`` into ``outdir``.

    With ``obc_depth_control`` (default) the open-boundary depths are first
    set to their NEXT_OBC depths (:func:`apply_obc_depth_control`), so the
    written ``_dep.dat`` is exactly the bathymetry FVCOM will run with.

    Always writes ``_grd.dat``, ``_dep.dat``, ``_obc.dat``; ``_cor.dat``
    when ``cor`` is given; ``_spg.dat`` when ``sponge`` is given or
    ``write_empty_spg`` is set; ``.2dm`` unless ``twodm=False``.
    Returns the mapping of file kind to written path.
    """
    # A mesh read by read_fvcom_case carries the type its file declared, and
    # that is the right default: making every caller remember to pass it is
    # how a base declaring type 2 got written back out as type 1 with nothing
    # said (fourth review).
    if obc_type is None:
        obc_type = getattr(mesh, "obc_type", 1)
    # every input is checked before any file is opened, and the set is
    # written to a temporary directory and moved into place only when all of
    # it is written: a bad optional input no longer leaves a half-replaced
    # case (review of the extend tools, round 7 F5)
    _validate_for_export(mesh)
    if cor is not None:
        _check_cor(mesh, cor)
    # normalised once and passed on: an iterator would be spent by the check
    # and write an empty sponge (review round 8 F5)
    sponge_rows = _check_sponge(mesh, sponge) if sponge is not None else None
    if obc_depth_control and mesh.open_boundaries:
        mesh, _ = apply_obc_depth_control(mesh)
    outdir = Path(outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    import shutil
    import tempfile

    stage = Path(tempfile.mkdtemp(dir=outdir, prefix=f".{casename}.export."))
    # The stage directory is removed only after a known outcome: nothing
    # replaced yet, everything published, or everything put back. Any other
    # exit -- a failed restore, or an interrupt in the middle of one -- keeps
    # it, with the previous files in .previous (review rounds 8 F2, 9 F3).
    state = "staging"
    try:
        staged: dict[str, Path] = {
            "grd": write_grd(mesh, stage / f"{casename}_grd.dat"),
            "dep": write_dep(mesh, stage / f"{casename}_dep.dat"),
            "obc": write_obc(mesh, stage / f"{casename}_obc.dat", obc_type=obc_type),
        }
        if cor is not None:
            staged["cor"] = write_cor(mesh, stage / f"{casename}_cor.dat", cor)
        if sponge is not None or write_empty_spg:
            staged["spg"] = write_spg(mesh, stage / f"{casename}_spg.dat", sponge_rows)
        if twodm:
            staged["2dm"] = write_2dm(
                mesh, stage / f"{casename}.2dm", z_convention=z_convention,
            )
        keep = stage / ".previous"
        keep.mkdir()
        written: dict[str, Path] = {}
        moved: list[Path] = []
        state = "publishing"
        try:
            for kind, p in staged.items():
                dst = outdir / p.name
                if dst.exists():
                    shutil.copy2(dst, keep / p.name)
                # recorded before the move, so an interrupt right after it is
                # still rolled back (review round 10 F6)
                moved.append(dst)
                p.replace(dst)
                written[kind] = dst
        except BaseException:
            state = "restoring"
            unrestored = []
            for dst in moved:
                try:
                    if (keep / dst.name).exists():
                        (keep / dst.name).replace(dst)
                    else:
                        dst.unlink(missing_ok=True)
                except OSError:
                    unrestored.append(dst.name)
            if unrestored:
                raise OSError(f"export into {outdir} failed and {unrestored} could not be "
                              f"restored; the previous files are in {keep}") from None
            state = "restored"
            raise
        state = "done"
    finally:
        if state in ("staging", "done", "restored"):
            shutil.rmtree(stage, ignore_errors=True)
    return written


__all__ = [
    "apply_obc_depth_control",
    "export_fvcom_case",
    "fvcom_next_obc",
    "write_2dm",
    "write_cor",
    "write_dep",
    "write_grd",
    "write_obc",
    "write_spg",
]


# ---------------------------------------------------------------------------
# Readers
# ---------------------------------------------------------------------------
#
# Added because local refinement takes a FINISHED FVCOM case as its base: the
# mesh and the depth file are inputs, not something to rebuild (owner,
# 2026-09-22).  The production case is `_grd.dat` + `_dep.dat` + `_obc.dat`,
# so those have to be readable, and the two files disagree on purpose -- the
# grd carries a depth column from whenever it was made, and the dep file is
# the one the baseline's bathymetry tag names.  On goto2023 node 1 that is
# 4.312072 m in the grd and 7.161207 m in
# `TokyoBay_dep_m7001tp_rfac0p2_cap300.dat`.  The dep file wins; the grd's
# column is read only to be ignored.


def _header_count(line: str, label: str) -> int:
    if "=" not in line:
        raise ValueError(f"expected a '{label} = N' header, got {line!r}")
    key, value = line.split("=", 1)
    if key.strip().lower() != label.lower():
        raise ValueError(f"expected header {label!r}, got {key.strip()!r}")
    return int(value.strip().split()[0])


def read_grd(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """Read ``casename_grd.dat``; returns ``(nodes (N,2), elements (M,3))``.

    Element rows are ``CELL# N1 N2 N3`` with any trailing columns ignored
    (OceanMesh2D writes the cell number again), node rows ``NODE# X Y`` with
    an optional depth column that is NOT returned.  Node ids are 1-indexed in
    the file and 0-indexed here.
    """
    path = Path(path).resolve()
    with path.open() as f:
        lines = [ln for ln in (x.strip() for x in f) if ln]
    n_nodes = _header_count(lines[0], "Node Number")
    n_cells = _header_count(lines[1], "Cell Number")
    # exactly the declared records, each with its fields: FVCOM reads to the
    # end and refuses surplus nodes, and a short row is no coordinate pair
    # (review of the extend tools, round 20 F7)
    if len(lines) != 2 + n_cells + n_nodes:
        raise ValueError(
            f"{path.name}: {len(lines) - 2} rows for {n_cells} cells + "
            f"{n_nodes} nodes")
    if any(len(ln.split()) < 4 for ln in lines[2:2 + n_cells]):
        raise ValueError(f"{path.name}: a cell row has fewer than 4 fields")
    if any(len(ln.split()) < 3 for ln in lines[2 + n_cells:]):
        raise ValueError(f"{path.name}: a node row has fewer than 3 fields")
    # record ids 1..NE and 1..NP in order, as FVCOM reads them: it finds the
    # node block by its first id (review of the extend tools, round 19 F6)
    cell_ids = [int(ln.split()[0]) for ln in lines[2:2 + n_cells]]
    if cell_ids != list(range(1, n_cells + 1)):
        raise ValueError(f"{path.name}: cell ids must run 1..{n_cells} in order")
    node_ids = [int(ln.split()[0]) for ln in lines[2 + n_cells:2 + n_cells + n_nodes]]
    if node_ids != list(range(1, n_nodes + 1)):
        raise ValueError(f"{path.name}: node ids must run 1..{n_nodes} in order")
    elements = np.array(
        [[int(w) for w in ln.split()[1:4]] for ln in lines[2:2 + n_cells]],
        dtype=np.int64) - 1
    rows = lines[2 + n_cells:2 + n_cells + n_nodes]
    nodes = np.array([[float(w) for w in ln.split()[1:3]] for ln in rows],
                     dtype=float)
    # A NaN coordinate reaches FVCOM as a mesh it cannot build and reaches
    # every geometric test here as a comparison that is quietly False; the
    # dep file was checked for this and the grd file was not (fourth review).
    if not np.isfinite(nodes).all():
        raise ValueError(
            f"{path.name}: {int((~np.isfinite(nodes)).any(axis=1).sum())} "
            "nodes have a non-finite coordinate")
    if elements.min() < 0 or elements.max() >= n_nodes:
        raise ValueError(f"{path.name}: connectivity references a node outside "
                         f"1..{n_nodes}")
    return nodes, elements


def read_dep(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """Read ``casename_dep.dat``; returns ``(xy (N,2), depths (N,))``."""
    path = Path(path).resolve()
    with path.open() as f:
        lines = [ln for ln in (x.strip() for x in f) if ln]
    n_nodes = _header_count(lines[0], "Node Number")
    data = np.array([[float(w) for w in ln.split()[:3]]
                     for ln in lines[1:1 + n_nodes]], dtype=float)
    if data.shape != (n_nodes, 3):
        raise ValueError(f"{path.name}: expected {n_nodes} rows of 'X Y H'")
    return data[:, :2], data[:, 2]


def read_obc(path: str | Path, with_types: bool = False):
    """Read ``casename_obc.dat``; returns the 0-indexed node ids in file order.

    File order is along-boundary order for every file this project writes or
    consumes, and the refinement driver relies on it: it splits the outer
    boundary loop at the arc's two ends.

    ``with_types`` also returns the FVCOM boundary type of each row.  The type
    is part of the model input -- odd is elevation only, even adds nonlinear
    flux (``mod_obcs.F``) -- and dropping it on read means a case written back
    out gets the writer's default of 1.  A base declaring type 2 came back
    as type 1 with nothing said (fourth review).
    """
    path = Path(path).resolve()
    with path.open() as f:
        lines = [ln for ln in (x.strip() for x in f) if ln]
    n = _header_count(lines[0], "OBC Node Number")
    rows = [ln.split() for ln in lines[1:1 + n]]
    if len(rows) != n:
        raise ValueError(f"{path.name}: expected {n} OBC rows")
    ids = np.array([int(r[1]) for r in rows], dtype=np.int64) - 1
    if not with_types:
        return ids
    types = np.array([int(r[2]) if len(r) > 2 else 1 for r in rows],
                     dtype=np.int64)
    return ids, types


def read_obc_types(path: str | Path) -> list[int]:
    """The distinct FVCOM boundary types in an ``_obc.dat``, in file order."""
    _, types = read_obc(path, with_types=True)
    seen: list[int] = []
    for t in types.tolist():
        if t not in seen:
            seen.append(int(t))
    return seen


def boundary_loops(elements: np.ndarray) -> list[np.ndarray]:
    """Ordered closed walks of the mesh boundary, one per loop."""
    tri = np.asarray(elements, dtype=np.int64)
    e = np.sort(np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]), axis=1)
    u, c = np.unique(e, axis=0, return_counts=True)
    b = u[c == 1]
    nbr: dict[int, list[int]] = {}
    for x, y in b.tolist():
        nbr.setdefault(x, []).append(y)
        nbr.setdefault(y, []).append(x)
    bad = [v for v, w in nbr.items() if len(w) != 2]
    if bad:
        raise ValueError(f"{len(bad)} boundary nodes are not on exactly two "
                         "boundary edges; the mesh boundary is not a set of loops")
    seen: set[int] = set()
    loops: list[np.ndarray] = []
    for start in nbr:
        if start in seen:
            continue
        walk, prev, cur = [start], None, start
        seen.add(start)
        while True:
            a, b2 = nbr[cur]
            nxt = a if a != prev else b2
            if nxt == start:
                break
            walk.append(nxt)
            seen.add(nxt)
            prev, cur = cur, nxt
        loops.append(np.asarray(walk, dtype=np.int64))
    return loops


def read_fvcom_case(
    grd: str | Path,
    dep: str | Path,
    obc: str | Path | None = None,
    *,
    title: str | None = None,
    coord_tol_m: float = 1e-3,
) -> Fort14Mesh:
    """Read a finished FVCOM case into a :class:`Fort14Mesh`.

    The depths come from ``dep``, never from the grd's own column. The node
    coordinates in the two files are checked against each other: a dep file
    built for a different mesh is otherwise silently accepted, and the
    refinement contract is about keeping those depths on those nodes.

    Land boundaries are derived, because FVCOM has no land-boundary file: the
    outer loop minus the open-boundary run is ibtype 20, every other loop is
    an island (ibtype 21), matching what :func:`read_fort14` would give.
    """
    if not np.isfinite(coord_tol_m) or coord_tol_m < 0:
        raise ValueError("coord_tol_m must be finite and non-negative; a NaN "
                         "tolerance passes every comparison")
    nodes, elements = read_grd(grd)
    dep_xy, depths = read_dep(dep)
    if dep_xy.shape[0] != nodes.shape[0]:
        raise ValueError(
            f"{Path(dep).name} has {dep_xy.shape[0]} nodes, "
            f"{Path(grd).name} has {nodes.shape[0]}")
    if not np.isfinite(dep_xy).all() or not np.isfinite(depths).all():
        raise ValueError(f"{Path(dep).name} contains non-finite values")
    off = np.linalg.norm(dep_xy - nodes, axis=1)
    # `off.max() > tol` is False for NaN, so a single NaN coordinate used to
    # turn this check off entirely rather than fail it (fourth review); the
    # finiteness test above is what catches that.
    if float(off.max()) > coord_tol_m:
        raise ValueError(
            f"{Path(dep).name} does not sit on {Path(grd).name}: worst node "
            f"offset {off.max():.3g} m")

    open_boundaries: list[np.ndarray] = []
    obc_types: list[int] = []
    if obc is not None:
        ids, types = read_obc(obc, with_types=True)
        if ids.size:
            open_boundaries = [ids]
            obc_types = sorted(set(types.tolist()))

    loops = boundary_loops(elements)
    area = [abs(float(np.dot(nodes[lp, 0], np.roll(nodes[lp, 1], -1))
                      - np.dot(nodes[lp, 1], np.roll(nodes[lp, 0], -1))) / 2)
            for lp in loops]
    outer = loops[int(np.argmax(area))]
    land: list[tuple[int, np.ndarray]] = []
    if open_boundaries:
        arc = set(open_boundaries[0].tolist())
        if not arc.issubset(set(outer.tolist())):
            raise ValueError("the open boundary is not on the outer loop")
        ring = np.roll(outer, -int(np.where(outer == open_boundaries[0][0])[0][0]))
        # the direction is the one the list takes from its first node, not a
        # guess from membership (both neighbours may be open: review round 2 F13)
        if len(open_boundaries[0]) > 1 and ring[1] != open_boundaries[0][1]:
            ring = np.roll(ring[::-1], 1)
        stop = int(np.where(ring == open_boundaries[0][-1])[0][0])
        # the list must BE the walk along the loop, in order and without
        # repeats, not only the same set of nodes (review of the extend tools,
        # round 1, F29: [0, 2, 1, 5, 8] passed on a 3x3 grid)
        if (len(open_boundaries[0]) != len(arc)
                or not np.array_equal(ring[:stop + 1], open_boundaries[0])):
            raise ValueError("the open-boundary nodes are not a consecutive walk along "
                             "the outer loop (repeated, out of order, or with gaps)")
        land.append((20, np.append(ring[stop:], ring[0])))
    else:
        land.append((20, np.append(outer, outer[0])))
    land += [(21, lp) for lp in loops if lp is not outer]

    if len(obc_types) > 1:
        raise ValueError(
            f"{Path(obc).name} mixes OBC types {obc_types}; this reader keeps "
            "one type per segment, and silently collapsing them would change "
            "the model input")
    return Fort14Mesh(
        title=title or Path(grd).stem,
        nodes=nodes, depths=depths, elements=elements,
        open_boundaries=open_boundaries, land_boundaries=land,
        obc_type=obc_types[0] if obc_types else 1)
