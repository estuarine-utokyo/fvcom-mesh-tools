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

from pathlib import Path
from typing import Sequence

import numpy as np

from fvcom_mesh_tools.io.fort14 import Fort14Mesh

_COORD_FMT = "{:.8f}"
_DEPTH_FMT = "{:.6f}"


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
            x = _COORD_FMT.format(mesh.nodes[i, 0])
            y = _COORD_FMT.format(mesh.nodes[i, 1])
            f.write(f"{i + 1} {x} {y}\n")
    return path


def write_dep(mesh: Fort14Mesh, path: str | Path) -> Path:
    """Write ``casename_dep.dat`` (rows ``X Y H``, node order)."""
    path = Path(path).resolve()
    with path.open("w") as f:
        f.write(f"Node Number = {mesh.n_nodes}\n")
        for i in range(mesh.n_nodes):
            x = _COORD_FMT.format(mesh.nodes[i, 0])
            y = _COORD_FMT.format(mesh.nodes[i, 1])
            f.write(f"{x} {y} {_DEPTH_FMT.format(mesh.depths[i])}\n")
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
    if isinstance(obc_type, int):
        seg_types = [obc_type] * n_segs
    else:
        seg_types = [int(t) for t in obc_type]
        if len(seg_types) != n_segs:
            raise ValueError(
                f"obc_type has {len(seg_types)} entries for {n_segs} open segments"
            )
    for t in seg_types:
        if not 1 <= t <= 10:
            raise ValueError(f"OBC type {t} outside FVCOM's valid range 1-10")

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
    cor = np.asarray(cor, dtype=np.float64)
    if cor.shape != (mesh.n_nodes,):
        raise ValueError(
            f"cor shape {cor.shape} does not match n_nodes = {mesh.n_nodes}"
        )
    path = Path(path).resolve()
    with path.open("w") as f:
        f.write(f"Node Number = {mesh.n_nodes}\n")
        for i in range(mesh.n_nodes):
            x = _COORD_FMT.format(mesh.nodes[i, 0])
            y = _COORD_FMT.format(mesh.nodes[i, 1])
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
    rows = list(sponge) if sponge else []
    for node, _r, _c in rows:
        if not 0 <= int(node) < mesh.n_nodes:
            raise ValueError(f"sponge node id {node} out of range")
    path = Path(path).resolve()
    with path.open("w") as f:
        f.write(f"Sponge Node Number = {len(rows)}\n")
        for node, radius, damping in rows:
            f.write(f"{int(node) + 1} {float(radius):.4f} {float(damping):.6f}\n")
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
            x = _COORD_FMT.format(mesh.nodes[i, 0])
            y = _COORD_FMT.format(mesh.nodes[i, 1])
            f.write(f"ND {i + 1} {x} {y} {_DEPTH_FMT.format(z[i])}\n")
        for seg in mesh.open_boundaries:
            ids = [int(v) + 1 for v in np.asarray(seg, dtype=np.int64)]
            if not ids:
                continue
            ids[-1] = -ids[-1]
            for k in range(0, len(ids), 10):
                chunk = " ".join(str(v) for v in ids[k : k + 10])
                f.write(f"NS {chunk}\n")
    return path


def export_fvcom_case(
    mesh: Fort14Mesh,
    outdir: str | Path,
    casename: str,
    *,
    obc_type: int | Sequence[int] = 1,
    cor: np.ndarray | Sequence[float] | None = None,
    sponge: Sequence[tuple[int, float, float]] | None = None,
    write_empty_spg: bool = False,
    twodm: bool = True,
    z_convention: str = "depth",
) -> dict[str, Path]:
    """Write the full FVCOM input set for ``casename`` into ``outdir``.

    Always writes ``_grd.dat``, ``_dep.dat``, ``_obc.dat``; ``_cor.dat``
    when ``cor`` is given; ``_spg.dat`` when ``sponge`` is given or
    ``write_empty_spg`` is set; ``.2dm`` unless ``twodm=False``.
    Returns the mapping of file kind to written path.
    """
    outdir = Path(outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {
        "grd": write_grd(mesh, outdir / f"{casename}_grd.dat"),
        "dep": write_dep(mesh, outdir / f"{casename}_dep.dat"),
        "obc": write_obc(mesh, outdir / f"{casename}_obc.dat", obc_type=obc_type),
    }
    if cor is not None:
        written["cor"] = write_cor(mesh, outdir / f"{casename}_cor.dat", cor)
    if sponge is not None or write_empty_spg:
        written["spg"] = write_spg(mesh, outdir / f"{casename}_spg.dat", sponge)
    if twodm:
        written["2dm"] = write_2dm(
            mesh, outdir / f"{casename}.2dm", z_convention=z_convention,
        )
    return written


__all__ = [
    "export_fvcom_case",
    "write_2dm",
    "write_cor",
    "write_dep",
    "write_grd",
    "write_obc",
    "write_spg",
]
