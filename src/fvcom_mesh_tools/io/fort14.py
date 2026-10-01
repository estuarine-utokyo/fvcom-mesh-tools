"""Read ADCIRC/FVCOM ``fort.14`` unstructured-mesh files."""

from __future__ import annotations

from dataclasses import dataclass
from io import TextIOBase
from pathlib import Path

import numpy as np


@dataclass
class Fort14Mesh:
    """In-memory representation of a fort.14 mesh.

    Coordinate columns are passed through unchanged (lon/lat or projected
    metres, depending on the source file). Node and element IDs in the file
    are 1-indexed; arrays returned by :func:`read_fort14` are 0-indexed so
    they can be used directly to index into ``nodes`` / ``depths``.

    Attributes
    ----------
    title:
        First line of the file.
    nodes:
        ``(NP, 2)`` float array of ``(x, y)`` coordinates in file order.
    depths:
        ``(NP,)`` float array of the fourth-column "z" / depth values.
    elements:
        ``(NE, 3)`` int array of 0-indexed node indices for each triangular
        element.
    open_boundaries:
        One int array per open-boundary segment, holding 0-indexed node
        indices in along-boundary order.
    land_boundaries:
        One ``(ibtype, ids)`` tuple per land/normal-flow boundary segment,
        where ``ibtype`` is the integer boundary-type code (0 = normal
        coast in the ADCIRC convention) and ``ids`` is a 0-indexed int
        array of node indices.
    obc_type:
        FVCOM's open-boundary condition type for the open boundary (1 in the
        goto2023 production case, 3 in others). fort.14 has no such column;
        it is carried here because ``dataclasses.replace`` -- which every
        transformation in this package uses -- drops an attribute that is not
        a field, and a depth-control pass silently turned a type 3 boundary
        back into type 1 (fourth review).
    """

    title: str
    nodes: np.ndarray
    depths: np.ndarray
    elements: np.ndarray
    open_boundaries: list[np.ndarray]
    land_boundaries: list[tuple[int, np.ndarray]]
    obc_type: int = 1

    @property
    def n_nodes(self) -> int:
        return int(self.nodes.shape[0])

    @property
    def n_elements(self) -> int:
        return int(self.elements.shape[0])

    @property
    def bbox(self) -> tuple[float, float, float, float]:
        x = self.nodes[:, 0]
        y = self.nodes[:, 1]
        return float(x.min()), float(y.min()), float(x.max()), float(y.max())


def _read_first_int(f: TextIOBase) -> int:
    line = f.readline()
    if not line:
        raise ValueError("unexpected EOF while reading an integer")
    return int(line.split()[0])


def _read_two_ints(f: TextIOBase) -> tuple[int, int]:
    line = f.readline()
    if not line:
        raise ValueError("unexpected EOF while reading a (count, ibtype) pair")
    parts = line.split()
    return int(parts[0]), int(parts[1])


def _read_node_ids(f: TextIOBase, n: int) -> np.ndarray:
    ids = np.empty(n, dtype=np.int64)
    for j in range(n):
        line = f.readline()
        if not line:
            raise ValueError(f"unexpected EOF while reading boundary node {j + 1}/{n}")
        ids[j] = int(line.split()[0]) - 1
    return ids


def read_fort14(path: str | Path) -> Fort14Mesh:
    """Parse a fort.14 file into a :class:`Fort14Mesh`.

    The standard layout is assumed: header line, ``NE NP`` line, ``NP``
    node rows ``id x y depth``, ``NE`` triangular-element rows
    ``id 3 n1 n2 n3``, then the boundary block (open boundaries followed
    by land/normal-flow boundaries). Node indices in the returned arrays
    are 0-indexed.
    """
    path = Path(path).resolve()
    with path.open("r") as f:
        title = f.readline().rstrip("\n")
        ne, np_ = _read_two_ints(f)

        # ndmin=2: a single row would otherwise come back one-dimensional
        # (review of the extend tools, round 9 F11)
        node_block = np.loadtxt(f, max_rows=np_, dtype=np.float64, ndmin=2)
        if node_block.shape != (np_, 4):
            raise ValueError(
                f"node block shape {node_block.shape} does not match expected ({np_}, 4); "
                f"check whether the header is 'NE NP' (ADCIRC convention)"
            )
        # The record ids are what the connectivity and boundaries refer to:
        # they must be 1..NP and 1..NE in order, and every element a
        # triangle (type 3) -- otherwise a reference would be read as another
        # node (round 9 F12).
        if not np.array_equal(node_block[:, 0], np.arange(1, np_ + 1)):
            raise ValueError("node ids must run 1..NP in order")
        nodes = node_block[:, 1:3].copy()
        depths = node_block[:, 3].copy()

        elem_block = np.loadtxt(f, max_rows=ne, dtype=np.int64, ndmin=2)
        if elem_block.shape != (ne, 5):
            raise ValueError(
                f"element block shape {elem_block.shape} does not match expected ({ne}, 5)"
            )
        if not np.array_equal(elem_block[:, 0], np.arange(1, ne + 1)):
            raise ValueError("element ids must run 1..NE in order")
        if ne and not (elem_block[:, 1] == 3).all():
            raise ValueError("every element must be a triangle (type 3)")
        elements = (elem_block[:, 2:5] - 1).copy()
        if ne and (elements.min() < 0 or elements.max() >= np_):
            raise ValueError(f"element node references outside 1..{np_}")

        # Boundary counts are not negative, their totals agree with NETA and
        # NVEL, and every id is a node 1..NP (review of the extend tools,
        # round 10 F9).
        nope = _read_first_int(f)
        neta = _read_first_int(f)
        if nope < 0 or neta < 0:
            raise ValueError(f"negative open-boundary counts: NOPE {nope}, NETA {neta}")
        open_boundaries: list[np.ndarray] = []
        for _ in range(nope):
            n = _read_first_int(f)
            if n < 0:
                raise ValueError(f"negative open-boundary length {n}")
            open_boundaries.append(_read_node_ids(f, n))
        if sum(len(b) for b in open_boundaries) != neta:
            raise ValueError(f"open boundaries hold {sum(len(b) for b in open_boundaries)} "
                             f"nodes, NETA says {neta}")

        nbou = _read_first_int(f)
        nvel = _read_first_int(f)
        if nbou < 0 or nvel < 0:
            raise ValueError(f"negative land-boundary counts: NBOU {nbou}, NVEL {nvel}")
        land_boundaries: list[tuple[int, np.ndarray]] = []
        for _ in range(nbou):
            n, ibtype = _read_two_ints(f)
            if n < 0:
                raise ValueError(f"negative land-boundary length {n}")
            land_boundaries.append((ibtype, _read_node_ids(f, n)))
        if sum(len(b) for _t, b in land_boundaries) != nvel:
            raise ValueError(f"land boundaries hold "
                             f"{sum(len(b) for _t, b in land_boundaries)} nodes, NVEL says {nvel}")
        for b in [*open_boundaries, *(b for _t, b in land_boundaries)]:
            if b.size and (b.min() < 0 or b.max() >= np_):
                raise ValueError(f"boundary node ids outside 1..{np_}")

    return Fort14Mesh(
        title=title,
        nodes=nodes,
        depths=depths,
        elements=elements,
        open_boundaries=open_boundaries,
        land_boundaries=land_boundaries,
    )


def write_fort14(mesh: Fort14Mesh, path: str | Path) -> None:
    """Write a :class:`Fort14Mesh` to ``path`` in standard ADCIRC fort.14 layout.

    The output is round-trip safe: ``read_fort14(write_fort14(m, p))`` recovers
    the same node coordinates, depths, element connectivity, and boundary
    structure as ``m``. The exact numeric formatting of the source file is not
    preserved; coordinates are written as the shortest text that reads back
    as the same double (``repr``), so every coordinate round-trips exactly
    and triangles with very small but positive signed area survive without
    being pancaked to zero.

    Element and boundary indices are checked before the file is opened:
    whole numbers in ``[0, NP)``, never truncated (review of the extend
    tools, round 10 F9).

    Depths carry ``.17g``, not the ``.10e`` this used to write. Eleven
    significant figures is enough for a depth that was itself read from a
    fort.14 and is not enough for one that was computed: 5.12345678912345 m
    came back 2.3e-11 m different, which is physically nothing and is still a
    round trip this docstring promises and did not deliver.
    """
    from fvcom_mesh_tools.io.fvcom_native import _indices

    path = Path(path).resolve()
    n_nodes = mesh.n_nodes
    n_elements = mesh.n_elements
    # shapes too, before the destination is touched (review round 11 F5)
    nodes_a, depths_a = np.asarray(mesh.nodes), np.asarray(mesh.depths)
    if nodes_a.ndim != 2 or nodes_a.shape[1] != 2 or depths_a.shape != (n_nodes,):
        raise ValueError(f"nodes must be (NP, 2) and depths (NP,), not {nodes_a.shape} "
                         f"and {depths_a.shape}")
    if mesh.n_elements:
        if np.asarray(mesh.elements).dtype.kind not in "iu":
            raise ValueError("elements must be integers")
        if np.asarray(mesh.elements).shape[1:] != (3,):
            raise ValueError(f"elements must be (NE, 3), not {np.asarray(mesh.elements).shape}")
        _indices(mesh.elements, n_nodes, "elements", ndim=2)
    for b in mesh.open_boundaries:
        _indices(b, n_nodes, "an open boundary")
    for _t, b in mesh.land_boundaries:
        _indices(b, n_nodes, "a land boundary")
    n_open_segs = len(mesh.open_boundaries)
    n_open_nodes = sum(len(b) for b in mesh.open_boundaries)
    n_land_segs = len(mesh.land_boundaries)
    n_land_nodes = sum(len(ids) for _, ids in mesh.land_boundaries)

    # written beside the destination and moved over it: a failure part-way
    # leaves an existing file as it was (round 11 F5)
    import os
    import tempfile

    fd, tmp_name = tempfile.mkstemp(dir=path.parent, prefix=path.name + ".", suffix=".tmp")
    tmp = Path(tmp_name)
    try:
        with os.fdopen(fd, "w") as f:
            _write_fort14_body(f, mesh, n_nodes, n_elements, n_open_segs, n_open_nodes,
                               n_land_segs, n_land_nodes)
        # mkstemp makes the file 0600; give it the mode a plain open would
        umask = os.umask(0)
        os.umask(umask)
        os.chmod(tmp, 0o666 & ~umask)
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def _write_fort14_body(f, mesh, n_nodes, n_elements, n_open_segs, n_open_nodes,
                       n_land_segs, n_land_nodes) -> None:
    f.write(f"{mesh.title}\n")
    f.write(f"{n_elements} {n_nodes}\n")

    for i in range(n_nodes):
        x, y = mesh.nodes[i]
        # shortest round-trip-exact text, as the depths already were:
        # .15f cut 0.12345678912345678 (review of the extend tools,
        # round 3 F8)
        f.write(f"{i + 1:>10d}  {float(x)!r}  {float(y)!r}  {mesh.depths[i]:.17g}\n")

    for i in range(n_elements):
        n0, n1, n2 = mesh.elements[i]
        f.write(f"{i + 1:>10d}  3  {n0 + 1:>10d}  {n1 + 1:>10d}  {n2 + 1:>10d}\n")

    f.write(f"{n_open_segs} = Number of open boundaries\n")
    f.write(f"{n_open_nodes} = Total number of open boundary nodes\n")
    for k, ids in enumerate(mesh.open_boundaries, start=1):
        f.write(f"{len(ids)} = Number of nodes for open boundary {k}\n")
        for node in ids:
            f.write(f"{int(node) + 1}\n")

    f.write(f"{n_land_segs} = Number of normal flow boundaries\n")
    f.write(f"{n_land_nodes} = Total number of land boundary nodes\n")
    for k, (ibtype, ids) in enumerate(mesh.land_boundaries, start=1):
        f.write(f"{len(ids)} {ibtype} = Number of nodes for land boundary {k}\n")
        for node in ids:
            f.write(f"{int(node) + 1}\n")
