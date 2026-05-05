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
    """

    title: str
    nodes: np.ndarray
    depths: np.ndarray
    elements: np.ndarray
    open_boundaries: list[np.ndarray]
    land_boundaries: list[tuple[int, np.ndarray]]

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

        node_block = np.loadtxt(f, max_rows=np_, dtype=np.float64)
        if node_block.shape != (np_, 4):
            raise ValueError(
                f"node block shape {node_block.shape} does not match expected ({np_}, 4); "
                f"check whether the header is 'NE NP' (ADCIRC convention)"
            )
        nodes = node_block[:, 1:3].copy()
        depths = node_block[:, 3].copy()

        elem_block = np.loadtxt(f, max_rows=ne, dtype=np.int64)
        if elem_block.shape != (ne, 5):
            raise ValueError(
                f"element block shape {elem_block.shape} does not match expected ({ne}, 5)"
            )
        elements = (elem_block[:, 2:5] - 1).copy()

        nope = _read_first_int(f)
        _ = f.readline()  # NETA: redundant total of open boundary nodes
        open_boundaries: list[np.ndarray] = []
        for _ in range(nope):
            n = _read_first_int(f)
            open_boundaries.append(_read_node_ids(f, n))

        nbou = _read_first_int(f)
        _ = f.readline()  # NVEL: redundant total of land boundary nodes
        land_boundaries: list[tuple[int, np.ndarray]] = []
        for _ in range(nbou):
            n, ibtype = _read_two_ints(f)
            land_boundaries.append((ibtype, _read_node_ids(f, n)))

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
    preserved; coordinates are written with 15 decimal digits (which exceeds
    float64's ~15-17 significant figures at coordinates ~140 deg) so that
    every writable double round-trips exactly and triangles with very small
    but positive signed area survive without being pancaked to zero.
    """
    path = Path(path).resolve()
    n_nodes = mesh.n_nodes
    n_elements = mesh.n_elements
    n_open_segs = len(mesh.open_boundaries)
    n_open_nodes = sum(len(b) for b in mesh.open_boundaries)
    n_land_segs = len(mesh.land_boundaries)
    n_land_nodes = sum(len(ids) for _, ids in mesh.land_boundaries)

    with path.open("w") as f:
        f.write(f"{mesh.title}\n")
        f.write(f"{n_elements} {n_nodes}\n")

        for i in range(n_nodes):
            x, y = mesh.nodes[i]
            f.write(f"{i + 1:>10d}  {x:.15f}  {y:.15f}  {mesh.depths[i]:.10e}\n")

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
