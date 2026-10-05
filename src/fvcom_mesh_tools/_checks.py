"""Argument checks shared by the mesh tools."""

from __future__ import annotations

import numbers

import numpy as np

_SCALARS = (numbers.Number, str, bytes, type(None), np.generic)


def _leaves(root):
    """Yield the leaves of ``root``: arrays of a plain dtype whole (a masked array as
    itself, so its mask can be seen), and the elements of lists, tuples and object arrays
    one by one, to any depth. Anything array-like (a memoryview, a pandas object) is
    brought to an array first. Iterative, each container visited once, and no array is
    built from a nested list, so shared or self-holding nesting costs nothing (review
    rounds 62-65)."""
    stack, visited = [root], set()
    while stack:
        obj = stack.pop()
        if isinstance(obj, _SCALARS):
            yield obj
            continue
        if isinstance(obj, (list, tuple)):
            if id(obj) in visited:
                continue
            visited.add(id(obj))
            items = obj
        else:
            arr = obj if isinstance(obj, np.ndarray) else np.asarray(obj)
            if arr.dtype.kind != "O":
                yield arr
                continue
            if id(obj) in visited:
                continue
            visited.add(id(obj))
            items = arr.ravel()
        for v in items:
            if isinstance(v, _SCALARS):
                yield v
            else:
                stack.append(v)


def no_masked(**arrays) -> None:
    """Refuse masked values wherever they sit, also in a masked array inside a list:
    a float conversion would use the data under the mask (review round 65 F3)."""
    for name, a in arrays.items():
        if a is None:
            continue
        for leaf in _leaves(a):
            if isinstance(leaf, np.ndarray) and np.ma.is_masked(leaf) \
                    or leaf is np.ma.masked:
                raise ValueError(f"{name} has masked values")


def no_complex(**arrays) -> None:
    """Refuse complex values, wherever they sit: a float conversion would keep only the
    real part and go on with altered data (review rounds 36 F2, 37-38, 64 F4)."""
    for name, a in arrays.items():
        if a is None:
            continue
        for leaf in _leaves(a):
            if isinstance(leaf, np.ndarray):
                if leaf.dtype.kind == "c":
                    raise ValueError(f"{name} has complex values")
            elif isinstance(leaf, (complex, np.complexfloating)):
                raise ValueError(f"{name} has complex values")


def real_scalar(value, name: str) -> float:
    """A finite real number (NumPy scalars too; not a bool, not complex) as a
    float (review round 33 F5)."""
    if isinstance(value, (bool, np.bool_)) or not isinstance(
            value, (int, float, np.integer, np.floating)):
        raise ValueError(f"{name} must be a real number, not {value!r}")
    v = float(value)
    if not np.isfinite(v):
        raise ValueError(f"{name} must be finite, not {value!r}")
    return v


def positive_whole(value, name: str) -> int:
    """A positive whole number, not 1.9 and not a bool (review round 33 F5)."""
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) \
            or value < 1:
        raise ValueError(f"{name} must be a positive whole number, not {value!r}")
    return int(value)


def checked_geometry(nodes, elements, what: str = "the mesh") -> tuple[np.ndarray, np.ndarray]:
    """``(xy, tri)``: finite float coordinates (N, >= 2) and whole in-range
    triangle indices (M, 3), refusing masks and complex values before any
    conversion (review round 39 F2, F3)."""
    from fvcom_mesh_tools.io.fvcom_native import _indices

    for name, a in (("nodes", nodes), ("elements", elements)):
        no_masked(**{f"{what}: {name}": a})
    no_complex(nodes=nodes, elements=elements)
    xy = np.asarray(nodes, float)
    if xy.ndim != 2 or xy.shape[1] < 2 or not np.isfinite(xy).all():
        raise ValueError(f"{what}: nodes must be finite (N, >= 2), not {xy.shape}")
    # lengths, areas and cross products of coordinates this large overflow
    # (review round 44 F2); metres and degrees are many orders below
    if np.abs(xy).max(initial=0.0) > 1e100:
        raise ValueError(f"{what}: coordinates beyond 1e100 are not supported")
    tri = _indices(elements, len(xy), f"{what}: elements", ndim=2)
    if tri.shape[1:] != (3,):
        raise ValueError(f"{what}: elements must be (M, 3), not {tri.shape}")
    return xy, tri


def checked_flags(flags, n: int, name: str) -> np.ndarray:
    """A boolean ``(n,)`` permission array; a mask or another dtype or shape is
    refused (a masked permission is no permission; review round 39 F1)."""
    if np.ma.is_masked(flags):
        raise ValueError(f"{name} has masked values")
    arr = np.asarray(flags)
    if arr.dtype != bool or arr.shape != (n,):
        raise ValueError(f"{name} must be a boolean ({n},) array, not {arr.dtype} {arr.shape}")
    return arr


def promoted_nodes(nodes, what: str = "nodes") -> np.ndarray:
    """``nodes`` as a float64 array, refusing masks and complex values first, so
    that differences and cross products are not computed in a narrow dtype and
    no unknown or imaginary part is dropped (review round 45)."""
    no_masked(**{what: nodes})
    no_complex(**{what: nodes})
    return np.asarray(nodes, dtype=np.float64)


def checked_planar(nodes, elements, what: str = "the mesh") -> tuple[np.ndarray, np.ndarray]:
    """``checked_geometry`` for the planar APIs: exactly two coordinate columns,
    so that no helper measures some lengths in two dimensions and others in
    three, and no column is silently dropped (review round 47 F2)."""
    xy, tri = checked_geometry(nodes, elements, what)
    if xy.shape[1] != 2:
        raise ValueError(f"{what}: nodes must be (N, 2), not {xy.shape}")
    return xy, tri


def no_bool(**arrays) -> None:
    """Refuse booleans, also inside object arrays and nested lists/arrays to any depth:
    ``True`` would become the number 1 in a float conversion (review rounds 55-64)."""
    for name, a in arrays.items():
        if a is None:
            continue
        for leaf in _leaves(a):
            if isinstance(leaf, np.ndarray):
                if leaf.dtype.kind == "b":
                    raise ValueError(f"{name} has boolean values, not numbers")
            elif isinstance(leaf, (bool, np.bool_)):
                raise ValueError(f"{name} has boolean values, not numbers")


def is_finite_real(v) -> bool:
    """A finite real number, also an integer too large for a float (which is refused
    rather than crashing the check; review round 58 F5). Not a bool."""
    import math

    if isinstance(v, (bool, np.bool_)):
        return False
    if isinstance(v, (int, np.integer)):
        return abs(int(v)) < 10 ** 300
    if isinstance(v, (float, np.floating)):
        return math.isfinite(v)
    return False
