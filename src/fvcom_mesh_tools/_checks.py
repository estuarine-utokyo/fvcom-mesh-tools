"""Argument checks shared by the mesh tools."""

from __future__ import annotations

import numpy as np


def no_complex(**arrays) -> None:
    """Refuse complex values: a float conversion would keep only the real
    part and go on with altered data (review rounds 36 F2, 37 F1-F2)."""
    for name, a in arrays.items():
        if a is None:
            continue
        arr = np.asarray(a)
        # an object array can hold complex numbers under a dtype that says
        # nothing about them (review round 38 F2)
        if arr.dtype.kind == "c" or (arr.dtype.kind == "O" and any(
                isinstance(v, (complex, np.complexfloating)) for v in arr.ravel())):
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
        if np.ma.is_masked(a):
            raise ValueError(f"{what}: {name} has masked values")
    no_complex(nodes=nodes, elements=elements)
    xy = np.asarray(nodes, float)
    if xy.ndim != 2 or xy.shape[1] < 2 or not np.isfinite(xy).all():
        raise ValueError(f"{what}: nodes must be finite (N, >= 2), not {xy.shape}")
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
