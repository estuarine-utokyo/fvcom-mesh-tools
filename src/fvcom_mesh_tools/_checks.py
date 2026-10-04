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
