"""Harmonic constants from ocean tide models, at arbitrary points.

Only the NAO.99 format is read for now (NAO.99Jb, regional around Japan at
1/12 deg, and NAO.99b, global at 0.5 deg; Matsumoto et al., 2000). The
files live under ``$DATA_DIR/tides/models`` (see the README there).

The NAO format, as checked against the coastline: seven header lines, then
for each latitude row from NORTH to south the row of amplitudes (0.01 cm)
followed by the row of phases (0.01 deg), ``(10i6)``, 999999 where there is
no ocean. Reading every amplitude first and every phase after gives a
striped, wrong field -- which is why the reader checks the counts.

Phases are Greenwich phase lags. Interpolation works on the complex
constant ``A exp(-i G)``, so a phase that wraps through 360 deg between two
grid nodes is not averaged to 180 deg. A point whose four surrounding nodes
are not all ocean takes the constant of the nearest ocean node within
``fill_cells`` cells: open boundaries end on the coast, where the model grid
is land.
"""

from __future__ import annotations

import re
from datetime import datetime
from pathlib import Path

import numpy as np

__all__ = [
    "NAO_CONSTITUENTS", "fvcom_spectral", "load_nao", "read_nao", "sample_constants",
    "spectral_text",
]

#: The 16 constituents of NAO.99 and their file stems.
NAO_CONSTITUENTS = {
    "M2": "m2", "S2": "s2", "N2": "n2", "K2": "k2", "K1": "k1", "O1": "o1",
    "P1": "p1", "Q1": "q1", "M1": "m1", "J1": "j1", "OO1": "oo1", "2N2": "2n2",
    "MU2": "mu2", "NU2": "nu2", "L2": "l2", "T2": "t2",
}
_MISSING = 999999


def read_nao(path) -> dict:
    """One NAO.99 constituent file: lon (ascending), lat (ascending), amp (m), phase (deg)."""
    lines = Path(path).read_text().splitlines()
    head = " ".join(lines[:7])

    def num(key):
        m = re.search(rf"{key}\s*=\s*([-0-9./ ]+?)(?=\s+[a-zA-Z]|$)", head)
        if not m:
            raise ValueError(f"{path}: no '{key}' in the header")
        return m.group(1).strip()

    def frac(text):
        a, _, b = text.partition("/")
        return float(a) / float(b) if b else float(a)

    xmin, ymax = float(num("xmin")), float(num("ymax"))
    dx, dy = frac(num("dx")), frac(num("dy"))
    nx, ny = int(num("mend")), int(num("nend"))
    content = re.search(r"Content\s*:\s*(\S+)", head)
    values: list[int] = []
    for line in lines[7:]:
        line = line.rstrip()
        values.extend(int(line[i:i + 6]) for i in range(0, len(line), 6) if line[i:i + 6].strip())
    if len(values) != 2 * nx * ny:
        raise ValueError(f"{path}: {len(values)} values, expected {2 * nx * ny}")
    rows = np.asarray(values, float).reshape(2 * ny, nx)
    amp, pha = rows[0::2], rows[1::2]          # interleaved per latitude row
    miss = (amp == _MISSING) | (pha == _MISSING)
    amp = np.where(miss, np.nan, amp / 1e4)   # 0.01 cm -> m
    pha = np.where(miss, np.nan, pha / 100)
    lat = ymax - np.arange(ny) * dy            # row 0 is the north
    return {
        "constituent": content.group(1).upper() if content else None,
        "lon": xmin + np.arange(nx) * dx,
        "lat": lat[::-1], "amp": amp[::-1], "phase": pha[::-1],
    }


def load_nao(directory, constituents, suffix="_j") -> dict[str, dict]:
    """Read several constituents from a NAO.99 directory (``_j`` = NAO.99Jb)."""
    out = {}
    for c in constituents:
        key = c.upper()
        if key not in NAO_CONSTITUENTS:
            raise ValueError(f"NAO.99 has no constituent {c}; it has {sorted(NAO_CONSTITUENTS)}")
        path = Path(directory) / f"{NAO_CONSTITUENTS[key]}{suffix}.nao"
        if not path.exists():
            raise FileNotFoundError(path)
        out[key] = read_nao(path)
    return out


def sample_constants(grid: dict, lon, lat, fill_cells: int = 3):
    """Amplitude (m) and Greenwich phase (deg) at points; NaN where no ocean is near.

    Returns ``(amp, phase, filled)``; ``filled`` marks points that took the
    nearest ocean node because their cell was not all ocean.
    """
    from scipy.spatial import cKDTree

    lon = np.asarray(lon, float)
    lat = np.asarray(lat, float)
    glon, glat = grid["lon"], grid["lat"]
    z = grid["amp"] * np.exp(-1j * np.radians(grid["phase"]))
    dx, dy = glon[1] - glon[0], glat[1] - glat[0]
    fi = (lon - glon[0]) / dx
    fj = (lat - glat[0]) / dy
    i0 = np.clip(np.floor(fi).astype(int), 0, len(glon) - 2)
    j0 = np.clip(np.floor(fj).astype(int), 0, len(glat) - 2)
    ti, tj = fi - i0, fj - j0
    inside = (fi >= 0) & (fi <= len(glon) - 1) & (fj >= 0) & (fj <= len(glat) - 1)
    val = (z[j0, i0] * (1 - ti) * (1 - tj) + z[j0, i0 + 1] * ti * (1 - tj)
           + z[j0 + 1, i0] * (1 - ti) * tj + z[j0 + 1, i0 + 1] * ti * tj)
    val = np.where(inside, val, np.nan)
    filled = ~np.isfinite(val) & inside
    if filled.any():
        jj, ii = np.nonzero(np.isfinite(z))
        tree = cKDTree(np.column_stack([ii, jj]))
        d, k = tree.query(np.column_stack([fi[filled], fj[filled]]))
        ok = d <= fill_cells
        got = np.full(int(filled.sum()), np.nan, complex)
        got[ok] = z[jj[k[ok]], ii[k[ok]]]
        val[filled] = got
        filled[np.flatnonzero(filled)[~ok]] = False
    return np.abs(val), np.degrees(-np.angle(val)) % 360, filled


def fvcom_spectral(names, amp, phase, start: datetime, mid: datetime, lat: float):
    """Greenwich constants -> FVCOM spectral forcing (period s, amplitude m, phase deg).

    FVCOM's spectral open-boundary forcing is ``A cos(2 pi t / T - phi)`` with
    ``t`` in seconds since the file's Time Origin, so the astronomical
    argument has to go into the phase: ``phi = G - V0(start) - u`` and
    ``A = f * amp``. ``V0`` is taken at ``start`` (the Time Origin); the
    nodal factors ``f`` and ``u`` are frozen at ``mid`` -- FVCOM cannot vary
    them, and over a few months they change by a few percent / degrees at
    most (K1, O1). Astronomy is utide's (``FUV``, Greenwich, exact nodal),
    so that a harmonic analysis with utide sees the same convention.

    ``amp`` and ``phase`` have shape ``(n_names, n_points)``; times are UTC.
    """
    from utide._ut_constants import ut_constants
    from utide.harmonics import FUV

    all_names = list(ut_constants.const.name)
    try:
        lind = np.array([all_names.index(n.upper()) for n in names])
    except ValueError as err:
        raise ValueError(f"utide does not know a constituent in {list(names)}") from err
    flags = np.array([0, 0, 0, 0])
    t0, tm = (np.array([d.toordinal() + (d - datetime(d.year, d.month, d.day)).total_seconds()
                        / 86400.0]) for d in (start, mid))
    _, _, v0 = FUV(t0, t0[0], lind, float(lat), flags)
    f, u, _ = FUV(tm, tm[0], lind, float(lat), flags)
    f, u, v0 = f.ravel(), u.ravel() * 360.0, v0.ravel() * 360.0
    period = 3600.0 / ut_constants.const.freq[lind]
    amp = np.asarray(amp, float)
    phase = np.asarray(phase, float)
    return (period, f[:, None] * amp,
            (phase - v0[:, None] - u[:, None]) % 360.0)


def spectral_text(names, period, amp, phase, origin: str) -> str:
    """An FVCOM non-Julian (spectral) tidal forcing file, several constituents.

    ``amp`` and ``phase`` have shape ``(n_names, n_obc)``, in open-boundary
    order; the first column of each row is the open-boundary ordinal, not the
    mesh node number.
    """
    amp = np.atleast_2d(np.asarray(amp, float))
    phase = np.atleast_2d(np.asarray(phase, float))
    if amp.shape != phase.shape or not amp.shape[0] == len(names) == len(period):
        raise ValueError("names, period, amp and phase disagree in shape")
    if not (np.isfinite(amp).all() and np.isfinite(phase).all()):
        raise ValueError("non-finite amplitude or phase")
    n = amp.shape[1]
    lines = [f"Tidal Component Number = {len(names)}"]
    lines += [f"{i} = {c} {p:.10f}" for i, (c, p) in enumerate(zip(names, period), 1)]
    lines += [f"Time Origin = {origin}", f"OBC Node Number = {n}"]
    for label, values in (("Amplitude", amp), ("Phase", phase), ("Eref", np.zeros((1, n)))):
        lines.append(label)
        lines += [f"{j + 1} " + " ".join(f"{v:.8f}" for v in values[:, j]) for j in range(n)]
        lines.append(label)
    return "\n".join(lines) + "\n"
