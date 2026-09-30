"""Harmonic constants from ocean tide models, at arbitrary points.

Three formats are read: NAO.99 (NAO.99Jb, regional around Japan at 1/12
deg, and NAO.99b, global at 0.5 deg; Matsumoto et al., 2000), the netCDF of
the OTIS/TPXO atlases (TPXO10-atlas-v2, global at 1/30 deg; Egbert and
Erofeeva, 2002), and the FES2022 netCDF (FES2022b from AVISO+, global at
1/30 deg, extrapolated onto land near the coast). The files live under
``$DATA_DIR/tides/models`` (see the README there). All come out as the same
grid dictionary -- ``lon`` and ``lat`` ascending, ``amp`` (m) and ``phase``
(Greenwich lag, deg) indexed ``[lat, lon]``, NaN where there is no ocean --
which ``sample_constants`` takes.

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
    "EQUILIBRIUM", "NAO_CONSTITUENTS", "astronomy", "fvcom_spectral", "load_nao", "load_tpxo",
    "load_fes", "read_fes", "read_nao", "read_tpxo", "sample_constants", "spectral_text",
]

#: The 16 constituents of NAO.99 and their file stems.
NAO_CONSTITUENTS = {
    "M2": "m2", "S2": "s2", "N2": "n2", "K2": "k2", "K1": "k1", "O1": "o1",
    "P1": "p1", "Q1": "q1", "M1": "m1", "J1": "j1", "OO1": "oo1", "2N2": "2n2",
    "MU2": "mu2", "NU2": "nu2", "L2": "l2", "T2": "t2",
}
_MISSING = 999999

#: Equilibrium tide per constituent: amplitude (m, Cartwright and Tayler,
#: 1971), the elasticity factor beta = 1 + k - h (Wahr, 1981), and the
#: species, as FVCOM's -DEQUI_TIDE reads them. The values are the ones
#: ADCIRC tables for its tidal potential.
EQUILIBRIUM = {
    "M2": (0.242334, 0.693, "SEMIDIURNAL"), "S2": (0.112841, 0.693, "SEMIDIURNAL"),
    "N2": (0.046398, 0.693, "SEMIDIURNAL"), "K2": (0.030704, 0.693, "SEMIDIURNAL"),
    "K1": (0.141565, 0.736, "DIURNAL"), "O1": (0.100514, 0.695, "DIURNAL"),
    "P1": (0.046843, 0.706, "DIURNAL"), "Q1": (0.019256, 0.695, "DIURNAL"),
}


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


def read_tpxo(path, window=None) -> dict:
    """One TPXO atlas elevation file (``h_<c>_*.nc``), optionally a lon/lat window.

    The file holds the complex amplitude ``hRe + i hIm`` in millimetres on
    ``(nx, ny)``; amplitude ``|h|`` and Greenwich phase ``atan2(-hIm, hRe)``,
    as its own attributes say. Land is exactly 0. ``window`` is
    ``(lon0, lon1, lat0, lat1)`` in degrees east (0-360 as the file) and
    north; reading only a window keeps a 1/30-deg global file cheap.
    """
    import netCDF4

    with netCDF4.Dataset(path) as ds:
        lon = np.asarray(ds["lon_z"][:], float)
        lat = np.asarray(ds["lat_z"][:], float)
        if window is None:
            ix, iy = slice(None), slice(None)
        else:
            lon0, lon1, lat0, lat1 = window
            sel_x = np.flatnonzero((lon >= lon0) & (lon <= lon1))
            sel_y = np.flatnonzero((lat >= lat0) & (lat <= lat1))
            if len(sel_x) < 2 or len(sel_y) < 2:
                raise ValueError(f"{path}: window {window} holds fewer than 2x2 nodes")
            ix = slice(sel_x[0], sel_x[-1] + 1)
            iy = slice(sel_y[0], sel_y[-1] + 1)
        re_ = np.asarray(ds["hRe"][ix, iy], float).T
        im_ = np.asarray(ds["hIm"][ix, iy], float).T
        name = bytes(np.asarray(ds["con"][:])).decode(errors="ignore").strip().upper() \
            if "con" in ds.variables else None
    lon, lat = lon[ix], lat[iy]
    if not (np.all(np.diff(lon) > 0) and np.all(np.diff(lat) > 0)):
        raise ValueError(f"{path}: coordinates are not ascending")
    land = (re_ == 0) & (im_ == 0)
    amp = np.where(land, np.nan, np.hypot(re_, im_) / 1000.0)
    pha = np.where(land, np.nan, np.degrees(np.arctan2(-im_, re_)) % 360)
    return {"constituent": name, "lon": lon, "lat": lat, "amp": amp, "phase": pha}


def load_tpxo(directory, constituents, window=None,
              pattern="h_{c}_tpxo10_atlas_30_v2.nc") -> dict[str, dict]:
    """Read several constituents from a TPXO atlas directory (``pattern`` names a file)."""
    out = {}
    for c in constituents:
        path = Path(directory) / pattern.format(c=c.lower())
        if not path.exists():
            raise FileNotFoundError(path)
        out[c.upper()] = read_tpxo(path, window)
    return out


def read_fes(path, window=None) -> dict:
    """One FES2022 file (``<c>_fes2022.nc``), optionally a lon/lat window.

    ``amplitude`` (cm) and ``phase`` (Greenwich lag, deg) on ``(lat, lon)``,
    NaN (or the fill value) where there is no ocean. ``window`` is
    ``(lon0, lon1, lat0, lat1)``, degrees east 0-360 as the file.
    """
    import netCDF4

    with netCDF4.Dataset(path) as ds:
        lon = np.asarray(ds["lon"][:], float)
        lat = np.asarray(ds["lat"][:], float)
        if window is None:
            ix, iy = slice(None), slice(None)
        else:
            lon0, lon1, lat0, lat1 = window
            sel_x = np.flatnonzero((lon >= lon0) & (lon <= lon1))
            sel_y = np.flatnonzero((lat >= lat0) & (lat <= lat1))
            if len(sel_x) < 2 or len(sel_y) < 2:
                raise ValueError(f"{path}: window {window} holds fewer than 2x2 nodes")
            ix = slice(sel_x[0], sel_x[-1] + 1)
            iy = slice(sel_y[0], sel_y[-1] + 1)
        amp = np.ma.filled(np.ma.asarray(ds["amplitude"][iy, ix], float), np.nan)
        pha = np.ma.filled(np.ma.asarray(ds["phase"][iy, ix], float), np.nan)
        units = getattr(ds["amplitude"], "units", "cm")
    lon, lat = lon[ix], lat[iy]
    if not (np.all(np.diff(lon) > 0) and np.all(np.diff(lat) > 0)):
        raise ValueError(f"{path}: coordinates are not ascending")
    if units != "cm":
        raise ValueError(f"{path}: amplitude in {units!r}, expected 'cm'")
    bad = ~np.isfinite(amp) | ~np.isfinite(pha) | (np.abs(amp) > 1e5)
    amp = np.where(bad, np.nan, amp / 100.0)
    pha = np.where(bad, np.nan, pha % 360)
    return {"constituent": Path(path).name.split("_")[0].upper(), "lon": lon, "lat": lat,
            "amp": amp, "phase": pha}


def load_fes(directory, constituents, window=None, pattern="{c}_fes2022.nc") -> dict[str, dict]:
    """Read several constituents from a FES2022 directory (``pattern`` names a file)."""
    out = {}
    for c in constituents:
        path = Path(directory) / pattern.format(c=c.lower())
        if not path.exists():
            raise FileNotFoundError(path)
        out[c.upper()] = read_fes(path, window)
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


def astronomy(names, start: datetime, mid: datetime, lat: float):
    """Period (s), nodal factor f, and V0 + u (deg) of each constituent.

    ``V0`` is the Greenwich astronomical argument at ``start``; ``f`` and
    ``u`` are frozen at ``mid``. Astronomy is utide's (``FUV``, Greenwich,
    exact nodal corrections), so a harmonic analysis with utide sees the same
    convention. Times are UTC.
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
    period = 3600.0 / ut_constants.const.freq[lind]
    return period, f.ravel(), (v0.ravel() + u.ravel()) * 360.0 % 360.0


def fvcom_spectral(names, amp, phase, start: datetime, mid: datetime, lat: float):
    """Greenwich constants -> FVCOM spectral forcing (period s, amplitude m, phase deg).

    FVCOM's spectral open-boundary forcing is ``A cos(2 pi t / T - phi)`` with
    ``t`` in seconds since the file's Time Origin, so the astronomical
    argument has to go into the phase: ``phi = G - V0(start) - u`` and
    ``A = f * amp``. ``V0`` is taken at ``start`` (the Time Origin); the
    nodal factors ``f`` and ``u`` are frozen at ``mid`` -- FVCOM cannot vary
    them, and over a few months they change by a few percent / degrees at
    most (K1, O1). See :func:`astronomy`.

    ``amp`` and ``phase`` have shape ``(n_names, n_points)``; times are UTC.
    """
    period, f, v0u = astronomy(names, start, mid, lat)
    amp = np.asarray(amp, float)
    phase = np.asarray(phase, float)
    return period, f[:, None] * amp, (phase - v0u[:, None]) % 360.0


def spectral_text(names, period, amp, phase, origin: str, equilibrium=None,
                  sal_beta: float | None = None) -> str:
    """An FVCOM non-Julian (spectral) tidal forcing file, several constituents.

    ``amp`` and ``phase`` have shape ``(n_names, n_obc)``, in open-boundary
    order; the first column of each row is the open-boundary ordinal, not the
    mesh node number.

    ``equilibrium=(f, v0u)`` adds, to each component line, what an FVCOM built
    with ``-DEQUI_TIDE`` reads: the equilibrium amplitude, beta and species
    (:data:`EQUILIBRIUM`), then the nodal factor and V0 + u (deg) at the Time
    Origin -- the owner's FVCOM extension, which makes the tidal potential use
    the same astronomy as the open boundary.

    ``sal_beta`` adds a ``SAL Beta`` line: the self-attraction and loading
    in the scalar approximation (``beta * zeta`` added to the equilibrium
    tide; typically 0.08-0.12), which the same extension reads. It needs
    ``equilibrium``.
    """
    amp = np.atleast_2d(np.asarray(amp, float))
    phase = np.atleast_2d(np.asarray(phase, float))
    if amp.shape != phase.shape or not amp.shape[0] == len(names) == len(period):
        raise ValueError("names, period, amp and phase disagree in shape")
    if not (np.isfinite(amp).all() and np.isfinite(phase).all()):
        raise ValueError("non-finite amplitude or phase")
    if sal_beta is not None and (equilibrium is None or not 0.0 <= sal_beta < 0.5):
        raise ValueError("sal_beta needs equilibrium and must be in [0, 0.5)")
    n = amp.shape[1]
    lines = [f"Tidal Component Number = {len(names)}"]
    for i, (c, p) in enumerate(zip(names, period), 1):
        line = f"{i} = {c} {p:.10f}"
        if equilibrium is not None:
            if c.upper() not in EQUILIBRIUM:
                raise ValueError(f"no equilibrium tide for {c}; known {sorted(EQUILIBRIUM)}")
            a_eq, beta, kind = EQUILIBRIUM[c.upper()]
            f, v0u = equilibrium[0][i - 1], equilibrium[1][i - 1]
            line += f" {a_eq:.6f} {beta:.3f} {kind} {f:.8f} {v0u:.6f}"
        lines.append(line)
    if sal_beta is not None:
        lines.append(f"SAL Beta = {sal_beta:.6f}")
    lines += [f"Time Origin = {origin}", f"OBC Node Number = {n}"]
    for label, values in (("Amplitude", amp), ("Phase", phase), ("Eref", np.zeros((1, n)))):
        lines.append(label)
        lines += [f"{j + 1} " + " ".join(f"{v:.8f}" for v in values[:, j]) for j in range(n)]
        lines.append(label)
    return "\n".join(lines) + "\n"
