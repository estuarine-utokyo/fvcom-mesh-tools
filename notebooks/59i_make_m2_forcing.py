"""PoC #59i (stage 1) — synthetic M2 OBC elevation forcing + tide case.

Builds ``outputs/fvcom_tide/`` for the tokyo_bay_v2 M2 test:

* ``tokyo_bay_v2_m2.nc`` — FVCOM "julian" time-series elevation
  forcing (format mirrors TB-FVCOM's production
  ``tb_julian_obc_tidegauge_2020.nc``): uniform M2,
  amplitude 0.40 m (≈ the Uraga/Sagami M2 amplitude), 10-min
  sampling, MJD time axis, obc_nodes = the v2 open-boundary list.
* the run namelist, cloned from the #59g smoke nml: forcing ON,
  3 simulated days (2020-01-01 .. 2020-01-04), 12 h ramp
  (IRAMP = 8640 internal steps at 5 s), hourly NetCDF output.
* ``input/`` populated from ``outputs/fvcom_inputs_v2``.

Success criteria (checked by 59i stage 3): run completes, no NaN,
bay-head M2 amplitude / mouth amplitude in a sane 1.1-1.4 band
(Tokyo vs Uraga), channel velocities O(1 m/s).
"""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

from fvcom_mesh_tools.io import read_fort14

REPO = Path(__file__).resolve().parents[1]
MESH = REPO / "outputs" / "59h_gate_passed.14"
V2_INPUTS = REPO / "outputs" / "fvcom_inputs_v2"
SMOKE_NML = REPO / "outputs" / "fvcom_smoke" / "tokyo_bay_v1_smoke_run.nml"
CASE_DIR = REPO / "outputs" / "fvcom_tide"
CASENAME = "tokyo_bay_v2_tide"

M2_AMP_M = 0.40
M2_PERIOD_H = 12.4206012
MJD_START = 58848.5   # 2019-12-31 12:00 UTC (half a day before START_DATE)
MJD_END = 58853.0     # 2020-01-05 00:00 UTC
DT_MIN = 10.0


def main() -> int:
    mesh = read_fort14(MESH)
    if len(mesh.open_boundaries) != 1:
        raise SystemExit(f"expected 1 open segment, got {len(mesh.open_boundaries)}")
    obc = np.asarray(mesh.open_boundaries[0], dtype=np.int64) + 1  # 1-indexed
    nobc = obc.size

    CASE_DIR.mkdir(parents=True, exist_ok=True)
    (CASE_DIR / "input").mkdir(exist_ok=True)
    (CASE_DIR / "output").mkdir(exist_ok=True)
    for f in V2_INPUTS.glob("tokyo_bay_v2_*.dat"):
        shutil.copy2(f, CASE_DIR / "input" / f.name)
    (CASE_DIR / "input" / "sigma.dat").write_text(
        "NUMBER OF SIGMA LEVELS = 11\nSIGMA COORDINATE TYPE = UNIFORM\n",
    )
    fvcom_bin = Path.home() / "Github" / "TB-FVCOM" / "hydro" / "bin" / "fvcom"
    link = CASE_DIR / "fvcom"
    if not link.exists():
        link.symlink_to(fvcom_bin)

    # --- M2 forcing NetCDF -------------------------------------------------
    t = np.arange(MJD_START, MJD_END + 1e-9, DT_MIN / 1440.0)
    phase = 2.0 * np.pi * ((t - t[0]) * 24.0) / M2_PERIOD_H
    elev = (M2_AMP_M * np.sin(phase))[:, None] * np.ones((1, nobc))

    nc_path = CASE_DIR / "input" / f"{CASENAME}_m2.nc"
    with Dataset(nc_path, "w", format="NETCDF4_CLASSIC") as ds:
        ds.createDimension("nobc", nobc)
        ds.createDimension("time", None)
        ds.createDimension("DateStrLen", 26)
        v = ds.createVariable("obc_nodes", "i4", ("nobc",))
        v.long_name = "Open Boundary Node Number"
        v.grid = "obc_grid"
        v[:] = obc.astype(np.int32)
        v = ds.createVariable("time", "f8", ("time",))
        v.long_name = "time"
        v.units = "days since 1858-11-17 00:00:00"
        v.format = "modified julian day (MJD)"
        v.time_zone = "UTC"
        v[:] = t
        v = ds.createVariable("elevation", "f4", ("time", "nobc"))
        v.long_name = "Open Boundary Elevation"
        v.units = "meters"
        v[:] = elev.astype(np.float32)
        ds.type = "FVCOM TIME SERIES ELEVATION FORCING FILE"
        ds.title = f"Synthetic M2 ({M2_AMP_M} m) for {CASENAME}"
        ds.history = "Created by notebooks/59i_make_m2_forcing.py"
    print(f"[59i] wrote {nc_path}  (nobc={nobc}, nt={t.size})", flush=True)

    # --- namelist ----------------------------------------------------------
    txt = SMOKE_NML.read_text()
    subs = [
        ("CASE_TITLE      = 'tokyo_bay_v1 grid-acceptance smoke (zero forcing)',",
         f"CASE_TITLE      = '{CASENAME} synthetic M2 tidal test',"),
        ("END_DATE        = '2020-01-01 00:02:00'",
         "END_DATE        = '2020-01-04 00:00:00'"),
        (" IRAMP           = 0,", " IRAMP           = 8640,"),
        (" NC_ON             = F,", " NC_ON             = T,"),
        ("NC_FIRST_OUT      = '2020-01-01 00:00:00',",
         "NC_FIRST_OUT      = '2020-01-01 00:00:00',"),
        ("NC_OUT_INTERVAL   = 'days=1.0',",
         "NC_OUT_INTERVAL   = 'seconds=3600.',"),
        ("OBC_ELEVATION_FORCING_ON   = F,", "OBC_ELEVATION_FORCING_ON   = T,"),
        ("OBC_ELEVATION_FILE         = 'none',",
         f"OBC_ELEVATION_FILE         = '{CASENAME}_m2.nc',"),
        ("OBC_NODE_LIST_FILE         = 'tokyo_bay_v1_obc.dat',",
         f"OBC_NODE_LIST_FILE         = '{CASENAME.replace('_tide', '')}_obc.dat',"),
        ("GRID_FILE            = 'tokyo_bay_v1_grd.dat',",
         "GRID_FILE            = 'tokyo_bay_v2_grd.dat',"),
        ("DEPTH_FILE           = 'tokyo_bay_v1_dep.dat',",
         "DEPTH_FILE           = 'tokyo_bay_v2_dep.dat',"),
        ("CORIOLIS_FILE        = 'tokyo_bay_v1_cor.dat',",
         "CORIOLIS_FILE        = 'tokyo_bay_v2_cor.dat',"),
        ("SPONGE_FILE          = 'tokyo_bay_v1_spg.dat'",
         "SPONGE_FILE          = 'tokyo_bay_v2_spg.dat'"),
    ]
    for old, new in subs:
        if old not in txt:
            raise SystemExit(f"[59i] nml pattern missing: {old!r}")
        txt = txt.replace(old, new)
    nml = CASE_DIR / f"{CASENAME}_run.nml"
    nml.write_text(txt)
    print(f"[59i] wrote {nml}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
