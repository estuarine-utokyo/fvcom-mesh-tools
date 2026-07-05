"""PoC #100 (stage 1) — M2 tide case for the v5 pipeline mesh.

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
MESH = REPO / "outputs" / "pipeline_v5" / "tokyo_bay_v5_final.14"
V2_INPUTS = REPO / "outputs" / "pipeline_v5" / "fvcom_inputs"
TIDE_NML = REPO / "outputs" / "fvcom_tide_v4" / "tokyo_bay_v4_tide_run.nml"
CASE_DIR = REPO / "outputs" / "fvcom_tide_v5"
CASENAME = "tokyo_bay_v5_tide"

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
    for f in V2_INPUTS.glob("tokyo_bay_v5_*.dat"):
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
    print(f"[100] wrote {nc_path}  (nobc={nobc}, nt={t.size})", flush=True)

    # --- namelist ----------------------------------------------------------
    txt = TIDE_NML.read_text()
    txt = txt.replace("tokyo_bay_v4_tide", CASENAME)
    txt = txt.replace("tokyo_bay_v4", "tokyo_bay_v4")
    nml = CASE_DIR / f"{CASENAME}_run.nml"
    nml.write_text(txt)
    print(f"[100] wrote {nml}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
