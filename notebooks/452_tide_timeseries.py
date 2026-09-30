"""Tide runs against tide gauges as time series, over a short window.

    python notebooks/452_tide_timeseries.py --run enshu_eq=$ROOT/enshu --run base_eq=$ROOT/base \\
        --out $ROOT/timeseries.png [--start 2021-03-01 --days 14] [--stations tokyo,yokosuka,mera]

The model's ``zeta`` is taken at each gauge as 450 takes it (linear in its
triangle, else the nearest node). The observation carries the mean sea level
on T.P. and the non-tidal signal (weather, seasonal), the model neither: both
are shown minus their own mean over the window, and the lower panel of each
gauge shows observation minus model.
"""

from __future__ import annotations

import argparse
import os
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
from pyproj import Transformer

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fvcom_mesh_tools.io.fvcom_native import read_fvcom_case  # noqa: E402

MESH_EPSG = 32654
MJD0 = datetime(1858, 11, 17)

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--run", action="append", required=True, help="LABEL=RUN_DIR")
p.add_argument("--out", type=Path, required=True)
p.add_argument("--start", default="2021-03-01")
p.add_argument("--days", type=float, default=14.0)
p.add_argument("--stations", default="tokyo,yokosuka,mera")
a = p.parse_args()
if not os.environ.get("DATA_DIR"):
    raise SystemExit("set DATA_DIR")
GAUGES = Path(os.environ["DATA_DIR"]) / "tides/nc_raw"

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402
import netCDF4  # noqa: E402

t0 = datetime.fromisoformat(a.start)
t1 = t0 + timedelta(days=a.days)
stations = a.stations.split(",")
to_utm = Transformer.from_crs(4326, MESH_EPSG, always_xy=True)


def observed(station):
    f = GAUGES / f"{station}_{t0.year}.nc"
    with netCDF4.Dataset(f) as ds:
        var = "sea_level" if "sea_level" in ds.variables else "ssh"
        h = np.asarray(ds["time"][:], float)
        z = np.ma.filled(np.ma.asarray(ds[var][:], float), np.nan)
        lon, lat, name = float(ds.longitude), float(ds.latitude), ds.getncattr("station_name")
    t = np.array([datetime(t0.year, 1, 1) + timedelta(hours=float(x)) for x in h])
    keep = (t >= t0) & (t <= t1)
    return t[keep], z[keep], lon, lat, name


obs = {s: observed(s) for s in stations}
model = {}
for item in a.run:
    label, _, path = item.partition("=")
    run = Path(path)
    m = read_fvcom_case(run / "input/m2_grd.dat", run / "input/m2_dep.dat",
                        run / "input/m2_obc.dat", title=label)
    finder = mtri.Triangulation(m.nodes[:, 0], m.nodes[:, 1], m.elements).get_trifinder()
    with netCDF4.Dataset(run / "output/m2_0001.nc") as ds:
        days = np.asarray(ds["time"][:], float)
        tm = np.array([MJD0 + timedelta(days=float(d)) for d in days])
        keep = (tm >= t0) & (tm <= t1)
        zeta = np.asarray(ds["zeta"][np.flatnonzero(keep), :], float)
    tm = tm[keep]
    for s in stations:
        _, _, lon, lat, _ = obs[s]
        x, y = to_utm.transform(lon, lat)
        e = int(finder(x, y))
        if e >= 0:
            nodes = m.elements[e]
            lam = np.linalg.solve(np.r_[m.nodes[nodes, :2].T, np.ones((1, 3))], [x, y, 1.0])
            z = zeta[:, nodes] @ lam
            where = "inside"
        else:
            d = np.hypot(m.nodes[:, 0] - x, m.nodes[:, 1] - y)
            z = zeta[:, int(d.argmin())]
            where = f"nearest node {d.min() / 1e3:.1f} km"
        model[(label, s)] = (tm, z, where)

labels = [item.partition("=")[0] for item in a.run]
colors = dict(zip(labels, ["tab:green", "tab:red", "tab:blue", "tab:orange"]))
fig, axes = plt.subplots(2 * len(stations), 1, figsize=(14, 3.6 * len(stations)), sharex=True,
                         gridspec_kw={"height_ratios": [3, 1] * len(stations)})
for i, s in enumerate(stations):
    ax, axr = axes[2 * i], axes[2 * i + 1]
    t, zo, _, _, name = obs[s]
    ok = np.isfinite(zo)
    ax.plot(t[ok], (zo[ok] - zo[ok].mean()) * 100, color="k", lw=1.6, label="observed")
    for lab in labels:
        tm, z, where = model[(lab, s)]
        zm = (z - z.mean()) * 100
        ax.plot(tm, zm, color=colors[lab], lw=1.0, label=f"{lab} ({where})")
        common, io, im = np.intersect1d(t[ok], tm, return_indices=True)
        r = (zo[ok][io] - zo[ok].mean()) * 100 - zm[im]
        axr.plot(common, r, color=colors[lab], lw=0.9,
                 label=f"{lab}: rms {np.sqrt(np.mean(r ** 2)):.1f} cm")
    ax.set_ylabel("sea level (cm)")
    ax.set_title(f"{name}: minus window mean", fontsize=10, loc="left")
    ax.legend(fontsize=8, ncol=len(labels) + 1, loc="upper right")
    axr.axhline(0, color="0.6", lw=0.6)
    axr.set_ylabel("obs - model (cm)", fontsize=8)
    axr.legend(fontsize=8, ncol=len(labels), loc="upper right")
    ax.grid(alpha=0.3)
    axr.grid(alpha=0.3)
axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
axes[-1].set_xlabel(f"{t0:%Y} (UTC)")
fig.suptitle(f"NAO.99Jb tide runs vs tide gauges, {t0:%Y-%m-%d} .. {t1:%Y-%m-%d}")
fig.tight_layout()
a.out.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(a.out, dpi=120)
print(f"[452] wrote {a.out}", flush=True)
