"""Harmonic constants of tide runs against tide gauges, same analysis on both.

    python notebooks/450_tide_gauge_compare.py --run enshu=$ROOT/enshu --run base=$ROOT/base \\
        --out $ROOT/compare [--spinup-days 15]

Each run is a directory staged by 449 (``manifest.json``, ``input/m2_grd.dat``,
``output/m2_0001.nc``). Observed hourly sea level comes from
``$DATA_DIR/tides/nc_raw/<station>_<year>.nc`` (T.P., UTC). The model's
``zeta`` is interpolated linearly in its triangle at the gauge (the nearest
node when the gauge falls outside the mesh -- gauges sit on the coast), taken
at the hours the gauge observed, and both series go through the same
``utide.solve`` (OLS, nodal corrections, the run's constituents). Reading
published constants instead would bring in another convention and another
analysis period.

Writes ``constants.csv`` (one row per run, station and constituent),
``summary.json`` and ``compare.png``.
"""

from __future__ import annotations

import argparse
import json
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
MJD0 = datetime(1858, 11, 17).toordinal()
STATIONS = ("tokyo", "chiba", "yokohama", "yokosuka", "kurihama", "aburatsubo", "mera")

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--run", action="append", required=True, help="LABEL=RUN_DIR")
p.add_argument("--out", type=Path, required=True)
p.add_argument("--spinup-days", type=float, default=15.0)
a = p.parse_args()
if not os.environ.get("DATA_DIR"):
    raise SystemExit("set DATA_DIR")
GAUGES = Path(os.environ["DATA_DIR"]) / "tides/nc_raw"
a.out.mkdir(parents=True, exist_ok=True)

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402
import netCDF4  # noqa: E402
import utide  # noqa: E402

runs = {}
for item in a.run:
    label, _, path = item.partition("=")
    runs[label] = Path(path).resolve()
manifests = {k: json.loads((v / "manifest.json").read_text()) for k, v in runs.items()}
names = manifests[next(iter(runs))]["constituents"]
if any(m["constituents"] != names or m["start"] != manifests[next(iter(runs))]["start"]
       for m in manifests.values()):
    raise SystemExit("the runs differ in constituents or start")
start = datetime.fromisoformat(manifests[next(iter(runs))]["start"])
t_lo = start + timedelta(days=a.spinup_days)
t_hi = min(datetime.fromisoformat(m["end"]) for m in manifests.values())
to_utm = Transformer.from_crs(4326, MESH_EPSG, always_xy=True)
print(f"[450] window {t_lo} .. {t_hi} ({(t_hi - t_lo).days} d), constituents {names}",
      flush=True)


def ordinal(dts):
    return np.array([d.toordinal() + (d - datetime(d.year, d.month, d.day)).total_seconds()
                     / 86400 for d in dts])


def observed(station):
    """Hourly observed sea level (m) and Python-ordinal times in the window."""
    t, z = [], []
    for year in range(t_lo.year, t_hi.year + 1):
        f = GAUGES / f"{station}_{year}.nc"
        if not f.exists():
            continue
        with netCDF4.Dataset(f) as ds:
            var = "sea_level" if "sea_level" in ds.variables else "ssh"
            h = np.asarray(ds["time"][:], float)
            v = np.ma.filled(np.ma.asarray(ds[var][:], float), np.nan)
            lat, lon = float(ds.latitude), float(ds.longitude)
        t.append(datetime(year, 1, 1).toordinal() + h / 24)
        z.append(v)
    if not t:
        return None
    t, z = np.concatenate(t), np.concatenate(z)
    lo, hi = ordinal([t_lo, t_hi])
    keep = np.isfinite(z) & (t >= lo) & (t <= hi)
    t, idx = np.unique(np.round(t[keep] * 24) / 24, return_index=True)
    return t, z[keep][idx], lon, lat


def analyse(t, z, lat):
    c = utide.solve(t, z, lat=lat, constit=names, method="ols", conf_int="linear",
                    nodal=True, trend=False, verbose=False, epoch="python")
    return {n: (float(A), float(g), float(e)) for n, A, g, e in zip(c.name, c.A, c.g, c.A_ci)}


rows = []
located = {}
obs_cache = {}
for label, run in runs.items():
    m = read_fvcom_case(run / "input/m2_grd.dat", run / "input/m2_dep.dat",
                        run / "input/m2_obc.dat", title=label)
    tri = mtri.Triangulation(m.nodes[:, 0], m.nodes[:, 1], m.elements)
    finder = tri.get_trifinder()
    with netCDF4.Dataset(run / "output/m2_0001.nc") as ds:
        mt = MJD0 + np.asarray(ds["time"][:], float)
        zeta = np.asarray(ds["zeta"][:], float)
    if not np.isfinite(zeta).all():
        raise SystemExit(f"{label}: non-finite zeta")
    for st in STATIONS:
        if st not in obs_cache:
            obs_cache[st] = observed(st)
        ob = obs_cache[st]
        if ob is None:
            continue
        t, zo, lon, lat = ob
        x, y = to_utm.transform(lon, lat)
        e = int(finder(x, y))
        if e >= 0:
            nodes = m.elements[e]
            lam = np.linalg.solve(np.r_[m.nodes[nodes, :2].T, np.ones((1, 3))], [x, y, 1.0])
            zm_all = zeta[:, nodes] @ lam
            how, dist = "inside", 0.0
        else:
            d = np.hypot(m.nodes[:, 0] - x, m.nodes[:, 1] - y)
            k = int(d.argmin())
            zm_all = zeta[:, k]
            how, dist = "nearest_node", float(d[k])
        located[(label, st)] = {"how": how, "distance_m": dist,
                                "depth_m": float(m.depths[m.elements[e]].mean() if e >= 0
                                                 else m.depths[k])}
        # the model at the hours the gauge observed (outputs are hourly, on the hour)
        pos = np.searchsorted(mt, t - 1e-6)
        ok = (pos < len(mt)) & (np.abs(mt[np.minimum(pos, len(mt) - 1)] - t) < 1e-4)
        tm, zm, zob = t[ok], zm_all[pos[ok]], zo[ok]
        if len(tm) < 24 * 30:
            print(f"[450] {label}/{st}: only {len(tm)} common hours; skipped", flush=True)
            continue
        cm, co = analyse(tm, zm, lat), analyse(tm, zob, lat)
        for n in names:
            Am, gm, _ = cm[n]
            Ao, go, eo = co[n]
            vd = abs(Ao * np.exp(-1j * np.radians(go)) - Am * np.exp(-1j * np.radians(gm)))
            rows.append({"run": label, "station": st, "constituent": n, "hours": int(len(tm)),
                         "obs_amp_m": Ao, "obs_g_deg": go, "obs_amp_ci_m": eo,
                         "model_amp_m": Am, "model_g_deg": gm,
                         "d_amp_m": Am - Ao, "d_g_deg": (gm - go + 180) % 360 - 180,
                         "vector_diff_m": float(vd), **located[(label, st)]})
        print(f"[450] {label}/{st}: {how} ({dist:.0f} m), {len(tm)} h, "
              f"M2 obs {co['M2'][0]:.3f}/{co['M2'][1]:.1f} model {cm['M2'][0]:.3f}/"
              f"{cm['M2'][1]:.1f}; K1 obs {co['K1'][0]:.3f}/{co['K1'][1]:.1f} model "
              f"{cm['K1'][0]:.3f}/{cm['K1'][1]:.1f}", flush=True)

if not rows:
    raise SystemExit("no station could be compared")
keys = list(rows[0])
with open(a.out / "constants.csv", "w") as fh:
    fh.write(",".join(keys) + "\n")
    for r in rows:
        fh.write(",".join(f"{r[k]:.5f}" if isinstance(r[k], float) else str(r[k])
                          for k in keys) + "\n")

summary = {"window": [t_lo.isoformat(), t_hi.isoformat()], "runs": {}}
for label in runs:
    s = {}
    for n in names:
        rr = [r for r in rows if r["run"] == label and r["constituent"] == n]
        s[n] = {"rms_vector_diff_m": float(np.sqrt(np.mean([r["vector_diff_m"] ** 2 for r in rr]))),
                "mean_d_amp_m": float(np.mean([r["d_amp_m"] for r in rr])),
                "mean_d_g_deg": float(np.mean([r["d_g_deg"] for r in rr])),
                "n_stations": len(rr)}
    summary["runs"][label] = {"manifest": manifests[label], "by_constituent": s}
(a.out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
for label in runs:
    print(f"[450] {label}: " + ", ".join(
        f"{n} rmsD {v['rms_vector_diff_m'] * 100:.1f} cm dA {v['mean_d_amp_m'] * 100:+.1f} cm "
        f"dg {v['mean_d_g_deg']:+.1f}"
        for n, v in summary["runs"][label]["by_constituent"].items()), flush=True)

# figure: amplitude and phase per station for the four largest constituents
show = [n for n in ("M2", "S2", "K1", "O1") if n in names]
sts = [s for s in STATIONS if any(r["station"] == s for r in rows)]
fig, ax = plt.subplots(2, len(show), figsize=(4.2 * len(show), 7), squeeze=False)
w = 0.8 / (len(runs) + 1)
for j, n in enumerate(show):
    for i, lab in enumerate(["obs", *runs]):
        amp, g = [], []
        for s in sts:
            r = next((r for r in rows if r["station"] == s and r["constituent"] == n
                      and (lab == "obs" or r["run"] == lab)), None)
            amp.append(np.nan if r is None else r["obs_amp_m" if lab == "obs" else "model_amp_m"])
            g.append(np.nan if r is None else r["obs_g_deg" if lab == "obs" else "model_g_deg"])
        xpos = np.arange(len(sts)) + (i - len(runs) / 2) * w
        ax[0, j].bar(xpos, np.array(amp) * 100, w, label=lab, color="k" if lab == "obs" else None)
        ax[1, j].bar(xpos, g, w, color="k" if lab == "obs" else None)
    ax[0, j].set_title(n)
    ax[0, j].set_ylabel("amplitude (cm)")
    ax[1, j].set_ylabel("phase lag G (deg, UTC)")
    gs = [r["obs_g_deg"] for r in rows if r["constituent"] == n]
    ax[1, j].set_ylim(min(gs) - 25, max(gs) + 25)
    for k in range(2):
        ax[k, j].set_xticks(range(len(sts)), sts, rotation=60, fontsize=8)
ax[0, 0].legend(fontsize=8)
fig.suptitle(f"NAO.99Jb tide runs vs gauges, {t_lo:%Y-%m-%d} .. {t_hi:%Y-%m-%d}")
fig.tight_layout()
fig.savefig(a.out / "compare.png", dpi=130)
print(f"[450] wrote {a.out}", flush=True)
