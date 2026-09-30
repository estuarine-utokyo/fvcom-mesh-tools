"""Where does a tide run depart from the tide model that forces it? Maps, every node.

    python notebooks/451_tide_model_map.py --run $ROOT/enshu --out $ROOT/map_enshu \\
        [--spinup-days 15] [--show M2,S2,K1,O1]

For a run staged by 449, fit the run's constituents to ``zeta`` at every
node by least squares over the run after spin-up, with the frequencies,
V0 and the frozen f and u that 449 used for the forcing -- so the fit
returns Greenwich constants directly comparable with NAO.99Jb sampled at
the same node. On the open boundary the two agree by construction (a check);
inside, the ratio model/NAO shows where the model amplifies or damps a
constituent. NAO.99Jb assimilated coastal tide gauges but is 1/12 deg: in
the bay it is a reference, not the truth.

Writes ``map_<C>.png`` per shown constituent, ``nodes.npz`` (all
constants) and ``summary.json`` (statistics on the open boundary, on any
other open-boundary line given with ``--line``, and at the gauges).
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
from pyproj import Transformer

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fvcom_mesh_tools.io.fvcom_native import read_fvcom_case  # noqa: E402
from fvcom_mesh_tools.tide_models import (  # noqa: E402
    fvcom_spectral,
    load_fes,
    load_nao,
    load_tpxo,
    sample_constants,
)

MESH_EPSG = 32654
GAUGES = {"tokyo": (139.77, 35.66), "chiba": (140.05, 35.57), "yokohama": (139.65, 35.45),
          "yokosuka": (139.651, 35.288), "kurihama": (139.72, 35.23),
          "aburatsubo": (139.614, 35.159), "mera": (139.825, 34.919)}

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--run", type=Path, required=True)
p.add_argument("--out", type=Path, required=True)
p.add_argument("--spinup-days", type=float, default=15.0)
p.add_argument("--show", default="M2,S2,K1,O1")
p.add_argument("--reference", choices=("forcing", "nao99jb", "tpxo10", "fes2022"),
               default="forcing",
               help="the tide model to compare with (default: the one forcing the run)")
p.add_argument("--line", action="append", default=[],
               help="LABEL=obc.dat of another mesh: report the run along that line")
a = p.parse_args()
a.out.mkdir(parents=True, exist_ok=True)

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402
import netCDF4  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402

man = json.loads((a.run / "manifest.json").read_text())
names = man["constituents"]
start = datetime.fromisoformat(man["start"])
mid = datetime.fromisoformat(man["nodal_f_u_at"])
LABEL = {"nao99jb": "NAO.99Jb", "tpxo10": "TPXO10-atlas-v2", "fes2022": "FES2022b"}
DEFAULT_DIR = {"nao99jb": "tides/models/NAO.99Jb/ocean", "tpxo10": "tides/models/TPXO10_atlas_v2",
               "fes2022": "tides/models/FES2022b/ocean_tide_extrapolated"}
forcing = man.get("tide_model_name", "nao99jb")
model_name = forcing if a.reference == "forcing" else a.reference
if model_name == forcing:
    nao_dir = Path(man["tide_model"])
elif os.environ.get("DATA_DIR"):
    nao_dir = Path(os.environ["DATA_DIR"]) / DEFAULT_DIR[model_name]
else:
    raise SystemExit("set DATA_DIR to compare with another tide model")
if not nao_dir.exists():
    raise SystemExit(f"tide model directory {nao_dir} is missing")
REF = LABEL[model_name]

m = read_fvcom_case(a.run / "input/m2_grd.dat", a.run / "input/m2_dep.dat",
                    a.run / "input/m2_obc.dat", title="run")
to_ll = Transformer.from_crs(MESH_EPSG, 4326, always_xy=True)
lon, lat = (np.asarray(v) for v in to_ll.transform(m.nodes[:, 0], m.nodes[:, 1]))
with netCDF4.Dataset(a.run / "output/m2_0001.nc") as ds:
    days = np.asarray(ds["time"][:], float)            # MJD
    zeta = np.asarray(ds["zeta"][:], float)
t = (days - (start - datetime(1858, 11, 17)).total_seconds() / 86400) * 86400  # s since start
keep = t >= a.spinup_days * 86400
t, zeta = t[keep], zeta[keep]

# least squares with the forcing's own frequencies; FVCOM form A cos(w t - phi)
ob = np.asarray(m.open_boundaries[0])
one = np.ones((len(names), 1))
period, f_amp, f_pha = fvcom_spectral(names, one, 0 * one, start, mid, float(lat[ob].mean()))
f = f_amp[:, 0]
v0u = (-f_pha[:, 0]) % 360          # V0 + u, since phi = G - V0 - u with G = 0
w = 2 * np.pi / period
X = np.column_stack([np.ones_like(t)] + [np.cos(wk * t) for wk in w]
                    + [np.sin(wk * t) for wk in w])
coef, *_ = np.linalg.lstsq(X, zeta, rcond=None)
ca, sa = coef[1:1 + len(w)], coef[1 + len(w):]
amp_fv = np.hypot(ca, sa)
phi = np.degrees(np.arctan2(sa, ca)) % 360
A_mod = amp_fv / f[:, None]
G_mod = (phi + v0u[:, None]) % 360
resid = zeta - X @ coef

if model_name == "nao99jb":
    grids = load_nao(nao_dir, names)
else:
    load = load_tpxo if model_name == "tpxo10" else load_fes
    grids = load(nao_dir, names, window=(float(lon.min()) - 0.5, float(lon.max()) + 0.5,
                                        float(lat.min()) - 0.5, float(lat.max()) + 0.5))
A_nao = np.empty_like(A_mod)
G_nao = np.empty_like(G_mod)
for k, c in enumerate(names):
    A_nao[k], G_nao[k], _ = sample_constants(grids[c], lon, lat)
np.savez_compressed(a.out / "nodes.npz", names=np.array(names), lon=lon, lat=lat,
                    A_model=A_mod, G_model=G_mod, A_nao=A_nao, G_nao=G_nao,
                    resid_rms=np.sqrt((resid ** 2).mean(axis=0)))


def dphase(g1, g0):
    return (g1 - g0 + 180) % 360 - 180


def stats(idx):
    out = {}
    for k, c in enumerate(names):
        ok = np.isfinite(A_nao[k, idx])
        r = A_mod[k, idx][ok] / A_nao[k, idx][ok]
        dg = dphase(G_mod[k, idx][ok], G_nao[k, idx][ok])
        out[c] = {"model_amp_m": float(A_mod[k, idx].mean()),
                  "nao_amp_m": float(np.nanmean(A_nao[k, idx])),
                  "ratio_mean": float(r.mean()) if len(r) else None,
                  "ratio_range": [float(r.min()), float(r.max())] if len(r) else None,
                  "dg_mean_deg": float(dg.mean()) if len(dg) else None,
                  "n": int(ok.sum())}
    return out


tri = mtri.Triangulation(m.nodes[:, 0], m.nodes[:, 1], m.elements)
finder = tri.get_trifinder()
to_utm = Transformer.from_crs(4326, MESH_EPSG, always_xy=True)
summary = {"run": str(a.run), "window_days": float((t[-1] - t[0]) / 86400),
           "resid_rms_m_max": float(np.sqrt((resid ** 2).mean(axis=0)).max()),
           "open_boundary": stats(ob), "lines": {}, "gauges": {}}
for item in a.line:
    label, _, path = item.partition("=")
    xy = np.loadtxt(path, skiprows=1, dtype=int)[:, 1] - 1
    other = read_fvcom_case(Path(path).with_name(Path(path).name.replace("_obc", "_grd")),
                            Path(path).with_name(Path(path).name.replace("_obc", "_dep")),
                            path, title=label)
    pts = other.nodes[xy, :2]
    nearest = np.array([int(np.argmin(np.hypot(*(m.nodes[:, :2] - q).T))) for q in pts])
    summary["lines"][label] = stats(np.unique(nearest))
for g, (glon, glat) in GAUGES.items():
    x, y = to_utm.transform(glon, glat)
    k = int(np.argmin(np.hypot(m.nodes[:, 0] - x, m.nodes[:, 1] - y)))
    summary["gauges"][g] = {
        c: {"model_amp_m": float(A_mod[i, k]), "model_g_deg": float(G_mod[i, k]),
            "nao_amp_m": float(A_nao[i, k]), "nao_g_deg": float(G_nao[i, k])}
        for i, c in enumerate(names)}
(a.out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
ob_short = {c: (round(v["ratio_mean"], 3), round(v["dg_mean_deg"], 2))
            for c, v in summary["open_boundary"].items()}
print("[451] open boundary: " + json.dumps(ob_short), flush=True)
for label, s in summary["lines"].items():
    short = {c: (round(v["model_amp_m"], 4), round(v["nao_amp_m"], 4),
                 round(v["ratio_mean"], 3), round(v["dg_mean_deg"], 2)) for c, v in s.items()}
    print(f"[451] line {label}: " + json.dumps(short), flush=True)
for g, s in summary["gauges"].items():
    print(f"[451] {g}: " + ", ".join(f"{c} model {v['model_amp_m']:.3f}/{v['model_g_deg']:.1f} "
                                     f"{REF} {v['nao_amp_m']:.3f}/{v['nao_g_deg']:.1f}"
                                     for c, v in s.items() if c in ("M2", "K1", "O1")), flush=True)

# maps: coast (land boundary) black, open boundary red
edges = np.sort(np.vstack([m.elements[:, [0, 1]], m.elements[:, [1, 2]], m.elements[:, [2, 0]]]),
                axis=1)
u, cnt = np.unique(edges, axis=0, return_counts=True)
bnd = u[cnt == 1]
obset = set(ob.tolist())
is_open = np.array([a_ in obset and b_ in obset for a_, b_ in bnd])
xk, yk = m.nodes[:, 0] / 1e3, m.nodes[:, 1] / 1e3
trik = mtri.Triangulation(xk, yk, m.elements)
for c in a.show.split(","):
    k = names.index(c)
    fig, ax = plt.subplots(1, 3, figsize=(17, 5.6), constrained_layout=True)
    panels = ((A_mod[k] * 100, "viridis", None, f"{c} model amplitude (cm)"),
              (A_mod[k] / A_nao[k], "RdBu_r", (0.8, 1.2), f"{c} amplitude model / {REF}"),
              (dphase(G_mod[k], G_nao[k]), "RdBu_r", (-10, 10),
               f"{c} phase model - {REF} (deg)"))
    for axi, (val, cmap, lim, title) in zip(ax, panels):
        v = np.where(np.isfinite(val), val, np.nan)
        kw = {"vmin": lim[0], "vmax": lim[1]} if lim else {}
        pc = axi.tripcolor(trik, np.nan_to_num(v, nan=np.nanmean(v)), cmap=cmap,
                           shading="gouraud", **kw)
        segs = np.stack([np.c_[xk[bnd[:, 0]], yk[bnd[:, 0]]], np.c_[xk[bnd[:, 1]], yk[bnd[:, 1]]]],
                        axis=1)
        axi.add_collection(LineCollection(segs[~is_open], colors="k", linewidths=0.5))
        axi.add_collection(LineCollection(segs[is_open], colors="r", linewidths=1.0))
        axi.set_aspect("equal")
        axi.set_title(title, fontsize=10)
        axi.set_xlabel("UTM54 x (km)")
        fig.colorbar(pc, ax=axi, shrink=0.8)
    fig.savefig(a.out / f"map_{c}.png", dpi=120)
    plt.close(fig)
print(f"[451] wrote {a.out}", flush=True)
