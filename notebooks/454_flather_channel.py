"""An idealised channel for checking FVCOM open boundaries against theory.

    python notebooks/454_flather_channel.py stage   --root DIR --mode clamped|flather|flather0
    python notebooks/454_flather_channel.py analyse --root DIR [--root DIR ...] --out PNG

A flat-bottomed channel, 200 km long and 20 km wide, 100 m deep, open at
x = 0 and closed at x = L, no rotation, almost no friction, forced by M2 of
amplitude A = 0.5 m at the open boundary. Linear theory gives the standing
wave

    eta(x, t) = A cos(k (L - x)) / cos(k L) cos(w t),
    u_n(0, t) = A sqrt(g/H) tan(k L) sin(w t)   (outward at x = 0),

k = w / sqrt(g H). Modes:

- ``clamped``: eta = eta_T on the boundary (FVCOM's usual condition);
- ``flather``: Flather with the theoretical u_n (must reproduce the same);
- ``flather0``: Flather with u_T = 0 -- an inconsistent velocity on purpose,
  which a sound Flather boundary must survive (it then lets part of the
  wave leave);
- ``eq_given`` / ``eq_legacy``: clamped, plus the M2 tidal potential (FVCOM
  -DEQUI_TIDE), with f and V0+u from utide in the file, or without them so
  that FVCOM uses its own monthly astronomy -- the path that crashed before
  FVCOM 2026-10-01. The two must agree.

``analyse`` fits M2 along the channel's centre line over the last 3 days
and compares amplitude and phase with the theory.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fvcom_mesh_tools.io.fort14 import Fort14Mesh  # noqa: E402
from fvcom_mesh_tools.io.fvcom_native import export_fvcom_case  # noqa: E402
from fvcom_mesh_tools.tide_models import astronomy, spectral_text  # noqa: E402

G = 9.81
L, W, H, DX = 200e3, 20e3, 100.0, 2e3
A = 0.5
PERIOD = 44714.1643133778
X0, Y0 = 400e3, 3900e3          # anywhere in UTM zone 54; only the geometry matters
DAYS = 5.0


def theory(x):
    w = 2 * np.pi / PERIOD
    k = w / np.sqrt(G * H)
    return A * np.abs(np.cos(k * (L - x)) / np.cos(k * L)), k


def _load(name, file):
    spec = importlib.util.spec_from_file_location(name, ROOT / "notebooks" / file)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def build_mesh():
    nx, ny = int(L / DX) + 1, int(W / DX) + 1
    xs, ys = np.meshgrid(np.arange(nx) * DX, np.arange(ny) * DX, indexing="ij")
    nodes = np.column_stack([xs.ravel() + X0, ys.ravel() + Y0])
    idx = np.arange(nx * ny).reshape(nx, ny)
    tri = []
    for i in range(nx - 1):
        for j in range(ny - 1):
            a, b, c, d = idx[i, j], idx[i + 1, j], idx[i + 1, j + 1], idx[i, j + 1]
            # alternate the diagonal so the mesh has no preferred direction
            tri += ([[a, b, c], [a, c, d]] if (i + j) % 2 == 0 else [[a, b, d], [b, c, d]])
    tri = np.asarray(tri)
    # counter-clockwise
    p = nodes[tri]
    area = (p[:, 1, 0] - p[:, 0, 0]) * (p[:, 2, 1] - p[:, 0, 1]) \
        - (p[:, 2, 0] - p[:, 0, 0]) * (p[:, 1, 1] - p[:, 0, 1])
    tri[area < 0] = tri[area < 0][:, [0, 2, 1]]
    obc = idx[0, :]                       # x = 0
    return Fort14Mesh("channel", nodes, np.full(len(nodes), H), tri, [obc], []), obc


def stage(a):
    M383 = _load("m2_383", "383_m2_case_prep.py")
    mesh, obc = build_mesh()
    root = a.root.resolve()
    if (root / "STAGED").exists():
        raise SystemExit(f"{root} already holds a run")
    inp, out = root / "input", root / "output"
    for d in (inp, out):
        if len(f"{d}/") > 80:
            raise SystemExit(f"{d}/ exceeds FVCOM's 80 characters")
        d.mkdir(parents=True, exist_ok=True)
    export_fvcom_case(mesh, inp, "m2", cor=np.zeros(mesh.n_nodes), twodm=False,
                      obc_depth_control=False)
    (inp / "sigma.dat").write_text("NUMBER OF SIGMA LEVELS = 6\nSIGMA COORDINATE TYPE = UNIFORM\n")
    (inp / "m2_spg.dat").write_text("Sponge Node Number = 0\n")
    n = len(obc)
    amp = np.full((1, n), A)
    pha = np.zeros((1, n))                # eta_T = A cos(w t)
    normal = None
    if a.mode != "clamped":
        un = A * np.sqrt(G / H) * np.tan(theory(0.0)[1] * L) if a.mode == "flather" else 0.0
        # A sqrt(g/H) tan(kL) sin(w t) = ... cos(w t - 90 deg)
        normal = (np.full((1, n), abs(un)), np.full((1, n), 90.0 if un >= 0 else 270.0))
    M383.START = "2021-01-01 00:00:00"
    M383.END = f"2021-01-{1 + int(DAYS):02d} 00:00:00"
    M383.NC_OUT_INTERVAL_SECONDS = 600.0
    M383.DTE = 10.0
    M383.RAMP_SECONDS = 43200.0
    tide = spectral_text(["M2"], [PERIOD], amp, pha, M383.START, normal_velocity=normal)
    # an FVCOM built with -DEQUI_TIDE needs the equilibrium items on the
    # component line: a zero equilibrium tide (no potential in the channel),
    # with f and V0+u given -- without them FVCOM takes its own monthly
    # astronomy, which crashed (SIGSEGV in ELEVATION_EQUI, 2026-09-30)
    if a.mode.startswith("eq_"):
        from datetime import datetime

        from pyproj import Transformer
        lat = Transformer.from_crs(32654, 4326, always_xy=True).transform(X0, Y0)[1]
        _, f, v0u = astronomy(["M2"], datetime(2021, 1, 1), datetime(2021, 1, 3, 12), lat)
        extra = " 0.242334 0.693 SEMIDIURNAL"
        if a.mode == "eq_given":
            extra += f" {f[0]:.8f} {v0u[0]:.6f}"
    else:
        extra = " 0.0 0.693 SEMIDIURNAL 1.0 0.0"
    tide = tide.replace(f"1 = M2 {PERIOD:.10f}", f"1 = M2 {PERIOD:.10f}{extra}")
    (inp / "m2_tide.dat").write_text(tide)
    text = M383.namelist(inp, out)
    for key, value in (("CASE_TITLE", f"'454 channel {a.mode}'"), ("NC_VELOCITY", "F"),
                       ("NC_SALT_TEMP", "F"), ("NC_VERTICAL_VEL", "F"),
                       ("BOTTOM_ROUGHNESS_MINIMUM", "0.0001"),
                       ("BOTTOM_ROUGHNESS_LENGTHSCALE", "0.00001"),
                       ("SPONGE_FILE", "'m2_spg.dat'")):
        text, c = re.subn(rf"(?im)^(\s*{key}\s*=)[^\n]*", rf"\g<1> {value},", text)
        if c != 1:
            raise SystemExit(f"namelist key {key}: found {c}")
    (root / "m2_run.nml").write_text(text)
    (root / "manifest.json").write_text(json.dumps({
        "mode": a.mode, "L_m": L, "W_m": W, "H_m": H, "dx_m": DX, "A_m": A,
        "period_s": PERIOD, "kL_rad": float(theory(0.0)[1] * L),
        "u_T_amp_m_s": None if normal is None else float(normal[0][0, 0]),
        "n_nodes": mesh.n_nodes, "n_elements": mesh.n_elements, "days": DAYS,
        "start": M383.START, "end": M383.END, "dte_seconds": M383.DTE,
    }, indent=1) + "\n")
    (root / "STAGED").write_text("ok\n")
    print(f"[454] staged {a.mode} in {root}", flush=True)


def analyse(a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import netCDF4

    fig, ax = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    xs = np.linspace(0, L, 200)
    th, _ = theory(xs)
    ax[0].plot(xs / 1e3, th, "k", lw=2, label="theory")
    ax[1].axhline(0, color="k", lw=2)
    for root in a.root:
        man = json.loads((root / "manifest.json").read_text())
        mesh, _ = build_mesh()
        with netCDF4.Dataset(root / "output/m2_0001.nc") as ds:
            t = np.asarray(ds["time"][:], float) * 86400
            z = np.asarray(ds["zeta"][:], float)
        t = t - t[0]
        keep = t >= (DAYS - 3) * 86400
        w = 2 * np.pi / PERIOD
        X = np.column_stack([np.ones(keep.sum()), np.cos(w * t[keep]), np.sin(w * t[keep])])
        c, *_ = np.linalg.lstsq(X, z[keep], rcond=None)
        amp = np.hypot(c[1], c[2])
        pha = np.degrees(np.arctan2(c[2], c[1]))          # A cos(w t - pha)
        x = mesh.nodes[:, 0] - X0
        mid = np.isclose(mesh.nodes[:, 1] - Y0, W / 2)
        o = np.argsort(x[mid])
        xm, am, pm = x[mid][o], amp[mid][o], pha[mid][o]
        tm, _ = theory(xm)
        err = am - tm
        # theory: phase 0 everywhere (standing wave, no node inside for kL < pi/2)
        print(f"[454] {man['mode']}: amplitude error max {np.abs(err).max() * 100:.2f} cm "
              f"(rms {np.sqrt(np.mean(err ** 2)) * 100:.2f} cm), phase {pm.min():+.2f} .. "
              f"{pm.max():+.2f} deg, closed-end amp {am[-1]:.4f} m (theory {tm[-1]:.4f}), "
              f"max |zeta| {np.abs(z).max():.3f} m, mean {c[0].mean() * 100:+.2f} cm", flush=True)
        ax[0].plot(xm / 1e3, am, "o-", ms=3, label=man["mode"])
        ax[1].plot(xm / 1e3, (pm + 180) % 360 - 180, "o-", ms=3, label=man["mode"])
    ax[0].set_ylabel("M2 amplitude (m)")
    ax[1].set_ylabel("M2 phase (deg)")
    ax[1].set_xlabel("distance from the open boundary (km)")
    ax[0].legend()
    ax[0].set_title("Channel, 100 m deep, closed at 200 km: FVCOM against linear theory")
    fig.tight_layout()
    fig.savefig(a.out, dpi=120)
    print(f"[454] wrote {a.out}", flush=True)


p = argparse.ArgumentParser(description=__doc__)
sub = p.add_subparsers(dest="cmd", required=True)
s = sub.add_parser("stage")
s.add_argument("--root", type=Path, required=True)
s.add_argument("--mode", choices=("clamped", "flather", "flather0", "eq_given", "eq_legacy"),
               required=True)
n = sub.add_parser("analyse")
n.add_argument("--root", type=Path, action="append", required=True)
n.add_argument("--out", type=Path, required=True)
a = p.parse_args()
stage(a) if a.cmd == "stage" else analyse(a)
