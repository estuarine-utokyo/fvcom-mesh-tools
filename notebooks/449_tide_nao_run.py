"""Stage an astronomical-tide run of one case, a tide model on its open boundary.

    python notebooks/449_tide_nao_run.py --case outputs/extend_tokyo_bay_enshu/TokyoBayEnshu \\
        --root $WORK_DIR/scratch/tide_<stamp>/enshu [--days 200] [--start 2021-01-01]

A tide-only hindcast: no wind, no rivers, no density. Each open-boundary
node gets the eight major constituents of the tide model (``--tide-model``:
NAO.99Jb, Matsumoto et al., 2000, by default; TPXO10-atlas-v2; or FES2022b) at its
position, converted to FVCOM's spectral form by
``tide_models.fvcom_spectral`` (V0 at the start, f and u frozen at the
middle of the run). The result is meant to be compared with tide gauges by a
harmonic analysis that uses the same astronomy (utide) on both series.

200 days by default: K1/P1 and S2/K2 need 183 days to be told apart
(Rayleigh), and the first 15 days are spin-up.

With ``--equilibrium`` the forcing file also carries the tidal potential
(equilibrium tide, ``tide_models.EQUILIBRIUM``, with the same f and V0+u as
the open boundary); it needs an FVCOM built with -DEQUI_TIDE and the owner's
extension that reads f and V0+u (FVCOM ``octopus/build_fvcom_equi.sh``).

Everything else -- namelist, sigma, the external step rule, the sponge -- is
448's, which is 383/414's.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import re
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
from pyproj import Transformer

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fvcom_mesh_tools.io.fvcom_native import (  # noqa: E402
    apply_obc_depth_control,
    export_fvcom_case,
    fvcom_next_obc,
    read_fvcom_case,
)
from fvcom_mesh_tools.tide_models import (  # noqa: E402
    astronomy,
    fvcom_spectral,
    load_fes,
    load_nao,
    load_tpxo,
    load_tpxo_transport,
    sample_constants,
    spectral_text,
)


def _load(name, file):
    spec = importlib.util.spec_from_file_location(name, ROOT / "notebooks" / file)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


M383 = _load("m2_383", "383_m2_case_prep.py")
M414 = _load("m2_414", "414_refine_m2_prep.py")
MESH_EPSG = 32654
CONSTITUENTS = ("M2", "S2", "N2", "K2", "K1", "O1", "P1", "Q1")
OUT_INTERVAL_S = 3600.0

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--case", type=Path, required=True, help="case prefix (…/<CASE>)")
p.add_argument("--root", type=Path, required=True)
p.add_argument("--start", default="2021-01-01")
p.add_argument("--days", type=float, default=200.0)
p.add_argument("--constituents", default=",".join(CONSTITUENTS),
               help="comma-separated (utide names, e.g. LDA2 for lambda2); default the major 8")
p.add_argument("--tide-model", choices=("nao99jb", "tpxo10", "fes2022"), default="nao99jb")
p.add_argument("--nao", type=Path, default=None,
               help="the tide model's directory (default under $DATA_DIR/tides/models)")
p.add_argument("--z0", type=float, default=None,
               help="bottom roughness length (m); default 383's 0.002693138")
p.add_argument("--cd-min", type=float, default=None,
               help="minimum bottom drag coefficient; default 383's 0.003")
p.add_argument("--flather", action="store_true",
               help="Flather boundary: also give the tide model's normal velocity "
                    "(TPXO transports / model depth); needs --tide-model=tpxo10")
p.add_argument("--flather-alpha", type=float, default=None,
               help="relaxation (0, 1] of the Flather correction term")
p.add_argument("--flather-min-depth", type=float, default=None,
               help="clamp the boundary nodes shallower than this (m); Flather elsewhere")
p.add_argument("--no-obc-depth-control", action="store_true",
               help="keep the open-boundary depths as built (FVCOM's OBC_DEPTH_CONTROL_ON = F)")
p.add_argument("--itd-coefficient", type=float, default=None,
               help="topographic (internal-tide) drag, 1/s: a linear bottom stress "
                    "C*H*|grad H|**2*u below --itd-min-depth (owner's FVCOM)")
p.add_argument("--itd-min-depth", type=float, default=200.0)
p.add_argument("--sal-beta", type=float, default=None,
               help="self-attraction and loading, scalar approximation (needs --equilibrium)")
p.add_argument("--equilibrium", action="store_true",
               help="add the tidal potential (needs an FVCOM built with -DEQUI_TIDE)")
a = p.parse_args()
CONSTITUENTS = tuple(c.strip().upper() for c in a.constituents.split(",") if c.strip())
if a.flather and a.tide_model != "tpxo10":
    raise SystemExit("--flather needs --tide-model=tpxo10 (the only model with transports here)")
DEFAULT_DIR = {"nao99jb": "tides/models/NAO.99Jb/ocean", "tpxo10": "tides/models/TPXO10_atlas_v2",
               "fes2022": "tides/models/FES2022b/ocean_tide_extrapolated"}
if a.nao is None:
    if not os.environ.get("DATA_DIR"):
        raise SystemExit("set DATA_DIR or pass --nao")
    a.nao = Path(os.environ["DATA_DIR"]) / DEFAULT_DIR[a.tide_model]

grd, dep, obc = (Path(f"{a.case}_{k}.dat") for k in ("grd", "dep", "obc"))
mesh = read_fvcom_case(grd, dep, obc, title="tide")
if len(mesh.open_boundaries) != 1:
    raise SystemExit("one open boundary is needed")
start = datetime.fromisoformat(a.start)
end = start + timedelta(days=a.days)
M383.START = start.strftime("%Y-%m-%d %H:%M:%S")
M383.END = end.strftime("%Y-%m-%d %H:%M:%S")
M383.NC_OUT_INTERVAL_SECONDS = OUT_INTERVAL_S
M383.DTE = M414.dividing_step(OUT_INTERVAL_S, M383.ISPLIT, M414.external_step(mesh))

if a.no_obc_depth_control:
    change = np.zeros(len(np.asarray(mesh.open_boundaries[0])))
else:
    mesh, change = apply_obc_depth_control(mesh)
case = a.root
# FVCOM keeps INPUT_DIR / OUTPUT_DIR in CHARACTER(LEN=80): a longer path is cut
# and the output lands beside the directory (a 81-character root, 2026-09-30)
for d in (case / "input", case / "output"):
    if len(f"{d.resolve()}/") > 80:
        raise SystemExit(f"{d.resolve()}/ exceeds FVCOM's 80 characters; use a shorter --root")
if (case / "STAGED").exists() or any((case / "output").glob("*.nc")):
    raise SystemExit(f"{case} already holds a run; give a new --root")
inp, out = case / "input", case / "output"
inp.mkdir(parents=True, exist_ok=True)
out.mkdir(exist_ok=True)
to_ll = Transformer.from_crs(MESH_EPSG, 4326, always_xy=True)
lon, lat = (np.asarray(v) for v in to_ll.transform(mesh.nodes[:, 0], mesh.nodes[:, 1]))
export_fvcom_case(mesh, inp, "m2", cor=lat, twodm=False, obc_depth_control=False)
(inp / "sigma.dat").write_text("NUMBER OF SIGMA LEVELS = 6\nSIGMA COORDINATE TYPE = UNIFORM\n")

# sponge: 448's (383's production values)
ob = np.asarray(mesh.open_boundaries[0])
xy = mesh.nodes[ob, :2]
e = np.linalg.norm(np.diff(xy, axis=0), axis=1)
radius = np.r_[e[0], 0.5 * (e[1:] + e[:-1]), e[-1]]
(inp / "m2_spg.dat").write_text(
    f"Sponge Node Number = {len(ob)}\n"
    + "".join(f"{n + 1} {r:.6f} {0.001:.6f}\n" for n, r in zip(ob, radius)))

# tide: the tide model at every open-boundary node
if a.tide_model == "nao99jb":
    grids = load_nao(a.nao, CONSTITUENTS)
else:
    window = (float(lon[ob].min()) - 0.5, float(lon[ob].max()) + 0.5,
              float(lat[ob].min()) - 0.5, float(lat[ob].max()) + 0.5)
    load = load_tpxo if a.tide_model == "tpxo10" else load_fes
    grids = load(a.nao, CONSTITUENTS, window=window)
amp = np.empty((len(CONSTITUENTS), len(ob)))
pha = np.empty_like(amp)
filled = np.zeros(len(ob), bool)
for k, c in enumerate(CONSTITUENTS):
    amp[k], pha[k], f = sample_constants(grids[c], lon[ob], lat[ob])
    filled |= f
if not (np.isfinite(amp).all() and np.isfinite(pha).all()):
    bad = np.flatnonzero(~np.isfinite(amp).all(axis=0))
    raise SystemExit(f"{a.tide_model} has no ocean near open-boundary node(s) {bad.tolist()}")
mid = start + (end - start) / 2
period, famp, fpha = fvcom_spectral(CONSTITUENTS, amp, pha, start, mid, float(lat[ob].mean()))
_, f_nodal, v0u = astronomy(CONSTITUENTS, start, mid, float(lat[ob].mean()))
normal = None
if a.flather:
    # outward normal as FVCOM takes it: from NEXT_OBC to the node, expressed
    # east/north through the two points' lon/lat (the UTM grid is rotated)
    nxt, _ = fvcom_next_obc(mesh.nodes[:, :2], mesh.elements, ob)
    lat0 = np.radians(lat[ob])
    de = (lon[ob] - lon[nxt]) * np.cos(lat0)
    dn = lat[ob] - lat[nxt]
    norm = np.hypot(de, dn)
    if (norm == 0).any():
        raise SystemExit("an open-boundary node coincides with its NEXT_OBC")
    ne, nn = de / norm, dn / norm
    tr = load_tpxo_transport(a.nao, CONSTITUENTS, window=window)
    h_obc = mesh.depths[ob]                     # the depth FVCOM runs with
    un_amp = np.empty_like(amp)
    un_pha = np.empty_like(amp)
    for k, c in enumerate(CONSTITUENTS):
        au, pu, _ = sample_constants(tr[c]["u"], lon[ob], lat[ob])
        av, pv, _ = sample_constants(tr[c]["v"], lon[ob], lat[ob])
        zn = (au * np.exp(-1j * np.radians(pu)) * ne
              + av * np.exp(-1j * np.radians(pv)) * nn) / h_obc
        if not np.isfinite(zn).all():
            raise SystemExit(f"no TPXO transport near open-boundary node(s) for {c}")
        un_amp[k], un_pha[k] = np.abs(zn), np.degrees(-np.angle(zn)) % 360
    _, un_famp, un_fpha = fvcom_spectral(CONSTITUENTS, un_amp, un_pha, start, mid,
                                         float(lat[ob].mean()))
    normal = (un_famp, un_fpha)
(inp / "m2_tide.dat").write_text(spectral_text(
    CONSTITUENTS, period, famp, fpha, M383.START,
    equilibrium=(f_nodal, v0u) if a.equilibrium else None, sal_beta=a.sal_beta,
    normal_velocity=normal, flather_alpha=a.flather_alpha,
    flather_min_depth=a.flather_min_depth, zero_unknown_equilibrium=True))
with open(case / "obc_constants.csv", "w") as fh:
    fh.write("obc,node,lon,lat,filled,"
             + ",".join(f"{c}_amp_m,{c}_g_deg" for c in CONSTITUENTS) + "\n")
    for j, n in enumerate(ob):
        fh.write(f"{j + 1},{n + 1},{lon[n]:.6f},{lat[n]:.6f},{int(filled[j])},"
                 + ",".join(f"{amp[k, j]:.5f},{pha[k, j]:.3f}" for k in range(len(CONSTITUENTS)))
                 + "\n")

# namelist: 383's, with tide-only output (zeta and the depth-mean velocity)
text = M383.namelist(inp, out)
changes = [("CASE_TITLE", f"'449 {a.tide_model} tide'"), ("NC_VELOCITY", "F"),
           ("NC_SALT_TEMP", "F"), ("NC_VERTICAL_VEL", "F")]
if a.z0 is not None:
    changes.append(("BOTTOM_ROUGHNESS_LENGTHSCALE", repr(a.z0)))
if a.cd_min is not None:
    changes.append(("BOTTOM_ROUGHNESS_MINIMUM", repr(a.cd_min)))
for key, value in changes:
    text, count = re.subn(rf"(?im)^(\s*{key}\s*=)[^\n]*", rf"\g<1> {value},", text)
    if count != 1:
        raise SystemExit(f"namelist key {key}: expected once, found {count}")
if a.itd_coefficient is not None:
    # the keys are not in the template (FVCOM's default is off): add them
    text, count = re.subn(
        r"(?im)^(\s*BOTTOM_ROUGHNESS_MINIMUM\s*=[^\n]*\n)",
        rf"\g<1> BOTTOM_ITD_COEFFICIENT = {a.itd_coefficient!r},\n"
        rf" BOTTOM_ITD_MIN_DEPTH = {a.itd_min_depth!r},\n", text)
    if count != 1:
        raise SystemExit(f"namelist key BOTTOM_ROUGHNESS_MINIMUM: expected once, found {count}")
if a.no_obc_depth_control:
    # the key is not in the template (FVCOM's default is T): add it to its block
    text, count = re.subn(r"(?im)^(\s*OBC_ON\s*=[^\n]*\n)", r"\g<1> OBC_DEPTH_CONTROL_ON = F,\n",
                          text)
    if count != 1:
        raise SystemExit(f"namelist key OBC_ON: expected once, found {count}")
(case / "m2_run.nml").write_text(text)

manifest = {
    "purpose": f"astronomical tide hindcast, {a.tide_model} on the open boundary",
    "case": str(a.case), "start": M383.START, "end": M383.END, "days": a.days,
    "nodal_f_u_at": mid.isoformat(), "constituents": list(CONSTITUENTS),
    "tide_model_name": a.tide_model, "tide_model": str(a.nao),
    "equilibrium_tide": bool(a.equilibrium),
    "bottom_z0_m": a.z0, "bottom_cd_min": a.cd_min,
    "obc_depth_control": not a.no_obc_depth_control, "sal_beta": a.sal_beta,
    "flather": bool(a.flather), "flather_alpha": a.flather_alpha,
    "flather_min_depth_m": a.flather_min_depth,
    "obc_normal_velocity_m2_max_m_s": float(normal[0][0].max()) if normal else None,
    "itd_coefficient_per_s": a.itd_coefficient, "itd_min_depth_m": a.itd_min_depth,
    "n_obc_filled_from_nearest": int(filled.sum()),
    "obc_amp_range_m": {c: [float(amp[k].min()), float(amp[k].max())]
                        for k, c in enumerate(CONSTITUENTS)},
    "dte_seconds": M383.DTE, "isplit": M383.ISPLIT,
    "max_dte_allowed_seconds": M414.external_step(mesh),
    "out_interval_seconds": OUT_INTERVAL_S, "n_obc_nodes": int(len(ob)),
    "sponge_radius_m": [float(radius.min()), float(radius.max())], "sponge_coef": 0.001,
    "obc_depth_control_change_m_max": float(np.max(np.abs(change))) if len(change) else 0.0,
    "n_nodes": mesh.n_nodes, "n_elements": mesh.n_elements,
    "depth_m": [float(mesh.depths.min()), float(mesh.depths.max())],
}
(case / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
(case / "STAGED").write_text(datetime.now().isoformat() + "\n")
print("[449] " + json.dumps(manifest), flush=True)
