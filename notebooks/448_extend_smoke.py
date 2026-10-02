"""Stage a short FVCOM run of one case: does the mesh integrate at all?

    python notebooks/448_extend_smoke.py --case outputs/extend_tokyo_bay_enshu/TokyoBayEnshu \\
        --root $WORK_DIR/scratch/smoke_<stamp> [--days 2] [--gauge MERA]

A STABILITY test for a mesh whose open boundary has no tidal forcing yet
(the wide Tokyo Bay mesh: its forcing will come from JCOPE-T DA). The whole
open boundary gets one M2 amplitude and phase, those of a published gauge
facing the open sea (MERA by default), so the model is driven by a tide of
the right size; the answer is not meant to be right, only finite.

Everything else -- namelist, sigma, tide file format, the external step
rule -- is 383/414's, imported. The sponge follows 383's production values:
coefficient 0.001 at every open-boundary node, radius = the node's own
boundary spacing (goto2023 uses about one spacing, 1.0-1.9 km on 1.7 km).
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
from pyproj import Transformer

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from fvcom_mesh_tools.io.fvcom_namelist import end_after  # noqa: E402
from fvcom_mesh_tools.io.fvcom_native import (  # noqa: E402
    apply_obc_depth_control,
    export_fvcom_case,
    read_fvcom_case,
)
from fvcom_mesh_tools.outdir import reserve  # noqa: E402


def _load(name, file):
    spec = importlib.util.spec_from_file_location(name, ROOT / "notebooks" / file)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


M383 = _load("m2_383", "383_m2_case_prep.py")
M414 = _load("m2_414", "414_refine_m2_prep.py")
MESH_EPSG = 32654

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--case", type=Path, required=True, help="case prefix (…/<CASE>)")
p.add_argument("--root", type=Path, required=True)
p.add_argument("--days", type=float, default=2.0)
p.add_argument("--gauge", default="MERA")
a = p.parse_args()
# absolute paths: the namelist is read from inside the case directory (review F24)
a.case, a.root = a.case.resolve(), a.root.resolve()
# a reused root could leave an old history beside new inputs (review F22);
# reserved atomically, before any work (review r2 F1)
# FVCOM would cut a longer run directory (round 4 F13): checked before the
# root is taken
M383.check_fvcom_dirs(a.root / "extended" / "input", a.root / "extended" / "output")
# a finite positive run, checked before the root is taken (review round 12 F5)
try:
    M383.END = end_after(M383.START, a.days)
except ValueError as err:
    raise SystemExit(str(err)) from None
a.root = reserve(a.root)

grd, dep, obc = (Path(f"{a.case}_{k}.dat") for k in ("grd", "dep", "obc"))
mesh = read_fvcom_case(grd, dep, obc, title="smoke")
if len(mesh.open_boundaries) != 1:
    raise SystemExit("one open boundary is needed")
# every depth change first, then the time step from the depths that run
# (review F23)
mesh, change = apply_obc_depth_control(mesh)
M383.DTE = M414.dividing_step(M383.NC_OUT_INTERVAL_SECONDS, M383.ISPLIT,
                              M414.external_step(mesh))
period, gauges = M383.tide_constants() if a.gauge in ("ABURATUBO", "MERA") else (None, None)
if gauges is None:
    raise SystemExit("--gauge must be ABURATUBO or MERA (the gauges 383 reads)")
g = gauges[a.gauge]

case = a.root / "extended"
inp, out = case / "input", case / "output"
inp.mkdir(parents=True, exist_ok=True)
out.mkdir(exist_ok=True)
_, lat = Transformer.from_crs(MESH_EPSG, 4326, always_xy=True).transform(
    mesh.nodes[:, 0], mesh.nodes[:, 1])
# the smoke test forces M2 on the open boundary: the boundary must be type 1
# (elevation), whatever the case declared; type 3 is a zero-elevation clamp
# that ignores the forcing (review round 18 F1). Read back to make sure.
written = export_fvcom_case(mesh, inp, "m2", cor=lat, twodm=False, obc_depth_control=False,
                            obc_type=1)
_types = {line.split()[2] for line in written["obc"].read_text().splitlines()[1:] if line.strip()}
if _types != {"1"}:
    raise SystemExit(f"the staged open boundary has types {_types}, not 1")
(inp / "sigma.dat").write_text("NUMBER OF SIGMA LEVELS = 6\nSIGMA COORDINATE TYPE = UNIFORM\n")
ob = np.asarray(mesh.open_boundaries[0])
xy = mesh.nodes[ob, :2]
e = np.linalg.norm(np.diff(xy, axis=0), axis=1)
radius = np.r_[e[0], 0.5 * (e[1:] + e[:-1]), e[-1]]
(inp / "m2_spg.dat").write_text(
    f"Sponge Node Number = {len(ob)}\n"
    + "".join(f"{n + 1} {r:.6f} {0.001:.6f}\n" for n, r in zip(ob, radius)))
amp = np.full(len(ob), g["amplitude_m"])
phase = np.full(len(ob), g["phase_deg"])
(inp / "m2_tide.dat").write_text(M383.spectral_text(period, amp, phase))
(case / "m2_run.nml").write_text(M383.namelist(inp, out))
manifest = {
    "purpose": "stability smoke test; uniform M2 on the open boundary, not a tidal hindcast",
    "case": str(a.case), "days": a.days, "start": M383.START, "end": M383.END,
    "dte_seconds": M383.DTE, "isplit": M383.ISPLIT,
    "max_dte_allowed_seconds": M414.external_step(mesh),
    "gauge": a.gauge, "amplitude_m": g["amplitude_m"], "phase_deg": g["phase_deg"],
    "period_seconds": period, "n_obc_nodes": int(len(ob)),
    "obc_type": 1, "case_obc_type": int(getattr(mesh, "obc_type", 1)),
    "sponge_radius_m": [float(radius.min()), float(radius.max())], "sponge_coef": 0.001,
    "obc_depth_control_change_m_max": float(np.max(np.abs(change))) if len(change) else 0.0,
    "n_nodes": mesh.n_nodes, "n_elements": mesh.n_elements,
    "depth_m": [float(mesh.depths.min()), float(mesh.depths.max())],
}
(a.root / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
(a.root / "STAGED").write_text(datetime.now().isoformat() + "\n")
print("[448] " + json.dumps(manifest), flush=True)
