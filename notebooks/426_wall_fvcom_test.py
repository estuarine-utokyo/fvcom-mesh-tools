# Does a split wall stop water in FVCOM, and does an unsplit one not?
#
# A 6 km x 2 km channel, 100 m right triangles, 10 m deep, an M2 tide of
# 0.5 m at the west end.  Three structures, all laid along grid lines so that
# the edges exist in BOTH meshes:
#
#   pier        x = 3000 m, from the south shore to mid-channel  (root + tip)
#   breakwater  x = 4500 m, detached, y = 600..1400 m            (two tips)
#   seal        x = 5500 m, shore to shore                       (closes a basin)
#
# The seal is the decisive part.  Behind a wall that is really boundary, the
# basin east of x = 5500 m has no open boundary and no connection, so its
# water level must stay at ZERO while the tide runs outside.  The control
# case is the same grid with the same edges and no split: there the claim in
# docs/linear_structures_design.md §2 says the tide goes straight through.
#
#   python notebooks/426_wall_fvcom_test.py prep <root>
#   python notebooks/426_wall_fvcom_test.py analyze <root>
import importlib.util
import json
import re
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from fvcom_mesh_tools.io.fort14 import Fort14Mesh  # noqa: E402
from fvcom_mesh_tools.io.fvcom_native import export_fvcom_case  # noqa: E402
from fvcom_mesh_tools.walls import split_along_walls, wall_edges_from_path  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "m2_383", ROOT / "notebooks" / "383_m2_case_prep.py")
m383 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m383)

H = 100.0
NX, NY = 61, 21
DEPTH = 10.0
AMP = 0.5
DAYS = 3.0
DTE, ISPLIT = 2.0, 10
SEAL_X = 5500.0
WALLS = {
    "pier": [(30, j) for j in range(0, 11)],
    "breakwater": [(45, j) for j in range(6, 15)],
    "seal": [(55, j) for j in range(0, NY)],
}


def node(i, j):
    return j * NX + i


def grid():
    gx, gy = np.meshgrid(np.arange(NX) * H, np.arange(NY) * H)
    xy = np.column_stack([gx.ravel(), gy.ravel()])
    tri = []
    for j in range(NY - 1):
        for i in range(NX - 1):
            a = node(i, j)
            tri += [[a, a + 1, a + NX + 1], [a, a + NX + 1, a + NX]]
    return xy, np.asarray(tri, dtype=np.int64)


def prep(root: Path):
    xy, tri = grid()
    walls = np.vstack([wall_edges_from_path([node(i, j) for i, j in path])
                       for path in WALLS.values()])
    sxy, stri, copy_of, rep = split_along_walls(xy, tri, walls)
    print(f"[426] split: {json.dumps({k: v for k, v in rep.items() if k != 'pairs'})}")
    obc = np.array([node(0, j) for j in range(NY)])
    period, _ = m383.tide_constants()
    start = datetime.fromisoformat(m383.START)
    end = (start + timedelta(days=DAYS)).strftime("%Y-%m-%d %H:%M:%S")
    for case, (p, t) in {"split": (sxy, stri), "edges": (xy, tri)}.items():
        mesh = Fort14Mesh(title=f"426 {case}", nodes=p, depths=np.full(len(p), DEPTH),
                          elements=t, open_boundaries=[obc], land_boundaries=[])
        inp, out = root / case / "input", root / case / "output"
        inp.mkdir(parents=True, exist_ok=True)
        out.mkdir(exist_ok=True)
        export_fvcom_case(mesh, inp, "m2", cor=np.full(len(p), 35.0),
                          write_empty_spg=True, twodm=False, obc_depth_control=False)
        (inp / "sigma.dat").write_text(
            "NUMBER OF SIGMA LEVELS = 6\nSIGMA COORDINATE TYPE = UNIFORM\n")
        (inp / "m2_tide.dat").write_text(
            m383.spectral_text(period, np.full(len(obc), AMP), np.zeros(len(obc))))
        nml = m383.namelist(inp, out)
        for key, value in {"END_DATE": f"'{end}'", "EXTSTEP_SECONDS": str(DTE),
                           "ISPLIT": str(ISPLIT),
                           "IRAMP": str(round(86400.0 / (DTE * ISPLIT))),
                           "CASE_TITLE": f"'426 wall test, {case}'"}.items():
            nml, n = re.subn(rf"(?im)^(\s*{key}\s*=)[^\n]*", rf"\g<1> {value},", nml)
            if n != 1:
                raise ValueError(f"namelist key {key}: found {n}")
        (root / case / "m2_run.nml").write_text(nml)
        print(f"[426] {case}: NP={len(p):,} NE={len(t):,} -> {root / case}")
    (root / "walls.json").write_text(json.dumps(
        {"walls": WALLS, "split": {k: v for k, v in rep.items() if k != "pairs"},
         "period_s": period, "end": end}, indent=1))


def m2_fit(t, z, period):
    w = 2 * np.pi / period
    a = np.column_stack([np.ones_like(t), np.cos(w * t), np.sin(w * t)])
    coef, *_ = np.linalg.lstsq(a, z, rcond=None)
    return np.hypot(coef[1], coef[2])


def analyze(root: Path):
    import netCDF4

    meta = json.loads((root / "walls.json").read_text())
    period = meta["period_s"]
    res = {}
    for case in ("split", "edges"):
        with netCDF4.Dataset(root / case / "output" / "m2_0001.nc") as ds:
            t = np.asarray(ds["time"][:], float) * 86400.0
            z = np.asarray(ds["zeta"][:], float)
            x = np.asarray(ds["x"][:], float)
            u = np.asarray(ds["ua"][:], float)
            v = np.asarray(ds["va"][:], float)
            nv = np.asarray(ds["nv"][:], int).T - 1
        keep = t >= t[0] + 86400.0            # after the one-day ramp
        amp = m2_fit(t[keep], z[keep], period)
        basin = x > SEAL_X + 1.0
        outside = (x > SEAL_X - 600.0) & (x < SEAL_X - 1.0)
        speed = np.hypot(u, v)[keep]
        r = {
            "finite": bool(np.isfinite(z).all() and np.isfinite(speed).all()),
            "zeta_max_abs_m": float(np.abs(z).max()),
            "basin_zeta_max_abs_m": float(np.abs(z[:, basin]).max()),
            "basin_m2_amp_m": float(amp[basin].max()),
            "outside_m2_amp_m": float(np.median(amp[outside])),
            "speed_p99_ms": float(np.percentile(speed, 99)),
            "speed_max_ms": float(speed.max()),
        }
        if case == "split":
            tips = [node(30, 10), node(45, 6), node(45, 14)]
            at_tip = np.isin(nv, tips).any(axis=1)
            r["tip_speed_max_ms"] = float(speed[:, at_tip].max())
        res[case] = r
        print(f"[426] {case}: {json.dumps(r)}")
    s, e = res["split"], res["edges"]
    verdict = {
        "split_basin_is_still": s["basin_m2_amp_m"] < 1e-4 and s["basin_zeta_max_abs_m"] < 1e-3,
        "unsplit_basin_sees_the_tide": e["basin_m2_amp_m"] > 0.5 * e["outside_m2_amp_m"],
        "split_run_is_finite": s["finite"],
        "tips_are_not_a_hotspot": s.get("tip_speed_max_ms", 0) <= 3 * s["speed_p99_ms"],
    }
    print(f"[426] VERDICT {json.dumps(verdict)}")
    (root / "result.json").write_text(json.dumps({"cases": res, "verdict": verdict}, indent=1))


if __name__ == "__main__":
    {"prep": prep, "analyze": analyze}[sys.argv[1]](Path(sys.argv[2]).resolve())
