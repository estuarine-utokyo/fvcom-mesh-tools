"""Stage the 383 M2 experiment; default scratch writes occur only when submitted.

Use --root <repository-local-directory> for a lightweight staging check.
ASCII spectral format: FVCOM mod_input.F READ_JULIAN_OBC/Parse_tide;
mod_force.F TIDAL_ELEVATION; mod_obcs.F uses A*cos(omega*t-phase).
No source inputs are modified. Repeated preparation preserves model output.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
from pyproj import Transformer

sys.dont_write_bytecode = True  # In particular, never write into the shared tide product.
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from fvcom_mesh_tools.dem.m7001 import production_depths  # noqa: E402
from fvcom_mesh_tools.io.fort14 import Fort14Mesh, read_fort14  # noqa: E402
from fvcom_mesh_tools.io.fvcom_native import apply_obc_depth_control, write_dep  # noqa: E402

BASE = Path("/octfs/work/G16445/v61021").resolve()
TB = BASE / "Github/TB-FVCOM/input/goto2023"
FV = BASE / "Github/FVCOM"
TIDES = Path("/octfs/work/G16445/share/Data/tides").resolve()
DEFAULT_RUN_ROOT = BASE / "scratch/m2_383"
DEP = TB / "grid/TokyoBay_dep_m7001tp_rfac0p2_cap300.dat"
DEFAULT_B_MESH = ROOT / "outputs/sample_repro/sample_repro_final.14"
# The mesh under test.  Overridable so a candidate can be validated without
# first promoting it into the repository (--mesh).
B_MESH = DEFAULT_B_MESH
TEMPLATE = FV / "Tests/PrecipEvap/run/out_pe_rain/pe_rain_run.nml"
# The bay's M2 response approaches its steady state with an e-folding time of
# about 3.8 days (measured on job 115305: the amplitude at HARUMI was still
# climbing 0.406 -> 0.447 m between day 6 and day 10).  The earlier 11-day run
# therefore reported a harmonic fit over a GROWING transient, not an M2
# constant.  20 days leaves a residual of exp(-20/3.8) = 0.5 %, and the
# analysis fits the last 5 days.
START = "2021-01-01 00:00:00"
END = "2021-01-21 00:00:00"
DTE = 5.0  # TB-FVCOM production EXTSTEP_SECONDS; implied dt >= 15.4 s for all meshes
ISPLIT = 10
# The output interval, in SECONDS, declared rather than written into the
# namelist by hand: FVCOM refuses a run whose NC_OUT_INTERVAL is not a whole
# number of internal steps, so whoever chooses the step has to know it.
NC_OUT_INTERVAL_SECONDS = 1800.0
# The tanh ramp, in SECONDS. FVCOM's IRAMP counts INTERNAL steps
# (DTI = DTE * ISPLIT), so a hard-coded IRAMP means the ramp changes whenever
# the external step does: 8640 was 5 days at DTE = 5 and would be 1.5 days at
# DTE = 1.5, while the manifest said 86,400 s either way (fourth review).
RAMP_SECONDS = 86400.0
SPINUP = 15 * 86400
DURATION = 20 * 86400
STATIONS = ("TOKYO-SIBAURA", "HARUMI", "TIBA-TIBA LIGHT", "SINKO", "YOKOSUKA")


def load_a():
    lines = (TB / "grid/TokyoBay_grd.dat").read_text().splitlines()
    nn, ne = [int(line.split("=")[1]) for line in lines[:2]]
    tri = np.array([s.split()[1:4] for s in lines[2 : 2 + ne]], int) - 1
    xy = np.array([s.split()[1:3] for s in lines[2 + ne : 2 + ne + nn]], float)
    dep = np.loadtxt(DEP, skiprows=1)
    if not np.allclose(dep[:, :2], xy, atol=0.02, rtol=0):
        raise ValueError("A depth coordinates do not match its grid")
    obc = np.loadtxt(TB / "grid/TokyoBay_obc.dat", skiprows=1, dtype=int)[:, 1] - 1
    return Fort14Mesh("goto2023 production", xy, dep[:, 2], tri, [obc], [])


def arc_position(points, arc):
    """Nearest projection onto an ordered polyline; return distance and offset, metres."""
    delta = np.diff(arc, axis=0)
    length = np.linalg.norm(delta, axis=1)
    if np.any(length == 0):
        raise ValueError("Repeated arc vertices")
    fraction = np.clip(np.einsum("nsi,si->ns", points[:, None] - arc[:-1], delta) / length**2, 0, 1)
    closest = arc[:-1] + fraction[..., None] * delta
    distance = np.linalg.norm(points[:, None] - closest, axis=2)
    idx = distance.argmin(axis=1)
    row = np.arange(len(points))
    return np.r_[0, np.cumsum(length)][idx] + fraction[row, idx] * length[idx], distance[row, idx]


def tide_constants():
    """Freeze the product's nodal factors/astronomy at Jan 1, 2021 UTC.

    Reproduce ptide.predict's January 1 formula, evaluated at 09:00 JST.
    The +omega*9 cancels its -omega*Ztm (all stations are in Japan).
    This is an idealized single-frequency experiment, not a daily nodal prediction.
    """
    source = TIDES / "ptide/ptide.py"
    spec = importlib.util.spec_from_file_location("ptide_383", source)
    pt = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pt)
    constants = TIDES / "ptide/Tide-ToKYOWAN.txt"
    iy, dl = 21, int((2021 + 3) / 4) - 500
    s = (211.728 + 129.38471 * iy + 13.176396 * dl) % 360
    h = (279.974 - 0.23871 * iy + 0.985647 * dl) % 360
    p = (83.298 + 40.66229 * iy + 0.111404 * dl) % 360
    n = (125.071 - 19.32812 * iy - 0.052954 * dl) % 360
    f, u = pt.nodal_factors(s, h, p, n)
    k = 30  # ptide constituent ordering: M2
    speed = (15 * pt.A1 + 0.54901652 * pt.A2 + 0.04106864 * pt.A3 + 0.00464181 * pt.A4)[k]
    result = {}
    for name in ("ABURATUBO", "MERA", *STATIONS):
        amp, lag, lon, lat, _ = pt.read_constants(constants, name)
        argument = (pt.A2 * s + pt.A3 * h + pt.A4 * p + pt.C + u)[k] + pt.A1[k] * lon
        result[name] = dict(
            lon=lon,
            lat=lat,
            amplitude_m=float(f[k] * amp[k] / 100),
            phase_deg=float((lag[k] - argument) % 360),
            published_amplitude_m=float(amp[k] / 100),
            published_kappa_deg=float(lag[k]),
            nodal_factor=float(f[k]),
        )
    return float(360 * 3600 / speed), result


def spectral_text(period, amp, phase):
    lines = [
        "Tidal Component Number = 1",
        f"1 = M2 {period:.10f}",
        f"Time Origin = {START}",
        f"OBC Node Number = {len(amp)}",
    ]
    # The first column is OBC ordinal, NOT global mesh node ID.
    for label, values in [("Amplitude", amp), ("Phase", phase), ("Eref", np.zeros(len(amp)))]:
        lines += [label] + [f"{i} {v:.10f}" for i, v in enumerate(values, 1)] + [label]
    return "\n".join(lines) + "\n"


def namelist(input_dir, output_dir):
    text = TEMPLATE.read_text()
    text = text[text.index("&NML_CASE") :]
    changes = dict(
        CASE_TITLE="'383 M2 mesh comparison'",
        END_DATE=f"'{END}'",
        INPUT_DIR=f"'{input_dir}/'",
        OUTPUT_DIR=f"'{output_dir}/'",
        EXTSTEP_SECONDS=str(DTE),
        ISPLIT=str(ISPLIT),
        IRAMP=str(max(1, round(RAMP_SECONDS / (DTE * ISPLIT)))),
        MIN_DEPTH="0.1",
        IREPORT="360",
        NC_OUT_INTERVAL=f"'seconds = {NC_OUT_INTERVAL_SECONDS:.1f}'",
        NC_EVAP_PRECIP="F",
        PRECIPITATION_ON="F",
        PRECIPITATION_PRC="0.0",
        OBC_ON="T",
        OBC_NODE_LIST_FILE="'m2_obc.dat'",
        OBC_ELEVATION_FORCING_ON="T",
        OBC_ELEVATION_FILE="'m2_tide.dat'",
        GRID_FILE="'m2_grd.dat'",
        DEPTH_FILE="'m2_dep.dat'",
        CORIOLIS_FILE="'m2_cor.dat'",
        SPONGE_FILE="'m2_spg.dat'",
        SIGMA_LEVELS_FILE="'sigma.dat'",
        TEMPERATURE_ACTIVE="F",
        SALINITY_ACTIVE="F",
        BOTTOM_ROUGHNESS_LENGTHSCALE="0.002693138",
        BOTTOM_ROUGHNESS_MINIMUM="0.003",
    )
    for key, value in changes.items():
        text, count = re.subn(rf"(?im)^(\s*{key}\s*=)[^\n]*", rf"\g<1> {value},", text)
        if count != 1:
            raise ValueError(f"Template key {key}: expected once, found {count}")
    return text


def mesh_metrics(mesh):
    xy = mesh.nodes[mesh.elements]
    edge = np.linalg.norm(xy - np.roll(xy, 1, axis=1), axis=2).min(axis=1)
    dt = edge / np.sqrt(9.81 * mesh.depths[mesh.elements].max(axis=1))
    if np.any(mesh.depths <= 0) or not np.isfinite(mesh.depths).all():
        raise ValueError("Invalid positive-down bathymetry")
    if dt.min() < 2 * DTE:
        raise ValueError(f"Insufficient external-step margin: implied dt={dt.min()}")
    return dict(
        nodes=mesh.n_nodes,
        elements=mesh.n_elements,
        depth_min_m=float(mesh.depths.min()),
        depth_max_m=float(mesh.depths.max()),
        implied_dt_s=float(dt.min()),
    )


def prepare(run_root):
    a, b = load_a(), read_fort14(B_MESH)
    # The experiment only holds together if case B carries the SAME open
    # boundary as case A -- one segment of 13 nodes on the goto2023 arc, which
    # is an input to both.  Node and element counts are properties of the
    # candidate, so they are pinned only for the default mesh, where an
    # unannounced change would be a regression.
    if list(map(len, b.open_boundaries)) != [13]:
        raise ValueError(
            f"case B must have one 13-node open boundary, got "
            f"{list(map(len, b.open_boundaries))}")
    if B_MESH == DEFAULT_B_MESH and (b.n_nodes, b.n_elements) != (3393, 5849):
        raise ValueError(
            f"the repository mesh changed: {b.n_nodes} nodes, {b.n_elements} "
            "elements, expected 3393/5849")
    period, gauges = tide_constants()
    arc = a.nodes[a.open_boundaries[0]]
    ll_to_xy = Transformer.from_crs(4326, 32654, always_xy=True)
    anchors = [gauges[n] for n in ("ABURATUBO", "MERA")]
    anchor_xy = np.array([ll_to_xy.transform(g["lon"], g["lat"]) for g in anchors])
    anchor_s, anchor_offset = arc_position(anchor_xy, arc)
    order = np.argsort(anchor_s)
    if np.ptp(anchor_s) < 1000:
        raise ValueError("Tidal anchors do not span the mouth")
    # Case B_m7001 carries depths built the way A's were: the M7001 survey on
    # the T.P. datum, floored at 3 m, r-factor smoothed to r <= 0.2, capped at
    # 300 m (TB-FVCOM MESH.md / build_bathy_variants.py).  Interpolating A's
    # FINISHED node depths instead would carry A's own smoothing and A's own
    # coastline into B, which is the thing the comparison is trying to hold
    # apart.  Validated on A's mesh: median -0.64 m, MAE 1.14 m against the
    # production file, same 3-300 m range.
    xy_to_ll = Transformer.from_crs(32654, 4326, always_xy=True)
    b_lon, b_lat = xy_to_ll.transform(b.nodes[:, 0], b.nodes[:, 1])
    ba_depth, ba_report = production_depths(b_lon, b_lat, b.elements)
    print("B_m7001 depth:", json.dumps(ba_report), flush=True)
    ba = Fort14Mesh(
        "B M7001 production-recipe depth",
        b.nodes,
        ba_depth,
        b.elements,
        b.open_boundaries,
        b.land_boundaries,
    )
    # FVCOM overwrites each OBC node's depth with its NEXT_OBC depth at start-up
    # (OBC_DEPTH_CONTROL_ON). Apply it here so the staged bathymetry is exactly
    # what FVCOM runs with and 384 can require output h == staged depth.
    obc_depth_change = {}
    a, obc_depth_change["A"] = apply_obc_depth_control(a)
    b, obc_depth_change["B_own"] = apply_obc_depth_control(b)
    ba, obc_depth_change["B_m7001"] = apply_obc_depth_control(ba)
    run_root.mkdir(parents=True, exist_ok=True)
    metadata = dict(
        epoch_utc=START,
        end_utc=END,
        spinup_seconds=SPINUP,
        duration_seconds=DURATION,
        period_seconds=period,
        dte_seconds=DTE,
        isplit=ISPLIT,
        sigma_layers=5,
        roughness_length_m=0.002693138,
        minimum_drag_coefficient=0.003,
        ramp_tanh_timescale_seconds=RAMP_SECONDS,
        physics=(
            "3D constant T=20 C/S=30; inactive scalars; closure mixing; "
            "production sponge mapped by arc position"
        ),
        ranks=8,
        phase_convention="eta=A*cos(2*pi*seconds_since_epoch/period-phase)",
        observation_kind="Published JCG/JMA harmonics; frozen 2021 nodal corrections",
        forcing_source=str(TIDES / "ptide/Tide-ToKYOWAN.txt"),
        gauges=gauges,
        anchor_names=["ABURATUBO", "MERA"],
        anchor_arc_m=anchor_s.tolist(),
        anchor_offset_m=anchor_offset.tolist(),
        bathymetry_A=str(DEP),
        bathymetry_B="SRTM15 Kanto, min 2 m; notebooks/325",
        bathymetry_B_m7001="M7001 T.P. survey, min 3 m, r-factor <= 0.2, cap 300 m "
                           "(A's own recipe, applied to B's nodes)",
        bathymetry_B_m7001_report=ba_report,
        runs={},
    )
    # Production sponge (TokyoBay_spg.dat: per-OBC-node radius and damping).
    # Without it every run, the production mesh A included, blew up within
    # ~1 km of the open boundary after 3-8 days (job 115180). The same
    # radius/coefficient profile is mapped to each mesh by along-arc position.
    spg_lines = (TB / "grid/TokyoBay_spg.dat").read_text().split("\n")[1:]
    spg = {int(t[0]) - 1: (float(t[1]), float(t[2])) for t in (ln.split() for ln in spg_lines) if t}
    a_obc = a.open_boundaries[0]
    if sorted(spg) != sorted(a_obc.tolist()):
        raise ValueError("TokyoBay_spg.dat nodes differ from the goto2023 OBC nodes")
    a_pos, _ = arc_position(a.nodes[a_obc], arc)
    a_order = np.argsort(a_pos)
    spg_radius = np.array([spg[n][0] for n in a_obc])[a_order]
    spg_coef = np.array([spg[n][1] for n in a_obc])[a_order]
    for label, mesh in [("A", a), ("B_own", b), ("B_m7001", ba)]:
        case = run_root / label
        inp, out = case / "input", case / "output"
        inp.mkdir(parents=True, exist_ok=True)
        out.mkdir(exist_ok=True)
        if label == "A":
            for kind in ("grd", "cor", "obc"):
                shutil.copyfile(TB / f"grid/TokyoBay_{kind}.dat", inp / f"m2_{kind}.dat")
            write_dep(a, inp / "m2_dep.dat")
        else:
            subprocess.run(
                [
                    "fmesh-export-fvcom",
                    str(B_MESH),
                    "--outdir",
                    str(inp),
                    "--casename",
                    "m2",
                    "--cor",
                    "crs",
                    "--crs",
                    "EPSG:32654",
                    "--write-empty-spg",
                    "--no-2dm",
                ],
                check=True,
            )
            if label == "B_m7001":
                write_dep(mesh, inp / "m2_dep.dat")
        (inp / "sigma.dat").write_text(
            "NUMBER OF SIGMA LEVELS = 6\nSIGMA COORDINATE TYPE = UNIFORM\n"
        )
        obc = mesh.open_boundaries[0]
        pos, offset = arc_position(mesh.nodes[obc], arc)
        if offset.max() > 100:
            raise ValueError(f"{label} boundary is not on the reference arc: {offset.max()} m")
        radius = np.interp(pos, a_pos[a_order], spg_radius)
        coef = np.interp(pos, a_pos[a_order], spg_coef)
        (inp / "m2_spg.dat").write_text(
            f"Sponge Node Number = {len(obc)}\n"
            + "".join(f"{n + 1} {r:.6f} {c:.6f}\n" for n, r, c in zip(obc, radius, coef))
        )
        amp = np.interp(pos, anchor_s[order], np.array([g["amplitude_m"] for g in anchors])[order])
        # Unwrap before interpolation; endpoint extension is explicitly recorded.
        lag = np.unwrap(np.deg2rad([anchors[i]["phase_deg"] for i in order]))
        phase = np.rad2deg(np.interp(pos, anchor_s[order], lag)) % 360
        (inp / "m2_tide.dat").write_text(spectral_text(period, amp, phase))
        (case / "m2_run.nml").write_text(namelist(inp, out))
        metrics = mesh_metrics(mesh)
        metrics.update(
            obc_depth_control_change_m=obc_depth_change[label].tolist(),
            obc_node_ids=(obc + 1).tolist(),
            obc_arc_m=pos.tolist(),
            obc_offset_m=offset.tolist(),
            obc_amplitude_m=amp.tolist(),
            obc_phase_deg=phase.tolist(),
            anchor_endpoint_extension_nodes=int(
                ((pos < min(anchor_s)) | (pos > max(anchor_s))).sum()
            ),
        )
        metadata["runs"][label] = metrics
        np.savez(case / "mesh.npz", xy=mesh.nodes, tri=mesh.elements, depth=mesh.depths, obc=obc)
        print(label, json.dumps(metrics), flush=True)
    paths = [
        DEP,
        B_MESH,
        TEMPLATE,
        TB / "grid/TokyoBay_grd.dat",
        TB / "grid/TokyoBay_obc.dat",
        TB / "grid/TokyoBay_cor.dat",
        TIDES / "ptide/Tide-ToKYOWAN.txt",
        TIDES / "ptide/ptide.py",
        Path(__file__).resolve(),
    ]
    metadata["sha256"] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (run_root / "manifest.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return metadata


def main():
    global B_MESH
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_RUN_ROOT)
    parser.add_argument("--mesh", type=Path, default=DEFAULT_B_MESH,
                        help="fort.14 for case B (default: the repository mesh)")
    args = parser.parse_args()
    B_MESH = args.mesh.resolve()
    if not B_MESH.exists():
        raise SystemExit(f"mesh not found: {B_MESH}")
    print(f"[383] case B mesh = {B_MESH}", flush=True)
    prepare(args.root.resolve())


if __name__ == "__main__":
    main()
