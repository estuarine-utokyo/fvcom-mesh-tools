"""Stage an M2 test of a local refinement: the base, and the base refined.

Two integrations that differ in ONE thing. 383 compares meshes built by
different pipelines and has to hold bathymetry, coastline and open boundary
apart; here the refined case IS the base case with a patch in it, so the
depths outside the patch are bit-identical, the coastline is identical, the
open boundary is the same node list in the same order, and any difference in
the answer is the patch.

Both cases run at the SAME external step, set by the refined mesh: its
shortest-edge CFL is 3.34 s against the base's 14.39 s, so a comparison at
each mesh's own step would be a comparison of two time steps as well as two
meshes.

Everything else -- tide, sponge, sigma, namelist -- is 383's, imported
rather than copied so the two experiments cannot drift apart.

    python notebooks/414_refine_m2_prep.py \\
        --root /octfs/work/.../scratch/m2r_<stamp> \\
        --base outputs/base_tool/TokyoBayTool \\
        --refined outputs/refine_futtsu_nori_tool/fvcom/futtsu_nori_tool
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from pyproj import Transformer

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from fvcom_mesh_tools.io.fvcom_native import (  # noqa: E402
    apply_obc_depth_control,
    export_fvcom_case,
    read_fvcom_case,
)

_spec = importlib.util.spec_from_file_location(
    "m2_383", ROOT / "notebooks" / "383_m2_case_prep.py")
M383 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(M383)

MESH_EPSG = 32654


def case_paths(prefix: Path) -> tuple[Path, Path, Path]:
    return (Path(f"{prefix}_grd.dat"), Path(f"{prefix}_dep.dat"),
            Path(f"{prefix}_obc.dat"))


def external_step(mesh, safety: float = 2.0) -> float:
    """The largest external step this mesh allows, by the shortest edge.

    383's rule, kept: ``implied dt >= safety * DTE``. The reported number is
    the shortest-edge measure, not the minimum altitude -- the altitude is the
    stricter of the two and is what the refinement report quotes, so this is
    deliberately the same convention 383 and FVCOM's own CFL use.
    """
    xy = mesh.nodes[mesh.elements]
    edge = np.linalg.norm(xy - np.roll(xy, 1, axis=1), axis=2).min(axis=1)
    dt = edge / np.sqrt(9.81 * mesh.depths[mesh.elements].max(axis=1))
    return float(dt.min()) / safety


def prepare(run_root: Path, cases: dict[str, Path], dte: float | None) -> dict:
    meshes = {}
    for label, prefix in cases.items():
        grd, dep, obc = case_paths(prefix)
        for p in (grd, dep, obc):
            if not p.exists():
                raise SystemExit(f"{label}: missing {p}")
        meshes[label] = read_fvcom_case(grd, dep, obc, title=f"m2 {label}")
        print(f"[414] {label}: NP={meshes[label].n_nodes:,} "
              f"NE={meshes[label].n_elements:,}, depth "
              f"{meshes[label].depths.min():.3f}-{meshes[label].depths.max():.3f} m, "
              f"max DTE {external_step(meshes[label]):.2f} s", flush=True)

    ref, cand = list(cases)
    a, b = meshes[ref], meshes[cand]
    if len(a.open_boundaries) != 1 or len(b.open_boundaries) != 1:
        raise SystemExit("both cases need exactly one open boundary")
    if not np.array_equal(a.nodes[a.open_boundaries[0]],
                          b.nodes[b.open_boundaries[0]]):
        raise SystemExit("the two cases do not share the same open-boundary arc; "
                         "the comparison would confound the patch with the forcing")

    # One step for both, and it is the refined mesh's.
    M383.DTE = float(dte) if dte else min(external_step(m) for m in meshes.values())
    M383.DTE = float(f"{M383.DTE:.3g}")
    print(f"[414] EXTSTEP_SECONDS = {M383.DTE} for both cases "
          f"(ISPLIT {M383.ISPLIT}, DTI {M383.DTE * M383.ISPLIT:g} s)", flush=True)

    period, gauges = M383.tide_constants()
    arc = a.nodes[a.open_boundaries[0]]
    ll_to_xy = Transformer.from_crs(4326, MESH_EPSG, always_xy=True)
    anchors = [gauges[n] for n in ("ABURATUBO", "MERA")]
    anchor_xy = np.array([ll_to_xy.transform(g["lon"], g["lat"]) for g in anchors])
    anchor_s, anchor_offset = M383.arc_position(anchor_xy, arc)
    order = np.argsort(anchor_s)

    # The production sponge, mapped by along-arc position exactly as 383 does.
    # Without it every mesh blew up within ~1 km of the open boundary.
    spg_lines = (M383.TB / "grid/TokyoBay_spg.dat").read_text().split("\n")[1:]
    spg = {int(t[0]) - 1: (float(t[1]), float(t[2]))
           for t in (ln.split() for ln in spg_lines) if t}
    g_grd = M383.TB / "grid/TokyoBay_grd.dat"
    lines = g_grd.read_text().splitlines()
    nn, ne = (int(ln.split("=")[1]) for ln in lines[:2])
    g_xy = np.array([s.split()[1:3] for s in lines[2 + ne: 2 + ne + nn]], float)
    g_obc = np.array([int(ln.split()[1]) - 1 for ln in
                      (M383.TB / "grid/TokyoBay_obc.dat").read_text().splitlines()[1:]])
    if sorted(spg) != sorted(g_obc.tolist()):
        raise SystemExit("TokyoBay_spg.dat nodes differ from the goto2023 OBC nodes")
    g_pos, _ = M383.arc_position(g_xy[g_obc], arc)
    g_order = np.argsort(g_pos)
    spg_radius = np.array([spg[n][0] for n in g_obc])[g_order]
    spg_coef = np.array([spg[n][1] for n in g_obc])[g_order]

    run_root.mkdir(parents=True, exist_ok=True)
    to_ll = Transformer.from_crs(MESH_EPSG, 4326, always_xy=True)
    metadata = dict(
        epoch_utc=M383.START, end_utc=M383.END,
        spinup_seconds=M383.SPINUP, duration_seconds=M383.DURATION,
        period_seconds=period, dte_seconds=M383.DTE, isplit=M383.ISPLIT,
        sigma_layers=5, roughness_length_m=0.002693138,
        minimum_drag_coefficient=0.003, ramp_tanh_timescale_seconds=86400,
        physics="3D constant T=20 C/S=30; inactive scalars; closure mixing; "
                "production sponge mapped by arc position",
        ranks=64,
        phase_convention="eta=A*cos(2*pi*seconds_since_epoch/period-phase)",
        observation_kind="Published JCG/JMA harmonics; frozen 2021 nodal corrections",
        forcing_source=str(M383.TIDES / "ptide/Tide-ToKYOWAN.txt"),
        gauges=gauges, anchor_names=["ABURATUBO", "MERA"],
        anchor_arc_m=anchor_s.tolist(), anchor_offset_m=anchor_offset.tolist(),
        experiment="local refinement: the base, and the base refined",
        cases={k: str(v) for k, v in cases.items()},
        reference_case=ref,
        runs={},
    )

    for label in cases:
        mesh = meshes[label]
        mesh, obc_change = apply_obc_depth_control(mesh)
        case = run_root / label
        inp, out = case / "input", case / "output"
        inp.mkdir(parents=True, exist_ok=True)
        out.mkdir(exist_ok=True)
        _, lat = to_ll.transform(mesh.nodes[:, 0], mesh.nodes[:, 1])
        export_fvcom_case(mesh, inp, "m2", cor=lat, twodm=False,
                          obc_depth_control=False)
        (inp / "sigma.dat").write_text(
            "NUMBER OF SIGMA LEVELS = 6\nSIGMA COORDINATE TYPE = UNIFORM\n")
        obc = mesh.open_boundaries[0]
        pos, offset = M383.arc_position(mesh.nodes[obc], arc)
        if offset.max() > 100:
            raise SystemExit(f"{label} boundary is not on the reference arc: "
                             f"{offset.max()} m")
        radius = np.interp(pos, g_pos[g_order], spg_radius)
        coef = np.interp(pos, g_pos[g_order], spg_coef)
        (inp / "m2_spg.dat").write_text(
            f"Sponge Node Number = {len(obc)}\n"
            + "".join(f"{n + 1} {r:.6f} {c:.6f}\n"
                      for n, r, c in zip(obc, radius, coef)))
        amp = np.interp(pos, anchor_s[order],
                        np.array([g["amplitude_m"] for g in anchors])[order])
        lag = np.unwrap(np.deg2rad([anchors[i]["phase_deg"] for i in order]))
        phase = np.rad2deg(np.interp(pos, anchor_s[order], lag)) % 360
        (inp / "m2_tide.dat").write_text(M383.spectral_text(period, amp, phase))
        (case / "m2_run.nml").write_text(M383.namelist(inp, out))
        metrics = M383.mesh_metrics(mesh)
        metrics.update(
            obc_depth_control_change_m=obc_change.tolist(),
            obc_node_ids=(obc + 1).tolist(),
            obc_arc_m=pos.tolist(), obc_offset_m=offset.tolist(),
            obc_amplitude_m=amp.tolist(), obc_phase_deg=phase.tolist(),
            anchor_endpoint_extension_nodes=int(
                ((pos < min(anchor_s)) | (pos > max(anchor_s))).sum()),
        )
        metadata["runs"][label] = metrics
        np.savez(case / "mesh.npz", xy=mesh.nodes, tri=mesh.elements,
                 depth=mesh.depths, obc=obc)
        print(f"[414] {label}: " + json.dumps(
            {k: v for k, v in metrics.items() if not isinstance(v, list)}), flush=True)

    files = [p for prefix in cases.values() for p in case_paths(prefix)]
    files += [M383.TEMPLATE, M383.TB / "grid/TokyoBay_spg.dat",
              M383.TIDES / "ptide/Tide-ToKYOWAN.txt",
              Path(__file__).resolve(), ROOT / "notebooks/383_m2_case_prep.py"]
    metadata["sha256"] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in files}
    (run_root / "manifest.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--base", type=Path, required=True,
                        help="prefix of the base FVCOM case (without _grd.dat)")
    parser.add_argument("--refined", type=Path, required=True,
                        help="prefix of the refined FVCOM case")
    parser.add_argument("--dte", type=float, default=None,
                        help="external step in seconds; default is the "
                             "smallest mesh's allowance")
    args = parser.parse_args()
    prepare(args.root.resolve(),
            {"base": args.base.resolve(), "refined": args.refined.resolve()},
            args.dte)


if __name__ == "__main__":
    main()
