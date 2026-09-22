"""Analyze 383 FVCOM output; run in the batch job after all three integrations.

Phases are cosine lags relative to the manifest UTC epoch, not Greenwich
harmonic constants. Volume change includes real exchange through the open
boundary: the detided trend is a drift diagnostic, not a conservation residual.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import netCDF4
import numpy as np
from pyproj import Transformer
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from fvcom_mesh_tools.plotting import add_atlas_grid, use_readable_style  # noqa: E402

DEFAULT_RUN_ROOT = Path("/octfs/work/G16445/v61021/scratch/m2_383").resolve()
# Set from the manifest at run time. The cases are a property of the
# experiment, not of this script: 383 stages A / B_own / B_m7001, and 414
# stages base / refined on one base. The first label is the reference every
# other one is differenced against.
LABELS: tuple[str, ...] = ("A", "B_own", "B_m7001")
STATIONS = ("TOKYO-SIBAURA", "HARUMI", "TIBA-TIBA LIGHT", "SINKO", "YOKOSUKA")


def design(t, period, trend=False):
    t = np.asarray(t, float)
    cols = [np.ones(len(t))]
    for multiple in (1, 2, 3):  # M2, M4, M6; no independent MS4 in M2-only forcing.
        arg = multiple * 2 * np.pi * t / period
        cols.extend([np.cos(arg), np.sin(arg)])
    if trend:
        cols.append((t - t.mean()) / 86400)
    return np.column_stack(cols)


def harmonic_fit(t, z, period, trend=False):
    if not np.isfinite(z).all():
        raise ValueError("Nonfinite harmonic input")
    matrix = design(t, period, trend)
    if len(t) < 2 * matrix.shape[1] or np.linalg.matrix_rank(matrix) < matrix.shape[1]:
        raise ValueError("Insufficient independent harmonic samples")
    coef = np.linalg.lstsq(matrix, z, rcond=None)[0]
    amp = np.hypot(coef[1], coef[2])
    phase = np.degrees(np.arctan2(coef[2], coef[1])) % 360
    return amp, phase, coef


def station_weights(xy, tri, station_xy, wet):
    """Barycentric weights of each station inside its containing element.

    The nearest WET NODE is not a fair sample: two meshes put their nearest
    node in different places -- up to 1,062 m from the gauge in the Tokyo Bay
    comparison -- so part of the difference between cases was just the
    difference between two sampling points.  Interpolating inside the element
    that contains the gauge samples both meshes at the SAME position.

    Returns ``(node_ids, weights, distance_m)`` per station, where distance is
    0 for a station inside the mesh and the nearest-node distance otherwise
    (a station outside every element falls back to that node, weight 1).
    """
    p = xy[tri]
    v0 = p[:, 1] - p[:, 0]
    v1 = p[:, 2] - p[:, 0]
    den = v0[:, 0] * v1[:, 1] - v1[:, 0] * v0[:, 1]
    ids = np.zeros((len(station_xy), 3), dtype=int)
    w = np.zeros((len(station_xy), 3))
    dist = np.zeros(len(station_xy))
    candidates = np.flatnonzero(wet)
    tree = cKDTree(xy[candidates])
    for k, sxy in enumerate(np.asarray(station_xy, float)):
        v2 = sxy - p[:, 0]
        b1 = (v2[:, 0] * v1[:, 1] - v1[:, 0] * v2[:, 1]) / den
        b2 = (v0[:, 0] * v2[:, 1] - v2[:, 0] * v0[:, 1]) / den
        b0 = 1.0 - b1 - b2
        inside = (b0 >= -1e-9) & (b1 >= -1e-9) & (b2 >= -1e-9)
        hit = np.flatnonzero(inside & wet[tri].all(axis=1))
        if hit.size:
            e = int(hit[0])
            ids[k] = tri[e]
            w[k] = [b0[e], b1[e], b2[e]]
            dist[k] = 0.0
        else:
            d, local = tree.query(sxy)
            node = int(candidates[local])
            ids[k] = node
            w[k] = [1.0, 0.0, 0.0]
            dist[k] = float(d)
    return ids, w, dist


def phase_difference(a, b):
    return (np.asarray(a) - np.asarray(b) + 180) % 360 - 180


def node_areas(xy, tri):
    p = xy[tri]
    d, e = p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]
    area = np.abs(d[:, 0] * e[:, 1] - d[:, 1] * e[:, 0]) / 2
    if np.any(area <= 0):
        raise ValueError("Degenerate triangles")
    return np.bincount(tri.ravel(), weights=np.repeat(area / 3, 3), minlength=len(xy))


def numeric(v):
    return np.issubdtype(v.dtype, np.number)


def read_run(case, manifest):
    mesh = dict(np.load(case / "mesh.npz"))
    files = sorted((case / "output").glob("m2_[0-9][0-9][0-9][0-9].nc"))
    if not files:
        raise ValueError(f"No FVCOM output in {case}")
    ts, zs, us, vs = [], [], [], []
    bad = {}
    wet = np.ones(len(mesh["xy"]), bool)
    for path in files:
        with netCDF4.Dataset(path) as ds:
            # Check all emitted numeric fields, including turbulence/scalars/vertical velocity.
            # Scan by time record so large velocity arrays do not accumulate in memory.
            for name, var in ds.variables.items():
                if not numeric(var):
                    continue
                slices = range(var.shape[0]) if var.dimensions[:1] == ("time",) else [None]
                for i in slices:
                    data = np.ma.asarray(var[:] if i is None else var[i])
                    count = int(np.ma.count_masked(data) + (~np.isfinite(data.compressed())).sum())
                    if count:
                        bad[name] = bad.get(name, 0) + count
            time = ds["time"]
            dates = netCDF4.num2date(
                time[:], time.units, calendar=getattr(time, "calendar", "standard")
            )
            t = netCDF4.date2num(dates, "seconds since " + manifest["epoch_utc"])
            ts.extend(t)
            z = np.ma.asarray(ds["zeta"][:]).filled(np.nan)
            if ds["zeta"].dimensions != ("time", "node"):
                raise ValueError("Unexpected zeta dimensions")
            if z.shape[1] != len(mesh["xy"]):
                raise ValueError("Output mesh does not match prepared inputs")
            # Compare topology by sorted vertex IDs, allowing FVCOM's CW reversal.
            output_tri = np.asarray(ds["nv"][:]).T.astype(int) - 1
            if not np.array_equal(np.sort(output_tri, axis=1), np.sort(mesh["tri"], axis=1)):
                raise ValueError("Output connectivity differs from prepared mesh")
            # Every node, OBC included: the prep applies FVCOM's OBC depth
            # control up front, so FVCOM must not change any depth.
            if not np.allclose(ds["h"][:], mesh["depth"], atol=1e-3, rtol=0):
                raise ValueError("Output bathymetry differs from prepared mesh")
            zs.extend(z)
            wet &= np.all(z + mesh["depth"] > 0.1, axis=0)
            if "wet_nodes" in ds.variables:
                wet &= np.all(np.asarray(ds["wet_nodes"][:]) > 0, axis=0)
            for i in range(len(t)):
                u = np.ma.asarray(ds["u"][i]).filled(np.nan)
                v = np.ma.asarray(ds["v"][i]).filled(np.nan)
                us.append(float(np.max(np.abs(u))))
                vs.append(float(np.max(np.hypot(u, v))))
    t, z = np.asarray(ts), np.asarray(zs)
    order = np.argsort(t)
    t, z = t[order], z[order]
    health = dict(
        files=[str(p) for p in files],
        records=len(t),
        nonfinite_or_masked=bad,
        finite_pass=not bad,
        max_abs_zeta_m=float(np.max(np.abs(z))),
        max_abs_u_m_s=float(np.max(us)),
        max_horizontal_speed_m_s=float(np.max(vs)),
    )
    if bad:
        return mesh, t, z, wet, health
    if len(t) < 2 or np.any(np.diff(t) <= 0) or np.max(np.diff(t)) > 1801:
        raise ValueError(f"{case.name}: duplicate or missing output times")
    if abs(t[0]) > 1 or abs(t[-1] - manifest["duration_seconds"]) > 1:
        raise ValueError(f"{case.name}: incomplete run: {t[0]} .. {t[-1]} seconds")
    select = t >= manifest["spinup_seconds"]
    # What the fit needs is enough M2 cycles to separate M2, M4 and M6, not a
    # fixed number of days.  A fixed 7-day floor also forces a LONGER window
    # than is wanted once the spin-up is the long part of the run: averaging
    # over more of a still-converging signal biases the constants low.  Four
    # cycles (about 2.1 days) is the bar; 20 days with 15 of spin-up leaves
    # 5 days, or 9.7 cycles.
    MIN_CYCLES = 4.0
    cycles = np.ptp(t[select]) / manifest["period_seconds"]
    if cycles < MIN_CYCLES:
        raise ValueError(
            f"analysis window is {cycles:.1f} M2 cycles, need {MIN_CYCLES:g}")
    weights = node_areas(mesh["xy"], mesh["tri"])
    volume = (z + mesh["depth"]) @ weights
    _, _, coef = harmonic_fit(t[select], volume[select], manifest["period_seconds"], trend=True)
    drift = float(coef[-1])
    health.update(
        complete_pass=True,
        always_wet_nodes=int(wet.sum()),
        initial_volume_m3=float(volume[0]),
        final_volume_m3=float(volume[-1]),
        volume_change_m3=float(volume[-1] - volume[0]),
        volume_change_fraction=float((volume[-1] - volume[0]) / volume[0]),
        volume_range_m3=[float(volume.min()), float(volume.max())],
        detided_analysis_drift_m3_per_day=drift,
        detided_analysis_drift_fraction_per_day=float(drift / volume[0]),
        mass_change_kg_at_rho1025=float(1025 * (volume[-1] - volume[0])),
        volume_note="Open-domain storage, not flux-budget conservation error; rho=1025 kg/m3",
        extrema_note="30-minute snapshots; also inspect the FVCOM log for intermediate failures",
    )
    return mesh, t, z, wet, health


def analyze(run_root, output, figure):
    global LABELS

    manifest = json.loads((run_root / "manifest.json").read_text())
    LABELS = tuple(manifest["runs"])
    print(f"[384] cases: {', '.join(LABELS)} (reference = {LABELS[0]})", flush=True)
    output.mkdir(parents=True, exist_ok=True)
    rows, health, maps = {}, {}, {}
    tr = Transformer.from_crs(4326, 32654, always_xy=True)
    ll = Transformer.from_crs(32654, 4326, always_xy=True)
    station_xy = np.array(
        [tr.transform(manifest["gauges"][s]["lon"], manifest["gauges"][s]["lat"]) for s in STATIONS]
    )
    for label in LABELS:
        try:
            mesh, t, z, wet, health[label] = read_run(run_root / label, manifest)
            if not health[label]["finite_pass"]:
                continue
            select = t >= manifest["spinup_seconds"]
            amp, phase, coef = harmonic_fit(t[select], z[select], manifest["period_seconds"])
            # Independent halves show whether the nominal spin-up was sufficient.
            half = (t[select][0] + t[-1]) / 2
            a1, p1, _ = harmonic_fit(
                t[select & (t < half)], z[select & (t < half)], manifest["period_seconds"]
            )
            a2, p2, _ = harmonic_fit(t[t >= half], z[t >= half], manifest["period_seconds"])
            if not np.flatnonzero(wet).size:
                raise ValueError("No nodes remain wet throughout the run")
            sids, sw, sdist = station_weights(mesh["xy"], mesh["tri"], station_xy, wet)
            # Interpolate the harmonic COEFFICIENTS, not amplitude and phase:
            # averaging phases across a node triple is meaningless near 0/360.
            c_cos = (coef[1][sids] * sw).sum(axis=1)
            c_sin = (coef[2][sids] * sw).sum(axis=1)
            s_amp = np.hypot(c_cos, c_sin)
            s_phase = np.degrees(np.arctan2(c_sin, c_cos)) % 360
            s_mean = (coef[0][sids] * sw).sum(axis=1)
            h1c = (a1[sids] * sw).sum(axis=1)
            h2c = (a2[sids] * sw).sum(axis=1)
            hp1 = (p1[sids] * sw).sum(axis=1)
            hp2 = (p2[sids] * sw).sum(axis=1)
            lon, lat = ll.transform(mesh["xy"][:, 0], mesh["xy"][:, 1])
            maps[label] = (np.column_stack([lon, lat]), mesh["tri"], amp, phase)
            for k, (station, distance) in enumerate(zip(STATIONS, sdist)):
                row = rows.setdefault(station, {"observed": manifest["gauges"][station]})
                obs = row["observed"]
                row[label] = dict(
                    node_id=int(sids[k, 0] + 1),
                    element_nodes=[int(n + 1) for n in sids[k]],
                    element_weights=[float(v) for v in sw[k]],
                    interpolated=bool(distance == 0.0),
                    distance_m=float(distance),
                    node_lon=float(lon[sids[k, 0]]),
                    node_lat=float(lat[sids[k, 0]]),
                    amplitude_m=float(s_amp[k]),
                    phase_deg=float(s_phase[k]),
                    amplitude_minus_observed_m=float(s_amp[k] - obs["amplitude_m"]),
                    phase_minus_observed_deg=float(
                        phase_difference(s_phase[k], obs["phase_deg"])),
                    half_window_amplitude_change_m=float(h2c[k] - h1c[k]),
                    half_window_phase_change_deg=float(phase_difference(hp2[k], hp1[k])),
                    mean_m=float(s_mean[k]),
                )
        except (ValueError, OSError, KeyError) as exc:
            health.setdefault(label, {}).update(error=str(exc), complete_pass=False)
    # SPIN-UP GATE.  The M2 response approaches steady state slowly (e-folding
    # ~3.8 days in Tokyo Bay), so a harmonic fit over a window that is still
    # growing reports a transient average, not an M2 constant -- which is how
    # job 115305's 11-day run came to understate every amplitude.  The two
    # halves of the analysis window must agree.
    TOL_M = 0.002
    convergence = {}
    for label in maps:
        changes = [abs(rows[s2][label]["half_window_amplitude_change_m"])
                   for s2 in rows if label in rows[s2]]
        worst = max(changes) if changes else float("nan")
        converged = bool(changes) and worst < TOL_M
        convergence[label] = {
            "max_half_window_amplitude_change_m": worst,
            "tolerance_m": TOL_M,
            "converged": converged,
        }
        print(f"[384] {label}: half-window amplitude change {worst * 1000:.2f} mm "
              f"({'converged' if converged else 'NOT CONVERGED -- integrate longer'})",
              flush=True)

    report = dict(
        design=manifest,
        health=health,
        spinup=convergence,
        stations=rows,
        caveats=[
            "Published harmonics are observational references, not contemporaneous measurements.",
            "No equivalence tolerance supplied; report differences without a match verdict.",
            "B_m7001: depths rebuilt from M7001 by A's recipe, not interpolated from A.",
            "Stations are interpolated inside the containing element, so every case "
            "is sampled at the gauge position rather than at its own nearest node.",
            "Trust the constants only where spinup.converged is true.",
        ],
    )
    if len(maps) == len(LABELS):
        ref = LABELS[0]
        pairs = [(b, ref) for b in LABELS[1:]]
        pairs += [(LABELS[1], LABELS[2])] if len(LABELS) == 3 else []
        for row in rows.values():
            row["differences"] = {}
            for b, a in pairs:
                row["differences"][f"{b}_minus_{a}"] = dict(
                    amplitude_m=row[b]["amplitude_m"] - row[a]["amplitude_m"],
                    phase_deg=float(phase_difference(row[b]["phase_deg"], row[a]["phase_deg"])),
                )

    # Preserve health diagnostics even when a run failed; never emit NaN JSON tokens.
    def clean(value):
        if isinstance(value, float) and not np.isfinite(value):
            return None
        if isinstance(value, dict):
            return {k: clean(v) for k, v in value.items()}
        if isinstance(value, list):
            return [clean(v) for v in value]
        return value

    (output / "comparison.json").write_text(
        json.dumps(clean(report), indent=2, allow_nan=False) + "\n"
    )
    print(json.dumps(clean(health), indent=2, allow_nan=False), flush=True)
    if len(maps) != 3:
        raise RuntimeError("Incomplete or unhealthy run(s); see comparison.json and FVCOM logs")
    # A gauge that sits OUTSIDE the mesh gets the nearest wet node instead, and
    # two meshes put that node in different places -- so at those stations the
    # A/B difference is partly the difference between two sampling points, not
    # between two meshes.  Say so rather than letting the number stand alone.
    outside = sorted({s2 for s2 in rows for label in LABELS
                      if label in rows[s2] and not rows[s2][label]["interpolated"]})
    if outside:
        print(f"[384] {len(outside)} of {len(rows)} gauges lie OUTSIDE the mesh and are "
              f"sampled at a nearby node, not at the gauge: {', '.join(outside)}",
              flush=True)
        for s2 in outside:
            d = {label: rows[s2][label]["distance_m"] for label in LABELS if label in rows[s2]}
            print("[384]   " + s2 + ": " + ", ".join(f"{k} {v:.0f} m" for k, v in d.items()),
                  flush=True)

    table = ["station,case,node,distance_m,in_mesh,M2_m,phase_deg,delta_obs_m,delta_obs_deg"]
    for s, row in rows.items():
        obs = row["observed"]
        table.append(f"{s},observed,,,,{obs['amplitude_m']:.6f},{obs['phase_deg']:.6f},0,0")
        for label in LABELS:
            r = row[label]
            table.append(
                f"{s},{label},{r['node_id']},{r['distance_m']:.1f},"
                f"{'yes' if r['interpolated'] else 'no'},"
                f"{r['amplitude_m']:.6f},{r['phase_deg']:.6f},"
                f"{r['amplitude_minus_observed_m']:.6f},{r['phase_minus_observed_deg']:.6f}"
            )
    text = "\n".join(table) + "\n"
    print(text)
    print(
        "PAIRWISE DIFFERENCES", json.dumps({s: r["differences"] for s, r in rows.items()}, indent=2)
    )
    (output / "stations.csv").write_text(text)
    plot(maps, rows, figure)


def plot(maps, rows, path):
    use_readable_style()
    fig = plt.figure(figsize=(19, 17), layout="constrained")
    grid = fig.add_gridspec(3, 3, height_ratios=(1, 1, 0.65))
    max_amp = max(v[2].max() for v in maps.values())
    bounds = np.concatenate([v[0] for v in maps.values()])
    for col, (label, (xy, tri, amp, phase)) in enumerate(maps.items()):
        triang = mtri.Triangulation(xy[:, 0], xy[:, 1], tri)
        for row, values in enumerate((amp, phase)):
            ax = fig.add_subplot(grid[row, col])
            ax.set_facecolor("0.90")
            artist = ax.tripcolor(
                triang,
                values,
                shading="gouraud",
                cmap="viridis" if row == 0 else "twilight",
                vmin=0,
                vmax=max_amp if row == 0 else 360,
            )
            # Phase contours are masked at the 0/360 seam.
            if row == 0:
                seam = np.ptp(phase[tri], axis=1) > 180
                contour_tri = mtri.Triangulation(xy[:, 0], xy[:, 1], tri, mask=seam)
                cs = ax.tricontour(
                    contour_tri, phase, levels=np.arange(0, 361, 15), colors="white", linewidths=0.5
                )
                ax.clabel(cs, fontsize=7, fmt="%d°")
            for i, s in enumerate(STATIONS):
                g = rows[s]["observed"]
                ax.plot(g["lon"], g["lat"], "r.", markersize=5)
                ax.annotate(str(i + 1), (g["lon"], g["lat"]), fontsize=8)
            ax.set_xlim(bounds[:, 0].min() - 0.01, bounds[:, 0].max() + 0.01)
            ax.set_ylim(bounds[:, 1].min() - 0.01, bounds[:, 1].max() + 0.01)
            ax.set_aspect(1 / np.cos(np.deg2rad(35.35)))
            add_atlas_grid(ax, crs="EPSG:4326")
            ax.set_title(f"{label}: M2 {'amplitude' if row == 0 else 'phase lag'}")
            fig.colorbar(artist, ax=ax, label="m" if row == 0 else "degrees from UTC epoch")
    axa, axp = fig.add_subplot(grid[2, :2]), fig.add_subplot(grid[2, 2])
    x = np.arange(len(STATIONS))
    for label in (*LABELS, "observed"):
        axa.plot(x, [rows[s][label]["amplitude_m"] for s in STATIONS], "o-", label=label)
        # Plot model-observed phase differences, avoiding wrapped phase chart artifacts.
        if label != "observed":
            axp.plot(
                x, [rows[s][label]["phase_minus_observed_deg"] for s in STATIONS], "o-", label=label
            )
    for ax in (axa, axp):
        ax.set_xticks(x, [f"{i + 1}. {s}" for i, s in enumerate(STATIONS)], rotation=25, ha="right")
        ax.legend(fontsize=9)
        ax.grid(alpha=0.25)
    axa.set_ylabel("M2 amplitude (m)")
    axp.set_ylabel("Model − observed phase (degrees)")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"Saved {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_RUN_ROOT)
    parser.add_argument("--output", type=Path, default=ROOT / "outputs/m2_383")
    parser.add_argument(
        "--figure", type=Path, default=ROOT / "outputs/figures/383_m2_comparison.png"
    )
    args = parser.parse_args()
    analyze(args.root.resolve(), args.output.resolve(), args.figure.resolve())


if __name__ == "__main__":
    main()
