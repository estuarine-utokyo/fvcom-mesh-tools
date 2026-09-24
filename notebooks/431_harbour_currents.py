# The harbour with and without walls: the tide's currents around the
# breakwaters, from the 20-day M2 runs of the two port meshes.
#
# The meshes differ, so nothing is interpolated from one to the other: each
# run is drawn on its own mesh.  Every solid boundary -- coast, quay, wall --
# is black (plotting.draw_mesh); a split wall is still counted, as a segment
# that is boundary TWICE, once from each side.
#
#   FMESH_NOWALL=<run dir> FMESH_WALL=<run dir> \
#       python notebooks/431_harbour_currents.py <output dir>
import json
import os
import sys
from pathlib import Path

import numpy as np

os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt  # noqa: E402
import netCDF4  # noqa: E402
from matplotlib.tri import Triangulation  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from fvcom_mesh_tools.plotting import draw_mesh  # noqa: E402

BOX = (391900.0, 394100.0, 3908350.0, 3910600.0)   # the Kimitsu port, UTM 54N
PERIOD = 12.4206012 * 3600.0


def load(run: Path):
    with netCDF4.Dataset(run / "output" / "m2_0001.nc") as ds:
        x = np.asarray(ds["x"][:], float)
        y = np.asarray(ds["y"][:], float)
        nv = np.asarray(ds["nv"][:], int).T - 1
        t = np.asarray(ds["time"][:], float) * 86400.0
        keep = t >= t[-1] - 5 * 86400.0            # the last five days
        ua = np.asarray(ds["ua"][keep], float)
        va = np.asarray(ds["va"][keep], float)
        zeta = np.asarray(ds["zeta"][keep], float)
        t = t[keep]
    w = 2 * np.pi / PERIOD
    a = np.column_stack([np.ones_like(t), np.cos(w * t), np.sin(w * t)])
    coef, *_ = np.linalg.lstsq(a, zeta, rcond=None)
    return {"x": x, "y": y, "nv": nv, "speed": np.hypot(ua, va).max(axis=0),
            "amp": np.hypot(coef[1], coef[2])}


def boundary_segments(x, y, nv):
    e = np.sort(np.vstack([nv[:, [0, 1]], nv[:, [1, 2]], nv[:, [2, 0]]]), axis=1)
    u, c = np.unique(e, axis=0, return_counts=True)
    b = u[c == 1]
    seg = np.stack([np.column_stack([x[b[:, 0]], y[b[:, 0]]]),
                    np.column_stack([x[b[:, 1]], y[b[:, 1]]])], axis=1)
    r = np.round(seg, 3)
    swap = (r[:, 0, 0] > r[:, 1, 0]) | ((r[:, 0, 0] == r[:, 1, 0]) & (r[:, 0, 1] > r[:, 1, 1]))
    r[swap] = r[swap][:, ::-1]
    key = r.reshape(-1, 4)
    _, inv, cnt = np.unique(key, axis=0, return_inverse=True, return_counts=True)
    twice = cnt[inv] == 2
    return seg[~twice], seg[twice]


def main(out: Path):
    runs = {"no walls": Path(os.environ["FMESH_NOWALL"]),
            "walls (split)": Path(os.environ["FMESH_WALL"])}
    data = {k: load(v) for k, v in runs.items()}
    x0, x1, y0, y1 = BOX
    fig, axes = plt.subplots(2, 2, figsize=(14, 13), constrained_layout=True)
    summary = {}
    for col, (name, d) in enumerate(data.items()):
        tri = Triangulation(d["x"], d["y"], d["nv"])
        coast, walls = boundary_segments(d["x"], d["y"], d["nv"])
        inside = (d["x"] > x0) & (d["x"] < x1) & (d["y"] > y0) & (d["y"] < y1)
        ec = d["nv"][:, 0]
        in_el = inside[ec]
        summary[name] = {"n_wall_segments": int(len(walls)),
                         "amp_m_in_box": [float(d["amp"][inside].min()),
                                          float(d["amp"][inside].max())],
                         "max_speed_ms_in_box": float(d["speed"][in_el].max())}
        for row, (field, label, kw) in enumerate((
                ("speed", "max depth-averaged speed, last 5 days (m/s)",
                 {"vmin": 0, "vmax": 0.3, "cmap": "viridis"}),
                ("amp", "M2 amplitude (m)", {"cmap": "plasma"}))):
            ax = axes[row, col]
            if field == "speed":
                pc = ax.tripcolor(tri, facecolors=d["speed"], **kw)
            else:
                v = d["amp"][inside]
                pc = ax.tripcolor(tri, d["amp"], shading="gouraud",
                                  vmin=float(np.percentile(v, 1)),
                                  vmax=float(np.percentile(v, 99)), cmap=kw["cmap"])
            draw_mesh(ax, np.column_stack([d["x"], d["y"]]), d["nv"], mesh_color="w", mesh_lw=0.15,
                      boundary_lw=1.2, zorder=3)
            ax.set_xlim(x0, x1)
            ax.set_ylim(y0, y1)
            ax.set_aspect("equal")
            ax.set_title(f"{name}: {label}")
            fig.colorbar(pc, ax=ax, shrink=0.8)
    png = out / "harbour_currents.png"
    fig.savefig(png, dpi=130)
    (out / "harbour_currents.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))
    print(f"wrote {png}")


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve())
