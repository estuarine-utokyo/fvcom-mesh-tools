# Build a base case from THIS project's mesh, with the production bathymetry.
#
#   python notebooks/422_tool_base.py \
#       outputs/verify_409.115302/fit/sample_repro_final.14 outputs/base_tool
#
# Local refinement takes a finished FVCOM case as its base. The goto2023 case
# is one; this makes the other, so the two can be refined the same way and
# compared: the mesh comes from this project's pipeline (notebook 325 ->
# 331, coast-fitted), and the depths are rebuilt by the SAME recipe the
# `current` hydro baseline names -- m7001tp_rfac0p2_cap300, which is
# `dem.m7001.production_depths`: M7001 on the T.P. datum, a 3 m floor,
# r-factor smoothing to 0.2, then a 300 m cap.
#
# The SRTM15 depths the mesh arrives with are discarded. They are what
# notebook 325 interpolated to get a mesh at all, not a bathymetry anyone
# intends to run: min 2.00 m against the 3 m floor, max 735.3 m against the
# 300 m cap, 2,258 edges above r = 0.2.
import json
import sys
from pathlib import Path

import numpy as np
from pyproj import Transformer

from fvcom_mesh_tools.dem.m7001 import production_depths
from fvcom_mesh_tools.io.fort14 import read_fort14, write_fort14
from fvcom_mesh_tools.io.fvcom_native import export_fvcom_case, read_fvcom_case
from fvcom_mesh_tools.qa import run_qa

MESH_EPSG = 32654
src = Path(sys.argv[1] if len(sys.argv) > 1
           else "outputs/verify_409.115302/fit/sample_repro_final.14").resolve()
OUT = Path(sys.argv[2] if len(sys.argv) > 2 else "outputs/base_tool").resolve()
CASE = sys.argv[3] if len(sys.argv) > 3 else "TokyoBayTool"
OUT.mkdir(parents=True, exist_ok=True)


def say(msg):
    print(f"[base] {msg}", flush=True)


mesh = read_fort14(src)
say(f"{src.name}: NP={mesh.n_nodes:,} NE={mesh.n_elements:,}, "
    f"{len(mesh.open_boundaries)} open boundary")
say(f"  depths as built (SRTM15, from notebook 325): "
    f"{mesh.depths.min():.2f}-{mesh.depths.max():.2f} m")

to_ll = Transformer.from_crs(f"EPSG:{MESH_EPSG}", "EPSG:4326", always_xy=True)
lon, lat = to_ll.transform(mesh.nodes[:, 0], mesh.nodes[:, 1])
depth, report = production_depths(lon, lat, mesh.elements)
say("production_depths: " + json.dumps(report))

# FVCOM's _dep.dat carries six decimals, so the file cannot hold a raw
# double. Rounding BEFORE writing is what makes the file and the object the
# same thing -- and the file is the base, so everything downstream inherits
# exactly what the model will read. goto2023's own dep file is six decimals
# too. The change is at most 5e-7 m on depths of 3 m and up.
depth = np.round(depth, 6)
mesh.depths = depth
e = np.unique(np.sort(np.vstack([mesh.elements[:, [0, 1]], mesh.elements[:, [1, 2]],
                                 mesh.elements[:, [2, 0]]]), axis=1), axis=0)
r = np.abs(depth[e[:, 0]] - depth[e[:, 1]]) / (depth[e[:, 0]] + depth[e[:, 1]])
say(f"depth {depth.min():.3f}-{depth.max():.3f} m, r-factor max {r.max():.4f} "
    f"({int((r > 0.2 + 1e-9).sum())} edges above 0.2)")
if r.max() > 0.2 + 1e-6:
    # The cap is applied after the smoothing, so it can reintroduce a jump at
    # the 300 m boundary. The goto2023 product does not have one; if this
    # mesh does, it has to be said rather than shipped.
    say(f"WARNING: the 300 m cap reintroduced r = {r.max():.4f}")

# Written as an FVCOM case, the same three files the refinement driver reads
# for goto2023, so the two bases differ only in the mesh and its depths.
# Depth control at the open boundary is NOT applied: it rewrites OBC depths,
# and the point of this file is that its depths are the recipe's.
written = export_fvcom_case(mesh, OUT, CASE, cor=lat, obc_depth_control=False)
write_fort14(mesh, OUT / f"{CASE}.14")
back = read_fvcom_case(written["grd"], written["dep"], written["obc"])
if not np.array_equal(back.depths, depth):
    raise SystemExit("the written case does not carry the depths that were built")
say(f"wrote {', '.join(sorted(written))} + {CASE}.14 in {OUT}")

qa = run_qa(back, name=CASE, path=written["grd"], max_offenders=10_000)
(OUT / f"{CASE}_qa.json").write_text(json.dumps(qa.to_dict(), indent=1, default=float))
say(f"QA {qa.n_gate_total - qa.n_gate_failed}/{qa.n_gate_total}")
for c in qa.checks:
    if c.status == "fail":
        say(f"  FAIL {c.check_id} {c.requirement} | {c.observed}")
(OUT / "report.json").write_text(json.dumps(
    {"source_mesh": str(src), "case": CASE, "depths": report,
     "r_after_cap": float(r.max()),
     "qa": {"n_gate_total": qa.n_gate_total, "n_gate_failed": qa.n_gate_failed}},
    indent=1, default=float))
say("done")
