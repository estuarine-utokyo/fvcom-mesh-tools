# Stage-2 M3 acceptance demo: Tokyo-port channel refine directive
# (150 m) on the v6 all-pass mesh; directives -> auto heal ->
# detectors + figure.
import os
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
from fvcom_mesh_tools.io import read_fort14, write_fort14
from fvcom_mesh_tools.finishing import (apply_directives,
    detect_violations, execute_patches, plan_patches)

m = read_fort14('outputs/pipeline_v6r/tokyo_bay_v6_final.14')
n0 = len(m.elements)
d = [{"polygon": [[139.80, 35.44], [139.86, 35.44],
                  [139.86, 35.49], [139.80, 35.49]],
      "target_h_m": 150.0}]
m, led = apply_directives(m, d, utm_epsg=32654)
print(led, flush=True)
assert led[0]["outcome"].startswith("applied")
obc = [int(v) for b in m.open_boundaries for v in b]
det = detect_violations(m.nodes, m.elements)
patches = plan_patches(m.nodes, m.elements, det, obc_nodes=obc)
print(f"[demo] auto patches after directive: {len(patches)}", flush=True)
if patches:
    m.nodes, led2 = execute_patches(m.nodes, m.elements, patches,
                                    obc_nodes=obc)
det2 = detect_violations(m.nodes, m.elements)
print("[demo] final:", {c: len(det2[c]["elements"])
      for c in ("c1", "c2", "c4", "c5", "pinch")}, flush=True)
write_fort14(m, 'outputs/pipeline_v6r/tokyo_bay_v6_directive_demo.14')
import matplotlib.pyplot as plt
from pyproj import Transformer
tr = Transformer.from_crs("EPSG:32654", "EPSG:4326", always_xy=True)
lon, lat = tr.transform(m.nodes[:, 0], m.nodes[:, 1])
fig, ax = plt.subplots(figsize=(10, 9))
ax.triplot(lon, lat, m.elements, lw=0.3, color="steelblue")
ax.set_xlim(139.78, 139.88); ax.set_ylim(35.42, 35.51)
ax.set_aspect(1 / np.cos(np.deg2rad(35.46)))
ax.set_title(f"M3 directive demo: Tokyo port 150 m "
             f"(NE {n0:,} -> {len(m.elements):,})")
fig.savefig('outputs/figures/m3_directive_demo.png', dpi=200,
            bbox_inches="tight")
print("[demo] figure saved", flush=True)
