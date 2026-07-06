import numpy as np
from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.autofinish import apply_directives, detect_violations

m = read_fort14('outputs/pipeline_v6r/tokyo_bay_v6_final.14')
base_nodes = m.nodes.copy()
d = [{"polygon": [[139.80, 35.44], [139.86, 35.44],
                  [139.86, 35.49], [139.80, 35.49]],
      "target_h_m": 150.0}]
m, led = apply_directives(m, d, utm_epsg=32654)
print(led, flush=True)
det = detect_violations(m.nodes, m.elements)
print("[iso] post-directive detectors:",
      {c: len(det[c]["elements"]) for c in ("c1","c2","c4","c5")},
      flush=True)
cen = m.nodes[m.elements].mean(axis=1)
for e in det["c4"]["elements"]:
    print("  c4 at", cen[e].round(0), flush=True)
c = np.array([396750.0, 3950650.0])
from scipy.spatial import cKDTree
selb = np.linalg.norm(base_nodes - c, axis=1) < 2000
dd, _ = cKDTree(m.nodes).query(base_nodes[selb])
print("[iso] Tokyo-port zone displacement max:", round(float(dd.max()), 3),
      flush=True)
