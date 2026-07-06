import numpy as np, shapely
from scipy.spatial import cKDTree
from shapely.geometry import Polygon
from pyproj import Transformer
from fvcom_mesh_tools.io import read_fort14
from oceanmesh.mesh_merge import _fix, _boundary_polygons, _reconstruct_sizing
from oceanmesh.mesh_merge import cat_meshes

m = read_fort14('outputs/pipeline_v6r/tokyo_bay_v6_final.14')
P0, T0 = m.nodes.copy(), m.elements.copy()
zone = np.linalg.norm(P0 - np.array([396750.0, 3950650.0]), axis=1) < 2000
Z = P0[zone]
def disp(nodes, tag):
    dd, _ = cKDTree(nodes).query(Z)
    print(f"[probe] {tag}: zone disp max = {dd.max():.4f}", flush=True)

tr = Transformer.from_crs("EPSG:4326","EPSG:32654",always_xy=True)
ring = np.array([[139.80,35.44],[139.86,35.44],[139.86,35.49],[139.80,35.49],[139.80,35.44]])
x, y = tr.transform(ring[:,0], ring[:,1])
pg = Polygon(np.column_stack([x, y]))
cen = P0[T0].mean(axis=1)
inside = shapely.contains_xy(pg, cen[:,0], cen[:,1])
sub_p, sub_t = _fix(P0, T0[inside])
hole_p, hole_t = _fix(P0, T0[~inside])
disp(hole_p, "after hole _fix")

import oceanmesh as om
new_pt = om.remesh_patch(P0, T0, np.column_stack([x, y])[:-1], target_h=150.0, seed=42)
disp(new_pt[0], "after full remesh_patch")

cp, ct = cat_meshes(sub_p, sub_t, hole_p, hole_t, tol=1e-3)
disp(cp, "after cat(sub, hole) alone")
