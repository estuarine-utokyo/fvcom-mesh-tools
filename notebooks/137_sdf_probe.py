import json
import numpy as np
import oceanmesh as om
from oceanmesh import Region
from oceanmesh.signed_distance_function import (
    multiscale_signed_distance_function)
from pyproj import Transformer

DEG = 1.0 / 111194.9266
BBOX = (139.40, 34.90, 140.30, 35.90)
POLY = [[139.58, 35.20], [139.58, 35.74], [140.16, 35.74],
        [140.16, 35.32], [139.90, 35.20]]
prep = 'outputs/pipeline_v6r/prep/'
region = Region((BBOX[0], BBOX[2], BBOX[1], BBOX[3]), 4326)
shore_o = om.Shoreline(prep + 'land_outer.shp', region.bbox, 1000.0 * DEG)
sdf_o = om.signed_distance_function(shore_o)
poly = np.asarray(POLY + [POLY[0]], float)
shore_i = om.Shoreline(prep + 'land_opened.shp', poly, 100.0 * DEG)
sdf_i = om.signed_distance_function(shore_i)
union, nests = multiscale_signed_distance_function([sdf_o, sdf_i])

qa = json.load(open('outputs/pipeline_v6r/qa_report.json'))
lo = next(c for c in qa['checks'] if c['check_id'] == 'land_overlap')
tr = Transformer.from_crs("EPSG:32654", "EPSG:4326", always_xy=True)
pts = np.array([tr.transform(o['x'], o['y']) for o in lo['offenders']])
for name, f in (("sdf_o", sdf_o), ("sdf_i", sdf_i),
                ("nest0", nests[0]), ("nest1", nests[1]),
                ("union", union)):
    v = np.asarray(f.eval(pts))
    print(f"{name}: neg(water)={int((v<0).sum())}/{len(v)} "
          f"min={v.min():.5f} max={v.max():.5f}", flush=True)
cov = np.asarray(nests[0].domain[1].eval(pts))
print("covering(inner) at offenders:", np.round(cov, 5).tolist(), flush=True)
print("manual nest0 = max(sdf_o, -cov):",
      np.round(np.maximum(np.asarray(sdf_o.eval(pts)), -cov), 5).tolist(),
      flush=True)
print("boubox head:", np.round(shore_i.boubox[:8], 3).tolist(), flush=True)
print("done", flush=True)
