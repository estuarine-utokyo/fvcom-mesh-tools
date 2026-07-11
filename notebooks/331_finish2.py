# Finishing chain on the UTM sample-repro mesh: the packaged
# fvcom_mesh_tools.algorithms.obc_finish.finish_obc_mesh chain
# (perp-local -> phase_h frozen -> compact -> perp/R4 flips ->
# phase_h -> compact -> C4 flips, OBC-line displacement verified).
import numpy as np
from fvcom_mesh_tools.algorithms.obc_finish import finish_obc_mesh
from fvcom_mesh_tools.io import read_fort14, write_fort14

SRC = "outputs/sample_repro/sample_repro_utm.14"
DST = "outputs/sample_repro/sample_repro_final.14"

mesh = read_fort14(SRC)
mesh, info = finish_obc_mesh(mesh, seed=42)
for k, v in info.items():
    print(f"[fin] {k}: {v}", flush=True)
write_fort14(mesh, DST)
print(f"[fin] wrote {DST} NP={mesh.n_nodes:,} NE={len(mesh.elements):,}",
      flush=True)
