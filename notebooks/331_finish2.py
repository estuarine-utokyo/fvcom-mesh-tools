# Finishing chain on the UTM sample-repro mesh: the packaged
# fvcom_mesh_tools.algorithms.obc_finish.finish_obc_mesh chain
# (perp-local -> phase_h frozen -> compact -> perp/R4 flips ->
# phase_h -> compact -> C4 flips, OBC-line displacement verified).
import numpy as np
from fvcom_mesh_tools.algorithms.obc_finish import finish_obc_mesh
from fvcom_mesh_tools.channel_policy import resolve_narrow_channels
from fvcom_mesh_tools.io import read_fort14, write_fort14

SRC = "outputs/sample_repro/sample_repro_utm.14"
DST = "outputs/sample_repro/sample_repro_final.14"

mesh = read_fort14(SRC)
# narrow-channel policy (owner): through/big-port -> widen,
# dead-end/small-port -> delete basin and all
# delete-only here: through-canal widening is done properly by the
# two-pass refinement in the generator (fans violate C1 in narrow
# canals and get undone by phase_h -- churn)
mesh, cinfo = resolve_narrow_channels(mesh, min_basin_elements=6,
                                      apply_widen=False,
                                      small_cluster_delete=0)
print(f"[fin] channel policy: flagged={cinfo['n_flagged']} "
      f"widened={cinfo['n_widened']} "
      f"deleted={cinfo['n_deleted_elements']}", flush=True)
for cl in cinfo.get("clusters", []):
    print(f"[fin]   cluster n={cl['n_members']} -> {cl['action']} "
          f"(neighbor basins {cl['neighbor_sizes']})", flush=True)
mesh, info = finish_obc_mesh(mesh, seed=42)
for k, v in info.items():
    print(f"[fin] {k}: {v}", flush=True)
write_fort14(mesh, DST)
print(f"[fin] wrote {DST} NP={mesh.n_nodes:,} NE={len(mesh.elements):,}",
      flush=True)
