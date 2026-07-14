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
# narrow-channel policy (owner 2026-07-13: "do not create
# one-mesh-wide channels"): STRICT mode -- the ledger's confirmed
# one-wide criterion joins the w/h flag, throats into small
# appendixes are pruned WITH the appendix (measured census run
# 6195614: 8/14 chokes guarded pockets of 1-23 elements; the
# sample meshes none of them), iterated until no tail is left.
# Loop/through chokes are NEVER deleted (severance/detour class,
# Keihin precedent) -- they stay in the one-wide ledger.
# min_basin_elements=25 deliberately outranks the pre-mesh keep
# bar (min_basin_cells=6 ~ 12 elements): a pocket only reachable
# through a one-wide neck is not usable water at this h.
mesh, cinfo = resolve_narrow_channels(mesh, min_basin_elements=25,
                                      apply_widen=False,
                                      small_cluster_delete=0,
                                      strict_boundary_flag=True,
                                      max_rounds=8)
print(f"[fin] channel policy: flagged={cinfo['n_flagged']} "
      f"widened={cinfo['n_widened']} "
      f"deleted={cinfo['n_deleted_elements']}", flush=True)
for cl in cinfo.get("clusters", []):
    print(f"[fin]   cluster n={cl['n_members']} -> {cl['action']} "
          f"(neighbor basins {cl['neighbor_sizes']})", flush=True)
# meshing land in mesh CRS: the widen-then-split choke operator
# needs it for the wall-thickness guard
import geopandas as _gpd
import json as _json
from shapely.ops import unary_union as _uu

_land_utm = _uu(list(_gpd.read_file(
    "outputs/sample_repro/land_channel_adj.shp")
    .to_crs(32654).geometry))
mesh, info = finish_obc_mesh(mesh, seed=42, land_union=_land_utm)
_wops = (info.get("choke_widen") or {}).get("ops", [])
with open("outputs/sample_repro/widen_ops.json", "w") as _f:
    _json.dump(_wops, _f)
for k, v in info.items():
    print(f"[fin] {k}: {v}", flush=True)
write_fort14(mesh, DST)
print(f"[fin] wrote {DST} NP={mesh.n_nodes:,} NE={len(mesh.elements):,}",
      flush=True)
