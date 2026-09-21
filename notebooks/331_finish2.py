# Finishing chain on the UTM sample-repro mesh: the packaged
# fvcom_mesh_tools.algorithms.obc_finish.finish_obc_mesh chain
# (perp-local -> phase_h frozen -> compact -> perp/R4 flips ->
# phase_h -> compact -> C4 flips, OBC-line displacement verified).
import os

from fvcom_mesh_tools.one_wide import finishing_one_wide
from fvcom_mesh_tools.sizing import load_sizing

_recipe = load_sizing(os.environ["SR_SIZING"]) if os.environ.get("SR_SIZING") else None
ONE_WIDE = finishing_one_wide("outputs/sample_repro/channel_policy.json", _recipe)
print(f"[fin] one_wide={ONE_WIDE}", flush=True)
from fvcom_mesh_tools.algorithms.obc_finish import finish_obc_mesh
from fvcom_mesh_tools.channel_policy import resolve_narrow_channels
from fvcom_mesh_tools.io import read_fort14, write_fort14

SRC = "outputs/sample_repro/sample_repro_utm.14"
DST = "outputs/sample_repro/sample_repro_final.14"

mesh = read_fort14(SRC)
# The selected policy controls ACTION; strict detection still reports in both modes.
# Legacy forbid behavior below (owner 2026-07-13: "do not create
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
                                      max_rounds=8, one_wide=ONE_WIDE)
print(f"[fin] channel policy: flagged={cinfo['n_flagged']} "
      f"widened={cinfo['n_widened']} "
      f"deleted={cinfo['n_deleted_elements']}", flush=True)
for cl in cinfo.get("clusters", []):
    print(f"[fin]   cluster n={cl['n_members']} -> {cl['action']} "
          f"(neighbor basins {cl['neighbor_sizes']})", flush=True)
# meshing land in mesh CRS: the widen-then-split choke operator
# needs it for the wall-thickness guard
import json as _json

import geopandas as _gpd
from shapely.ops import unary_union as _uu

_land_utm = _uu(list(_gpd.read_file(
    "outputs/sample_repro/land_channel_adj.shp")
    .to_crs(32654).geometry))
# The stochastic local repair takes a seed; a handful of stubborn
# elements can depend on it, so it is selectable for sweeps.
mesh, info = finish_obc_mesh(
    mesh, seed=int(os.environ.get("SR_FIN_SEED", 42)),
    land_union=_land_utm, one_wide=ONE_WIDE)
_wops = (info.get("choke_widen") or {}).get("ops", [])
_wops += (info.get("choke_widen_2") or {}).get("ops", [])
# HUMAN-JUDGMENT mesh edits (owner 2026-07-15): coordinate-
# addressed, gated, loud-fail ledger applied after the automatic
# chain -- see algorithms/mesh_edits.py
from pathlib import Path as _P

from pyproj import Transformer as _Tr

from fvcom_mesh_tools.algorithms.mesh_edits import apply_mesh_edits

_tr43 = _Tr.from_crs(4326, 32654, always_xy=True)
_eds = []
for _mf in sorted(_P("recipes/mesh_edits/sample_repro").glob("*.json")):
    _md = _json.loads(_mf.read_text())
    for _o in _md.get("ops", []):
        _o.setdefault("id", _md.get("id", _mf.stem))
        _eds.append(_o)
if _eds:
    mesh, _minfo = apply_mesh_edits(
        mesh, _eds, _land_utm,
        lambda lo, la: _tr43.transform(lo, la))
    _wops += _minfo.get("widen_ops", [])
    print(f"[fin] mesh_edits: applied {_minfo['applied']}, "
          f"FAILED {_minfo['failed']}", flush=True)
    for _r in _minfo["results"]:
        print(f"[fin]   {_r.get('id')}: {_r.get('op')} "
              f"nodes {_r.get('nodes_f14')} -> "
              f"{'OK ' + str({k: v for k, v in _r.items() if k in ('new_edge', 'min_angle', 'split_frac', 'local_worst_deg')}) if _r.get('ok') else 'FAILED: ' + str(_r.get('reason'))}",
              flush=True)
with open("outputs/sample_repro/widen_ops.json", "w") as _f:
    _json.dump(_wops, _f)
# COASTLINE FIT (owner 2026-09-21: "the mesh shoreline is visibly off the OSM
# coastline where a human could place it by eye").  DistMesh only projects a
# point back when it has stepped OUTSIDE the domain, so a boundary node that
# stopped just inside is never pulled to the coast; the mesh shoreline then
# zig-zags across the real one (measured: median 10 m, p90 84 m, 32 % of the
# boundary nodes on the land side).  This pass moves those nodes onto the
# polygon the generator was given, under the C1/C2/C4 gates plus a time-step
# floor, and changes no topology.
if os.environ.get("SR_COAST_FIT", "on") != "off":
    import numpy as _np

    from fvcom_mesh_tools.coast_fit import fit_boundary_to_coast

    _obc = (_np.concatenate([_np.asarray(s, int).ravel()
                             for s in mesh.open_boundaries])
            if mesh.open_boundaries else _np.empty(0, int))
    # The width guard measures against the SAME polygon the fit targets, so
    # the two agree about where the water is.  Measuring against raw OSM
    # instead treats every intended widening as forbidden ground and costs
    # most of the gain (p90 offset 14.0 m -> 41.6 m) for no measured benefit:
    # with 364's width read from four samples per element rather than one,
    # both references leave the gate passing.
    _cf = fit_boundary_to_coast(mesh.nodes, mesh.elements, _land_utm,
                                fixed=_obc, depths=mesh.depths)
    mesh.nodes[:, :2] = _cf.nodes[:, :2]
    print(f"[fin] {_cf.summary()}", flush=True)
    with open("outputs/sample_repro/coast_fit.json", "w") as _f:
        _json.dump(_cf.to_dict(), _f, indent=1)
else:
    print("[fin] coast fit: OFF", flush=True)
for k, v in info.items():
    print(f"[fin] {k}: {v}", flush=True)
write_fort14(mesh, DST)
print(f"[fin] wrote {DST} NP={mesh.n_nodes:,} NE={len(mesh.elements):,}",
      flush=True)
