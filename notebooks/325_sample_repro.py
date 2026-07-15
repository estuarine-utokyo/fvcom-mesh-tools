# Reproduce the goto2023 SAMPLE (TokyoBay_grd.dat) with the
# faithful OM2D stack: bay-only domain, the sample's OBC arc as a
# CONSTRAINED input line (pfix+egfix), uniform-Courant sizing
# (h = sqrt(gH)*dt/Cr clipped to [hmin, maxel]), OSM coastline.
import os, sys, logging, time
import faulthandler
faulthandler.enable()
faulthandler.dump_traceback_later(600, repeat=True)
import numpy as np
logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                    format="%(levelname)s %(name)s: %(message)s")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
import oceanmesh as om
from oceanmesh import DEM, Region, Shoreline

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
OUT = Path("outputs/sample_repro")
OUT.mkdir(parents=True, exist_ok=True)
DEG = 1.0 / 111e3
t0 = time.time()

# knobs (single-pass targeted tuning against the sample bands)
H0 = float(os.environ.get("SR_H0", 290.0))          # m
MAXEL = float(os.environ.get("SR_MAXEL", 1400.0))    # m
GRADE = float(os.environ.get("SR_GRADE", 0.165))
DT = float(os.environ.get("SR_DT", 15.0))            # s
CRMIN = float(os.environ.get("SR_CRMIN", 0.45))
FS = float(os.environ.get("SR_FS", 3.0))

# The OBC is an INPUT: the goto2023 sample's 13-node smooth arc
# (TokyoBay_obc.dat + TokyoBay_grd.dat, EPSG:32654 -> 4326),
# NW (Miura/Otsu coast) -> SE (mid-channel corner). Mesh nodes are
# CONSTRAINED onto this line via pfix+egfix.
OBC_ARC = np.array([
    [139.6713, 35.1396], [139.6737, 35.1288], [139.6772, 35.1168],
    [139.6816, 35.1031], [139.6871, 35.0877], [139.6946, 35.0705],
    [139.7000, 35.0576], [139.7069, 35.0445], [139.7134, 35.0327],
    [139.7216, 35.0184], [139.7289, 35.0047], [139.7373, 34.9916],
    [139.7497, 34.9750]])
OBC_SEG = np.column_stack([np.arange(12), np.arange(1, 13)])

# bay-only domain POLYGON matching the goto2023 sample geometry:
# the SW crossing follows the OBC arc, extended into Miura land at
# the NW end; the southern closure runs east at lat~34.973 to the
# Boso coast (the sample carries this as an artificial land line).
poly = np.array(
    [[139.83, 34.973], [140.12, 34.973], [140.12, 35.75],
     [139.60, 35.75], [139.60, 35.20], [139.6642, 35.1546]]
    + OBC_ARC.tolist()
    + [[139.83, 34.973]])
bbox = (139.60, 140.12, 34.96, 35.75)
reg = Region(bbox, 4326)

# geometry-stage narrow-channel policy (owner 2026-07-12): decide
# channel fates on the SHORELINE, before meshing -- deterministic
# and minimum-mesh-size preserving. Through / big-basin channels
# get their banks pushed into land (width -> 1.8x h_mesh) so two
# STANDARD-size rows fit; dead-ends and small basins are closed.
import geopandas as gpd
from shapely.geometry import Polygon as _Poly
from shapely.ops import unary_union as _uu
from fvcom_mesh_tools.prep.channel_policy_geom import (
    apply_channel_policy_to_land,
)

# ROLLED BACK to default OFF (owner 2026-07-12): the connectivity
# comparator (342) measured 214 sample-water closures + 342
# on-land elements from this stage -- worse than the ~19 one-wide
# cells it was meant to fix. Re-enable only per-site with the
# comparator as gate (SR_CH_POLICY=on).
CH_SHP = OUT / "land_channel_adj.shp"
_land_g = gpd.read_file("outputs/tb_varres_3r/land_osm_wide.shp")
_dom = _Poly(poly)
_cosw = float(np.cos(np.deg2rad(35.35)))
if os.environ.get("SR_CH_POLICY", "off") != "on":
    chinfo = {"widened": [], "closed": [], "n_narrow": 0}
    _new_land = _uu(list(_land_g.geometry))
    print("[sr] channel policy (geometry stage): OFF (rollback)",
          flush=True)
else:
    _new_land, chinfo = apply_channel_policy_to_land(
        _uu(list(_land_g.geometry)), _dom,
        h_mesh_m=H0 * 1.2, obc_point=tuple(OBC_ARC[6]),
        metric_scale=(111e3 * _cosw, 111e3))
if os.environ.get("SR_CH_POLICY", "off") == "on":
    print(f"[sr] channel policy (geometry stage): "
          f"{len(chinfo['widened'])} widened, "
          f"{len(chinfo['closed'])} closed "
          f"(of {chinfo['n_narrow']} narrow corridors)", flush=True)
for r in chinfo["widened"] + chinfo["closed"]:
    print(f"[sr]   {r['action']}: ({r['center'][0]:.3f}, "
          f"{r['center'][1]:.3f}) area={r['area_cells']:.1f} "
          f"extent={r['extent_cells']:.1f} cells "
          f"members={r['n_members']}", flush=True)

# ARC-BASED CHANNEL EDITS (owner-approved 2026-07-12): each JSON in
# recipes/edits/sample_repro/ is ONE waterway, given as its along-
# channel arc + width: {"id", "note", "arc": [[lon,lat],...],
# "width_m", "min_gap_m"}. Applied in filename order — a
# VERSION-CONTROLLED, deterministic manual-edit log.
# carve_channel_corridor is barrier-safe: it raises instead of
# piercing land that separates other water.
import json as _json
from fvcom_mesh_tools.channel_arcs import carve_channel_corridor
CH_PF, CH_EG = [], []      # bank pfix/egfix accumulated per edit
# A/B lever (explicit opt-in, loud): comma-separated edit STEMS to
# skip -- used to validate that an automatic stage reproduces a
# manual edit before retiring it. Never silent: every skip prints.
_SR_EXCL = {s.strip() for s in os.environ.get(
    "SR_EDITS_EXCLUDE", "").split(",") if s.strip()}
for _ef in sorted(Path("recipes/edits/sample_repro").glob("*.json")):
    if _ef.stem in _SR_EXCL:
        print(f"[sr] channel edit {_ef.stem}: EXCLUDED "
              f"(SR_EDITS_EXCLUDE)", flush=True)
        continue
    _ed = _json.loads(_ef.read_text())
    if _ed.get("type") == "land_patch":
        # data-crack correction: force a region to LAND (e.g. the
        # sliver between the artificial west domain edge and the
        # OSM land data that meshed as 16-deg sliver cells)
        import shapely.geometry as _sg
        _patch = _sg.shape(_ed["geometry"])
        _new_land = _uu([_new_land, _patch])
        print(f"[sr] channel edit {_ed.get('id', _ef.stem)} "
              f"(land_patch): +{_patch.area * (111e3 * _cosw) * 111e3 / 1e4:.1f} ha",
              flush=True)
        continue
    if _ed.get("type") == "water_patch":
        # owner-approved water footprint (e.g. the sample's own
        # meshed water): land inside the polygon becomes water
        import shapely.geometry as _sg
        _patch = _sg.shape(_ed["geometry"])
        _opened = _patch.intersection(_new_land)
        _new_land = _new_land.difference(_patch)
        print(f"[sr] channel edit {_ed.get('id', _ef.stem)} "
              f"(water_patch rev {_ed.get('rev')}): opened "
              f"{_opened.area * (111e3 * _cosw) * 111e3 / 1e4:.1f}"
              f" ha of land", flush=True)
        continue
    _tol = _ed.get("arc_on_land_tol_m")   # explicit opt-in only
    _w = (np.asarray(_ed["widths_m"], float)
          if "widths_m" in _ed else float(_ed["width_m"]))
    _new_land, _einfo = carve_channel_corridor(
        _new_land, np.asarray(_ed["arc"], float), _w,
        min_gap_m=float(_ed.get("min_gap_m", 150.0)),
        metric_scale=(111e3 * _cosw, 111e3), domain_poly=_dom,
        arc_on_land_tol_m=None if _tol is None else float(_tol))
    if _ed.get("bank_pfix"):
        _off = sum(len(a) for a in CH_PF)
        CH_PF.append(np.asarray(_ed["bank_pfix"], float))
        CH_EG.append(np.asarray(_ed["bank_egfix"], int) + _off)
    print(f"[sr] channel edit {_ed.get('id', _ef.stem)}: "
          f"len={_einfo['arc_length_m']:.0f} m "
          f"width_max={_einfo['width_max_m']:.0f} m "
          f"arc_on_land={_einfo['arc_on_land_m']:.0f} m "
          f"land_removed={_einfo['land_removed_m2']/1e4:.2f} ha "
          f"banks={'yes' if _ed.get('bank_pfix') else 'no'}",
          flush=True)

# MATHEMATICAL waterway policy (owner 2026-07-12): detect
# sub-two-row waterways from the OSM geometry ITSELF (no
# reference mesh), then keep/close by the operational rules --
# through & big-port: widen to two standard rows along the arc
# (banks pushed into land, barrier-safe); dead ends & small
# basins: fill. Runs AFTER manual edits (data corrections).
if os.environ.get("SR_WATERWAYS", "on") == "on":
    from fvcom_mesh_tools.waterways import (
        apply_waterway_policy,
        detect_waterways,
    )
    _land_pre_ww = _new_land   # post-edit land: pass-2 base
    _recs = detect_waterways(
        _land_pre_ww, _dom, h_mesh_m=1.2 * H0,
        obc_point=tuple(OBC_ARC[6]),
        metric_scale=(111e3 * _cosw, 111e3))
    # OSM waterway CENTRELINES authorize bridge-gap opening
    # (bridge vs levee is undecidable from the land polygon
    # alone). Source: the archived Geofabrik Kanto dump.
    _dd = os.environ.get("DATA_DIR")
    if not _dd:
        raise RuntimeError(
            "DATA_DIR is not set -- required for the OSM "
            "waterway-centreline layer (geodata/OSM/"
            "geofabrik_kanto/gis_osm_waterways_free_1.shp)")
    _wl = gpd.read_file(
        f"{_dd}/geodata/OSM/geofabrik_kanto/"
        "gis_osm_waterways_free_1.shp",
        bbox=(139.60, 34.96, 140.12, 35.75))
    _wl = _wl[_wl["fclass"].isin(["river", "canal", "stream"])]
    print(f"[sr] OSM waterway centrelines: {len(_wl)} "
          f"(river/canal/stream)", flush=True)
    # SAMPLE-CONSERVATIVE widths (owner verdict 2026-07-13: the
    # v2 escalation -- widen_factor 1.0 + 1.7 h bar -- meshed
    # narrow water the sample rightly leaves out; the sample is
    # the quality bar). The artifact FIXES stay on regardless:
    # verified dup-skip, crumb cleanup, severance override,
    # boundary short-edge collapse.
    _new_land, _winfo = apply_waterway_policy(
        _land_pre_ww, _dom, _recs, h_mesh_m=1.2 * H0,
        metric_scale=(111e3 * _cosw, 111e3),
        h_grade_per_m=1.2 * GRADE,
        open_bridges="auto",
        waterway_lines=list(_wl.geometry),
        widen_factor=float(os.environ.get(
            "SR_WIDEN_FACTOR", "0.875")),
        attain_bar_h=float(os.environ.get(
            "SR_ATTAIN_BAR", "1.5")),
        force_two_rows=(os.environ.get(
            "SR_FORCE2ROWS", "off") == "on"))
    # NORMALIZATION (owner rule 2026-07-15: water we decided not
    # to resolve is LAND for later geometry decisions). Two fixed
    # passes, never a loop: pass 1 above learns which corridors
    # are kept; normalization converts every basin left behind a
    # sub-floor passage into land; pass 2 re-detects and re-carves
    # on the normalized land, so corridors widen SYMMETRICALLY and
    # arcs re-snap onto the real channel axis (the edit_004 +
    # edit_005 mechanism, automated).
    if os.environ.get("SR_NORMALIZE", "off") == "on":
        from fvcom_mesh_tools.waterways import (
            normalize_unresolved_water,
        )
        _nfills, _ninfo = normalize_unresolved_water(
            _new_land, _dom, h_mesh_m=1.2 * H0,
            obc_point=tuple(OBC_ARC[6]),
            metric_scale=(111e3 * _cosw, 111e3),
            keep_tubes=_winfo["refine_arcs"])
        print(f"[sr] normalize: components "
              f"{_ninfo['components_filled']}, basins "
              f"{_ninfo['basin_parts_filled']}, fringes "
              f"{_ninfo['fringes_filled']}, neck stubs left "
              f"{_ninfo['neck_stubs']}, filled "
              f"{_ninfo['area_filled_ha']:.1f} ha", flush=True)
        for _f in _ninfo["fills"]:
            if _f["area_ha"] >= 1.0:
                print(f"[sr]   fill {_f['area_ha']:8.1f} ha at "
                      f"({_f['center'][0]:.4f}, "
                      f"{_f['center'][1]:.4f})  {_f['why']}",
                      flush=True)
        (OUT / "normalize.json").write_text(_json.dumps(
            {k: v for k, v in _ninfo.items()}))
        if _nfills:
            _land_pre2 = _uu([_land_pre_ww, *_nfills])
            _recs = detect_waterways(
                _land_pre2, _dom, h_mesh_m=1.2 * H0,
                obc_point=tuple(OBC_ARC[6]),
                metric_scale=(111e3 * _cosw, 111e3))
            _new_land, _winfo = apply_waterway_policy(
                _land_pre2, _dom, _recs, h_mesh_m=1.2 * H0,
                metric_scale=(111e3 * _cosw, 111e3),
                h_grade_per_m=1.2 * GRADE,
                open_bridges="auto",
                waterway_lines=list(_wl.geometry),
                widen_factor=float(os.environ.get(
                    "SR_WIDEN_FACTOR", "0.875")),
                attain_bar_h=float(os.environ.get(
                    "SR_ATTAIN_BAR", "1.5")),
                force_two_rows=(os.environ.get(
                    "SR_FORCE2ROWS", "off") == "on"))
            print("[sr] normalize: pass 2 (re-detect + symmetric "
                  "re-carve on normalized land) done", flush=True)
    # forced two-row ladder constraints from marginal kept
    # branches join the constrained-node set (same plumbing as
    # the manual-edit bank chains)
    if _winfo.get("band_n"):
        _off = sum(len(a) for a in CH_PF)
        CH_PF.extend(np.asarray(x, float)
                     for x in _winfo["band_pfix"])
        CH_EG.extend(np.asarray(x, int) + _off
                     for x in _winfo["band_egfix"])
        print(f"[sr] two-row ladder constraints: "
              f"+{_winfo['band_n']} pfix nodes", flush=True)
        # persist band geometry: QA violations get correlated
        # against ladder extents when judging the forced rows
        (OUT / "ladder_bands.json").write_text(_json.dumps(
            [np.asarray(x, float).round(6).tolist()
             for x in _winfo["band_pfix"]]))
    print(f"[sr] waterways (OSM-native): kept {_winfo['kept']}, "
          f"closed {_winfo['closed']}, "
          f"ignored {_winfo['ignored']}, "
          f"blocked {len(_winfo['blocked'])} "
          f"(retried ok {_winfo['retried']}, "
          f"marginal-connectivity {_winfo['marginal_kept']}, "
          f"dup-skipped {_winfo['dup_skipped']}, "
          f"crumbs {_winfo.get('crumbs_dropped', 0)}), "
          f"land_removed "
          f"{_winfo['land_removed_m2']/1e4:.1f} ha", flush=True)
    # persist the records so the one-wide checker sweeps EVERY
    # kept waterway's arc (not just manual edits) and blocked
    # ones stay visible in the ledger
    _wjson = []
    for _r in _recs:
        _c = _r["geometry"].representative_point()
        _wjson.append({
            "action": _r["action"], "kind": _r["kind"],
            "extent_cells": _r["extent_cells"],
            "basin_cells": _r["basin_cells"],
            "center": [round(float(_c.x), 5),
                       round(float(_c.y), 5)],
            "retry": _r.get("retry"),
            "reason": _r.get("reason"),
            "arc": (np.asarray(_r["arc"], float).round(6).tolist()
                    if _r.get("arc") is not None else None),
            "widths_m": (np.asarray(_r["width_m"], float)
                         .round(1).tolist()
                         if _r.get("width_m") is not None
                         else None),
            "marginal_branches": _r.get("marginal_branches"),
            # EVERY carved branch arc (not just the longest):
            # side branches carry their own widths, and the
            # one-wide sweep + choke diagnosis must see them
            "arcs_done": ([
                [np.asarray(_a2, float).round(6).tolist(),
                 np.asarray(_w2, float).round(1).tolist()]
                for _a2, _w2 in _r["arcs_done"]]
                if _r.get("arcs_done") else None),
        })
    (OUT / "waterways.json").write_text(_json.dumps(_wjson))
    print(f"[sr] waterways records -> {OUT / 'waterways.json'}",
          flush=True)
    for _r in _recs:
        if _r["action"] == "ignore":
            continue
        _c = _r["geometry"].representative_point()
        print(f"[sr]   {_r['action']:7s} {_r['kind']:8s} "
              f"ext={_r['extent_cells']:6.1f} "
              f"basin={_r['basin_cells']:6.1f} "
              f"at ({_c.x:.4f}, {_c.y:.4f})"
              + (f"  REASON: {_r.get('reason', '')[:90]}"
                 if _r["action"] == "blocked" else ""),
              flush=True)
else:
    print("[sr] waterways (OSM-native): OFF", flush=True)
_geoms = list(_new_land.geoms) if hasattr(_new_land, "geoms")     else [_new_land]
gpd.GeoDataFrame(geometry=_geoms, crs=_land_g.crs).to_file(CH_SHP)

sh = Shoreline(str(CH_SHP), poly, H0 * DEG)
sdf = om.signed_distance_function(sh)
dem = DEM(str(OM2D / "datasets/TokyoBay/dem/SRTM15_kanto_15s.nc"),
          bbox=reg, nc_reader="coords")
print(f"[sr] inputs +{time.time()-t0:.0f}s", flush=True)

f = om.feature_sizing_function(
    sh, sdf, r=FS, max_edge_length=MAXEL * DEG,
    lattice_anchor=(dem.bbox[0], dem.bbox[2]))
# SR_MODE=courant: the sample's measured recipe (anatomy 2026-07-11,
# rev.3 after the E/W coast split):
#   The sample has a TILTED, UNIFORM-BASE mouth zone -- coastal
#   base ~930 (mesh) everywhere south of the line through
#   (139.695, 35.215)-(139.795, 35.155)  [W coast coarse to 35.215,
#   E coast only to ~35.155; same latitude, different sides ->
#   not a latitude function], and ~447 (mesh) in the bay. Inside
#   each zone:  h = min(base + g*d_coast, cap),
#   cap = 1030 (bay) / MAXEL (mouth);  limgrad(g) then smooths the
#   zone seam over ~3 km exactly like the sample's 437-665
#   transition bands. Field = mesh/1.2 (DistMesh calibration).
ZBASE = float(os.environ.get("SR_ZBASE", 630.0))   # mouth-zone base, m
ZW = float(os.environ.get("SR_ZW_LAT", 35.215))    # zone lat at 139.695
ZE = float(os.environ.get("SR_ZE_LAT", 35.17))     # zone lat at 139.795
DMAX = float(os.environ.get("SR_DMAX", 1030.0))    # bay far-field cap, m
if os.environ.get("SR_MODE", "courant") == "courant":
    fdst = om.distance_sizing_function(sh, rate=GRADE,
                                       max_edge_length=None)
    vals = np.ma.filled(np.ma.asarray(fdst.values), MAXEL * DEG)
    d_m = (vals / DEG - H0) / GRADE            # metric coast distance
    lon_g, lat_g = fdst.create_grid()
    lat_line = ZW + (ZE - ZW) / 0.10 * (lon_g - 139.695)
    in_zone = lat_g < lat_line
    base_m = np.where(in_zone, ZBASE, H0)
    cap_m = np.where(in_zone, MAXEL, DMAX)
    h_m = np.minimum(base_m + GRADE * d_m, cap_m)
    fdst.values = h_m * DEG
    fdst.build_interpolant()
    f = fdst
    print(f"[sr] mode=courant (two-zone distance field: bay {H0:.0f}"
          f"/{DMAX:.0f}, mouth {ZBASE:.0f}/{MAXEL:.0f}, seam "
          f"({139.695},{ZW})-({139.795},{ZE}), g={GRADE})", flush=True)
g, dt_eff = om.finalize_sizing(
    [f], dem=dem, shoreline=sh, hmin=H0,
    max_edge_length=MAXEL, gradation=GRADE, courant=None)
print(f"[sr] sizing done +{time.time()-t0:.0f}s", flush=True)

# GRADED, DEPTH-DILATED CFL FLOOR (owner 2026-07-14: the CFL
# condition had become stricter than the sample -- implied dt
# 13.06 vs 15.99 s, tail over DEEP water far from any channel).
# The plain OM2D `dt` limiter (post-limgrad raise) FAILED twice
# in one run (6197174): cells straddling the Uraga slope carry
# deep nodes their lattice-point floor never saw (dt 12.68 s at
# 36->117 m node depths), and the un-graded plateau edge made
# C1x3/C4x2. Fix: floor = dt*sqrt(g*(H_dilated+1))/Cr_max with
# H_dilated = max depth within ~500 m (one cell), then a GRADED
# EXPANSION (outward slope GRADE) before taking the pointwise
# max with the limgrad-feasible field -- the result is
# gradation-feasible by construction, so no re-limgrad and no
# steps. No-op for H <= ~21 m: channels stay untouched.
_cfl_dt = os.environ.get("SR_CFL_DT", "18")
if _cfl_dt != "off":
    from scipy.ndimage import maximum_filter
    _dtt = float(_cfl_dt)
    # Cr 0.75 was tried for realization-undershoot margin and
    # REVERTED (run 6202126: the stronger global coarsening
    # perturbed transitions -- dt 15.26, C1x2, two new chokes).
    # 0.9 + the finish-operator CFL gates is the measured
    # optimum; the residual -0.4 s vs the sample is ONE cell's
    # DistMesh undershoot on a correctly-floored broad area.
    _crm = float(os.environ.get("SR_CFL_CRMAX", "0.9"))
    _lon_gr, _lat_gr = g.create_grid()
    _vals_m = np.asarray(np.ma.filled(np.ma.asarray(g.values),
                                      MAXEL * DEG), float) / DEG
    _zf = np.asarray(dem.eval(np.column_stack(
        [_lon_gr.ravel(), _lat_gr.ravel()])))
    _Hf = np.clip(-_zf.reshape(_lon_gr.shape), 0.0, None)
    # lattice spacings in metres (axis 0 / axis 1)
    _sp0 = float(np.hypot(
        np.diff(_lon_gr, axis=0).mean() * _cosw * 111e3,
        np.diff(_lat_gr, axis=0).mean() * 111e3))
    _sp1 = float(np.hypot(
        np.diff(_lon_gr, axis=1).mean() * _cosw * 111e3,
        np.diff(_lat_gr, axis=1).mean() * 111e3))
    _npx = max(1, int(np.ceil(500.0 / max(min(_sp0, _sp1), 1.0))))
    # OPENING first (run 6197233: the narrow dredged Chiba
    # approach, a few hundred m of 30-60 m water crossing shallow
    # surroundings, raised a thin coarse ridge that left C1/C4
    # slivers at L9-d3 twice): deep STRIPS narrower than ~800 m
    # are excluded from the floor -- their Courant tail stays and
    # is accepted; only BROAD deep water coarsens.
    from scipy.ndimage import minimum_filter
    # 400 -> 250 m (owner 2026-07-15: the 400 m opening also
    # dropped the BROAD Yokohama approach strip -- 284 m cells
    # over 33 m water put the RAW-mesh floor at 15.60 s, below
    # the sample. Only the narrow Chiba strip (the L9-d3 C1/C4
    # ridge, twice) needs excluding; strips wider than ~500 m
    # keep the floor.)
    _kop = max(1, int(np.ceil(250.0 / max(min(_sp0, _sp1), 1.0))))
    _Hop = maximum_filter(
        minimum_filter(_Hf, size=2 * _kop + 1),
        size=2 * _kop + 1)
    _Hd_raw = maximum_filter(_Hf, size=2 * _npx + 1)
    _Hd = maximum_filter(_Hop, size=2 * _npx + 1)
    # SHORE-ADJACENT deep water keeps the full floor (run
    # 6197284: the opening also dropped the Yokosuka naval
    # basins, 41-73 m AGAINST the shore, dt back to 11.5 s).
    # Those are safe to coarsen -- one flank is land, so no thin
    # coarse ridge inside fine water can form (6197233 had zero
    # violations there); only OFFSHORE narrow strips (Chiba
    # approach, shallow on both sides) stay excluded.
    _dsp = (_Hd_raw > 21.0) & (_Hd < _Hd_raw - 1.0)
    if bool(_dsp.any()):
        import shapely
        from shapely.strtree import STRtree as _ST2
        _lnd2 = list(gpd.read_file(CH_SHP).geometry)
        _st2 = _ST2(_lnd2)
        _pi = np.nonzero(_dsp)
        _pp = shapely.points(_lon_gr[_pi], _lat_gr[_pi])
        _ni = _st2.nearest(_pp)
        _dl2 = shapely.distance(
            _pp, np.array(_lnd2, dtype=object)[_ni]) * 111e3
        _nsh = _dl2 < 800.0
        _Hd[_pi[0][_nsh], _pi[1][_nsh]] = \
            _Hd_raw[_pi[0][_nsh], _pi[1][_nsh]]
        print(f"[sr] CFL floor: {int(_nsh.sum())} shore-adjacent "
              f"deep cells kept, {int((~_nsh).sum())} offshore "
              f"strip cells excluded", flush=True)
    _F = np.minimum(_dtt * np.sqrt(9.81 * (_Hd + 1.0)) / _crm,
                    MAXEL)
    # graded expansion: F(x) >= F(y) - GRADE*dist(x,y)
    for _it in range(300):
        _Fn = _F
        _Fn = np.maximum(_Fn, np.pad(
            _F[1:, :], ((0, 1), (0, 0))) - GRADE * _sp0)
        _Fn = np.maximum(_Fn, np.pad(
            _F[:-1, :], ((1, 0), (0, 0))) - GRADE * _sp0)
        _Fn = np.maximum(_Fn, np.pad(
            _F[:, 1:], ((0, 0), (0, 1))) - GRADE * _sp1)
        _Fn = np.maximum(_Fn, np.pad(
            _F[:, :-1], ((0, 0), (1, 0))) - GRADE * _sp1)
        if float(np.max(_Fn - _F)) < 0.5:
            _F = _Fn
            break
        _F = _Fn
    _sel = _F > _vals_m + 0.5
    _nup = int(_sel.sum())
    _vals_m[_sel] = _F[_sel]
    g.values = _vals_m * DEG
    g.build_interpolant()
    print(f"[sr] graded CFL floor: dt={_dtt}s Cr<={_crm}, "
          f"H dilation {_npx} px, raised {_nup} lattice cells "
          f"({_it + 1} sweeps)", flush=True)

# CHANNEL REFINEMENT corridors, v2 (owner 2026-07-13: finer cells
# are allowed wherever they do NOT tighten the CFL condition --
# more nodes at unchanged dt are nearly free). Along every carved
# branch, lower the size field so an INTEGER number of rows fits
# the ACHIEVED width: n = ceil(W/(1.2*H0) - 0.15) rows (min 2,
# max 4), field target W/(1.2*n). v1 lessons (run 6190761):
# (a) a flat H<=30 m CFL cap was violated by deep dredged
# channels (implied dt 12.0 -> 10.1 s) -- the floor now uses the
# REAL DEM depth per station, dt >= 15 s with realization margin;
# (b) the two-row-only trigger (W < 2.2 h) missed the awkward
# 2.5-3 h band where DistMesh oscillates between 2 and 3 rows
# (N7 12-cell cluster) -- the row-aware target covers up to 3.2 h.
# Applied BEFORE the OBC corridor (boundary-priority), with a
# re-gradation pass so transitions honour GRADE.
# OWNER VERDICT 2026-07-13: default OFF. The refined field cannot
# be confined to the kept corridors -- the gradation halo lowers
# neighbouring water too, so sub-cell creeks the policy left
# unresolved suddenly mesh (590 beyond-sample elements in narrow
# original water, run 6191458), and a signed-distance field can
# never represent a land wall thinner than the LOCAL cell, so
# finer cells punch through thin levees. The sample stays better;
# keep the machinery only as a measured negative result.
if (os.environ.get("SR_CH_REFINE", "off") == "on"
        and "_winfo" in dir() and _winfo.get("refine_arcs")):
    _vals = np.asarray(np.ma.filled(np.ma.asarray(g.values),
                                    MAXEL * DEG), dtype=float)
    _vpre = _vals.copy()    # pre-refinement field: the CFL guard
    #                         must never COARSEN the baseline
    _lon_gr, _lat_gr = g.create_grid()
    _DT_FLOOR = 15.0        # global implied dt is 12.0 s; the
    #                         margin absorbs sub-target edges
    _n_low = 0
    for _ra, _rw in _winfo["refine_arcs"]:
        _ra = np.asarray(_ra, float)
        _rw = np.asarray(_rw, float)
        _nrow = np.clip(np.ceil(_rw / (1.2 * H0) - 0.15), 2, 4)
        _tf = np.clip(_rw / (1.2 * _nrow), 0.5 * H0, H0)
        _do = (_tf < H0 - 1.0) & (_rw < 3.2 * 1.2 * H0)
        if not bool(_do.any()):
            continue
        # branch-window subgrid, then per-station corridor discs
        _rmax = 0.75 * float(_rw[_do].max()) / 111e3
        _bb = ((_lon_gr >= _ra[:, 0].min() - 1.5 * _rmax)
               & (_lon_gr <= _ra[:, 0].max() + 1.5 * _rmax)
               & (_lat_gr >= _ra[:, 1].min() - 1.5 * _rmax)
               & (_lat_gr <= _ra[:, 1].max() + 1.5 * _rmax))
        _bi = np.nonzero(_bb)
        if len(_bi[0]) == 0:
            continue
        _blon, _blat = _lon_gr[_bi], _lat_gr[_bi]
        _bval = _vals[_bi]
        # CFL floor from the MAX DEM depth around each disc: a
        # station on a shallow shelf must not refine cells whose
        # triangles will straddle a steep slope into deep water
        # (Uraga west side, dt 8.03 s run 6191155 -- even
        # per-lattice-point depth missed the slope the CELL
        # spans). Conservative by design: a disc that touches
        # deep water simply does not refine.
        _bz = np.asarray(dem.eval(
            np.column_stack([_blon, _blat]))).ravel()
        _bH = np.clip(-_bz, 2.0, None)
        for _k in np.nonzero(_do)[0]:
            _px, _py = _ra[_k]
            _r = 0.75 * _rw[_k] / 111e3
            _dd = np.hypot((_blon - _px) * _cosw, _blat - _py)
            _m = _dd < _r
            _near = _dd < _r + 400.0 / 111e3
            if not bool(_near.any()):
                continue
            _Hd = float(_bH[_near].max())
            _flo = max(0.5 * H0,
                       _DT_FLOOR * np.sqrt(9.81 * _Hd) / 1.2)
            _tv = max(_tf[_k], _flo) * DEG
            _sel = _m & (_bval > _tv)
            _n_low += int(_sel.sum())
            _bval[_sel] = _tv
        _vals[_bi] = _bval
    _n_guard = 0
    if _n_low:
        g.values = _vals
        g = om.enforce_mesh_gradation(g, gradation=GRADE)
        # GLOBAL depth-aware CFL guard (run 6191278: the
        # gradation pass spreads refined values ~1 km outward,
        # and the halo crossed the Uraga slope -- dt 11.48 s at
        # 635 m from the disc). Every lattice cell is floored at
        # its DILATED-depth CFL bound, but never above the
        # PRE-refinement field: the baseline mesh (dt 12.0 s) is
        # left exactly as it was.
        from scipy.ndimage import maximum_filter
        _vg = np.asarray(np.ma.filled(np.ma.asarray(g.values),
                                      MAXEL * DEG), dtype=float)
        _zf = np.asarray(dem.eval(np.column_stack(
            [_lon_gr.ravel(), _lat_gr.ravel()])))
        _Hf = np.clip(-_zf.reshape(_lon_gr.shape), 2.0, None)
        _sp_m = min(
            float(np.abs(np.diff(_lon_gr, axis=0)).max()
                  + np.abs(np.diff(_lon_gr, axis=1)).max())
            * _cosw * 111e3,
            float(np.abs(np.diff(_lat_gr, axis=0)).max()
                  + np.abs(np.diff(_lat_gr, axis=1)).max())
            * 111e3)
        _npx = max(1, int(np.ceil(450.0 / max(_sp_m, 1.0))))
        _Hd = maximum_filter(_Hf, size=2 * _npx + 1)
        _guard = np.minimum(
            np.maximum(0.5 * H0,
                       _DT_FLOOR * np.sqrt(9.81 * _Hd) / 1.2)
            * DEG,
            _vpre)
        _low2 = _vg < _guard
        _n_guard = int(_low2.sum())
        _vg[_low2] = _guard[_low2]
        g.values = _vg
        g.build_interpolant()
    print(f"[sr] channel refinement v2: lowered {_n_low} lattice "
          f"cells over {len(_winfo['refine_arcs'])} branch arcs "
          f"(row-aware n=2-4, DEM-depth CFL floor dt>="
          f"{_DT_FLOOR:.0f} s; depth guard re-raised {_n_guard} "
          f"cells)", flush=True)

# OBC boundary-band (fvcom_mesh_tools.obc_band, default ON):
# smooth inner guide line at k*local-target-size offsets + a
# boundary-priority size corridor along the whole southern
# crossing, applied AFTER limgrad. General rule (sample-calibrated
# K=1.25); no sample-specific numbers remain.
if os.environ.get("SR_OBC_LADDER", "on") == "on":
    from fvcom_mesh_tools.obc_band import (
        apply_corridor, build_obc_band, corridor_targets)
    h_arc_m = np.asarray(g.eval(OBC_ARC)).ravel() / DEG * 1.2
    # end sizes: a field eval AT an arc end that abuts a coast or
    # an artificial closure is contaminated by the coastal halo
    # (TB SE corner: ~820 m where the boundary-adjacent water is
    # 1680 m class). No silent guessing -- override explicitly.
    if os.environ.get("SR_OBC_H0"):
        h_arc_m[0] = float(os.environ["SR_OBC_H0"])
    if os.environ.get("SR_OBC_H1"):
        h_arc_m[-1] = float(os.environ["SR_OBC_H1"])
    band = build_obc_band(OBC_ARC, h_arc_m, k_offset=1.25,
                          skip_ends=1)
    PFIX, SEGS = band["pfix"], band["egfix"]
    print(f"[sr] OBC band: offsets "
          f"{band['offsets_m'].min():.0f}-"
          f"{band['offsets_m'].max():.0f} m (K=1.25 x local size)",
          flush=True)
    closure = np.array([[139.7497, 34.9750], [139.83, 34.973]])
    h_coast_m = float(np.asarray(
        g.eval(np.array([[139.82, 34.974]]))).ravel()[0]) / DEG * 1.2
    # corridor target = the BAND size (K x local), not the ambient
    # local size -- the boundary band must stay one class coarser
    pts_m, tgt_m = corridor_targets(
        OBC_ARC, band["offsets_m"], closure_ll=closure,
        h_closure_end_m=h_coast_m)
    lon_g2, lat_g2 = g.create_grid()
    g.values, n_up = apply_corridor(
        lon_g2, lat_g2, np.asarray(g.values, dtype=float),
        pts_m, tgt_m, grade=GRADE,
        arc_mean_lat=float(OBC_ARC[:, 1].mean()))
    g.build_interpolant()
    print(f"[sr] boundary corridor: raised {n_up} lattice cells "
          f"(post-limgrad, boundary-priority)", flush=True)
    # ladder-band size override: the forced two-row rows only
    # stay clean when the ambient sizing matches their spacing
    # (96 bare pfix nodes -> 27 C1 violations, run 6188830)
    if "_winfo" in dir() and _winfo.get("band_n"):
        _bp = np.vstack(_winfo["band_size_pts"])
        _bt = np.concatenate(_winfo["band_size_tgt"])
        _bpm = np.column_stack([
            _bp[:, 0] * np.cos(np.deg2rad(35.35)) * 111e3,
            _bp[:, 1] * 111e3])
        g.values, _n3 = apply_corridor(
            lon_g2, lat_g2, np.asarray(g.values, dtype=float),
            _bpm, _bt, grade=GRADE, arc_mean_lat=35.35)
        g.build_interpolant()
        print(f"[sr] ladder-band corridor: raised {_n3} lattice "
              f"cells", flush=True)
else:
    PFIX, SEGS = OBC_ARC, OBC_SEG

# channel-bank constraints from the applied edits: the same
# pfix+egfix primitive as the OBC ladder, holding both banks of a
# sub-cell-width channel so the 1-row band is meshed, not bridged
if CH_PF:
    _chp = np.vstack(CH_PF)
    _che = np.vstack(CH_EG) + len(PFIX)
    PFIX = np.vstack([np.asarray(PFIX, dtype=float), _chp])
    SEGS = np.vstack([np.asarray(SEGS, dtype=int), _che])
    print(f"[sr] channel bank constraints: +{len(_chp)} pfix, "
          f"+{len(_che)} egfix", flush=True)

# the built-in msh.clean('default') runs with pfix nodes pinned and
# (fork feature) egfix-carrying faces excluded from the boundary
# deletion loop, so the constrained OBC line survives the clean
from scipy.spatial import cKDTree
import shapely
from fvcom_mesh_tools.algorithms.obc_finish import (
    prune_one_wide_protected,
)

arc_line = shapely.LineString(OBC_ARC)
TOL = 1e-6  # deg (~0.1 m): constrained nodes are exact


def _loop_area(pts, lp):
    x, y = pts[lp, 0], pts[lp, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1))
                       - np.dot(y, np.roll(x, -1)))


def build_mesh(g):
    """One full generate->clean->bc pass on the current sizing."""
    p, t = om.generate_mesh(sdf, g, max_iter=60, seed=0,
                            pfix=PFIX, egfix=SEGS)
    ne0 = len(t)
    p, t = prune_one_wide_protected(p, t, PFIX)
    p, t = om.make_mesh_boundaries_traversable(p, t)
    print(f"[sr] 1-wide pruning (pfix-protected): NE {ne0:,} -> "
          f"{len(t):,}", flush=True)
    print(f"[sr] mesh NP={len(p):,} NE={len(t):,} "
          f"+{time.time()-t0:.0f}s", flush=True)
    d_arc, arc_idx = cKDTree(p).query(OBC_ARC)
    if (d_arc > 1e-8).any():
        bad = np.where(d_arc > 1e-8)[0]
        raise RuntimeError(
            f"OBC arc nodes lost during cleanup: arc rows "
            f"{bad.tolist()} (offsets {d_arc[bad]} deg).")
    b = om.interp_bathymetry(p, t, dem, method="cell-averaging",
                             min_depth=2.0)
    loops = om.boundary_loops(t)
    outer = max(loops, key=lambda lp: abs(_loop_area(p, lp)))
    arc_set = set(arc_idx.tolist())
    if not arc_set.issubset(set(outer.tolist())):
        raise RuntimeError("OBC arc nodes are not all on the outer "
                           "boundary loop -- inspect the mesh.")
    start = int(np.where(outer == arc_idx[0])[0][0])
    ring = np.roll(outer, -start)
    if shapely.distance(shapely.Point(p[ring[1]]), arc_line) > TOL:
        ring = np.roll(ring[::-1], 1)
    stop = int(np.where(ring == arc_idx[-1])[0][0])
    open_str = ring[:stop + 1]
    off = np.array([shapely.distance(shapely.Point(p[v]), arc_line)
                    for v in open_str])
    if (off > TOL).any() or not arc_set.issubset(
            set(open_str.tolist())):
        raise RuntimeError(
            "outer-loop run between the OBC arc ends leaves the "
            f"constrained line (max offset {off.max():.2e} deg).")
    land_str = np.append(ring[stop:], ring[0])
    bc = {"open": [open_str], "land": [land_str],
          "island": [lp for lp in loops if lp is not outer]}
    print(f"[sr] bc: open={len(bc['open'])} land={len(bc['land'])} "
          f"island={len(bc['island'])}", flush=True)
    for kk, s2 in enumerate(bc["open"]):
        print(f"[sr]  open[{kk}]: {len(s2)} nodes at "
              f"({p[s2, 0].min():.3f}-{p[s2, 0].max():.3f}, "
              f"{p[s2, 1].min():.3f}-{p[s2, 1].max():.3f})",
              flush=True)
    return p, t, b, bc


p, t, b, bc = build_mesh(g)

om.write_fort14(str(OUT / "sample_repro.14"), p, t, depth=b,
                boundaries=bc)
np.save(OUT / "p.npy", p); np.save(OUT / "t.npy", t)
# kickoff: FVCOM production mesh is cartesian UTM54N (EPSG:32654)
from pyproj import Transformer as _T
_tr = _T.from_crs("EPSG:4326", "EPSG:32654", always_xy=True)
_xu, _yu = _tr.transform(p[:, 0], p[:, 1])
om.write_fort14(str(OUT / "sample_repro_utm.14"),
                np.column_stack([_xu, _yu]), t, depth=b,
                boundaries=bc)
print("[sr] wrote UTM54N fort.14", flush=True)

# per-node Courant at DT (msh.CalcCFL port) -- the design target is
# Cr <= 0.5 everywhere with dt = DT
cr = om.calc_cfl(p, t, b, dt=DT)
print(f"[sr] Cr(dt={DT}) p50/p90/p99/max = "
      f"{np.percentile(cr, [50, 90, 99]).round(3)} {cr.max():.3f} "
      f"(n>0.5: {(cr > 0.5).sum()}, n>0.6: {(cr > 0.6).sum()})",
      flush=True)

from oceanmesh.fix_mesh import simp_qual
from pyproj import Transformer
tr = Transformer.from_crs("EPSG:4326", "EPSG:32654", always_xy=True)
xx, yy = tr.transform(p[:, 0], p[:, 1])
q = simp_qual(np.column_stack([xx, yy]), t)
print(f"[sr] qual min={q.min():.3f} mean={q.mean():.3f}", flush=True)

# band table vs the sample
e = np.vstack([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]])
e = np.unique(np.sort(e, axis=1), axis=0)
lat = 0.5 * (p[e[:, 0], 1] + p[e[:, 1], 1])
L = np.hypot((p[e[:, 0], 0] - p[e[:, 1], 0])
             * np.cos(np.deg2rad(lat)),
             p[e[:, 0], 1] - p[e[:, 1], 1]) / DEG
print("[sr] bands (sample: 1130/1146/541/547/534/518):", flush=True)
for lo, hi in [(34.95, 35.10), (35.10, 35.20), (35.20, 35.30),
               (35.30, 35.40), (35.40, 35.55), (35.55, 35.70)]:
    s = (lat >= lo) & (lat < hi)
    if s.sum():
        print(f"[sr]  {lo:.2f}-{hi:.2f}: p10/50/90 = "
              f"{np.percentile(L[s], [10, 50, 90]).round(0)}",
              flush=True)
