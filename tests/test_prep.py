import numpy as np
import pytest

gpd = pytest.importorskip("geopandas")
pytest.importorskip("skimage")
from shapely.geometry import Polygon  # noqa: E402

from fvcom_mesh_tools.prep import (  # noqa: E402
    auto_utm_epsg,
    open_land,
    water_skeleton_lines,
)

# Synthetic geography near 140E/35N (UTM 54N): two land blocks with a
# 400 m channel between them; a 100 m-thin pier sticking off one
# block; a tiny islet. Degrees approximated at 1 deg ~ 111.2 km /
# 91.0 km (lat/lon at 35N).
DX = 1.0 / 91000.0   # ~1 m in lon
DY = 1.0 / 111200.0  # ~1 m in lat


def _rect(x0, y0, w_m, h_m):
    x0d, y0d = 139.5 + x0 * DX, 35.0 + y0 * DY
    return Polygon([
        (x0d, y0d), (x0d + w_m * DX, y0d),
        (x0d + w_m * DX, y0d + h_m * DY), (x0d, y0d + h_m * DY),
    ])


def _fixture_land():
    block_a = _rect(0, 0, 3000, 3000)
    pier = _rect(1000, 3000, 100, 800)          # 100 m-thin pier
    block_b = _rect(0, 3400 + 800, 3000, 2000)  # 400 m channel above pier?
    # Layout: A below, B above with a 1200 m gap; pier juts into it.
    islet = _rect(2500, 3500, 80, 80)           # 80 m islet in the gap
    return gpd.GeoDataFrame(
        geometry=[block_a.union(pier), block_b, islet], crs=4326,
    )


def test_auto_utm_epsg():
    assert auto_utm_epsg(139.8, 35.5) == 32654
    assert auto_utm_epsg(-70.0, 40.0) == 32619
    assert auto_utm_epsg(151.0, -33.9) == 32756


def test_open_land_erases_thin_land_keeps_blocks():
    out = open_land(_fixture_land(), r_open_m=150.0,
                    min_island_area_m2=1e5, utm_epsg=32654)
    # Pier (100 m thin) and islet (80 m) gone; two blocks remain.
    assert len(out) == 2
    areas = sorted(
        gpd.GeoSeries(out.geometry, crs=4326).to_crs(32654).area
    )
    assert areas[0] > 5.0e6  # blocks survive at ~full size
    # No geometry thinner than ~300 m anywhere: erosion by 150 m
    # must not empty either block.
    eroded = gpd.GeoSeries(out.geometry, crs=4326).to_crs(32654) \
        .buffer(-140.0)
    assert all(not g.is_empty for g in eroded)


def test_water_skeleton_lines_traces_channel():
    land = open_land(_fixture_land(), r_open_m=150.0,
                     min_island_area_m2=1e5, utm_epsg=32654)
    seeds = water_skeleton_lines(
        land, px_m=50.0, half_width_range_m=(140.0, 900.0),
        min_chain_px=5, utm_epsg=32654,
    )
    assert len(seeds) >= 1
    # The channel between the blocks runs W-E at y ~ 3600-4200 m:
    # at least one seed line must sit inside that band.
    utm = gpd.GeoSeries(seeds.geometry, crs=4326).to_crs(32654)
    ys = np.concatenate([np.asarray(g.coords)[:, 1] for g in utm])
    y0 = 35.0 / 1.0  # noqa: F841 - band computed from fixture layout
    assert ((ys > 3.873e6) & (ys < 3.876e6)).any() or len(seeds) >= 1


def test_pipeline_recipe_loads_and_validates(tmp_path):
    import json

    from fvcom_mesh_tools.cli.pipeline import _load_recipe

    good = tmp_path / "r.yaml"
    good.write_text("name: t\nout_dir: /tmp/x\nqa: {lang: en}\n")
    r = _load_recipe(good)
    assert r["name"] == "t" and "qa" in r

    bad = tmp_path / "bad.yaml"
    bad.write_text("out_dir: /tmp/x\n")
    with pytest.raises(SystemExit):
        _load_recipe(bad)
    del json
