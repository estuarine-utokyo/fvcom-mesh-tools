"""The depth ladder, against the products themselves.

These need `DATA_DIR` and the three grids, so they skip where those are not
mounted. They are worth having anyway: the numbers in
`docs/refine_coast_and_bathy_design.md` §2.2 are the reason the design reports
sounding distance instead of grid coverage, and a design whose evidence no
test re-measures is a design nobody can check later.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

from fvcom_mesh_tools.dem import tokyo_bay as tb

pytestmark = pytest.mark.skipif(
    not (os.environ.get("DATA_DIR")
         and (Path(os.environ["DATA_DIR"]) / tb._PRODUCTS["m7001"][0]).exists()),
    reason="DATA_DIR and the Tokyo Bay depth products are not mounted")

#: The Futtsu nori ground as `futtsu_nori.yaml` declares it.
FUTTSU = (139.7881, 35.3228)

#: The declared polygons, which is what the design measured.  A disc dropped
#: near a fishery is NOT the same thing: a first version of this test put one
#: at (139.9, 35.40) for Banzu and found M7001 missing over much of it, which
#: says where that guess landed and nothing about the fishery.
GEOMETRY = Path(__file__).resolve().parents[1] / "recipes" / "refine" / "geometry"


def _disc(lon, lat, r_m=300.0, n=40):
    k = 111_000.0 * float(np.cos(np.radians(lat)))
    a = np.linspace(0, 2 * np.pi, n, endpoint=False)
    rr = np.linspace(0.1, 1.0, 6)[:, None]
    return (lon + (rr * np.cos(a) * r_m / k).ravel(),
            lat + (rr * np.sin(a) * r_m / 111_000.0).ravel())


def _polygon_samples(path: Path, name: str, k: int = 40):
    import json

    import shapely

    for f in json.loads(path.read_text())["features"]:
        if f["properties"].get("NAME") == name:
            g = shapely.geometry.shape(f["geometry"])
            break
    else:                                                   # pragma: no cover
        pytest.skip(f"{name} is not in {path.name}")
    w, s, e, n = g.bounds
    gx, gy = np.meshgrid(np.linspace(w, e, k), np.linspace(s, n, k))
    m = shapely.contains(g, shapely.points(gx.ravel(), gy.ravel()))
    return gx.ravel()[m], gy.ravel()[m]


def test_m7001_answers_inside_the_fisheries_and_the_rung_says_so():
    """The fallback does not engage here, which is the point of §2.2.

    A finite grid value is not a sounding: the product interpolates across its
    own gaps, so "M7001 covers this" is not evidence of a survey nearby.
    """
    cases = [_disc(*FUTTSU),
             _polygon_samples(GEOMETRY / "banzu_nori_demo.geojson", "banzu_nori"),
             _polygon_samples(GEOMETRY / "futtsu_two_beds_demo.geojson",
                              "futtsu_nori_east")]
    for x, y in cases:
        depth, rung, dist = tb.sample(x, y)
        assert np.isfinite(depth).all()
        assert (rung == 0).mean() > 0.98, "M7001 should answer almost everywhere"
        assert (dist == 0).all(), "nothing here should have been extrapolated"


def test_the_ladder_falls_through_and_then_extrapolates():
    """Far out in the Pacific, M7001 and the 30 m grid have nothing."""
    x, y = np.array([142.5, 143.0]), np.array([33.0, 32.5])
    depth, rung, dist = tb.sample(x, y)
    assert np.isfinite(depth).all()
    assert (rung != 0).all(), "the M7001 rectangle does not reach here"
    assert set(np.unique(rung).tolist()) <= {2, tb.EXTRAPOLATED}


def test_refusing_to_extrapolate_is_offered_and_says_so():
    x, y = np.array([100.0, 101.0]), np.array([10.0, 11.0])   # the Indian Ocean
    depth, rung, dist = tb.sample(x, y, extrapolate=False)
    assert not np.isfinite(depth).any()
    assert (rung == tb.EXTRAPOLATED).all()
    depth2, _, dist2 = tb.sample(x, y, extrapolate=True)
    assert np.isfinite(depth2).all() and (dist2 > 0).all(), (
        "an extrapolated depth has to carry the distance it was invented over")


def test_the_provenance_report_counts_what_it_says():
    x, y = _disc(*FUTTSU)
    _, rung, dist = tb.sample(x, y)
    rep = tb.provenance_report(rung, dist)
    assert rep["n_nodes"] == len(x)
    assert sum(rep[f"n_from_{n}"] for n in tb.LADDER) + rep["n_extrapolated"] \
        == rep["n_nodes"]


def test_a_sounding_is_further_away_than_the_grid_is_fine():
    """The measurement the design rests on, re-made.

    Inside these fisheries the 181 m M7001 grid is complete while its
    soundings are 128-188 m apart, so the depth a 30 m node receives is the
    product's own interpolation. The report is the distance, not a verdict.
    """
    pytest.importorskip("pandas")
    x, y = _disc(*FUTTSU)
    d = tb.sounding_distance(x, y)
    assert np.isfinite(d).all()
    assert np.median(d) > 30.0, (
        "if the soundings were denser than the target, the design's central "
        "caveat would be wrong and should be rewritten rather than asserted")


def test_non_finite_coordinates_are_refused():
    with pytest.raises(ValueError, match="non-finite"):
        tb.sample([np.nan], [35.0])


def test_empty_input_is_empty_output():
    depth, rung, dist = tb.sample([], [])
    assert depth.size == rung.size == dist.size == 0
