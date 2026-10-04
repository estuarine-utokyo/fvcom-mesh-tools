"""Open-boundary design: coast-normal ends, rounded corners, spacing."""

import numpy as np
import pytest
import shapely

from fvcom_mesh_tools.obc_design import (
    bearing_vector,
    coast_normal,
    fillet,
    ray_intersection,
    resample,
)


def test_bearing_vector_is_compass():
    assert np.allclose(bearing_vector(0), [0, 1])
    assert np.allclose(bearing_vector(90), [1, 0])


def test_coast_normal_points_to_sea_and_lands_on_the_coast():
    land = shapely.box(-50_000, 0, 50_000, 40_000)          # sea to the south
    b, q = coast_normal(land, 1_000, -2_000)
    assert b == pytest.approx(180.0) and np.allclose(q, [1_000, 0])
    land = shapely.box(0, -50_000, 40_000, 50_000)           # sea to the west
    b, _ = coast_normal(land, -500, 0)
    assert b == pytest.approx(270.0)


def test_ray_intersection_and_its_refusals():
    p = ray_intersection([0, 0], [1, 0], [10, -5], [0, 1])
    assert np.allclose(p, [10, 0])
    with pytest.raises(ValueError, match="parallel"):
        ray_intersection([0, 0], [1, 0], [0, 1], [1, 0])
    with pytest.raises(ValueError, match="behind"):
        ray_intersection([0, 0], [1, 0], [-10, -5], [0, 1])


def test_fillet_keeps_the_sides_straight_and_the_arc_tangent():
    v = np.array([[0, 0], [100_000, 0], [100_000, 100_000]], float)
    line = fillet(v, [20_000])
    assert np.allclose(line[0], v[0]) and np.allclose(line[-1], v[-1])
    # the arc starts 20 km before the corner and ends 20 km after it
    assert np.allclose(line[1], [80_000, 0], atol=1e-6)
    assert np.allclose(line[-2], [100_000, 20_000], atol=1e-6)
    # every arc point is 20 km from the centre (80 km, 20 km)
    arc = line[1:-1]
    assert np.allclose(np.hypot(arc[:, 0] - 80_000, arc[:, 1] - 20_000), 20_000)


def test_fillet_refuses_a_radius_that_does_not_fit():
    v = np.array([[0, 0], [10_000, 0], [10_000, 10_000]], float)
    with pytest.raises(ValueError, match="does not fit"):
        fillet(v, [20_000])
    with pytest.raises(ValueError, match="radii"):
        fillet(v, [])


def test_resample_never_goes_below_the_spacing():
    """Review F5: the remainder is merged into the last step, dropping preceding
    nodes where needed so that no step falls below the spacing."""
    line = np.array([[0, 0], [10_000, 0]], float)
    n = resample(line, 3_000)
    d = np.diff(n[:, 0])
    assert np.allclose(n[0], [0, 0]) and np.allclose(n[-1], [10_000, 0])
    assert d.min() >= 3_000 - 1e-6 and d.tolist() == pytest.approx([3_000, 3_000, 4_000])
    n = resample(np.array([[0, 0], [8_000, 0]], float), 3_000)      # was 3000, 3000, 2000
    assert np.diff(n[:, 0]).min() >= 3_000 - 1e-6
    # variable spacing: each step is at least the spacing at both of its ends
    f = lambda p: np.where(p[:, 0] < 5_000, 1_000.0, 2_500.0)       # noqa: E731
    v = resample(line, f)
    d = np.diff(v[:, 0])
    ends = np.maximum(f(v[:-1]), f(v[1:]))
    assert np.all(d >= ends - 1e-6)
    with pytest.raises(ValueError, match="positive"):
        resample(line, 0)
    with pytest.raises(ValueError, match="shorter than its spacing"):
        resample(line, 20_000)


def test_resample_refuses_non_finite_spacing_and_coordinates():
    """Review F19: NaN spacing used to loop for ever."""
    line = np.array([[0, 0], [10_000, 0]], float)
    with pytest.raises(ValueError, match="finite"):
        resample(line, np.nan)
    with pytest.raises(ValueError, match="non-finite"):
        resample(np.array([[0, 0], [np.nan, 0]], float), 1_000)


def test_coast_normal_does_not_cross_a_thin_strip_of_land():
    """Review F7: a far test point beyond a 1 km strip used to pass."""
    strip = shapely.box(-50_000, 0, 50_000, 1_000)           # sea on both sides
    for land in (strip, shapely.Polygon(list(strip.exterior.coords)[::-1])):
        b, q = coast_normal(land, 0, -10)
        assert b == pytest.approx(180.0)                     # away from the strip


def test_resample_floor_holds_on_chords_and_both_ends():
    """Round 2 F4: variable spacing and bends; F12: a step must advance."""
    line = np.array([[0, 0], [10_900, 0]], float)
    f = lambda p: np.where(p[:, 0] < 5_200, 1_000.0, 2_500.0)       # noqa: E731
    v = resample(line, f)
    d = np.linalg.norm(np.diff(v, axis=0), axis=1)
    assert np.all(d >= np.maximum(f(v[:-1]), f(v[1:])) - 1e-6)
    th = np.linspace(0, np.pi / 2, 200)
    arc = np.column_stack([20_000 * np.cos(th), 20_000 * np.sin(th)])
    v = resample(arc, 3_000)
    assert np.linalg.norm(np.diff(v, axis=0), axis=1).min() >= 3_000 * (1 - 1e-9)


def test_fillet_refuses_bad_radii():
    """Round 2 F20: a negative radius made a path outside the corner."""
    v = np.array([[0, 0], [10_000, 0], [10_000, 10_000]], float)
    for r in (-1_000, 0, float("nan")):
        with pytest.raises(ValueError, match="finite and positive"):
            fillet(v, [r])


def test_coast_normal_does_not_jump_a_50_m_strip():
    """Round 2 F5: probes 60 m apart jumped a 50 m strip."""
    strip = shapely.box(-50_000, 0, 50_000, 50)
    for land in (strip, shapely.Polygon(list(strip.exterior.coords)[::-1])):
        b, _ = coast_normal(land, 0, -10)
        assert b == pytest.approx(180.0)


def test_fillet_refuses_a_reversal_and_passes_a_straight_continuation():
    """Round 3 F13: a 180-degree reversal came back unchanged."""
    with pytest.raises(ValueError, match="reverses"):
        fillet(np.array([[0, 0], [10_000, 0], [0, 0]], float), [1_000])
    line = fillet(np.array([[0, 0], [10_000, 0], [20_000, 0]], float), [1_000])
    assert np.allclose(line[:, 1], 0)


def test_a_masked_line_is_refused():
    """Review of the extend tools, round 25 F4."""
    import numpy as np
    import pytest

    from fvcom_mesh_tools.obc_design import resample

    line = np.ma.array([[0.0, 0.0], [10000.0, 0.0]], mask=[[False, False], [True, False]])
    with pytest.raises(ValueError, match="masked"):
        resample(line, 1000.0)


def test_ray_intersection_needs_known_finite_rays():
    """Review of the extend tools, round 26 F6."""
    import numpy as np
    import pytest

    from fvcom_mesh_tools.obc_design import ray_intersection

    with pytest.raises(ValueError):
        ray_intersection(np.ma.array([0.0, 0.0], mask=[True, False]), [1, 0], [1, 1], [0, -1])
    with pytest.raises(ValueError):
        ray_intersection([np.nan, 0], [1, 0], [1, 1], [0, -1])
    assert np.allclose(ray_intersection([0, 0], [1, 0], [1, 1], [0, -1]), [1, 0])


def test_ray_intersection_does_not_depend_on_direction_length():
    """Review of the extend tools, round 27 F5."""
    import numpy as np
    import pytest

    from fvcom_mesh_tools.obc_design import ray_intersection

    assert np.allclose(ray_intersection([0, 0], [1e-7, 0], [1, 1], [0, -1e-7]), [1, 0])
    with pytest.raises(ValueError, match="not finite"):
        ray_intersection([1.6e308, 0], [1, 1], [1.6e308, 1e308], [1, -1])


def test_ray_intersection_does_not_depend_on_the_direction_magnitude():
    # review round 28 F4: the norm overflowed at 1e200 and underflowed at 1e-200
    for scale in (1.0, 1e200, 1e-200):
        pt = ray_intersection([0, 0], [scale, 0], [1, 1], [0, -scale])
        assert pt.tolist() == pytest.approx([1.0, 0.0])
