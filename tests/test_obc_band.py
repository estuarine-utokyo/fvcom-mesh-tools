"""Tests for the open-boundary band construction (obc_band)."""

import numpy as np
import pytest

from fvcom_mesh_tools.obc_band import (
    apply_corridor,
    build_obc_band,
    corridor_targets,
)


def _straight_arc(n=7, lat0=35.0, lon0=139.7, dlat=-0.01):
    """North-to-south straight arc; domain to the EAST (left of the
    walk direction for a southward walk is east)."""
    lats = lat0 + dlat * np.arange(n)
    lons = np.full(n, lon0)
    return np.column_stack([lons, lats])


class TestBuildObcBand:
    def test_offsets_follow_local_size(self):
        arc = _straight_arc()
        h = np.full(7, 1000.0)
        band = build_obc_band(arc, h, k_offset=1.25, skip_ends=1)
        inner = band["inner_ll"]
        assert len(inner) == 5
        # offset distance ~ 1250 m east of the arc
        cosw = np.cos(np.deg2rad(arc[:, 1].mean()))
        dx_m = (inner[:, 0] - arc[1:-1, 0]) * cosw * 111e3
        assert np.allclose(dx_m, 1250.0, rtol=0.02)
        # inward = east (left of a southward walk)
        assert (dx_m > 0).all()

    def test_egfix_indexes_pfix(self):
        arc = _straight_arc()
        band = build_obc_band(arc, np.full(7, 800.0))
        pfix, egfix = band["pfix"], band["egfix"]
        assert egfix.max() == len(pfix) - 1
        assert len(egfix) == (7 - 1) + (len(band["inner_ll"]) - 1)

    def test_smoothing_damps_noisy_sizing(self):
        arc = _straight_arc(n=9)
        h = np.full(9, 1000.0)
        h[4] = 3000.0  # single-node spike
        band = build_obc_band(arc, h, smooth_passes=2, taper="local")
        off = band["offsets_m"]
        assert off[4] < 1.25 * 3000.0 * 0.75  # spike damped

    def test_too_short_arc_raises(self):
        with pytest.raises(ValueError, match=">= 3 nodes"):
            build_obc_band(_straight_arc(n=2), np.full(2, 1000.0))

    def test_bad_sizes_raise(self):
        arc = _straight_arc()
        with pytest.raises(ValueError, match="positive"):
            build_obc_band(arc, np.zeros(7))
        with pytest.raises(ValueError, match="shape"):
            build_obc_band(arc, np.full(6, 1000.0))

    def test_skip_ends_too_large_raises(self):
        arc = _straight_arc(n=5)
        with pytest.raises(ValueError, match="skip_ends"):
            build_obc_band(arc, np.full(5, 1000.0), skip_ends=2)

    def test_linear_ends_taper_is_monotone(self):
        arc = _straight_arc(n=9)
        h = np.full(9, 1000.0)
        h[0], h[-1] = 800.0, 1600.0
        h[4] = 3000.0  # interior spike must NOT affect the taper
        band = build_obc_band(arc, h)
        off = band["offsets_m"]
        assert np.isclose(off[0], 1000.0)
        assert np.isclose(off[-1], 2000.0)
        assert (np.diff(off) > 0).all()

    def test_bad_taper_raises(self):
        arc = _straight_arc()
        with pytest.raises(ValueError, match="taper"):
            build_obc_band(arc, np.full(7, 1000.0), taper="cubic")


class TestCorridor:
    def test_targets_are_field_scale(self):
        arc = _straight_arc()
        h = np.full(7, 1200.0)
        pts, tgt = corridor_targets(arc, h, mesh_factor=1.2)
        assert np.allclose(tgt, 1000.0, rtol=0.01)
        assert len(pts) == len(tgt) > 7

    def test_closure_tapers_to_coastal_size(self):
        arc = _straight_arc()
        h = np.full(7, 1200.0)
        closure = np.array([[139.7, 34.94], [139.78, 34.94]])
        pts, tgt = corridor_targets(
            arc, h, closure_ll=closure, h_closure_end_m=600.0)
        assert np.isclose(tgt[-1], 500.0, rtol=0.02)   # 600/1.2
        assert np.isclose(tgt.max(), 1000.0, rtol=0.02)

    def test_closure_without_end_size_raises(self):
        arc = _straight_arc()
        with pytest.raises(ValueError, match="h_closure_end_m"):
            corridor_targets(arc, np.full(7, 1200.0),
                             closure_ll=np.array([[139.7, 34.9],
                                                  [139.8, 34.9]]))

    def test_apply_corridor_boundary_priority(self):
        arc = _straight_arc()
        h = np.full(7, 1200.0)
        pts, tgt = corridor_targets(arc, h)
        lon_g, lat_g = np.meshgrid(
            np.linspace(139.6, 139.8, 40),
            np.linspace(34.9, 35.05, 40), indexing="ij")
        vals = np.full(lon_g.shape, 300.0 / 111e3)  # fine field
        out, n_up = apply_corridor(
            lon_g, lat_g, vals, pts, tgt, grade=0.2,
            arc_mean_lat=float(arc[:, 1].mean()))
        assert n_up > 0
        # cells ON the arc raised to ~1000 m (field scale)
        on_arc = (np.abs(lon_g - 139.7) < 0.002) \
            & (lat_g < 35.0) & (lat_g > 34.95)
        assert np.all(out[on_arc] * 111e3 > 900.0)
        # far cells untouched
        far = lon_g > 139.79
        assert np.allclose(out[far], 300.0 / 111e3)
