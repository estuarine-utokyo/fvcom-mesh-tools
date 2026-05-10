"""Tests for :func:`fvcom_mesh_tools.mesh_engine.oceanmesh.courant_sizing_function`.

The function applies the linear-wave-theory characteristic-velocity
formula (long-wave celerity + particle velocity, with an overland
surrogate) to convert a depth field into a per-cell maximum element
size such that the approximate Courant number stays below the supplied
target.

These tests stub the DEM interface so they don't require a real
raster file, and verify:

  1. The deep-water branch exactly matches ``(nu * sqrt(g/h) +
     sqrt(g*h)) * dt / C``.
  2. The overland branch clips correctly to ``2 * sqrt(g*nu) * dt / C``.
  3. ``min_edgelength`` / ``max_edge_length`` clamps work in degrees.
  4. Invalid inputs raise ``ValueError`` with a clear message.
"""
from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("oceanmesh")
from fvcom_mesh_tools.mesh_engine.oceanmesh import courant_sizing_function

GRAVITY = 9.81


class _StubDEM:
    """Minimal duck-typed DEM that ``courant_sizing_function`` can
    consume. Provides ``create_grid()``, ``eval()``, and the ``bbox`` /
    ``dx`` / ``dy`` attributes the function reads.
    """

    def __init__(self, depths_2d: np.ndarray, *, bbox, dx, dy):
        self._depths = np.asarray(depths_2d, dtype=float)
        self.bbox = bbox
        self.dx = dx
        self.dy = dy

    def create_grid(self):
        # Return any shape-consistent (lon, lat) the function evaluates;
        # courant_sizing_function only uses .eval() output, not these.
        return np.zeros_like(self._depths), np.zeros_like(self._depths)

    def eval(self, _xy):
        return self._depths


def _expected_h_max_m(depth_m: float, target_C: float, dt: float, nu: float
                     ) -> float:
    h = abs(depth_m) if abs(depth_m) >= 1.0 else 1.0
    if h > nu:
        char_vel = nu * np.sqrt(GRAVITY / h) + np.sqrt(GRAVITY * h)
    else:
        char_vel = 2.0 * np.sqrt(GRAVITY * nu)
    return char_vel * dt / target_C


def test_courant_sizing_deep_water_branch_matches_formula() -> None:
    depths = -np.array([[10.0, 50.0, 100.0]])  # all deep, |h| >> 2 m
    dem = _StubDEM(
        depths,
        bbox=(139.0, 140.0, 35.0, 35.5),  # om CRS=4326 (lon0, lon1, lat0, lat1)
        dx=0.001, dy=0.001,
    )
    grid = courant_sizing_function(
        dem, target_courant=0.7, timestep_s=5.0,
        wave_amplitude_m=2.0, crs=4326,
    )
    # Convert grid values (degrees) back to metres and compare cell-by-cell.
    # oceanmesh's upstream sizing functions feed the mean latitude in
    # degrees directly into np.cos (i.e. they don't deg→rad first).
    # We match that convention so the metre / degree scaling lines up
    # with feature / gradient / wavelength sizing.
    mean_lat = float(np.mean(dem.bbox[2:]))
    meters_per_degree = (
        111132.92
        - 559.82 * np.cos(2 * mean_lat)
        + 1.175 * np.cos(4 * mean_lat)
        - 0.0023 * np.cos(6 * mean_lat)
    )
    h_max_m = grid.values * meters_per_degree
    expected = np.array([
        _expected_h_max_m(d, 0.7, 5.0, 2.0) for d in depths.ravel()
    ]).reshape(depths.shape)
    np.testing.assert_allclose(h_max_m, expected, rtol=1e-9, atol=1e-9)


def test_courant_sizing_overland_branch_uses_amplitude_floor() -> None:
    """Cells with |depth| <= wave_amplitude (or land) clamp to the
    overland surrogate ``2 * sqrt(g * nu) * dt / C``."""
    nu = 2.0
    # Tiny ocean (|h| < 1 → safety-floored to 1) AND truly shallow (|h| ≤ nu).
    depths = -np.array([[0.5, 1.0, 1.5, 2.0]])
    dem = _StubDEM(
        depths, bbox=(139.0, 140.0, 35.0, 35.5), dx=0.001, dy=0.001,
    )
    grid = courant_sizing_function(
        dem, target_courant=0.7, timestep_s=5.0, wave_amplitude_m=nu,
        crs=4326,
    )
    # oceanmesh's upstream sizing functions feed the mean latitude in
    # degrees directly into np.cos (i.e. they don't deg→rad first).
    # We match that convention so the metre / degree scaling lines up
    # with feature / gradient / wavelength sizing.
    mean_lat = float(np.mean(dem.bbox[2:]))
    meters_per_degree = (
        111132.92
        - 559.82 * np.cos(2 * mean_lat)
        + 1.175 * np.cos(4 * mean_lat)
        - 0.0023 * np.cos(6 * mean_lat)
    )
    h_max_m = grid.values * meters_per_degree
    expected = np.array(
        [_expected_h_max_m(d, 0.7, 5.0, nu) for d in depths.ravel()]
    ).reshape(depths.shape)
    np.testing.assert_allclose(h_max_m, expected, rtol=1e-9, atol=1e-9)


def test_courant_sizing_min_edge_length_floor_in_degrees() -> None:
    """``min_edgelength`` clamps the output (in degrees) — useful
    when very deep cells would otherwise produce values larger than
    the target ``hmax_deg``."""
    depths = -np.full((1, 3), 1000.0)  # very deep → very large h_max
    dem = _StubDEM(
        depths, bbox=(139.0, 140.0, 35.0, 35.5), dx=0.001, dy=0.001,
    )
    grid = courant_sizing_function(
        dem, target_courant=0.7, timestep_s=5.0,
        max_edge_length=0.001, crs=4326,
    )
    assert (grid.values <= 0.001 + 1e-12).all()


def test_courant_sizing_invalid_target_courant_raises() -> None:
    dem = _StubDEM(
        -np.ones((1, 1)) * 10.0,
        bbox=(139.0, 140.0, 35.0, 35.5), dx=0.001, dy=0.001,
    )
    with pytest.raises(ValueError, match="target_courant"):
        courant_sizing_function(dem, target_courant=0.0, timestep_s=5.0)


def test_courant_sizing_invalid_timestep_raises() -> None:
    dem = _StubDEM(
        -np.ones((1, 1)) * 10.0,
        bbox=(139.0, 140.0, 35.0, 35.5), dx=0.001, dy=0.001,
    )
    with pytest.raises(ValueError, match="timestep"):
        courant_sizing_function(dem, target_courant=0.7, timestep_s=0.0)


def test_courant_sizing_invalid_wave_amplitude_raises() -> None:
    dem = _StubDEM(
        -np.ones((1, 1)) * 10.0,
        bbox=(139.0, 140.0, 35.0, 35.5), dx=0.001, dy=0.001,
    )
    with pytest.raises(ValueError, match="wave_amplitude"):
        courant_sizing_function(
            dem, target_courant=0.7, timestep_s=5.0, wave_amplitude_m=0.0,
        )
