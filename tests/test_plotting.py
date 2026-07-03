"""Tests for the §10 figure system (plotting.py, fmesh-plot-mesh)."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest

from fvcom_mesh_tools.cli import meshplot
from fvcom_mesh_tools.io import Fort14Mesh, write_fort14
from fvcom_mesh_tools.plotting import (
    ReferenceGrid,
    _col_index,
    _col_label,
    plot_mesh_overview,
)


def _mesh() -> Fort14Mesh:
    # 3x3-node right-triangle grid, 12 km spacing -> spans 24 km so a
    # 5 km reference grid has several columns/rows.
    n = 3
    nodes = np.array(
        [[i * 12_000.0, j * 12_000.0] for j in range(n) for i in range(n)],
    )
    elements = []
    for j in range(n - 1):
        for i in range(n - 1):
            a, b = j * n + i, j * n + i + 1
            c, d = (j + 1) * n + i + 1, (j + 1) * n + i
            elements.append([a, b, c])
            elements.append([a, c, d])
    return Fort14Mesh(
        title="plot-test",
        nodes=nodes,
        depths=np.full(n * n, 5.0),
        elements=np.asarray(elements),
        open_boundaries=[np.array([3, 6])],
        land_boundaries=[(20, np.array([6, 7, 8, 5, 2, 1, 0, 3]))],
    )


# ---------------------------------------------------------------------------
# Column labels + grid math
# ---------------------------------------------------------------------------


def test_col_labels_roundtrip():
    for i in (0, 1, 25, 26, 27, 51, 52, 701, 702):
        assert _col_index(_col_label(i)) == i
    assert _col_label(0) == "A"
    assert _col_label(25) == "Z"
    assert _col_label(26) == "AA"
    with pytest.raises(ValueError):
        _col_index("A1")


def test_reference_grid_from_bbox():
    g = ReferenceGrid.from_bbox((1000.0, 2000.0, 26_000.0, 23_000.0),
                                cell_m=5000.0)
    assert g.x0 == 0.0 and g.y0 == 25_000.0
    assert g.n_cols == 6 and g.n_rows == 5
    # A1 = north-west cell.
    assert g.cell_bbox("A1") == (0.0, 20_000.0, 5000.0, 25_000.0)
    # Rows count southward.
    assert g.cell_bbox("A5") == (0.0, 0.0, 5000.0, 5000.0)
    assert g.cell_label(2, 3) == "C4"


def test_parse_ref_cell_range_alias():
    g = ReferenceGrid.from_bbox((0.0, 0.0, 30_000.0, 25_000.0),
                                cell_m=5000.0,
                                aliases={"Banzu": (1.0, 2.0, 3.0, 4.0)})
    c4 = g.parse_ref("C4")
    assert c4 == g.cell_bbox("C4")
    rng = g.parse_ref("C4-D5")
    assert rng[0] == g.cell_bbox("C4")[0]
    assert rng[1] == g.cell_bbox("D5")[1]
    assert rng[2] == g.cell_bbox("D5")[2]
    assert rng[3] == g.cell_bbox("C4")[3]
    assert g.parse_ref("C4–D5") == rng  # en-dash accepted
    assert g.parse_ref("Banzu") == (1.0, 2.0, 3.0, 4.0)
    with pytest.raises(ValueError):
        g.parse_ref("Z99")
    with pytest.raises(ValueError):
        g.parse_ref("4C")


# ---------------------------------------------------------------------------
# Figures (smoke, no xcoast)
# ---------------------------------------------------------------------------


def test_plot_mesh_overview_writes_png(tmp_path):
    p = plot_mesh_overview(_mesh(), tmp_path / "overview.png", coast=None)
    assert p.exists() and p.stat().st_size > 10_000


def test_plot_mesh_overview_zoom_by_ref(tmp_path):
    p = plot_mesh_overview(
        _mesh(), tmp_path / "zoom.png", coast=None, zoom="B2", dpi=80,
    )
    assert p.exists()
    # The zoom canvas must stay at figsize x dpi — off-view grid
    # labels must not inflate the saved area (the 1.3-gigapixel bug).
    from PIL import Image

    w, h = Image.open(p).size
    assert w <= 11 * 80 * 1.2 and h <= 10 * 80 * 1.2
    with pytest.raises(ValueError):
        plot_mesh_overview(
            _mesh(), tmp_path / "bad.png", coast=None, cell_m=None, zoom="B2",
        )


def test_cli_plot(tmp_path, capsys):
    f14 = tmp_path / "m.14"
    write_fort14(_mesh(), f14)
    rc = meshplot.main([
        str(f14), "--zoom", "B2",
        "--alias", "Spot=1000,1000,9000,9000", "--zoom", "Spot",
        "--dpi", "80",
    ])
    assert rc == 0
    assert (tmp_path / "m_overview.png").exists()
    assert (tmp_path / "m_B2.png").exists()
    assert (tmp_path / "m_Spot.png").exists()
    out = capsys.readouterr().out
    assert "overview" in out


def test_cli_plot_bad_alias(tmp_path):
    f14 = tmp_path / "m.14"
    write_fort14(_mesh(), f14)
    assert meshplot.main([str(f14), "--alias", "broken"]) == 2
