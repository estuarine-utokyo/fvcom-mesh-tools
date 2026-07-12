"""Mesh figures with land rendering and a labeled reference grid.

Implements the kickoff §10 visualization contract:

* **Whole-domain figure** — the triangulation with land rendered by
  ``xcoast`` behind it (optional; degrades gracefully when xcoast or
  its cached data are unavailable) plus a **rectangular reference
  grid** in the mesh's projected CRS (true squares, default 5 km).
  Columns are labeled A, B, C, … (west → east; spreadsheet-style AA
  after Z), rows 1, 2, 3, … (north → south) — a cell is e.g. "C4".
* **Regional zoom panels** per grid cell or named-region alias.
* **Prose addressing** — :meth:`ReferenceGrid.parse_ref` converts
  "C4", a range "C4-D5", or a registered alias into a CRS bbox, so
  chat instructions like "raise resolution in C4" map to polygons.

matplotlib is imported lazily; callers that only need the grid math
never pull it in. Set ``MPLBACKEND=Agg`` for headless use.
"""

from __future__ import annotations

import math
import re
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from fvcom_mesh_tools.io import Fort14Mesh

MESH_PNG_DPI: int = 600
"""Default raster DPI for ``*.png`` outputs that visualise a mesh.

Mesh triangulations stay readable when zoomed in only at high DPI; 600 dpi
is the project-wide default for any figure showing the triangulation
itself (``triplot`` / ``tripcolor`` plots, boundary maps, side-by-side
mesh comparisons). Histograms and other non-spatial plots can keep the
matplotlib default. Pass this value explicitly to
``fig.savefig(..., dpi=MESH_PNG_DPI)`` so that the choice is visible at
the call site and consistent across notebooks and scripts.
"""

_REF_RE = re.compile(r"^([A-Z]+)([1-9][0-9]*)$")


def use_readable_style() -> None:
    """Project-wide figure style: text must stay READABLE in the
    saved PNG (owner 2026-07-12: default matplotlib sizes were
    illegible in report figures). Call this once at the top of
    every figure script, BEFORE creating figures, and do not pass
    smaller explicit ``fontsize=`` values at call sites.
    """
    import matplotlib as mpl

    mpl.rcParams.update({
        "font.size": 14,
        "axes.titlesize": 15,
        "axes.labelsize": 13,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "legend.fontsize": 12,
        "figure.titlesize": 18,
    })


def _col_label(i: int) -> str:
    """0 -> A, 25 -> Z, 26 -> AA (spreadsheet style)."""
    if i < 0:
        raise ValueError(f"column index must be >= 0, got {i}")
    out = ""
    i += 1
    while i:
        i, rem = divmod(i - 1, 26)
        out = chr(ord("A") + rem) + out
    return out


def _col_index(label: str) -> int:
    """A -> 0, Z -> 25, AA -> 26."""
    if not label or not label.isalpha():
        raise ValueError(f"invalid column label {label!r}")
    i = 0
    for ch in label.upper():
        i = i * 26 + (ord(ch) - ord("A") + 1)
    return i - 1


@dataclass
class ReferenceGrid:
    """A labeled square reference grid in projected coordinates.

    The origin ``(x0, y0)`` is the north-west corner; columns advance
    east, rows advance SOUTH (row 1 is the northernmost band, matching
    the kickoff's map-reading convention).
    """

    x0: float
    y0: float
    cell_m: float
    n_cols: int
    n_rows: int
    aliases: dict[str, tuple[float, float, float, float]] = field(default_factory=dict)

    @classmethod
    def from_bbox(
        cls,
        bbox: tuple[float, float, float, float],
        *,
        cell_m: float = 5000.0,
        aliases: dict[str, tuple[float, float, float, float]] | None = None,
    ) -> ReferenceGrid:
        """Grid covering ``bbox = (xmin, ymin, xmax, ymax)``, origin
        snapped outward to whole-cell multiples for clean coordinates."""
        xmin, ymin, xmax, ymax = bbox
        if not (xmax > xmin and ymax > ymin and cell_m > 0):
            raise ValueError("bbox must be non-degenerate and cell_m > 0")
        x0 = math.floor(xmin / cell_m) * cell_m
        y0 = math.ceil(ymax / cell_m) * cell_m
        return cls(
            x0=x0,
            y0=y0,
            cell_m=float(cell_m),
            n_cols=int(math.ceil((xmax - x0) / cell_m)),
            n_rows=int(math.ceil((y0 - ymin) / cell_m)),
            aliases=dict(aliases or {}),
        )

    def cell_label(self, i_col: int, j_row: int) -> str:
        return f"{_col_label(i_col)}{j_row + 1}"

    def cell_bbox(self, label: str) -> tuple[float, float, float, float]:
        """Bbox of one cell like ``"C4"`` (case-insensitive)."""
        m = _REF_RE.match(label.strip().upper())
        if not m:
            raise ValueError(f"invalid cell reference {label!r}")
        i = _col_index(m.group(1))
        j = int(m.group(2)) - 1
        if not (0 <= i < self.n_cols and 0 <= j < self.n_rows):
            raise ValueError(
                f"cell {label!r} outside grid "
                f"(cols A-{_col_label(self.n_cols - 1)}, rows 1-{self.n_rows})"
            )
        xmin = self.x0 + i * self.cell_m
        ymax = self.y0 - j * self.cell_m
        return (xmin, ymax - self.cell_m, xmin + self.cell_m, ymax)

    def parse_ref(self, ref: str) -> tuple[float, float, float, float]:
        """Resolve ``"C4"``, a range ``"C4-D5"`` (also en-dash), or a
        registered alias into a bbox (union for ranges)."""
        key = ref.strip()
        if key in self.aliases:
            return self.aliases[key]
        parts = re.split(r"[-–]", key.upper())
        if len(parts) == 1:
            return self.cell_bbox(parts[0])
        if len(parts) == 2:
            a = self.cell_bbox(parts[0])
            b = self.cell_bbox(parts[1])
            return (min(a[0], b[0]), min(a[1], b[1]),
                    max(a[2], b[2]), max(a[3], b[3]))
        raise ValueError(f"invalid region reference {ref!r}")

    def iter_cells(self):
        for j in range(self.n_rows):
            for i in range(self.n_cols):
                yield self.cell_label(i, j), self.cell_bbox(self.cell_label(i, j))


def _add_reference_grid(ax, grid: ReferenceGrid, *, zorder: float = 3.0) -> None:
    x1 = grid.x0 + grid.n_cols * grid.cell_m
    y1 = grid.y0 - grid.n_rows * grid.cell_m
    style = {"color": "0.45", "lw": 0.5, "zorder": zorder, "alpha": 0.8}
    for i in range(grid.n_cols + 1):
        x = grid.x0 + i * grid.cell_m
        ax.plot([x, x], [y1, grid.y0], **style)
    for j in range(grid.n_rows + 1):
        y = grid.y0 - j * grid.cell_m
        ax.plot([grid.x0, x1], [y, y], **style)
    # clip_on: labels outside the current view (zoom panels) must not
    # be drawn NOR counted by bbox_inches="tight" — an off-view text
    # would otherwise blow the saved figure up to the full grid extent.
    for i in range(grid.n_cols):
        ax.text(
            grid.x0 + (i + 0.5) * grid.cell_m, grid.y0 + 0.1 * grid.cell_m,
            _col_label(i), ha="center", va="bottom",
            fontsize=7, color="0.25", zorder=zorder, clip_on=True,
        )
    for j in range(grid.n_rows):
        ax.text(
            grid.x0 - 0.1 * grid.cell_m, grid.y0 - (j + 0.5) * grid.cell_m,
            str(j + 1), ha="right", va="center",
            fontsize=7, color="0.25", zorder=zorder, clip_on=True,
        )
    for name, (axmin, aymin, axmax, aymax) in grid.aliases.items():
        ax.plot(
            [axmin, axmax, axmax, axmin, axmin],
            [aymin, aymin, aymax, aymax, aymin],
            color="tab:purple", lw=1.0, zorder=zorder,
        )
        ax.text(
            0.5 * (axmin + axmax), aymax, name, ha="center", va="bottom",
            fontsize=7, color="tab:purple", zorder=zorder, clip_on=True,
        )


def _add_coast(ax, coast, crs: str, coast_config=None, *, zorder: float = 0.5) -> bool:
    """Render xcoast land behind the mesh. Returns True on success;
    warns and returns False when xcoast or its data are unavailable."""
    try:
        import xcoast
        from xcoast.render import add_land_to_plain_axes

        mask = xcoast.load(coast, config=coast_config)
        land = mask.land_gdf.to_crs(crs)
        add_land_to_plain_axes(
            ax, land, facecolor="0.88", edgecolor="0.55",
            linewidth=0.3, zorder=zorder,
        )
        return True
    except Exception as e:  # pragma: no cover - environment dependent
        warnings.warn(f"land rendering skipped ({e})", stacklevel=2)
        return False


def add_atlas_grid(ax, crs: str = "EPSG:4326", grid=None,
                   color: str = "crimson", labels: bool = True):
    """Overlay the atlas reference grid (user-AI location pointing:
    column letters west->east, row numbers north->south, quadrant
    a/b/c/d). Default grid = gridref.TOKYO_BAY_GRID. Draws only
    within the current axis limits; call AFTER setting xlim/ylim."""
    import numpy as np
    from pyproj import Transformer

    from fvcom_mesh_tools.gridref import TOKYO_BAY_GRID

    g = grid or TOKYO_BAY_GRID
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    if crs in ("EPSG:4326", "lonlat", "4326"):
        def tf(lon, lat):
            return lon, lat
        inv = tf
    else:
        _t = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
        _i = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)
        tf, inv = _t.transform, _i.transform
    lo0, la0 = inv(x0, y0)
    lo1, la1 = inv(x1, y1)
    for i in range(g.ncol + 1):
        gx = g.lon0 + i * g.dlon
        if lo0 - g.dlon < gx < lo1 + g.dlon:
            xs, ys = tf(np.full(50, gx),
                        np.linspace(la0, la1, 50))
            ax.plot(xs, ys, color=color, lw=0.5, alpha=0.5,
                    zorder=4)
    for j in range(g.nrow + 1):
        gy = g.lat1 - j * g.dlat
        if la0 - g.dlat < gy < la1 + g.dlat:
            xs, ys = tf(np.linspace(lo0, lo1, 50),
                        np.full(50, gy))
            ax.plot(xs, ys, color=color, lw=0.5, alpha=0.5,
                    zorder=4)
    if not labels:
        return
    for i in range(g.ncol):
        gx = g.lon0 + (i + 0.5) * g.dlon
        if lo0 < gx < lo1:
            for la in (la0 + 0.012 * (la1 - la0),
                       la1 - 0.02 * (la1 - la0)):
                px, py = tf(gx, la)
                ax.text(px, py, g.col_letter(i), color=color,
                        ha="center", fontsize=11,
                        fontweight="bold", zorder=5)
    for j in range(1, g.nrow + 1):
        gy = g.lat1 - (j - 0.5) * g.dlat
        if la0 < gy < la1:
            for lo in (lo0 + 0.008 * (lo1 - lo0),
                       lo1 - 0.02 * (lo1 - lo0)):
                px, py = tf(lo, gy)
                ax.text(px, py, str(j), color=color, va="center",
                        fontsize=11, fontweight="bold", zorder=5)


def plot_mesh_overview(
    mesh: Fort14Mesh,
    out_png: str | Path,
    *,
    crs: str = "EPSG:32654",
    cell_m: float | None = 5000.0,
    aliases: dict[str, tuple[float, float, float, float]] | None = None,
    coast: str | tuple[float, float, float, float] | None = None,
    coast_config=None,
    zoom: str | tuple[float, float, float, float] | None = None,
    zoom_margin_frac: float = 0.05,
    dpi: int = 200,
    title: str | None = None,
    atlas: bool = True,
) -> Path:
    """Whole-domain (or zoomed) mesh figure per kickoff §10.

    ``coast`` is an xcoast preset name or a lon/lat bbox (rendered
    behind the mesh, reprojected to ``crs``); ``None`` disables land.
    ``zoom`` is a grid/alias reference ("C4", "C4-D5") or a CRS bbox;
    ``cell_m=None`` disables the reference grid. Returns the PNG path.
    """
    import matplotlib.pyplot as plt

    out_png = Path(out_png).resolve()
    bbox = mesh.bbox
    grid = (
        ReferenceGrid.from_bbox(bbox, cell_m=cell_m, aliases=aliases)
        if cell_m else None
    )

    fig, ax = plt.subplots(figsize=(11, 10))
    if coast is not None:
        _add_coast(ax, coast, crs, coast_config)
    ax.triplot(
        mesh.nodes[:, 0], mesh.nodes[:, 1], mesh.elements,
        color="#33658a", lw=0.12, zorder=2,
    )
    for seg in mesh.open_boundaries:
        seg = np.asarray(seg, dtype=np.int64)
        ax.plot(
            mesh.nodes[seg, 0], mesh.nodes[seg, 1],
            color="red", lw=1.8, zorder=4, label="open boundary",
        )
    if grid is not None:
        _add_reference_grid(ax, grid)

    if zoom is not None:
        if isinstance(zoom, str):
            if grid is None:
                raise ValueError("zoom by reference requires cell_m (grid enabled)")
            zb = grid.parse_ref(zoom)
        else:
            zb = zoom
        mx = zoom_margin_frac * max(zb[2] - zb[0], zb[3] - zb[1])
        ax.set_xlim(zb[0] - mx, zb[2] + mx)
        ax.set_ylim(zb[1] - mx, zb[3] + mx)
    else:
        mx = 0.02 * max(bbox[2] - bbox[0], bbox[3] - bbox[1])
        ax.set_xlim(bbox[0] - mx, bbox[2] + mx)
        ax.set_ylim(bbox[1] - mx, bbox[3] + mx)

    ax.set_aspect("equal")
    ax.set_xlabel(f"x (m, {crs})")
    ax.set_ylabel(f"y (m, {crs})")
    if atlas:
        try:
            add_atlas_grid(ax, crs=crs)
        except Exception as exc:  # never block a figure on the grid
            import warnings

            warnings.warn(f"atlas grid overlay failed: {exc}")
    zoom_note = f"  zoom={zoom}" if isinstance(zoom, str) else ""
    ax.set_title(
        title
        or f"{mesh.title}  NP={mesh.n_nodes:,} NE={mesh.n_elements:,}{zoom_note}"
    )
    if mesh.open_boundaries:
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(handles[:1], labels[:1], loc="lower right", fontsize=8)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    # bbox_inches="tight" only for the whole-domain view: on zoom
    # panels it would expand the canvas to every off-view artist.
    fig.savefig(
        out_png, dpi=dpi,
        bbox_inches="tight" if zoom is None else None,
    )
    plt.close(fig)
    return out_png


__all__ = [
    "MESH_PNG_DPI",
    "ReferenceGrid",
    "plot_mesh_overview",
]
