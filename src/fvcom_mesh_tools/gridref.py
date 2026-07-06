"""Atlas-style grid references for pointing at mesh locations.

Users reference locations as ``"F12"`` (column letter west->east,
row number north->south) with an optional quadrant suffix
(``a``=NW, ``b``=NE, ``c``=SW, ``d``=SE): ``"F12c"``. The same
string resolves to a lon/lat polygon, so recipe directives can say
``cell: F12`` instead of coordinates.
"""

from __future__ import annotations

import string

__all__ = ["GridRef", "TOKYO_BAY_GRID"]


class GridRef:
    def __init__(self, lon0, lat0, lon1, lat1, dlon=0.05, dlat=0.05):
        self.lon0, self.lat0 = float(lon0), float(lat0)
        self.lon1, self.lat1 = float(lon1), float(lat1)
        self.dlon, self.dlat = float(dlon), float(dlat)
        self.ncol = int(round((self.lon1 - self.lon0) / self.dlon))
        self.nrow = int(round((self.lat1 - self.lat0) / self.dlat))
        if self.ncol > 26:
            raise ValueError("more than 26 columns; enlarge dlon")

    def col_letter(self, i):
        return string.ascii_uppercase[i]

    def cell_bounds(self, ref):
        """'F12' or 'F12c' -> (lon_min, lat_min, lon_max, lat_max)."""
        ref = ref.strip().upper()
        quad = None
        if ref and ref[-1] in "ABCD" and ref[:-1] and ref[-2].isdigit():
            quad = ref[-1].lower()
            ref = ref[:-1]
        col = string.ascii_uppercase.index(ref[0])
        row = int(ref[1:])
        if not (0 <= col < self.ncol and 1 <= row <= self.nrow):
            raise ValueError(f"cell {ref} outside grid "
                             f"(A1..{self.col_letter(self.ncol-1)}"
                             f"{self.nrow})")
        x0 = self.lon0 + col * self.dlon
        # rows count from the NORTH edge (map-atlas convention)
        y1 = self.lat1 - (row - 1) * self.dlat
        y0 = y1 - self.dlat
        x1 = x0 + self.dlon
        if quad:
            xm, ym = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
            x0, x1 = (x0, xm) if quad in "ac" else (xm, x1)
            y0, y1 = (ym, y1) if quad in "ab" else (y0, ym)
        return (x0, y0, x1, y1)

    def cell_polygon(self, ref):
        x0, y0, x1, y1 = self.cell_bounds(ref)
        return [[x0, y0], [x1, y0], [x1, y1], [x0, y1]]

    def point_to_cell(self, lon, lat, quadrant=False):
        col = int((lon - self.lon0) / self.dlon)
        row = int((self.lat1 - lat) / self.dlat) + 1
        if not (0 <= col < self.ncol and 1 <= row <= self.nrow):
            raise ValueError("point outside grid")
        ref = f"{self.col_letter(col)}{row}"
        if quadrant:
            x0, y0, x1, y1 = self.cell_bounds(ref)
            qx = "a" if lon < 0.5 * (x0 + x1) else "b"
            ref += ("a" if qx == "a" else "b") if lat >= 0.5 * (
                y0 + y1) else ("c" if qx == "a" else "d")
        return ref


# canonical Tokyo Bay grid = the v6 prep bbox at 0.05 deg
TOKYO_BAY_GRID = GridRef(139.40, 34.90, 140.30, 35.90)
