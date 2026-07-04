# DATA_INVENTORY — input datasets for Tokyo Bay FVCOM meshing

Audited 2026-07-03 (header-level inspection on GENKAI: `ncdump -h`, `ogrinfo -so`,
`gdalinfo`, file listings; no full-data reads). Root: `$DATA_DIR`
(= group-shared `Data/` directory; never hardcode the absolute path).

Project fixed settings (see `fvcom_mesh_kickoff.md` §5): horizontal CRS **UTM 54N**
(all sources reprojected), vertical datum **T.P.**, minimum depth clip **2 m**.

Sign convention used below: "depth" is positive down; "elevation" positive up.

---

## Precedence (resolved)

| Rank | Bathymetry (within footprint) | Coastline |
|---|---|---|
| 1 | `Futtsu_JFA_2023` | **`OSM` (land polygons 2026-03)** |
| 2 | `M7001` / `tokyo_bay` (see overlap note) | `coastline/tokyo_bay` (MLIT C23) |
| 3 | `mesh500` | `GSHHS` (fallback only) |
| 4 | `GEBCO` (fallback only) | — |

Coastline precedence REVISED 2026-07-04 (user decision, superseding the
kickoff §4 order): OSM is the default — the MLIT C23 extract (2021) is
older than the OSM land polygons (2026-03), and OSM applies worldwide,
matching the general-purpose-toolkit goal. In practice the meshing
shoreline is the xcoast "true land" product (land polygons minus
inland water); the same polygons serve as the snap/QA conformity
reference so mesh-vs-land figures are exact conformity checks.

Overlap resolution within rank 2: `M7001` is the survey-authority source (JHA
soundings/contours) but has **no intertidal data** (marks start ≥ 1 m below chart
datum); the `tokyo_bay` 30 m grid (`depth_0030-*.nc`) fills tidal flats (e.g.
Sanbanze) inside its footprint. Rule: **M7001 wins where both have data; the
tokyo_bay 30 m grid wins in the intertidal/shallow gap (< ~1 m below CD) and
anywhere M7001 has no soundings**.

---

## Bathymetry — `$DATA_DIR/bathymetry/`

### 1. Futtsu_JFA_2023 — precedence #1 (15 GB, 43 tiles)

| Field | Value |
|---|---|
| CRS | **EPSG:6677** (JGD2011 / Japan Plane Rectangular CS IX), meters — `.prj` present; CSV headers confirm zone 9 |
| Vertical datum | **T.P.** (directory suffix `_TP`; elevations signed, negative below datum) |
| Resolution | **0.5 m grid** (800 × 600 points per 400 m × 300 m tile); `ORIGINAL_TP/` holds the denser raw point cloud (~14 GB) |
| Coverage | Futtsu area only (sample tile corner 35°19′45.75″N, 139°48′40.80″E). Full mosaic bbox: needs batch job (trivial 43-header scan) |
| Source | Japan Fisheries Agency survey, 測量年 2022 (files dated 2022-12); "JFA 2023" per directory name |
| License | **Undocumented on disk** — record provenance from external correspondence |
| Formats | `05MCSV_TP/` gridded CSV-like `.txt` (`id, Y_easting_m, X_northing_m, elev_m, flag`); `MESH_TP/` byte-identical copies; `05MDTM/` JPGIS `.lem` binary + Shift-JIS headers; `SHP/` 43 shapefiles |
| Known issues | Only source in a projected CRS — reproject EPSG:6677 → UTM 54N (not the geographic→UTM path used for the others) |

### 2. M7001 — precedence #2 (409 MB, JHA vector product)

| Field | Value |
|---|---|
| CRS | WGS84 geographic (datum flag `1 = WGS84` per bundled JHA format spec; shapefiles have **no `.prj`** — coordinates are lon/lat degrees) |
| Vertical datum | **Chart datum (基本水準面 ≈ LLW), NOT T.P.** — conversion **`depth_TP = depth_CD + 1.13 m`** (documented in the dataset README) is **mandatory**; the blanket "assume MSL = T.P." rule must NOT be applied here |
| Resolution | Vector soundings/contours (3.95 M points; in the Tokyo Bay domain 411,170 depth points, node-to-nearest-sounding p50 0.075 km / p95 0.51 km) |
| Coverage | 137.93–140.88 E, 33.65–35.88 N (coastline layer; contours reach 141.56 E) — southern Kanto |
| Source | JHA (日本水路協会) M7001 関東南部 Ver. 2.4 (2022-10), derived from JCG charts |
| License | Commercial product — **not redistributable**; keep outside the repo |
| Formats | J-BIRD fixed-column ASCII (324 MB) + 3 shapefiles: 海岸線 coastline (5,329 LineString), 低潮線 low-tide line (3,317), 等深線 depth contours (38,954) |
| Known issues | **No intertidal data** (depths start ≥ 1 m below CD) — tidal flats must come from `tokyo_bay` 30 m grid / Futtsu_JFA within their footprints |

### 3. tokyo_bay — precedence #2, intertidal gap-filler (131 MB, regridded NetCDF)

| Field | Value |
|---|---|
| CRS | **EPSG:4326** (WKT embedded in every file) |
| Vertical datum | Not stated in headers — **assume MSL = T.P. 0 per project rule** (supported by the M7001 README noting these grids are T.P.-referenced) |
| Resolution / coverage | `depth_0030-11+12+13+14+15.nc` ~30 m, inner bay (139.565–140.172 E, 35.10–35.86 N); `depth_0090-05+06+07.nc` ~90 m (139.05–140.51 E); `depth_0270-03.nc` ~270 m (138.84–141.08 E); `depth_0810-01.nc` ~810 m (129.92–143.27 E, 29.46–38.75 N) |
| Source | Regridded set from a prior study (GDAL-translated NetCDF, Sep 2021; duplicated under `coastline/tokyo_bay/Wang/`) |
| License | Derived research product; undocumented — treat as internal |
| Known issues | 16 `dem_*.nc` files are **scenario DEMs** (reclamation/SLR variants), not survey data — exclude from the depth-merge stack |

### 4. mesh500 — precedence #3 (12 MB, ASCII grid points)

| Field | Value |
|---|---|
| CRS | Geographic lat/lon (datum unstated; assume WGS84) |
| Vertical datum | Unstated — **assume MSL = T.P. 0 per project rule** |
| Resolution | ~500 m grid (consistent with JODC J-EGG500; source inferred, not documented on disk) |
| Coverage | 33–37 N, 138–142 E (14 tiles, 443,670 points; Kanto coastal seas incl. Tokyo, Sagami, Kashima) |
| Format | 4-column ASCII: `flag lat lon depth_m` |
| License / issues | No README/license on disk — provenance undocumented; role = Japanese coast **outside** Tokyo Bay |

### 5. GEBCO — precedence #4, fallback only (4.0 GB)

| Field | Value |
|---|---|
| Status | **BLOCKED: `GEBCO_2024.nc` is actually a ZIP archive** (contains the real 7.47 GB NetCDF + grid documentation + terms-of-use PDFs). Unusable until extracted — **needs a batch job** (disk: ~7.5 GB) |
| Nominal spec | GEBCO 2024, global 15 arc-sec, EPSG:4326, elevation rel. MSL — confirm from headers after extraction |
| Known issues | Low accuracy near Japan — fallback only (outer domain far-field) |

### Other subdirs (not in the precedence stack)

- `SRTM15plus/` — `SRTM15+.nc` 6.1 GB global 15″ topo-bathy; alternative far-field
  fallback (readable now, unlike GEBCO). The MATLAB reference scripts used an
  SRTM15+ slice for the outer region.
- `CUDEM/` — 3 small test-fixture NetCDFs (not Tokyo Bay).
- `osaka_bay/` — SRTM15+ GeoTIFF crop (used by the Osaka Bay generality PoC).

---

## Coastline — `$DATA_DIR/coastline/`

### 1. tokyo_bay (MLIT C23) — precedence #1 (68 MB)

| Field | Value |
|---|---|
| CRS | WGS84 geographic (`_INNER`/`_OUTER` have `.prj`; the merged `C23-06_TOKYOBAY.shp` **lacks `.prj`** — coordinates are degrees, treat as WGS84) |
| Content | `MLIT_C23/`: `C23-06_TOKYOBAY.shp` (1,325 LineString), `_INNER` (162), `_OUTER` (1,094), + DXF. Extent 139.142–140.394 E, 34.898–35.860 N |
| Source | MLIT National Land Numerical Information C23 coastline (bay-merged custom extract; Sep 2021 originals preserved in `Wang/`) |
| License | MLIT NLNI free-use terms (no license file on disk — flag) |
| Datum note | Planimetric; C23 shoreline is nominally spring-high-water line on topo maps — acceptable under the no-wetting/drying assumption |
| Known issues | **`Futtsu/` polygon shapefiles (`Futtsu_coastline.shp`, `coastline_2.shp`) have continental-scale extents (e.g. 72.9–169.3 E) — suspect stray vertices; verify geometry before any use** (these are the shapes the MATLAB reference scripts consume). `SMS/` holds legacy manual SMS projects — historical reference only (SMS phase-out policy) |

### 2. OSM — precedence #2 (~6 GB, shapefiles only)

| Field | Value |
|---|---|
| CRS | EPSG:4326 |
| Content | `land-polygons-split-4326/` — osmdata.openstreetmap.de land polygons (862,652 features, global, data date **2026-03-21**): the canonical OSM coastline-as-polygons product. `geofabrik_kanto/` — 18-layer Geofabrik Kanto extract (data date 2026-03-19): `water_a`, `waterways`, `natural` are the layers relevant to river-bank/water geometry |
| License | ODbL 1.0 (READMEs on disk) |
| Known issues | Shapefile only (no PBF/GeoJSON); river banks come from `water_a`/`waterways`, not the land polygons |

### 3. GSHHS — precedence #3, fallback only (311 MB)

| Field | Value |
|---|---|
| CRS | WGS84 (`.prj` present) |
| Content | GSHHG v2.3.7 (2017-06-15), 5 resolutions (c/l/i/h/f) × levels L1–L6 |
| License | LGPL (per README) |
| Known issues | 2017 vintage, generalized — misses recent Tokyo Bay reclamation; last-resort fallback (outer-domain far field only) |

---

## Reference meshes — `$DATA_DIR/mesh/`

- `fvcom/tokyo_bay/` — empty (destination for produced FVCOM inputs, presumably).
- `reference/tokyo_bay/` — legacy OceanMesh2D-era Futtsu meshes (`tb_futtsu.14`
  2021-12, `tb_futtsu20220311.14`, `tokyobay_futtsu5.0.14` + `.mat` +
  `genMesh_futtsu*.m`). Useful as comparison baselines; unrelated to
  `bathymetry/mesh500` despite the name.

---

## Action items surfaced by this audit

1. **GEBCO**: extract the zip via batch job (or decide SRTM15+ suffices as far-field
   fallback and drop GEBCO from the stack).
2. **Futtsu coastline shapefiles**: geometry check (stray-vertex removal) before use.
3. **M7001 datum shift**: encode `+1.13 m` CD→T.P. in the depth-merge recipe as a
   per-source `datum_shift_m` field — never a global assumption.
4. **Missing `.prj`s** (`C23-06_TOKYOBAY.shp`, M7001 shapefiles): assign WGS84
   explicitly at load time in the pipeline; do not rely on GDAL defaults.
5. **Provenance/licenses** for Futtsu_JFA_2023, mesh500, MLIT extracts: record from
   external sources; none are documented on disk.
