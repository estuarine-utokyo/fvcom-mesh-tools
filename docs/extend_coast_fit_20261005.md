# Coastline fit of the extended Tokyo Bay mesh (2026-10-05)

Question: does the extended mesh (`tokyo_bay_enshu`, 14,673 nodes, 27,011 elements, GENKAI
build, grd sha256 `966f8d05...`) follow the coastline well enough? The build's own report
said 285 boundary nodes on land and a maximum offset of 709 m. Script
`notebooks/455_coast_fit_analysis.py`, job `jobs/genkai/455_coast_analysis.sh`
(job 7003166; figures were written beside the numbers, in `scratch/coastfit2/`).
All distances are against the land polygon the build itself used (`land_with_base.shp`).

## Node view (the 959 new land-boundary nodes)

| | value |
|---|---|
| median \|offset\| | 0.003 m |
| p90 / p99 / max \|offset\| | 108 / 422 / 709 m |
| p90 / max of \|offset\| / local edge length h | 0.19 / 0.73 |
| in water / on land | 676 / 283 |
| on land by more than 1 m / 10 m / 0.1 h | 86 / 69 / 50 |
| in water by more than 1 m / 0.1 h | 181 / 82 |
| median h | 647 m |

The "285 nodes on land" of the build report is mostly sign noise: the fit puts nodes ON the
coastline (median 0.003 m), and half of them land on the land side of it by millimetres.
The nodes that are really on land are 69 (7 %) by more than 10 m, and none
by more than 0.73 h.

## Coast view (coast within 2 km of a new node, not the base's own)

14118 samples at 100 m (1412 km): distance to the nearest boundary edge of the new
mesh, median 91 m, p90 481 m (coast_h is 600 m), p99 1827 m;
64 % within 150 m; 7.0 % (about 100 km) farther than 600 m.

## What the departures are

* The largest omitted stretch, about 17 km at 138.50E 35.00N, is the **Shimizu port basin**
  (Suruga Bay): a harbour with a narrow entrance, left as land by the resolution policy
  (water below the element size is land). Not a defect.
* The other omitted samples are small coves, river mouths and harbours of the Izu peninsula,
  the Boso coast and the Izu islands, by the same policy.
* The nodes on land by more than 10 m sit on the islands (Oshima, Niijima, Miyakejima, ...)
  and capes, where the island is only a few elements across: a boundary node whose move
  onto the coast would take an element under the 30 degree minimum angle stays put, and
  the chord between two such nodes cuts the corner by up to 0.7 h.

Verdict: no systematic fault. The mesh follows the coastline to 0.2 h (p90); the exceptions
are islands and capes at the resolution limit, and harbours treated as land. A tighter fit
would cost minimum angle on the islands; it is not worth it unless an island's own
dynamics matters (then give that island its own coast_h in the recipe).
