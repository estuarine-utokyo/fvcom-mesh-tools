"""Correlate final-QA violations with forced two-row ladder bands.

Reads outputs/sample_repro/sample_repro_final_qa.json (violation
coordinates in UTM 54N) and outputs/sample_repro/ladder_bands.json
(band pfix nodes in lon/lat), prints per-violation distance to the
nearest ladder node and a per-band violation tally.  Diagnostic
only -- the accept/reject decision stays with the QA gates.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from pyproj import Transformer

OUT = Path(__file__).resolve().parents[1] / "outputs/sample_repro"

qa = json.loads((OUT / "sample_repro_final_qa.json").read_text())
bands_ll = json.loads((OUT / "ladder_bands.json").read_text())
tr = Transformer.from_crs(4326, 32654, always_xy=True)

band_xy = []
for b in bands_ll:
    a = np.asarray(b, float)
    x, y = tr.transform(a[:, 0], a[:, 1])
    band_xy.append(np.column_stack([x, y]))
allb = np.vstack(band_xy) if band_xy else np.zeros((0, 2))
owner = np.concatenate([np.full(len(b), i)
                        for i, b in enumerate(band_xy)]) \
    if band_xy else np.zeros(0, int)
print(f"[361] bands: {len(band_xy)}, ladder nodes: {len(allb)}")

WATCH = {"c1_min_angle", "c2_max_angle", "c4_area_change",
         "c5_valence", "delaunay_local", "channel_wh"}
tally: dict[int, int] = {}
for ch in qa["checks"]:
    if ch["check_id"] not in WATCH or not ch["offenders"]:
        continue
    print(f"[361] {ch['check_id']}: {ch['n_violations']} "
          f"({ch['status']})")
    for off in ch["offenders"]:
        if "x" not in off or len(allb) == 0:
            continue
        d = np.hypot(allb[:, 0] - off["x"], allb[:, 1] - off["y"])
        j = int(d.argmin())
        near = d[j] < 800.0
        if near:
            tally[int(owner[j])] = tally.get(int(owner[j]), 0) + 1
        print(f"[361]   {off.get('kind', '?')} {off.get('id', '?')}"
              f" at ({off['x']:.0f},{off['y']:.0f})"
              f" -> band {owner[j]} dist {d[j]:.0f} m"
              f" {'NEAR-LADDER' if near else 'elsewhere'}")
print(f"[361] per-band violation tally (<800 m): "
      f"{dict(sorted(tally.items())) or 'none near any band'}")
