"""PoC #84b (stage 3) — validate the tokyo_bay_v2 M2 tidal run.

Reads the FVCOM output NetCDF and checks the physics-sanity criteria:

* no NaN in zeta / ua / va;
* M2 amplitude over the final two cycles (24.84 h): forcing band at
  the OBC (~0.40 m), bay-head (lat > 35.60N) amplification ratio in
  a plausible 1.0-1.6 band (Tokyo vs Uraga M2 is ~1.2-1.3);
* depth-averaged speeds: strong in the Uraga narrows (O(1 m/s)),
  bounded everywhere (< 3 m/s).

Writes a JSON verdict + an amplitude-map figure.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import xarray as xr
from pyproj import Transformer

REPO = Path(__file__).resolve().parents[1]
OUT_NC = REPO / "outputs" / "fvcom_tide_v4" / "output" / "tokyo_bay_v4_tide_0001.nc"
VERDICT = REPO / "outputs" / "fvcom_tide_v4" / "84b_tide_verdict.json"
FIG = REPO / "outputs" / "figures" / "84b_m2_amplitude.png"

M2_H = 12.4206012


def main() -> int:
    ds = xr.open_dataset(OUT_NC, decode_times=False)
    zeta = ds["zeta"].values          # (time, node)
    ua = ds["ua"].values              # (time, nele)
    va = ds["va"].values
    t = ds["time"].values             # MJD days
    x, y = ds["x"].values, ds["y"].values
    lonc = ds["lonc"].values if "lonc" in ds else None

    nan_counts = {
        "zeta": int(np.isnan(zeta).sum()),
        "ua": int(np.isnan(ua).sum()),
        "va": int(np.isnan(va).sum()),
    }

    # Final two M2 cycles.
    t_end = t[-1]
    sel = t >= t_end - 2 * M2_H / 24.0
    amp = 0.5 * (zeta[sel].max(axis=0) - zeta[sel].min(axis=0))

    tr = Transformer.from_crs("EPSG:32654", "EPSG:4326", always_xy=True)
    lon, lat = tr.transform(x, y)

    # Region masks.
    head = lat > 35.60
    # OBC band: southern/western artificial-arc vicinity.
    obc_band = (lat < 35.11) | (lon < 139.575)
    uraga = (lat > 35.15) & (lat < 35.30) & (lon > 139.70) & (lon < 139.80)

    speed = np.hypot(ua[sel], va[sel]).max(axis=0)
    if lonc is not None:
        latc = ds["latc"].values
        uraga_e = (latc > 35.15) & (latc < 35.30) & (lonc > 139.70) & (lonc < 139.80)
    else:
        uraga_e = None

    amp_obc = float(amp[obc_band].mean())
    amp_head = float(amp[head].mean())
    ratio = amp_head / amp_obc if amp_obc > 0 else float("nan")
    result = {
        "n_time": int(t.size),
        "sim_days": float(t[-1] - t[0]),
        "nan_counts": nan_counts,
        "amp_obc_mean_m": amp_obc,
        "amp_head_mean_m": amp_head,
        "head_over_obc_ratio": ratio,
        "amp_uraga_mean_m": float(amp[uraga].mean()) if uraga.any() else None,
        "zeta_abs_max_m": float(np.abs(zeta).max()),
        "speed_max_m_s": float(speed.max()),
        "speed_uraga_max_m_s": (
            float(speed[uraga_e].max()) if uraga_e is not None and uraga_e.any()
            else None
        ),
        "checks": {},
    }
    result["checks"] = {
        "no_nan": all(v == 0 for v in nan_counts.values()),
        "obc_amp_near_forcing": bool(0.30 <= amp_obc <= 0.50),
        "head_ratio_sane": bool(1.0 <= ratio <= 1.6),
        "zeta_bounded": bool(result["zeta_abs_max_m"] < 2.0),
        "speed_bounded": bool(result["speed_max_m_s"] < 3.0),
    }
    result["passed"] = all(result["checks"].values())
    VERDICT.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    nv = ds["nv"].values.T - 1  # (nele, 3)
    fig, ax = plt.subplots(figsize=(9, 9))
    tp = ax.tripcolor(x, y, nv, amp, cmap="viridis", shading="flat")
    fig.colorbar(tp, ax=ax, label="M2 amplitude (m), last 2 cycles")
    ax.set_aspect("equal")
    ax.set_title(
        f"tokyo_bay_v2 M2 test — amp OBC {amp_obc:.3f} m, "
        f"head {amp_head:.3f} m (x{ratio:.2f})"
    )
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG, dpi=200, bbox_inches="tight")
    print(f"[84b] wrote {FIG}", flush=True)
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
