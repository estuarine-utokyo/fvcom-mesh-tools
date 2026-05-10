"""PoC #38: Phase E Stage 2 (medial-axis CDT) end-to-end validation.

PoC #37 closed the Go/No-Go on Stage 2 and identified
``--under-resolved-min-channel-elements 10`` as the production sweet
spot (51 components / 1,070 flagged elements / mean ``long_axis_m /
h_local_median ≈ 6.5`` on the cleaned PoC #19 Tokyo Bay mesh, with
medial-axis estimate 0.66× the centroid-widen cost).

This PoC validates the implementation end-to-end: applies the new
``mode='medial'`` to the same cleaned PoC #19 input, compares
quality / topology / detector-6 residuals against the existing
centroid-widen Phase E (PoC #29), and reports skip-reason histogram
for any components rejected during retriangulation.

The medial mode is only meaningful for components above the filter
threshold; the comparison is therefore "medial at min_n=10" vs
"widen at min_n=10" — both run on the same input. Since the existing
PoC #29 widen ran with the unfiltered detector, we additionally run
a fresh widen at min_n=10 here so the two are directly comparable.

Outputs:
    outputs/38_phase_e_medial.14
    outputs/38_phase_e_widen_at_min_n_10.14
    outputs/38_phase_e_medial_summary.txt
    outputs/38_phase_e_medial_summary.json
"""
from __future__ import annotations

import json
import time
from pathlib import Path

from fvcom_mesh_tools.diagnostics import under_resolved_channels_flag
from fvcom_mesh_tools.io import read_fort14, write_fort14
from fvcom_mesh_tools.mesh_clean import repair_under_resolved_channels
from fvcom_mesh_tools.quality import compute_metrics

REPO = Path(__file__).resolve().parent.parent
OUT_DIR = REPO / "outputs"
INPUT = OUT_DIR / "19_tokyo_bay_oceanmesh_cleaned.14"
SUMMARY_TXT = OUT_DIR / "38_phase_e_medial_summary.txt"
SUMMARY_JSON = OUT_DIR / "38_phase_e_medial_summary.json"
MEDIAL_OUT = OUT_DIR / "38_phase_e_medial.14"
WIDEN_OUT = OUT_DIR / "38_phase_e_widen_at_min_n_10.14"

MIN_W_H = 3.0
MIN_CHANNEL_ELEMENTS = 10


def _detector6_residual(mesh) -> dict:
    flag, _info = under_resolved_channels_flag(
        mesh, min_w_h=MIN_W_H, min_channel_elements=1,
    )
    flag_filt, _ = under_resolved_channels_flag(
        mesh, min_w_h=MIN_W_H, min_channel_elements=MIN_CHANNEL_ELEMENTS,
    )
    return {
        "n_flagged_unfiltered": int(flag.sum()),
        "n_flagged_filtered_min_n_10": int(flag_filt.sum()),
    }


def _summary_block(label: str, mesh, info: dict, wall_s: float) -> dict:
    metrics = compute_metrics(mesh)
    detector = _detector6_residual(mesh)
    return {
        "label": label,
        "wall_s": wall_s,
        "info": info,
        "n_nodes": int(mesh.n_nodes),
        "n_elements": int(mesh.n_elements),
        "metrics": metrics,
        "detector6": detector,
    }


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if not INPUT.exists():
        raise SystemExit(f"input missing: {INPUT}")

    mesh_in = read_fort14(INPUT)
    print(
        f"[38] input: NP={mesh_in.n_nodes:,}  NE={mesh_in.n_elements:,}"
    )
    in_block = _summary_block("input", mesh_in, {}, 0.0)

    print("[38] running mode=widen min_channel_elements=10 (baseline)")
    t0 = time.perf_counter()
    widen_mesh, widen_info = repair_under_resolved_channels(
        mesh_in, mode="widen", min_w_h=MIN_W_H,
        min_channel_elements=MIN_CHANNEL_ELEMENTS,
    )
    widen_block = _summary_block(
        "widen_min_n_10", widen_mesh, widen_info, time.perf_counter() - t0,
    )
    write_fort14(widen_mesh, WIDEN_OUT)
    print(
        f"[38]   NP={widen_mesh.n_nodes:,}  NE={widen_mesh.n_elements:,}  "
        f"alpha={widen_block['metrics']['alpha_mean']:.4f}"
    )

    print("[38] running mode=medial min_channel_elements=10 (Stage 2)")
    t0 = time.perf_counter()
    medial_mesh, medial_info = repair_under_resolved_channels(
        mesh_in, mode="medial", min_w_h=MIN_W_H,
        min_channel_elements=MIN_CHANNEL_ELEMENTS,
    )
    medial_block = _summary_block(
        "medial_min_n_10", medial_mesh, medial_info,
        time.perf_counter() - t0,
    )
    write_fort14(medial_mesh, MEDIAL_OUT)
    print(
        f"[38]   NP={medial_mesh.n_nodes:,}  NE={medial_mesh.n_elements:,}  "
        f"alpha={medial_block['metrics']['alpha_mean']:.4f}"
    )

    payload = {
        "min_w_h": MIN_W_H,
        "min_channel_elements": MIN_CHANNEL_ELEMENTS,
        "input": in_block,
        "widen": widen_block,
        "medial": medial_block,
    }
    SUMMARY_JSON.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    out_lines = [
        "PoC #38 Phase E Stage 2 (medial-axis CDT) validation",
        f"input: {INPUT.name}",
        f"min_w_h={MIN_W_H}  min_channel_elements={MIN_CHANNEL_ELEMENTS}",
        "",
        "                                 NP            NE      alpha     "
        "frac<20°    n_flag(filt)   wall_s",
        "  -----------------------  ----------    ----------  ---------  "
        "----------  --------------  -------",
    ]
    for blk in (in_block, widen_block, medial_block):
        m = blk["metrics"]
        det = blk["detector6"]
        out_lines.append(
            f"  {blk['label']:<23}  {blk['n_nodes']:>10,}    "
            f"{blk['n_elements']:>10,}  {m['alpha_mean']:>9.4f}  "
            f"{m['frac_lt_20deg'] * 100:>9.4f}%  "
            f"{det['n_flagged_filtered_min_n_10']:>14,}  "
            f"{blk['wall_s']:>7.1f}"
        )
    out_lines.append("")
    info = medial_block["info"]
    out_lines.append("medial mode info:")
    for k in (
        "n_components", "n_components_replaced", "n_components_skipped",
        "n_nodes_inserted", "n_elements_removed", "n_elements_inserted",
    ):
        out_lines.append(f"  {k}: {info.get(k, 0):,}")
    sr = info.get("skip_reasons", {})
    if sr:
        out_lines.append("  skip_reasons:")
        for k, v in sorted(sr.items(), key=lambda kv: -kv[1]):
            out_lines.append(f"    {k}: {v}")

    SUMMARY_TXT.write_text("\n".join(out_lines), encoding="utf-8")
    print()
    print("\n".join(out_lines))
    print(f"\nwrote {MEDIAL_OUT}")
    print(f"wrote {WIDEN_OUT}")
    print(f"wrote {SUMMARY_TXT}")
    print(f"wrote {SUMMARY_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
