"""``fmesh-mesh-pipeline`` CLI: progressive `clean → quality → repeat` loop.

Three cumulative *rungs* of `fmesh-mesh-clean` phases are applied:

    * **rung 0** — A + B + C: drop disjoint components, trim dead-ends,
      widen 1-cell channels. The conservative default of
      ``fmesh-mesh-clean``.
    * **rung 1** — rung 0 + D + F + G: balance over-connected nodes
      (Lawson edge swap), delete extreme-skew triangles, run
      Laplacian smoothing (with the flipped-triangle safety net).
    * **rung 2** — rung 1 + E: widen elements flagged as
      under-resolved by the medial-axis channel-width detector.
      Adds ~3× new elements per flagged element so this is the most
      destructive rung; reach for it only when the threshold gate
      still fails after rung 1.

After each rung the toolkit's unified quality metrics
(:func:`fvcom_mesh_tools.quality.compute_metrics`) are computed and
optionally checked against user-supplied thresholds:
``--min-alpha``, ``--max-frac-lt-20deg``, ``--max-valence``,
``--max-overconnected``, ``--max-flipped``,
``--max-disjoint-elems``. The loop stops at the first rung that
passes all thresholds; the corresponding mesh is written to
``OUTPUT``. If no rung passes (or thresholds were not supplied),
the *last* rung's mesh is written.

Exit code: 0 if every supplied threshold passes (or none were
supplied); 1 if any threshold fails after every attempted rung.

JSON output (``--summary PATH``, default
``<output_stem>_pipeline_summary.json``) records per-rung metrics,
the diff vs the input, and the threshold-check results so the
caller can audit which rung achieved which numbers.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from fvcom_mesh_tools.cli.meshclean import _infer_bbox
from fvcom_mesh_tools.diagnostics import (
    DEFAULT_MAX_NBR_ELEM,
    DEFAULT_MIN_W_H,
)
from fvcom_mesh_tools.io import read_fort14, write_fort14
from fvcom_mesh_tools.mesh_clean import (
    DEFAULT_BBOX_TOL_M,
    DEFAULT_SKEWED_MAX_ANGLE_DEG,
    DEFAULT_SKEWED_MIN_ANGLE_DEG,
    DEFAULT_SMOOTH_LAPLACIAN_ITERS,
    DEFAULT_SMOOTH_LAPLACIAN_TOL,
    clean_mesh,
)
from fvcom_mesh_tools.quality import (
    check_thresholds,
    compute_metrics,
    format_comparison_table,
    format_threshold_table,
)


def _build_rung_overlays(args: argparse.Namespace) -> list[tuple[str, dict]]:
    """Return the cumulative kwargs added at each rung, in order.

    Rung 0 contributes nothing beyond ``base_kwargs``; rung N's
    overlay is merged on top of rung N-1's accumulated dict.
    """
    return [
        ("rung0:A+B+C", {}),
        (
            "rung1:+D+F+G",
            {
                # Phase D
                "repair_overconnected_iters": int(args.overconnected_iters),
                "max_nbr_elem": int(args.max_nbr_elem),
                "overconn_min_angle_floor_deg": float(args.overconn_min_angle_floor),
                # Phase F
                "repair_skewed": True,
                "repair_skewed_min_angle_deg":
                    float(args.repair_skewed_min_angle_deg),
                "repair_skewed_max_angle_deg":
                    float(args.repair_skewed_max_angle_deg),
                # Phase G
                "smooth_laplacian": True,
                "smooth_laplacian_iters": int(args.smooth_laplacian_iters),
                "smooth_laplacian_tol": float(args.smooth_laplacian_tol),
                "smooth_repair_flipped": True,
            },
        ),
        (
            "rung2:+E",
            {
                "under_resolved_mode": args.under_resolved_mode,
                "under_resolved_min_w_h": float(args.under_resolved_min_w_h),
                "under_resolved_min_channel_elements":
                    int(args.under_resolved_min_channel_elements),
            },
        ),
    ]


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="fmesh-mesh-pipeline",
        description=(
            "Progressively clean a fort.14 through three cumulative "
            "rungs of fmesh-mesh-clean phases (A+B+C → +D+F+G → +E), "
            "evaluating fmesh-mesh-quality thresholds after each. "
            "Stops at the first passing rung; exits 1 on threshold "
            "failure when thresholds were supplied. Boundaries are "
            "re-derived via DEM-bbox proximity, matching "
            "fmesh-buildmesh."
        ),
    )
    p.add_argument("input", type=Path, help="Input fort.14.")
    p.add_argument("output", type=Path, help="Output fort.14.")

    p.add_argument(
        "--max-iters", type=int, default=3,
        help=(
            "Maximum number of rungs to attempt (1, 2, or 3). "
            "Rung indices: 0=A+B+C, 1=+D+F+G, 2=+E. Default 3 "
            "(attempt all)."
        ),
    )
    p.add_argument(
        "--best-rung", action="store_true",
        help=(
            "Disable the first-passing-rung early-stop. Attempt every "
            "rung up to --max-iters and pick the one that maximises "
            "alpha_mean among the gate-passing rungs (ties broken in "
            "favour of the lower rung index — the lighter repair). If "
            "no rung passes the gate (or no thresholds are supplied), "
            "the rung with the highest alpha_mean overall is chosen. "
            "More compute (every rung runs) but quality is monotonic in "
            "rung depth less often than one would hope, so this option "
            "is the way to find that out."
        ),
    )

    g_bbox = p.add_argument_group("boundary classification")
    g_bbox.add_argument(
        "--bbox", type=float, nargs=4, default=None,
        metavar=("XMIN", "YMIN", "XMAX", "YMAX"),
        help=(
            "DEM bbox driving open / land classification. Default: "
            "infer from input mesh's existing open-boundary nodes."
        ),
    )
    g_bbox.add_argument(
        "--bbox-tol-m", type=float, default=DEFAULT_BBOX_TOL_M,
        help=f"Tolerance in metres. Default {DEFAULT_BBOX_TOL_M:g}.",
    )
    g_bbox.add_argument(
        "--land-ibtype", type=int, default=20,
        help="ibtype for re-derived land segments. Default 20.",
    )
    g_bbox.add_argument(
        "--open-merge-coast-gap", type=int, default=0,
        help="Bridge a short land run between two open runs. Default 0.",
    )

    g_rung1 = p.add_argument_group("rung 1 tuning (Phase D / F / G)")
    g_rung1.add_argument(
        "--overconnected-iters", type=int, default=20,
        help="Phase D iter cap when rung 1+ is reached. Default 20.",
    )
    g_rung1.add_argument(
        "--max-nbr-elem", type=int, default=DEFAULT_MAX_NBR_ELEM,
        help=f"Phase D MAX_NBR_ELEM cap. Default {DEFAULT_MAX_NBR_ELEM}.",
    )
    g_rung1.add_argument(
        "--overconn-min-angle-floor", type=float, default=0.0,
        help="Phase D min-angle floor (deg). Default 0.",
    )
    g_rung1.add_argument(
        "--repair-skewed-min-angle-deg", type=float,
        default=DEFAULT_SKEWED_MIN_ANGLE_DEG,
        help=f"Phase F lower angle bound. Default {DEFAULT_SKEWED_MIN_ANGLE_DEG:g}°.",
    )
    g_rung1.add_argument(
        "--repair-skewed-max-angle-deg", type=float,
        default=DEFAULT_SKEWED_MAX_ANGLE_DEG,
        help=f"Phase F upper angle bound. Default {DEFAULT_SKEWED_MAX_ANGLE_DEG:g}°.",
    )
    g_rung1.add_argument(
        "--smooth-laplacian-iters", type=int,
        default=DEFAULT_SMOOTH_LAPLACIAN_ITERS,
        help=f"Phase G iter cap. Default {DEFAULT_SMOOTH_LAPLACIAN_ITERS}.",
    )
    g_rung1.add_argument(
        "--smooth-laplacian-tol", type=float,
        default=DEFAULT_SMOOTH_LAPLACIAN_TOL,
        help=f"Phase G tolerance. Default {DEFAULT_SMOOTH_LAPLACIAN_TOL:g}.",
    )

    g_rung2 = p.add_argument_group("rung 2 tuning (Phase E)")
    g_rung2.add_argument(
        "--under-resolved-mode",
        choices=["widen", "medial"], default="widen",
        help=(
            "Phase E policy in rung 2. 'widen' (default) inserts a "
            "centroid in every flagged element (cheap, lifts w/h to "
            "~1.73× the original). 'medial' replaces each face-face "
            "channel of >= --under-resolved-min-channel-elements "
            "members with the Stage 2 medial-axis CDT path (PoC #37 "
            "production sweet spot is min-channel-elements 10)."
        ),
    )
    g_rung2.add_argument(
        "--under-resolved-min-w-h", type=float, default=DEFAULT_MIN_W_H,
        help=f"Phase E w/h threshold. Default {DEFAULT_MIN_W_H:g}.",
    )
    g_rung2.add_argument(
        "--under-resolved-min-channel-elements", type=int, default=1,
        help=(
            "Phase E: ignore detector-6 flags whose channel component "
            "has fewer than N flagged elements. Default 1 (no filter); "
            "raise to skip the small isolated clusters PoC #35 "
            "characterised. Particularly relevant for "
            "--under-resolved-mode medial — PoC #37 recommends 10."
        ),
    )

    g_thresh = p.add_argument_group("quality thresholds (any → CI gate)")
    g_thresh.add_argument("--min-alpha", type=float, default=None)
    g_thresh.add_argument("--max-frac-lt-20deg", type=float, default=None)
    g_thresh.add_argument("--max-valence", type=int, default=None)
    g_thresh.add_argument("--max-overconnected", type=int, default=None)
    g_thresh.add_argument("--max-flipped", type=int, default=None)
    g_thresh.add_argument("--max-disjoint-elems", type=int, default=None)

    p.add_argument(
        "--summary", type=Path, default=None,
        help=(
            "Path for the JSON summary. Default: "
            "<output_stem>_pipeline_summary.json next to OUTPUT."
        ),
    )
    p.add_argument(
        "--quiet", action="store_true",
        help="Suppress the per-rung table and final report on stdout.",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if not args.input.exists():
        print(f"input not found: {args.input}", file=sys.stderr)
        return 2
    if not 1 <= args.max_iters <= 3:
        print("--max-iters must be in [1, 3].", file=sys.stderr)
        return 2
    if args.bbox_tol_m <= 0:
        print("--bbox-tol-m must be > 0.", file=sys.stderr)
        return 2
    if args.max_nbr_elem < 3:
        print("--max-nbr-elem must be >= 3.", file=sys.stderr)
        return 2
    if args.under_resolved_min_channel_elements < 1:
        print("--under-resolved-min-channel-elements must be >= 1.",
              file=sys.stderr)
        return 2
    if args.under_resolved_min_w_h <= 0:
        print("--under-resolved-min-w-h must be > 0.", file=sys.stderr)
        return 2
    if args.smooth_laplacian_iters < 1:
        print("--smooth-laplacian-iters must be >= 1.", file=sys.stderr)
        return 2
    if args.smooth_laplacian_tol <= 0:
        print("--smooth-laplacian-tol must be > 0.", file=sys.stderr)
        return 2
    if args.repair_skewed_min_angle_deg >= args.repair_skewed_max_angle_deg:
        print("--repair-skewed-min-angle-deg must be < --repair-skewed-max-angle-deg.",
              file=sys.stderr)
        return 2
    args.output.parent.mkdir(parents=True, exist_ok=True)

    mesh_in = read_fort14(args.input)
    if args.bbox is None:
        bbox, bbox_source = _infer_bbox(mesh_in)
    else:
        bbox = tuple(args.bbox)
        bbox_source = "user-supplied"

    base_kwargs: dict = dict(
        bbox=bbox,
        bbox_tol_m=float(args.bbox_tol_m),
        land_ibtype=int(args.land_ibtype),
        open_merge_coast_gap=int(args.open_merge_coast_gap),
    )

    thresholds = {
        k: v for k, v in {
            "min_alpha_mean": args.min_alpha,
            "max_frac_lt_20deg": args.max_frac_lt_20deg,
            "max_valence": args.max_valence,
            "max_overconnected": args.max_overconnected,
            "max_flipped": args.max_flipped,
            "max_disjoint_elems": args.max_disjoint_elems,
        }.items() if v is not None
    }

    rung_overlays = _build_rung_overlays(args)[: args.max_iters]
    initial_metrics = compute_metrics(mesh_in, max_nbr_elem=args.max_nbr_elem)

    # Each entry: history-dict (JSON-friendly) + the cleaned mesh
    # itself, which we keep in memory so --best-rung can pick a
    # non-final rung as the output.
    rung_results: list[tuple[dict, object]] = []

    cumulative_kwargs = dict(base_kwargs)
    for rung_idx, (rung_label, overlay) in enumerate(rung_overlays):
        cumulative_kwargs.update(overlay)
        cleaned, clean_info = clean_mesh(mesh_in, **cumulative_kwargs)
        metrics = compute_metrics(cleaned, max_nbr_elem=args.max_nbr_elem)
        passed, checks = check_thresholds(metrics, **thresholds)
        rung_results.append((
            {
                "rung_index": rung_idx,
                "rung_label": rung_label,
                "kwargs_added": overlay,
                "phases_run": [p["name"] for p in clean_info.get("phases", [])],
                "metrics": metrics,
                "passed": bool(passed),
                "checks": [c.to_dict() for c in checks],
            },
            cleaned,
        ))
        # First-passing-rung early-stop, unless --best-rung asked us
        # to explore every rung.
        if not args.best_rung and thresholds and passed:
            break

    history = [entry for entry, _ in rung_results]

    chosen_idx, selection_reason = _select_rung(
        rung_results, thresholds=bool(thresholds), best_rung=args.best_rung,
    )
    final_entry, final_mesh = rung_results[chosen_idx]
    final_metrics = final_entry["metrics"]
    final_passed = final_entry["passed"]
    final_rung_label = final_entry["rung_label"]

    write_fort14(final_mesh, args.output)
    summary_path = args.summary or args.output.with_name(
        args.output.stem + "_pipeline_summary.json"
    )
    payload = {
        "input_path": str(args.input.resolve()),
        "output_path": str(args.output.resolve()),
        "bbox": list(bbox),
        "bbox_source": bbox_source,
        "thresholds": thresholds,
        "max_iters": args.max_iters,
        "best_rung_mode": bool(args.best_rung),
        "initial_metrics": initial_metrics,
        "history": history,
        "final": {
            "rung_label": final_rung_label,
            "rung_index": final_entry["rung_index"],
            "metrics": final_metrics,
            "passed": (bool(final_passed) if thresholds else None),
            "selection_reason": selection_reason,
        },
    }
    summary_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    if not args.quiet:
        print(f"input:  {args.input}")
        print(f"output: {args.output}")
        print(f"bbox:   {bbox}  ({bbox_source})")
        print()
        rows = [("input", initial_metrics)]
        for h in history:
            rows.append((h["rung_label"], h["metrics"]))
        print(format_comparison_table(rows))
        print()
        print(f"final rung: {final_rung_label}  "
              f"(phases: {final_entry['phases_run']})")
        print(f"selection: {selection_reason}")
        if thresholds:
            print("\nthresholds (against final):")
            print(format_threshold_table(
                [_dict_to_check(c) for c in final_entry["checks"]]
            ))
            print(f"\noverall: {'PASS' if final_passed else 'FAIL'}")
        else:
            print("(no thresholds supplied — gate not evaluated)")
        print(f"\nwrote {summary_path}")

    if thresholds:
        return 0 if final_passed else 1
    return 0


def _select_rung(
    rung_results: list,
    *,
    thresholds: bool,
    best_rung: bool,
) -> tuple[int, str]:
    """Pick which rung's output to write as the final fort.14.

    Returns ``(index, reason)``. Selection rules:

    * **default (``best_rung=False``)**: the last rung that ran. With
      thresholds set, the loop already early-stops at the first
      passing rung, so "last" is "first passing" or "last attempted"
      on failure. Without thresholds, "last" is the deepest rung.
    * **``best_rung=True``**: among the rungs that pass the gate
      (or all rungs if no thresholds), pick the one with the highest
      ``alpha_mean``. Ties broken in favour of the lower rung index
      (lighter repair). If no rung passes, fall back to the rung
      with the highest ``alpha_mean`` overall.

    A NaN ``alpha_mean`` (e.g. an empty mesh) is ranked below any
    finite alpha so it never wins the tie.
    """
    if not rung_results:
        raise ValueError("rung_results must not be empty")

    if not best_rung:
        return len(rung_results) - 1, "first-passing-rung-or-last"

    def _alpha(entry: dict) -> float:
        v = entry["metrics"].get("alpha_mean")
        if v is None:
            return float("-inf")
        if isinstance(v, float) and v != v:  # NaN
            return float("-inf")
        return float(v)

    # Stable order: lower rung_idx first → max() returns the first
    # max → ties broken in favour of the lighter repair.
    passing = [(i, e) for i, (e, _) in enumerate(rung_results)
               if (not thresholds) or e["passed"]]
    if passing:
        idx, _ = max(passing, key=lambda pair: _alpha(pair[1]))
        return idx, "best-alpha-mean (gate-passing)"
    # No passing rung: pick the best alpha across all.
    idx, _ = max(
        ((i, e) for i, (e, _) in enumerate(rung_results)),
        key=lambda pair: _alpha(pair[1]),
    )
    return idx, "best-alpha-mean (no rung passed gate)"


def _dict_to_check(d: dict):
    """Reconstruct a ThresholdCheck from its dict form (for the
    formatter, which expects the dataclass)."""
    from fvcom_mesh_tools.quality import ThresholdCheck

    return ThresholdCheck(
        metric=d["metric"], op=d["op"],
        threshold=d["threshold"], actual=d["actual"],
        passed=bool(d["passed"]),
    )


if __name__ == "__main__":
    sys.exit(main())
