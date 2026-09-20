"""Measure a lon/lat fort.14 mesh against a recipe; emit one CSV row per element.

Usage: python notebooks/391_sizing_report.py mesh.14 recipes/sizing/tokyo_bay.yaml

Constraint labels are a reconstruction, not historical builder provenance:
coastal_target + grade*distance to mesh land boundary, capped at maxel, raised
by the recipe CFL floor, overridden by regions, then graded on the element
adjacency graph. The legacy 325 zones, dilated floor and boundary corridors
cannot be recovered from a mesh and recipe. Measured dt uses minimum triangle
altitude and maximum vertex depth (positive-down metres). Run large meshes in
an NQSV batch job. CSV goes to stdout; caveats/global dt go to stderr.
"""
from __future__ import annotations

import argparse
import csv
import heapq
import sys
from pathlib import Path

import numpy as np
import shapely
from shapely.geometry import LineString
from shapely.ops import unary_union

from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.sizing import apply_sizing_regions, load_sizing


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mesh", type=Path)
    parser.add_argument("recipe", type=Path)
    args = parser.parse_args()
    mesh, cfg = read_fort14(args.mesh.resolve()), load_sizing(args.recipe.resolve())
    ll, tri = mesh.nodes, mesh.elements
    if not len(tri) or np.any(abs(ll[:, 0]) > 180) or np.any(abs(ll[:, 1]) >= 90):
        parser.error("mesh must contain elements and use geographic lon/lat coordinates")
    if not np.isfinite(mesh.depths).all() or np.any(mesh.depths < 0):
        parser.error("mesh depths must be finite positive-down metres")
    xy = (ll - ll.mean(axis=0))*[111000*np.cos(np.deg2rad(ll[:, 1].mean())), 111000]
    centres = ll[tri].mean(axis=1)
    cxy = xy[tri].mean(axis=1)
    lines = [LineString(xy[ids]) for _, ids in mesh.land_boundaries if len(ids) >= 2]
    if not lines:
        parser.error("land boundary segments are required for coastal distance")
    distance = shapely.distance(shapely.points(cxy), unary_union(lines))
    coast = cfg["coastal_target_m"] + cfg["gradation"]*distance
    cap = cfg["max_edge_length_m"]
    base = np.minimum(coast, cap)
    depth = mesh.depths[tri].max(axis=1)
    speed = np.sqrt(9.81*depth)
    floor = speed*cfg["cfl"]["dt_s"]/cfg["cfl"]["cr"]
    labels = np.where(coast < cap, "coastal_target", "max_edge_length").astype(object)
    labels[floor > base] = "cfl_floor"
    base = np.maximum(base, floor)
    # A huge slope disables the lattice limiter for these unstructured centroids;
    # mesh adjacency below supplies the physically meaningful gradation graph.
    target, _ = apply_sizing_regions(
        base[:, None], centres[:, 0, None], centres[:, 1, None], cfg["regions"],
        gradation=1e100, hmin_m=cfg.get("hmin_m"))
    target = target[:, 0]
    labels[target != base] = "region_target"
    if "hmin_m" in cfg:
        labels[target == cfg["hmin_m"]] = "hmin"
    neighbours = [[] for _ in tri]
    edges = {}
    for i, nodes in enumerate(tri):
        for a, b in ((nodes[0], nodes[1]), (nodes[1], nodes[2]), (nodes[2], nodes[0])):
            key = tuple(sorted((int(a), int(b))))
            if key in edges:
                j = edges[key]
                cost = cfg["gradation"]*float(np.linalg.norm(cxy[i]-cxy[j]))
                neighbours[i].append((j, cost))
                neighbours[j].append((i, cost))
            else:
                edges[key] = i
    queue = [(float(v), i) for i, v in enumerate(target)]
    heapq.heapify(queue)
    while queue:
        value, i = heapq.heappop(queue)
        if value > target[i]:
            continue
        for j, cost in neighbours[i]:
            if value+cost < target[j]:
                target[j] = value+cost
                labels[j] = "gradation_limit"
                heapq.heappush(queue, (target[j], j))
    vertices = xy[tri]
    lengths = np.linalg.norm(vertices-np.roll(vertices, 1, axis=1), axis=2)
    a, b = vertices[:, 1]-vertices[:, 0], vertices[:, 2]-vertices[:, 0]
    altitude = np.abs(a[:, 0]*b[:, 1]-a[:, 1]*b[:, 0])/lengths.max(axis=1)
    measured_dt = np.divide(cfg["cfl"]["cr"]*altitude, speed,
                            out=np.full(len(tri), np.inf), where=speed > 0)
    print("Binding labels are reconstructed constraints, not historical provenance. "
          "Legacy zones/dilated floor/corridors are unavailable. "
          f"Measured global dt at recipe Cr: {measured_dt.min():.6g} s", file=sys.stderr)
    writer = csv.writer(sys.stdout)
    writer.writerow(["element", "reconstructed_binding", "coastal_target_m", "cfl_floor_m",
                     "reconstructed_size_m", "min_edge_m", "min_altitude_m", "depth_max_m",
                     "measured_dt_s", "measured_dt_cr1_s"])
    for i in range(len(tri)):
        writer.writerow([i+1, labels[i], coast[i], floor[i], target[i], lengths[i].min(),
                         altitude[i], depth[i], measured_dt[i], measured_dt[i]/cfg["cfl"]["cr"]])


if __name__ == "__main__":
    main()
