#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=8
#PBS -l memsz_job=64GB
#PBS -l elapstim_req=01:00:00
#PBS -N fmesh_386
#PBS -j o
#PBS -o logs/386_issue_map.pbs.log
#PBS -r n
#============================================================================
# Re-measure the CERTIFIED mesh and map every flagged site.
#
# The 385 ablation left the per-run registries (one_wide_cells.json,
# ow_registry.json, normalize.json) describing the edit-free mesh while the
# .14 files were restored to the certified ones. This job reruns the
# comparators on the certified mesh, so the registries and the log match it
# again, then draws the issue map and zoom panels.
#
# Usage (from the repository root): qsub jobs/octopus/386_issue_map.sh
#============================================================================
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
. jobs/octopus/common.sh 386_issue_map 8

export SR_NORMALIZE=on
export SR_OBC_H1=1680

echo "=== connectivity comparator (certified mesh) ==="
python notebooks/342_connectivity_check.py || true
echo "=== one-wide ledger ==="
python notebooks/346_one_wide_flags.py || true
echo "=== over-resolution / wall integrity ==="
python notebooks/364_over_resolution.py || true
echo "=== issue map ==="
python notebooks/386_issue_map.py
echo "end=$(date -Is)"
