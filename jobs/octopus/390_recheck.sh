#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=16
#PBS -l memsz_job=64GB
#PBS -l elapstim_req=01:00:00
#PBS -N fmesh_390
#PBS -j o
#PBS -o logs/390_recheck.pbs.log
#PBS -r n
#============================================================================
# Re-measure the certified mesh with the CORRECTED detectors and draw both
# maps: the flagged-site map (386) and the new resolvable-unmeshed-water map
# (389). Assumes outputs/sample_repro/ holds the self-consistent certified
# set (job 115215 regenerated it).
#
# Usage (from the repository root): qsub jobs/octopus/390_recheck.sh
#============================================================================
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
. jobs/octopus/common.sh 390_recheck 16

export SR_NORMALIZE=on SR_OBC_H1=1680
echo "=== 364 (corrected WALL criterion) ==="
python notebooks/364_over_resolution.py || true
export FMESH_OVERWRITE=1
echo "=== 389 resolvable unmeshed water (new detector) ==="
python notebooks/389_coverage_gaps.py || true
echo "=== 386 issue map ==="
python notebooks/386_issue_map.py
echo "end=$(date -Is)"
