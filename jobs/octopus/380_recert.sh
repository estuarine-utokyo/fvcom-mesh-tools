#!/bin/bash
#PBS -q OCT-S
#PBS --group=G16445
#PBS -l cpunum_job=8
#PBS -l memsz_job=64GB
#PBS -l elapstim_req=02:00:00
#PBS -N fmesh_380
#PBS -j o
#PBS -o logs/380_recert.pbs.log
#PBS -r n
#============================================================================
# OCTOPUS port of notebooks/380_junction_recert.pjsub (GENKAI).
#
# Regenerates and finishes the sample-repro mesh, then runs the QA gate,
# connectivity comparator, one-wide ledger, over-resolution/wall and
# OW05 checks and the standing comparison figures. The "342 on PREVIOUS
# mesh" gate-validation step of the GENKAI job is omitted.
#
# Migration check: the last GENKAI result (run 6218996: 3393 nodes,
# 5849 elements, QA 21/21, implied dt 16.32 s) is kept in
# outputs/sample_repro/ref_genkai_6218996/ and compared at the end.
#
# Usage (from the repository root):
#   qsub jobs/octopus/380_recert.sh
# Log: logs/380_recert.<jobid>.log
#============================================================================
set -euo pipefail
cd "${PBS_O_WORKDIR:-$PWD}"
. jobs/octopus/common.sh 380_recert 4

export SR_NORMALIZE=on
export SR_OBC_H1=1680
REF=outputs/sample_repro/ref_genkai_6218996

echo "=== unit tests ==="
python -m pytest -p no:cacheprovider tests/test_waterways.py \
    tests/test_waterways_normalize.py tests/test_channel_arcs.py -q
echo "=== generate (junction bridges + normalize on, edits kept) ==="
python notebooks/325_sample_repro.py
echo "=== finish ==="
python notebooks/331_finish2.py
echo "=== QA final (non-fatal for measurement) ==="
fmesh-mesh-qa outputs/sample_repro/sample_repro_final.14 || true
echo "=== connectivity comparator (incl. mesh-level probe) ==="
python notebooks/342_connectivity_check.py || true
echo "=== one-wide ledger ==="
python notebooks/346_one_wide_flags.py || true
echo "=== over-resolution / wall integrity ==="
python notebooks/364_over_resolution.py || true
echo "=== OW05 axis-centering check ==="
python notebooks/372_axis_check.py || true
echo "=== comparison figures (standing deliverable) ==="
python notebooks/369_ow05w.py || true
python notebooks/367_ow05_detail.py || true
python notebooks/376_problem_sites_3way.py || true
python notebooks/378_full_domain_3way.py || true

echo "=== migration check: GENKAI 6218996 vs this run ==="
for f in sample_repro.14 sample_repro_final.14; do
    if cmp -s "${REF}/${f}" "outputs/sample_repro/${f}"; then
        echo "${f}: byte-identical to GENKAI"
    else
        echo "${f}: DIFFERS from GENKAI"
    fi
done
fmesh-mesh-quality "${REF}/sample_repro_final.14" \
    outputs/sample_repro/sample_repro_final.14 || true

echo "end=$(date -Is)"
