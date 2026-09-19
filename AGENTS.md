# AGENTS.md

Instructions for AI coding agents other than Claude Code (e.g. Codex)
working in this repository. Read `CLAUDE.md` as well: its project
overview, license policy, code organization and conventions apply to
you unchanged. This file adds the rules that `CLAUDE.md` inherits from
the owner's global settings, plus the rules for delegated work.

## Language

- Source-code comments, docstrings, documentation, commit messages and
  job scripts: English.
- Your final report is read by the Claude Code session that delegated
  the task; write it in English unless the task says otherwise.

## Environment (OCTOPUS, Osaka University)

- Batch system: NQSV. Submit with `qsub`, check with `qstat`. Job
  scripts live in `jobs/octopus/` and source `jobs/octopus/common.sh`
  (logs go to `logs/<name>.<jobid>.log`). Copy an existing script as a
  template. The `notebooks/*.pjsub` files are for GENKAI (`pjsub`) and
  do not run here.
- The login node is shared: anything heavier than a few seconds of CPU
  (mesh generation, finishing, sweeps, figure batches) must run as a
  batch job, never directly. Unit tests (`pytest -q`, ~30 s) and
  `ruff check` are fine on the login node.
- Inside the Codex sandbox, `qsub`/`qstat` are blocked (the NQSV
  wrapper is not executable there). Do not try to work around it:
  write or edit the job script, and put the exact `qsub` command in
  your report. The delegating session submits it and reads the log.
- Compute nodes have no network access.
- Conda env: `oceanmesh-bench` under
  `/octfs/work/G16445/v61021/miniforge3`. Activate with
  `. /octfs/work/G16445/v61021/miniforge3/etc/profile.d/conda.sh &&
  conda activate oceanmesh-bench`.
- Dependencies come from conda-forge only (`mamba install -c
  conda-forge ...`). Do not use `pip install` except
  `pip install -e . --no-deps --no-build-isolation` for the owner's own
  local repositories. Do not add dependencies without saying so in
  your report.
- Shared data: `$DATA_DIR` = `/octfs/work/G16445/share/Data`
  (read-only for the task; never write or delete there).
- Sibling repositories under `/octfs/work/G16445/v61021/Github/`
  (`oceanmesh` fork, `xcoast`, `TB-FVCOM`, ...) are read-only unless
  the task names them.

## Rules for delegated tasks

- Stay inside the scope of the task. If the task turns out to need a
  change outside it, stop and report instead of doing it.
- Do not commit, push, create branches, or rewrite git history. Leave
  your changes uncommitted in the working tree; the owner approves
  commits.
- Do not delete or overwrite files under `outputs/` that you did not
  create in this task.
- Verify with measurement: when a change affects a mesh, state the
  numbers (QA gates, node/element counts, implied dt) from an actual
  run, or say explicitly that it was not run.
- Final report: what you changed (files), how you verified it
  (commands and results, including failures), and anything left open.
