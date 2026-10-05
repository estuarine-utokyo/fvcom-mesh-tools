# FVCOM run directories longer than 80 bytes (2026-10-06)

FVCOM kept `INPUT_DIR` and `OUTPUT_DIR` in `CHARACTER(LEN=80)`, and the buffers that join them to a file
name were 80..500 characters. A longer run directory was cut without a word and the run failed with
"FILE ... NOT FOUND" on the cut path, so every tool here carried an 80-byte guard and runs needed a short
root. Fortran has no such limit; the 80 was FVCOM's own declaration.

## What changed

FVCOM (`uk-fabm/v5.1.0-dev`, commits 52165c4e, 7aed2187, 04e7064e, ee0c8193, 2fcf9ee4):

- the path variables are `CHARACTER(LEN=1024)`: `INPUT_DIR`, `OUTPUT_DIR`, the joined-path buffers
  (`pathnfile`, `fname`, ...), `NCF%FNAME`, `GRID%NAME` (the file identifiers of history, average and
  surface output are full paths), `INFOFILE` and the command-line buffer (an absolute `--logfile`), and the
  absolute-path entries `FABM_YAML_FILE`, `MS_*_FILE`, `MUD_INITIAL_TEMP_FILE`;
- `Setup_Sed` / `Setup_Fluid_Mud` take `CHARACTER(LEN=*)`;
- `FIND_FILE_BYNAME` first looks for the whole name or the last path component, then falls back to the old
  substring match (a file registered as `D/pressure_wind.nc` was returned for `wind.nc`).

This repository: `FVCOM_DIR_MAX` is 800 (a 160-character `*_FILE` entry or the generated names must still
fit in 1024); notebooks 427, 449 and 454 use the shared `check_fvcom_dirs`. **A directory over 80 bytes needs
an FVCOM built from these commits**; an older binary still cuts at 80, so keep short roots where one is in use
(OCTOPUS until its FVCOM is updated).

## Verification

`jobs/genkai/fvcom_longpath_test.sh` runs one 2-minute case with the old and the new binary in a 71-byte and
a 426-byte directory (absolute `INPUT_DIR`/`OUTPUT_DIR`, history, average, surface and restart output, an
absolute 426-byte `--logfile`). The comparison covers every variable of the five files (schema, shape, type,
mask, NaN, scalars; only `file_date` is exempt): max |diff| 0 in both directories; the old binary fails in the
long one. GENKAI job 7005133.

## Review (gpt-6.1-sol, five rounds; `scratch/reviews/fvcom-pathlen/`, not tracked)

Fixed: file identifiers sharing a 160-character prefix, the directory test file, sediment dummies widened
against an 80-character actual, FABM/subdomain/Lagrangian/nesting buffers, absolute `--logfile`, direct
absolute paths in FABM/multi-sigma/mud temperature, the substring lookup, the guard margin, and the test's
comparator and argument handling.

Known limits, not changed (independent of the buffer length, or outside the build and use here):
EnKF/RRK/ETKF and `SEMI_IMPLICIT` buffers; `"./"//INPUT_DIR` joins in optional readers (dam, par, vegetation,
sediment, bbl, wqm: an absolute directory turns into `.//abs`); subdomain output file names (average ignores
`OUTPUT_DIR`, history and surface share `_0001.nc`); NetCDF error-message buffers (120); a restart base name
over 80 characters; two forcing files with the same base name in different sub-directories are still
ambiguous for `FIND_FILE_BYNAME`; relocation of a case onto a destination with symlinked children.
