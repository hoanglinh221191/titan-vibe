# Reproducibility record for 1.1.0

The executable and the host source under `source/` come from TITAN commit
`dcb2bb5ffbcc9cf27d843a55f53d31f765b01794`. The excluded HoloLift paths are listed in
`source/SOURCE_PROVENANCE.json`. The executable hash and the imported DLLs are
recorded in `SOURCE_AUDIT.json`.

The embedded host version is `0.2.0-dev.450+gdcb2bb5ffbcc.dirty`. Its suffix records the
HoloLift-excluded export and the pinned metadata. The TITAN-VIBE distribution
version 1.1.0 does not relabel that identity.

This is the executable for the VIBE calculations of the revised methods
article. `validation/` records the checks shipped with the package:
- `QUICKSTART_RESULT.txt`: the quick-start examples against their
  `EXPECTED.json`.
- `HELDOUT480_REPLAY.json`: the 480 held-out inputs against the references of
  this release.
- `REAL367_REPLAY.json` and `FIRST_FRAME5_REPLAY.json`: the five molecular frame
  sets, compared with the deposited coordinates and with the report fields of
  this release.

The changes relative to 1.0.3 are counted in `SOURCE_AUDIT.json`. The original
held-out endpoint table is kept in
`examples/heldout480/ARCHIVED_CASE_ENDPOINTS.csv`.

The designed-case oracle, the MDVWhole and GROMACS comparisons, and the timing
protocol are archived in the data deposit (DOI 10.5281/zenodo.22710653 and its
later versions). A replay establishes that outputs agree. It is not
independent physical truth or a new held-out study.
