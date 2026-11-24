# Reproducibility record for 1.2.0

The executable and the source under `source/` come from TITAN commit
`994ffa0ad9cd69488643ec8a1313e816ff4ef862`. `source/SOURCE_PROVENANCE.json` records how the VIBE-only source
tree was exported; every exported file equals the file of that commit, apart
from line endings and the pinned build metadata in `source/Makefile`. The
executable hash and the imported DLLs are recorded in `SOURCE_AUDIT.json`.

The embedded TITAN version is `0.2.0-dev.467+g994ffa0ad9cd.dirty`. Its suffix records the
source export and the pinned metadata. The TITAN-VIBE distribution version
1.2.0 does not relabel that identity.

`validation/` records the checks shipped with the package:
- `QUICKSTART_RESULT.txt`: the quick-start examples against their
  `EXPECTED.json`.
- `HELDOUT480_REPLAY.json`: the 480 held-out inputs against the references
  (byte-identical PDBs and equal report fields).
- `REAL367_REPLAY.json` and `FIRST_FRAME5_REPLAY.json`: the five molecular frame
  sets, compared with the deposited coordinates and with the report fields of
  `validation/REAL367_FIELDS_v1.1.0.csv`.
- `LINUX_*`: the same checks for a static Linux build of `source/` with GCC
  13.3.0 (`LINUX_BUILD.txt`); the held-out replay ran with
  `--coordinate-tolerance 0.001`.

The original held-out endpoint table is kept in
`examples/heldout480/ARCHIVED_CASE_ENDPOINTS.csv`.

The designed-case oracle, the MDVWhole and GROMACS comparisons, and the timing
protocol are archived in the data deposit (DOI 10.5281/zenodo.22710653 and its
later versions). A replay establishes that outputs agree. It is not
independent physical truth or a new held-out study.
