# Reproducibility record for 1.0.3

The executable, the host source under `source/`, and the focused library under
`src/` all come from TITAN commit `50ed41edecfe7b9cc5b72fec1ed5d38827d4c6ac`.
Source normalization, adapted host headers, the executable hash, and imported
DLLs are recorded in `SOURCE_AUDIT.json`.

The embedded host version is `0.2.0-dev.304+g50ed41edecfe.dirty`. The suffix
records the HoloLift-excluded export and pinned metadata. The TITAN-VIBE
distribution version 1.0.3 does not relabel that identity.

Every VIBE calculation of the methods article was run with this executable.
`validation/` records the checks shipped with the package. Replay of the 480
held-out inputs and of the five molecular frame sets compares outputs with the
reference outputs of this release. The original held-out endpoint table is kept
in `examples/heldout480/ARCHIVED_CASE_ENDPOINTS.csv`.

The designed-case oracle, the MDVWhole and GROMACS comparisons, and the timing
protocol are archived in the data deposit (DOI 10.5281/zenodo.22710653 and its
later versions). A replay establishes output agreement; it is not independent
physical truth or a new held-out study.
