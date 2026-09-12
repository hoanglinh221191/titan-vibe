# TITAN-VIBE 1.0.2 with executable and host source

This companion archive provides `titan.exe`, the complete buildable TITAN host
source needed for VIBE, and runnable reconstruction checks. It supplements the
VIBE v1.0.1 static-library release at https://doi.org/10.5281/zenodo.22669272.
Version 1.0.1 and its DOI remain available as the preceding static-library release.

The host is exported from TITAN commit
`6a8385cf1c63c42046b8a38053f003311e85e633`, the source used for the September 11
submission controls. All `src/pbctopo` files, `solv_pbc_pbctopo.cpp`, and
`mol_topology.cpp` match the published VIBE algorithm after line-ending
normalization. The public library has five adapted host headers; the executable
uses their complete TITAN counterparts. HoloLift source and its standalone runner
are omitted, and no HoloLift objects are linked. Other TITAN host functions remain.

The binary reports `0.2.0-dev.242+g6a8385cf1c63.dirty`: `dirty` identifies the
documented source-export/build changes, not a change to the VIBE algorithm. This
is the embedded TITAN host identity; the TITAN-VIBE distribution version is 1.0.2.
The two version identifiers describe different scopes. See `SOURCE_AUDIT.json`
and `titan_build_manifest.json` for exact identities.

## Run without a compiler

Use Windows x86-64 on an AVX2/FMA-capable processor. The binary was built with
Clang 22.1.4 and the portable profile; it depends only on Windows/UCRT system
DLLs. Python 3.10 or newer is used for the replay drivers; they require only the
standard library. Paths with spaces are supported. Run commands from the
extracted package directory in PowerShell:

```powershell
.\titan.exe --version
python .\verify_package.py
python .\replay.py --output .\my_smoke
python .\replay.py --all --output .\my_480
```

The first reconstruction command selects 20 cases spanning every available
stratum/topology/component-count combination. `--all` runs all 480 held-out
inputs in separate directories. Four independent one-thread processes run by
default; use `--jobs 1` for serial execution. The driver refuses an existing
output directory. It compares each output PDB SHA-256 and 18 scientific report
fields with the archived results and returns a nonzero exit code on disagreement.
`REPLAY_RESULT.json` records the binary hash, counts, and per-case findings.

This is an implementation replay of fixed inputs, not a newly held-out experiment
or a new independent physical ground-truth assessment. The original scientific
endpoint table is retained as `examples/heldout480/ARCHIVED_CASE_ENDPOINTS.csv`.

## Reproduce the five molecular systems from the deposited data

Obtain and extract the September 11 data update from
https://doi.org/10.5281/zenodo.22710653. Point `--data-root` at the directory
`vibe_data_update_20260911` containing `INPUT_PATH_MAP.json`,
`input_dependencies`, and `jcim_revision_controls_20260911`.

```powershell
python .\replay_real.py --data-root 'D:\review\vibe_data_update_20260911' --output .\my_five_frames
python .\replay_real.py --data-root 'D:\review\vibe_data_update_20260911' --all-frames --output .\my_367_frames
```

The first command runs the first frame of each of the five systems. The second
runs all 367 frames used in the comparison: 102 tetramer, 101 hexamer, 102
protein-glycan, 31 peptide-pair, and 31 phase-separated frames. The driver uses
the deposited full-atom frame subsets, source PDBs, index selections, and topology
include trees; no complete upstream trajectory or GROMACS installation is needed
for this reconstruction stage. It adapts input paths into a new `arg.in` without
editing archived inputs or outputs. It compares atom records/coordinates at PDB
precision and the same 18 scientific report fields for every frame. See
`REAL_REPLAY_RESULT.json` for completion counts and comparisons.

These commands reproduce the VIBE reconstruction endpoints. Reproducing
MDVWhole/GROMACS results, original timing measurements, input extraction, or
the independent oracle analysis additionally uses the corresponding archived
protocols and their own software requirements. Historical results and timing
tables are not replaced by the companion-build checks.

For other inputs, edit a copy of an argument file and run:

```powershell
.\titan.exe -a 'D:\review\new_run\arg.in'
```

The VIBE option is `-pbctopo=true`; the prescribed phase-separated mode is
`-pbctopo-mode=llps-local`. Output files are written in the working directory.
Use a new working directory for each run.

## Build the supplied host source

The tested Windows toolchain is MSYS2 UCRT64 Clang 22.1.4 with GNU make, LLVM LLD,
GNU windres, Python, and the toolchain's OpenMP/static runtime libraries. Eigen
and trajectory I/O source are included. In PowerShell, with that toolchain at
`C:\msys64\ucrt64\bin`:

```powershell
$env:PATH = 'C:\msys64\ucrt64\bin;' + $env:PATH
mingw32-make -C source PROFILE=portable -j 6 all
.\source\titan.exe --version
python .\replay.py --binary .\source\titan.exe --all --output .\rebuilt_480
```

The Makefile pins source identity and marks the HoloLift-excluded build as
modified. Its Linux build rules are also included, but this companion's native
Linux executable has not been built or validated. Do not assume byte-identical
binaries across toolchains; use the reconstruction checks above.

## Interpretation of the eight-contact cap

In the released algorithm, candidate contacts are sorted and truncated to eight
per molecular-unit pair in `src/pbctopo/anchor_graph.cpp` before
`build_pbctopo_contact_hypotheses` receives the retained edges. Consequently,
contact support and ambiguity are conditional on that retained set. A relative
translation represented only by discarded contacts need not be detected. These
counts are not interface sizes; they are also distinct from the separate
all-heavy-atom post hoc counts in the molecular-system comparison.

## Provenance and licenses

`SHA256SUMS.txt` covers the packaged files. `validation` contains the companion
build's actual validation records. Original protocols, results, and public DOI
records remain identifiable separately. `source/LICENSE` is the TITAN license;
third-party notices in the included source and force-field files are retained.
