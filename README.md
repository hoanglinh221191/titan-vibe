# TITAN-VIBE 1.0.3 with executable and host source

This archive provides `titan.exe`, the complete buildable TITAN host source
needed for VIBE, the focused VIBE static library, and runnable reconstruction
checks. It is the software used for the VIBE methods article.

The host is exported from TITAN commit
`50ed41edecfe7b9cc5b72fec1ed5d38827d4c6ac`. The files in `src/pbctopo`,
`src/solv_pbc_pbctopo.cpp`, and `src/mol_topology.cpp` of the focused library
match that commit after line-ending normalization. The library keeps five
adapted host headers; the executable uses their complete TITAN counterparts.
HoloLift source and its runners are omitted, and no HoloLift objects are
linked. Other TITAN host functions remain.

The binary reports `0.2.0-dev.304+g50ed41edecfe.dirty`. `dirty` marks the
documented source export (HoloLift excluded, build metadata pinned), not a
change to the VIBE algorithm. This is the embedded TITAN host identity; the
TITAN-VIBE distribution version is 1.0.3. See `SOURCE_AUDIT.json` and
`titan_build_manifest.json` for exact identities, and `RELEASE_NOTES_v1.0.3.md`
for what changed since 1.0.2.

## Run without a compiler

Use Windows x86-64 on an AVX2/FMA-capable processor. The binary was built with
Clang 22.1.4 and the portable profile; it depends only on Windows/UCRT system
DLLs. The replay drivers need Python 3.10 or newer and only the standard
library. Run commands from the extracted package directory in PowerShell:

```powershell
.\titan.exe --version
python .\verify_package.py
python .\replay.py --output .\my_smoke
python .\replay.py --all --output .\my_480
```

The first reconstruction command selects 20 cases spanning every stratum,
topology source, and component count; `--all` runs all 480 held-out inputs in
separate directories, four one-thread processes at a time by default (`--jobs 1`
for serial execution). The driver refuses an existing output directory. It
compares each output PDB SHA-256 and the scientific report fields with the
reference outputs of this release and returns a nonzero exit code on
disagreement; the name of the search route (`mode`) is recorded but is not a
pass criterion. `REPLAY_RESULT.json` records the binary hash, counts, and
per-case findings.

This is an implementation replay of fixed inputs, not a new held-out experiment
or an independent physical ground-truth assessment. The original endpoint table
of the held-out confirmation is kept as
`examples/heldout480/ARCHIVED_CASE_ENDPOINTS.csv`.

## Reproduce the five molecular systems from the deposited data

Obtain and extract the September 11 data update from
https://doi.org/10.5281/zenodo.22710653 and point `--data-root` at the directory
`vibe_data_update_20260911` containing `INPUT_PATH_MAP.json`,
`input_dependencies`, and `jcim_revision_controls_20260911`.

```powershell
python .\replay_real.py --data-root 'D:\review\vibe_data_update_20260911' --output .\my_five_frames
python .\replay_real.py --data-root 'D:\review\vibe_data_update_20260911' --all-frames --output .\my_367_frames
```

The first command runs the first frame of each of the five systems; the second
runs all 367 comparison frames (102 tetramer, 101 hexamer, 102 protein-glycan,
31 peptide-pair, and 31 phase-separated frames). The driver uses the deposited
full-atom frame subsets, source PDBs, index selections, and topology include
trees; no upstream trajectory or GROMACS installation is needed for this stage.
It writes a new `arg.in` for each run without editing archived files, and
compares atom records at PDB precision and the scientific report fields for
every frame.

These commands reproduce the VIBE reconstructions. Reproducing the MDVWhole and
GROMACS results, timing measurements, input extraction, or the independent
oracle analysis uses the archived protocols in the data deposit and their own
software.

For other inputs, edit a copy of an argument file and run:

```powershell
.\titan.exe -a 'D:\review\new_run\arg.in'
```

The VIBE option is `-pbctopo=true`; the phase-separated mode is
`-pbctopo-mode=llps-local`. Output files are written in the working directory;
use a new working directory for each run.

## Build

The focused library builds with CMake 3.20+ and a C++23 compiler:

```powershell
cmake -S . -B build -G Ninja
cmake --build build
```

The host was built with MSYS2 UCRT64 Clang 22.1.4, GNU make, LLVM LLD, GNU
windres, Python, and the toolchain's OpenMP and static runtime libraries. Eigen
and trajectory I/O source are included. In PowerShell, with that toolchain at
`C:\msys64\ucrt64\bin`:

```powershell
$env:PATH = 'C:\msys64\ucrt64\bin;' + $env:PATH
mingw32-make -C source PROFILE=portable -j 6 all
.\source\titan.exe --version
python .\replay.py --binary .\source\titan.exe --all --output .\rebuilt_480
```

The Makefile pins the source identity and marks the HoloLift-excluded build as
modified. Its Linux rules are included, but no Linux executable of this release
has been built or validated. Do not assume byte-identical binaries across
toolchains; use the reconstruction checks above.

## Retained contacts

Candidate contacts are sorted and truncated to eight per pair of molecular
units in `src/pbctopo/anchor_graph.cpp` before translation hypotheses are built,
so contact support and ambiguity refer to that retained set. In the article's
tests, builds retaining 32 or 1,024 pairs gave the same placements and contact
interpretations. The retained counts are not interface sizes.

## Provenance and licenses

`SHA256SUMS.txt` covers the packaged files, and `validation` holds this
release's validation records. `source/LICENSE` is the TITAN license;
third-party notices in the included source and force-field files are retained.
