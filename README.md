# TITAN-VIBE 1.2.0

TITAN-VIBE reconstructs molecules and multi-chain assemblies across the
periodic boundaries of molecular dynamics frames, in orthorhombic and
triclinic cells, and reports for every frame whether the stored frame fixes
the relative placement of the parts. The VIBE methods article uses this
software.

This package provides:
- `titan-vibe.exe`, a Windows x86-64 executable;
- the source of TITAN-VIBE under `source/`, which builds `titan-vibe` on
  Windows and on Linux;
- small examples with several snapshots each;
- the 480 held-out inputs with their reference outputs, and replay drivers.

## What TITAN-VIBE is

VIBE is implemented in the TITAN program (`source/src/pbctopo` and its host
files). TITAN-VIBE is the VIBE-only build of TITAN: `source/` holds only the
570 files that this build reads (154 C++ source files, about 8 MB):
the structure, index, GROMACS topology and trajectory (xtc, trr) readers, the
force-field library that types atoms when no topology is given, VIBE, the
TITAN engine without its energy terms, and the PDB writer. None of TITAN's
energy, solvation, kinetics or analysis modules is included.
`source/src/vibe_only/README.md` explains how they are left out.

- `titan-vibe` accepts the VIBE options and the input, output, frame and
  thread options. Any other TITAN option is refused with a message, so no run
  can reach a module that is not built.
- `-h` and `-vibe-help` print the options and the meaning of every label.
- Every run writes coordinates, so `-nocal=true` is required.

The source is exported from TITAN commit `994ffa0ad9cd69488643ec8a1313e816ff4ef862`, and the executable
reports `0.2.0-dev.467+g994ffa0ad9cd.dirty`. Here `dirty` marks the documented source export (only
the VIBE-only files, build metadata pinned in `source/Makefile`), not a change
to the VIBE algorithm. This is the embedded TITAN identity; the TITAN-VIBE
distribution version is 1.2.0. Exact identities are in `SOURCE_AUDIT.json`
and `titan_build_manifest.json`. `RELEASE_NOTES_v1.2.0.md` lists what
changed since 1.1.0. On every check shipped here (the quick-start examples,
the 480 held-out inputs and the 367 frames of five molecular systems) the
outputs equal those of 1.1.0.

## Quick start

Write one option per line in an argument file, then run it. Each example
folder under `examples/quickstart` has its argument files; run them from that
folder, and the outputs are written there:

```powershell
cd examples\quickstart\01_peptide_pair
..\..\..\titan-vibe.exe -a arg.in
..\..\..\titan-vibe.exe -vibe-help      # the options below and the meaning of every label
```

On Linux, after building (see Build), the executable is
`../../../source/titan-vibe`.

### 1. A trajectory of a small assembly

`examples/quickstart/01_peptide_pair` holds 31 frames of two amyloid-beta
peptides in a triclinic cell (from the data deposit). `arg.in`:

```
-ff=charmm36
-path-a=../../../forcefields
-path-n=../../../forcefields
-s=reference.pdb
-traj=frames.xtc
-n=index.ndx
-g=0-0
-common_topo=topol.top
-pbctopo=true
-1pdb=true
-nocal=true
-no-dssp=true
-noprint-side_chain_contact=true
-opdb=placed.pdb
```

`placed.pdb` has one model per frame with both peptides whole and in
contact. `pbctopo_frame_report.csv` has one row per frame; in 28 frames
`evidence_state` is `supported_connected`, and in frames 0, 4 and 28 it is
`extended_periodic` with `contact_cycle_rank` 1: there a peptide touches two
copies of its partner (and in frame 28 also its own copy), so no finite
placement keeps every contact.

### 2. What the labels mean, on six designed cases

`examples/quickstart/02_labels` holds six small designed objects, built so
that the answer is known. Each has four wrapped frames at different positions
in the cell, and every frame gets the same label:

| case | object | evidence_state | contact_cycle_rank |
|---|---|---|---|
| finite_line | four molecules in a row | supported_connected | 0 |
| stack_rank1 | molecules that close a row through +a | extended_periodic | 1 |
| sheet_rank2 | molecules that tile the ab plane | extended_periodic | 2 |
| rod_self | one rod 4 A from its own image | extended_periodic | 1 |
| rod_pair | two rods, two gaps across the boundary | extended_periodic | 1 |
| ring_bonded | a bonded ring through the boundary | not_evaluated (fallback_reason bonded_winding, bonded_winding_rank 1) | - |

### 3. Two glycosylated chains

`examples/quickstart/03_glycoprotein_pair` holds five frames (10 ns apart) of
two protein chains with their N-glycans (14,393 atoms). The glycans are
bonded to their chains, so each chain with its glycans is one part. In every
frame both parts are written whole and in contact: `certified`,
`supported_connected`, `contact_cycle_rank` 0.

### 4. A condensate, local and global

`examples/quickstart/04_condensate_martini` holds six frames (60 ns apart) of
128 Martini 3 peptides, about 100 of them in one droplet.

- **Local mode.** `arg_local.in` uses `-pbctopo-mode=llps-local`, the mode for
  condensates and membranes. It makes every peptide whole where it is and does
  not place peptides relative to each other. The frame report still gives the
  rank of the contact network of the whole molecules. In frame 2 that rank is
  1 along a: the droplet touches its own image across the 202 Å cell. In the
  other frames it is 0. (A bilayer that spans the cell has rank 2.)
- **Global mode.** `arg_global.in` runs global VIBE on the same frames. Frame
  2 is labelled `extended_periodic`. The others are `weak_ambiguous`: with
  21-26 contact components, a frame does not fix where the free peptides sit
  relative to the droplet.

### Checking the examples

```powershell
python run_quickstart.py --output quickstart_runs
```

`run_quickstart.py` copies every example into a new directory, runs it, and
checks two things against the example's `EXPECTED.json`:
- the frame-report labels, frame by frame;
- the centroid of every molecule in the output PDB, within 0.01 Å.

### Reading the frame report

| evidence_state | meaning |
|---|---|
| supported_connected | contacts connect every part and every cycle of contacts closes: the frame fixes the placement up to one common lattice translation |
| weak_disconnected | some parts are joined to the rest by no contact path: their relative images are not fixed by the frame |
| extended_periodic | the contact network closes through the boundary: no finite placement keeps every contact; `contact_cycle_rank` and `contact_cycle_b*` give the periodic directions |
| weak_ambiguous | the contacts close, but the selected placement leaves a relation unsupported |
| contradicted | the selected placement loses an observed contact |

`input_overlap_pairs` counts atom pairs of two parts that are closer than
lambda (r_i + r_j) in the input. Under periodic boundary conditions such a
pair overlaps in the simulated system whatever placement is written, so it
flags the input, not the placement.

## Run without a compiler (Windows)

Use Windows x86-64 on an AVX2/FMA-capable processor. The executable
(7.8 MB) was built with Clang 22.1.4 and the portable profile, and it
depends only on Windows/UCRT system DLLs. The Python drivers need Python 3.10
or newer and only the standard library. Run these commands from the extracted
package directory in PowerShell:

```powershell
.\titan-vibe.exe --version
python .\verify_package.py
python .\run_quickstart.py --output .\my_quickstart
python .\replay.py --output .\my_smoke
python .\replay.py --all --output .\my_480
```

- `verify_package.py` checks every file against `SHA256SUMS.txt`.
- `run_quickstart.py` runs the examples above.
- `replay.py` runs 20 held-out cases that span every stratum, topology source
  and component count. With `--all` it runs all 480 held-out inputs, each in
  its own directory, four one-thread processes at a time by default (`--jobs 1`
  runs them one after another). It compares each output PDB SHA-256 and the
  scientific report fields with the reference outputs, and returns a nonzero
  exit code on any disagreement.
- `REPLAY_RESULT.json` records the binary hash, the counts and the findings for
  every case.
- Every driver refuses an existing output directory, and every driver takes
  `--binary` for another executable.

This is an implementation replay of fixed inputs. It is not a new held-out
experiment or an independent physical ground-truth assessment. The original
endpoint table of the held-out confirmation is kept as
`examples/heldout480/ARCHIVED_CASE_ENDPOINTS.csv`. The references are the
outputs of release 1.1.0, which this release reproduces; the 1.0.3 references
are kept in `examples/heldout480/manifest_v1.0.3.json`.

## Build

The source needs a C++23 compiler with OpenMP and GNU make; Eigen and the
trajectory I/O source are included. `source/Makefile` builds `titan-vibe` by
default (`VIBE_ONLY ?= 1`) and pins the source identity.

### Windows

The executable of this package was built with MSYS2 UCRT64 Clang 22.1.4, GNU
make, LLVM LLD, GNU windres, Python, and the toolchain's OpenMP and static
runtime libraries. In PowerShell, with that toolchain at
`C:\msys64\ucrt64\bin`:

```powershell
$env:PATH = 'C:\msys64\ucrt64\bin;' + $env:PATH
mingw32-make -C source -j 6 all
.\source\titan-vibe.exe --version
python .\run_quickstart.py --binary .\source\titan-vibe.exe --output .\rebuilt_quickstart
python .\replay.py --binary .\source\titan-vibe.exe --all --output .\rebuilt_480
```

### Linux

Tested with GCC 13.3.0 on Ubuntu 24.04.1 LTS under WSL2 (x86-64), statically linked:

```bash
make -C source CXX=g++-13 LDFLAGS="-static -fopenmp" -j 8 all
./source/titan-vibe --version
python3 run_quickstart.py --output rebuilt_quickstart
python3 replay.py --all --coordinate-tolerance 0.001 --output rebuilt_480
```

`make -C source -j 8 all` builds with the system `g++` (GCC 13 or newer) an
executable that loads only the C library at run time; it was also tested on
the quick-start examples. On Linux the drivers use
`source/titan-vibe` by default.

GCC and Clang can round a coordinate differently in its last bit, and a PDB
prints coordinates to 0.001 Å. The Linux build gives the quick-start results
of `EXPECTED.json` and the 367 frames of the five molecular systems at PDB
precision (`replay_real.py`, below). Of the 480 held-out outputs,
459 are byte-identical to the references and 21
differ in some coordinates by 0.001 Å, one unit of the last printed digit,
with the same labels and report fields. `--coordinate-tolerance 0.001`
accepts exactly this: the same records in the same order, coordinates within
0.001 Å of `examples/heldout480/<case>/expected_oracle_v2_pbctopo.pdb`, and
every report field equal except `wrapped_coordinate_hash`, a checksum of the
wrapped input coordinates at full precision, which differs between the two
compilers in every case. Do not assume byte-identical outputs across
toolchains; use these checks.

## Reproduce the five molecular systems from the deposited data

Obtain the September 11 data update from
https://doi.org/10.5281/zenodo.22710653 and extract it. Point `--data-root` at
the directory `vibe_data_update_20260911`, which contains
`INPUT_PATH_MAP.json`, `input_dependencies` and
`jcim_revision_controls_20260911`:

```powershell
python .\replay_real.py --data-root 'D:\review\vibe_data_update_20260911' --output .\my_five_frames
python .\replay_real.py --data-root 'D:\review\vibe_data_update_20260911' --all-frames --output .\my_367_frames
```

- The first command runs the first frame of each of the five systems.
- The second runs all 367 comparison frames: 102 tetramer, 101 hexamer, 102
  protein-glycan, 31 peptide-pair and 31 phase-separated frames.
- The driver uses the deposited frame subsets, source PDBs, index selections
  and topology include trees. No upstream trajectory or GROMACS installation
  is needed.
- It writes a new `arg.in` for each run, without editing archived files.
- It compares the atom records with the deposited outputs at PDB precision,
  and the scientific report fields with
  `validation/REAL367_FIELDS_v1.1.0.csv`. It also counts the fields that
  differ from the deposited 1.0.3 reports.

These commands reproduce the VIBE reconstructions. The MDVWhole and GROMACS
results, timing measurements, input extraction and the independent oracle
analysis use the archived protocols in the data deposit and their own software.

For other inputs, edit a copy of an argument file and run:

```powershell
.\titan-vibe.exe -a 'D:\review\new_run\arg.in'
```

Output files are written to the working directory, so use a new working
directory for each run. `.\titan-vibe.exe -vibe-help` lists the options and the
labels.

## Retained contacts

Candidate contacts are sorted and truncated to eight per pair of molecular
units in `source/src/pbctopo/anchor_graph.cpp` before translation hypotheses
are built, so contact support and ambiguity refer to that retained set. In the
article's tests, builds retaining 32 or 1,024 pairs gave the same placements
and contact interpretations. The retained counts are not interface sizes. The
contact-cycle lattice (`contact_cycle_*`) is computed from every heavy-atom
pair within the contact cutoff, and the input-overlap counts (`input_overlap_*`)
from every atom pair below the overlap threshold. Neither uses the retained
set.

## Provenance and licenses

- `SHA256SUMS.txt` covers the packaged files.
- `validation/` holds the validation records of this release, for the Windows
  executable and for the Linux build.
- `source/LICENSE` is the TITAN license. Third-party notices in the included
  source and force-field files are retained.
- `source/SOURCE_PROVENANCE.json` records how the source tree was exported.
- The example topologies are reduced to the sections VIBE reads (see each
  example's README), so they are not simulation topologies.
