# TITAN-VIBE

TITAN-VIBE is the reference implementation of Virtual Intact Block Extraction
(VIBE), a framewise method for reconstructing molecular assemblies across
periodic boundaries. It supports orthorhombic and triclinic cells, uses an
explicit molecular topology when available, and reports ambiguity or an
incomplete bounded search instead of silently presenting every placement as a
unique solution.

This first public release is intentionally limited to VIBE. It is a focused
source extraction from TITAN, not a release of the complete TITAN program.
HoloLift, Poisson-Boltzmann solvers, boundary-element solvers, PBCN/FixBox, and
unrelated analysis modules are not included.

## Method scope

VIBE determines relative periodic images for connected molecular parts in one
simulation frame. It also reconstructs within-molecule continuity needed to
define those parts. The method does not infer thermodynamic stability and does
not claim an unlimited search over every lattice translation.

The implementation combines:

1. topology bonds, or chain and residue metadata when topology is unavailable;
2. nearest-copy calculations in a triclinic lattice;
3. virtual-box and component-placement searches over a declared finite domain;
4. continuity, steric, and required-contact checks after placement;
5. explicit reports for supported, ambiguous, weakly observed, and incomplete
   outcomes.

Exactness is conditional on the reported search domain and completed numerical
checks. Large component assignments use a deterministic restricted search and
are not labeled exhaustive.

## Repository layout

- `src/pbctopo/`: VIBE geometry, search, scoring, certification, and reporting;
- `src/solv_pbc_pbctopo.cpp`: frame-level adapter used by TITAN;
- `src/mol_topology.*`: GROMACS topology graph construction;
- `src/MOL.h`, `src/Multi_MOL.h`, and related headers: minimal host declarations
  required to compile the published adapter;
- `test/`: independent algorithm-level smoke tests;
- `examples/vibe.arg.in`: example settings for integration in TITAN.

See `docs/SCOPE.md` for the exact inclusion boundary.

## Build

Requirements are CMake 3.20 or newer, a C++23 compiler, and an x86-64 processor
with AVX2 support. GCC 13 or newer and recent Clang releases are supported. The
AVX2 requirement matches the native build used for the reported calculations;
the individual mathematical tests retain scalar paths where applicable.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 8
```

The build creates a static reference library named `titan_vibe`. The adapter
retains calls supplied by the host application for trajectory input, output,
and logging; this repository therefore does not create a standalone `titan`
executable.

## Tests

On Windows PowerShell 7 with `g++` available:

```powershell
./test/run_core_smokes.ps1 -Compiler g++.exe
```

Each smoke compiles the required VIBE source files directly and checks one or
more mathematical or reporting invariants. Large trajectories and confirmatory
publication data are kept outside Git and identified by manifests and hashes.

## Reproducibility

See `docs/REPRODUCIBILITY.md` for release checks, interpretation of reported
states, and the boundary between source tests and manuscript evidence.

## Citation and license

Citation metadata are provided in `CITATION.cff`. TITAN-VIBE is distributed
under the MIT License; see `LICENSE`.
