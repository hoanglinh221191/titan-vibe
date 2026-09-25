# Release scope

Version 1.0.3 includes the focused VIBE static reference library in `src/`,
the runnable complete TITAN host under `source/`, and `titan.exe` at the
repository root. The host supplies input, topology parsing, output, logging,
PBCI fallback and other TITAN analysis functions, including the
Poisson-Boltzmann and other source needed by its build. HoloLift implementation
files and runners are excluded from the host source and the linked executable.
The optional VIBE H9 diagnostic export does not include the HoloLift solver.

The eight-contact cap is applied before retained-contact translation hypotheses
are built. Contact support and ambiguity are conditional on that retained set,
not on every contact that could be found in the input.

The top-level CMake target builds the focused static library; the Makefile under
`source/` builds the runnable host. The README gives commands for prebuilt
execution and for rebuilding the host. Large molecular frame subsets are in the
separate data deposit; the 480 compact held-out inputs are included here.
