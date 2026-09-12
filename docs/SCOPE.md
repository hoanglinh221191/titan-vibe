# Release scope

Version 1.0.2 includes both the focused VIBE static reference library in `src/`
and the runnable complete TITAN host under `source/`, with `titan.exe` at the
repository root. The host supplies input, topology parsing, output, logging,
PBCI fallback and existing TITAN analysis functions. Unlike 1.0.1, this release
includes the host Poisson-Boltzmann and other analysis source needed by its build.
HoloLift implementation files and its standalone runner are excluded from the
host source and linked executable. Optional VIBE H9 diagnostic exports do not
include the HoloLift solver.

The VIBE algorithm files have not changed from 1.0.1. The eight-contact cap is
applied before retained-contact translation hypotheses are constructed. Contact
support and ambiguity are conditional on that retained set, not on every contact
that could be found in the input.

The top-level CMake target still builds the focused static library. The Makefile
under `source/` builds the runnable host. The README provides commands for both
prebuilt execution and rebuilding the host. Large molecular frame subsets remain
in the separate data deposit; 480 compact held-out input cases are included.
