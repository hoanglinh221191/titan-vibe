# Release scope

Version 1.1.0 ships three things:
- `titan.exe` at the repository root;
- the runnable, complete TITAN host under `source/`, where VIBE is implemented
  in `source/src/pbctopo`;
- small examples under `examples/quickstart`.

The host supplies input, topology parsing, output, logging, PBCI fallback and
other TITAN analysis functions, including the Poisson-Boltzmann code and the
other source its build needs. HoloLift implementation files, tests and runners
are excluded from the host source and from the linked executable. The focused
VIBE static library of the 1.0.x releases is not shipped, because the VIBE
solver now spans several host directories.

The eight-contact cap is applied before the retained-contact translation
hypotheses are built. Contact support and ambiguity are therefore conditional
on that retained set, not on every contact in the input. The contact-cycle
lattice and the input-overlap counts use every pair within their cutoffs.

The Makefile under `source/` builds the host, and the README gives commands
for prebuilt execution and for rebuilding. The large molecular frame subsets
are in the separate data deposit. The 480 compact held-out inputs and the
quick-start snapshots are included here. The quick-start snapshots are
selections, with water and ions removed, from simulations of the methods
article; their topologies are reduced to the sections VIBE reads.
