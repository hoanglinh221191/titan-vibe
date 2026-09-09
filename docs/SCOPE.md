# Release scope

This repository contains the source code needed to inspect, compile, and test
the VIBE periodic-reconstruction method. It is a focused extraction from the
TITAN research codebase, not a distribution of the complete TITAN program.

Included:

- the VIBE algorithms in `src/pbctopo/`;
- the frame-level TITAN adapter in `src/solv_pbc_pbctopo.cpp`;
- the molecular-topology graph builder;
- minimal host data declarations required to compile the adapter;
- independent algorithm-level tests.

The source-identical frame adapter includes optional H9 diagnostic exports
(disabled by default). These exports do not include the HoloLift solver.

Excluded:

- HoloLift and all trajectory-wide lifting code;
- Poisson-Boltzmann, boundary-element, and TABI solvers;
- unrelated TITAN analysis modules and the private TITAN executable;
- molecular trajectories and unpublished validation archives.

The static library preserves unresolved host calls in the adapter by design.
Applications embedding VIBE must supply coordinate storage, trajectory I/O,
logging, and output writing through their host integration.
