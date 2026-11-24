# Release scope

Version 1.2.0 ships three things:
- `titan-vibe.exe` at the repository root, for Windows x86-64;
- the source of TITAN-VIBE under `source/`: the files that the VIBE-only
  build of TITAN reads, with a Makefile for Windows and Linux; VIBE is
  implemented in `source/src/pbctopo` and its host files;
- small examples under `examples/quickstart`, and the 480 held-out inputs with
  their reference outputs under `examples/heldout480`.

The source supplies the structure, index, topology and trajectory readers, the
force-field library, VIBE, the TITAN engine without its energy terms, and the
PDB writer. TITAN's energy, solvation, kinetics, analysis and HoloLift code is
not included; `source/src/vibe_only/` supplies the entry points the engine
names for those modules, and `titan-vibe` refuses their options.

The eight-contact cap is applied before the retained-contact translation
hypotheses are built. Contact support and ambiguity are therefore conditional
on that retained set, not on every contact in the input. The contact-cycle
lattice and the input-overlap counts use every pair within their cutoffs.

The large molecular frame subsets are in the separate data deposit. The
quick-start snapshots are selections, with water and ions removed, from
simulations of the methods article; their topologies are reduced to the
sections VIBE reads.
