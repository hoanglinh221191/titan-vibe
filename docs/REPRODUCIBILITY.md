# Reproducibility record

## Source provenance

This repository is a focused extraction of the VIBE implementation used for
the associated methods study. The extraction was prepared on 20 August 2026
from the TITAN development tree and contains the frame solver, all VIBE
submodules, the topology graph builder, and only the host declarations required
to compile them.

The release deliberately omits unrelated TITAN implementations. In particular,
it contains no HoloLift, Poisson-Boltzmann, boundary-element, TABI, PBCN/FixBox,
or general trajectory-analysis source. The public Git repository starts with a
new root commit so omitted code is not recoverable from its history.

## Release checks

Before creating a release tag:

1. configure and build the C++23 AVX2 static library with eight parallel jobs;
2. run every test listed by `test/run_core_smokes.ps1`;
3. scan every tracked path and file for excluded module implementations;
4. record the public commit, annotated release tag, and source-archive SHA256;
5. keep trajectory inputs and manuscript evidence in a separate immutable data
   archive with its own manifest and checksums.

The tests in this repository establish implementation invariants and regression
behavior. They do not replace the designed-case and trajectory-level
experiments reported in the methods article.

## Scientific interpretation

- `certified` means that the stated finite search and post-selection checks
  completed successfully.
- `bounded_exhaustive` does not mean exhaustive over every lattice translation.
- `weak_observation` means the wrapped frame lacks enough information for a
  uniquely supported relative placement.
- `llps-local` reconnects atoms within each molecule or chain and intentionally
  preserves the intermolecular spatial distribution.

## Host integration

The static library compiles the same VIBE frame solver used by TITAN. A host
program must provide trajectory reading, coordinate storage, output writing,
and the small logging/executor interfaces referenced by the adapter. The
algorithm-level tests are self-contained and do not require the private TITAN
executable.

## External data

The publication data archive should contain the exact argument files, input
manifests, checksums, source tables, figure-generation scripts, and validation
summaries. Its DOI or accession will be added here and to the manuscript after
deposition.
