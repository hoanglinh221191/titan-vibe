# Reproducibility record

## Source provenance

Version 1.0.1 refreshes the focused VIBE extraction on 9 September 2026 from
clean TITAN commit `6a8385cf1c63c42046b8a38053f003311e85e633`.
All VIBE submodules, the frame adapter and topology implementation are copied
from that commit. Five minimal host declaration headers are retained from
the public extraction and adjusted only when required for compilation;
they are not byte-identical copies of the complete TITAN host headers.
The file-level source mapping is recorded in `SOURCE_PROVENANCE.json`.

The separately rebuilt complete TITAN executable retains embedded identity
`0.2.0-dev.242+g6a8385cf1c63` (clean development source). The public source
package version `1.0.1` does not change that executable identity and this
repository still does not build a standalone TITAN executable.

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
