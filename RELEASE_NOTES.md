# TITAN-VIBE 1.2.0

TITAN-VIBE is now built from a VIBE-only TITAN source tree (`titan-vibe`), on Windows and on Linux; reconstruction is faster, with the same outputs as 1.1.0. See `RELEASE_NOTES_v1.2.0.md`.

## Version 1.1.0

Steric overlaps are reported instead of choosing a placement; a new evidence state, extended_periodic, reports contact networks that close through the periodic boundary; the common-topology reader follows [ molecules ]; quick-start examples with several snapshots. The focused library is no longer shipped. See `RELEASE_NOTES_v1.1.0.md`.

## Version 1.0.3

Faster reconstruction; candidates are checked in the lattice gauge that placed them; status, warm start and rescue paths follow the contact assessment of the selected placement. See `RELEASE_NOTES_v1.0.3.md`.

## Version 1.0.2

Executable distribution with buildable host source; see `RELEASE_NOTES_v1.0.2.md`.

## Version 1.0.1

Refresh VIBE implementation files from clean TITAN commit
`6a8385cf1c63c42046b8a38053f003311e85e633`.
The public source package and the separately rebuilt full TITAN executable
have distinct version identities; see `docs/REPRODUCIBILITY.md`.
The existing `titan-v1.0.0` tag is preserved.

## Version 1.0.0

Initial public reference release of VIBE.

This is a VIBE-only source release. Other TITAN modules and the private TITAN
executable are outside its scope.

- Framewise reconstruction in orthorhombic and triclinic cells.
- Topology-aware and metadata-based molecular connectivity.
- Certified nearest-copy geometry with fail-closed unresolved states.
- Virtual-box, component graph, exhaustive small-system, and restricted
  large-system search paths.
- Post-selection continuity and steric checks.
- Explicit identifiability, ambiguity, and finite-domain completeness reports.
- Local molecule/chain unwrapping mode for phase-separated systems.

This release does not claim trajectory-wide temporal lifting, unlimited search
completeness, thermodynamic stability, or a runtime advantage over other tools.
