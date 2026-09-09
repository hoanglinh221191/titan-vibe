# TITAN-VIBE 1.0.1

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
