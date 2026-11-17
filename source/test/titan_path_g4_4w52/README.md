# 4W52 TITAN-Path gate fixture

This fixture validates the external stock-GROMACS/TITAN-Path loop on the
T4-lysozyme L99A–benzene system.

## Route provenance

- `opes_route_a.tpath` is retained only as a rejected candidate. The public
  `eq_structure*.gro` files are independent bound-state starting structures,
  not ordered unbinding frames. `derive_opes_route.py` makes each molecule
  whole and proves that the apparent wrapped route is invalid. See
  `source_manifest.json`; never use this route for MD.
- `hj_literature_route_a.tpath` is the accepted pilot input. It is explicitly
  a literature-seeded geometric hypothesis for the H/J-cyan channel, not an
  observed OPES trajectory or CAVER centerline. Its deterministic construction
  is verified by `derive_hj_literature_route.py` and bound by
  `source_manifest_hj.json`.

## Controller rule that matters

At a direction change, the new GROMACS coordinate is projected from columns
3-5 of the final instantaneous `pullx` row (`pull-print-components=yes`). Do
not use TITAN/PDB COMs: their PBC image differs from the GROMACS pull-group
image on this system. The rejected `run_001` mechanical pilot documents the
resulting 0.0500824 nm extension discontinuity. Corrected `run_003` reduces the
observed residual to 3.0e-7 nm.

## Reproduction outline

1. Checkout the pinned public input repository recorded in the manifests.
2. Build the Linux `titan_wsl` binary from this source tree.
3. Create the initial TPR with `smd_bootstrap.mdp` for the 2-ps mechanical gate
   or `smd_bootstrap_20ps.mdp` for the production-interval gate.
4. Validate the selected controller JSON with `--preflight`.
5. Run the controller; `INCOMPLETE_MAX_SEGMENTS` is the expected two-segment
   pilot stop.
6. Run `pilot_report.py` on the work directory. It checks observed restart
   extension, algebraic invariants, direction cap, registration, clearance,
   output time grids, hashes, and native `gmx check` results.

The pilot evidence is a software/mechanical gate. It is not a publication
replica and cannot support pathway accuracy or affinity claims by itself.
