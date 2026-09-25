# TITAN-VIBE 1.0.3

Reconstruction is faster. Candidates are checked in the lattice gauge that placed
them, and every code path now takes the frame status from the contact assessment
of the selected placement.

- **Speed.**
  - A nearest-copy query now takes an exact shortcut when the rounded
    reduced-basis image is provably the unique minimum.
  - The steric scan walks a sorted cell array.
  - Candidates in a complete uncertainty band are checked in parallel on the TITAN
    executor.
  - None of these steps changes the selected placements or the coordinates.
- **Consistency test in the placing gauge.** Offsets recovered from placed
  coordinates are shifted by one common lattice vector before they are checked.
  A candidate is therefore no longer rejected for a whole-assembly translation.
- **Status and warm start follow the selected evidence.** The warm-state commit
  and the rescue paths now use the same assessment that 1.0.2 already applied when
  writing the frame report.
- **Evidence export.** The new option `-pbctopo-evidence-export` writes the
  per-candidate contact relations and a connectivity margin next to the frame
  report.
- **Masses.** Atomic masses are read from the topology.
- **Triclinic cells.** The host builds the lattice vectors of general triclinic
  cells exactly from the cell lengths and angles. This can change the last
  binary digits of the cell matrix and therefore the input checksum
  (`wrapped_coordinate_hash`) of a frame.
- **Rejected candidates.** Rejected evaluations of an assignment go to
  `pbctopo_rejected_candidates.csv`, not to the certified band.
- **Options.** Every command-line option of 1.0.2 is kept, including
  `-pbctopo-component-shell-radius` and `-pbctopo-h9-orbit`.

**Reproduction against 1.0.2.**
- *480 held-out inputs:* the placements, statuses, evidence states and consistency
  results are the same as in the archived analysis. 80 reports differ only in the
  name of the search route. One coordinate differs by 0.001 Å from PDB rounding.
- *Molecular frames:* the five representative frames and the 367 comparison frames
  give byte-identical coordinate files. In 2 of the 31 peptide-pair frames the
  input checksum changes because of the exact triclinic cell matrix.
- *Designed frames:* the 3,600 generated and 4,320 noisy frames give the same
  statuses and evidence states as 1.0.2.

As in 1.0.2, the host source omits HoloLift. The data deposit is
https://doi.org/10.5281/zenodo.22710653.
