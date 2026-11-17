# TITAN-VIBE 1.1.0

Steric overlaps no longer choose a placement, a new evidence state reports
contact networks that close through the periodic boundary, and a short VIBE
reference is printed by `titan.exe -vibe-help`.

## Steric overlaps are reported, not used to choose

Under periodic boundary conditions, an atom pair of two parts that overlaps at
the relative image a placement uses is also an overlapping pair of the
simulated system, whatever placement is written. A steric test on placed
coordinates can therefore only remove a placement that the contacts support.
In 1.0.x it did exactly that on noisy designed frames.

- `-pbctopo-steric=report` (default): overlaps do not reject candidates and do
  not enter the pair-shift score. The frame report counts the overlapping atom
  pairs between parts in the input (`input_overlap_pairs`,
  `input_overlap_unit_pairs`, `input_overlap_min_distance_A`,
  `input_overlap_min_ratio`), with the 1.0.x thresholds
  (lambda (r_i + r_j), lambda = `-pbctopo-vdw-clash-scale`, default 0.55;
  `-pbctopo-clash-cutoff` for atoms without radii).
- `-pbctopo-steric=filter` restores the 1.0.x candidate rejection.
- In metadata mode, intra-component overlaps that are already in the wrapped
  input no longer enter `hard_min_steric_margin_A`.
- In `llps-local` mode `nonclashing_candidate_count` is `na`.

## Contact networks that close through the boundary

The observed contacts define a graph on the parts: an edge (u, v, s) for every
heavy-atom pair of parts u and v within the contact distance, with s the
relative lattice image that pair requires, and a loop (u, u, s), s != 0, when
a part touches its own image. A finite placement that keeps every contact
exists exactly when the image vectors sum to zero around every cycle. The
cycle sums span a sublattice of Z^3, and its rank is now reported.

- New evidence state `extended_periodic`, with `fallback_reason`
  `evidence_periodic`: the cycle lattice has rank 1 (a filament), 2 (a sheet)
  or 3. This includes a pair in contact at two images (the 1.0.x
  `weak_ambiguous` case), a part touching its own image (labelled
  `supported_connected` in 1.0.x), and longer rings. `weak_ambiguous` now
  means a rank-0 network whose selected placement leaves a relation
  unsupported.
- New columns: `contact_cycle_rank`, the generators `contact_cycle_b{0,1,2}_*`,
  `contact_graph_components`, `contact_periodic_graph_components`,
  `contact_max_component_rank`, `self_image_units`, `self_image_atom_pairs`,
  `self_image_min_distance_A`. They are also written in `llps-local` mode,
  whose units are the whole molecules (a bilayer has rank 2).
- A bonded cycle that runs through the boundary (no whole lift exists) is
  reported as `fallback_reason` `bonded_winding`, with `bonded_winding_rank`
  and `topological_winding_rank` set.

## Topology input

- The `-common_topo` reader builds the system in `[ molecules ]` order from
  every `[ moleculetype ]` of every included file. Before, it read one
  molecule type per file, in include order, and only `[ atoms ]` rows with a
  mass column. Topologies with several molecule types in one file, an include
  order different from `[ molecules ]`, a molecule name listed twice, or
  7-column `[ atoms ]` rows (mass from `[ atomtypes ]`, as in GROMACS) are
  now read in the right atom order. Topologies of one molecule type per file,
  included in `[ molecules ]` order, give the same result as before.
- A `[ atomtypes ]` section starts only at its header; a comment that
  mentions the word no longer starts it.

## Interface

- `titan.exe -vibe-help` prints the options a reconstruction needs and the
  meaning of every evidence state.
- The frame-report schema is version 8 (columns appended to version 7);
  readers of this release accept versions 7 and 8.

## Examples

`examples/quickstart` holds small inputs with several snapshots each:
- `01_peptide_pair`: 31 frames, two peptides, triclinic cell.
- `02_labels`: six designed objects, one for each label.
- `03_glycoprotein_pair`: five frames, two glycosylated chains.
- `04_condensate_martini`: six frames, 128 Martini 3 peptides, local and
  global mode.

`run_quickstart.py` runs them and compares the frame-report labels and the
molecule centroids of the output with `EXPECTED.json`.

## Reproduction against 1.0.3

**Held-out inputs (480).** 416 cases give the same coordinates and report
fields as the 1.0.3 references.
- The 60 cases of the `invalid_stress` stratum are designs with overlapping
  atoms. In 1.0.3 the steric test removed the designed placement and they were
  labelled `weak_ambiguous`. In 1.1.0 all 60 are placed as designed and
  labelled `certified` / `supported_connected`, and all 60 report the overlap
  in `input_overlap_pairs`. No case of the other strata reports an overlap.
- In 20 of these 60 cases, `feasible_assignments_within_domain` changes from
  124 to 125.
- In 4 other cases, one or two atom coordinates differ in the last printed
  digit (0.001 Å). The placement is the same.
- The 1.0.3 references are kept in `examples/heldout480/manifest_v1.0.3.json`.

**Molecular systems (367 comparison frames and the two Figure 1 frames).** The
reconstructed atom records are identical to the deposited 1.0.3 outputs at PDB
precision. The report labels change in four frames, all to
`extended_periodic`:
- tetramer frame 20 (`supported_connected` in 1.0.3): one chain touches its
  own periodic image at 4.09 Å;
- peptide-pair frames 0, 4 and 28 (`weak_ambiguous` in 1.0.3): `contact_cycle_rank`
  is 1, and in frame 28 a peptide also touches its own image at 3.47 Å.

Every other report field compared by `replay_real.py` is unchanged. The fields
of this release are in `validation/REAL367_FIELDS_v1.1.0.csv`.
