# 04 A peptide condensate (Martini 3), six snapshots

The system is 128 amylin peptides: 83 beads each, 10,624 beads in all, in
Martini 3, with water and ions removed. The cell is 202 × 202 × 607 Å. The
six frames are 60 ns apart (22.55-22.85 µs), late in a phase-separation run.
About 100 peptides form one droplet; the others are free or sit in clusters
of two or three.

| file | content |
|---|---|
| `reference.pdb` | the peptides, with CONECT records |
| `frames.xtc` | the six frames, copied from the trajectory without re-encoding |
| `index.ndx` | one group, every bead |
| `topol.top`, `molecule_0.itp` | the peptide topology |
| `martini_types.itp` | `[ defaults ]` and the `[ atomtypes ]` / `[ nonbond_params ]` rows of the 19 bead types the peptide uses, cut from `martini_v3.0.0.itp`; the full file gives the same output |
| `arg_local.in` | `-pbctopo-mode=llps-local`: every peptide made whole where it is, without placing peptides relative to each other (`local.pdb`) |
| `arg_global.in` | global VIBE on the same frames (`global.pdb`) |

Run both from this folder. Each run overwrites the frame report, so copy
`pbctopo_frame_report.csv` aside between the two runs:

```powershell
..\..\..\titan.exe -a arg_local.in
..\..\..\titan.exe -a arg_global.in
```

Contacts are bead pairs of two peptides within 8 Å (the Martini default).

| frame | time (µs) | `contact_graph_components` | `contact_cycle_rank` | global `evidence_state` |
|---|---|---|---|---|
| 0 | 22.546 | 26 | 0 | `weak_ambiguous` |
| 1 | 22.606 | 22 | 0 | `weak_ambiguous` |
| 2 | 22.666 | 23 | 1, along a | `extended_periodic` |
| 3 | 22.726 | 26 | 0 | `weak_ambiguous` |
| 4 | 22.786 | 24 | 0 | `weak_ambiguous` |
| 5 | 22.846 | 21 | 0 | `weak_ambiguous` |

**Local mode.** `status` is `local_unwrap` and `evidence_state` is
`not_applicable` in every frame. The contact network of the whole peptides is
still analysed.

**Frame 2.** The contacts of the droplet close through the a boundary
(`contact_cycle_b0` = (1, 0, 0)). The droplet touches its own image across
the 202 Å cell, so no finite placement keeps all of its contacts.

**Global mode, other frames.** These frames are `weak_ambiguous`
(`fallback_reason` `evidence_underdetermined`). With 21-26 contact
components, the frame does not fix where most free peptides sit relative to
the droplet.

For a condensate, use the local mode. The global run shows why.
