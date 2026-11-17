# 03 Two glycosylated chains, five snapshots

Chains A and B of PDB entry 7WBL with their N-glycans: 14,393 atoms from an
AMBER ff19SB simulation, with water and ions removed. The five frames are
10 ns apart (0-40 ns) in a cubic cell of 133.2 Å.

| file | content |
|---|---|
| `reference.pdb` | frame 0 of the selection, with CONECT records |
| `frames.xtc` | the five frames, copied from the trajectory without re-encoding |
| `index.ndx` | one group, every atom |
| `topol.top`, `toppar/` | the force-field atom types and the `[ atoms ]` and `[ bonds ]` of both chains (reduced for the example; not a simulation topology) |
| `arg.in` | the VIBE run |

Run it from this folder:

```powershell
..\..\..\titan.exe -a arg.in
```

`placed.pdb` has five models. In every frame of `pbctopo_frame_report.csv`:
- `status` is `certified`;
- `evidence_state` is `supported_connected`;
- `component_count` is 2;
- `contact_cycle_rank` is 0.

The glycans are bonded to their chains in the topology, so each chain with its
glycans is one part. Both parts are written whole and in contact. Each frame
takes under a second.
