# TITAN-VIBE 1.2.0

TITAN-VIBE is now the VIBE-only build of TITAN, built on Windows and on
Linux, and reconstruction is faster. On every check shipped with the package
the outputs equal those of 1.1.0.

## A VIBE-only build

- `source/` holds only the files that the VIBE-only build of TITAN reads
  (`VIBE_ONLY=1`, the default in this tree): 570 files, 154 of them C++
  source files, about 8 MB. Release 1.1.0 shipped the full TITAN host
  (47 MB). The tree contains the structure, index, topology and trajectory
  readers, VIBE, the engine without its energy terms, and the PDB writer.
- The executable is `titan-vibe` (`titan-vibe.exe`, 7.8 MB; the 1.1.0
  `titan.exe` was 22.6 MB).
- `titan-vibe` accepts the VIBE options and the input, output, frame and
  thread options; any other option is refused with a message before it is
  parsed. `-h` prints the VIBE reference, as `-vibe-help` does. `-nocal=true`
  is required, and interactive mode is refused.
- The engine is unchanged. The entry points of the modules left out are
  supplied by `source/src/vibe_only/` (see its README): those the engine calls
  on every run do nothing, and the others stop with an internal error naming
  the function.

## Linux

- The source builds on Linux with GCC 13 (tested: GCC 13.3.0,
  Ubuntu 24.04.1 LTS under WSL2: the static build with every check below, and the default
  build, which loads only the C library, on the quick-start examples). The
  README gives the commands.
- The Linux build reproduces the quick-start results and the 367 frames of the
  five molecular systems at PDB precision. Of the 480 held-out outputs,
  21 differ from the references in some coordinates by 0.001 Å,
  one unit of the last printed digit; labels and report fields are the same.
- `replay.py --coordinate-tolerance 0.001` checks such a build. The reference
  PDBs it needs are now in `examples/heldout480/<case>/`.

## Faster reconstruction, same output

The VIBE solver keeps its results bit for bit while doing less work:
- when the components are placed rigidly, the minimum spanning tree of the
  anchors is computed exactly on a sparse edge set (otherwise the dense scan
  runs);
- beam continuity terms are computed once per depth and lattice shift;
- trials of the refinement's first stage are scored in parallel and committed
  in serial order;
- the beam states of a depth are expanded in parallel;
- for 256 anchors or more, the dense anchor MST fuses relaxation and selection
  into one AVX2 pass;
- the final beam rescoring is speculative and parallel, committed in order;
- certification skips atoms that cannot meet another component, and the
  nearest images of the edges are computed once per frame.

Wall time for 21 frames (median of 5 rounds after one warm-up, Intel Core
Ultra 7 255H, Windows 11; both releases in the same rounds, in rotated order;
1.2.0 is a build of this source):

| System | 1.1.0, 8 threads | 1.2.0, 8 threads | 1.1.0, 1 thread | 1.2.0, 1 thread |
|---|---|---|---|---|
| protein-glycan complex | 14.0 s | 8.3 s | 34.0 s | 12.8 s |
| peptide pair | 1.65 s | 1.21 s | 2.99 s | 1.43 s |
| hexamer | 12.0 s | 10.1 s | 12.1 s | 10.0 s |
| tetramer | 6.9 s | 6.5 s | 6.9 s | 6.6 s |
| condensate | 3.3 s | 3.3 s | 3.6 s | 3.3 s |

The tetramer, hexamer and condensate times are dominated by reading the
solvated frames. Peak memory is unchanged.

## Replay drivers

- Every driver runs `titan-vibe.exe` on Windows and `source/titan-vibe` on
  Linux by default, and takes `--binary`.
- The held-out references and `validation/REAL367_FIELDS_v1.1.0.csv` are
  unchanged from 1.1.0.
