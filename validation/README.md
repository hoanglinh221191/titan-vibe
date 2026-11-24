Records of checks run on this package's titan-vibe.exe and on a Linux build of source/.

- QUICKSTART_RESULT.txt runs the examples of examples/quickstart against their EXPECTED.json (frame-report
  labels and molecule centroids).
- HELDOUT480_REPLAY.json replays the 480 held-out inputs against the reference outputs
  (examples/heldout480/manifest.json: PDB SHA-256 and report fields); these references are the outputs of
  release 1.1.0, and the 1.0.3 references are kept in manifest_v1.0.3.json.
- FIRST_FRAME5_REPLAY.json and REAL367_REPLAY.json replay the first frame and all 367 comparison frames of
  the five molecular systems, including the two Figure 1 frames: atom records are compared with the outputs
  deposited on 11 September 2026 at PDB precision, and the scientific report fields with
  REAL367_FIELDS_v1.1.0.csv; fields_changed_from_1_0_3 counts the report fields that differ from the
  deposited 1.0.3 reports.
- LINUX_BUILD.txt describes the Linux executable (GCC, static); LINUX_QUICKSTART_RESULT.txt,
  LINUX_HELDOUT480_REPLAY.json (with --coordinate-tolerance 0.001) and LINUX_REAL367_REPLAY.json are the
  same checks for it.

These records establish implementation agreement at the stated precision, not independent physical truth.
