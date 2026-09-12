# Reproducibility record for 1.0.2

Version 1.0.2 adds the executable distribution and complete buildable host with
HoloLift omitted. See the root README for exact commands. The focused VIBE source
remains identical to 1.0.1; the host snapshot is based on TITAN commit
`6a8385cf1c63c42046b8a38053f003311e85e633`. Source normalization, adapted host
headers, executable hash, and imported DLLs are recorded in SOURCE_AUDIT.json.

The embedded host version is `0.2.0-dev.242+g6a8385cf1c63.dirty`. The suffix records
the HoloLift-excluded build and pinned metadata. TITAN-VIBE distribution version
1.0.2 does not relabel the underlying research-source identity.

Validation already completed for this exact executable: 480/480 held-out fixed
input PDB hashes and 18 scientific report fields match archived outputs; the
first frame of each of the five representative systems matches coordinates at
PDB precision and all 18 report fields. Two additional Figure 1 frames match.
A longer replay was stopped at the user's request to avoid rerunning every test;
its partial results are retained and are not described as full completion.

No VIBE algorithm was changed, no new MD simulation or timing campaign was run,
and the full library test suite was not rerun just to change the distribution
version. Existing tests remain available. The release commit intentionally skips
automatic CI for this packaging-only change; validation records are local checks,
not a claim of a fresh remote CI run.

The original designed-case oracle and the original MDVWhole/GROMACS and timing
protocols remain in data DOI 10.5281/zenodo.22710653. A replay establishes output
agreement, not independent physical truth or a new held-out study.
