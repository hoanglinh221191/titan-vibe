# TITAN-VIBE 1.0.2

The previous release built only a static VIBE library, so readers could not run
the deposited frame-reconstruction workflow directly. Version 1.0.2 supplies a
Windows x86-64 `titan.exe`, full buildable TITAN host source with HoloLift omitted,
and standard-library Python drivers for the 480 compact held-out inputs and
five deposited molecular frame sets. The existing focused library is retained.

The VIBE algorithm is unchanged from 1.0.1. Source audit confirms the eight-contact
cap precedes translation-hypothesis grouping; ambiguity is conditional on retained
contacts. The executable's embedded TITAN snapshot identity is documented separately
from this distribution version.

Validation: 480 designed-input reconstructions match archived PDB hashes and 18
scientific report fields; five representative first frames and two additional
figure frames match. Package checksums, no-HoloLift symbol checks and Windows DLL
imports are recorded. No new MD or complete-suite rerun was performed for release.

The earlier 1.0.1 release and DOI remain available. The data deposit is
https://doi.org/10.5281/zenodo.22710653.
