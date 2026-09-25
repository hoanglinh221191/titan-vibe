"""Match the engine's Eigen ABI when compiling a linked scientific driver.

Do not link an AVX2-only driver to native AVX-512 objects: Eigen's dynamic
allocator and inline packet loads then disagree about 32/64-byte alignment.
Native objects and their drivers must be built on the same host. Unknown or
missing profiles are rejected instead of guessing an ABI.
"""
import json
from pathlib import Path

def engine_arch_flags(objects_root):
    manifest = Path(objects_root).resolve().parent / 'generated/titan_build_manifest.json'
    profile = json.loads(manifest.read_text(encoding='utf-8'))['profile']
    if profile == 'native':
        return ['-march=native']
    if profile == 'portable':
        return ['-mavx2', '-mfma', '-mtune=generic']
    raise ValueError(f'Unsupported engine profile {profile!r} in {manifest}')
