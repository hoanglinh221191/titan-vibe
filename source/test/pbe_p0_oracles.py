"""Compile production PB kernels and run independent oracles on Linux/Windows.

screen and k4 are the P0 oracles (with --verify-regression, the original
formula must be rejected). far is the treecode far-field oracle: the radial
Taylor kernels against a long-double reference, the fast and four-lane far
kernels against both, timing skipped. ses is the SES self-intersection
regression: true overlaps still found, shared-vertex-only contacts no longer
counted. near is the curved-cell near quadrature across narrow gaps against
an order-256 reference. p1 is the degree of the p1 far field: the error of
every kernel K1-K4 falls as h^2 when a source cell shrinks by h. parity is
the FMM near-field branches on a small PN sphere under a leaf-1 tree: p1
packed and scalar, p4 packed and scalar agree to rounding, converge with
theta in every block, and reproduce the direct operator when every pair is
direct."""
import argparse
import hashlib
import json
import os
from pathlib import Path
from science_build_profile import engine_arch_flags
import subprocess

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--compiler', default='clang++' if os.name == 'nt' else 'g++')
parser.add_argument('--work-root', type=Path, default=ROOT / 'test/tmp/pbe_p0_oracles')
parser.add_argument('--case', choices=['screen', 'k4', 'far', 'ses', 'near',
                                       'p1', 'parity'])
parser.add_argument('--objects-root', type=Path, default=ROOT / 'build/native/src')
parser.add_argument('--verify-regression', action='store_true',
                    help='Also verify each oracle rejects the original P0 formula')
args = parser.parse_args()
args.work_root.mkdir(parents=True, exist_ok=True)

def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

common = [args.compiler, '-std=c++23', '-O2', '-fno-access-control',
          *engine_arch_flags(args.objects_root),
          '-ffunction-sections', '-fdata-sections', '-Isrc', '-Isrc/eigen',
          '-Isrc/gmxtraj/include', '-Isrc/tabi', '-pthread']
if os.name == 'nt':
    common += ['-fuse-ld=lld']
sources = {
    'screen': ['src/pbe/grid/linear_solver.cpp'],
    'k4': ['src/tabi/' + name + '.cpp' for name in (
        'tabi_near', 'tabi_treecode_series', 'tabi_treecode_internal',
        'tabi_treecode_leaf', 'tabi_fmm_m2l_cartesian', 'tabi_fmm_m2l_local',
        'tabi_fmm_m2l_prototype', 'tabi_preconditioner')],
    'far': ['src/tabi/tabi_treecode_series.cpp'],
    'ses': ['src/tabi/tabi_calderon_aux_surface.cpp'],
    'near': ['src/tabi/tabi_node_patch_curved.cpp'],
    'p1': [],
    'parity': [],
}
drivers = {'far': 'test/tabi_treecode_radial_oracle.cpp',
           'ses': 'test/bie_surface_intersection_regression.cpp',
           'near': 'test/bie_curved_near_quadrature.cpp',
           'p1': 'test/tabi_fmm_p1_degree_probe.cpp',
           'parity': 'test/tabi_fmm_branch_parity.cpp'}
# A driver that includes an engine source links without that source's object.
included = {'p1': ['tabi/tabi_fmm_m2l_prototype.o']}
run_args = {'far': ['--no-timing']}
for name, files in sources.items():
    if args.case and args.case != name:
        continue
    exe = args.work_root.resolve() / (name + ('.exe' if os.name == 'nt' else ''))
    driver = drivers.get(name, 'test/titan_pbe_' + name + '_oracle.cpp')
    cmd = common + [driver] + files
    replaced = {str(Path(f).with_suffix('.o').relative_to('src')).lower() for f in files}
    replaced |= {str(Path(o)).lower() for o in included.get(name, [])}
    objects = [p.resolve() for p in args.objects_root.rglob('*.o')
               if p.name.lower() != 'main.o' and
               str(p.relative_to(args.objects_root)).lower() not in replaced]
    if not objects:
        raise RuntimeError('Build TITAN first; no objects under ' + str(args.objects_root))
    rsp = args.work_root.resolve() / (name + '_objects.rsp')
    rsp.write_text('\n'.join('"' + p.as_posix() + '"' for p in objects))
    cmd += ['@' + str(rsp), '-fopenmp', '-o', str(exe)]
    subprocess.run(cmd, cwd=ROOT, check=True)
    result = subprocess.run([str(exe)] + run_args.get(name, []), cwd=args.work_root,
                            capture_output=True, text=True)
    (args.work_root / (name + '.log')).write_text(result.stdout + result.stderr)
    print(result.stdout + result.stderr, end='', flush=True)
    result.check_returncode()
    if args.verify_regression and name in ('screen', 'k4'):
        source = files[0]
        text = (ROOT / source).read_text()
        old, new = (
            ('return eps_out * kappa * kappa * shape.h * shape.h;',
             'return kappa * kappa * shape.h * shape.h;') if name == 'screen' else
            ('return -fsecond * nxu * nyu - (fprime * inv_r) * (nxny - nxu * nyu);',
             'return fsecond * nxu * nyu + (fprime * inv_r) * (nxny - nxu * nyu);'))
        if text.count(old) != 1:
            raise RuntimeError('Regression mutation anchor changed: ' + source)
        mutant = args.work_root.resolve() / (name + '_old_formula.cpp')
        mutant.write_text(text.replace(old, new))
        regression_exe = exe.with_name(exe.stem + '_regression' + exe.suffix)
        regression_cmd = [str(mutant) if c == source else
                          str(regression_exe) if c == str(exe) else c for c in cmd]
        subprocess.run(regression_cmd, cwd=ROOT, check=True)
        result = subprocess.run([str(regression_exe)], cwd=args.work_root,
                                capture_output=True, text=True)
        (args.work_root / (name + '_regression.log')).write_text(result.stdout + result.stderr)
        if result.returncode != 1 or 'FAIL ' not in result.stderr:
            raise RuntimeError('Oracle did not reject the old formula: ' + name)
        print('PASS regression rejected original ' + name + ' formula', flush=True)
    record = {
        'status': 'PASS', 'command': cmd,
        'compiler': subprocess.check_output([args.compiler, '--version'], text=True),
        'compiled_sources_sha256': {f: sha256(ROOT / f) for f in files + [driver]},
        'linked_objects_sha256': {p.as_posix(): sha256(p) for p in objects},
        'binary_sha256': sha256(exe),
        'original_formula_rejected': args.verify_regression and name in ('screen', 'k4'),
    }
    (args.work_root / (name + '_verification.json')).write_text(json.dumps(record, indent=2))
