"""Oracles that call engine functions on synthetic inputs, linked against
the objects of a finished build:

  kin16  DHAM frame weights after window matching (run_dham_postproc)
  kin17  MBAR demux map resolved by time
         (TitanMBAR::run_temperature_ladder_reweighting)
  kin18  MBAR NPT: reduced potential beta (U + p V), analytic Gamma ensemble
  kin19  DHAM trajectory breaks with a repeated t = 0 frame
  kin20  MBAR block bootstrap of AR(1) energies against the analytic SE
  kin21  energy grid coarser than the trajectory grid (one frame per row)
  kin14  transfer entropy with Miller-Madow (chi^2 bound, spurious edges)
  kin06  MSM implied timescales: NaN with rank kept, multiplicity of 1
  kin22  -entropy (Duan-Zhang) on Gaussian energies, untrimmed, O(N) series
  kin11  2-D FES: empty bins NaN/"nan", basins next to empty bins (PCA/TICA)
  kin23  -dham/-mbar temperature files and demux maps: comments skipped,
         other text refused with file and line
  gbsa   fused GBSA kernel (cal_energy_gbsa_avx2): intramolecular Delta-GB
         (ENE-15) and cross G_polar (ENE-14) against the scalar long-double
         sums, and the engine's tiling against one call
  gbsa_timing  (not run by default) in-process timing of the kernel parts
  smoke_qnative, smoke_s2, smoke_msd, smoke_membrane, smoke_contact_map,
  smoke_dihedral  the module unit smokes (test/titan_*_smoke.cpp), built with
         -DTITAN_SMOKE_LINK_ENGINE so their stubs give way to the engine's
  path01 PATH-01: -de-pathway-radius keeps a recovery band of non-zero width
         and explicit -pathway-recovery-* options win in any order

  python test/engine_linked_oracles.py --work-root test/tmp/engine_linked
  python test/engine_linked_oracles.py --source-root C:/titan_kinbase --work-root test/tmp/engine_linked_base

--source-root selects the tree whose headers and build/<profile>/src objects
are used (a baseline export shows the old failures). Exit 1 if any oracle fails.
"""
import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from science_build_profile import engine_arch_flags  # noqa: E402

ROOT = HERE.parent
DRIVERS = {
    'kin16': 'titan_kin16_dham_weights_oracle.cpp',
    'kin17': 'titan_kin17_mbar_demux_oracle.cpp',
    'kin18': 'titan_kin18_mbar_npt_pv_oracle.cpp',
    'kin19': 'titan_kin19_dham_break_dt_oracle.cpp',
    'kin20': 'titan_kin20_mbar_block_bootstrap_oracle.cpp',
    'kin21': 'titan_kin21_energy_grid_oracle.cpp',
    'kin14': 'titan_kin14_transfer_entropy_oracle.cpp',
    'kin06': 'titan_kin06_msm_its_oracle.cpp',
    'kin22': 'titan_kin22_interaction_entropy_oracle.cpp',
    'kin11': 'titan_kin11_fes_empty_bins_oracle.cpp',
    'kin23': 'titan_kin23_input_parsers_oracle.cpp',
    'gbsa': 'titan_gbsa_kernel_oracle.cpp',
    # timing probe, not an oracle: prints min/median ms of the kernel parts
    'gbsa_timing': 'titan_gbsa_kernel_timing.cpp',
    'smoke_qnative': 'titan_qnative_smoke.cpp',
    'smoke_s2': 'titan_s2_smoke.cpp',
    'smoke_msd': 'titan_msd_smoke.cpp',
    'smoke_membrane': 'titan_membrane_smoke.cpp',
    'smoke_contact_map': 'titan_contact_map_smoke.cpp',
    'smoke_dihedral': 'titan_dihedral_smoke.cpp',
    'path01': 'titan_path01_alias_oracle.cpp',
}
DEFAULT_CASES = ('kin16,kin17,kin18,kin19,kin20,kin21,kin22,kin14,kin06,kin11,kin23,gbsa,'
                 'smoke_qnative,smoke_s2,smoke_msd,smoke_membrane,'
                 'smoke_contact_map,smoke_dihedral,path01')


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--compiler', default='clang++' if os.name == 'nt' else 'g++')
    p.add_argument('--source-root', type=Path, default=ROOT)
    p.add_argument('--objects-root', type=Path, default=None,
                   help='default: <source-root>/build/native/src')
    p.add_argument('--work-root', type=Path, default=ROOT / 'test/tmp/engine_linked')
    p.add_argument('--cases', default=DEFAULT_CASES)
    a = p.parse_args()
    src_root = a.source_root.resolve()
    objects_root = (a.objects_root or src_root / 'build/native/src').resolve()
    work = a.work_root.resolve()
    work.mkdir(parents=True, exist_ok=True)
    objects = sorted(x for x in objects_root.rglob('*.o') if x.name.lower() != 'main.o')
    if not any(x.name == 'titan_mbar.o' for x in objects):
        raise SystemExit(f'no engine objects under {objects_root}: build TITAN first')
    rsp = work / 'objects.rsp'
    rsp.write_text('\n'.join('"' + x.as_posix() + '"' for x in objects), encoding='utf-8')
    suffix = '.exe' if os.name == 'nt' else ''
    # the drivers need the compiler's runtime DLLs (MSYS2 UCRT64 on Windows)
    env = os.environ.copy()
    compiler_path = shutil.which(a.compiler)
    if compiler_path:
        env['PATH'] = str(Path(compiler_path).parent) + os.pathsep + env.get('PATH', '')
    status = {}
    for case in [c.strip() for c in a.cases.split(',') if c.strip()]:
        exe = work / (case + suffix)
        cmd = [a.compiler, '-std=c++23', '-O2', *engine_arch_flags(objects_root),
               '-I' + str(src_root / 'src'), '-I' + str(src_root / 'src/eigen'),
               '-I' + str(src_root / 'src/gmxtraj/include'),
               '-I' + str(objects_root.parent / 'generated'),
               '-DTITAN_SMOKE_LINK_ENGINE',
               '-pthread', '-fopenmp', str(HERE / DRIVERS[case]), '@' + str(rsp),
               '-o', str(exe)]
        if os.name == 'nt':
            cmd.insert(1, '-fuse-ld=lld')
        build = subprocess.run(cmd, capture_output=True, text=True)
        (work / (case + '_build.log')).write_text(build.stdout + build.stderr, encoding='utf-8')
        if build.returncode:
            print(build.stdout + build.stderr)
            status[case] = 'BUILD FAILED'
            continue
        case_dir = work / case
        case_dir.mkdir(exist_ok=True)
        run = subprocess.run([str(exe)], cwd=case_dir, capture_output=True, text=True,
                             encoding='utf-8', errors='replace', env=env)
        (work / (case + '_run.log')).write_text(run.stdout + run.stderr, encoding='utf-8')
        lines = [line for line in run.stdout.splitlines()
                 if line.startswith(('PASS', 'FAIL', 'case ', '  '))
                 and not any(w in line for w in ('K-Means', 'Converged in'))]
        print('\n'.join(lines))
        if run.returncode and not lines:
            print(f'{case}: exit code {run.returncode} without output')
        status[case] = 'PASS' if run.returncode == 0 else 'FAIL'
    print(' '.join(f'{k}={v}' for k, v in status.items()))
    return 0 if status and all(v == 'PASS' for v in status.values()) else 1


if __name__ == '__main__':
    sys.exit(main())
