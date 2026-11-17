"""ENE-14 oracle: the Debye-screened Coulomb energy of the AVX2 kernels must
match a double-precision reference to the printing precision.

AMBER99SB-ILDN variant of the H10 dimer (test/fixtures/e2e_h10/
topology_variants), lib_ff, one frame, -elec-screening at several salt
concentrations, with and without -denfo (nofilter and intrin kernels).
Reference: sum over chain A x chain B of 138.935457839 q_i q_j exp(-kappa r)/r
with kappa = 3.291 sqrt(I) nm^-1 (titan_electrostatics.h) and the per-atom
charges of the run's TITAN_PARAMETER_AUDIT.

Before ENE-14 the kernels used f_exp_avx (relative error up to 3.2e-6 per
exponential), which shows up as a deviation of the summed energy that is
several times the printing resolution of cal.xvg (7 significant digits).

  python test/screening_exp_oracle.py --binary titan.exe --work-root test/tmp/screening_exp
"""
import argparse
import math
import os
import shutil
import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))

from topology_oracles import lj_reference as ljr  # noqa: E402
from topology_oracles.gromacs_reference import index_groups, pdb_coordinates, variant_dir  # noqa: E402
from topology_oracles.titan_run import load_audit, run, xvg_rows  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
FF = 'amber99sb-ildn'
SALTS = [0.15, 1.0, 4.0]
CUTOFF_A = 12.0
TOL_REL = 6e-7   # half a unit in the 7th significant digit, plus summation order


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--binary', default='titan.exe' if os.name == 'nt' else 'titan')
    ap.add_argument('--work-root', default=str(ROOT / 'test' / 'tmp' / 'screening_exp'))
    a = ap.parse_args()
    binary = Path(a.binary).resolve() if Path(a.binary).exists() else Path(a.binary)
    work = Path(a.work_root).resolve()
    if work.exists():
        shutil.rmtree(work)
    vdir = variant_dir(FF)
    libff = (vdir / 'lib_ff_subset').as_posix() + '/'
    groups = list(index_groups(vdir / 'index.ndx').values())
    xyz = pdb_coordinates(vdir / 'conf.pdb')
    dist = ljr.distances([xyz[i - 1] for i in groups[1]], [xyz[i - 1] for i in groups[2]])
    failures, charges = [], None
    for salt in SALTS:
        kappa = 3.291 * math.sqrt(salt)
        for cutoff in (None, CUTOFF_A):
            extra = ['-elec-screening', f'-saltcon={salt}'] + ([f'-denfo={cutoff:g}'] if cutoff else [])
            name = f'I={salt} M' + (f', cutoff {cutoff:g} A' if cutoff else '')
            run_dir = work / f'salt{salt}_{"cut" if cutoff else "nocut"}'
            args = [f'-ff={FF}', f'-f={(vdir / "conf.pdb").as_posix()}',
                    f'-n={(vdir / "index.ndx").as_posix()}', '-g=1-2', '-eu=kJ/mol',
                    '-noprint-side_chain_contact=true', '-nt=1', f'-path-a={libff}',
                    f'-path-n={libff}'] + extra
            ok, _ = run(binary, run_dir, args)
            if not ok:
                failures.append(f'{name}: run did not finish')
                continue
            if charges is None:
                audit = load_audit(run_dir / 'params.tsv')
                charges = ([float(r['charge_e']) for r in audit if r['molecule'] == 'mol1'],
                           [float(r['charge_e']) for r in audit if r['molecule'] == 'mol2'])
            qa, qb = charges
            ref, k = 0.0, 0
            for q1 in qa:
                for q2 in qb:
                    r = dist[k]
                    k += 1
                    if cutoff is None or r <= cutoff * 0.1:
                        ref += 138.935457839 * q1 * q2 * math.exp(-kappa * r) / r
            got = xvg_rows(run_dir / 'cal.xvg')[0][1]
            rel = abs(got - ref) / abs(ref)
            status = 'ok' if rel <= TOL_REL else 'FAIL'
            if rel > TOL_REL:
                failures.append(f'{name}: elec {got:.7g} vs reference {ref:.9g} (rel {rel:.2e})')
            print(f'  {name:26s} elec {got:14.7g}  reference {ref:16.9g}  rel {rel:.2e}  {status}')
    if failures:
        print('FAIL')
        for f in failures:
            print('  -', f)
        return 1
    print('PASS: screened Coulomb energies match the double-precision reference')
    return 0


if __name__ == '__main__':
    sys.exit(main())
