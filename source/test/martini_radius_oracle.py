"""ENE-13 oracle: MARTINI bead radii must be r_min/2 in Angstrom,
5 * 2^(1/6) * sigma_ii [nm], the unit every other force field uses
(take_para) and that the pbctopo steric test and -pqr expect.

Runs titan -ff=martini -pqr on test/fixtures/martini_ene13 (two copies of a
7-bead peptide, chain A and chain B) and checks every bead of the PQR files:

  radius   5 * 2^(1/6) * sigma of the bead type's diagonal pair, taken from
           the [ nonbond_params ] row "t t" of martini_subset.itp, else from
           [ atomtypes ] (MARTINI 2 types carry sigma = 0 there)
  charge   the [ atoms ] charge of cgpep.itp

Before ENE-13 the radius was sigma * 2^(1/6) / 2 in nm, ten times too small
(0.264 instead of 2.638 for sigma = 0.47 nm).

  python test/martini_radius_oracle.py --binary titan.exe --work-root test/tmp/martini_radius

Exit code 0 when every check passes, 1 otherwise.
"""
import argparse
import os
import shutil
import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))

from topology_oracles.titan_run import run  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / 'test' / 'fixtures' / 'martini_ene13'
TWO_SIXTH = 2.0 ** (1.0 / 6.0)
TOL_RADIUS = 1.5e-3   # the PQR prints three decimals
TOL_CHARGE = 1.5e-3


def rows_of(text, wanted):
    out, current = [], None
    for line in text.splitlines():
        s = line.split(';', 1)[0].strip()
        if s.startswith('['):
            current = s.strip('[] ').strip().lower()
        elif current == wanted and s:
            out.append(s.split())
    return out


def reference():
    """[(bead, type, charge, radius_A)] of one molecule, in [ atoms ] order."""
    ff = (FIXTURE / 'martini_subset.itp').read_text(encoding='utf-8')
    sigma = {c[0]: float(c[-2]) for c in rows_of(ff, 'atomtypes')}
    for c in rows_of(ff, 'nonbond_params'):
        if c[0] == c[1]:
            sigma[c[0]] = float(c[3])
    beads = []
    for c in rows_of((FIXTURE / 'cgpep.itp').read_text(encoding='utf-8'), 'atoms'):
        beads.append((c[4], c[1], float(c[6]), 5.0 * TWO_SIXTH * sigma[c[1]]))
    return beads


def pqr_atoms(path):
    """[(name, charge, radius)] of the ATOM records."""
    out = []
    for line in Path(path).read_text(encoding='utf-8', errors='replace').splitlines():
        if line.startswith('ATOM'):
            cols = line.split()
            out.append((cols[2], float(cols[-2]), float(cols[-1])))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--binary', default='titan.exe' if os.name == 'nt' else 'titan')
    ap.add_argument('--work-root', default=str(ROOT / 'test' / 'tmp' / 'martini_radius'))
    a = ap.parse_args()
    binary = Path(a.binary).resolve() if Path(a.binary).exists() else Path(a.binary)
    work = Path(a.work_root).resolve()
    if work.exists():
        shutil.rmtree(work)
    args = ['-ff=martini', f'-common_topo={(FIXTURE / "topol.top").as_posix()}',
            f'-f={(FIXTURE / "conf.pdb").as_posix()}', f'-n={(FIXTURE / "index.ndx").as_posix()}',
            '-g=1-2', '-pqr', '-eu=kJ/mol', '-nt=1', '-noprint-side_chain_contact=true']
    ok, _ = run(binary, work, args, audit=False)
    failures = []
    if not ok:
        failures.append(f'titan did not finish (see {work / "stdout.log"})')
    pqrs = sorted(work.rglob('*.pqr'))
    if ok and len(pqrs) != 2:
        failures.append(f'expected a receptor and a ligand PQR file, found {[p.name for p in pqrs]}')
    ref = reference()
    worst_r = worst_q = 0.0
    for pqr in pqrs:
        atoms = pqr_atoms(pqr)
        if [name for name, _, _ in atoms] != [b[0] for b in ref]:
            failures.append(f'{pqr.name}: beads {[n for n, _, _ in atoms]}, expected {[b[0] for b in ref]}')
            continue
        for (name, q, r), (_, btype, q_ref, r_ref) in zip(atoms, ref):
            worst_r, worst_q = max(worst_r, abs(r - r_ref)), max(worst_q, abs(q - q_ref))
            if abs(r - r_ref) > TOL_RADIUS:
                failures.append(f'{pqr.name} {name} ({btype}): radius {r:.3f} A, expected {r_ref:.3f} A')
            if abs(q - q_ref) > TOL_CHARGE:
                failures.append(f'{pqr.name} {name} ({btype}): charge {q:.3f} e, expected {q_ref:.3f} e')
    print(f'{len(pqrs)} PQR files, {len(ref)} beads each; largest |difference| radius '
          f'{worst_r:.2e} A, charge {worst_q:.2e} e; expected radii '
          + ', '.join(sorted({f'{t} {r:.3f} A' for _, t, _, r in ref})))
    if failures:
        print('FAIL')
        for f in failures[:20]:
            print('  -', f)
        return 1
    print('PASS: MARTINI bead radii are r_min/2 in Angstrom')
    return 0


if __name__ == '__main__':
    sys.exit(main())
