"""ENE-22 oracle: energy groups joined by a covalent bond.

Two ALA residues, chain A (group 1, N CA C O CB) and chain B (group 2,
N H CA C O CB), on a synthetic monolithic topology (no local force-field
files needed). In the bonded geometry chain B's N sits 1.33 A from chain A's
C and the topology lists that bond, so with nrexcl = 3 every pair within
three bonds across the junction is excluded.

  gbsa_bonded       -common_topo -gbsa: MM/GBSA must apply the topology's
                    exclusions. E_elec and E_vdW (cal.xvg) must equal the
                    pair sum over receptor x ligand without the excluded
                    pairs. The first ENE-22 fix refused this run; before
                    that it ran and counted the bonded pairs at 1.3-2.5 A.
  gbsa_bonded_14    as above, with two cross-junction 1-4 pairs listed in
                    [ pairs ] (gen-pairs: epsilon x fudgeLJ, charge product
                    x fudgeQQ, plain Coulomb), which must be added.
  topology_bonded   the same geometry on the -common_topo route without
                    -gbsa (topology energies): the same E_elec/E_vdW.
  gbsa_hbond        chains apart, with only an N-H...O hydrogen bond between
                    them (H...O 1.8 A, N...O 2.8 A): no exclusion applies,
                    the plain pair sum, and no invalid_bonded_groups refusal.

  python test/bonded_groups_oracle.py --binary titan.exe --work-root <dir>
"""
import argparse
import math
import subprocess
import sys
from collections import deque
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument('--binary', type=Path, default=ROOT / 'titan.exe')
p.add_argument('--work-root', type=Path, default=ROOT / 'test/tmp/bonded_groups_oracle')
a = p.parse_args()
binary = a.binary.resolve()
work = a.work_root.resolve()
work.mkdir(parents=True, exist_ok=True)

K_COULOMB = 138.935457839  # kJ mol^-1 nm e^-2
FUDGE_LJ, FUDGE_QQ = 0.5, 0.833333333333
LJ = {'X': (0.3, 0.2), 'Y': (0.1, 0.05)}  # sigma nm, epsilon kJ/mol
CHARGE = {'N': -0.3, 'H': 0.3, 'O': -0.5, 'C': 0.5}

RES_A = [('N', 0.000, 0.000, 0.000), ('CA', 1.450, 0.000, 0.000), ('C', 2.050, 1.350, 0.000),
         ('O', 1.450, 2.350, 0.000), ('CB', 1.950, -0.750, 1.250)]


def residue_b(shift, h_x=-0.500):
    """Residue B atoms N H CA C O CB, translated by `shift`; the N-H bond
    points along (h_x, -0.85)."""
    local = [('N', 0.000, 0.000, 0.000), ('H', h_x, -0.850, 0.000), ('CA', 1.450, 0.000, 0.000),
             ('C', 2.050, 1.350, 0.000), ('O', 1.450, 2.350, 0.000), ('CB', 1.950, -0.750, 1.250)]
    return [(n, x + shift[0], y + shift[1], z + shift[2]) for n, x, y, z in local]


def bonds_of(n_a, bonded):
    b0 = n_a
    bonds = [(1, 2), (2, 3), (3, 4), (2, 5),
             (b0 + 1, b0 + 2), (b0 + 1, b0 + 3), (b0 + 3, b0 + 4), (b0 + 4, b0 + 5), (b0 + 3, b0 + 6)]
    if bonded:
        bonds.append((3, b0 + 1))  # peptide C(A) - N(B)
    return bonds


def excluded_pairs(n_atoms, bonds, nrexcl=3):
    graph = {i: set() for i in range(1, n_atoms + 1)}
    for i, j in bonds:
        graph[i].add(j)
        graph[j].add(i)
    out = set()
    for start in graph:
        depth = {start: 0}
        queue = deque([start])
        while queue:
            node = queue.popleft()
            if depth[node] == nrexcl:
                continue
            for nxt in graph[node]:
                if nxt not in depth:
                    depth[nxt] = depth[node] + 1
                    queue.append(nxt)
        out.update(frozenset((start, k)) for k in depth if k != start)
    return out


def expected_energies(atoms, n_a, bonds, listed):
    """Receptor x ligand Coulomb and LJ (kJ/mol) with nrexcl = 3 exclusions
    and the listed 1-4 pairs, from first principles."""
    excl = excluded_pairs(len(atoms), bonds)
    elec = vdw = 0.0
    for i in range(1, n_a + 1):
        for j in range(n_a + 1, len(atoms) + 1):
            ni, xi, yi, zi = atoms[i - 1]
            nj, xj, yj, zj = atoms[j - 1]
            r = 0.1 * math.dist((xi, yi, zi), (xj, yj, zj))
            ti, tj = ('Y' if ni == 'H' else 'X'), ('Y' if nj == 'H' else 'X')
            sigma = 0.5 * (LJ[ti][0] + LJ[tj][0])
            eps = math.sqrt(LJ[ti][1] * LJ[tj][1])
            qq = CHARGE.get(ni, 0.0) * CHARGE.get(nj, 0.0)
            sr6 = (sigma / r) ** 6
            if frozenset((i, j)) not in excl:
                elec += K_COULOMB * qq / r
                vdw += 4 * eps * (sr6 * sr6 - sr6)
            if (i, j) in listed or (j, i) in listed:
                elec += K_COULOMB * FUDGE_QQ * qq / r
                vdw += 4 * FUDGE_LJ * eps * (sr6 * sr6 - sr6)
    return elec, vdw


def write_case(case, b_atoms, extra, bonded, listed):
    case.mkdir(parents=True, exist_ok=True)
    lines = ['CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1           1',
             'MODEL        1']
    serial = 0
    for chain, resnr, atoms in (('A', 1, RES_A), ('B', 2, b_atoms)):
        for name, x, y, z in atoms:
            serial += 1
            lines.append(f'ATOM  {serial:5d} {name:^4s} ALA {chain}{resnr:4d}    '
                         f'{x + 10:8.3f}{y + 10:8.3f}{z + 10:8.3f}  1.00  0.00           {name[0]}')
        lines.append('TER')
    lines += ['ENDMDL', 'END']
    (case / 'fixture.pdb').write_text('\n'.join(lines) + '\n', encoding='ascii')
    n_a, n_b = len(RES_A), len(b_atoms)
    (case / 'fixture.ndx').write_text(
        '[ System ]\n' + ' '.join(str(i) for i in range(1, n_a + n_b + 1)) +
        '\n[ A ]\n' + ' '.join(str(i) for i in range(1, n_a + 1)) +
        '\n[ B ]\n' + ' '.join(str(i) for i in range(n_a + 1, n_a + n_b + 1)) + '\n', encoding='ascii')
    names = [n for n, *_ in RES_A] + [n for n, *_ in b_atoms]
    top = ['[ defaults ]', f'1 2 yes {FUDGE_LJ} {FUDGE_QQ}',
           '[ atomtypes ]', 'X 6 12.0 0.0 A 0.3 0.2', 'Y 1 1.008 0.0 A 0.1 0.05',
           '[ moleculetype ]', '; name nrexcl', 'Fixture 3', '[ atoms ]']
    for i, name in enumerate(names, 1):
        resnr = 1 if i <= n_a else 2
        top.append(f'{i} {"Y" if name == "H" else "X"} {resnr} ALA {name} {i} '
                   f'{CHARGE.get(name, 0.0)} {1.008 if name == "H" else 12.0}')
    bonds = bonds_of(n_a, bonded)
    top += ['[ bonds ]'] + [f'{i} {j} 1' for i, j in bonds]
    if listed:
        top += ['[ pairs ]'] + [f'{i} {j} 1' for i, j in sorted(listed)]
    top += ['[ system ]', 'ENE-22 oracle', '[ molecules ]', 'Fixture 1']
    (case / 'fixture.top').write_text('\n'.join(top) + '\n', encoding='ascii')
    options = [f'-f={(case / "fixture.pdb").as_posix()}', '-1pdb=true', '-ff=AMBER99SB',
               f'-n={(case / "fixture.ndx").as_posix()}', '-g=1-2',
               '-mol-1-type=protein', '-mol-2-type=lig',
               f'-path-a={case.as_posix()}/', f'-path-n={case.as_posix()}/',
               f'-common_topo={(case / "fixture.top").as_posix()}', '-nt=1'] + extra
    (case / 'arg.in').write_text('\n'.join(options) + '\n', encoding='utf-8')
    run = subprocess.run([str(binary), '-a', str(case / 'arg.in')], cwd=case, capture_output=True,
                         text=True, encoding='utf-8', errors='replace', timeout=120)
    out = run.stdout + run.stderr
    (case / 'result.log').write_text(out, encoding='utf-8')
    energies = None
    xvg = case / 'cal.xvg'
    if xvg.exists():
        for line in xvg.read_text(encoding='utf-8').splitlines():
            if line.strip() and not line.lstrip().startswith(('#', '@')):
                cols = line.split()
                energies = (float(cols[1]), float(cols[2]))
                break
    # the coordinates as written to the PDB (%8.3f), not the unrounded ones
    atoms = [(n, round(x + 10, 3), round(y + 10, 3), round(z + 10, 3)) for n, x, y, z in RES_A + b_atoms]
    return run.returncode, out, energies, expected_energies(atoms, n_a, bonds, listed)


# peptide geometry: B's N 1.33 A from A's C, its H pointing away from the C
bonded_b = residue_b((2.050 + 1.330, 1.350, 0.0), h_x=0.500)
# hydrogen-bonded geometry: B's H 1.8 A from A's O, N...O 2.8 A (N-H...O linear)
o = (1.450, 2.350, 0.000)
hbond_b = residue_b((o[0] + 0.500 * 2.8 / 0.986, o[1] + 0.850 * 2.8 / 0.986, 0.0))
# cross-junction 1-4 pairs: CA(A)-C(A)-N(B)-CA(B) and O(A)-C(A)-N(B)-H(B)
LISTED_14 = {(2, 8), (4, 7)}


def close(x, y):
    return abs(x - y) <= 1e-5 * max(1.0, abs(y))


failures = []
for name, b_atoms, extra, bonded, listed in (
        ('gbsa_bonded', bonded_b, ['-gbsa=true'], True, set()),
        ('gbsa_bonded_14', bonded_b, ['-gbsa=true'], True, LISTED_14),
        ('topology_bonded', bonded_b, [], True, set()),
        ('gbsa_hbond', hbond_b, ['-gbsa=true'], False, set())):
    code, out, got, want = write_case(work / name, b_atoms, extra, bonded, listed)
    refused = 'invalid_bonded_groups' in out
    ok = code == 0 and not refused and got is not None and close(got[0], want[0]) and close(got[1], want[1])
    got_text = 'none' if got is None else f'{got[0]:.6f}/{got[1]:.6f}'
    print(f'{"PASS" if ok else "FAIL"} {name}: exit {code}, refused {refused}, '
          f'E_elec/E_vdW {got_text} (expected {want[0]:.6f}/{want[1]:.6f})')
    if not ok:
        failures.append(name)
        print(out[-1500:])
sys.exit(1 if failures else 0)
