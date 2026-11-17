"""Build test/fixtures/martini_ene13/: two copies of a 7-bead MARTINI 2
peptide (LYS-PHE), the smallest system that exercises the MARTINI topology
reader of TITAN (src/topology/martini_*.cpp) for test/martini_radius_oracle.py.

  martini_subset.itp   [ defaults ] (comb-rule 2), the [ atomtypes ] rows of
                       P5, Qd, C3 and SC4 and the [ nonbond_params ] among them,
                       copied from martini_v2.1-dna.itp (example/5, not
                       tracked). As in that file, the atom types carry
                       sigma = epsilon = 0 and every LJ pair, the diagonal
                       included, comes from [ nonbond_params ].
  cgpep.itp            the molecule (charge +1 on the Qd bead of LYS)
  topol.top            both files and [ molecules ] CGPEP 2
  conf.pdb, index.ndx  chain A and chain B, 0.8 nm apart along x

  python test/fixtures/martini_ene13/make_fixture.py [--martini example/5/martini_v2.1-dna.itp]
"""
import argparse
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
TYPES = ['P5', 'Qd', 'C3', 'SC4']
# residue, bead, type, charge, offset (nm) of each bead within one copy
BEADS = [(1, 'LYS', 'BB', 'P5', 0.0, (0.00, 0.00, 0.00)),
         (1, 'LYS', 'SC1', 'C3', 0.0, (0.00, 0.33, 0.00)),
         (1, 'LYS', 'SC2', 'Qd', 1.0, (0.00, 0.62, 0.05)),
         (2, 'PHE', 'BB', 'P5', 0.0, (0.35, 0.00, 0.00)),
         (2, 'PHE', 'SC1', 'SC4', 0.0, (0.35, 0.31, 0.00)),
         (2, 'PHE', 'SC2', 'SC4', 0.0, (0.35, 0.50, 0.21)),
         (2, 'PHE', 'SC3', 'SC4', 0.0, (0.35, 0.58, -0.05))]
BONDS = [(1, 2), (2, 3), (1, 4), (4, 5), (5, 6), (6, 7), (5, 7)]
SHIFT = (0.80, 0.05, 0.10)   # copy B relative to copy A, nm
NOTE = '; written by test/fixtures/martini_ene13/make_fixture.py\n'


def section_lines(text, wanted):
    out, current = [], None
    for line in text.splitlines():
        s = line.split(';', 1)[0].strip()
        if s.startswith('['):
            current = s.strip('[] ').strip().lower()
            continue
        if current == wanted and s:
            out.append(line.rstrip())
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--martini', default=str(ROOT / 'example' / '5' / 'martini_v2.1-dna.itp'))
    a = ap.parse_args()
    text = Path(a.martini).read_text(encoding='utf-8', errors='replace')

    atomtypes = [l for l in section_lines(text, 'atomtypes') if l.split()[0] in TYPES]
    pairs = [l for l in section_lines(text, 'nonbond_params')
             if l.split()[0] in TYPES and l.split()[1] in TYPES]
    if len(atomtypes) != len(TYPES):
        raise SystemExit(f'atom types found: {[l.split()[0] for l in atomtypes]}')
    subset = [NOTE.rstrip(), f'; MARTINI 2.1 subset ({", ".join(TYPES)}) of {Path(a.martini).name}',
              '', '[ defaults ]', '; nbfunct  comb-rule', '  1        2', '',
              '[ atomtypes ]'] + atomtypes + ['', '[ nonbond_params ]'] + pairs
    (HERE / 'martini_subset.itp').write_text('\n'.join(subset) + '\n', encoding='utf-8', newline='\n')

    itp = [NOTE.rstrip(), '[ moleculetype ]', '; name nrexcl', 'CGPEP 1', '', '[ atoms ]',
           '; nr type resnr residue atom cgnr charge mass']
    for k, (resnr, res, bead, btype, q, _) in enumerate(BEADS, start=1):
        itp.append(f'{k:4d} {btype:>5s} {resnr:4d} {res:>5s} {bead:>4s} {k:4d} {q:7.3f} 72.0')
    itp += ['', '[ bonds ]', '; i j funct length force']
    itp += [f'{i:4d} {j:4d} 1 0.350 1250' for i, j in BONDS]
    (HERE / 'cgpep.itp').write_text('\n'.join(itp) + '\n', encoding='utf-8', newline='\n')

    top = [NOTE.rstrip(), '#include "martini_subset.itp"', '#include "cgpep.itp"', '',
           '[ system ]', 'CGPEP dimer for the ENE-13 oracle', '', '[ molecules ]', 'CGPEP 2']
    (HERE / 'topol.top').write_text('\n'.join(top) + '\n', encoding='utf-8', newline='\n')

    pdb, serial = ['REMARK    written by test/fixtures/martini_ene13/make_fixture.py'], 0
    groups = {}
    for chain, shift in (('A', (0.0, 0.0, 0.0)), ('B', SHIFT)):
        for resnr, res, bead, _, _, xyz in BEADS:
            serial += 1
            x, y, z = (10.0 * (c + s) for c, s in zip(xyz, shift))
            pdb.append(f'ATOM  {serial:5d} {bead:<4s} {res:3s} {chain}{resnr:4d}    '
                       f'{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00')
            groups.setdefault(chain, []).append(serial)
    pdb.append('END')
    (HERE / 'conf.pdb').write_text('\n'.join(pdb) + '\n', encoding='utf-8', newline='\n')
    every = [s for c in groups.values() for s in c]
    ndx = [f'[ System ]', ' '.join(map(str, every))]
    for chain, serials in groups.items():
        ndx += [f'[ chain_{chain} ]', ' '.join(map(str, serials))]
    (HERE / 'index.ndx').write_text('\n'.join(ndx) + '\n', encoding='utf-8', newline='\n')
    print(f'{len(atomtypes)} atom types, {len(pairs)} pairs, {serial} beads -> {HERE}')


if __name__ == '__main__':
    main()
