"""The GROMACS side of a topology variant (test/fixtures/e2e_h10/
topology_variants/<ff>): per-atom type, charge and mass from the itp files
pdb2gmx wrote, LJ parameters of each type and pair from ffnonbonded.itp, and
the coordinates of conf.pdb."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / 'test' / 'fixtures' / 'e2e_h10'
VARIANTS = FIXTURE / 'topology_variants'
sys.path.insert(0, str(FIXTURE / 'tools'))

from gmx_files import itp_atoms, section_rows  # noqa: E402

# GROMACS combination rule: 1 = C6/C12 geometric, 2 = Lorentz-Berthelot,
# 3 = geometric sigma and epsilon.
COMB_RULE = {'gromos54a7': 1, 'amber99sb-ildn': 2, 'oplsaa': 3, 'charmm36': 2}
LIBFF_NAME = {'gromos54a7': 'GROMOS54A7', 'oplsaa': 'OPLSAA',
              'amber99sb-ildn': 'AMBER99SB-ILDN', 'charmm36': 'CHARMM36'}
CHAINS = ['A', 'B']


def variant_dir(ff):
    return VARIANTS / ff


def topology_atoms(ff):
    """[ atoms ] of chain A then chain B (the order of conf.pdb)."""
    atoms = []
    for chain in CHAINS:
        text = (variant_dir(ff) / f'topol_Protein_chain_{chain}.itp').read_text(encoding='utf-8')
        atoms.append(itp_atoms(text))
    return atoms


def lj_types(ffnonbonded_text):
    """{type: (v, w)}: the last two columns of [ atomtypes ] (sigma/epsilon,
    or C6/C12 for GROMOS)."""
    return {cols[0]: (float(cols[-2]), float(cols[-1]))
            for cols in section_rows(ffnonbonded_text.splitlines(), 'atomtypes')}


def lj_pairs(text):
    """{(type_i, type_j): (v, w)} of [ nonbond_params ], both orders."""
    pairs = {}
    for cols in section_rows(text.splitlines(), 'nonbond_params'):
        v, w = float(cols[3]), float(cols[4])
        pairs[(cols[0], cols[1])] = (v, w)
        pairs[(cols[1], cols[0])] = (v, w)
    return pairs


def pdb_coordinates(pdb_path):
    """[(x, y, z)] in nm of the ATOM/HETATM records, in file order."""
    xyz = []
    for line in Path(pdb_path).read_text(encoding='utf-8', errors='replace').splitlines():
        if line.startswith(('ATOM', 'HETATM')):
            xyz.append((float(line[30:38]) * 0.1, float(line[38:46]) * 0.1,
                        float(line[46:54]) * 0.1))
    return xyz


def index_groups(ndx_path):
    """{group name: [1-based atom numbers]} in file order."""
    groups, current = {}, None
    for line in Path(ndx_path).read_text(encoding='utf-8').splitlines():
        s = line.strip()
        if s.startswith('['):
            current = s.strip('[] ').strip()
            groups[current] = []
        elif s and current is not None:
            groups[current] += [int(x) for x in s.split()]
    return groups
