"""ENE-11 oracle: pair-specific Lennard-Jones parameters ([ nonbond_params ],
NBFIX) must reach every energy path, for every force field.

Systems: the H10 dimer as pdb2gmx variants (GROMOS 54A7, OPLS-AA,
AMBER99SB-ILDN; test/fixtures/e2e_h10/topology_variants) and as CHARMM36
(test/fixtures/e2e_h10), one frame each (-f <pdb>), chain A - chain B.

For each force field the oracle writes two lib_ff tables: "plain" without
[ nonbond_params ] and "override" with them. GROMOS takes the GROMACS rows of
the fixture; the others get three pair overrides chosen among the type pairs
that contribute most to the chain A - chain B vdW energy (sigma x 1.05 and
epsilon x 1.6 of the combination-rule value), CHARMM36 also the NBFIX rows of
charmm36_subset.ff/nbfix.itp. Each energy path then runs with both tables:

  default    lib_ff, no cutoff           *_intrin_nofilter kernels
  cutoff     -denfo 10 (1 nm)            *_intrin kernels
  screening  -elec-screening (AMBER/CHARMM), with and without the cutoff
  force      -cal-force                  cal_energy_force_* (force.cpp)
  edecomp    -edecomp residue            titan_energy_decomposition.cpp
  gbsa       -gbsa -igb 5                gbsa_intrin.cpp
  topology   -common_topo (titan_nb), and -common_topo -gbsa, with the
             [ nonbond_params ] in the subset ffnonbonded.itp (variants)

and its vdW energy must equal the Python pair sum of lj_reference.py with
the same per-atom parameters (the run's TITAN_PARAMETER_AUDIT) and the same
pair table, both without and with the overrides; the change the overrides
make must be at least 1 kJ/mol, so that a path that ignores them fails.

  python test/nonbond_params_oracle.py --binary titan.exe --work-root test/tmp/nonbond_params

Exit code 0 when every check passes, 1 otherwise.
"""
import argparse
import os
import shutil
import sys
from pathlib import Path

sys.dont_write_bytecode = True  # fixture tools are imported from a tracked folder
sys.path.insert(0, str(Path(__file__).resolve().parent))

from topology_oracles import lj_reference as ljr  # noqa: E402
from topology_oracles.gromacs_reference import (COMB_RULE, FIXTURE, LIBFF_NAME,  # noqa: E402
                                                index_groups, lj_pairs, lj_types,
                                                pdb_coordinates, variant_dir)
from topology_oracles.titan_run import load_audit, matrix_total, run, vdw_energies  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
FORCE_FIELDS = ['amber99sb-ildn', 'oplsaa', 'gromos54a7', 'charmm36']
CUTOFF_A = 10.0
KERNELS = {
    'default': [],
    'cutoff': [f'-denfo={CUTOFF_A:g}'],
    'screening': ['-elec-screening', '-saltcon=0.15'],
    'screening_cutoff': ['-elec-screening', '-saltcon=0.15', f'-denfo={CUTOFF_A:g}'],
    'force': ['-cal-force'],
    'edecomp': ['-edecomp=true', '-edecomp-level=residue'],
    'gbsa': ['-gbsa=true', '-igb=5'],
}
SCREENING_FF = {'amber99sb-ildn', 'charmm36'}
# Paths whose vdW is wrong for a reason other than [ nonbond_params ]; they are
# reported, not run, until those issues are fixed.
KNOWN = {}   # ENE-16 and ENE-19 fixed (22/10/2026): every path runs
TOL_ABS, TOL_REL = 2e-3, 2e-6   # kJ/mol; cal.xvg prints 7 significant digits
MIN_EFFECT = 1.0                # kJ/mol


def system(ff):
    if ff == 'charmm36':
        return dict(pdb=FIXTURE / 'protein_ref.pdb', ndx=FIXTURE / 'prot.ndx',
                    libff=FIXTURE / 'lib_ff_charmm36', topology=None)
    vdir = variant_dir(ff)
    return dict(pdb=vdir / 'conf.pdb', ndx=vdir / 'index.ndx', libff=vdir / 'lib_ff_subset',
                topology=vdir)


def base_args(ff, sysd, extra):
    return [f'-ff={ff}', f'-f={sysd["pdb"].as_posix()}', f'-n={sysd["ndx"].as_posix()}',
            '-g=1-2', '-eu=kJ/mol', '-noprint-side_chain_contact=true', '-nt=1'] + extra


def libff_args(libff_dir):
    path = Path(libff_dir).as_posix() + '/'
    return [f'-path-a={path}', f'-path-n={path}']


def without_nonbond_params(text):
    out, skip = [], False
    for line in text.splitlines():
        s = line.split(';', 1)[0].strip()
        if s.startswith('['):
            skip = s.strip('[] ').strip().lower() == 'nonbond_params'
        if not skip:
            out.append(line)
    return '\n'.join(out).rstrip() + '\n'


def write_table(path, text, rows):
    text = without_nonbond_params(text)
    if rows:
        text += '\n[ nonbond_params ]\n' + ''.join(f'{i} {j} 1 {v:.12g} {w:.12g}\n'
                                                   for i, j, v, w in rows)
    Path(path).write_text(text, encoding='utf-8', newline='\n')


def write_libff(src, dest, ff, rows):
    dest.mkdir(parents=True)
    name = LIBFF_NAME[ff]
    shutil.copyfile(src / f'{name}_aminoacids.rtp', dest / f'{name}_aminoacids.rtp')
    dat = f'ff{name}nonbonded.dat'
    write_table(dest / dat, (src / dat).read_text(encoding='utf-8'), rows)


def write_topology(src, dest, ff, rows):
    shutil.copytree(src, dest, ignore=shutil.ignore_patterns('lib_ff_subset'))
    itp = dest / f'{ff}_subset.ff' / 'ffnonbonded.itp'
    text = itp.read_text(encoding='utf-8')
    defined = set(lj_types(text))   # grompp rejects pairs of undefined types
    write_table(itp, text, [r for r in rows if r[0] in defined and r[1] in defined])
    # The AMBER forcefield.itp keeps GROMACS's "*****" banner before the first
    # directive: titan_nb must ignore it as grompp does (ENE-18).


def atoms_of(audit, molecule):
    return [(r['atom_type'], float(r['sigma_nm']), float(r['epsilon_kjmol']))
            for r in audit if r['molecule'] == molecule]


def choose_overrides(ff, sysd, rule, atoms_a, atoms_b, dist):
    """(type_i, type_j, v, w) rows of the override table."""
    if ff == 'gromos54a7':
        text = (sysd['libff'] / f'ff{LIBFF_NAME[ff]}nonbonded.dat').read_text(encoding='utf-8')
        pairs = lj_pairs(text)
        return sorted({(min(i, j), max(i, j), v, w) for (i, j), (v, w) in pairs.items()})
    params = {}
    for t, v, w in atoms_a + atoms_b:
        params[t] = (v, w)
    contrib = ljr.type_pair_energies(atoms_a, atoms_b, dist, rule)
    chosen = sorted((e, k) for k, e in contrib.items()
                    if params[k[0]][1] > 0 and params[k[1]][1] > 0)[:3]
    rows = []
    for _, (ti, tj) in chosen:
        v, w = ljr.mix(rule, *params[ti], *params[tj])
        rows.append((ti, tj, v * 1.05, w * 1.6))
    if ff == 'charmm36':
        nbfix = (FIXTURE / 'charmm36_subset.ff' / 'nbfix.itp').read_text(encoding='utf-8')
        for (i, j), (v, w) in sorted(lj_pairs(nbfix).items()):
            if i <= j and i in params and j in params:
                rows.append((i, j, v, w))
    return rows


def energy_of(kernel, run_dir):
    if kernel == 'edecomp':
        # one PDB frame: no _<frame> suffix
        found = sorted(Path(run_dir).glob('edecomp_mm_residue_matrix_vdw.csv*'))
        return matrix_total(found[0]) if found else None
    values = vdw_energies(run_dir)
    return values[0] if values else None


def close(a, b):
    return abs(a - b) <= TOL_ABS + TOL_REL * max(abs(a), abs(b))


def check_ff(binary, work, ff, failures):
    sysd = system(ff)
    rule = COMB_RULE[ff]
    groups = list(index_groups(sysd['ndx']).values())
    xyz = pdb_coordinates(sysd['pdb'])
    xyz_a = [xyz[i - 1] for i in groups[1]]
    xyz_b = [xyz[i - 1] for i in groups[2]]
    dist = ljr.distances(xyz_a, xyz_b)

    plain_lib = work / 'lib_plain'
    write_libff(sysd['libff'], plain_lib, ff, [])
    ok, _ = run(binary, work / 'probe', base_args(ff, sysd, libff_args(plain_lib)))
    if not ok:
        failures.append(f'{ff}: lib_ff run did not finish ({work / "probe"})')
        return
    audit = load_audit(work / 'probe' / 'params.tsv')
    atoms_a, atoms_b = atoms_of(audit, 'mol1'), atoms_of(audit, 'mol2')
    rows = choose_overrides(ff, sysd, rule, atoms_a, atoms_b, dist)
    over_lib = work / 'lib_override'
    write_libff(sysd['libff'], over_lib, ff, rows)
    table = {}
    for i, j, v, w in rows:
        table[(i, j)] = table[(j, i)] = (v, w)

    refs = {}
    for cut in (None, CUTOFF_A * 0.1):
        refs[(cut, False)] = ljr.interaction(atoms_a, atoms_b, dist, rule, None, cut)
        refs[(cut, True)] = ljr.interaction(atoms_a, atoms_b, dist, rule, table, cut)
    print(f'{ff}: {len(rows)} override rows; reference vdW {refs[(None, False)]:.4f} -> '
          f'{refs[(None, True)]:.4f} kJ/mol')
    if abs(refs[(None, True)] - refs[(None, False)]) < MIN_EFFECT:
        failures.append(f'{ff}: the overrides change the reference by less than {MIN_EFFECT} kJ/mol')

    cases = [k for k in KERNELS if not k.startswith('screening') or ff in SCREENING_FF]
    routes = [('lib_ff', k) for k in cases]
    if sysd['topology'] is not None:
        routes += [('topology', 'default'), ('topology', 'gbsa')]
        write_topology(sysd['topology'], work / 'topo_plain', ff, [])
        write_topology(sysd['topology'], work / 'topo_override', ff, rows)
    for route, kernel in routes:
        if (ff, route, kernel) in KNOWN:
            print(f'  {route:8s} {kernel:16s} not run: {KNOWN[(ff, route, kernel)]}')
            continue
        cut = CUTOFF_A * 0.1 if 'cutoff' in kernel else None
        got = {}
        for override in (False, True):
            if route == 'lib_ff':
                extra = libff_args(over_lib if override else plain_lib)
            else:
                top = work / ('topo_override' if override else 'topo_plain') / 'topol.top'
                extra = [f'-common_topo={top.as_posix()}']
            run_dir = work / f'{route}_{kernel}_{"override" if override else "plain"}'
            ok, _ = run(binary, run_dir, base_args(ff, sysd, KERNELS[kernel] + extra))
            got[override] = energy_of(kernel, run_dir) if ok else None
        name = f'{ff} {route} {kernel}'
        if got[False] is None or got[True] is None:
            failures.append(f'{name}: run did not finish or wrote no vdW energy')
            continue
        ref0, ref1 = refs[(cut, False)], refs[(cut, True)]
        status = []
        if not close(got[False], ref0):
            status.append(f'plain {got[False]:.4f} vs reference {ref0:.4f}')
        if not close(got[True], ref1):
            status.append(f'override {got[True]:.4f} vs reference {ref1:.4f}')
        if status:
            failures.append(f'{name}: ' + '; '.join(status))
        print(f'  {route:8s} {kernel:16s} plain {got[False]:12.4f}  override {got[True]:12.4f}  '
              f'(reference {ref0:12.4f} / {ref1:12.4f})  {"FAIL" if status else "ok"}')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--binary', default='titan.exe' if os.name == 'nt' else 'titan')
    ap.add_argument('--work-root', default=str(ROOT / 'test' / 'tmp' / 'nonbond_params'))
    ap.add_argument('--ff', nargs='+', default=FORCE_FIELDS, choices=FORCE_FIELDS)
    a = ap.parse_args()
    binary = Path(a.binary).resolve() if Path(a.binary).exists() else Path(a.binary)
    work = Path(a.work_root).resolve()
    if work.exists():
        shutil.rmtree(work)
    failures = []
    for ff in a.ff:
        check_ff(binary, work / ff, ff, failures)
    if failures:
        print('FAIL')
        for f in failures:
            print('  -', f)
        return 1
    print('PASS: [ nonbond_params ] reach every energy path')
    return 0


if __name__ == '__main__':
    sys.exit(main())
