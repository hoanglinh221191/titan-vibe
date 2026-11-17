"""ENE-08 / ENE-09 oracle: the lib_ff parameter path (no -common_topo) must
reproduce the force-field topology atom by atom.

Part 1, CHARMM36 (test/fixtures/e2e_h10, 627 + 627 atoms, three frames).
Two runs export the initialised per-atom parameters through
TITAN_PARAMETER_AUDIT:

  topology  -common_topo topol.top          charges/types from the topology
  lib_ff    -path-a/-path-n lib_ff_charmm36  charges/types from the rtp + .dat

The lib_ff run must give the same atom type, charge, mass, sigma, epsilon and
radius as the topology run for every atom (ENE-08: no type may have been
shadowed by a [ pairtypes ] row; ENE-09: the N-terminal NH3+ and C-terminal
COO- atoms must carry the CHARMM patch charges), no atom may stay UNASSIGNED,
and the total charge of each chain must be an integer. A third run renames the
histidines to an AMBER protonation name (HID) to exercise the assigner's
histidine fallback, which must yield the same parameters.

Part 2, the pdb2gmx variants (test/fixtures/e2e_h10/topology_variants/<ff>:
GROMOS 54A7, OPLS-AA, AMBER99SB-ILDN; tools/make_variants.py). A lib_ff run
on conf.pdb with the variant's lib_ff_subset must give, for every atom, the
type, charge and mass that pdb2gmx wrote to the itp, and the LJ parameters of
ffnonbonded.itp for that type (ENE-09: GROMOS and OPLS-AA terminal patches
incl. the GLY/PRO forms, the O1/O2 and OC1/OC2 carboxylate names, the
histidine names of each family, the AMBER N-/C-prefixed residues).

  python test/fflib_parameter_audit_oracle.py --binary titan.exe --work-root test/tmp/fflib_audit
  python test/fflib_parameter_audit_oracle.py --binary titan.exe --libff lib_ff   # full CHARMM36 library

Exit code 0 when every check passes, 1 otherwise.
"""
import argparse
import os
import shutil
import sys
from pathlib import Path

sys.dont_write_bytecode = True  # fixture tools are imported from a tracked folder
sys.path.insert(0, str(Path(__file__).resolve().parent))

from topology_oracles.gromacs_reference import (COMB_RULE, lj_types, topology_atoms,  # noqa: E402
                                                variant_dir)
from topology_oracles.titan_run import load_audit, run  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / 'test' / 'fixtures' / 'e2e_h10'
VARIANTS = ['gromos54a7', 'oplsaa', 'amber99sb-ildn']

BASE = ['-ff=charmm36', '-s={pdb}', '-n={fx}/prot.ndx', '-eu=kCal/mol',
        '-noprint-side_chain_contact=true', '-nt=1', '-traj={fx}/h10_window.xtc',
        '-g=1-2', '-b=8500', '-e=8600', '-dt=50']
# Charges, sigma and radius must agree to rounding. epsilon and mass are printed
# with fewer digits in the lib_ff tables than in ffnonbonded.itp (0.19246 vs
# 0.192464 kJ/mol, 15.9994 vs 15.999 u), so those two carry the file precision.
TOL = dict(charge_e=1e-6, mass_u=1e-3, sigma_nm=1e-9, epsilon_kjmol=1e-5, radius_A=1e-6)
TWO_SIXTH = 2.0 ** (1.0 / 6.0)
# Masses pdb2gmx takes from a .tdb that differ from the atom type: the OPLS-AA
# NH3+ patch of GROMACS 2026.3 (aminoacids.n.tdb) gives opls_287 14.0027 u,
# its [ atomtypes ] row 14.00670 u. TITAN uses the atom-type mass.
TDB_MASSES = {('oplsaa', 'opls_287'): 14.0067}


def label(r):
    return f"{r['molecule']} {r['residue']}{r['residue_number']}/{r['atom']}"


def common_checks(name, rows, failures):
    unassigned = [r for r in rows if r['atom_type'] == 'UNASSIGNED']
    if unassigned:
        failures.append(f'{name}: {len(unassigned)} UNASSIGNED atoms, e.g. '
                        + ', '.join(label(r) for r in unassigned[:6]))
    sums = {}
    for r in rows:
        sums[r['molecule']] = sums.get(r['molecule'], 0.0) + float(r['charge_e'])
    for mol, q in sums.items():
        if abs(q - round(q)) > 1e-6:
            failures.append(f'{name}: total charge of {mol} is {q:.6f} e, not an integer')
    return sums


# ---------------------------------------------------------------- part 1

def charmm_args(args, pdb, fixture, libff):
    return [a.replace('{fx}', fixture.as_posix()).replace('{pdb}', pdb.as_posix())
             .replace('{libff}', libff.as_posix().rstrip('/') + '/') for a in args]


def compare_runs(name, ref, test, failures):
    if len(ref) != len(test):
        failures.append(f'{name}: {len(test)} atoms, topology run has {len(ref)}')
        return
    zero = [r for r in test if float(r['sigma_nm']) == 0.0 or float(r['epsilon_kjmol']) == 0.0]
    if zero:
        failures.append(f'{name}: {len(zero)} atoms with sigma or epsilon = 0 (ENE-08), e.g. '
                        + ', '.join(f"{label(r)}:{r['atom_type']}" for r in zero[:6]))
    worst = {k: 0.0 for k in TOL}
    bad = []
    for a, b in zip(ref, test):
        if a['atom'] != b['atom'] or a['molecule'] != b['molecule']:
            bad.append(f"atom order differs at {label(a)} vs {b['atom']}")
            continue
        # residue names are not compared: the PDB reader renames HIS to the
        # protonation name only on some paths, the parameters are what matters
        if a['atom_type'] != b['atom_type']:
            bad.append(f"{label(a)}: type {b['atom_type']} vs {a['atom_type']}")
        for key, tol in TOL.items():
            d = abs(float(a[key]) - float(b[key]))
            worst[key] = max(worst[key], d)
            if d > tol:
                bad.append(f"{label(a)}: {key} {b[key]} vs {a[key]}")
    if bad:
        failures.append(f'{name}: {len(bad)} parameter differences vs the topology run, e.g. '
                        + '; '.join(bad[:8]))
    print(f'  {name}: largest |difference| vs topology: ' +
          ', '.join(f'{k}={v:.2e}' for k, v in worst.items()))


def check_charmm(binary, work, libff, failures):
    fixture = work / 'fixture'
    shutil.copytree(FIXTURE, fixture, ignore=shutil.ignore_patterns(
        'h10_long.xtc', 'lib_ff_charmm36', 'topology_variants', 'tools'))
    pdb = fixture / 'protein_ref.pdb'

    print('CHARMM36 topology run ...')
    ok, _ = run(binary, work / 'topology',
                charmm_args(BASE + ['-common_topo={fx}/topol.top'], pdb, fixture, libff))
    if not ok:
        failures.append('CHARMM36 topology run did not finish')
        return
    topo = load_audit(work / 'topology' / 'params.tsv')
    topo_q = common_checks('topology', topo, failures)

    print('CHARMM36 lib_ff run ...')
    ok, stdout = run(binary, work / 'lib_ff',
                     charmm_args(BASE + ['-path-a={libff}', '-path-n={libff}'], pdb, fixture, libff))
    if not ok:
        failures.append('CHARMM36 lib_ff run did not finish (see lib_ff/stdout.log, stderr.log)')
    else:
        lib = load_audit(work / 'lib_ff' / 'params.tsv')
        if 'rows outside [ atomtypes ]' not in stdout:
            failures.append('lib_ff run did not report skipped [ pairtypes ] rows: the fixture no '
                            'longer reproduces the ENE-08 table layout')
        compare_runs('lib_ff', topo, lib, failures)
        for mol, q in common_checks('lib_ff', lib, failures).items():
            if abs(q - topo_q[mol]) > 1e-6:
                failures.append(f'lib_ff: total charge of {mol} is {q:.6f} e '
                                f'(topology {topo_q[mol]:.6f} e)')
            print(f'  lib_ff: total charge {mol} = {q:+.3f} e (topology {topo_q[mol]:+.3f} e)')

    print('CHARMM36 lib_ff run with histidines renamed to HID ...')
    text = pdb.read_text(encoding='utf-8', errors='replace')
    his_pdb = fixture / 'protein_ref_hid.pdb'
    his_pdb.write_text(''.join(line[:17] + 'HID' + line[20:]
                               if line.startswith(('ATOM', 'HETATM')) and line[17:20] == 'HIS'
                               else line for line in text.splitlines(keepends=True)),
                       encoding='utf-8')
    ok, _ = run(binary, work / 'lib_ff_hid',
                charmm_args(BASE + ['-path-a={libff}', '-path-n={libff}'], his_pdb, fixture, libff))
    if not ok:
        failures.append('lib_ff run with HID residue names did not finish')
    else:
        compare_runs('lib_ff (HID names)', topo, load_audit(work / 'lib_ff_hid' / 'params.tsv'),
                     failures)


# ---------------------------------------------------------------- part 2

def reference_radius(rule, v, w):
    """take_para: r_min/2 in Angstrom; GROMOS from C6/C12, 0 without C6."""
    if rule == 1:
        return 0.0 if v < 1e-13 else 5.0 * (2.0 * w / v) ** (1.0 / 6.0)
    return 5.0 * TWO_SIXTH * v


def check_variant(binary, work, ff, failures):
    """The lib_ff route (ENE-09) and the -common_topo route (ENE-16: OPLS-AA
    8-column [ atomtypes ], radius by sigma except for GROMOS) of one variant."""
    vdir = variant_dir(ff)
    libff = (vdir / 'lib_ff_subset').as_posix() + '/'
    base = [f'-ff={ff}', f'-f={(vdir / "conf.pdb").as_posix()}',
            f'-n={(vdir / "index.ndx").as_posix()}', '-g=1-2', '-eu=kJ/mol',
            '-noprint-side_chain_contact=true', '-nt=1']
    routes = [('lib_ff', base + [f'-path-a={libff}', f'-path-n={libff}']),
              ('topology', base + [f'-common_topo={(vdir / "topol.top").as_posix()}'])]
    for route, args in routes:
        name = f'{ff} {route}'
        print(f'{name} run ...')
        ok, stdout = run(binary, work / ff / route, args)
        if not ok:
            failures.append(f'{name}: run did not finish (see {ff}/{route}/stdout.log)')
            continue
        warnings = ['not fully assigned', 'not an integer', 'has no [ nonbond_params ]',
                    'could not be read', 'names differ']
        for warning in warnings:
            if warning in stdout:
                failures.append(f'{name}: run printed "{warning}"')
        compare_variant(name, ff, route, load_audit(work / ff / route / 'params.tsv'), failures)


def compare_variant(name, ff, route, rows, failures):
    vdir = variant_dir(ff)
    chains = topology_atoms(ff)
    ref = [a for chain in chains for a in chain]
    lj = lj_types((vdir / f'{ff}_subset.ff' / 'ffnonbonded.itp').read_text(encoding='utf-8'))
    rule = COMB_RULE[ff]
    if len(rows) != len(ref):
        failures.append(f'{name}: {len(rows)} atoms, the itp files have {len(ref)}')
        return
    per_mol = [sum(1 for r in rows if r['molecule'] == m) for m in ('mol1', 'mol2')]
    if per_mol != [len(c) for c in chains]:
        failures.append(f'{name}: atoms per molecule {per_mol}, chains {[len(c) for c in chains]}')
    bad, worst = [], dict(charge_e=0.0, mass_u=0.0, lj=0.0, radius_A=0.0)
    for r, a in zip(rows, ref):
        where = f"{r['molecule']} {a['residue']}{a['resnr']}/{a['atom']}"
        if r['atom'] != a['atom']:
            bad.append(f"{where}: atom {r['atom']}")
            continue
        if r['atom_type'] != a['type']:
            bad.append(f"{where}: type {r['atom_type']} vs {a['type']}")
            continue
        v, w = lj[a['type']]
        # the topology route takes the itp mass, the lib_ff route the atom type's
        mass = a['mass'] if route == 'topology' else TDB_MASSES.get((ff, a['type']), a['mass'])
        diffs = dict(charge_e=abs(float(r['charge_e']) - a['charge']),
                     mass_u=abs(float(r['mass_u']) - mass),
                     lj=max(abs(float(r['sigma_nm']) - v) / max(abs(v), 1e-30),
                            abs(float(r['epsilon_kjmol']) - w) / max(abs(w), 1e-30)),
                     radius_A=abs(float(r['radius_A']) - reference_radius(rule, v, w)))
        limits = dict(charge_e=1e-6, mass_u=1e-3, lj=1e-12, radius_A=1e-9)
        for key, d in diffs.items():
            worst[key] = max(worst[key], d)
            if d > limits[key]:
                bad.append(f'{where} ({a["type"]}): {key} differs by {d:.3g}')
    if bad:
        failures.append(f'{name}: {len(bad)} differences from the pdb2gmx topology, e.g. '
                        + '; '.join(bad[:8]))
    sums = common_checks(name, rows, failures)
    print(f'  {name}: largest |difference| vs pdb2gmx: '
          + ', '.join(f'{k}={v:.2e}' for k, v in worst.items())
          + '; charges ' + ', '.join(f'{m} {q:+.3f} e' for m, q in sums.items()))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--binary', default='titan.exe' if os.name == 'nt' else 'titan')
    ap.add_argument('--work-root', default=str(ROOT / 'test' / 'tmp' / 'fflib_audit'))
    ap.add_argument('--libff', default=str(FIXTURE / 'lib_ff_charmm36'),
                    help='CHARMM36 lib_ff directory (default: the subset of the H10 fixture)')
    ap.add_argument('--variants', nargs='*', default=VARIANTS, choices=VARIANTS,
                    help='pdb2gmx variants of part 2 (none: part 1 only)')
    a = ap.parse_args()
    binary = Path(a.binary).resolve() if Path(a.binary).exists() else Path(a.binary)
    work = Path(a.work_root).resolve()
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    failures = []

    check_charmm(binary, work, Path(a.libff).resolve(), failures)
    for ff in a.variants:
        check_variant(binary, work / 'variants', ff, failures)

    if failures:
        print('FAIL')
        for f in failures:
            print('  -', f)
        return 1
    print('PASS: lib_ff parameters reproduce the topology for every atom')
    return 0


if __name__ == '__main__':
    sys.exit(main())
