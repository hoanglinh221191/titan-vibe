"""STR-32 oracle: backbone torsions of frames whose box is not a valid cell.

The dihedral module checked only the box lengths. A CRYST1 with positive
a/b/c but degenerate angles (here gamma = 0) built a NaN cell matrix, and
every torsion of the frame came out NaN. Such frames must now be measured
with plain coordinates, like frames without a box.

Three copies of the same 3-frame peptide hairpin (test/qnative_oracle.py):
  no_box        no CRYST1
  degenerate    CRYST1 100 100 100, angles 90 90 0
  cubic         CRYST1 100 100 100, angles 90 90 90 (molecule whole, far
                from the faces: the minimum image changes nothing)
All three must give the same finite phi/psi, equal to the torsions computed
here from the coordinates.

  python test/dihedral_box_oracle.py --binary titan.exe --work-root <dir>
"""
import argparse
import csv
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import md_oracle_fixture as fx  # noqa: E402
import qnative_oracle as qo  # noqa: E402

CRYST = {
    'no_box': None,
    'degenerate': 'CRYST1  100.000  100.000  100.000  90.00  90.00   0.00 P 1           1',
    'cubic': 'CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1           1',
}


def write_trajectory(path, frames, cryst):
    lines = []
    for index, (time_ps, atoms) in enumerate(frames):
        if cryst:
            lines.append(cryst)
        lines.append(f'MODEL     {index + 1:4d}')
        lines.append(f'TITLE     t= {time_ps:.5f} step= {index}')
        for serial, (chain, resnr, resname, atom, x, y, z) in enumerate(atoms, 1):
            atom_field = f' {atom:<3s}' if len(atom) < 4 else atom
            lines.append(f'ATOM  {serial:5d} {atom_field} {resname:<4s}{chain}{resnr:4d}    '
                         f'{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {fx.element(atom)}')
        lines.append('TER')
        lines.append('ENDMDL')
    Path(path).write_text('\n'.join(lines) + '\nEND\n')


def torsion(p0, p1, p2, p3):
    b0 = [p0[k] - p1[k] for k in range(3)]
    b1 = [p2[k] - p1[k] for k in range(3)]
    b2 = [p3[k] - p2[k] for k in range(3)]
    n1 = math.sqrt(sum(v * v for v in b1))
    b1 = [v / n1 for v in b1]
    dot = lambda a, b: sum(x * y for x, y in zip(a, b))
    v = [b0[k] - dot(b0, b1) * b1[k] for k in range(3)]
    w = [b2[k] - dot(b2, b1) * b1[k] for k in range(3)]
    cross = [b1[1] * v[2] - b1[2] * v[1], b1[2] * v[0] - b1[0] * v[2], b1[0] * v[1] - b1[1] * v[0]]
    return math.degrees(math.atan2(dot(cross, w), dot(v, w)))


def expected_phi_psi(atoms):
    pos = {(r, n): (x, y, z) for _, r, _, n, x, y, z in atoms}
    out = {}
    for r in range(2, 10):
        if all(k in pos for k in ((r - 1, 'C'), (r, 'N'), (r, 'CA'), (r, 'C'), (r + 1, 'N'))):
            out[(r, 'phi')] = torsion(pos[(r - 1, 'C')], pos[(r, 'N')], pos[(r, 'CA')], pos[(r, 'C')])
            out[(r, 'psi')] = torsion(pos[(r, 'N')], pos[(r, 'CA')], pos[(r, 'C')], pos[(r + 1, 'N')])
    return out


def read_angles(case):
    """{(frame, resnr, kind): degrees} from the backbone time series."""
    files = sorted(Path(case).glob('*dihedral*timeseries*.csv'))
    files = [f for f in files if 'sidechain' not in f.name]
    if not files:
        return None, []
    rows = list(csv.DictReader(files[0].open()))
    return files[0].name, rows


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--binary', type=Path, required=True)
    p.add_argument('--work-root', type=Path, required=True)
    a = p.parse_args()
    a.binary = a.binary.resolve()
    work = a.work_root.resolve()
    work.mkdir(parents=True, exist_ok=True)
    frames = [(t, qo.hairpin(s)) for t, s in ((0.0, 0.0), (10.0, 1.6), (20.0, 3.0))]
    atoms0 = frames[0][1]
    top, ndx = work / 'hairpin.top', work / 'hairpin.ndx'
    fx.write_topology(top, atoms0)
    fx.write_index(ndx, [('System', list(range(1, len(atoms0) + 1)))])
    checks = fx.Checks()
    results = {}
    for name, cryst in CRYST.items():
        case = work / name
        case.mkdir(exist_ok=True)
        traj = case / 'hairpin.pdb'
        write_trajectory(traj, frames, cryst)
        code, out = fx.run_titan(a.binary, case, fx.base_args(work, traj, top, ndx) +
                                 ['-dihedral=true', '-dihedral-kind=backbone'])
        fname, rows = read_angles(case)
        checks.check(code == 0 and fname is not None, f'{name}: run and time series', f'exit {code}')
        values = []
        for row in rows:
            for key, value in row.items():
                if key and ('phi' in key.lower() or 'psi' in key.lower() or key.lower() in ('angle_deg', 'value_deg', 'degrees')):
                    try:
                        values.append(float(value))
                    except (TypeError, ValueError):
                        pass
        results[name] = (rows, values)
        # phi of the first residue and psi/omega of the last are undefined.
        last = max((int(r['resnr']) for r in rows), default=0)
        bad = 0
        for r in rows:
            for kind in ('phi_deg', 'psi_deg', 'omega_deg'):
                terminal = (kind == 'phi_deg' and int(r['resnr']) == 1) or \
                           (kind != 'phi_deg' and int(r['resnr']) == last)
                if not terminal and not math.isfinite(float(r[kind])):
                    bad += 1
        checks.check(rows and bad == 0, f'{name}: every torsion away from the termini finite',
                     f'{len(rows)} rows, {bad} NaN')
    ref_rows, ref_values = results['no_box']
    for name in ('degenerate', 'cubic'):
        rows, values = results[name]
        same = len(values) == len(ref_values) and all(
            abs(x - y) <= 1e-6 or (not math.isfinite(x) and not math.isfinite(y))
            for x, y in zip(values, ref_values))
        checks.check(same, f'{name}: torsions equal to the no-box run', f'{len(values)} vs {len(ref_values)}')
    expected = expected_phi_psi(atoms0)
    finite_ref = [v for v in ref_values if math.isfinite(v)]
    hits = sum(1 for e in expected.values() if any(abs(e - v) < 1e-3 for v in finite_ref))
    checks.check(expected and hits == len(expected),
                 'no_box: frame-0 phi/psi equal to the torsions computed from the coordinates',
                 f'{hits}/{len(expected)} matched')
    return checks.exit_code()


if __name__ == '__main__':
    sys.exit(main())
