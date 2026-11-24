"""Replay the 480 held-out inputs; compare coordinates and scientific report fields with the references.

  python replay.py --output <new dir>                 # 20 stratified smoke cases
  python replay.py --all --output <new dir>           # all 480 cases
  python replay.py --all --coordinate-tolerance 0.001 --output <new dir>   # a build made with another compiler

By default the output PDB of every case must equal the reference byte for byte (its SHA-256 is in
examples/heldout480/manifest.json). Builds made with another compiler (for example GCC on Linux) can round
a coordinate differently in the last printed digit, 0.001 A. With --coordinate-tolerance T > 0, a PDB that
differs from the reference must still have the same records in the same order, with every coordinate
within T of the reference (examples/heldout480/<case>/expected_oracle_v2_pbctopo.pdb), and the
full-precision checksum of the wrapped input coordinates (wrapped_coordinate_hash) is recorded instead of
compared. Every other report field is compared exactly in both modes.
"""
import argparse, csv, hashlib, json, os, subprocess, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parent
# Recorded, not compared: the search route name (equal outputs can be reached by different routes).
RECORDED = ('mode',)
# Also recorded, not compared, with --coordinate-tolerance: a checksum over the last bits of the input.
FULL_PRECISION = ('wrapped_coordinate_hash',)


def default_binary():
    names = ['titan-vibe.exe'] if os.name == 'nt' else ['titan-vibe']
    candidates = [ROOT / n for n in names] + [ROOT / 'source' / n for n in names]
    if os.name != 'nt':
        candidates.reverse()  # a Linux build is made in source/
    return next((c for c in candidates if c.is_file()), candidates[0])


def compare_pdb(got, want, tol):
    """(None, worst, n) if the records agree within tol; otherwise the first difference."""
    a, b = got.read_text().splitlines(), want.read_text().splitlines()
    if len(a) != len(b):
        return {'reason': f'{len(a)} lines, reference has {len(b)}'}, 0.0, 0
    worst, atoms = 0.0, 0
    for i, (x, y) in enumerate(zip(a, b)):
        if x == y:
            continue
        if not (x.startswith(('ATOM  ', 'HETATM')) and y.startswith(('ATOM  ', 'HETATM'))) \
                or x[:30] != y[:30] or x[54:] != y[54:]:
            return {'reason': f'line {i + 1} differs outside the coordinates', 'actual': x, 'expected': y}, worst, atoms
        d = max(abs(float(x[k:k + 8]) - float(y[k:k + 8])) for k in (30, 38, 46))
        worst = max(worst, d)
        atoms += 1
        if d > tol + 1e-9:
            return {'reason': f'line {i + 1}: coordinate differs by {d:.4f} A', 'actual': x, 'expected': y}, worst, atoms
    return None, worst, atoms


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--all', action='store_true', help='Run all 480 held-out inputs; default: 20 stratified smoke cases')
    p.add_argument('--output', required=True, type=Path, help='New output directory; existing directories are refused')
    p.add_argument('--binary', type=Path, default=default_binary(),
                   help='Default: titan-vibe.exe (Windows) or source/titan-vibe (Linux) in this package')
    p.add_argument('--jobs', type=int, default=4, help='Independent one-thread processes (default 4)')
    p.add_argument('--coordinate-tolerance', type=float, default=0.0, metavar='A',
                   help='0 (default): output PDBs must equal the references byte for byte; '
                        'T > 0: coordinates may differ by at most T Angstrom (0.001 = one unit of the last digit)')
    a = p.parse_args()
    if not a.binary.is_file():
        raise SystemExit(f'{a.binary}: no such executable (give --binary)')
    tol = a.coordinate_tolerance
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    cases = json.loads((ROOT / 'examples/heldout480/manifest.json').read_text())
    skip = RECORDED + (FULL_PRECISION if tol > 0 else ())
    results = []
    binary_hash = hashlib.sha256(a.binary.read_bytes()).hexdigest()

    def run_case(case):
        target = out / case['case_id']
        target.mkdir(parents=True)
        inp = ROOT / 'examples/heldout480' / case['case_id']
        text = (inp / 'arg.template').read_text().replace('@FF@', (ROOT / 'forcefields').as_posix()) \
            .replace('@CASE@', inp.as_posix())
        (target / 'arg.in').write_text(text)
        t = time.perf_counter()
        with (target / 'stdout.txt').open('wb') as so, (target / 'stderr.txt').open('wb') as se:
            run = subprocess.run([str(a.binary.resolve()), '-a', 'arg.in'], cwd=target, stdout=so, stderr=se,
                                 timeout=180)
        result = {'case_id': case['case_id'], 'exit_code': run.returncode, 'elapsed_s': time.perf_counter() - t}
        pdb = target / 'oracle_v2_pbctopo.pdb'
        report = target / 'pbctopo_frame_report.csv'
        result['coordinate_hash_matches'] = pdb.is_file() and \
            hashlib.sha256(pdb.read_bytes()).hexdigest() == case['expected_pdb_sha256']
        coordinates_ok = result['coordinate_hash_matches']
        if pdb.is_file() and not coordinates_ok and tol > 0:
            problem, worst, atoms = compare_pdb(pdb, inp / 'expected_oracle_v2_pbctopo.pdb', tol)
            result['max_coordinate_difference_A'] = round(worst, 6)
            result['atoms_with_coordinate_difference'] = atoms
            if problem:
                result['coordinate_problem'] = problem
            coordinates_ok = problem is None
        result['coordinates_match'] = coordinates_ok
        if report.is_file():
            current = next(csv.DictReader(report.open()))
            result['field_mismatches'] = {k: {'expected': v, 'actual': current.get(k)}
                                          for k, v in case['expected_fields'].items()
                                          if k not in skip and current.get(k) != v}
            result['route'] = {'expected': case['expected_fields'].get('mode'), 'actual': current.get('mode')}
            if tol > 0:
                result['wrapped_coordinate_hash_matches'] = \
                    current.get('wrapped_coordinate_hash') == case['expected_fields'].get('wrapped_coordinate_hash')
        else:
            result['field_mismatches'] = {'report': 'missing'}
        result['pass'] = run.returncode == 0 and coordinates_ok and not result['field_mismatches']
        return result

    selected = [c for c in cases if a.all or c['smoke']]
    summary = {}
    with ThreadPoolExecutor(max_workers=max(1, a.jobs)) as pool:
        pending = [pool.submit(run_case, c) for c in selected]
        for task in as_completed(pending):
            result = task.result()
            results.append(result)
            summary = {'byte_identical_pdbs': sum(r['coordinate_hash_matches'] for r in results)}
            if tol > 0:
                summary['pdbs_within_tolerance_only'] = sum(r['coordinates_match'] and not r['coordinate_hash_matches']
                                                            for r in results)
                summary['max_coordinate_difference_A'] = max(
                    [r.get('max_coordinate_difference_A', 0.0) for r in results])
                summary['wrapped_coordinate_hash_differs'] = sum(
                    r.get('wrapped_coordinate_hash_matches') is False for r in results)
            (out / 'REPLAY_RESULT.json').write_text(json.dumps(
                {'purpose': 'Implementation replay against archived coordinates and scientific report fields; '
                            'not a new independent truth assessment',
                 'binary_sha256': binary_hash, 'coordinate_tolerance_A': tol,
                 'all_pass': all(x['pass'] for x in results), 'completed': len(results), 'expected': len(selected),
                 'summary': summary, 'results': sorted(results, key=lambda r: r['case_id'])}, indent=2))
            print(len(results), result['case_id'], 'PASS' if result['pass'] else 'FAIL', flush=True)
    print('RESULT:', 'PASS' if all(x['pass'] for x in results) else 'FAIL', summary)
    return 0 if all(x['pass'] for x in results) else 1


if __name__ == '__main__':
    raise SystemExit(main())
