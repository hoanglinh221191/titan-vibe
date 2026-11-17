"""Run the titan binary on an arg.in and read back what the oracles compare:
the per-atom parameter audit (TITAN_PARAMETER_AUDIT) and the energies."""
import csv
import gzip
import os
import subprocess
from pathlib import Path


def run(binary, run_dir, args, audit=True):
    """(ok, stdout) of `titan -a arg.in` in run_dir; params.tsv when audit."""
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / 'arg.in').write_text('\n'.join(args) + '\n', encoding='utf-8')
    env = dict(os.environ)
    env.pop('TITAN_PARAMETER_AUDIT', None)
    if audit:
        env['TITAN_PARAMETER_AUDIT'] = str(run_dir / 'params.tsv')
    with open(run_dir / 'stdout.log', 'wb') as so, open(run_dir / 'stderr.log', 'wb') as se:
        rc = subprocess.run([str(binary), '-a', 'arg.in'], cwd=run_dir, stdout=so, stderr=se,
                            stdin=subprocess.DEVNULL, env=env).returncode
    stdout = (run_dir / 'stdout.log').read_text(errors='replace')
    ok = rc == 0 and 'TERMINATED SUCCESSFULLY' in stdout
    if audit and not (run_dir / 'params.tsv').exists():
        ok = False
    return ok, stdout


def load_audit(path):
    with open(path, encoding='utf-8') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def xvg_rows(path):
    rows = []
    for line in Path(path).read_text(encoding='utf-8', errors='replace').splitlines():
        s = line.strip()
        if s and not s.startswith(('#', '@')):
            rows.append([float(x) for x in s.split()])
    return rows


def vdw_energies(run_dir):
    """VdW column of cal.xvg (time, elec, vdw, total), one value per frame."""
    return [row[2] for row in xvg_rows(Path(run_dir) / 'cal.xvg')]


def matrix_total(path, component='vdw'):
    """Sum of signed_mean over the rows of one component of an -edecomp
    residue matrix csv(.gz) (long format: row, column, component, means)."""
    opener = gzip.open if str(path).endswith('.gz') else open
    total = 0.0
    with opener(path, 'rt', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            if row.get('component') == component:
                total += float(row['signed_mean'])
    return total
