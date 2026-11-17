"""ENE-12 oracle: topology parameters are assigned to the structure by
position only when the atom names agree.

  unit     test/topology_name_match_oracle.cpp: the hashed name matcher equals
           the former linear search on random structures; name rules
  same     AMBER99SB-ILDN variant (test/fixtures/e2e_h10/topology_variants),
           -common_topo -gbsa: no name warning, run finishes
  renamed  two atom names changed in each chain: warning naming them, run
           finishes (2 of 627 atoms, under the 5 % limit)
  rotated  the atom order inside every residue rotated by one: refused with
           "not in the atom order of the topology", non-zero exit (before
           ENE-12 every atom silently took its neighbour's charge)

  python test/topology_name_check_oracle.py --binary titan.exe --compiler clang++

Exit code 0 when every check passes, 1 otherwise.
"""
import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))

from topology_oracles.gromacs_reference import variant_dir  # noqa: E402
from topology_oracles.titan_run import run  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
FF = 'amber99sb-ildn'


def unit(compiler, work):
    exe = work / ('name_match' + ('.exe' if os.name == 'nt' else ''))
    build = subprocess.run([compiler, '-std=c++23', '-O1', '-Isrc', 'test/topology_name_match_oracle.cpp',
                            'src/topology/gmx_topology_match.cpp', '-o', str(exe)],
                           cwd=ROOT, capture_output=True, text=True)
    if build.returncode:
        return False, build.stdout + build.stderr
    result = subprocess.run([str(exe)], capture_output=True, text=True)
    return result.returncode == 0, result.stdout + result.stderr


def edited_pdb(src, dest, mode):
    lines = src.read_text(encoding='utf-8').splitlines()
    atom_rows = [k for k, l in enumerate(lines) if l.startswith(('ATOM', 'HETATM'))]
    if mode == 'renamed':
        seen = {}
        for k in atom_rows:
            chain, name = lines[k][21], lines[k][12:16].strip()
            if name == 'CB' and seen.get(chain, 0) < 2:   # first two CB of each chain
                lines[k] = lines[k][:12] + ' CX ' + lines[k][16:]
                seen[chain] = seen.get(chain, 0) + 1
    elif mode == 'rotated':
        groups, current = [], []
        for k in atom_rows:
            key = lines[k][21:27]
            if current and lines[current[0]][21:27] != key:
                groups.append(current)
                current = []
            current.append(k)
        groups.append(current)
        new = list(lines)
        for g in groups:
            names = [lines[k][12:16] for k in g]
            names = names[1:] + names[:1]
            for k, n in zip(g, names):
                new[k] = lines[k][:12] + n + lines[k][16:]
        lines = new
    dest.write_text('\n'.join(lines) + '\n', encoding='utf-8', newline='\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--binary', default='titan.exe' if os.name == 'nt' else 'titan')
    ap.add_argument('--compiler', default='clang++' if os.name == 'nt' else 'g++')
    ap.add_argument('--work-root', default=str(ROOT / 'test' / 'tmp' / 'topology_name_check'))
    a = ap.parse_args()
    binary = Path(a.binary).resolve() if Path(a.binary).exists() else Path(a.binary)
    work = Path(a.work_root).resolve()
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    failures = []

    ok, text = unit(a.compiler, work)
    print(text.strip())
    if not ok:
        failures.append('unit oracle failed')

    vdir = variant_dir(FF)
    for mode in ('same', 'renamed', 'rotated'):
        pdb = work / f'conf_{mode}.pdb'
        edited_pdb(vdir / 'conf.pdb', pdb, mode)
        # -gbsa: the plain -common_topo energy route binds titan_nb, which
        # already refuses any name mismatch; GBSA (like -nocal analyses)
        # relies on the positional assignment alone.
        args = [f'-ff={FF}', f'-f={pdb.as_posix()}', f'-n={(vdir / "index.ndx").as_posix()}',
                '-g=1-2', '-eu=kJ/mol', '-noprint-side_chain_contact=true', '-nt=1',
                '-gbsa=true', '-igb=5', f'-common_topo={(vdir / "topol.top").as_posix()}']
        finished, stdout = run(binary, work / mode, args, audit=False)
        out = stdout + (work / mode / 'stderr.log').read_text(errors='replace')
        warned = 'names differ' in out
        refused = 'not in the atom order of the topology' in out
        print(f'  {mode:8s} finished={finished} warned={warned} refused={refused}')
        expected = {'same': (True, False, False), 'renamed': (True, True, False),
                    'rotated': (False, True, True)}[mode]
        if (finished, warned, refused) != expected:
            failures.append(f'{mode}: finished/warned/refused = {(finished, warned, refused)}, '
                            f'expected {expected}')
        if mode == 'renamed' and 'atom names at 2 ' not in out:
            failures.append('renamed: the warning does not report 2 mismatched atoms per chain')

    if failures:
        print('FAIL')
        for f in failures:
            print('  -', f)
        return 1
    print('PASS: topology/structure names are checked before positional assignment')
    return 0


if __name__ == '__main__':
    sys.exit(main())
