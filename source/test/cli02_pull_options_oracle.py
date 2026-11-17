"""CLI-02 oracle: numeric values of the pull/DE/pathway options.

The block used bare std::stoi/stoul/stoull: "-5" for a count became
2^64 - 5, "12x" was read as 12, and a seed above the unsigned-int range was
truncated, all without a message. Each bad value must now stop the run with
an error that names the option; well-formed values must not trigger it.

  python test/cli02_pull_options_oracle.py --binary titan.exe --work-root <dir>
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument('--binary', type=Path, default=ROOT / 'titan.exe')
p.add_argument('--work-root', type=Path, default=ROOT / 'test/tmp/cli02_pull_options_oracle')
a = p.parse_args()
binary = a.binary.resolve()
work = a.work_root.resolve()
work.mkdir(parents=True, exist_ok=True)
(work / 'tiny.pdb').write_text(
    'ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n'
    'ATOM      2  CA  ALA B   2       5.000   0.000   0.000  1.00  0.00           C\nEND\n',
    encoding='ascii')

BAD = [('-de-memory-size', '-5'), ('-de-min-pop-num', '-3'), ('-pull-max-evals', '-1'),
       ('-de-pop-num', '12x'), ('-de-ite-num', '1e3'), ('-de-stagnation-limit', 'ten'),
       ('-de-pathway-group', '-1'), ('-pull-seed', '99999999999'),
       ('-pathway-corridor-radial-samples', '-16')]
GOOD = [('-de-memory-size', '7'), ('-de-pop-num', '12'), ('-pull-seed', '4294967295')]

failures = 0
for flag, value in BAD + GOOD:
    bad = (flag, value) in BAD
    run = subprocess.run([str(binary), '-f', str(work / 'tiny.pdb'), '-1pdb', flag, value],
                         cwd=work, capture_output=True, text=True, encoding='utf-8',
                         errors='replace', timeout=60, env=dict(os.environ, TITAN_NO_PAUSE='1'))
    out = run.stdout + run.stderr
    named = flag in out and ('must be a non-negative integer' in out or 'must be an integer' in out)
    ok = (run.returncode != 0 and named) if bad else not named
    print(f'{"PASS" if ok else "FAIL"} {flag} {value}: exit {run.returncode}, flagged {named}')
    failures += 0 if ok else 1
print(f'{"FAIL" if failures else "PASS"} CLI-02 pull/DE/pathway numeric options')
sys.exit(1 if failures else 0)
