#!/usr/bin/env python3
"""TITAN and APBS timed side by side (RSOS revision, cost comparisons).

This machine drifts: the same TITAN binary on the same protein took 7.9 s and
18 s an hour apart (clock 3 GHz against a 5 GHz boost), so a TITAN/APBS
ratio is only sound when both programs run under one machine state.  Each
case here runs its TITAN settings and its APBS grids back to back, and every
repeat reverses the order, so drift falls on both sides alike; the spread
between repeats is the A/A check.

Sets:
  proteins27  the 27 CHARMM proteins at 0.15 M: TITAN protocol and the
              balanced preset q7 against APBS (paper mode) at 0.5 A
  binding     6LZG and 1BRS cycles: TITAN against the three APBS runs at 0.5 A
  lowsalt5    the 5 sweep proteins at 1e-5 M: TITAN protocol at h = 1.25 and
              1.0 A against APBS (lowsalt mode) at 0.5, 0.3, 0.25 and 0.2 A,
              for the equal-accuracy comparison

    python test/tabi_rsos_interleaved_timing.py --titan bin/titan.exe \
        --work-root D:/data_MD/tabi/rsos_v06/interleaved --sets proteins27,binding,lowsalt5
"""

import argparse
import csv
import statistics
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tabi_rsos_apbs_sweep as apbs  # noqa: E402
import tabi_rsos_campaign as camp  # noqa: E402
import tabi_rsos_runner as runner  # noqa: E402

Q7 = camp.PRESETS['q7_alpha05_remesh3_theta060']
LOWSALT = ['-pbe-salt=1e-5']


def proteins27():
    for name in camp.FAGBI27_IDS:
        pqr = camp.charmm(name)
        yield name, [('titan', 'protocol', pqr, 1.25, camp.ANALYTIC),
                     ('titan', 'q7', pqr, 1.25, Q7),
                     ('apbs', 'paper_g0.5', [pqr], 0.5, 'paper')]


def binding():
    lzg = (camp.SIXLZG / '6lzg_complex_charmm_ph7.pqr',
           camp.SIXLZG / '6lzg_ace2_chain_A_charmm_ph7.pqr',
           camp.SIXLZG / '6lzg_rbd_chain_B_charmm_ph7.pqr')
    bb = (camp.BINDING_JPCB / '1brs_complex_amber_ph7.pqr',
          camp.BINDING_JPCB / '1brs_barnase_chain_A_amber_ph7.pqr',
          camp.BINDING_JPCB / '1brs_barstar_chain_D_amber_ph7.pqr')
    for name, files in (('6lzg', lzg), ('1brs', bb)):
        yield name, [('titan', 'protocol', files[0], 1.25,
                      camp.ANALYTIC + camp.binding_args(*files)),
                     ('apbs', 'paper_g0.5', list(files), 0.5, 'paper')]


def lowsalt5():
    for name in camp.HSWEEP_IDS:
        pqr = camp.charmm(name)
        entries = [('titan', f'protocol_h{h:g}', pqr, h, camp.ANALYTIC + LOWSALT)
                   for h in (1.25, 1.0)]
        entries += [('apbs', f'lowsalt_g{g:g}', [pqr], g, 'lowsalt')
                    for g in (0.5, 0.3, 0.25, 0.2)]
        yield name, entries


SETS = {'proteins27': proteins27, 'binding': binding, 'lowsalt5': lowsalt5}


def run_entry(args, root, entry):
    program, setting, inputs, level, extra = entry
    if program == 'titan':
        r = runner.run_case(args.titan, root, inputs, level, extra, {}, args.threads,
                            args.timeout)
        cycle = any(a.startswith('-pbe-binding-cycle=') for a in extra)
        energy = r['binding_dgpol_kcalmol'] if cycle else r['energy_kcalmol']
        return {'wall_s': r['wall_s'], 'peak_rss_mb': r['peak_rss_mb'],
                'energy_kcalmol': energy, 'status': r['status']}
    walls, peaks, energies, status = [], [], [], 'ok'
    for pqr in inputs:
        r = apbs.run_case(args.apbs, root / Path(pqr).stem, pqr, level, extra,
                          args.threads, args.timeout)
        walls.append(r['wall_s'])
        peaks.append(r['peak_rss_mb'])
        energies.append(r['energy_kcalmol'])
        status = status if r['exit'] == 0 else f"exit_{r['exit']}"
    # A binding cycle is the complex minus its two components.
    energy = energies[0] - sum(energies[1:]) if len(energies) == 3 else energies[0]
    return {'wall_s': sum(walls), 'peak_rss_mb': max(peaks), 'energy_kcalmol': energy,
            'status': status}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--titan', type=Path, required=True)
    ap.add_argument('--apbs', type=Path, default=apbs.DEFAULT_APBS)
    ap.add_argument('--work-root', type=Path, required=True)
    ap.add_argument('--sets', default=','.join(SETS))
    ap.add_argument('--repeats', type=int, default=2)
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--timeout', type=float, default=7200)
    args = ap.parse_args()
    args.work_root.mkdir(parents=True, exist_ok=True)
    out = args.work_root / 'interleaved.csv'
    fields = ['set', 'case', 'repeat', 'position', 'program', 'setting', 'wall_s',
              'peak_rss_mb', 'energy_kcalmol', 'status']
    rows = []
    for set_name in args.sets.split(','):
        for case, entries in SETS[set_name]():
            for repeat in range(args.repeats):
                order = entries if repeat % 2 == 0 else list(reversed(entries))
                for position, entry in enumerate(order):
                    program, setting = entry[0], entry[1]
                    root = args.work_root / set_name / f'{case}_{setting}_r{repeat}'
                    row = {'set': set_name, 'case': case, 'repeat': repeat,
                           'position': position, 'program': program, 'setting': setting}
                    row.update(run_entry(args, root, entry))
                    rows.append(row)
                    print(f"{set_name:10s} {case:5s} r{repeat} {program:5s} {setting:16s} "
                          f"{row['wall_s']:8.1f} s {row['peak_rss_mb']:7.0f} MB "
                          f"E={row['energy_kcalmol']:.4f} {row['status']}", flush=True)
                    with open(out, 'w', newline='') as f:
                        w = csv.DictWriter(f, fieldnames=fields)
                        w.writeheader()
                        w.writerows(rows)
    # Ratios within one repeat of one case, then medians over cases.
    by = defaultdict(dict)
    for r in rows:
        by[(r['set'], r['case'], r['repeat'])][r['setting']] = r
    for set_name in args.sets.split(','):
        keys = [k for k in by if k[0] == set_name]
        settings = sorted({s for k in keys for s in by[k]})
        titan = [s for s in settings if not s.startswith(('paper', 'lowsalt'))]
        grids = [s for s in settings if s.startswith(('paper', 'lowsalt'))]
        for t in titan:
            for g in grids:
                wall = [by[k][t]['wall_s'] / by[k][g]['wall_s'] for k in keys
                        if t in by[k] and g in by[k]]
                rss = [by[k][t]['peak_rss_mb'] / by[k][g]['peak_rss_mb'] for k in keys
                       if t in by[k] and g in by[k]]
                if wall:
                    print(f'{set_name}: TITAN {t} / APBS {g}: wall median {statistics.median(wall):.2f}'
                          f' (faster in {sum(w < 1 for w in wall)}/{len(wall)}), RSS median'
                          f' {statistics.median(rss):.2f}')
        # A/A: the same setting in two repeats.
        for s in settings:
            pairs = defaultdict(list)
            for k in keys:
                if s in by[k]:
                    pairs[k[1]].append(by[k][s]['wall_s'])
            spread = [max(v) / min(v) for v in pairs.values() if len(v) > 1 and min(v) > 0]
            if spread:
                print(f'{set_name}: {s} repeat spread median {statistics.median(spread):.2f}'
                      f' max {max(spread):.2f}')
    print(f'wrote {out}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
