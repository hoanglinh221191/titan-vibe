#!/usr/bin/env python3
"""Accuracy-matched cost of TITAN BIE-PB and APBS (RSOS revision, reviewer 2).

Reads the TITAN spacing sweep (tabi_rsos_runner.py matrix.csv) and the APBS
grid sweeps (tabi_rsos_apbs_sweep.py apbs_<mode>.csv) and, per protein and
method, measures each level against that method's own converged value:

- APBS: E(g) = E0 + C g^p fitted through the three finest grids (p from
  their two differences; the finest value when they are not monotone).
- TITAN: the finest spacing, h = 0.5 A.  Its sweep moves by a few 1e-4
  between h = 1.25 and 0.5 without a monotone trend, so no extrapolation is
  attempted; the spread of h <= 0.75 is reported as its uncertainty.

Writes per-level and per-target tables (the cheapest level within 1, 0.5,
0.2% of its own limit, with wall time and peak working set) and a
work-precision figure.

    python test/tabi_rsos_workprecision.py --titan hsweep/matrix.csv \
        --apbs apbs_sweep/apbs_matched.csv --out analysis
"""

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

TARGETS = (1e-2, 5e-3, 2e-3)


def load(path, level_key):
    data = defaultdict(list)
    for row in csv.DictReader(open(path, newline='')):
        energy = float(row['energy_kcalmol'])
        if not math.isfinite(energy) or str(row.get('exit', '0')) not in ('0',):
            continue
        case = row['case'].split('_')[0]
        data[case].append({'level': float(row[level_key]), 'energy': energy,
                           'wall': float(row['wall_s']),
                           'rss': float(row['peak_rss_mb'])})
    for rows in data.values():
        rows.sort(key=lambda r: -r['level'])
    return data


def extrapolate(rows):
    """E0 and order p from the three finest levels of a sweep."""
    (g1, e1), (g2, e2), (g3, e3) = [(r['level'], r['energy']) for r in rows[-3:]]
    d12, d23 = e1 - e2, e2 - e3
    if d12 == 0 or d23 == 0 or d12 * d23 < 0:
        return e3, math.nan
    ratio = d12 / d23

    def f(p):
        return (g1 ** p - g2 ** p) / (g2 ** p - g3 ** p) - ratio

    lo, hi = 0.2, 6.0
    if f(lo) * f(hi) > 0:
        return e3, math.nan
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(lo) * f(mid) <= 0:
            hi = mid
        else:
            lo = mid
    p = 0.5 * (lo + hi)
    c = d23 / (g2 ** p - g3 ** p)
    return e3 - c * g3 ** p, p


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--titan', type=Path, required=True)
    parser.add_argument('--apbs', type=Path, nargs='+', required=True)
    parser.add_argument('--titan-max-h', type=float, default=1.25,
                        help='coarsest TITAN spacing taken as usable')
    parser.add_argument('--common-reference', action='store_true',
                        help='measure every method against the finest TITAN value; '
                             'meaningful only when both solve one boundary-value '
                             'problem (the low-salt sweeps)')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    titan = load(args.titan, 'h_A')
    methods = {}
    titan_ref = {}
    for case, rows in titan.items():
        rows = [r for r in rows if r['level'] <= args.titan_max_h + 1e-12]
        ref = rows[-1]['energy']
        titan_ref[case] = ref
        spread = max(abs(r['energy'] - ref) for r in rows if r['level'] <= 0.75)
        methods.setdefault('TITAN BIE', {})[case] = (rows, ref, math.nan, spread)
    for path in args.apbs:
        mode = Path(path).stem.replace('apbs_', '')
        for case, rows in load(path, 'spacing_A').items():
            e0, p = extrapolate(rows)
            if args.common_reference and case in titan_ref:
                e0 = titan_ref[case]
            methods.setdefault(f'APBS {mode}', {})[case] = (rows, e0, p, math.nan)

    level_rows, target_rows = [], []
    for method, cases in methods.items():
        for case, (rows, ref, p, spread) in sorted(cases.items()):
            for r in rows:
                err = abs(r['energy'] - ref) / abs(ref)
                level_rows.append({'method': method, 'case': case,
                                   'level_A': r['level'], 'energy_kcalmol': r['energy'],
                                   'reference_kcalmol': ref, 'rel_error': err,
                                   'order_p': p, 'titan_spread_kcalmol': spread,
                                   'wall_s': r['wall'], 'peak_rss_mb': r['rss']})
            for target in TARGETS:
                ok = [r for r in rows if abs(r['energy'] - ref) / abs(ref) <= target]
                best = min(ok, key=lambda r: r['wall']) if ok else None
                target_rows.append({
                    'method': method, 'case': case, 'target_rel_error': target,
                    'level_A': best['level'] if best else '',
                    'wall_s': f"{best['wall']:.2f}" if best else '',
                    'peak_rss_mb': f"{best['rss']:.0f}" if best else '',
                    'reached': bool(best)})

    with open(args.out / 'workprecision_levels.csv', 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(level_rows[0]))
        writer.writeheader()
        writer.writerows(level_rows)
    with open(args.out / 'workprecision_targets.csv', 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(target_rows[0]))
        writer.writeheader()
        writer.writerows(target_rows)

    # Converged values side by side.
    with open(args.out / 'converged_values.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        names = sorted(methods)
        writer.writerow(['case'] + [f'{m} (kcal/mol)' for m in names])
        for case in sorted(titan):
            writer.writerow([case] + [f"{methods[m][case][1]:.3f}" if case in methods[m] else ''
                                      for m in names])

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    cases = sorted(titan)
    fig, axes = plt.subplots(2, len(cases), figsize=(3.2 * len(cases), 6.2),
                             sharey='row', constrained_layout=True)
    styles = {'TITAN BIE': ('o-', '#1f77b4'), 'APBS matched': ('s-', '#d62728'),
              'APBS paper': ('^--', '#7f7f7f')}
    for col, case in enumerate(cases):
        for method, cases_m in methods.items():
            if case not in cases_m:
                continue
            rows, ref, _, _ = cases_m[case]
            pts = [(r['wall'], r['rss'], abs(r['energy'] - ref) / abs(ref))
                   for r in rows if r['energy'] != ref]
            if not pts:
                continue
            marker, color = styles.get(method, ('x-', 'k'))
            axes[0, col].loglog([p[0] for p in pts], [p[2] for p in pts], marker,
                                color=color, label=method, ms=4)
            axes[1, col].loglog([p[1] for p in pts], [p[2] for p in pts], marker,
                                color=color, label=method, ms=4)
        axes[0, col].set_title(case.upper())
        axes[0, col].set_xlabel('wall time (s, 4 threads)')
        axes[1, col].set_xlabel('peak memory (MB)')
        for row in (0, 1):
            for t in TARGETS:
                axes[row, col].axhline(t, color='0.85', lw=0.8, zorder=0)
    ylabel = ('relative error vs TITAN at h = 0.5 A' if args.common_reference
              else 'relative error vs own converged value')
    axes[0, 0].set_ylabel(ylabel)
    axes[1, 0].set_ylabel(ylabel)
    axes[0, 0].legend(fontsize=7)
    fig.savefig(args.out / 'workprecision.png', dpi=200)
    fig.savefig(args.out / 'workprecision.svg')
    print('wrote', args.out)


if __name__ == '__main__':
    main()
