#!/usr/bin/env python3
"""Node-patch against the screened Kirkwood sphere, through the real CLI.

One neutral sphere of radius R carries a unit charge at a fraction t of R
from its centre (a single ion when t = 0).  For one sphere the SES is the
sphere itself, so the BIE model and the analytic linearised-PB model are the
same boundary-value problem and every difference is discretisation error.

Cases run through tabi_rsos_runner.run_case, i.e. the paper protocol plus a
named variant, with a tighter GMRES tolerance so the solve does not limit the
measured error.  One CSV row per (variant, offset, h) carries the relative
error and the observed order between successive spacings.

    python test/tabi_node_patch_kirkwood.py --variants paper_np_flat,paper_np_curved
"""

import argparse
import decimal
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tabi_rsos_runner as runner  # noqa: E402

# The solve must not limit the measured discretisation error.
TIGHT = ['-pbe-rms=1e-6', '-pbe-dphi=1e-6', '-pbe-bie-matvec-audit=off']


def screened_kirkwood(t, radius, charge=1.0, epsin=1.0, epsout=80.0,
                      salt=0.15, temperature=298.15):
    """Linearised-PB Kirkwood energy (kcal/mol) of one charge in a sphere.

    Same series as test/sphere_bie_step1.py: exterior response
    gamma_l = l - z r_l with r_l = k_{l+1}(z)/k_l(z), r_{l+1} = (2l+3)/z + 1/r_l.
    kappa follows TITAN's debye_kappa_inv_a, including the eps_out scaling.
    """
    kappa = (0.3291 * math.sqrt(salt) * math.sqrt(298.15 / temperature)
             * math.sqrt(78.5 / epsout)) if salt > 0 else 0.0
    with decimal.localcontext() as ctx:
        ctx.prec = 70
        D = decimal.Decimal
        T, R, q, ei, eo = (D(repr(v)) for v in (t, radius, charge, epsin, epsout))
        total, power, ratio = D(0), D(1), T * T
        if kappa == 0.0:
            for l in range(10000):
                total += (ei - eo) * (l + 1) * power / (ei * l + eo * (l + 1))
                power *= ratio
                if ratio == 0 or power / (1 - ratio) < D('1e-40'):
                    break
        else:
            z = D(repr(kappa)) * R
            r = 1 + 1 / z
            for l in range(10000):
                gamma = D(l) - z * r
                total += (eo * gamma + ei * (l + 1)) * power / (ei * D(l) - eo * gamma)
                power *= ratio
                r = (2 * D(l) + 3) / z + 1 / r
                if ratio == 0 or power / (1 - ratio) < D('1e-40'):
                    break
        return float(D('0.5') * D('332.063713299') * q * q / (ei * R) * total)


def write_pqr(path, radius, offset_fraction, charge=1.0):
    lines = []
    if offset_fraction == 0.0:
        lines.append('ATOM      1  ION ION     1       0.000   0.000   0.000'
                     f'{charge:8.4f}{radius:8.4f}       X')
    else:
        # A neutral shell sphere and a small buried charge carrier.  The
        # carrier lies wholly inside the shell, so the SES is the shell.
        x = offset_fraction * radius
        if x + 0.2 >= radius:
            raise ValueError('charge carrier would reach the sphere surface')
        lines.append('ATOM      1  SHL SPH     1       0.000   0.000   0.000'
                     f'{0.0:8.4f}{radius:8.4f}       X')
        lines.append(f'ATOM      2  Q1  ION     2    {x:8.3f}{0.0:8.3f}{0.0:8.3f}'
                     f'{charge:8.4f}{0.2:8.4f}       X')
    lines.append('END')
    path.write_text('\n'.join(lines) + '\n', encoding='ascii')


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--titan', type=Path, default=runner.ROOT / 'titan.exe')
    parser.add_argument('--work-root', type=Path,
                        default=runner.ROOT / 'test/tmp/tabi_node_patch_kirkwood')
    parser.add_argument('--radius', type=float, default=2.0)
    parser.add_argument('--offsets', default='0,0.5,0.8')
    parser.add_argument('--h', default='1.0,0.75,0.5,0.35,0.25')
    parser.add_argument('--variants', default='paper_np_flat,paper_centroid')
    parser.add_argument('--salt', type=float, default=0.15)
    parser.add_argument('--nt', type=int, default=4)
    parser.add_argument('--timeout', type=float, default=1800)
    args = parser.parse_args()

    for name in args.variants.split(','):
        if name not in runner.VARIANTS:
            sys.exit(f'unknown variant {name}; known: {", ".join(runner.VARIANTS)}')
    args.work_root.mkdir(parents=True, exist_ok=True)
    out = args.work_root / 'kirkwood.csv'
    fields = ['offset_fraction', 'exact_kcalmol', 'rel_error',
              'observed_order'] + runner.FIELDS
    rows = []
    for name in args.variants.split(','):
        extra, env = runner.VARIANTS[name]
        for t in (float(v) for v in args.offsets.split(',')):
            exact = screened_kirkwood(t, args.radius, salt=args.salt)
            pqr = args.work_root / f'kirkwood_R{args.radius:g}_t{t:g}.pqr'
            write_pqr(pqr, args.radius, t)
            previous = None
            for h in (float(v) for v in args.h.split(',')):
                tag = f'{name}_t{t:g}_h{h:g}'.replace('.', 'p')
                row = {'case': pqr.stem, 'variant': name, 'h_A': h,
                       'offset_fraction': t, 'exact_kcalmol': exact}
                # The variant follows TIGHT so that its own options (the
                # guarded p4 settings) take precedence.
                row.update(runner.run_case(args.titan, args.work_root / tag, pqr, h,
                                           TIGHT + list(extra)
                                           + [f'-pbe-salt={args.salt}'],
                                           env, args.nt,
                                           args.timeout))
                energy = row['energy_kcalmol']
                rel = (energy - exact) / abs(exact) if math.isfinite(energy) else math.nan
                order = ''
                if previous and math.isfinite(rel) and math.isfinite(previous[1]) \
                        and rel != 0 and previous[1] != 0:
                    order = f'{math.log(abs(previous[1]) / abs(rel)) / math.log(previous[0] / h):.3f}'
                previous = (h, rel)
                row['rel_error'] = rel
                row['observed_order'] = order
                rows.append(row)
                print(f"{name:16s} t={t:<4g} h={h:<5g} n={row['panels']:<7d} "
                      f"it={row['iterations']:<4d} E={energy:12.5f} exact={exact:12.5f} "
                      f"rel={rel:+.3e} p={order} [{row['status']}|{row['exit']}] "
                      f"{row['wall_s']:.1f}s {row['node_patch_representation']}",
                      flush=True)
                runner.write_rows(out, rows, fields)
    print(f'wrote {out}')


if __name__ == '__main__':
    main()
