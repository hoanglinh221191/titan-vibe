"""Verify a retained node-patch spherical-dual diagnostic campaign.

Run the C++ companion as:
  DRIVER points offset salt direction far_order near_order [solve]
Run DRIVER without arguments for kernel/trace controls. `solve` assembles
four dense experimental matrices; it is not a production solver replacement.
This verifier reads evidence and never retries or overwrites simulations.
"""
import argparse
import json
import math
from pathlib import Path

from sphere_bie_step1 import screened_oracle


def verify(root):
    rows = json.loads((root / 'rows.json').read_text())
    checks = []

    def check(name, passed, value=None):
        checks.append(dict(name=name, passed=bool(passed), value=value))

    check('matrix_36_rows', len(rows) == 36, len(rows))
    hi = [x for x in rows if x['far_order'] == 8]
    lo = {x['name'].split('_q')[0]: x for x in rows if x['far_order'] == 4}
    for x in hi:
        name = x['name']
        reference, _ = screened_oracle(x['offset'], salt=x['salt'])
        check(name + '/decimal_energy_oracle', abs(x['energy_exact'] - reference) < 1e-9)
        check(name + '/partition_area', abs(x['area_relative_error']) < 1e-11)
        check(name + '/boundary_constant_mode', x['winding_max'] < 1e-10)
        for eq in (1, 2):
            value = x[f'curved_exact_eq{eq}_relative_rms']
            check(name + f'/exact_trace_eq{eq}', value < 1e-10, value)
        value = abs(x['projection_curved_exact'] / reference - 1)
        check(name + '/exact_trace_energy_projection', value < 1e-10, value)
        value = abs((x['projection_curved_constant'] - lo[name.split('_q')[0]]['projection_curved_constant']) / reference)
        check(name + '/projection_quadrature_convergence', value < 1e-6, value)
        check(name + '/actual_unknown_count', x['scalar_unknowns'] == x['unknown_nodes'] * (2 if x['salt'] > 0 else 1))

    prototypes = json.loads((root / 'prototype_rows.json').read_text())
    extended = json.loads((root / 'extended_rows.json').read_text())
    check('prototype_10_cases', len(prototypes) == 10, len(prototypes))
    check('extended_5_cases', len(extended) == 5, len(extended))
    for x in prototypes + extended:
        check(x['name'] + '/four_ablations', len(x['solves']) == 4)
        for s in x['solves']:
            check(x['name'] + f"/mode{s['mode']}_dense_residual", math.isfinite(s['residual']) and s['residual'] < 1e-10, s['residual'])
        if '_o0_' in x['name']:
            value = abs(x['solves'][3]['energy_curved_projection'] / x['energy_exact'] - 1)
            check(x['name'] + '/center_full_curved', value < 1e-8, value)

    # Paired quadrature refinement of the full matrix, not just its energy map.
    for points in (96, 384):
        pair = [x for x in prototypes if x['name'].startswith(f'solve_p{points}_o0.8_s0.15_')]
        check(f'p{points}/quadrature_pair_exists', len(pair) == 2)
        if len(pair) == 2:
            value = abs((pair[0]['solves'][3]['energy_curved_projection'] - pair[1]['solves'][3]['energy_curved_projection']) / pair[0]['energy_exact'])
            check(f'p{points}/full_solve_quadrature_convergence', value < 1e-6, value)
    report = dict(checks=checks, passed=sum(x['passed'] for x in checks), total=len(checks),
                  production_modified=False, molecular_validation=False,
                  residual_scope='16-17 sampled rows; full dense solves for ablations',
                  note='No universal 1% gate: the coarse off-axis constant-trace discretization is inaccurate.')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('campaign', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    report = verify(args.campaign)
    if args.output:
        if args.output.exists():
            raise SystemExit('Refusing to overwrite an existing validation report')
        args.output.write_text(json.dumps(report, indent=2))
    print(f"{report['passed']}/{report['total']} checks passed")
    for x in report['checks']:
        if not x['passed']:
            print('FAIL', x)
    raise SystemExit(0 if report['passed'] == report['total'] else 1)
