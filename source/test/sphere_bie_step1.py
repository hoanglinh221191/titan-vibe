"""Bounded sphere BIE diagnosis (no production changes).

python test/sphere_bie_step1.py --driver DRIVER_EXE --output NEW_DIRECTORY
Requires the native driver titan_sphere_bie_step1.cpp linked with matching
engine architecture flags. Raw evidence is retained; no automatic retries.
"""
import argparse
import csv
import decimal
import hashlib
import json
import math
import os
import pathlib
import subprocess
import time


def oracle(offset, radius=1.7, charge=1, epsin=1, epsout=78.5):
    """Independent high precision scalar Kirkwood sum, with a tail bound.

    Coefficient absolute value is <= 1 for 0 < epsin <= epsout.
    At the charge location the geometric ratio is offset**2. The omitted
    dimensionless tail is <= ratio**(l+1)/(1-ratio).
    """
    with decimal.localcontext() as ctx:
        ctx.prec = 70
        D = decimal.Decimal
        t, R, q, ei, eo = map(lambda v: D(str(v)), (offset, radius, charge, epsin, epsout))
        total, power, ratio = D(0), D(1), t*t
        for l in range(10000):
            total += (ei-eo)*(l+1)*power/(ei*l+eo*(l+1))
            power *= ratio
            tail = power/(1-ratio)
            if tail < D('1e-40'):
                scale = D('0.5')*D('332.063713299')*q*q/(ei*R)
                return float(scale*total), float(abs(scale)*tail)
        raise RuntimeError('Kirkwood oracle failed to converge')


def screened_oracle(offset, radius=1.7, charge=1, epsin=1, epsout=78.5, salt=0.0,
                    temperature=298.15):
    """Linearised-PB Kirkwood sum, independent of the driver's own copy.

    Matching potential and dielectric flux on the sphere gives the exterior
    response gamma_l = kappa*a*k_l'(kappa*a)/k_l(kappa*a) = l - z*r_l, with
    r_l = k_{l+1}(z)/k_l(z) obeying r_{l+1} = (2l+3)/z + 1/r_l, r_0 = 1+1/z.
    As z -> 0, gamma_l -> -(l+1) and this reduces to oracle() term by term.
    """
    if salt <= 0:
        return oracle(offset, radius, charge, epsin, epsout)
    kappa = 0.3291*math.sqrt(salt)*math.sqrt(298.15/temperature)*math.sqrt(78.5/epsout)
    with decimal.localcontext() as ctx:
        ctx.prec = 70
        D = decimal.Decimal
        z = D(str(kappa))*D(str(radius))
        t, R, q, ei, eo = map(lambda v: D(str(v)), (offset, radius, charge, epsin, epsout))
        total, power, ratio, r = D(0), D(1), t*t, 1+1/z
        for l in range(10000):
            gamma = D(l)-z*r
            total += (eo*gamma+ei*(l+1))*power/(ei*D(l)-eo*gamma)
            power *= ratio
            r = (2*D(l)+3)/z+1/r
            if ratio == 0 or power/(1-ratio) < D('1e-40'):
                scale = D('0.5')*D('332.063713299')*q*q/(ei*R)
                return float(scale*total), float(abs(scale)*power/(1-ratio) if ratio else 0.0)
        raise RuntimeError('screened Kirkwood oracle failed to converge')


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--driver',required=True,type=pathlib.Path)
    p.add_argument('--output',required=True,type=pathlib.Path)
    p.add_argument('--baseline',type=pathlib.Path,
                   help='Pre-fix rows.json; enables corrected Poisson regression gates')
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    old_rows=json.loads(a.baseline.read_text()) if a.baseline else []
    # Declared before executing the matrix. Accuracy is reported separately
    # from reproducibility/oracle/convergence checks: a small GMRES residual
    # does not make a coarse mesh pass the 1% physical-accuracy gate.
    thresholds=dict(solver_residual=2e-10,oracle_relative=2e-13,
                    geometry_area_relative=1e-11,radial_relative=1e-12,
                    centered_projection_relative=1e-11,
                    physical_energy_relative=0.01,
                    self_quadrature_absolute=2e-5,
                    curved_constant_identity_absolute=2e-5,
                    timeout_seconds=90)
    (a.output/'thresholds.json').write_text(json.dumps(thresholds,indent=2))
    rows=[];checks=[]
    def check(name,passed): checks.append(dict(name=name,passed=bool(passed)))

    def run(mesh,points,offset,near='analytic_singular',diagonal='solid-angle',
            tol=1e-10,radius=1.7,charge=1,epsin=1,epsout=78.5,shift=0,direction=0,tag='matrix',
            disc='centroid',salt=0.0):
        name=f'{len(rows):02d}_{tag}_{mesh}_{points}_{offset}'
        args=[mesh,points,offset,near,diagonal,tol,radius,charge,epsin,epsout,shift,direction]
        env=dict(os.environ,TITAN_STEP1_DISCRETIZATION=disc,TITAN_STEP1_SALT=repr(salt))
        (a.output/(name+'.command.json')).write_text(json.dumps(
            dict(argv=[str(a.driver.resolve()),*map(str,args)],discretization=disc,salt_m=salt),indent=2))
        begin=time.perf_counter()
        try:
            proc=subprocess.run([str(a.driver.resolve()),*map(str,args)],capture_output=True,text=True,
                                timeout=thresholds['timeout_seconds'],env=env)
        except subprocess.TimeoutExpired as e:
            (a.output/(name+'.log')).write_bytes((e.stdout or b'')+(e.stderr or b''))
            raise
        (a.output/(name+'.log')).write_text(proc.stdout+'\n'+proc.stderr)
        if proc.returncode:
            raise RuntimeError(f'{name} failed: {proc.stderr}')
        result=json.loads(next(x[7:] for x in proc.stdout.splitlines() if x.startswith('RESULT ')))
        ref,tail=screened_oracle(offset,radius,charge,epsin,epsout,salt)
        row=dict(name=name,tag=tag,mesh=mesh,points=points,offset=offset,near=near,
                 diagonal=diagonal,tol=tol,radius=radius,charge=charge,epsin=epsin,epsout=epsout,
                 shift=shift,direction=direction,disc=disc,salt=salt,
                 oracle=ref,oracle_tail_bound=tail,
                 wall_seconds=time.perf_counter()-begin,**result)
        scale=max(abs(ref),1.)
        row['energy_relative_error']=abs(row['energy']-ref)/scale
        row['projection_relative_error']=abs(row['exact_trace_projection_energy']-ref)/scale
        row['physical_1_percent_pass']=row['energy_relative_error']<=thresholds['physical_energy_relative']
        rows.append(row)
        (a.output/'rows.json').write_text(json.dumps(rows,indent=2))
        check(name+': independent oracle',abs(row['reference']-ref)/scale<thresholds['oracle_relative'])
        check(name+': solver converged',row['residual']<max(2*tol,thresholds['solver_residual']))
        check(name+': native energy',row['energy_source']=='tabi_boundary_integral')
        check(name+': spherical boundary',abs(row['area']/(4*math.pi*radius**2)-1)<thresholds['geometry_area_relative'] and
              row['radial_error']/radius<thresholds['radial_relative'])
        # The operator sampler models a spherical triangle at zero salt; it is
        # switched off for node-patch cells and for screened kernels, so its
        # gates only apply where the driver reports that it ran.
        if row.get('sampled',True):
            check(name+': quadrature refinement',row['self_quadrature_64_128_max_abs']<thresholds['self_quadrature_absolute'])
            check(name+': independent constant identity',row['curved_constant_identity_max_abs']<thresholds['curved_constant_identity_absolute'])
            if offset==0:
                check(name+': exact-trace energy projection',row['projection_relative_error']<thresholds['centered_projection_relative'])
        print(f"{name}: N={row['panels']}, E={row['energy']:.9f}, error={100*row['energy_relative_error']:.4f}%, residual={row['residual']:.2g}",flush=True)
        return row

    base={}
    for mesh,points in [('icosphere',[80,320,1280,4096]),('latlong',[16,64,289,1225])]:
        for offset in [0,.5,.8]:
            series=[run(mesh,n,offset) for n in points]
            for r in series:base[(mesh,r['points'],offset)]=r
            errors=[r['energy_relative_error'] for r in series]
            if offset==0:
                check(mesh+': centered energy converges monotonically',all(x>y for x,y in zip(errors,errors[1:])))
            check(f'{mesh}/{offset}: finest improves preceding level',errors[-1]<errors[-2])
            for prev,r in zip(series,series[1:]):
                r['observed_order_h']=math.log(prev['energy_relative_error']/r['energy_relative_error'])/math.log(math.sqrt(r['panels']/prev['panels']))
    for offset in [0,.8]:
        baseline=base[('icosphere',1280,offset)]
        for near in ['centroid','duffy_adaptive']:
            r=run('icosphere',1280,offset,near=near,tag='near_control')
            reference=baseline
            if a.baseline and near=='centroid':
                reference=next(x for x in old_rows if x['tag']=='matrix' and
                    x['mesh']=='icosphere' and x['points']==1280 and x['offset']==offset)
            check(r['name']+': near mode matches its specified reference',abs(r['energy']-reference['energy'])<1e-10)
        for diagonal in ['flat','off']:
            r=run('icosphere',1280,offset,diagonal=diagonal,tag='diagonal_control')
            check(r['name']+': JCP smooth jump unaffected by metadata',abs(r['energy']-baseline['energy'])<1e-10)
    # Salted and node-patch coverage.  The matrix above only ever exercised
    # centroid at zero salt, which is why a correction gated on salt_m <= 0
    # could sit in the tree unnoticed while every ionic-strength run used the
    # uncorrected operator.
    for offset in [0,.5,.8]:
        series=[run('icosphere',n,offset,salt=.15,tag='salted') for n in [320,1280,5120]]
        errors=[r['energy_relative_error'] for r in series]
        check(f'salted/{offset}: finest improves preceding level',errors[-1]<errors[-2])
        # Coarse meshes miss 1 percent at 0.8R with or without salt, exactly as
        # the unscreened matrix already records; gate the resolved level only.
        check(series[-1]['name']+': salted 1 percent gate',series[-1]['physical_1_percent_pass'])
    # Physics invariant, and the specific regression that was missed: an
    # infinitesimal ionic strength must not change the answer, yet it flips
    # coupled = kappa > 1e-14 and therefore the whole kernel path.
    for offset in [0,.8]:
        trace=run('icosphere',1280,offset,salt=1e-14,tag='salt_limit_control')
        # kappa = 3.3e-8 1/A is a real but tiny shift, about 7e-10 relative;
        # the defect this guards against was 0.69 kcal/mol, five orders larger.
        check(trace['name']+': salt -> 0 limit matches the unscreened solve',
              abs(trace['energy']-base[('icosphere',1280,offset)]['energy'])<1e-6)
    # Node-patch is the manuscript's discretization and had no oracle coverage.
    # These rows are recorded, not gated at 1 percent: the self term still
    # carries only its leading curvature correction.
    # A centred charge in a sphere has a constant boundary trace, so offset 0
    # exercises only the l=0 mode.  Gating node-patch on it alone tested almost
    # nothing; the displaced rows below carry the real convergence evidence.
    for offset in [0,.5,.8]:
        node=[run('icosphere',n,offset,salt=.15,disc='node-patch',tag='node_patch')
              for n in [320,1280,5120]]
        errors=[r['energy_relative_error'] for r in node]
        check(f'node-patch/{offset}: finest within 1 percent',errors[-1]<0.01)
        if offset==0:
            # The constant mode is the discrete Gauss identity.  A curved
            # representation reproduces it exactly; a flat one cannot.
            check('node-patch: constant mode is exact',errors[-1]<1e-8)
        if offset==.5:
            check('node-patch: monotone under refinement',
                  all(x>y for x,y in zip(errors,errors[1:])))
            order=math.log(errors[0]/errors[-1])/math.log(math.sqrt(
                node[-1]['panels']/node[0]['panels']))
            node[-1]['observed_order_h']=order
            # Second order is the point of the curved dual cell; centroid is
            # first order on the same series.  Guard the qualitative gain.
            check('node-patch: displaced trace converges faster than first order',
                  order>1.5)
    r=run('icosphere',1280,.8,tol=1e-12,tag='tolerance_control')
    check('tight tolerance has negligible energy effect',abs(r['energy']-base[('icosphere',1280,.8)]['energy'])<1e-7)
    shifted=run('icosphere',320,.5,shift=100,tag='translation_control')
    check('translation energy invariant',abs(shifted['energy']-base[('icosphere',320,.5)]['energy'])<1e-8)
    scaled=run('icosphere',320,.5,radius=3.4,charge=2,tag='scaling_control')
    check('q squared / radius energy scaling',abs(scaled['energy']-2*base[('icosphere',320,.5)]['energy'])<1e-8)
    run('icosphere',1280,.5,epsin=4,tag='dielectric_control')
    homogeneous=run('icosphere',320,.5,epsin=78.5,tag='homogeneous_control')
    check('equal dielectric has zero reaction energy',abs(homogeneous['energy'])<1e-12)
    run('icosphere',1280,.8,direction=1,tag='rotation_sensitivity')
    centered=base[('icosphere',1280,0)]
    if a.baseline:
        old=next(x for x in old_rows if x['tag']=='matrix' and x['mesh']=='icosphere' and x['points']==1280 and x['offset']==0)
        check('corrected operator defect improves by at least 20x',
              centered['sample_operator_rel']<.05*old['sample_operator_rel'])
        check('curved self contribution is present',centered['zero_self_samples']==0 and centered['near_duffy_pairs']>0)
        for offset in [0,.8]:
            r=run('icosphere',1280,offset,near='auto',tag='auto_control')
            check(r['name']+': auto uses corrected kernel',abs(r['energy']-base[('icosphere',1280,offset)]['energy'])<1e-10)
        check('finest icosphere passes 1% at all offsets',all(base[('icosphere',4096,t)]['physical_1_percent_pass'] for t in [0,.5,.8]))
    else:
        check('self omission dominates sampled centered operator defect',
              centered['sample_self_added_rel'] < .15*centered['sample_operator_rel'])
    check('curved quadrature removes centered operator defect',centered['sample_curved_quadrature_rel']<1e-7)
    check('GMRES error cannot explain centered operator defect',centered['sample_operator_rel']>1e6*centered['residual'])
    with (a.output/'results.csv').open('w',newline='') as f:
        fields=list(dict.fromkeys(k for r in rows for k in r))
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    (a.output/'rows.json').write_text(json.dumps(rows,indent=2))
    summary=dict(diagnostic_checks=checks,diagnostic_pass=all(x['passed'] for x in checks),
                 physical_gate_relative=thresholds['physical_energy_relative'],
                 default_icosphere_320_all_offsets_pass=all(base[('icosphere',320,t)]['physical_1_percent_pass'] for t in [0,.5,.8]),
                 finest_icosphere_5120_all_offsets_pass=all(base[('icosphere',4096,t)]['physical_1_percent_pass'] for t in [0,.5,.8]),
                 production_code_changed=bool(a.baseline))
    (a.output/'summary.json').write_text(json.dumps(summary,indent=2))
    paths=[a.driver,pathlib.Path(__file__),pathlib.Path('test/titan_sphere_bie_step1.cpp')]
    paths += list(a.output.glob('*.json'))+list(a.output.glob('*.log'))+[a.output/'results.csv']
    (a.output/'sha256.json').write_text(json.dumps({str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in paths},indent=2))
    print(json.dumps({k:v for k,v in summary.items() if k!='diagnostic_checks'}),flush=True)
    return not summary['diagnostic_pass']


if __name__=='__main__':raise SystemExit(main())
