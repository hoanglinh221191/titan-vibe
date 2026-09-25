"""Validate compiled TITAN GBn2/SES-R6 endpoints against frozen Amber results.

GBn2 uses the existing 0.01 kJ/mol + 1e-4 relative gate. SES-R6 compares
different quadratures, with 0.5 kJ/mol + 0.2% energy and 0.01 A radius gates.
These are numerical model checks, not validation of binding affinities.
"""
import argparse, hashlib, json, math, subprocess
from pathlib import Path

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--reference',type=Path,required=True)
p.add_argument('--driver',type=Path,required=True)
p.add_argument('--work-root',type=Path,required=True)
a=p.parse_args();w=a.work_root.resolve();w.mkdir(parents=True,exist_ok=False)
sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()
ref=json.loads(a.reference.read_text());base=a.reference.resolve().parent
for path,digest in ref['files_sha256'].items():
 if sha(base/path)!=digest:raise ValueError('reference hash mismatch: '+path)

def run(fixture,label):
 r=subprocess.run([str(a.driver.resolve()),str(fixture.resolve())],capture_output=True,text=True)
 (w/(label+'.log')).write_text(r.stdout+r.stderr)
 r.check_returncode()
 energy=list(map(float,next(l for l in r.stdout.splitlines() if l.startswith('ENERGY ')).split()[1:]))
 radii=[list(map(float,l.split()[1:])) for l in r.stdout.splitlines() if l.startswith('RADII')]
 if len(energy)!=2 or len(radii)!=4 or not all(math.isfinite(v) and v>0 for row in radii for v in row):
  raise ValueError('invalid driver output')
 return energy,radii

rows=[]
for case in ref['records']:
 values,radii=run(base/case['fixture'],case['case'])
 checks=[]
 for k,v in zip(['elec','gb'],values):
  expected=case['delta_kj'][k]
  tol=(.5+.002*abs(expected)) if case['igb']==66 and k=='gb' else (.01+.0001*abs(expected))
  err=abs(v-expected)
  checks.append(dict(component=k,titan_kj=v,amber_kj=expected,error=err,tolerance=tol,passed=err<=tol))
 if case['igb']==66:
  expected=[case['endpoints_kj'][e]['radii_A'] for e in ['rec','lig','complex']]
  actual=[radii[0],radii[1],radii[2]+radii[3]]
  for endpoint,t,v in zip(['rec','lig','complex'],expected,actual):
   assert len(t)==len(v)
   err=max(abs(x-y) for x,y in zip(t,v))
   checks.append(dict(component=endpoint+'_radii_A',error=err,tolerance=.01,passed=err<=.01))
 rows.append(dict(case=case['case'],checks=checks))
 print(case['case'],'PASS' if all(c['passed'] for c in checks) else 'FAIL',flush=True)

# Independent exact centered-sphere radii. Check decreasing error on three meshes.
sphere=next(c for c in ref['records'] if c['case'].startswith('spheres_igb66_salt0.0'))
source=(base/sphere['fixture']).read_text().splitlines();head=source[0].split();convergence=[]
for spacing in [.5,.25,.125]:
 head[-1]=str(spacing);f=w/f'sphere_{spacing}.txt';f.write_text(' '.join(head)+'\n'+'\n'.join(source[1:])+'\n')
 _,radii=run(f,f'sphere_{spacing}')
 errors=[abs(radii[i][0]-r/(1+.028*r)) for i,r in enumerate([1.7,1.5])]
 convergence.append(dict(spacing_A=spacing,max_radius_error_A=max(errors)))
converged=all(convergence[i+1]['max_radius_error_A']<convergence[i]['max_radius_error_A'] for i in range(2)) and convergence[-1]['max_radius_error_A']<.001
ok=converged and all(c['passed'] for row in rows for c in row['checks'])
report={'status':'PASS' if ok else 'FAIL','scope':__doc__,'records':rows,
 'sphere_convergence':convergence,'sphere_convergence_passed':converged,
 'reference_sha256':sha(a.reference),'driver_sha256':sha(a.driver),'script_sha256':sha(__file__)}
(w/'comparison.json').write_text(json.dumps(report,indent=2)+'\n')
raise SystemExit(not ok)
