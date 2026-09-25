"""Optional independent APBS grid/phase sweep on exactly the supplied balls.

Use srfm mol and srad 0, not an SES/spline boundary. Solvent and vacuum
calculations share each grid so the singular grid self energy cancels.
This records grid sensitivity; a finite-grid result is not exact ground truth.
"""
import argparse
import hashlib
import json
import os
import pathlib
import re
import subprocess
import time

p=argparse.ArgumentParser()
p.add_argument('--apbs',type=pathlib.Path,required=True)
p.add_argument('--atoms',type=pathlib.Path,required=True)
p.add_argument('--output',type=pathlib.Path,required=True)
a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
atoms=[list(map(float,line.split())) for line in a.atoms.read_text().splitlines()[1:] if line.strip()]
lo=[min(b[d]-b[3] for b in atoms) for d in range(3)]
hi=[max(b[d]+b[3] for b in atoms) for d in range(3)]
assert max(hi[d]-lo[d] for d in range(3))<24
center=[(x+y)/2 for x,y in zip(lo,hi)]
env=dict(os.environ);env['PATH']=str(a.apbs.resolve().parent)+os.pathsep+env.get('PATH','')
env['OMP_NUM_THREADS']='2';env['OPENBLAS_NUM_THREADS']='2'
meta=dict(apbs=str(a.apbs.resolve()),apbs_sha256=hashlib.sha256(a.apbs.read_bytes()).hexdigest(),
          atoms_sha256=hashlib.sha256(a.atoms.read_bytes()).hexdigest(),
          surface='mol',srad=0,epsin=1,epsout=78.5,salt=0,fglen=32,cglen=64,
          omp_threads=2,energy_conversion_kj_to_kcal=4.184)
(a.output/'protocol.json').write_text(json.dumps(meta,indent=2))
rows=[]
for n,phase in [(65,0),(129,0),(129,.5),(257,0),(257,.5)]:
    name=f'n{n}_phase{phase}';out=a.output/name;out.mkdir()
    pqr=''.join(f'ATOM {i:6d} C LIG 1 '+ ' '.join(format(v,'.17g') for v in [b[0],b[1],b[2],b[4],b[3]])+'\n' for i,b in enumerate(atoms,1))
    (out/'molecule.pqr').write_text(pqr)
    h=32/(n-1);fine=[c+phase*h for c in center]
    def block(name,eps):
        return f'''elec name {name}
  mg-auto
  dime {n} {n} {n}
  cglen 64 64 64
  fglen 32 32 32
  cgcent {' '.join(map(str,center))}
  fgcent {' '.join(map(str,fine))}
  mol 1
  lpbe
  bcfl sdh
  pdie 1
  sdie {eps}
  srfm mol
  srad 0
  sdens 10
  swin 0.3
  chgm spl2
  temp 298.15
  calcenergy total
  calcforce no
end
'''
    text='read\n  mol pqr molecule.pqr\nend\n'+block('solvent',78.5)+block('vacuum',1)+'print elecEnergy solvent - vacuum end\nquit\n'
    (out/'apbs.in').write_text(text)
    cmd=[str(a.apbs.resolve()),'apbs.in']
    (out/'command.json').write_text(json.dumps(cmd,indent=2))
    row=dict(n=n,spacing=h,phase=phase);start=time.perf_counter()
    print('START',name,flush=True)
    with (out/'run.log').open('wb') as log:
        try:row['returncode']=subprocess.run(cmd,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,
            timeout=300,creationflags=0x08000000 if os.name=='nt' else 0).returncode
        except subprocess.TimeoutExpired:row['timeout']=True
    matches=re.findall(r'Global net ELEC energy\s*=\s*([-+0-9.eE]+)',(out/'run.log').read_text(errors='replace'))
    if matches:row['energy_kcal_mol']=float(matches[-1])/4.184
    row['seconds']=time.perf_counter()-start;rows.append(row)
    (a.output/'rows.json').write_text(json.dumps(rows,indent=2));print('DONE',json.dumps(row),flush=True)
raise SystemExit(0 if all(r.get('returncode')==0 and 'energy_kcal_mol' in r for r in rows) else 1)
