"""GIST CLI routing/water ingestion parity using rounded reference frames."""
import argparse,hashlib,json,subprocess,math,random
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--driver',type=Path,required=True);p.add_argument('--reference',type=Path);p.add_argument('--work-root',type=Path,required=True);a=p.parse_args();w=a.work_root.resolve();w.mkdir(parents=True,exist_ok=False)
nf=32;nw=4
if a.reference:
 lines=(a.reference/'titan_fixture.txt').read_text().splitlines();coords=[[list(map(float,line.split())) for line in lines[1+f*12:1+(f+1)*12]] for f in range(nf)]
else:
 rng=random.Random(87512);coords=[]
 for f in range(nf):
  frame=[]
  for c in [[2.2,2.2,.2],[-2.2,2.2,.2],[2.2,-2.2,.2],[-2.2,-2.2,.2]]:
   q=[rng.gauss(0,1) for _ in range(4)];norm=math.sqrt(sum(v*v for v in q));qw,x,y,z=[v/norm for v in q];r=[[1-2*(y*y+z*z),2*(x*y-z*qw),2*(x*z+y*qw)],[2*(x*y+z*qw),1-2*(x*x+z*z),2*(y*z-x*qw)],[2*(x*z-y*qw),2*(y*z+x*qw),1-2*(x*x+y*y)]];frame.append(c)
   for h in [[.9572,0,0],[.9572*math.cos(math.radians(104.52)),.9572*math.sin(math.radians(104.52)),0]]:frame.append([c[i]+sum(r[i][j]*h[j] for j in range(3)) for i in range(3)])
  coords.append(frame)
pdb=['CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1           1'];rounded=[]
for f,water in enumerate(coords):
 pdb+=['TITLE Numerical GIST fixture t= '+str(f),f'MODEL     {f+1:4d}'];points=[[-1,-1,0],[1,-1,0],[0,1,0],[0,0,0]]+water
 for i,xyz in enumerate(points):
  if i<3:atom,res,chain,resid='CA','ALA','A',i+1
  elif i==3:atom,res,chain,resid='C','LIG','B',4
  else:atom,res,chain,resid=['OW','HW1','HW2'][(i-4)%3],'SOL','W',5+(i-4)//3
  x,y,z=[float(f'{v:.3f}') for v in xyz]
  if i>=4:rounded.append((x,y,z))
  pdb.append(f'ATOM  {i+1:5d} {atom:^4s} {res:3s} {chain}{resid:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {"H" if atom.startswith("H") else "O" if atom=="OW" else "C"}')
  if i in [2,3]:pdb.append('TER')
 pdb.append('ENDMDL')
pdb.append('END');(w/'fixture.pdb').write_text('\n'.join(pdb)+'\n')
top=['[ defaults ]','1 2 yes 0.5 0.8333333333','[ atomtypes ]','DU 6 12.0 0.0 A 0.3 0.0','OW 8 15.9994 -0.834 A 0.31507 0.636','HW 1 1.008 0.417 A 0 0','[ moleculetype ]','Fixture 0','[ atoms ]']
for i in range(16):
 if i<4:typ,resid,res,atom,q,mass='DU',i+1,'ALA' if i<3 else 'LIG','CA' if i<3 else 'C',0,12
 else:typ='OW' if (i-4)%3==0 else 'HW';resid=5+(i-4)//3;res='SOL';atom=['OW','HW1','HW2'][(i-4)%3];q=-.834 if typ=='OW' else .417;mass=15.9994 if typ=='OW' else 1.008
 top.append(f'{i+1} {typ} {resid} {res} {atom} {i+1} {q} {mass}')
top+=['[ system ]','Numerical GIST fixture','[ molecules ]','Fixture 1'];(w/'fixture.top').write_text('\n'.join(top)+'\n');(w/'groups.ndx').write_text('[ System ]\n'+' '.join(map(str,range(1,17)))+'\n[ Receptor ]\n1 2 3\n[ Ligand ]\n4\n[ Water ]\n'+' '.join(map(str,range(5,17)))+'\n')
case=w/'cli';case.mkdir();opts=[f'-f={(w/"fixture.pdb").as_posix()}','-1pdb=true','-ff=AMBER99SB',f'-n={(w/"groups.ndx").as_posix()}','-g=1-2','-mol-1-type=protein','-mol-2-type=lig',f'-path-a={w.as_posix()}/',f'-path-n={w.as_posix()}/',f'-common_topo={(w/"fixture.top").as_posix()}','-nt=2','-nocal=true','-gist=true','-gist-padding=6','-gist-resolution=0.5','-gist-temperature=300','-gist-cutoff=9','-gist-min-count=2']
opts.remove('-nocal=true') # GIST uses the standard analysis path, not I/O-only mode.
(case/'arg.in').write_text('\n'.join(opts)+'\n');r=subprocess.run([str(a.binary.resolve()),'-a',str(case/'arg.in')],cwd=case,capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=120);(case/'run.log').write_text(r.stdout+r.stderr,encoding='utf-8');r.check_returncode()
reports=list(case.glob('*_observables.json'));assert len(reports)==1,reports;status=json.loads(reports[0].read_text());assert status['processed_frames']==nf and status['orientation_samples']==nf*nw and status['E_ww_available'] and status['orientational_entropy_available'],status
fixture=w/'rounded.txt';fixture.write_text(f'{nf} {nw}\n'+'\n'.join(' '.join(map(str,xyz)) for xyz in rounded)+'\n');direct=w/'direct';direct.mkdir();r=subprocess.run([str(a.driver.resolve()),str(fixture),str(direct/'gist')],cwd=direct,capture_output=True,text=True,encoding='utf-8',errors='replace');(direct/'run.log').write_text(r.stdout+r.stderr,encoding='utf-8');r.check_returncode()
def dx(path):
 lines=path.read_text().splitlines();i=next(i for i,l in enumerate(lines) if 'data follows' in l);n=int(lines[i].split('items')[1].split()[0]);values=[]
 for line in lines[i+1:]:
  if len(values)>=n:break
  values.extend(map(float,line.split()))
 assert len(values)==n;return values
checks={}
for suffix in ['density','E_sw','E_ww_unref','minus_TdS_orient_nn']:
 actual=list(case.glob('*_'+suffix+'.dx'));assert len(actual)==1
 x=dx(actual[0]);y=dx(direct/('gist_'+suffix+'.dx'));assert len(x)==len(y);errors=[]
 for u,v in zip(x,y):
  assert math.isnan(u)==math.isnan(v)
  if math.isfinite(u):errors.append(abs(u-v))
 assert errors and max(errors)<1e-8,(suffix,max(errors));checks[suffix]={'max_error':max(errors),'finite_voxels':len(errors)}
variants=[]
def variant(name,changes,expected_error=None,index_text=None):
 folder=w/name;folder.mkdir();options=opts.copy()
 for flag,value in changes.items():options=[o for o in options if not o.startswith(flag+'=')]+[flag+'='+value]
 if index_text is not None:
  index=folder/'groups.ndx';index.write_text(index_text);options=[o for o in options if not o.startswith('-n=')]+['-n='+index.as_posix()]
 (folder/'arg.in').write_text('\n'.join(options)+'\n');r=subprocess.run([str(a.binary.resolve()),'-a',str(folder/'arg.in')],cwd=folder,capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=120);text=r.stdout+r.stderr;(folder/'run.log').write_text(text,encoding='utf-8')
 if expected_error:
  assert r.returncode!=0 and expected_error in text,(name,r.returncode,text[-1600:]);assert not list(folder.glob('*_observables.json'))
 else:r.check_returncode()
 variants.append({'name':name,'status':'PASS','expected_rejection':bool(expected_error)});return folder
explicit=variant('explicit_water',{'-n_wat':'3'});report=json.loads(next(explicit.glob('*_observables.json')).read_text());assert report['orientation_samples']==128
settings=variant('settings',{'-gist-padding':'5','-gist-resolution':'1','-gist-temperature':'315','-gist-cutoff':'7','-gist-water-model':'SPCE','-gist-min-count':'3'});report=json.loads(next(settings.glob('*_observables.json')).read_text());assert report['temperature_K']==315 and report['oxygen_cutoff_A']==7 and report['water_model']=='SPCE';assert len(dx(next(settings.glob('*_density.dx'))))==1000
minimum=variant('min_count',{'-gist-min-count':'33'});assert json.loads(next(minimum.glob('*_observables.json')).read_text())['occupied_voxels']==0
disabled=variant('disabled',{'-gist':'false'});assert not list(disabled.glob('*_observables.json'))
variant('invalid_numeric',{'-gist-padding':'5garbage'},'-gist-padding must be numeric')
variant('invalid_count',{'-gist-min-count':'2.5'},'-gist-min-count must be an integer')
variant('invalid_water_index',{'-n_wat':'3garbage'},'-n_wat must be an integer')
variant('system_as_water',{'-n_wat':'0'},'GIST water group')
variant('solute_as_water',{'-n_wat':'1'},'GIST water group overlaps solute')
variant('nocal',{'-nocal':'true'},'-gist cannot be combined with -nocal')
original_index=(w/'groups.ndx').read_text();variant('ambiguous_water',{},'GIST water group is ambiguous',original_index+'[ WAT ]\n5 6 7\n')
alias=variant('identical_water_alias',{},index_text=original_index+'[ WAT ]\n'+' '.join(map(str,range(5,17)))+'\n');assert json.loads(next(alias.glob('*_observables.json')).read_text())['orientation_samples']==128
variant('missing_water',{},'GIST requires -n_wat',original_index.replace('[ Water ]','[ Arbitrary ]'))
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
(w/'verification.json').write_text(json.dumps({'status':'PASS','frames':nf,'waters':nw,'checks':checks,'variants':variants,'binary_sha256':sha(a.binary),'driver_sha256':sha(a.driver),'reference_sha256':sha(a.reference/'reference.json') if a.reference else None,'script_sha256':sha(__file__)},indent=2)+'\n');print('PASS CLI reference routing',checks,variants)
