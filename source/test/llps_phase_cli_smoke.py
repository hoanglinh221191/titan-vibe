"""End-to-end CLI parsing and explicit LLPS phase observation windows."""
import argparse,hashlib,json,math,subprocess,csv
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary',type=Path,required=True);p.add_argument('--work-root',type=Path,required=True)
a=p.parse_args();w=a.work_root.resolve();w.mkdir(parents=True,exist_ok=False);binary=a.binary.resolve()
coords=[(0,0,2),(0,0,-2),(0,0,25),(0,0,70),(90,90,90)]
top=['[ defaults ]','1 2 yes 0.5 0.8333333333','[ atomtypes ]','X 6 12.0 0.0 A 0.3 0.2','[ moleculetype ]','Fixture 0','[ atoms ]']
for i in range(5):top.append(f'{i+1} X {i+1} ALA {"CA" if i<4 else "X"} {i+1} 0 12')
top+=['[ system ]','Numerical phase fixture','[ molecules ]','Fixture 1']
(w/'fixture.top').write_text('\n'.join(top)+'\n');(w/'fixture.ndx').write_text('[ System ]\n1 2 3 4 5\n[ Receptor ]\n1 2 3 4\n[ Ligand ]\n5\n')
pdb=['CRYST1  200.000  200.000  200.000  80.00  85.00  65.00 P 1           1']
for frame in [1,2]:
 pdb += [f'TITLE Synthetic phase t= {frame-1}',f'MODEL     {frame:4d}']
 for i,(x,y,z) in enumerate(coords,1):
  atom='CA' if i<5 else 'X';pdb.append(f'ATOM  {i:5d} {atom:^4s} ALA {chr(64+i)}{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C');pdb.append('TER')
 pdb.append('ENDMDL')
pdb.append('END');(w/'fixture.pdb').write_text('\n'.join(pdb)+'\n')
rows=[]
for name in ['droplet','slab','invalid_window','dbscan_enabled','paper_guard','paper_droplet','paper_slab']:
 case=w/name;case.mkdir();geometry='slab' if name=='slab' else 'droplet'
 opts=[f'-f={(w/"fixture.pdb").as_posix()}','-1pdb=true','-ff=AMBER99SB',f'-n={(w/"fixture.ndx").as_posix()}','-g=1-2',
  '-mol-1-type=protein','-mol-2-type=lig',f'-path-a={w.as_posix()}/',f'-path-n={w.as_posix()}/',f'-common_topo={(w/"fixture.top").as_posix()}',
  '-nt=1','-llps=true',f'-llps-dbscan={"true" if name in ["dbscan_enabled","paper_guard","paper_droplet","paper_slab"] else "false"}','-llps-density=off','-llps-dbscan-eps=7','-llps-dbscan-minpts=2',
  f'-llps-phase-geometry={geometry}','-llps-phase-axis=c',f'-llps-phase-dense-nm={3 if name=="invalid_window" else 1}',
  '-llps-phase-dilute-inner-nm=2','-llps-phase-dilute-outer-nm=3']
 if name=='paper_guard':opts+=['-llps-paper=true']
 if name in ['paper_droplet','paper_slab']:opts+=['-llps-paper=true',f'-llps-paper-geometry={name[6:]}']
 (case/'arg.in').write_text('\n'.join(opts)+'\n')
 r=subprocess.run([str(binary),'-a',str(case/'arg.in')],cwd=case,capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=120)
 output=r.stdout+r.stderr;(case/'result.log').write_text(output,encoding='utf-8')
 if name in ['invalid_window','paper_guard']:
  expected='phase windows require' if name=='invalid_window' else 'paper profiles require'
  if r.returncode==0 or expected not in output:raise RuntimeError(name+' incorrect rejection: '+output[-1600:])
 else:
  if r.returncode:raise RuntimeError(name+' failed: '+output[-1600:])
  files=[x for x in case.rglob('*') if x.is_file() and 'timeseries' in x.name and 'llps' in x.name.lower()]
  tables=[]
  for f in files:
   text=f.read_text(errors='replace')
   if 'phase_status' in text:tables.append((f,text))
  if len(tables)!=1:raise RuntimeError('missing phase timeseries '+str(files))
  data=[line.split() for line in tables[0][1].splitlines() if line and not line.startswith('#')]
  if len(data)!=2 or any(row[-1]!='ok' or int(row[8])!=2 or int(row[9])!=1 for row in data):raise RuntimeError('wrong phase counts '+str(data))
  if name=='dbscan_enabled':
   cluster_rows=list(csv.DictReader((case/'llps_dbscan_timeseries.csv').open()))
   if len(cluster_rows)!=2 or any(int(row['NumCondensates'])!=1 or int(row['GiantCluster_ProChains'])!=2 or int(row['NoiseBeads'])!=2 for row in cluster_rows):raise RuntimeError('wrong DBSCAN CLI clusters '+str(cluster_rows))
  if name in ['paper_droplet','paper_slab']:
   profile=case/('llps_slab_profiles.csv' if name=='paper_slab' else 'llps_radial_profiles.csv')
   bins=list(csv.DictReader(profile.open()));assert bins and all(int(row['Frames'])==2 for row in bins)
   volume=sum(float(row['ShellVolumeMean_nm3']) for row in bins)
   counts=sum(float(row['ProteinDensity'])*float(row['ShellVolumeMean_nm3']) for row in bins)
   assert abs(counts-4)<2e-3,counts
   if name=='paper_slab':
    ca,cb,cg=[math.cos(math.radians(x)) for x in [80,85,65]];expected=8000*math.sqrt(1-ca*ca-cb*cb-cg*cg+2*ca*cb*cg)
    assert abs(volume-expected)<1e-4,(volume,expected)
 rows.append({'case':name,'exit_code':r.returncode,'status':'PASS'});print('PASS LLPS CLI',name,flush=True)
(w/'verification.json').write_text(json.dumps(dict(binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),cases=rows),indent=2)+'\n')
