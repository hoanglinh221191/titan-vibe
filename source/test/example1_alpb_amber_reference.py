"""ALPB reference using AmberTools elsize -det and explicit Sander arad.

Original simulation topology is unavailable. Reuses the qualified reconstructed
topology and coordinates; this does not certify historical example/1 energies.
"""
import argparse,hashlib,json,re,subprocess
from pathlib import Path
import parmed as pmd

p=argparse.ArgumentParser(description=__doc__);p.add_argument('--reference-root',type=Path,required=True);p.add_argument('--amber-bin',type=Path,required=True);p.add_argument('--work-root',type=Path,required=True);a=p.parse_args();w=a.work_root.resolve();w.mkdir(parents=True,exist_ok=False)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
base=json.loads((a.reference_root/'reference.json').read_text());records=[];inputs={str(a.reference_root/'reference.json'):sha(a.reference_root/'reference.json')}
def elsize(path):
 r=subprocess.run([str(a.amber_bin/'elsize'),str(path),'-det'],capture_output=True,text=True);r.check_returncode();(path.with_suffix('.elsize.log')).write_text(r.stdout+r.stderr)
 return float(r.stdout.strip().split()[-1])
sphere=w/'sphere.pqr';sphere.write_text('ATOM      1  X   NUM     1        0.00000000     0.00000000     0.00000000  0.00000000 2.00000000\n');sphere_a=elsize(sphere);assert abs(sphere_a-2)<1e-5
for model in [0,25,51]:
 case=w/f'model_{model}';case.mkdir();endpoints={};sizes={}
 for name in ['receptor','ligand','complex']:
  top=(a.reference_root/f'{name}.prmtop').resolve();coords=(a.reference_root/f'model_{model}/{name}.rst7').resolve();inputs[str(top)]=sha(top);inputs[str(coords)]=sha(coords)
  parm=pmd.load_file(str(top),xyz=str(coords));pqr=case/f'{name}.pqr'
  # Match the declared TITAN A_det size convention (LJ Rmin/2 radii), while GB
  # intrinsic radii in the unchanged prmtop remain mbondi2.
  with pqr.open('w') as f:
   for i,(atom,xyz) in enumerate(zip(parm.atoms,parm.coordinates),1):
    f.write(f'ATOM {i:6d} {atom.name:4s} NUM {i:6d} {xyz[0]:.10f} {xyz[1]:.10f} {xyz[2]:.10f} {atom.charge:.10f} {atom.rmin:.10f}\n')
  size=elsize(pqr);sizes[name]=size
  for salt in [0.,.15]:
   stem=case/f'{name}_salt_{salt}';inp=stem.with_suffix('.in');out=stem.with_suffix('.out')
   # Use string suffixes: Path.with_suffix would collapse the decimal salt.
   inp=Path(str(stem)+'.in');out=Path(str(stem)+'.out')
   inp.write_text(f'ALPB reconstructed single point\n&cntrl\n imin=1,maxcyc=0,ncyc=0,ntb=0,igb=5,alpb=1,arad={size:.10f},\n cut=999,rgbmax=999,intdiel=1,extdiel=78.5,saltcon={salt},gbsa=0,\n ntpr=1,ntc=1,ntf=1,\n/\n')
   r=subprocess.run([str(a.amber_bin/'sander'),'-O','-i',str(inp),'-o',str(out),'-p',str(top),'-c',str(coords),'-r',str(stem)+'.rst7','-inf',str(stem)+'.mdinfo'],cwd=case,capture_output=True,text=True,timeout=300)
   Path(str(stem)+'.run.log').write_text(r.stdout+r.stderr);r.check_returncode()
   text=out.read_text();energy_lines=[line for line in text.splitlines() if line.strip().startswith('VDWAALS')]
   assert energy_lines,'No final energy record'
   values={k:float(v)*4.184 for k,v in re.findall(r'(VDWAALS|EEL|EGB)\s*=\s*([+\-\d.Ee]+)',energy_lines[-1])};assert len(values)==3
   prior=next(r for r in base['records'] if r['model']==model and r['salt']==salt)['endpoints_kj'][name]
   # Sander's ALPB branch also scales direct EEL by 1/(1+alpha*beta).
   # AmberClassic egb.F90 documents this explicitly at intdieli/temp1.
   # Verify that convention against the independently computed non-ALPB EEL;
   # retain raw EEL and report a separate unscaled vacuum-MM convention.
   scale=1+.571412/78.5
   assert abs(values['EEL']-prior['elec']/scale)<.001,(model,name,'ALPB EEL scale',values['EEL'],prior['elec']/scale)
   assert abs(values['VDWAALS']-prior['vdw'])<.001
   values['EEL_vacuum_convention']=values['EEL']*scale
   endpoints.setdefault(str(salt),{})[name]=values
 for salt in [0.,.15]:
  en=endpoints[str(salt)];delta={k:en['complex'][k]-en['receptor'][k]-en['ligand'][k] for k in ['EEL','EEL_vacuum_convention','VDWAALS','EGB']}
  records.append({'model':model,'salt':salt,'arad_A':sizes,'endpoints_kj_mol':en,'delta_kj_mol':delta});print('ALPB reference',model,salt,delta,sizes,flush=True)
(w/'reference.json').write_text(json.dumps({'status':'REFERENCE_COMPLETED','scope':__doc__,'igb':5,'alpb':1,'intdiel':1,'extdiel':78.5,'arad_method':'AmberTools elsize -det, PQR LJ Rmin/2; explicitly set independently for each endpoint/frame','arad_sphere_oracle_A':sphere_a,'eel_convention':'Raw Sander ALPB EEL scales direct Coulomb by 1/(1+0.571412/78.5). EEL_vacuum_convention explicitly undoes it and is checked against non-ALPB Sander. Raw total energy parity is not asserted.','source_explanation':'https://raw.githubusercontent.com/Amber-MD/AmberClassic/main/src/msander/egb.F90','records':records,'inputs_sha256':inputs,'script_sha256':sha(__file__)},indent=2)+'\n')
