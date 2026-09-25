"""Independent single-point Amber references. No MD and no fitted parameters.

Run with AmberTools Python. Known mbondi3 radii are explicit in the small
fixtures so this validates TITAN's assignment as well as the GB calculation.
"""
import argparse, hashlib, json, os, re, subprocess
from pathlib import Path
import numpy as np
import parmed as pmd
import sander

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--work-root',type=Path,required=True)
p.add_argument('--include-examples',action='store_true')
a=p.parse_args(); w=a.work_root.resolve(); w.mkdir(parents=True,exist_ok=False)
Z={'H':1,'C':6,'N':7,'O':8,'P':15,'S':16,'F':9}
# name, element/type, residue, residue number, x,y,z,q, radius
fixtures={
 'protein': ([['C','C','ALA',1,0.,0.,0.,.2,1.7],['N','N','ALA',1,1.4,0.,0.,-.3,1.55],
              ['H','H','ALA',1,2.35,0.,0.,.1,1.3],['S','S','CYS',2,-1.8,.3,.2,-.2,1.8]],
             [['C','C','LIG',1,4.,1.,0.,.4,1.7],['O','O','LIG',1,5.2,1.,0.,-.4,1.5]]),
 'mbondi3': ([['OD1','O','ASP',1,0.,0.,0.,-.5,1.4],['OD2','O','ASP',1,0.,2.2,0.,-.5,1.4],
              ['HE','H','ARG',2,3.7,0.,0.,.4,1.17],['NE','N','ARG',2,2.7,0.,0.,-.4,1.55]],
             [['O','O','ALA',1,5.5,1.,0.,-.4,1.4],['OXT','O','ALA',1,5.5,3.,0.,-.6,1.4]]),
 'rna': ([['P','P','A',1,0.,0.,0.,1.,1.85],['O1P','O','A',1,1.5,0.,0.,-.7,1.5],
         ['C1','C','A',1,0.,2.5,0.,.2,1.7],['N1','N','A',1,0.,4.,0.,-.3,1.55],
         ['H1','H','A',1,0.,5.,0.,.2,1.3]],
        [['C1','C','DA5',1,5.,2.,0.,.1,1.7],['N1','N','DA5',1,6.4,2.,0.,-.1,1.55],
         ['P','P','DA5',1,7.,4.,0.,.5,1.85]]),
 'spheres': ([['C','C','LIG',1,0.,0.,0.,1.,1.7]],
             [['O','O','LIG',1,4.,0.,0.,-1.,1.5]]),
 'neck_surface': ([['C','C','LIG',1,0.,0.,0.,.5,1.7],
                   ['C2','C','LIG',1,1.5,0.,0.,.5,1.7]],
                  [['O','O','LIG',1,3.8,.4,0.,-1.,1.5]]),
}

if a.include_examples:
 for number,files in [(1,['rec_00.pqr','lg_00.pqr']),(4,['rec_0.pqr','lg_0.pqr'])]:
  groups=[]
  for filename in files:
   path=Path(__file__).resolve().parents[1]/'example'/str(number)/filename
   atoms=[]
   for line in path.read_text().splitlines():
    if not line.startswith(('ATOM','HETATM')): continue
    f=line.split();name=f[2];e=next(c for c in name if c.isalpha());res=f[3]
    atoms.append([name,e,res,int(f[-6]),*map(float,f[-5:-2]),float(f[-2]),0.])
   xyz=np.array([r[4:7] for r in atoms]);heavy=[i for i,r in enumerate(atoms) if r[1]!='H']
   for row in atoms:
    name,e,res,resnr=row[:4];r={'C':1.7,'N':1.55,'O':1.5,'F':1.5,'P':1.85,'S':1.8}.get(e,1.2)
    if e=='H':
     ds=np.sum((xyz[heavy]-row[4:7])**2,axis=1);j=int(np.argmin(ds))
     r=1.3 if ds[j]<1.55**2 and atoms[heavy[j]][1]=='N' else 1.2
     if res=='ARG' and name.startswith(('HH','HE')):r=1.17
    if e=='O' and ((res=='ASP' and name.startswith('OD')) or (res=='GLU' and name.startswith('OE')) or name in ['OXT','OT1','OT2']):r=1.4
    if e=='O' and name=='O' and any(t[3]==resnr and t[0]=='OXT' for t in atoms):r=1.4
    row[-1]=r
   groups.append(atoms)
  fixtures[f'example{number}']=tuple(groups)

def topology(atoms, folder, endpoint):
 s=pmd.Structure(); types={}
 for row in atoms:
  name,e,res,resnr,x,y,z,q,r=row
  if e not in types:
   t=pmd.topologyobjects.AtomType(e,len(types)+1,float(Z[e]*2),Z[e]);t.set_lj_params(.1/4.184,3.*2**(1/6)/2);types[e]=t
  atom=pmd.Atom(name=name,type=e,atomic_number=Z[e],mass=float(Z[e]*2),charge=q,solvent_radius=r,screen=.8)
  # TITAN accepts CHARMM/PDB RNA names. Amber's isnucat needs canonical names.
  amber_res={'ADE':'A','CYT':'C','GUA':'G','URA':'U','THY':'DT'}.get(res,res)
  atom.atom_type=types[e];s.add_atom(atom,amber_res,resnr)
 s.coordinates=np.array([row[4:7] for row in atoms])
 parm=pmd.amber.AmberParm.from_structure(s)
 parm.parm_data['RADII']=[row[-1] for row in atoms]
 top=folder/(endpoint+'.prmtop'); crd=folder/(endpoint+'.rst7')
 parm.save(str(top));parm.save(str(crd),format='rst7')
 return top,crd,s.coordinates

records=[]
for name,(rec,lig) in fixtures.items():
 for model in ([66] if name in ['spheres','neck_surface'] else [8]):
  for salt in [0.,.15]:
   for cutoff in ([0.] if model==66 else [25.,999.]):
    key=f'{name}_igb{model}_salt{salt}_cut{cutoff}'
    folder=w/key;folder.mkdir()
    fixture=folder/'fixture.txt'
    fixture.write_text(f'{len(rec)} {len(lig)} {model} {salt} {cutoff if cutoff<900 else 0} 0.25\n'+
        ''.join(' '.join(map(str,row))+'\n' for row in rec+lig))
    endpoints={}
    for endpoint,atoms in [('rec',rec),('lig',lig),('complex',rec+lig)]:
     top,crd,xyz=topology(atoms,folder,endpoint)
     if model==8:
      opts=sander.gas_input(8);opts.cut=999;opts.rgbmax=cutoff;opts.saltcon=salt
      opts.intdiel=1;opts.extdiel=78.5;opts.gbsa=0
      with sander.setup(str(top),xyz,None,opts):
       energy,_=sander.energy_forces()
       endpoints[endpoint]={'elec':energy.elec*4.184,'gb':energy.gb*4.184}
     else:
      mdin=folder/(endpoint+'.in');out=folder/(endpoint+'.out')
      mdin.write_text('GBNSR6 reference\n &cntrl\n inp=1\n /\n &gb\n '+
        # gbnsr6 uses mM and fixes T=300 K; compare kappa at TITAN's 298.15 K.
        f'epsin=1.,epsout=78.5,istrng={salt*1000.*300./298.15},space=0.125,arcres=0.1,dprob=1.4,B=0.028,alpb=1,rbornstat=1,chagb=0\n /\n')
      command=[str(Path(os.environ['AMBERHOME'])/'bin/gbnsr6')]
      # Amber uses fixed-width path buffers; filenames relative to cwd avoid truncation.
      r=subprocess.run(command+['-i',mdin.name,'-o',out.name,'-p',top.name,'-c',crd.name],cwd=folder,capture_output=True,text=True)
      (folder/(endpoint+'.log')).write_text(r.stdout+r.stderr);r.check_returncode()
      text=out.read_text();vals={k:float(v)*4.184 for k,v in re.findall(r'\b(EELEC|EGB)\s*=\s*([-+\d.Ee]+)',text)}
      assert set(vals)=={'EELEC','EGB'},text[-2000:]
      endpoints[endpoint]={'elec':vals['EELEC'],'gb':vals['EGB'],
         'radii_A':[1./float(v) for v in re.findall(r'rinv\s+\d+\s+([-+\d.Ee]+)',text)]}
    records.append({'case':key,'igb':model,'fixture':str(fixture.relative_to(w)),
       'endpoints_kj':endpoints,'delta_kj':{k:endpoints['complex'][k]-endpoints['rec'][k]-endpoints['lig'][k] for k in ['elec','gb']}})
    print(key,records[-1]['delta_kj'],flush=True)
manifest={'scope':__doc__,'records':records,'parmed_version':pmd.__version__,
 'nucleic_residue_map':{'ADE':'A','CYT':'C','GUA':'G','URA':'U','THY':'DT'},
 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'files_sha256':{str(f.relative_to(w)):hashlib.sha256(f.read_bytes()).hexdigest() for f in w.rglob('*') if f.is_file()}}
(w/'reference.json').write_text(json.dumps(manifest,indent=2)+'\n')
