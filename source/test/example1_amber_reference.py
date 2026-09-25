"""Independent AmberTools comparison fixture from example/1 PQR coordinates/charges.

This is a numerical MM+GB fixture, not reconstruction of the original MD force
field. LJ is explicitly assigned by element; mbondi/mbondi2 radii and screening
are recorded per atom. No dynamics. No SASA comparison or affinity claim.
"""
import argparse,hashlib,json,math
from pathlib import Path
import numpy as np
import parmed as pmd
import sander
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--work-root',type=Path,required=True)
p.add_argument('--example',type=Path,default=Path(__file__).resolve().parents[1]/'example/1')
p.add_argument('--rgbmax',type=float,default=999.0,help='Amber egb rgbmax in Angstrom; 999 approximates the uncut sum, 25 is the sander default')
p.add_argument('--receptor',default='rec_00.pqr',help='receptor PQR inside --example')
p.add_argument('--ligand',default='lg_00.pqr',help='ligand PQR inside --example')
a=p.parse_args();w=a.work_root.resolve();w.mkdir(parents=True,exist_ok=False)
records=[]
for filename in [a.receptor,a.ligand]:
 atoms=[]
 for line in (a.example/filename).read_text().splitlines():
  if not line.startswith(('ATOM','HETATM')):continue
  fields=line.split();name=fields[2];element=''.join(c for c in name if c.isalpha())[0]
  xyz=list(map(float,fields[-5:-2]));charge=float(fields[-2])
  assert element in 'HCNOSPF',element
  atoms.append(dict(name=name,element=element,xyz=xyz,charge=charge))
 records.append(atoms)
assert all(records)
screen={'H':.85,'C':.72,'N':.79,'O':.85,'S':.96,'P':.86,'F':.88}
radii={'C':1.7,'N':1.55,'O':1.5,'S':1.8,'P':1.85,'F':1.5}
number={'H':1,'C':6,'N':7,'O':8,'F':9,'P':15,'S':16}
refs=[]
for igb in (1,2,5):
 case=w/('igb'+str(igb));case.mkdir()
 for atoms in records:
  coords=np.array([x['xyz'] for x in atoms]);heavy=[i for i,x in enumerate(atoms) if x['element']!='H']
  for atom in atoms:
   e=atom['element'];r=radii.get(e,1.2)
   if e=='H':
    distances=np.sum((coords[heavy]-atom['xyz'])**2,axis=1);j=int(np.argmin(distances));partner=atoms[heavy[j]]['element'] if distances[j]<1.55**2 else 'N'
    r=(1.3 if partner in 'CN' else .8 if partner in 'OS' else 1.2) if igb==1 else (1.3 if partner=='N' else 1.2)
   atom['radius']=r;atom['screen']=screen[e]
   atom['sigma_nm']=.1 if e=='H' else .3;atom['epsilon_kj']=.001 if e=='H' else .1
 with (case/'fixture.txt').open('w') as f:
  f.write(f'{len(records[0])} {len(records[1])} {igb}\n')
  for atoms in records:
   for atom in atoms:f.write(' '.join(map(str,[atom['name'],atom['element'],*atom['xyz'],atom['charge'],atom['sigma_nm'],atom['epsilon_kj'],atom['radius'],atom['screen']]))+'\n')
 energies={}
 for endpoint,atoms in [('receptor',records[0]),('ligand',records[1]),('complex',records[0]+records[1])]:
  structure=pmd.Structure()
  types={}
  for atom in atoms:
   e=atom['element']
   if e not in types:
    typ=pmd.topologyobjects.AtomType(e,len(types)+1,float(number[e]*2),number[e]);typ.set_lj_params(atom['epsilon_kj']/4.184,atom['sigma_nm']*10*2**(1/6)/2);types[e]=typ
   item=pmd.Atom(name=atom['name'],type=e,atomic_number=number[e],mass=float(number[e]*2),charge=atom['charge'],solvent_radius=atom['radius'],screen=atom['screen'])
   item.atom_type=types[e];structure.add_atom(item,'TST',1)
  structure.coordinates=np.array([atom['xyz'] for atom in atoms])
  parm=pmd.amber.AmberParm.from_structure(structure)
  parm.parm_data['RADII']=[atom['radius'] for atom in atoms];parm.parm_data['SCREEN']=[atom['screen'] for atom in atoms]
  topology=case/(endpoint+'.prmtop');parm.save(str(topology),overwrite=False)
  opts=sander.gas_input(igb);opts.cut=999;opts.rgbmax=a.rgbmax;opts.saltcon=0;opts.intdiel=1;opts.extdiel=78.5;opts.gbsa=0
  with sander.setup(str(topology),structure.coordinates,None,opts):
   terms,_=sander.energy_forces()
   energies[endpoint]={name:float(getattr(terms,name))*4.184 for name in ['elec','vdw','gb']}
 refs.append(dict(igb=igb,endpoints_kj=energies,delta_kj={k:energies['complex'][k]-energies['receptor'][k]-energies['ligand'][k] for k in ['elec','vdw','gb']}))
 print('PASS Amber reference',igb,refs[-1]['delta_kj'],flush=True)
manifest={'scope':__doc__,'parmed_version':pmd.__version__,'n_receptor':len(records[0]),'n_ligand':len(records[1]),'salt_molar':0,'int_diel':1,'ext_diel':78.5,'rgbmax_angstrom':a.rgbmax,'sasa_enabled':False,'fixed_tolerance_kj':.01,'fixed_relative_tolerance':.0001,'references':refs,'example_dir':str(a.example),'receptor':a.receptor,'ligand':a.ligand,'inputs_sha256':{name:hashlib.sha256((a.example/name).read_bytes()).hexdigest() for name in [a.receptor,a.ligand]}}
(w/'reference.json').write_text(json.dumps(manifest,indent=2)+'\n')
