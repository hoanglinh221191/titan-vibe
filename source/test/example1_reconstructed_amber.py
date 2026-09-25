"""Sander MM/GB and analytical CPPTRAJ SASA on an explicitly reconstructed topology.
Original topology is unavailable. This does not certify historical example/1 energies.
"""
import argparse,hashlib,json,re,shutil,subprocess,math
from pathlib import Path
import numpy as np
import parmed as pmd
from parmed.tools.actions import changeRadii
import sander
from sasa_reference_support import exposed_sphere_mask
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--prepared',type=Path,required=True);p.add_argument('--work-root',type=Path,required=True);a=p.parse_args()
root=Path(__file__).resolve().parents[1];w=a.work_root.resolve();w.mkdir(parents=True,exist_ok=False)
ff=Path('/opt/gromacs/2026_3/share/gromacs/top/amber99sb-ildn.ff');shutil.copytree(ff,w/ff.name)
shutil.copy2(root/'example/1/MOL.top',w/'MOL.top')
text=(a.prepared/'protein.top').read_text();text=text.replace('[ moleculetype ]','#include "MOL.top"\n\n[ moleculetype ]',1);text+='\nMOL 1\n';(w/'combined.top').write_text(text)
top=pmd.gromacs.GromacsTopologyFile(str(w/'combined.top'))
nr=len(top.atoms)-43;assert nr==6063
models=[];boxes=[];current=[];current_box=None
for line in (root/'example/1/complex_2.pdb').read_text().splitlines():
 if line.startswith('CRYST1'):current_box=[float(line[i:j]) for i,j in [(6,15),(15,24),(24,33),(33,40),(40,47),(47,54)]]
 elif line.startswith('MODEL'):current=[]
 elif line.startswith(('ATOM','HETATM')):current.append((line[12:16].strip(),line[17:20].strip(),int(line[22:26]),[float(line[i:i+8]) for i in [30,38,46]]))
 elif line.startswith('ENDMDL'):models.append(current);boxes.append(current_box)
assert len(models)==52
normalize=lambda name:name[1:]+name[0] if name[:1].isdigit() else name
assert [atom.name for atom in top.atoms]==[normalize(row[0]) for row in models[0]],'Reconstructed topology/input atom order differs after standard leading-digit hydrogen alias normalization'
xyz0=np.array([row[3] for row in models[0]])
bond_distances=[float(np.linalg.norm(xyz0[b.atom1.idx]-xyz0[b.atom2.idx])) for b in top.bonds]
(w/'bond_geometry_audit.json').write_text(json.dumps({'bonds':len(bond_distances),'max_length_a':max(bond_distances),'over_2_5_a':sum(d>2.5 for d in bond_distances)},indent=2)+'\n')
assert max(bond_distances)<=2.5,'Reconstructed covalent topology contains long bonds in input coordinates'
receptor=top[:nr];ligand=top[nr:];endpoints={'receptor':receptor,'ligand':ligand,'complex':top}
parms={};sasa_radii={}
for name,st in endpoints.items():
 parm=pmd.amber.AmberParm.from_structure(st);changeRadii(parm,'mbondi2').execute();parms[name]=parm
 sasa_radii[name]=[atom.rmin for atom in parm.atoms]
 parm.save(str(w/(name+'.prmtop')))
 sasa=pmd.load_file(str(w/(name+'.prmtop')))
 for atom,radius in zip(sasa.atoms,sasa_radii[name]):atom.solvent_radius=radius
 sasa.parm_data['RADII']=sasa_radii[name];sasa.save(str(w/(name+'_sasa.prmtop')))
records=[]
for model in [0,25,51]:
 coords=np.array([row[3] for row in models[model]]);assert len(coords)==len(top.atoms)
 case=w/f'model_{model}';case.mkdir()
 box=boxes[model];assert box and all(abs(angle-90)<1e-8 for angle in box[3:]),'This fixture preprocessor requires orthogonal input cells'
 lengths=np.array(box[:3]);translation=-np.rint((coords[nr:].mean(axis=0)-coords[:nr].mean(axis=0))/lengths)*lengths
 coords[nr:]+=translation
 (case/'coordinate_preparation.json').write_text(json.dumps({'method':'rigid ligand lattice translation to nearest receptor center image; applied identically in both engines','input_cell_a':box,'ligand_translation_a':translation.tolist(),'coordinate_sha256':hashlib.sha256(coords.tobytes()).hexdigest()},indent=2)+'\n')
 assert max(np.linalg.norm(coords[b.atom1.idx]-coords[b.atom2.idx]) for b in top.bonds)<=2.5,'Frame has an unwrapped covalent-bond defect'
 areas={};energy_by_salt={}
 for name,parm in parms.items():
  xyz=coords if name=='complex' else coords[:nr] if name=='receptor' else coords[nr:]
  parm.coordinates=xyz;crd=case/(name+'.rst7');parm.save(str(crd))
  script=case/(name+'_sasa.in');out=case/(name+'_sasa.dat')
  mask,removed=exposed_sphere_mask(xyz,sasa_radii[name])
  (case/(name+'_sasa_domain.json')).write_text(json.dumps({'mask':mask,'removed_fully_contained_spheres':removed,'operation':'exact union-preserving dominance removal'},indent=2)+'\n')
  script.write_text(f'parm {w/(name+"_sasa.prmtop")}\ntrajin {crd}\nmolsurf S {mask} out {out} probe 0.0 offset 1.4 radii gb\nrun\n')
  r=subprocess.run(['/home/linh/micromamba/envs/ambertools26/bin/cpptraj','-i',str(script)],cwd=case,capture_output=True,text=True);(case/(name+'_sasa.log')).write_text(r.stdout+r.stderr);r.check_returncode()
  data=[line.split() for line in out.read_text().splitlines() if line and not line.startswith('#')];assert len(data)==1
  areas[name]=float(data[0][1]);assert math.isfinite(areas[name]) and areas[name]>0
  for salt in [0.,.15]:
   opts=sander.gas_input(5);opts.cut=999;opts.rgbmax=999;opts.saltcon=salt;opts.intdiel=1;opts.extdiel=78.5;opts.gbsa=0
   with sander.setup(str(w/(name+'.prmtop')),xyz,None,opts):
    e,_=sander.energy_forces();values={key:float(getattr(e,key))*4.184 for key in ['elec','vdw','gb','elec_14','vdw_14','bond','angle','dihedral']}
   energy_by_salt.setdefault(str(salt),{})[name]=values
 for salt in [0.,.15]:
  fixture=case/f'fixture_salt_{salt}.txt'
  with fixture.open('w') as f:
   f.write(f'{nr} 43 5 {salt}\n')
   for atom,xyz,radius in zip(parms['complex'].atoms,coords,sasa_radii['complex']):
    symbol=pmd.periodic_table.Element[atom.atomic_number]
    f.write(' '.join(map(str,[atom.name,atom.type,*xyz,atom.charge,atom.sigma*.1,atom.epsilon*4.184,radius,atom.solvent_radius,atom.screen]))+'\n')
  energy=energy_by_salt[str(salt)];delta={k:energy['complex'][k]-energy['receptor'][k]-energy['ligand'][k] for k in energy['complex']}
  ds=areas['complex']-areas['receptor']-areas['ligand'];delta['sasa_a2']=ds;delta['nonpolar']=ds*.0072*4.184
  delta['total_mmgbsa']=sum(delta[k] for k in ['elec','vdw','gb','nonpolar'])
  records.append(dict(model=model,salt=salt,fixture=str(fixture.relative_to(w)),endpoints_kj=energy,endpoint_sasa_a2=areas,delta=delta));print('REFERENCE',model,salt,delta,flush=True)
manifest=dict(scope=__doc__,topology='Reconstructed AMBER99SB-ILDN protein + retained MOL.top ligand; original atom order verified; original PDB coordinates used, not regenerated hydrogen positions.',models=[0,25,51],gb_radii='mbondi2',sasa_radii='LJ Rmin/2 copied explicitly to separate SASA prmtop RADII',sasa_reference='CPPTRAJ molsurf, probe=0, radius offset=1.4: union boundary of inflated spheres',igb=5,extdiel=78.5,intdiel=1,rgbmax_a=999,cut_a=999,gamma_kcal_a2=.0072,mmgb_tolerance_abs_kj=.01,mmgb_tolerance_rel=1e-4,sasa_delta_tolerance_abs_a2=5,sasa_delta_tolerance_rel=.01,records=records,source_sha256={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [root/'example/1/MOL.top',root/'example/1/complex_2.pdb',a.prepared.resolve()/'protein.top']})
(w/'reference.json').write_text(json.dumps(manifest,indent=2)+'\n')
