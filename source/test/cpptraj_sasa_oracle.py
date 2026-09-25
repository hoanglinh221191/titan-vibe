"""Validate CPPTRAJ's inflated-sphere/probe-zero SASA reference analytically."""
import argparse,math,subprocess,json
from pathlib import Path
import parmed as pmd,numpy as np
from sasa_reference_support import exposed_sphere_mask
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cpptraj',required=True);p.add_argument('--work-root',type=Path,required=True);a=p.parse_args();w=a.work_root.resolve();w.mkdir(parents=True,exist_ok=False)
cases=[('isolated',[1.7],[[0,0,0]],4*math.pi*3.1**2),('overlap',[1.7,1.7],[[0,0,0],[4,0,0]],4*math.pi*3.1**2+2*math.pi*3.1*4),('contained',[1.7,.6],[[0,0,0],[.2,0,0]],4*math.pi*3.1**2)]
results=[]
for name,radii,coords,expected in cases:
 st=pmd.Structure();ty=pmd.topologyobjects.AtomType('C',1,12.,6);ty.set_lj_params(.1,1.7)
 for i,r in enumerate(radii):
  atom=pmd.Atom(name='C',type='C',atomic_number=6,mass=12.,charge=0,solvent_radius=r,screen=.72);atom.atom_type=ty;st.add_atom(atom,'TST',i+1)
 st.coordinates=np.array(coords);parm=pmd.amber.AmberParm.from_structure(st);top=w/(name+'.prmtop');crd=w/(name+'.rst7');parm.save(str(top));parm.coordinates=np.array(coords);parm.save(str(crd))
 mask,removed=exposed_sphere_mask(coords,radii)
 out=w/(name+'.dat');inp=w/(name+'.in');inp.write_text(f'parm {top}\ntrajin {crd}\nmolsurf S {mask} out {out} probe 0.0 offset 1.4 radii gb\nrun\n')
 r=subprocess.run([a.cpptraj,'-i',str(inp)],cwd=w,capture_output=True,text=True);(w/(name+'.log')).write_text(r.stdout+r.stderr);r.check_returncode();value=float([x for x in out.read_text().splitlines() if x and not x.startswith('#')][0].split()[1]);error=abs(value-expected)
 results.append({'case':name,'expected_a2':expected,'cpptraj_a2':value,'absolute_error_a2':error,'status':'PASS' if error<.0005 else 'FAIL'});print(results[-1],flush=True)
with (w/'verification.json').open('w') as f:json.dump({'cases':results},f,indent=2);f.write('\n')
assert all(r['status']=='PASS' for r in results)
