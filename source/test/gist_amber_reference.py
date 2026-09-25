"""Independent CPPTRAJ GIST reference on a declared synthetic TIP3P fixture.

Four noninteracting alignment/ligand ghost atoms plus four waters. No MD,
bulk-reference subtraction, cutoff mismatch or historical topology claim.
All water pairs lie inside the TITAN 9 A oxygen cutoff.
"""
import argparse,hashlib,json,subprocess
from pathlib import Path
import numpy as np
import parmed as pmd
from parmed.topologyobjects import Atom,AtomType,Bond,BondType
from parmed.amber import AmberParm,NetCDFTraj
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--work-root',type=Path,required=True);p.add_argument('--cpptraj',type=Path,required=True);a=p.parse_args();w=a.work_root.resolve();w.mkdir(parents=True,exist_ok=False)
 st=pmd.Structure();types=[]
 for i,(name,mass,anum,eps,rmin) in enumerate([('DU',12.,6,0.,0.),('OW',15.9994,8,.636/4.184,3.1507*2**(1/6)/2),('HW',1.008,1,0.,0.)]):
  t=AtomType(name,i+1,mass,anum);t.set_lj_params(eps,rmin);types.append(t)
 for i,name in enumerate(['CA','CA','CA','C']):
  atom=Atom(name=name,type='DU',atomic_number=6,mass=12.,charge=0.);atom.atom_type=types[0];st.add_atom(atom,'DUM',i+1)
 bond=BondType(450.,.9572);st.bond_types.append(bond)
 for water in range(4):
  start=len(st.atoms)
  for name,typeid,charge in [('O',1,-.834),('H1',2,.417),('H2',2,.417)]:
   t=types[typeid];atom=Atom(name=name,type=t.name,atomic_number=t.atomic_number,mass=t.mass,charge=charge);atom.atom_type=t;st.add_atom(atom,'WAT',5+water)
  st.bonds.extend([Bond(st.atoms[start],st.atoms[start+1],type=bond),Bond(st.atoms[start],st.atoms[start+2],type=bond)])
 st.box=[50.,50.,50.,90.,90.,90.];st.coordinates=np.zeros((16,3));parm=AmberParm.from_structure(st);parm.save(str(w/'fixture.prmtop'),overwrite=False)
 nframes=1024;rng=np.random.default_rng(8675309);frames=[];centers=np.array([[2.2,2.2,.2],[-2.2,2.2,.2],[2.2,-2.2,.2],[-2.2,-2.2,.2]])
 h=np.array([[.9572,0,0],[.9572*np.cos(np.deg2rad(104.52)),.9572*np.sin(np.deg2rad(104.52)),0]])
 quat=rng.normal(size=(nframes,4,4));quat/=np.linalg.norm(quat,axis=2)[:,:,None]
 def matrix(q):
  w,x,y,z=q;return np.array([[1-2*(y*y+z*z),2*(x*y-z*w),2*(x*z+y*w)],[2*(x*y+z*w),1-2*(x*x+z*z),2*(y*z-x*w)],[2*(x*z-y*w),2*(y*z+x*w),1-2*(x*x+y*y)]])
 for f in range(nframes):
  coords=np.zeros((16,3));coords[:4]=[[-1,-1,0],[1,-1,0],[0,1,0],[0,0,0]]
  for j,c in enumerate(centers):coords[4+3*j]=c;coords[5+3*j:7+3*j]=c+h@matrix(quat[f,j]).T
  frames.append(coords)
 frames=np.array(frames);np.save(w/'coordinates.npy',frames)
 # Batch-write standard AMBER NetCDF. The per-frame ParmEd scipy fallback
 # repeatedly rewrites the file and is prohibitively slow on mounted storage.
 from scipy.io import netcdf_file
 with netcdf_file(str(w/'fixture.nc'),'w') as traj:
  traj.Conventions='AMBER';traj.ConventionVersion='1.0';traj.application='TITAN independent numerical fixture';traj.program='scipy';traj.programVersion='fixture-1'
  for name,size in [('frame',nframes),('atom',16),('spatial',3),('cell_spatial',3),('cell_angular',3),('label',5)]:traj.createDimension(name,size)
  xyz=traj.createVariable('coordinates','f',('frame','atom','spatial'));xyz.units='angstrom';xyz[:]=frames
  time=traj.createVariable('time','f',('frame',));time.units='picosecond';time[:]=np.arange(nframes)
  lengths=traj.createVariable('cell_lengths','d',('frame','cell_spatial'));lengths.units='angstrom';lengths[:]=np.full((nframes,3),50.)
  angles=traj.createVariable('cell_angles','d',('frame','cell_angular'));angles.units='degree';angles[:]=np.full((nframes,3),90.)
 # Read back the NetCDF (typically float32) and give exactly those coordinates
 # to TITAN and the independent numerical oracle.
 traj=NetCDFTraj.open_old(str(w/'fixture.nc'));actual=np.array(traj.coordinates);traj.close()
 with (w/'titan_fixture.txt').open('w') as out:
  out.write(f'{nframes} 4\n')
  for frame in actual:
   for xyz in frame[4:]:out.write(' '.join(format(float(v),'.17g') for v in xyz)+'\n')
 cmd=f'parm {w}/fixture.prmtop\ntrajin {w}/fixture.nc\ngist gridcntr 0 0 0 griddim 24 24 24 gridspacn 0.5 temp 300 refdens 0.0334 nopme nocom prefix {w}/amber floatfmt scientific floatprec 12\nrun\n'
 (w/'cpptraj.in').write_text(cmd);r=subprocess.run([str(a.cpptraj),'-i',str(w/'cpptraj.in')],capture_output=True,text=True,env=__import__('os').environ|{'OMP_NUM_THREADS':'2'});(w/'cpptraj.log').write_text(r.stdout+r.stderr);r.check_returncode()
 # Independent reference: all six water pairs, explicit nine charge pairs;
 # entropy from matrix-relative trace rather than production quaternion KD.
 qcharge=np.array([-.834,.417,.417]);energies=np.zeros((nframes,4));orientation=[]
 for i in range(4):
  waters=actual[:,4+3*i:7+3*i];x=waters[:,1]-waters[:,0];x/=np.linalg.norm(x,axis=1)[:,None];y=waters[:,2]-waters[:,0];y-=np.sum(x*y,axis=1)[:,None]*x;y/=np.linalg.norm(y,axis=1)[:,None];orientation.append(np.stack([x,y,np.cross(x,y)],axis=2))
 for i in range(4):
  for j in range(i+1,4):
   wi=actual[:,4+3*i:7+3*i].astype(float);wj=actual[:,4+3*j:7+3*j].astype(float);d=np.linalg.norm(wi[:,:,None,:]-wj[:,None,:,:],axis=3);sr6=(3.1507/d[:,0,0])**6;e=1389.3546*np.sum(qcharge[:,None]*qcharge[None,:]/d,axis=(1,2))+4*.636*(sr6**2-sr6);energies[:,i]+=.5*e;energies[:,j]+=.5*e
 records=[]
 for j,mat in enumerate(orientation):
  # Trace(R_i^T R_j) gives 1+2 cos(theta); float64 recomputation avoids
  # float32 orientation arithmetic degrading the external comparison.
  waters=actual[:,4+3*j:7+3*j].astype(float);x=waters[:,1]-waters[:,0];x/=np.linalg.norm(x,axis=1)[:,None];y=waters[:,2]-waters[:,0];y-=np.sum(x*y,axis=1)[:,None]*x;y/=np.linalg.norm(y,axis=1)[:,None];mat=np.stack([x,y,np.cross(x,y)],axis=2)
  dots=mat.reshape(nframes,9)@mat.reshape(nframes,9).T;np.fill_diagonal(dots,-1);theta=np.arccos(np.clip((dots.max(axis=1)-1)/2,-1,1));s=np.mean(np.log(nframes*(theta-np.sin(theta))/np.pi))+np.euler_gamma
  records.append({'water':j,'center_A':centers[j].tolist(),'Eww_kj':float(energies[:,j].mean()),'minus_TSorient_kj':float(-.008314462*300*s)})
 result={'status':'REFERENCE_COMPLETE','scope':__doc__,'frames':nframes,'waters':4,'numerical_reference':records,'inputs_sha256':{p.name:sha(p) for p in [w/'fixture.prmtop',w/'fixture.nc',w/'titan_fixture.txt',w/'cpptraj.in']},'cpptraj_sha256':sha(a.cpptraj),'script_sha256':sha(__file__)}
 (w/'reference.json').write_text(json.dumps(result,indent=2)+'\n');print('REFERENCE_COMPLETE',flush=True)
if __name__=='__main__':main()
