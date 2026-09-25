"""Exact removal of spheres that cannot contribute to an inflated-sphere union boundary."""
import numpy as np
def exposed_sphere_mask(coordinates,radii,probe=1.4):
 xyz=np.asarray(coordinates,dtype=float);r=np.asarray(radii,dtype=float)+probe
 if xyz.shape!=(len(r),3) or not np.isfinite(xyz).all() or not np.isfinite(r).all() or (r<=0).any():raise ValueError('Invalid SASA sphere input')
 keep=[]
 for i in range(len(r)):
  bigger=np.flatnonzero((r>r[i]+1e-10)|((abs(r-r[i])<=1e-10)&(np.arange(len(r))<i)))
  if len(bigger):
   d=np.linalg.norm(xyz[bigger]-xyz[i],axis=1)
   if np.any(d+r[i]<=r[bigger]+1e-10):continue
  keep.append(i)
 if not keep:raise ValueError('SASA dominance filter removed every sphere')
 # Compact one-based CPPTRAJ atom mask.
 ranges=[];start=last=keep[0]+1
 for value in [i+1 for i in keep[1:]]:
  if value==last+1:last=value
  else:ranges.append(str(start) if start==last else f'{start}-{last}');start=last=value
 ranges.append(str(start) if start==last else f'{start}-{last}')
 return '@'+','.join(ranges),len(r)-len(keep)
