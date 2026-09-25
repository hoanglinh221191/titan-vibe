"""Independent union-of-balls area: integrate exposed azimuth, no triangles."""
import argparse
import json
import math
import numpy as np
from numpy.polynomial.legendre import leggauss

RULES=[leggauss(n) for n in (16,32)]


def integrate(f,lo,hi,tolerance=1e-10,depth=0):
    values=[]
    for x,w in RULES:
        values.append((hi-lo)/2*sum(float(wi)*f((lo+hi)/2+(hi-lo)/2*float(xi))
                                   for xi,wi in zip(x,w)))
    error=abs(values[1]-values[0])
    if error<tolerance:return values[1],error
    if depth==16:raise RuntimeError('independent area quadrature did not converge')
    mid=(lo+hi)/2
    a,ea=integrate(f,lo,mid,tolerance/2,depth+1)
    b,eb=integrate(f,mid,hi,tolerance/2,depth+1)
    return a+b,ea+eb


def union_area(atoms):
    total=error=0.
    for i,a in enumerate(atoms):
        center=np.asarray(a[:3]); radius=a[3]; caps=[]; covered=False
        for j,b in enumerate(atoms):
            if i==j:continue
            delta=np.asarray(b[:3])-center;d=np.linalg.norm(delta);r=b[3]
            if d==0 and r==radius:
                if j<i:covered=True
                continue
            if d+radius<=r:covered=True;break
            if d>=radius+r or d+r<=radius:continue
            caps.append((delta/d,(d*d+radius*radius-r*r)/(2*d*radius)))
        if covered:continue
        if not caps:total+=4*math.pi*radius*radius;continue
        points=[-1.,1.]
        for n,c in caps:
            span=math.sqrt(max(0.,(1-n[2]*n[2])*(1-c*c)))
            points.extend([float(n[2]*c-span),float(n[2]*c+span)])
        # Circle-circle intersections delimit changes of interval ordering.
        for j,(n,c) in enumerate(caps):
            for m,k in caps[:j]:
                dot=np.dot(n,m);den=1-dot*dot
                if den<1e-14:continue
                v=((c-dot*k)*n+(k-dot*c)*m)/den
                norm2=np.dot(v,v)
                if norm2<=1:
                    axis=np.cross(n,m)/math.sqrt(den)
                    dz=math.sqrt(max(0.,1-norm2))*axis[2]
                    points.extend([float(v[2]-dz),float(v[2]+dz)])
        points=sorted(set(max(-1.,min(1.,x)) for x in points))

        def exposed(z):
            intervals=[]; radial=math.sqrt(max(0.,1-z*z))
            for n,c in caps:
                amp=radial*math.hypot(n[0],n[1]);rhs=c-n[2]*z
                if amp<1e-16:
                    if rhs<0:return 0.
                    continue
                x=rhs/amp
                if x<=-1:return 0.
                if x>=1:continue
                phi=math.atan2(n[1],n[0])%(2*math.pi);half=math.acos(x)
                lo=phi-half;hi=phi+half
                if lo<0:intervals.extend([(lo+2*math.pi,2*math.pi),(0,hi)])
                elif hi>2*math.pi:intervals.extend([(lo,2*math.pi),(0,hi-2*math.pi)])
                else:intervals.append((lo,hi))
            end=covered_angle=0.
            for lo,hi in sorted(intervals):
                covered_angle+=max(0.,hi-max(lo,end));end=max(end,hi)
            return max(0.,2*math.pi-covered_angle)
        value=err=0.
        for lo,hi in zip(points,points[1:]):
            if hi-lo<1e-14:continue
            # The sine-square substitution smooths square-root tangencies at
            # both ends of every interval before Gauss integration.
            def transformed(t):
                return exposed(lo+(hi-lo)*math.sin(t)**2)*(hi-lo)*math.sin(2*t)
            part,epart=integrate(transformed,0,math.pi/2,1e-10/max(1,len(points)))
            value+=part;err+=epart
        total+=radius*radius*value;error+=radius*radius*err
    return dict(area=total,estimated_quadrature_error=error)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('input');a=p.parse_args()
    rows=[list(map(float,s.split())) for s in open(a.input).read().splitlines()[1:] if s.strip()]
    check=union_area([[0,0,0,1.7],[2.5,0,0,1.7]])
    assert abs(check['area']-(4*math.pi*1.7**2+2*math.pi*1.7*2.5))<1e-8
    print(json.dumps(union_area(rows),indent=2))
