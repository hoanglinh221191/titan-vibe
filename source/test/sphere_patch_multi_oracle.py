"""Independent axisymmetric dielectric-sphere oracle, no TITAN surface/kernel.

Centers must be collinear, spheres equal-radius and disjoint, charges centered.
Interior regular and exterior decaying Legendre expansions enforce potential
and dielectric-flux continuity at Gauss-Legendre collocation points. Increase
the expansion degree independently of the TITAN panel count.
"""
import json
import numpy as np
from numpy.polynomial.legendre import leggauss, legvander


def energy(z, charges, radius=1.7, degree=32, epsin=1., epsout=78.5):
    z=np.asarray(z)/radius
    q=np.asarray(charges)
    count=len(z); m=degree+1
    if count>1 and min(abs(z[i]-z[j]) for i in range(count) for j in range(i))<=2:
        raise ValueError('oracle requires disjoint spheres')
    mu,_=leggauss(m); P=legvander(mu,degree); l=np.arange(m)
    A=np.zeros((2*count*m,2*count*m)); rhs=np.zeros(2*count*m)
    for i in range(count):
        pot=slice(2*i*m,(2*i+1)*m); flux=slice((2*i+1)*m,(2*i+2)*m)
        ai=slice(2*i*m,(2*i+1)*m)
        A[pot,ai]=P; A[flux,ai]=(epsin/epsout)*P*l
        rhs[pot]=-q[i]/epsin; rhs[flux]=q[i]/epsout
        for j in range(count):
            bj=slice((2*j+1)*m,(2*j+2)*m)
            dz=z[i]-z[j]; r=np.sqrt(1+dz*dz+2*dz*mu)
            nu=(dz+mu)/r; V=legvander(nu,degree)
            dV=np.zeros_like(V)
            for k in range(1,m):dV[:,k]=k*(V[:,k-1]-nu*V[:,k])/(1-nu*nu)
            dr=(1+dz*mu)/r; dnu=(mu-nu*dr)/r
            power=r[:,None]**(-l-1)
            derivative=power*(dV*dnu[:,None]-V*(l+1)*dr[:,None]/r[:,None])
            A[pot,bj]=-power*V; A[flux,bj]=-derivative
    coefficients=np.linalg.solve(A,rhs)
    residual=float(np.linalg.norm(A@coefficients-rhs)/np.linalg.norm(rhs))
    # a_i contains the total regular field in cavity i, including the other
    # charges. Reaction potential subtracts their bare Coulomb field too.
    value=332.063713299/(2*radius)*sum(q[i]*(coefficients[2*i*m]-
        sum(q[j]/(epsin*abs(z[i]-z[j])) for j in range(count) if i!=j))
        for i in range(count))
    return dict(energy=float(value),residual=residual,degree=degree)


if __name__=='__main__':
    single=energy([0],[1]); reference=-332.063713299/(2*1.7)*(1-1/78.5)
    assert abs(single['energy']/reference-1)<1e-12
    assert abs(energy([0,5],[1,-1],epsout=1)['energy'])<1e-10
    results={name:[energy(z,q,degree=d) for d in (16,32,64)] for name,z,q in (
        ('disjoint2',[0,5],[1,-1]),('disjoint3',[0,5,10],[1,-.5,.3]))}
    for values in results.values():
        assert abs(values[-1]['energy']-values[-2]['energy'])<1e-9
        assert max(v['residual'] for v in values)<1e-12
    print(json.dumps(results,indent=2))
