"""Reproduce the optional experiments; these numerical checks are not proofs.

Python and NumPy. Outputs JSON under the selected output directory.
No external data, randomness, files to download, or network access.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
from pathlib import Path
import json
import numpy as np

PI = np.pi
ALPHA, DELTA, ETA = 1/3, 1/4, 1/2
TIMES = np.linspace(0, 0.01, 101)
OBS_TIMES = np.array([1e-3,5e-4,2e-4,1e-4,5e-5])
ALL_TIMES = np.unique(np.concatenate((TIMES,OBS_TIMES)))

def i2(z):
    term=z*z/8
    out=term.copy()
    for m in range(1,20):
        term=term*(z*z/4)/(m*(m+2))
        out+=term
    return out

def rk4(rhs,y,times,dt):
    t=0.
    values=[]
    for end in times:
        steps=int(np.ceil(max(0.,end-t)/dt-1e-10))
        if steps:
            h=(end-t)/steps
            for _ in range(steps):
                k1=rhs(t,y)
                k2=rhs(t+h/2,y+h*k1/2)
                k3=rhs(t+h/2,y+h*k2/2)
                k4=rhs(t+h,y+h*k3)
                y=y+h*(k1+2*k2+2*k3+k4)/6
                t+=h
        t=end
        values.append(y.copy())
    return np.stack(values,axis=1)

class SampledSolution:
    def __init__(self,values):
        self.values=values
    def sol(self,times):
        ids=np.searchsorted(ALL_TIMES,times)
        assert np.max(np.abs(ALL_TIMES[ids]-times))<1e-14
        return self.values[:,ids]

def initial(n):
    a = np.arange(n)/n
    x = a.copy()
    for _ in range(12):
        x -= (x + ALPHA*np.sin(2*PI*x)/(2*PI)-a)/(1+ALPHA*np.cos(2*PI*x))
    q = 1/(1+ALPHA*np.cos(2*PI*x))
    return x, q

def fields(x, weight):
    c, s = np.cos(2*PI*x), np.sin(2*PI*x)
    c0, s0 = c.mean(), s.mean()
    ch, sh = (c*weight).mean(), (s*weight).mean()
    force = 2*PI*(s*c0-c*s0)
    fx = (2*PI)**2*(c*c0+s*s0)
    mobility = 1+weight.mean()+ETA*(c*ch+s*sh)
    mx = 2*PI*ETA*(-s*ch+c*sh)
    return -mobility*force, -mx*force-mobility*fx

def continuum(n, k, tight=False):
    x, q = initial(n)
    def rhs(t,y):
        xx, qq = y[:n], y[n:]
        h = np.exp(-qq)*np.i0(qq*DELTA*np.exp(-4*PI**2*k*k*t))
        u,b = fields(xx,h)
        return np.concatenate((u,b*qq))
    return SampledSolution(rk4(rhs,np.concatenate((x,q)),ALL_TIMES,
        1e-5 if tight else 2e-5))

def particles(n,tight=False):
    count = n**3
    x0,q0 = initial(count)
    g = q0*(1+DELTA*np.cos(2*PI*np.arange(count)/n))
    g /= g.mean()
    x = np.concatenate(([0.],np.cumsum(g[:-1])/count))
    def rhs(t,xx):
        gap = count*np.diff(np.concatenate((xx,[xx[0]+1])))
        u,_ = fields(xx,np.exp(-gap))
        return u+(n*n/count)*(gap-np.roll(gap,1))
    return rk4(rhs,x,TIMES,1e-5 if tight else 2e-5)

def main(output_dir):
    ref=continuum(4096,1,True)
    ref_double=continuum(8192,1,True)
    vals=ref.sol(TIMES)
    doubled=ref_double.sol(TIMES)
    spatial=max(np.max(np.abs(vals[:4096]-doubled[:8192:2])),
                np.max(np.abs(vals[4096:]-doubled[8192::2])))
    rows=[]
    for n in [4,8,16]:
        count=n**3
        base=particles(n)
        xx=particles(n,True)
        gap=count*np.diff(np.vstack((xx,xx[:1]+1)),axis=0)
        basegap=count*np.diff(np.vstack((base,base[:1]+1)),axis=0)
        step=4096//count
        X,Q=vals[:4096:step],vals[4096::step]
        G=Q*(1+DELTA*np.cos(2*PI*np.arange(count)/n)[:,None]
            *np.exp(-4*PI**2*TIMES)[None,:])
        ee=np.max(np.abs(gap-G),axis=0)
        dd=np.max(np.abs(xx-X),axis=0)
        error=float(np.max(ee+dd))
        rows.append(dict(n=n,N=count,scale=count**(-2/3),
            max_mesh_E_plus_D=error,scaled_error=error*n*n,
            time_refinement_difference=float(np.max(np.abs(gap-basegap))+
                                            np.max(np.abs(xx-base)))))
    # Independently evaluate the analytic integral C_2 using Bessel identities.
    x,q=initial(8192)
    z=q*DELTA
    diss=(2*PI*q*DELTA)**2*np.exp(-q)*(np.i0(z)-i2(z))/2
    c,s=np.cos(2*PI*x),np.sin(2*PI*x)
    conv=diss.mean()+ETA*(c*(c*diss).mean()+s*(s*diss).mean())
    force=PI*ALPHA*s
    C2=float(3*np.mean(force**2*conv))
    ks=continuum(8192,2,True)
    times=OBS_TIMES
    x1=ref_double.sol(times)[:8192]
    x2=ks.sol(times)[:8192]
    difference=ALPHA*np.mean(np.sin(PI*(x1+x2))*np.sin(PI*(x2-x1)),axis=0)
    # Repeat with half label count and looser time tolerances.
    b1=continuum(4096,1)
    b2=continuum(4096,2)
    y1,y2=b1.sol(times)[:4096],b2.sol(times)[:4096]
    difference_base=ALPHA*np.mean(np.sin(PI*(y1+y2))*np.sin(PI*(y2-y1)),axis=0)
    result=dict(parameters=dict(alpha=ALPHA,delta=DELTA,eta=ETA,beta=0,
        m0=1,gamma=1,kappa=1,lambda_=1,T=.01,time_mesh_points=101),
        method='classical RK4 dt=1e-5 (refinement from 2e-5); analytic affine heat cell; periodic label quadrature',
        reference_spatial_refinement_max=float(spatial),particle_rows=rows,
        C2=C2,C2_over_2=C2/2,
        observable_rows=[dict(t=float(t),ratio=float(d/t**2),
            refinement_difference=float(abs(d-db)/t**2))
            for t,d,db in zip(times,difference,difference_base)])
    assert spatial<1e-9
    assert all(r['time_refinement_difference']<1e-7 for r in rows)
    assert abs(result['observable_rows'][-1]['ratio']-C2/2)/(C2/2)<0.01

    path = Path(output_dir) / 'data' / 'sanity_results.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2), flush=True)
    return result
