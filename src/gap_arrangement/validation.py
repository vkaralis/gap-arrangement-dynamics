"""Independent equation, integration, transport and asymptotic checks.
Writes numerical audit_results.json and audit_curves.npz, not an editorial report.
"""
import os,tempfile
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('MPLCONFIGDIR',os.path.join(tempfile.gettempdir(),'gk-matplotlib'))
from pathlib import Path
import json,time
import numpy as np
from scipy.integrate import solve_ivp,quad
from scipy.optimize import linprog
from scipy.signal import resample
from . import simulation as base
P=2*np.pi
BETA=.5
TS=base.TIMES

def log(s):print(s,flush=True)

def init(n,k=1,m=None):
    a=np.arange(n)/n; x=a.copy()
    for _ in range(15):x-=(x+np.sin(P*x)/(3*P)-a)/(1+np.cos(P*x)/3)
    q=1/(1+np.cos(P*x)/3)
    if m is None:return x,q
    th=np.arange(m)/m
    return np.c_[x,q[:,None]*(1+.25*np.cos(P*k*th))]

def complex_fields(x,h):
    z=np.exp(1j*P*x);zbar=np.mean(z.conj());hbar=np.mean(h*z.conj())
    f=P*(z*zbar).imag;fx=P*P*(z*zbar).real
    mobility=1+np.mean(h)+.5*(z*hbar).real
    return -mobility*f,(.5*P*(z*hbar).imag)*f-mobility*fx

def dense_fields(x,h):
    d=P*(x[:,None]-x[None,:]);f=P*np.sin(d).mean(1)
    fx=P*P*np.cos(d).mean(1);m=1+((1+.5*np.cos(d))*h).mean(1)
    mx=(-.5*P*np.sin(d)*h).mean(1)
    return -m*f,-mx*f-m*fx

def solve(rhs,y,times=TS,rtol=2e-12,atol=2e-14,max_step=np.inf):
    # A diffusion-scale cap also controls stiff high-frequency dense-output modes.
    # Initial trial steps are capped to avoid overflow in rejected RK stages.
    shape=y.shape
    sol=solve_ivp(lambda t,v:rhs(t,v.reshape(shape)).ravel(),(0,float(times[-1])),y.ravel(),method='DOP853',t_eval=times,rtol=rtol,atol=atol,
        max_step=max_step,first_step=min(1e-6,max_step,float(times[-1])))
    assert sol.success and np.isfinite(sol.y).all(),sol.message
    return sol.y.T.reshape((len(times),)+shape)

def independent_particles(n,k,lam):
    N=n**3;x,q=init(N);g=q*(1+.25*np.cos(P*k*np.arange(N)/n));g/=g.mean()
    prefix=np.r_[0,np.cumsum(g[:-1])/N]
    y0=np.r_[g,prefix.mean()]
    def reconstruct(y):
        pre=np.r_[0,np.cumsum(y[:-2])/N]
        return pre+y[-1]-pre.mean()
    def rhs(t,y):
        gg=y[:-1];xx=reconstruct(y);u,_=complex_fields(xx,np.exp(-gg));a=gg+.5*gg**2
        return np.r_[N*(np.roll(u,-1)-u)+lam*n*n*(np.roll(a,-1)-2*a+np.roll(a,1)),u.mean()]
    yy=solve(rhs,y0,max_step=1/(4*lam*n*n) if lam else np.inf)
    xx=np.array([reconstruct(y) for y in yy])
    return xx,yy[:,:-1]

def cell_rhs(na,m,lam=1.,beta=.5,fd=False):
    omega=P*np.fft.fftfreq(m,d=1/m)
    def rhs(t,y):
        x,g=y[:,0],y[:,1:];u,b=complex_fields(x,np.exp(-g).mean(1));a=g+beta*g*g
        if fd:diff=m*m*(np.roll(a,-1,axis=1)-2*a+np.roll(a,1,axis=1))
        else:diff=np.fft.ifft(-omega**2*np.fft.fft(a,axis=1),axis=1).real
        return np.c_[u,b[:,None]*g+lam*diff]
    return rhs

def independent_cell(na,m,k=1,lam=1.,beta=.5,fd=False,times=TS,rtol=2e-12,atol=2e-14):
    return solve(cell_rhs(na,m,lam,beta,fd),init(na,k,m),times,rtol,atol,
        max_step=1/(lam*(np.pi*m)**2) if lam else np.inf)

def ref_w1(x1,x2,n=8192):
    a=np.arange(x1.shape[1])/x1.shape[1];aa=np.arange(n)/n
    y1=resample(x1-a,n,axis=1)+aa;y2=resample(x2-a,n,axis=1)+aa
    return np.array([base.circle_w1(a,b) for a,b in zip(y1,y2)])

def coefficient_quad(beta):
    rho=lambda x:1+np.cos(P*x)/3
    def D(x):
        q=1/rho(x)
        return quad(lambda th:np.exp(-q*(1+.25*np.cos(P*th)))*(1+2*beta*q*(1+.25*np.cos(P*th)))*(P*q*.25*np.sin(P*th))**2,0,1,epsabs=2e-12,epsrel=2e-12)[0]
    moments=[quad(lambda x:rho(x)*D(x)*(1 if j==0 else np.cos(P*x)),0,1,epsabs=2e-12,epsrel=2e-12)[0] for j in [0,1]]
    conv=lambda x:moments[0]+.5*np.cos(P*x)*moments[1]
    c=3*quad(lambda x:rho(x)*(np.pi/3*np.sin(P*x))**2*conv(x),0,1,epsabs=2e-12,epsrel=2e-12)[0]
    # For these symmetric data the circular tangent median is zero.
    # Half the integral of rho*|acceleration difference| reduces to this half-cell integral.
    cw=3*quad(lambda x:rho(x)*abs(np.pi/3*np.sin(P*x))*conv(x),0,.5,epsabs=2e-12,epsrel=2e-12)[0]
    return c,cw

def transport_checks():
    rng=np.random.default_rng(7713);errors=[];rot=[];sym=[]
    for it in range(60):
        n=2+it%8;m=2+(3*it)%9
        x,y=rng.random(n),rng.random(m)
        if it%5==0:x[:2]=0
        if it%7==0:y[0]=1
        d=abs(x[:,None]-y);cost=np.minimum(d,1-d)
        constraints=np.r_[np.kron(np.eye(n),np.ones((1,m))),np.kron(np.ones((1,n)),np.eye(m))]
        lp=linprog(cost.ravel(),A_eq=constraints,b_eq=np.r_[np.full(n,1/n),np.full(m,1/m)],bounds=(0,None),method='highs')
        assert lp.success
        w=base.circle_w1(x,y);errors.append(abs(w-lp.fun))
        z=rng.random();rot.append(abs(w-base.circle_w1(x+z,y+z)));sym.append(abs(w-base.circle_w1(y,x)))
    return dict(cases=60,max_LP_error=max(errors),max_rotation_error=max(rot),max_symmetry_error=max(sym))

def cell_diagnostics(y,lam):
    na=y.shape[1];m=y.shape[2]-1;g=y[:,:,1:];q=g.mean(2);r=g/q[:,:,None]
    ka=P*np.fft.fftfreq(na,d=1/na);kt=P*np.fft.fftfreq(m,d=1/m)
    dx=1+np.fft.ifft(1j*ka[None,:]*np.fft.fft(y[:,:,0]-np.arange(na)/na,axis=1),axis=1).real
    var=((r-1)**2).mean(2)
    bound=.25**2/2*np.exp(-8*np.pi**2*lam*TS[:,None])
    lawerr=[];varerr=[]
    for ti in [0,len(TS)//2,-1]:
        x=y[ti,:,0];gg=g[ti];qq=q[ti];rr=r[ti]
        gt=cell_rhs(na,m,lam)(TS[ti],y[ti])[:,1:]
        _,b=complex_fields(x,np.exp(-gg).mean(1))
        grad=np.fft.ifft(1j*kt*np.fft.fft(gg,axis=1),axis=1).real
        lhs=(-np.exp(-gg)*gt).mean(1)
        rhs=-b*(gg*np.exp(-gg)).mean(1)-lam*(np.exp(-gg)*(1+gg)*grad**2).mean(1)
        lawerr.append(float(np.max(abs(lhs-rhs))))
        rt=(gt-b[:,None]*gg)/qq[:,None]
        vlhs=2*((rr-1)*rt).mean(1);vrhs=-2*lam*((1+gg)*(grad/qq[:,None])**2).mean(1)
        varerr.append(float(np.max(abs(vlhs-vrhs))))
    return dict(min_G=float(g.min()),max_mean_mass_error=float(np.max(abs(q.mean(1)-1))),compatibility_error=float(np.max(abs(dx-q))),max_law_identity_residual=max(lawerr),max_variance_identity_residual=max(varerr),max_variance_bound_violation=float(np.max(var-bound)),max_variance_increase=float(np.max(np.diff(var,axis=0))))

def main(output_dir):
    OUT = Path(output_dir) / "data"
    OUT.mkdir(parents=True, exist_ok=True)
    start=time.monotonic();result={};curves={'times':TS}
    log('Independent fields and circular transport')
    x,q=init(47);h=np.exp(-q*(1+.1*np.cos(P*np.arange(47)/47)))
    direct=dense_fields(x,h);production=base.fields(x,h);complex_=complex_fields(x,h)
    result['fields']=dict(max_moment_vs_dense=float(np.max(abs(np.array(direct)-production))),max_complex_vs_dense=float(np.max(abs(np.array(direct)-complex_))))
    result['transport']=transport_checks()
    log('Independent physical-coordinate adaptive quadrature')
    result['coefficients']={}
    for beta in [0.,.5]:
        c,cw=coefficient_quad(beta)
        result['coefficients'][str(beta)]=dict(C2=c,C2_over_2=c/2,W1_leading_coefficient=cw,error_vs_label_quadrature=abs(c-base.coefficient(beta)))
    log('Independent conservative gap ODE with DOP853')
    prows=[]
    for n in [8,16,32]:
        for lam in [0.,1.]:
            for k in [1,2]:
                log(f'Particle audit N={n**3}, lambda={lam}, k={k}')
                x,g=independent_particles(n,k,lam);ref=base.particles(n,k,lam)
                gref=n**3*np.diff(np.c_[ref,ref[:,:1]+1],axis=1)
                prows.append(dict(N=n**3,k=k,lambda_=lam,max_X_difference=float(np.max(abs(x-ref))),max_G_difference=float(np.max(abs(g-gref))),max_mass_error=float(np.max(abs(g.mean(1)-1))),min_G=float(g.min())))
    result['independent_particles']=prows
    log('Cell solver: independent DOP853, separate spatial and time checks')
    cells={}
    for k in [1,2]:
        y=independent_cell(64,64,k);cells[k]=y
        production=base.cell(64,64,k)
        result[f'cell_solver_k{k}']=dict(max_X_difference=float(np.max(abs(y[:,:,0]-production[:,:,0]))),max_G_difference=float(np.max(abs(y[:,:,1:]-production[:,:,1:]))))
        result[f'cell_diagnostics_k{k}']=cell_diagnostics(y,1.)
    ref=cells[2]
    densea=independent_cell(128,64,2)
    denseth=independent_cell(64,128,2)
    loose=independent_cell(64,64,2,rtol=2e-10,atol=2e-12)
    result['separate_refinement']=dict(label_error=float(np.max(abs(densea[:,::2]-ref))),cell_error=float(np.max(abs(np.concatenate([denseth[:,:,:1],denseth[:,:,1::2]],axis=2)-ref))),time_tolerance_error=float(np.max(abs(loose-ref))))
    log('Independent centered finite-difference cell discretization')
    fdrows=[];ref32=cells[2][:,::2]
    for m in [16,32,64,128]:
        y=independent_cell(32,m,2,fd=True)
        gx=resample(ref32[:,:,1:],m,axis=2)
        fdrows.append(dict(cell_points=m,max_X_error=float(np.max(abs(y[:,:,0]-ref32[:,:,0]))),max_G_error=float(np.max(abs(y[:,:,1:]-gx)))))
    result['finite_difference']=fdrows
    log('Affine exact solution and genuine zero-diffusion cell tests')
    aff=independent_cell(64,64,2,beta=0);aexact=base.affine(64,2)
    expected=aexact[:,:,1,None]*(1+.25*np.exp(-4*np.pi**2*4*TS[:,None,None])*np.cos(2*P*np.arange(64)[None,None,:]/64))
    result['affine_formula']=dict(max_X_error=float(np.max(abs(aff[:,:,0]-aexact[:,:,0]))),max_G_error=float(np.max(abs(aff[:,:,1:]-expected))))
    z1=independent_cell(64,65,1,lam=0);z2=independent_cell(64,65,2,lam=0)
    result['zero_cell']=dict(max_X_difference=float(np.max(abs(z1[:,:,0]-z2[:,:,0]))),max_mean_G_difference=float(np.max(abs(z1[:,:,1:].mean(2)-z2[:,:,1:].mean(2)))))
    log('Early-time coefficients and W1 quadrature')
    small=np.array([0.,5e-6,1e-5,2e-5,5e-5,1e-4,2e-4,5e-4,1e-3])
    sy=[independent_cell(128,32,k,times=small) for k in [1,2]]
    J=base.observable(sy[0][:,:,0],sy[1][:,:,0]);w=ref_w1(sy[0][:,:,0],sy[1][:,:,0])
    result['early_time']=[dict(t=float(t),J_over_t2=float(j/t**2),W1_over_t2=float(v/t**2)) for t,j,v in zip(small[1:],J[1:],w[1:])]
    w512=ref_w1(cells[1][:,:,0],cells[2][:,:,0],512);w8192=ref_w1(cells[1][:,:,0],cells[2][:,:,0],8192);w16384=ref_w1(cells[1][:,:,0],cells[2][:,:,0],16384)
    result['W1_reference']=dict(max_512_to_8192=float(np.max(abs(w512-w8192))),max_8192_to_16384=float(np.max(abs(w8192-w16384))),final_W1_8192=float(w8192[-1]))
    log('Full trajectory comparison with independently evolved reduced model')
    def reduced(na):
        x,q=init(na)
        def rhs(t,y):
            u,b=complex_fields(y[:,0],np.exp(-y[:,1]));return np.c_[u,b*y[:,1]]
        return solve(rhs,np.c_[x,q])
    red=reduced(64);lrows=[]
    for lam in [1.,2.,4.,8.,16.,32.]:
        log(f'Relaxation lambda={lam}')
        y=independent_cell(64,32,1,lam=lam)
        e=np.max(abs(y[:,:,0]-red[:,:,0]),axis=1)+np.max(abs(y[:,:,1:].mean(2)-red[:,:,1]),axis=1)
        q=y[:,:,1:].mean(2);v=((y[:,:,1:]/q[:,:,None]-1)**2).mean(2).max(1)
        lrows.append(dict(lambda_=lam,max_trajectory_error=float(e.max()),scaled_error=float(lam*e.max()),final_error=float(e[-1])))
        curves[f'relaxation_error_{int(lam)}']=e;curves[f'variance_{int(lam)}']=v
    tight=independent_cell(64,64,1,lam=32)
    e=np.max(abs(tight[:,:,0]-red[:,:,0]),axis=1)+np.max(abs(tight[:,:,1:].mean(2)-red[:,:,1]),axis=1)
    result['relaxation']=lrows;result['relaxation_refinement_error']=float(np.max(abs(e-curves['relaxation_error_32'])))
    result['elapsed_seconds']=time.monotonic()-start
    validate_results(result)
    (OUT/'audit_results.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savez_compressed(OUT/'audit_curves.npz',**curves)
    log(json.dumps(result,indent=2))
    return result

def validate_results(r):
    assert r['fields']['max_moment_vs_dense']<1e-12
    assert max(r['transport'][k] for k in ['max_LP_error','max_rotation_error','max_symmetry_error'])<1e-11
    assert all(v['error_vs_label_quadrature']<1e-10 for v in r['coefficients'].values())
    assert all(p['max_X_difference']<1e-9 and p['max_G_difference']<1e-7 and p['min_G']>0 and p['max_mass_error']<1e-12 for p in r['independent_particles'])
    for k in [1,2]:
        assert r[f'cell_solver_k{k}']['max_G_difference']<1e-8
        d=r[f'cell_diagnostics_k{k}']
        assert d['min_G']>0 and d['compatibility_error']<1e-8 and d['max_mean_mass_error']<1e-10
        assert d['max_law_identity_residual']<1e-8 and d['max_variance_identity_residual']<1e-8
        assert d['max_variance_bound_violation']<1e-12 and d['max_variance_increase']<1e-12
    assert max(r['separate_refinement'].values())<1e-7
    assert all(r['finite_difference'][i]['max_G_error']/r['finite_difference'][i+1]['max_G_error']>3.8 for i in range(3))
    assert max(r['affine_formula'].values())<1e-8
    assert max(r['zero_cell'].values())<1e-11
    assert r['W1_reference']['max_8192_to_16384']<1e-10
    assert r['relaxation_refinement_error']<1e-8
    rr=r['relaxation']
    assert all(rr[i]['max_trajectory_error']>rr[i+1]['max_trajectory_error']>0 for i in range(len(rr)-1))
    assert 1.9<rr[-2]['max_trajectory_error']/rr[-1]['max_trajectory_error']<2.1
    assert abs(r['early_time'][0]['J_over_t2']/r['coefficients']['0.5']['C2_over_2']-1)<.002
