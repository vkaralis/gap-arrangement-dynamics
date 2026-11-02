"""Deterministic numerical illustrations of gap arrangement; no external data.
Run from any directory. Numerical refinement is not a proof or interval bound.
"""
import os, tempfile
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('MPLCONFIGDIR',os.path.join(tempfile.gettempdir(),'gk-matplotlib'))
from pathlib import Path
import json, sys, time
import numpy as np
import scipy
from scipy.interpolate import CubicSpline
from scipy.optimize import linprog
from scipy.signal import resample
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PI=np.pi
ALPHA,DELTA,ETA,BETA=1/3,1/4,1/2,1/2
T=.02
GRID=np.linspace(0,T,101)
SMALL=np.array([.001,.0005,.0002,.0001,.00005])
TIMES=np.unique(np.r_[GRID,SMALL])

def initial(n):
    a=np.arange(n)/n
    x=a.copy()
    for _ in range(12):
        x-=(x+ALPHA*np.sin(2*PI*x)/(2*PI)-a)/(1+ALPHA*np.cos(2*PI*x))
    return x,1/(1+ALPHA*np.cos(2*PI*x))

def fields(x,h):
    c,s=np.cos(2*PI*x),np.sin(2*PI*x)
    c0,s0=c.mean(),s.mean()
    ch,sh=(c*h).mean(),(s*h).mean()
    f=2*PI*(s*c0-c*s0)
    fx=(2*PI)**2*(c*c0+s*s0)
    m=1+h.mean()+ETA*(c*ch+s*sh)
    mx=2*PI*ETA*(-s*ch+c*sh)
    return -m*f,-mx*f-m*fx

def rk4(rhs,y,dt,times=TIMES):
    t=0.; saved=[]
    for end in times:
        count=int(np.ceil(max(0.,end-t)/dt-1e-9))
        if count:
            h=(end-t)/count
            for _ in range(count):
                k1=rhs(t,y); k2=rhs(t+h/2,y+h*k1/2)
                k3=rhs(t+h/2,y+h*k2/2); k4=rhs(t+h,y+h*k3)
                y=y+h/6*(k1+2*k2+2*k3+k4); t+=h
        t=end; saved.append(y.copy())
    return np.array(saved)

def cell(na,ntheta,k,lam=1.,beta=BETA,dt=5e-6):
    x,q=initial(na); theta=np.arange(ntheta)/ntheta
    g=q[:,None]*(1+DELTA*np.cos(2*PI*k*theta))
    omega2=(2*PI*np.fft.rfftfreq(ntheta,d=1/ntheta))**2
    def rhs(t,y):
        xx,gg=y[:,0],y[:,1:]
        u,b=fields(xx,np.exp(-gg).mean(axis=1))
        diffusion=np.fft.irfft(-omega2*np.fft.rfft(gg+beta*gg*gg,axis=1),n=ntheta,axis=1)
        return np.c_[u,b[:,None]*gg+lam*diffusion]
    return rk4(rhs,np.c_[x,g],dt)

def affine(na,k,lam=1.,dt=1e-5):
    x,q=initial(na)
    def rhs(t,y):
        xx,qq=y[:,0],y[:,1]
        h=np.exp(-qq)*np.i0(qq*DELTA*np.exp(-4*PI**2*k*k*lam*t))
        u,b=fields(xx,h)
        return np.c_[u,b*qq]
    return rk4(rhs,np.c_[x,q],dt)

def particles(n,k,lam=1.,dt=1e-5):
    N=n**3; x,q=initial(N)
    g=q*(1+DELTA*np.cos(2*PI*k*np.arange(N)/n)); g/=g.mean()
    x=np.r_[0,np.cumsum(g[:-1])/N]
    def rhs(t,xx):
        gap=N*np.diff(np.r_[xx,xx[0]+1])
        u,_=fields(xx,np.exp(-gap))
        ag=gap+BETA*gap*gap
        return u+(lam*n*n/N)*(ag-np.roll(ag,1))
    return rk4(rhs,x,dt)

def circle_w1(x,y):
    """Exact W1 for equal-weight atomic measures and circular geodesic cost.
    Integral |CDF_x-CDF_y-c| minimized at a length-weighted median c.
    """
    x=np.mod(x,1); y=np.mod(y,1)
    points=np.r_[0,x,y,1]
    jumps=np.r_[0,np.full(len(x),1/len(x)),np.full(len(y),-1/len(y)),0]
    order=np.argsort(points,kind='stable'); points=points[order]
    levels=np.cumsum(jumps[order])[:-1]; lengths=np.diff(points)
    order=np.argsort(levels); pos=np.searchsorted(np.cumsum(lengths[order]),.5)
    median=levels[order[min(pos,len(order)-1)]]
    return float(np.dot(lengths,np.abs(levels-median)))

def reference_w1(x1,x2,points=8192):
    """Refine smooth periodic displacement maps before atomic W1 quadrature."""
    a=np.arange(x1.shape[1])/x1.shape[1]; target=np.arange(points)/points
    y1=resample(x1-a,points,axis=1)+target
    y2=resample(x2-a,points,axis=1)+target
    return np.array([circle_w1(u,v) for u,v in zip(y1,y2)])

def verify_w1():
    rng=np.random.default_rng(20261102)
    for n in [3,5,7]:
        x,y=rng.random(n),rng.random(n)
        d=np.abs(x[:,None]-y); costs=np.minimum(d,1-d)
        constraints=np.r_[np.kron(np.eye(n),np.ones((1,n))),np.kron(np.ones((1,n)),np.eye(n))]
        lp=linprog(costs.ravel(),A_eq=constraints,b_eq=np.full(2*n,1/n),bounds=(0,None),method='highs')
        assert lp.success and abs(circle_w1(x,y)-lp.fun)<1e-11
    assert abs(circle_w1(np.array([.99]),np.array([.01]))-.02)<1e-14

def sample_reference(state,N,p):
    na,cols=state.shape; m=cols-1; a=np.arange(na)/na; target=np.arange(N)/N
    displacement=state[:,0]-a
    xs=CubicSpline(np.r_[a,1],np.r_[displacement,displacement[0]],bc_type='periodic')(target)+target
    gs=np.empty(N)
    assert m%p==0
    for j in range(p):
        v=state[:,1+j*(m//p)]
        gs[j::p]=CubicSpline(np.r_[a,1],np.r_[v,v[0]],bc_type='periodic')(target[j::p])
    return xs,gs

def coefficient(beta,na=2048,m=256):
    x,q=initial(na); th=np.arange(m)/m
    g=q[:,None]*(1+DELTA*np.cos(2*PI*th))
    gt=-2*PI*q[:,None]*DELTA*np.sin(2*PI*th)
    diss=(np.exp(-g)*(1+2*beta*g)*gt*gt).mean(axis=1)
    c,s=np.cos(2*PI*x),np.sin(2*PI*x)
    conv=diss.mean()+ETA*(c*(c*diss).mean()+s*(s*diss).mean())
    return float(3*np.mean((PI*ALPHA*s)**2*conv))

def observable(x1,x2):
    # Stable trig difference avoids subtraction of nearly equal cosines.
    return ALPHA*np.mean(np.sin(PI*(x1+x2))*np.sin(PI*(x2-x1)),axis=1)

def log(msg): print(msg,flush=True)

def main(output_dir):
    output_dir = Path(output_dir)
    OUT = output_dir / "data"
    OUT.mkdir(parents=True, exist_ok=True)
    require_audit(OUT)
    start=time.monotonic(); verify_w1()
    refs={}; refinement={}; coarse_maps={}
    for k in [1,2]:
        log(f'Nonlinear cell k={k}: 512 labels, 64 cell points')
        refs[k]=cell(512,64,k)
        log(f'Coarser cell k={k}: 256 labels, 32 cell points')
        coarse=cell(256,32,k,dt=2e-5)
        coarse_maps[k]=coarse[:,:,0]
        fine_sub=np.concatenate([refs[k][:,::2,:1],refs[k][:,::2,1::2]],axis=2)
        refinement[str(k)]=dict(max_X=float(np.max(abs(coarse[:,:,0]-fine_sub[:,:,0]))),max_G=float(np.max(abs(coarse[:,:,1:]-fine_sub[:,:,1:]))))
    x1,x2=refs[1][:,:,0],refs[2][:,:,0]
    W=reference_w1(x1,x2)
    obs=observable(x1,x2)
    coarse_w=reference_w1(coarse_maps[1],coarse_maps[2])
    w1_cell_refinement=float(np.max(abs(W-coarse_w)))
    w1_quadrature_refinement=float(np.max(abs(W-reference_w1(x1,x2,16384))))
    small_ids=np.searchsorted(TIMES,SMALL)
    ratio_refinement=float(np.max(abs(observable(coarse_maps[1],coarse_maps[2])[small_ids]-obs[small_ids])/SMALL**2))
    aff1,aff2=affine(1024,1),affine(1024,2)
    obs_aff=observable(aff1[:,:,0],aff2[:,:,0])
    # Independent zero-diffusion implementation uses the exact frozen-cell law.
    z1,z2=affine(1024,1,lam=0),affine(1024,2,lam=0)
    zero=float(np.max(abs(z1-z2)))
    rows=[]; curves={}; zero_curves={}; arrays={'times':TIMES,'cell_w1':W,'cell_observable':obs,'affine_observable':obs_aff}
    grid_ids=np.searchsorted(TIMES,GRID)
    for n in [8,16,32]:
        log(f'Particles p={n}, N={n**3}, lambda=1, k=1,2')
        xx1,xx2=particles(n,1),particles(n,2)
        w=np.array([circle_w1(a,b) for a,b in zip(xx1,xx2)])
        curves[n]=w; arrays[f'w1_n{n}']=w
        max_error=0.; min_gap=float('inf'); max_gap=0.
        for ti in grid_ids:
            xr,gr=sample_reference(refs[1][ti],n**3,n)
            gap=n**3*np.diff(np.r_[xx1[ti],xx1[ti,0]+1])
            max_error=max(max_error,float(np.max(abs(gap-gr))+np.max(abs(xx1[ti]-xr))))
            min_gap=min(min_gap,float(gap.min()));max_gap=max(max_gap,float(gap.max()))
        log(f'Particle step refinement p={n}')
        tight1=particles(n,1,dt=5e-6); tight2=particles(n,2,dt=5e-6)
        tg=n**3*np.diff(np.c_[tight1,tight1[:,:1]+1],axis=1)
        bg=n**3*np.diff(np.c_[xx1,xx1[:,:1]+1],axis=1)
        dt_error=float(np.max(abs(tg-bg))+np.max(abs(tight1-xx1)))
        w_tight=np.array([circle_w1(a,b) for a,b in zip(tight1,tight2)])
        log(f'Zero-diffusion particles p={n}')
        zz1,zz2=particles(n,1,lam=0),particles(n,2,lam=0)
        zw=np.array([circle_w1(a,b) for a,b in zip(zz1,zz2)])
        zero_curves[n]=zw; arrays[f'zero_w1_n{n}']=zw
        rows.append(dict(p=n,N=n**3,scale=n**-2,max_mesh_E_plus_D=max_error,scaled_error=max_error*n*n,
            time_refinement_E_plus_D=dt_error,time_refinement_W1=float(np.max(abs(w_tight-w))),
            W1_initial=float(w[0]),W1_final=float(w[-1]),W1_final_error=float(abs(w[-1]-W[-1])),
            zero_W1_initial=float(zw[0]),zero_W1_final=float(zw[-1]),min_gap=min_gap,max_gap=max_gap))
        del xx1,xx2,tight1,tight2,tg,bg,zz1,zz2
    ratios=[]
    c2=coefficient(BETA); c2_aff=coefficient(0)
    for t in SMALL:
        i=np.searchsorted(TIMES,t)
        ratios.append(dict(t=float(t),affine_J_over_t2=float(obs_aff[i]/t**2),nonlinear_J_over_t2=float(obs[i]/t**2),nonlinear_W1_over_t2=float(W[i]/t**2)))
    result=dict(parameters=dict(alpha=ALPHA,delta=DELTA,eta=ETA,beta=BETA,lambda_=1,T=T),
        environment=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__),
        method='RK4; Fourier collocation in theta; periodic trapezoidal a quadrature; cubic periodic reference sampling; exact circular atomic W1; continuum displacement maps resampled at 8192 labels',
        cell_refinement=refinement,W1_reference_refinement=w1_quadrature_refinement,W1_cell_solver_refinement=w1_cell_refinement,W1_reference_labels=8192,observable_ratio_refinement=ratio_refinement,zero_continuum_difference=zero,C2=c2,C2_over_2=c2/2,affine_C2_over_2=c2_aff/2,
        final_cell_W1=float(W[-1]),particle_rows=rows,observable_rows=ratios,elapsed_seconds=time.monotonic()-start)
    validate_results(result)
    (OUT/'simulation_results.json').write_text(json.dumps(result,indent=2)+'\n', encoding='utf-8')
    np.savez_compressed(OUT/'simulation_curves.npz',**arrays)
    make_outputs(result, arrays, OUT, output_dir)
    log(json.dumps(result,indent=2))
    return result


def validate_results(result):
    """Acceptance thresholds for the fixed experiment settings."""
    assert result['zero_continuum_difference'] < 1e-13
    assert all(r['min_gap'] > 0 and r['time_refinement_E_plus_D'] < 1e-7
               for r in result['particle_rows'])
    assert all(v['max_G'] < 1e-6 for v in result['cell_refinement'].values())
    assert abs(result['observable_rows'][-1]['nonlinear_J_over_t2'] /
               result['C2_over_2'] - 1) < .02
    assert result['W1_reference_refinement'] < 1e-10
    assert result['W1_cell_solver_refinement'] < 1e-10


def require_audit(data_dir):
    from .validation import validate_results as validate_audit
    for name in ['audit_results.json', 'audit_curves.npz']:
        if not (data_dir / name).is_file():
            raise FileNotFoundError(
                f'Missing {name}: run "validate" first, or use "reproduce".')
    validate_audit(json.loads((data_dir / 'audit_results.json').read_text()))


def make_outputs(result, a, data_dir, output_dir):
    """Render both figures and tables from saved results, without integration."""
    from .exports import save_figure
    OUT = Path(data_dir)
    output_dir = Path(output_dir)
    require_audit(OUT)
    (output_dir / "figures").mkdir(parents=True, exist_ok=True)
    (output_dir / "tables").mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Serif','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':160,'savefig.dpi':300,'mathtext.fontset':'dejavuserif'})
    colors=['#0072B2','#D55E00','#009E73']
    fig,axes=plt.subplots(2,2,figsize=(7.05,5.4),layout='constrained')
    ax=axes.ravel()
    ns=np.array([r['N'] for r in result['particle_rows']]);err=np.array([r['max_mesh_E_plus_D'] for r in result['particle_rows']])
    ax[0].loglog(ns,err,'o-',color=colors[0],label='Nonlinear particle error')
    ax[0].loglog(ns,err[-1]*(ns/ns[-1])**(-2/3),'--',color='.35',label=r'$N^{-2/3}$ guide')
    ax[0].set(xlabel=r'$N$',ylabel=r'$\max_{t_j}(E_N+D_N)$',title='(a) Particle approximation');ax[0].legend(fontsize=8)
    ids=np.searchsorted(a['times'],SMALL)
    ax[1].semilogx(SMALL,a['affine_observable'][ids]/SMALL**2,'o-',color=colors[0],label=r'$\beta=0$')
    ax[1].semilogx(SMALL,a['cell_observable'][ids]/SMALL**2,'s-',color=colors[1],label=r'$\beta=1/2$')
    ax[1].axhline(result['affine_C2_over_2'],color=colors[0],ls=':',lw=1)
    ax[1].axhline(result['C2_over_2'],color=colors[1],ls=':',lw=1)
    ax[1].set(xlabel=r'$t$',ylabel=r'$(J_2-J_1)/t^2$',title='(b) Second-order response');ax[1].legend(fontsize=8)
    audit=json.loads((OUT/'audit_results.json').read_text())
    with np.load(OUT/'audit_curves.npz') as curves:
        lam=np.array([r['lambda_'] for r in audit['relaxation']])
        error=np.array([r['max_trajectory_error'] for r in audit['relaxation']])
        ax[2].loglog(lam,error,'o-',color=colors[0],label='Full vs. reduced trajectory')
        ax[2].loglog(lam,error[-1]*lam[-1]/lam,'--',color='.35',label=r'$\lambda^{-1}$ guide')
        ax[2].set(xlabel=r'$\lambda$',ylabel=r'$\max_{t_j}e_\lambda(t_j)$',title='(c) Fast-relaxation error')
        ax[2].legend(fontsize=7.5)
        for ll,cc in [(1,colors[0]),(4,colors[1])]:
            ax[3].semilogy(curves['times'],curves[f'variance_{ll}'],color=cc,label=rf'$\lambda={ll}$ computed')
            ax[3].semilogy(curves['times'],DELTA**2/2*np.exp(-8*PI**2*ll*curves['times']),ls=':',color=cc,label=rf'$\lambda={ll}$ bound')
        ax[3].set(xlabel=r'$t$',ylabel=r'$\max_a V(t,a)$',title='(d) Cell variance decay')
        ax[3].legend(fontsize=7.5)
    save_figure(fig, output_dir / 'figures/figure1_approximation')
    plt.close(fig)
    fig,ax=plt.subplots(1,2,figsize=(7.05,3.1),layout='constrained')
    for c,n in zip(colors,[8,16,32]):
        ax[0].plot(a['times'],a[f'w1_n{n}'],color=c,label=f'$N={n**3}$')
        ax[1].plot(a['times'],a[f'zero_w1_n{n}'],color=c,label=f'$N={n**3}$')
    ax[0].plot(a['times'],a['cell_w1'],'k--',lw=1.6,label='Cell reference')
    ax[1].axhline(0,color='k',ls='--',label='Exact cell limit')
    for j in range(2):
        ax[j].set(xlabel=r'$t$',ylabel=r'$W_1$ on the circle')
        ax[j].ticklabel_format(axis='y',style='sci',scilimits=(0,0));ax[j].legend(fontsize=7.5)
    inset=ax[0].inset_axes([.40,.38,.56,.36])
    for c,n in zip(colors,[8,16,32]): inset.plot(a['times'],a[f'w1_n{n}'],color=c,lw=1)
    inset.plot(a['times'],a['cell_w1'],'k--',lw=1.2)
    inset.set(xlim=(.008,.020),ylim=(0,6e-5));inset.tick_params(labelsize=6)
    inset.ticklabel_format(axis='y',style='sci',scilimits=(0,0));inset.yaxis.get_offset_text().set_fontsize(6)
    ax[0].set_title(r'(a) Nonlinear diffusion, $\lambda=1$')
    ax[1].set_title(r'(b) Control, $\lambda=0$')
    save_figure(fig, output_dir / 'figures/figure2_separation')
    plt.close(fig)
    rows=result['particle_rows']
    table=r'''\begin{table}[htbp]
\centering\small
\caption{Nonlinear diffusion $A(s)=s+s^2/2$: particle approximation and circular separation at $T=0.02$. The maximum uses 101 output times. The last column is the zero-diffusion control.}\label{tab:nonlinear}
\begin{tabular}{@{}rrrrr@{}}\toprule
$N$ & $\max(E_N+D_N)$ & $N^{2/3}\max(E_N+D_N)$ & $W_{1,N}(T)$ & $W_{1,N}^{\lambda=0}(T)$ \\\midrule
'''
    for r in rows:
        table+=f"{r['N']} & {r['max_mesh_E_plus_D']:.6f} & {r['scaled_error']:.4f} & {r['W1_final']:.7f} & {r['zero_W1_final']:.7f} \\\\\n"
    table+='\\bottomrule\\end{tabular}\n\\end{table}\n'
    (output_dir/'tables/table1.tex').write_text(table)
    table=r'''\begin{table}[htbp]
\centering\small
\caption{Small-time response for affine and nonlinear diffusion, with $\lambda=1$. The last column uses the nonlinear cell reference. The $t=0$ row gives the independently integrated theoretical observable coefficients, not a time-step extrapolation.}\label{tab:response}
\begin{tabular}{@{}rrrr@{}}\toprule
$t$ & $(J_2-J_1)/t^2$, $\beta=0$ & $(J_2-J_1)/t^2$, $\beta=1/2$ & $W_1/t^2$, $\beta=1/2$ \\\midrule
'''
    for r in result['observable_rows']:
        table+=f"{r['t']:.5f} & {r['affine_J_over_t2']:.6f} & {r['nonlinear_J_over_t2']:.6f} & {r['nonlinear_W1_over_t2']:.6f} \\\\\n"
    table+=f"$0$ (coefficient) & {result['affine_C2_over_2']:.6f} & {result['C2_over_2']:.6f} & --- \\\\\n"
    table+='\\bottomrule\\end{tabular}\n\\end{table}\n'
    (output_dir/'tables/table2.tex').write_text(table)

    # CSV exports retain the full numerical precision behind the displayed tables.
    import csv
    for name, rows in [('table1_particles.csv', result['particle_rows']),
                       ('table2_response.csv', result['observable_rows'])]:
        with (output_dir / 'tables' / name).open('w', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    with (output_dir / 'tables' / 'theoretical_coefficients.csv').open('w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle)
        writer.writerow(['beta', 'C2_over_2'])
        writer.writerow([0, result['affine_C2_over_2']])
        writer.writerow([BETA, result['C2_over_2']])
