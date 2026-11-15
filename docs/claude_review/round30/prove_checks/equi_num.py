# numerical search: equidistant clusters (code classes from equi_codes.py), D, M, c=0
import json, sys, numpy as np, itertools
from scipy.optimize import least_squares
D=int(sys.argv[1]); M=int(sys.argv[2]); starts=int(sys.argv[3]) if len(sys.argv)>3 else 200
L=D-1-D//2
cls=json.load(open('equi_%d_%d.json'%(D,M)))
rng=np.random.default_rng(0)
def build(cl):
    cols=[]
    for t,n in cl: cols+= [t]*n
    S=np.array(cols).T  # (M-1) x D, entries +-1
    return (S==-1).astype(float)
def resid(x,F):
    rho=x[:D]+1j*x[D:]
    P=np.array([rho**j for j in range(1,L+1)]).imag  # L x D
    r=(F@P.T).ravel()
    # gauge: sum Re rho=0, sum |rho|^2 = D
    g=[x[:D].sum(), (np.abs(rho)**2).sum()-D]
    return np.concatenate([r,g])
def nondeg(x):
    rho=x[:D]+1j*x[D:]
    a=np.min(np.abs(rho.imag))
    b=min(min(abs(rho[k]-rho[l]),abs(rho[k]-np.conj(rho[l]))) for k in range(D) for l in range(k+1,D))
    return min(a,b)
for ci,cl in enumerate(cls):
    F=build(cl)
    best=None; good=[]
    for s in range(starts):
        x0=rng.normal(size=2*D)
        sol=least_squares(resid,x0,args=(F,),xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=4000)
        r=np.linalg.norm(sol.fun); nd=nondeg(sol.x)
        if r<1e-10 and nd>1e-3: good.append((nd,sol.x))
        if best is None or (r<1e-10 and nd>best[0]): best=(nd if r<1e-10 else -r, sol.x)
    if good:
        nd,x=max(good,key=lambda z:z[0])
        J=least_squares(resid,x,args=(F,)).jac
        sv=np.linalg.svd(J,compute_uv=False)
        rank=int((sv>1e-8*sv[0]).sum())
        print('class',ci,'FOUND',len(good),'nondeg',round(nd,4),'local dim ~',2*D-rank, flush=True)
    else:
        print('class',ci,'none; best',best[0],flush=True)
