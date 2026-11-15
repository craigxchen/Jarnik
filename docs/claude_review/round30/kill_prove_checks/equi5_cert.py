# Independent: find and certify (Krawczyk) a nondegenerate 5-point equidistant (rank-4) cluster at D=8, c=0.
import itertools, numpy as np, mpmath as mp, sys, pickle
from scipy.optimize import least_squares
exec(open('equi_lattice.py').read().split("# short integer vectors")[0].replace("M=int(sys.argv[1])","M=5"))
D=8;L=3
keys=list(classes.keys())
rng=np.random.default_rng(int(sys.argv[1]) if len(sys.argv)>1 else 0)
def F(x,dl):
    rho=x[:D]+1j*x[D:]
    P=np.array([rho**m for m in range(1,L+1)]).imag
    return (dl@P.T).ravel()
def nd(x):
    rho=x[:D]+1j*x[D:]; allr=np.concatenate([rho,rho.conj()])
    return min(np.abs(allr[i]-allr[j]) for i in range(2*D) for j in range(i+1,2*D))
found=None
for ci,key in enumerate(keys):
    c=classes[key]; pts=np.array([ones]+[W[j] for j in c]); dl=(pts[1:]-pts[0])//2
    best=None
    for st in range(300):
        x0=rng.normal(size=2*D)
        def res(z):
            x=np.concatenate([[0.0],z[:D-1],[1.0],z[D-1:]])
            return F(x,dl)
        sol=least_squares(res,x0[:2*D-2],xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=3000)
        x=np.concatenate([[0.0],sol.x[:D-1],[1.0],sol.x[D-1:]])
        if np.linalg.norm(sol.fun)<1e-11:
            q=nd(x)/np.sqrt((x**2).sum())
            if best is None or q>best[0]: best=(q,x)
    print('class',ci,'best normalised separation',None if best is None else best[0],flush=True)
    if best and best[0]>1e-3 and found is None: found=(ci,dl,pts,best[1])
pickle.dump(found,open('equi5_found.pkl','wb'))
