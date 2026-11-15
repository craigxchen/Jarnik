# Numerical realisation of cluster types in the FF model (c arbitrary).
# Type = list of sign vectors eps_i in {+-1}^D, normalised: +1 means the factor (t-rho_k), Im rho_k>0.
# Equations: sum_k (eps_i-eps_1)_k/2 * Im(rho_k^m) = 0, m=1..L, i=2..M.
import numpy as np
from scipy.optimize import least_squares
import itertools
def LofD(D,c): return D-1-(D+2*c)//2
def residual(z, E, L):
    D=E.shape[1]; x=z[:D]; u=z[D:]; y=np.exp(u)          # Im rho = e^u > 0
    rho=x+1j*y
    P=np.array([np.imag(rho**m) for m in range(1,L+1)])   # L x D
    Dl=(E[1:]-E[0])/2.0                                   # (M-1) x D
    res=(Dl@P.T).ravel()
    # gauge: sum x = 0, sum log y = 0  (real affine maps: translation and scaling)
    g=[x.sum(), u.sum()]
    return np.concatenate([res, g])
def degeneracy(rho):
    D=len(rho); allr=np.concatenate([rho,np.conj(rho)])
    dmin=min(abs(allr[i]-allr[j]) for i in range(2*D) for j in range(i+1,2*D))
    return dmin/ max(1,np.max(abs(rho)))
def solve(E, L, starts=30, seed=0, tol=1e-11, scale=1.0):
    E=np.array(E,float); D=E.shape[1]; rng=np.random.default_rng(seed); best=None
    for s in range(starts):
        z0=np.concatenate([rng.normal(0,scale,D), rng.normal(0,0.7,D)])
        try:
            r=least_squares(residual,z0,args=(E,L),method='lm',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=4000)
        except Exception: continue
        rho=r.x[:D]+1j*np.exp(r.x[D:])
        cost=np.linalg.norm(r.fun)
        dg=degeneracy(rho)
        if best is None or (cost<tol and dg>1e-4 and (best[0]>=tol or dg>best[2])) or (best[0]>=tol and cost<best[0]):
            best=(cost,rho,dg)
        if cost<tol and dg>1e-2: break
    return best
def cluster_through(rho, L, base, tol=1e-8):
    """all sign vectors e with sum (e-base)/2 Im rho^m = 0; brute force 2^D"""
    D=len(rho); P=np.array([np.imag(rho**m) for m in range(1,L+1)])
    out=[]
    for bits in itertools.product([1,-1],repeat=D):
        e=np.array(bits); v=P@((e-base)/2)
        if np.linalg.norm(v)<tol*max(1,np.abs(P).max()): out.append(e)
    return np.array(out)
def rank_of(E):
    E=np.array(E,float); return np.linalg.matrix_rank(E[1:]-E[0]) if len(E)>1 else 0
