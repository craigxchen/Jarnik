# s-form with general index set K of real e_k, sign-pattern parametrisation
import itertools, numpy as np, sys
from scipy.optimize import least_squares
SUBS=[frozenset(c) for k in range(1,5) for c in itertools.combinations(range(1,5),k)]
IDX={T:i for i,T in enumerate(SUBS)}; MEM=[[IDX[T] for T in SUBS if j in T] for j in range(1,5)]
def unpack(x,eps): return x[:15]+1j*np.array(eps)*np.exp(x[15:30])
def eqs(v,K):
    out=[]
    for j in range(4):
        c=np.poly(v[MEM[j]]); sc=np.max(np.abs(v[MEM[j]]))
        out+=[c[k].imag/ (1+np.abs(c[k])) for k in K]
    return np.array(out)
def diag(v):
    allv=np.concatenate([v,np.conj(v)]); iu=np.triu_indices(30,1)
    D=np.abs(allv[:,None]-allv[None,:])[iu]
    return dict(imin=np.min(np.abs(v.imag)),sep=np.min(D),rmax=np.max(np.abs(v)),rmin=np.min(np.abs(v)))
def solve(K,eps,rng,nslice,x0=None,nfev=3000):
    L=rng.normal(size=(nslice,30)); c0=rng.normal(size=nslice)
    def F(x):
        v=unpack(x,eps)
        return np.concatenate([eqs(v,K),[np.log(abs(v[14]))],0.1*(L@x-c0)])
    if x0 is None: x0=np.concatenate([rng.normal(size=15)*rng.choice([0.5,1,2]), rng.normal(size=15)*0.7])
    sol=least_squares(F,x0,method='lm',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=nfev)
    return sol.x, np.linalg.norm(eqs(unpack(sol.x,eps),K))
if __name__=='__main__':
    from signpat import *
