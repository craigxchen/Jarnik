# general k: coordinates 1..k-1 (+ k normalised), blocks T subset [k-1], branch pair +-1, contact mu = 2^(k-3) each
import itertools, numpy as np
from scipy.optimize import least_squares
class Sys:
    def __init__(s, k):
        s.k=k; n=k-1; s.n=n
        s.BIG=[frozenset(c) for kk in range(2,n+1) for c in itertools.combinations(range(1,n+1),kk)]
        s.nb=len(s.BIG); s.BIDX={T:i for i,T in enumerate(s.BIG)}
        s.SMEM=[[s.BIDX[T] for T in s.BIG if j in T] for j in range(1,n+1)]
        s.deg=2**(n-1); s.mu=s.deg//2
        s.Wc=(np.poly1d([1,0,-1])**s.mu).coeffs.astype(float)
        s.nx=2*s.nb-1
    def unpack(s,x):
        r=np.zeros(s.nb,dtype=complex); r[:-1]=x[:s.nb-1]+1j*x[s.nb-1:2*s.nb-2]; r[-1]=1j*x[-1]; return r
    def solve_j(s,r,j):
        P=np.array([1.0+0j])
        for kk in s.SMEM[j]: P=np.convolve(P,[1.0,-r[kk]])
        P9=np.concatenate([[0],P]); tP=np.concatenate([P,[0]])
        V=np.stack([P9.real,P9.imag,tP.real,tP.imag],axis=1)
        coef,*_=np.linalg.lstsq(V,s.Wc,rcond=None)
        mu=coef[1]+1j*coef[0]; lam=coef[3]+1j*coef[2]
        return s.Wc-V@coef, lam, mu, lam*tP+mu*P9
    def residual(s,x):
        r=s.unpack(x); return np.concatenate([s.solve_j(r,j)[0] for j in range(s.n)])
    def roots(s,x):
        r=s.unpack(x); out=list(r); gs=[]
        for j in range(s.n):
            _,lam,mu,g=s.solve_j(r,j); out.append(-mu/lam); gs.append(g)
        return np.array(out), gs
    def diag(s,x):
        v,gs=s.roots(x)
        iu=np.triu_indices(len(v),1)
        sep=np.min(np.abs(v[:,None]-v[None,:])[iu]); csep=np.min(np.abs(v[:,None]-np.conj(v)[None,:]))
        tt=0.3
        us=[np.polyval(g,tt)/np.conj(np.polyval(g,tt)) for g in gs]
        du=min([abs(a-b) for a,b in itertools.combinations(us,2)]+[abs(u-1) for u in us])
        return dict(imin=np.min(np.abs(v.imag)),sep=sep,csep=csep,du=du,rmax=np.max(np.abs(v)))
def run(k,seed,N,out,eps=0.0,nfev=800):
    S=Sys(k); rng=np.random.default_rng(seed); nc=nd=0
    iu=None
    def F(x):
        res=S.residual(x)
        if eps>0:
            v,_=S.roots(x); m=len(v); iu=np.triu_indices(m,1)
            D=np.abs(v[:,None]-v[None,:])[iu]; C=np.abs(v[:,None]-np.conj(v)[None,:])
            res=res*(1+eps*(np.sum(1/D**2)+np.sum(1/C**2)/2))
        return res
    f=open(out,'w')
    for it in range(N):
        x0=rng.normal(size=S.nx)*rng.choice([0.5,1,2])
        try: sol=least_squares(F,x0,method='lm',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=nfev)
        except Exception: continue
        x=sol.x; rn=np.linalg.norm(S.residual(x))
        if not np.isfinite(rn) or rn>1e-10: continue
        nc+=1; d=S.diag(x)
        ok=d['imin']>1e-3 and d['sep']>1e-3 and d['csep']>1e-3 and d['du']>1e-5 and d['rmax']<1e4
        nd+=ok
        print(it,f"res={rn:.1e}",' '.join(f"{a}={b:.3e}" for a,b in d.items()),'NONDEG' if ok else '',list(x) if ok else '',file=f,flush=True)
    print('converged',nc,'nondeg',nd,file=f); f.close()
if __name__=='__main__':
    import sys
    run(int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]),sys.argv[4],float(sys.argv[5]),int(sys.argv[6]) if len(sys.argv)>6 else 800)
