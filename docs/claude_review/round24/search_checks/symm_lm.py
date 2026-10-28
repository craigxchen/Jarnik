import numpy as np, sys, itertools
from scipy.optimize import least_squares
from symm import *
def peqs(rho):
    out=[]
    for j in range(4):
        w=rho[MEM[j]]
        out+=[np.sum(w**k).imag for k in (1,2,3)]+[np.sum(w**(-k)).imag for k in (1,2,3)]+[np.sin(np.sum(np.angle(w)))]
    return np.array(out)
def ndiag(rho):
    v=rho; iu=np.triu_indices(30,1); D=np.abs(v[:,None]-v[None,:])[iu]
    kap=[np.poly(rho[MEM[j]])[4].imag for j in range(4)]
    return dict(imin=np.min(np.abs(v.imag)),sep=np.min(D),kap=np.min(np.abs(kap)),rmax=np.max(np.abs(v)),rmin=np.min(np.abs(v)))
reps={'id':(1,2,3,4,5),'(12)':(2,1,3,4,5),'(15)':(5,2,3,4,1),'(12)(34)':(2,1,4,3,5),'(12)(35)':(2,1,5,4,3)}
if __name__=='__main__':
    N=int(sys.argv[1]); out=open(sys.argv[2],'w'); rng=np.random.default_rng(int(sys.argv[3]))
    cases=[(g,pn,e) for g in ['-s','1/s','-1/s'] for pn in reps for e in (1,-1)]
    half=int(sys.argv[4]); cases=cases[half::2]
    for g,pn,e in cases:
        r=analyse([sym_elt(g,reps[pn],e),KAPPA],rng)
        if not isinstance(r,dict) or r['npar']==0: continue
        build=r['build']; npar=r['npar']; nc=nd=0; best=None
        for it in range(N):
            L=rng.normal(size=npar); c0=rng.normal()
            def F(p):
                rho=build(p); res=peqs(rho)
                extra=[0.1*(L@p-c0)] if r['gauge'] else []
                return np.concatenate([res,extra])
            p0=rng.normal(size=npar)*rng.choice([0.5,1,2])
            try: sol=least_squares(F,p0,method='lm',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=1000)
            except Exception: continue
            rho=build(sol.x); rn=np.linalg.norm(peqs(rho))
            if not np.isfinite(rn) or rn>1e-10: continue
            nc+=1; d=ndiag(rho)
            ok=d['imin']>1e-3 and d['sep']>1e-3 and d['kap']>1e-6 and d['rmax']<1e4 and d['rmin']>1e-4
            if ok:
                nd+=1; print('  NONDEG',g,pn,e,' '.join(f"{a}={b:.3e}" for a,b in d.items()),list(sol.x),file=out,flush=True)
        print(g,pn,e,'npar',npar,'gauge',r['gauge'],'converged',nc,'nondeg',nd,file=out,flush=True)
