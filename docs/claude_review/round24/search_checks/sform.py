# s-form: rho_T (15 complex), e_k(R_j) real for k in {1,2,3,5,6,7,8}; gauge |rho_1234| fixed by rho_1234 = i*exp(x)?? we fix rho_1234 = cos(phi)+i sin(phi) (unit circle)
import itertools, numpy as np, sys
from scipy.optimize import least_squares
SUBS=[frozenset(c) for k in range(1,5) for c in itertools.combinations(range(1,5),k)]
IDX={T:i for i,T in enumerate(SUBS)}; MEM=[[IDX[T] for T in SUBS if j in T] for j in range(1,5)]
KS=[1,2,3,5,6,7,8]
def unpack(x):
    r=x[:14]+1j*x[14:28]; return np.concatenate([r,[np.exp(1j*x[28])]])
def residual(x):
    r=unpack(x); out=[]
    for j in range(4):
        v=r[MEM[j]]
        out.append(np.array([np.sum(v**k).imag for k in (1,2,3)]+[np.sum(v**(-k)).imag for k in (1,2,3)]+[np.sin(np.sum(np.angle(v)))]))
    return np.concatenate(out)
def diag(x):
    v=unpack(x); allv=np.concatenate([v,np.conj(v)])
    iu=np.triu_indices(30,1); D=np.abs(allv[:,None]-allv[None,:])[iu]
    kap=[np.poly(v[MEM[j]])[4].imag for j in range(4)]
    return dict(imin=np.min(np.abs(v.imag)), sep=np.min(D), kap=np.min(np.abs(kap)), rmax=np.max(np.abs(v)), rmin=np.min(np.abs(v)))
def run(seed,N,out,eps,nfev=1500):
    rng=np.random.default_rng(seed); nc=nd=0; f=open(out,'w')
    def F(x):
        res=np.concatenate([residual(x),[0.1*(L@x-c0)]])
        if eps>0:
            v=unpack(x); allv=np.concatenate([v,np.conj(v)]); iu=np.triu_indices(30,1)
            D=np.abs(allv[:,None]-allv[None,:])[iu]
            kap=np.array([np.poly(v[MEM[j]])[4].imag for j in range(4)])
            res=res*(1+eps*(np.sum(1/D**2)+np.sum(1/kap**2)+np.sum(1/np.abs(v)**2)))
        return res
    for it in range(N):
        L=rng.normal(size=29); c0=rng.normal()
        x0=np.concatenate([rng.normal(size=28)*rng.choice([0.5,1,2]),[rng.uniform(0,2*np.pi)]])
        try: sol=least_squares(F,x0,method='lm',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=nfev)
        except Exception: continue
        x=sol.x; rn=np.linalg.norm(residual(x))
        if not np.isfinite(rn) or rn>1e-10: continue
        nc+=1; d=diag(x)
        ok=d['imin']>1e-3 and d['sep']>1e-3 and d['kap']>1e-6 and d['rmax']<1e4 and d['rmin']>1e-4
        nd+=ok
        print(it,f"res={rn:.1e}",' '.join(f"{a}={b:.3e}" for a,b in d.items()),'NONDEG' if ok else '',list(x) if ok else '',file=f,flush=True)
    print('converged',nc,'nondeg',nd,file=f); f.close()
if __name__=='__main__':
    run(int(sys.argv[1]),int(sys.argv[2]),sys.argv[3],float(sys.argv[4]))
