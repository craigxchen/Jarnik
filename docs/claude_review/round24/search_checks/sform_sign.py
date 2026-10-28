# s-form with fixed half-plane sign pattern: rho_T = x_T + i eps_T exp(eta_T)
import itertools, numpy as np, sys
from scipy.optimize import least_squares
SUBS=[frozenset(c) for k in range(1,5) for c in itertools.combinations(range(1,5),k)]
IDX={T:i for i,T in enumerate(SUBS)}; MEM=[[IDX[T] for T in SUBS if j in T] for j in range(1,5)]
exec(open('signpat.py').read().split('print(len(pats))')[0])
def reps():
    def act(p,perm,fl):
        d={frozenset(perm[x-1] for x in T):b*fl for b,T in zip(p,SUBS)}
        return tuple(d[T] for T in SUBS)
    seen=set(); out=[]
    for p in pats:
        if p in seen: continue
        seen|={act(p,perm,fl) for perm in itertools.permutations(range(1,5)) for fl in (1,-1)}; out.append(p)
    return out
REPS=reps()
def unpack(x,eps): return x[:15]+1j*np.array(eps)*np.exp(x[15:30])
def eqs(v):
    out=[]
    for j in range(4):
        w=v[MEM[j]]
        out+= [np.sum(w**k).imag for k in (1,2,3)]+[np.sum(w**(-k)).imag for k in (1,2,3)]+[np.sin(np.sum(np.angle(w)))]
    return np.array(out)
def diag(v):
    allv=np.concatenate([v,np.conj(v)]); iu=np.triu_indices(30,1)
    D=np.abs(allv[:,None]-allv[None,:])[iu]
    kap=[np.poly(v[MEM[j]])[4].imag for j in range(4)]
    return dict(imin=np.min(np.abs(v.imag)),sep=np.min(D),kap=np.min(np.abs(kap)),rmax=np.max(np.abs(v)),rmin=np.min(np.abs(v)))
if __name__=='__main__':
    pid=int(sys.argv[1]); N=int(sys.argv[2]); seed=int(sys.argv[3]); out=open(sys.argv[4],'w')
    eps=REPS[pid]; rng=np.random.default_rng(seed); nc=nd=0
    for it in range(N):
        L=rng.normal(size=30); c0=rng.normal()
        def F(x):
            v=unpack(x,eps)
            return np.concatenate([eqs(v),[np.log(abs(v[14])), 0.1*(L@x-c0)]])
        x0=np.concatenate([rng.normal(size=15)*rng.choice([0.5,1,2]), rng.normal(size=15)*0.7])
        try: sol=least_squares(F,x0,method='lm',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=3000)
        except Exception: continue
        x=sol.x; v=unpack(x,eps); rn=np.linalg.norm(eqs(v))
        if not np.isfinite(rn) or rn>1e-10: continue
        nc+=1; d=diag(v)
        ok=d['imin']>1e-3 and d['sep']>1e-3 and d['kap']>1e-6 and d['rmax']<1e4 and d['rmin']>1e-4
        nd+=ok
        print(it,f"res={rn:.1e}",' '.join(f"{a}={b:.3e}" for a,b in d.items()),'NONDEG' if ok else '',list(v) if ok else '',file=out,flush=True)
    print('pattern',pid,''.join('+' if b>0 else '-' for b in eps),'converged',nc,'nondeg',nd,file=out)
