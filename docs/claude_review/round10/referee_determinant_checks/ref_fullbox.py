# Referee Proposition R1 check (full box, no slab cut): base locus {delta: f_delta|_B in Poly_{<m}|_B}
# equals the sparse strata { #{x<n : delta_x != delta_n} * 2T < m }.  Exact rational test of the
# Newton-basis coefficient criterion on random rational zeta-1 values is replaced here by a direct
# numerical projection (evidence) plus the exact sup-norm of the Newton dual vectors (<= 2^{|alpha|}).
import numpy as np, itertools, math, sys
rng=np.random.default_rng(3)
out=[]
def pr(s): print(s,flush=True); out.append(s)
for (n,T,m) in [(4,1,3),(4,1,5),(5,1,3),(5,1,5),(4,2,5),(4,2,7)]:
    d=n-1
    B=np.array(list(itertools.product(range(-T,T+1),repeat=d)),dtype=float)
    mons=[a for t in range(m) for a in itertools.product(range(t+1),repeat=d) if sum(a)==t]
    P=np.stack([np.prod(B**np.array(a),axis=1) for a in mons],axis=1)
    Uu,sv,_=np.linalg.svd(P,full_matrices=False); rk=int((sv>sv[0]*1e-10).sum()); Q=Uu[:,:rk]
    res={}
    for s in range(0,d+1):
        vals=[]
        for _ in range(20):
            dl=np.zeros(d); S=rng.choice(d,size=s,replace=False); dl[S]=rng.uniform(0.3,2.5,size=s)
            f=np.exp(1j*(B@dl)); r=f-Q@(Q.T@f); vals.append(np.linalg.norm(r)/np.linalg.norm(f))
        res[s]=max(vals) if s*2*T<m else min(vals)
    pred={s:(s*2*T<m) for s in res}
    pr(f"n={n} T={T} m={m}: support size s -> (max residual if predicted zero / min residual if predicted nonzero): "
       + ", ".join(f"s={s}:{'Z' if pred[s] else 'NZ'} {res[s]:.1e}" for s in res))
open(sys.argv[1],"w").write("\n".join(out)+"\n")
