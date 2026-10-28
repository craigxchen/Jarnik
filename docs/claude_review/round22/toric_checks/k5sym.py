# k=5, d=1 extremal-curve system under the symmetry ansatz r_{gT} = -conj(r_T), g an involution of [4]
# (the reflection s -> -s is a real affine reparametrisation, so the ansatz is consistent; c_{gj} = +-c_j).
# Full equation set (redundant) solved by least squares; nondegeneracy checked a posteriori.
import numpy as np, itertools, sys
from scipy.optimize import least_squares
Ts=[T for r in range(1,5) for T in itertools.combinations(range(4),r)]
Mem=np.array([[1.0 if j in T else 0.0 for T in Ts] for j in range(4)])
def make(g):
    gT=lambda T: tuple(sorted(g[x] for x in T))
    fixed=[T for T in Ts if gT(T)==T]; reps=[]
    seen=set()
    for T in Ts:
        if T in seen or T in fixed: continue
        reps.append(T); seen.add(T); seen.add(gT(T))
    return gT,fixed,reps
def build(g):
    gT,fixed,reps=make(g)
    nf,nr=len(fixed),len(reps)
    # params: imag parts of fixed (first fixed has Im=1 for scale) ; complex reps
    def unpack(z):
        r={}
        r[fixed[0]]=1j
        for q,T in enumerate(fixed[1:]): r[T]=1j*z[q]
        off=nf-1
        for q,T in enumerate(reps):
            v=z[off+2*q]+1j*z[off+2*q+1]; r[T]=v; r[gT(T)]=-np.conj(v)
        return np.array([r[T] for T in Ts])
    nv=nf-1+2*nr
    def F(z):
        r=unpack(z); pw=np.array([r**n for n in range(1,8)])
        return (Mem@pw.imag.T).T.ravel()
    return unpack,nv,F,fixed,reps
g=[int(c) for c in sys.argv[1]]  # e.g. 1032 for (12)(34) in 0-based
ntry=int(sys.argv[2]); seed=int(sys.argv[3])
unpack,nv,F,fixed,reps=build(g)
print("g=",g,"fixed",fixed,"pairs",len(reps),"unknowns",nv)
rng=np.random.default_rng(seed); good=[]; conv=0
for it in range(ntry):
    z0=rng.normal(size=nv)*rng.choice([0.3,1,3])
    sol=least_squares(F,z0,method='lm',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=5000)
    m=np.max(np.abs(sol.fun)); r=unpack(sol.x)
    if m<1e-10:
        conv+=1
        c=np.array([np.prod(-r[Mem[j]>0]).imag for j in range(4)])
        nd=min(np.min(np.abs(r.imag)), np.min(np.abs(c)), min(abs(c[a]-c[b]) for a in range(4) for b in range(a+1,4)))
        if nd>1e-5:
            good.append(sol.x); print(f"try {it} res {m:.1e} nd {nd:.2e} c={np.round(c,5)}\n r={np.round(r,6)}",flush=True)
print("converged",conv,"nondegenerate",len(good))
np.save(f"k5sym_{sys.argv[1]}_{seed}.npy",np.array(good) if good else np.zeros((0,nv)))
