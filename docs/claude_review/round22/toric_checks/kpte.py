# General k: balanced block-degree-1 extremal curves <=> Im p_n({r_T: T∋j}) = 0, n=1..2^{k-2}-1, j=1..k-1,
# T ranging over nonempty subsets of [k-1]. Unknowns r_T (2^{k-1}-1 complex). Gauge: sum Re r = 0, Im r_top = 1.
# Nondegeneracy checked a posteriori (min |Im r_T|, |c_j| distinct nonzero). Calibration at k=4 (3-dim family exists).
import numpy as np, itertools, sys
from scipy.optimize import least_squares
k=int(sys.argv[1]); ntry=int(sys.argv[2]); seed=int(sys.argv[3])
K=k-1; Ts=[T for r in range(1,K+1) for T in itertools.combinations(range(K),r)]; nT=len(Ts)
deg=2**(K-1); nmax=deg-1
iTop=Ts.index(tuple(range(K))); idx=[i for i in range(nT) if i!=iTop]
Mem=np.array([[1.0 if j in T else 0.0 for T in Ts] for j in range(K)])
nv=2*nT-2
Da=np.zeros((nT,nv)); Da[:nT-1,:nT-1]=np.eye(nT-1); Da[nT-1,:nT-1]=-1
Db=np.zeros((nT,nv))
for q,i in enumerate(idx): Db[i,nT-1+q]=1
def unpack(z):
    a=np.empty(nT); b=np.empty(nT); a[:nT-1]=z[:nT-1]; a[nT-1]=-z[:nT-1].sum(); b[idx]=z[nT-1:]; b[iTop]=1.0
    return a+1j*b
def F(z):
    r=unpack(z); pw=np.array([r**n for n in range(1,nmax+1)])
    return (Mem@pw.imag.T).T.ravel()
def J(z):
    r=unpack(z); dr=np.array([n*r**(n-1) for n in range(1,nmax+1)])
    rows=[]
    for n in range(nmax):
        for j in range(K):
            rows.append((Mem[j]*dr[n].imag)@Da+(Mem[j]*dr[n].real)@Db)
    return np.array(rows)
rng=np.random.default_rng(seed); good=0; best=[]
for it in range(ntry):
    sc=rng.choice([0.3,1,3]); z0=rng.normal(size=nv)*sc
    sol=least_squares(F,z0,jac=J,method='trf',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=2000)
    m=np.max(np.abs(sol.fun)); r=unpack(sol.x)
    c=np.array([np.prod(-r[Mem[j]>0]).imag for j in range(K)])
    nd=min(np.min(np.abs(r.imag)), np.min(np.abs(c)), min([abs(c[a]-c[b]) for a in range(K) for b in range(a+1,K)]+[9]))
    best.append((m,nd))
    if m<1e-10 and nd>1e-4:
        good+=1
        if good<=5: print(f"try {it}: res {m:.1e} nondeg {nd:.3e} r={np.round(r,4)}")
print(f"k={k}: {ntry} tries, nondegenerate solutions {good}; converged (res<1e-10) {sum(1 for m,_ in best if m<1e-10)}")
ms=sorted(best); print("smallest residuals among nondegenerate-ish (nd>1e-3):", sorted([m for m,nd in best if nd>1e-3])[:8])
