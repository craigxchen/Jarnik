# Calibration: k=4 system cut to expected dimension 0 by 3 random affine constraints on the unknowns
# (fix z[0], z[1], z[nT]) -> how often does LM from random starts reach a NONDEGENERATE isolated solution?
import numpy as np, sys, itertools
from scipy.optimize import least_squares
src=open('kpte.py').read().split("rng=np.random")[0].replace("k=int(sys.argv[1]); ntry=int(sys.argv[2]); seed=int(sys.argv[3])","k=4")
exec(src)
rng=np.random.default_rng(7); stats=[]
for inst in range(5):
    vals=rng.normal(size=3); fixidx=[0,1,nT]
    free=[i for i in range(nv) if i not in fixidx]
    def full(x):
        z=np.empty(nv); z[fixidx]=vals; z[free]=x; return z
    Fr=lambda x: F(full(x)); Jr=lambda x: J(full(x))[:,free]
    good=0; conv=0
    for it in range(100):
        x0=rng.normal(size=len(free))*rng.choice([0.3,1,3])
        sol=least_squares(Fr,x0,jac=Jr,method='lm',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=3000)
        if np.max(np.abs(sol.fun))<1e-10:
            conv+=1; r=unpack(full(sol.x)); c=np.array([np.prod(-r[Mem[j]>0]).imag for j in range(K)])
            nd=min(np.min(np.abs(r.imag)),np.min(np.abs(c)),min(abs(c[a]-c[b]) for a in range(K) for b in range(a+1,K)))
            if nd>1e-5: good+=1
    print(f"instance {inst}: converged {conv}/100, nondegenerate {good}/100",flush=True)
