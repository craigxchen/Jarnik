import numpy as np, itertools
from scipy.optimize import least_squares
exec(open('kpte.py').read().split("rng=np.random")[0].replace("k=int(sys.argv[1]); ntry=int(sys.argv[2]); seed=int(sys.argv[3])","k=5"))
rng=np.random.default_rng(5); from collections import Counter; kinds=Counter()
for it in range(60):
    z0=rng.normal(size=nv)*rng.choice([0.3,1,3])
    sol=least_squares(F,z0,jac=J,method='trf',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=2000)
    if np.max(np.abs(sol.fun))>1e-10: continue
    r=unpack(sol.x)
    realT=tuple(Ts[i] for i in range(nT) if abs(r[i].imag)<1e-6)
    c=np.array([np.prod(-r[Mem[j]>0]).imag for j in range(K)])
    # conjugate pairs among roots
    pairs=[(Ts[a],Ts[b]) for a in range(nT) for b in range(a+1,nT) if abs(r[a]-np.conj(r[b]))<1e-6]
    kinds[(len(realT), tuple(np.abs(c)<1e-6), len(pairs))]+=1
    if it<6: print("real roots at",realT,"| c~0:",np.abs(c)<1e-6,"| conj pairs",pairs)
print(kinds)
