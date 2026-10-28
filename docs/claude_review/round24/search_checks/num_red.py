import sys, numpy as np
from scipy.optimize import least_squares
from red import *
seed=int(sys.argv[1]); N=int(sys.argv[2]); out=open(sys.argv[3],'w')
rng=np.random.default_rng(seed); nconv=0; nnd=0
for it in range(N):
    x0 = rng.normal(size=21)*rng.choice([0.5,1,2])
    try:
        sol = least_squares(residual, x0, method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=600)
    except Exception as e:
        continue
    x=sol.x; rn=np.linalg.norm(residual(x))
    if rn>1e-10: continue
    nconv+=1
    d=diagnostics(x)
    nd = d['imin']>1e-3 and d['sep']>1e-3 and d['csep']>1e-3 and d['du']>1e-5 and d['rmax']<1e4
    nnd+=nd
    print(it, f"res={rn:.1e}", ' '.join(f"{k}={v:.3e}" for k,v in d.items()), 'NONDEG' if nd else '', list(x) if nd else '', file=out, flush=True)
print('converged',nconv,'nondeg',nnd,file=out)
