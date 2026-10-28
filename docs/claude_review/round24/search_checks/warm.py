import sys, numpy as np
from scipy.optimize import least_squares
from kset import *
exec(open('signpat.py').read().split('print(len(pats))')[0])
K0=[int(c) for c in sys.argv[1].split(',')]; K1=[int(c) for c in sys.argv[2].split(',')]
N=int(sys.argv[3]); rng=np.random.default_rng(int(sys.argv[4]))
for it in range(N):
    eps=pats[rng.integers(len(pats))]
    x,rn=solve(K0,eps,rng,max(1,29-4*len(K0)))
    if rn>1e-10: continue
    d0=diag(unpack(x,eps))
    if not (d0['imin']>1e-3 and d0['sep']>1e-3): continue
    def F(y):
        v=unpack(y,eps); return np.concatenate([eqs(v,K1),[np.log(abs(v[14]))]])
    sol=least_squares(F,x,method='trf',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=5000)
    v=unpack(sol.x,eps); r1=np.linalg.norm(eqs(v,K1)); d=diag(v)
    extra=[k for k in K1 if k not in K0]
    resid_extra=[ [np.poly(v[MEM[j]])[k].imag for j in range(4)] for k in extra]
    print(it,f"start-res(K0)={rn:.1e} end-res(K1)={r1:.2e}",' '.join(f"{a}={b:.2e}" for a,b in d.items()), 'extra',np.round(np.array(resid_extra),4).tolist(),flush=True)
