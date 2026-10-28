import sys, numpy as np, itertools
from scipy.optimize import least_squares
from kset import *
exec(open('signpat.py').read().split('print(len(pats))')[0])
FULL=frozenset(range(1,6)); PLACES=[T for T in SUBS]+[FULL-T for T in SUBS]
def lab(A): return ''.join(str(a) for a in sorted(A))
K0=[1,2,3]; K1=[int(c) for c in sys.argv[1].split(',')]
rng=np.random.default_rng(int(sys.argv[2]))
for it in range(int(sys.argv[3])):
    eps=pats[rng.integers(len(pats))]
    x,rn=solve(K0,eps,rng,17)
    v=unpack(x,eps); d=diag(v)
    if rn>1e-10 or d['imin']<1e-3 or d['sep']<1e-3: continue
    def F(y):
        v=unpack(y,eps); return np.concatenate([eqs(v,K1),[np.log(abs(v[14]))]])
    sol=least_squares(F,x,method='trf',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=8000)
    v=unpack(sol.x,eps); r=np.linalg.norm(eqs(v,K1)); d=diag(v)
    allv=np.concatenate([v,np.conj(v)])
    coinc=[(lab(PLACES[a]),lab(PLACES[b])) for a in range(30) for b in range(a+1,30) if abs(allv[a]-allv[b])<1e-4]
    nearreal=[lab(PLACES[a]) for a in range(15) if abs(v[a].imag)<1e-3]
    print(it,f"res={r:.1e} imin={d['imin']:.1e} sep={d['sep']:.1e} rmax={d['rmax']:.1e} rmin={d['rmin']:.1e}",'coinc',coinc[:12],'real',nearreal,flush=True)
