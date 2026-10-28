# like staged.py but report which places coincide when a stage collapses
import sys, numpy as np, itertools
from scipy.optimize import least_squares
from kset import *
exec(open('signpat.py').read().split('print(len(pats))')[0])
FULL=frozenset(range(1,6))
PLACES=[T for T in SUBS]+[FULL-T for T in SUBS]
STAGES=[[1,2,3],[1,2,3,8],[1,2,3,7,8],[1,2,3,6,7,8],[1,2,3,5,6,7,8]]
def lab(A): return ''.join(str(a) for a in sorted(A))
N=int(sys.argv[1]); rng=np.random.default_rng(int(sys.argv[2])); out=open(sys.argv[3],'w')
for it in range(N):
    eps=pats[rng.integers(len(pats))]
    try: x,rn=solve(STAGES[0],eps,rng,17)
    except Exception: continue
    if rn>1e-10: continue
    d=diag(unpack(x,eps))
    if not (d['imin']>1e-3 and d['sep']>1e-3): continue
    line=f"{it} {''.join('+' if e>0 else '-' for e in eps)}"
    for st in range(1,5):
        K=STAGES[st]
        def F(y):
            v=unpack(y,eps); return np.concatenate([eqs(v,K),[np.log(abs(v[14]))]])
        sol=least_squares(F,x,method='trf',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=8000)
        v=unpack(sol.x,eps); r=np.linalg.norm(eqs(v,K)); d=diag(v)
        allv=np.concatenate([v,np.conj(v)])
        coinc=[(lab(PLACES[a]),lab(PLACES[b])) for a in range(30) for b in range(a+1,30) if abs(allv[a]-allv[b])<1e-5]
        line+=f" | s{st} res={r:.0e} imin={d['imin']:.0e} coinc={coinc}"
        if r>1e-9 or d['imin']<1e-3 or coinc: break
        x=sol.x
        if st==4: line+=" REACHED44 "+repr(list(v))
    print(line,file=out,flush=True)
