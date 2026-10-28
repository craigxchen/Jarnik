# staged continuation: K = {1,2,3} -> +8 -> +7 -> +6 -> +5 (contact (4,0) -> (4,4)), nondegeneracy tracked
import sys, numpy as np
from scipy.optimize import least_squares
from kset import *
exec(open('signpat.py').read().split('print(len(pats))')[0])
STAGES=[[1,2,3],[1,2,3,8],[1,2,3,7,8],[1,2,3,6,7,8],[1,2,3,5,6,7,8]]
N=int(sys.argv[1]); rng=np.random.default_rng(int(sys.argv[2])); out=open(sys.argv[3],'w')
reach=[0]*5
for it in range(N):
    eps=pats[rng.integers(len(pats))]
    try: x,rn=solve(STAGES[0],eps,rng,17)
    except Exception: continue
    if rn>1e-10: continue
    d=diag(unpack(x,eps))
    if not (d['imin']>1e-3 and d['sep']>1e-3): continue
    reach[0]+=1; line=f"{it} {''.join('+' if e>0 else '-' for e in eps)} s0 ok"
    for st in range(1,5):
        K=STAGES[st]
        def F(y):
            v=unpack(y,eps); return np.concatenate([eqs(v,K),[np.log(abs(v[14]))]])
        try: sol=least_squares(F,x,method='trf',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=8000)
        except Exception: break
        v=unpack(sol.x,eps); r=np.linalg.norm(eqs(v,K)); d=diag(v)
        good = r<1e-10 and d['imin']>1e-3 and d['sep']>1e-3 and d['rmax']<1e4 and d['rmin']>1e-4
        line+=f" | s{st} res={r:.1e} imin={d['imin']:.1e} sep={d['sep']:.1e} rmax={d['rmax']:.1e}"
        if not good: break
        reach[st]+=1; x=sol.x
        if st==4: line+=" REACHED44 "+repr(list(v))
    print(line,file=out,flush=True)
print('reach',reach,file=out)
