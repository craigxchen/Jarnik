import sys, numpy as np, itertools
from kset import *
exec(open('signpat.py').read().split('print(len(pats))')[0])
Kname=sys.argv[1]; K={'33':[1,2,6,7,8],'43':[1,2,3,6,7,8],'44':[1,2,3,5,6,7,8],'40':[1,2,3],'30':[1,2],'50':[1,2,3,4],'41':[1,2,3,8],'42':[1,2,3,7,8]}[Kname]
N=int(sys.argv[2]); rng=np.random.default_rng(int(sys.argv[3])); out=open(sys.argv[4],'w')
nslice=30-1-4*len(K)  # make square-ish: unknowns 30, eqs 4|K| + 1 gauge + slices
nslice=max(nslice,0)+0
stats={}
for it in range(N):
    eps=pats[rng.integers(len(pats))]
    try: x,rn=solve(K,eps,rng,max(1,nslice))
    except Exception as e: continue
    if not np.isfinite(rn) or rn>1e-10: continue
    d=diag(unpack(x,eps))
    ok=d['imin']>1e-3 and d['sep']>1e-3 and d['rmax']<1e4 and d['rmin']>1e-4
    stats[ok]=stats.get(ok,0)+1
    print(it,f"res={rn:.1e}",' '.join(f"{a}={b:.3e}" for a,b in d.items()),'NONDEG' if ok else '', list(x) if ok else '', ''.join('+' if e>0 else '-' for e in eps),file=out,flush=True)
print('stats',stats,file=out)
