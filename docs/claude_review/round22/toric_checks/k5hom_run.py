# Sample Newton-homotopy endpoints; classify stalled/ended endpoints by nondegeneracy measures.
import numpy as np, sys, warnings
warnings.filterwarnings('ignore')
from k5complex import F,J,unpack,cvals
from k5homotopy import track
n=int(sys.argv[1]); rng=np.random.default_rng(int(sys.argv[2])); out=[]
for k in range(n):
    z0=(rng.normal(size=28)+1j*rng.normal(size=28))*rng.choice([0.5,1,2])
    z,t,st=track(z0)
    r,rp=unpack(z); c=cvals(r,rp)
    allr=np.concatenate([r,rp]); dmin=min(abs(allr[a]-allr[b]) for a in range(30) for b in range(a+1,30))
    sv=np.linalg.svd(J(z),compute_uv=False) if np.all(np.isfinite(z)) else [0,1]
    out.append((st,t,float(np.min(np.abs(c))),float(min(abs(c[a]-c[b]) for a in range(4) for b in range(a+1,4))),float(dmin),float(sv[-1]/sv[0]),float(np.linalg.norm(z))))
    print(*out[-1],flush=True)
