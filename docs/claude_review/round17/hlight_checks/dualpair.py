# sigma/r of 'dual pairs' f=(g_b1 +- g_b2)/2 + corrections at b1,b2; compared with best 4-term (sparse) chars.
import numpy as np, itertools
from hlib import *
def primes3(lo,hi):
    out=[]
    for p in range(lo,hi):
        if p%4==3 and all(p%d for d in range(2,int(p**.5)+1)): out.append(p)
    return out
for q in [7,11,19,23,31,43,47,59,67,71,79,83,103,199,499,1019,2003,3251]:
    ch=chi_vec(q); best=None
    for sgn in [1,-1]:
      for b2 in [1,2,3]:
        g1=np.roll(ch,0); g2=np.array([ch[(x-b2)%q] for x in range(q)])
        g1=np.array([ch[x%q] for x in range(q)])
        base=(g1+sgn*g2)
        for e1,e2 in itertools.product([-1,0,1],repeat=2):
            f=base.copy()
            f[0]=0; f[b2]=0
            if np.any(f%2): continue
            f//=2; f[0]=e1; f[b2]=e2
            if np.abs(f).sum()<=2: continue
            s=sigma_formula(f,q,ch)
            if best is None or s<best[0]: best=(s,sgn,b2,e1,e2,int(np.abs(f).sum()),int(f.sum()))
    # 4-term: min over e_0+e_x-e_y-e_z (by AGL symmetry fix 0, scan small set)
    print(f"q={q} dualpair min sigma={best[0]} sigma/r={best[0]/(2*q):.4f} (sgn,b2,e1,e2,||f||1,F)={best[1:]}",flush=True)
