# Exact: (S_min) fails even when the dependence is a single CIRCUIT (full-support relation, rank deficiency 1).
# Even cycle C_{2k}: vertices v=0..2k-1, g_v=T+v+i; coordinate u_x=g_x*g_{x+1 mod 2k}; prod_{x even} u_x = prod_{x odd} u_x.
import sys; sys.path.insert(0,'.')
from gauss import *
from fractions import Fraction
def rank(rows):
    M=[[Fraction(v) for v in r] for r in rows]; rk=0
    for c in range(len(M[0])):
        piv=next((i for i in range(rk,len(M)) if M[i][c]!=0),None)
        if piv is None: continue
        M[rk],M[piv]=M[piv],M[rk]
        for i in range(len(M)):
            if i!=rk and M[i][c]!=0:
                f=M[i][c]/M[rk][c]; M[i]=[a-f*b for a,b in zip(M[i],M[rk])]
        rk+=1
    return rk
for T in [10**12,10**30]:
  for k in [2,3,4,6,10,15]:
    m=2*k; g=[(T+v,1) for v in range(m)]
    u=[mul(g[x],g[(x+1)%m]) for x in range(m)]
    pe=(1,0); po=(1,0)
    for x in range(m):
        if x%2==0: pe=mul(pe,u[x])
        else: po=mul(po,u[x])
    assert pe==po and nonreal_pairs(u)
    # exponent rank of the ratios u_x/u_0 in the free basis g_v (as formal symbols): edge vectors minus edge 0
    E=[[1 if v in (x,(x+1)%m) else 0 for v in range(m)] for x in range(m)]
    rk=rank([[a-b for a,b in zip(E[x],E[0])] for x in range(1,m)])
    n=m-1; hl,G=h_line(u); s,a=minanchor(u)
    print(f"T=1e{len(str(T))-1} cycle 2k={m}: n={n} formal rank={rk} (n-1={n-1}) h(ell)={hl:.3f} 2logT={2*log(T):.3f} G={G} "
          f"minanchor={s:.3f} ratio/n={s/hl/n:.4f} S-(n/3)h(ell)={s-n/3*hl:.1f}",flush=True)
