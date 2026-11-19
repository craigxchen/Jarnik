import itertools, math
from f2lib import *
from irrlib import IRR, vP, iota
def irr_upto(e): return [P for P in IRR if deg(P)<=e]
for r in range(2,6):
    m=2**r; d=r+1; Q=2**d; L=(d-2)*2**(d-1); D=(m-1)*L
    masks=[a for a in range(Q) if a&1]
    H={a:[v for v in range(Q) if bin(a&v).count('1')%2==0] for a in masks}
    tot=0
    for a,b in itertools.combinations(masks,2):
        J=set(H[a])&set(H[b]); PJ1=1
        for v in J: PJ1=mul(PJ1,1^v)
        tot+=iota(PJ1)
    s=m; closed=3/8*s*(s-1)*(int(math.log2(s))//2)
    print(f"m={m} D={D}  (4/s)*sum iota={4*tot/s:.1f}  closed-form bound={closed:.1f}  ratio D/(m^2 log2 m)={D/(m*m*math.log2(m)):.3f}")
