from pell4 import *
from itertools import combinations
def F(k):
    a,b=0,1
    for _ in range(k): a,b=b,a+b
    return a
def primary(z): return z[0]%2==1 and z[1]%2==0 and (z[0]+z[1])%4==1
beta=(-1,-2)
for d in [(0,1,2),(0,1,3),(0,2,5),(1,2,4),(0,3,4)]:
    print("offsets",d)
    for n in (1,2,4,8,12):
        bl=[(F(6*(n+dj)+1),F(6*(n+dj))) for dj in d]
        assert all(primary(b) and primary(cj(b)) for b in bl) and primary(beta) and primary(cj(beta))
        pts=words(bl,pre0=beta,pre=cj(beta))
        assert all(primary(p) for p in pts)
        g,Hp,e,K,q=analyse(pts)
        gp = g if primary(g) else None
        # normalise gcd to primary (odd) and check quotients primary
        units=[(1,0),(0,1),(-1,0),(0,-1)]
        gpr=[mul(u,g) for u in units if primary(mul(u,g))]
        qp=[gdiv(p,gpr[0]) for p in pts] if gpr else None
        allprim = qp is not None and all(primary(z) for z in qp)
        print(f"  n={n} N(gcd)={N(g)} bits={Hp.bit_length()} logD/logH={float(e):.5f} arc/R'^(1/3)={float(K[0]):.4f} quotients primary={allprim}")
