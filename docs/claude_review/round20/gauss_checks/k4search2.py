"""Filtered search: nondegenerate rational k=4 degree-1 families (rho=i). Requires Im sigma_1 != 0 for each of the
three triples, singleton roots non-real, all 7 roots distinct and pairwise non-conjugate."""
from fractions import Fraction as F
import itertools, sys
from k4search import phi, neighbours
def sing(t):
    us=[x[0] for x in t]; vs=[x[1] for x in t]
    # complex arithmetic with Fractions
    def mul(p,q): return (p[0]*q[0]-p[1]*q[1], p[0]*q[1]+p[1]*q[0])
    s1=(sum(us),sum(vs)); s2=(F(0),F(0))
    for x,y in itertools.combinations(t,2):
        m=mul(x,y); s2=(s2[0]+m[0],s2[1]+m[1])
    if s1[1]==0: return None
    x=(s1[0]*s1[1]-s2[1])/s1[1]; return (x, -s1[1])
def good(a,b,c):
    rho=(F(0),F(1)); roots=[rho,a,b,c]
    for t in [(rho,a,b),(rho,a,c),(rho,b,c)]:
        r=sing(t)
        if r is None or r[1]==0: return None
        roots.append(r)
    if any(r[1]==0 for r in roots): return None
    S=set(roots)
    if len(S)<7: return None
    if any((r[0],-r[1]) in S for r in roots): return None
    return roots
H = int(sys.argv[1]) if len(sys.argv)>1 else 3
vals = sorted({F(p,q) for q in range(1,H+1) for p in range(-3*H,3*H+1)} - {F(0), F(-1)})
ugrid = sorted({F(p,q) for q in range(1,H+1) for p in range(-3*H,3*H+1)})
out=[]
for a in itertools.product(ugrid, vals):
    N=[b for b in neighbours(a, vals) if b!=a]
    for b,c in itertools.combinations(N,2):
        if b!=c and phi(b,c)==0:
            g=good(a,b,c)
            if g: out.append(g)
print("nondegenerate solutions:", len(out))
for g in out[:12]: print([(str(x[0]),str(x[1])) for x in g])
