# Exact check of (II.4): J(u) = prod_{i<j} t_ij / H^3 for u_i = P_i/conj(P_i),
# and of the beta-exponent formula in the uniqueness calculation, on random Gaussian blocks.
import random
from fractions import Fraction as F
from itertools import combinations
random.seed(31)
class G:  # Gaussian rational
    __slots__=("a","b")
    def __init__(s,a,b=0): s.a=F(a); s.b=F(b)
    def __add__(s,o): o=o if isinstance(o,G) else G(o); return G(s.a+o.a,s.b+o.b)
    def __sub__(s,o): o=o if isinstance(o,G) else G(o); return G(s.a-o.a,s.b-o.b)
    def __mul__(s,o): o=o if isinstance(o,G) else G(o); return G(s.a*o.a-s.b*o.b,s.a*o.b+s.b*o.a)
    def conj(s): return G(s.a,-s.b)
    def nrm(s): return s.a*s.a+s.b*s.b
    def inv(s): n=s.nrm(); return G(s.a/n,-s.b/n)
    def __truediv__(s,o): o=o if isinstance(o,G) else G(o); return s*o.inv()
    def __pow__(s,k):
        r=G(1)
        base=s if k>=0 else s.inv()
        for _ in range(abs(k)): r=r*base
        return r
    def __eq__(s,o): o=o if isinstance(o,G) else G(o); return s.a==o.a and s.b==o.b
def prod(xs):
    r=G(1)
    for x in xs: r=r*x
    return r
marks=range(5); pairs=list(combinations(marks,2))
ntrial=0; nunif=0
for trial in range(300):
    while True:
        beta={e:G(random.randint(-60,60),random.randint(-60,60)) for e in pairs}
        if all(b.nrm()!=0 for b in beta.values()): break
    n={e:beta[e].nrm() for e in pairs}
    P=[prod(beta[e] for e in pairs if i in e) for i in marks]
    u=[p/p.conj() for p in P]
    t={}
    ok=True
    for (i,j) in pairs:
        det=(P[i].conj()*P[j]).b    # det(P_i,P_j)=Im(conj(P_i)P_j)
        if det==0: ok=False;break
        assert det % n[(i,j)] ==0 if det.denominator==1 else False
        t[(i,j)]=det/n[(i,j)]
    if not ok: continue
    H=1
    for e in pairs: H*=n[e]
    J=prod(u[i]-u[j] for i,j in pairs)/(G(0,2)**10)/(prod(u)**2)
    assert J==G(prod(G(t[e]) for e in pairs).a/F(H)**3), trial
    ntrial+=1
    # uniqueness exponent formula with random exponents
    a=[random.randint(-3,3) for _ in marks]
    v={e:random.randint(-3,3) for e in pairs}
    a[4]-= sum(a)+sum(v.values())   # homogeneity
    s=[sum(v[e] for e in pairs if i in e) for i in marks]
    lhs=prod(u[i]**a[i] for i in marks)*prod((u[i]-u[j])**v[(i,j)] for i,j in pairs)
    rhs=G(0,-2)**sum(v.values())
    for (i,j) in pairs:
        rhs=rhs*G(t[(i,j)])**v[(i,j)]*beta[(i,j)]**(v[(i,j)]+a[i]+a[j])*beta[(i,j)].conj()**(v[(i,j)]-s[i]-s[j]-a[i]-a[j])
    assert lhs==rhs, trial
    nunif+=1
print("II.4 identity verified exactly on",ntrial,"random block sets; exponent formula verified on",nunif)
