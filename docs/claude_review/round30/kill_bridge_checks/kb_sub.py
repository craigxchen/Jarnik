# Sub-cluster constants C_inf^(k) for all subsets of size 5,6,7 of the cube cluster (independent code).
import itertools, math, sympy as sp
from fractions import Fraction as Fr
from kb_gauss import *
from kb_local import badprimes
from kb_local2 import C_inf
y=[Fr(-9,4),Fr(-2,3),Fr(1,4),Fr(8,3)]
rts=roots_cube(y)
P=badprimes(rts)
# gammas from exact sympy polys
t=sp.symbols('t')
def poly(m):
    rr=[]
    for k,(a,b) in enumerate(rts):
        bb=-b if (m>>k)&1 else b
        rr.append(sp.Rational(a.numerator,a.denominator)+sp.I*sp.Rational(bb.numerator,bb.denominator))
    return sp.Poly(sp.expand(sp.prod([t-x for x in rr])),t,domain='QQ_I')
polys={m:poly(m) for m in MASKS}
def gam(a,b):
    d=(polys[a]-polys[b]); co=d.all_coeffs()[::-1]
    return abs(complex(sp.N(co[4]))) if len(co)>4 else 0.0
best={}
for k in (5,6,7):
    for S in itertools.combinations(MASKS,k):
        # common columns -> polynomial gcd degree
        common=[j for j in range(8) if len(set((m>>j)&1 for m in S))==1]
        assert not common
        Gam=max(gam(a,b) for a,b in itertools.combinations(S,2))
        c,_=C_inf(rts,list(S),Gam,P=P)
        if k not in best or c<best[k][0]: best[k]=(c,S)
        print(k,[format(m,'08b')[::-1] for m in S],round(Gam,4),c,flush=True)
print("BEST",{k:(v[0],[format(m,'08b')[::-1] for m in v[1]]) for k,v in best.items()})
