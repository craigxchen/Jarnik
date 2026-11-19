import sympy as sp
from sympy import I, Rational as Q, Poly, expand, Matrix, gcd
from fam import *
seedsets = [(Q(1,2),Q(-2),Q(7,3),Q(4)), (Q(-1,3),Q(5,2),Q(6),Q(-7)), (Q(0),Q(1),Q(3,5),Q(-11,4))]
for n in (1,2):
  for sd in seedsets:
    G,Gb,zs = marks(n, sd)
    D = 8*n+1
    degs = [Poly(z,t).degree() for z in zs]
    norms = [expand(z*expand(sp.conjugate(z))) for z in zs]
    eqn = all(expand(norms[0]-m)==0 for m in norms)
    pd = sorted(set(Poly(expand(zs[a]-zs[b]),t).degree() for a in range(5) for b in range(a+1,5)))
    g = zs[0]
    for z in zs[1:]:
        g = sp.gcd(Poly(g,t,domain='QQ_I'), Poly(z,t,domain='QQ_I'))
    lc = [Poly(z,t).LC() for z in zs]
    print(n, sd, 'degs',set(degs),'D',D,'equal norms',eqn,'pair-diff degs',pd,'gcd deg',Poly(g,t).degree(),'LC',set(lc))
