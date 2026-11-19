import random, sys
sys.path.insert(0,'.')
from gf2 import *
from sympy import Poly, symbols, GF
from sympy import factor_list
t=symbols('t')
rng=random.Random(5)
for trial in range(60):
    f = rng.getrandbits(rng.randint(2,90)) | (1<<rng.randint(1,90))
    # random products incl. powers
    g = mul(f, mul(f, rng.getrandbits(30)|1)) if trial%3==0 else f
    if deg(g)<1: continue
    fac = factor(g, seed=trial)
    # compare with sympy
    coeffs=[ (g>>i)&1 for i in range(deg(g),-1,-1)]
    P=Poly(coeffs,t,modulus=2)
    _,fl=P.factor_list()
    sym = sorted((sum(int(c)%2<<i for i,c in enumerate(reversed(q.all_coeffs()))), e) for q,e in fl)
    mine = sorted(fac.items())
    assert sym==mine, (g, sym, mine)
print("factor agrees with sympy on 60 random polys")
# deriv check: d(fg)=f'g+fg'
for _ in range(100):
    a=rng.getrandbits(40); b=rng.getrandbits(40)
    assert deriv(mul(a,b)) == mul(deriv(a),b)^mul(a,deriv(b))
print("Leibniz ok")
