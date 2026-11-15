# 5-subsets of the cube: gcd degree, mu0, alpha, rank, cones (exact) and EXACT Jacobian rank of sweep map.
import itertools, sympy as sp
from collections import Counter
t,y0,y1,y2=sp.symbols('t y0 y1 y2')
I=sp.I
ys=[y0,y1,y2,-(y0+y1+y2)]
eps=(1,-1,-1,1)
roots=[]
for e,v in zip(eps,ys): r=e/v+I*v; roots+=[r,-r]
def conjr(x): return sp.conjugate(x).subs({sp.conjugate(y0):y0,sp.conjugate(y1):y1,sp.conjugate(y2):y2})
MASKS=[0b00000000,0b00001111,0b00110011,0b01010101,0b10101010,0b11001100,0b11110000,0b11111111]
def z(m): return sp.prod([t-(conjr(roots[k]) if (m>>k)&1 else roots[k]) for k in range(8)])
Z={m:z(m) for m in MASKS}
pt={y0:sp.Rational(-9,4),y1:sp.Rational(-2,3),y2:sp.Rational(1,4)}
ZJ={m:sp.Poly(sp.expand(Z[m].subs(pt)),t) for m in MASKS}
bits=lambda m:[(m>>k)&1 for k in range(8)]
F=y0*y1*y2*(-(y0+y1+y2))-1
grad=[sp.diff(F,v).subs(pt) for v in (y0,y1,y2)]
tang=sp.Matrix([grad]).nullspace()   # 2 tangent vectors
stats=Counter(); cones=0; jr=Counter()
t0=sp.Rational(1,3)
for S in itertools.combinations(MASKS,5):
    g=sum(1 for k in range(8) if len(set(bits(m)[k] for m in S))==1)
    mu0=8-max((ZJ[a]-ZJ[b]).degree() for a,b in itertools.combinations(S,2))
    rk=sp.Matrix([[bits(m)[k]-bits(S[0])[k] for k in range(8)] for m in S[1:]]).rank()
    stats[(g,mu0,sp.Rational(2*(8-g),mu0),rk)]+=1
    for a in range(5):
        o=[j for j in range(5) if j!=a]
        for (b,c) in [(o[0],o[1]),(o[0],o[2]),(o[0],o[3])]:
            d,e=[j for j in o if j not in (b,c)]
            A=lambda j:ZJ[S[j]]-ZJ[S[a]]
            if (A(b)*A(c)-A(d)*A(e)).is_zero: cones+=1
    # Jacobian of (z_j/z_S0) j=1..4 w.r.t. (tangent1,tangent2,t) at (yJ,t0), exact
    rows=[]
    for m in S[1:]:
        u=Z[m]/Z[S[0]]
        gy=[sp.diff(u,v).subs(pt).subs(t,t0) for v in (y0,y1,y2)]
        row=[sp.nsimplify(sp.expand(sum(gy[i]*tv[i] for i in range(3)))) for tv in tang]+[sp.diff(u,t).subs(pt).subs(t,t0)]
        rows.append([sp.simplify(x) for x in row])
    jr[(rk,sp.Matrix(rows).rank())]+=1
print("(g,mu0,alpha,rank):",dict(stats)); print("cone identities:",cones,"of",56*15)
print("(code rank, exact Jacobian rank):",dict(jr))
