# find rational y=(a,b,c,d): a+b+c+d=0, abcd=1, no pair-sum zero; build the 8-root cluster over Q(i) and verify EXACTLY.
from fractions import Fraction as Fr
import itertools, sympy as sp
from math import isqrt
def issq(fr):
    if fr<0: return None
    n,d=fr.numerator,fr.denominator; a,b=isqrt(n),isqrt(d)
    return Fr(a,b) if a*a==n and b*b==d else None
sols=[]
vals=sorted({Fr(p,q) for p in range(-12,13) for q in range(1,7) if p!=0})
for a,b in itertools.combinations(vals,2):
    ab=a*b; disc=ab*ab*(a+b)**2-4*ab
    s=issq(disc)
    if s is None: continue
    for sg in (1,-1):
        c=(-ab*(a+b)+sg*s)/(2*ab)
        if c==0: continue
        d=-(a+b+c)
        y=[a,b,c,d]
        if any(y[i]+y[j]==0 for i in range(4) for j in range(i+1,4)): continue
        if len(set(y))<4: continue
        sols.append(tuple(sorted(y)))
sols=sorted(set(sols),key=lambda y:max(max(abs(v.numerator),v.denominator) for v in y))
print(len(sols)); print(sols[:5])
t=sp.symbols('t')
I=sp.I
def cj(p):
    P=sp.Poly(sp.expand(p),t)
    return sp.expand(sum(sp.conjugate(c)*t**m[0] for m,c in P.terms()))
def build(y):
    eps=[1,-1,-1,1]
    r=[sp.Rational(e)/sp.Rational(v)+I*sp.Rational(v) for e,v in zip(eps,y)]  # x=eps/y
    roots={'000':r[0],'010':-r[0],'001':r[1],'011':-r[1],'100':r[2],'110':-r[2],'101':r[3],'111':-r[3]}
    return roots
def verify(y):
    roots=build(y)
    Q=sp.expand(sp.prod([t-v for v in roots.values()]))
    pts=[]
    for F in itertools.product([0,1],repeat=8):
        pass
    keys=list(roots)
    cl=[]
    for mask in range(256):
        F=[keys[k] for k in range(8) if (mask>>k)&1]
        z=sp.expand(sp.prod([t-sp.conjugate(roots[k]) if k in F else t-roots[k] for k in keys]))
        dz=sp.Poly(sp.expand(z-Q),t)
        if dz.is_zero or dz.degree()<=4: cl.append((mask,z))
    N=sp.expand(Q*cj(Q))
    ok=all(sp.expand(z*cj(z)-N)==0 for m,z in cl)
    pair=max(sp.Poly(sp.expand(z1-z2),t).degree() for (m1,z1),(m2,z2) in itertools.combinations(cl,2))
    return len(cl),ok,pair,Q,N
for y in sols[:3]:
    n,ok,pair,Q,N=verify(y)
    print("y=",y," cluster size",n," norms equal",ok," max deg diff",pair, " deg N",sp.Poly(N,t).degree())
    print("  Q=",Q)
    print("  N squarefree:",sp.Poly(N,t).sqf_part().degree()==16)
