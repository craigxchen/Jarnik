"""Exact check of the building blocks of Prop 3.1: V (in H^0(3B)), V_abc (H^0(B)), D_ab (H^0(B)).
For each: support in S_m, iota-parity (f(x^{-1}) = +-f(x)), exact order at e (all monomial moments
sum_s c_s s^alpha vanish for |alpha| < t, some moment of degree t nonzero), and sigma_X-behaviour
(Q-rational iff c_{-s} = conj(c_s)). Laurent polynomials as dicts exponent-tuple -> integer coefficient."""
import itertools
from collections import defaultdict
from sympy import symbols, expand, Poly
def mul(f,g):
    h=defaultdict(int)
    for a,ca in f.items():
        for b,cb in g.items():
            h[tuple(x+y for x,y in zip(a,b))]+=ca*cb
    return {k:v for k,v in h.items() if v}
def e(i,c=1):
    t=[0]*5; t[i]=c; return tuple(t)
def mono(v): return {tuple(v):1}
def diff(a,b): # x_a - x_b
    return {e(a):1, e(b):-1}
def V_abc(a,b,c):
    f=mul(mul(diff(a,b),diff(b,c)),diff(a,c))
    sh=[0]*5; sh[a]=sh[b]=sh[c]=-1
    return mul(f,mono(sh))
def V():
    f={tuple([0]*5):1}
    for a,b in itertools.combinations(range(5),2): f=mul(f,diff(a,b))
    return mul(f,mono([-2]*5))
def D(a,b):
    v=[0]*5; v[a]=1; v[b]=-1; w=[-x for x in v]
    return {tuple(v):1, tuple(w):1, tuple([0]*5):-2}
def inS(f,m): return all(sum(s)==0 and sum(abs(x) for x in s)<=2*m for s in f)
def parity(f):
    g={tuple(-x for x in s):c for s,c in f.items()}
    if g==f: return +1
    if g=={s:-c for s,c in f.items()}: return -1
    return 0
y=symbols('y0:4')
def order(f,maxd=20):
    for j in range(maxd+1):
        P=expand(sum(c*(sum(s[i]*y[i] for i in range(4)))**j for s,c in f.items()))
        if P!=0: return j
    return None
for name,f,m in [("V",V(),3),("V_123",V_abc(0,1,2),1),("D_12",D(0,1),1),
                 ("V_123*V_145",mul(V_abc(0,1,2),V_abc(0,3,4)),2),
                 ("V_123*V_145*V_245",mul(mul(V_abc(0,1,2),V_abc(0,3,4)),V_abc(1,3,4)),3),
                 ("V_123*D_45",mul(V_abc(0,1,2),D(3,4)),2)]:
    par=parity(f)
    # Q-rationality for real integer coefficients: c_{-s} = conj(c_s) = c_s  <=> even
    print(f"{name}: in S_{m}: {inS(f,m)}, iota-parity {par:+d}, ord_e = {order(f)}, "
          f"Q-rational as is: {par==1}, i*f Q-rational: {par==-1}, #terms {len(f)}")
