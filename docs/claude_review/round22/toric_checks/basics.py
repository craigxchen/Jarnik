# Exact toric bookkeeping for the Losev-Manin (permutohedral) variety LM_k, k=3,4,5.
# N = Z^k/Z1, rays 1_A; P_{-K} = {sum x = 0, <x,1_A> <= 1 for all proper A}.
import itertools
from fractions import Fraction as F
from math import comb, factorial
import sympy as sp

def subsets(k):
    return [A for r in range(1,k) for A in itertools.combinations(range(k),r)]

def check(k):
    print(f"=== k={k}")
    rays = subsets(k)
    # P_{-K} = {sum x = 0, sum_{A} x <= 1}; claim = {sum x=0, sum x_+ <= 1} = {||x||_1 <= 2}
    # vertices of max cones (perm sigma): <m,1_{A_t}> = -1 for A_t = {s1..st}
    verts=set()
    for s in itertools.permutations(range(k)):
        m=[0]*k
        # partial sums -1 for t=1..k-1, total 0
        prev=0
        for t in range(k):
            val = -1 if t<k-1 else 0
            m[s[t]] = val-prev; prev=val
        verts.add(tuple(m))
    roots = {tuple((1 if j==a else -1 if j==b else 0) for j in range(k)) for a in range(k) for b in range(k) if a!=b}
    print("vertices of max cones = roots e_a-e_b:", verts==roots, "#distinct", len(verts), "#maxcones", factorial(k))
    # lattice points of mR, Ehrhart
    def count(m):
        c=0
        for x in itertools.product(range(-m,m+1),repeat=k-1):
            last=-sum(x); y=list(x)+[last]
            if sum(v for v in y if v>0)<=m: c+=1
        return c
    ms=list(range(0,k+1)); vals=[count(m) for m in ms]
    t=sp.symbols('t'); poly=sp.interpolate(list(zip(ms,vals)),t)
    print("Ehrhart L(mR):", sp.expand(poly), "check m=k+1:", poly.subs(t,k+1)==count(k+1))
    lead=sp.Poly(poly,t).LC(); vol_norm=lead*factorial(k-1)
    print("(-K)^%d = normalized volume ="%(k-1), vol_norm, " binom(2k-2,k-1) =", comb(2*k-2,k-1))
    # widths of R in direction v = 2(max v - min v): check on random v
    import random
    for _ in range(5):
        v=[random.randint(-3,3) for _ in range(k)]
        w=max(sum(a*b for a,b in zip(r,v)) for r in roots)-min(sum(a*b for a,b in zip(r,v)) for r in roots)
        assert w==2*(max(v)-min(v))
    print("width_v(R)=2(max v-min v): ok on samples")
    # mean-width bound Omega for R: 2^{1-k} sum_{P subset [k]} width(1_P) = 2^{1-k}*(2^k-2)*2
    Om=F(2*(2**k-2),2**(k-1))
    vand_m = max(sum(sorted(range(k))[-r:]) - r*F(k-1,2) for r in range(1,k))  # Vandermonde support in direction 1_A
    print("Omega(R) =",Om, " Vandermonde ord",comb(k,2)," needs m =",vand_m," ratio",F(comb(k,2))/vand_m)
    # balanced class gamma: E_A.gamma=1/2; -K.gamma = (2^k-2)/2
    print("-K.gamma =",F(2**k-2,2), " Liouville-critical mu (per gamma) =", F(2**k-2,2)/Om*1, "(= ell/ord)")
for k in (3,4,5): check(k)
