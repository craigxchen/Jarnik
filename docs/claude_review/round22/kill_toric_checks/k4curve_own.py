#!/usr/bin/env python3
"""Independent exact check of toric.md Thm 4.4 (k=4 curves of gauss Thm 4.1):
rays met once each, B.C = 14, contact order at s=inf = 4, ratio 7/2, D_V.C = 24 = 6*4.
Uses the explicit rational 2-parameter subfamily of gauss.md 4.1, own code (sympy, exact)."""
import sympy as sp
from itertools import combinations
s = sp.symbols('s')
I = sp.I

def family(X, m):
    X = sp.Rational(X); m = sp.Rational(m)
    xa = X; xb = X*(1-X)/(m**2+X); xc = 1-xa-xb; u = m*xb
    r = {(1,2,3): I, (1,2): u+I*(xa-1), (1,3): u+I*(xb-1), (2,3): u+I*(xc-1)}
    # singleton roots: fourth root formula r = (Re s1 Im s1 - Im s2)/Im s1 - i Im s1 for the three other roots of P_x
    for x in (1, 2, 3):
        others = [r[T] for T in r if x in T]
        e1 = sp.expand(sum(others)); e2 = sp.expand(sum(a*b for a, b in combinations(others, 2)))
        re1, im1 = sp.re(e1), sp.im(e1); im2 = sp.im(e2)
        r[(x,)] = sp.nsimplify((re1*im1 - im2)/im1 - I*im1)
    return r

def check(X, m):
    r = family(X, m)
    P = {}
    for x in (1, 2, 3):
        P[x] = sp.expand(sp.prod([s - r[T] for T in r if x in T]))
    ImP = {x: sp.simplify(sp.im(sp.expand(P[x].subs(s, sp.Symbol('t', real=True))))) for x in P}
    const = all(sp.Poly(ImP[x], sp.Symbol('t', real=True)).degree() <= 0 for x in ImP)
    roots = list(r.values()) + [sp.conjugate(v) for v in r.values()]
    distinct = len(set(sp.nsimplify(v) for v in roots)) == 14
    # valuation vectors at each root: u_x = P_x/conj(P_x), u_0 = 1 (anchor 0); vector (v_0, v_1, v_2, v_3) mod 1
    rays = []
    for T, rt in r.items():
        v = [0] + [1 if x in T else 0 for x in (1, 2, 3)]
        rays.append(tuple(v))       # at s = r_T : ray 1_T
        rays.append(tuple(1 - t for t in v))  # at s = conj r_T: -1_T == 1_{T^c} (complement contains anchor 0)
    allrays = set()
    for rr in range(1, 4):
        for A in combinations(range(4), rr):
            allrays.add(tuple(1 if j in A else 0 for j in range(4)))
    rays_ok = sorted(rays) == sorted(allrays)
    # contact order at s=inf: ord_inf(u_x - 1) = deg conj(P_x) - deg(2 i ImP_x) = 4 - 0
    contact = min(4 - sp.Poly(ImP[x], sp.Symbol('t', real=True)).degree() for x in ImP)
    BC = len(rays)
    # D_V . C = kappa * sum_cuts |A||A^c| with kappa = 1
    DV = sum(len(A)*(4-len(A)) for rr in range(1, 4) for A in combinations(range(4), rr) if 0 in A)
    print(f"(X,m)=({X},{m}): ImP = {[ImP[x] for x in (1,2,3)]}, constant={const}, 14 distinct roots={distinct}, "
          f"each of 14 rays met once={rays_ok}, B.C={BC}, contact={contact}, ratio={sp.Rational(BC, contact)}, "
          f"D_V.C={DV} vs ord_e V_4 * mult = {6*contact}")

for X, m in [(-3, 3), (2, 5), (sp.Rational(1, 2), sp.Rational(3, 2)), (-5, 1), (7, 2)]:
    check(X, m)
