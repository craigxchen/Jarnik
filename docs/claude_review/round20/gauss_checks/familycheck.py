"""(a) Are the (L_5)/(L_6) refuting families circle (Gaussian) data, read with the isotropic points at +-i?
Circle data require, for every visible block T and every row j in T:  n_T | N(P_j) = X_j^2 + Y_j^2  (isotropy:
X_j = rho_T Y_j mod n_T with rho_T^2 = -1).  We report, per member, the fraction
  phi = log gcd(N(P_j), prod_{T∋j} n_T) / log prod_{T∋j} n_T   (min over j);   circle data need phi = 1.
(1) ptolemy.md C1 (k=5): exact resultants Res_s(n_T(s), X_i(s)^2+1) (nonzero => n_T(s) | X_i^2+1 only for |n_T(s)| <= |Res|).
(2) threshold/L8.md E_t (k=6) at t0=(2,3,-5,1), points N=1..14 of item 109 / Thm B (scripts imported read-only)."""
import sys, os, math
from fractions import Fraction as Fr
from math import gcd, log
from itertools import combinations
import sympy as sp
s = sp.symbols('s')
a = [1,2,3,4]; tau = [4,6,7,0,2]
X = {i+1: sp.expand(a[i]*(s-tau[0])*sp.prod([s-tau[l] for l in range(1,5) if l != i+1])) for i in range(4)}
blocks = {}
blocks[(1,2,3,4)] = s - tau[0]
for l in range(1,5): blocks[tuple(j for j in range(1,5) if j != l)] = s - tau[l]
for i, j in combinations(range(1,5), 2):
    ai, aj = a[i-1], a[j-1]; blocks[(i,j)] = sp.expand((ai-aj)*s - (ai*tau[j]-aj*tau[i]))
print("(1) C1 family: visible blocks", {k: str(v) for k,v in blocks.items() if 2 <= len(k) <= 3})
mx = 0
for T, nT in blocks.items():
    if not (2 <= len(T) <= 3): continue
    for i in T:
        r = sp.resultant(nT, X[i]**2 + 1, s); assert r != 0; mx = max(mx, abs(int(r)))
print("    all resultants Res(n_T, X_i^2+1) nonzero; max |Res| =", mx,
      "=> for every s with some |n_T(s)| > %d the member is not circle data with rho = i" % mx)
# also: Y_i = 1 forces N(P_i) = X_i^2+1 ; check phi on the explicit class s = 52346140003 + 85900454640 q
for q in (0, 1, 7):
    sv = 52346140003 + 85900454640*q
    Xv = {i: int(X[i].subs(s, sv)) for i in X}; nv = {T: abs(int(v.subs(s, sv))) for T, v in blocks.items() if 2 <= len(T) <= 3}
    ph = min(log(gcd(Xv[j]**2+1, math.prod(nv[T] for T in nv if j in T))) / log(math.prod(nv[T] for T in nv if j in T)) for j in X)
    print("    C1 member s=%d: min_j phi = %.4f (circle data need 1)" % (sv, ph))
# (2) E_t
sys.dont_write_bytecode = True
sys.path.insert(0, '/home/user/Jarnik/docs/claude_review/threshold/L8/k6fam')
from fam import family
from ell import quotient_curve
from lift import lift_T
from profile import profile, sides, mobius_to_std, row, det
def anchored_P(xs, pt):
    d, b, Dm = profile(xs, pt); x1, x2, x3 = xs[:3]; M = mobius_to_std(x1, x2, x3)
    pts = [None, Fr(0), Fr(1)] + [M(z) for z in pt]; V = [row(p) for p in pts]
    Dl = 1
    for j in range(1, 6):
        y = abs(V[j][1]); Dl = Dl * y // gcd(Dl, y)
    P = [(1, 0)]
    for j in range(1, 6):
        Xj, Yj = V[j][0] * Dl, V[j][1]; g = gcd(Xj, Yj); P.append((Xj // g, Yj // g))
    n = {}
    for C in sides():
        T = C if 0 not in C else frozenset(range(6)) - C
        n[T] = n.get(T, 1) * d[C]
    return P, n
xs = [Fr(2), Fr(3), Fr(-5), Fr(1)]
F, C, g, pts = quotient_curve(xs); k = g[3]
Q1, Q2 = pts[0], pts[1]; Dl = C.add(Q2, (Q1[0], -Q1[1])); R = Q1
print("(2) E_t, t0 = (2,3,-5,1):")
for N in range(1, 15):
    R = C.add(R, Dl); L = lift_T(F, R[0] / k)
    try: P, n = anchored_P(xs, L[0])
    except ZeroDivisionError: continue
    vis = {T: v for T, v in n.items() if 2 <= len(T) <= 4}
    phis = []
    for j in range(1, 6):
        pr = math.prod(v for T, v in vis.items() if j in T)
        if pr > 1: phis.append(log(gcd(P[j][0]**2 + P[j][1]**2, pr)) / log(pr))
    w = sum(log(v) for v in vis.values()) / len(vis)
    print("    N=%2d w=%.1f  min_j phi = %.4f  max_j phi = %.4f (circle data need 1)" % (N, w, min(phis), max(phis)))
