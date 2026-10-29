"""K2: 2-D example S = S' u (S'+(7,3)) of interp Prop 3.3.  Exact over Q.  Own code.
(1) Sidon check; (2) Hilbert function and (SHIFT_e) for all difference shifts and for {+-t};
(3) complete enumeration of rational forms +-Y^a prod_i Q_i(Y^{v_i}) (Q_i products of cyclotomic polynomials)
    supported in S, via centrally symmetric convex lattice polygons with vertices in S (every Newton polygon of such
    a product is a zonogon whose vertices carry nonzero coefficients)."""
import itertools, sys, math
sys.path.insert(0, '.')
from klib import *
import sympy as sp
from math import gcd
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
Sp = [(0,5),(1,3),(1,5),(2,2),(3,4),(4,0),(4,5)]
t = (7,3)
S = sorted(set(Sp + [(a+t[0], b+t[1]) for a, b in Sp]))
pr("|S| =", len(S))
D = [(a[0]-b[0], a[1]-b[1]) for a in Sp for b in Sp if a != b]
pr("S' Sidon (all %d nonzero differences distinct): %s" % (len(D), len(set(D)) == len(D)))
H = hilbert(S, 6, exact=True)
pr("H_S =", H, "sum", sum(H))
Lall = sorted(set((a[0]-b[0], a[1]-b[1]) for a in S for b in S if a != b))
for nm, L in (("all %d differences" % len(Lall), Lall), ("{+-t}", [t, (-t[0], -t[1])])):
    pr("L =", nm, [(e,) + shift_dims_Q(S, L, e) for e in range(5)])
# (3) enumeration
y1, y2 = sp.symbols('y1 y2')
X = sp.symbols('X')
Sset = set(S)
# cyclotomic products of degree l (as polys in X)
ks = [k for k in range(1, 100) if sp.totient(k) <= 11]
def cycprods(l):
    res = []
    def rec(i, deg, cur, m1):
        if deg == l:
            res.append((cur, m1)); return
        for jj in range(i, len(ks)):
            k = ks[jj]; ph = int(sp.totient(k))
            if deg + ph <= l:
                rec(jj, deg + ph, sp.expand(cur * sp.cyclotomic_poly(k, X)), m1 + (k == 1))
    rec(0, 0, sp.Integer(1), 0)
    return res
CP = {l: cycprods(l) for l in range(1, 12)}
def convex_position_sorted(V):
    # V centrally symmetric point list; return vertices in ccw order if in strictly convex position else None
    cx = sum(p[0] for p in V) / len(V); cy = sum(p[1] for p in V) / len(V)
    Vs = sorted(V, key=lambda p: math.atan2(p[1]-cy, p[0]-cx))
    n = len(Vs)
    for i in range(n):
        a, b, c = Vs[i], Vs[(i+1) % n], Vs[(i+2) % n]
        cr = (b[0]-a[0])*(c[1]-b[1]) - (b[1]-a[1])*(c[0]-b[0])
        if cr <= 0:
            return None
    return Vs
best = 0; bestex = []; npoly = 0; nforms = 0
sums = {}
for p in S:
    for q in S:
        sums.setdefault((p[0]+q[0], p[1]+q[1]), set()).add(p)
seen = set()
for c2, pts in sums.items():
    pts = sorted(pts)
    pairs = sorted(set(tuple(sorted([p, (c2[0]-p[0], c2[1]-p[1])])) for p in pts if p != (c2[0]-p[0], c2[1]-p[1])))
    for r in range(1, len(pairs) + 1):
        for sub in itertools.combinations(pairs, r):
            V = [p for pr_ in sub for p in pr_]
            if r == 1:
                Vs = list(sub[0]); gens = [(Vs[1][0]-Vs[0][0], Vs[1][1]-Vs[0][1])]
                a0 = Vs[0]
            else:
                Vs = convex_position_sorted(V)
                if Vs is None:
                    continue
                n = len(Vs); g = n // 2
                gens = [(Vs[(i+1) % n][0]-Vs[i][0], Vs[(i+1) % n][1]-Vs[i][1]) for i in range(g)]
                a0 = None
            key = tuple(sorted(V))
            if key in seen:
                continue
            seen.add(key); npoly += 1
            # generators: (length l, primitive v)
            gl = []
            for (gx, gy) in gens:
                l = gcd(abs(gx), abs(gy)); gl.append((l, (gx // l, gy // l)))
            if any(l > 11 for l, _ in gl):
                continue
            for choice in itertools.product(*[CP[l] for l, _ in gl]):
                Ppoly = sp.Integer(1); order = 0
                for (Qx, m1), (l, v) in zip(choice, gl):
                    # Q(Y^v) with Y^v = y1^vx y2^vy (Laurent)
                    Ppoly = sp.expand(Ppoly * Qx.subs(X, y1**v[0] * y2**v[1]))
                    order += m1
                nforms += 1
                # support as exponent set (Laurent): translate so that it fits in S
                Pp = sp.Poly(sp.expand(Ppoly * y1**30 * y2**30), y1, y2)
                supp = [(e[0]-30, e[1]-30) for e in Pp.monoms()]
                # all translates a with supp + a subset S: anchor min point
                ok = False
                for s0 in S:
                    a = (s0[0]-supp[0][0], s0[1]-supp[0][1])
                    if all((u[0]+a[0], u[1]+a[1]) in Sset for u in supp):
                        ok = True; break
                if ok:
                    if order > best:
                        best = order; bestex = [(sp.factor(Ppoly), a)]
                    elif order == best and len(bestex) < 4:
                        bestex.append((sp.factor(Ppoly), a))
pr("centrally symmetric vertex sets examined: %d; candidate forms: %d" % (npoly, nforms))
pr("max order at e of a rational cyclotomic-binomial product supported in S:", best)
for ex in bestex:
    pr("   example:", ex)
open("k2_twodim_out.txt", "w").write("\n".join(out) + "\n")
