"""Exact contact profile of actual rational six-point configurations on E_t, by gcd/Moebius inversion
(no factoring).  Anchors x1,x2,x3 -> infinity,0,1; rows primitive; sides C with <=1 anchor."""
from fractions import Fraction as Fr
from math import gcd, log
from itertools import combinations, product
from lift import lift_T
from ell import quotient_curve

LAB = ['A1', 'A2', 'A3', 'a', 'b', 'c']

def mobius_to_std(x1, x2, x3):
    # M(z) = ((z - x2)(x3 - x1)) / ((z - x1)(x3 - x2)) : x1->inf, x2->0, x3->1
    def M(z):
        if z == x1:
            return None
        return ((z - x2) * (x3 - x1)) / ((z - x1) * (x3 - x2))
    return M

def row(q):
    if q is None:
        return (1, 0)
    q = Fr(q)
    return (q.numerator, q.denominator)  # (X, Y) with point X/Y; primitive

def det(u, v):
    return u[0] * v[1] - u[1] * v[0]

def sides():
    mov = [3, 4, 5]
    S = []
    for j in range(3):
        for r in range(1, 4):
            for s in combinations(mov, r):
                S.append(frozenset((j,) + s))
    for r in (2, 3):
        for s in combinations(mov, r):
            S.append(frozenset(s))
    return S

def gcd_list(xs):
    g = 0
    for x in xs:
        g = gcd(g, x)
    return g

def profile(xs, pt):
    x1, x2, x3 = xs[:3]
    M = mobius_to_std(x1, x2, x3)
    pts = [None, Fr(0), Fr(1)] + [M(z) for z in pt]
    V = [row(p) for p in pts]
    D = {frozenset((i, j)): abs(det(V[i], V[j])) for i, j in combinations(range(6), 2)}
    SC = sides()
    d = {}
    for C in sorted(SC, key=len, reverse=True):
        G = gcd_list([D[frozenset(p)] for p in combinations(sorted(C), 2)])
        sup = 1
        for C2 in SC:
            if C2 > C:
                sup *= d[C2]
        assert G % sup == 0, (C, G, sup)
        d[C] = G // sup
    b = {}
    for p, val in D.items():
        prodC = 1
        for C in SC:
            if p <= C:
                prodC *= d[C]
        assert val % prodC == 0
        b[p] = val // prodC
    return d, b, D

if __name__ == '__main__':
    import itertools
    xs = [Fr(2), Fr(3), Fr(-5), Fr(1)]
    F, C, g, pts = quotient_curve(xs)
    k = g[3]
    def comb(v):
        R = None
        for n, P in zip(v, pts):
            if n < 0:
                P = (P[0], -P[1]); n = -n
            for _ in range(n):
                R = C.add(R, P)
        return R
    rows = []
    for v in [(-2, -1, 0, 0), (-2, 0, 0, -1), (-3, 0, 0, 0), (-4, 0, -1, 0), (-4, -1, -1, -1), (-4, -2, -1, -1), (-4, -3, -2, -1)]:
        R = comb(v)
        if R is None: continue
        L = lift_T(F, R[0] / k)
        if not L: continue
        pt = L[0]
        d, b, D = profile(xs, pt)
        logs = sorted(log(x) for x in d.values())
        maxb = max(b.values())
        pg = max(gcd(d[A], d[B]) for A, B in combinations(d, 2))
        rows.append((v, logs[0], logs[-1], maxb, pg))
        print(v, 'log d_C range [%.2f, %.2f]' % (logs[0], logs[-1]), ' mean %.2f' % (sum(logs) / len(logs)),
              ' max b_ij =', maxb, ' max pairwise gcd(d_C,d_D) =', pg)
