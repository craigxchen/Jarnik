"""Ansatz u_i(s) = f_s(theta_i), f_s = A(s) th^3 + B(s) th^2 + C(s) th: triple collision of T <=> [A:B:C] = p_T,
p_T = [1 : -e1(T) : e2(T)].  Balanced curve at k points <=> plane curve of degree d=2^{k-4}... through the
needed triple points.  Exact rank of the evaluation matrix of degree-d plane monomials at the points."""
from fractions import Fraction as Fr
from itertools import combinations
import random, sys

def rank(rows):
    rows = [list(r) for r in rows]; rk = 0; ncol = len(rows[0])
    for c in range(ncol):
        piv = next((i for i in range(rk, len(rows)) if rows[i][c] != 0), None)
        if piv is None: continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        for i in range(len(rows)):
            if i != rk and rows[i][c] != 0:
                f = rows[i][c] / rows[rk][c]
                rows[i] = [a - f * b for a, b in zip(rows[i], rows[rk])]
        rk += 1
    return rk

def monos(d):
    return [(a, b, d - a - b) for a in range(d + 1) for b in range(d + 1 - a)]

random.seed(1)
for k, d in ((6, 3), (7, 6)):
    th = [Fr(random.randint(-40, 40), random.randint(1, 5)) for _ in range(k)]
    pts = []
    for T in combinations(range(k), 3):
        t = [th[i] for i in T]
        e1 = sum(t); e2 = t[0] * t[1] + t[0] * t[2] + t[1] * t[2]
        pts.append((Fr(1), -e1, e2))
    M = [[p[0] ** a * p[1] ** b * p[2] ** c for (a, b, c) in monos(d)] for p in pts]
    print('k=%d: %d triple points, degree-%d forms: %d monomials, evaluation rank %d' % (k, len(pts), d, len(monos(d)), rank(M)))
