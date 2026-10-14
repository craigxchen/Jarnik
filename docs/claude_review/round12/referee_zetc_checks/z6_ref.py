import itertools, numpy as np
from fractions import Fraction
from rlib import Exact, slabbox
from zonocheck import zono
from scipy.optimize import linprog
cols = [[1,2,3,4,5],[0,1,2,3,4,5],[0,3,4],[1,4,5],[1,2,3,4],[0,1,2,3,4,5],[0,2,4,5]]
Ls = [2,2,1,2,2,1,1]
cuts = [(tuple(c), L) for c, L in zip(cols, Ls)]
A = slabbox(6, cuts)
E = Exact(A)
print("|A|", len(A), "H", E.H, "rho", E.rho, "dim V(A,rho) =", len(A) - E.H[-2])
# twin check: sign matrix points (0..5, base=6) x distinct columns
dcols = sorted(set(tuple(c) for c in cols))
rows = [tuple(1 if x in c else 0 for c in dcols) for x in range(6)] + [tuple(0 for c in dcols)]
print("distinct columns", len(dcols), "twin-free:", len(set(rows)) == 7)
z, edges, k0 = zono(6, cuts, A); print("zonotope length (verified)", z, edges)
# LP upper bound on zeta: functionals f in coords + column sums; constraint over primitive v in A-A
Aset = np.array(A)
diffs = set()
for a in A:
    for b in A:
        v = tuple(x - y for x, y in zip(a, b))
        if any(v):
            g = np.gcd.reduce(np.abs(v)); v = tuple(x // g for x in v)
            if next(x for x in v if x) < 0: v = tuple(-x for x in v)
            diffs.add(v)
diffs = sorted(diffs)
funcs = [tuple(1 if y == x else 0 for y in range(6)) for x in range(6)] + [tuple(1 if y in c else 0 for y in range(6)) for c in dcols]
width = [max(sum(f[i]*k[i] for i in range(6)) for k in A) - min(sum(f[i]*k[i] for i in range(6)) for k in A) for f in funcs]
Amat = np.array([[abs(sum(f[i]*v[i] for i in range(6))) for f in funcs] for v in diffs], float)
res = linprog(np.array(width, float), A_ub=-Amat, b_ub=-np.ones(len(diffs)), bounds=[(0, None)]*len(funcs), method="highs")
print("LP zeta bound (float)", res.fun)
# rationalise and verify exactly
y = [Fraction(x).limit_denominator(60) for x in res.x]
ok = all(sum(yy * int(a) for yy, a in zip(y, row)) >= 1 for row in Amat)
print("exact certificate valid:", ok, "value", sum(yy * w for yy, w in zip(y, width)))
