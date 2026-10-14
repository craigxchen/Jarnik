# exact top forms of the 8-point sections of the Paley-12 (L=2) box, and of M8 for comparison; factor over Q
import itertools, sys
import flint, sympy as sp
from tclib import box_points, cube_monomials, mon_eval
def kernel_forms(pts, m, d):
    mons = [a for a in cube_monomials(d) if sum(a) < m]
    M = flint.fmpz_mat(len(mons), len(pts), [mon_eval(k, a) for a in mons for k in pts])
    X, nullity = M.nullspace()
    basis = []
    for j in range(nullity):
        c = [flint.fmpq(X[i, j]) for i in range(len(pts))]
        basis.append(c)
    return basis
def form(pts, c, m, us):
    expr = 0
    for k, ck in zip(pts, c):
        if ck != 0:
            expr += sp.Rational(int(ck.p), int(ck.q)) * (sum(ki*ui for ki, ui in zip(k, us)))**m
    return sp.expand(expr)
q = 11; D = [0,1,3,4,5,9]; L = 2
cols = [sorted((i + x) % q for x in D) for i in range(q)]
us = sp.symbols('u0:7')
reps = [(0,1,2,3,4,5,6), (0,1,2,3,4,5,7), (0,1,2,3,4,5,9)]
for Y in reps:
    idx = {y: i for i, y in enumerate(Y)}
    rc = [[idx[x] for x in c if x in idx] for c in cols]
    rc = [c for c in rc if c]
    pts = box_points(7, [(c, L) for c in rc])
    B = kernel_forms(pts, 8, 7)
    print('Y=%s |A_Y|=%d dimV(8)=%d' % (Y, len(pts), len(B)), flush=True)
    for c in B:
        F = form(pts, c, 8, us)
        fac = sp.factor_list(F)
        print('  factors:', [(sp.Poly(f, *us).total_degree(), e, (str(f) if sp.Poly(f,*us).total_degree() <= 2 else '...')) for f, e in fac[1]], 'const', fac[0], flush=True)
