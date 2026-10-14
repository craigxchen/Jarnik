"""Exact test of ZE-TC^rat on rational planes.

For a model A and a rational plane M = span(v1, v2) (integer vectors), N = nu_A(M) (generic value; computed
exactly as nu at v1 + T v2 for a large integer T, then confirmed by the gcd being a nonzero polynomial).
The points u = v1 + t v2 (t real) with nu_A(u) < N are the common real roots of
    g_c(t) = Phi^{(N)}_c(v1 + t v2) = sum_k c_k (<k,v1> + t <k,v2>)^N,   c in V(A, N).
If the gcd of the g_c has an irreducible factor of degree >= 2 with a real root, then that u has rational hull M
but nu(u) < nu(M): ZE-TC^rat fails at order nu(u) + 1 (exactly).  Rational roots are harmless (u on a rational line).
All arithmetic exact (fmpz nullspace, fmpq_poly gcd/factor).
"""
import sys, random, itertools
import flint
from tclib import box_points, cube_monomials, mon_eval, nu_exact


def V_basis(pts, N, d):
    mons = [a for a in cube_monomials(d) if sum(a) < N]
    M = flint.fmpz_mat(len(mons), len(pts), [mon_eval(k, a) for a in mons for k in pts])
    X, nullity = M.nullspace()
    return [[int(X[i, j]) for i in range(len(pts))] for j in range(nullity)]


def plane_gcd(pts, v1, v2, N, basis):
    a = [sum(x * y for x, y in zip(k, v1)) for k in pts]
    b = [sum(x * y for x, y in zip(k, v2)) for k in pts]
    g = None
    for c in basis:
        poly = flint.fmpq_poly([0])
        for ck, ak, bk in zip(c, a, b):
            if ck:
                poly += ck * flint.fmpq_poly([ak, bk]) ** N
        if poly == 0:
            continue
        g = poly if g is None else g.gcd(poly)
        if g.degree() == 0:
            break
    return g


def analyse(pts, d, planes, label):
    out = []
    cache = {}
    for v1, v2 in planes:
        T = 1000003
        N, _ = nu_exact(pts, [x + T * y for x, y in zip(v1, v2)])
        if N == 0:
            continue
        if N not in cache:
            cache[N] = V_basis(pts, N, d)
        g = plane_gcd(pts, v1, v2, N, cache[N])
        bad = []
        if g is not None and g.degree() > 0:
            for f, e in g.factor()[1]:
                if f.degree() >= 2:
                    # real root?  use sympy for real root count
                    import sympy as sp
                    t = sp.symbols('t')
                    fs = sp.Poly([sp.Rational(int(f[i].p), int(f[i].q)) for i in range(f.degree(), -1, -1)], t)
                    nreal = len(sp.real_roots(fs))
                    bad.append((f.degree(), nreal))
        out.append((v1, v2, N, g.degree() if g is not None else None, bad))
        flag = 'VIOLATION' if any(nr > 0 for _, nr in bad) else ''
        print(label, 'v1=%s v2=%s nu(M)=%d gcd_deg=%s nonlinear_factors(deg,#real)=%s %s' % (v1, v2, N, out[-1][3], bad, flag), flush=True)
    return out


if __name__ == '__main__':
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
    nplanes = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    models = []
    # two-cut blocks (crossing), d = 5, 6
    for (a, b, c, L1, L2) in [(2, 1, 2, 1, 1), (2, 2, 2, 1, 1), (2, 2, 1, 2, 1), (1, 3, 2, 2, 1), (2, 1, 2, 2, 2)]:
        d = a + b + c
        S1 = list(range(a + b)); S2 = list(range(a, d))
        models.append(('twocut(%d,%d,%d;%d,%d)' % (a, b, c, L1, L2), d, box_points(d, [(S1, L1), (S2, L2)])))
    # M8 (validation: J4 says only rational roots can occur)
    fano = [[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]]
    cols = [[x for x in range(7) if x not in l] for l in fano]
    models.append(('M8', 7, box_points(7, [(c, 2) for c in cols])))
    for label, d, pts in models:
        planes = []
        for _ in range(nplanes):
            v1 = [rng.randint(-3, 3) for _ in range(d)]
            v2 = [rng.randint(-3, 3) for _ in range(d)]
            planes.append((v1, v2))
        analyse(pts, d, planes, label)
