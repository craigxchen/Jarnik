"""T2 negative control at k=4 (gauss Thm 4.1 instance X=-3, m=3).
(1) 4-point Cox form: gamma_T(s) = s - r_T (T subset [3]), Pi_x = prod_{T ni x} gamma_T, Im Pi_x(s) constant (exact, sympy).
(2) Along integer s: primitive reduction, residues t_ab = Im Q_ab bounded, singleton decoupling (gamma_x is the unique
    solution of Im(gamma R_x) = t_x near the box), balance.
(3) Pencil test: the tight sections at k=4 with the support of V_4 (Vandermonde, order 6) and of pair products:
    dim A_t(S) = 1 for S = all lattice points of the Newton polytope (Theorem R at n=3) -> no tight pencil exists,
    so the pencil principle (Theorem S) says nothing at k=4, consistent with the dense family.
"""
import sympy as sp
from itertools import combinations, product
from gauss_lib import *
import rigidity

s = sp.symbols('s')
I = sp.I
roots = {frozenset([1, 2]): -6-4*I, frozenset([1, 3]): -6-3*I, frozenset([2, 3]): -6+5*I,
         frozenset([1]): -7+6*I, frozenset([2]): -3-2*I, frozenset([3]): -4-3*I, frozenset([1, 2, 3]): I}
CUTS4 = list(roots.keys())
print('(1) Im Pi_x(s):')
for x in [1, 2, 3]:
    P = sp.expand(sp.prod([s - roots[T] for T in CUTS4 if x in T]))
    print('   x=%d  Im P_x = %s' % (x, sp.simplify(sp.im(sp.expand(P.subs(s, sp.Symbol('u', real=True)))))))

def gi(z):
    z = complex(z)
    return (int(round(z.real)), int(round(z.imag)))

print('(2) along s:')
maxt = 0; nrows = 0; dec_ok = True
for sv in [10**k + 7 for k in range(2, 13)] + [123456789012345678901]:
    g = {T: gsub((sv, 0), gi(roots[T])) for T in CUTS4}
    zs = []
    for x in [1, 2, 3]:
        z = (1, 0)
        for T in CUTS4:
            z = gmul(z, g[T] if x in T else gconj(g[T]))
        zs.append(z)
    z0 = gprod([gconj(g[T]) for T in CUTS4])
    zs = zs + [z0]
    G = ggcd_list(zs)
    red = [gexactdiv(z, G) for z in zs]
    N = gnorm(red[0])
    # residues |t_ab| from |zeta_a - zeta_b|^2 = 4 N t^2 / e^{d_ab}; e^{d_ab} = |c_ab|^2/(4 t^2) ... use c_ab:
    ts = []
    for a, b in combinations(range(4), 2):
        c = gexactdiv(gsub(red[a], red[b]), ggcd(red[a], red[b]))
        ts.append(gnorm(c))  # = 4 t_ab^2 within one unit class, else 2*(odd)
    # singleton decoupling: Im(gamma_x R_x) = t_x with R_x = prod_{T ni x, |T|>=2} gamma_T
    for x in [1, 2, 3]:
        R = gprod([g[T] for T in CUTS4 if x in T and len(T) >= 2])
        gx = g[frozenset([x])]
        tx = gmul(gx, R)[1]
        # all solutions (u,v) of u*Im R + v*Re R = tx lie on a line with spacing |R|; the one within |R|/2 of 0 is unique
        if not (abs(gx[0]) + abs(gx[1]) < (gnorm(R))**0.5 / 2):
            dec_ok = False
    nrows += 1
    maxt = max(maxt, max(ts))
    print('   s=%-22d gcd norm=%-10d max|c_ab|^2=%-6d log N/log s=%.3f' % (sv, gnorm(G), max(ts), sp.log(N, 10).evalf()/sp.log(sv, 10).evalf()))
print('   singleton decoupling (gamma_x unique short solution): %s' % dec_ok)

print('(3) tight-section rigidity at k=4 (n=3 coordinates after dropping the base coordinate):')
def lattice_points_of_zonotope(gens, n):
    pts = {tuple([0]*n)}
    for gv in gens:
        pts |= {tuple(p[i]+gv[i] for i in range(n)) for p in pts}
    # fill: all lattice points of the convex hull (box filter + LP-free test via sympy convex hull is heavy; use
    # the sumset plus integer points of the bounding box that are in the hull by a half-space test from scipy)
    import numpy as np
    from scipy.spatial import ConvexHull
    P = np.array(sorted(pts))
    if len(P) <= n or np.linalg.matrix_rank(P - P[0]) < n:
        return sorted(pts)
    H = ConvexHull(P)
    lo = P.min(0); hi = P.max(0)
    out = []
    for q in product(*[range(lo[i], hi[i]+1) for i in range(n)]):
        if all(np.dot(eq[:-1], q) + eq[-1] <= 1e-9 for eq in H.equations):
            out.append(q)
    return out
def rv(a, b, k=4):
    v = [0]*k; v[a] += 1; v[b] -= 1
    return tuple(v[:k-1])
pairs4 = [(a, b) for a in range(4) for b in range(a+1, 4)]
tests = {'V_4 (all 6 pairs)': [rv(a, b) for a, b in pairs4],
         'V_4 squared': [rv(a, b) for a, b in pairs4]*2,
         '(12)(34)': [rv(0, 1), rv(2, 3)], '(12)(13)(23)': [rv(0, 1), rv(0, 2), rv(1, 2)]}
for name, gens in tests.items():
    S = lattice_points_of_zonotope(gens, 3)
    K = rigidity.MW(S, 3)
    print('   %-16s |S|=%-4d MW=%s reg=%d dim A_MW=%d' % (name, len(S), K, rigidity.reg(S, 3), rigidity.dimA(S, 3, int(K))))
