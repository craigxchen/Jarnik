"""Exact checks for the critical-arc example and identities (II.3)–(II.9)."""

from functools import reduce
from itertools import combinations
from math import prod

import sympy as sp

x = sp.symbols("x", real=True)
marks = range(4)
pairs = list(combinations(marks, 2))
offsets = (0, 1, 3, 7)
linear = [x + c + sp.I for c in offsets]
stars = {
    pair: sp.Poly(prod(f if j in pair else sp.conjugate(f)
                       for j, f in enumerate(linear)), x, extension=sp.I)
    for pair in pairs[:-1]
}
base = next(iter(stars.values()))
assert all((v - base).degree() <= 2 for v in stars.values())
assert len({sp.im(v.nth(2)) for v in stars.values()}) == 5
assert reduce(sp.gcd, stars.values()).degree() == 0
assert stars[0, 2]*stars[1, 3] == stars[0, 3]*stars[1, 2]
for v in stars.values():
    assert sp.expand_complex(v.as_expr()*sp.conjugate(v.as_expr())
                             - base.as_expr()*sp.conjugate(base.as_expr())).expand() == 0
print("Five quartics have equal norm, critical separation, gcd one and an exact relation.")

marks = range(5)
pairs = list(combinations(marks, 2))
blocks = [sp.Integer(a) + sp.I*b for a, b in
          [(2, 1), (3, 2), (4, 1), (5, 2), (6, 1),
           (5, 4), (7, 2), (6, 5), (8, 3), (8, 5)]]
n = {edge: int(sp.expand_complex(v*sp.conjugate(v))) for edge, v in zip(pairs, blocks)}
P = [sp.expand(prod(v for edge, v in zip(pairs, blocks) if i in edge)) for i in marks]
H = prod(n.values())
u = [sp.cancel(v/sp.conjugate(v), extension=sp.I) for v in P]
def norm(i, j):
    return n[tuple(sorted((i, j)))]

s, t = {}, {}
for i in marks:
    for j in marks:
        if i != j:
            word = sp.expand(sp.conjugate(P[i])*P[j])/norm(i, j)
            s[i, j], t[i, j] = map(sp.simplify, word.as_real_imag())
            assert s[i, j].is_Integer and t[i, j].is_Integer and t[i, j] != 0

# Factorwise cancellation keeps the Vandermonde product small.
J = sp.S.One
for i, j in pairs:
    J = sp.cancel(J*(u[i]-u[j])/(2*sp.I), extension=sp.I)
for v in u:
    J = sp.cancel(J/v**2, extension=sp.I)
assert J == sp.Rational(prod(t[i, j] for i, j in pairs), H**3)

for i, j in pairs:
    D = prod(norm(i, k)*norm(j, k) for k in marks if k not in (i, j))
    assert D == s[i, j]**2 + t[i, j]**2
    for k in set(marks) - {i, j}:
        l, m = sorted(set(marks) - {i, j, k})
        A = (t[j, k]**2*norm(j, k)*norm(i, l)*norm(i, m)
             + t[i, k]**2*norm(i, k)*norm(j, l)*norm(j, m)
             - t[i, j]**2*norm(i, j)*norm(k, l)*norm(k, m))
        assert A == 2*t[i, k]*t[j, k]*s[i, j]
        assert 4*t[i, k]**2*t[j, k]**2*D - A**2 == 4*(t[i, j]*t[i, k]*t[j, k])**2
print("Cubic-cover identity, six-norm identity and all thirty sextic identities verified.")
adj = sp.Matrix([[int(e != f and bool(set(e) & set(f))) for e in pairs] for f in pairs])
assert adj.det() != 0 and all(c.is_Integer for c in 6*adj.inv())
print("The edge adjacency determinant is", adj.det(), "and 6 times its inverse is integral.")

# A genuine ten-term prime progression. Exact trial division suffices here.
from math import isqrt

ap_norms = [34913 + 2100*r for r in range(10)]
ap_blocks = []
for n_e in ap_norms:
    assert n_e % 4 == 1 and all(n_e % j for j in range(2, isqrt(n_e) + 1))
    x_e = next(j for j in range(1, isqrt(n_e)+1)
               if isqrt(n_e-j*j)**2 == n_e-j*j)
    ap_blocks.append(x_e + sp.I*isqrt(n_e-x_e*x_e))
ap_stars = [sp.expand(prod(z for edge, z in zip(pairs, ap_blocks) if i in edge))
            for i in marks]
for edge, n_e in zip(pairs, ap_norms):
    i, j = edge
    determinant = sp.im(sp.expand(sp.conjugate(ap_stars[i])*ap_stars[j]))
    assert determinant != 0 and determinant % n_e == 0
assert sp.Matrix([(1, 1, 0), (1, 5, 6), (1, 9, 20)]).det() == 32
print("Verified a squarefree ten-term progression and the generic quartet determinant.")
