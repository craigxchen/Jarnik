"""Exact QQ fixed-ratio K5 reduction.

For x01=1, x02=2, this computes a completed SymPy Groebner basis for the
four vertex cube-difference cubics in variables (a,b,c,d), using QQ and
grevlex order, then reduces P9 and vertex square-sum differences. This is
only a fixed-ratio result and does not address arbitrary edge ratios.
"""

import sympy as s
from itertools import combinations

A, B, C, D = s.symbols('a b c d')
xs = {
    (0, 1): s.Integer(1), (0, 2): s.Integer(2),
    (0, 3): -3 + A + C + D,
    (0, 4): -2 - A + 2*B + C + D,
    (1, 2): -4 + B + C + 2*D,
    (1, 3): 1 - A + B + C,
    (1, 4): A, (2, 3): B, (2, 4): C, (3, 4): D,
}
edges = list(combinations(range(5), 2))
# Verify the displayed substitution and completeness of the linear solve.
vertex_sums = [sum(xs[e] for e in edges if v in e) for v in range(5)]
assert all(s.expand(vertex_sums[v] - vertex_sums[0]) == 0 for v in range(1, 5))
# Pivot variables are (x03,x04,x12,x13); free variables are (x14,x23,x24,x34).
# Derive the coefficients from graph incidence, including cancellation
# when an edge meets both compared vertices.
pivot_edges = ((0, 3), (0, 4), (1, 2), (1, 3))
pivot_matrix = s.Matrix([
    [int(v in edge) - int(0 in edge) for edge in pivot_edges]
    for v in range(1, 5)
])
assert pivot_matrix.det() == 1
print('linear_differences_zero=True pivot_det=', pivot_matrix.det(), flush=True)
F = []
for v in range(1, 5):
    F.append(s.expand(sum(xs[e]**3 for e in edges if v in e)
                     - sum(xs[e]**3 for e in edges if 0 in e)))

print('domain=QQ order=grevlex variables=(a,b,c,d) method=f5b', flush=True)
G = s.groebner(F, A, B, C, D, order='grevlex', domain=s.QQ, method='f5b')
print('zero_dimensional=', G.is_zero_dimensional, 'basis_len=', len(G.polys),
      'domain=', G.domain, 'order=', G.order, flush=True)
assert G.is_zero_dimensional and len(G.polys) == 28

P9 = s.Integer(1)
for j in range(1, 10):
    P9 *= xs[edges[0]]**2 - xs[edges[j]]**2
_, rem = G.reduce(s.Poly(s.expand(P9), A, B, C, D, domain=s.QQ).as_expr())
print('P9 remainder=', s.factor(rem), 'zero=', s.expand(rem) == 0, flush=True)
assert s.expand(rem) == 0

for v in range(1, 5):
    sq = s.expand(sum(xs[e]**2 for e in edges if v in e)
                  - sum(xs[e]**2 for e in edges if 0 in e))
    _, rem = G.reduce(sq)
    print('square_diff_vertex', v, 'zero=', s.expand(rem) == 0,
          'remainder=', s.factor(rem), flush=True)
