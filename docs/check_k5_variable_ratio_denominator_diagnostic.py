"""K5 adjacent-ratio computation over the rational function field QQ(t).

Fix x01=1 and x02=t. The four linear vertex-sum differences are solved
directly from the incidence equations, with free variables (a,b,c,d). A
completed grevlex Groebner basis over QQ(t) is then used to reduce the nine
successive factors of P9 = product(j=1..9)(x01^2-x_j^2). The generic-field
identity P9=0 is exact over QQ(t); specialization at poles of the rational
certificate requires separate checking. The quotient-denominator output is a
pole diagnostic only: Buchberger transformation coefficients back to the
original generators are not tracked here.
"""

import sympy as s
from itertools import combinations

t = s.symbols('t')
a, b, c, d = s.symbols('a b c d')
p03, p04, p12, p13 = s.symbols('p03 p04 p12 p13')
pivots = (p03, p04, p12, p13)
x = {
    (0, 1): s.Integer(1), (0, 2): t,
    (0, 3): p03, (0, 4): p04, (1, 2): p12, (1, 3): p13,
    (1, 4): a, (2, 3): b, (2, 4): c, (3, 4): d,
}
edges = list(combinations(range(5), 2))
linear = [sum(x[e] for e in edges if v in e)
          - sum(x[e] for e in edges if 0 in e) for v in range(1, 5)]
M, _ = s.linear_eq_to_matrix(linear, pivots)
assert M.det() == 1
sol = s.solve(linear, pivots, dict=True)[0]
for e in x:
    x[e] = s.expand(x[e].subs(sol))
vertex_sums = [sum(x[e] for e in edges if v in e) for v in range(5)]
assert all(s.expand(vertex_sums[v] - vertex_sums[0]) == 0 for v in range(1, 5))
assert x[(0, 3)] == a + c + d - t - 1
assert x[(0, 4)] == -a + 2*b + c + d - 2
assert x[(1, 2)] == b + c + 2*d - t - 2
assert x[(1, 3)] == -a + b + c + t - 1

K = s.QQ.frac_field(t)
F = [s.Poly(s.expand(sum(x[e]**3 for e in edges if v in e)
             - sum(x[e]**3 for e in edges if 0 in e)),
             a, b, c, d, domain=K) for v in range(1, 5)]
print('linear_det=', M.det(), 'domain=QQ(t) order=grevlex', flush=True)
G = s.groebner(F, a, b, c, d, order='grevlex', domain=K, method='f5b')
print('completed=', G.is_zero_dimensional, 'basis_len=', len(G.polys),
      'domain=', G.domain, flush=True)
assert G.is_zero_dimensional and len(G.polys) == 28

Qden = s.Integer(1)
prod = s.Integer(1)
for j in range(1, 10):
    prod = s.expand(prod * (x[edges[0]]**2 - x[edges[j]]**2))
    quot, rem = G.reduce(s.Poly(prod, a, b, c, d, domain=K).as_expr())
    for qi in quot:
        qexpr = qi.as_expr() if hasattr(qi, 'as_expr') else qi
        for cc in s.Poly(qexpr, a, b, c, d, domain=K).coeffs():
            Qden = s.lcm(Qden, s.denom(s.cancel(cc)))
    print('prefix', j, 'zero=', rem == 0,
          'terms=', len(rem.terms()) if rem != 0 else 0, flush=True)
    assert (j < 9 and rem != 0) or (j == 9 and rem == 0)
    if rem != 0:
        prod = rem.as_expr()
print('quotient_denominator=', s.factor(Qden), 'degree=', s.degree(Qden,t), flush=True)
