"""Independent check of the 'sixteen-point' two-moment claim used in VI.5.

Claim: if ten complex numbers r_e (e an edge of K5) have equal first and equal
second star power sums at all five vertices, then they take at most 4 distinct values.

Method (independent of Astra's Jacobian certificate):
 1. normalize first star sums to 0 (shift by S/4): r lies in V = ker(incidence), dim 5.
 2. four homogeneous quadrics Q_i - Q_5 = 0 on V ~ C^5.
 3. Groebner basis (sympy, grevlex) of the homogeneous ideal; Hilbert polynomial via
    counting standard monomials in degrees up to 12: if constant = 16 then the projective
    scheme is zero-dimensional of degree 16 (no positive-dimensional component).
 4. affine chart via a random linear form = 1, lex basis -> solve numerically, count
    distinct values among the ten r_e for each solution.
"""
from itertools import combinations
import random, sys
import sympy as sp

edges = list(combinations(range(5), 2))
M = sp.Matrix([[int(i in e) for e in edges] for i in range(5)])
K = M.nullspace()
assert len(K) == 5
xs = sp.symbols('x0:5')
r = [sum(xs[j]*K[j][e] for j in range(5)) for e in range(10)]
r = [sp.nsimplify(sp.expand(v)) for v in r]
Q = [sp.expand(sum(r[k]**2 for k, e in enumerate(edges) if i in e)) for i in range(5)]
eqs = [sp.expand(Q[i]-Q[4]) for i in range(4)]
G = sp.groebner(eqs, *xs, order='grevlex')
lead = [sp.Poly(g, *xs).monoms(order='grevlex')[0] for g in G.exprs]
print('grevlex basis size', len(G.exprs), 'max deg', max(sum(m) for m in lead))

def divisible(m, l):
    return all(a >= b for a, b in zip(m, l))

def monos(deg, n=5):
    if n == 1:
        yield (deg,)
        return
    for a in range(deg+1):
        for rest in monos(deg-a, n-1):
            yield (a,)+rest

for d in range(0, 13):
    cnt = sum(1 for m in monos(d) if not any(divisible(m, l) for l in lead))
    print('degree', d, 'standard monomials', cnt)

# affine chart
random.seed(1)
c = [random.randint(1, 50) for _ in range(5)]
chart = sum(ci*xi for ci, xi in zip(c, xs)) - 1
GL = sp.groebner(eqs+[chart], *xs, order='lex')
print('lex basis size', len(GL.exprs))
last = sp.Poly(GL.exprs[-1], xs[-1])
print('univariate degree in x4:', last.degree())
sols = sp.solve_poly_system(eqs+[chart], *xs)
print('number of affine solutions', len(sols))
maxdistinct = 0
for s in sols:
    vals = [complex(sp.N(v.subs(dict(zip(xs, s))), 30)) for v in r]
    dist = []
    for v in vals:
        if not any(abs(v-w) < 1e-12 for w in dist):
            dist.append(v)
    maxdistinct = max(maxdistinct, len(dist))
    print(len(dist), end=' ')
print('\nmax distinct values over all solutions:', maxdistinct)
