"""Referee checks (independent of Astra's ii_grid_check.py):
 (1) circulant identity det[1,z_i..z_{i+4}] = +- prod_{k=1..5} sum_j xi^{kj} z_j, generic symbols;
 (2) summation-by-parts identity;
 (3) 2x2 block determinant formula, generic symbols, d=1..5;
 (4) full modified simplex design determinant = predicted product, symbolic for d=1, exact random rationals d=1..3;
 (5) coefficient of z_j^5 in the Fourier product is nonzero.
"""
import random
from fractions import Fraction
from itertools import product
import sympy as sp

# (1)
z = sp.symbols('z0:6')
xi = sp.Rational(1, 2) + sp.sqrt(3) * sp.I / 2
M = sp.Matrix([[1] + [z[(i + h) % 6] for h in range(5)] for i in range(6)])
det = sp.expand(M.det(method='berkowitz'))
prod_f = sp.expand(sp.prod([sum(sp.expand(xi ** (k * j)) * z[j] for j in range(6)) for k in range(1, 6)]))
prod_f = sp.expand(sp.nsimplify(prod_f))
d1 = sp.simplify(sp.expand(det - prod_f))
d2 = sp.simplify(sp.expand(det + prod_f))
print('(1) det - prod == 0:', d1 == 0, '; det + prod == 0:', d2 == 0)
P = sp.Poly(prod_f, *z)
print('(5) coefficients of z_j^5:', [sp.simplify(P.coeff_monomial(z[j] ** 5)) for j in range(6)])
print('    number of monomials in det:', len(sp.Poly(det, *z).terms()))

# (2)
h = sp.symbols('h')
lhs = (1 - h) * sum(h ** j * z[j] for j in range(6))
rhs = sum((z[j] - z[j - 1]) * (h ** j - 1) for j in range(1, 6))
print('(2) LHS-RHS reduced mod h^6-1:', sp.rem(sp.expand(lhs - rhs), h ** 6 - 1, h))

# (3) and (4): Newton-basis modified simplex design
def total_degree(k, d):
    if k == 0:
        yield ()
        return
    for a in range(d + 1):
        for t in total_degree(k - 1, d - a):
            yield (a,) + t

def design(d, node, u, v):
    A = list(total_degree(5, d))
    special = [(d, 0, 0, 0, 0), (d - 1, 1, 0, 0, 0)]
    rows = [tuple(node[hh][b[hh]] for hh in range(5)) for b in A if b not in special]
    tail = (node[2][0], node[3][0], node[4][0])
    rows += [(u, v) + tail, (v, u) + tail]
    return A, rows

def predicted(d, node, u, v):
    A = list(total_degree(5, d))
    special = [(d, 0, 0, 0, 0), (d - 1, 1, 0, 0, 0)]
    val = 1
    npair = 0
    for b in A:
        if b in special:
            continue
        for hh in range(5):
            for j in range(b[hh]):
                val *= (node[hh][b[hh]] - node[hh][j]); npair += 1
    Pu = 1; Pv = 1
    for j in range(d - 1):
        Pu *= (u - node[0][j]); Pv *= (v - node[0][j]); npair += 2
    val *= Pu * Pv * (u - v) * (u + v - node[0][d - 1] - node[1][0])
    npair += 1
    return val, npair

for d in range(1, 6):
    a, b, u, v = sp.symbols('a b u v')
    s = sp.symbols('s0:%d' % max(d - 1, 1))
    Pf = lambda t: sp.prod([t - s[j] for j in range(d - 1)])
    # entries: col (d,0..): P(t)(t-a) at x1=t ; col (d-1,1,..): P(x1)(x2-b)
    B = sp.Matrix([[Pf(u) * (u - a), Pf(u) * (v - b)], [Pf(v) * (v - a), Pf(v) * (u - b)]])
    claim = Pf(u) * Pf(v) * (u - v) * (u + v - a - b)
    print('(3) d=%d block det - claim =' % d, sp.expand(B.det() - claim), ' (+claim:', sp.expand(B.det() + claim) == 0, ')')

def monomial_det(A, rows):
    import sympy
    Mx = sympy.Matrix([[sympy.prod([r[hh] ** al[hh] for hh in range(5)]) for al in A] for r in rows])
    return Mx

# symbolic d=1
d = 1
node = [[sp.Symbol('s%d_%d' % (hh, j)) for j in range(d + 1)] for hh in range(5)]
u, v = sp.symbols('u v')
A, rows = design(d, node, u, v)
D = sp.factor(monomial_det(A, rows).det())
pv, npair = predicted(d, node, u, v)
print('(4) d=1 symbolic factorization:', D)
print('    equals +-predicted:', sp.expand(D - pv) == 0 or sp.expand(D + pv) == 0, ' V=', sum(map(sum, A)), 'pair+1 factors=', npair + 1)

# exact random rational tests d=1..3
def frac_det(Mx):
    n = len(Mx); Mx = [row[:] for row in Mx]; det = Fraction(1)
    for c in range(n):
        piv = next((r for r in range(c, n) if Mx[r][c] != 0), None)
        if piv is None:
            return Fraction(0)
        if piv != c:
            Mx[c], Mx[piv] = Mx[piv], Mx[c]; det = -det
        det *= Mx[c][c]
        for r in range(c + 1, n):
            f = Mx[r][c] / Mx[c][c]
            if f:
                Mx[r] = [x - f * y for x, y in zip(Mx[r], Mx[c])]
    return det

random.seed(1)
for d in (1, 2, 3):
    for trial in range(3):
        vals = random.sample(range(-500, 500), 5 * (d + 1) + 2)
        node = [[Fraction(vals[5 * j + hh]) for j in range(d + 1)] for hh in range(5)]
        u, v = Fraction(vals[-2]), Fraction(vals[-1])
        A, rows = design(d, node, u, v)
        Mx = [[Fraction(1) * eval_mon(r, al) if False else
               __import__('math').prod([r[hh] ** al[hh] for hh in range(5)]) for al in A] for r in rows]
        D = frac_det(Mx)
        pv, npair = predicted(d, node, u, v)
        V = sum(map(sum, A))
        print('(4) d=%d trial %d: n=%d V=%d  det/pred = %s  linear factors=%d' % (d, trial, len(A), V, D / pv, npair + 1))
