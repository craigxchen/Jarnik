"""Exact toy check of Lemma 4.1: on Y = {k in {-1,0,1}^d : |k(X_i)| <= c}, the polynomial
Q = l_u^m - sum_{i<o} binom(m,i) eta^i R_i(k) l_v(k)^i (R_i = tensor-Lagrange interpolant of l_y^{m-i} in the
statistics t_i = k(X_i)) has leading form l_u^m (deg of subtracted part < m) and equals sum_{i>=o} binom(m,i) eta^i l_y^{m-i} l_v^i on Y.
Toy: d = 8, q = 2 bit-sets, c = 1, rho = 2cq = 4, m = 7, o = 3; exact rationals."""
import itertools, math
from fractions import Fraction as F
import sympy as sp
d, q, c, m = 8, 2, 1, 7
rho = 2 * c * q; o = m - rho
k = sp.symbols('k0:%d' % d)
X = [[0, 1, 2, 3], [0, 1, 4, 5]]           # two overlapping 'bit' sets
lam = [F(1, 3), F(2, 3)]                    # y = sum lam_i 1_{X_i}
y = [sum(lam[i] for i in range(q) if x in X[i]) for x in range(d)]
v = [F(x + 1, 10) for x in range(d)]
eta = F(1, 7)
ly = sum(sp.Rational(y[x].numerator, y[x].denominator) * k[x] for x in range(d))
lv = sum(sp.Rational(v[x].numerator, v[x].denominator) * k[x] for x in range(d))
t = [sum(k[x] for x in Xi) for Xi in X]
vals = list(range(-c, c + 1))
def lag(j, s):  # Lagrange basis on vals, as polynomial in expression s
    return sp.prod([(s - a) / sp.Integer(vals[j] - a) for a in vals if a != vals[j]])
def interp(e):  # interpolate F(t1,t2) = (lam.t)^e on the box
    tot = 0
    for idx in itertools.product(range(len(vals)), repeat=q):
        pt = [vals[i] for i in idx]
        Fv = sum(sp.Rational(lam[i].numerator, lam[i].denominator) * pt[i] for i in range(q)) ** e
        tot += Fv * sp.prod([lag(idx[i], t[i]) for i in range(q)])
    return sp.expand(tot)
etaS = sp.Rational(1, 7)
lu = ly + etaS * lv
Q = sp.expand(lu ** m - sum(sp.binomial(m, i) * etaS ** i * interp(m - i) * lv ** i for i in range(o)))
P = sp.Poly(Q - sp.expand(lu ** m), *k)
print("degree of Q - l_u^m:", P.total_degree(), "(< m =", m, ")")
from fractions import Fraction
def ev(expr_poly, pt):
    tot = Fraction(0)
    for mon, coef in expr_poly.terms():
        t = Fraction(int(coef.p), int(coef.q))
        for e, xv in zip(mon, pt):
            if e: t *= xv ** e
        tot += t
    return tot
QP = sp.Poly(Q, *k)
et = Fraction(1, 7)
ok = 0; tot = 0
for kk in itertools.product([-1, 0, 1], repeat=d):
    tv = [sum(kk[x] for x in Xi) for Xi in X]
    if max(abs(a) for a in tv) > c: continue
    tot += 1
    lyv = sum(y[x] * kk[x] for x in range(d)); lvv = sum(v[x] * kk[x] for x in range(d))
    rhs = sum(math.comb(m, i) * et ** i * lyv ** (m - i) * lvv ** i for i in range(o, m + 1))
    ok += (ev(QP, kk) == rhs)
print("identity holds on", ok, "of", tot, "points of Y (exact rationals)")
