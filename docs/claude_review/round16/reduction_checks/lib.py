"""Core exact routines for round 16 (reduction of MIN^pair to L2 / anti-concentration form).

Conventions: points k' in {-1,0,1}^d (base point n has u_n = 0, k_n = -sum k' unconstrained at T = 1).
A box is a list of tuples. Psi^Sigma_A(u) = min_{deg P < m} sum_{k in A} (<k,u>^m - P(k))^2  (counting measure).
Everything exact over Q via python-flint.
"""
import itertools
from fractions import Fraction
from math import comb, factorial
import flint


def monomials(d, maxdeg):
    """exponent tuples alpha in N^d with |alpha| <= maxdeg"""
    out = []
    def rec(i, rem, cur):
        if i == d:
            out.append(tuple(cur)); return
        for a in range(rem + 1):
            cur.append(a); rec(i + 1, rem - a, cur); cur.pop()
    rec(0, maxdeg, [])
    return out


def mono_eval(k, alpha):
    v = 1
    for kx, a in zip(k, alpha):
        if a:
            v *= kx ** a
    return v


def poly_basis_matrix(A, m):
    """flint fmpq_mat with rows = points, columns = independent monomials of degree < m on A (pivot set)."""
    d = len(A[0])
    mons = [a for a in monomials(d, m - 1)]
    # reduce: on ternary points k^a depends only on a mod: a>=3 -> a-2 parity pattern; keep a_x<=2 suffices
    mons = [a for a in mons if max(a, default=0) <= 2]
    M = flint.fmpq_mat(len(mons), len(A), [mono_eval(k, a) for a in mons for k in A])
    R, rank = M.rref()
    # pivot columns of M^T  == independent monomials: compute rref of M^T instead
    MT = M.transpose()
    R2, rank2 = MT.rref()
    piv = []
    r = 0
    for c in range(MT.ncols()):
        if r < rank2 and R2[r, c] != 0:
            piv.append(c); r += 1
    B = flint.fmpq_mat(len(A), len(piv), [mono_eval(k, mons[c]) for k in A for c in piv])
    return B, [mons[c] for c in piv]


class Box:
    def __init__(self, A, m):
        self.A = [tuple(k) for k in A]
        self.m = m
        self.d = len(self.A[0])
        self.B, self.piv = poly_basis_matrix(self.A, m)
        Bt = self.B.transpose()
        self.G = Bt * self.B
        self.Ginv = self.G.inv()
        self.r0 = len(self.piv)

    def g(self, u):
        m = self.m
        vals = [sum(Fraction(ux) * kx for ux, kx in zip(u, k)) ** m for k in self.A]
        return flint.fmpq_mat(len(self.A), 1, [flint.fmpq(v.numerator, v.denominator) for v in vals])

    def residual(self, u):
        g = self.g(u)
        b = self.B.transpose() * g
        coef = self.Ginv * b
        return g - self.B * coef

    def psi(self, u):
        """counting-measure Psi^Sigma_A(u) (exact Fraction)"""
        res = self.residual(u)
        s = (res.transpose() * res)[0, 0]
        return Fraction(int(s.p), int(s.q))


def cube(d):
    return list(itertools.product((-1, 0, 1), repeat=d))


def one_cut(d, S, L):
    return [k for k in cube(d) if abs(sum(k[x] for x in S)) <= L]


def psi_cube_formula(u, m):
    """counting-measure Psi on {-1,0,1}^d: 3^d sum_{alpha<=2,|alpha|=m} (m!/alpha!)^2 (2/3)^#1 (2/9)^#2 u^{2alpha}"""
    d = len(u)
    tot = Fraction(0)
    for alpha in itertools.product((0, 1, 2), repeat=d):
        if sum(alpha) != m:
            continue
        c = Fraction(factorial(m))
        w = Fraction(1)
        for a, ux in zip(alpha, u):
            c /= factorial(a)
            if a == 1:
                w *= Fraction(2, 3) * Fraction(ux) ** 2
            elif a == 2:
                w *= Fraction(2, 9) * Fraction(ux) ** 4
        tot += c * c * w
    return tot * 3 ** d
