# Referee's own exact library (independent of the author's lib.py).
# Psi^Sigma_A(u) = min_{deg P < m} sum_{k in A} (<k,u>^m - P(k))^2, exact over Q (python-flint fmpq_mat).
import itertools
from fractions import Fraction
import flint

def monos(d, maxdeg):
    """exponent vectors alpha in {0,1,2}^d with |alpha| <= maxdeg (enough on ternary sets)"""
    out = []
    for a in itertools.product(range(3), repeat=d):
        if sum(a) <= maxdeg:
            out.append(a)
    return out

def ev(k, a):
    v = 1
    for ki, ai in zip(k, a):
        v *= ki ** ai
    return v

def indep_cols(rows_as_cols, U):
    """rows_as_cols: list of column vectors (length U) ints; return indices of a maximal independent subset"""
    chosen = []
    basis = []  # echelon rows over Q as dict pivot->vector
    piv = []
    for idx, col in enumerate(rows_as_cols):
        v = [flint.fmpq(x) for x in col]
        for (p, b) in zip(piv, basis):
            if v[p] != 0:
                f = v[p] / b[p]
                v = [vi - f * bi for vi, bi in zip(v, b)]
        nz = [i for i, x in enumerate(v) if x != 0]
        if nz:
            piv.append(nz[0]); basis.append(v); chosen.append(idx)
    return chosen

def poly_basis(A, m):
    d = len(A[0])
    ms = monos(d, m - 1)
    cols = [[ev(k, a) for k in A] for a in ms]
    ch = indep_cols(cols, len(A))
    return [cols[i] for i in ch]

def psi_sigma(A, m, u, B=None):
    """u: list of Fractions/ints. returns exact fmpq residual sum of squares."""
    if B is None:
        B = poly_basis(A, m)
    U = len(A)
    g = []
    for k in A:
        s = sum(flint.fmpq(Fraction(ui).numerator, Fraction(ui).denominator) * ki for ki, ui in zip(k, u))
        g.append(s ** m)
    r = len(B)
    if r == 0:
        return sum(x * x for x in g)
    Bm = flint.fmpq_mat(U, r, [flint.fmpq(B[j][i]) for i in range(U) for j in range(r)])
    gv = flint.fmpq_mat(U, 1, g)
    G = Bm.transpose() * Bm
    b = Bm.transpose() * gv
    x = G.solve(b)
    res = (gv.transpose() * gv)[0, 0] - (b.transpose() * x)[0, 0]
    return res

def cheb_zero(A, m, u, B=None):
    return psi_sigma(A, m, u, B) == 0
