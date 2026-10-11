"""Exact upper bound for kappa_1 = max{ ||lambda||_1 : lambda in R^M_0, N(lambda) <= 1 },
N(lambda) = sum_S om_S |lambda(S)|.

||lambda||_1 = max over s in {+-1}^M of s.lambda, and by LP duality, for each s,
   max{ s.lambda : N(lambda) <= 1, sum lambda = 0 } = min{ max_S |alpha_S|/om_S :
                                sum_S alpha_S 1_S(x) = s_x + c for all x, some c in R }.
For each s (s and -s give the same value; s_0 = +1 wlog) we take a float optimal alpha from HiGHS,
rationalise it, repair it to satisfy the equality system exactly (exact particular solution of the
residual system), and evaluate max|alpha_S|/om_S exactly: a valid upper bound for that s.
The maximum over s is an exact upper bound for kappa_1."""
import itertools
from fractions import Fraction
import numpy as np
from scipy.optimize import linprog


def _exact_particular(Aint, b):
    """Exact solution x (Fractions) of A x = b (A integer m x n, full row rank m), or None."""
    m = len(Aint); n = len(Aint[0])
    R = [[Fraction(v) for v in row] + [Fraction(bb)] for row, bb in zip(Aint, b)]
    piv = []
    r = 0
    for c in range(n):
        p = next((i for i in range(r, m) if R[i][c] != 0), None)
        if p is None:
            continue
        R[r], R[p] = R[p], R[r]
        pv = R[r][c]
        R[r] = [v / pv for v in R[r]]
        for i in range(m):
            if i != r and R[i][c] != 0:
                f = R[i][c]
                R[i] = [a - f * bq for a, bq in zip(R[i], R[r])]
        piv.append(c)
        r += 1
        if r == m:
            break
    for i in range(r, m):
        if R[i][-1] != 0:
            return None
    x = [Fraction(0)] * n
    for i, c in enumerate(piv):
        x[c] = R[i][-1]
    return x


def kappa1_bound(M, T, om):
    om = [Fraction(w) for w in om]
    nT = len(T)
    # system matrix over unknowns (alpha_S..., c):  sum_S alpha_S 1_S(x) - c = s_x
    A = [[(T[j] >> x) & 1 for j in range(nT)] + [-1] for x in range(M)]
    Af = np.array(A, dtype=float)
    best = Fraction(0)
    for tail in itertools.product((1, -1), repeat=M - 1):
        s = (1,) + tail
        if sum(s) in (M, -M):
            continue          # constant s: s.lambda = 0
        # LP: min t s.t. -t om_S <= alpha_S <= t om_S, A (alpha, c) = s
        nv = nT + 2
        c = np.zeros(nv); c[-1] = 1
        Aub = []
        for j in range(nT):
            r1 = np.zeros(nv); r1[j] = 1; r1[-1] = -float(om[j]); Aub.append(r1)
            r2 = np.zeros(nv); r2[j] = -1; r2[-1] = -float(om[j]); Aub.append(r2)
        Aeq = np.hstack([Af, np.zeros((M, 1))])
        res = linprog(c, A_ub=np.array(Aub), b_ub=np.zeros(2 * nT), A_eq=Aeq, b_eq=np.array(s, dtype=float),
                      bounds=[(None, None)] * (nT + 1) + [(0, None)], method='highs')
        if res.status != 0:
            return None
        al = [Fraction(float(v)).limit_denominator(10 ** 6) for v in res.x[:nT + 1]]
        resid = [Fraction(s[x]) - sum(A[x][j] * al[j] for j in range(nT + 1)) for x in range(M)]
        corr = _exact_particular(A, resid)
        if corr is None:
            return None
        al = [a + d for a, d in zip(al, corr)]
        assert all(sum(A[x][j] * al[j] for j in range(nT + 1)) == s[x] for x in range(M))
        val = max(abs(al[j]) / om[j] for j in range(nT))
        if val > best:
            best = val
    return best
