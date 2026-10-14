"""fatpoint library (round 13).  T = 1 slab boxes, tangent spaces, product (generalised zonotope) certificates.

Conventions (as round12/fake_checks/flib.py, copied, not imported, since the repository is read-only):
rows (points) 0..M-1, base point 0, coordinates of k' are rows 1..M-1 (index x-1).  A sign column is a pattern
S (set of rows); S and its complement give the same slab; relative pattern P omits the base row, and on Z^M_0 the
slab reads |k(P)| <= L.  A = {k' in {-1,0,1}^d : |k(P)| <= L for all patterns}.

Tangent space.  For c in V(A,m) (moments of degree < m vanish) Phi_c(u) = sum_k c_k <k,u>^m.  Its coefficient on u^a
is multinom(m;a) * mu_a(c), mu_a(c) = sum_k c_k k^a.  Only reduced exponents a in {0,1,2}^d occur.  T_m = {mu(c)}.
Product certificate: P(u) = prod_i <z_i,u> lies in span{Phi_c} iff p_a/multinom(m;a) (a reduced, and p_a = 0 for
non-reduced a) lies in T_m.
"""
import itertools
from math import factorial
import numpy as np
from flint import nmod_mat, fmpz_mat, fmpq_mat, fmpq

PR = 2147475331  # prime used in sibling round-13 scripts


def legendre(a, q):
    a %= q
    return 0 if a == 0 else (1 if pow(a, (q - 1) // 2, q) == 1 else -1)


def paley_I(q):
    idx = ['inf'] + list(range(q))
    M = q + 1
    H = np.zeros((M, M), dtype=int)
    for a, x in enumerate(idx):
        for b, y in enumerate(idx):
            if x == 'inf':
                H[a, b] = 1
            elif y == 'inf':
                H[a, b] = -1
            else:
                H[a, b] = (1 if x == y else 0) + legendre(y - x, q)
    return H


def canon(M, S):
    S = frozenset(S)
    if 0 in S:
        S = frozenset(range(M)) - S
    return tuple(sorted(S))


def sdc_patterns(q):
    K = paley_I(q)
    M = q + 1
    if not ((K + K.T) == 2 * np.eye(M, dtype=int)).all():
        K = K.T
    assert ((K + K.T) == 2 * np.eye(M, dtype=int)).all()
    Dg = np.diag(K[:, 0])
    K = Dg @ K @ Dg
    pats = set()
    for a in range(1, M):
        S = [x for x in range(M) if K[x, a] == -1]
        Sp = sorted(set(S) ^ {0, a})
        pats.add(canon(M, S))
        pats.add(canon(M, Sp))
    return M, sorted(pats)


def pal_patterns(q):
    H = paley_I(q)
    M = q + 1
    Hn = H.copy()
    Hn[1:, :] *= -1
    lab = {0: 1}
    for i in range(q):
        lab[1 + i] = 1 + i
    pats = set()
    for a in range(1, M):
        S = tuple(x for x in range(M) if Hn[x, a] == 1)
        if 0 < len(S) < M:
            pats.add(canon(M, S))
    for x in range(M):
        c = Hn[:, lab[x]].copy()
        c[x] *= -1
        S = tuple(y for y in range(M) if c[y] == 1)
        if 0 < len(S) < M:
            pats.add(canon(M, S))
    return M, sorted(pats)


def sylvester_patterns(k):
    M = 2 ** k
    pats = set()
    for a in range(1, M):
        S = [x for x in range(M) if bin(a & x).count('1') % 2 == 1]
        pats.add(canon(M, S))
    return M, sorted(pats)


def paley12_patterns():
    """Paley-12 of zetc §5.5 / loccompute: coordinates Z_11, base infinity, columns i + D, D = {0,1,3,4,5,9}.
    Rows: 0 = infinity (base), row 1+i = i in Z_11."""
    D = [0, 1, 3, 4, 5, 9]
    M = 12
    pats = set()
    for i in range(11):
        S = [1 + ((i + t) % 11) for t in D]
        pats.add(canon(M, S))
    return M, sorted(pats)


MODELS = {
    'sdc7': lambda: sdc_patterns(7), 'pal7': lambda: pal_patterns(7), 'syl8': lambda: sylvester_patterns(3),
    'sdc11': lambda: sdc_patterns(11), 'pal11': lambda: pal_patterns(11), 'paley12': paley12_patterns,
    'syl16': lambda: sylvester_patterns(4),
}


def rel(P):
    return [x - 1 for x in P]


def box(M, pats, L):
    d = M - 1
    cols = [rel(P) for P in pats]
    inc = [[j for j, C in enumerate(cols) if x in C] for x in range(d)]
    rem = [[sum(1 for y in C if y >= x) for C in cols] for x in range(d + 1)]
    out = []
    k = [0] * d
    sums = [0] * len(cols)

    def rec(x):
        if x == d:
            out.append(tuple(k))
            return
        for v in (-1, 0, 1):
            ok = True
            for j in inc[x]:
                sums[j] += v
            for j in range(len(cols)):
                if abs(sums[j]) - rem[x + 1][j] > L:
                    ok = False
                    break
            if ok:
                k[x] = v
                rec(x + 1)
            for j in inc[x]:
                sums[j] -= v
        k[x] = 0

    rec(0)
    return out


def reduced_monomials(d, deg):
    res = []

    def rec(x, left, cur):
        if x == d:
            if left == 0:
                res.append(tuple(cur))
            return
        if left > 2 * (d - x):
            return
        for e in (0, 1, 2):
            if e <= left:
                cur.append(e)
                rec(x + 1, left - e, cur)
                cur.pop()

    rec(0, deg, [])
    return res


def mono_val(k, a):
    v = 1
    for kx, ax in zip(k, a):
        if ax:
            if kx == 0:
                return 0
            if ax == 1:
                v *= kx
            # ax == 2: kx^2 = 1
    return v


def eval_matrix(A, monos):
    """rows = monomials, cols = points; entries in {-1,0,1}."""
    E = np.zeros((len(monos), len(A)), dtype=np.int64)
    Aarr = np.array(A, dtype=np.int64)
    for r, a in enumerate(monos):
        col = np.ones(len(A), dtype=np.int64)
        for x, ax in enumerate(a):
            if ax == 1:
                col *= Aarr[:, x]
            elif ax == 2:
                col *= Aarr[:, x] ** 2
        E[r] = col
    return E


def to_nmod(Mnp, p=PR):
    r, c = Mnp.shape
    return nmod_mat(r, c, [int(x) % p for x in Mnp.flatten()], p)


def nullspace_mod(Mnp, p=PR):
    """basis of right kernel of Mnp mod p, as numpy int64 array (rows = basis vectors)."""
    N = to_nmod(Mnp, p)
    X, nullity = N.nullspace()
    cols = [[int(X[i, j]) for i in range(X.nrows())] for j in range(nullity)]
    return np.array(cols, dtype=np.int64).reshape(nullity, Mnp.shape[1])


def rank_mod(Mnp, p=PR):
    return to_nmod(Mnp, p).rank()


def multinom(m, a):
    v = factorial(m)
    for ax in a:
        v //= factorial(ax)
    return v


def tangent_space_mod(A, m, p=PR):
    """returns (monos_m, Tbasis mod p) with Tbasis rows spanning T_m = {mu(c) : c in V(A,m)} mod p.
    NOTE: mod p this is the image of the mod-p kernel; its dimension is >= the rational H_A(m) only if ranks agree;
    used for search, never as a certificate."""
    d = len(A[0])
    low = []
    for e in range(m):
        low += reduced_monomials(d, e)
    monos_m = reduced_monomials(d, m)
    Elow = eval_matrix(A, low)
    V = nullspace_mod(Elow, p)  # rows: c vectors
    Em = eval_matrix(A, monos_m)
    # Em has entries in {-1,0,1}; V entries < p < 2^31, so V @ Em.T has |entries| < |A| 2^31: no int64 overflow
    T = (V % p) @ Em.T % p if V.shape[0] else np.zeros((0, len(monos_m)), dtype=np.int64)
    assert len(A) < 2 ** 31
    # row-reduce T
    Tn = to_nmod(T, p)
    Tn = Tn.rref()[0]
    rk = Tn.rank()
    Tb = np.array([[int(Tn[i, j]) for j in range(Tn.ncols())] for i in range(rk)], dtype=np.int64)
    return monos_m, Tb


def product_coeffs(Z, m):
    """Z: list of m integer vectors (length d).  Returns dict exponent-tuple -> integer coefficient of prod <z_i,u>."""
    d = len(Z[0])
    poly = {tuple([0] * d): 1}
    for z in Z:
        new = {}
        nz = [(x, zx) for x, zx in enumerate(z) if zx]
        for e, c in poly.items():
            for x, zx in nz:
                e2 = list(e)
                e2[x] += 1
                e2 = tuple(e2)
                new[e2] = new.get(e2, 0) + c * zx
        poly = {e: c for e, c in new.items() if c}
    return poly
