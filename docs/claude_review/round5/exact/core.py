"""Exact core for round5/exact.md: the normalised (leading-order) angle-model problem.

Conventions
-----------
* Rows x = 0..M-1.  A *pattern* is a nonempty subset S of {1..M-1}; it stands for the binary column
  1_S and also for its complement (identical constraints, since characters are zero-sum).  Patterns are
  bitmasks with bit x <-> row x (bit 0 always clear).
* A character is lambda in Q^M with sum 0, stored as mu = (lambda_1..lambda_{M-1}) (lambda_0 = -sum mu).
  lambda(S) = sum_{x in S} mu_x.  ||lambda||_1 = sum|mu_x| + |sum mu|.
* Support T (list of distinct patterns), weights om (Fractions >= 1 in the normalised problem).
  N(lambda) = sum_S om_S |lambda(S)|, Omega = sum om, sigma_lambda = N(lambda) - Omega/2.
* L_T = {lambda : lambda(S) in Z, S in T};  G_T = L_T / Z^M_0, finite iff the M-1 columns of the
  |T| x (M-1) matrix Q (Q[S][x] = [x in S]) are independent (rows affinely independent on T).
* Windings: k in Z^M with k_0 = 0 act by h(lambda) = lambda.k = mu.kappa mod 1, kappa = (k_1..k_{M-1}).
  A character is *relevant* for (T, k) iff lambda in L_T \\ {0} and lambda.k in Z.
* Requirement r(lambda) = rho0 + rho1 * ||lambda||_1 / 2.  (T, om, k) is *admissible at (rho0, rho1)*
  iff sigma_lambda >= r(lambda) for every relevant lambda ('strict': >).

Everything here is exact (Python ints / Fractions; numpy int64 only on integer data whose size is
checked).  Floating point appears only in `left_inverse_lp`, whose output is corrected to an exact
left inverse before use.
"""
import itertools, math
from fractions import Fraction
import numpy as np


# ------------------------------------------------------------------------------------------------
# patterns

def all_patterns(M):
    return list(range(2, 1 << M, 2))


def pattern_of(M, rows):
    """Canonical pattern of a subset of rows (complemented if it contains row 0)."""
    s = 0
    for x in rows:
        s |= 1 << x
    if s & 1:
        s ^= (1 << M) - 1
    return s


def rows_of(s, M):
    return [x for x in range(M) if (s >> x) & 1]


def Qmat(M, T):
    """|T| x (M-1) 0/1 int matrix, column j <-> row x = j+1."""
    Q = np.zeros((len(T), M - 1), dtype=np.int64)
    for i, s in enumerate(T):
        for x in range(1, M):
            if (s >> x) & 1:
                Q[i, x - 1] = 1
    return Q


# ------------------------------------------------------------------------------------------------
# exact linear algebra helpers

def frac_solve_left_inverse(Q):
    """Exact left inverse A (list of lists of Fractions, (M-1) x |T|) with A Q = I, built from a
    maximal independent set of rows of Q.  Returns None if rank < number of columns."""
    rows = [list(map(int, r)) for r in Q]
    m = len(rows)
    n = len(rows[0])
    # Gaussian elimination to select independent rows
    chosen = []
    basis = []  # reduced rows (Fractions)
    piv_cols = []
    for i, r in enumerate(rows):
        v = [Fraction(a) for a in r]
        for b, pc in zip(basis, piv_cols):
            if v[pc] != 0:
                f = v[pc] / b[pc]
                v = [a - f * c for a, c in zip(v, b)]
        nz = [j for j in range(n) if v[j] != 0]
        if nz:
            chosen.append(i)
            basis.append(v)
            piv_cols.append(nz[0])
            if len(chosen) == n:
                break
    if len(chosen) < n:
        return None
    # invert the square submatrix B = Q[chosen]
    B = [[Fraction(rows[i][j]) for j in range(n)] for i in chosen]
    Aug = [B[i] + [Fraction(int(i == j)) for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if Aug[r][c] != 0)
        Aug[c], Aug[p] = Aug[p], Aug[c]
        pv = Aug[c][c]
        Aug[c] = [a / pv for a in Aug[c]]
        for r in range(n):
            if r != c and Aug[r][c] != 0:
                f = Aug[r][c]
                Aug[r] = [a - f * b for a, b in zip(Aug[r], Aug[c])]
    Binv = [row[n:] for row in Aug]           # B^{-1}: n x n,  B^{-1} B = I
    A = [[Fraction(0)] * m for _ in range(n)]
    for jj, i in enumerate(chosen):
        for x in range(n):
            A[x][i] = Binv[x][jj]
    return A


def int_rank(Q):
    return np.linalg.matrix_rank(np.array(Q, dtype=float))


# ------------------------------------------------------------------------------------------------
# Smith-type diagonalisation with column transform:  U Q V = D

def diagonalise(Q):
    """Returns (d, V): V unimodular (n x n, Python ints) and d (length n, ints) with
    {mu in Q^n : Q mu in Z^m} = V . prod_i (1/d_i) Z  (d_i = 0 means the column space is deficient).
    Algorithm: integer row/column reduction (rows transforms discarded)."""
    A = [list(map(int, r)) for r in Q]
    m = len(A)
    n = len(A[0]) if m else 0
    V = [[int(i == j) for j in range(n)] for i in range(n)]

    def col_swap(i, j):
        for r in A:
            r[i], r[j] = r[j], r[i]
        for r in V:
            r[i], r[j] = r[j], r[i]

    def col_add(src, dst, f):          # col_dst += f * col_src
        if f:
            for r in A:
                r[dst] += f * r[src]
            for r in V:
                r[dst] += f * r[src]

    d = []
    t = 0
    while t < n:
        # pivot: smallest nonzero |entry| in the remaining block
        piv = None
        for i in range(t, m):
            for j in range(t, n):
                if A[i][j] and (piv is None or abs(A[i][j]) < abs(A[piv[0]][piv[1]])):
                    piv = (i, j)
        if piv is None:
            d.extend([0] * (n - t))
            break
        A[t], A[piv[0]] = A[piv[0]], A[t]
        col_swap(t, piv[1])
        while True:
            p = A[t][t]
            done = True
            for i in range(t + 1, m):           # row ops (allowed: they do not affect the lattice)
                if A[i][t]:
                    f = A[i][t] // p
                    A[i] = [a - f * b for a, b in zip(A[i], A[t])]
                    if A[i][t]:
                        done = False
            for j in range(t + 1, n):
                if A[t][j]:
                    col_add(t, j, -(A[t][j] // p))
                    if A[t][j]:
                        done = False
            if done:
                break
            # move smallest nonzero of row t / column t to pivot
            best = (abs(A[t][t]), 'r', t)
            for i in range(t + 1, m):
                if A[i][t] and abs(A[i][t]) < best[0]:
                    best = (abs(A[i][t]), 'r', i)
            for j in range(t + 1, n):
                if A[t][j] and abs(A[t][j]) < best[0]:
                    best = (abs(A[t][j]), 'c', j)
            if best[1] == 'r' and best[2] != t:
                A[t], A[best[2]] = A[best[2]], A[t]
            elif best[1] == 'c':
                col_swap(t, best[2])
        d.append(abs(A[t][t]))
        t += 1
    return d, V


def int_mat_inv(V):
    n = len(V)
    Aug = [[Fraction(V[i][j]) for j in range(n)] + [Fraction(int(i == j)) for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if Aug[r][c] != 0)
        Aug[c], Aug[p] = Aug[p], Aug[c]
        pv = Aug[c][c]
        Aug[c] = [a / pv for a in Aug[c]]
        for r in range(n):
            if r != c and Aug[r][c] != 0:
                f = Aug[r][c]
                Aug[r] = [a - f * b for a, b in zip(Aug[r], Aug[c])]
    inv = [row[n:] for row in Aug]
    for row in inv:
        for v in row:
            assert v.denominator == 1, 'V not unimodular'
    return [[int(v) for v in row] for row in inv]


class Group:
    """G_T = L_T / Z^{M-1} (mu-coordinates)."""

    def __init__(self, M, T):
        self.M = M
        self.T = list(T)
        Q = Qmat(M, T)
        self.Q = Q
        if len(T) < M - 1:
            self.full_rank = False
            self.d_all = None
            return
        d, V = diagonalise(Q)
        self.full_rank = all(x != 0 for x in d)
        self.d_all = d
        if not self.full_rank:
            return
        self.V = V
        self.keep = [i for i, x in enumerate(d) if x > 1]
        self.invs = [d[i] for i in self.keep]
        self.order = 1
        for x in self.invs:
            self.order *= x
        n = M - 1
        self.gens = [[Fraction(V[r][i], d[i]) for r in range(n)] for i in self.keep]
        Vinv = int_mat_inv(V)
        self.VinvT = [[Vinv[j][i] for j in range(n)] for i in range(n)]   # (V^{-1})^T

    def element(self, nvec):
        mu = [Fraction(0)] * (self.M - 1)
        for ni, g in zip(nvec, self.gens):
            if ni:
                for x in range(self.M - 1):
                    mu[x] += ni * g[x]
        return mu

    def elements(self):
        return itertools.product(*[range(x) for x in self.invs])

    def hom_value(self, eta, nvec):
        return sum(Fraction(e * ni, dd) for e, ni, dd in zip(eta, nvec, self.invs)) % 1

    def kappa_of(self, eta):
        """kappa in Z^{M-1} with mu.kappa = h_eta(mu) mod 1 for mu in L_T."""
        n = self.M - 1
        full = [0] * n
        for j, i in enumerate(self.keep):
            full[i] = eta[j]
        return [sum(self.VinvT[r][i] * full[i] for i in range(n)) for r in range(n)]

    def eta_of_kappa(self, kappa):
        """Inverse: eta_j = (V^T kappa)_i mod d_i."""
        n = self.M - 1
        out = []
        for j, i in enumerate(self.keep):
            out.append(sum(self.V[r][i] * kappa[r] for r in range(n)) % self.invs[j])
        return out

    def kernel(self, eta):
        return [nv for nv in self.elements() if any(nv) and self.hom_value(eta, nv) == 0]

    def minimal_kernels(self, max_homs=None):
        """Inclusion-minimal kernels of homomorphisms (as frozensets of nonzero elements), with an eta
        realising each."""
        elems = [nv for nv in self.elements() if any(nv)]
        kers = {}
        count = 0
        for eta in itertools.product(*[range(x) for x in self.invs]):
            count += 1
            if max_homs and count > max_homs:
                break
            K = frozenset(nv for nv in elems if self.hom_value(eta, nv) == 0)
            if K not in kers:
                kers[K] = eta
        ks = list(kers)
        minimal = [K for K in ks if not any((K2 < K) for K2 in ks)]
        return [(kers[K], K) for K in minimal]


def frac_part(mu):
    return [v - math.floor(v) for v in mu]


# ------------------------------------------------------------------------------------------------
# coercivity certificate

def left_inverse_lp(Q, om_float, target='l1'):
    """Float LP: left inverse A (n x m, A Q = I) minimising max_S c_S(A)/om_S, where
    c_S = sum_x |A_xS| + |sum_x A_xS|  (target 'l1').  Returns float array or None."""
    from scipy.optimize import linprog
    m, n = Q.shape
    # variables: A (n*m), P (n*m) >= |A|, s (m) >= |sum_x A_xS|, t
    nA = n * m
    nv = 2 * nA + m + 1
    c = np.zeros(nv); c[-1] = 1
    rows = []; rhs = []
    def idx(x, S):
        return x * m + S
    for x in range(n):
        for S in range(m):
            r = np.zeros(nv); r[idx(x, S)] = 1; r[nA + idx(x, S)] = -1; rows.append(r); rhs.append(0)
            r = np.zeros(nv); r[idx(x, S)] = -1; r[nA + idx(x, S)] = -1; rows.append(r); rhs.append(0)
    for S in range(m):
        r = np.zeros(nv); r[[idx(x, S) for x in range(n)]] = 1; r[2 * nA + S] = -1; rows.append(r); rhs.append(0)
        r = np.zeros(nv); r[[idx(x, S) for x in range(n)]] = -1; r[2 * nA + S] = -1; rows.append(r); rhs.append(0)
        r = np.zeros(nv); r[[nA + idx(x, S) for x in range(n)]] = 1; r[2 * nA + S] = 1; r[-1] = -om_float[S]
        rows.append(r); rhs.append(0)
    Aeq = np.zeros((n * n, nv)); beq = np.zeros(n * n)
    for x in range(n):
        for y in range(n):
            for S in range(m):
                Aeq[x * n + y, idx(x, S)] = Q[S, y]
            beq[x * n + y] = 1.0 if x == y else 0.0
    bounds = [(None, None)] * nA + [(0, None)] * (nA + m + 1)
    res = linprog(c, A_ub=np.array(rows), b_ub=np.array(rhs), A_eq=Aeq, b_eq=beq, bounds=bounds,
                  method='highs')
    if res.status != 0:
        return None
    return res.x[:nA].reshape(n, m)


def exact_left_inverse(Q, A_float=None, den=10 ** 6):
    """Exact A with A Q = I: rationalise A_float and correct with an exact left inverse A0:
    A = A_r + (I - A_r Q) A0."""
    A0 = frac_solve_left_inverse(Q)
    if A0 is None:
        return None
    if A_float is None:
        return A0
    m, n = Q.shape
    Ar = [[Fraction(float(A_float[x, S])).limit_denominator(den) for S in range(m)] for x in range(n)]
    # E = I - Ar Q   (n x n)
    E = [[Fraction(int(x == y)) - sum(Ar[x][S] * int(Q[S, y]) for S in range(m)) for y in range(n)]
         for x in range(n)]
    A = [[Ar[x][S] + sum(E[x][y] * A0[y][S] for y in range(n)) for S in range(m)] for x in range(n)]
    # exact check
    for x in range(n):
        for y in range(n):
            assert sum(A[x][S] * int(Q[S, y]) for S in range(m)) == int(x == y)
    return A


def coercivity(Q, om, A):
    """kappa_inf, kappa_1 (Fractions) with ||lambda||_inf <= kappa_inf N(lambda), ||lambda||_1 <= kappa_1 N(lambda)
    for every lambda in Q^M_0, from the exact left inverse A."""
    m, n = Q.shape
    k_inf = Fraction(0)
    k_1 = Fraction(0)
    for S in range(m):
        col = [A[x][S] for x in range(n)]
        s_abs = sum(abs(v) for v in col)
        s_sum = abs(sum(col))
        k_1 = max(k_1, (s_abs + s_sum) / om[S])
        k_inf = max(k_inf, max(max(abs(v) for v in col), s_sum) / om[S])
    return k_inf, k_1


# ------------------------------------------------------------------------------------------------
# enumeration of lattice points of a class in an l1/linf box (exact integers)

def enum_class(F, D, Binf_num, B1_num, chunk_limit=4_000_000):
    """All integer vectors v in Z^{n} with v = F (mod D) componentwise (F in [0, D)^n integers),
    |v_x| <= Binf_num, |sum v| <= Binf_num, sum|v_x| + |sum v| <= B1_num, v != 0.
    Returned as an int64 array (rows); these are D*mu of the lifts mu.  Yields blocks."""
    n = len(F)
    # candidate values per coordinate
    cands = []
    for x in range(n):
        lo = -Binf_num
        hi = Binf_num
        # smallest value >= lo congruent to F[x] mod D
        a0 = lo + ((F[x] - lo) % D)
        vals = list(range(a0, hi + 1, D))
        cands.append(np.array(vals, dtype=np.int64))
    # recursive vectorised build with l1 pruning
    partial = np.zeros((1, 0), dtype=np.int64)
    l1 = np.zeros(1, dtype=np.int64)
    sm = np.zeros(1, dtype=np.int64)
    for x in range(n):
        c = cands[x]
        if len(c) == 0:
            return
        rep = np.repeat(np.arange(partial.shape[0]), len(c))
        vv = np.tile(c, partial.shape[0])
        nl1 = l1[rep] + np.abs(vv)
        nsm = sm[rep] + vv
        # remaining coordinates contribute >= 0; also final |sum| >= |nsm| - (remaining max)
        ok = nl1 <= B1_num
        rep = rep[ok]; vv = vv[ok]; nl1 = nl1[ok]; nsm = nsm[ok]
        partial = np.hstack([partial[rep], vv[:, None]])
        l1 = nl1
        sm = nsm
        if partial.shape[0] > chunk_limit * 8:
            raise MemoryError('enumeration too large: %d partial vectors' % partial.shape[0])
    tot = l1 + np.abs(sm)
    ok = (tot <= B1_num) & (np.abs(sm) <= Binf_num) & np.any(partial != 0, axis=1)
    out = partial[ok]
    for i in range(0, out.shape[0], chunk_limit):
        yield out[i:i + chunk_limit]


def lcm_list(xs):
    out = 1
    for v in xs:
        out = out * v // math.gcd(out, v)
    return out
