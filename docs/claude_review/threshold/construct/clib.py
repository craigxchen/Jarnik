"""Independent exact tools for construct.md (written for this round; does not import other agents' code).

* hilb(pts, q): Hilbert function of a finite point set A in Z^n modulo the prime q,
  h(d) = rank of {polynomials of degree <= d} restricted to A.  The basis used is the
  binomial (Newton) basis prod_j binom(x_j - min_j, a_j), sum a <= d, which spans the same space
  as the monomials of degree <= d.  rank mod q <= rank over Q, so "h(d) = |A| mod q" proves
  reg(A) <= d over Q.  A lower bound reg(A) >= m is certified separately (explicit measure).
* reg2(pts): reg with two primes; raises if they disagree.
* cut widths, the half-separating LP (exact rational simplex, Bland's rule).
"""
import itertools
from fractions import Fraction
import numpy as np

Q1, Q2 = 2097143, 1999993          # primes < 2^21 so that int64 products of 3 residues are safe


def _isprime(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


assert _isprime(Q1) and _isprime(Q2)


def compositions(n, d):
    """all a in N^n with sum a == d"""
    if n == 0:
        if d == 0:
            yield ()
        return
    if n == 1:
        yield (d,)
        return
    for a in range(d + 1):
        for rest in compositions(n - 1, d - a):
            yield (a,) + rest


def _binom_table(X, q, amax):
    """B[j][a] = binom(X[:,j], a) mod q as int64 vectors (X >= 0)."""
    npts, n = X.shape
    inv = [0, 1] + [0] * amax
    for a in range(2, amax + 1):
        inv[a] = pow(a, q - 2, q)
    B = []
    for j in range(n):
        col = [np.ones(npts, dtype=np.int64)]
        for a in range(1, amax + 1):
            # binom(x, a) = binom(x, a-1) * (x - a + 1) / a
            v = (col[-1] * ((X[:, j] - a + 1) % q)) % q
            v = (v * inv[a]) % q
            col.append(v)
        B.append(col)
    return B


def _rref(V, q):
    V = V.copy()
    rows, cols = V.shape
    r = 0
    piv = []
    for c in range(cols):
        if r == rows:
            break
        nz = np.nonzero(V[r:, c])[0]
        if len(nz) == 0:
            continue
        i = r + nz[0]
        if i != r:
            V[[r, i]] = V[[i, r]]
        V[r] = (V[r] * pow(int(V[r, c]), q - 2, q)) % q
        colv = V[:, c].copy()
        colv[r] = 0
        nzr = np.nonzero(colv)[0]
        if len(nzr):
            V[nzr] = (V[nzr] - (colv[nzr, None] * V[r][None, :]) % q) % q
        piv.append(c)
        r += 1
    return V[:r], piv


def hilb(pts, q=Q1, dmax=None):
    """Hilbert function list h[0..reg] (stops when h = |A|) of the point set pts (tuples)."""
    A = np.array(pts, dtype=np.int64)
    if A.ndim == 1:
        A = A[:, None]
    A = A - A.min(axis=0)
    size, n = A.shape
    amax = int(A.max()) + 2 if size else 2
    amax = min(amax, size + 2)
    Bt = _binom_table(A, q, amax)
    basis = np.zeros((0, size), dtype=np.int64)
    piv = []
    h = []
    d = 0
    while True:
        rows = []
        for a in compositions(n, d):
            if max(a) > amax:
                continue
            v = np.ones(size, dtype=np.int64)
            for j, aj in enumerate(a):
                if aj:
                    v = (v * Bt[j][aj]) % q
            rows.append(v)
        if rows:
            V = np.array(rows, dtype=np.int64)
            if piv:
                V = (V - (V[:, piv] @ basis) % q) % q
            R, p2 = _rref(V, q)
            if p2:
                if piv:
                    basis = (basis - (basis[:, p2] @ R) % q) % q
                basis = np.vstack([basis, R])
                piv = piv + p2
        h.append(len(piv))
        if len(piv) == size:
            return h
        d += 1
        if dmax is not None and d > dmax:
            return h


def reg2(pts):
    """regularity index, computed modulo two primes (must agree)."""
    h1 = hilb(pts, Q1)
    h2 = hilb(pts, Q2)
    if h1 != h2:
        raise RuntimeError("primes disagree: %s %s" % (h1, h2))
    return len(h1) - 1


def proj(pts):
    """drop the last coordinate (affine bijection on a slice sum k = s)."""
    return [tuple(k[:-1]) for k in pts]


# ---------------------------------------------------------------- cuts and widths
def cuts(M):
    """nontrivial cuts up to complement: bitmasks P with element M-1 not in P, P != 0."""
    return [P for P in range(1, 1 << (M - 1))]


def popc(x):
    return bin(x).count("1")


def widths(pts, M):
    """w[P] = max_{k} k(P) - min_{k} k(P) over pts, for every cut P (bitmask, M-1 not in P)."""
    A = np.array(pts, dtype=np.int64)
    w = {}
    for P in cuts(M):
        idx = [j for j in range(M) if (P >> j) & 1]
        v = A[:, idx].sum(axis=1)
        w[P] = int(v.max() - v.min())
    return w


def separates(P, i, j):
    return ((P >> i) & 1) != ((P >> j) & 1)


# ---------------------------------------------------------------- exact LP (maximize)
def simplex_max(c, A_ub, b_ub, A_eq, b_eq):
    """maximize c.x s.t. A_ub x <= b_ub, A_eq x = b_eq, x >= 0, exact Fractions.
    Two-phase tableau simplex with Bland's rule.  Returns (value, x)."""
    F = Fraction
    m1, m2 = len(A_ub), len(A_eq)
    n = len(c)
    # rows: ub rows with slack, eq rows; make rhs >= 0
    rows = []
    rhs = []
    nslack = m1
    for i in range(m1):
        row = [F(v) for v in A_ub[i]] + [F(0)] * nslack
        row[n + i] = F(1)
        r = F(b_ub[i])
        if r < 0:
            row = [-v for v in row]
            r = -r
        rows.append(row)
        rhs.append(r)
    for i in range(m2):
        row = [F(v) for v in A_eq[i]] + [F(0)] * nslack
        r = F(b_eq[i])
        if r < 0:
            row = [-v for v in row]
            r = -r
        rows.append(row)
        rhs.append(r)
    m = len(rows)
    N = n + nslack
    # artificial variables for every row (simple, robust)
    T = [rows[i] + [F(1) if k == i else F(0) for k in range(m)] + [rhs[i]] for i in range(m)]
    basis = [N + i for i in range(m)]
    tot = N + m

    def pivot(r, col):
        pv = T[r][col]
        T[r] = [v / pv for v in T[r]]
        for i in range(m):
            if i != r and T[i][col] != 0:
                f = T[i][col]
                T[i] = [a - f * b for a, b in zip(T[i], T[r])]
        basis[r] = col

    def run(obj, allowed):
        # obj: list of length tot (maximize obj.x)
        while True:
            # reduced costs
            cb = [obj[b] for b in basis]
            best = None
            for j in range(tot):
                if j not in allowed or j in basis:
                    continue
                rc = obj[j] - sum(cb[i] * T[i][j] for i in range(m))
                if rc > 0:
                    best = j
                    break  # Bland
            if best is None:
                return
            ratios = []
            for i in range(m):
                if T[i][best] > 0:
                    ratios.append((T[i][-1] / T[i][best], basis[i], i))
            if not ratios:
                raise RuntimeError("unbounded")
            ratios.sort()
            pivot(ratios[0][2], best)

    # phase 1: maximize -sum artificials
    obj1 = [F(0)] * N + [F(-1)] * m
    run(obj1, set(range(tot)))
    if any(basis[i] >= N and T[i][-1] != 0 for i in range(m)):
        raise RuntimeError("infeasible")
    # drive artificials out of basis where possible
    for i in range(m):
        if basis[i] >= N:
            for j in range(N):
                if T[i][j] != 0:
                    pivot(i, j)
                    break
    obj2 = [F(v) for v in c] + [F(0)] * nslack + [F(0)] * m
    run(obj2, set(range(N)))
    x = [F(0)] * N
    for i in range(m):
        if basis[i] < N:
            x[basis[i]] = T[i][-1]
    val = sum(F(c[j]) * x[j] for j in range(n))
    return val, x[:n]


def adversary(w, M):
    """max over half-separating cut measures mu of E_mu w_P (exact). Returns (value, mu dict)."""
    P_list = cuts(M)
    c = [w[P] for P in P_list]
    A_ub, b_ub = [], []
    for i in range(M):
        for j in range(i + 1, M):
            A_ub.append([-1 if separates(P, i, j) else 0 for P in P_list])
            b_ub.append(Fraction(-1, 2))
    A_eq = [[1] * len(P_list)]
    b_eq = [1]
    val, x = simplex_max(c, A_ub, b_ub, A_eq, b_eq)
    mu = {P: x[t] for t, P in enumerate(P_list) if x[t] != 0}
    return val, mu


def sym_adversary(wsize, M):
    """For S_M-symmetric width profiles: wsize[j] = width of cuts of size j (1 <= j <= M//2).
    The optimal mu may be taken symmetric; pi_j = Pr[uniform size-j cut separates a fixed pair]."""
    best = None
    pis = {j: Fraction(2 * j * (M - j), M * (M - 1)) for j in range(1, M // 2 + 1)}
    half = Fraction(1, 2)
    js = list(pis)
    # optimum of a 2-constraint LP is a mixture of at most two sizes
    for a in js:
        if pis[a] >= half:
            v = Fraction(wsize[a])
            if best is None or v > best[0]:
                best = (v, {a: Fraction(1)})
        for b in js:
            if pis[a] < half <= pis[b]:
                # theta*pi_a + (1-theta)*pi_b = 1/2
                th = (pis[b] - half) / (pis[b] - pis[a])
                v = th * wsize[a] + (1 - th) * wsize[b]
                if best is None or v > best[0]:
                    best = (v, {a: th, b: 1 - th})
    return best
