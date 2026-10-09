# Exact / modular tools for the tau(M) problem.
# Point sets live in Z^M (k_1..k_M); polynomials are taken in the coordinates k_2..k_M
# (k_1 = s - sum is determined on a slice sum k = s).
import itertools, math
from fractions import Fraction
import numpy as np

P1 = 2147483647          # 2^31 - 1
P2 = 2147483629          # 2^31 - 19
P3 = 2147483587          # 2^31 - 61

def is_prime(n):
    if n < 2: return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0: return n == q
    d, r = n - 1, 0
    while d % 2 == 0: d //= 2; r += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1: break
        else:
            return False
    return True

assert is_prime(P1) and is_prime(P2) and is_prime(P3)

def slice_points(M, s, D, shell_only):
    """k in Z^M with sum k = s and |k|_1 = D (shell) or |k|_1 <= D, |k|_1 = D mod 2 (ball)."""
    out = []
    def rec(prefix, rem):
        j = len(prefix)
        if j == M - 1:
            last = s - sum(prefix)
            l1 = D - rem + abs(last)
            if abs(last) <= rem and (l1 - D) % 2 == 0 and (not shell_only or l1 == D):
                out.append(tuple(prefix) + (last,))
            return
        for a in range(-rem, rem + 1):
            rec(prefix + [a], rem - abs(a))
    rec([], D)
    return out

def shell(M, s, t): return slice_points(M, s, t, True)
def ball(M, s, D): return slice_points(M, s, D, False)

def monomials(n, d):
    """exponent tuples in n variables of total degree exactly d"""
    out = []
    for c in itertools.combinations_with_replacement(range(n), d):
        e = [0] * n
        for i in c: e[i] += 1
        out.append(tuple(e))
    return out

def eval_matrix_mod(pts, dmax, p):
    """columns = monomials in k_2..k_M ordered by degree; returns (T, degs) with
    T[i, j] = monomial_j(pts_i) mod p and degs[j] = degree of monomial j."""
    X = np.array([k[1:] for k in pts], dtype=np.int64) % p
    npts, n = X.shape
    cols, degs = [], []
    frontier = {(): np.ones(npts, dtype=np.int64)}
    cols.append(frontier[()]); degs.append(0)
    for d in range(1, dmax + 1):
        newf = {}
        for key, v in frontier.items():
            start = key[-1] if key else 0
            for j in range(start, n):
                newf[key + (j,)] = (v * X[:, j]) % p
        frontier = newf
        for v in newf.values():
            cols.append(v); degs.append(d)
    return np.array(cols, dtype=np.int64).T.copy(), degs

def pivot_columns_mod(T, p, stop_at=None):
    """Gaussian elimination left to right; returns list of pivot column indices."""
    A = T.copy() % p
    nr, nc = A.shape
    r = 0
    piv = []
    for c in range(nc):
        if r == nr: break
        nz = np.nonzero(A[r:, c])[0]
        if len(nz) == 0: continue
        i = r + nz[0]
        if i != r: A[[r, i]] = A[[i, r]]
        inv = pow(int(A[r, c]), p - 2, p)
        A[r, c:] = (A[r, c:] * inv) % p
        col = A[:, c].copy(); col[r] = 0
        nzr = np.nonzero(col)[0]
        if len(nzr):
            A[np.ix_(nzr, np.arange(c, nc))] = (A[np.ix_(nzr, np.arange(c, nc))]
                                                - (col[nzr, None] * A[r, c:][None, :]) % p) % p
        piv.append(c); r += 1
    return piv

def hilbert_mod(pts, p, dmax=None):
    """h_A(d) = rank of polys of degree <= d on A, for d = 0..until full (mod p).
    Returns list h[0..reg] with h[reg] = |A| (if reached within dmax)."""
    n = len(pts)
    if dmax is None: dmax = n
    d = 0
    # grow dmax adaptively
    dm = 4
    while True:
        T, degs = eval_matrix_mod(pts, dm, p)
        piv = pivot_columns_mod(T, p)
        h = []
        for d in range(dm + 1):
            h.append(sum(1 for c in piv if degs[c] <= d))
        if h[-1] == n or dm >= dmax:
            # truncate at first full
            out = []
            for v in h:
                out.append(v)
                if v == n: break
            return out
        dm = min(2 * dm, dmax)

def reg_mod(pts, p):
    """reg_p(A) = min d with h_A(d) = |A| mod p.  reg_Q(A) <= reg_p(A)."""
    h = hilbert_mod(pts, p)
    assert h[-1] == len(pts)
    return len(h) - 1, h

# ---------- exact certification of a lower bound reg_Q(A) >= m ----------
def kernel_vector_exact(pts, d, primes=(P1, P2, P3)):
    """Try to find a nonzero integer vector c on pts with sum_k c_k k^alpha = 0 for all |alpha| <= d
    (alpha over monomials in k_2..k_M).  Uses RREF mod several primes, CRT, rational
    reconstruction, then verifies EXACTLY with Python integers. Returns c or None."""
    n = len(pts)
    mons = [()]
    allm = []
    for dd in range(d + 1): allm += monomials(len(pts[0]) - 1, dd)
    # exact integer matrix rows = monomials, cols = points
    def mon_val(k, e):
        v = 1
        for ki, ei in zip(k[1:], e):
            if ei: v *= ki ** ei
        return v
    Vex = [[mon_val(k, e) for k in pts] for e in allm]
    images = []
    for p in primes:
        A = np.array([[x % p for x in row] for row in Vex], dtype=np.int64)
        piv = pivot_columns_mod(A, p)
        # full RREF to read the kernel
        R = rref_mod(A, p)
        images.append((p, piv, R))
    pivsets = [tuple(x[1]) for x in images]
    if len(set(pivsets)) != 1: return None
    piv = list(pivsets[0])
    free = [j for j in range(n) if j not in piv]
    if not free: return None
    f = free[0]
    # kernel vector with x_f = 1, other free = 0, x_piv[i] = -R[i, f]
    vecs = []
    for p, _, R in images:
        x = [0] * n; x[f] = 1
        for i, c in enumerate(piv): x[c] = (-int(R[i, f])) % p
        vecs.append(x)
    # CRT
    Mod = 1
    for p in primes: Mod *= p
    comb = []
    for j in range(n):
        r = 0
        for (p, _, _), x in zip(images, vecs):
            Mp = Mod // p
            r = (r + x[j] * Mp * pow(Mp, -1, p)) % Mod
        comb.append(r)
    fr = []
    for r in comb:
        q = ratrecon(r, Mod)
        if q is None: return None
        fr.append(q)
    L = 1
    for q in fr: L = L * q.denominator // math.gcd(L, q.denominator)
    c = [int(q * L) for q in fr]
    g = 0
    for v in c: g = math.gcd(g, v)
    c = [v // g for v in c]
    # exact verification
    for row in Vex:
        if sum(a * b for a, b in zip(row, c)) != 0: return None
    if all(v == 0 for v in c): return None
    return c

def rref_mod(A, p):
    A = A.copy() % p
    nr, nc = A.shape
    r = 0
    for c in range(nc):
        if r == nr: break
        nz = np.nonzero(A[r:, c])[0]
        if len(nz) == 0: continue
        i = r + nz[0]
        if i != r: A[[r, i]] = A[[i, r]]
        inv = pow(int(A[r, c]), p - 2, p)
        A[r] = (A[r] * inv) % p
        col = A[:, c].copy(); col[r] = 0
        nzr = np.nonzero(col)[0]
        if len(nzr):
            A[nzr] = (A[nzr] - (col[nzr, None] * A[r][None, :]) % p) % p
        r += 1
    return A[:r]

def ratrecon(a, m):
    """rational reconstruction of a mod m with |num|, den <= sqrt(m/2)"""
    bound = math.isqrt(m // 2)
    r0, r1 = m, a % m
    s0, s1 = 0, 1
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound: return None
    return Fraction(r1, s1) if s1 > 0 else Fraction(-r1, -s1)

def order_exact(pts, c, maxdeg):
    """exact order of the measure c on pts: max m such that c annihilates all monomials of degree < m
    (monomials in k_2..k_M); searches up to maxdeg."""
    n = len(pts[0]) - 1
    for d in range(maxdeg + 1):
        for e in monomials(n, d):
            tot = 0
            for k, ck in zip(pts, c):
                v = ck
                for ki, ei in zip(k[1:], e):
                    if ei: v *= ki ** ei
                tot += v
            if tot != 0: return d
    return None  # >= maxdeg + 1
