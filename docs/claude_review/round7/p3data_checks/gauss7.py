"""Exact helpers for round 7 (p3data): primes, Gaussian integers, norm-one groups T_Q, Paley profiles.

Everything here is integer arithmetic.  Conventions follow round6/sharp.md:
  * Q = q^a odd prime power, T_Q = norm-one subgroup of (Z[i]/Q)^*, cyclic of order
    m_Q = (q - chi(q)) q^(a-1), chi(q) = +1 for q = 1 mod 4, -1 for q = 3 mod 4.
  * For a Gaussian prime pi of norm p (p prime to Q), t = pi/conj(pi) = pi^2/p in T_Q.
"""
import random
import numpy as np


# ---------------------------------------------------------------- primes
def is_prime(n):
    """Deterministic Miller-Rabin for n < 3.3e24."""
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41)
    for p in small:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2; s += 1
    for a in small:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0:2] = b'\x00\x00'
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def odd_prime_powers_upto(X):
    """List of (q, a, Q=q^a) for odd prime powers Q <= X, sorted by Q."""
    out = []
    for q in primes_upto(X):
        if q == 2:
            continue
        Q, a = q, 1
        while Q <= X:
            out.append((q, a, Q))
            Q *= q; a += 1
    return sorted(out, key=lambda t: t[2])


def mQ(q, a):
    return (q - 1 if q % 4 == 1 else q + 1) * q ** (a - 1)


def two_squares(p):
    """p = 1 mod 4 prime -> (a, b) with a^2 + b^2 = p, a odd > 0, b even > 0."""
    assert p % 4 == 1
    c = 2
    while pow(c, (p - 1) // 2, p) != p - 1:
        c += 1
    x = pow(c, (p - 1) // 4, p)
    a, b = p, x
    lim = int(p ** 0.5)
    while b > lim:
        a, b = b, a % b
    u = b
    v2 = p - u * u
    v = int(round(v2 ** 0.5))
    while v * v > v2:
        v -= 1
    while (v + 1) * (v + 1) <= v2:
        v += 1
    assert u * u + v * v == p
    if u % 2 == 0:
        u, v = v, u
    return abs(u), abs(v)


# ---------------------------------------------------------------- Gaussian arithmetic mod Q
def gmul(z, w, Q):
    return ((z[0] * w[0] - z[1] * w[1]) % Q, (z[0] * w[1] + z[1] * w[0]) % Q)


def gpow(z, k, Q):
    r = (1 % Q, 0)
    while k:
        if k & 1:
            r = gmul(r, z, Q)
        z = gmul(z, z, Q)
        k >>= 1
    return r


def gconj(z, Q):
    return (z[0] % Q, (-z[1]) % Q)


def factor_small(n):
    f, d = [], 2
    while d * d <= n:
        if n % d == 0:
            f.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        f.append(n)
    return f


class TGroup:
    """Norm-one group T_Q of Z[i]/Q with a generator g and a full discrete-log table.
    dlog[re*Q + im] = k with g^k = re + i*im, or -1 if not in T_Q."""

    def __init__(self, q, a, seed=0):
        self.q, self.a = q, a
        Q = q ** a
        self.Q = Q
        m = mQ(q, a)
        self.m = m
        rng = random.Random(1000003 * q + a + seed)
        pf = factor_small(m)
        while True:
            z = (rng.randrange(Q), rng.randrange(Q))
            nz = (z[0] * z[0] + z[1] * z[1]) % Q
            if nz % q == 0:
                continue
            # t = z / conj(z) = z^2 / N(z)
            ninv = pow(nz, -1, Q)
            z2 = gmul(z, z, Q)
            t = (z2[0] * ninv % Q, z2[1] * ninv % Q)
            assert gpow(t, m, Q) == (1 % Q, 0)
            if all(gpow(t, m // l, Q) != (1 % Q, 0) for l in pf):
                g = t
                break
        self.g = g
        tab = -np.ones(Q * Q, dtype=np.int64)
        x = (1 % Q, 0)
        for k in range(m):
            idx = x[0] * Q + x[1]
            assert tab[idx] == -1
            tab[idx] = k
            x = gmul(x, g, Q)
        assert x == (1 % Q, 0)
        self.tab = tab
        self.iota = int(tab[0 * Q + 1])      # i = g^iota
        assert self.iota in (m // 4, 3 * m // 4)

    def dlog_t(self, A, B, P):
        """Vectorised: for Gaussian primes pi = A + iB of norm P (int64 arrays), return
        ell = dlog(pi / conj pi) = dlog(pi^2 / P) mod m."""
        Q = self.Q
        A = np.asarray(A, dtype=object) % Q
        B = np.asarray(B, dtype=object) % Q
        P = np.asarray(P, dtype=object) % Q
        Pinv = np.array([pow(int(p), -1, Q) for p in P], dtype=object)
        re = ((A * A - B * B) * Pinv) % Q
        im = ((2 * A * B) * Pinv) % Q
        idx = (re * Q + im).astype(np.int64)
        ell = self.tab[idx]
        assert (ell >= 0).all()
        return ell


# ---------------------------------------------------------------- Paley profiles
def legendre(a, q):
    a %= q
    if a == 0:
        return 0
    return 1 if pow(a, (q - 1) // 2, q) == 1 else -1


def paley_normalised(q):
    """Normalised Paley-I Hadamard matrix of order q+1 (q = 3 mod 4 prime): first column all ones.
    Row 0 = infinity, rows 1..q = F_q; columns 0 (constant), 1..q = labels F_q.
    This is round5/checks/lib.py normalize(paley1(q))."""
    assert q % 4 == 3
    n = q + 1
    S = np.zeros((n, n), dtype=np.int64)
    S[0, 1:] = 1
    S[1:, 0] = -1
    for i in range(q):
        for j in range(q):
            S[1 + i, 1 + j] = legendre(j - i, q)
    H = np.eye(n, dtype=np.int64) + S
    H = H * H[:, [0]]
    assert (H.T @ H == n * np.eye(n, dtype=np.int64)).all()
    assert (H[:, 0] == 1).all()
    return H


def profile_construct(q, b):
    """round5/construct.md section 2: b unflipped copies of every label 1..M-1 plus one flipped copy
    per row x of label a(x) (a(0) = a* = 1, a(x) = x for x >= 1).  r = b(M-1) + M.  Returns S (+-1)."""
    H = paley_normalised(q)
    M = q + 1
    assign = [1] + list(range(1, M))
    cols = []
    for a in range(1, M):
        for _ in range(b):
            cols.append(H[:, a].copy())
    for x in range(M):
        col = H[:, assign[x]].copy()
        col[x] = -col[x]
        cols.append(col)
    return np.array(cols, dtype=np.int64).T


def profile_lower(q, b):
    """round5/lower.md section 3.1 / sharp.md section 4, canonical assignment: b physical columns per
    label; row x >= 1 flipped in copy 0 of label x; row 0 (infinity) flipped in copy 2 of label a* = 1.
    r = b(M-1).  Returns S (+-1)."""
    assert b >= 4
    H = paley_normalised(q)
    M = q + 1
    cols = []
    for a in range(1, M):
        for c in range(b):
            col = H[:, a].copy()
            if c == 0:
                col[a] = -col[a]          # row x = a flipped in physical column (a, 0)
            if a == 1 and c == 2:
                col[0] = -col[0]          # infinity row flipped in (a*, 2)
            cols.append(col)
    return np.array(cols, dtype=np.int64).T


def check_profile(S):
    """rank M, rows distinct, every column nonconstant."""
    M = S.shape[0]
    assert np.linalg.matrix_rank(S.astype(float)) == M
    assert len({tuple(r) for r in S}) == M
    assert all(len(set(S[:, j])) == 2 for j in range(S.shape[1]))


# ---------------------------------------------------------------- prime choices
def primes_1mod4_window(start, count):
    out = []
    n = start + (1 - start) % 4
    while len(out) < count:
        if is_prime(n):
            out.append(n)
        n += 4
    return out


def gaussian_primes(ps, rng, randomize_conj=True):
    """For each p: (A, B) with A^2 + B^2 = p, A odd > 0, B even; B negated (conjugate) at random."""
    A, B = [], []
    for p in ps:
        a, b = two_squares(p)
        if randomize_conj and rng.random() < 0.5:
            b = -b
        A.append(a); B.append(b)
    return A, B
