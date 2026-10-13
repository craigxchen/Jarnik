"""Referee's independent library (does not import the author's lib.py).

Torus T_{q^a} = norm-1 elements of (Z[i]/q^a)^*, cyclic of order m_q q^(a-1), m_q = q - chi(q).
Labels L_j = log_g(pi_j/conj pi_j); iota = log_g(i).
"""
import math, itertools, random
from functools import lru_cache

def is_prime(n):
    if n < 2: return False
    if n % 2 == 0: return n == 2
    f = 3
    while f * f <= n:
        if n % f == 0: return False
        f += 2
    return True

def primes_upto(n):
    return [p for p in range(2, n + 1) if is_prime(p)]

def chi4(q):
    return 1 if q % 4 == 1 else -1

def two_squares(p):
    for a in range(1, int(math.isqrt(p)) + 1):
        b2 = p - a * a
        b = math.isqrt(b2)
        if b * b == b2:
            return (a, b)
    raise ValueError(p)

def gmul(z, w, mod=None):
    a, b = z; c, d = w
    re, im = a * c - b * d, a * d + b * c
    if mod is not None:
        re %= mod; im %= mod
    return (re, im)

def gconj(z):
    return (z[0], -z[1])

def gpow(z, e, mod):
    r = (1, 0); b = (z[0] % mod, z[1] % mod)
    while e:
        if e & 1: r = gmul(r, b, mod)
        b = gmul(b, b, mod); e >>= 1
    return r

def ginv(z, mod):
    # inverse of z mod mod: z^{-1} = conj z / N(z)
    a, b = z
    nz = (a * a + b * b) % mod
    inv = pow(nz, -1, mod)
    return ((a * inv) % mod, (-b * inv) % mod)

def legendre(a, q):
    a %= q
    if a == 0: return 0
    return 1 if pow(a, (q - 1) // 2, q) == 1 else -1

def factor(n):
    f = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1; n //= d
        d += 1
    if n > 1: f[n] = f.get(n, 0) + 1
    return f

class Torus:
    """Norm-1 torus mod q^a with a fixed generator and a dlog table (for small orders)."""
    def __init__(self, q, a=1, seed=1):
        self.q = q; self.a = a; self.mod = q ** a
        self.m1 = q - chi4(q)
        self.order = self.m1 * q ** (a - 1)
        rng = random.Random(seed * 1000003 + q * 7 + a)
        fac = factor(self.order)
        while True:
            z = (rng.randrange(self.mod), rng.randrange(self.mod))
            if (z[0] * z[0] + z[1] * z[1]) % q == 0: continue
            g = gmul(z, ginv(gconj(z), self.mod), self.mod)  # z / conj z
            ok = all(gpow(g, self.order // p, self.mod) != (1, 0) for p in fac)
            if ok and gpow(g, self.order, self.mod) == (1, 0):
                self.g = g; break
        self.table = None
        if self.order <= 2_000_000:
            t = {}; x = (1, 0)
            for k in range(self.order):
                t[x] = k; x = gmul(x, self.g, self.mod)
            assert len(t) == self.order
            self.table = t
    def dlog(self, z):
        z = (z[0] % self.mod, z[1] % self.mod)
        return self.table[z]
    def ratio_log(self, z):
        """log_g (z/conj z) for z a unit mod q."""
        w = gmul(z, ginv(gconj(z), self.mod), self.mod)
        return self.dlog(w)
    def iota(self):
        return self.dlog((0, 1))

def gaussian_primes(ps):
    return [two_squares(p) for p in ps]

def genuine_labels(ps, q, a=1):
    T = Torus(q, a)
    pis = gaussian_primes(ps)
    return T, [T.ratio_log(pi) for pi in pis]

def nu(ps, q):
    return [1 if legendre(p, q) == -1 else 0 for p in ps]

def Zv(pis, v):
    z = (1, 0)
    for pi, e in zip(pis, v):
        if e > 0:
            for _ in range(e): z = gmul(z, pi)
        elif e < 0:
            for _ in range(-e): z = gmul(z, gconj(pi))
    return z

def Ds(z, s):
    X, Y = z
    return [Y, X - Y, X, X + Y][s]

def VK(r, K):
    """all v in Z^r with 1<=||v||_1<=K"""
    out = []
    def rec(i, rem, cur):
        if i == r:
            if any(cur): out.append(tuple(cur))
            return
        for x in range(-rem, rem + 1):
            cur.append(x); rec(i + 1, rem - abs(x), cur); cur.pop()
    rec(0, K, [])
    return out

def VK_mod_sign(r, K):
    out = []
    for v in VK(r, K):
        # canonical: first nonzero positive
        for x in v:
            if x != 0:
                if x > 0: out.append(v)
                break
    return out

def w_of(ps, v):
    return sum(abs(x) * math.log(p) for x, p in zip(v, ps))

def budget(ps, v, s):
    return w_of(ps, v) / 2 + (math.log(2) / 2 if s % 2 else 0.0)

def split_primes_upto(n):
    return [p for p in primes_upto(n) if p % 4 == 1]
