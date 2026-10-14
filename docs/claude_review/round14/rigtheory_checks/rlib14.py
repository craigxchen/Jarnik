"""Round-14 library (written for this round; logic re-derived from round10 L4K.md §1, no imports from the repo).
Torus T_Q = norm-one elements of (Z[i]/Q)^*, cyclic of order m_Q = (q - chi(q)) q^(a-1).
Genuine labels L_j = log_g(pi_j / conj pi_j), normalised so that iota = log_g(i) = m/4.
"""
import math, random, itertools
from functools import lru_cache

def is_prime(n):
    if n < 2: return False
    if n % 2 == 0: return n == 2
    f = 3
    while f * f <= n:
        if n % f == 0: return False
        f += 2
    return True

def primes_upto(n): return [p for p in range(2, n + 1) if is_prime(p)]
def split_primes_upto(n): return [p for p in primes_upto(n) if p % 4 == 1]
def chi4(q): return 0 if q % 2 == 0 else (1 if q % 4 == 1 else -1)

def two_squares(p):
    for a in range(1, math.isqrt(p) + 1):
        b2 = p - a * a; b = math.isqrt(b2)
        if b * b == b2: return (a, b)
    raise ValueError(p)

def gmul(z, w, mod=None):
    a, b = z; c, d = w
    re, im = a * c - b * d, a * d + b * c
    if mod is not None: re %= mod; im %= mod
    return (re, im)
def gconj(z): return (z[0], -z[1])
def gpow(z, e, mod):
    r = (1, 0); b = (z[0] % mod, z[1] % mod)
    while e:
        if e & 1: r = gmul(r, b, mod)
        b = gmul(b, b, mod); e >>= 1
    return r
def ginv(z, mod):
    a, b = z; nz = (a * a + b * b) % mod; inv = pow(nz, -1, mod)
    return ((a * inv) % mod, (-b * inv) % mod)
def legendre(a, q):
    a %= q
    if a == 0: return 0
    return 1 if pow(a, (q - 1) // 2, q) == 1 else -1
def factor(n):
    f = {}; d = 2
    while d * d <= n:
        while n % d == 0: f[d] = f.get(d, 0) + 1; n //= d
        d += 1
    if n > 1: f[n] = f.get(n, 0) + 1
    return f

class Torus:
    def __init__(self, q, a=1, seed=7):
        self.q, self.a, self.mod = q, a, q ** a
        self.m = (q - chi4(q)) * q ** (a - 1)
        rng = random.Random(seed * 1000003 + 31 * q + a)
        fac = factor(self.m)
        while True:
            z = (rng.randrange(self.mod), rng.randrange(self.mod))
            if (z[0] * z[0] + z[1] * z[1]) % q == 0: continue
            g = gmul(z, ginv(gconj(z), self.mod), self.mod)
            if all(gpow(g, self.m // p, self.mod) != (1, 0) for p in fac):
                break
        t = {}; x = (1, 0)
        for k in range(self.m):
            t[x] = k; x = gmul(x, g, self.mod)
        assert len(t) == self.m
        self.g, self.table = g, t
        io = t[(0, 1)]
        assert io in (self.m // 4, 3 * self.m // 4)
        self.flip = (io != self.m // 4)   # if so, use generator g^{-1}: logs negate
    def dlog(self, z):
        k = self.table[(z[0] % self.mod, z[1] % self.mod)]
        return (-k) % self.m if self.flip else k
    def ratio_log(self, z):
        return self.dlog(gmul(z, ginv(gconj(z), self.mod), self.mod))

def gaussian_primes(ps): return [two_squares(p) for p in ps]

def genuine_labels(ps, q, a=1):
    T = Torus(q, a); pis = gaussian_primes(ps)
    L = [T.ratio_log(pi) for pi in pis]
    assert T.dlog((0, 1)) == T.m // 4
    return T.m, L

def nu(ps, q): return [1 if legendre(p, q) == -1 else 0 for p in ps]

def Zv(pis, v):
    z = (1, 0)
    for pi, e in zip(pis, v):
        b = pi if e > 0 else gconj(pi)
        for _ in range(abs(e)): z = gmul(z, b)
    return z

def Ds(z, s):
    X, Y = z
    return [Y, X - Y, X, X + Y][s % 4]

def VK(r, K, lo=1):
    out = []
    def rec(i, rem, cur):
        if i == r:
            if lo <= sum(abs(x) for x in cur): out.append(tuple(cur))
            return
        for x in range(-rem, rem + 1):
            cur.append(x); rec(i + 1, rem - abs(x), cur); cur.pop()
    rec(0, K, [])
    return out

def canon(v):
    for x in v:
        if x != 0: return x > 0
    return False

def VK_mod_sign(r, K, lo=1): return [v for v in VK(r, K, lo) if canon(v)]

def vq(x, q):
    x = abs(x); e = 0
    if x == 0: return 10 ** 9
    while x % q == 0: x //= q; e += 1
    return e
