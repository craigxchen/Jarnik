"""Referee's own exact library (independent of L4_checks/).

Gaussian integers as Python int pairs. T_Q = norm-one subgroup of (Z[i]/Q)^*, Q = q^A, q odd prime.
Discrete logs by Pohlig-Hellman + baby-step giant-step (exact).
"""
import math, random

def sieve(n):
    s = bytearray([1]) * (n + 1)
    s[0:2] = b'\x00\x00'
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]

def is_prime(n):
    if n < 2:
        return False
    small = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    for p in small:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2; s += 1
    for a in small:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True

def factor(n):
    f = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f

def legendre(a, q):
    a %= q
    if a == 0:
        return 0
    return 1 if pow(a, (q - 1) // 2, q) == 1 else -1

def jacobi(a, n):
    assert n > 0 and n % 2 == 1
    a %= n
    res = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                res = -res
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            res = -res
        a %= n
    return res if n == 1 else 0

def two_squares(p):
    """x^2 + y^2 = p for prime p = 1 mod 4 (x odd > 0, y even > 0)."""
    for x in range(1, int(math.isqrt(p)) + 1):
        y2 = p - x * x
        y = math.isqrt(y2)
        if y * y == y2:
            if x % 2 == 0:
                x, y = y, x
            return (x, y)
    raise ValueError(p)

# Gaussian arithmetic (exact, unbounded)
def gmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])

def gconj(a):
    return (a[0], -a[1])

def gpow(a, e):
    r = (1, 0)
    b = a
    while e:
        if e & 1:
            r = gmul(r, b)
        b = gmul(b, b)
        e >>= 1
    return r

# arithmetic mod Q
def mmul(a, b, Q):
    return ((a[0] * b[0] - a[1] * b[1]) % Q, (a[0] * b[1] + a[1] * b[0]) % Q)

def mpow(a, e, Q):
    r = (1, 0)
    b = (a[0] % Q, a[1] % Q)
    while e:
        if e & 1:
            r = mmul(r, b, Q)
        b = mmul(b, b, Q)
        e >>= 1
    return r

def minv(a, Q):
    n = (a[0] * a[0] + a[1] * a[1]) % Q
    ni = pow(n, -1, Q)
    return ((a[0] * ni) % Q, (-a[1] * ni) % Q)

def ratio(a, Q):
    """a / conj(a) mod Q (a a unit mod q)."""
    return mmul(a, minv(gconj(a), Q), Q)

class Torus:
    """T_Q for Q = q^A, with generator g and discrete log."""
    def __init__(self, q, A, rng=None):
        self.q, self.A = q, A
        self.Q = q ** A
        chi = 1 if q % 4 == 1 else -1
        self.chi = chi
        self.m = (q - chi) * q ** (A - 1)
        self.fm = factor(self.m)
        rng = rng or random.Random(12345 + q)
        while True:
            w = (rng.randrange(self.Q), rng.randrange(self.Q))
            if (w[0] * w[0] + w[1] * w[1]) % q == 0:
                continue
            g = ratio(w, self.Q)
            if all(mpow(g, self.m // l, self.Q) != (1, 0) for l in self.fm):
                self.g = g
                break
        self.iota = self.log((0, 1))

    def log(self, h):
        """discrete log of h (in T_Q) to base g, Pohlig-Hellman."""
        Q, m = self.Q, self.m
        h = (h[0] % Q, h[1] % Q)
        residues, moduli = [], []
        for l, e in self.fm.items():
            le = l ** e
            g1 = mpow(self.g, m // le, Q)
            h1 = mpow(h, m // le, Q)
            # solve g1^x = h1, x mod l^e, digit by digit
            gl = mpow(g1, l ** (e - 1), Q)  # order l
            # baby-step table for order-l subgroup
            table = {}
            cur = (1, 0)
            for k in range(l):
                table[cur] = k
                cur = mmul(cur, gl, Q)
            x = 0
            for k in range(e):
                t = mmul(h1, mpow(minv(g1, Q), x, Q), Q)  # g1^{-x}: inverse of norm-one = conj
                t = mpow(t, l ** (e - 1 - k), Q)
                d = table[t]
                x += d * l ** k
            residues.append(x % le)
            moduli.append(le)
        # CRT
        x, M = 0, 1
        for r, mm in zip(residues, moduli):
            # solve x' = x mod M, x' = r mod mm
            t = ((r - x) * pow(M, -1, mm)) % mm
            x = x + M * t
            M *= mm
        assert mpow(self.g, x % m, Q) == h, "dlog failed"
        return x % m
