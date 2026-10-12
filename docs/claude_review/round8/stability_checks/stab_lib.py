"""Shared helpers for the round-8 stability checks (pure Python integers, mpmath for logs).

Gaussian integers are tuples (a, b) = a + b i.  Everything that is claimed to be *exact* in
stability.md is done with Python integers; logarithms are evaluated with mpmath at 60 digits
only where an inequality between logarithms is checked (the printed slacks are far from 0).
"""
from math import comb, isqrt, atan2, pi
import mpmath

mpmath.mp.dps = 60


# ---------------------------------------------------------------- Gaussian arithmetic
def mul(z, w):
    return (z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0])


def sub(z, w):
    return (z[0] - w[0], z[1] - w[1])


def conj(z):
    return (z[0], -z[1])


def norm(z):
    return z[0] * z[0] + z[1] * z[1]


def _rdiv(a, n):
    return (2 * a + n) // (2 * n)


def gdivmod(z, w):
    n = norm(w)
    t = mul(z, conj(w))
    q = (_rdiv(t[0], n), _rdiv(t[1], n))
    return q, sub(z, mul(q, w))


def ggcd(z, w):
    while w != (0, 0):
        _, r = gdivmod(z, w)
        z, w = w, r
    return z


def gpow(z, k):
    out = (1, 0)
    for _ in range(k):
        out = mul(out, z)
    return out


def exact_div(z, w):
    n = norm(w)
    t = mul(z, conj(w))
    assert t[0] % n == 0 and t[1] % n == 0
    return (t[0] // n, t[1] // n)


# ---------------------------------------------------------------- integers
def factor(n):
    out = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, isqrt(n) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


_GP = {}


def gauss_prime(p):
    """Gaussian prime pi with N(pi) = p, for p = 1 mod 4 (cached)."""
    if p in _GP:
        return _GP[p]
    for a in range(2, p):
        t = pow(a, (p - 1) // 4, p)
        if t * t % p == p - 1:
            break
    g = ggcd((p, 0), (t, 1))
    assert norm(g) == p
    _GP[p] = g
    return g


def points_of_norm(n, fac=None):
    """All Gaussian integers of norm n, sorted by argument in [0, 2 pi)."""
    if fac is None:
        fac = factor(n)
    base = [(1, 0)]
    for q, e in fac.items():
        if q == 2:
            base = [mul(z, gpow((1, 1), e)) for z in base]
        elif q % 4 == 3:
            if e % 2:
                return []
            base = [mul(z, (q ** (e // 2), 0)) for z in base]
        else:
            pq = gauss_prime(q)
            pb = conj(pq)
            new = []
            for z in base:
                for a in range(e + 1):
                    new.append(mul(z, mul(gpow(pq, a), gpow(pb, e - a))))
            base = new
    pts = []
    for z in base:
        w = z
        for _ in range(4):
            pts.append(w)
            w = (-w[1], w[0])
    pts = sorted(set(pts), key=lambda z: atan2(z[1], z[0]) % (2 * pi))
    return pts


# ---------------------------------------------------------------- F*(M) of two_thirds Thm 5.2
def E(n, q):
    if n <= 0:
        return 0
    b, r = divmod(n, q)
    return q * comb(b, 2) + r * b


def depth_sum(n, q0, p):
    tot = 0
    q = q0
    while q < n:
        tot += E(n, q)
        q *= p
    return tot


def psi8(M, p):
    w = [8 * depth_sum(n, p - 1, p) for n in range(M + 1)]
    f = [None] * (M + 1)
    for S in range(1, M):
        f[S] = w[S] + (2 * S - M) ** 2
    for S in range(1, M):
        for S2 in range(S + 1, M):
            v = f[S] + w[S2 - S] + (2 * S2 - M) ** 2
            if v < f[S2]:
                f[S2] = v
    return min(f[S] + w[M - S] for S in range(1, M))


_FS = {}


def Fstar(M):
    """F*(M) as an mpmath number (sum of rational multiples of logs of primes)."""
    if M in _FS:
        return _FS[M]
    if M < 2:
        _FS[M] = mpmath.mpf(0)
        return _FS[M]
    tot = mpmath.mpf(comb(M, 2)) * mpmath.log(2) / 2
    for q in primes_upto(M + 2):
        if q == 2:
            continue
        if q % 4 == 3:
            if q + 1 >= M:
                continue
            tot += depth_sum(M, q + 1, q) * mpmath.log(q)
        else:
            if q - 1 >= M:
                continue
            free8 = 8 * depth_sum(M, q - 1, q)
            g8 = min(free8, psi8(M, q))
            tot += mpmath.mpf(g8) / 8 * mpmath.log(q)
    _FS[M] = tot
    return tot


def bfun(s, W, logC):
    """b(s) = W/s + (4(s-1)/s) log C - 8 F*(s)/s^2 : the sub-cluster gcd bound."""
    return W / s + 4 * mpmath.mpf(s - 1) / s * logC - 8 * Fstar(s) / (s * s)


K0 = mpmath.mpf('13.9092')  # two_thirds Theorem 7.5 constant (upper bound)
