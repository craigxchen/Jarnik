"""Exact Gaussian-integer helpers for the round-9 L1/L2 checks.

Gaussian integers are tuples (a, b) = a + b i, Python ints throughout.  Floating point is used
only for angles (ordering / arc constants) and for printing logarithms.
"""
from math import isqrt, atan2, pi, log, gcd


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


def exact_div(z, w):
    n = norm(w)
    t = mul(z, conj(w))
    assert t[0] % n == 0 and t[1] % n == 0, (z, w)
    return (t[0] // n, t[1] // n)


def gpow(z, k):
    out = (1, 0)
    for _ in range(k):
        out = mul(out, z)
    return out


UNITS = [(1, 0), (0, 1), (-1, 0), (0, -1)]


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


def points_of_norm(fac):
    """All Gaussian integers of norm prod p^e (fac: dict p->e, all p = 1 mod 4)."""
    base = [(1, 0)]
    for q, e in fac.items():
        assert q % 4 == 1
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


def reduced_chord_data(zx, zy):
    """Return (g, beta, s, B) with zx - zy = g*beta, beta = i^{kx} gamma_s B (B rational int).

    g = gcd(zx, zy).  alpha = zx/g, alpha' = zy/g.  We find the unit eps with alpha' = eps conj(alpha)
    (coprime, equal norm => exists for split-only norms) and return s with eps = i^s * (unit of alpha)
    normalisation as in the notes: writing alpha = i^{kx} A, alpha' = i^{ky} conj(A) we have
    beta = i^{kx}(A - i^{s} conj A), s = ky - kx.  We take A = alpha (kx = 0), so eps = i^{s}.
    """
    g = ggcd(zx, zy)
    al = exact_div(zx, g)
    al2 = exact_div(zy, g)
    beta = sub(al, al2)
    ca = conj(al)
    s = None
    for k, u in enumerate(UNITS):
        if mul(u, ca) == al2:
            s = k
            break
    assert s is not None, (zx, zy, al, al2)
    X, Y = al
    D = [Y, X - Y, X, X + Y][s]
    gam = [(0, 2), (1, -1), (2, 0), (1, 1)][s]
    assert mul(gam, (D, 0)) == beta, (beta, gam, D)
    return g, beta, s, D
