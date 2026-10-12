"""Shared helpers for round-8 checks (exact Gaussian arithmetic, profiles, ranks).

Gaussian integers are pairs (a, b) meaning a + b i.  All arithmetic is exact
(Python integers).  Angles use mpmath with 60 digits.
"""
from fractions import Fraction
from math import gcd
import mpmath
import sympy

mpmath.mp.dps = 60


def gmul(u, v):
    return (u[0] * v[0] - u[1] * v[1], u[0] * v[1] + u[1] * v[0])


def gconj(u):
    return (u[0], -u[1])


def gnorm(u):
    return u[0] * u[0] + u[1] * u[1]


def gpow(u, e):
    r = (1, 0)
    for _ in range(e):
        r = gmul(r, u)
    return r


def gdivexact(u, v):
    """u / v if v divides u in Z[i], else None."""
    n = gnorm(v)
    w = gmul(u, gconj(v))
    if w[0] % n or w[1] % n:
        return None
    return (w[0] // n, w[1] // n)


def ggcd(u, v):
    while v != (0, 0):
        n = gnorm(v)
        w = gmul(u, gconj(v))
        # rounded quotient
        q = ((2 * w[0] + n) // (2 * n), (2 * w[1] + n) // (2 * n))
        r = (u[0] - gmul(q, v)[0], u[1] - gmul(q, v)[1])
        u, v = v, r
    return u


def ggcd_list(zs):
    g = (0, 0)
    for z in zs:
        g = ggcd(g, z) if g != (0, 0) else z
    return g


def two_squares(p):
    """x, y > 0 with x^2 + y^2 = p for a prime p = 1 mod 4 (or p = 2)."""
    if p == 2:
        return (1, 1)
    # sqrt(-1) mod p
    for c in range(2, p):
        t = pow(c, (p - 1) // 4, p)
        if (t * t) % p == p - 1:
            break
    a, b = p, t
    lim = int(p ** 0.5)
    while b > lim:
        a, b = b, a % b
    x = b
    y2 = p - x * x
    y = sympy.integer_nthroot(y2, 2)[0]
    assert x * x + y * y == p
    return (x, y)


def gauss_prime(p):
    """Fixed Gaussian prime above a split prime p."""
    x, y = two_squares(p)
    return (x, y)


def gval(z, pi):
    """Valuation of z at the Gaussian prime pi."""
    v = 0
    while True:
        w = gdivexact(z, pi)
        if w is None:
            return v
        z = w
        v += 1


def angle(z):
    return mpmath.atan2(z[1], z[0])


def int_rank(rows):
    """Rank over Q of a list of integer vectors."""
    if not rows:
        return 0
    m = [[Fraction(x) for x in r] for r in rows]
    rk = 0
    ncol = len(m[0])
    for c in range(ncol):
        piv = None
        for i in range(rk, len(m)):
            if m[i][c] != 0:
                piv = i
                break
        if piv is None:
            continue
        m[rk], m[piv] = m[piv], m[rk]
        for i in range(len(m)):
            if i != rk and m[i][c] != 0:
                f = m[i][c] / m[rk][c]
                m[i] = [a - f * b for a, b in zip(m[i], m[rk])]
        rk += 1
        if rk == len(m):
            break
    return rk


def naive_height_logs(p, pi):
    """log H(pi/conj(pi)) where H = max(p, 2|Re(pi^2)|) (primitive minimal polynomial)."""
    re2 = pi[0] * pi[0] - pi[1] * pi[1]
    return mpmath.log(max(p, 2 * abs(re2)))


def profile(points):
    """Normalise a list of Gaussian integers of one norm and return
    (N, primes, pis, exponent vectors, g) with exponents at the fixed pi_p."""
    g = ggcd_list(points)
    pts = [gdivexact(z, g) for z in points]
    N = gnorm(pts[0])
    assert all(gnorm(z) == N for z in pts)
    fac = sympy.factorint(N)
    primes = sorted(fac)
    pis = {}
    for p in primes:
        assert p % 4 == 1, (p, N)
        pis[p] = gauss_prime(p)
    vecs = []
    for z in pts:
        vecs.append([gval(z, pis[p]) for p in primes])
    return N, fac, primes, pis, vecs, pts


def slw_exponent(sub_pts, sub_vecs, primes, pis):
    """gamma = log(1/max|Lambda_i|) / sum_{j in U} log H(alpha_j) for the forms
    Lambda_i = Log(z_i/z_0) (i >= 1); U = primes where the exponents vary."""
    th = [angle(z) for z in sub_pts]
    t0 = th[0]
    lam = max(abs(mpmath.atan2(mpmath.sin(t - t0), mpmath.cos(t - t0))) for t in th[1:])
    U = [j for j in range(len(primes)) if len({v[j] for v in sub_vecs}) > 1]
    denom = sum(naive_height_logs(primes[j], pis[primes[j]]) for j in U)
    return float(mpmath.log(1 / lam) / denom), float(lam), U
